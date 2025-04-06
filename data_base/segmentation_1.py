import cv2
import numpy as np
from PIL import Image
import os

class Segmentation():
    def segmentation(img, black_lavel, type='dark'):

        """
        Накладывает маску на полученное изображения, используя его контрастность. Готовит файл
        к отправке в базу данных, конвертируя jpg/png в BLOB
    
        Параметры
        ----------
        img : Путь до изображения : str 
    
        black_lavel : уровень черного, по которому будет производиться сегментация. 
                      Число от 0 до 255, где 0 - полностью черное, 255 - полностью белое : int
        
        type: тип сегментации, dark - выделение темных объектов. != dark - выделение светлых: str 
            
        Результат
        -------
        width : ширина выделенного объекта в пикселях : numpy.int32
        
        blob_data : исходное изображение с наложенной маской сегментации, 
                    переведенное в формат BLOB : bytes
        """
            
        # Чтение изображения
        try:
            image = cv2.imread(img)
        except:
            image = img

        sup = Image.fromarray(image)

        # Преобразование в оттенки серого
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

        # Бинаризация
        if type == 'dark':
            _, thresh = cv2.threshold(
                gray, black_lavel, 255, cv2.THRESH_BINARY_INV)
        else:
            _, thresh = cv2.threshold(
                gray, black_lavel, 255, cv2.THRESH_BINARY)

        sup = Image.fromarray(thresh)
        sup.save('sup_image.jpg')

        # Поиск контуров
        contours, _ = cv2.findContours(
            thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

        # Поиск самого темного объекта
        def darkest_object_func(contours):
            darkest_object = None
            min_avg_intensity = 256
            for cnt in contours:
                mask = np.zeros(gray.shape[:2], dtype=np.uint8)
                cv2.drawContours(mask, [cnt], -1, 255, -1)
                mean_intensity = cv2.mean(gray, mask=mask)[0]
                if mean_intensity < min_avg_intensity:
                    darkest_object = cnt
                    min_avg_intensity = mean_intensity
            return darkest_object

        # Поиск самого светлого объекта

        def brightest_object_func(contours):
            brightest_object = None
            min_avg_intensity = -1
            for cnt in contours:
                mask = np.zeros(gray.shape[:2], dtype=np.uint8)
                cv2.drawContours(mask, [cnt], -1, 255, -1)
                mean_intensity = cv2.mean(gray, mask=mask)[0]
                if mean_intensity > min_avg_intensity:
                    brightest_object = cnt
                    min_avg_intensity = mean_intensity
            return brightest_object

        # Нахождение крайних точек самого темного объекта
        if type == 'dark':
            # print('USE DARK TYPE')
            segment_object = darkest_object_func(contours)
        else:
            # print('USE BRIGHT TYPE')
            segment_object = brightest_object_func(contours)

        if segment_object is not None:
            leftmost = tuple(
                segment_object[segment_object[:, :, 0].argmin()][0])
            rightmost = tuple(
                segment_object[segment_object[:, :, 0].argmax()][0])
            topmost = tuple(
                segment_object[segment_object[:, :, 1].argmin()][0])
            bottommost = tuple(
                segment_object[segment_object[:, :, 1].argmax()][0])

            lenght = np.round(0.53836990 * (rightmost[0] - leftmost[0]), 5)
            width = np.round(0.53836990 * (bottommost[1] - topmost[1]), 5)

        else:
            print(f"Самый {type} объект не найден.")
            return (0, 0)

        # Рисуем контур самого темного объекта
        if segment_object is not None:
            cv2.drawContours(image, [segment_object], -1, (0, 255, 0), 2)

            # Добавляем текст с размерами
            cv2.putText(image, f'Width: {width} micrometers', (20, 30),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 0), 2)
            cv2.putText(image, f'Lenght: {lenght} micrometers', (20, 60),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 0), 2)

        # Преобразуем изображение из BGR в RGB
        image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

        # Преобразуем в объект PIL
        image_pil = Image.fromarray(image_rgb)

        # Сохраняем и показываем изображение
        image_pil.save('test_pic_2.jpg')

        def image_to_blob(image_path):
            with open(image_path, 'rb') as file:
                blob_data = file.read()
            return blob_data

        return (width, image_pil)

    def calculate_percentile_brightness(image_path, procentage: int, type='dark'):
        # Открываем изображение и конвертируем в черно-белое
        image = Image.open(image_path).convert('L')
        pixels = np.array(image)
        if type != 'dark':
            procentage = abs(100-procentage)
        percentile = np.percentile(pixels, procentage)

        return percentile

    def crop_center_square(image_path):
        # Открываем изображение
        img = Image.open(image_path)

        # Получаем размеры изображения
        width, height = img.size
        size = height // 2
        # Вычисляем координаты для обрезки центрального квадрата
        left = width // 2
        top = (height - size) // 2
        right = width
        bottom = (height + size) // 2

        # Обрезаем изображение
        cropped_img = img.crop((left, top, right, bottom))
        # image_pil = Image.fromarray(cropped_img)
        return np.array(cropped_img)    
