import cv2
import numpy as np
from PIL import Image

def segmentation(img, black_lavel):
    """
    Накладывает маску на полученное изображения, используя его контрастность. Готовит файл
    к отправке в базу данных, конвертируя jpg/png в BLOB

    Параметры
    ----------
    img : Путь до изображения : str 

    black_lavel : уровень черного, по которому будет производиться сегментация. 
                  Число от 0 до 255, где 0 - полностью черное, 255 - полностью белое : int

    Результат
    -------
    width : ширина выделенного объекта в пикселях : numpy.int32
    
    blob_data : исходное изображение с наложенной маской сегментации, 
                переведенное в формат BLOB : bytes
    """
    # Чтение изображения
    image = cv2.imread(img)

    # Преобразование в оттенки серого
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    # Бинаризация
    _, thresh = cv2.threshold(gray, black_lavel, 255, cv2.THRESH_BINARY_INV)

    # Поиск контуров
    contours, _ = cv2.findContours(
        thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

    # Поиск самого темного объекта
    darkest_object = None
    min_avg_intensity = 256
    for cnt in contours:
        mask = np.zeros(gray.shape[:2], dtype=np.uint8)
        cv2.drawContours(mask, [cnt], -1, 255, -1)
        mean_intensity = cv2.mean(gray, mask=mask)[0]
        if mean_intensity < min_avg_intensity:
            darkest_object = cnt
            min_avg_intensity = mean_intensity

    # Нахождение крайних точек самого темного объекта
    if darkest_object is not None:
        leftmost = tuple(darkest_object[darkest_object[:, :, 0].argmin()][0])
        rightmost = tuple(darkest_object[darkest_object[:, :, 0].argmax()][0])
        topmost = tuple(darkest_object[darkest_object[:, :, 1].argmin()][0])
        bottommost = tuple(darkest_object[darkest_object[:, :, 1].argmax()][0])

        width = rightmost[0] - leftmost[0]
        height = bottommost[1] - topmost[1]

        print(f'Длина объекта: {width} пикселей')
        print(f'Высота объекта: {height} пикселей')
    else:
        print("Самый темный объект не найден.")

    # Рисуем контур самого темного объекта
    if darkest_object is not None:
        cv2.drawContours(image, [darkest_object], -1, (0, 255, 0), 2)

        # Добавляем текст с размерами
        cv2.putText(image, f'Width: {width}px', (20, 30),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 0), 2)
        cv2.putText(image, f'Height: {height}px', (20, 60),
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

    return (width, image_to_blob('test_pic_2.jpg'))
