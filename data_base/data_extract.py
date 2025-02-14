import sqlite3

def text_extract():
    # Подключение к базе данных (или создание новой)
    conn = sqlite3.connect('ILT_data_base.db')
    cursor = conn.cursor()

    # SQL-запрос для выборки данных
    select_query = "SELECT * FROM laser_modes"  # Напишите ваш запрос на языке SQLlite

    # Выполнение запроса
    cursor.execute(select_query)

    # Получение всех результатов
    rows = cursor.fetchall()

    # Вывод результатов
    for row in rows:
        print(row)

    # Закрытие соединения
    conn.close()


def jpeg_extract():
    # Подключение к базе данных (или создание новой)
    conn = sqlite3.connect('ILT_data_base.db')
    cursor = conn.cursor()

    # SQL-запрос для выборки микрофотографий
    select_query = "SELECT micro_photo FROM microscope_results;"   # Напишите ваш запрос на языке SQLlite
    cursor.execute(select_query)
    rows = cursor.fetchall()
    for row in rows:
        blob_data = row[0]

    # Закрытие соединения
    conn.close()

    #Преобразование изображения из blob в файл jpg
    with open('изображение_из_базы.jpg', 'wb') as file:
        file.write(blob_data)
     