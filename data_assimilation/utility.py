import time
from typing import TextIO

from tqdm import tqdm

from data_assimilation import FieldReader


def wait_artss(wait_time: float, artss_data_path: str):
    t_cur = FieldReader.get_t_current(path=artss_data_path)
    pbar = tqdm(total=wait_time)
    pbar.update(t_cur)
    while t_cur <= wait_time:
        time.sleep(1)
        t_new = FieldReader.get_t_current(path=artss_data_path)
        pbar.update(t_new - t_cur)
        t_cur = t_new

    pbar.close()


def write_da_data(file_da: TextIO, parameters: dict):
    """
    write sensor data and ARTSS data to file
    :param file_da: name of file to write to
    :param parameters:
    """
    for key in parameters:
        file_da.write(f';{key}:{parameters[key]}')
    file_da.write('\n')


def log(message: str, file_debug: TextIO):
    print(message)
    file_debug.write(f'{message}\n')


def kelvin_to_celsius(kelvin):
    return kelvin - 273.5
