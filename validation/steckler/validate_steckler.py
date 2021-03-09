import numpy as np
import steckler
import ARTSS
import plot_helper as ph
import util as util

if __name__ == '__main__':
    steckler16 = steckler.Steckler()
    steckler16.read_test_case("test_16")
    artss = ARTSS.ARTSS()
    i = 2000
    artss.add_data_file(f'steckler_num_{i:06d}.csv')
    artss.add_xml('steckler.xml')
    artss.print_info()

    # steckler
    steckler_door_vel = np.flip(util.create_2D_array_from_dict(steckler16.opening_vel_dict), axis=0)
    ph.plot_colormap_vel(steckler_door_vel, title='Steckler (door)')

    steckler_door_temp = np.flip(util.create_2D_array_from_dict(steckler16.opening_temp_dict), axis=0)
    ph.plot_colormap_temp(steckler_door_temp, title='Steckler (door)')

    # artss
    offset_x = artss.domain.get_obstacle('left wall')[0]['pox2']
    offset_z = artss.domain.get_obstacle('front wall')[0]['poz2']
    opening_x, opening_y, opening_z = artss.get_indices_from_artss(steckler16.coordinates_dict['sensors opening y'], steckler16.coordinates_dict['sensors opening z'], steckler16.coordinates_dict['sensors opening x'], offset_x=offset_x, offset_z=offset_z)
    room_x, room_y, room_z = artss.get_indices_from_artss(steckler16.coordinates_dict['sensors room y'], steckler16.coordinates_dict['sensors room z'], steckler16.coordinates_dict['sensors room x'], offset_x=offset_x, offset_z=offset_z)

    artss_vel_x = artss.filter_data_on_indices(opening_x, opening_y, opening_z, data_type='vel_x')
    ph.plot_colormap_vel(artss_vel_x, title='ARTSS (door)')

    artss_temperature_kelvin = artss.filter_data_on_indices(opening_x, opening_y, opening_z, data_type='temperature')
    artss_temperature_celsius = ARTSS.to_celsius(artss_temperature_kelvin)
    ph.plot_colormap_temp(artss_temperature_celsius, title='ARTSS (door)')

    # artss whole door
    interval_x = [opening_x[0], opening_x[-1]]
    interval_y = [opening_y[0], opening_y[-1]]
    interval_z = [opening_z[0], opening_z[-1]]
    artss_temperature_kelvin_door = artss.clip_data(interval_x, interval_y, interval_z, data_type='temperature')
    artss_temperature_celsius_door = ARTSS.to_celsius(artss_temperature_kelvin_door)
    ph.plot_colormap_temp(artss_temperature_celsius_door, title='ARTSS (door)')
