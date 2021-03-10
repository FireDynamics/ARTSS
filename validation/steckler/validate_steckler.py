import numpy as np
import steckler
import ARTSS
import plot_helper as ph
import util as util

if __name__ == '__main__':
    # steckler
    steckler16 = steckler.Steckler()
    steckler16.read_test_case("test_16")  # specify test case

    # ARTSS
    artss = ARTSS.ARTSS()
    i = 200  # number of csv file
    xml_name='Test_NavierStokesTempTurb_Steckler'  # name of xml file (without file ending)
    artss.add_xml(xml_name+'.xml')  # add xml

    artss.add_data_file(f'{xml_name}_num_{i:06d}.csv')  # add csv file
    artss.print_info()  # print info of ARTSS simluation
    
    # calculation of indices
    offset_x = artss.domain.get_obstacle('left wall')[0]['pox2']
    offset_z = artss.domain.get_obstacle('front wall')[0]['poz2']
    opening_x, opening_y, opening_z = artss.get_indices_from_artss(steckler16.coordinates_dict['sensors opening y'], steckler16.coordinates_dict['sensors opening z'], steckler16.coordinates_dict['sensors opening x'], offset_x=offset_x, offset_z=offset_z)
    room_x, room_y, room_z = artss.get_indices_from_artss(steckler16.coordinates_dict['sensors room y'], steckler16.coordinates_dict['sensors room z'], steckler16.coordinates_dict['sensors room x'], offset_x=offset_x, offset_z=offset_z)

    artss_vel_x = artss.filter_data_on_indices(opening_x, opening_y, opening_z, data_type='vel_x')  # get velocity data of door opening
    artss_temperature_kelvin = artss.filter_data_on_indices(opening_x, opening_y, opening_z, data_type='temperature')  # temperature data in kelvin
    artss_temperature_celsius = ARTSS.to_celsius(artss_temperature_kelvin)  # get temperature data of door opening in celsius
    ph.plot_colormap_vel(artss_vel_x, title='ARTSS (door)',savefig=True,filetitle='artss_vel')
    ph.plot_colormap_temp(artss_temperature_celsius, title='ARTSS (door)',savefig=True,filetitle='artss_temp')

    # artss whole door
    interval_x = [opening_x[0], opening_x[-1]]
    interval_y = [opening_y[0], opening_y[-1]]
    interval_z = [opening_z[0], opening_z[-1]]
    artss_temperature_kelvin_door = artss.clip_data(interval_x, interval_y, interval_z, data_type='temperature')
    artss_temperature_celsius_door = ARTSS.to_celsius(artss_temperature_kelvin_door)
    ph.plot_colormap_temp(artss_temperature_celsius_door, title='ARTSS (door)',savefig=True,filetitle='artss_temp_original_resolution')

    # steckler
    steckler_door_vel = np.flip(util.create_2D_array_from_dict(steckler16.opening_vel_dict), axis=0)  # get velocity data of door opening
    steckler_door_temp = np.flip(util.create_2D_array_from_dict(steckler16.opening_temp_dict), axis=0)  # get temperature data of door opening
    ph.plot_colormap_vel(steckler_door_vel, title='Steckler (door)',savefig=True,filetitle='steckler_vel')  # plot velocity data
    ph.plot_colormap_temp(steckler_door_temp, title='Steckler (door)',savefig=True,filetitle='steckler_temp')  # plot temperature data

    ph.plot_colormap_vel(steckler_door_vel-artss_temperature_celsius,title='difference',savefig=True,filetitle='difference')
