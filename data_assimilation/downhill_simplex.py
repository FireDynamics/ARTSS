import math
import os
import time
from pprint import pprint
from typing import TextIO, Tuple, Dict

import numpy as np
import pandas
import pandas as pd
import scipy.optimize as op

import TCP_client
import data_assimilation
import fds_utility
from ARTSS import XML, Domain
from data_assimilation import FieldReader, DAFile
from utility import wait_artss, log, write_da_data, kelvin_to_celsius


def nelder_mead_special(client: TCP_client,
                        file_da: TextIO, file_debug: TextIO,
                        sensor_times: pandas.Index, devc_info: dict,
                        cur: dict, delta: dict, keys: list,
                        artss_data_path: str, artss: XML, domain: Domain,
                        fds_data: pandas.DataFrame,
                        heat_source: dict,
                        n_iterations: int,
                        parallel=True,
                        precondition=False):
    t_artss = 0.0
    t_revert = 0.0
    t_sensor = 0.0

    def f(x):
        print(x)
        print(keys)
        new_para = {}
        for i in range(len(keys)):
            new_para[keys[i]] = x[i]
        diff_cur, _ = do_rollback(client=client,
                                  t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                  new_para=new_para, heat_source=heat_source.values(),
                                  sub_file_name='scipy',
                                  artss_data_path=artss_data_path,
                                  fds_data=fds_data,
                                  devc_info=devc_info,
                                  file_da=file_da, file_debug=file_debug)
        log(f'diff: {diff_cur["T"]}', file_debug)
        return diff_cur['T']

    def call_back(xk) -> bool:
        print(f'xk: {xk}')
        file_debug.write(f'xk: {xk}')
        return True

    for count, t_sensor in enumerate(sensor_times):
        pprint(cur)

        wait_artss(t_sensor, artss_data_path)
        t_artss, t_revert = FieldReader.get_time_step_artss(t_sensor,
                                                            artss_data_path,
                                                            dt=artss.get_dt(),
                                                            time_back=6)
        wait_artss(t_artss, artss_data_path)

        start_time = time.time()
        log(f't_sensor: {t_sensor} t_artss: {t_artss} t_revert: {t_revert}', file_debug)
        field_reader = FieldReader(t_artss, path=artss_data_path)
        file_da.write(f'original data;time_sensor:{t_sensor};time_artss:{t_artss}\n')
        write_da_data(file_da=file_da, parameters=cur)
        diff_orig, minima_x = comparison_sensor_simulation_data(devc_info, fds_data, field_reader, t_sensor, file_da)

        file_da.write(f'original: t_artss:{t_artss};t_sensor:{t_sensor};differenceT:{diff_orig["T"]};HRR:{cur["HRR"]};x0:{cur["x0"]}\n')

        log(f'org: {diff_orig["T"]}', file_debug)
        log(f't_revert: {t_revert}', file_debug)

        if diff_orig['T'] < 1e-5:
            log(f'skip, difference: {diff_orig["T"]}', file_debug)
            continue

        if precondition:
            log('precondition for x0', file_debug)
            precondition = False
            safe_keys = keys
            key_precon = ['x0']
            keys = key_precon
            x0 = [cur[x] for x in key_precon]
            initial_simplex = np.array([x0] * (len(key_precon) + 1))
            for i in range(len(key_precon)):
                initial_simplex[i, i] += delta[key_precon[i]]
            interim_start = time.time()
            res = op.minimize(f,
                              x0=np.array(x0),
                              callback=call_back,
                              tol=1e-5,
                              method='Nelder-Mead',
                              options={
                                  'maxiter': 5,
                                  'disp': True,
                                  'initial_simplex': initial_simplex
                              })
            interim_end = time.time()
            log(f't_artss: {t_artss}; interim time: {interim_end - interim_start}', file_debug)
            log(f'org: {diff_orig}', file_debug)
            log(f'res: {res}', file_debug)
            for i in range(len(keys)):
                cur[keys[i]] = res.x[i]
            log(f'res_x (new parameter): {list(res.x)}\n', file_debug)
            log(f'res_sim (final simplex): {res.final_simplex}\n', file_debug)
            log(f'res_fun: {res.fun}\n', file_debug)
            log(f'res_n (number of iterations): {res.nit}\n', file_debug)
            log(f'res_nfev: {res.nfev}\n', file_debug)
            log(f'res_suc (exit successfully): {res.success}\n', file_debug)
            log(f'res_msg (message if exit was unsuccessful): {res.message}\n', file_debug)
            log(f'final rollback with {cur}', file_debug)
            diff_cur, _ = do_rollback(client=client,
                                      t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                      new_para=cur, heat_source=heat_source.values(),
                                      sub_file_name='final',
                                      artss_data_path=artss_data_path,
                                      fds_data=fds_data,
                                      devc_info=devc_info,
                                      file_da=file_da, file_debug=file_debug)
            file_da.write(f'final: t_artss:{t_artss};t_sensor:{t_sensor};differenceT:{diff_cur["T"]};HRR:{cur["HRR"]};x0:{cur["x0"]}\n')
            file_da.flush()
            file_debug.flush()
            end_time = time.time()
            log(f'time: {end_time - start_time}', file_debug)
            keys = safe_keys

        x0 = [cur[x] for x in keys]
        delta['HRR'] = cur['HRR'] * 0.05,
        initial_simplex = np.array([x0] * (len(keys) + 1))
        for i in range(len(keys)):
            initial_simplex[i, i] += delta[keys[i]]
        interim_start = time.time()
        res = op.minimize(f,
                          x0=np.array(x0),
                          callback=call_back,
                          tol=1e-5,
                          method='Nelder-Mead',
                          options={
                              'maxiter': n_iterations,
                              'disp': True,
                              'initial_simplex': initial_simplex
                          })
        interim_end = time.time()
        log(f't_artss: {t_artss}; interim time: {interim_end - interim_start}', file_debug)
        log(f'org: {diff_orig}', file_debug)
        log(f'res: {res}', file_debug)
        for i in range(len(keys)):
            cur[keys[i]] = res.x[i]
        log(f'res_x (new parameter): {list(res.x)}\n', file_debug)
        log(f'res_sim (final simplex): {res.final_simplex}\n', file_debug)
        log(f'res_fun: {res.fun}\n', file_debug)
        log(f'res_n (number of iterations): {res.nit}\n', file_debug)
        log(f'res_nfev: {res.nfev}\n', file_debug)
        log(f'res_suc (exit successfully): {res.success}\n', file_debug)
        log(f'res_msg (message if exit was unsuccessful): {res.message}\n', file_debug)
        log(f'final rollback with {cur}', file_debug)
        diff_cur, _ = do_rollback(client=client,
                                  t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                  new_para=cur, heat_source=heat_source.values(),
                                  sub_file_name='final',
                                  artss_data_path=artss_data_path,
                                  fds_data=fds_data,
                                  devc_info=devc_info,
                                  file_da=file_da, file_debug=file_debug)
        file_da.write(f'final: t_artss:{t_artss};t_sensor:{t_sensor};differenceT:{diff_cur["T"]};HRR:{cur["HRR"]};x0:{cur["x0"]}\n')
        file_da.flush()
        file_debug.flush()
        end_time = time.time()
        log(f'time: {end_time - start_time}', file_debug)


def opt_scipy(client: TCP_client,
              file_da: TextIO, file_debug: TextIO,
              sensor_times: pandas.Index, devc_info: dict,
              cur: dict, delta: dict, keys: list,
              artss_data_path: str, artss: XML, domain: Domain,
              fds_data: pandas.DataFrame,
              heat_source: dict,
              n_iterations: int,
              parallel=True):
    t_artss = 0.0
    t_revert = 0.0
    t_sensor = 0.0

    bounds = op.Bounds(
        lb=np.asarray([0, domain.domain_param['X1']]),
        ub=np.asarray([np.inf, domain.domain_param['X2']])
    )

    def f(x):
        print(x)
        print(keys)
        new_para = {}
        for i in range(len(keys)):
            new_para[keys[i]] = x[i]
        diff_cur, _ = do_rollback(client=client,
                                  t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                  new_para=new_para, heat_source=heat_source.values(),
                                  sub_file_name='scipy',
                                  artss_data_path=artss_data_path,
                                  fds_data=fds_data,
                                  devc_info=devc_info,
                                  file_da=file_da, file_debug=file_debug)
        log(f'diff: {diff_cur["T"]}', file_debug)
        return diff_cur['T']

    def call_back(xk) -> bool:
        print(f'xk: {xk}')
        file_debug.write(f'xk: {xk}')
        return True

    iterations = np.ones(len(sensor_times)) * n_iterations
    # iterations[0] = 15
    for count, t_sensor in enumerate(sensor_times):
        pprint(cur)

        wait_artss(t_sensor, artss_data_path)
        t_artss, t_revert = FieldReader.get_time_step_artss(t_sensor,
                                                            artss_data_path,
                                                            dt=artss.get_dt(),
                                                            time_back=6)
        wait_artss(t_artss, artss_data_path)
        if count == 0:
            t_revert = 0.2
        start_time = time.time()
        log(f't_sensor: {t_sensor} t_artss: {t_artss} t_revert: {t_revert}', file_debug)
        field_reader = FieldReader(t_artss, path=artss_data_path)
        file_da.write(f'original data;time_sensor:{t_sensor};time_artss:{t_artss}\n')
        write_da_data(file_da=file_da, parameters=cur)
        diff_orig, minima_x = comparison_sensor_simulation_data(devc_info, fds_data, field_reader, t_sensor, file_da)

        file_da.write(f'original: t_artss:{t_artss};t_sensor:{t_sensor};differenceT:{diff_orig["T"]};HRR:{cur["HRR"]};x0:{cur["x0"]};z0:{cur["z0"]}\n')

        log(f'org: {diff_orig["T"]}', file_debug)
        log(f't_revert: {t_revert}', file_debug)

        if diff_orig['T'] < 1e-5:
            log(f'skip, difference: {diff_orig["T"]}', file_debug)
            continue
        x0 = [cur[x] for x in keys]
        delta['HRR'] = cur['HRR'] * 0.05,
        initial_simplex = np.array([x0] * (len(keys) + 1))
        for i in range(len(keys)):
            initial_simplex[i, i] += delta[keys[i]]
        interim_start = time.time()
        res = op.minimize(f,
                          x0=np.array(x0),
                          callback=call_back,
                          tol=1e-5,
                          method='Nelder-Mead',
                          bounds=bounds,
                          options={
                              'maxiter': iterations[count],
                              'disp': True,
                              'initial_simplex': initial_simplex
                          })
        interim_end = time.time()
        log(f't_artss: {t_artss}; interim time: {interim_end - interim_start}', file_debug)
        log(f'org: {diff_orig}', file_debug)
        log(f'res: {res}', file_debug)
        for i in range(len(keys)):
            cur[keys[i]] = res.x[i]
        log(f'res_x (new parameter): {list(res.x)}\n', file_debug)
        log(f'res_sim (final simplex): {res.final_simplex}\n', file_debug)
        log(f'res_fun: {res.fun}\n', file_debug)
        log(f'res_n (number of iterations): {res.nit}\n', file_debug)
        log(f'res_nfev: {res.nfev}\n', file_debug)
        log(f'res_suc (exit successfully): {res.success}\n', file_debug)
        log(f'res_msg (message if exit was unsuccessful): {res.message}\n', file_debug)
        log(f'final rollback with {cur}', file_debug)
        diff_cur, _ = do_rollback(client=client,
                                  t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                  new_para=cur, heat_source=heat_source.values(),
                                  sub_file_name='final',
                                  artss_data_path=artss_data_path,
                                  fds_data=fds_data,
                                  devc_info=devc_info,
                                  file_da=file_da, file_debug=file_debug)
        file_da.write(f'final: t_artss:{t_artss};t_sensor:{t_sensor};differenceT:{diff_cur["T"]};HRR:{cur["HRR"]};x0:{cur["x0"]};z0:{cur["z0"]}\n')
        file_da.flush()
        file_debug.flush()
        end_time = time.time()
        log(f'time: {end_time - start_time}', file_debug)


def do_rollback(client: TCP_client,
                t_sensor, t_artss, t_revert,
                new_para: dict, heat_source: [dict, dict, dict],
                sub_file_name: str,
                artss_data_path: str,
                fds_data,
                devc_info: dict,
                file_da: TextIO, file_debug: TextIO) -> [dict, float]:
    config_file_name, _ = write_changes_xml(
        new_para,
        heat_source,
        f'{t_artss}_{sub_file_name}',
        file_da=file_da,
        path=artss_data_path
    )
    print(f'make adjustment with {new_para}')
    file_debug.write(f'make adjustment with: {new_para}\n')
    write_da_data(file_da=file_da, parameters=new_para)
    client.send_message(data_assimilation.create_message(t_revert, config_file_name))
    wait_artss(t_artss, artss_data_path)

    field_reader = FieldReader(t_artss, path=artss_data_path)
    diff_cur, min_pos_x = comparison_sensor_simulation_data(
        devc_info,
        fds_data,
        field_reader,
        t_sensor,
        file_da
    )
    return diff_cur, min_pos_x


def comparison_sensor_simulation_data(devc_info: dict, sensor_data: pd.DataFrame, field_reader: FieldReader,
                                      t_sensor: float, file_da: TextIO) -> Tuple[Dict[str, float], float]:
    diff: dict = {'T': [], 'C': []}
    fields_sim = field_reader.read_field_data()
    min_pos_x: float = 0
    min_sensor_val = math.inf
    for key in devc_info:
        # sensor data
        type_sensor: str = devc_info[key]['type']
        index_sensor: int = devc_info[key]['index']
        value_sensor: float = sensor_data[key][t_sensor]

        value_sim = fields_sim[type_sensor][index_sensor]
        difference = kelvin_to_celsius(value_sim) - value_sensor
        if abs(difference) > 2:
            diff[type_sensor].append(difference)
        file_da.write(
            f'sensor:{key};time_sensor:{t_sensor};time_artss:{field_reader.t};sensor_val:{value_sensor};artss_val:{value_sim};diff:{difference};considered:{abs(difference) > 2};position:{devc_info[key]["XYZ"]}\n')
        file_da.flush()

        if difference < min_sensor_val:
            min_sensor_val = difference
            min_pos_x = devc_info[key]['XYZ'][0]
    print(f'diff: {diff["T"]}')
    result = {}
    for key in diff:
        if len(diff[key]) == 0:
            result[key] = 0
            file_da.write(f'differences:{key}:{result[key]};no_of_sensors:{len(diff[key])}\n')
            continue
        result[key] = np.sqrt(sum(np.array(diff[key]) ** 2))
        file_da.write(f'differences:{key}:{result[key]};no_of_sensors:{len(diff[key])}\n')
    return result, min_pos_x


def write_changes_xml(change: dict, source: list, file_name: str, file_da: TextIO, path='.') -> [str, str]:
    source_type, temperature_source, random = data_assimilation.change_heat_source(*source,
                                                                                   changes={'source_type': {},
                                                                                            'temperature_source': change,
                                                                                            'random': {}})
    write_da_data(file_da=file_da, parameters=temperature_source)
    da = DAFile()
    da.create_temperature_source_changes(
        source_type=source_type,
        temperature_source=temperature_source,
        random=random)

    config_file_name = f'config_{file_name}.xml'
    config_file_path = os.path.join(path, config_file_name)
    da.write_xml(config_file_path)
    return config_file_name, config_file_path


def start(fds_data_path: str, fds_input_file_name: str, artss_data_path: str, parallel=True, port=7777):
    artss_data_path = os.path.abspath(artss_data_path)
    client = TCP_client.TCPClient()
    client.set_server_address('localhost', port)
    # client.set_server_address('172.17.0.2', 7777)
    client.connect()

    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()
    domain = Domain(domain_param=xml.domain, obstacles=xml.obstacles,
                    enable_computational_domain=xml.computational_domain)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = fds_utility.read_fds_data(fds_data_path,
                                                                                        fds_input_file_name, domain)
    print('devices:')
    pprint(devc_info_temperature)

    sensor_times = fds_data.index[3:]
    print('sensor times:', fds_data.index)

    heat_source = xml.get_temperature_source()

    file_da = open(os.path.join(artss_data_path, 'da_details.csv'), 'w')
    file_debug = open(os.path.join(artss_data_path, 'da_debug_details.dat'), 'w')
    pprint(devc_info_temperature, stream=file_debug)

    delta = {
        'HRR': float(heat_source['temperature_source']['HRR']) * 0.05,
        'x0': 1.5,
        'y0': domain.domain_param['dy'] * 1,
        'z0': 0.5
    }

    cur = {
        'HRR': float(heat_source['temperature_source']['HRR']),
        'x0': float(heat_source['temperature_source']['x0']),
        'y0': float(heat_source['temperature_source']['y0']),
        'z0': float(heat_source['temperature_source']['z0'])
    }

    keys = ['HRR', 'x0', 'z0']

    log(f'cur: {cur}', file_debug)
    log(f'delta: {delta}', file_debug)
    log(f'keys: {keys}', file_debug)
    #nelder_mead_special(client=client, file_da=file_da, file_debug=file_debug,
    #                    sensor_times=sensor_times,
    #                    devc_info=devc_info_temperature, fds_data=fds_data,
    #                    artss_data_path=artss_data_path,
    #                    domain=domain, heat_source=heat_source,
    #                    cur=cur, delta=delta, keys=keys, n_iterations=5, artss=xml,
    #                    parallel=parallel,
    #                    precondition=True)
    opt_scipy(client=client, file_da=file_da, file_debug=file_debug,
              sensor_times=sensor_times,
              devc_info=devc_info_temperature, fds_data=fds_data,
              artss_data_path=artss_data_path,
              domain=domain, heat_source=heat_source,
              cur=cur, delta=delta, keys=keys, n_iterations=5, artss=xml,
              parallel=parallel)
