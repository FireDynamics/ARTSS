import math
import os
import time
from typing import Dict, TextIO

import fdsreader
import multiprocessing
import numpy as np
import pandas
import pandas as pd
import matplotlib.pyplot as plt

from copy import copy
from pprint import pprint

from tqdm import tqdm

import ARTSS
import TCP_client
import data_assimilation
from ARTSS import Domain, XML, DAFile
from data_assimilation import FieldReader
from main import create_message


def write_da_data(file_da: TextIO, parameters: dict):
    for key in parameters:
        file_da.write(f';{key}:{parameters[key]}')
    file_da.write('\n')


def map_minima(client, artss_data_path, cur, delta, sensor_times, devc_info_thermocouple, devc_info_temperature,
               fds_data, source_type, temperature_source, random, file_da, cwd, artss):
    def f(t_end, param):
        t_artss, t_revert = get_time_step_artss(t_end, artss_data_path)
        config_file_name, _ = write_changes_xml(
            param,
            [source_type, temperature_source, random],
            f'{t_artss}_f',
            file_da=file_da,
            path=cwd
        )
        client.send_message(create_message(t_revert, config_file_name))
        wait_artss(t_end, artss_data_path, artss)

        field_reader = FieldReader(t_artss, path=artss_data_path)
        ret, _ = comparison_sensor_simulation_data(
            devc_info_thermocouple,
            fds_data,
            field_reader,
            t_artss,
            t_end
        )

        return ret['T']

    print('org: {f(t_artss, t_revert, dict())}')
    t = sensor_times[2]
    pprint(cur)
    wait_artss(t, artss_data_path, artss)

    t_artss, t_revert = get_time_step_artss(t, artss_data_path)
    pprint((t, t_artss, t_revert))
    field_reader = FieldReader(t_artss, path=artss_data_path)
    diff_orig, _ = comparison_sensor_simulation_data(devc_info_temperature, fds_data, field_reader, t_artss, t)
    print(f'org: {diff_orig["T"]}')
    print(t_revert)

    with open('map_2.csv', 'w') as out:
        for x in range(20, 85, 10):
            for y in [delta['y0'] * t for t in range(8)]:
                ret = f(t, {'x0': x, 'y0': y})
                print(f'map: {x},{y},{ret}')
                out.write(f'{x},{y},{ret}\n')

    return


def start_new_instance(output_file: str, directory: str, artss_exe_path: str,
                       artss_exe: str = 'artss_data_assimilation_serial'):
    cwd = os.getcwd()
    os.chdir(directory)
    exe_command = f'mpirun --np 2 {os.path.join(artss_exe_path, artss_exe)} {output_file}'
    print(exe_command)
    os.system(exe_command)
    os.chdir(cwd)


def create_start_xml(change: dict, input_file: str, output_file: str, t_revert: float, file: str, t_end: float,
                     directory: str, artss_path: str, file_debug: TextIO, dir_name: str = '') -> [str, str]:
    f = os.path.join('..', file)
    f = f.replace("/", "\/")

    command = ''
    name = ''
    for key in change:
        command += f' --{key} {change[key]}'
        name += f'{key}_'
    name = name[:-1]

    if not dir_name == '':
        name = dir_name
    output_dir = os.path.join(directory, '.vis', name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    out_file = os.path.join(output_dir, output_file)
    log(f'out_file: {out_file}', file_debug)
    log(f'out_dir: {output_dir}', file_debug)
    command = f'{os.path.join(artss_path, "tools", "change_xml.sh")} -i {input_file} -o {out_file} --da {t_revert} "{f}" 7779 --tend {t_end} --loglevel off' + command
    log(f'{command}', file_debug)
    os.system(command)
    return out_file, output_dir


def calc_diff(t_artss: float, t_sensor: float,
              data_path: str, fds_data: pd.DataFrame, devc_info: dict,
              file_da: TextIO) -> [Dict[str, float], float]:
    field_reader = FieldReader(t_artss, path=data_path)
    file_da.write(f'time_sensor:{t_sensor};time_artss:{t_artss}\n')
    diff_cur, min_pos_x = comparison_sensor_simulation_data(
        devc_info,
        fds_data,
        field_reader,
        t_sensor,
        file_da
    )
    return diff_cur, min_pos_x


def do_rollback(client: TCP_client,
                t_sensor, t_artss, t_revert,
                new_para: dict, heat_source: [dict, dict, dict],
                sub_file_name: str,
                artss_data_path: str, artss: XML,
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
    client.send_message(create_message(t_revert, config_file_name))
    wait_artss(t_sensor, artss_data_path, artss)

    field_reader = FieldReader(t_artss, path=artss_data_path)
    diff_cur, min_pos_x = comparison_sensor_simulation_data(
        devc_info,
        fds_data,
        field_reader,
        t_sensor,
        file_da
    )
    return diff_cur, min_pos_x


def log(message: str, file_debug: TextIO):
    print(message)
    file_debug.write(f'{message}\n')


def continuous_gradient_parallel(client: TCP_client,
                                 file_da: TextIO, file_debug: TextIO,
                                 sensor_times: pandas.Index, devc_info: dict,
                                 cur: dict, delta: dict, keys: list,
                                 artss_path: str, artss_data_path: str,
                                 artss: XML, domain: Domain,
                                 fds_data: pandas.DataFrame,
                                 heat_source: dict,
                                 n_iterations: int,
                                 precondition: bool = False):
    no_cpus = multiprocessing.cpu_count()
    artss_exe_path = os.path.join(artss_path, 'build', 'bin')
    xml_input_file = os.path.join(artss_data_path, FieldReader.get_xml_file_name(path=artss_data_path))
    n = max(1, n_iterations)
    for t_sensor in sensor_times[3:]:
        pprint(cur)
        wait_artss(t_sensor, artss_data_path, artss)

        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=artss.get_dt(), time_back=6)

        log(f't_sensor: {t_sensor} t_artss: {t_artss} t_revert: {t_revert}', file_debug)
        field_reader = FieldReader(t_artss, path=artss_data_path)
        file_da.write(f'original data;time_sensor:{t_sensor};time_artss:{t_artss}\n')
        write_da_data(file_da=file_da, parameters=cur)
        diff_orig, minima_x = comparison_sensor_simulation_data(devc_info, fds_data, field_reader, t_sensor, file_da)

        log(f'org: {diff_orig["T"]}', file_debug)
        log(f't revert: {t_revert}', file_debug)

        if sum(diff_orig.values()) < 1e-5:
            continue

        if precondition:
            precondition = False
            cur['x0'] = minima_x
            heat_source['temperature_source']['x0'] = cur['x0']
            file_da.write(f'preconditioning;')
            diff_orig, _ = do_rollback(client=client,
                                       t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                       new_para={'x0': cur['x0']}, heat_source=heat_source.values(),
                                       sub_file_name='precondition',
                                       artss_data_path=artss_data_path, artss=artss,
                                       fds_data=fds_data,
                                       devc_info=devc_info,
                                       file_da=file_da, file_debug=file_debug)

        nabla = {}

        jobs = {}
        directories = {}
        counter = 1
        for p in keys:
            output_file, output_dir = create_start_xml(change={p: cur[p] + delta[p]},
                                                       input_file=xml_input_file, output_file=f'config_{p}.xml',
                                                       t_revert=t_revert, t_end=t_artss,
                                                       directory=artss_data_path,
                                                       file=f'{t_revert:.5e}',
                                                       artss_path=artss_path,
                                                       file_debug=file_debug)
            directories[p] = output_dir
            job = multiprocessing.Process(target=start_new_instance, args=(
                output_file, output_dir, artss_exe_path, 'artss_data_assimilation_serial'))
            job.start()
            counter += 1
            jobs[p] = job
            if counter >= no_cpus:
                # wait if you want to start more multi processes than cpus available
                j = jobs.pop(list(jobs.keys())[0])
                j.join()
        print('calc nabla parallel')
        for p in keys:
            if p in jobs.keys():
                job = jobs[p]
                print(f'wait for {job.name} [{job}]')
                job.join()
            file_da.write(f'calc_nabla;{p}:{cur[p] + delta[p]}')
            diff_cur, _ = calc_diff(t_artss=t_artss, t_sensor=t_sensor,
                                    data_path=directories[p], fds_data=fds_data, devc_info=devc_info,
                                    file_da=file_da)
            print(p, diff_cur['T'])
            nabla[p] = (diff_cur['T'] - diff_orig['T']) / delta[p]
        log(f'nabla: {nabla}', file_debug)
        file_da.write(f'nabla without restrictions;')
        write_da_data(file_da, nabla)

        log(f'nabla without restrictions: {nabla}', file_debug)
        restrictions = {}
        if 'HRR' in keys:
            restrictions['HRR'] = 1

        if 'x0' in keys:
            if nabla['x0'] < 0:
                restrictions['x0'] = (domain.domain_param['X2'] - cur['x0']) / -nabla['x0']
            elif nabla['x0'] > 0:
                restrictions['x0'] = (cur['x0'] - domain.domain_param['X1']) / nabla['x0']
            else:
                restrictions['x0'] = 0
        if 'y0' in keys:
            if nabla['y0'] < 0:
                restrictions['y0'] = (domain.domain_param['Y2'] - cur['y0']) / -nabla['y0']
            elif nabla['y0'] > 0:
                restrictions['y0'] = (cur['y0'] - domain.domain_param['Y1']) / nabla['y0']
            else:
                restrictions['y0'] = 0
        if 'z0' in keys:
            if nabla['z0'] < 0:
                restrictions['z0'] = (domain.domain_param['Z2'] - cur['z0']) / -nabla['z0']
            elif nabla['z0'] > 0:
                restrictions['z0'] = (cur['z0'] - domain.domain_param['Z1']) / nabla['z0']
            else:
                restrictions['z0'] = 0
        log(f'restrictions: {restrictions}', file_debug)

        # search direction
        nabla = np.asarray([nabla[x] for x in keys])
        D = np.eye(len(keys))  # -> hesse matrix
        d = np.dot(-D, nabla)

        # calc alpha
        sigma = 0.01

        alpha_min = 1
        for k in keys:
            alpha_min = min(alpha_min, restrictions[k])
        alpha = 0.9 * alpha_min
        log(f'mins: 1.0, {restrictions.values()}', file_debug)
        log(f'alpha: {alpha}', file_debug)
        while n > 0:
            x = np.asarray([cur[x] for x in keys])
            new_para = x + (alpha * d)
            new_para = {
                p: new_para[i]
                for i, p in enumerate(keys)
            }
            pprint(new_para)

            file_da.write(f'iterate:{n};')
            diff_cur, _ = do_rollback(client=client,
                                      t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                      new_para=new_para, heat_source=heat_source.values(),
                                      sub_file_name=n,
                                      artss_data_path=artss_data_path, artss=artss,
                                      fds_data=fds_data,
                                      devc_info=devc_info,
                                      file_da=file_da, file_debug=file_debug)
            log(f'conditional statement {diff_orig["T"] + np.dot(alpha * sigma * nabla, d)}', file_debug)
            if diff_cur['T'] < diff_orig['T'] + np.dot(alpha * sigma * nabla, d):
                log(f'found alpha: {alpha}', file_debug)
                log(f'found x_k+1: {new_para}', file_debug)
                log(f"al1: {np.dot(alpha * sigma * nabla, d)}", file_debug)
                log(f"als: {diff_cur['T']} < {diff_orig['T'] + np.dot(alpha * sigma * nabla, d)}", file_debug)
                cur = copy(new_para)
                break

            alpha = 0.7 * alpha
            n = n - 1
            log(f'ls: {diff_cur["T"]}', file_debug)

        log(f'final nabla: {nabla}', file_debug)
        file_da.write(f'final;')
        write_da_data(file_da=file_da, parameters=cur)
        for key in cur:
            heat_source['temperature_source'][key] = cur[key]

        log(f'alpha: {alpha}', file_debug)
        log(f'using cur: {cur}', file_debug)


def continuous_gradient(client: TCP_client,
                        file_da: TextIO, file_debug: TextIO,
                        sensor_times: pandas.Index, devc_info: dict,
                        cur: dict, delta: dict, keys: list,
                        artss_data_path: str, artss: XML, domain: Domain,
                        fds_data: pandas.DataFrame,
                        heat_source: dict,
                        n_iterations: int,
                        precondition: bool = False):
    n = max(1, n_iterations)
    for t_sensor in sensor_times[2:]:
        pprint(cur)
        wait_artss(t_sensor, artss_data_path, artss)

        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=artss.get_dt(), time_back=6)

        log(f't_sensor: {t_sensor} t_artss: {t_artss} t_revert: {t_revert}', file_debug)
        field_reader = FieldReader(t_artss, path=artss_data_path)
        file_da.write(f'original data;time_sensor:{t_sensor};time_artss:{t_artss}\n')
        write_da_data(file_da=file_da, parameters=cur)
        diff_orig, minima_x = comparison_sensor_simulation_data(devc_info, fds_data, field_reader, t_sensor, file_da)

        log(f'org: {diff_orig["T"]}', file_debug)
        log(f't revert: {t_revert}', file_debug)

        if sum(diff_orig.values()) < 1e-5:
            continue

        if precondition:
            precondition = False
            cur['x0'] = minima_x
            heat_source['temperature_source']['x0'] = cur['x0']
            file_da.write(f'preconditioning;')
            diff_orig, _ = do_rollback(client=client,
                                       t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                       new_para={'x0': cur['x0']}, heat_source=heat_source.values(),
                                       sub_file_name='precondition',
                                       artss_data_path=artss_data_path, artss=artss,
                                       fds_data=fds_data,
                                       devc_info=devc_info,
                                       file_da=file_da, file_debug=file_debug)

        nabla = {}

        for p in keys:
            file_da.write(f'calc_nabla;')
            diff_cur, _ = do_rollback(client=client,
                                      t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                      new_para={p: cur[p] + delta[p]}, heat_source=heat_source.values(),
                                      sub_file_name=p,
                                      artss_data_path=artss_data_path, artss=artss,
                                      fds_data=fds_data,
                                      devc_info=devc_info,
                                      file_da=file_da, file_debug=file_debug)
            # calc new nabla
            nabla[p] = (diff_cur['T'] - diff_orig['T']) / delta[p]
            log(f'cur: {diff_cur["T"]}', file_debug)
            log(f'nabla[{p}]: {nabla[p]}', file_debug)

        log(f'nabla without restrictions: {nabla}', file_debug)
        restrictions = {}
        if 'HRR' in keys:
            restrictions['HRR'] = 1

        if 'x0' in keys:
            if nabla['x0'] < 0:
                restrictions['x0'] = (domain.domain_param['X2'] - cur['x0']) / -nabla['x0']
            elif nabla['x0'] > 0:
                restrictions['x0'] = (cur['x0'] - domain.domain_param['X1']) / nabla['x0']
            else:
                restrictions['x0'] = 0
        if 'y0' in keys:
            if nabla['y0'] < 0:
                restrictions['y0'] = (domain.domain_param['Y2'] - cur['y0']) / -nabla['y0']
            elif nabla['y0'] > 0:
                restrictions['y0'] = (cur['y0'] - domain.domain_param['Y1']) / nabla['y0']
            else:
                restrictions['y0'] = 0
        if 'z0' in keys:
            if nabla['z0'] < 0:
                restrictions['z0'] = (domain.domain_param['Z2'] - cur['z0']) / -nabla['z0']
            elif nabla['z0'] > 0:
                restrictions['z0'] = (cur['z0'] - domain.domain_param['Z1']) / nabla['z0']
            else:
                restrictions['z0'] = 0
        log(f'restrictions: {restrictions}', file_debug)

        # search direction
        nabla = np.asarray([nabla[x] for x in keys])
        D = np.eye(len(keys))  # -> hesse matrix
        d = np.dot(-D, nabla)

        # calc alpha
        sigma = 0.01

        alpha_min = 1
        for k in keys:
            alpha_min = min(alpha_min, restrictions[k])
        alpha = 0.9 * alpha_min
        log(f'mins: 1.0, {restrictions.values()}', file_debug)
        log(f'alpha: {alpha}', file_debug)
        while n > 0:
            x = np.asarray([cur[x] for x in keys])
            new_para = x + (alpha * d)
            new_para = {
                p: new_para[i]
                for i, p in enumerate(keys)
            }
            pprint(new_para)

            file_da.write(f'iterate:{n};')
            diff_cur, _ = do_rollback(client=client,
                                      t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                      new_para=new_para, heat_source=heat_source.values(),
                                      sub_file_name=n,
                                      artss_data_path=artss_data_path, artss=artss,
                                      fds_data=fds_data,
                                      devc_info=devc_info,
                                      file_da=file_da, file_debug=file_debug)
            log(f'conditional statement {diff_orig["T"] + np.dot(alpha * sigma * nabla, d)}', file_debug)
            if diff_cur['T'] < diff_orig['T'] + np.dot(alpha * sigma * nabla, d):
                log(f'found alpha: {alpha}', file_debug)
                log(f'found x_k+1: {new_para}', file_debug)
                log(f"al1: {np.dot(alpha * sigma * nabla, d)}", file_debug)
                log(f"als: {diff_cur['T']} < {diff_orig['T'] + np.dot(alpha * sigma * nabla, d)}", file_debug)
                cur = copy(new_para)
                break

            alpha = 0.7 * alpha
            n = n - 1
            log(f'ls: {diff_cur["T"]}', file_debug)

        file_da.write(f'final;')
        write_da_data(file_da=file_da, parameters=cur)
        for key in cur:
            heat_source['temperature_source'][key] = cur[key]

        log(f'alpha: {alpha}', file_debug)
        log(f'using cur: {cur}', file_debug)


def one_time_gradient_parallel(client: TCP_client,
                               file_da: TextIO, file_debug: TextIO,
                               sensor_times: pandas.Index, devc_info: dict,
                               cur: dict, delta: dict, keys: list,
                               artss_path: str, artss_data_path: str,
                               artss: XML, domain: Domain,
                               fds_data: pandas.DataFrame,
                               heat_source: dict,
                               n_iterations: int,
                               precondition: bool = False):
    no_cpus = multiprocessing.cpu_count()
    artss_exe_path = os.path.join(artss_path, 'build', 'bin')
    xml_input_file = os.path.join(artss_data_path, FieldReader.get_xml_file_name(path=artss_data_path))
    for t_sensor in sensor_times[2:]:
        pprint(cur)
        wait_artss(t_sensor, artss_data_path, artss)

        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=artss.get_dt(), time_back=6)

        log(f't_sensor: {t_sensor} t_artss: {t_artss} t_revert: {t_revert}', file_debug)
        field_reader = FieldReader(t_artss, path=artss_data_path)
        file_da.write(f'original data;time_sensor:{t_sensor};time_artss:{t_artss}\n')
        write_da_data(file_da=file_da, parameters=cur)
        diff_orig, minima_x = comparison_sensor_simulation_data(devc_info, fds_data, field_reader, t_sensor, file_da)

        log(f'org: {diff_orig["T"]}', file_debug)
        log(f't revert: {t_revert}', file_debug)

        if sum(diff_orig.values()) < 1e-5:
            continue

        if precondition:
            precondition = False
            cur['x0'] = minima_x
            heat_source['temperature_source']['x0'] = cur['x0']
            file_da.write(f'preconditioning;')
            diff_orig, _ = do_rollback(client=client,
                                       t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                       new_para={'x0': cur['x0']}, heat_source=heat_source.values(),
                                       sub_file_name='precondition',
                                       artss_data_path=artss_data_path, artss=artss,
                                       fds_data=fds_data,
                                       devc_info=devc_info,
                                       file_da=file_da, file_debug=file_debug)

        nabla = {}

        jobs = {}
        directories = {}
        counter = 1
        for p in keys:
            file_da.write(f'calc_nabla;')
            output_file, output_dir = create_start_xml(change={p: cur[p] + delta[p]},
                                                       input_file=xml_input_file, output_file=f'config_{p}.xml',
                                                       t_revert=t_revert, t_end=t_artss,
                                                       directory=artss_data_path,
                                                       file=f'{t_revert:.5e}',
                                                       artss_path=artss_path,
                                                       file_debug=file_debug)
            directories[p] = output_dir
            job = multiprocessing.Process(target=start_new_instance, args=(
                output_file, output_dir, artss_exe_path, 'artss_data_assimilation_serial'))
            job.start()
            counter += 1
            print(f'after start job {p} with name {job.name}, {job}')
            jobs[p] = job
            print(f'added job {p}')
            if counter >= no_cpus:
                # wait if you want to start more multi processes than cpus available
                j = jobs.pop(list(jobs.keys())[0])
                j.join()
        print('calc nabla parallel')
        new_para = {}
        for p in keys:
            if p in jobs.keys():
                job = jobs[p]
                print(f'wait for {job.name} [{job}]')
                job.join()
            diff_cur, _ = calc_diff(t_artss=t_artss, t_sensor=t_sensor,
                                    data_path=directories[p], fds_data=fds_data, devc_info=devc_info,
                                    file_da=file_da)
            print(p, diff_cur['T'])
            nabla[p] = (diff_cur['T'] - diff_orig['T']) / delta[p]
            log(f'nabla[{p}]={nabla[p]}', file_debug)
            if nabla[p] > 0:
                new_para[p] = cur[p] + delta[p]

        log(f'nabla without restrictions: {nabla}', file_debug)
        restrictions = {}
        if 'HRR' in keys:
            restrictions['HRR'] = 1

        if 'x0' in keys:
            if nabla['x0'] < 0:
                restrictions['x0'] = (domain.domain_param['X2'] - cur['x0']) / -nabla['x0']
            elif nabla['x0'] > 0:
                restrictions['x0'] = (cur['x0'] - domain.domain_param['X1']) / nabla['x0']
            else:
                restrictions['x0'] = 0
        if 'y0' in keys:
            if nabla['y0'] < 0:
                restrictions['y0'] = (domain.domain_param['Y2'] - cur['y0']) / -nabla['y0']
            elif nabla['y0'] > 0:
                restrictions['y0'] = (cur['y0'] - domain.domain_param['Y1']) / nabla['y0']
            else:
                restrictions['y0'] = 0
        if 'z0' in keys:
            if nabla['z0'] < 0:
                restrictions['z0'] = (domain.domain_param['Z2'] - cur['z0']) / -nabla['z0']
            elif nabla['z0'] > 0:
                restrictions['z0'] = (cur['z0'] - domain.domain_param['Z1']) / nabla['z0']
            else:
                restrictions['z0'] = 0
        log(f'restrictions: {restrictions}', file_debug)
        # search direction
        nabla = np.asarray([nabla[x] for x in keys])
        D = np.eye(len(keys))  # -> hesse matrix
        d = np.dot(-D, nabla)

        # calc alpha
        sigma = 0.01

        alpha_min = 1
        for k in keys:
            alpha_min = min(alpha_min, restrictions[k])
        alpha = 0.9 * alpha_min
        log(f'mins: 1.0, {restrictions.values()}', file_debug)
        log(f'alpha: {alpha}', file_debug)
        n = n_iterations
        while n > 0:
            x = np.asarray([cur[x] for x in keys])
            new_para = x + (alpha * d)
            new_para = {
                p: new_para[i]
                for i, p in enumerate(keys)
            }
            pprint(new_para)

            file_da.write(f'iterate:{n};')
            diff_cur, _ = do_rollback(client=client,
                                      t_sensor=t_sensor, t_artss=t_artss, t_revert=t_revert,
                                      new_para=new_para, heat_source=heat_source.values(),
                                      sub_file_name=n,
                                      artss_data_path=artss_data_path, artss=artss,
                                      fds_data=fds_data,
                                      devc_info=devc_info,
                                      file_da=file_da, file_debug=file_debug)
            log(f'conditional statement {diff_orig["T"] + np.dot(alpha * sigma * nabla, d)}', file_debug)
            if diff_cur['T'] < diff_orig['T'] + np.dot(alpha * sigma * nabla, d):
                log(f'found alpha: {alpha}', file_debug)
                log(f'found x_k+1: {new_para}', file_debug)
                log(f"al1: {np.dot(alpha * sigma * nabla, d)}", file_debug)
                log(f"als: {diff_cur['T']} < {diff_orig['T'] + np.dot(alpha * sigma * nabla, d)}", file_debug)
                cur = copy(new_para)
                break

            alpha = 0.7 * alpha
            n = n - 1
            log(f'ls: {diff_cur["T"]}', file_debug)


def start(fds_data_path: str, fds_input_file_name: str, artss_data_path: str, artss_path: str, parallel=False):
    artss_path = os.path.abspath(artss_path)
    artss_data_path = os.path.abspath(artss_data_path)
    client = TCP_client.TCPClient()
    # client.set_server_address('172.17.0.2', 7777)
    client.connect()

    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index

    heat_source = xml.get_temperature_source()

    file_da = open('da_details.csv', 'w')
    file_debug = open('da_debug_details.dat', 'w')

    delta = {
        'HRR': float(heat_source['temperature_source']['HRR']) * 0.5,
        'x0': domain.domain_param['dx'] * 10,
        'y0': domain.domain_param['dy'] * 1,
        'z0': domain.domain_param['dz'] * 1
    }

    cur = {
        'HRR': float(heat_source['temperature_source']['HRR']),
        'x0': float(heat_source['temperature_source']['x0']),
        'y0': float(heat_source['temperature_source']['y0']),
        'z0': float(heat_source['temperature_source']['z0'])
    }

    keys = ['HRR', 'x0']
    if parallel:
        continuous_gradient_parallel(client=client, file_da=file_da, file_debug=file_debug,
                                     sensor_times=sensor_times,
                                     devc_info=devc_info_temperature, fds_data=fds_data,
                                     artss_data_path=artss_data_path,
                                     artss_path=artss_path,
                                     domain=domain, heat_source=heat_source,
                                     cur=cur, delta=delta, keys=keys, n_iterations=5, artss=xml,
                                     precondition=False)
    else:
        continuous_gradient(client=client, file_da=file_da, file_debug=file_debug,
                            sensor_times=sensor_times,
                            devc_info=devc_info_temperature, fds_data=fds_data,
                            artss_data_path=artss_data_path,
                            domain=domain, heat_source=heat_source,
                            cur=cur, delta=delta, keys=keys, n_iterations=5, artss=xml,
                            precondition=False)

    # map_minima(client, artss_data_path, cur, delta, sensor_times, devc_info_thermocouple, devc_info_temperature, fds_data, source_type, temperature_source, random, file_da, cwd, xml)

    file_debug.close()
    file_da.close()


def wait_artss(t_sensor, artss_data_path, artss: ARTSS.XML):
    time_step = artss.get_output_time_step()
    if t_sensor % time_step == 0:
        t = t_sensor
    else:
        t = (t_sensor // time_step + 1) * time_step
    # time.sleep(30)  # needed for artss to rewrite meta file
    t_cur = FieldReader.get_t_current(path=artss_data_path)
    pbar = tqdm(total=t)
    pbar.update(t_cur)
    while t_cur <= t:
        time.sleep(1)
        t_new = FieldReader.get_t_current(path=artss_data_path)
        pbar.update(t_new - t_cur)
        t_cur = t_new

    pbar.close()


def write_changes_xml(change: dict, source: list, file_name: str, file_da: TextIO, path='.') -> [str, str]:
    source_type, temperature_source, random = data_assimilation.change_heat_source(*source,
                                                                                   changes={'source_type': {},
                                                                                            'temperature_source': change,
                                                                                            'random': {}})
    write_da_data(file_da=file_da, parameters=temperature_source)
    da = DAFile()
    da.create_config({'u': False, 'v': False, 'w': False, 'p': False, 'T': False, 'C': False})
    da.create_temperature_source_changes(
        source_type=source_type,
        temperature_source=temperature_source,
        random=random)

    config_file_name = f'config_{file_name}.xml'
    config_file_path = os.path.join(path, config_file_name)
    da.write_xml(config_file_path)
    return config_file_name, config_file_path


def get_time_step_artss(t_sensor: float, artss_data_path: str, dt: float, time_back: float = 10) -> [float, float]:
    files = FieldReader.get_all_time_steps(path=artss_data_path)
    times = list(filter(lambda x: x >= t_sensor, files))
    t_revert = ([dt] + list(filter(lambda x: x <= max(0., t_sensor - time_back), files)))[-1]
    return times[0], t_revert


def read_fds_data(input_path: str, input_file_name: str, artss: Domain) -> [dict, dict, pd.DataFrame]:
    full_name = os.path.join(input_path, input_file_name)
    str_devc = read_devc_from_input_file(full_name + '.fds')
    dict_devc = parse_devc(str_devc)
    chid_file_name = get_chid_from_input_file(full_name + '.fds')

    fds_data = read_devc_from_csv_file(os.path.join(input_path, chid_file_name + '_devc.csv'))
    # fds_data = read_fds_file(input_path, artss)

    devc_temperature = {}
    devc_thermocouple = {}
    for d in dict_devc:
        dict_devc[d]['type']: str = 'T'

        ijk = artss.get_ijk_from_xyz(dict_devc[d]['XYZ'][0], dict_devc[d]['XYZ'][2], dict_devc[d]['XYZ'][1])
        dict_devc[d]['index']: int = artss.get_index(*ijk)
        if d.startswith('Temperatur'):
            devc_temperature[d] = dict_devc[d]
        elif d.startswith('Thermocouple'):
            devc_thermocouple[d] = dict_devc[d]

    return devc_temperature, devc_thermocouple, fds_data


def read_devc_from_csv_file(file_name: str) -> pd.DataFrame:
    df = pd.read_csv(file_name, skiprows=1)
    df.index = df.iloc[:, 0]
    df.drop('Time', axis=1, inplace=True)
    return df


def get_chid_from_input_file(input_file: str) -> str:
    with open(input_file, 'r') as inp:
        for line in inp:
            if line.startswith('&HEAD'):
                line = line[len('&HEAD'):].strip()
                tmp_arr = line.split('=')
                index = tmp_arr.index('CHID')
                chid_str = tmp_arr[index + 1].strip()
                index_end = chid_str[1:].index("'")
                chid = chid_str[1:index_end + 1]
                break
    return chid


def parse_devc(str_devc: list, start='&DEVC') -> dict:
    devc = {}
    for line in str_devc:
        line = line[len(start):-1]
        line.strip()
        tmp_arr = line.split(',')
        parameters = []
        for i in tmp_arr:
            if type(i) == str and '=' in i:
                parameters.append(i)
            else:
                parameters[-1] += f',{i}'

        param_dict = {}
        for p in parameters:
            tmp_arr = p.split('=')
            param_dict[tmp_arr[0].strip()] = tmp_arr[1].strip().strip("'").strip('"')
        id = param_dict.pop('ID')
        if 'XYZ' in param_dict.keys():
            param_dict['XYZ'] = [float(i) for i in param_dict['XYZ'].split(',')]
        devc[id] = param_dict
    return devc


def read_devc_from_input_file(fds_file_name: str, start="&DEVC", end="/") -> list:
    """
    Reads FDS input out of a list of string and collects namelist definitions 
    of a given type, e.g. "&DEVC ".
    Reads all definitions line by line, without processing.

    :param fds_file_name: List of string, of a single namelist definition.
    :param start: String, indicate the namelist, default "&DEVC ".
    :param end: String, indicate the end of the namelist, default "/".
    :return: list with lines
    inspired by Tristan Hehnen
    """

    namelist_lines = []

    collect_line = False
    fds_input_content = open(fds_file_name, 'r')
    for line in fds_input_content:
        line = line.strip()
        if collect_line is True:
            namelist_lines[-1] += f' {line}'
        elif line.startswith(start):
            namelist_lines.append(line)
            collect_line = True
        collect_line = collect_line and (end not in line)
    fds_input_content.close()
    return namelist_lines


def read_fds_file(data_path: str, artss: Domain) -> pd.DataFrame:
    data = fdsreader.Simulation(data_path)
    devc = data.devices
    return_val: pd.DataFrame = pd.DataFrame([], index=devc['Time'].data)
    return_val.index.names = ['Time']
    return_val.columns.names = ['Position of Thermocouple']
    for key in devc:
        if key.startswith('Thermocouple'):
            if devc[key] is not None:
                print(key, devc[key].xyz, *devc[key].xyz, devc[key].xyz[1])
                # be careful, in FDS the third parameter indicates the height, but in ARTSS it is the second
                i, j, k = artss.get_ijk_from_xyz(devc[key].xyz[0], devc[key].xyz[2], devc[key].xyz[1])
                index = artss.get_index(i, j, k)
                print(f'index: {index} of {i}|{j}|{k}')
                if index in return_val:
                    print(f'index  {index} already exists, take average')
                    tmp = (devc[key].data + return_val[index]) / 2.
                    return_val[index] = tmp
                else:
                    return_val[index] = devc[key].data
            else:
                print(f'{key} is None, skipped')
    return return_val


def comparison_sensor_simulation_data(devc_info: dict, sensor_data: pd.DataFrame, field_reader: FieldReader,
                                      t_sensor: float, file_da: TextIO) -> (Dict[str, float], float):
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
        if abs(difference) > 0.5:
            diff[type_sensor].append(difference)
        file_da.write(
            f'sensor:{key};time_sensor:{t_sensor};time_artss:{field_reader.t};sensor_val:{value_sensor};artss_val:{value_sim};diff:{difference};considered:{abs(difference) > 0.5};position:{devc_info[key]["XYZ"]}\n')
        file_da.flush()

        if difference < min_sensor_val:
            min_sensor_val = difference
            min_pos_x = devc_info[key]['XYZ'][0]
    print(f'diff: {diff["T"]}')
    result = {}
    for key in diff:
        if len(diff[key]) == 0:
            continue
        result[key] = np.sqrt(sum(np.array(diff[key]) ** 2))
        file_da.write(f'differences:{key}:{result[key]};no_of_sensors:{len(diff[key])}\n')
    return result, min_pos_x


def calc_single_differences(devc_info: dict, sensor_data: pd.DataFrame, field_reader: FieldReader, t_artss: float,
                            t_sensor: float) -> dict:
    diff: dict = {}
    fields_sim = field_reader.read_field_data()
    tmp_x = []
    tmp_y = []
    tmp_y2 = []
    index = []
    for key in devc_info:
        type_sensor: str = devc_info[key]['type']
        index_sensor: int = devc_info[key]['index']
        value_sensor: float = sensor_data[key][t_sensor]

        value_sim = fields_sim[type_sensor][index_sensor]
        diff[key] = (kelvin_to_celsius(value_sim) - value_sensor)
        tmp_x.append(key)
        index.append(index_sensor)
        tmp_y.append(value_sim)
        tmp_y2.append(value_sensor)

    # print(f' calc single diff {t_artss}, {t_sensor}')
    # print(index)
    # print(tmp_x)
    # print(tmp_y)
    # print(tmp_y2)
    return diff


def kelvin_to_celsius(kelvin):
    return kelvin - 273.5


def plot_differences(fds_data_path: str, fds_input_file_name: str, artss_data_path: str):
    cwd = os.getcwd()
    a_path = os.path.join(artss_data_path, '30')
    xml = XML(FieldReader.get_xml_file_name(path=a_path), path=a_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index
    print(sensor_times)
    for t in range(1, 8):
        t_sensor = sensor_times[t]
        fig, axs = plt.subplots(2, 1)
        for i in [20, 30, 40, 50, 60, 70, 80]:
            t_artss, t_revert = get_time_step_artss(t_sensor, os.path.join(artss_data_path, str(i)))
            field_reader = FieldReader(t_artss, path=os.path.join(artss_data_path, str(i)))
            fields = field_reader.read_field_data()
            diff = calc_single_differences(devc_info_temperature, fds_data, field_reader, t_artss, t_sensor)
            x = []
            y = []
            y_fds = []
            y_artss = []
            tmp_dict = {}
            for k in devc_info_temperature:
                x.append(devc_info_temperature[k]['XYZ'][0])
                y.append(diff[k])
                y_fds.append(fds_data[k][t_sensor])
                if x[-1] in tmp_dict:
                    tmp_dict[x[-1]].append(y[-1])
                else:
                    tmp_dict[x[-1]] = [y[-1]]
                index = domain.get_index(
                    *domain.get_ijk_from_xyz(devc_info_temperature[k]['XYZ'][0], devc_info_temperature[k]['XYZ'][2],
                                             devc_info_temperature[k]['XYZ'][1]))
                y_artss.append(kelvin_to_celsius(fields['T'][index]))
            line_x = []
            line_y = []
            for k in tmp_dict:
                line_x.append(k)
                line_y.append(sum(tmp_dict[k]) / len(tmp_dict[k]))
            axs[0].plot(line_x, line_y, label=f'{i} avg: {sum(y) / len(y):.2e}')
            axs[0].scatter(x, y, label=i)
            axs[1].plot(x, y_artss, label=f'{i}')
        axs[0].set_title('difference')
        axs[1].set_title('sensor data')
        axs[1].plot(x, y_fds, label='fds')
        plt.xlabel('position of temperature sensor')
        # axs[0].ylabel('difference [C]')
        # axs[1].ylabel('T [C]')
        axs[0].legend()
        axs[1].legend()
        plt.savefig(f'diff_{t_artss}.pdf')
        plt.close()


def process_data(devc_info_temperature: dict, artss_times: list, artss_data_path: str):
    f = open('artss_data.dat', 'w', buffering=1)
    artss_data = {}
    for path in ['with_da', 'without_da']:
        artss_data[path] = {}
        for sensor in devc_info_temperature:
            artss_data[path][sensor] = []
        for t_artss in artss_times:
            print('t:', t_artss)
            field_reader = FieldReader(t_artss, path=os.path.join(artss_data_path, path))
            fields_sim = field_reader.read_field_data()
            for sensor in devc_info_temperature:
                index_sensor: int = devc_info_temperature[sensor]['index']
                value_artss = kelvin_to_celsius(fields_sim['T'][index_sensor])
                artss_data[path][sensor].append(value_artss)
        for sensor in devc_info_temperature:
            f.write(f'sensor:{sensor}')
            f.write(f'x:{artss_times}\n')
            f.write(f'y:{artss_data[path][sensor]}\n')
            f.write(f'label:{path}\n')
    f.close()
    return artss_data


def plot_comparison_da(fds_data_path: str, fds_input_file_name: str, artss_data_path: str):
    cwd = os.getcwd()
    a_path = os.path.join(artss_data_path, 'with_da')
    xml = XML(FieldReader.get_xml_file_name(path=a_path), path=a_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index[:31]
    print(sensor_times)
    print('artss path', artss_data_path)
    artss_times = FieldReader.get_all_time_steps(a_path)[:501]
    print('artss times', artss_times)
    print('sensor times', sensor_times)
    print('first row fds', fds_data.columns)

    artss_data = process_data(devc_info_temperature, artss_times, artss_data_path)

    for sensor in devc_info_temperature:
        for path in artss_data:
            plt.plot(artss_times, artss_data[path][sensor], label=path.replace('_', ' '))
            print(type(artss_times))
            print(type(artss_data[path][sensor]))
        fds_values = []
        for t in sensor_times:
            fds_values.append(fds_data[sensor.replace('Temperature', 'Thermocouple')][t])
        plt.plot(sensor_times, fds_values, label='fds thermocouple')
        print(type(fds_values))

        fds_values = []
        for t in sensor_times:
            fds_values.append(fds_data[sensor][t])
        plt.plot(sensor_times, fds_values, label='fds temperature')
        print(type(fds_values))
        plt.legend()
        plt.ylabel('temperature [°C]')
        plt.xlabel('time [s]')
        sensor_xpos = devc_info_temperature[sensor]['XYZ'][0]
        plt.savefig(f'comparison_sensor_{sensor_xpos}.pdf')
        plt.close()


def plot_sensor_data(fds_data_path: str, fds_input_file_name: str, artss_data_path: str):
    cwd = os.getcwd()
    a_path = os.path.join(artss_data_path, '30')
    xml = XML(FieldReader.get_xml_file_name(path=a_path), path=a_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index[:13]
    print(sensor_times)
    print('artss path', artss_data_path)
    artss_times = FieldReader.get_all_time_steps(a_path)
    print('artss times', artss_times)
    print('sensor times', sensor_times)

    tmp_dict = {}
    axis_name = set()
    for d in devc_info_temperature:
        xpos = devc_info_temperature[d]['XYZ'][0]
        axis_name.add(xpos)
        if xpos in tmp_dict:
            tmp_dict[xpos].append(d)
        else:
            tmp_dict[xpos] = [d]

    for key in tmp_dict:
        print(key, tmp_dict[key])

    axis_name = list(axis_name)
    axis_name.sort()

    artss_data = {}
    for x_pos in [20, 30, 40, 50, 60, 70, 80]:
        artss_data[x_pos] = {}
        for sensor_x in axis_name:
            artss_data[x_pos][sensor_x] = [[], []]
        for t_artss in artss_times:
            field_reader = FieldReader(t_artss, path=os.path.join(artss_data_path, str(x_pos)))
            fields_sim = field_reader.read_field_data()
            for sensor_x in axis_name:
                for i in range(len(tmp_dict[sensor_x])):
                    index_sensor: int = devc_info_temperature[tmp_dict[sensor_x][i]]['index']
                    value_artss = kelvin_to_celsius(fields_sim['T'][index_sensor])
                    artss_data[x_pos][sensor_x][i].append(value_artss)

    for index, sensor_x in enumerate(axis_name):
        fig, axs = plt.subplots(2, 1)
        for i in range(len(tmp_dict[sensor_x])):
            key = tmp_dict[sensor_x][i]
            key_temp = key.replace('Thermocouple', 'Temperature')
            for x_pos in [30, 40, 50, 60, 70, 80]:
                axs[i].plot(artss_times, artss_data[x_pos][sensor_x][i], label=f'artss xpos {x_pos}')

            fds_values = []
            for t in sensor_times:
                fds_values.append(fds_data[key][t])
            axs[i].plot(sensor_times, fds_values, label='fds thermocouple')

            fds_values = []
            for t in sensor_times:
                fds_values.append(fds_data[key_temp][t])
            axs[i].plot(sensor_times, fds_values, label='fds temperature')

        plt.xlabel('t [s]')
        plt.ylabel('T [C]')
        plt.legend()
        plt.savefig(f'verlauf_sensor_{sensor_x}.pdf')
        plt.close()


def compare_distance(fds_data_path: str, fds_input_file_name: str, a_data_path: str, artss_path: str):
    artss_data_path = os.path.abspath(os.path.join(a_data_path, 'distance'))
    artss_exe_path = os.path.abspath(os.path.join(artss_path, 'build', 'bin'))
    xml_input_file = os.path.join(artss_data_path, FieldReader.get_xml_file_name(path=artss_data_path))

    file_da = open(os.path.join(artss_data_path, 'da_details.csv'), 'w')
    file_plot = open(os.path.join(artss_data_path, 'da_plot_details.csv'), 'w')
    file_plot.write(f'time;x pos;diff\n')
    file_debug = open(os.path.join(artss_data_path, 'da_debug_details.dat'), 'w')

    xml = XML(FieldReader.get_xml_file_name(path=artss_data_path), path=artss_data_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    log(f'devc info: {devc_info_temperature}', file_debug)

    x_pos_orig = 52.5
    sensor_times = fds_data.index[3:5]
    log(sensor_times, file_debug)
    log(f'artss path {artss_data_path}', file_debug)
    artss_times = FieldReader.get_all_time_steps(artss_data_path)
    log(f'artss times {artss_times}', file_debug)
    log(f'sensor times {sensor_times}', file_debug)

    heat_source = xml.get_temperature_source()
    x0 = float(heat_source['temperature_source']['x0'])
    dx = float(domain.domain_param['dx'])
    sigma_x = float(heat_source['temperature_source']['sigma_x'])

    steps = np.array(range(-20, 20))
    delta_dx = steps * dx
    delta_sigma_x = steps * sigma_x
    log(f'sigma x: {delta_sigma_x}', file_debug)
    log(f'x: {delta_dx}', file_debug)

    no_cpus = multiprocessing.cpu_count()
    devc_info = devc_info_temperature

    for t_sensor in sensor_times:
        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(), time_back=6)

        log(f't_sensor: {t_sensor} t_artss: {t_artss} t_revert: {t_revert}', file_debug)
        field_reader = FieldReader(t_artss, path=artss_data_path)
        diff_orig, minima_x = comparison_sensor_simulation_data(devc_info, fds_data, field_reader, t_sensor, file_da)
        file_plot.write(f'{t_sensor};{x_pos_orig};{diff_orig}\n')
        log(f'org: {diff_orig["T"]}', file_debug)
        log(f't revert: {t_revert}', file_debug)

        jobs = {}
        directories = {}
        counter = 1

        l_diff = []
        l_xpos = []
        l_xpos_r = []

        for delta in delta_dx:
            x_pos = x0 + delta
            log(f'xpos: {x_pos}', file_debug)
            output_file, output_dir = create_start_xml(change={'x0': x_pos},
                                                       input_file=xml_input_file, output_file=f'config_{delta:.5e}.xml',
                                                       t_revert=t_revert, t_end=t_artss,
                                                       directory=artss_data_path,
                                                       file=f'{t_revert:.5e}',
                                                       artss_path=artss_path,
                                                       file_debug=file_debug,
                                                       dir_name=f'x_{x_pos}')
            directories[delta] = output_dir
            job = multiprocessing.Process(target=start_new_instance, args=(
                output_file, output_dir, artss_exe_path, 'artss_data_assimilation_serial'))
            job.start()
            counter += 1
            jobs[delta] = job
            if counter >= no_cpus:
                # wait if you want to start more multi processes than cpus available
                j = jobs.pop(list(jobs.keys())[0])
                j.join()

        for delta in delta_sigma_x:
            x_pos = x0 + delta
            if delta in jobs.keys():
                job = jobs[delta]
                print(f'wait for {job.name} [{job}]')
                job.join()
            file_da.write(f'sigma_x;{x_pos};{delta}\n')
            diff_cur, _ = calc_diff(t_artss=t_artss, t_sensor=t_sensor,
                                    data_path=directories[delta], fds_data=fds_data, devc_info=devc_info,
                                    file_da=file_da)
            l_diff.append(diff_cur['T'])
            l_xpos.append(x_pos)
            l_xpos_r.append(delta)
            file_plot.write(f'{t_artss};{x_pos};{diff_cur["T"]}\n')

        file_plot.flush()
        file_da.flush()
        log(f'x_pos: {l_xpos}\ndiff: {l_diff}', file_debug)
        file_debug.flush()

        tmp_int = int(t_artss)
        plt.plot(l_xpos, l_diff)
        plt.xlabel('x position of heat source')
        plt.ylabel('difference to FDS data')
        plt.title('Move heat source based on cell width dx')
        plt.savefig(f'distance_dx_{tmp_int:02d}.pdf', format='pdf')
        plt.close()

        plt.plot(l_xpos_r, l_diff)
        plt.title('Move heat source based on cell width dx')
        plt.xlabel('relative x position of heat source')
        plt.ylabel('difference to FDS data')
        plt.savefig(f'distance_dx_{tmp_int:02d}_relative.pdf', format='pdf')
        plt.close()
