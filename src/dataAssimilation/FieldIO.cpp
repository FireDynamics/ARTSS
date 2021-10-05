/// \file       FieldIO.cpp
/// \brief      
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.

#include "FieldIO.h"
#include <chrono>
#include <ctime>
#include <algorithm>
#include <iomanip>

#include <fmt/format.h>
#include "../Domain.h"
#include "DataAssimilation.h"
#include "Zero.h"
#include "HRRChanger.h"


FieldIO::FieldIO(const SolverController &solver_controller): m_solver_controller(solver_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    Parameters *params = Parameters::getInstance();
    std::string init = params->get("data_assimilation/class_name");
    if (init == AssimilationMethods::None) {
        m_func = new Zero();
    } else if (init == AssimilationMethods::HRRChanger) {
        m_func = new HRRChanger(m_solver_controller.get_temperature_source_function());
    } else {
#ifndef BENCHMARKING
        m_logger->critical("Data Assimilation class {} is not defined", init);
#endif
        std::exit(1);
        // TODO Error Handling
    }

    std::string dt = Parameters::getInstance()->get("physical_parameters/dt");
    m_dt = std::stod(dt);
    std::string t_end = Parameters::getInstance()->get("physical_parameters/t_end");

    //std::vector<std::string> t_end_parts = Utility::split(Utility::trim(t_end), '.');
    //std::vector<std::string> dt_parts = Utility::split(Utility::trim(dt), '.');

    int decimal_number_digits;
    int whole_number_digits;
    Utility::calculate_whole_and_decimal_digits(t_end, dt, &decimal_number_digits, &whole_number_digits);

    m_length_time_stamp = decimal_number_digits + whole_number_digits + 1;
    m_format = "{:0" + std::to_string(m_length_time_stamp) + "." + std::to_string(decimal_number_digits) + "f}";

    size_t n = static_cast<size_t>(std::stod(t_end) / m_dt) + 1;
    m_positions = new long[n];

    std::string header = create_header();
    std::string header_extra = m_func->add_header_data();
    header += header_extra +  "\n";
    m_positions[0] = static_cast<long>(header.length()) + 1;

    std::ofstream output_file(m_filename, std::ios_base::out);
    output_file.write(header.c_str(), m_positions[0]);
    output_file.close();
}


// ========================================== write ================================================
// *************************************************************************************************
/// \brief  write fields u, v, w, p, T and C into the predefined file
/// \param  t_cur   time step from which the fields should be written out
/// \param  u       data of field u to be written out
/// \param  v       data of field v to be written out
/// \param  w       data of field w to be written out
/// \param  p       data of field p to be written out
/// \param  T       data of field T to be written out
/// \param  C       data of field C to be written out
// *************************************************************************************************
void FieldIO::write(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C) {
    std::string output = fmt::format(m_format + "\n", t_cur);
    Field fields[] = {u, v, w, p, T, C};
    size_t size = u.get_size();
    for (Field &f: fields){
        for (size_t i = 0; i < size - 1; i++) {
            output.append(fmt::format(("{};"), f[i]));
        }
        output.append(fmt::format(("{}\n"), f[size - 1]));
    }
    size_t n = static_cast<size_t>(t_cur / m_dt) - 1;
    long length = static_cast<long>(output.length());
    m_positions[n + 1] = m_positions[n] + length;

    std::fstream output_file(m_filename);
    // write field at position dependent on time step
    output_file.seekp(m_positions[n], std::ios_base::beg);
    output_file.write(output.c_str(), length);

    std::string body_extra = m_func->update_body_data();
    output_file.write(body_extra.c_str(), static_cast<long>(body_extra.length()));

    // overwrite current time step
    output_file.seekp(m_pos_time_step, std::ios_base::beg);
    output_file.write(fmt::format(m_format, t_cur).c_str(), m_length_time_stamp);

    output_file.seekp(m_pos_header_extra, std::ios_base::beg);
    std::string header_extra = m_func->update_header_data();
    output_file.write(header_extra.c_str(), static_cast<long>(header_extra.length()));

    output_file.close();
}

// ========================================== read =================================================
// *************************************************************************************************
/// \brief  read fields u, v, w, p, T and C from the predefined file
/// \param  t_cur   time step which should be read
/// \param  u       field u to store the read data
/// \param  v       field v to store the read data
/// \param  w       field w to store the read data
/// \param  p       field p to store the read data
/// \param  T       field T to store the read data
/// \param  C       field C to store the read data
// *************************************************************************************************
void FieldIO::read(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C) {
    std::ifstream input_file(m_filename, std::ifstream::binary);
    size_t n = static_cast<size_t>(t_cur / m_dt) - 1;
    long pos = m_positions[n];
    std::string line;
    input_file.seekg(pos);

    Field fields[] = {u, v, w, p, T, C};
    for (Field &f: fields) {
        getline(input_file, line);
        std::vector<std::string> splitted_string = Utility::split(line, ';');
        size_t counter = 0;
        for (const std::string &part: splitted_string) {
            f.data[counter] = std::stod(part);
            counter++;
        }
    }
}

// ========================================== read =================================================
// *************************************************************************************************
/// \brief  read fields u, v, w, p, T and C from the given file
/// \param  file_name   file in which the data is stored
/// \param  u           field u to store the read data
/// \param  v           field v to store the read data
/// \param  w           field w to store the read data
/// \param  p           field p to store the read data
/// \param  T           field T to store the read data
/// \param  C           field C to store the read data
// *************************************************************************************************
void FieldIO::read(std::string &file_name, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C) {
    std::ifstream input_file(file_name, std::ifstream::binary);
    if (input_file.is_open()) {
        std::string string_time_step;
        getline(input_file, string_time_step);
        getline(input_file, string_time_step);
        getline(input_file, string_time_step);
        getline(input_file, string_time_step);
        std::string line;
        Field fields[] = {u, v, w, p, T, C};
        for (Field &f: fields) {
            getline(input_file, line);
            std::vector<std::string> splitted_string = Utility::split(line, ';');
            size_t counter = 0;
            for (const std::string &part: splitted_string) {
                f.data[counter] = static_cast<real>(std::stold(part));
                counter++;
            }
        }
    } else {
#ifndef BENCHMARKING
        m_logger->warn("File {} not found, skipping assimilation", file_name);
#endif
    }
}

//================================= write header ===================================================
// *************************************************************************************************
/// \brief  write header for field storage. includes essential information: Nx, Ny, Nz
/// \details header format:
/// ###DOMAIN;<Nx>;<Ny>;<Nz>
/// ###FIELDS;u;v;w;p;T;concentration
/// ###DATE:<date>;XML:<XML>
// *************************************************************************************************
std::string FieldIO::create_header() {
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();
    auto Nz = domain->get_Nz();

    std::string string_t_cur_text = "Current time step;";
    m_pos_time_step = static_cast<long>(string_t_cur_text.length());
    std::string header = fmt::format(string_t_cur_text + m_format + ";dt;{}\n", 0.0, m_dt);
    header.append(fmt::format("###DOMAIN;{};{};{}\n", Nx, Ny, Nz));
    header.append(fmt::format("###FIELDS;u;v;w;p;T;concentration\n"));
    header.append(fmt::format("###DATE;{}", std::ctime(&end_time)));
    header.append(fmt::format("###XML;{}\n", Parameters::getInstance()->get_filename()));
    return header;
}

