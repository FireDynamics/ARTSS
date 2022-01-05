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
#include <fstream>
#include "../domain/DomainData.h"


FieldIO::FieldIO() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    auto domain_data = DomainData::getInstance();
    real t_end = domain_data->get_physical_parameters().t_end;

    size_t n = static_cast<size_t>(t_end / domain_data->get_physical_parameters().dt) + 1;
    m_positions = new long[n];

    std::string header = create_header();
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
void FieldIO::write_fields(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C) {
    std::string output = fmt::format(m_format + "\n", t_cur);
    Field fields[] = {u, v, w, p, T, C};
    size_t size = u.get_size();
    for (Field &f: fields){
        for (size_t i = 0; i < size - 1; i++) {
            output.append(fmt::format(("{};"), f[i]));
        }
        output.append(fmt::format(("{}\n"), f[size - 1]));
    }
    size_t n = static_cast<size_t>(t_cur / DomainData::getInstance()->get_physical_parameters().dt) - 1;
    long length = static_cast<long>(output.length());
    m_positions[n + 1] = m_positions[n] + length;

    std::fstream output_file(m_filename);
    // write field at position dependent on time step
    output_file.seekp(m_positions[n], std::ios_base::beg);
    output_file.write(output.c_str(), length);

    // overwrite current time step
    output_file.seekp(m_pos_time_step, std::ios_base::beg);
    output_file.write(fmt::format(m_format, t_cur).c_str(), m_length_time_stamp);

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
void FieldIO::read_fields(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C) {
    std::ifstream input_file(m_filename, std::ifstream::binary);
    size_t n = static_cast<size_t>(t_cur / DomainData::getInstance()->get_physical_parameters().dt) - 1;
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
void FieldIO::read_fields(const std::string &file_name, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C) {
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
        m_logger->warn("File {} not found, skipping data_assimilation", file_name);
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

    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_number_of_cells(CoordinateAxis::X);
    size_t Ny = domain_data->get_number_of_cells(CoordinateAxis::Y);
    size_t Nz = domain_data->get_number_of_cells(CoordinateAxis::Z);

    std::string string_t_cur_text = "Current time step;";
    m_pos_time_step = static_cast<long>(string_t_cur_text.length());
    std::string header = fmt::format(string_t_cur_text + m_format + ";dt;{}\n", 0.0, DomainData::getInstance()->get_physical_parameters().dt);
    header.append(fmt::format("###DOMAIN;{};{};{}\n", Nx, Ny, Nz));
    header.append(fmt::format("###FIELDS;u;v;w;p;T;concentration\n"));
    header.append(fmt::format("###DATE;{}", std::ctime(&end_time)));
    return header;
}

void FieldIO::read_fields(const std::string &file_name,
                          const real t_cur,
                          std::vector<FieldType> fields,
                          Field &u, Field &v, Field &w,
                          Field &p, Field &T, Field &C) {

    if (!fields.empty()) {
        // TODO read onl the provided files
        read_fields(file_name, u, v, w, p, T, C);
    } else {
        read_fields(t_cur, u, v, w, p, T, C);
    }

}

