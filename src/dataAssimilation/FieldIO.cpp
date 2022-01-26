/// \file       FieldIO.cpp
/// \brief      Class for reading/writing the raw data of the fields u, v, w, p, T, C
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.

#include "FieldIO.h"
#include <chrono>
#include <ctime>
#include <algorithm>
#include <iomanip>

#include <fmt/compile.h>
#include <fstream>
#include "../domain/DomainData.h"


FieldIO::FieldIO(const std::string &xml_file_name, const std::string &output_file_name) :
        m_file_name(output_file_name) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    auto domain_data = DomainData::getInstance();
    real t_end = domain_data->get_physical_parameters().t_end;

    size_t n = static_cast<size_t>(std::round(t_end / domain_data->get_physical_parameters().dt)) + 1;
    m_positions = new long[n];

    std::string header = create_header(xml_file_name);
    m_positions[0] = static_cast<long>(header.length()) + 1;

    std::ofstream output_file(m_file_name, std::ios_base::out);
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
    std::string output = fmt::format("{:.5e}\n", t_cur);
    Field fields[] = {u, v, w, p, T, C};
    size_t size = u.get_size();
    for (Field &f: fields) {
        for (size_t i = 0; i < size - 1; i++) {
            output.append(fmt::format(("{};"), f[i]));
        }
        output.append(fmt::format(("{}\n"), f[size - 1]));
    }
    size_t n = static_cast<size_t>(std::round(t_cur / DomainData::getInstance()->get_physical_parameters().dt)) - 1;
    long length = static_cast<long>(output.length());
    m_positions[n + 1] = m_positions[n] + length;

    std::fstream output_file(m_file_name);
    // write field at position dependent on time step
    m_logger->debug("times: {:>10d} write to: {:>20d}", n, m_positions[n]);
    output_file.seekp(m_positions[n], std::ios_base::beg);
    output_file.write(output.c_str(), length);

    // overwrite current time step
    output_file.seekp(m_pos_time_step, std::ios_base::beg);
    output_file.write(fmt::format("{:.5e}", t_cur).c_str(), 11);

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
    m_logger->debug("read original data");
    std::ifstream input_file(m_file_name, std::ifstream::binary);
    size_t n = static_cast<size_t>(std::round(t_cur / DomainData::getInstance()->get_physical_parameters().dt)) - 1;
    long pos = m_positions[n];
    m_logger->debug("times: {:>10d} read from: {:>20d}", n, m_positions[n]);
    std::string line;
    input_file.seekg(pos);

    std::getline(input_file, line);
    m_logger->debug("read time step {}", line);

    read_field(input_file, u);
    read_field(input_file, v);
    read_field(input_file, w);
    read_field(input_file, p);
    read_field(input_file, T);
    read_field(input_file, C);
}

// ================================= write header ==================================================
// *************************************************************************************************
/// \brief  write header for field storage. includes essential debugger notation: Nx, Ny, Nz
/// \details header format:
/// ###DOMAIN;<Nx>;<Ny>;<Nz>
/// ###FIELDS;u;v;w;p;T;concentration
/// ###DATE:<date>;XML:<XML>
// *************************************************************************************************
std::string FieldIO::create_header(const std::string &xml_file_name) {
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_number_of_cells(CoordinateAxis::X);
    size_t Ny = domain_data->get_number_of_cells(CoordinateAxis::Y);
    size_t Nz = domain_data->get_number_of_cells(CoordinateAxis::Z);

    std::string string_t_cur_text = "Current time step;";
    m_pos_time_step = static_cast<long>(string_t_cur_text.length());
    std::string header = fmt::format("{}{:.5e};dt;{}\n", string_t_cur_text, 0.0,
                                     DomainData::getInstance()->get_physical_parameters().dt);
    header.append(fmt::format("###DOMAIN;{};{};{}\n", Nx, Ny, Nz));
    header.append(fmt::format("###FIELDS;u;v;w;p;T;concentration\n"));
    header.append(fmt::format("###DATE;{}", std::ctime(&end_time)));
    header.append(fmt::format("###XML;{}\n", xml_file_name));
    return header;
}

void FieldIO::read_field(std::ifstream &file_stream, Field &field) {
    std::string line;
    auto to_double = [](const auto &v) { return std::stod(v); };

    std::getline(file_stream, line);
    std::vector<std::string> divided_string = Utility::split(line, ';');
    std::transform(divided_string.begin(), divided_string.end(), field.data, to_double);
}

void FieldIO::read_fields(const Settings::data_assimilation::field_changes &field_changes,
                          Field &u, Field &v, Field &w,
                          Field &p, Field &T, Field &C) {
    std::string line;
    if (!field_changes.changed) {
        return;
    }
    // no changes -> read original file
    std::ifstream file_changes(field_changes.file_name, std::ifstream::binary);
    if (file_changes.is_open()) {  // could not open file -> read original file + warning
        if (field_changes.u_changed) {
            read_field(file_changes, u);
            m_logger->debug("read changed u Field");
            std::getline(file_changes, line);
        }
        if (field_changes.v_changed) {
            read_field(file_changes, v);
            m_logger->debug("read changed v Field");
        }
        if (field_changes.w_changed) {
            read_field(file_changes, w);
            m_logger->debug("read changed w Field");
        }
        if (field_changes.p_changed) {
            read_field(file_changes, p);
            m_logger->debug("read changed p Field");
        }
        if (field_changes.T_changed) {
            read_field(file_changes, T);
            m_logger->debug("read changed T Field");
        }
        if (field_changes.C_changed) {
            read_field(file_changes, C);
            m_logger->debug("read changed C Field");
        }
    } else {
        m_logger->warn(fmt::format("File '{}' could not be opened. No changes will be applied.", field_changes.file_name));
    }
}

void FieldIO::read_fields(const real t_cur,
                          const Settings::data_assimilation::field_changes &field_changes,
                          Field &u, Field &v, Field &w,
                          Field &p, Field &T, Field &C) {
    std::ifstream file_original(m_file_name, std::ifstream::binary);
    size_t n = static_cast<size_t>(std::round(t_cur / DomainData::getInstance()->get_physical_parameters().dt)) - 1;
    file_original.seekg(m_positions[n]);

    std::string line;
    std::getline(file_original, line);
    m_logger->debug("read time step {}", line);

    if (field_changes.changed) {  // no changes -> read original file
        std::ifstream file_changes(field_changes.file_name, std::ifstream::binary);
        if (file_changes.is_open()) {  // could not open file -> read original file + warning
            if (field_changes.u_changed) {
                read_field(file_changes, u);
                std::getline(file_original, line);
                m_logger->debug("read changed u Field");
            } else {
                read_field(file_original, u);
                std::getline(file_changes, line);
            }
            if (field_changes.v_changed) {
                read_field(file_changes, v);
                std::getline(file_original, line);
                m_logger->debug("read changed v Field");
            } else {
                read_field(file_original, v);
                std::getline(file_changes, line);
            }
            if (field_changes.w_changed) {
                read_field(file_changes, w);
                std::getline(file_original, line);
                m_logger->debug("read changed w Field");
            } else {
                read_field(file_original, w);
                std::getline(file_changes, line);
            }
            if (field_changes.p_changed) {
                read_field(file_changes, p);
                std::getline(file_original, line);
                m_logger->debug("read changed p Field");
            } else {
                read_field(file_original, p);
                std::getline(file_changes, line);
            }
            if (field_changes.T_changed) {
                read_field(file_changes, T);
                std::getline(file_original, line);
                m_logger->debug("read changed T Field");
            } else {
                read_field(file_original, T);
                std::getline(file_changes, line);
            }
            if (field_changes.C_changed) {
                read_field(file_changes, C);
                std::getline(file_original, line);
                m_logger->debug("read changed C Field");
            } else {
                read_field(file_original, C);
                std::getline(file_changes, line);
            }
        } else {
            m_logger->warn(fmt::format("File '{}' could not be opened, original data will be loaded", field_changes.file_name));
            read_fields(t_cur, u, v, w, p, T, C);
        }
    } else {
        m_logger->debug("no field changes");
        read_fields(t_cur, u, v, w, p, T, C);
    }
}

