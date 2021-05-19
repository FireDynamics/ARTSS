/// \file       FieldIO.cpp
/// \brief      
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.

#include "FieldIO.h"
#include <chrono>
#include <ctime>

#include "../Domain.h"


FieldIO::FieldIO(FieldController *field_controller) : m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_dt = Parameters::getInstance()->get_real("physical_parameters/dt");
    real t_end = Parameters::getInstance()->get_real("physical_parameters/t_end");
    size_t n = static_cast<size_t>(t_end / m_dt) + 1;
    m_positions = new size_t[n];
    const char* header = create_header().c_str();
    std::ofstream output_file(m_filename, std::ios_base::out);
    output_file.write(header, m_positions[0]);
}

void FieldIO::write_out(real t_cur) {
    std::string output;

    auto u = m_field_controller->get_field_u_data();
    auto v = m_field_controller->get_field_v_data();
    auto w = m_field_controller->get_field_w_data();
    auto p = m_field_controller->get_field_w_data();
    auto T = m_field_controller->get_field_T_data();
    auto C = m_field_controller->get_field_concentration_data();

    int t = 0;
    real *fields[] = {u, v, w, p, T, C};
    auto size = Domain::getInstance()->get_size();
    for (auto f: fields){
        for (size_t i = 0; i < size - 1; i++) {
            output.append(fmt::format(("{};"), f[i]));
        }
        output.append(fmt::format(("{}\n"), f[size - 1]));
    }

    size_t n = static_cast<size_t>(t_cur / m_dt) - 1;
    size_t length = output.length();
    m_positions[n] = m_positions[n] + length;

    std::ofstream output_file(m_filename, std::ios_base::app);
    output_file.seekp(m_positions[n]);
    output_file.write(output.c_str(), length);
}

void FieldIO::read(real t_cur, Field *u, Field *v, Field *w, Field *p, Field *T, Field *C) {
    std::ifstream input_file(m_filename, std::ifstream::binary);
    size_t n = static_cast<size_t>(t_cur / m_dt) - 1;
    size_t pos = m_positions[n];
    size_t length = m_positions[n + 1] - m_positions[n];
    char *buffer = new char [length];
    std::string line;
    input_file.seekg(pos);

    Field *fields[] = {u, v, w, p, T, C};
    for (Field *f: fields) {
        getline(input_file, line);
        std::strcpy(buffer, line.c_str());
        char *token_delimiter = strtok(buffer, ";");
        size_t counter = 0;
        while (token_delimiter != nullptr) {
            f->data[counter] = atof(token_delimiter);
            token_delimiter = strtok(nullptr, ";");
            counter++;
        }
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

    std::string header = fmt::format("###DOMAIN;{};{};{}\n", Nx, Ny, Nz);
    header.append(fmt::format(("###FIELDS;u;v;w;p;T;concentration\n")));
    header.append(fmt::format(("###DATE:{};XML:{}\n"),
                             std::ctime(&end_time), Parameters::getInstance()->get_filename()));
    m_positions[0] = header.length() + 1;
    return header;
}
