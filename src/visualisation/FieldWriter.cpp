/// \file       FieldWriter.cpp
/// \brief      
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.

#include "FieldWriter.h"

#include <chrono>
#include <ctime>
#include <fmt/compile.h>

#include "../Domain.h"


FieldWriter::FieldWriter(FieldController &field_controller, real dt) :
        m_field_controller(field_controller), m_dt(dt) {
    std::string header = create_header();
    std::ofstream output_file(m_filename, std::ofstream::binary);
    output_file << header;
}

void FieldWriter::write_out(real t_cur) {
    std::string output;
    auto n = static_cast<int>(t_cur / m_dt);

    auto u = m_field_controller.get_field_u_data();
    auto v = m_field_controller.get_field_v_data();
    auto w = m_field_controller.get_field_w_data();
    auto p = m_field_controller.get_field_w_data();
    auto T = m_field_controller.get_field_T_data();
    auto C = m_field_controller.get_field_concentration_data();

    auto fields = {u, v, w, p, T, C};
    auto size = Domain::getInstance()->get_size();
    for (auto f: fields){
        for (size_t i = 0; i < size - 1; i++) {
            output.append(fmt::format(FMT_COMPILE("{};"), f[i]));
        }
        output.append(fmt::format(FMT_COMPILE("{}\n"), f[size - 1]));
    }
}

//================================= write header ===================================================
// *************************************************************************************************
/// \brief  write header for field storage. includes essential information: Nx, Ny, Nz
/// \details header format:
/// ###DOMAIN;<Nx>;<Ny>;<Nz>
/// ###FIELDS;u;v;w;p;T;concentration
/// ###DATE:<date>;XML:<XML>
/// <field temperature ambient>
// *************************************************************************************************
std::string FieldWriter::create_header() {
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();
    auto Nz = domain->get_Nz();

    std::string header = fmt::format(FMT_COMPILE("###DOMAIN;{};{};{}\n"), Nx, Ny, Nz);
    header.append(fmt::format(FMT_COMPILE("###FIELDS;u;v;w;p;T;concentration\n")));
    header.append(fmt::format(FMT_COMPILE("###DATE:{};XML:{}\n"),
                              std::ctime(&end_time), Parameters::getInstance()->get_filename()));

    auto size = domain->get_size();
    auto T_ambient = m_field_controller.get_field_T_ambient_data();
    for (size_t i = 0; i < size - 1; i++) {
        header.append(fmt::format(FMT_COMPILE("{};"), T_ambient[i]));
    }
    header.append(fmt::format(FMT_COMPILE("{}\n"), T_ambient[size - 1]));
    return header;
}
