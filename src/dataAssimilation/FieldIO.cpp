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
#include "../utility/Mapping.h"


FieldIO::FieldIO(const std::string &xml_file_name, const std::string &output_file_name) :
        m_file_name(output_file_name),
        m_logger(Utility::create_logger(typeid(this).name())) {

    auto domain_data = DomainData::getInstance();
    real t_end = domain_data->get_physical_parameters().t_end;

    size_t n = static_cast<size_t>(std::round(t_end / domain_data->get_physical_parameters().dt)) + 1;
    m_positions = new long[n];

    HighFive::File out_file(m_file_name, HighFive::File::ReadWrite | HighFive::File::Create);
    create_header(out_file, xml_file_name);
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
    HighFive::File out_file(m_file_name, HighFive::File::ReadWrite);
    auto group_name = fmt::format("{:.5e}", t_cur);
    m_logger->error("attempt to write @t:{}", t_cur);

    if (out_file.exist(group_name)) {
        out_file.unlink(group_name);
    }
    
    HighFive::Group t_group = out_file.createGroup(group_name);
    // Field fields[] = {u, v, w, p, T, C};
    Field fields[] = {u, v, w, T, C};
    size_t size = u.get_size();
    std::vector<size_t> dims{1, size};
    for (Field &f: fields) {
        auto field_name = Mapping::get_field_type_name(f.get_type());
        m_logger->error("attempt to write @t:{}:{}", t_cur, field_name);
        t_group.createDataSet<real>(field_name, HighFive::DataSpace(dims));
    }
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
    HighFive::File input_file(m_file_name, HighFive::File::ReadOnly);

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
void FieldIO::create_header(HighFive::File &file, const std::string &xml_file_name) {
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    if (file.exist("metadata")) {
        file.unlink("metadata");
    }
    HighFive::Group meta_group = file.createGroup("metadata");

    auto domain_data = DomainData::getInstance();

    auto dt = domain_data->get_physical_parameters().dt;
    HighFive::DataSet dt_set = meta_group.createDataSet<real>("dt", HighFive::DataSpace::From(dt));
    dt_set.write(dt);

    size_t Nx = domain_data->get_number_of_cells(CoordinateAxis::X);
    size_t Ny = domain_data->get_number_of_cells(CoordinateAxis::Y);
    size_t Nz = domain_data->get_number_of_cells(CoordinateAxis::Z);

    size_t domain[3] = {Nx, Ny, Nz};
    HighFive::DataSet domain_set = meta_group.createDataSet<size_t>("domain", HighFive::DataSpace::From(domain));
    domain_set.write(domain);

    std::vector<std::string> fields{"u" , "v", "w", "p", "T", "concentration"};
    HighFive::DataSet fields_set = meta_group.createDataSet<std::string>("fields", HighFive::DataSpace::From(fields));
    fields_set.write(fields);

    std::string date[1] = {std::ctime(&end_time)};
    HighFive::DataSet date_set = meta_group.createDataSet<std::string>("date", HighFive::DataSpace::From(date));
    date_set.write(date);

    std::string xml[1] = {xml_file_name};
    HighFive::DataSet xml_set = meta_group.createDataSet<std::string>("xml", HighFive::DataSpace::From(date));
    xml_set.write(xml);
}

void FieldIO::read_field(HighFive::File &file, Field &field) {
    auto field_name = Mapping::get_field_type_name(field.get_type());
    auto ds = file.getDataSet(field_name);
    ds.read(field.data);
}

void FieldIO::read_fields(const Settings::data_assimilation::field_changes &field_changes,
                          Field &u, Field &v, Field &w,
                          Field &p, Field &T, Field &C) {
    int n;
    std::string line;
    if (!field_changes.changed) {
        return;
    }

    // no changes -> read original file
    m_logger->info("opening dat file {}", field_changes.file_name);
    try {
        HighFive::File file_changes(field_changes.file_name, HighFive::File::ReadOnly);

        if (field_changes.u_changed) {
            read_field(file_changes, u);
            m_logger->debug("read changed u Field");
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
    } catch (const std::exception &ex) {
        m_logger->warn(fmt::format("File '{}' could not be opened. No changes will be applied.", field_changes.file_name));
        m_logger->warn(fmt::format("Exception during reading {}", ex.what()));
    }
}

void FieldIO::read_fields(const real t_cur,
                          const Settings::data_assimilation::field_changes &field_changes,
                          Field &u, Field &v, Field &w,
                          Field &p, Field &T, Field &C) {
    m_logger->debug("read time step {}", t_cur);
    if (field_changes.changed) {  // no changes -> read original file
        m_logger->debug("no field changes");
        read_fields(t_cur, u, v, w, p, T, C);
        return;
    }

    try {
        HighFive::File org_file(m_file_name, HighFive::File::ReadOnly);
        HighFive::File new_file(field_changes.file_name, HighFive::File::ReadOnly);

        if (field_changes.u_changed) {
            read_field(new_file, u);
            m_logger->debug("read changed u Field");
        } else {
            read_field(org_file, u);
        }
        if (field_changes.v_changed) {
            read_field(new_file, v);
            m_logger->debug("read changed v Field");
        } else {
            read_field(org_file, v);
        }
        if (field_changes.w_changed) {
            read_field(new_file, w);
            m_logger->debug("read changed w Field");
        } else {
            read_field(org_file, w);
        }
        if (field_changes.p_changed) {
            read_field(new_file, p);
            m_logger->debug("read changed p Field");
        } else {
            read_field(org_file, p);
        }
        if (field_changes.T_changed) {
            read_field(new_file, T);
            m_logger->debug("read changed T Field");
        } else {
            read_field(org_file, T);
        }
        if (field_changes.C_changed) {
            read_field(new_file, C);
            m_logger->debug("read changed C Field");
        } else {
            read_field(org_file, C);
        }
    } catch (const std::exception &ex) {
        m_logger->warn(fmt::format("File '{}' could not be opened, original data will be loaded", field_changes.file_name));
        m_logger->warn(fmt::format("Exception during reading {}", ex.what()));
        read_fields(t_cur, u, v, w, p, T, C);
    }
}
