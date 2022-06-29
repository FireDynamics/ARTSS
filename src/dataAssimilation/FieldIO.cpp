/// \file       FieldIO.cpp
/// \brief      Class for reading/writing the raw data of the fields u, v, w, p, T, C
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.

#include "FieldIO.h"

#include <ctime>
#include <chrono>
#include <iomanip>

#include <fmt/compile.h>
#include <fstream>

#include "../domain/DomainData.h"
#include "../utility/Mapping.h"


FieldIO::FieldIO(const std::string &xml_file_name, const std::string &output_dir) :
        m_path(output_dir),
        m_xml_filename(xml_file_name),
        m_logger(Utility::create_logger(typeid(this).name())) {
    namespace fs = std::filesystem;
    fs::create_directories(m_path);
    m_meta_path = m_path;
    m_meta_path /= "meta";

    create_meta_file(0);
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
void FieldIO::write_fields(real t_current, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C) {
    auto tstr = fmt::format("{:.5e}", t_current);
    std::filesystem::path file_name = m_path;
    file_name /= tstr;

    HighFive::File out_file(file_name, HighFive::File::ReadWrite | HighFive::File::Create);
    m_logger->debug("attempt to write @t:{}", t_current);

    create_header(out_file);

    if (out_file.exist(tstr)) {
        out_file.unlink(tstr);
    }

    HighFive::Group t_group = out_file.createGroup(tstr);
    Field fields[] = {u, v, w, p, T, C};
    size_t size = u.get_size();
    std::vector<size_t> dims{1, size};
    for (Field &f: fields) {
        auto field_name = Mapping::get_field_type_name(f.get_type());
        m_logger->debug("attempt to write @t:{}:{}", t_current, field_name);
        HighFive::DataSet dsf = t_group.createDataSet<real>(field_name, HighFive::DataSpace(dims));
        dsf.write(f.get_data());
        m_logger->debug("sum of {}: {}", Mapping::get_field_type_name(f.get_type()), f.get_sum());
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
    m_logger->debug("read original data @t: {}", t_cur);
    auto t_cur_str = fmt::format("{:.5e}", t_cur);
    std::filesystem::path file_name = m_path;
    file_name /= t_cur_str;
    HighFive::File input_file(file_name, HighFive::File::ReadOnly);

    read_vis_field(input_file, u, t_cur);
    read_vis_field(input_file, v, t_cur);
    read_vis_field(input_file, w, t_cur);
    read_vis_field(input_file, p, t_cur);
    read_vis_field(input_file, T, t_cur);
    read_vis_field(input_file, C, t_cur);
}

// ================================= write header ==================================================
// *************************************************************************************************
/// \brief  write header for field storage. includes essential debugger notation: xml file name,
///         which fields are written, date, Nx, Ny, Nz,
// *************************************************************************************************
void FieldIO::create_header(HighFive::File &file) {
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

    std::string xml[1] = {m_xml_filename};
    HighFive::DataSet xml_set = meta_group.createDataSet<std::string>("xml", HighFive::DataSpace::From(date));
    xml_set.write(xml);

    real t_cur[1] = {0};
    HighFive::DataSet t_cur_set = meta_group.createDataSet<real>("t_cur", HighFive::DataSpace::From(t_cur));
    t_cur_set.write(t_cur);
}

void FieldIO::read_field(HighFive::File &file, Field &field) {
    m_logger->debug("read changed field data of {} in {}", Mapping::get_field_type_name(field.get_type()), file.getName());
    auto field_name = Mapping::get_field_type_name(field.get_type());
    auto ds = file.getDataSet(field_name);
    ds.read(field.data);
}

void FieldIO::read_vis_field(HighFive::File &file, Field &field, const real t) {
    auto field_name = Mapping::get_field_type_name(field.get_type());
    auto full_name = fmt::format("/{:.5e}/{}", t, field_name);
    m_logger->info("opening vis file at {}", full_name);
    auto ds = file.getDataSet(full_name);
    ds.read(field.data);

    m_logger->debug("read @t {} sum of {}: {}", t, Mapping::get_field_type_name(field.get_type()), field.get_sum());
}

void FieldIO::read_changed_fields(const Settings::data_assimilation::field_changes &field_changes,
                                  Field &u, Field &v, Field &w,
                                  Field &p, Field &T, Field &C) {
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
    create_meta_file(t_cur);
    if (!field_changes.changed) {  // no changes -> read original file
        m_logger->debug("no field changes");
        read_fields(t_cur, u, v, w, p, T, C);
        return;
    }

    try {
        auto t_cur_str = fmt::format("{:.5e}", t_cur);
        std::filesystem::path file_name = m_path;
        file_name /= t_cur_str;
        HighFive::File org_file(file_name, HighFive::File::ReadOnly);
        HighFive::File new_file(field_changes.file_name, HighFive::File::ReadOnly);

        if (field_changes.u_changed) {
            read_field(new_file, u);
            m_logger->debug("read changed u Field");
        } else {
            read_vis_field(org_file, u, t_cur);
        }
        if (field_changes.v_changed) {
            read_field(new_file, v);
            m_logger->debug("read changed v Field");
        } else {
            read_vis_field(org_file, v, t_cur);
        }
        if (field_changes.w_changed) {
            read_field(new_file, w);
            m_logger->debug("read changed w Field");
        } else {
            read_vis_field(org_file, w, t_cur);
        }
        if (field_changes.p_changed) {
            read_field(new_file, p);
            m_logger->debug("read changed p Field");
        } else {
            read_vis_field(org_file, p, t_cur);
        }
        if (field_changes.T_changed) {
            read_field(new_file, T);
            m_logger->debug("read changed T Field");
        } else {
            read_vis_field(org_file, T, t_cur);
        }
        if (field_changes.C_changed) {
            read_field(new_file, C);
            m_logger->debug("read changed C Field");
        } else {
            read_vis_field(org_file, C, t_cur);
        }
    } catch (const std::exception &ex) {
        m_logger->warn(fmt::format("File '{}' could not be opened, original data will be loaded", field_changes.file_name));
        m_logger->warn(fmt::format("Exception during reading {}", ex.what()));
        read_fields(t_cur, u, v, w, p, T, C);
    }
}

void FieldIO::create_meta_file(real t_cur) {
    std::ofstream meta_file;
    meta_file.open(m_meta_path);
    meta_file << fmt::format("t:{}\n", t_cur);
    meta_file << fmt::format("xml_name:{}\n", m_xml_filename);
    meta_file.close();
}
