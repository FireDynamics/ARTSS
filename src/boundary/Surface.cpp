//
// Created by linh on 01.10.19.
//

#include "Surface.h"


// TODO(issue 15): surface implementing
//  - underscores instead of camel case
//  - create file description
//  - remove surface cells from domain domain boundary cells
//  - develop a concept for boundary conditions
Surface::Surface(tinyxml2::XMLElement *element) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->debug("################ SURFACE ################");
#endif
    m_bc_values = new real[number_of_field_types];
    m_boundary_conditions = new BoundaryCondition[number_of_field_types];

    m_name = element->Attribute("name");
    m_patch = PatchObject::match_patch(element->Attribute("patch"));
#ifndef BENCHMARKING
    m_logger->debug("read surface '{}' for patch {}", m_name, PatchObject::get_patch_name(m_patch));
#endif

    auto cur_elem = element->FirstChildElement();
    while (cur_elem) {
        std::string node_name = cur_elem->Value();
        if (node_name == "boundary") {
            create_boundary(cur_elem);
        } else if (node_name == "geometry") {
            create_geometry(cur_elem);
        } else {
#ifndef BENCHMARKING
            m_logger->warn("Ignoring unknown node {}", node_name);
#endif
        }
        cur_elem = cur_elem->NextSiblingElement();
    }

    DomainData *domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();
    createSurface(Nx, Ny);
    print();
#ifndef BENCHMARKING
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::Surface(const std::string &name, Patch patch, Coordinate &start, Coordinate &end, size_t level) :
    m_patch(patch), m_name(name), m_start(start), m_end(end) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->info("################ SURFACE ################");
#endif

    DomainData *domain = DomainData::getInstance();
    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);

    init(Nx, Ny);
#ifndef BENCHMARKING
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::~Surface() {
    for (BoundaryData *bd: dataList) {
        delete (bd);
    }
    delete (m_surface_list);
}

void Surface::print() {
#ifndef BENCHMARKING
    m_logger->info("Surface name {}", m_name);
    m_logger->info("strides: X: {}, Y: {}, Z:{}",
                   get_stride_z(),
                   get_stride_y(),
                   get_stride_x());
    m_logger->info("size of Surface: {}", m_size_surfaceList);
    m_logger->info("coords: ({}|{}) ({}|{}) ({}|{})",
                   m_start[X], m_end[X],
                   m_start[Y], m_end[Y],
                   m_start[Z], m_end[Z]);
#endif
}

void Surface::init(size_t Nx, size_t Ny) {
    m_size_surfaceList = get_stride_x() * get_stride_y() * get_stride_z();
    m_surface_list = new size_t[m_size_surfaceList];

    createSurface(Nx, Ny);
    print();
}

void Surface::createSurface(size_t Nx, size_t Ny) {
    size_t counter = 0;
#ifndef BENCHMARKING
    m_logger->info("list size of sList: {}", m_size_surfaceList);
#endif
    // fill sList with corresponding indices
    for (size_t k = m_start[Z]; k < m_end[Z]; ++k) {
        for (size_t j = m_start[Y]; j < m_end[Y]; ++j) {
            for (size_t i = m_start[X]; i < m_end[X]; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                m_surface_list[counter++] = idx;
            }
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("control create Surface ({}|{})", counter,
                   m_size_surfaceList);
    m_logger->debug("end of creating sList");
#endif
}

void Surface::apply_boundary_conditions(Field &field, FieldType field_type, bool sync) {
    // m_bdc_boundary->apply_boundary_condition(dataField, indexFields, patch_starts, patch_ends, fieldType, level, sync);
    //TODO(issue 5)
}

size_t Surface::get_matching_index(real surface_coordinate, real spacing, real start_coordinate) {
    return static_cast<int>(round((-start_coordinate + surface_coordinate) / spacing));
}

void Surface::create_boundary(tinyxml2::XMLElement *element) {
    std::vector<std::string> field_strings = Utility::split(element->Attribute("field"), ',');
    BoundaryCondition boundary_condition = BoundaryData::match_boundary_condition(element->Attribute("type"));
    auto value = element->DoubleAttribute("value");

    for (const std::string &f: field_strings) {
        FieldType field_type = Field::match_field(f);
        m_boundary_conditions[field_type] = boundary_condition;
        m_bc_values[field_type] = value;
    }
}

void Surface::create_geometry(tinyxml2::XMLElement *cur_elem) {
    DomainData *domain_data = DomainData::getInstance();
    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real sx1, sx2;
    real sy1, sy2;
    real sz1, sz2;

    size_t i1, i2, j1, j2, k1, k2;
    if (m_patch == Patch::FRONT || m_patch == Patch::BACK) {
        sx1 = cur_elem->DoubleAttribute("sx1");
        sx2 = cur_elem->DoubleAttribute("sx2");
        i1 = get_matching_index(sx1, dx, X1) + 1;
        i2 = get_matching_index(sx2, dx, X1);

        sy1 = cur_elem->DoubleAttribute("sy1");
        sy2 = cur_elem->DoubleAttribute("sy2");
        j1 = get_matching_index(sy1, dy, Y1) + 1;
        j2 = get_matching_index(sy2, dy, Y1);

        if (m_patch == Patch::FRONT) {
            k1 = 0;
        } else {
            k1 = domain_data->get_Nz() - 1;
        }
        k2 = k1;

        m_start.set_coordinate(i1, j1, k1);
        m_end.set_coordinate(i2, j2, k2);
    } else if (m_patch == Patch::BOTTOM || m_patch == Patch::TOP) {
        sx1 = cur_elem->DoubleAttribute("sx1");
        sx2 = cur_elem->DoubleAttribute("sx2");
        i1 = get_matching_index(sx1, dx, X1) + 1;
        i2 = get_matching_index(sx2, dx, X1);

        sz1 = cur_elem->DoubleAttribute("sz1");
        sz2 = cur_elem->DoubleAttribute("sz2");
        k1 = get_matching_index(sz1, dz, Z1) + 1;
        k2 = get_matching_index(sz2, dz, Z1);

        if (m_patch == Patch::BOTTOM) {
            j1 = 0;
        } else {
            j1 = domain_data->get_Ny() - 1;
        }
        j2 = j1;

        m_start.set_coordinate(i1, j1, k1);
        m_end.set_coordinate(i2, j2, k2);
    } else if (m_patch == Patch::LEFT || m_patch == Patch::RIGHT) {
        sy1 = cur_elem->DoubleAttribute("sy1");
        sy2 = cur_elem->DoubleAttribute("sy2");
        j1 = get_matching_index(sy1, dy, Y1) + 1;
        j2 = get_matching_index(sy2, dy, Y1);

        sz1 = cur_elem->DoubleAttribute("sz1");
        sz2 = cur_elem->DoubleAttribute("sz2");
        k1 = get_matching_index(sz1, dz, Z1) + 1;
        k2 = get_matching_index(sz2, dz, Z1);

        if (m_patch == Patch::LEFT) {
            i1 = 0;
        } else {
            i1 = domain_data->get_Nx() - 1;
        }
        i2 = i1;

        m_start.set_coordinate(i1, j1, k1);
        m_end.set_coordinate(i2, j2, k2);
    }
}
