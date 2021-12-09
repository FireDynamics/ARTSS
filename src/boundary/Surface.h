//
// Created by linh on 01.10.19.
//

#ifndef ARTSS_BOUNDARY_SURFACE_H_
#define ARTSS_BOUNDARY_SURFACE_H_

#include <vector>
#include "BoundaryData.h"

#include "../DomainData.h"
#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"

#include "BoundaryDataController.h"
#include "Coordinate.h"

class Surface {
 public:
    Surface(Settings::Settings const &settings, Settings::SurfaceSetting const &surface_setting);
    Surface(Settings::Settings const &settings,
            const std::string &name,
            Patch patch,
            Coordinate &start, Coordinate &end,
            size_t level);
    ~Surface();

    size_t * get_surface_list() { return m_surface_list; }
    size_t get_size_surface_list() const { return m_size_surfaceList; }

    Patch get_patch() { return m_patch; }

    size_t get_stride_x() { return m_end[X] - m_start[X] + 1;}
    size_t get_stride_y() { return m_end[Y] - m_start[Y] + 1;}
    size_t get_stride_z() { return m_end[Z] - m_start[Z] + 1;}

    std::string get_name() { return m_name; }
    size_t get_id() { return m_id; }
    void set_id(size_t id) { m_id = id; }

    void apply_boundary_conditions(Field &field, bool sync);

    void print();

    Coordinate & get_start_coordinates() { return m_start; }
    Coordinate & get_end_coordinates() { return m_end; }

private:
    Settings::Settings const &m_settings;
    size_t m_level = 0;
    size_t m_id;
    Patch m_patch;
    std::string m_name;

    Coordinate m_start;
    Coordinate m_end;

    size_t *m_surface_list;  // indices of surface
    size_t m_size_surfaceList;

    std::vector<BoundaryData*> dataList;

    void init(size_t Nx, size_t Ny);
    void createSurface(size_t Nx, size_t Ny);

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    size_t get_matching_index(real surface_coordinate, real spacing, real start_coordinate);

    void create_boundary(tinyxml2::XMLElement *element);

    void create_geometry(tinyxml2::XMLElement *element);
};


#endif /* ARTSS_BOUNDARY_SURFACE_H_ */
