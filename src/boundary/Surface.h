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
    Surface(real x1, real x2,
            real y1, real y2,
            real z1, real z2,
            const std::string &name,
            Patch patch);
    Surface(const std::string &name,
            Patch patch,
            Coordinate<size_t> &start, Coordinate<size_t> &end,
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

    Coordinate<size_t> & get_start_coordinates() { return m_start; }
    Coordinate<size_t> & get_end_coordinates() { return m_end; }

private:
    size_t m_level = 0;
    size_t m_id;
    Patch m_patch;
    std::string m_name;

    Coordinate<size_t> m_start;
    Coordinate<size_t> m_end;

    size_t *m_surface_list;  // indices of surface
    size_t m_size_surfaceList;

    std::vector<BoundaryData*> dataList;

    void init(size_t Nx, size_t Ny);
    void createSurface(size_t Nx, size_t Ny);

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    size_t get_matching_index(real surface_coordinate, real spacing, real start_coordinate);
};


#endif /* ARTSS_BOUNDARY_SURFACE_H_ */
