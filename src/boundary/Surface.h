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
    Surface(Settings::Settings const &settings,
            real x1, real x2,
            real y1, real y2,
            real z1, real z2,
            const std::string &name);
    Surface(Settings::Settings const &settings,
            Coordinate<size_t> &coords_start, Coordinate<size_t> &coords_end,
            size_t level,
            const std::string &name,
            Patch patch);
    ~Surface();

    size_t * get_surface_list() { return m_surface_list; }
    size_t get_size_surface_list() const { return m_size_surface_list; }

    Patch get_patch() { return m_patch; }


    std::string get_name() { return m_name; }

    void print();

    size_t get_stride(CoordinateAxis axis) { return m_end[axis] - m_start[axis] + 1; };
    size_t get_start_index(CoordinateAxis axis) { return m_start[axis]; }
    size_t get_end_index(CoordinateAxis axis) { return m_end[axis]; }

private:
    Settings::Settings const &m_settings;
    Patch m_patch;
    std::string m_name;

    Coordinate<size_t> m_start;
    Coordinate<size_t> m_end;

    size_t *m_surface_list;  // indices of surface
    size_t m_size_surface_list;

    void init();

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    size_t get_matching_index(real surface_coordinate, real spacing, real start_coordinate);

    size_t m_level;
};


#endif /* ARTSS_BOUNDARY_SURFACE_H_ */
