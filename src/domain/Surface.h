/// \file       Surface.h
/// \brief      Data class of surface object
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_SURFACE_H_
#define ARTSS_BOUNDARY_SURFACE_H_

#include <vector>


#include "BoundaryData.h"
#include "BoundaryDataController.h"
#include "Coordinate.h"
#include "DomainData.h"
#include "../utility/Utility.h"
#include "../utility/Mapping.h"
#include "../utility/settings/Settings.h"

class Surface {
 public:
    Surface(real x1, real x2,
            real y1, real y2,
            real z1, real z2,
            const std::string &name,
            Patch patch);
    Surface(Coordinate<size_t> &coords_start, Coordinate<size_t> &coords_end,
            size_t level,
            const std::string &name,
            Patch patch);
    ~Surface() = default;

    const size_t * get_surface_list() { return m_surface_list.data(); }
    size_t get_size_surface_list() const { return m_size_surface_list; }

    Patch get_patch() { return m_patch; }
    size_t get_level() { return m_level; }
    std::string get_name() { return m_name; }

    void print();

    size_t get_stride(CoordinateAxis axis) { return m_strides[axis]; };
    [[nodiscard]] const Coordinate<size_t> & get_start_coordinates() const { return m_start; }
    [[nodiscard]] const Coordinate<size_t> & get_end_coordinates() const { return m_end; }
    size_t get_start_index(CoordinateAxis axis) { return m_start[axis]; }
    size_t get_end_index(CoordinateAxis axis) { return m_end[axis]; }

private:
    Patch m_patch;
    std::string m_name;

    Coordinate<size_t> m_start;
    Coordinate<size_t> m_end;
    Coordinate<size_t> m_strides;

    std::vector<size_t> m_surface_list;  // indices of surface
    size_t m_size_surface_list = 0;

    void init();

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    size_t get_matching_index(real surface_coordinate, real spacing, real start_coordinate);

    size_t m_level;
};


#endif /* ARTSS_BOUNDARY_SURFACE_H_ */
