/// \file       DomainData.h
/// \brief      stores information about the domain
/// \date       Jul 16, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DOMAIN_DOMAINDATA_H_
#define ARTSS_DOMAIN_DOMAINDATA_H_

#include <cmath>
#ifdef _OPENACC
#include <accelmath.h>
#endif
#include <cstddef>

#include "Coordinate.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"

class DomainData {
 public:
    explicit DomainData(const Settings::physical_parameters &physical_params,
                        const Settings::domain_parameters &domain_params,
                        size_t multigrid_level);
    ~DomainData();
    static void init(const Settings::physical_parameters &physical_params,
                     const Settings::domain_parameters &domain_params,
                     size_t multigrid_level) {
        single = std::make_unique<DomainData>(physical_params, domain_params, multigrid_level);
    }
    static void reset() { single.reset(); }
    static DomainData *getInstance() { return single.get(); }

    // getter
    [[deprecated("Replaced by get_number_of_inner_cells")]]
    size_t inline get_nx(size_t level = 0) const { return get_number_of_inner_cells(CoordinateAxis::X, level); }
    [[deprecated("Replaced by get_number_of_inner_cells")]]
    size_t inline get_ny(size_t level = 0) const { return get_number_of_inner_cells(CoordinateAxis::Y, level); }
    [[deprecated("Replaced by get_number_of_inner_cells")]]
    size_t inline get_nz(size_t level = 0) const { return get_number_of_inner_cells(CoordinateAxis::Z, level); }
    size_t inline get_number_of_inner_cells(CoordinateAxis axis, size_t level = 0) const {
        return this->number_of_inner_cells[level][axis];
    }

    [[deprecated("Replaced by get_number_of_cells")]]
    size_t inline get_Nx(size_t level = 0) const { return get_number_of_cells(CoordinateAxis::X, level); }
    [[deprecated("Replaced by get_number_of_cells")]]
    size_t inline get_Ny(size_t level = 0) const { return get_number_of_cells(CoordinateAxis::Y, level); }
    [[deprecated("Replaced by get_number_of_cells")]]
    size_t inline get_Nz(size_t level = 0) const { return get_number_of_cells(CoordinateAxis::Z, level); }
    size_t inline get_number_of_cells(CoordinateAxis axis, size_t level = 0) const {
        return this->number_of_inner_cells[level][axis] + 2;
    }

    [[deprecated("Replaced by get_start_coord_CD")]]
    real inline get_x1() const { return get_start_coord_CD(CoordinateAxis::X); }
    [[deprecated("Replaced by get_start_coord_CD")]]
    real inline get_y1() const { return get_start_coord_CD(CoordinateAxis::Y); }
    [[deprecated("Replaced by get_start_coord_CD")]]
    real inline get_z1() const { return get_start_coord_CD(CoordinateAxis::Z); }
    real inline get_start_coord_CD(CoordinateAxis axis) const {
        return start_coords_CD[axis];
    }

    [[deprecated("Replaced by get_end_coord_CD")]]
    real inline get_x2() const { return get_end_coord_CD(CoordinateAxis::X); }
    [[deprecated("Replaced by get_end_coord_CD")]]
    real inline get_y2() const { return get_end_coord_CD(CoordinateAxis::Y); }
    [[deprecated("Replaced by get_end_coord_CD")]]
    real inline get_z2() const { return get_end_coord_CD(CoordinateAxis::Z); }
    real inline get_end_coord_CD(CoordinateAxis axis) const {
        return end_coords_CD[axis];
    }

    [[deprecated("Replaced by get_start_coord_PD")]]
    real inline get_X1() const { return get_start_coord_PD(CoordinateAxis::X); }
    [[deprecated("Replaced by get_start_coord_PD")]]
    real inline get_Y1() const { return get_start_coord_PD(CoordinateAxis::Y); }
    [[deprecated("Replaced by get_start_coord_PD")]]
    real inline get_Z1() const { return get_start_coord_PD(CoordinateAxis::Z); }
    real inline get_start_coord_PD(CoordinateAxis axis) const {
        return start_coords_PD[axis];
    }

    [[deprecated("Replaced by get_end_coord_PD")]]
    real inline get_X2() const { return get_end_coord_PD(CoordinateAxis::X); }
    [[deprecated("Replaced by get_end_coord_PD")]]
    real inline get_Y2() const { return get_end_coord_PD(CoordinateAxis::Y); }
    [[deprecated("Replaced by get_end_coord_PD")]]
    real inline get_Z2() const { return get_end_coord_PD(CoordinateAxis::Z); }
    real inline get_end_coord_PD(CoordinateAxis axis) const {
        return end_coords_PD[axis];
    }

    [[deprecated("Replaced by get_length_PD")]]
    real inline get_Lx() const { return get_length_PD(CoordinateAxis::X); }
    [[deprecated("Replaced by get_length_PD")]]
    real inline get_Ly() const { return get_length_PD(CoordinateAxis::Y); }
    [[deprecated("Replaced by get_length_PD")]]
    real inline get_Lz() const { return get_length_PD(CoordinateAxis::Z); }
    real inline get_length_PD(CoordinateAxis axis) const { return length_PD[axis]; }

    [[deprecated("Replaced by get_length_CD")]]
    real inline get_lx() const { return get_length_CD(CoordinateAxis::X); }
    [[deprecated("Replaced by get_length_CD")]]
    real inline get_ly() const { return get_length_CD(CoordinateAxis::Y); }
    [[deprecated("Replaced by get_length_CD")]]
    real inline get_lz() const { return get_length_CD(CoordinateAxis::Z); }
    real inline get_length_CD(CoordinateAxis axis) const {
        return fabs(end_coords_CD[axis] - start_coords_CD[axis]);
    }

    [[deprecated("Replaced by get_spacing")]]
    real inline get_dx(size_t level = 0) const { return get_spacing(CoordinateAxis::X, level); }
    [[deprecated("Replaced by get_spacing")]]
    real inline get_dy(size_t level = 0) const { return get_spacing(CoordinateAxis::Y, level); }
    [[deprecated("Replaced by get_spacing")]]
    real inline get_dz(size_t level = 0) const { return get_spacing(CoordinateAxis::Z, level); }
    real inline get_spacing(CoordinateAxis axis, size_t level = 0) const {
        return get_length_CD(CoordinateAxis(axis)) / static_cast<real>(get_number_of_inner_cells(CoordinateAxis(axis), level));
    }

    // start index of computational domain without ghost cells
    [[deprecated("Replaced by get_start_index_CD")]]
    size_t inline get_index_x1(size_t level = 0) const { return get_start_index_CD(CoordinateAxis::X, level); }
    [[deprecated("Replaced by get_start_index_CD")]]
    size_t inline get_index_y1(size_t level = 0) const { return get_start_index_CD(CoordinateAxis::Y, level); }
    [[deprecated("Replaced by get_start_index_CD")]]
    size_t inline get_index_z1(size_t level = 0) const { return get_start_index_CD(CoordinateAxis::Z, level); }
    size_t inline get_start_index_CD(CoordinateAxis axis, size_t level = 0) const {
        return static_cast<size_t> (std::round((start_coords_CD[axis] - start_coords_PD[axis]) / get_spacing(axis, level))) + 1;
    }

    // end index of computational domain without ghost cells
    [[deprecated("Replaced by get_end_index_CD")]]
    size_t inline get_index_x2(size_t level = 0) const { return get_end_index_CD(CoordinateAxis::X, level); }
    [[deprecated("Replaced by get_end_index_CD")]]
    size_t inline get_index_y2(size_t level = 0) const { return get_end_index_CD(CoordinateAxis::Y, level); }
    [[deprecated("Replaced by get_end_index_CD")]]
    size_t inline get_index_z2(size_t level = 0) const { return get_end_index_CD(CoordinateAxis::Z, level); }
    size_t inline get_end_index_CD(CoordinateAxis axis, size_t level = 0) const {
        return static_cast<size_t> (std::round((end_coords_CD[axis] - start_coords_PD[axis]) / get_spacing(axis, level)));
    }

    size_t inline get_levels() const { return m_levels; }

    size_t get_size(size_t level = 0) const;
    const Settings::physical_parameters& get_physical_parameters() const { return m_physical_parameters; }

    bool resize(const Coordinate<long>& shift_start, const Coordinate<long>& shift_end);

    void print();
    void print_details();

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    static std::unique_ptr<DomainData> single;  // Singleton
    void calc_MG_values();

    static real calc_new_coord(real oldCoord, long shift, real cell_width);

    static bool set_new_value(long shift, real start_coord_p, real end_coord_p, real old_coord, real cell_width, real *new_coord);

    // CD = computational domain
    Coordinate<size_t> *number_of_inner_cells;  // nx/ny/nz
    Coordinate<real> start_coords_CD;  // x1/y1/z1
    Coordinate<real> end_coords_CD;  // x2/y2/z2

    // PD = physical domain (non changing)
    Coordinate<real> length_PD;  // Lx/Ly/Lz
    Coordinate<real> start_coords_PD;  // X1/Y1/Z1
    Coordinate<real> end_coords_PD;  // X2/Y2/Z2

    size_t m_levels = 0;

    const struct Settings::physical_parameters m_physical_parameters;

    void control();
};

#endif /* ARTSS_DOMAIN_DOMAINDATA_H_ */
