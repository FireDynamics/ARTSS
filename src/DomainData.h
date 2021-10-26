/// \file       Domain.h
/// \brief      XML Domain parameters to variables
/// \date       Jul 16, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DOMAINDATA_H_
#define ARTSS_DOMAINDATA_H_

#include <cmath>
#ifdef _OPENACC
#include <accelmath.h>
#endif
#include <cstddef>

#include "utility/GlobalMacrosTypes.h"
#include "utility/Parameters.h"
#include "utility/Utility.h"

class DomainData {
 public:
    DomainData();

    static DomainData *getInstance();

    // getter
    size_t inline get_nx(size_t level = 0) const { return this->m_nx[level]; }
    size_t inline get_ny(size_t level = 0) const { return this->m_ny[level]; }
    size_t inline get_nz(size_t level = 0) const { return this->m_nz[level]; }

    size_t inline get_Nx(size_t level = 0) const { return get_nx(level) + 2; }
    size_t inline get_Ny(size_t level = 0) const { return get_ny(level) + 2; }
    size_t inline get_Nz(size_t level = 0) const { return get_nz(level) + 2; }

    real inline get_x1() const { return this->m_x1; }
    real inline get_x2() const { return this->m_x2; }
    real inline get_y1() const { return this->m_y1; }
    real inline get_y2() const { return this->m_y2; }
    real inline get_z1() const { return this->m_z1; }
    real inline get_z2() const { return this->m_z2; }

    real inline get_X1() const { return this->m_X1; }
    real inline get_X2() const { return this->m_X2; }
    real inline get_Y1() const { return this->m_Y1; }
    real inline get_Y2() const { return this->m_Y2; }
    real inline get_Z1() const { return this->m_Z1; }
    real inline get_Z2() const { return this->m_Z2; }

    real inline get_Lx() const { return fabs(m_X2 - m_X1); }
    real inline get_Ly() const { return fabs(m_Y2 - m_Y1); }
    real inline get_Lz() const { return fabs(m_Z2 - m_Z1); }
    real inline get_lx() const { return fabs(m_x2 - m_x1); }
    real inline get_ly() const { return fabs(m_y2 - m_y1); }
    real inline get_lz() const { return fabs(m_z2 - m_z1); }

    real inline get_dx(size_t level = 0) const { return this->get_lx() / static_cast<double>(m_nx[level]); }
    real inline get_dy(size_t level = 0) const { return this->get_ly() / static_cast<double>(m_ny[level]); }
    real inline get_dz(size_t level = 0) const { return this->get_lz() / static_cast<double>(m_nz[level]); }

    // start and end index of computational domain without ghost cells
    size_t inline get_index_x1(size_t level = 0) const { return static_cast<size_t> (std::round((m_x1 - m_X1) / get_dx(level))) + 1; }
    size_t inline get_index_x2(size_t level = 0) const { return static_cast<size_t> (std::round((m_x2 - m_X1) / get_dx(level))); }
    size_t inline get_index_y1(size_t level = 0) const { return static_cast<size_t> (std::round((m_y1 - m_Y1) / get_dy(level))) + 1; }
    size_t inline get_index_y2(size_t level = 0) const { return static_cast<size_t> (std::round((m_y2 - m_Y1) / get_dy(level))); }
    size_t inline get_index_z1(size_t level = 0) const { return static_cast<size_t> (std::round((m_z1 - m_Z1) / get_dz(level))) + 1; }
    size_t inline get_index_z2(size_t level = 0) const { return static_cast<size_t> (std::round((m_z2 - m_Z1) / get_dz(level))); }

    size_t inline get_levels() const { return m_levels; }

    size_t inline get_size(size_t level = 0) const { return get_Nx(level) * get_Ny(level) * get_Nz(level); }

    bool resize(long shift_x1, long shift_x2, long shift_y1, long shift_y2, long shift_z1, long shift_z2);

    void print();
    void printDetails();

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    static DomainData *single; // Singleton
    void calc_MG_values();

    static real calc_new_coord(real oldCoord, long shift, real cell_width);

    static bool set_new_value(long shift, real startCoord_p, real endCoord_p, real oldCoord, real cell_width, real *newCoord);

    size_t *m_nx, *m_ny, *m_nz;
    real m_x1, m_x2, m_y1, m_y2, m_z1, m_z2;
    real m_X1, m_X2, m_Y1, m_Y2, m_Z1, m_Z2;
    size_t m_levels = 0;

    void control();
};

#endif /* ARTSS_DOMAINDATA_H_ */
