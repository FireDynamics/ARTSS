/// \file       Domain.cpp
/// \brief      XML Domain parameters to variables
/// \date       July 16, 2018
/// \author   My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DOMAIN_H
#define ARTSS_DOMAIN_H


#include <cstddef>
#include <cstdint>
#include <cmath>
#include "utility/GlobalMacrosTypes.h"

class Domain {
 public:
    Domain();

    static Domain* getInstance();

    // getter
    size_t inline Getnx() {return this->m_nx[0];}
    size_t inline Getny() {return this->m_ny[0];}
    size_t inline Getnz() {return this->m_nz[0];}
    size_t inline Getnx(size_t level) {return this->m_nx[level];}
    size_t inline Getny(size_t level) {return this->m_ny[level];}
    size_t inline Getnz(size_t level) {return this->m_nz[level];}

    size_t inline GetNx() {
        return static_cast<size_t> (std::round(GetLx() / Getdx() + 2));
    }
    size_t inline GetNy() {
        return static_cast<size_t> (std::round(GetLy() / Getdy() + 2));
    }
    size_t inline GetNz() {
        return static_cast<size_t> (std::round(GetLz() / Getdz() + 2));
    }
    size_t inline GetNx(size_t level) {
        return static_cast<size_t> (std::round(GetLx() / Getdx(level) + 2));
    }
    size_t inline GetNy(size_t level) {
        return static_cast<size_t> (std::round(GetLy() / Getdy(level) + 2));
    }
    size_t inline GetNz(size_t level) {
        return static_cast<size_t> (std::round(GetLz() / Getdz(level) + 2));
    }

    real inline Getx1() {return this->m_x1;}
    real inline Getx2() {return this->m_x2;}
    real inline Gety1() {return this->m_y1;}
    real inline Gety2() {return this->m_y2;}
    real inline Getz1() {return this->m_z1;}
    real inline Getz2() {return this->m_z2;}

    real inline GetX1() {return this->m_X1;}
    real inline GetX2() {return this->m_X2;}
    real inline GetY1() {return this->m_Y1;}
    real inline GetY2() {return this->m_Y2;}
    real inline GetZ1() {return this->m_Z1;}
    real inline GetZ2() {return this->m_Z2;}

    real inline GetLx() {return fabs(m_X2-m_X1);}
    real inline GetLy() {return fabs(m_Y2-m_Y1);}
    real inline GetLz() {return fabs(m_Z2-m_Z1);}
    real inline Getlx() {return fabs(m_x2-m_x1);}
    real inline Getly() {return fabs(m_y2-m_y1);}
    real inline Getlz() {return fabs(m_z2-m_z1);}

    real inline Getdx() {
        return this->Getlx()/static_cast<real>(m_nx[0]-2);
    }
    real inline Getdy() {
        return this->Getly()/static_cast<real>(m_ny[0]-2);
    }
    real inline Getdz() {
        return this->Getlz()/static_cast<real>(m_nz[0]-2);
    }
    real inline Getdx(size_t level) {
        return this->Getlx()/static_cast<real>(m_nx[level]-2);
    }
    real inline Getdy(size_t level) {
        return this->Getly()/static_cast<real>(m_ny[level]-2);
    }
    real inline Getdz(size_t level) {
        return this->Getlz()/static_cast<real>(m_nz[level]-2);
    }

    // start and end index of computational domain without ghost cells
    size_t inline GetIndexx1() { return GetIndexx1(0); }
    size_t inline GetIndexx2() { return GetIndexx2(0); }
    size_t inline GetIndexy1() { return GetIndexy1(0); }
    size_t inline GetIndexy2() { return GetIndexy2(0); }
    size_t inline GetIndexz1() { return GetIndexz1(0); }
    size_t inline GetIndexz2() { return GetIndexz2(0); }
    size_t inline GetIndexx1(size_t level) {
        return static_cast<size_t> (std::round((m_x1 - m_X1) / Getdx(level)))+1;
    }
    size_t inline GetIndexx2(size_t level) {
        return static_cast<size_t> (std::round((m_x2 - m_X1) / Getdx(level)));
    }
    size_t inline GetIndexy1(size_t level) {
        return static_cast<size_t> (std::round((m_y1 - m_Y1) / Getdy(level)))+1;
    }
    size_t inline GetIndexy2(size_t level) {
        return static_cast<size_t> (std::round((m_y2 - m_Y1) / Getdy(level)));
    }
    size_t inline GetIndexz1(size_t level) {
        return static_cast<size_t> (std::round((m_z1 - m_Z1) / Getdz(level)))+1;
    }
    size_t inline GetIndexz2(size_t level) {
        return static_cast<size_t> (std::round((m_z2 - m_Z1) / Getdz(level)));
    }

    size_t inline GetLevels() {return m_levels;}

    size_t inline GetSize() {
        return GetNx()*GetNy()*GetNz();
    }
    size_t inline GetSize(size_t level) {
        return GetNx(level)*GetNy(level)*GetNz(level);
    }

    bool Resize(
            int64_t shift_x1, int64_t shift_x2,
            int64_t shift_y1, int64_t shift_y2,
            int64_t shift_z1, int64_t shift_z2);

    void print();
    void printDetails();

 private:
    static Domain* single;  // Singleton
    void calcMGValues();
    static real calcNewCoord(real oldCoord , int64_t shift, real cellwidth);
    static bool setNewValue(
            int64_t shift, real startCoord_p, real endCoord_p,
            real oldCoord, real cellwidth, real *newCoord);

    size_t *m_nx, *m_ny, *m_nz;
    real m_x1, m_x2, m_y1, m_y2, m_z1, m_z2;
    real m_X1, m_X2, m_Y1, m_Y2, m_Z1, m_Z2;
    size_t m_levels = 0;
};

#endif //ARTSS_DOMAIN_H
