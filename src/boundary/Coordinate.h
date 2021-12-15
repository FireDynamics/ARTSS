/// \file       Coordinate.h
/// \brief      
/// \date       Oct 26, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_COORDINATE_H_
#define ARTSS_COORDINATE_H_

#include <cstdlib>
#include <string>
#include <vector>


#include "../utility/GlobalMacrosTypes.h"

inline static const std::vector<std::string> axis_names = {"X", "Y", "Z"};

const size_t number_of_axis = 3;
enum CoordinateAxis : int {
    UNKNOWN_AXIS = -1,
    X = 0,
    Y = 1,
    Z = 2
};

class Coordinate {
  public:
    Coordinate(size_t i, size_t j, size_t k) {
        m_coordinates = new size_t[number_of_axis];
        m_coordinates[CoordinateAxis::X] = i;
        m_coordinates[CoordinateAxis::Y] = j;
        m_coordinates[CoordinateAxis::Z] = k;
    };
    Coordinate() {
        m_coordinates = new size_t[number_of_axis];
    };

    static std::string get_axis_name(size_t axis) {
        return axis_names[axis];
    }

    ~Coordinate() { delete[] m_coordinates; }
    inline size_t &operator[](size_t i) const { return m_coordinates[i]; }

    Coordinate &operator+=(const Coordinate &rhs) {
        auto rhs_patches = rhs.m_coordinates;
        for (size_t i = 0; i < number_of_axis; ++i) {
            this->m_coordinates[i] += rhs_patches[i];
        }
        return *this;
    }

    Coordinate &operator+=(const size_t x) {
        for (size_t i = 0; i < number_of_axis; ++i) {
            this->m_coordinates[i] += x;
        }
        return *this;
    }

    Coordinate &operator*=(const Coordinate &rhs) {
        auto rhs_patches = rhs.m_coordinates;
        for (size_t i = 0; i < number_of_axis; ++i) {
            this->m_coordinates[i] *= rhs_patches[i];
        }
        return *this;
    }
    Coordinate &operator*=(const size_t x) {
        for (size_t i = 0; i < number_of_axis; ++i) {
            this->m_coordinates[i] *= x;
        }
        return *this;
    }

    Coordinate(const Coordinate &original) {
        m_coordinates = new size_t[number_of_axis];
        for (size_t axis = 0; axis < number_of_axis; axis++) {
            m_coordinates[axis] = original.m_coordinates[axis];
        }
    }

    void set_coordinate(size_t i, size_t j, size_t k) { m_coordinates[X] = i; m_coordinates[Y] = j, m_coordinates[Z] = k; }
    size_t get_index(size_t Nx, size_t Ny) { return IX(m_coordinates[X], m_coordinates[Y], m_coordinates[Z], Nx, Ny); }
  private:
    size_t *m_coordinates;
};


#endif /* ARTSS_COORDINATE_H_ */
