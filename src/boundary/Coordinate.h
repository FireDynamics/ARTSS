/// \file       Coordinate.h
/// \brief
/// \date       Oct 26, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_COORDINATE_H_
#define ARTSS_COORDINATE_H_

#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>
#include "../utility/GlobalMacrosTypes.h"
#ifndef BENCHMARKING
//#include <fmt/core.h>
//#include <fmt/compile.h>
#endif

inline static const std::vector<std::string> axis_names = {"X", "Y", "Z"};

const size_t number_of_axes = 3;
enum CoordinateAxis : int {
    UNKNOWN_AXIS = -1,
    X = 0,
    Y = 1,
    Z = 2
};

template<class numeral>    // TODO(cvm) restrict to numerical values
class Coordinate {
public:
    Coordinate(numeral x, numeral y, numeral z) {
        m_coordinates = new numeral[number_of_axes];
        m_coordinates[CoordinateAxis::X] = x;
        m_coordinates[CoordinateAxis::Y] = y;
        m_coordinates[CoordinateAxis::Z] = z;
    };

    Coordinate() {
        m_coordinates = new numeral[number_of_axes];
        std::fill(m_coordinates, m_coordinates + number_of_axes, 0);
    };

    Coordinate(const Coordinate &original) {
        m_coordinates = new numeral[number_of_axes];
        for (size_t axis = 0; axis < number_of_axes; axis++) {
            m_coordinates[axis] = original.m_coordinates[axis];
        }
    }

    static std::string get_axis_name(CoordinateAxis axis) {
        return axis_names[axis];
    }

    ~Coordinate() { delete[] m_coordinates; }

    inline numeral &operator[](size_t i) const { return m_coordinates[i]; }

    Coordinate &operator+=(const Coordinate &rhs) {
        auto rhs_patches = rhs.m_coordinates;
        for (size_t i = 0; i < number_of_axes; ++i) {
            this->m_coordinates[i] += rhs_patches[i];
        }
        return *this;
    }

    Coordinate &operator+=(const numeral x) {
        for (size_t i = 0; i < number_of_axes; ++i) {
            this->m_coordinates[i] += x;
        }
        return *this;
    }

    Coordinate &operator*=(const Coordinate &rhs) {
        auto rhs_patches = rhs.m_coordinates;
        for (size_t i = 0; i < number_of_axes; ++i) {
            this->m_coordinates[i] *= rhs_patches[i];
        }
        return *this;
    }

    Coordinate &operator*=(const numeral x) {
        for (size_t i = 0; i < number_of_axes; ++i) {
            this->m_coordinates[i] *= x;
        }
        return *this;
    }


    void set_coordinate(numeral x, numeral y, numeral z) {
        m_coordinates[X] = x;
        m_coordinates[Y] = y,
                m_coordinates[Z] = z;
    }

    //TODO(cvm) restrict to numeral=size_t, partial specialisation?
    size_t get_index(size_t Nx, size_t Ny) { return IX(m_coordinates[X], m_coordinates[Y], m_coordinates[Z], Nx, Ny); }

    void copy(Coordinate<numeral> original) {
        for (size_t axis = 0; axis < number_of_axes; axis++) {
            m_coordinates[axis] = original.m_coordinates[axis];
        }
    }
    CoordinateAxis static match_axis(const std::string &string) {
        std::string upper_case;
        //std::transform(string.begin(), string.end(), upper_case.begin(), ::toupper);
        for (size_t an = 0; an < axis_names.size(); an++) {
            if (axis_names[an] == string) return (CoordinateAxis) an;
        }
        return UNKNOWN_AXIS;
    }

private:
    numeral *m_coordinates;
};

#ifndef BENCHMARKING
//template <typename T> struct fmt::formatter<Coordinate<T>> {
//    // Parses the format specifier, if needed (in my case, only return an iterator to the context)
//    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
//
//    // Actual formatting. The second parameter is the format specifier and the next parameters are the actual values from my custom type
//    template <typename FormatContext>
//    auto format(const Coordinate<T> &coord, FormatContext &ctx) -> decltype(ctx.out()) {
//        return format_to(
//                ctx.out(),
//                "({}|{}|{})", coord[0], coord[1], coord[2]);
//    }
//};
#endif

#endif /* ARTSS_COORDINATE_H_ */
