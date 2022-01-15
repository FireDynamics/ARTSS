/// \file       Coordinate.h
/// \brief
/// \date       Oct 26, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_DOMAIN_COORDINATE_H_
#define ARTSS_DOMAIN_COORDINATE_H_


#include <cstdlib>
#include <algorithm>
#include <string>
#include <stdexcept>

#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Mapping.h"

#ifndef BENCHMARKING
#include <fmt/format.h>
#endif

template<class numeral>
class Coordinate {
public:
    Coordinate(numeral numeral_x, numeral numeral_y, numeral numeral_z) {
        x = numeral_x;
        y = numeral_y;
        z = numeral_z;
    };

    Coordinate() = default;

    Coordinate(const Coordinate &original) {
        x = original.x;
        y = original.y;
        z = original.z;
    }

    inline numeral &operator[](size_t i) {
        if (i == CoordinateAxis::X) {
            return x;
        }
        if (i == CoordinateAxis::Y) {
            return y;
        }
        if (i == CoordinateAxis::Z) {
            return z;
        }
        throw std::runtime_error("unknown axis");
    }

    const inline numeral &operator[](size_t i) const {
        if (i == CoordinateAxis::X) {
            return x;
        }
        if (i == CoordinateAxis::Y) {
            return y;
        }
        if (i == CoordinateAxis::Z) {
            return z;
        }
        throw std::runtime_error("unknown axis");
    }

    Coordinate &operator+=(const Coordinate &rhs) {
        auto rhs_patches = rhs.m_coordinates;
        for (size_t i = 0; i < number_of_axes; ++i) {
            this->m_coordinates[i] += rhs_patches[i];
        }
        return *this;
    }

    // TODO
    //bool operator<=>(const Coordinate &rhs) const {
    //    return x <=> rhs.x && y <=> rhs.y && z <=> rhs.z;
    //}
    bool operator!=(const Coordinate &rhs) const {
        return x != rhs.x || y != rhs.y || z != rhs.z;
    }


    Coordinate &operator+=(const numeral n) {
        x += n;
        y += n;
        z += n;
        return *this;
    }

    Coordinate &operator*=(const Coordinate &rhs) {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        return *this;
    }

    Coordinate &operator*=(const numeral n) {
        for (size_t i = 0; i < number_of_axes; ++i) {
            this->m_coordinates[i] *= n;
        }
        return *this;
    }


    void set_coordinate(numeral numeral_x, numeral numeral_y, numeral numeral_z) {
        x = numeral_x;
        y = numeral_y,
        z = numeral_z;
    }

    size_t get_index(size_t Nx, size_t Ny) { return IX(x, y, z, Nx, Ny); }

    void copy(Coordinate<numeral> original) {
        x = original.x;
        y = original.y;
        z = original.z;
    }

private:
    numeral x = 0;
    numeral y = 0;
    numeral z = 0;
};

#ifndef BENCHMARKING
template <typename T> struct fmt::formatter<Coordinate<T>> {
    // Parses the format specifier, if needed (in my case, only return an iterator to the context)
    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

    // Actual formatting. The second parameter is the format specifier and the next parameters are the actual values from my custom type
    template <typename FormatContext>
    auto format(const Coordinate<T> &coord, FormatContext &ctx) -> decltype(ctx.out()) {
        return format_to(
                ctx.out(),
                "({}|{}|{})", coord[X], coord[Y], coord[Z]);
    }
};
#endif

#endif /* ARTSS_DOMAIN_COORDINATE_H_ */
