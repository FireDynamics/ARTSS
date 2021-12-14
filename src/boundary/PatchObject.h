/// \file       PatchObject.h
/// \brief      
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_PATCHOBJECT_H_
#define ARTSS_PATCHOBJECT_H_

#include <cstdlib>
#include <string>
#include "Coordinate.h"


const size_t number_of_patches = 6;
enum Patch : int {
    UNKNOWN_PATCH = -1,
    LEFT = 0,
    RIGHT = 1,
    BOTTOM = 2,
    TOP = 3,
    FRONT = 4,
    BACK = 5
};

class PatchObject {
  public:
    PatchObject();
    ~PatchObject();

    inline size_t &operator[](size_t i) const { return m_patches[i]; }

    PatchObject &operator+=(const PatchObject &rhs) {
        auto rhs_patches = rhs.m_patches;
        for (size_t i = 0; i < number_of_patches; ++i) {
            this->m_patches[i] += rhs_patches[i];
        }
        return *this;
    }

    size_t get_sum();
    void add_value(size_t patch, size_t value) { m_patches[patch] += value; }

    static std::string get_patch_name(size_t p);
    static Patch match_patch(const std::string &string);
    /// \brief get patch of the given axis.
    /// \details e.g. Axis X (0) start = false, results in Patch Right (1)
    static Patch to_patch(CoordinateAxis axis, bool start) {
        size_t p = axis * 2;
        if (!start) {
            p++;
        }
        return Patch(p);
    }
  private:
    size_t *m_patches;
};


#endif /* ARTSS_PATCHOBJECT_H_ */
