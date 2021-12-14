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
#include "../utility/Mapping.h"


class PatchObject {
  public:
    PatchObject() = default;
    ~PatchObject() = default;

    inline size_t &operator[](size_t i) { return m_patches[i]; }  // r/w
    inline size_t const &operator[](size_t i) const { return m_patches[i]; }  // read only

    PatchObject &operator+=(const PatchObject &rhs) {
        auto rhs_patches = rhs.m_patches;
        for (size_t i = 0; i < number_of_patches; ++i) {
            this->m_patches[i] += rhs_patches[i];
        }
        return *this;
    }

    size_t get_sum() {
        size_t sum = 0;
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            sum += m_patches[patch];
        }
        return sum;
    }
    void add_value(size_t patch, size_t value) { m_patches[patch] += value; }

  private:
    std::array<size_t, number_of_patches> m_patches;
};


#endif /* ARTSS_PATCHOBJECT_H_ */
