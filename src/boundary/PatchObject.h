/// \file       PatchObject.h
/// \brief      
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_PATCHOBJECT_H_
#define ARTSS_PATCHOBJECT_H_


#include <cstdlib>
#include "BoundaryData.h"

class PatchObject {
  public:
    PatchObject();
    ~PatchObject();
    size_t *m_patches;

    inline size_t &operator[](size_t i) const { return m_patches[i]; }

    PatchObject &operator+=(const PatchObject &rhs) {
        auto rhs_patches = rhs.m_patches;
        for (size_t i = 0; i < number_of_patches; ++i) {
            this->m_patches[i] += rhs_patches[i];
        }
        return *this;
    }
};


#endif /* ARTSS_PATCHOBJECT_H_ */
