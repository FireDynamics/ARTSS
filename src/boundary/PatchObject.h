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

const size_t number_of_patches = 6;
enum Patch : int {
    UNKNOWN_PATCH = -1,
    FRONT = 0,
    BACK = 1,
    BOTTOM = 2,
    TOP = 3,
    LEFT = 4,
    RIGHT = 5
};

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

    static std::string get_patch_name(size_t p);
    static Patch match_patch(const std::string &string);
};


#endif /* ARTSS_PATCHOBJECT_H_ */
