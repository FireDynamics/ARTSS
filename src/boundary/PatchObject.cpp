/// \file       PatchObject.cpp
/// \brief      
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include "PatchObject.h"
#include <vector>

inline static const std::vector<std::string> patch_names = {"front", "back", "bottom", "top", "left", "right"};

PatchObject::PatchObject() {
    m_patches = new size_t[number_of_patches];
    std::fill(m_patches, m_patches + number_of_patches, 0);
}

PatchObject::~PatchObject() {
    delete[] m_patches;
}

std::string PatchObject::get_patch_name(size_t p) {
    return patch_names[p];
}

// *************************************************************************************************
/// \brief  matches string to patch_names
/// \param  string           string to be matched
// *************************************************************************************************
Patch PatchObject::match_patch(const std::string &string) {
    for (size_t pn = 0; pn < patch_names.size(); pn++) {
        if (patch_names[pn] == string) return (Patch) pn;
    }
    return UNKNOWN_PATCH;
}

size_t PatchObject::get_sum() {
    size_t sum = 0;
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        sum += m_patches[patch];
    }
    return sum;
}
