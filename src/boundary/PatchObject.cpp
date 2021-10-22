/// \file       PatchObject.cpp
/// \brief      
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "PatchObject.h"

PatchObject::PatchObject() {
    m_patches = new size_t[number_of_patches];
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_patches[patch] = 0;
    }
}
