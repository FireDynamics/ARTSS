/// \file       Mapping.cpp
/// \brief      
/// \date       Dec 14, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "Mapping.h"

#include <algorithm>

namespace Mapping {

// *************************************************************************************************
/// \brief  matches string to axis_names
/// \param  string           string to be matched
// *************************************************************************************************
CoordinateAxis match_axis(std::string string) {
    std::transform(string.begin(), string.end(), string.begin(), ::toupper);
    for (size_t an = 0; an < axis_names.size(); an++) {
        if (axis_names[an] == string) return (CoordinateAxis) an;
    }
    return UNKNOWN_AXIS;
}

std::string get_axis_name(CoordinateAxis axis) {
    return axis_names[axis];
}

CoordinateAxis to_axis(Patch patch) {
    return CoordinateAxis(static_cast<size_t>(patch / 2));
}

// *************************************************************************************************
/// \brief  matches string to patch_names
/// \param  string           string to be matched
// *************************************************************************************************
Patch match_patch(const std::string &string) {
    for (size_t pn = 0; pn < patch_names.size(); pn++) {
        if (patch_names[pn] == string) return (Patch) pn;
    }
    return UNKNOWN_PATCH;
}

std::string get_patch_name(Patch p) {
    return patch_names[p];
}

/// \brief get patch of the given axis.
/// \details e.g. Axis X (0) start = false, results in Patch Right (1)
Patch to_patch(CoordinateAxis axis, bool start) {
    size_t p = axis * 2;
    if (!start) {
        p++;
    }
    return Patch(p);
}

//====================================== Matches ===================================================
// *************************************************************************************************
/// \brief  matches string to field_type_names
/// \param  string           string to be matched
// *************************************************************************************************
FieldType match_field(const std::string &string) {
    for (size_t fn = 0; fn < field_type_names.size(); fn++) {
        if (field_type_names[fn] == string) return (FieldType) fn;
    }
    return UNKNOWN_FIELD;
}

std::string get_field_type_name(FieldType f) {
    if (f != FieldType::UNKNOWN_FIELD) {
        return field_type_names[f];
    } else {
        return "UNKNOWN FIELD";
    }
}

// *************************************************************************************************
/// \brief  matches string to boundary_condition_names
/// \param  string           string to be matched
// *************************************************************************************************
BoundaryCondition match_boundary_condition(const std::string &string) {
    for (size_t tn = 0; tn < boundary_condition_names.size(); tn++) {
        if (boundary_condition_names[tn] == string) return static_cast<BoundaryCondition>(tn);
    }
    return UNKNOWN_CONDITION;
}

std::string get_boundary_condition_name(BoundaryCondition bc) {
    return boundary_condition_names[bc];
}

std::vector<CoordinateAxis> get_axes(Patch patch) {
    std::vector<CoordinateAxis> axes;
    axes.resize(2);
    if (patch == Patch::LEFT || patch == Patch::RIGHT) {
        axes[0] = Y;
        axes[1] = Z;
    } else if (patch == Patch::BOTTOM || patch == Patch::TOP) {
        axes[0] = X;
        axes[1] = Z;
    } else {  // p == Patch::LEFT || p == Patch::RIGHT
        axes[0] = X;
        axes[1] = Y;
    }
    return axes;
}

std::array<Patch, 2> get_patches(CoordinateAxis axis) {
    return {Patch(axis * 2), Patch(axis * 2 + 1)};
}
}  // end namespace Mapping
