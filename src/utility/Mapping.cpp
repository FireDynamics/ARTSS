/// \file       Mapping.cpp
/// \brief      
/// \date       Dec 14, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "Mapping.h"

#include <algorithm>

namespace Mapping {

    std::string get_axis_name(CoordinateAxis axis) {
        return axis_names[axis];
    }

    CoordinateAxis match_axis(const std::string &string) {
        std::string upper_case;
        std::transform(string.begin(), string.end(), upper_case.begin(), ::toupper);
        for (size_t an = 0; an < axis_names.size(); an++) {
            if (axis_names[an] == string) return (CoordinateAxis) an;
        }
        return UNKNOWN_AXIS;
    }

    CoordinateAxis to_axis(Patch patch) {
        return CoordinateAxis(static_cast<size_t>(patch / 2));
    }

    std::string get_patch_name(Patch p) {
        return patch_names[p];
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
}