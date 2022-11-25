/// \file       Mapping.h
/// \brief      
/// \date       Dec 14, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_MAPPING_H_
#define ARTSS_MAPPING_H_

#include <array>
#include <vector>
#include <string>

#include "../utility/GlobalMacrosTypes.h"

inline static const std::vector<std::string> state_name = {"unmodified", "modified", "new", "deleted"};
constexpr size_t number_of_states = 4;
enum State : int {
    UNKNOWN_STATE = -1,
    UNMODIFIED = 0,
    MODIFIED = 1,
    NEW = 2,
    DELETED = 3
};

inline static const std::vector<std::string> axis_names = {"X", "Y", "Z"};
constexpr size_t number_of_axes = 3;
enum CoordinateAxis : int {
    UNKNOWN_AXIS = -1,
    X = 0,
    Y = 1,
    Z = 2
};
inline static const std::vector<CoordinateAxis> all_coordinate_axes = {CoordinateAxis::X, CoordinateAxis::Y, CoordinateAxis::Z};

inline static const std::vector<std::string> patch_names = {"left", "right", "bottom", "top", "front", "back"};
constexpr size_t number_of_patches = 6;
enum Patch : int {
    UNKNOWN_PATCH = -1,
    LEFT = 0,
    RIGHT = 1,
    BOTTOM = 2,
    TOP = 3,
    FRONT = 4,
    BACK = 5
};
inline static const std::vector<Patch> all_patches = {Patch::LEFT, Patch::RIGHT, Patch::BOTTOM, Patch::TOP, Patch::FRONT, Patch::BACK};

inline static const std::vector<std::string> field_type_names = {"rho", "u", "v", "w", "p", "T", "nu"};
constexpr size_t number_of_field_types = 7;
enum FieldType : int {
    UNKNOWN_FIELD = -1, RHO = 0, U = 1, V = 2, W = 3, P = 4, T = 5, NU = 6
};

inline static const std::vector<std::string> boundary_condition_names = {"neumann", "dirichlet", "periodic"};
constexpr size_t number_of_boundary_conditions = 3;
enum BoundaryCondition : int {
    UNKNOWN_CONDITION = -1,
    NEUMANN = 0,
    DIRICHLET = 1,
    PERIODIC = 2
};

namespace Mapping {
    std::string get_axis_name(CoordinateAxis axis);
    //std::string get_axis_name(size_t axis) { return get_axis_name(CoordinateAxis(axis)); }
    CoordinateAxis match_axis(std::string string);
    CoordinateAxis to_axis(Patch patch);
    //CoordinateAxis to_axis(size_t patch) { return to_axis(Patch(patch)); }

    std::vector<CoordinateAxis> get_axes(Patch patch);
    std::string get_patch_name(Patch p);
    //std::string get_patch_name(size_t p) { return get_patch_name(Patch(p)); }
    Patch match_patch(const std::string &string);
    Patch to_patch(CoordinateAxis axis, bool start);
    //Patch to_patch(size_t axis, bool start) { return to_patch(CoordinateAxis(axis), start); }
    std::array<Patch, 2> get_patches(CoordinateAxis axis);

    std::string get_field_type_name(FieldType f);
    //std::string get_field_type_name(size_t f) { return get_field_type_name(FieldType(f)); };
    FieldType match_field(const std::string& string);

    std::string get_boundary_condition_name(BoundaryCondition bc);
    //std::string get_boundary_condition_name(size_t bc) { return get_boundary_condition_name(BoundaryCondition(bc)); };
    BoundaryCondition match_boundary_condition(const std::string &string);

    std::string get_state_name(State state);
    State match_state(const std::string &string);
}

#endif /* ARTSS_MAPPING_H_ */
