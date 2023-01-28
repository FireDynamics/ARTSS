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

inline static const std::vector<std::string> state_name = {"XML", "unmodified", "modified", "new", "deleted"};
constexpr size_t number_of_states = 5;
/// \brief 5 different states for obstacles, negligible if no data assimilation is used. State
///        always refers to the difference between the current ARTSS state and the config file read in.
/// \details
///  XML: state of an obstacle from the original XML. (replaces current obstacle with the data from
///       the XML, not meant for later usage)\n
///  UNMODIFIED: obstacle already exist in ARTSS, no changes necessary. Usage especially for
///              obstacles which were at least once modified/newly created\n
///  MODIFIED: obstacle is to be changed\n
///  NEW: obstacle doesn't exist in ARTSS and has to be newly created. This also counts for
///       obstacles which were defined in the XML, got deleted and now need to be newly created\n
///  DELETED: obstacle exists in ARTSS and has to be deleted\n
/// \example At the start of ARTSS, obstacles have the state XML. The first obstacle changes happen:
/// an obstacle is deleted (state: deleted) and another one changed (state: modified). With the
/// second obstacle change the deleted obstacle is to be restored. The deleted obstacle has to be
/// newly created (state:new) and and the previously changed obstacle gets the status unmodified
/// because it is not changed any further.
enum State : int {
    UNKNOWN_STATE = -1,
    XML = 0,
    UNMODIFIED = 1,
    MODIFIED = 2,
    NEW = 3,
    DELETED = 4
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
