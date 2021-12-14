/// \file       Coordinate.cpp
/// \brief      
/// \date       Dec 14, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include "Coordinate.h"

namespace Axis {
CoordinateAxis match_axis(const std::string &string) {
    std::string upper_case = string;
    std::transform(upper_case.begin(), upper_case.end(), upper_case.begin(), ::toupper);
    for (size_t an = 0; an < axis_names.size(); an++) {
        if (axis_names[an] == upper_case) return (CoordinateAxis) an;
    }
    return UNKNOWN_AXIS;
}

std::string get_axis_name(CoordinateAxis axis) {
    return axis_names[axis];
}

//CoordinateAxis to_axis(Patch patch) {
//    return CoordinateAxis(static_cast<size_t>(patch / 2));
//}
}