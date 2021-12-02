/// \file       IAdaptionFunction.h
/// \brief      Interface for adaption methods
/// \date       December 04, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACE_IADAPTIONFUNCTION_H
#define ARTSS_INTERFACE_IADAPTIONFUNCTION_H

#include "../boundary/Coordinate.h"

class IAdaptionFunction {
public:
    virtual void apply_changes(Coordinate<long> *shift_start, Coordinate<long> *shift_end) = 0;
    virtual bool update(Coordinate<long> *shift_start, Coordinate<long> *shift_end) = 0;
    virtual bool has_reduction() = 0;
};

#endif /* ARTSS_INTERFACE_IADAPTIONFUNCTION_H */
