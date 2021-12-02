/// \file       IAdaptionFunction.h
/// \brief      Interface for adaption methods
/// \date       December 04, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACE_IADAPTIONFUNCTION_H
#define ARTSS_INTERFACE_IADAPTIONFUNCTION_H

class IAdaptionFunction {
public:
    virtual void apply_changes(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) = 0;
    virtual bool update(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) = 0;
    virtual bool has_reduction() = 0;
};

#endif /* ARTSS_INTERFACE_IADAPTIONFUNCTION_H */
