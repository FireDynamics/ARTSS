/// \file       IAdaptionFunction.h
/// \brief      Interface for adaption methods
/// \date       December 04, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACE_IADAPTIONFUNCTION_H
#define ARTSS_INTERFACE_IADAPTIONFUNCTION_H

class IAdaptionFunction {
public:
    virtual void applyChanges() = 0;
    virtual bool update() = 0;
    virtual bool hasReduction() = 0;
};

#endif /* ARTSS_INTERFACE_IADAPTIONFUNCTION_H */
