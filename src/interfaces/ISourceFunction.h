/// \file       ISourceFunction.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_ISOURCEFUNCTION_H
#define ARTSS_ISOURCEFUNCTION_H

#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"

struct SourceFunctionNames{
    inline static const std::string RampTanh = "RampTanh";
};

class ISourceFunction {
public:
    virtual void update_source(Field *out, real t_cur) = 0;
};


#endif /* ARTSS_ISOURCEFUNCTION_H */
