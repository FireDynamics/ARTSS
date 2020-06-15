/// \file 		AdapationFunction.h
/// \brief 		Interface for Adapation method
/// \date 		December 04, 2018
/// \author 	My Linh Wuerzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACE_ADPATIONFUNCTIONI_H
#define ARTSS_INTERFACE_ADPATIONFUNCTIONI_H


class IAdaptionFunction {
public:
    virtual void applyChanges() = 0;
    virtual bool update() = 0;
    virtual bool hasReduction()=0;
private:
};


#endif /* ARTSS_INTERFACE_ADPATIONFUNCTIONI_H */
