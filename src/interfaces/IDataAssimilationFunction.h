/// \file      IDataAssimilationFunction.h
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_SRC_INTERFACES_IDATAASSIMILATIONFUNCTION_H
#define ARTSS_SRC_INTERFACES_IDATAASSIMILATIONFUNCTION_H


class IDataAssimilationFunction {
 public:
    virtual bool control(Field *u, Field *v, Field *w, Field *p, Field *T, Field *C) = 0;
    virtual void assimilate(Field *u, Field *v, Field *w, Field *p, Field *T, Field *C) = 0;
};
#endif /* ARTSS_SRC_INTERFACES_IDATAASSIMILATIONFUNCTION_H */
