/// \file       Cube.h
/// \brief      
/// \date       Sep 29, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_SOURCE_CUBE_H_
#define ARTSS_SOURCE_CUBE_H_


#include "../interfaces/ISourceFunction.h"

class Cube: public ISourceFunction {
 public:
    explicit Cube(const Settings::solver::sources::cube &cube);
    void update_source(Field &out, real t_cur) override;
    void update_obstacle_change() override { }

 private:
    const Settings::solver::sources::cube m_settings;
    Field m_source_field;
    void set_up();
};


#endif /* ARTSS_SOURCE_CUBE_H_ */
