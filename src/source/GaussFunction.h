/// \file       GaussFunction.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_GAUSSFUNCTION_H
#define ARTSS_GAUSSFUNCTION_H

#include <math.h>
#include <memory>

#include "../field/Field.h"
#include "../boundary/Obstacle.h"
#include "../interfaces/ISourceFunction.h"

class GaussFunction: public ISourceFunction {
 public:
    GaussFunction(real HRR, real cp,
            real x0, real y0, real z0,
            real sigma_x, real sigma_y, real sigma_z, real tau);

    GaussFunction(real HRR, real cp,
            real x0, real y0, real z0,
            real sigma_x, real sigma_y, real sigma_z, real tau,
            std::shared_ptr<spdlog::logger> logger);

    ~GaussFunction();

    void update_source(Field *out, real t_cur) override;

    bool test_obstacles_blocks(int level,
            int i0, int j0, int k0,
            int i, int j, int k,
            Obstacle** obst_list, int obst_id);
    bool test_obstacle_blocks(int i0, int j0, int k0,
            int i, int j, int k,
            const Obstacle &obst);

    void create_spatial_values(real HRR, real cp,
            real x0, real y0, real z0,
            real sigma_x, real sigma_y, real sigma_z);

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    Field *m_field_spatial_values;
    real m_tau;
    real get_time_value(real t_cur);
};


#endif /* ARTSS_GAUSSFUNCTION_H */
