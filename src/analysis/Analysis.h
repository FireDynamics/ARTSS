/// \file       Analysis.h
/// \brief      Calculates residual, compares analytical and numerical solutions, saves variables
/// \date       July 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ANALYSIS_ANALYSIS_H_
#define ARTSS_ANALYSIS_ANALYSIS_H_

#include "../interfaces/ISolver.h"
#include "../utility/GlobalMacrosTypes.h"
#include "Solution.h"

class Analysis {
public:
    explicit Analysis(Solution *solution);

    void analyse(ISolver *solver, real t);

    void calc_L2_norm_mid_point(ISolver *solver, real t, real *sum);
    void calc_RMS_error(real sum_u, real sum_p, real sum_T);

    bool check_time_step_VN(Field *u, real dt);
    bool check_time_step_CFL(Field *u, Field *v, Field *w, real dt);

    real set_DT_with_CFL(Field *u, Field *v, Field *w);

    void save_variables_in_file(ISolver *solv);

private:
    real m_tol = 1e-7;

    bool compare_solutions(read_ptr num, read_ptr ana, FieldType type, real t);

    real calc_absolute_spatial_error(read_ptr num, read_ptr ana);
    real calc_relative_spatial_error(read_ptr num, read_ptr ana);

    static void write_file(const real *field, const std::string& filename, size_t *inner_list, size_t size_inner_list, size_t *boundary_list, size_t size_boundary_list, size_t *obstacle_list,
                    size_t size_obstacle_list);

    bool has_analytic_solution = false;
    Solution *m_solution;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_ANALYSIS_ANALYSIS_H_ */
