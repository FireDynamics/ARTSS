/// \file       Analysis.h
/// \brief      Calculates residual, compares analytical and numerical solutions, saves variables
/// \date       July 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ANALYSIS_ANALYSIS_H_
#define ARTSS_ANALYSIS_ANALYSIS_H_

#include "Solution.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"
#include "../field/Field.h"
#include "../field/FieldController.h"

class Analysis {
 public:
    Analysis(const Settings::solver::solution &solution_settings, Solution &solution);

    void analyse(FieldController *solver, real t);

    void calc_L2_norm_mid_point(FieldController *solver, real t, real *sum);
    void calc_RMS_error(real sum_u, real sum_p, real sum_T);
    real calc_CFL(Field const &u, Field const &v, Field const &w, real dt) const;

    bool check_time_step_VN(real dt);

    void save_variables_in_file(FieldController *field_controller);

 private:
    const Settings::solver::solution &m_solution_settings;

    bool compare_solutions(read_ptr num, read_ptr ana, FieldType type, real t);

    real calc_absolute_spatial_error(read_ptr num, read_ptr ana);
    real calc_relative_spatial_error(read_ptr num, read_ptr ana);

    static void write_file(const Field &field, const std::string& filename);
    static void write_obstacles(const Field &field, const std::string &filename);

    Solution &m_solution;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_ANALYSIS_ANALYSIS_H_ */
