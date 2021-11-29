/// \file       Solution.cpp
/// \brief      Calculates analytical solution
/// \details    This class calculates the analytical solution of different test cases
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Solution.h"
#include "../Functions.h"


Solution::Solution(Settings::Settings const &settings, const std::string &initial_condition, bool has_analytical_solution) :
        m_settings(settings),
        m_u_analytical_solution(Field(FieldType::U)),
        m_v_analytical_solution(Field(FieldType::V)),
        m_w_analytical_solution(Field(FieldType::W)),
        m_p_analytical_solution(Field(FieldType::P)),
        m_T_analytical_solution(Field(FieldType::T)),
        m_has_analytical_solution(has_analytical_solution) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(m_settings, typeid(this).name());
#endif

    // set function pointer to chosen initial condition
    if (initial_condition == FunctionNames::gauss_bubble) {
        m_init_function = &Solution::gauss_bubble;
    } else if (initial_condition == FunctionNames::exp_sinus_prod) {
        m_init_function = &Solution::exp_sinus_prod;
    } else if (initial_condition == FunctionNames::exp_sinus_sum) {
        m_init_function = &Solution::exp_sinus_sum;
    } else if (initial_condition == FunctionNames::hat) {
        m_init_function = &Solution::hat;
    } else if (initial_condition == FunctionNames::sin_sin_sin) {
        m_init_function = &Solution::sin_sin_sin;
    } else if (initial_condition == FunctionNames::mcdermott) {
        m_init_function = &Solution::mcDermott;
    } else if (initial_condition == FunctionNames::vortex) {
        m_init_function = &Solution::vortex;
    } else if (initial_condition == FunctionNames::vortex_y) {
        m_init_function = &Solution::vortex_y;
    } else if (initial_condition == FunctionNames::beltrami) {
        m_init_function = &Solution::beltrami;
    } else if (initial_condition == FunctionNames::buoyancy_mms) {
        m_init_function = &Solution::buoyancy_mms;
    } else {
#ifndef BENCHMARKING
        m_logger->info("Analytical solution set to zero!");
#endif
        m_init_function = &Solution::zero;

        m_u_analytical_solution.set_value(0);
        m_v_analytical_solution.set_value(0);
        m_w_analytical_solution.set_value(0);
        m_p_analytical_solution.set_value(0);
        m_T_analytical_solution.set_value(0);
    }
}

void Solution::gauss_bubble(const real t) {
    // Advection test case
    real u_lin = m_settings.get_real("initial_conditions/u_lin");
    real v_lin = m_settings.get_real("initial_conditions/v_lin");
    real w_lin = m_settings.get_real("initial_conditions/w_lin");
    real x_shift = m_settings.get_real("initial_conditions/x_shift");
    real y_shift = m_settings.get_real("initial_conditions/y_shift");
    real z_shift = m_settings.get_real("initial_conditions/z_shift");
    real l = m_settings.get_real("initial_conditions/l");

    Functions::gauss_bubble(m_u_analytical_solution, t,
            u_lin, v_lin, w_lin, x_shift, y_shift, z_shift, l);
    Functions::gauss_bubble(m_v_analytical_solution, t,
            u_lin, v_lin, w_lin, x_shift, y_shift, z_shift, l);
    Functions::gauss_bubble(m_w_analytical_solution, t,
            u_lin, v_lin, w_lin, x_shift, y_shift, z_shift, l);
}

void Solution::exp_sinus_prod(const real t) {
    // Diffusion test case
    real nu = m_settings.get_real("physical_parameters/nu");
    real l = m_settings.get_real("initial_conditions/l");

    Functions::exp_sinus_prod(m_u_analytical_solution, t, nu, l);
    Functions::exp_sinus_prod(m_v_analytical_solution, t, nu, l);
    Functions::exp_sinus_prod(m_w_analytical_solution, t, nu, l);
}

void Solution::exp_sinus_sum(const real t) {
    // Diffusion test case
    real nu = m_settings.get_real("physical_parameters/nu");

    Functions::exp_sinus_sum(m_u_analytical_solution,
                           m_v_analytical_solution,
                           m_w_analytical_solution, t, nu);
}

void Solution::hat(const real) {
    // Diffusion test case
    real start_x = m_settings.get_real("initial_conditions/x1");
    real end_x = m_settings.get_real("initial_conditions/x2");
    real start_y = m_settings.get_real("initial_conditions/y1");
    real end_y = m_settings.get_real("initial_conditions/y2");
    real start_z = m_settings.get_real("initial_conditions/z1");
    real end_z = m_settings.get_real("initial_conditions/z2");
    real val_in = m_settings.get_real("initial_conditions/val_in");
    real val_out = m_settings.get_real("initial_conditions/val_out");

    Functions::hat(m_u_analytical_solution,
            start_x, end_x, start_y, end_y, start_z, end_z, val_in, val_out);  // TODO time dependency?
    Functions::hat(m_v_analytical_solution,
            start_x, end_x, start_y, end_y, start_z, end_z, val_in, val_out);
    Functions::hat(m_w_analytical_solution,
            start_x, end_x, start_y, end_y, start_z, end_z, val_in, val_out);
}

void Solution::sin_sin_sin(const real) {
    // Pressure test case
    real l = m_settings.get_real("initial_conditions/l");
    Functions::fac_sin_sin_sin(m_p_analytical_solution, l);  // TODO time dependency?
}

void Solution::mcDermott(const real t) {
    // NavierStokes test case
    real nu = m_settings.get_real("physical_parameters/nu");
    real A = m_settings.get_real("initial_conditions/A");

    Functions::mcdermott(m_u_analytical_solution,
                         m_v_analytical_solution,
                         m_w_analytical_solution,
                         m_p_analytical_solution,
                         t, nu, A);
}

void Solution::vortex(const real) {
    real u_lin = m_settings.get_real("initial_conditions/u_lin");
    real v_lin = m_settings.get_real("initial_conditions/v_lin");
    real pa = m_settings.get_real("initial_conditions/pa");
    real rhoa = m_settings.get_real("initial_conditions/rhoa");

    Functions::vortex(m_u_analytical_solution,
                      m_v_analytical_solution,
                      m_w_analytical_solution,
                      m_p_analytical_solution,
                      u_lin, v_lin, pa, rhoa);  // TODO time dependency
}

void Solution::vortex_y(const real) {
    real u_lin = m_settings.get_real("initial_conditions/u_lin");
    real v_lin = m_settings.get_real("initial_conditions/v_lin");
    real pa = m_settings.get_real("initial_conditions/pa");
    real rhoa = m_settings.get_real("initial_conditions/rhoa");

    Functions::vortex_y(m_u_analytical_solution,
                       m_v_analytical_solution,
                       m_w_analytical_solution,
                       m_p_analytical_solution,
                       u_lin, v_lin, pa, rhoa);  // TODO time dependency
}

void Solution::beltrami(const real t) {
    real a = m_settings.get_real("initial_conditions/a");  // 0.25 * M_PI;
    real d = m_settings.get_real("initial_conditions/d");  // 0.5 * M_PI;
    real nu = m_settings.get_real("physical_parameters/nu");  // 1;
    Functions::beltrami(m_u_analytical_solution,
                        m_v_analytical_solution,
                        m_w_analytical_solution,
                        m_p_analytical_solution,
                        t, a, d, nu);
}

void Solution::zero(const real) {
    // do nothing
}

void Solution::buoyancy_mms(const real t) {
    real g = m_settings.get_real("physical_parameters/g");  // 1;
    real nu = m_settings.get_real("physical_parameters/nu");  // 1;
    real beta = m_settings.get_real("physical_parameters/beta");
    real rhoa = m_settings.get_real("initial_conditions/rhoa");

    Functions::buoyancy_mms(m_u_analytical_solution,
                           m_v_analytical_solution,
                           m_w_analytical_solution,
                           m_p_analytical_solution,
                           m_T_analytical_solution,
                           t, nu, beta, g, rhoa);
}

// =================== Calculate analytical solution based on test case ==================
// ***************************************************************************************
/// \brief  calculates analytical solution
/// \param  t   time
// ***************************************************************************************
void Solution::calc_analytical_solution(const real t) {
    if (m_current_time_step == t) {
        return;
    }

    if (!m_has_analytical_solution) {
        return;
    }

    m_current_time_step = t;
    (*this.*m_init_function)(t);
}
