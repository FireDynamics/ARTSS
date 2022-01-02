/// \file       Solution.cpp
/// \brief      Calculates analytical solution
/// \details    This class calculates the analytical solution of different test cases
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Solution.h"
#include "../Functions.h"


Solution::Solution(const Settings::initial_conditions_parameters &ic_parameters,
                   const Settings::solver::solution &solution_parameters) :
        m_ic_settings(ic_parameters),
        m_solution_settings(solution_parameters),
        m_u_analytical_solution(Field(FieldType::U)),
        m_v_analytical_solution(Field(FieldType::V)),
        m_w_analytical_solution(Field(FieldType::W)),
        m_p_analytical_solution(Field(FieldType::P)),
        m_T_analytical_solution(Field(FieldType::T)) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    // set function pointer to chosen initial condition
    std::string initial_condition = m_ic_settings.usr_fct;
    if (initial_condition == FunctionNames::gauss_bubble) {
        m_init_function2 = [this](real t){this->gauss_bubble(t);};
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
    auto gauss = std::get<Settings::initial_conditions::gauss_bubble>(m_ic_settings.ic.value());
    Functions::gauss_bubble(m_u_analytical_solution, t,gauss);
    Functions::gauss_bubble(m_v_analytical_solution, t,gauss);
    Functions::gauss_bubble(m_w_analytical_solution, t,gauss);
}

void Solution::exp_sinus_prod(const real t) {
    // Diffusion test case
    auto exp_sinus_prod = std::get<Settings::initial_conditions::exp_sinus_prod>(m_ic_settings.ic.value());
    Functions::exp_sinus_prod(m_u_analytical_solution, t, exp_sinus_prod);
    Functions::exp_sinus_prod(m_v_analytical_solution, t, exp_sinus_prod);
    Functions::exp_sinus_prod(m_w_analytical_solution, t, exp_sinus_prod);
}

void Solution::exp_sinus_sum(const real t) {
    // Diffusion test case
    Functions::exp_sinus_sum(m_u_analytical_solution,
                           m_v_analytical_solution,
                           m_w_analytical_solution, t);
}

void Solution::hat(const real) {
    // Diffusion test case
    auto hat = std::get<Settings::initial_conditions::hat>(m_ic_settings.ic.value());
    Functions::hat(m_u_analytical_solution, hat);
    Functions::hat(m_v_analytical_solution, hat);
    Functions::hat(m_w_analytical_solution, hat);
}

void Solution::sin_sin_sin(const real) {
    // Pressure test case
    auto sin_sin_sin = std::get<Settings::initial_conditions::sin_sin_sin>(m_ic_settings.ic.value());
    Functions::fac_sin_sin_sin(m_p_analytical_solution, sin_sin_sin);  // TODO time dependency?
}

void Solution::mcDermott(const real t) {
    // NavierStokes test case
    auto mc_dermott = std::get<Settings::initial_conditions::mc_dermott>(m_ic_settings.ic.value());
    Functions::mcdermott(m_u_analytical_solution,
                         m_v_analytical_solution,
                         m_w_analytical_solution,
                         m_p_analytical_solution,
                         t, mc_dermott);
}

void Solution::vortex(const real) {
    auto vortex = std::get<Settings::initial_conditions::vortex>(m_ic_settings.ic.value());
    Functions::vortex(m_u_analytical_solution,
                      m_v_analytical_solution,
                      m_w_analytical_solution,
                      m_p_analytical_solution,
                      vortex);  // TODO time dependency
}

void Solution::vortex_y(const real) {
    auto vortex = std::get<Settings::initial_conditions::vortex>(m_ic_settings.ic.value());
    Functions::vortex_y(m_u_analytical_solution,
                       m_v_analytical_solution,
                       m_w_analytical_solution,
                       m_p_analytical_solution,
                       vortex);  // TODO time dependency
}

void Solution::beltrami(const real t) {
    auto beltrami = std::get<Settings::initial_conditions::beltrami>(m_ic_settings.ic.value());
    Functions::beltrami(m_u_analytical_solution,
                        m_v_analytical_solution,
                        m_w_analytical_solution,
                        m_p_analytical_solution,
                        t, beltrami);
}

void Solution::zero(const real) {
    // do nothing
}

void Solution::buoyancy_mms(const real t) {
    Functions::buoyancy_mms(m_u_analytical_solution,
                           m_v_analytical_solution,
                           m_w_analytical_solution,
                           m_p_analytical_solution,
                           m_T_analytical_solution,
                           t);
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

    if (!m_solution_settings.analytical_solution) {
        return;
    }

    m_current_time_step = t;
    (*this.*m_init_function)(t);
}
