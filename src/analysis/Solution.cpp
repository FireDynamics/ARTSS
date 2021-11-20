/// \file       Solution.cpp
/// \brief      Calculates analytical solution
/// \details    This class calculates the analytical solution of different test cases
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Solution.h"
#include "../Functions.h"


Solution::Solution(Settings const &settings, const std::string &initial_condition, bool has_analytical_solution) :
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
    if (initial_condition == FunctionNames::GaussBubble) {
        m_init_function = &Solution::gauss_bubble;
    } else if (initial_condition == FunctionNames::ExpSinusProd) {
        m_init_function = &Solution::exp_sinus_prod;
    } else if (initial_condition == FunctionNames::ExpSinusSum) {
        m_init_function = &Solution::exp_sinus_sum;
    } else if (initial_condition == FunctionNames::Hat) {
        m_init_function = &Solution::hat;
    } else if (initial_condition == FunctionNames::SinSinSin) {
        m_init_function = &Solution::sin_sin_sin;
    } else if (initial_condition == FunctionNames::McDermott) {
        m_init_function = &Solution::mcDermott;
    } else if (initial_condition == FunctionNames::Vortex) {
        m_init_function = &Solution::vortex;
    } else if (initial_condition == FunctionNames::VortexY) {
        m_init_function = &Solution::vortex_y;
    } else if (initial_condition == FunctionNames::Beltrami) {
        m_init_function = &Solution::beltrami;
    } else if (initial_condition == FunctionNames::BuoyancyMMS) {
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
    Functions::GaussBubble(m_settings, m_u_analytical_solution, t);
    Functions::GaussBubble(m_settings, m_v_analytical_solution, t);
    Functions::GaussBubble(m_settings, m_w_analytical_solution, t);
}

void Solution::exp_sinus_prod(const real t) {
    // Diffusion test case
    Functions::ExpSinusProd(m_settings, m_u_analytical_solution, t);
    Functions::ExpSinusProd(m_settings, m_v_analytical_solution, t);
    Functions::ExpSinusProd(m_settings, m_w_analytical_solution, t);
}

void Solution::exp_sinus_sum(const real t) {
    // Diffusion test case
    Functions::ExpSinusSum(m_settings,
                           m_u_analytical_solution,
                           m_v_analytical_solution,
                           m_w_analytical_solution, t);
}

void Solution::hat(const real) {
    // Diffusion test case
    Functions::Hat(m_settings, m_u_analytical_solution);  // TODO time dependency?
    Functions::Hat(m_settings, m_v_analytical_solution);
    Functions::Hat(m_settings, m_w_analytical_solution);
}

void Solution::sin_sin_sin(const real) {
// Pressure test case
    Functions::FacSinSinSin(m_settings, m_p_analytical_solution);  // TODO time dependency?
}

void Solution::mcDermott(const real t) {
// NavierStokes test case
    Functions::McDermott(m_settings,
                         m_u_analytical_solution,
                         m_v_analytical_solution,
                         m_w_analytical_solution,
                         m_p_analytical_solution, t);
}

void Solution::vortex(const real) {
    Functions::Vortex(m_settings,
                      m_u_analytical_solution,
                      m_v_analytical_solution,
                      m_w_analytical_solution,
                      m_p_analytical_solution);  // TODO time dependency
}

void Solution::vortex_y(const real) {
    Functions::VortexY(m_settings,
                       m_u_analytical_solution,
                       m_v_analytical_solution,
                       m_w_analytical_solution,
                       m_p_analytical_solution);  // TODO time dependency
}

void Solution::beltrami(const real t) {
    Functions::Beltrami(m_settings,
                        m_u_analytical_solution,
                        m_v_analytical_solution,
                        m_w_analytical_solution,
                        m_p_analytical_solution,
                        t);
}

void Solution::zero(const real) {
    // do nothing
}

void Solution::buoyancy_mms(const real t) {
    Functions::BuoyancyMMS(m_settings,
                           m_u_analytical_solution,
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

    if (!m_has_analytical_solution) {
        return;
    }

    m_current_time_step = t;
    (*this.*m_init_function)(t);
}
