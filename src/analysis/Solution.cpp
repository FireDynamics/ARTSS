/// \file       Solution.cpp
/// \brief      Calculates analytical solution
/// \details    This class calculates the analytical solution of different test cases
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Solution.h"


Solution::Solution(std::string &initial_condition) :
    u_a(Field(FieldType::U)),
    v_a(Field(FieldType::V)),
    w_a(Field(FieldType::W)),
    p_a(Field(FieldType::P)),
    T_a(Field(FieldType::T))  {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    init(initial_condition);

    auto params = Parameters::getInstance();
    m_has_analytical_solution = (params->get("solver/solution/available") == "Yes");
}

void Solution::init(std::string &initial_condition) {
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
    }
}

void Solution::gauss_bubble(const real t) {
    // Advection test case
    Functions::GaussBubble(u_a, t);
    Functions::GaussBubble(v_a, t);
    Functions::GaussBubble(w_a, t);
}

void Solution::exp_sinus_prod(const real t) {
    // Diffusion test case
    Functions::ExpSinusProd(u_a, t);
    Functions::ExpSinusProd(v_a, t);
    Functions::ExpSinusProd(w_a, t);
}

void Solution::exp_sinus_sum(const real t) {
    // Diffusion test case
    Functions::ExpSinusSum(u_a, v_a, w_a, t);
}

void Solution::hat(const real) {
    // Diffusion test case
    Functions::Hat(u_a);  // TODO time dependency?
    Functions::Hat(v_a);
    Functions::Hat(w_a);
}

void Solution::sin_sin_sin(const real) {
// Pressure test case
    Functions::FacSinSinSin(p_a);  // TODO time dependency?
}

void Solution::mcDermott(const real t) {
// NavierStokes test case
    Functions::McDermott(u_a, v_a, w_a, p_a, t);
}

void Solution::vortex(const real) {
    Functions::Vortex(u_a, v_a, w_a, p_a);  // TODO time dependency
}

void Solution::vortex_y(const real) {
    Functions::VortexY(u_a, v_a, w_a, p_a);  // TODO time dependency
}

void Solution::beltrami(const real t) {
    Functions::Beltrami(u_a, v_a, w_a, p_a, t);
}

void Solution::zero(const real) {
    // do nothing
}

void Solution::buoyancy_mms(const real t) {
    Functions::BuoyancyMMS(u_a, v_a, w_a, p_a, T_a, t);
}

// =================== Calculate analytical solution based on test case ==================
// ***************************************************************************************
/// \brief  calculates analytical solution
/// \param  t   time
// ***************************************************************************************
void Solution::calc_analytical_solution(real t) {
    if (m_current_time_step == t) {
        return;
    }

    if (!m_has_analytical_solution) {
        return;
    }

    m_current_time_step = t;
    (*this.*m_init_function)(t);
}
