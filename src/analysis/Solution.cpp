/// \file       Solution.cpp
/// \brief      Calculates analytical solution
/// \details    This class calculates the analytical solution of different test cases
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Solution.h"
#include "../utility/Utility.h"
#include "../utility/Parameters.h"
#include "../Functions.h"

Solution::Solution() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    u_a = new Field(FieldType::U, 0.0);
    v_a = new Field(FieldType::V, 0.0);
    w_a = new Field(FieldType::W, 0.0);
    p_a = new Field(FieldType::P, 0.0);
    T_a = new Field(FieldType::T, 0.0);

    init();
}

Solution::~Solution() {
    delete u_a;
    delete v_a;
    delete w_a;
    delete p_a;
    delete T_a;
}

void Solution::init() {
    auto params = Parameters::getInstance();
    std::string initialCondition = params->get("initial_conditions/usr_fct");

    // set function pointer to chosen initial condition
    if (initialCondition == FunctionNames::GaussBubble) {
        m_init_function = &Solution::gauss_bubble;
    } else if (initialCondition == FunctionNames::ExpSinusProd) {
        m_init_function = &Solution::exp_sinus_prod;
    } else if (initialCondition == FunctionNames::ExpSinusSum) {
        m_init_function = &Solution::exp_sinus_sum;
    } else if (initialCondition == FunctionNames::Hat) {
        m_init_function = &Solution::hat;
    } else if (initialCondition == FunctionNames::SinSinSin) {
        m_init_function = &Solution::sin_sin_sin;
    } else if (initialCondition == FunctionNames::McDermott) {
        m_init_function = &Solution::mcDermott;
    } else if (initialCondition == FunctionNames::Vortex) {
        m_init_function = &Solution::vortex;
    } else if (initialCondition == FunctionNames::VortexY) {
        m_init_function = &Solution::vortex_y;
    } else if (initialCondition == FunctionNames::Beltrami) {
        m_init_function = &Solution::beltrami;
    } else if (initialCondition == FunctionNames::BuoyancyMMS) {
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

void Solution::hat(const real t) {
    // Diffusion test case
    Functions::Hat(u_a); // TODO time dependency?
    Functions::Hat(v_a);
    Functions::Hat(w_a);
}

void Solution::sin_sin_sin(real t) {
// Pressure test case
    Functions::FacSinSinSin(p_a); // TODO time dependency?
}

void Solution::mcDermott(const real t) {
// NavierStokes test case
    Functions::McDermott(u_a, v_a, w_a, p_a, t);
}

void Solution::vortex(const real t) {
    Functions::Vortex(u_a, v_a, w_a, p_a); // TODO time dependency
}

void Solution::vortex_y(const real t) {
    Functions::VortexY(u_a, v_a, w_a, p_a); // TODO time dependency
}

void Solution::beltrami(const real t) {
    Functions::Beltrami(u_a, v_a, w_a, p_a, t);
}

void Solution::zero(const real t) {
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
    m_current_time_step = t;
    (*this.*m_init_function)(t);
}
