/// \file 		Solution.cpp
/// \brief 		Calculates analytical solution
/// \details	This class calculates the analytical solution of different test cases
/// \date 		Jul 11, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "Solution.h"
#include "../utility/Parameters.h"
#include "../Functions.h"

Solution::Solution() {
    ua = new Field(FieldType::U, 0.0);
    va = new Field(FieldType::V, 0.0);
    wa = new Field(FieldType::W, 0.0);
    pa = new Field(FieldType::P, 0.0);
    Ta = new Field(FieldType::T, 0.0);

    init();
}

Solution::~Solution() {
    delete ua;
    delete va;
    delete wa;
    delete pa;
    delete Ta;
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
        m_init_function = &Solution::fac_sin_sin_sin;
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
        std::cout << "Analytical solution set to zero!" << std::endl;
        m_init_function = &Solution::zero;
    }
}

void Solution::gauss_bubble(const real t) {
    // Advection test case
    Functions::GaussBubble(ua, t);
    Functions::GaussBubble(va, t);
    Functions::GaussBubble(wa, t);
}

void Solution::exp_sinus_prod(const real t) {
    // Diffusion test case
    Functions::ExpSinusProd(ua, t);
    Functions::ExpSinusProd(va, t);
    Functions::ExpSinusProd(wa, t);
}

void Solution::exp_sinus_sum(const real t) {
    // Diffusion test case
    Functions::ExpSinusSum(ua, va, wa, t);
}

void Solution::hat(const real t) {
    // Diffusion test case
    Functions::Hat(ua); // TODO time dependency?
    Functions::Hat(va);
    Functions::Hat(wa);
}

void Solution::fac_sin_sin_sin(const real t) {
// Pressure test case
    Functions::FacSinSinSin(pa); // TODO time dependency?
}

void Solution::mcDermott(const real t) {
// NavierStokes test case
    Functions::McDermott(ua, va, wa, pa, t);
}

void Solution::vortex(const real t) {
    Functions::Vortex(ua, va, wa, pa); // TODO time dependency
}

void Solution::vortex_y(const real t) {
    Functions::VortexY(ua, va, wa, pa); // TODO time dependency
}

void Solution::beltrami(const real t) {
    Functions::Beltrami(ua, va, wa, pa, t);
}

void Solution::zero(const real t) {
    // do nothing
}

void Solution::buoyancy_mms(const real t) {
    Functions::BuoyancyMMS(ua, va, wa, pa, Ta, t);
}

// =================== Calculate analytical solution based on test case ==================
// ***************************************************************************************
/// \brief  calculates analytical solution
/// \param	t		time
// ***************************************************************************************
void Solution::CalcAnalyticalSolution(const real t) {
    if (m_current_timestep == t) {
        return;
    }
    m_current_timestep = t;
    (*this.*m_init_function)(t);
}
