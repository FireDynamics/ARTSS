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

    CalcAnalyticalSolution(0.);
}

Solution::~Solution() {
    delete ua;
    delete va;
    delete wa;
    delete pa;
    delete Ta;
}

// =================== Calculate analytical solution based on test case ==================
// ***************************************************************************************
/// \brief  calculates analytical solution
/// \param	t		time
// ***************************************************************************************
void Solution::CalcAnalyticalSolution(const real t) {
//TODO not every function has a full analytic solution with time parameter
    auto params = Parameters::getInstance();

    // if analytical solution available
    if (params->get("solver/solution/available") == "Yes") {
        std::string initialCondition = params->get("initial_conditions/usr_fct");
        // Advection test case
        if (initialCondition == FunctionNames::GaussBubble) {
            Functions::GaussBubble(ua, t);
            Functions::GaussBubble(va, t);
            Functions::GaussBubble(wa, t);
        }

            // Diffusion test case
        else if (initialCondition == FunctionNames::ExpSinusProd) {

            Functions::ExpSinusProd(ua, t);
            Functions::ExpSinusProd(va, t);
            Functions::ExpSinusProd(wa, t);
        }
            // Diffusion test case
        else if (initialCondition == FunctionNames::ExpSinusSum) {

            Functions::ExpSinusSum(ua, va, wa, t);
        }
            // Diffusion test case
        else if (initialCondition == FunctionNames::Hat) {
            Functions::Hat(ua);
            Functions::Hat(va);
            Functions::Hat(wa);
        }
            // Pressure test case
        else if (initialCondition == FunctionNames::SinSinSin) {
            Functions::FacSinSinSin(pa);
        }

            // NavierStokes test cases
        else if (initialCondition == FunctionNames::McDermott) {
            Functions::McDermott(ua, va, wa, pa, t);
        } else if (initialCondition == FunctionNames::Vortex) {
            Functions::Vortex(ua, va, wa, pa);
        } else if (initialCondition == FunctionNames::VortexY) {
            Functions::VortexY(ua, va, wa, pa);
        } else if (initialCondition == FunctionNames::Beltrami) {
            Functions::Beltrami(ua, va, wa, pa, t);
        } else if (initialCondition == FunctionNames::BuoyancyMMS) {
            Functions::BuoyancyMMS(ua, va, wa, pa, Ta, t);
        } else {
            std::cout << "Analytical solution set to zero!" << std::endl;
        }
    }
}
