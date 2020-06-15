/// \file 		Solution.cpp
/// \brief 		Calculates analytical solution
/// \details	This class calculates the analytical solution of different test cases
/// \date 		Jul 11, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#include "Solution.h"
#include "../utility/Parameters.h"
#include "../Functions.h"
#include "../boundary/BoundaryController.h"
#include "../utility/Utility.h"

Solution::Solution() {
    m_logger = Utility::createLogger(typeid(this).name());
    ua = new Field(FieldType::U, 0.0);
    va = new Field(FieldType::V, 0.0);
    wa = new Field(FieldType::W, 0.0);
    u0a = new Field(FieldType::U, 0.0);
    v0a = new Field(FieldType::V, 0.0);
    w0a = new Field(FieldType::W, 0.0);
    pa = new Field(FieldType::P, 0.0);
    p0a = new Field(FieldType::P, 0.0);
    Ta = new Field(FieldType::T, 0.0);
    T0a = new Field(FieldType::T, 0.0);

    SetUp();
}

Solution::~Solution() {
    delete ua;
    delete va;
    delete wa;
    delete u0a;
    delete v0a;
    delete w0a;
    delete pa;
    delete p0a;
    delete Ta;
    delete T0a;
}

// =================== Calculate analytical solution based on test case ==================
// ***************************************************************************************
/// \brief  calculates analytical solution
/// \param	t		time
// ***************************************************************************************
void Solution::CalcAnalyticalSolution(const real t) {
//TODO not every function has a full analytic solution with time parameter
//TODO put each if block into a function and use a function pointer to prevent if queries each time. set function pointer in setUp function. also prevents multiple outputs of analytical solution set to zero
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
            m_logger->info("Analytical solution set to zero!");
        }
    }

    if (t == 0.) {

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t size_bList = boundary->getSize_boundaryList();

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            u0a->data[idx] = ua->data[idx];
            v0a->data[idx] = va->data[idx];
            w0a->data[idx] = wa->data[idx];
            p0a->data[idx] = pa->data[idx];
            T0a->data[idx] = Ta->data[idx];
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            u0a->data[idx] = ua->data[idx];
            v0a->data[idx] = va->data[idx];
            w0a->data[idx] = wa->data[idx];
            p0a->data[idx] = pa->data[idx];
            T0a->data[idx] = Ta->data[idx];
        }
    }
}

// ============ Calculate analytical solution at t=0 based on test case ==================
// ***************************************************************************************
/// \brief  calculates analytical solution at t=0
// ***************************************************************************************
void Solution::SetUp() {
    CalcAnalyticalSolution(0.);
}
