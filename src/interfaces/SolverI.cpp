/// \file 		SolverI.h
/// \brief 		Interface holds solvers for solving governing equations
/// \date 		May 20, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>
#include "SolverI.h"
#include "../utility/Parameters.h"
#include "../Functions.h"
#include "../source/ExplicitEulerSource.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"
#include "../solver/SolverSelection.h"
#include "../utility/Utility.h"

SolverI::SolverI() {
#ifndef PROFILING
    m_logger = Utility::createLogger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();

    // Variables
    // Velocities
    u = new Field(FieldType::U, 0.0);
    u0 = new Field(FieldType::U, 0.0);
    u_tmp = new Field(FieldType::U, 0.0);

    v = new Field(FieldType::V, 0.0);
    v0 = new Field(FieldType::V, 0.0);
    v_tmp = new Field(FieldType::V, 0.0);

    w = new Field(FieldType::W, 0.0);
    w0 = new Field(FieldType::W, 0.0);
    w_tmp = new Field(FieldType::W, 0.0);

    // Turbulent diffusivities
    nu_t = new Field(FieldType::U, 0.0);
    kappa_t = new Field(FieldType::T, 0.0);
    gamma_t = new Field(FieldType::RHO, 0.0);

    // Pressure
    p = new Field(FieldType::P, 0.0);
    p0 = new Field(FieldType::P, 0.0);
    rhs = new Field(FieldType::P, 0.0);

    // Temperature
    T = new Field(FieldType::T, 0.0);
    T0 = new Field(FieldType::T, 0.0);
    T_tmp = new Field(FieldType::T, 0.0);
    T_a = new Field(FieldType::T, 0.0);

    // Concentration
    C = new Field(FieldType::RHO, 0.0);
    C0 = new Field(FieldType::RHO, 0.0);
    C_tmp = new Field(FieldType::RHO, 0.0);

    // Forces
    f_x = new Field(FieldType::U, 0.0);
    f_y = new Field(FieldType::V, 0.0);
    f_z = new Field(FieldType::W, 0.0);

    // Sources
    S_T = new Field(FieldType::T, 0.0);
    S_C = new Field(FieldType::RHO, 0.0);

    // Fields for sight of boundaries
    sight = new Field(FieldType::RHO, 1.0);

    // Source term for momentum/temperature/concentration
    m_string_solver = params->get("solver/description");
    if (m_string_solver == SolverTypes::NSSolver or \
        m_string_solver == SolverTypes::NSTempSolver or \
        m_string_solver == SolverTypes::NSTempTurbSolver or \
        m_string_solver == SolverTypes::NSTempConSolver or \
        m_string_solver == SolverTypes::NSTempTurbConSolver or \
        m_string_solver == SolverTypes::NSTurbSolver) {
        if (params->get("solver/source/type") == SourceMethods::ExplicitEuler) {
            this->sou_vel = new ExplicitEulerSource();
        } else {
            m_logger->critical("Source method not yet implemented! Simulation stopped!");
            std::exit(1);
            //TODO Error handling
        }
    }
    if (m_string_solver == SolverTypes::NSTempSolver or \
        m_string_solver == SolverTypes::NSTempConSolver or \
        m_string_solver == SolverTypes::NSTempTurbConSolver or \
        m_string_solver == SolverTypes::NSTempTurbSolver) {

        // Source of temperature
        if (params->get("solver/temperature/source/type") == SourceMethods::ExplicitEuler) {
            this->sou_temp = new ExplicitEulerSource();
        } else {
            m_logger->critical("Source method not yet implemented! Simulation stopped!");
            std::exit(1);
            //TODO Error handling
        }
    }
    if (m_string_solver == SolverTypes::NSTempConSolver or \
        m_string_solver == SolverTypes::NSTempTurbConSolver) {

        // Source of concentration
        if (params->get("solver/concentration/source/type") == SourceMethods::ExplicitEuler) {
            this->sou_con = new ExplicitEulerSource();
        } else {
            m_logger->critical("Source method not yet implemented! Simulation stopped!");
            std::exit(1);
            //TODO Error handling
        }
    }

    SetUp();
}

// ==================================== Destructor ====================================
// ***************************************************************************************
SolverI::~SolverI() {
    delete u;
    delete v;
    delete w;
    delete u0;
    delete v0;
    delete w0;
    delete u_tmp;
    delete v_tmp;
    delete w_tmp;

    delete nu_t;
    delete kappa_t;
    delete gamma_t;

    delete p;
    delete p0;
    delete rhs;

    delete T;
    delete T0;
    delete T_tmp;
    delete T_a;

    delete C;
    delete C0;
    delete C_tmp;

    delete f_x;
    delete f_y;
    delete f_z;
    delete S_T;
    delete S_C;

    delete sight;

    if (m_string_solver == SolverTypes::NSSolver or \
        m_string_solver == SolverTypes::NSTempSolver or \
        m_string_solver == SolverTypes::NSTempConSolver or \
        m_string_solver == SolverTypes::NSTempTurbConSolver or \
        m_string_solver == SolverTypes::NSTempTurbSolver or \
        m_string_solver == SolverTypes::NSTurbSolver) {
        delete sou_vel;
    }
    if (m_string_solver == SolverTypes::NSTempSolver or \
        m_string_solver == SolverTypes::NSTempTurbSolver) {
        delete sou_temp;
    }
    if (m_string_solver == SolverTypes::NSTempConSolver or \
        m_string_solver == SolverTypes::NSTempTurbConSolver) {
        delete sou_con;
    }
}

// ==================================== Set Up ========================================
// ***************************************************************************************
/// \brief  initializes numerical and temporary solution
// ***************************************************************************************
void SolverI::SetUp() {
    m_logger->info("Start initializing....");

    // Initialization of variables at time t=0
    Init();

    // Initialization of temp and new variables
    Init_f(u_tmp, u0->data);
    Init_f(v_tmp, v0->data);
    Init_f(w_tmp, w0->data);
    Init_f(u, u0->data);
    Init_f(v, v0->data);
    Init_f(w, w0->data);
    Init_f(p, p0->data);
    Init_f(T, T0->data);
    Init_f(T_tmp, T0->data);
    Init_f(T_a, T0->data);
    Init_f(C, C0->data);
    Init_f(C_tmp, C0->data);
}

//================================== 0) Initialization ===================================
// ***************************************************************************************
/// \brief  initializes variables \a u, \a v, \a p, \a d, \a T at \f$ t=0 \f$
// ***************************************************************************************
void SolverI::Init() {

    auto params = Parameters::getInstance();
    std::string string_init_usr_fct = params->get("initial_conditions/usr_fct");

    if (string_init_usr_fct == FunctionNames::GaussBubble) {
        if (m_string_solver == SolverTypes::AdvectionSolver) {
            Functions::GaussBubble(u0, 0.);
            Functions::GaussBubble(v0, 0.);
            Functions::GaussBubble(w0, 0.);
        }
    } else if (string_init_usr_fct == FunctionNames::Drift) {
        if (m_string_solver == SolverTypes::AdvectionSolver or \
            m_string_solver == SolverTypes::NSSolver or \
            m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver or \
            m_string_solver == SolverTypes::NSTurbSolver) {
            Functions::Drift(u0, v0, w0, p0);
        }
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == FunctionNames::ExpSinusProd) {
        // Diffusion test case
        if (m_string_solver == SolverTypes::DiffusionSolver or \
            m_string_solver == SolverTypes::DiffusionTurbSolver) {
            Functions::ExpSinusProd(u0, 0.);
            Functions::ExpSinusProd(v0, 0.);
            Functions::ExpSinusProd(w0, 0.);
        }
    } else if (string_init_usr_fct == FunctionNames::Hat) {
        if (m_string_solver == SolverTypes::DiffusionSolver or \
            m_string_solver == SolverTypes::DiffusionTurbSolver) {
            Functions::Hat(u0);
            Functions::Hat(v0);
            Functions::Hat(w0);
        }
    } else if (string_init_usr_fct == FunctionNames::ExpSinusSum) {
        // Burgers (=nonlinear Advection + Diffusion) test case
        if (m_string_solver == SolverTypes::AdvectionDiffusionSolver) {
            Functions::ExpSinusSum(u0, v0, w0, 0.);
        }
    } else if (string_init_usr_fct == FunctionNames::SinSinSin) {
        if (m_string_solver == SolverTypes::PressureSolver) {
            // Pressure test case
            Functions::SinSinSin(rhs);
        }
    } else if (string_init_usr_fct == FunctionNames::McDermott) {
        // NavierStokes test case: McDermott (no force, no temperature) 2D
        if (m_string_solver == SolverTypes::NSSolver or \
            m_string_solver == SolverTypes::NSTurbSolver or \
            m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            Functions::McDermott(u0, v0, w0, p0, 0.);
            Init_c(p0, 0.);
        }
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == FunctionNames::Vortex) {
        // NavierStokes test case: Vortex (no force, no temperature) 2D
        if (m_string_solver == SolverTypes::NSSolver or \
            m_string_solver == SolverTypes::NSTurbSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver or \
            m_string_solver == SolverTypes::NSTempSolver) {
            Functions::Vortex(u0, v0, w0, p0);
            Init_c(p0, 0.);
        }
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == FunctionNames::VortexY) {
        // NavierStokes test case: Vortex (no force, no temperature) 2D
        if (m_string_solver == SolverTypes::NSSolver or \
            m_string_solver == SolverTypes::NSTurbSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver or \
            m_string_solver == SolverTypes::NSTempSolver) {
            Functions::VortexY(u0, v0, w0, p0);
            Init_c(p0, 0.);
        }
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == FunctionNames::Beltrami) {
        // NavierStokes test case: Beltrami  (no force, no temperature) 3D
        if (m_string_solver == SolverTypes::NSSolver or \
            m_string_solver == SolverTypes::NSTurbSolver) {
            Functions::Beltrami(u0, v0, w0, p0, 0.);
            //Init_c(p0, 0.);
        }
    } else if (string_init_usr_fct == FunctionNames::BuoyancyMMS) {
        // NavierStokesTemp test case
        // MMS
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            Functions::BuoyancyMMS(u0, v0, w0, p0, T0, 0.0);
            Init_c(p0, 0.);
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == FunctionNames::Uniform) {
// Uniform Temperature unequal to zero
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            real val = params->getReal("initial_conditions/val");
            Functions::Uniform(T0, val);
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == FunctionNames::RandomT) {
        //Random
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            // Random temperature
            real Ta = params->getReal("initial_conditions/Ta");        // ambient temperature in KELVIN!
            real A = params->getReal("initial_conditions/A");            // amplitude
            size_t range = static_cast<size_t>(params->getInt("initial_conditions/range")); // range of random numbers (0-range)
            Functions::Random(T0, Ta, A, range);
        }
        //Random concentration
        if ((m_string_solver == SolverTypes::NSTempConSolver or \
             m_string_solver == SolverTypes::NSTempTurbConSolver)
            and params->get("initial_conditions/con_fct") == "RandomC") {

            real Ca = params->getReal("initial_conditions/Ca");        //ambient concentration
            real A_C = params->getReal("initial_conditions/A_C");        //amplitude
            size_t range_C = static_cast<size_t>(params->getInt("initial_conditions/range_C")); //range of random numbers (0-range)

            Functions::Random(C0, Ca, A_C, range_C);
        }
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == "LayersT") {
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            Functions::Layers(T0);
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == "LayersT,RandomT") {
        if (m_string_solver == SolverTypes::NSTempSolver or \
            m_string_solver == SolverTypes::NSTempConSolver or \
            m_string_solver == SolverTypes::NSTempTurbConSolver or \
            m_string_solver == SolverTypes::NSTempTurbSolver) {
            Functions::Layers(T0);
            real A = params->getReal("initial_conditions/A");        //amplitude
            size_t range = static_cast<size_t>(params->getReal("initial_conditions/range"));    //range of random numbers (0-range)
            Functions::Random(T0, T0, A, range);
            ForceSource();
            TemperatureSource();
        }
    } else if (string_init_usr_fct == FunctionNames::Zero) {
        // NavierStokes test case: Channel Flow (with uniform force in x-direction)
        if ((m_string_solver == SolverTypes::NSSolver or m_string_solver == SolverTypes::NSTurbSolver)
            and params->get("solver/source/force_fct") == FunctionNames::Uniform) {
            real val_x = 0;
            real val_y = 0;
            real val_z = 0;
            val_x = params->getReal("solver/source/val_x");
            val_y = params->getReal("solver/source/val_y");
            val_z = params->getReal("solver/source/val_z");
            std::string dir = params->get("solver/source/dir");

            if (dir.find('x') != std::string::npos) {
                Functions::Uniform(f_x, val_x);
            }
            if (dir.find('y') != std::string::npos) {
                Functions::Uniform(f_y, val_y);
            }
            if (dir.find('z') != std::string::npos) {
                Functions::Uniform(f_z, val_z);
            }
        } else {
#ifndef PROFILING
            m_logger->info("Initial values all set to zero!");
#endif
        }
    }

// Sight of boundaries
    auto boundary = BoundaryController::getInstance();
    size_t *iList = boundary->get_innerList_level_joined();
    size_t size_iList = boundary->getSize_innerList();

    for (size_t i = 0; i < size_iList; i++) {
        size_t idx = iList[i];
        sight->data[idx] = 0.0;
    }
}

// ========================================== Init =======================================
// ***************************************************************************************
/// \brief  initializes variable with another field
/// \param	out		output pointer
/// \param	in		input pointer
// ***************************************************************************************
void SolverI::Init_f(Field *out, const real *in) {

    auto boundary = BoundaryController::getInstance();

    size_t *iList = boundary->get_innerList_level_joined();
    size_t size_iList = boundary->getSize_innerList();
    size_t *bList = boundary->get_boundaryList_level_joined();
    size_t size_bList = boundary->getSize_boundaryList();
    size_t *oList = boundary->get_obstacleList();
    size_t size_oList = boundary->getSize_obstacleList();

    // inner cells
    for (size_t i = 0; i < size_iList; i++) {
        size_t idx = iList[i];
        out->data[idx] = in[idx];
    }
    // boundary
    for (size_t i = 0; i < size_bList; i++) {
        size_t idx = bList[i];
        out->data[idx] = in[idx];
    }
    // obstacle
    for (size_t i = 0; i < size_oList; i++) {
        size_t idx = oList[i];
        out->data[idx] = in[idx];
    }
}

// ========================================== Init =======================================
// ***************************************************************************************
/// \brief  initializes variable with constant
/// \param	out		output pointer
/// \param	in		constant real value
// ***************************************************************************************
void SolverI::Init_c(Field *out, const real in) {

    auto boundary = BoundaryController::getInstance();

    size_t *iList = boundary->get_innerList_level_joined();
    size_t size_iList = boundary->getSize_innerList();
    size_t *bList = boundary->get_boundaryList_level_joined();
    size_t size_bList = boundary->getSize_boundaryList();
    size_t *oList = boundary->get_obstacleList();
    size_t size_oList = boundary->getSize_obstacleList();

    //inner cells
    for (size_t i = 0; i < size_iList; i++) {
        size_t idx = iList[i];
        out->data[idx] = in;
    }
    //boundary
    for (size_t i = 0; i < size_bList; i++) {
        size_t idx = bList[i];
        out->data[idx] = in;
    }
    //obstacle
    for (size_t i = 0; i < size_oList; i++) {
        size_t idx = oList[i];
        out->data[idx] = in;
    }
}

// ========================================== Set up boundary =======================================
// ***************************************************************************************
/// \brief  initializes boundary cells
/// \param  sync	synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SolverI::SetUpBoundary(bool sync) {

    auto boundary = BoundaryController::getInstance();

    boundary->applyBoundary(u0->data, u0->GetType(), sync);
    boundary->applyBoundary(v0->data, v0->GetType(), sync);
    boundary->applyBoundary(w0->data, w0->GetType(), sync);
    boundary->applyBoundary(T0->data, T0->GetType(), sync);
    boundary->applyBoundary(p0->data, p0->GetType(), sync);
    boundary->applyBoundary(C0->data, C0->GetType(), sync);

    boundary->applyBoundary(u->data, u->GetType(), sync);
    boundary->applyBoundary(v->data, v->GetType(), sync);
    boundary->applyBoundary(w->data, w->GetType(), sync);
    boundary->applyBoundary(p->data, p->GetType(), sync);
    boundary->applyBoundary(T->data, T->GetType(), sync);
    boundary->applyBoundary(C->data, C->GetType(), sync);

    boundary->applyBoundary(u_tmp->data, u_tmp->GetType(), sync);
    boundary->applyBoundary(v_tmp->data, v_tmp->GetType(), sync);
    boundary->applyBoundary(w_tmp->data, w_tmp->GetType(), sync);
    boundary->applyBoundary(T_tmp->data, T_tmp->GetType(), sync);
    boundary->applyBoundary(T_a->data, T0->GetType(), sync);
    boundary->applyBoundary(C_tmp->data, C_tmp->GetType(), sync);
}

//======================================== Couple velocity ====================================
// ***************************************************************************************
/// \brief  couples vector (sets tmp and zero-th variables to current numerical solution)
/// \param  a		current field in x- direction (const)
/// \param  a0		zero-th field in x- direction
/// \param  a_tmp	temporal field in x- direction
/// \param  b		current field in y- direction (const)
/// \param  b0		zero-th field in y- direction
/// \param  b_tmp	temporal field in y- direction
/// \param  c		current field in z- direction (const)
/// \param  c0		zero-th field in z- direction
/// \param  c_tmp	temporal field in z- direction
/// \param  sync	synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SolverI::CoupleVector(const Field *a, Field *a0, Field *a_tmp, const Field *b, Field *b0, Field *b_tmp, const Field *c, Field *c0, Field *c_tmp, bool sync) {

    // local variables and parameters for GPU
    auto d_a = a->data;
    auto d_a0 = a0->data;
    auto d_a_tmp = a_tmp->data;
    auto d_b = b->data;
    auto d_b0 = b0->data;
    auto d_b_tmp = b_tmp->data;
    auto d_c = c->data;
    auto d_c0 = c0->data;
    auto d_c_tmp = c_tmp->data;

    auto size = Domain::getInstance()->GetSize(a0->GetLevel());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(    d_a[:size], d_a0[:size], d_a_tmp[:size], d_b[:size], d_b0[:size], d_b_tmp[:size], \
                            d_c[:size], d_c0[:size], d_c_tmp[:size], d_iList[:bsize_i], d_bList[:bsize_b], d_oList[:bsize_o])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_a0[i] = d_a[i];
            d_b0[i] = d_b[i];
            d_c0[i] = d_c[i];
            d_a_tmp[i] = d_a[i];
            d_b_tmp[i] = d_b[i];
            d_c_tmp[i] = d_c[i];
        }
        // boundary
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_bList[j];
            d_a0[i] = d_a[i];
            d_b0[i] = d_b[i];
            d_c0[i] = d_c[i];
            d_a_tmp[i] = d_a[i];
            d_b_tmp[i] = d_b[i];
            d_c_tmp[i] = d_c[i];
        }
        // obstacles
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t i = d_oList[j];
            d_a0[i] = d_a[i];
            d_b0[i] = d_b[i];
            d_c0[i] = d_c[i];
            d_a_tmp[i] = d_a[i];
            d_b_tmp[i] = d_b[i];
            d_c_tmp[i] = d_c[i];
        }

        if (sync) {
#pragma acc wait
        }
    }//end data region
}

//======================================= Couple scalar ==================================
// ***************************************************************************************
/// \brief  couples vector (sets tmp and zero-th variables to current numerical solution)
/// \param  a		current field in x- direction (const)
/// \param  a0		zero-th field in x- direction
/// \param  a_tmp	temporal field in x- direction
/// \param  sync	synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SolverI::CoupleScalar(const Field *a, Field *a0, Field *a_tmp, bool sync) {

    // local variables and parameters for GPU
    auto d_a = a->data;
    auto d_a0 = a0->data;
    auto d_a_tmp = a_tmp->data;

    auto size = Domain::getInstance()->GetSize(a0->GetLevel());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(d_a[:size], d_a0[:size], d_a_tmp[:size], d_iList[:bsize_i], d_bList[:bsize_b], d_oList[:bsize_o])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_a0[i] = d_a[i];
            d_a_tmp[i] = d_a[i];
        }
        // boundary
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_bList[j];
            d_a0[i] = d_a[i];
            d_a_tmp[i] = d_a[i];
        }
        // obstacles
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t i = d_oList[j];
            d_a0[i] = d_a[i];
            d_a_tmp[i] = d_a[i];
        }
        if (sync) {
#pragma acc wait
        }
    }//end data region
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates variables for the next iteration step or time dependent parameters such as temperature source function
/// \param  sync	synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SolverI::UpdateData(bool sync) {

    // local variables and parameters for GPU
    auto bsize = Domain::getInstance()->GetSize();

    const auto d_u = u->data;                        //due to const correctness
    const auto d_v = v->data;
    const auto d_w = w->data;
    auto d_u0 = u0->data;
    auto d_v0 = v0->data;
    auto d_w0 = w0->data;
    auto d_u_tmp = u_tmp->data;
    auto d_v_tmp = v_tmp->data;
    auto d_w_tmp = w_tmp->data;
    const auto d_p = p->data;
    auto d_p0 = p0->data;
    const auto d_T = T->data;
    auto d_T0 = T0->data;
    auto d_T_tmp = T_tmp->data;
    const auto d_C = C->data;
    auto d_C0 = C0->data;
    auto d_C_tmp = C_tmp->data;

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(    d_iList[:bsize_i], d_bList[:bsize_b], d_oList[:bsize_o], d_u[:bsize], d_v[:bsize], d_w[:bsize], \
                            d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], \
                            d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize])
    {
        // inner
#pragma acc parallel loop independent present(    d_iList[:bsize_i], d_u[:bsize], d_v[:bsize], d_w[:bsize], \
                                                d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], \
                                                d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], \
                                                d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize]) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t idx = d_iList[j];
            d_u0[idx] = d_u[idx];
            d_v0[idx] = d_v[idx];
            d_w0[idx] = d_w[idx];
            d_u_tmp[idx] = d_u[idx];
            d_v_tmp[idx] = d_v[idx];
            d_w_tmp[idx] = d_w[idx];
            d_p0[idx] = d_p[idx];
            d_T0[idx] = d_T[idx];
            d_T_tmp[idx] = d_T[idx];
            d_C0[idx] = d_C[idx];
            d_C_tmp[idx] = d_C[idx];
        }
        // boundary
#pragma acc parallel loop independent present(d_bList[:bsize_b], d_u[:bsize], d_v[:bsize], d_w[:bsize], d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize]) async
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t idx = d_bList[j];
            d_u0[idx] = d_u[idx];
            d_v0[idx] = d_v[idx];
            d_w0[idx] = d_w[idx];
            d_u_tmp[idx] = d_u[idx];
            d_v_tmp[idx] = d_v[idx];
            d_w_tmp[idx] = d_w[idx];
            d_p0[idx] = d_p[idx];
            d_T0[idx] = d_T[idx];
            d_T_tmp[idx] = d_T[idx];
            d_C0[idx] = d_C[idx];
            d_C_tmp[idx] = d_C[idx];
        }
        // obstacles
#pragma acc parallel loop independent present(d_oList[:bsize_o], d_u[:bsize], d_v[:bsize], d_w[:bsize], d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize]) async
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t idx = d_oList[j];
            d_u0[idx] = d_u[idx];
            d_v0[idx] = d_v[idx];
            d_w0[idx] = d_w[idx];
            d_u_tmp[idx] = d_u[idx];
            d_v_tmp[idx] = d_v[idx];
            d_w_tmp[idx] = d_w[idx];
            d_p0[idx] = d_p[idx];
            d_T0[idx] = d_T[idx];
            d_T_tmp[idx] = d_T[idx];
            d_C0[idx] = d_C[idx];
            d_C_tmp[idx] = d_C[idx];
        }

        if (sync) {
#pragma acc wait
        }
    } //end data region
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters such as momentum/temperature/concentration source functions
/// \param	t		time
/// \param  sync	synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SolverI::UpdateSources(real t, bool sync) {

    auto params = Parameters::getInstance();

// Momentum source
    if (m_string_solver == SolverTypes::NSSolver or \
        m_string_solver == SolverTypes::NSTempSolver or \
        m_string_solver == SolverTypes::NSTempConSolver or \
        m_string_solver == SolverTypes::NSTempTurbConSolver or \
        m_string_solver == SolverTypes::NSTempTurbSolver or \
        m_string_solver == SolverTypes::NSTurbSolver) {
        std::string forceFct = params->get("solver/source/force_fct");
        if (forceFct == SourceMethods::Zero or \
            forceFct == SourceMethods::Uniform) {

        } else if (forceFct == SourceMethods::Buoyancy) {
#ifndef PROFILING
            m_logger->info("Update f(T) ...");
#endif
            MomentumSource();
        } else {
            m_logger->critical("Source function not yet implemented! Simulation stopped!");
            std::exit(1);
            //TODO Error handling
        }
    }

// Temperature source
    if (m_string_solver == SolverTypes::NSTempSolver or \
        m_string_solver == SolverTypes::NSTempConSolver or \
        m_string_solver == SolverTypes::NSTempTurbConSolver or \
        m_string_solver == SolverTypes::NSTempTurbSolver) {

        std::string tempFct = params->get("solver/temperature/source/temp_fct");
        if (tempFct == SourceMethods::Zero) {

        } else if (tempFct == SourceMethods::BuoyancyST_MMS) {
            sou_vel->BuoyancyST_MMS(S_T, t, sync);
        } else if (tempFct == SourceMethods::GaussST and params->get("solver/temperature/source/ramp_fct") == FunctionNames::RampTanh) {
            // get parameters for Gauss function
            real HRR = params->getReal("solver/temperature/source/HRR");    // heat release rate in [kW]
            real cp = params->getReal("solver/temperature/source/cp");        // specific heat capacity in [kJ/ kg K]
            real x0 = params->getReal("solver/temperature/source/x0");
            real y0 = params->getReal("solver/temperature/source/y0");
            real z0 = params->getReal("solver/temperature/source/z0");
            real sigmax = params->getReal("solver/temperature/source/sigmax");
            real sigmay = params->getReal("solver/temperature/source/sigmay");
            real sigmaz = params->getReal("solver/temperature/source/sigmaz");

            sou_vel->Gauss(S_T, HRR, cp, x0, y0, z0, sigmax, sigmay, sigmaz, sync);

            auto d_S_T = S_T->data;
            auto bsize = Domain::getInstance()->GetSize();

            auto boundary = BoundaryController::getInstance();
            size_t *d_iList = boundary->get_innerList_level_joined();
            auto bsize_i = boundary->getSize_innerList();

            real t_ramp = Functions::RampTanh(t);

            auto d_T = T->data;

            // ramp-up
#pragma acc parallel loop independent present(d_iList[:bsize_i], d_S_T[:bsize]) async
            for (size_t l = 0; l < bsize_i; ++l) {
                const size_t idx = d_iList[l];
                d_S_T[idx] *= t_ramp;
            }
        } else {
            m_logger->critical("Source function not yet implemented! Simulation stopped!");
            std::exit(1);
            //TODO Error handling
        }
    }

// Concentration source
    if (m_string_solver == SolverTypes::NSTempConSolver or \
        m_string_solver == SolverTypes::NSTempTurbConSolver) {
        std::string conFct = params->get("solver/concentration/source/con_fct");
        if (conFct == SourceMethods::Zero) {

        } else if (conFct == SourceMethods::GaussSC \
 and params->get("solver/concentration/source/ramp_fct") == FunctionNames::RampTanh) {
            //get parameters for Gauss function
            real HRR = params->getReal("solver/concentration/source/HRR");       // heat release rate in [kW]
            real Hc = params->getReal("solver/concentration/source/Hc");        // heating value in [kJ/kg]
            real Ys = params->getReal("solver/concentration/source/Ys");        // soot yield in [g/g]
            real YsHRR = Ys * HRR;
            real x0 = params->getReal("solver/concentration/source/x0");
            real y0 = params->getReal("solver/concentration/source/y0");
            real z0 = params->getReal("solver/concentration/source/z0");
            real sigmax = params->getReal("solver/concentration/source/sigmax");
            real sigmay = params->getReal("solver/concentration/source/sigmay");
            real sigmaz = params->getReal("solver/concentration/source/sigmaz");

            sou_con->Gauss(S_C, YsHRR, Hc, x0, y0, z0, sigmax, sigmay, sigmaz, sync);
            auto d_S_C = S_C->data;
            auto bsize = Domain::getInstance()->GetSize();

            auto boundary = BoundaryController::getInstance();
            size_t *d_iList = boundary->get_innerList_level_joined();
            auto bsize_i = boundary->getSize_innerList();

            real t_ramp = Functions::RampTanh(t);

            // ramp-up
#pragma acc parallel loop independent present(d_iList[:bsize_i], d_S_C[:bsize]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                size_t idx = d_iList[j];
                d_S_C[idx] *= t_ramp;
            }
        } else {
            m_logger->critical("Source function not yet implemented! Simulation stopped!");
            std::exit(1);
            //TODO Error handling
        }
    }
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters temperature source functions
// ***************************************************************************************
void SolverI::TemperatureSource() {
// Temperature source
    if (Parameters::getInstance()->get("solver/temperature/source/temp_fct") == FunctionNames::BuoyancyST_MMS) {
        Functions::BuoyancyST_MMS(S_T, 0.);
    }
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters force source functions
// ***************************************************************************************
void SolverI::ForceSource() {
    auto params = Parameters::getInstance();
    // Force
    if (params->get("solver/source/force_fct") == SourceMethods::Buoyancy) {
        std::string dir = params->get("solver/source/dir");
        Init_f(T_a, T0->data);

        if (dir.find('x') != std::string::npos) {
            Functions::BuoyancyForce(f_x, T0, T_a);
        }
        if (dir.find('y') != std::string::npos) {
            Functions::BuoyancyForce(f_y, T0, T_a);
        }
        if (dir.find('z') != std::string::npos) {
            Functions::BuoyancyForce(f_z, T0, T_a);
        }
    }
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters momentum source functions
// ***************************************************************************************
void SolverI::MomentumSource() {
    //Momentum source
    auto params = Parameters::getInstance();
    std::string dir_vel = params->get("solver/source/dir");
    if (dir_vel.find('x') != std::string::npos) {
        sou_vel->BuoyancyForce(f_x, T, T_a);
    }
    if (dir_vel.find('y') != std::string::npos) {
        sou_vel->BuoyancyForce(f_y, T, T_a);
    }
    if (dir_vel.find('z') != std::string::npos) {
        sou_vel->BuoyancyForce(f_z, T, T_a);
    }
}
