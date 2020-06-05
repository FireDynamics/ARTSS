/// \file 		Analysis.cpp
/// \brief 		Calculates residual, compares analytical and numerical solutions, saves variables
/// \date 		July 11, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "Analysis.h"
#include "../boundary/BoundaryController.h"
#include "Solution.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../utility/Utility.h"

Analysis::Analysis() {
    m_logger = Utility::createLogger(typeid(this).name());
    auto params = Parameters::getInstance();
    hasAnalyticSolution = params->get("solver/solution/available") == "Yes";
    if (hasAnalyticSolution) {
        m_tol = params->getReal("solver/solution/tol");
    } else {
        m_logger->info("No analytical solution available!\n");
    }
}

// ===================================== Start analysis ==================================
// ***************************************************************************************
/// \brief  starts analysis to compare numerical and analytical solutions
/// \param	solver		pointer to solver
/// \param	t			current time
// ***************************************************************************************
void Analysis::Analyse(SolverI *solver, const real t) {
    m_logger = Utility::createLogger(typeid(this).name());
    //TODO statement t == 0.
    Solution solution;

    auto params = Parameters::getInstance();

    if (hasAnalyticSolution) {

        tinyxml2::XMLElement* rootElement = params->getRootElement();
        tinyxml2::XMLElement *xmlParameter = rootElement->FirstChildElement("boundaries");

        auto curElem = xmlParameter->FirstChildElement();

        solution.CalcAnalyticalSolution(t);
        m_logger->info("Compare to analytical solution:");

        while (curElem) {
            std::string nodeName(curElem->Value());

            if (nodeName == "boundary") {
                std::string field = curElem->Attribute("field");
                if (field.find(BoundaryData::getFieldTypeName(FieldType::U)) != std::string::npos) {
                    if (t == 0.) {
                        CompareSolutions(solver->GetU0(), solution.GetU0(), solver->u0->GetType(), 0.);
                    } else {
                        CompareSolutions(solver->GetU(), solution.GetU(), solver->u->GetType(), t);
                    }
                }
                if (field.find(BoundaryData::getFieldTypeName(FieldType::V)) != std::string::npos) {
                    if (t == 0.) {
                        CompareSolutions(solver->GetV0(), solution.GetV0(), solver->v0->GetType(), 0.);
                    } else {
                        CompareSolutions(solver->GetV(), solution.GetV(), solver->v->GetType(), t);
                    }
                }
                if (field.find(BoundaryData::getFieldTypeName(FieldType::W)) != std::string::npos) {
                    if (t == 0.) {
                        CompareSolutions(solver->GetW0(), solution.GetW0(), solver->w0->GetType(), 0.);
                    } else {
                        CompareSolutions(solver->GetW(), solution.GetW(), solver->w->GetType(), t);
                    }
                }
                if (field.find(BoundaryData::getFieldTypeName(FieldType::P)) != std::string::npos) {
                    if (t == 0.) {
                        CompareSolutions(solver->GetP0(), solution.GetP0(), solver->p0->GetType(), 0.);
                    } else {
                        CompareSolutions(solver->GetP(), solution.GetP(), solver->p->GetType(), t);
                    }
                }
                if (field.find(BoundaryData::getFieldTypeName(FieldType::T)) != std::string::npos) {
                    if (t == 0.) {
                        CompareSolutions(solver->GetT0(), solution.GetT0(), solver->T0->GetType(), 0.);
                    } else {
                        CompareSolutions(solver->GetT(), solution.GetT(), solver->T->GetType(), t);
                    }
                }
            }//end if
            curElem = curElem->NextSiblingElement();
        }//end while
    }
}

// ======================= Compare analytical and numerical solution =====================
// ***************************************************************************************
/// \brief  compares analytical solution and numerical solution, returns true when verification passed
/// \param	num		numerical solution
/// \param	ana		analytical solution
/// \param	type	type of variable
/// \param	t		current time
// ***************************************************************************************
bool Analysis::CompareSolutions(read_ptr num, read_ptr ana, const FieldType type, const real t) {

    bool verification = false;

    real res;
// Choose absolute or relative based error calculation
    res = CalcAbsoluteSpatialError(num, ana);
    //res = CalcRelativeSpatialError(num, ana);

    if (res <= m_tol) {
        m_logger->info("{} PASSED Test at time {} with error e = {}",
                BoundaryData::getFieldTypeName(type), t, res);
        verification = true;
    } else {
        m_logger->warn("{} FAILED Test at time {} with error e = {}",
                BoundaryData::getFieldTypeName(type), t, res);
    }
    return verification;
}

// ================================== Calculate absolute error ===========================
// ***************************************************************************************
/// \brief  calculates absolute spatial error based on L2-norm
/// \param	num		numerical solution
/// \param	ana		analytical solution
// ***************************************************************************************
real Analysis::CalcAbsoluteSpatialError(read_ptr num, read_ptr ana) {

    real sum = 0.;
    real r;

    auto boundary = BoundaryController::getInstance();
    size_t *innerList = boundary->get_innerList_level_joined();
    size_t size_iList = boundary->getSize_innerList();

    // weighted 2-norm
    // absolute error
    for (size_t i = 0; i < size_iList; i++) {
        size_t idx = innerList[i];
        r = fabs(num[idx] - ana[idx]);
        sum += r * r;
    }

    //weight
    real nr = static_cast<real>(size_iList);
    real eps = sqrt(1. / nr * sum);

    m_logger->info("Absolute error ||e|| = {}", eps);
    //std::cout << "num =" << num[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)] 		<< std::endl;
    //std::cout << "ana =" << ana[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)] 		<< std::endl;
    //std::cout << "num =" << num[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]	<< std::endl;
    //std::cout << "ana =" << ana[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]	<< std::endl;
    return eps;
}

// ================================== Calculate relative error ===========================
// ***************************************************************************************
/// \brief  calculates relative spatial error based on L2-norm
/// \param	num		numerical solution
/// \param	ana		analytical solution
// ***************************************************************************************
real Analysis::CalcRelativeSpatialError(read_ptr num, read_ptr ana) {

    real sumr = 0.;
    real rr;

    auto boundary = BoundaryController::getInstance();
    size_t *innerList = boundary->get_innerList_level_joined();
    size_t size_iList = boundary->getSize_innerList();

    // relative part with norm of analytical solution as denominator
    for (size_t i = 0; i < size_iList; i++) {
        rr = ana[innerList[i]];
        sumr += rr * rr;
    }

    //weight
    real nr = static_cast<real>(size_iList);
    real adenom = sqrt(1. / nr * sumr);

    real eps;
    real zero_tol = 10e-20;
    real epsa = CalcAbsoluteSpatialError(num, ana);

    // zero absolute error => zero relative error
    if (epsa <= zero_tol) {
        eps = 0.0;

        // zero denominator => take 2-norm of numerical solution as denominator
    } else if (adenom <= zero_tol) {
        sumr = 0.;

        // relative part with norm of numerical solution as quotient
        for (size_t i = 0; i < size_iList; i++) {
            rr = num[innerList[i]];
            sumr += rr * rr;
        }

        real ndenom = sqrt(1. / nr * sumr);

        eps = epsa / ndenom;

    } else {
        eps = epsa / adenom;
    }

    m_logger->info("Relative error ||e|| = {}", eps);
    /*std::cout << "num =" << num[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)] 		<< std::endl;
    std::cout << "ana =" << ana[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)] 		<< std::endl;
    std::cout << "num =" << num[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]	<< std::endl;
    std::cout << "ana =" << ana[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]	<< std::endl;*/

    return eps;
}

// ============= Calculate absolute error at center to be averaged over time ==============
// ***************************************************************************************
/// \brief  calculates absolute spatial error at time t at midpoint based on L2-norm
/// \param	solver		pointer to solver
/// \param	t			current time
/// \param	sum			pointer to sum for (u,p,T results)
// ***************************************************************************************
void Analysis::CalcL2NormMidPoint(SolverI *solver, real t, real *sum) {
    Solution solution;

    auto boundary = BoundaryController::getInstance();
    size_t *iList = boundary->get_innerList_level_joined();

    //take median of indices in iList to get center point ix
    //std::nth_element(iList.begin(), iList.begin() + iList.size()/2, iList.end());
    //size_t ix = iList[iList.size()/2];

    size_t ix = iList[boundary->getSize_innerList() / 2];

    auto params = Parameters::getInstance();
    if (params->get("solver/solution/available") == "Yes") {
        solution.CalcAnalyticalSolution(t);

        // local variables and parameters
        auto d_ua = solution.GetU();
        auto d_pa = solution.GetP();
        auto d_Ta = solution.GetT();

        auto d_u = solver->GetU();
        auto d_p = solver->GetP();
        auto d_T = solver->GetT();

        real ru = fabs((d_u[ix] - d_ua[ix]));
        real rp = fabs((d_p[ix] - d_pa[ix]));
        real rT = fabs((d_T[ix] - d_Ta[ix]));
        sum[0] += ru * ru;
        sum[1] += rp * rp;
        sum[2] += rT * rT;
    }
}

// ================================= Calculate RMS error ==================================
// ***************************************************************************************
/// \brief  calculates absolute spatial error at time t at midpoint based on L2-norm
/// \param	solver		pointer to solver
/// \param	t			current time
/// \param	sum			pointer to sum for (u,p,T results)
// ***************************************************************************************
void Analysis::CalcRMSError(real sumu, real sump, real sumT) {

    auto params = Parameters::getInstance();

    if (params->get("solver/solution/available") == "Yes") {
        // local variables and parameters
        real dt = params->getReal("physical_parameters/dt");
        real t_end = params->getReal("physical_parameters/t_end");
        auto Nt = static_cast<size_t>(std::round(t_end / dt));
        real rNt = 1. / static_cast<real>(Nt);
        real epsu = sqrt(rNt * sumu);
        real epsp = sqrt(rNt * sump);
        real epsT = sqrt(rNt * sumT);

        m_logger->info("RMS error of u at domain center is e_RMS = {}", epsu);
        m_logger->info("RMS error of p at domain center is e_RMS = {}", epsp);
        m_logger->info("RMS error of T at domain center is e_RMS = {}", epsT);
    }
}

// =============================== Check Von Neumann condition ===========================
// ***************************************************************************************
/// \brief  checks Von Neumann condition on time step (returns true or false)
/// \param	u			x-velocity field
/// \param	dt			time step size
// ***************************************************************************************
bool Analysis::CheckTimeStepVN(Field *u, real dt) {

    bool VN_check;

    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();

    // local variables and parameters
    real nu = params->getReal("physical_parameters/nu");

    real dx = domain->Getdx(u->GetLevel());
    real dy = domain->Getdy(u->GetLevel());
    real dz = domain->Getdz(u->GetLevel());

    real dx2sum = (dx * dx + dy * dy + dz * dz);
    real rdx2 = 1. / dx2sum;

    real VN = dt * nu * rdx2;

    VN_check = VN < 0.5;

    std::cout << "VN = " << VN << std::endl;

    return VN_check;
}

// ================================= Check CFL condition ==================================
// ***************************************************************************************
/// \brief  checks CFL condition on time step (returns true or false)
/// \param	u			x-velocity field
/// \param	v			y-velocity field
/// \param	w			z-velocity field
/// \param	dt			time step size
// ***************************************************************************************
bool Analysis::CheckTimeStepCFL(Field *u, Field *v, Field *w, real dt) {

    bool CFL_check;

    auto boundary = BoundaryController::getInstance();
    auto domain = Domain::getInstance();

    // local variables and parameters
    size_t *innerList = boundary->get_innerList_level_joined();
    size_t sizei = boundary->getSize_innerList();

    real dx = domain->Getdx(u->GetLevel());
    real dy = domain->Getdy(u->GetLevel());
    real dz = domain->Getdz(u->GetLevel());

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;

    real *max_vel = new real[sizei];
    real uvrdx, uvwrdx, maxvelrdx;

    //TODO correct?
    for (size_t i = 0; i < sizei; i++) {
        size_t idx = innerList[i];
        uvrdx = std::max(fabs(d_u[idx]) / dx, fabs(d_v[idx]) / dy);
        uvwrdx = std::max(uvrdx, fabs(d_w[idx]) / dz);

        max_vel[i] = uvwrdx;
        ++i;
    }

    maxvelrdx = *(std::max_element(max_vel, max_vel + sizei));

    real CFL = dt * maxvelrdx;

    CFL_check = CFL < 1.;

    m_logger->info("CFL = {}", CFL);

    return CFL_check;
}

// =============================== Set dt based on CFL condition ========================
// ***************************************************************************************
/// \brief  sets time step size based on CFL=0.8 (returns dt)
/// \param	u			x-velocity field
/// \param	v			y-velocity field
/// \param	w			z-velocity field
// ***************************************************************************************
real Analysis::SetDTwithCFL(Field *u, Field *v, Field *w) {

    auto boundary = BoundaryController::getInstance();
    auto domain = Domain::getInstance();

    // local variables and parameters
    size_t *innerList = boundary->get_innerList_level_joined();
    size_t sizei = boundary->getSize_innerList();

    real dx = domain->Getdx(u->GetLevel());
    real dy = domain->Getdy(u->GetLevel());
    real dz = domain->Getdz(u->GetLevel());

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;

    real *max_vel = new real[sizei];
    real uvrdx, uvwrdx, maxvelrdx;

    //TODO correct?
    for (size_t i = 0; i < sizei; i++) {
        size_t idx = innerList[i];
        uvrdx = std::max(fabs(d_u[idx]) / dx, fabs(d_v[idx]) / dy);
        uvwrdx = std::max(uvrdx, fabs(d_w[idx]) / dz);
        max_vel[i] = uvwrdx;
        ++i;
    }

    maxvelrdx = *(std::max_element(max_vel, max_vel + sizei));

    real CFL = 0.8;

    real DT = CFL / maxvelrdx;

    return DT;
}

// ==================================== Save variables ===================================
// ***************************************************************************************
/// \brief  saves variables in .dat files
/// \param	solv		pointer to solver
// ***************************************************************************************
void Analysis::SaveVariablesInFile(SolverI *solv) {
    //TODO do not write field out if not used
    auto boundary = BoundaryController::getInstance();
    size_t *innerList = boundary->get_innerList_level_joined();
    size_t size_innerList = boundary->getSize_innerList();
    size_t *boundaryList = boundary->get_boundaryList_level_joined();
    size_t size_boundaryList = boundary->getSize_boundaryList();
    size_t *obstacleList = boundary->get_obstacleList();
    size_t size_obstacleList = boundary->getSize_obstacleList();

    const real *dataField[numberOfFieldTypes-1];
    dataField[FieldType::U-1] = solv->GetU();
    dataField[FieldType::V-1] = solv->GetV();
    dataField[FieldType::W-1] = solv->GetW();
    dataField[FieldType::P-1] = solv->GetP();
    dataField[FieldType::T-1] = solv->GetT();

    for (size_t i = 0; i < numberOfFieldTypes -1; i++) {
        writeFile(dataField[i], BoundaryData::getFieldTypeName(static_cast<FieldType>(i+1)), innerList, size_innerList, boundaryList, size_boundaryList, obstacleList, size_obstacleList);
    }
}

void Analysis::writeFile(const real *field, std::string filename, size_t *innerList, size_t size_innerList, size_t *boundaryList, size_t size_boundaryList, size_t *obstacleList, size_t size_obstacleList) {

    std::ofstream out;
    out.open(filename + ".dat", std::ofstream::out);

    std::ofstream out_inner;
    out_inner.open(filename + "_inner.dat", std::ofstream::out);
    for (size_t idx = 0; idx < size_innerList; idx++) {
        out_inner << innerList[idx] << "|" << field[innerList[idx]] << std::endl;
        out << field[innerList[idx]] << std::endl;
    }
    out_inner.close();

    std::ofstream out_obstacle;
    out_obstacle.open(filename + "_obstacle.dat", std::ofstream::out);
    for (size_t idx = 0; idx < size_obstacleList; idx++) {
        out_obstacle << obstacleList[idx] << "|" << field[obstacleList[idx]] << std::endl;
        out << field[obstacleList[idx]] << std::endl;
    }
    out_obstacle.close();

    std::ofstream out_boundary;
    out_boundary.open(filename + "_boundary.dat", std::ofstream::out);
    for (size_t idx = 0; idx < size_boundaryList; idx++) {
        out_boundary << boundaryList[idx] << "|" << field[boundaryList[idx]] << std::endl;
        out << field[boundaryList[idx]] << std::endl;
    }
    out_boundary.close();

    out.close();
}
