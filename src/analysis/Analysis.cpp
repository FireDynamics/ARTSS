/// \file       Analysis.cpp
/// \brief      Calculates residual, compares analytical and numerical
/// solutions, saves variables \date       Jul 11, 2016 \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

#include "../Domain.h"
#include "../boundary/BoundaryController.h"
#include "../utility/Parameters.h"
#include "../utility/Utility.h"
#include "Analysis.h"

Analysis::Analysis(Solution *solution) {
#ifndef BENCHMARKING
  m_logger = Utility::create_logger(typeid(this).name());
#endif
  auto params = Parameters::getInstance();
  has_analytic_solution = params->get("solver/solution/available") == "Yes";
  if (has_analytic_solution) {
    m_tol = params->get_real("solver/solution/tol");
  } else {
#ifndef BENCHMARKING
    m_logger->info("No analytical solution available!");
#endif
  }
  m_solution = solution;
}

// ================================ Start analysis =============================
// *****************************************************************************
/// \brief  starts analysis to compare numerical and analytical solutions
/// \param  field_controller    pointer to solver
/// \param  t     current time
// ***************************************************************************************
void Analysis::analyse(FieldController *field_controller, real t) {
  if (has_analytic_solution) {
    m_solution->calc_analytical_solution(t);

    auto params = Parameters::getInstance();

    tinyxml2::XMLElement *xmlParameter = params->get_first_child("boundaries");
    auto curElem = xmlParameter->FirstChildElement();

#ifndef BENCHMARKING
    m_logger->info("Compare to analytical solution:");
#endif

    while (curElem) {
      std::string nodeName(curElem->Value());

      if (nodeName == "boundary") {
        std::string field = curElem->Attribute("field");
        if (field.find(BoundaryData::getFieldTypeName(FieldType::U)) !=
            std::string::npos) {
          compare_solutions(field_controller->get_field_u_data(),
                            m_solution->GetU(), FieldType::U, t);
        }
        if (field.find(BoundaryData::getFieldTypeName(FieldType::V)) !=
            std::string::npos) {
          compare_solutions(field_controller->get_field_v_data(),
                            m_solution->GetV(), FieldType::V, t);
        }
        if (field.find(BoundaryData::getFieldTypeName(FieldType::W)) !=
            std::string::npos) {
          compare_solutions(field_controller->get_field_w_data(),
                            m_solution->GetW(), FieldType::W, t);
        }
        if (field.find(BoundaryData::getFieldTypeName(FieldType::P)) !=
            std::string::npos) {
          compare_solutions(field_controller->get_field_p_data(),
                            m_solution->GetP(), FieldType::P, t);
        }
        if (field.find(BoundaryData::getFieldTypeName(FieldType::T)) !=
            std::string::npos) {
          compare_solutions(field_controller->get_field_T_data(),
                            m_solution->GetT(), FieldType::T, t);
        }
      }  // end if
      curElem = curElem->NextSiblingElement();
    }  // end while
  }
}

// ================== Compare analytical and numerical solution ================
// *****************************************************************************
/// \brief  compares analytical solution and numerical solution, returns true
/// when verification passed \param  num   numerical solution \param  ana
/// analytical solution \param  type  type of variable \param  t   current time
// ***************************************************************************************
bool Analysis::compare_solutions(
        read_ptr num, read_ptr ana,
        FieldType type, real t) {
  bool verification = false;

  // Choose absolute or relative based error calculation
  real res = calc_absolute_spatial_error(num, ana);
  // real res = calc_relative_spatial_error(num, ana);

  if (res <= m_tol) {
#ifndef BENCHMARKING
    m_logger->info("{} PASSED Test at time {} with error e = {}",
                   BoundaryData::getFieldTypeName(type), t, res);
#endif
    verification = true;
  } else {
#ifndef BENCHMARKING
    m_logger->warn("{} FAILED Test at time {} with error e = {}",
                   BoundaryData::getFieldTypeName(type), t, res);
#endif
  }
  return verification;
}

// ============================= Calculate absolute error ======================
// *****************************************************************************
/// \brief  calculates absolute spatial error based on L2-norm
/// \param  num   numerical solution
/// \param  ana   analytical solution
// ***************************************************************************************
real Analysis::calc_absolute_spatial_error(read_ptr num, read_ptr ana) {
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

  // weight
  real nr = static_cast<real>(size_iList);
  real eps = sqrt(1. / nr * sum);

  // TODO(n16h7) add. output
#ifndef BENCHMARKING
  m_logger->info("Absolute error ||e|| = {}", eps);
#endif
  // std::cout << "num =" << num[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)] <<
  // std::endl; std::cout << "ana =" << ana[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx,
  // m_ny)]        << std::endl; std::cout << "num =" << num[IX((m_nx-2)/2 + 1,
  // (m_ny-2)/2, 1, m_nx, m_ny)]    << std::endl; std::cout << "ana =" <<
  // ana[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]    << std::endl;
  return eps;
}

// ============================= Calculate relative error ======================
// *****************************************************************************
/// \brief  calculates relative spatial error based on L2-norm
/// \param  num   numerical solution
/// \param  ana   analytical solution
// ***************************************************************************************
real Analysis::calc_relative_spatial_error(read_ptr num, read_ptr ana) {
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

  // weight
  real nr = static_cast<real>(size_iList);
  real adenom = sqrt(1. / nr * sumr);

  real eps;
  real zero_tol = 10e-20;
  real epsa = calc_absolute_spatial_error(num, ana);

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

  // TODO(n16h7) add. output
#ifndef BENCHMARKING
  m_logger->info("Relative error ||e|| = {}", eps);
#endif
  // std::cout << "num =" << num[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)] <<
  // std::endl; std::cout << "ana =" << ana[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx,
  // m_ny)]      << std::endl; std::cout << "num =" << num[IX((m_nx-2)/2 + 1,
  // (m_ny-2)/2, 1, m_nx, m_ny)]  << std::endl; std::cout << "ana =" <<
  // ana[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]  << std::endl;

  return eps;
}

// ======== Calculate absolute error at center to be averaged over time ========
// *****************************************************************************
/// \brief  calculates absolute spatial error at time t at midpoint based on
/// L2-norm \param  field_controller    pointer to field_controller \param  t
/// current time \param  sum     pointer to sum for (u,p,T results)
// ***************************************************************************************
void Analysis::calc_L2_norm_mid_point(FieldController *field_controller, real t,
                                      real *sum) {
  auto boundary = BoundaryController::getInstance();
  size_t *iList = boundary->get_innerList_level_joined();
  size_t ix = iList[boundary->getSize_innerList() / 2];

  if (has_analytic_solution) {
    m_solution->calc_analytical_solution(t);

    // local variables and parameters
    auto d_ua = m_solution->GetU();
    auto d_pa = m_solution->GetP();
    auto d_Ta = m_solution->GetT();

    auto d_u = field_controller->get_field_u_data();
    auto d_p = field_controller->get_field_p_data();
    auto d_T = field_controller->get_field_T_data();

    real ru = fabs((d_u[ix] - d_ua[ix]));
    real rp = fabs((d_p[ix] - d_pa[ix]));
    real rT = fabs((d_T[ix] - d_Ta[ix]));
    sum[0] += ru * ru;
    sum[1] += rp * rp;
    sum[2] += rT * rT;
  }
}

// =========================== Calculate RMS error ============================
// *****************************************************************************
/// \brief  calculates absolute spatial error at time t at midpoint based on
/// L2-norm \param  solver    pointer to solver \param  t     current time
/// \param  sum     pointer to sum for (u,p,T results)
// ***************************************************************************************
void Analysis::calc_RMS_error(real sum_u, real sum_p, real sum_T) {
  auto params = Parameters::getInstance();

  if (has_analytic_solution) {
    // local variables and parameters
    real dt = params->get_real("physical_parameters/dt");
    real t_end = params->get_real("physical_parameters/t_end");
    auto Nt = static_cast<size_t>(std::round(t_end / dt));
    real rNt = 1. / static_cast<real>(Nt);
    real epsu = sqrt(rNt * sum_u);
    real epsp = sqrt(rNt * sum_p);
    real epsT = sqrt(rNt * sum_T);

#ifndef BENCHMARKING
    m_logger->info("RMS error of u at domain center is e_RMS = {}", epsu);
    m_logger->info("RMS error of p at domain center is e_RMS = {}", epsp);
    m_logger->info("RMS error of T at domain center is e_RMS = {}", epsT);
#endif
  }
}

// ========================== Check Von Neumann condition ======================
// *****************************************************************************
/// \brief  checks Von Neumann condition on time step (returns true or false)
/// \param  u     x-velocity field
/// \param  dt      time step size
// ***************************************************************************************
bool Analysis::check_time_step_VN(Field *u, real dt) {
  bool VN_check;

  auto params = Parameters::getInstance();
  auto domain = Domain::getInstance();

  // local variables and parameters
  real nu = params->get_real("physical_parameters/nu");

  real dx = domain->get_dx(u->get_level());
  real dy = domain->get_dy(u->get_level());
  real dz = domain->get_dz(u->get_level());

  real dx2sum = (dx * dx + dy * dy + dz * dz);
  real rdx2 = 1. / dx2sum;

  real VN = dt * nu * rdx2;

  VN_check = VN < 0.5;

#ifndef BENCHMARKING
  m_logger->info("VN = {}", VN);
#endif

  return VN_check;
}

// ================================ Calc CFL ===================================
// *****************************************************************************
/// \brief  Returns the max CFL overall cells
/// \param  u     x-velocity field
/// \param  v     y-velocity field
/// \param  w     z-velocity field
/// \param  dt      time step size
// *****************************************************************************
real Analysis::calc_CFL(Field const &u, Field const &v, Field const &w, real dt) const {
  real cfl_max = 0;  // highest seen C. C is always positiv, so 0 is a lower
                     // bound
  real cfl_local;    // C in the local cell

  auto boundary = BoundaryController::getInstance();
  auto domain = Domain::getInstance();

  // local variables and parameters
  size_t *innerList = boundary->get_innerList_level_joined();
  size_t sizei = boundary->getSize_innerList();

  real dx = domain->get_dx();
  real dy = domain->get_dy();
  real dz = domain->get_dz();

  // calc C for every cell and get the maximum
#pragma acc data present(u, v, w)
#pragma acc parallel loop reduction(max : cfl_max)
  for (size_t i = 0; i < sizei; i++) {
    size_t idx = innerList[i];
    // \frac{C}{\Delta t} = \frac{\Delta u}{\Delta x} +
    //                      \frac{\Delta v}{\Delta y} +
    //                      \frac{\Delta w}{\Delta z} +
    // https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition#The_two_and_general_n-dimensional_case
    cfl_local = std::fabs(u[idx]) / dx + std::fabs(v[idx]) / dy +
                std::fabs(w[idx]) / dz;
    cfl_max = std::max(cfl_max, cfl_local);
  }

  return dt * cfl_max;
}

// =============================== Save variables ==============================
// *****************************************************************************
/// \brief  saves variables in .dat files
/// \param  field_controller    pointer to solver
// ***************************************************************************************
void Analysis::save_variables_in_file(FieldController *field_controller) {
  // TODO do not write field out if not used
  auto boundary = BoundaryController::getInstance();
  size_t *innerList = boundary->get_innerList_level_joined();
  size_t size_innerList = boundary->getSize_innerList();
  size_t *boundaryList = boundary->get_boundaryList_level_joined();
  size_t size_boundaryList = boundary->getSize_boundaryList();
  size_t *obstacleList = boundary->get_obstacleList();
  size_t size_obstacleList = boundary->getSize_obstacleList();

  std::vector<FieldType> v_fields = boundary->get_used_fields();

  const real *dataField[numberOfFieldTypes];
  dataField[FieldType::RHO] = field_controller->get_field_concentration_data();
  dataField[FieldType::U] = field_controller->get_field_u_data();
  dataField[FieldType::V] = field_controller->get_field_v_data();
  dataField[FieldType::W] = field_controller->get_field_w_data();
  dataField[FieldType::P] = field_controller->get_field_p_data();
  dataField[FieldType::T] = field_controller->get_field_T_data();

  for (auto &v_field : v_fields) {
    write_file(dataField[v_field], BoundaryData::getFieldTypeName(v_field),
               innerList, size_innerList, boundaryList, size_boundaryList,
               obstacleList, size_obstacleList);
  }
}

void Analysis::write_file(const real *field, const std::string &filename,
                          size_t *inner_list, size_t size_inner_list,
                          size_t *boundary_list, size_t size_boundary_list,
                          size_t *obstacle_list, size_t size_obstacle_list) {
  std::ofstream out;
  out.open(filename + ".dat", std::ofstream::out);

  std::ofstream out_inner;
  out_inner.open(filename + "_inner.dat", std::ofstream::out);
  for (size_t idx = 0; idx < size_inner_list; idx++) {
    out_inner << inner_list[idx] << ";" << field[inner_list[idx]] << std::endl;
    out << field[inner_list[idx]] << std::endl;
  }
  out_inner.close();

  std::ofstream out_obstacle;
  out_obstacle.open(filename + "_obstacle.dat", std::ofstream::out);
  for (size_t idx = 0; idx < size_obstacle_list; idx++) {
    out_obstacle << obstacle_list[idx] << ";" << field[obstacle_list[idx]]
                 << std::endl;
    out << field[obstacle_list[idx]] << std::endl;
  }
  out_obstacle.close();

  std::ofstream out_boundary;
  out_boundary.open(filename + "_boundary.dat", std::ofstream::out);
  for (size_t idx = 0; idx < size_boundary_list; idx++) {
    out_boundary << boundary_list[idx] << ";" << field[boundary_list[idx]]
                 << std::endl;
    out << field[boundary_list[idx]] << std::endl;
  }
  out_boundary.close();

  out.close();
}
