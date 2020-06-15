/// \file       Domain.cpp
/// \brief      XML Domain parameters to variables
/// \date       July 16, 2018
/// \author   My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include <spdlog/spdlog.h>

#include <iostream>
#include <spdlog/spdlog.h>

#include "Domain.h"
#include "utility/Parameters.h"

Domain *Domain::single = nullptr; //Singleton

Domain::Domain() {
    auto params = Parameters::getInstance();
    auto solver = params->get("solver/description");
    if ( solver.find("NS") != std::string::npos ||  solver.find("Pressure") != std::string::npos){
      m_levels = static_cast<size_t> (params->getInt("solver/pressure/n_level"));
    }
    m_nx = new size_t[m_levels + 1];
    m_ny = new size_t[m_levels + 1];
    m_nz = new size_t[m_levels + 1];

    m_nx[0] = static_cast<size_t> (params->getInt("domain_parameters/nx")) + 2;
    m_ny[0] = static_cast<size_t> (params->getInt("domain_parameters/ny")) + 2;
    m_nz[0] = static_cast<size_t> (params->getInt("domain_parameters/nz")) + 2;

    m_x1 = params->getReal("domain_parameters/x1");
    m_x2 = params->getReal("domain_parameters/x2");
    m_y1 = params->getReal("domain_parameters/y1");
    m_y2 = params->getReal("domain_parameters/y2");
    m_z1 = params->getReal("domain_parameters/z1");
    m_z2 = params->getReal("domain_parameters/z2");
    m_X1 = params->getReal("domain_parameters/X1");
    m_X2 = params->getReal("domain_parameters/X2");
    m_Y1 = params->getReal("domain_parameters/Y1");
    m_Y2 = params->getReal("domain_parameters/Y2");
    m_Z1 = params->getReal("domain_parameters/Z1");
    m_Z2 = params->getReal("domain_parameters/Z2");

    calcMGValues();
#ifndef PROFILING
    //printDetails();
#endif
}

// =============================== Calculation of MultiGrid Values ========================
// ***************************************************************************************
/// \brief  Calculates amount of Cells in XYZ direction for each multigrid level
// ***************************************************************************************
void Domain::calcMGValues() {

    for (size_t l = 1; l < m_levels + 1; ++l) {
        m_nx[l] = (m_nx[l - 1] == 3) ? 3 : static_cast<size_t> (std::round((m_nx[l - 1] - 2) / 2 + 2));
        m_ny[l] = (m_ny[l - 1] == 3) ? 3 : static_cast<size_t> (std::round((m_ny[l - 1] - 2) / 2 + 2));
        m_nz[l] = (m_nz[l - 1] == 3) ? 3 : static_cast<size_t> (std::round((m_nz[l - 1] - 2) / 2 + 2));
    }
}

Domain *Domain::getInstance() {
    if (single == nullptr) {
        single = new Domain();
    }
    return single;
}

// =============================== Resize computational Domain ========================
// ***************************************************************************************
/// \brief  calculates the size of the new calculation domain
/// \params shift_x1  Amount of cells which will be resized in x direction
/// \params shift_x2  Amount of cells which will be resized in x direction
/// \params shift_y1  Amount of cells which will be resized in y direction
/// \params shift_y2  Amount of cells which will be resized in y direction
/// \params shift_z1  Amount of cells which will be resized in z direction
/// \params shift_z2  Amount of cells which will be resized in z direction
/// \return bool True if computational domain has been resized, False if not
// ***************************************************************************************
bool Domain::Resize(long shift_x1, long shift_x2, long shift_y1, long shift_y2, long shift_z1, long shift_z2) {
    auto dx = Getdx();
    auto dy = Getdy();
    auto dz = Getdz();

    bool update = false;
    real tmp;
    if (setNewValue(shift_x1, m_X1, m_X2, m_x1, dx, &tmp)) {
      m_x1 = tmp;
      update = true;
    }
    if (setNewValue(shift_x2, m_X1, m_X2, m_x2, dx, &tmp)) {
      m_x2 = tmp;
      update = true;
    }
    if (setNewValue(shift_y1, m_Y1, m_Y2, m_y1, dy, &tmp)) {
      m_y1 = tmp;
      update = true;
    }
    if (setNewValue(shift_y2, m_Y1, m_Y2, m_y2, dy, &tmp)) {
      m_y2 = tmp;
      update = true;
    }
    if (setNewValue(shift_z1, m_Z1, m_Z2, m_z1, dz, &tmp)) {
      m_z1 = tmp;
      update = true;
    }
    if (setNewValue(shift_z2, m_Z1, m_Z2, m_z2, dz, &tmp)) {
      m_z2 = tmp;
      update = true;
    }
#pragma acc wait
    if (update) {
#ifndef PROFILING
        spdlog::info("Resize domain: {}|{} {}|{} {}|{}", shift_x1, shift_x2,
                                                         shift_y1, shift_y2,
                                                         shift_z1, shift_z2);
#endif
      m_nx[0] = static_cast<size_t> (std::round(Getlx() / dx + 2));
      m_ny[0] = static_cast<size_t> (std::round(Getly() / dy + 2));
      m_nz[0] = static_cast<size_t> (std::round(Getlz() / dz + 2));

      calcMGValues();
    }
    return update;
}

// =============================== Calculation of new Coordinates  ========================
// ***************************************************************************************
/// \brief Calculates new coordinates depending on shift variable
/// \params oldCoord  coordinate of old computational domain
/// \params shift  Amount of cells which will be resized
/// \params cellWidth  dx,dy,dz
/// \return coordinate of new computational domain
// ***************************************************************************************
real Domain::calcNewCoord(real oldCoord, long shift, real cellwidth) {
    return oldCoord + static_cast<real>(shift) * cellwidth;
}

// =============================== Setting coordinates of new computational domain  ========================
// ***************************************************************************************
/// \brief If resizing is valid set new coordinates
/// \params shift  Amount of cells which will be resized
/// \params startCoord_p Start coordinate of physical domain
/// \params endCoord_p End coordinate of physical domain
/// \params oldCoord Old coordiante of computational domain
/// \params cellWidth  dx,dy,dz
/// \params newCoord New coordiante of computational domain
/// \return bool True if resize is valid
// ***************************************************************************************
bool Domain::setNewValue(long shift, real startCoord_p, real endCoord_p, real oldCoord, real cellwidth, real *newCoord) {
    bool update = false;
    if (shift != 0) {
        *newCoord = calcNewCoord(oldCoord, shift, cellwidth);
        if (*newCoord >= startCoord_p && *newCoord <= endCoord_p) {
            update = true;
        }
    }
    return update;
}
void Domain::print(){
    std::cout << "-- Domain" << std::endl;
    std::cout << "\t Domain size inner cells: (" << Getnx()-2 << "|" << Getny()-2 << "|" << Getnz()-2 << ")" << std::endl;
    std::cout << "\t step size (x|y|z): (" << Getdx() << "|" << Getdy() << "|" << Getdz() << ")" << std::endl;
}

void Domain::printDetails(){
  std::cout << "############### Domain Parameter ###############" << std::endl;
  for (size_t level = 0; level < m_levels+1; level++) {
    std::cout << "For Level " << level << " Nx: " << GetNx(level) << " Ny: " << GetNy(level) << " Nz: " << GetNz(level) << std::endl;
  }
  std::cout << std::endl;
  for (size_t level = 0; level < m_levels+1; level++) {
    std::cout << "For Level " << level << " nx: " << Getnx(level) << " ny: " << Getny(level) << " nz: " << Getnz(level) << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X: (" << GetX1() << "|" << GetX2() << ") x: (" << Getx1() << "|" << Getx2() << ")"  <<std::endl;
  std::cout << "Y: (" << GetY1() << "|" << GetY2() << ") y: (" << Gety1() << "|" << Gety2() << ")"  <<std::endl;
  std::cout << "Z: (" << GetZ1() << "|" << GetZ2() << ") z: (" << Getz1() << "|" << Getz2() << ")"  <<std::endl;

  std::cout << std::endl;
  std::cout << "Lx: " << GetLx() << " Ly:" << GetLy() << " Lz: " << GetLz() << std::endl;
  std::cout << "lx: " << Getlx() << " ly:" << Getly() << " lz: " << Getlz() << std::endl;
  std::cout << std::endl;
  for (size_t level = 0; level < m_levels+1; level++) {
    std::cout << "For Level " << level << " dx: " << Getdx(level) << " dy:" << Getdy(level) << " dz: " << Getdz(level) << std::endl;
  }
  std::cout << std::endl;
  for (size_t level = 0; level < m_levels+1; level++) {
    std::cout << "for Level " << level << " X: (" << GetIndexx1(level) << "|" << GetIndexx2(level) << ") Y: (" << GetIndexy1(level) << "|" << GetIndexy2(level) << ") Z: (" << GetIndexz1(level) << "|" << GetIndexz2(level) << ")" <<std::endl;
  }
  for (size_t level = 0; level < m_levels+1; level++) {
    std::cout << " For Level "  << level << " domain size: " << GetSize(level) << std::endl;
  }
  std::cout << "--------------- Domain Parameter end ---------------" << std::endl;
}
