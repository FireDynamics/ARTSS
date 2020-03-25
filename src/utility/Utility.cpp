/// \file 		Utility.cpp
/// \brief 		Offers some tools
/// \date 		October 01, 2019
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cstring>
#include <iterator>
#include <sstream>
#include "Utility.h"
#include "GlobalMacrosTypes.h"

// ================================= Calc i,j,k ==========================================
// ***************************************************************************************
/// \brief  calculates indices <i,j,k> from linear index idx
/// \param  idx     linear (global) index
/// \param  Nx      number of cells in x-direction of physical domain
/// \param  Ny      number of cells in y-direction of physical domain
// ***************************************************************************************
std::vector<size_t> Utility::coordinateFromLinearIndex(size_t idx, size_t Nx, size_t Ny) {
    std::vector<size_t> coord;
    size_t k = getCoordinateK(idx, Nx, Ny);
    size_t j = getCoordinateJ(idx, Nx, Ny, k);
    size_t i = getCoordinateI(idx, Nx, Ny, j, k);

    coord.push_back(i);
    coord.push_back(j);
    coord.push_back(k);

    return coord;
}

// ================================= Calculate coordinate i ==========================================
// ***************************************************************************************
/// \brief  Calculates the i coordinate
/// \param  idx     linear (global) index
/// \param  Nx      number of cells in x-direction of physical domain
/// \param  Ny      number of cells in y-direction of physical domain
/// \param  j       index of <i,j,k>
/// \param  k       index of <i,j,k>
// ***************************************************************************************
//size_t Utility::getCoordinateI(size_t idx, size_t Nx, size_t Ny, size_t j, size_t k) {
//    return idx - k * Nx * Ny - j * Nx;
//}

// ================================= Calculate coordinate j ==========================================
// ***************************************************************************************
/// \brief  Calculates the j coordinate
/// \param  idx     linear (global) index
/// \param  Nx      number of cells in x-direction of physical domain
/// \param  Ny      number of cells in y-direction of physical domain
/// \param  k       index of <i,j,k>
// ***************************************************************************************
//size_t Utility::getCoordinateJ(size_t idx, size_t Nx, size_t Ny, size_t k) {
//    return (idx - k * Nx * Ny) / Nx;
//}

// ================================= Calculate coordinate k ==========================================
// ***************************************************************************************
/// \brief  Calculates the k coordinate
/// \param  idx     linear (global) index
/// \param  Nx      number of cells in x-direction of physical domain
/// \param  Ny      number of cells in y-direction of physical domain
// ***************************************************************************************
//size_t Utility::getCoordinateK(size_t idx, size_t Nx, size_t Ny) {
//    return idx / (Nx * Ny);
//}

// ================================= Split string at character ==========================================
// ***************************************************************************************
/// \brief  Splits a string at a defined char
/// \param  text     String
/// \param  delimiter Where to split
// ***************************************************************************************
std::vector<std::string> Utility::split(const std::string &text, char delimiter) {
    return split(text.c_str(), delimiter);
}

// ================================= Split string at character ==========================================
// ***************************************************************************************
/// \brief  Splits a string at a defined char
/// \param  text     String
/// \param  delimiter Where to split
// ***************************************************************************************
std::vector<std::string> Utility::split(const char *text, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(text);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
