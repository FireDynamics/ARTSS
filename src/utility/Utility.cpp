/// \file      Utility.cpp
/// \brief     Offers some tools
/// \date      October 01, 2019
/// \author    My Linh Würzburger
/// \copyright <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

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

std::vector<size_t> Utility::mergeSortedListsToUniqueList(size_t *list1, size_t size_list1, size_t *list2, size_t size_list2) {
    std::vector<size_t> result;
    size_t counter1 = 0;
    size_t counter2 = 0;

    if (list1[counter1] < list2[counter2]) {
        result.push_back(list1[counter1]);
        counter1++;
    } else {
        result.push_back(list2[counter2]);
        counter2++;
    }
    while (counter1 < size_list1 && counter2 < size_list2) {
        if (list1[counter1] == result[result.size()-1]){
            counter1++;
            continue;
        }
        if (list2[counter2] == result[result.size()-1]){
            counter2++;
            continue;
        }
        if (list1[counter1] < list2[counter2]) {
            result.push_back(list1[counter1]);
            counter1++;
        } else {
            result.push_back(list2[counter2]);
            counter2++;
        }
    }
    if (counter1 < size_list1){
        for (size_t c = counter1; c < size_list1; c++){
            result.push_back(list1[c]);
        }
    }else{
        for (size_t c = counter2; c < size_list2; c++){
            result.push_back(list2[c]);
        }
    }
    return result;
}
