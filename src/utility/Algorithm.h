/// \file       Algorithm.h
/// \brief      
/// \date       Oct 26, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_ALGORITHM_H_
#define ARTSS_ALGORITHM_H_

#include <string>
#include <vector>
#include "settings/Settings.h"

namespace Algorithm {
    std::vector<size_t> merge_sort(const size_t *list1, const size_t *list2, size_t size_list1, size_t size_list2);
    std::vector<size_t> merge_sort_with_duplicates(const size_t *list1, size_t size_list1, const size_t *list2, size_t size_list2);
};


#endif /* ARTSS_ALGORITHM_H_ */
