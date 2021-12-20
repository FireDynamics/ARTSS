/// \file       Algorithm.cpp
/// \brief      
/// \date       Oct 26, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "Algorithm.h"
#include "Utility.h"

namespace Algorithm {
static std::string class_name = "Algorithm";

//======================================== merge sort ==============================================
// *************************************************************************************************
/// \brief merge two !!sorted!! lists. merge sort algorithm
/// \param list1
/// \param list2
/// \param size_list1 size of list1
/// \param size_list2 size of list2
/// \param merged_list result of merge sort. size should be size_list1 + size_list2
// *************************************************************************************************
std::vector<size_t> merge_sort(const size_t *list1, const size_t *list2,
                               const size_t size_list1, const size_t size_list2) {
    std::vector<size_t> merged_list;
    merged_list.resize(size_list1 + size_list2);
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(class_name);
#endif
    size_t counter_list1 = 0, counter_list2 = 0, counter_merged_list = 0;
    while (counter_list1 < size_list1 && counter_list2 < size_list2) {
        if (list1[counter_list1] < list2[counter_list2]) {
            merged_list[counter_merged_list++] = list1[counter_list1++];
        } else {
            merged_list[counter_merged_list++] = list2[counter_list2++];
        }
    }
#ifndef BENCHMARKING
    logger->debug("merge sort: end of one list counter_list1: {}|{} counter_list2: {}|{}",
                    counter_list1, size_list1, counter_list2, size_list2);
#endif

    while (counter_list1 < size_list1) {
        merged_list[counter_merged_list++] = list1[counter_list1++];
    }

    while (counter_list2 < size_list2) {
        merged_list[counter_merged_list++] = list2[counter_list2++];
    }
#ifndef BENCHMARKING
    logger->debug("merge sort: end of merge sort counter_list1: {}|{} counter_list2: {}|{} counter_merged_list: {}|{}",
                    counter_list1, size_list1, counter_list2, size_list2, counter_merged_list, size_list1 + size_list2);
#endif
    return merged_list;
}

//======================================== merge sort ==============================================
// *************************************************************************************************
/// \brief merge two !!sorted!! lists with duplicates. merge sort algorithm
/// \param list1
/// \param list2
/// \param size_list1 size of list1
/// \param size_list2 size of list2
/// \return vector with merged list1 and list2. entries are unique
// *************************************************************************************************
std::vector<size_t> merge_sort_with_duplicates(const size_t *list1, const size_t size_list1,
                                                const size_t *list2, const size_t size_list2) {
    std::vector<size_t> result;
    result.reserve(size_list1 + size_list2);
    size_t counter_list1 = 0;
    size_t counter_list2 = 0;

    if (list1[counter_list1] < list2[counter_list2]) {
        result.push_back(list1[counter_list1]);
        counter_list1++;
    } else if (list1[counter_list1] > list2[counter_list2]){
        result.push_back(list2[counter_list2]);
        counter_list2++;
    } else {
        result.push_back(list1[counter_list1]);
        counter_list1++;
        counter_list2++;
    }
    while (counter_list1 < size_list1 && counter_list2 < size_list2) {
        if (list1[counter_list1] < list2[counter_list2]) {
            result.push_back(list1[counter_list1]);
            counter_list1++;
        } else if (list1[counter_list1] > list2[counter_list2]){
            result.push_back(list2[counter_list2]);
            counter_list2++;
        } else {
            result.push_back(list1[counter_list1]);
            counter_list1++;
            counter_list2++;
        }
    }
    if (counter_list1 < size_list1) {
        if (list1[counter_list1] == result[result.size() - 1]) {
            counter_list1++;
        }
        for (size_t c = counter_list1; c < size_list1; c++) {
            result.push_back(list1[c]);
        }
    } else if (counter_list2 < size_list2) {
        if (list2[counter_list2] == result[result.size() - 1]) {
            counter_list2++;
        }
        for (size_t c = counter_list2; c < size_list2; c++) {
            result.push_back(list2[c]);
        }
    }
    result.shrink_to_fit();
    return result;
}
}
