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
    void merge_sort(Settings::Settings const &settings,
                    const size_t *list1, const size_t *list2,
                    const size_t size_list1, const size_t size_list2,
                    size_t *merged_list) {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(settings, class_name);
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
        std::vector <size_t> result;
        result.reserve(size_list1 + size_list2);
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
            if (list1[counter1] == result[result.size() - 1]) {
                counter1++;
                continue;
            }
            if (list2[counter2] == result[result.size() - 1]) {
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
        if (counter1 < size_list1) {
            if (list1[counter1] == result[result.size() - 1]) {
                counter1++;
            }
            for (size_t c = counter1; c < size_list1; c++) {
                result.push_back(list1[c]);
            }
        } else {
            if (list2[counter2] == result[result.size() - 1]) {
                counter2++;
            }
            for (size_t c = counter2; c < size_list2; c++) {
                result.push_back(list2[c]);
            }
        }
        result.shrink_to_fit();
        return result;
    }
}
