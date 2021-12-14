/// \file       SingleJoinedList.cpp
/// \brief      
/// \date       Oct 27, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include "src/joinedLists/SingleJoinedList.h"


class SingleJoinedListTest : public testing::Test {
};

TEST_F(SingleJoinedListTest, constructor_level0) {
#define BENCHMARKING
    size_t mg_level = 0;
    SingleJoinedList *jl = new SingleJoinedList(mg_level);
    for (size_t level = 0; level < mg_level + 1; level++) {
        ASSERT_EQ(jl->get_first_index(level), 0);
        ASSERT_EQ(jl->get_last_index(level), 0);
        ASSERT_EQ(jl->get_slice_size(level), 0);
    }
    ASSERT_EQ(jl->get_first_index(mg_level + 1), 0);
    ASSERT_EQ(jl->get_size(), 0);
#undef BENCHMARKING
}

TEST_F(SingleJoinedListTest, constructor_level2) {
    size_t mg_level = 2;
    SingleJoinedList *jl = new SingleJoinedList(mg_level);
    for (size_t level = 0; level < mg_level + 1; level++) {
        ASSERT_EQ(jl->get_first_index(level), 0);
        ASSERT_EQ(jl->get_last_index(level), 0);
        ASSERT_EQ(jl->get_slice_size(level), 0);
    }
    ASSERT_EQ(jl->get_first_index(mg_level + 1), 0);
    ASSERT_EQ(jl->get_size(), 0);
}