/// \file       PatchObject.cpp
/// \brief      
/// \date       Oct 27, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include <src/boundary/Coordinate.h>
#include "src/boundary/PatchObject.h"

TEST(PatchObjectTest, constructorVal) {
    PatchObject po;
    for (size_t i = 0; i < number_of_patches; i++) {
        EXPECT_EQ(po[i], 0);
    }
}

TEST(PatchObjectTest, constructorVal2) {
    PatchObject *po = new PatchObject();
    for (size_t i = 0; i < 6; i++) {
        EXPECT_EQ((*po)[i], 0);
    }
}

TEST(PatchObjectTest, sumInit) {
    PatchObject *po = new PatchObject();
    EXPECT_EQ(po->get_sum(), 0);
}

TEST(PatchObjectTest, sum) {
    PatchObject *po = new PatchObject();

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        po->add_value(patch, patch);
    }
    EXPECT_EQ(po->get_sum(), ((number_of_patches - 1) * number_of_patches) / 2);
}

TEST(PatchObjectTest, add) {
    PatchObject *po1 = new PatchObject();
    PatchObject *po2 = new PatchObject();

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        po1->add_value(patch, patch);
    }
    po2->add_value(Patch::FRONT, 6);
    po2->add_value(Patch::BOTTOM, 3);

    *po1 += *po2;

    EXPECT_EQ((*po1)[Patch::LEFT], 0);
    EXPECT_EQ((*po1)[Patch::RIGHT], 1);
    EXPECT_EQ((*po1)[Patch::BOTTOM], 5);
    EXPECT_EQ((*po1)[Patch::TOP], 3);
    EXPECT_EQ((*po1)[Patch::FRONT], 10);
    EXPECT_EQ((*po1)[Patch::BACK], 5);

    EXPECT_EQ((*po2)[Patch::LEFT], 0);
    EXPECT_EQ((*po2)[Patch::RIGHT], 0);
    EXPECT_EQ((*po2)[Patch::BOTTOM], 3);
    EXPECT_EQ((*po2)[Patch::TOP], 0);
    EXPECT_EQ((*po2)[Patch::FRONT], 6);
    EXPECT_EQ((*po2)[Patch::BACK], 0);
}
