/// \file       PatchObject.cpp
/// \brief      
/// \date       Oct 27, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include <src/boundary/Coordinate.h>
#include "src/boundary/PatchObject.h"

class PatchObjectTest : public testing::Test {
};

TEST_F(PatchObjectTest, constructor_val) {
    PatchObject *po = new PatchObject();
    for (size_t i = 0; i < 6; i++) {
        ASSERT_EQ((*po)[0], 0);
    }
}

TEST_F(PatchObjectTest, sum_init) {
    PatchObject *po = new PatchObject();
    ASSERT_EQ(po->get_sum(), 0);
}

TEST_F(PatchObjectTest, match_patch) {
    Patch patch_front = PatchObject::match_patch("front");
    ASSERT_EQ(patch_front, Patch::FRONT);
    Patch patch_back = PatchObject::match_patch("back");
    ASSERT_EQ(patch_back, Patch::BACK);
    Patch patch_bottom = PatchObject::match_patch("bottom");
    ASSERT_EQ(patch_bottom, Patch::BOTTOM);
    Patch patch_top = PatchObject::match_patch("top");
    ASSERT_EQ(patch_top, Patch::TOP);
    Patch patch_left = PatchObject::match_patch("left");
    ASSERT_EQ(patch_left, Patch::LEFT);
    Patch patch_right = PatchObject::match_patch("right");
    ASSERT_EQ(patch_right, Patch::RIGHT);
    Patch patch_unknown = PatchObject::match_patch("a");
    ASSERT_EQ(patch_unknown, Patch::UNKNOWN_PATCH);
}

TEST_F(PatchObjectTest, get_boundary_patch_name) {
    std::string front = PatchObject::get_patch_name(Patch::FRONT);
    ASSERT_EQ("front", front);
    std::string back = PatchObject::get_patch_name(Patch::BACK);
    ASSERT_EQ("back", back);
    std::string bottom = PatchObject::get_patch_name(Patch::BOTTOM);
    ASSERT_EQ("bottom", bottom);
    std::string top = PatchObject::get_patch_name(Patch::TOP);
    ASSERT_EQ("top", top);
    std::string left = PatchObject::get_patch_name(Patch::LEFT);
    ASSERT_EQ("left", left);
    std::string right = PatchObject::get_patch_name(Patch::RIGHT);
    ASSERT_EQ("right", right);
}

TEST_F(PatchObjectTest, sum) {
    PatchObject *po = new PatchObject();

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        po->add_value(patch, patch);
    }
    ASSERT_EQ(po->get_sum(), ((number_of_patches - 1) * number_of_patches) / 2);
}

TEST_F(PatchObjectTest, add) {
    PatchObject *po1 = new PatchObject();
    PatchObject *po2 = new PatchObject();

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        po1->add_value(patch, patch);
    }
    po2->add_value(Patch::LEFT, 6);
    po2->add_value(Patch::BOTTOM, 3);

    *po1 += *po2;

    ASSERT_EQ((*po1)[Patch::FRONT], 0);
    ASSERT_EQ((*po1)[Patch::BACK], 1);
    ASSERT_EQ((*po1)[Patch::BOTTOM], 5);
    ASSERT_EQ((*po1)[Patch::TOP], 3);
    ASSERT_EQ((*po1)[Patch::LEFT], 10);
    ASSERT_EQ((*po1)[Patch::RIGHT], 5);

    ASSERT_EQ((*po2)[Patch::FRONT], 0);
    ASSERT_EQ((*po2)[Patch::BACK], 0);
    ASSERT_EQ((*po2)[Patch::BOTTOM], 3);
    ASSERT_EQ((*po2)[Patch::TOP], 0);
    ASSERT_EQ((*po2)[Patch::LEFT], 6);
    ASSERT_EQ((*po2)[Patch::RIGHT], 0);
}

TEST_F(PatchObjectTest, coordinate_axis) {
    ASSERT_EQ(CoordinateAxis::X * 2, Patch::FRONT);
    ASSERT_EQ(CoordinateAxis::X * 2 + 1, Patch::BACK);
    ASSERT_EQ(CoordinateAxis::Y * 2, Patch::BOTTOM);
    ASSERT_EQ(CoordinateAxis::Y * 2 + 1, Patch::TOP);
    ASSERT_EQ(CoordinateAxis::Z * 2, Patch::LEFT);
    ASSERT_EQ(CoordinateAxis::Z * 2 + 1, Patch::RIGHT);
}