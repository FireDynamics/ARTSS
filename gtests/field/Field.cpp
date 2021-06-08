
#include <iostream>
#include <filesystem>

#include <gtest/gtest.h>

#include "src/Domain.h"
#include "src/field/Field.h"
#include "src/utility/Parameters.h"
#include "src/utility/GlobalMacrosTypes.h"


class FieldTest : public testing::Test {
};

TEST_F(FieldTest, constructor_val) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.5, 0, size);

    for (auto i=0; i < size; ++i) {
        ASSERT_EQ(a.data[i], 0.0);
        ASSERT_EQ(b.data[i], 0.5);
    }
}

TEST_F(FieldTest, set_val) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);

    a.set_value(0.5);
    for (auto i=0; i < size; ++i) {
        ASSERT_EQ(a.data[i], 0.5);
    }
}

TEST_F(FieldTest, copy_data) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.5, 0, size);

    a.copy_data(b);

    for (auto i=0; i < size; ++i) {
        ASSERT_EQ(a.data[i], b.data[i]);
    }
}

TEST_F(FieldTest, swap_field) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i=0; i < size; ++i) {
        a.data[i] = x + 0.0;
        b.data[i] = x + 0.5;
        x += 1.0;
    }

    Field::swap(a, b);

    x = 0.0;
    for (auto i=0; i < size; ++i) {
        ASSERT_EQ(a.data[i], x + 0.5);
        ASSERT_EQ(b.data[i], x + 0.0);
        x += 1.0;
    }
}
