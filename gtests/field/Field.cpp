
#include <iostream>
#include <filesystem>

#include <gtest/gtest.h>

#include "src/Domain.h"
#include "src/field/Field.h"
#include "src/utility/Parameters.h"
#include "src/utility/GlobalMacrosTypes.h"


class FieldTest : public testing::Test {
 protected:
};

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

    Field::swap(&a, &b);

    x = 0.0;
    for (auto i=0; i < size; ++i) {
        ASSERT_EQ(a.data[i], x + 0.5);
        ASSERT_EQ(b.data[i], x + 0.0);
        x += 1.0;
    }
}
