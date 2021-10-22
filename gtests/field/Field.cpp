
#include <iostream>
#include <filesystem>

#include <gtest/gtest.h>

#include "src/DomainData.h"
#include "src/field/Field.h"
#include "src/utility/Parameters.h"
#include "src/utility/GlobalMacrosTypes.h"


class FieldTest : public testing::Test {
};

TEST_F(FieldTest, constructor_val) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.5, 0, size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], 0.0);
        EXPECT_EQ(b[i], 0.5);
    }
}

TEST_F(FieldTest, set_val) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);

    a.set_value(0.5);
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], 0.5);
    }
}

TEST_F(FieldTest, stress_set_val) {
    size_t size = 100000;

    Field a(UNKNOWN_FIELD, 0.0, 0, size);

    for (int i = 0; i < 100000; ++i) {
        a.set_value(0.5);
    }

    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], 0.5);
    }
}


TEST_F(FieldTest, copy_data) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.5, 0, size);

    a.copy_data(b);

    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], b[i]);
    }
}

TEST_F(FieldTest, stress_copy_data) {
    size_t size = 100000;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.5, 0, size);

    for (int i = 0; i < 100000; ++i) {
        a.copy_data(b);
    }

    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], b[i]);
    }
}

TEST_F(FieldTest, swap_field) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x + 0.0;
        b[i] = x + 0.5;
        x += 1.0;
    }

    Field::swap(a, b);

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], x + 0.5);
        EXPECT_EQ(b[i], x + 0.0);
        x += 1.0;
    }
}

TEST_F(FieldTest, stress_swap_field) {
    size_t size = 100000;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x + 0.0;
        b[i] = x + 0.5;
        x += 1.0;
    }

    for (int i = 0; i < 100001; ++i) {
        Field::swap(a, b);
    }

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], x + 0.5);
        EXPECT_EQ(b[i], x + 0.0);
        x += 1.0;
    }
}

TEST_F(FieldTest, add_two_fields) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x + 0.0;
        b[i] = x + 0.5;
        x += 1.0;
    }

    a += b;

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], 2 * x + 0.5);
        x += 1.0;
    }
}

TEST_F(FieldTest, stress_add_two_fields) {
    size_t size = 100000;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 1.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x;
        x += 1.0;
    }

    for (int i = 0; i <= 100000; ++i) {
        a += b;
    }

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], x + 100001.);
        x += 1.0;
    }
}

TEST_F(FieldTest, mul_two_fields) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x + 0.0;
        b[i] = x + 0.5;
        x += 1.0;
    }

    a *= b;

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], b[i] * x);
        x += 1.0;
    }
}

TEST_F(FieldTest, stress_mul_two_fields) {
    size_t size = 100000;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);
    Field b(UNKNOWN_FIELD, 1.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x + 0.0;
        x += 1.0;
    }

    for (int i = 0; i < 100000; ++i) {
        a *= b;
    }

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], x);
        x += 1.0;
    }
}

TEST_F(FieldTest, add_scalar) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x;
        x += 1.0;
    }

    a += 0.5;

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], x + 0.5);
        x += 1.0;
    }
}

TEST_F(FieldTest, stress_add_scalar) {
    size_t size = 100000;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x;
        x += 1.0;
    }

    for (int i = 0; i < 100000; ++i) {
        a += 0.5;
    }

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], x + 100000. / 2);
        x += 1.0;
    }
}

TEST_F(FieldTest, mul_scalar) {
    size_t size = 100;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x;
        x += 1.0;
    }

    a *= 0.5;

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], x * 0.5);
        x += 1.0;
    }
}

TEST_F(FieldTest, stress_mul_scalar) {
    size_t size = 100000;
    Field a(UNKNOWN_FIELD, 0.0, 0, size);

    real x = 0.0;
    for (auto i = 0; i < size; ++i) {
        a[i] = x;
        x += 1.0;
    }

    for (int i = 0; i < 100000; ++i) {
        a *= 2.0;
        a *= 0.5;
    }

    x = 0.0;
    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], x);
        x += 1.0;
    }
}
