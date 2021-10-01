
#include <gtest/gtest.h>

#include "src/randomField/UniformRandom.h"


#define EPS 10e-5


class UniformRandomFieldTest : public testing::Test {
};

TEST_F(UniformRandomFieldTest, size_test) {
    size_t size = 10;

    UniformRandom noise_maker(1.0, 1.0, 0);
    Field a = noise_maker.random_field(size);

    EXPECT_EQ(a.get_size(), 10);
}

TEST_F(UniformRandomFieldTest, size_test_2000) {
    size_t size = 2000;

    UniformRandom noise_maker(1.0, 1.0, 0);
    Field a = noise_maker.random_field(size);

    EXPECT_EQ(a.get_size(), 2000);
}

TEST_F(UniformRandomFieldTest, reprod) {
    int range = 10;
    size_t size = 2 * range * 100;

    UniformRandom noise_maker(range, 1.0, 0);
    Field a = noise_maker.random_field(size);

    UniformRandom noise_maker2(range, 1.0, 0);
    Field b = noise_maker2.random_field(size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], b[i]);
    }
}

TEST_F(UniformRandomFieldTest, reprod2) {
    int range = 10;
    size_t size = 2 * range * 100;

    UniformRandom noise_maker(range, 1.0, 1336);
    Field a = noise_maker.random_field(size);

    UniformRandom noise_maker2(range, 1.0, 1336);
    Field b = noise_maker2.random_field(size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], b[i]);
    }
}

TEST_F(UniformRandomFieldTest, uniq) {
    int range = 10;
    size_t size = 2 * range * 10000;

    UniformRandom noise_maker(range, 1.0, 0);
    Field a = noise_maker.random_field(size);

    UniformRandom noise_maker2(range, 1.0, 1);
    Field b = noise_maker2.random_field(size);

    int count = 0;
    for (auto i = 0; i < size; ++i) {
        count += abs(a[i] - b[i]) >= EPS;
    }

    // TODO calculate probability of failure
    EXPECT_GT(count, 1000);
}

TEST_F(UniformRandomFieldTest, uniq2) {
    int range = 10;
    size_t size = 2 * range * 10000;

    UniformRandom noise_maker(range, 1.0, 1337);
    Field a = noise_maker.random_field(size);

    UniformRandom noise_maker2(range, 1.0, 1338);
    Field b = noise_maker2.random_field(size);

    int count = 0;
    for (auto i = 0; i < size; ++i) {
        count += abs(a[i] - b[i]) >= EPS;
    }

    // TODO calculate probability of failure
    EXPECT_GT(count, 1000);
}

TEST_F(UniformRandomFieldTest, success_difference) {
    int range = 10;
    size_t size = 2 * range * 10000;

    UniformRandom noise_maker(range, 1.0, 1337);
    Field a = noise_maker.random_field(size);
    Field b = noise_maker.random_field(size);

    int count = 0;
    for (auto i = 0; i < size; ++i) {
        count += abs(a[i] - b[i]) >= EPS;
    }

    // TODO calculate probability of failure
    EXPECT_GT(count, 1000);
}

TEST_F(UniformRandomFieldTest, range_10_ge_le) {
    int range = 10;
    size_t size = 2 * range * 100;

    UniformRandom noise_maker(range, 1.0, 0);
    Field a = noise_maker.random_field(size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_GE(a[i], -range);
        EXPECT_LE(a[i], range);
    }
}

TEST_F(UniformRandomFieldTest, range_1000_ge_le) {
    int range = 1000;
    size_t size = 2 * range * 100;

    UniformRandom noise_maker(range, 1.0, 1337);
    Field a = noise_maker.random_field(size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_GE(a[i], -range);
        EXPECT_LE(a[i], range);
    }
}

TEST_F(UniformRandomFieldTest, range_10_at_least_one) {
    int range = 10;
    size_t size = 2 * range * 1000;
    real step_size = 1.0;

    UniformRandom noise_maker(range, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -range; j <= range; ++j) {
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (abs(a[i] - j) < EPS) {
                found = true;
                break;
            }
        }
        // could fail in 1 : 10E-45 Cases
        EXPECT_TRUE(found) << "didn't found: " << j;
    }
}

TEST_F(UniformRandomFieldTest, range_10_at_least_one2) {
    int range = 10;
    int step_size = 2;
    size_t size = 2 * range * 1000;

    UniformRandom noise_maker(range, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = 0; j <= 2 * range / step_size; ++j) {
        real no = -range + j * step_size;
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (abs(a[i] - no) < EPS) {
                found = true;
                break;
            }
        }
        // could fail in 1 : 10E-45 Cases
        EXPECT_TRUE(found) << "didn't found: " << no;
    }
}

TEST_F(UniformRandomFieldTest, range_10_at_least_one_h) {
    int range = 10;
    real step_size = 0.5;
    size_t size = 2 * range * 1000;

    UniformRandom noise_maker(range, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = 0; j <= 2 * range / step_size; ++j) {
        real no = -range + j * step_size;
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (abs(a[i] - no) < EPS) {
                found = true;
                break;
            }
        }
        // could fail in 1 : 10E-45 Cases
        EXPECT_TRUE(found) << "didn't found: " << no;
    }
}

TEST_F(UniformRandomFieldTest, range_100_at_least_one) {
    int range = 100;
    size_t size = 2 * range * 10000;
    real step_size = 1.0;

    UniformRandom noise_maker(range, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = 0; j <= 2 * range / step_size; ++j) {
        real no = -range + j * step_size;
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (abs(a[i] == no) < EPS) {
                found = true;
                break;
            }
        }
        // could fail in 1 : 10E-88 Cases
        EXPECT_TRUE(found) << "didn't found: " << no;
    }
}

TEST_F(UniformRandomFieldTest, range_100_at_least_one2) {
    int range = 100;
    int step_size = 2;
    size_t size = 2 * range * 10000;

    UniformRandom noise_maker(range, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = 0; j <= 2 * range / step_size; ++j) {
        real no = -range + j * step_size;
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (abs(a[i] - no) < EPS) {
                found = true;
                break;
            }
        }
        // could fail in 1 : 10E-45 Cases
        EXPECT_TRUE(found) << "didn't found: " << no;
    }
}

TEST_F(UniformRandomFieldTest, range_100_at_least_one_h) {
    int range = 100;
    real step_size = 0.5;
    size_t size = 2 * range * 10000;

    UniformRandom noise_maker(range, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = 0; j <= 2 * range / step_size; ++j) {
        real no = -range + j * step_size;
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (abs(a[i] - no) < EPS) {
                found = true;
                break;
            }
        }
        // could fail in 1 : 10E-45 cases
        EXPECT_TRUE(found) << "didn't found: " << no;
    }
}

TEST_F(UniformRandomFieldTest, stress_test) {
    int range = 100;
    int step_size = 2;
    size_t size = 10000;

    for (int i = 0; i < 10000; ++i) {
        UniformRandom noise_maker(range, step_size, 0);
        Field a = noise_maker.random_field(size);
        EXPECT_GE(a[a[0]], -range);
        EXPECT_LE(a[a[0]], range);
    }
}

TEST_F(UniformRandomFieldTest, stress_test2) {
    int range = 100;
    int step_size = 2;
    size_t size = 20000;

    for (int i = 0; i < 10000; ++i) {
        UniformRandom noise_maker(range, step_size, 0);
        Field a = noise_maker.random_field(size);
        EXPECT_GE(a[a[0]], -range);
        EXPECT_LE(a[a[0]], range);
    }
}

TEST_F(UniformRandomFieldTest, step_size) {
    int range = 10;
    real step_size = 0.1;
    size_t size = 1000;

    UniformRandom noise_maker(range, step_size, 0);
    Field a = noise_maker.random_field(size);
    for (int i = 0; i < size; ++i) {
        real tmp = a[i] * 10;
        int integer = static_cast<int>(tmp);
        ASSERT_LE(tmp - integer, EPS);
    }
}

TEST_F(UniformRandomFieldTest, negative_positive) {
    int range = 10;
    real step_size = 0.1;
    size_t size = 10000;

    UniformRandom noise_maker(range, step_size, 0);
    Field a = noise_maker.random_field(size);
    int counter_positive = 0;
    int counter_negative = 0;
    for (int i = 0; i < size; ++i) {
        if (a[i] < 0) {
            counter_positive++;
        } else if (a[i] > 0) {
            counter_negative++;
        }
    }
    ASSERT_GE(counter_positive, size * 0.4);
    ASSERT_GE(counter_negative, size * 0.4);
}
