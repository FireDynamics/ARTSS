
#include <gtest/gtest.h>

#include "src/randomField/UniformRandom.h"


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
    int steps = 10;
    size_t size = 2 * steps * 100;

    UniformRandom noise_maker(steps, 1.0, 0);
    Field a = noise_maker.random_field(size);

    UniformRandom noise_maker2(steps, 1.0, 0);
    Field b = noise_maker2.random_field(size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], b[i]);
    }
}

TEST_F(UniformRandomFieldTest, reprod2) {
    int steps = 10;
    size_t size = 2 * steps * 100;

    UniformRandom noise_maker(steps, 1.0, 1336);
    Field a = noise_maker.random_field(size);

    UniformRandom noise_maker2(steps, 1.0, 1336);
    Field b = noise_maker2.random_field(size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_EQ(a[i], b[i]);
    }
}

TEST_F(UniformRandomFieldTest, uniq) {
    int steps = 10;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, 1.0, 0);
    Field a = noise_maker.random_field(size);

    UniformRandom noise_maker2(steps, 1.0, 1);
    Field b = noise_maker2.random_field(size);

    int count = 0;
    for (auto i = 0; i < size; ++i) {
        count += a[i] == b[i];
    }

    // TODO calculate probability of failure
    EXPECT_GT(count, 1000);
}

TEST_F(UniformRandomFieldTest, uniq2) {
    int steps = 10;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, 1.0, 1337);
    Field a = noise_maker.random_field(size);

    UniformRandom noise_maker2(steps, 1.0, 1338);
    Field b = noise_maker2.random_field(size);

    int count = 0;
    for (auto i = 0; i < size; ++i) {
        count += a[i] == b[i];
    }

    // TODO calculate probability of failure
    EXPECT_GT(count, 1000);
}

TEST_F(UniformRandomFieldTest, success_difference) {
    int steps = 10;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, 1.0, 1337);
    Field a = noise_maker.random_field(size);
    Field b = noise_maker.random_field(size);

    int count = 0;
    for (auto i = 0; i < size; ++i) {
        count += a[i] == b[i];
    }

    // TODO calculate probability of failure
    EXPECT_GT(count, 1000);
}

TEST_F(UniformRandomFieldTest, steps_10_ge_le) {
    int steps = 10;
    size_t size = 2 * steps * 100;

    UniformRandom noise_maker(steps, 1.0, 0);
    Field a = noise_maker.random_field(size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_GE(a[i], -steps);
        EXPECT_LE(a[i], steps);
    }
}

TEST_F(UniformRandomFieldTest, steps_1000_ge_le) {
    int steps = 1000;
    size_t size = 2 * steps * 100;

    UniformRandom noise_maker(steps, 1.0, 1337);
    Field a = noise_maker.random_field(size);

    for (auto i = 0; i < size; ++i) {
        EXPECT_GE(a[i], -steps);
        EXPECT_LE(a[i], steps);
    }
}

TEST_F(UniformRandomFieldTest, steps_10_at_least_one) {
    int steps = 10;
    size_t size = 2 * steps * 1000;

    UniformRandom noise_maker(steps, 1.0, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps; j <= steps; ++j) {
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (a[i] == j) {
                found = true;
            }
        }
        // could fail in 1 : 10E-45 Cases
        EXPECT_TRUE(found) << "didn't found: " << j;
    }
}

TEST_F(UniformRandomFieldTest, steps_10_at_least_one2) {
    int steps = 10;
    int step_size = 2;
    size_t size = 2 * steps * 1000;

    UniformRandom noise_maker(steps, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps / step_size; j <= steps / step_size; ++j) {
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (a[i] == j) {
                found = true;
            }
        }
        // could fail in 1 : 10E-45 Cases
        EXPECT_TRUE(found) << "didn't found: " << j;
    }
}

TEST_F(UniformRandomFieldTest, steps_10_at_least_one_h) {
    int steps = 10;
    real step_size = 0.5;
    size_t size = 2 * steps * 1000;

    UniformRandom noise_maker(steps, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps / step_size; j <= steps / step_size; ++j) {
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (a[i] == j) {
                found = true;
            }
        }
        // could fail in 1 : 10E-45 Cases
        EXPECT_TRUE(found) << "didn't found: " << j;
    }
}

TEST_F(UniformRandomFieldTest, steps_100_at_least_one) {
    int steps = 100;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, 1.0, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps; j <= steps; ++j) {
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (a[i] == j) {
                found = true;
            }
        }
        // could fail in 1 : 10E-88 Cases
        EXPECT_TRUE(found) << "didn't found: " << j;
    }
}

TEST_F(UniformRandomFieldTest, steps_100_at_least_one2) {
    int steps = 100;
    int step_size = 2;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps / step_size; j <= steps / step_size; ++j) {
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (a[i] == j) {
                found = true;
            }
        }
        // could fail in 1 : 10E-45 Cases
        EXPECT_TRUE(found) << "didn't found: " << j;
    }
}

TEST_F(UniformRandomFieldTest, steps_100_at_least_one_h) {
    int steps = 100;
    real step_size = 0.5;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, step_size, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps / step_size; j <= steps / step_size; ++j) {
        bool found = false;
        for (size_t i = 0; i < size; ++i) {
            if (a[i] == j) {
                found = true;
            }
        }
        // could fail in 1 : 10E-45 cases
        EXPECT_TRUE(found) << "didn't found: " << j;
    }
}

TEST_F(UniformRandomFieldTest, stress_test) {
    int steps = 100;
    int step_size = 2;
    size_t size = 10000;

    for (int i = 0; i < 10000; ++i) {
        UniformRandom noise_maker(steps, step_size, 0);
        Field a = noise_maker.random_field(size);
        EXPECT_GE(a[a[0]], -steps);
        EXPECT_LE(a[a[0]], steps);
    }
}

TEST_F(UniformRandomFieldTest, stress_test2) {
    int steps = 100;
    int step_size = 2;
    size_t size = 20000;

    for (int i = 0; i < 10000; ++i) {
        UniformRandom noise_maker(steps, step_size, 0);
        Field a = noise_maker.random_field(size);
        EXPECT_GE(a[a[0]], -steps);
        EXPECT_LE(a[a[0]], steps);
    }
}

TEST_F(UniformRandomFieldTest, step_size) {
    int steps = 20;
    real step_size = 0.1;
    size_t size = 1000;

    UniformRandom noise_maker(steps, step_size, 0);
    Field a = noise_maker.random_field(size);
    for (int i = 0; i < size; ++i) {
        real tmp = a[i] * 10;
        int integer = static_cast<int>(tmp);
        ASSERT_GE(tmp - integer, 0);
    }
}
