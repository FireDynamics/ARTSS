
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

    // TODO calculate probality of failure
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

    // TODO calculate probality of failure
    EXPECT_GT(count, 1000);
}

TEST_F(UniformRandomFieldTest, sucdiff) {
    int steps = 10;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, 1.0, 1337);
    Field a = noise_maker.random_field(size);
    Field b = noise_maker.random_field(size);

    int count = 0;
    for (auto i = 0; i < size; ++i) {
        count += a[i] == b[i];
    }

    // TODO calculate probality of failure
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

TEST_F(UniformRandomFieldTest, steps_10_atleatone) {
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

TEST_F(UniformRandomFieldTest, steps_10_atleatone2) {
    int steps = 10;
    int ssize = 2;
    size_t size = 2 * steps * 1000;

    UniformRandom noise_maker(steps, ssize, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps / ssize; j <= steps / ssize; ++j) {
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

TEST_F(UniformRandomFieldTest, steps_10_atleatoneh) {
    int steps = 10;
    real ssize = 0.5;
    size_t size = 2 * steps * 1000;

    UniformRandom noise_maker(steps, ssize, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps / ssize; j <= steps / ssize; ++j) {
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

TEST_F(UniformRandomFieldTest, steps_100_atleatone) {
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

TEST_F(UniformRandomFieldTest, steps_100_atleatone2) {
    int steps = 100;
    int ssize = 2;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, ssize, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps / ssize; j <= steps / ssize; ++j) {
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

TEST_F(UniformRandomFieldTest, steps_100_atleatoneh) {
    int steps = 100;
    real ssize = 0.5;
    size_t size = 2 * steps * 10000;

    UniformRandom noise_maker(steps, ssize, 0);
    Field a = noise_maker.random_field(size);

    for (int j = -steps / ssize; j <= steps / ssize; ++j) {
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

TEST_F(UniformRandomFieldTest, stess_test) {
    int steps = 100;
    int ssize = 2;
    size_t size = 10000;

    for (int i = 0; i < 10000; ++i) {
        UniformRandom noise_maker(steps, ssize, 0);
        Field a = noise_maker.random_field(size);
        EXPECT_GE(a[a[0]], -steps);
        EXPECT_LE(a[a[0]], steps);
    }
}

TEST_F(UniformRandomFieldTest, stess_test2) {
    int steps = 100;
    int ssize = 2;
    size_t size = 20000;

    for (int i = 0; i < 10000; ++i) {
        UniformRandom noise_maker(steps, ssize, 0);
        Field a = noise_maker.random_field(size);
        EXPECT_GE(a[a[0]], -steps);
        EXPECT_LE(a[a[0]], steps);
    }
}
