
#include <gtest/gtest.h>


#include "../../src/Domain.h"
#include "../../src/utility/Utility.h"
#include "../../src/source/GaussFunction.h"

class GaussFunctionTest : public testing::Test {
};

TEST_F(GaussFunctionTest, test_single_obst_corner_0) {
    auto logger = Utility::create_logger("test_single_obst_corner_0", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 5, 5, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_corner_1) {
    auto logger = Utility::create_logger("test_single_obst_corner_1", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 5, 59, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_corner_2) {
    auto logger = Utility::create_logger("test_single_obst_corner_2", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 59, 5, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_corner_3) {
    auto logger = Utility::create_logger("test_single_obst_corner_3", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 59, 59, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_corner_4) {
    auto logger = Utility::create_logger("test_single_obst_corner_4", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 59, 5, 5, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_corner_5) {
    auto logger = Utility::create_logger("test_single_obst_corner_5", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 59, 5, 59, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_corner_6) {
    auto logger = Utility::create_logger("test_single_obst_corner_6", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 59, 59, 5, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_corner_7) {
    auto logger = Utility::create_logger("test_single_obst_corner_7", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 59, 59, 59, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obstsurf__0) {
    auto logger = Utility::create_logger("test_single_obst_surf_000", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 17, 17, 5, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_surf_1) {
    auto logger = Utility::create_logger("test_single_obst_surf_1", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 17, 5, 17, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_surf_2) {
    auto logger = Utility::create_logger("test_single_obst_surf_2", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 17, 5, 5, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_surf_3) {
    auto logger = Utility::create_logger("test_single_obst_surf_3", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 17, 17, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_surf_4) {
    auto logger = Utility::create_logger("test_single_obst_surf_4", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 17, 5, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_surf_5) {
    auto logger = Utility::create_logger("test_single_obst_surf_5", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 5, 17, obst);
    EXPECT_EQ(ret, true);
}

TEST_F(GaussFunctionTest, test_single_obst_field) {
    auto logger = Utility::create_logger("test_single_obst", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, cube.get_size());

    // GaussFunction gauss(25000, 1.023415823, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 5);
    // GaussFunction::update_source(&out, 0.0);
}

TEST_F(GaussFunctionTest, test_single_obst_field_2d_1) {
    auto logger = Utility::create_logger("test_single_obst", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    // |ooooooo|
    // |ooooooo|
    // |ooooooo|
    // |oooXooo|
    // |.BBoooo|
    // |.BBoooo|
    // |...oooo|

    Obstacle obst(30, 30, 0, 31, 31, 64, 0, logger, cube);
    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, cube.get_size());

    int X0 = 32;
    int Y0 = 32;
    int x, y;

    // for (x=0; x < 64; ++x) {
    //     for (y=0; y < 64; ++y) {
    //         std::cerr << "XXX" << x << " " << y << " " << GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst) << std::endl;
    //     }
    // }
    for (x=0, y=0; x < X0 / 2; ++x) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }
    for (x=X0/2+1, y=0; x < 64; ++x) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x=0, y=Y0/2; x < X0/4*3; ++x) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }
    for (x=X0/4*3+1, y=Y0/2; x < 64; ++x) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x=0, y=0; y < Y0 / 2; ++y) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x=X0/2, y=0; y < Y0/4*3; ++y) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x=X0/2, y=Y0/4*3+1; y < 64; ++y) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x=0; x < 64; ++x) {
        for(y=0; y < 64; ++y) {
            auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
            if (x > X0 || y > Y0) {
                EXPECT_FALSE(ret);
            }
        }
    }
}

TEST_F(GaussFunctionTest, test_single_obst_field_2d_2) {
    auto logger = Utility::create_logger("test_single_obst", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);

    Obstacle obst(33, 33, 0, 34, 34, 64, 0, logger, cube);
    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, cube.get_size());

    int X0 = 32;
    int Y0 = 32;
    int x, y;

    for (x=0, y=64; x < X0/2*3; ++x) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }
    for (x=X0/2*3+1, y=64; x < 64; ++x) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x=64, y=0; y < Y0/2*3; ++y) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }
    for (x=64, y=X0/2*3+1; y < 64; ++y) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x=X0/2, y=Y0/4*3+1; y < 64; ++y) {
        auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x=0; x < 64; ++x) {
        for(y=0; y < 64; ++y) {
            auto ret = GaussFunction::test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
            if (x < X0 || y < Y0) {
                EXPECT_FALSE(ret);
            }
        }
    }
}

TEST_F(GaussFunctionTest, stress_test_single_obst_field1) {
    auto logger = Utility::create_logger("test_single_obst", "warn", "tmp.out");
    Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64*8, 64*8, 64*8, 1);

    Obstacle obst(30*8, 30*8, 0*8, 31*8, 31*8, 64*8, 0*8, logger, cube);
    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, cube.get_size());

    int X0 = 32*8;
    int Y0 = 32*8;
    int x, y, z, count;

    for (x=0; x < 64*8; ++x) {
        for (y=0; y < 64*8; ++y) {
            for (z=0; z < 64*8; ++z) {
                count += GaussFunction::test_obstacle_blocks(X0, Y0, 32*8, x, y, z, obst);
            }
        }
    }

    EXPECT_EQ(count, 1101004554);
}
