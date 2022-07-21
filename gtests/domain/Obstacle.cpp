#include <gtest/gtest.h>


#include "../../src/domain/DomainData.h"
#include "../../src/utility/Utility.h"
#include "../../src/domain/Obstacle.h"

class ObstacleTest : public testing::Test {
    void SetUp() override {
        DomainData::reset();
        Settings::domain_parameters domain_parameters{ };
        domain_parameters.enable_computational_domain = false;
        domain_parameters.number_of_inner_cells.set_coordinate(64, 64, 64);
        domain_parameters.start_coords_PD.set_coordinate(0, 0, 0);
        domain_parameters.end_coords_PD.set_coordinate(1, 1, 1);
        domain_parameters.start_coords_CD.copy(domain_parameters.start_coords_PD);
        domain_parameters.end_coords_CD.copy(domain_parameters.end_coords_PD);
        // config
        DomainData::init({ }, domain_parameters, 1);
    }

    void TearDown() override {
        DomainData::reset();
    }
};

TEST_F(ObstacleTest, test_single_obst_corner_0) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 5, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_corner_1) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 5, 59);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_corner_2) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 59, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_corner_3) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 59, 59);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_corner_4) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(59, 5, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_corner_5) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(59, 5, 59);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_corner_6) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(59, 59, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_corner_7) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(59, 59, 59);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obstsurf__0) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(17, 17, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_surf_1) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(17, 5, 17);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_surf_2) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(17, 5, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_surf_3) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 17, 17);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_surf_4) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 17, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_surf_5) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 5, 17);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, test_single_obst_field) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, DomainData::getInstance()->get_size());

// Obstacle gauss(25000, 1.023415823, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 5);
// Obstacle::update_source(&out, 0.0);
}

TEST_F(ObstacleTest, test_single_obst_field_2d_1) {
// |ooooooo|
// |ooooooo|
// |ooooooo|
// |oooXooo|
// |.BBoooo|
// |.BBoooo|
// |...oooo|

    Coordinate<size_t> obst_start(30, 30, 0);
    Coordinate<size_t> obst_end(31, 31, 64);
    Obstacle obst(obst_start, obst_end, 0, "cube");
    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, DomainData::getInstance()->get_size());

    int X0 = 32;
    int Y0 = 32;
    int x, y;

// for (x=0; x < 64; ++x) {
//     for (y=0; y < 64; ++y) {
//         std::cerr << "XXX" << x << " " << y << " " << obst.line_crosses(X0, Y0, 0, x, y, 0) << std::endl;
//     }
// }
    for (x = 0, y = 0; x < X0 / 2; ++x) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }
    for (x = X0 / 2 + 1, y = 0; x < 64; ++x) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x = 0, y = Y0 / 2; x < X0 / 4 * 3; ++x) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }
    for (x = X0 / 4 * 3 + 1, y = Y0 / 2; x < 64; ++x) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x = 0, y = 0; y < Y0 / 2; ++y) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x = X0 / 2, y = 0; y < Y0 / 4 * 3; ++y) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x = X0 / 2, y = Y0 / 4 * 3 + 1; y < 64; ++y) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x = 0; x < 64; ++x) {
        for (y = 0; y < 64; ++y) {
            Coordinate<size_t> start(X0, Y0, 0);
            Coordinate<size_t> end(x, y, 0);
            auto ret = obst.line_crosses(start, end);
            if (x > X0 || y > Y0) {
                EXPECT_FALSE(ret);
            }
        }
    }
}

TEST_F(ObstacleTest, test_single_obst_field_2d_2) {
    Coordinate<size_t> obst_start(33, 33, 0);
    Coordinate<size_t> obst_end(34, 34, 64);
    Obstacle obst(obst_start, obst_end, 0, "cube");
    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, DomainData::getInstance()->get_size());

    int X0 = 32;
    int Y0 = 32;
    int x, y;

    for (x = 0, y = 64; x < X0 / 2 * 3; ++x) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }
    for (x = X0 / 2 * 3 + 1, y = 64; x < 64; ++x) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x = 64, y = 0; y < Y0 / 2 * 3; ++y) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }
    for (x = 64, y = X0 / 2 * 3 + 1; y < 64; ++y) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_TRUE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x = X0 / 2, y = Y0 / 4 * 3 + 1; y < 64; ++y) {
        Coordinate<size_t> start(X0, Y0, 0);
        Coordinate<size_t> end(x, y, 0);
        auto ret = obst.line_crosses(start, end);
        EXPECT_FALSE(ret) << "Failed at (" << x << ", " << y << ")";
    }

    for (x = 0; x < 64; ++x) {
        for (y = 0; y < 64; ++y) {
            Coordinate<size_t> start(X0, Y0, 0);
            Coordinate<size_t> end(x, y, 0);
            auto ret = obst.line_crosses(start, end);
            if (x < X0 || y < Y0) {
                EXPECT_FALSE(ret);
            }
        }
    }
}

TEST_F(ObstacleTest, stress_test_single_obst_field1) {
    DomainData::reset();
    Settings::domain_parameters domain_parameters{ };
    domain_parameters.enable_computational_domain = false;
    domain_parameters.number_of_inner_cells.set_coordinate(64 * 8, 64 * 8, 64 * 8);
    domain_parameters.start_coords_PD.set_coordinate(0, 0, 0);
    domain_parameters.end_coords_PD.set_coordinate(1, 1, 1);
    domain_parameters.start_coords_CD.copy(domain_parameters.start_coords_PD);
    domain_parameters.end_coords_CD.copy(domain_parameters.end_coords_PD);
    // config
    DomainData::init({ }, domain_parameters, 1);

    Coordinate<size_t> obst_start(30 * 8, 30 * 8, 0 * 8);
    Coordinate<size_t> obst_end(31 * 8, 31 * 8, 64 * 8);
    Obstacle obst(obst_start, obst_end, 0, "cube");
    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, DomainData::getInstance()->get_size());

    int X0 = 32 * 8;
    int Y0 = 32 * 8;
    int x, y, z, count;

    for (x = 0; x < 64 * 8; ++x) {
        for (y = 0; y < 64 * 8; ++y) {
            for (z = 0; z < 64 * 8; ++z) {
                Coordinate<size_t> start(X0, Y0, 32 * 8);
                Coordinate<size_t> end(x, y, z);
                count += obst.line_crosses(start, end);
            }
        }
    }
}
