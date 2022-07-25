#include <gtest/gtest.h>


#include "../../src/domain/DomainData.h"
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

TEST_F(ObstacleTest, testSimple) {
    Coordinate<size_t> obst_start(25, 16, 25);
    Coordinate<size_t> obst_end(39, 18, 39);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 1, 32);

    for (size_t j = 1; j < 16; j++) {
        Coordinate<size_t> end(32, j, 32);
        bool ret = obst.line_crosses(start, end);
        EXPECT_FALSE(ret) << "Failed at y = " << j;
    }
    for (size_t j = 16; j <= 64; j++) {
        Coordinate<size_t> end(32, j, 32);
        bool ret = obst.line_crosses(start, end);
        EXPECT_TRUE(ret) << "Failed at y = " << j;
    }
}

TEST_F(ObstacleTest, testRoom) {
    std::vector<Obstacle> room;
    room.reserve(5);

    Coordinate<size_t> obst_wall_ceiling_start(22, 41, 22);
    Coordinate<size_t> obst_wall_ceiling_end(42, 42, 42);
    room.emplace_back(obst_wall_ceiling_start, obst_wall_ceiling_end, 0, "ceiling");

    Coordinate<size_t> obst_wall_left_start(22, 0, 22);
    Coordinate<size_t> obst_wall_left_end(25, 40, 42);
    room.emplace_back(obst_wall_left_start, obst_wall_left_end, 0, "left");

    Coordinate<size_t> obst_wall_right_start(39, 0, 22);
    Coordinate<size_t> obst_wall_right_end(42, 40, 42);
    room.emplace_back(obst_wall_right_start, obst_wall_right_end, 0, "right");

    Coordinate<size_t> obst_wall_front_start(26, 0, 22);
    Coordinate<size_t> obst_wall_front_end(38, 40, 25);
    room.emplace_back(obst_wall_front_start, obst_wall_front_end, 0, "front");

    Coordinate<size_t> obst_wall_back_start(26, 0, 39);
    Coordinate<size_t> obst_wall_back_end(38, 40, 42);
    room.emplace_back(obst_wall_back_start, obst_wall_back_end, 0, "back");

    // heat source
    Coordinate<size_t> start(32, 1, 32);
    auto no_inner_cells = DomainData::getInstance()->get_number_of_inner_cells();
    for (size_t i = 1; i <= no_inner_cells[CoordinateAxis::X]; i++) {
        for (size_t j = 1; j <= no_inner_cells[CoordinateAxis::Y]; j++) {
            for (size_t k = 1; k <= no_inner_cells[CoordinateAxis::Z]; k++) {
                Coordinate<size_t> end(i, j, k);
                bool ret = false;
                for (Obstacle &obst: room) {
                    bool tmp = obst.line_crosses(start, end);
                    ret = ret || tmp;
                    if (i > 25 && i < 39 && j < 41 && k > 25 && k < 39) {
                        EXPECT_FALSE(tmp) << "Failed for obst "<< obst.get_name() << " at (" << i << "|" << j << "|" << k << ")";
                    }
                }
                if (i > 25 && i < 39 && j < 41 && k > 25 && k < 39) {
                    EXPECT_FALSE(ret) << "Failed for (" << i << "|" << j << "|" << k << ")";
                } else {
                    EXPECT_TRUE(ret) << "Failed for (" << i << "|" << j << "|" << k << ")";
                }
            }
        }
    }
}

TEST_F(ObstacleTest, testSingleObstCorner0) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 5, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstCorner1) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 5, 59);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstCorner2) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 59, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstCorner3) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 59, 59);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstCorner4) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(59, 5, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstCorner5) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(59, 5, 59);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstCorner6) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(59, 59, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstCorner7) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(59, 59, 59);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstSurf0) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(17, 17, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstSurf1) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(17, 5, 17);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstSurf2) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(17, 5, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstSurf3) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 17, 17);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstSurf4) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 17, 5);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstSurf5) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Coordinate<size_t> start(32, 32, 32);
    Coordinate<size_t> end(5, 5, 17);
    auto ret = obst.line_crosses(start, end);
    EXPECT_EQ(ret, true);
}

TEST_F(ObstacleTest, testSingleObstField) {
    Coordinate<size_t> obst_start(16, 16, 16);
    Coordinate<size_t> obst_end(48, 48, 48);
    Obstacle obst(obst_start, obst_end, 0, "cube");

    Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, DomainData::getInstance()->get_size());

// Obstacle gauss(25000, 1.023415823, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 5);
// Obstacle::update_source(&out, 0.0);
}

TEST_F(ObstacleTest, testSingleObstField2d1) {
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

TEST_F(ObstacleTest, testSingleObstField2d2) {
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

TEST_F(ObstacleTest, stressTestSingleObstField1) {
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
