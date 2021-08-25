
#include <gtest/gtest.h>


#include "../../src/Domain.h"
#include "../../src/utility/Utility.h"
#include "../../src/source/GaussFunction.h"

class GaussFunctionTest : public testing::Test {
};

// TEST_F(GaussFunctionTest, test_single_obst_corner_0) {
//     auto logger = Utility::create_logger("test_single_obst_corner_0", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 5, 5, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_corner_1) {
//     auto logger = Utility::create_logger("test_single_obst_corner_1", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 5, 59, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_corner_2) {
//     auto logger = Utility::create_logger("test_single_obst_corner_2", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 59, 5, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_corner_3) {
//     auto logger = Utility::create_logger("test_single_obst_corner_3", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 59, 59, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_corner_4) {
//     auto logger = Utility::create_logger("test_single_obst_corner_4", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 59, 5, 5, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_corner_5) {
//     auto logger = Utility::create_logger("test_single_obst_corner_5", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 59, 5, 59, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_corner_6) {
//     auto logger = Utility::create_logger("test_single_obst_corner_6", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 59, 59, 5, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_corner_7) {
//     auto logger = Utility::create_logger("test_single_obst_corner_7", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 59, 59, 59, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obstsurf__0) {
//     auto logger = Utility::create_logger("test_single_obst_surf_000", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 17, 17, 5, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_surf_1) {
//     auto logger = Utility::create_logger("test_single_obst_surf_1", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 17, 5, 17, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_surf_2) {
//     auto logger = Utility::create_logger("test_single_obst_surf_2", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 17, 5, 5, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_surf_3) {
//     auto logger = Utility::create_logger("test_single_obst_surf_3", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 17, 17, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_surf_4) {
//     auto logger = Utility::create_logger("test_single_obst_surf_4", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 17, 5, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_surf_5) {
//     auto logger = Utility::create_logger("test_single_obst_surf_5", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     auto ret = GaussFunction::test_obstacle_blocks(32, 32, 32, 5, 5, 17, obst);
//     EXPECT_EQ(ret, true);
// }
// 
// TEST_F(GaussFunctionTest, test_single_obst_field) {
//     auto logger = Utility::create_logger("test_single_obst", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     Obstacle obst(16, 16, 16, 48, 48, 48, 0, logger, cube);
//     Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, cube.get_size());
// 
//     // GaussFunction gauss(25000, 1.023415823, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 5);
//     // gauss.update_source(&out, 0.0);
// }

TEST_F(GaussFunctionTest, test_single_obst_field_2d_1) {
    auto logger = Utility::create_logger("test_single_obst", "debug", "tmp.out");
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
    GaussFunction gauss(0, 0, 0, 0, 0, 0, 0, 0, 0, logger);

    int X0 = 32;
    int Y0 = 32;

    for (int x=0; x < 64; ++x) {
        for (int y=0; y < 64; ++y) {
            auto ret = gauss.test_obstacle_blocks(X0, Y0, 0, x, y, 0, obst);
            ASSERT_EQ(ret,
                    x < X0 && y < Y0) << "Failed at (" << x << ", " << y << ")";
        }
    }
}

// TEST_F(GaussFunctionTest, test_single_obst_field_2d_2) {
//     auto logger = Utility::create_logger("test_single_obst", "debug", "tmp.out");
//     Domain cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
//                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 64, 64, 64, 1);
// 
//     // |ooooooo|
//     // |ooooooo|
//     // |ooooooo|
//     // |oooXooo|
//     // |.BBoooo|
//     // |.BBoooo|
//     // |...oooo|
// 
//     Obstacle obst(33, 33, 0, 34, 34, 64, 0, logger, cube);
//     Field out(FieldType::UNKNOWN_FIELD, 0.0, 0, cube.get_size());
//     GaussFunction gauss(0, 0, 0, 0, 0, 0, 0, 0, 0, logger);
// 
//     int X0 = 32;
//     int Y0 = 32;
// 
//     for (int x=0; x < 64; ++x) {
//         for (int y=0; y < 64; ++y) {
//             auto ret = gauss.test_obstacle_blocks(X0, Y0, 1, x, y, 1, obst);
//             ASSERT_EQ(ret,
//                     x > X0 && y > Y0) << "Failed at (" << x << ", " << y << ")";
//         }
//     }
// }
