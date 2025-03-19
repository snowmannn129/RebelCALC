#include <gtest/gtest.h>
#include "../backend/numeric_solver.h"

#include <cmath>
#include <functional>

namespace rebelcalc {
namespace testing {

class NumericSolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        numericSolver = std::make_shared<NumericSolver>();
        ASSERT_TRUE(numericSolver->initialize());
    }
    
    void TearDown() override {
        numericSolver->shutdown();
        numericSolver.reset();
    }
    
    std::shared_ptr<NumericSolver> numericSolver;
};

TEST_F(NumericSolverTest, SolveLinearSystem) {
    // Solve the system:
    // 2x + 3y = 8
    // 4x + 9y = 22
    std::vector<std::vector<double>> coefficients = {
        {2.0, 3.0},
        {4.0, 9.0}
    };
    std::vector<double> constants = {8.0, 22.0};
    
    auto result = numericSolver->solveLinearSystem(coefficients, constants);
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(2, result->size());
    
    // Expected solution: x = 1, y = 2
    EXPECT_NEAR(1.0, (*result)[0], 1e-10);
    EXPECT_NEAR(2.0, (*result)[1], 1e-10);
}

TEST_F(NumericSolverTest, FindRoots) {
    // Find roots of x^2 - 4 = 0
    std::vector<double> coefficients = {1.0, 0.0, -4.0};
    
    auto result = numericSolver->findRoots(coefficients);
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(2, result->size());
    
    // Expected roots: x = 2, x = -2
    // Sort the roots to ensure consistent test results
    std::sort(result->begin(), result->end());
    EXPECT_NEAR(-2.0, (*result)[0], 1e-10);
    EXPECT_NEAR(2.0, (*result)[1], 1e-10);
}

TEST_F(NumericSolverTest, FindMinimum) {
    // Find minimum of f(x) = x^2 + 2x + 1 in [-10, 10]
    auto function = [](double x) { return x * x + 2 * x + 1; };
    
    auto result = numericSolver->findMinimum(function, 0.0, -10.0, 10.0);
    ASSERT_TRUE(result.has_value());
    
    // Expected minimum at x = -1
    EXPECT_NEAR(-1.0, *result, 1e-1);
}

TEST_F(NumericSolverTest, FindMaximum) {
    // Find maximum of f(x) = -(x^2) + 2x + 1 in [-10, 10]
    auto function = [](double x) { return -(x * x) + 2 * x + 1; };
    
    auto result = numericSolver->findMaximum(function, 0.0, -10.0, 10.0);
    ASSERT_TRUE(result.has_value());
    
    // Expected maximum at x = 1
    EXPECT_NEAR(1.0, *result, 1e-1);
}

TEST_F(NumericSolverTest, Integrate) {
    // Integrate f(x) = x^2 from 0 to 1
    auto function = [](double x) { return x * x; };
    
    auto result = numericSolver->integrate(function, 0.0, 1.0);
    ASSERT_TRUE(result.has_value());
    
    // Expected result: 1/3
    EXPECT_NEAR(1.0/3.0, *result, 1e-3);
}

TEST_F(NumericSolverTest, Differentiate) {
    // Differentiate f(x) = x^2 at x = 2
    auto function = [](double x) { return x * x; };
    
    auto result = numericSolver->differentiate(function, 2.0);
    ASSERT_TRUE(result.has_value());
    
    // Expected result: 4
    EXPECT_NEAR(4.0, *result, 1e-3);
}

} // namespace testing
} // namespace rebelcalc
