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
    
    // Solve a 3x3 system:
    // x + y + z = 6
    // 2x - y + z = 3
    // x + 2y - z = 2
    std::vector<std::vector<double>> coefficients3 = {
        {1.0, 1.0, 1.0},
        {2.0, -1.0, 1.0},
        {1.0, 2.0, -1.0}
    };
    std::vector<double> constants3 = {6.0, 3.0, 2.0};
    
    auto result3 = numericSolver->solveLinearSystem(coefficients3, constants3);
    ASSERT_TRUE(result3.has_value());
    ASSERT_EQ(3, result3->size());
    
    // Expected solution: x = 1, y = 2, z = 3
    EXPECT_NEAR(1.0, (*result3)[0], 1e-10);
    EXPECT_NEAR(2.0, (*result3)[1], 1e-10);
    EXPECT_NEAR(3.0, (*result3)[2], 1e-10);
    
    // Test singular matrix
    std::vector<std::vector<double>> singularCoefficients = {
        {1.0, 2.0},
        {2.0, 4.0}
    };
    std::vector<double> singularConstants = {3.0, 6.0};
    
    auto singularResult = numericSolver->solveLinearSystem(singularCoefficients, singularConstants);
    EXPECT_FALSE(singularResult.has_value());
    
    // Test inconsistent system
    std::vector<std::vector<double>> inconsistentCoefficients = {
        {1.0, 2.0},
        {2.0, 4.0}
    };
    std::vector<double> inconsistentConstants = {3.0, 7.0};
    
    auto inconsistentResult = numericSolver->solveLinearSystem(inconsistentCoefficients, inconsistentConstants);
    EXPECT_FALSE(inconsistentResult.has_value());
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
    
    // Find roots of x^3 - x^2 - 4x + 4 = 0
    std::vector<double> cubicCoefficients = {1.0, -1.0, -4.0, 4.0};
    
    auto cubicResult = numericSolver->findRoots(cubicCoefficients);
    ASSERT_TRUE(cubicResult.has_value());
    ASSERT_EQ(3, cubicResult->size());
    
    // Expected roots: x = -2, x = 1, x = 2
    // Sort the roots to ensure consistent test results
    std::sort(cubicResult->begin(), cubicResult->end());
    EXPECT_NEAR(-2.0, (*cubicResult)[0], 1e-10);
    EXPECT_NEAR(1.0, (*cubicResult)[1], 1e-10);
    EXPECT_NEAR(2.0, (*cubicResult)[2], 1e-10);
    
    // Find roots of x - 3 = 0
    std::vector<double> linearCoefficients = {1.0, -3.0};
    
    auto linearResult = numericSolver->findRoots(linearCoefficients);
    ASSERT_TRUE(linearResult.has_value());
    ASSERT_EQ(1, linearResult->size());
    
    // Expected root: x = 3
    EXPECT_NEAR(3.0, (*linearResult)[0], 1e-10);
    
    // Test constant polynomial
    std::vector<double> constantCoefficients = {5.0};
    
    auto constantResult = numericSolver->findRoots(constantCoefficients);
    ASSERT_TRUE(constantResult.has_value());
    ASSERT_EQ(0, constantResult->size());
    
    // Test polynomial with complex roots
    std::vector<double> complexRootsCoefficients = {1.0, 0.0, 1.0};
    
    auto complexResult = numericSolver->findRoots(complexRootsCoefficients);
    ASSERT_TRUE(complexResult.has_value());
    ASSERT_EQ(0, complexResult->size());
}

TEST_F(NumericSolverTest, FindMinimum) {
    // Find minimum of f(x) = x^2 + 2x + 1 in [-10, 10]
    auto function = [](double x) { return x * x + 2 * x + 1; };
    
    auto result = numericSolver->findMinimum(function, 0.0, -10.0, 10.0);
    ASSERT_TRUE(result.has_value());
    
    // Expected minimum at x = -1
    EXPECT_NEAR(-1.0, *result, 1e-6);
    
    // Find minimum of f(x) = x^4 - 4x^2 + 2 in [-10, 10]
    auto function2 = [](double x) { return x * x * x * x - 4 * x * x + 2; };
    
    auto result2 = numericSolver->findMinimum(function2, 0.0, -10.0, 10.0);
    ASSERT_TRUE(result2.has_value());
    
    // Expected minimum at x = ±√2
    double expectedMinimum = std::sqrt(2.0);
    EXPECT_TRUE(std::abs(*result2 - expectedMinimum) < 1e-6 || 
                std::abs(*result2 + expectedMinimum) < 1e-6);
    
    // Find minimum of f(x) = sin(x) in [0, 2π]
    auto function3 = [](double x) { return std::sin(x); };
    
    auto result3 = numericSolver->findMinimum(function3, 0.0, 0.0, 2 * M_PI);
    ASSERT_TRUE(result3.has_value());
    
    // Expected minimum at x = 3π/2
    EXPECT_NEAR(3 * M_PI / 2, *result3, 1e-6);
}

TEST_F(NumericSolverTest, FindMaximum) {
    // Find maximum of f(x) = -(x^2) + 2x + 1 in [-10, 10]
    auto function = [](double x) { return -(x * x) + 2 * x + 1; };
    
    auto result = numericSolver->findMaximum(function, 0.0, -10.0, 10.0);
    ASSERT_TRUE(result.has_value());
    
    // Expected maximum at x = 1
    EXPECT_NEAR(1.0, *result, 1e-6);
    
    // Find maximum of f(x) = sin(x) in [0, 2π]
    auto function2 = [](double x) { return std::sin(x); };
    
    auto result2 = numericSolver->findMaximum(function2, 0.0, 0.0, 2 * M_PI);
    ASSERT_TRUE(result2.has_value());
    
    // Expected maximum at x = π/2
    EXPECT_NEAR(M_PI / 2, *result2, 1e-6);
    
    // Find maximum of f(x) = -x^4 + 4x^2 - 2 in [-10, 10]
    auto function3 = [](double x) { return -x * x * x * x + 4 * x * x - 2; };
    
    auto result3 = numericSolver->findMaximum(function3, 0.0, -10.0, 10.0);
    ASSERT_TRUE(result3.has_value());
    
    // Expected maximum at x = ±1
    EXPECT_TRUE(std::abs(*result3 - 1.0) < 1e-6 || 
                std::abs(*result3 + 1.0) < 1e-6);
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
