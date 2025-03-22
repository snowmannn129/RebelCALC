#include <gtest/gtest.h>
#include "../backend/calculator.h"

namespace rebelcalc {
namespace testing {

class CalculatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        calculator = std::make_shared<Calculator>();
        ASSERT_TRUE(calculator->initialize());
    }
    
    void TearDown() override {
        calculator->shutdown();
        calculator.reset();
    }
    
    std::shared_ptr<Calculator> calculator;
};

TEST_F(CalculatorTest, EvaluateSimpleExpression) {
    auto result = calculator->evaluate("1+1");
    ASSERT_TRUE(result.has_value());
    
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(2.0, std::get<double>(*result));
}

TEST_F(CalculatorTest, EvaluateBasicArithmetic) {
    // Addition
    auto result = calculator->evaluate("3 + 4");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(7.0, std::get<double>(*result));
    
    // Subtraction
    result = calculator->evaluate("10 - 3");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(7.0, std::get<double>(*result));
    
    // Multiplication
    result = calculator->evaluate("3 * 4");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(12.0, std::get<double>(*result));
    
    // Division
    result = calculator->evaluate("12 / 4");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(3.0, std::get<double>(*result));
    
    // Exponentiation
    result = calculator->evaluate("2 ^ 3");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(8.0, std::get<double>(*result));
}

TEST_F(CalculatorTest, EvaluateWithParentheses) {
    // Simple grouping
    auto result = calculator->evaluate("(2 + 3) * 4");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(20.0, std::get<double>(*result));
    
    // Nested parentheses
    result = calculator->evaluate("2 * (3 + (4 - 2))");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(10.0, std::get<double>(*result));
}

TEST_F(CalculatorTest, EvaluateWithFunctions) {
    // Sin function
    auto result = calculator->evaluate("sin(0)");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_NEAR(0.0, std::get<double>(*result), 1e-10);
    
    // Cos function
    result = calculator->evaluate("cos(0)");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_NEAR(1.0, std::get<double>(*result), 1e-10);
    
    // Sqrt function
    result = calculator->evaluate("sqrt(16)");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(4.0, std::get<double>(*result));
    
    // Function with expression as argument
    result = calculator->evaluate("sqrt(3^2 + 4^2)");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(5.0, std::get<double>(*result));
    
    // Multiple argument function
    result = calculator->evaluate("pow(2, 3)");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(8.0, std::get<double>(*result));
}

TEST_F(CalculatorTest, EvaluateWithVariables) {
    // Set variables
    calculator->setVariable("x", 5.0);
    calculator->setVariable("y", 3.0);
    
    // Use variables in expressions
    auto result = calculator->evaluate("x + y");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(8.0, std::get<double>(*result));
    
    // Use built-in constants
    result = calculator->evaluate("pi * 2");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(M_PI * 2, std::get<double>(*result));
    
    // Mix variables and constants
    result = calculator->evaluate("x * pi");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(5.0 * M_PI, std::get<double>(*result));
}

TEST_F(CalculatorTest, EvaluateWithUnaryOperators) {
    // Unary minus
    auto result = calculator->evaluate("-5");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(-5.0, std::get<double>(*result));
    
    // Unary plus
    result = calculator->evaluate("+5");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(5.0, std::get<double>(*result));
    
    // Unary minus with expression
    result = calculator->evaluate("-(3 + 2)");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(-5.0, std::get<double>(*result));
}

TEST_F(CalculatorTest, EvaluateComplexExpressions) {
    // Set variables
    calculator->setVariable("x", 2.0);
    calculator->setVariable("y", 3.0);
    
    // Complex expression with multiple features
    auto result = calculator->evaluate("sin(pi/4) + cos(pi/4) + sqrt(x^2 + y^2)");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_NEAR(sin(M_PI/4) + cos(M_PI/4) + sqrt(2.0*2.0 + 3.0*3.0), 
                std::get<double>(*result), 1e-10);
    
    // Expression with multiple operations and precedence
    result = calculator->evaluate("2 + 3 * 4 - 5 / 2.5");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(2.0 + 3.0 * 4.0 - 5.0 / 2.5, std::get<double>(*result));
}

TEST_F(CalculatorTest, ErrorHandling) {
    // Division by zero
    auto result = calculator->evaluate("1 / 0");
    EXPECT_FALSE(result.has_value());
    
    // Unknown variable
    result = calculator->evaluate("unknown_var + 5");
    EXPECT_FALSE(result.has_value());
    
    // Unknown function
    result = calculator->evaluate("unknown_func(5)");
    EXPECT_FALSE(result.has_value());
    
    // Mismatched parentheses
    result = calculator->evaluate("(2 + 3");
    EXPECT_FALSE(result.has_value());
    
    // Invalid syntax
    result = calculator->evaluate("2 +* 3");
    EXPECT_FALSE(result.has_value());
}

TEST_F(CalculatorTest, GetBuiltInConstants) {
    auto pi = calculator->getVariable("pi");
    ASSERT_TRUE(pi.has_value());
    EXPECT_DOUBLE_EQ(M_PI, *pi);
    
    auto e = calculator->getVariable("e");
    ASSERT_TRUE(e.has_value());
    EXPECT_DOUBLE_EQ(M_E, *e);
}

TEST_F(CalculatorTest, SetAndGetVariable) {
    EXPECT_TRUE(calculator->setVariable("x", 42.0));
    
    auto x = calculator->getVariable("x");
    ASSERT_TRUE(x.has_value());
    EXPECT_DOUBLE_EQ(42.0, *x);
    
    // Invalid variable names
    EXPECT_FALSE(calculator->setVariable("", 10.0));
    EXPECT_FALSE(calculator->setVariable("123", 10.0));
    EXPECT_FALSE(calculator->setVariable("x@y", 10.0));
}

TEST_F(CalculatorTest, HasVariable) {
    EXPECT_FALSE(calculator->hasVariable("x"));
    
    calculator->setVariable("x", 42.0);
    EXPECT_TRUE(calculator->hasVariable("x"));
    
    // Built-in constants
    EXPECT_TRUE(calculator->hasVariable("pi"));
    EXPECT_TRUE(calculator->hasVariable("e"));
}

TEST_F(CalculatorTest, GetAllVariables) {
    calculator->setVariable("x", 1.0);
    calculator->setVariable("y", 2.0);
    
    auto variables = calculator->getAllVariables();
    EXPECT_EQ(4, variables.size()); // x, y, pi, e
    
    EXPECT_DOUBLE_EQ(1.0, variables["x"]);
    EXPECT_DOUBLE_EQ(2.0, variables["y"]);
    EXPECT_DOUBLE_EQ(M_PI, variables["pi"]);
    EXPECT_DOUBLE_EQ(M_E, variables["e"]);
}

TEST_F(CalculatorTest, ClearVariables) {
    calculator->setVariable("x", 42.0);
    calculator->clearVariables();
    
    auto x = calculator->getVariable("x");
    EXPECT_FALSE(x.has_value());
    
    // Built-in constants should still be available
    auto pi = calculator->getVariable("pi");
    ASSERT_TRUE(pi.has_value());
}

TEST_F(CalculatorTest, VariableAssignment) {
    // Simple assignment
    auto result = calculator->evaluate("x = 5");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(5.0, std::get<double>(*result));
    
    auto x = calculator->getVariable("x");
    ASSERT_TRUE(x.has_value());
    EXPECT_DOUBLE_EQ(5.0, *x);
    
    // Assignment with expression
    result = calculator->evaluate("y = 2 * x + 3");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(13.0, std::get<double>(*result));
    
    auto y = calculator->getVariable("y");
    ASSERT_TRUE(y.has_value());
    EXPECT_DOUBLE_EQ(13.0, *y);
}

TEST_F(CalculatorTest, CompoundAssignment) {
    // Initialize variables
    calculator->setVariable("x", 5.0);
    
    // Addition assignment
    auto result = calculator->evaluate("x += 3");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(8.0, std::get<double>(*result));
    
    auto x = calculator->getVariable("x");
    ASSERT_TRUE(x.has_value());
    EXPECT_DOUBLE_EQ(8.0, *x);
    
    // Subtraction assignment
    result = calculator->evaluate("x -= 2");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(6.0, std::get<double>(*result));
    
    x = calculator->getVariable("x");
    ASSERT_TRUE(x.has_value());
    EXPECT_DOUBLE_EQ(6.0, *x);
    
    // Multiplication assignment
    result = calculator->evaluate("x *= 2");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(12.0, std::get<double>(*result));
    
    x = calculator->getVariable("x");
    ASSERT_TRUE(x.has_value());
    EXPECT_DOUBLE_EQ(12.0, *x);
    
    // Division assignment
    result = calculator->evaluate("x /= 3");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(4.0, std::get<double>(*result));
    
    x = calculator->getVariable("x");
    ASSERT_TRUE(x.has_value());
    EXPECT_DOUBLE_EQ(4.0, *x);
    
    // Compound assignment with expression
    result = calculator->evaluate("x += 2 * 3");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<double>(*result));
    EXPECT_DOUBLE_EQ(10.0, std::get<double>(*result));
    
    x = calculator->getVariable("x");
    ASSERT_TRUE(x.has_value());
    EXPECT_DOUBLE_EQ(10.0, *x);
}

} // namespace testing
} // namespace rebelcalc
