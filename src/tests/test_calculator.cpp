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

TEST_F(CalculatorTest, GetBuiltInConstants) {
    auto pi = calculator->getVariable("pi");
    ASSERT_TRUE(pi.has_value());
    EXPECT_DOUBLE_EQ(M_PI, *pi);
    
    auto e = calculator->getVariable("e");
    ASSERT_TRUE(e.has_value());
    EXPECT_DOUBLE_EQ(M_E, *e);
}

TEST_F(CalculatorTest, SetAndGetVariable) {
    calculator->setVariable("x", 42.0);
    
    auto x = calculator->getVariable("x");
    ASSERT_TRUE(x.has_value());
    EXPECT_DOUBLE_EQ(42.0, *x);
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

} // namespace testing
} // namespace rebelcalc
