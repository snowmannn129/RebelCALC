#include <gtest/gtest.h>
#include "../backend/symbolic_engine.h"

namespace rebelcalc {
namespace testing {

class SymbolicEngineTest : public ::testing::Test {
protected:
    void SetUp() override {
        symbolicEngine = std::make_shared<SymbolicEngine>();
        ASSERT_TRUE(symbolicEngine->initialize());
    }
    
    void TearDown() override {
        symbolicEngine->shutdown();
        symbolicEngine.reset();
    }
    
    std::shared_ptr<SymbolicEngine> symbolicEngine;
};

TEST_F(SymbolicEngineTest, SimplifyExpression) {
    auto result = symbolicEngine->simplify("x + 0");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x", *result);
    
    result = symbolicEngine->simplify("x * 1");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x", *result);
    
    result = symbolicEngine->simplify("x * 0");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("0", *result);
}

TEST_F(SymbolicEngineTest, ExpandExpression) {
    auto result = symbolicEngine->expand("(x + y)^2");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x^2 + 2*x*y + y^2", *result);
    
    result = symbolicEngine->expand("(x + y)*(x - y)");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x^2 - y^2", *result);
}

TEST_F(SymbolicEngineTest, FactorExpression) {
    auto result = symbolicEngine->factor("x^2 - y^2");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("(x + y)*(x - y)", *result);
    
    result = symbolicEngine->factor("x^2 + 2*x*y + y^2");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("(x + y)^2", *result);
}

TEST_F(SymbolicEngineTest, SolveEquation) {
    auto result = symbolicEngine->solve("x + 5 = 10", "x");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("5", *result);
    
    result = symbolicEngine->solve("2*x = 10", "x");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("5", *result);
}

TEST_F(SymbolicEngineTest, DifferentiateExpression) {
    auto result = symbolicEngine->differentiate("x^2", "x");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("2*x", *result);
    
    result = symbolicEngine->differentiate("sin(x)", "x");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("cos(x)", *result);
}

TEST_F(SymbolicEngineTest, IntegrateExpression) {
    auto result = symbolicEngine->integrate("x", "x");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x^2/2", *result);
    
    result = symbolicEngine->integrate("x^2", "x");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x^3/3", *result);
}

TEST_F(SymbolicEngineTest, SubstituteExpression) {
    auto result = symbolicEngine->substitute("x + y", "x", "a + b");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("a + b + y", *result);
}

} // namespace testing
} // namespace rebelcalc
