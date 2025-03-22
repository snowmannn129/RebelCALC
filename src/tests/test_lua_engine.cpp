#include <gtest/gtest.h>
#include "../scripts/lua_engine.h"
#include "../backend/calculator.h"

namespace rebelcalc {
namespace testing {

class LuaEngineTest : public ::testing::Test {
protected:
    void SetUp() override {
        calculator = std::make_shared<Calculator>();
        ASSERT_TRUE(calculator->initialize());
        
        luaEngine = std::make_shared<LuaEngine>(calculator);
        ASSERT_TRUE(luaEngine->initialize());
    }
    
    void TearDown() override {
        luaEngine->shutdown();
        luaEngine.reset();
        
        calculator->shutdown();
        calculator.reset();
    }
    
    std::shared_ptr<Calculator> calculator;
    std::shared_ptr<LuaEngine> luaEngine;
};

TEST_F(LuaEngineTest, ExecuteSimpleScript) {
    auto result = luaEngine->executeScript("return 2 + 2");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("4", *result);
}

TEST_F(LuaEngineTest, ExecuteCalculatorFunction) {
    auto result = luaEngine->executeScript("return evaluate('1+1')");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("2", *result);
    
    result = luaEngine->executeScript("return solve('x + 5 = 10', 'x')");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("5", *result);
    
    result = luaEngine->executeScript("return differentiate('x^2', 'x')");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("2*x", *result);
    
    result = luaEngine->executeScript("return integrate('x', 'x')");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x^2/2", *result);
}

TEST_F(LuaEngineTest, ExecuteMultiLineScript) {
    std::string script = R"(
        local x = 10
        local y = 20
        return x + y
    )";
    
    auto result = luaEngine->executeScript(script);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("30", *result);
}

TEST_F(LuaEngineTest, ExecuteScriptWithVariables) {
    std::string script = R"(
        local result = evaluate('x^2 + 2*x + 1')
        return result
    )";
    
    calculator->setVariable("x", 2.0);
    
    auto result = luaEngine->executeScript(script);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("9", *result);
}

TEST_F(LuaEngineTest, ExecuteAdvancedFunctions) {
    auto result = luaEngine->executeScript("return simplify('x + 0')");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x", *result);
    
    result = luaEngine->executeScript("return expand('(x + y)^2')");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("x^2 + 2*x*y + y^2", *result);
    
    result = luaEngine->executeScript("return factor('x^2 - y^2')");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("(x + y)*(x - y)", *result);
    
    result = luaEngine->executeScript("return substitute('x + y', 'x', 'a + b')");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("a + b + y", *result);
}

TEST_F(LuaEngineTest, ExecuteMathLibraryFunctions) {
    auto result = luaEngine->executeScript("return math.factorial(5)");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("120", *result);
    
    result = luaEngine->executeScript("return math.gcd(12, 18)");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("6", *result);
    
    result = luaEngine->executeScript("return math.lcm(4, 6)");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("12", *result);
    
    result = luaEngine->executeScript("return math.pi");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(std::stod(*result) > 3.14 && std::stod(*result) < 3.15);
    
    result = luaEngine->executeScript("return math.e");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(std::stod(*result) > 2.71 && std::stod(*result) < 2.72);
}

TEST_F(LuaEngineTest, ExecuteCalculatorLibraryFunctions) {
    auto result = luaEngine->executeScript(R"(
        calculator.setVariable("x", 5)
        return calculator.getVariable("x")
    )");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("5", *result);
    
    result = luaEngine->executeScript(R"(
        return calculator.hasVariable("x")
    )");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("true", *result);
    
    result = luaEngine->executeScript(R"(
        calculator.clearVariables()
        return calculator.hasVariable("x")
    )");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ("false", *result);
}

TEST_F(LuaEngineTest, GetAvailableFunctions) {
    auto functions = luaEngine->getAvailableFunctions();
    ASSERT_FALSE(functions.empty());
    
    // Check for basic functions
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "evaluate") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "solve") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "differentiate") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "integrate") != functions.end());
    
    // Check for advanced functions
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "simplify") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "expand") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "factor") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "substitute") != functions.end());
    
    // Check for numeric functions
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "solveLinearSystem") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "findRoots") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "findMinimum") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "findMaximum") != functions.end());
}

TEST_F(LuaEngineTest, GetFunctionHelp) {
    auto help = luaEngine->getFunctionHelp("evaluate");
    ASSERT_TRUE(help.has_value());
    EXPECT_FALSE(help->empty());
}

} // namespace testing
} // namespace rebelcalc
