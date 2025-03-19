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

TEST_F(LuaEngineTest, GetAvailableFunctions) {
    auto functions = luaEngine->getAvailableFunctions();
    ASSERT_FALSE(functions.empty());
    
    // Check for some expected functions
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "evaluate") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "solve") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "differentiate") != functions.end());
    EXPECT_TRUE(std::find(functions.begin(), functions.end(), "integrate") != functions.end());
}

TEST_F(LuaEngineTest, GetFunctionHelp) {
    auto help = luaEngine->getFunctionHelp("evaluate");
    ASSERT_TRUE(help.has_value());
    EXPECT_FALSE(help->empty());
}

} // namespace testing
} // namespace rebelcalc
