#include "lua_engine.h"
#include "../backend/calculator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <functional>

// Include Lua headers
extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

namespace rebelcalc {

// Private implementation class
class LuaEngine::Impl {
public:
    Impl(std::shared_ptr<Calculator> calculator)
        : m_calculator(calculator), m_luaState(nullptr) {
    }
    
    ~Impl() {
        shutdown();
    }
    
    bool initialize() {
        // Create a new Lua state
        m_luaState = luaL_newstate();
        if (!m_luaState) {
            std::cerr << "Failed to create Lua state" << std::endl;
            return false;
        }
        
        // Load Lua libraries
        luaL_openlibs(m_luaState);
        
        // Set up the Lua environment
        if (!setupLuaEnvironment()) {
            std::cerr << "Failed to set up Lua environment" << std::endl;
            shutdown();
            return false;
        }
        
        return true;
    }
    
    void shutdown() {
        if (m_luaState) {
            lua_close(m_luaState);
            m_luaState = nullptr;
        }
    }
    
    std::optional<std::string> executeScript(const std::string& script) {
        if (!m_luaState) {
            std::cerr << "Lua state not initialized" << std::endl;
            return std::nullopt;
        }
        
        // Execute the script
        int status = luaL_dostring(m_luaState, script.c_str());
        if (status != LUA_OK) {
            std::string error = lua_tostring(m_luaState, -1);
            lua_pop(m_luaState, 1);
            std::cerr << "Lua error: " << error << std::endl;
            return std::nullopt;
        }
        
        // Check if there's a result on the stack
        if (lua_gettop(m_luaState) > 0) {
            // Convert the result to a string
            const char* result = lua_tostring(m_luaState, -1);
            lua_pop(m_luaState, 1);
            
            if (result) {
                return std::string(result);
            }
        }
        
        return "Script executed successfully";
    }
    
    bool loadScript(const std::string& filename) {
        if (!m_luaState) {
            std::cerr << "Lua state not initialized" << std::endl;
            return false;
        }
        
        // Load the script from file
        int status = luaL_loadfile(m_luaState, filename.c_str());
        if (status != LUA_OK) {
            std::string error = lua_tostring(m_luaState, -1);
            lua_pop(m_luaState, 1);
            std::cerr << "Lua error loading file: " << error << std::endl;
            return false;
        }
        
        // Execute the loaded script
        status = lua_pcall(m_luaState, 0, LUA_MULTRET, 0);
        if (status != LUA_OK) {
            std::string error = lua_tostring(m_luaState, -1);
            lua_pop(m_luaState, 1);
            std::cerr << "Lua error executing file: " << error << std::endl;
            return false;
        }
        
        return true;
    }
    
    void registerFunction(const std::string& name, 
                         std::function<std::string(const std::vector<std::string>&)> function) {
        m_functions[name] = function;
        
        // Register the function with Lua
        if (m_luaState) {
            // Create a closure that will call the C++ function
            lua_pushlightuserdata(m_luaState, this);
            lua_pushstring(m_luaState, name.c_str());
            lua_pushcclosure(m_luaState, &LuaEngine::Impl::luaFunctionDispatcher, 2);
            lua_setglobal(m_luaState, name.c_str());
        }
    }
    
    std::vector<std::string> getAvailableFunctions() const {
        std::vector<std::string> result;
        
        for (const auto& [name, _] : m_functions) {
            result.push_back(name);
        }
        
        return result;
    }
    
    std::optional<std::string> getFunctionHelp(const std::string& functionName) const {
        // This is a placeholder for a real help system
        // In a real implementation, this would look up help text from a database
        
        if (functionName == "evaluate") {
            return "evaluate(expression) - Evaluate a mathematical expression";
        } else if (functionName == "solve") {
            return "solve(equation, variable) - Solve an equation for a variable";
        } else if (functionName == "differentiate") {
            return "differentiate(expression, variable) - Differentiate an expression with respect to a variable";
        } else if (functionName == "integrate") {
            return "integrate(expression, variable) - Integrate an expression with respect to a variable";
        }
        
        return std::nullopt;
    }
    
    bool setupLuaEnvironment() {
        // Register built-in functions
        registerFunction("evaluate", [this](const std::vector<std::string>& args) {
            if (args.size() < 1) {
                return std::string("Usage: evaluate(expression)");
            }
            
            auto result = m_calculator->evaluate(args[0]);
            if (!result) {
                return std::string("Error evaluating expression");
            }
            
            std::stringstream ss;
            std::visit([&ss](auto&& arg) {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, double>) {
                    ss << arg;
                } else if constexpr (std::is_same_v<T, std::string>) {
                    ss << arg;
                } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                    for (size_t i = 0; i < arg.size(); ++i) {
                        ss << arg[i];
                        if (i < arg.size() - 1) {
                            ss << ", ";
                        }
                    }
                }
            }, *result);
            
            return ss.str();
        });
        
        registerFunction("solve", [this](const std::vector<std::string>& args) {
            if (args.size() < 2) {
                return std::string("Usage: solve(equation, variable)");
            }
            
            auto result = m_calculator->solve(args[0], args[1]);
            if (!result) {
                return std::string("Error solving equation");
            }
            
            std::stringstream ss;
            std::visit([&ss](auto&& arg) {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, double>) {
                    ss << arg;
                } else if constexpr (std::is_same_v<T, std::string>) {
                    ss << arg;
                } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                    for (size_t i = 0; i < arg.size(); ++i) {
                        ss << arg[i];
                        if (i < arg.size() - 1) {
                            ss << ", ";
                        }
                    }
                }
            }, *result);
            
            return ss.str();
        });
        
        registerFunction("differentiate", [this](const std::vector<std::string>& args) {
            if (args.size() < 2) {
                return std::string("Usage: differentiate(expression, variable)");
            }
            
            auto result = m_calculator->differentiate(args[0], args[1]);
            if (!result) {
                return std::string("Error differentiating expression");
            }
            
            return *result;
        });
        
        registerFunction("integrate", [this](const std::vector<std::string>& args) {
            if (args.size() < 2) {
                return std::string("Usage: integrate(expression, variable)");
            }
            
            auto result = m_calculator->integrate(args[0], args[1]);
            if (!result) {
                return std::string("Error integrating expression");
            }
            
            return *result;
        });
        
        // Add more advanced functions
        
        // Simplify expressions
        registerFunction("simplify", [this](const std::vector<std::string>& args) {
            if (args.size() < 1) {
                return std::string("Usage: simplify(expression)");
            }
            
            auto symbolicEngine = m_calculator->getSymbolicEngine();
            if (!symbolicEngine) {
                return std::string("Symbolic engine not available");
            }
            
            auto result = symbolicEngine->simplify(args[0]);
            if (!result) {
                return std::string("Error simplifying expression");
            }
            
            return *result;
        });
        
        // Expand expressions
        registerFunction("expand", [this](const std::vector<std::string>& args) {
            if (args.size() < 1) {
                return std::string("Usage: expand(expression)");
            }
            
            auto symbolicEngine = m_calculator->getSymbolicEngine();
            if (!symbolicEngine) {
                return std::string("Symbolic engine not available");
            }
            
            auto result = symbolicEngine->expand(args[0]);
            if (!result) {
                return std::string("Error expanding expression");
            }
            
            return *result;
        });
        
        // Factor expressions
        registerFunction("factor", [this](const std::vector<std::string>& args) {
            if (args.size() < 1) {
                return std::string("Usage: factor(expression)");
            }
            
            auto symbolicEngine = m_calculator->getSymbolicEngine();
            if (!symbolicEngine) {
                return std::string("Symbolic engine not available");
            }
            
            auto result = symbolicEngine->factor(args[0]);
            if (!result) {
                return std::string("Error factoring expression");
            }
            
            return *result;
        });
        
        // Substitute variables
        registerFunction("substitute", [this](const std::vector<std::string>& args) {
            if (args.size() < 3) {
                return std::string("Usage: substitute(expression, variable, replacement)");
            }
            
            auto symbolicEngine = m_calculator->getSymbolicEngine();
            if (!symbolicEngine) {
                return std::string("Symbolic engine not available");
            }
            
            auto result = symbolicEngine->substitute(args[0], args[1], args[2]);
            if (!result) {
                return std::string("Error substituting in expression");
            }
            
            return *result;
        });
        
        // Solve linear systems
        registerFunction("solveLinearSystem", [this](const std::vector<std::string>& args) {
            if (args.size() < 2) {
                return std::string("Usage: solveLinearSystem(coefficients, constants)");
            }
            
            auto numericSolver = m_calculator->getNumericSolver();
            if (!numericSolver) {
                return std::string("Numeric solver not available");
            }
            
            // Parse coefficients and constants from strings
            // This is a simplified implementation that expects specific formats
            // In a real implementation, this would be more robust
            
            // Example format for coefficients: "[[1,2],[3,4]]"
            // Example format for constants: "[5,6]"
            
            // For now, just return a placeholder result
            return std::string("Linear system solving not fully implemented in Lua yet");
        });
        
        // Find polynomial roots
        registerFunction("findRoots", [this](const std::vector<std::string>& args) {
            if (args.size() < 1) {
                return std::string("Usage: findRoots(coefficients)");
            }
            
            auto numericSolver = m_calculator->getNumericSolver();
            if (!numericSolver) {
                return std::string("Numeric solver not available");
            }
            
            // Parse coefficients from string
            // This is a simplified implementation that expects specific formats
            // In a real implementation, this would be more robust
            
            // Example format for coefficients: "[1,2,3]" for x^2 + 2x + 3
            
            // For now, just return a placeholder result
            return std::string("Polynomial root finding not fully implemented in Lua yet");
        });
        
        // Find minimum of a function
        registerFunction("findMinimum", [this](const std::vector<std::string>& args) {
            if (args.size() < 4) {
                return std::string("Usage: findMinimum(function, initialGuess, lowerBound, upperBound)");
            }
            
            auto numericSolver = m_calculator->getNumericSolver();
            if (!numericSolver) {
                return std::string("Numeric solver not available");
            }
            
            // Parse function and bounds from strings
            // This is a simplified implementation that expects specific formats
            // In a real implementation, this would be more robust
            
            // For now, just return a placeholder result
            return std::string("Function minimization not fully implemented in Lua yet");
        });
        
        // Find maximum of a function
        registerFunction("findMaximum", [this](const std::vector<std::string>& args) {
            if (args.size() < 4) {
                return std::string("Usage: findMaximum(function, initialGuess, lowerBound, upperBound)");
            }
            
            auto numericSolver = m_calculator->getNumericSolver();
            if (!numericSolver) {
                return std::string("Numeric solver not available");
            }
            
            // Parse function and bounds from strings
            // This is a simplified implementation that expects specific formats
            // In a real implementation, this would be more robust
            
            // For now, just return a placeholder result
            return std::string("Function maximization not fully implemented in Lua yet");
        });
        
        // Set up Lua standard libraries and environment
        
        // Set up the math library
        lua_getglobal(m_luaState, "math");
        if (lua_isnil(m_luaState, -1)) {
            lua_pop(m_luaState, 1);
            lua_newtable(m_luaState);
            lua_setglobal(m_luaState, "math");
            lua_getglobal(m_luaState, "math");
        }
        
        // Add constants to the math library
        lua_pushnumber(m_luaState, M_PI);
        lua_setfield(m_luaState, -2, "pi");
        
        lua_pushnumber(m_luaState, M_E);
        lua_setfield(m_luaState, -2, "e");
        
        // Add additional math functions
        
        // math.factorial
        lua_pushcfunction(m_luaState, [](lua_State* L) -> int {
            int n = luaL_checkinteger(L, 1);
            if (n < 0) {
                return luaL_error(L, "factorial of negative number");
            }
            
            lua_Integer result = 1;
            for (int i = 2; i <= n; ++i) {
                result *= i;
            }
            
            lua_pushinteger(L, result);
            return 1;
        });
        lua_setfield(m_luaState, -2, "factorial");
        
        // math.gcd
        lua_pushcfunction(m_luaState, [](lua_State* L) -> int {
            int a = luaL_checkinteger(L, 1);
            int b = luaL_checkinteger(L, 2);
            
            // Euclidean algorithm
            while (b != 0) {
                int temp = b;
                b = a % b;
                a = temp;
            }
            
            lua_pushinteger(L, a);
            return 1;
        });
        lua_setfield(m_luaState, -2, "gcd");
        
        // math.lcm
        lua_pushcfunction(m_luaState, [](lua_State* L) -> int {
            int a = luaL_checkinteger(L, 1);
            int b = luaL_checkinteger(L, 2);
            
            // Calculate GCD first
            int gcd = a;
            int tempB = b;
            while (tempB != 0) {
                int temp = tempB;
                tempB = gcd % tempB;
                gcd = temp;
            }
            
            // LCM = (a * b) / gcd
            lua_Integer lcm = (lua_Integer)a * b / gcd;
            
            lua_pushinteger(L, lcm);
            return 1;
        });
        lua_setfield(m_luaState, -2, "lcm");
        
        // Clean up the stack
        lua_pop(m_luaState, 1);
        
        // Create a calculator library
        lua_newtable(m_luaState);
        
        // Add calculator functions
        
        // calculator.setVariable
        lua_pushcfunction(m_luaState, [](lua_State* L) -> int {
            // Get the calculator instance
            lua_getfield(L, LUA_REGISTRYINDEX, "calculator");
            Calculator* calculator = static_cast<Calculator*>(lua_touserdata(L, -1));
            lua_pop(L, 1);
            
            // Get the arguments
            const char* name = luaL_checkstring(L, 1);
            double value = luaL_checknumber(L, 2);
            
            // Set the variable
            bool success = calculator->setVariable(name, value);
            
            // Return success or failure
            lua_pushboolean(L, success);
            return 1;
        });
        lua_setfield(m_luaState, -2, "setVariable");
        
        // calculator.getVariable
        lua_pushcfunction(m_luaState, [](lua_State* L) -> int {
            // Get the calculator instance
            lua_getfield(L, LUA_REGISTRYINDEX, "calculator");
            Calculator* calculator = static_cast<Calculator*>(lua_touserdata(L, -1));
            lua_pop(L, 1);
            
            // Get the arguments
            const char* name = luaL_checkstring(L, 1);
            
            // Get the variable
            auto value = calculator->getVariable(name);
            
            // Return the value or nil
            if (value) {
                lua_pushnumber(L, *value);
            } else {
                lua_pushnil(L);
            }
            
            return 1;
        });
        lua_setfield(m_luaState, -2, "getVariable");
        
        // calculator.hasVariable
        lua_pushcfunction(m_luaState, [](lua_State* L) -> int {
            // Get the calculator instance
            lua_getfield(L, LUA_REGISTRYINDEX, "calculator");
            Calculator* calculator = static_cast<Calculator*>(lua_touserdata(L, -1));
            lua_pop(L, 1);
            
            // Get the arguments
            const char* name = luaL_checkstring(L, 1);
            
            // Check if the variable exists
            bool exists = calculator->hasVariable(name);
            
            // Return the result
            lua_pushboolean(L, exists);
            return 1;
        });
        lua_setfield(m_luaState, -2, "hasVariable");
        
        // calculator.clearVariables
        lua_pushcfunction(m_luaState, [](lua_State* L) -> int {
            // Get the calculator instance
            lua_getfield(L, LUA_REGISTRYINDEX, "calculator");
            Calculator* calculator = static_cast<Calculator*>(lua_touserdata(L, -1));
            lua_pop(L, 1);
            
            // Clear the variables
            calculator->clearVariables();
            
            return 0;
        });
        lua_setfield(m_luaState, -2, "clearVariables");
        
        // Set the calculator library as a global
        lua_setglobal(m_luaState, "calculator");
        
        // Store the calculator instance in the registry
        lua_pushlightuserdata(m_luaState, m_calculator.get());
        lua_setfield(m_luaState, LUA_REGISTRYINDEX, "calculator");
        
        return true;
    }
    
    static int luaFunctionDispatcher(lua_State* L) {
        // Get the Impl instance and function name from the closure
        Impl* impl = static_cast<Impl*>(lua_touserdata(L, lua_upvalueindex(1)));
        const char* functionName = lua_tostring(L, lua_upvalueindex(2));
        
        // Get the arguments
        std::vector<std::string> args;
        int numArgs = lua_gettop(L);
        for (int i = 1; i <= numArgs; ++i) {
            const char* arg = lua_tostring(L, i);
            if (arg) {
                args.push_back(arg);
            } else {
                args.push_back("");
            }
        }
        
        // Call the C++ function
        auto it = impl->m_functions.find(functionName);
        if (it != impl->m_functions.end()) {
            std::string result = it->second(args);
            lua_pushstring(L, result.c_str());
            return 1;
        }
        
        // Function not found
        lua_pushstring(L, "Function not found");
        return 1;
    }

private:
    std::shared_ptr<Calculator> m_calculator;
    lua_State* m_luaState;
    std::unordered_map<std::string, std::function<std::string(const std::vector<std::string>&)>> m_functions;
};

// LuaEngine implementation

LuaEngine::LuaEngine(std::shared_ptr<Calculator> calculator)
    : m_calculator(calculator), m_impl(std::make_unique<Impl>(calculator)) {
}

LuaEngine::~LuaEngine() {
    shutdown();
}

bool LuaEngine::initialize() {
    try {
        return m_impl->initialize();
    } catch (const std::exception& e) {
        std::cerr << "Exception during Lua engine initialization: " << e.what() << std::endl;
        return false;
    }
}

void LuaEngine::shutdown() {
    if (m_impl) {
        m_impl->shutdown();
    }
}

std::optional<std::string> LuaEngine::executeScript(const std::string& script) {
    try {
        return m_impl->executeScript(script);
    } catch (const std::exception& e) {
        std::cerr << "Error executing Lua script: " << e.what() << std::endl;
        return std::nullopt;
    }
}

bool LuaEngine::loadScript(const std::string& filename) {
    try {
        return m_impl->loadScript(filename);
    } catch (const std::exception& e) {
        std::cerr << "Error loading Lua script: " << e.what() << std::endl;
        return false;
    }
}

void LuaEngine::registerFunction(const std::string& name, 
                               std::function<std::string(const std::vector<std::string>&)> function) {
    m_impl->registerFunction(name, function);
}

std::vector<std::string> LuaEngine::getAvailableFunctions() const {
    return m_impl->getAvailableFunctions();
}

std::optional<std::string> LuaEngine::getFunctionHelp(const std::string& functionName) const {
    return m_impl->getFunctionHelp(functionName);
}

bool LuaEngine::setupLuaEnvironment() {
    return m_impl->setupLuaEnvironment();
}

void LuaEngine::registerBuiltinFunctions() {
    // This is handled in the Impl class
}

std::string LuaEngine::formatLuaError(const std::string& error) const {
    // This is a placeholder for a real error formatter
    return "Lua error: " + error;
}

} // namespace rebelcalc
