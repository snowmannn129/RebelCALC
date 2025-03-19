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
