#pragma once

#include <string>
#include <memory>
#include <functional>
#include <optional>
#include <vector>
#include <unordered_map>

namespace rebelcalc {

// Forward declarations
class Calculator;

/**
 * Class for handling Lua scripting functionality
 */
class LuaEngine {
public:
    /**
     * Constructor
     * @param calculator The calculator instance
     */
    LuaEngine(std::shared_ptr<Calculator> calculator);
    
    /**
     * Destructor
     */
    ~LuaEngine();
    
    /**
     * Initialize the Lua engine
     * @return true if initialization was successful, false otherwise
     */
    bool initialize();
    
    /**
     * Shutdown the Lua engine
     */
    void shutdown();
    
    /**
     * Execute a Lua script
     * @param script The script to execute
     * @return The result of the script execution, or nullopt if execution failed
     */
    std::optional<std::string> executeScript(const std::string& script);
    
    /**
     * Load a Lua script from a file
     * @param filename The file to load
     * @return true if the script was loaded successfully, false otherwise
     */
    bool loadScript(const std::string& filename);
    
    /**
     * Register a C++ function to be callable from Lua
     * @param name The function name in Lua
     * @param function The C++ function to register
     */
    void registerFunction(const std::string& name, 
                         std::function<std::string(const std::vector<std::string>&)> function);
    
    /**
     * Get a list of available Lua functions
     * @return A vector of function names
     */
    std::vector<std::string> getAvailableFunctions() const;
    
    /**
     * Get help text for a Lua function
     * @param functionName The function name
     * @return The help text, or nullopt if the function doesn't exist
     */
    std::optional<std::string> getFunctionHelp(const std::string& functionName) const;

private:
    // Implementation details
    class Impl;
    std::unique_ptr<Impl> m_impl;
    
    // Components
    std::shared_ptr<Calculator> m_calculator;
    
    // Helper methods
    bool setupLuaEnvironment();
    void registerBuiltinFunctions();
    std::string formatLuaError(const std::string& error) const;
};

} // namespace rebelcalc
