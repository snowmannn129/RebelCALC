#pragma once

#include <string>
#include <memory>
#include <vector>
#include <functional>

namespace rebelcalc {

// Forward declarations
class Calculator;
class LuaEngine;

/**
 * Class for handling the terminal-based user interface
 */
class TerminalUI {
public:
    /**
     * Constructor
     * @param calculator The calculator instance
     * @param luaEngine The Lua scripting engine instance
     */
    TerminalUI(std::shared_ptr<Calculator> calculator, std::shared_ptr<LuaEngine> luaEngine);
    
    /**
     * Destructor
     */
    ~TerminalUI();
    
    /**
     * Initialize the terminal UI
     * @return true if initialization was successful, false otherwise
     */
    bool initialize();
    
    /**
     * Shutdown the terminal UI
     */
    void shutdown();
    
    /**
     * Run the main UI loop
     */
    void run();
    
    /**
     * Set the UI theme
     * @param themeName The name of the theme to set
     * @return true if the theme was set successfully, false otherwise
     */
    bool setTheme(const std::string& themeName);
    
    /**
     * Get the current UI theme
     * @return The name of the current theme
     */
    std::string getTheme() const;
    
    /**
     * Register a custom command
     * @param name The command name
     * @param handler The command handler function
     * @param helpText The help text for the command
     */
    void registerCommand(const std::string& name, 
                        std::function<bool(const std::vector<std::string>&)> handler,
                        const std::string& helpText);

private:
    // Implementation details
    class Impl;
    std::unique_ptr<Impl> m_impl;
    
    // Components
    std::shared_ptr<Calculator> m_calculator;
    std::shared_ptr<LuaEngine> m_luaEngine;
    
    // Helper methods
    bool processCommand(const std::string& command);
    void displayHelp() const;
    void displayResult(const std::string& result) const;
    void displayError(const std::string& error) const;
};

} // namespace rebelcalc
