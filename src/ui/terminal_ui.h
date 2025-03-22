#pragma once

#include <string>
#include <memory>
#include <vector>
#include <functional>
#include <unordered_map>
#include "../backend/matrix.h"

namespace rebelcalc {

// Forward declarations
class Calculator;
class LuaEngine;
class InputProcessor;
class ThemeManager;

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
    
    /**
     * Set whether to enable syntax highlighting
     * @param enable Whether to enable syntax highlighting
     */
    void enableSyntaxHighlighting(bool enable);
    
    /**
     * Get whether syntax highlighting is enabled
     * @return Whether syntax highlighting is enabled
     */
    bool isSyntaxHighlightingEnabled() const;
    
    /**
     * Set whether to enable autocompletion
     * @param enable Whether to enable autocompletion
     */
    void enableAutocompletion(bool enable);
    
    /**
     * Get whether autocompletion is enabled
     * @return Whether autocompletion is enabled
     */
    bool isAutocompletionEnabled() const;
    
    /**
     * Set whether to enable multi-line input
     * @param enable Whether to enable multi-line input
     */
    void enableMultiLineInput(bool enable);
    
    /**
     * Get whether multi-line input is enabled
     * @return Whether multi-line input is enabled
     */
    bool isMultiLineInputEnabled() const;
    
    /**
     * Create a new workspace
     * @param name The name of the workspace
     * @return true if the workspace was created successfully, false otherwise
     */
    bool createWorkspace(const std::string& name);
    
    /**
     * Switch to a workspace
     * @param name The name of the workspace to switch to
     * @return true if the workspace was switched to successfully, false otherwise
     */
    bool switchWorkspace(const std::string& name);
    
    /**
     * Get the current workspace
     * @return The name of the current workspace
     */
    std::string getCurrentWorkspace() const;
    
    /**
     * Get a list of available workspaces
     * @return A vector of workspace names
     */
    std::vector<std::string> getAvailableWorkspaces() const;
    
    /**
     * Split the current view
     * @param direction The direction to split (horizontal or vertical)
     * @return true if the view was split successfully, false otherwise
     */
    bool splitView(const std::string& direction);
    
    /**
     * Close the current split view
     * @return true if the view was closed successfully, false otherwise
     */
    bool closeSplitView();

private:
    // Implementation details
    class Impl;
    std::unique_ptr<Impl> m_impl;
    
    // Components
    std::shared_ptr<Calculator> m_calculator;
    std::shared_ptr<LuaEngine> m_luaEngine;
    std::shared_ptr<InputProcessor> m_inputProcessor;
    std::shared_ptr<ThemeManager> m_themeManager;
    
    // Helper methods
    bool processCommand(const std::string& command);
    void displayHelp() const;
    void displayResult(const std::string& result) const;
    void displayError(const std::string& error) const;
};

} // namespace rebelcalc
