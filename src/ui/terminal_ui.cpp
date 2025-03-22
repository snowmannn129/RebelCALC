#include "terminal_ui.h"
#include "input_processor.h"
#include "theme_manager.h"
#include "../backend/calculator.h"
#include "../scripts/lua_engine.h"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <regex>

namespace rebelcalc {

// Private implementation class
class TerminalUI::Impl {
private:
    // Components
    std::shared_ptr<Calculator> m_calculator;
    std::shared_ptr<LuaEngine> m_luaEngine;
    std::shared_ptr<InputProcessor> m_inputProcessor;
    std::shared_ptr<ThemeManager> m_themeManager;
    
    // State
    bool m_running;
    std::string m_theme;
    std::string m_currentWorkspace;
    std::vector<std::string> m_workspaces;
    bool m_splitViewActive;
    std::string m_splitViewDirection;
    
    // Commands
    struct CommandInfo {
        std::function<bool(const std::vector<std::string>&)> handler;
        std::string helpText;
    };
    
    std::unordered_map<std::string, CommandInfo> m_commands;
    
    // Helper methods
    void setupCompletionProviders();
    void setupSyntaxHighlightingRules();
    void displayVariables() const;
    void displayHelp() const;
    bool processInput(const std::string& input);
    bool processCommand(const std::string& commandLine);

public:
    Impl(std::shared_ptr<Calculator> calculator, std::shared_ptr<LuaEngine> luaEngine);
    ~Impl();
    
    bool initialize();
    void shutdown();
    void run();
    
    bool setTheme(const std::string& themeName);
    std::string getTheme() const;
    
    void enableSyntaxHighlighting(bool enable);
    bool isSyntaxHighlightingEnabled() const;
    
    void enableAutocompletion(bool enable);
    bool isAutocompletionEnabled() const;
    
    void enableMultiLineInput(bool enable);
    bool isMultiLineInputEnabled() const;
    
    bool createWorkspace(const std::string& name);
    bool switchWorkspace(const std::string& name);
    std::string getCurrentWorkspace() const;
    std::vector<std::string> getAvailableWorkspaces() const;
    
    bool splitView(const std::string& direction);
    bool closeSplitView();
    
    void registerCommand(const std::string& name, 
                        std::function<bool(const std::vector<std::string>&)> handler,
                        const std::string& helpText);
};

// TerminalUI implementation
TerminalUI::TerminalUI(std::shared_ptr<Calculator> calculator, std::shared_ptr<LuaEngine> luaEngine)
    : m_impl(std::make_unique<Impl>(calculator, luaEngine)),
      m_calculator(calculator), m_luaEngine(luaEngine) {
}

TerminalUI::~TerminalUI() = default;

bool TerminalUI::initialize() {
    return m_impl->initialize();
}

void TerminalUI::shutdown() {
    m_impl->shutdown();
}

void TerminalUI::run() {
    m_impl->run();
}

bool TerminalUI::setTheme(const std::string& themeName) {
    return m_impl->setTheme(themeName);
}

std::string TerminalUI::getTheme() const {
    return m_impl->getTheme();
}

void TerminalUI::enableSyntaxHighlighting(bool enable) {
    m_impl->enableSyntaxHighlighting(enable);
}

bool TerminalUI::isSyntaxHighlightingEnabled() const {
    return m_impl->isSyntaxHighlightingEnabled();
}

void TerminalUI::enableAutocompletion(bool enable) {
    m_impl->enableAutocompletion(enable);
}

bool TerminalUI::isAutocompletionEnabled() const {
    return m_impl->isAutocompletionEnabled();
}

void TerminalUI::enableMultiLineInput(bool enable) {
    m_impl->enableMultiLineInput(enable);
}

bool TerminalUI::isMultiLineInputEnabled() const {
    return m_impl->isMultiLineInputEnabled();
}

bool TerminalUI::createWorkspace(const std::string& name) {
    return m_impl->createWorkspace(name);
}

bool TerminalUI::switchWorkspace(const std::string& name) {
    return m_impl->switchWorkspace(name);
}

std::string TerminalUI::getCurrentWorkspace() const {
    return m_impl->getCurrentWorkspace();
}

std::vector<std::string> TerminalUI::getAvailableWorkspaces() const {
    return m_impl->getAvailableWorkspaces();
}

bool TerminalUI::splitView(const std::string& direction) {
    return m_impl->splitView(direction);
}

bool TerminalUI::closeSplitView() {
    return m_impl->closeSplitView();
}

void TerminalUI::registerCommand(const std::string& name, 
                               std::function<bool(const std::vector<std::string>&)> handler,
                               const std::string& helpText) {
    m_impl->registerCommand(name, handler, helpText);
}

bool TerminalUI::processCommand(const std::string& command) {
    // This is a placeholder. The actual implementation is in the Impl class.
    return false;
}

void TerminalUI::displayHelp() const {
    // This is a placeholder. The actual implementation is in the Impl class.
}

void TerminalUI::displayResult(const std::string& result) const {
    // This is a placeholder. The actual implementation is in the Impl class.
}

void TerminalUI::displayError(const std::string& error) const {
    // This is a placeholder. The actual implementation is in the Impl class.
}

// TerminalUI::Impl implementation

TerminalUI::Impl::Impl(std::shared_ptr<Calculator> calculator, std::shared_ptr<LuaEngine> luaEngine)
    : m_calculator(calculator), m_luaEngine(luaEngine), 
      m_inputProcessor(std::make_shared<InputProcessor>()), 
      m_themeManager(std::make_shared<ThemeManager>()),
      m_running(false), m_theme("default"), m_currentWorkspace("default"),
      m_splitViewActive(false), m_splitViewDirection("horizontal") {
    
    // Initialize workspaces
    m_workspaces.push_back("default");
    
    // Register built-in commands
    registerCommand("help", [this](const std::vector<std::string>& args) {
        displayHelp();
        return true;
    }, "Display help information");
    
    registerCommand("exit", [this](const std::vector<std::string>& args) {
        m_running = false;
        return true;
    }, "Exit the application");
    
    registerCommand("clear", [this](const std::vector<std::string>& args) {
        // Clear the screen (platform-dependent)
        #ifdef _WIN32
        system("cls");
        #else
        system("clear");
        #endif
        return true;
    }, "Clear the screen");
    
    registerCommand("theme", [this](const std::vector<std::string>& args) {
        if (args.size() < 1) {
            std::cout << "Current theme: " << m_theme << std::endl;
            std::cout << "Available themes: ";
            auto themes = m_themeManager->getAvailableThemes();
            for (size_t i = 0; i < themes.size(); ++i) {
                std::cout << themes[i];
                if (i < themes.size() - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << std::endl;
            return true;
        }
        
        return setTheme(args[0]);
    }, "Set or display the current theme");
    
    registerCommand("vars", [this](const std::vector<std::string>& args) {
        displayVariables();
        return true;
    }, "Display all defined variables");
    
    registerCommand("clear_vars", [this](const std::vector<std::string>& args) {
        m_calculator->clearVariables();
        std::cout << "Variables cleared" << std::endl;
        return true;
    }, "Clear all user-defined variables");
    
    registerCommand("syntax", [this](const std::vector<std::string>& args) {
        if (args.size() < 1) {
            std::cout << "Syntax highlighting is " << (m_inputProcessor->isSyntaxHighlightingEnabled() ? "enabled" : "disabled") << std::endl;
            return true;
        }
        
        if (args[0] == "on" || args[0] == "enable") {
            m_inputProcessor->enableSyntaxHighlighting(true);
            std::cout << "Syntax highlighting enabled" << std::endl;
        } else if (args[0] == "off" || args[0] == "disable") {
            m_inputProcessor->enableSyntaxHighlighting(false);
            std::cout << "Syntax highlighting disabled" << std::endl;
        } else {
            std::cout << "Usage: syntax [on|off]" << std::endl;
            return false;
        }
        
        return true;
    }, "Enable or disable syntax highlighting");
    
    registerCommand("autocomplete", [this](const std::vector<std::string>& args) {
        if (args.size() < 1) {
            std::cout << "Autocompletion is " << (m_inputProcessor->isAutocompletionEnabled() ? "enabled" : "disabled") << std::endl;
            return true;
        }
        
        if (args[0] == "on" || args[0] == "enable") {
            m_inputProcessor->enableAutocompletion(true);
            std::cout << "Autocompletion enabled" << std::endl;
        } else if (args[0] == "off" || args[0] == "disable") {
            m_inputProcessor->enableAutocompletion(false);
            std::cout << "Autocompletion disabled" << std::endl;
        } else {
            std::cout << "Usage: autocomplete [on|off]" << std::endl;
            return false;
        }
        
        return true;
    }, "Enable or disable autocompletion");
    
    registerCommand("multiline", [this](const std::vector<std::string>& args) {
        if (args.size() < 1) {
            std::cout << "Multi-line input is " << (m_inputProcessor->isMultiLineInputEnabled() ? "enabled" : "disabled") << std::endl;
            return true;
        }
        
        if (args[0] == "on" || args[0] == "enable") {
            m_inputProcessor->enableMultiLineInput(true);
            std::cout << "Multi-line input enabled" << std::endl;
        } else if (args[0] == "off" || args[0] == "disable") {
            m_inputProcessor->enableMultiLineInput(false);
            std::cout << "Multi-line input disabled" << std::endl;
        } else {
            std::cout << "Usage: multiline [on|off]" << std::endl;
            return false;
        }
        
        return true;
    }, "Enable or disable multi-line input");
    
    registerCommand("workspace", [this](const std::vector<std::string>& args) {
        if (args.size() < 1) {
            std::cout << "Current workspace: " << m_currentWorkspace << std::endl;
            std::cout << "Available workspaces: ";
            for (size_t i = 0; i < m_workspaces.size(); ++i) {
                std::cout << m_workspaces[i];
                if (i < m_workspaces.size() - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << std::endl;
            return true;
        }
        
        if (args[0] == "create") {
            if (args.size() < 2) {
                std::cout << "Usage: workspace create <name>" << std::endl;
                return false;
            }
            
            return createWorkspace(args[1]);
        } else if (args[0] == "switch") {
            if (args.size() < 2) {
                std::cout << "Usage: workspace switch <name>" << std::endl;
                return false;
            }
            
            return switchWorkspace(args[1]);
        } else {
            std::cout << "Usage: workspace [create|switch] <name>" << std::endl;
            return false;
        }
    }, "Manage workspaces");
    
    registerCommand("split", [this](const std::vector<std::string>& args) {
        if (args.size() < 1) {
            std::cout << "Usage: split [horizontal|vertical]" << std::endl;
            return false;
        }
        
        return splitView(args[0]);
    }, "Split the current view");
    
    registerCommand("close", [this](const std::vector<std::string>& args) {
        return closeSplitView();
    }, "Close the current split view");
    
    registerCommand("solve", [this](const std::vector<std::string>& args) {
        if (args.size() < 2) {
            std::cout << "Usage: solve <equation> <variable>" << std::endl;
            return false;
        }
        
        auto result = m_calculator->solve(args[0], args[1]);
        if (!result) {
            std::cout << "Failed to solve equation" << std::endl;
            return false;
        }
        
        std::visit([](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, double>) {
                std::cout << "Solution: " << arg << std::endl;
            } else if constexpr (std::is_same_v<T, std::string>) {
                std::cout << "Solution: " << arg << std::endl;
            } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                std::cout << "Solutions: ";
                for (size_t i = 0; i < arg.size(); ++i) {
                    std::cout << arg[i];
                    if (i < arg.size() - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << std::endl;
            }
        }, *result);
        
        return true;
    }, "Solve an equation for a variable");
    
    registerCommand("diff", [this](const std::vector<std::string>& args) {
        if (args.size() < 2) {
            std::cout << "Usage: diff <expression> <variable>" << std::endl;
            return false;
        }
        
        auto result = m_calculator->differentiate(args[0], args[1]);
        if (!result) {
            std::cout << "Failed to differentiate expression" << std::endl;
            return false;
        }
        
        std::cout << "Derivative: " << *result << std::endl;
        return true;
    }, "Differentiate an expression with respect to a variable");
    
    registerCommand("int", [this](const std::vector<std::string>& args) {
        if (args.size() < 2) {
            std::cout << "Usage: int <expression> <variable>" << std::endl;
            return false;
        }
        
        auto result = m_calculator->integrate(args[0], args[1]);
        if (!result) {
            std::cout << "Failed to integrate expression" << std::endl;
            return false;
        }
        
        std::cout << "Integral: " << *result << std::endl;
        return true;
    }, "Integrate an expression with respect to a variable");
    
    registerCommand("lua", [this](const std::vector<std::string>& args) {
        if (args.empty()) {
            std::cout << "Usage: lua <script>" << std::endl;
            return false;
        }
        
        std::string script = args[0];
        for (size_t i = 1; i < args.size(); ++i) {
            script += " " + args[i];
        }
        
        auto result = m_luaEngine->executeScript(script);
        if (!result) {
            std::cout << "Failed to execute Lua script" << std::endl;
            return false;
        }
        
        std::cout << *result << std::endl;
        return true;
    }, "Execute a Lua script");
    
    registerCommand("load", [this](const std::vector<std::string>& args) {
        if (args.size() < 1) {
            std::cout << "Usage: load <filename>" << std::endl;
            return false;
        }
        
        if (!m_luaEngine->loadScript(args[0])) {
            std::cout << "Failed to load script from file" << std::endl;
            return false;
        }
        
        std::cout << "Script loaded successfully" << std::endl;
        return true;
    }, "Load a Lua script from a file");
    
    // Matrix-specific commands
    registerCommand("matrix", [this](const std::vector<std::string>& args) {
        if (args.size() < 1) {
            std::cout << "Usage: matrix <subcommand> [args...]" << std::endl;
            std::cout << "Subcommands: create, identity, diagonal, random, transpose, inverse, det, trace, rank" << std::endl;
            return false;
        }
        
        std::string subcommand = args[0];
        
        if (subcommand == "create") {
            if (args.size() < 4) {
                std::cout << "Usage: matrix create <rows> <cols> <value>" << std::endl;
                return false;
            }
            
            try {
                size_t rows = std::stoul(args[1]);
                size_t cols = std::stoul(args[2]);
                double value = std::stod(args[3]);
                
                Matrix m(rows, cols, value);
                std::cout << std::endl << m << std::endl;
                return true;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
                return false;
            }
        } else if (subcommand == "identity") {
            if (args.size() < 2) {
                std::cout << "Usage: matrix identity <size>" << std::endl;
                return false;
            }
            
            try {
                size_t size = std::stoul(args[1]);
                Matrix m = Matrix::identity(size);
                std::cout << std::endl << m << std::endl;
                return true;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
                return false;
            }
        } else if (subcommand == "diagonal") {
            if (args.size() < 2) {
                std::cout << "Usage: matrix diagonal <value1> <value2> ..." << std::endl;
                return false;
            }
            
            try {
                std::vector<double> diagonal;
                for (size_t i = 1; i < args.size(); ++i) {
                    diagonal.push_back(std::stod(args[i]));
                }
                
                Matrix m = Matrix::diagonal(diagonal);
                std::cout << std::endl << m << std::endl;
                return true;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
                return false;
            }
        } else if (subcommand == "random") {
            if (args.size() < 4) {
                std::cout << "Usage: matrix random <rows> <cols> <min> <max>" << std::endl;
                return false;
            }
            
            try {
                size_t rows = std::stoul(args[1]);
                size_t cols = std::stoul(args[2]);
                double min = std::stod(args[3]);
                double max = (args.size() > 4) ? std::stod(args[4]) : 1.0;
                
                Matrix m = Matrix::random(rows, cols, min, max);
                std::cout << std::endl << m << std::endl;
                return true;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
                return false;
            }
        } else if (subcommand == "transpose") {
            std::cout << "Use the expression 'matrix.transpose()' to transpose a matrix" << std::endl;
            return true;
        } else if (subcommand == "inverse") {
            std::cout << "Use the expression 'matrix.inverse()' to invert a matrix" << std::endl;
            return true;
        } else if (subcommand == "det") {
            std::cout << "Use the expression 'matrix.determinant()' to calculate the determinant of a matrix" << std::endl;
            return true;
        } else if (subcommand == "trace") {
            std::cout << "Use the expression 'matrix.trace()' to calculate the trace of a matrix" << std::endl;
            return true;
        } else if (subcommand == "rank") {
            std::cout << "Use the expression 'matrix.rank()' to calculate the rank of a matrix" << std::endl;
            return true;
        } else {
            std::cout << "Unknown subcommand: " << subcommand << std::endl;
            return false;
        }
    }, "Matrix operations");
    
    // Set up completion providers
    setupCompletionProviders();
}

TerminalUI::Impl::~Impl() {
    shutdown();
}

bool TerminalUI::Impl::initialize() {
    // Initialize the terminal UI implementation
    if (!m_inputProcessor->initialize()) {
        std::cerr << "Failed to initialize input processor" << std::endl;
        return false;
    }
    
    if (!m_themeManager->initialize()) {
        std::cerr << "Failed to initialize theme manager" << std::endl;
        return false;
    }
    
    // Set up syntax highlighting rules
    setupSyntaxHighlightingRules();
    
    return true;
}

void TerminalUI::Impl::shutdown() {
    // Shutdown the terminal UI implementation
    m_inputProcessor->shutdown();
}

void TerminalUI::Impl::run() {
    m_running = true;
    
    while (m_running) {
        // Display prompt with theme
        std::string prompt = m_themeManager->formatText("> ", "prompt");
        
        // Read input using InputProcessor
        std::string input = m_inputProcessor->readLine(prompt);
        
        // Process input
        if (!input.empty()) {
            processInput(input);
        }
    }
}

bool TerminalUI::Impl::setTheme(const std::string& themeName) {
    // Set the theme using ThemeManager
    if (!m_themeManager->setCurrentTheme(themeName)) {
        std::cout << "Unknown theme: " << themeName << std::endl;
        return false;
    }
    
    m_theme = themeName;
    std::cout << "Theme set to " << themeName << std::endl;
    return true;
}

std::string TerminalUI::Impl::getTheme() const {
    return m_theme;
}

void TerminalUI::Impl::enableSyntaxHighlighting(bool enable) {
    m_inputProcessor->enableSyntaxHighlighting(enable);
}

bool TerminalUI::Impl::isSyntaxHighlightingEnabled() const {
    return m_inputProcessor->isSyntaxHighlightingEnabled();
}

void TerminalUI::Impl::enableAutocompletion(bool enable) {
    m_inputProcessor->enableAutocompletion(enable);
}

bool TerminalUI::Impl::isAutocompletionEnabled() const {
    return m_inputProcessor->isAutocompletionEnabled();
}

void TerminalUI::Impl::enableMultiLineInput(bool enable) {
    m_inputProcessor->enableMultiLineInput(enable);
}

bool TerminalUI::Impl::isMultiLineInputEnabled() const {
    return m_inputProcessor->isMultiLineInputEnabled();
}

std::string TerminalUI::Impl::getCurrentWorkspace() const {
    return m_currentWorkspace;
}

std::vector<std::string> TerminalUI::Impl::getAvailableWorkspaces() const {
    return m_workspaces;
}

void TerminalUI::Impl::registerCommand(const std::string& name, 
                                     std::function<bool(const std::vector<std::string>&)> handler,
                                     const std::string& helpText) {
    m_commands[name] = {handler, helpText};
}

bool TerminalUI::Impl::createWorkspace(const std::string& name) {
    // Check if the workspace already exists
    if (std::find(m_workspaces.begin(), m_workspaces.end(), name) != m_workspaces.end()) {
        std::cout << "Workspace '" << name << "' already exists" << std::endl;
        return false;
    }
    
    // Create the workspace
    m_workspaces.push_back(name);
    std::cout << "Workspace '" << name << "' created" << std::endl;
    
    // Switch to the new workspace
    return switchWorkspace(name);
}

bool TerminalUI::Impl::switchWorkspace(const std::string& name) {
    // Check if the workspace exists
    if (std::find(m_workspaces.begin(), m_workspaces.end(), name) == m_workspaces.end()) {
        std::cout << "Workspace '" << name << "' does not exist" << std::endl;
        return false;
    }
    
    // Switch to the workspace
    m_currentWorkspace = name;
    std::cout << "Switched to workspace '" << name << "'" << std::endl;
    return true;
}

bool TerminalUI::Impl::splitView(const std::string& direction) {
    // Check if the direction is valid
    if (direction != "horizontal" && direction != "vertical") {
        std::cout << "Invalid direction: " << direction << std::endl;
        std::cout << "Valid directions: horizontal, vertical" << std::endl;
        return false;
    }
    
    // Set the split view
    m_splitViewActive = true;
    m_splitViewDirection = direction;
    std::cout << "Split view " << direction << " activated" << std::endl;
    return true;
}

bool TerminalUI::Impl::closeSplitView() {
    // Check if a split view is active
    if (!m_splitViewActive) {
        std::cout << "No split view is active" << std::endl;
        return false;
    }
    
    // Close the split view
    m_splitViewActive = false;
    std::cout << "Split view closed" << std::endl;
    return true;
}

void TerminalUI::Impl::setupCompletionProviders() {
    // Command completion provider
    m_inputProcessor->addCompletionProvider("commands", [this](const std::string& input) {
        std::vector<std::string> completions;
        for (const auto& [name, _] : m_commands) {
            if (name.find(input) == 0) {
                completions.push_back(name);
            }
        }
        return completions;
    });
    
    // Variable completion provider
    m_inputProcessor->addCompletionProvider("variables", [this](const std::string& input) {
        std::vector<std::string> completions;
        auto variables = m_calculator->getAllVariables();
        for (const auto& [name, _] : variables) {
            if (name.find(input) == 0) {
                completions.push_back(name);
            }
        }
        return completions;
    });
    
    // Function completion provider
    m_inputProcessor->addCompletionProvider("functions", [](const std::string& input) {
        std::vector<std::string> functions = {
            "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh",
            "sqrt", "log", "log10", "exp", "abs", "floor", "ceil", "round"
        };
        
        std::vector<std::string> completions;
        for (const auto& func : functions) {
            if (func.find(input) == 0) {
                completions.push_back(func);
            }
        }
        return completions;
    });
}

void TerminalUI::Impl::setupSyntaxHighlightingRules() {
    // Clear existing rules
    m_inputProcessor->clearSyntaxHighlightingRules();
    
    // Add syntax highlighting rules
    m_inputProcessor->addSyntaxHighlightingRule("\\b(sin|cos|tan|asin|acos|atan|sinh|cosh|tanh|sqrt|log|log10|exp|abs|floor|ceil|round)\\b", "info");
    m_inputProcessor->addSyntaxHighlightingRule("\\b(pi|e)\\b", "info");
    m_inputProcessor->addSyntaxHighlightingRule("\\b(if|else|for|while|function|return)\\b", "prompt");
    m_inputProcessor->addSyntaxHighlightingRule("\\b(true|false|null)\\b", "warning");
    m_inputProcessor->addSyntaxHighlightingRule("\\b[0-9]+(\\.[0-9]+)?\\b", "result");
    m_inputProcessor->addSyntaxHighlightingRule("\\b[A-Za-z_][A-Za-z0-9_]*\\s*=", "prompt");
    m_inputProcessor->addSyntaxHighlightingRule("\"[^\"]*\"", "warning");
    m_inputProcessor->addSyntaxHighlightingRule("'[^']*'", "warning");
    m_inputProcessor->addSyntaxHighlightingRule("//.*$", "info");
    m_inputProcessor->addSyntaxHighlightingRule("/\\*.*?\\*/", "info");
}

void TerminalUI::Impl::displayVariables() const {
    std::cout << m_themeManager->formatText("Variables:", "info") << std::endl;
    
    // Get all variables
    auto variables = m_calculator->getAllVariables();
    
    // Built-in constants
    std::cout << "  Built-in constants:" << std::endl;
    std::cout << "    pi = " << m_themeManager->formatText(std::to_string(variables["pi"]), "result") << std::endl;
    std::cout << "    e = " << m_themeManager->formatText(std::to_string(variables["e"]), "result") << std::endl;
    
    // User-defined variables
    std::cout << "  User-defined variables:" << std::endl;
    
    bool hasUserVariables = false;
    for (const auto& [name, value] : variables) {
        if (name != "pi" && name != "e") {
            std::cout << "    " << name << " = " << m_themeManager->formatText(std::to_string(value), "result") << std::endl;
            hasUserVariables = true;
        }
    }
    
    if (!hasUserVariables) {
        std::cout << "    (No user-defined variables)" << std::endl;
    }
}

void TerminalUI::Impl::displayHelp() const {
    std::cout << m_themeManager->formatText("RebelCALC - Advanced Computational Engine for RebelSUITE", "info") << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;
    
    std::cout << m_themeManager->formatText("Available commands:", "info") << std::endl;
    std::cout << "------------------" << std::endl;
    
    // Calculate the maximum command length for alignment
    size_t maxCommandLength = 0;
    for (const auto& [name, _] : m_commands) {
        maxCommandLength = std::max(maxCommandLength, name.length());
    }
    
    // Display commands in alphabetical order
    std::vector<std::string> commandNames;
    for (const auto& [name, _] : m_commands) {
        commandNames.push_back(name);
    }
    
    std::sort(commandNames.begin(), commandNames.end());
    
    for (const auto& name : commandNames) {
        const auto& command = m_commands.at(name);
        std::cout << "  /" << std::left << std::setw(maxCommandLength + 2) << name
                  << command.helpText << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << m_themeManager->formatText("Expression syntax:", "info") << std::endl;
    std::cout << "-----------------" << std::endl;
    std::cout << "  Basic arithmetic: +, -, *, /, ^ (power)" << std::endl;
    std::cout << "  Parentheses: ( )" << std::endl;
    std::cout << "  Functions: sin, cos, tan, sqrt, log, etc." << std::endl;
    std::cout << "  Variables: x, y, etc." << std::endl;
    std::cout << "  Constants: pi, e" << std::endl;
    
    std::cout << std::endl;
    std::cout << m_themeManager->formatText("Variable assignments:", "info") << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "  Simple assignment: x = 5" << std::endl;
    std::cout << "  Compound assignments: x += 3, y *= 2, z -= 1, w /= 4" << std::endl;
    std::cout << "  Assignment with expressions: x = 2 * sin(pi/4) + 3" << std::endl;
    
    std::cout << std::endl;
    std::cout << m_themeManager->formatText("Matrix operations:", "info") << std::endl;
    std::cout << "-----------------" << std::endl;
    std::cout << "  Matrix creation:" << std::endl;
    std::cout << "    /matrix create <rows> <cols> <value>   # Create a matrix filled with a value" << std::endl;
    std::cout << "    /matrix identity <size>                # Create an identity matrix" << std::endl;
    std::cout << "    /matrix diagonal <val1> <val2> ...     # Create a diagonal matrix" << std::endl;
    std::cout << "    /matrix random <rows> <cols> <min> <max> # Create a random matrix" << std::endl;
    std::cout << std::endl;
    std::cout << "  Matrix operations in expressions:" << std::endl;
    std::cout << "    matrix1 + matrix2                      # Matrix addition" << std::endl;
    std::cout << "    matrix1 - matrix2                      # Matrix subtraction" << std::endl;
    std::cout << "    matrix1 * matrix2                      # Matrix multiplication" << std::endl;
    std::cout << "    matrix * scalar                        # Scalar multiplication" << std::endl;
    std::cout << "    matrix / scalar                        # Scalar division" << std::endl;
    std::cout << "    matrix.transpose()                     # Matrix transpose" << std::endl;
    std::cout << "    matrix.inverse()                       # Matrix inverse" << std::endl;
    std::cout << "    matrix.determinant()                   # Matrix determinant" << std::endl;
    std::cout << "    matrix.trace()                         # Matrix trace" << std::endl;
    std::cout << "    matrix.rank()                          # Matrix rank" << std::endl;
    
    std::cout << std::endl;
    std::cout << m_themeManager->formatText("Examples:", "info") << std::endl;
    std::cout << "---------" << std::endl;
    std::cout << "  > 2 + 2                   # Basic arithmetic" << std::endl;
    std::cout << "  > sin(pi/2)               # Function with constant" << std::endl;
    std::cout << "  > x = 5                   # Variable assignment" << std::endl;
    std::cout << "  > y = 2*x + 3             # Assignment with expression" << std::endl;
    std::cout << "  > x += 3                  # Compound assignment" << std::endl;
    std::cout << "  > (3 + 4) * 2             # Parentheses for grouping" << std::endl;
    std::cout << "  > sqrt(x^2 + y^2)         # Pythagorean theorem" << std::endl;
    std::cout << "  > /vars                   # Display all variables" << std::endl;
    std::cout << "  > /solve x^2 - 4 = 0 x    # Solve equation for x" << std::endl;
    std::cout << "  > /diff x^2 x             # Differentiate x^2 with respect to x" << std::endl;
    std::cout << "  > /matrix identity 3      # Create a 3x3 identity matrix" << std::endl;
    std::cout << "  > A = /matrix create 2 2 1 # Create a 2x2 matrix filled with 1" << std::endl;
    std::cout << "  > B = /matrix random 2 2 0 10 # Create a 2x2 random matrix with values between 0 and 10" << std::endl;
    std::cout << "  > A * B                   # Multiply matrices A and B" << std::endl;
    std::cout << "  > A.transpose()           # Transpose matrix A" << std::endl;
    std::cout << "  > A.determinant()         # Calculate the determinant of matrix A" << std::endl;
    
    std::cout << std::endl;
    std::cout << m_themeManager->formatText("UI Features:", "info") << std::endl;
    std::cout << "-----------" << std::endl;
    std::cout << "  /syntax [on|off]          # Enable or disable syntax highlighting" << std::endl;
    std::cout << "  /autocomplete [on|off]    # Enable or disable autocompletion" << std::endl;
    std::cout << "  /multiline [on|off]       # Enable or disable multi-line input" << std::endl;
    std::cout << "  /theme [name]             # Set or display the current theme" << std::endl;
    std::cout << "  /workspace [create|switch] <name> # Manage workspaces" << std::endl;
    std::cout << "  /split [horizontal|vertical] # Split the current view" << std::endl;
    std::cout << "  /close                    # Close the current split view" << std::endl;
}

bool TerminalUI::Impl::processInput(const std::string& input) {
    // Check if the input is a command
    if (input[0] == '/') {
        return processCommand(input.substr(1));
    }
    
    // Check if the input is empty
    if (input.empty()) {
        return true;
    }
    
    // Check if the input is a help request
    if (input == "help" || input == "?") {
        displayHelp();
        return true;
    }
    
    // Otherwise, treat it as an expression to evaluate
    auto result = m_calculator->evaluate(input);
    if (!result) {
        std::cout << m_themeManager->formatText("Error: Failed to evaluate expression", "error") << std::endl;
        std::cout << "Type '/help' for syntax help or '?' for a quick reference" << std::endl;
        return false;
    }
    
    // Check if the input was a variable assignment
    std::regex assignmentRegex(R"(^\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*([+\-*/]?=).+$)");
    std::smatch matches;
    bool isAssignment = std::regex_match(input, matches, assignmentRegex);
    
    // Display the result
    std::visit([this, isAssignment, &matches](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        
        if (isAssignment) {
            std::string variableName = matches[1].str();
            std::string assignmentOp = matches[2].str();
            
            if constexpr (std::is_same_v<T, double>) {
                std::cout << variableName << " " << assignmentOp << " " 
                          << m_themeManager->formatText(std::to_string(arg), "result") << std::endl;
            } else {
                std::cout << variableName << " " << assignmentOp << " ";
                
                if constexpr (std::is_same_v<T, std::string>) {
                    std::cout << m_themeManager->formatText(arg, "result") << std::endl;
                } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                    std::string result = "[";
                    for (size_t i = 0; i < arg.size(); ++i) {
                        result += std::to_string(arg[i]);
                        if (i < arg.size() - 1) {
                            result += ", ";
                        }
                    }
                    result += "]";
                    std::cout << m_themeManager->formatText(result, "result") << std::endl;
                } else if constexpr (std::is_same_v<T, Matrix>) {
                    std::cout << std::endl << arg << std::endl;
                }
            }
        } else {
            if constexpr (std::is_same_v<T, double>) {
                std::cout << m_themeManager->formatText(std::to_string(arg), "result") << std::endl;
            } else if constexpr (std::is_same_v<T, std::string>) {
                std::cout << m_themeManager->formatText(arg, "result") << std::endl;
            } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                std::string result = "[";
                for (size_t i = 0; i < arg.size(); ++i) {
                    result += std::to_string(arg[i]);
                    if (i < arg.size() - 1) {
                        result += ", ";
                    }
                }
                result += "]";
                std::cout << m_themeManager->formatText(result, "result") << std::endl;
            } else if constexpr (std::is_same_v<T, Matrix>) {
                std::cout << std::endl << arg << std::endl;
            }
        }
    }, *result);
    
    return true;
}

bool TerminalUI::Impl::processCommand(const std::string& commandLine) {
    // Parse the command and arguments
    std::istringstream iss(commandLine);
    std::string command;
    iss >> command;
    
    // Check if the command exists
    auto it = m_commands.find(command);
    if (it == m_commands.end()) {
        std::cout << "Unknown command: " << command << std::endl;
        std::cout << "Type '/help' for a list of available commands" << std::endl;
        return false;
    }
    
    // Parse arguments
    std::vector<std::string> args;
    std::string arg;
    while (iss >> arg) {
        args.push_back(arg);
    }
    
    // Execute the command
    return it->second.handler(args);
}

} // namespace rebelcalc
