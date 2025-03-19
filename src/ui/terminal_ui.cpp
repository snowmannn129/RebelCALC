#include "terminal_ui.h"
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

namespace rebelcalc {

// Private implementation class
class TerminalUI::Impl {
public:
    Impl(std::shared_ptr<Calculator> calculator, std::shared_ptr<LuaEngine> luaEngine)
        : m_calculator(calculator), m_luaEngine(luaEngine), m_running(false), m_theme("default") {
        
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
    }
    
    ~Impl() {
    }
    
    bool initialize() {
        // Initialize the terminal UI implementation
        return true;
    }
    
    void shutdown() {
        // Shutdown the terminal UI implementation
    }
    
    void run() {
        m_running = true;
        
        while (m_running) {
            // Display prompt
            std::cout << "> ";
            
            // Read input
            std::string input;
            std::getline(std::cin, input);
            
            // Process input
            if (!input.empty()) {
                processInput(input);
            }
        }
    }
    
    bool setTheme(const std::string& themeName) {
        // This is a placeholder for a real theme system
        // In a real implementation, this would load theme settings from a file
        
        if (themeName == "default" || themeName == "dark" || themeName == "light") {
            m_theme = themeName;
            std::cout << "Theme set to " << themeName << std::endl;
            return true;
        }
        
        std::cout << "Unknown theme: " << themeName << std::endl;
        return false;
    }
    
    std::string getTheme() const {
        return m_theme;
    }
    
    void registerCommand(const std::string& name, 
                        std::function<bool(const std::vector<std::string>&)> handler,
                        const std::string& helpText) {
        m_commands[name] = {handler, helpText};
    }
    
    bool processInput(const std::string& input) {
        // Check if the input is a command
        if (input[0] == '/') {
            return processCommand(input.substr(1));
        }
        
        // Otherwise, treat it as an expression to evaluate
        auto result = m_calculator->evaluate(input);
        if (!result) {
            std::cout << "Error evaluating expression" << std::endl;
            return false;
        }
        
        std::visit([](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, double>) {
                std::cout << "Result: " << arg << std::endl;
            } else if constexpr (std::is_same_v<T, std::string>) {
                std::cout << "Result: " << arg << std::endl;
            } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                std::cout << "Result: ";
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
    }
    
    bool processCommand(const std::string& commandLine) {
        // Parse the command and arguments
        std::istringstream iss(commandLine);
        std::string command;
        iss >> command;
        
        std::vector<std::string> args;
        std::string arg;
        while (iss >> arg) {
            args.push_back(arg);
        }
        
        // Look up the command
        auto it = m_commands.find(command);
        if (it == m_commands.end()) {
            std::cout << "Unknown command: " << command << std::endl;
            return false;
        }
        
        // Execute the command
        return it->second.handler(args);
    }
    
    void displayHelp() const {
        std::cout << "Available commands:" << std::endl;
        
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
        std::cout << "For expressions, just type them directly:" << std::endl;
        std::cout << "  > 2 + 2" << std::endl;
        std::cout << "  > sin(pi/2)" << std::endl;
        std::cout << "  > x^2 + 2*x + 1" << std::endl;
    }
    
    void displayVariables() const {
        std::cout << "Variables:" << std::endl;
        
        // Built-in constants
        std::cout << "  pi = " << *m_calculator->getVariable("pi") << " (built-in)" << std::endl;
        std::cout << "  e = " << *m_calculator->getVariable("e") << " (built-in)" << std::endl;
        
        // User-defined variables
        // This is a placeholder - in a real implementation, we'd iterate through all variables
        // For now, just show a message
        std::cout << "  (No user-defined variables)" << std::endl;
    }

private:
    std::shared_ptr<Calculator> m_calculator;
    std::shared_ptr<LuaEngine> m_luaEngine;
    bool m_running;
    std::string m_theme;
    
    struct CommandInfo {
        std::function<bool(const std::vector<std::string>&)> handler;
        std::string helpText;
    };
    
    std::unordered_map<std::string, CommandInfo> m_commands;
};

// TerminalUI implementation

TerminalUI::TerminalUI(std::shared_ptr<Calculator> calculator, std::shared_ptr<LuaEngine> luaEngine)
    : m_calculator(calculator), m_luaEngine(luaEngine), m_impl(std::make_unique<Impl>(calculator, luaEngine)) {
}

TerminalUI::~TerminalUI() {
    shutdown();
}

bool TerminalUI::initialize() {
    try {
        return m_impl->initialize();
    } catch (const std::exception& e) {
        std::cerr << "Exception during terminal UI initialization: " << e.what() << std::endl;
        return false;
    }
}

void TerminalUI::shutdown() {
    if (m_impl) {
        m_impl->shutdown();
    }
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

void TerminalUI::registerCommand(const std::string& name, 
                               std::function<bool(const std::vector<std::string>&)> handler,
                               const std::string& helpText) {
    m_impl->registerCommand(name, handler, helpText);
}

bool TerminalUI::processCommand(const std::string& command) {
    return m_impl->processCommand(command);
}

void TerminalUI::displayHelp() const {
    m_impl->displayHelp();
}

void TerminalUI::displayResult(const std::string& result) const {
    std::cout << result << std::endl;
}

void TerminalUI::displayError(const std::string& error) const {
    std::cerr << "Error: " << error << std::endl;
}

} // namespace rebelcalc
