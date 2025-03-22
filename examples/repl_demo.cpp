#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include "../src/ui/input_processor.h"
#include "../src/backend/calculator.h"
#include "../src/scripts/lua_engine.h"

using namespace rebelcalc;

// Function to provide command completions
std::vector<std::string> commandCompletionProvider(const std::string& input) {
    std::vector<std::string> commands = {
        "help", "exit", "clear", "vars", "clear_vars", "solve", "diff", "int", "lua", "load", "matrix"
    };
    
    std::vector<std::string> completions;
    for (const auto& cmd : commands) {
        if (cmd.find(input) == 0) {
            completions.push_back(cmd);
        }
    }
    
    return completions;
}

// Function to provide function completions
std::vector<std::string> functionCompletionProvider(const std::string& input) {
    std::vector<std::string> functions = {
        "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh",
        "log", "log10", "exp", "sqrt", "abs", "floor", "ceil", "round",
        "min", "max", "pow", "factorial", "gcd", "lcm"
    };
    
    std::vector<std::string> completions;
    for (const auto& func : functions) {
        if (func.find(input) == 0) {
            completions.push_back(func);
        }
    }
    
    return completions;
}

// Function to provide variable completions
std::vector<std::string> variableCompletionProvider(const std::string& input) {
    std::vector<std::string> variables = {
        "pi", "e", "x", "y", "z", "a", "b", "c", "result"
    };
    
    std::vector<std::string> completions;
    for (const auto& var : variables) {
        if (var.find(input) == 0) {
            completions.push_back(var);
        }
    }
    
    return completions;
}

int main() {
    std::cout << "RebelCALC Enhanced REPL Demo" << std::endl;
    std::cout << "===========================" << std::endl;
    std::cout << std::endl;
    
    // Create a calculator and Lua engine
    auto calculator = std::make_shared<Calculator>();
    auto luaEngine = std::make_shared<LuaEngine>(calculator);
    
    // Create an input processor
    InputProcessor inputProcessor;
    
    // Initialize the input processor
    if (!inputProcessor.initialize()) {
        std::cerr << "Failed to initialize input processor" << std::endl;
        return 1;
    }
    
    // Add completion providers
    inputProcessor.addCompletionProvider("commands", commandCompletionProvider);
    inputProcessor.addCompletionProvider("functions", functionCompletionProvider);
    inputProcessor.addCompletionProvider("variables", variableCompletionProvider);
    
    // Enable features
    inputProcessor.enableSyntaxHighlighting(true);
    inputProcessor.enableAutocompletion(true);
    inputProcessor.enableMultiLineInput(true);
    
    // Add custom syntax highlighting rules
    inputProcessor.addSyntaxHighlightingRule("\\b(help|exit|clear|vars|clear_vars|solve|diff|int|lua|load|matrix)\\b", "prompt");
    
    // Main loop
    bool running = true;
    while (running) {
        // Read a line of input
        std::string input = inputProcessor.readLine("> ");
        
        // Process the input
        if (input.empty()) {
            continue;
        }
        
        // Add to history
        inputProcessor.addToHistory(input);
        
        // Process commands
        if (input == "exit") {
            running = false;
        } else if (input == "help") {
            std::cout << "Available commands:" << std::endl;
            std::cout << "  help        - Display this help message" << std::endl;
            std::cout << "  exit        - Exit the program" << std::endl;
            std::cout << "  clear       - Clear the screen" << std::endl;
            std::cout << "  vars        - Display all variables" << std::endl;
            std::cout << "  clear_vars  - Clear all variables" << std::endl;
            std::cout << "  solve       - Solve an equation" << std::endl;
            std::cout << "  diff        - Differentiate an expression" << std::endl;
            std::cout << "  int         - Integrate an expression" << std::endl;
            std::cout << "  lua         - Execute a Lua script" << std::endl;
            std::cout << "  load        - Load a Lua script from a file" << std::endl;
            std::cout << "  matrix      - Matrix operations" << std::endl;
        } else if (input == "clear") {
            // Clear the screen
            std::cout << "\033[2J\033[H";
        } else {
            // Evaluate the expression
            try {
                auto result = calculator->evaluate(input);
                if (result) {
                    // Handle the variant type
                    std::visit([](auto&& arg) {
                        using T = std::decay_t<decltype(arg)>;
                        if constexpr (std::is_same_v<T, double>) {
                            std::cout << arg << std::endl;
                        } else if constexpr (std::is_same_v<T, std::string>) {
                            std::cout << arg << std::endl;
                        } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                            std::cout << "[";
                            for (size_t i = 0; i < arg.size(); ++i) {
                                std::cout << arg[i];
                                if (i < arg.size() - 1) {
                                    std::cout << ", ";
                                }
                            }
                            std::cout << "]" << std::endl;
                        } else if constexpr (std::is_same_v<T, Matrix>) {
                            std::cout << "Matrix result:" << std::endl;
                            // Display the matrix in a formatted way
                            for (size_t i = 0; i < arg.rows(); ++i) {
                                std::cout << "  [";
                                for (size_t j = 0; j < arg.cols(); ++j) {
                                    std::cout << arg(i, j);
                                    if (j < arg.cols() - 1) {
                                        std::cout << ", ";
                                    }
                                }
                                std::cout << "]" << std::endl;
                            }
                        } else {
                            std::cout << "Result of unknown type" << std::endl;
                        }
                    }, *result);
                } else {
                    std::cout << "Error: Failed to evaluate expression" << std::endl;
                }
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
        }
    }
    
    // Shutdown the input processor
    inputProcessor.shutdown();
    
    return 0;
}
