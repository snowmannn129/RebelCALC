# InputProcessor Usage Guide

The `InputProcessor` class provides an enhanced Read-Eval-Print Loop (REPL) interface for RebelCALC, with features like command history navigation, autocompletion, syntax highlighting, and multi-line input. This guide explains how to use the `InputProcessor` class in your applications.

## Features

- **Command History Navigation**: Use the up and down arrow keys to navigate through previously entered commands.
- **Autocompletion**: Press the Tab key to autocomplete commands, functions, and variables.
- **Syntax Highlighting**: Colorize different parts of the input based on customizable rules.
- **Multi-line Input**: Enter complex expressions or scripts that span multiple lines.
- **Cursor Movement and Editing**: Use arrow keys, Home, End, Backspace, and Delete for editing.
- **Customizable Appearance**: Define your own syntax highlighting rules and color schemes.

## Basic Usage

Here's a simple example of how to use the `InputProcessor` class:

```cpp
#include <iostream>
#include <string>
#include "ui/input_processor.h"

using namespace rebelcalc;

int main() {
    // Create an input processor
    InputProcessor inputProcessor;
    
    // Initialize the input processor
    if (!inputProcessor.initialize()) {
        std::cerr << "Failed to initialize input processor" << std::endl;
        return 1;
    }
    
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
        } else {
            std::cout << "You entered: " << input << std::endl;
        }
    }
    
    // Shutdown the input processor
    inputProcessor.shutdown();
    
    return 0;
}
```

## Advanced Features

### Command History

The `InputProcessor` class automatically maintains a history of commands entered by the user. You can also manually add commands to the history:

```cpp
// Add a command to the history
inputProcessor.addToHistory("2 + 2");

// Clear the history
inputProcessor.clearHistory();

// Save the history to a file
inputProcessor.saveHistory("history.txt");

// Load the history from a file
inputProcessor.loadHistory("history.txt");

// Get the history
const std::vector<std::string>& history = inputProcessor.getHistory();

// Set the maximum history size
inputProcessor.setMaxHistorySize(1000);

// Get the maximum history size
size_t maxHistorySize = inputProcessor.getMaxHistorySize();
```

### Autocompletion

You can add custom completion providers to suggest completions for different types of input:

```cpp
// Function to provide command completions
std::vector<std::string> commandCompletionProvider(const std::string& input) {
    std::vector<std::string> commands = {
        "help", "exit", "clear", "vars", "clear_vars", "solve", "diff", "int"
    };
    
    std::vector<std::string> completions;
    for (const auto& cmd : commands) {
        if (cmd.find(input) == 0) {
            completions.push_back(cmd);
        }
    }
    
    return completions;
}

// Add a completion provider
inputProcessor.addCompletionProvider("commands", commandCompletionProvider);

// Remove a completion provider
inputProcessor.removeCompletionProvider("commands");

// Get completions for a partial input
std::vector<std::string> completions = inputProcessor.getCompletions("he");
// Returns ["help"]

// Enable or disable autocompletion
inputProcessor.enableAutocompletion(true);

// Check if autocompletion is enabled
bool autocompletionEnabled = inputProcessor.isAutocompletionEnabled();
```

### Syntax Highlighting

You can add custom syntax highlighting rules to colorize different parts of the input:

```cpp
// Add a syntax highlighting rule
inputProcessor.addSyntaxHighlightingRule("\\b(sin|cos|tan)\\b", "info");

// Clear all syntax highlighting rules
inputProcessor.clearSyntaxHighlightingRules();

// Apply syntax highlighting to a string
std::string highlighted = inputProcessor.highlightSyntax("sin(x) + cos(y)");

// Enable or disable syntax highlighting
inputProcessor.enableSyntaxHighlighting(true);

// Check if syntax highlighting is enabled
bool syntaxHighlightingEnabled = inputProcessor.isSyntaxHighlightingEnabled();
```

The available color types are:
- `"foreground"`: Default foreground color
- `"background"`: Default background color
- `"prompt"`: Prompt color (cyan)
- `"result"`: Result color (green)
- `"error"`: Error color (red)
- `"warning"`: Warning color (yellow)
- `"info"`: Info color (blue)

### Multi-line Input

You can enable multi-line input for entering complex expressions or scripts:

```cpp
// Enable or disable multi-line input
inputProcessor.enableMultiLineInput(true);

// Check if multi-line input is enabled
bool multiLineInputEnabled = inputProcessor.isMultiLineInputEnabled();

// Set the indentation string for multi-line input
inputProcessor.setIndentation("  ");

// Get the indentation string for multi-line input
const std::string& indentation = inputProcessor.getIndentation();
```

When multi-line input is enabled, the user can press Enter to add a new line, and press Enter on an empty line to submit the input.

## Integration with Calculator

For a complete REPL experience, you can integrate the `InputProcessor` class with the `Calculator` class:

```cpp
#include <iostream>
#include <string>
#include "ui/input_processor.h"
#include "backend/calculator.h"

using namespace rebelcalc;

int main() {
    // Create a calculator
    auto calculator = std::make_shared<Calculator>();
    
    // Create an input processor
    InputProcessor inputProcessor;
    
    // Initialize the input processor
    if (!inputProcessor.initialize()) {
        std::cerr << "Failed to initialize input processor" << std::endl;
        return 1;
    }
    
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
```

## Example: REPL Demo

For a complete example of how to use the `InputProcessor` class, see the `repl_demo.cpp` file in the `examples` directory. This example demonstrates how to create a full-featured REPL interface with command history, autocompletion, syntax highlighting, and integration with the calculator.

## Keyboard Shortcuts

The `InputProcessor` class supports the following keyboard shortcuts:

- **Up Arrow**: Navigate to the previous command in history
- **Down Arrow**: Navigate to the next command in history
- **Left Arrow**: Move the cursor left
- **Right Arrow**: Move the cursor right
- **Home**: Move the cursor to the beginning of the line
- **End**: Move the cursor to the end of the line
- **Backspace**: Delete the character before the cursor
- **Delete**: Delete the character at the cursor
- **Tab**: Autocomplete the current input
- **Ctrl+C**: Cancel the current input
- **Ctrl+L**: Clear the screen
- **Escape**: Hide completions

## Platform Support

The `InputProcessor` class supports both Windows and Unix-like systems (Linux, macOS). It uses platform-specific code to handle terminal input and output, so it should work correctly on all supported platforms.

## Customization

The `InputProcessor` class is highly customizable. You can:

- Add your own completion providers for different types of input
- Define your own syntax highlighting rules for different parts of the input
- Set the maximum history size
- Enable or disable features like autocompletion, syntax highlighting, and multi-line input
- Set the indentation string for multi-line input

## Conclusion

The `InputProcessor` class provides a powerful and flexible REPL interface for RebelCALC. It enhances the user experience with features like command history, autocompletion, syntax highlighting, and multi-line input. By integrating it with the `Calculator` class, you can create a complete REPL environment for interactive calculations.
