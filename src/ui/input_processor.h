#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <memory>
#include <optional>

namespace rebelcalc {

/**
 * Class for processing user input in the terminal UI
 * Provides features like command history, autocompletion, and syntax highlighting
 */
class InputProcessor {
public:
    /**
     * Constructor
     */
    InputProcessor();
    
    /**
     * Destructor
     */
    ~InputProcessor();
    
    /**
     * Initialize the input processor
     * @return true if initialization was successful, false otherwise
     */
    bool initialize();
    
    /**
     * Shutdown the input processor
     */
    void shutdown();
    
    /**
     * Read a line of input from the user with enhanced features
     * @param prompt The prompt to display
     * @return The input string
     */
    std::string readLine(const std::string& prompt);
    
    /**
     * Add a command to the history
     * @param command The command to add
     */
    void addToHistory(const std::string& command);
    
    /**
     * Clear the command history
     */
    void clearHistory();
    
    /**
     * Save the command history to a file
     * @param filename The name of the file to save to
     * @return true if the history was saved successfully, false otherwise
     */
    bool saveHistory(const std::string& filename);
    
    /**
     * Load the command history from a file
     * @param filename The name of the file to load from
     * @return true if the history was loaded successfully, false otherwise
     */
    bool loadHistory(const std::string& filename);
    
    /**
     * Get the command history
     * @return The command history
     */
    const std::vector<std::string>& getHistory() const;
    
    /**
     * Set the maximum history size
     * @param size The maximum number of commands to store in history
     */
    void setMaxHistorySize(size_t size);
    
    /**
     * Get the maximum history size
     * @return The maximum number of commands stored in history
     */
    size_t getMaxHistorySize() const;
    
    /**
     * Add a completion provider
     * @param name The name of the provider
     * @param provider The completion provider function
     */
    void addCompletionProvider(const std::string& name, 
                              std::function<std::vector<std::string>(const std::string&)> provider);
    
    /**
     * Remove a completion provider
     * @param name The name of the provider to remove
     * @return true if the provider was removed successfully, false otherwise
     */
    bool removeCompletionProvider(const std::string& name);
    
    /**
     * Get completions for a partial input
     * @param input The partial input to complete
     * @return A vector of possible completions
     */
    std::vector<std::string> getCompletions(const std::string& input);
    
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
     * Set the indentation string for multi-line input
     * @param indentation The indentation string (e.g., "  " for two spaces)
     */
    void setIndentation(const std::string& indentation);
    
    /**
     * Get the indentation string for multi-line input
     * @return The indentation string
     */
    const std::string& getIndentation() const;
    
    /**
     * Add a syntax highlighting rule
     * @param pattern The regex pattern to match
     * @param colorType The type of color to use (foreground, background, prompt, result, error, warning, info)
     */
    void addSyntaxHighlightingRule(const std::string& pattern, const std::string& colorType);
    
    /**
     * Clear all syntax highlighting rules
     */
    void clearSyntaxHighlightingRules();
    
    /**
     * Apply syntax highlighting to a string
     * @param text The text to highlight
     * @return The highlighted text
     */
    std::string highlightSyntax(const std::string& text);

private:
    // Implementation details
    class Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace rebelcalc
