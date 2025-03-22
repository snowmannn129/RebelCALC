#include "input_processor.h"
#include "theme_manager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>
#include <filesystem>

#ifdef _WIN32
#include <conio.h>
#include <windows.h>
#else
#include <termios.h>
#include <unistd.h>
#include <sys/ioctl.h>
#endif

namespace rebelcalc {

// Key codes
enum KeyCode {
    KEY_BACKSPACE = 8,
    KEY_TAB = 9,
    KEY_ENTER = 13,
    KEY_ESCAPE = 27,
    KEY_SPACE = 32,
    KEY_DELETE = 127,
    KEY_UP = 1001,
    KEY_DOWN = 1002,
    KEY_RIGHT = 1003,
    KEY_LEFT = 1004,
    KEY_HOME = 1005,
    KEY_END = 1006,
    KEY_PAGEUP = 1007,
    KEY_PAGEDOWN = 1008,
    KEY_INSERT = 1009,
    KEY_F1 = 1010,
    KEY_F2 = 1011,
    KEY_F3 = 1012,
    KEY_F4 = 1013,
    KEY_F5 = 1014,
    KEY_F6 = 1015,
    KEY_F7 = 1016,
    KEY_F8 = 1017,
    KEY_F9 = 1018,
    KEY_F10 = 1019,
    KEY_F11 = 1020,
    KEY_F12 = 1021,
    KEY_CTRL_A = 1,
    KEY_CTRL_B = 2,
    KEY_CTRL_C = 3,
    KEY_CTRL_D = 4,
    KEY_CTRL_E = 5,
    KEY_CTRL_F = 6,
    KEY_CTRL_G = 7,
    KEY_CTRL_H = 8,
    KEY_CTRL_I = 9,
    KEY_CTRL_J = 10,
    KEY_CTRL_K = 11,
    KEY_CTRL_L = 12,
    KEY_CTRL_M = 13,
    KEY_CTRL_N = 14,
    KEY_CTRL_O = 15,
    KEY_CTRL_P = 16,
    KEY_CTRL_Q = 17,
    KEY_CTRL_R = 18,
    KEY_CTRL_S = 19,
    KEY_CTRL_T = 20,
    KEY_CTRL_U = 21,
    KEY_CTRL_V = 22,
    KEY_CTRL_W = 23,
    KEY_CTRL_X = 24,
    KEY_CTRL_Y = 25,
    KEY_CTRL_Z = 26
};

// Terminal utilities
namespace {

#ifdef _WIN32
// Windows implementation
HANDLE hStdin = GetStdHandle(STD_INPUT_HANDLE);
HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
DWORD originalConsoleMode;

void initTerminal() {
    GetConsoleMode(hStdin, &originalConsoleMode);
    DWORD mode = originalConsoleMode;
    mode &= ~(ENABLE_ECHO_INPUT | ENABLE_LINE_INPUT);
    mode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    SetConsoleMode(hStdin, mode);
}

void resetTerminal() {
    SetConsoleMode(hStdin, originalConsoleMode);
}

int getTerminalWidth() {
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(hStdout, &csbi);
    return csbi.srWindow.Right - csbi.srWindow.Left + 1;
}

int getTerminalHeight() {
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(hStdout, &csbi);
    return csbi.srWindow.Bottom - csbi.srWindow.Top + 1;
}

int readKey() {
    INPUT_RECORD record;
    DWORD numRead;
    
    while (true) {
        ReadConsoleInput(hStdin, &record, 1, &numRead);
        
        if (record.EventType == KEY_EVENT && record.Event.KeyEvent.bKeyDown) {
            WORD vk = record.Event.KeyEvent.wVirtualKeyCode;
            CHAR ch = record.Event.KeyEvent.uChar.AsciiChar;
            DWORD controlKeyState = record.Event.KeyEvent.dwControlKeyState;
            
            if (ch == 0) {
                // Special key
                switch (vk) {
                    case VK_UP: return KEY_UP;
                    case VK_DOWN: return KEY_DOWN;
                    case VK_RIGHT: return KEY_RIGHT;
                    case VK_LEFT: return KEY_LEFT;
                    case VK_HOME: return KEY_HOME;
                    case VK_END: return KEY_END;
                    case VK_PRIOR: return KEY_PAGEUP;
                    case VK_NEXT: return KEY_PAGEDOWN;
                    case VK_INSERT: return KEY_INSERT;
                    case VK_DELETE: return KEY_DELETE;
                    case VK_F1: return KEY_F1;
                    case VK_F2: return KEY_F2;
                    case VK_F3: return KEY_F3;
                    case VK_F4: return KEY_F4;
                    case VK_F5: return KEY_F5;
                    case VK_F6: return KEY_F6;
                    case VK_F7: return KEY_F7;
                    case VK_F8: return KEY_F8;
                    case VK_F9: return KEY_F9;
                    case VK_F10: return KEY_F10;
                    case VK_F11: return KEY_F11;
                    case VK_F12: return KEY_F12;
                }
            } else {
                // Regular key
                if (controlKeyState & (LEFT_CTRL_PRESSED | RIGHT_CTRL_PRESSED)) {
                    // Control key combination
                    if (ch >= 'a' && ch <= 'z') {
                        return ch - 'a' + 1; // Ctrl+A to Ctrl+Z
                    } else if (ch >= 'A' && ch <= 'Z') {
                        return ch - 'A' + 1; // Ctrl+A to Ctrl+Z
                    }
                }
                
                return ch;
            }
        }
    }
}

#else
// Unix implementation
struct termios originalTermios;

void initTerminal() {
    tcgetattr(STDIN_FILENO, &originalTermios);
    struct termios raw = originalTermios;
    raw.c_lflag &= ~(ECHO | ICANON);
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);
}

void resetTerminal() {
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &originalTermios);
}

int getTerminalWidth() {
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col;
}

int getTerminalHeight() {
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_row;
}

int readKey() {
    int nread;
    char c;
    char seq[3];
    
    while ((nread = read(STDIN_FILENO, &c, 1)) == 0);
    
    if (nread == -1) {
        return -1;
    }
    
    if (c == 27) {
        // Escape sequence
        if (read(STDIN_FILENO, &seq[0], 1) == 0) return KEY_ESCAPE;
        if (read(STDIN_FILENO, &seq[1], 1) == 0) return KEY_ESCAPE;
        
        if (seq[0] == '[') {
            if (seq[1] >= '0' && seq[1] <= '9') {
                if (read(STDIN_FILENO, &seq[2], 1) == 0) return KEY_ESCAPE;
                if (seq[2] == '~') {
                    switch (seq[1]) {
                        case '1': return KEY_HOME;
                        case '2': return KEY_INSERT;
                        case '3': return KEY_DELETE;
                        case '4': return KEY_END;
                        case '5': return KEY_PAGEUP;
                        case '6': return KEY_PAGEDOWN;
                        case '7': return KEY_HOME;
                        case '8': return KEY_END;
                    }
                }
            } else {
                switch (seq[1]) {
                    case 'A': return KEY_UP;
                    case 'B': return KEY_DOWN;
                    case 'C': return KEY_RIGHT;
                    case 'D': return KEY_LEFT;
                    case 'H': return KEY_HOME;
                    case 'F': return KEY_END;
                }
            }
        } else if (seq[0] == 'O') {
            switch (seq[1]) {
                case 'A': return KEY_UP;
                case 'B': return KEY_DOWN;
                case 'C': return KEY_RIGHT;
                case 'D': return KEY_LEFT;
                case 'H': return KEY_HOME;
                case 'F': return KEY_END;
                case 'P': return KEY_F1;
                case 'Q': return KEY_F2;
                case 'R': return KEY_F3;
                case 'S': return KEY_F4;
            }
        }
        
        return KEY_ESCAPE;
    } else if (c == 127) {
        return KEY_BACKSPACE;
    } else if (c == 13) {
        return KEY_ENTER;
    } else if (c == 9) {
        return KEY_TAB;
    } else if (c >= 1 && c <= 26) {
        return c; // Ctrl+A to Ctrl+Z
    }
    
    return c;
}
#endif

// Common utilities
void clearLine() {
    std::cout << "\r\033[K";
}

void moveCursorToColumn(int column) {
    std::cout << "\r\033[" << column << "C";
}

void moveCursorUp(int lines) {
    std::cout << "\033[" << lines << "A";
}

void moveCursorDown(int lines) {
    std::cout << "\033[" << lines << "B";
}

} // anonymous namespace

// InputProcessor implementation
class InputProcessor::Impl {
public:
    Impl()
        : m_maxHistorySize(1000),
          m_syntaxHighlightingEnabled(true),
          m_autocompletionEnabled(true),
          m_multiLineInputEnabled(false),
          m_indentation("  ") {
        
        // Add default syntax highlighting rules
        addSyntaxHighlightingRule("\\b(sin|cos|tan|sqrt|log|exp)\\b", "info");
        addSyntaxHighlightingRule("\\b(pi|e)\\b", "info");
        addSyntaxHighlightingRule("\\b(if|else|for|while|function|return)\\b", "prompt");
        addSyntaxHighlightingRule("\\b(true|false|null)\\b", "warning");
        addSyntaxHighlightingRule("\\b[0-9]+(\\.[0-9]+)?\\b", "result");
        addSyntaxHighlightingRule("\\b[A-Za-z_][A-Za-z0-9_]*\\s*=", "prompt");
        addSyntaxHighlightingRule("\"[^\"]*\"", "warning");
        addSyntaxHighlightingRule("'[^']*'", "warning");
        addSyntaxHighlightingRule("//.*$", "info");
        addSyntaxHighlightingRule("/\\*.*?\\*/", "info");
    }
    
    ~Impl() {
        shutdown();
    }
    
    bool initialize() {
        try {
            initTerminal();
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Exception during input processor initialization: " << e.what() << std::endl;
            return false;
        }
    }
    
    void shutdown() {
        resetTerminal();
    }
    
    std::string readLine(const std::string& prompt) {
        std::string input;
        std::vector<std::string> lines = {""};
        size_t currentLine = 0;
        size_t cursorPos = 0;
        size_t historyPos = m_history.size();
        std::vector<std::string> completions;
        size_t completionIndex = 0;
        bool showCompletions = false;
        
        // Display prompt
        std::cout << prompt;
        
        while (true) {
            // Display current line with syntax highlighting
            clearLine();
            std::cout << prompt;
            
            if (m_syntaxHighlightingEnabled) {
                std::cout << highlightSyntax(lines[currentLine]);
            } else {
                std::cout << lines[currentLine];
            }
            
            // Move cursor to the correct position
            moveCursorToColumn(prompt.length() + cursorPos);
            
            // Display completions if needed
            if (showCompletions && !completions.empty()) {
                std::cout << std::endl;
                
                int termWidth = getTerminalWidth();
                int maxCompletionLength = 0;
                for (const auto& completion : completions) {
                    maxCompletionLength = std::max(maxCompletionLength, static_cast<int>(completion.length()));
                }
                
                int completionsPerRow = std::max(1, termWidth / (maxCompletionLength + 2));
                int rows = (completions.size() + completionsPerRow - 1) / completionsPerRow;
                
                for (int row = 0; row < rows; ++row) {
                    for (int col = 0; col < completionsPerRow; ++col) {
                        int index = row + col * rows;
                        if (index < static_cast<int>(completions.size())) {
                            std::string completion = completions[index];
                            if (index == static_cast<int>(completionIndex)) {
                                // Highlight the selected completion
                                std::cout << "\033[7m" << completion << "\033[0m";
                            } else {
                                std::cout << completion;
                            }
                            
                            // Pad with spaces
                            int padding = maxCompletionLength - static_cast<int>(completion.length()) + 2;
                            for (int i = 0; i < padding; ++i) {
                                std::cout << " ";
                            }
                        }
                    }
                    
                    if (row < rows - 1) {
                        std::cout << std::endl;
                    }
                }
                
                // Move cursor back to the input line
                moveCursorUp(rows);
                moveCursorToColumn(prompt.length() + cursorPos);
            }
            
            // Read a key
            int key = readKey();
            
            if (key == KEY_ENTER) {
                // Enter key
                std::cout << std::endl;
                
                if (m_multiLineInputEnabled) {
                    // Check if this is the end of a multi-line input
                    if (lines[currentLine].empty() && currentLine > 0) {
                        // Empty line, end of multi-line input
                        break;
                    } else {
                        // Add a new line
                        lines.push_back("");
                        currentLine++;
                        cursorPos = 0;
                    }
                } else {
                    // Single-line input, we're done
                    break;
                }
            } else if (key == KEY_BACKSPACE) {
                // Backspace key
                if (cursorPos > 0) {
                    lines[currentLine].erase(cursorPos - 1, 1);
                    cursorPos--;
                } else if (currentLine > 0) {
                    // At the beginning of a line, merge with the previous line
                    cursorPos = lines[currentLine - 1].length();
                    lines[currentLine - 1] += lines[currentLine];
                    lines.erase(lines.begin() + currentLine);
                    currentLine--;
                }
            } else if (key == KEY_DELETE) {
                // Delete key
                if (cursorPos < lines[currentLine].length()) {
                    lines[currentLine].erase(cursorPos, 1);
                } else if (currentLine < lines.size() - 1) {
                    // At the end of a line, merge with the next line
                    lines[currentLine] += lines[currentLine + 1];
                    lines.erase(lines.begin() + currentLine + 1);
                }
            } else if (key == KEY_LEFT) {
                // Left arrow key
                if (cursorPos > 0) {
                    cursorPos--;
                } else if (currentLine > 0) {
                    // Move to the end of the previous line
                    currentLine--;
                    cursorPos = lines[currentLine].length();
                }
            } else if (key == KEY_RIGHT) {
                // Right arrow key
                if (cursorPos < lines[currentLine].length()) {
                    cursorPos++;
                } else if (currentLine < lines.size() - 1) {
                    // Move to the beginning of the next line
                    currentLine++;
                    cursorPos = 0;
                }
            } else if (key == KEY_UP) {
                // Up arrow key
                if (m_multiLineInputEnabled && currentLine > 0) {
                    // Move to the previous line
                    currentLine--;
                    cursorPos = std::min(cursorPos, lines[currentLine].length());
                } else if (historyPos > 0) {
                    // Navigate history
                    historyPos--;
                    lines[currentLine] = m_history[historyPos];
                    cursorPos = lines[currentLine].length();
                }
            } else if (key == KEY_DOWN) {
                // Down arrow key
                if (m_multiLineInputEnabled && currentLine < lines.size() - 1) {
                    // Move to the next line
                    currentLine++;
                    cursorPos = std::min(cursorPos, lines[currentLine].length());
                } else if (historyPos < m_history.size()) {
                    // Navigate history
                    historyPos++;
                    if (historyPos == m_history.size()) {
                        lines[currentLine] = "";
                    } else {
                        lines[currentLine] = m_history[historyPos];
                    }
                    cursorPos = lines[currentLine].length();
                }
            } else if (key == KEY_HOME) {
                // Home key
                cursorPos = 0;
            } else if (key == KEY_END) {
                // End key
                cursorPos = lines[currentLine].length();
            } else if (key == KEY_TAB) {
                // Tab key
                if (m_autocompletionEnabled) {
                    if (!showCompletions) {
                        // Get completions
                        completions = getCompletions(lines[currentLine].substr(0, cursorPos));
                        completionIndex = 0;
                        showCompletions = !completions.empty();
                    } else {
                        // Cycle through completions
                        completionIndex = (completionIndex + 1) % completions.size();
                    }
                    
                    if (!completions.empty()) {
                        // Apply the selected completion
                        std::string completion = completions[completionIndex];
                        
                        // Find the start of the word to complete
                        size_t wordStart = cursorPos;
                        while (wordStart > 0 && (std::isalnum(lines[currentLine][wordStart - 1]) || lines[currentLine][wordStart - 1] == '_')) {
                            wordStart--;
                        }
                        
                        // Replace the word with the completion
                        lines[currentLine].replace(wordStart, cursorPos - wordStart, completion);
                        cursorPos = wordStart + completion.length();
                    }
                } else if (m_multiLineInputEnabled) {
                    // Insert indentation
                    lines[currentLine].insert(cursorPos, m_indentation);
                    cursorPos += m_indentation.length();
                } else {
                    // Insert a tab character
                    lines[currentLine].insert(cursorPos, "\t");
                    cursorPos++;
                }
            } else if (key == KEY_CTRL_C) {
                // Ctrl+C, cancel input
                std::cout << "^C" << std::endl;
                return "";
            } else if (key == KEY_CTRL_L) {
                // Ctrl+L, clear screen
                std::cout << "\033[2J\033[H";
                std::cout << prompt << lines[currentLine];
            } else if (key == KEY_ESCAPE) {
                // Escape key, hide completions
                showCompletions = false;
            } else if (key >= 32 && key <= 126) {
                // Printable character
                lines[currentLine].insert(cursorPos, 1, static_cast<char>(key));
                cursorPos++;
                
                // Hide completions when typing
                showCompletions = false;
            }
        }
        
        // Combine lines into a single string
        std::string result;
        for (size_t i = 0; i < lines.size(); ++i) {
            result += lines[i];
            if (i < lines.size() - 1) {
                result += "\n";
            }
        }
        
        // Add to history if not empty
        if (!result.empty()) {
            addToHistory(result);
        }
        
        return result;
    }
    
    void addToHistory(const std::string& command) {
        // Don't add duplicates
        if (!m_history.empty() && m_history.back() == command) {
            return;
        }
        
        m_history.push_back(command);
        
        // Trim history if it exceeds the maximum size
        if (m_history.size() > m_maxHistorySize) {
            m_history.erase(m_history.begin());
        }
    }
    
    void clearHistory() {
        m_history.clear();
    }
    
    bool saveHistory(const std::string& filename) {
        try {
            std::ofstream file(filename);
            if (!file.is_open()) {
                return false;
            }
            
            for (const auto& command : m_history) {
                file << command << std::endl;
            }
            
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Exception during history save: " << e.what() << std::endl;
            return false;
        }
    }
    
    bool loadHistory(const std::string& filename) {
        try {
            std::ifstream file(filename);
            if (!file.is_open()) {
                return false;
            }
            
            m_history.clear();
            std::string line;
            while (std::getline(file, line)) {
                m_history.push_back(line);
            }
            
            // Trim history if it exceeds the maximum size
            if (m_history.size() > m_maxHistorySize) {
                m_history.erase(m_history.begin(), m_history.begin() + (m_history.size() - m_maxHistorySize));
            }
            
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Exception during history load: " << e.what() << std::endl;
            return false;
        }
    }
    
    const std::vector<std::string>& getHistory() const {
        return m_history;
    }
    
    void setMaxHistorySize(size_t size) {
        m_maxHistorySize = size;
        
        // Trim history if it exceeds the new maximum size
        if (m_history.size() > m_maxHistorySize) {
            m_history.erase(m_history.begin(), m_history.begin() + (m_history.size() - m_maxHistorySize));
        }
    }
    
    size_t getMaxHistorySize() const {
        return m_maxHistorySize;
    }
    
    void addCompletionProvider(const std::string& name, 
                              std::function<std::vector<std::string>(const std::string&)> provider) {
        m_completionProviders[name] = provider;
    }
    
    bool removeCompletionProvider(const std::string& name) {
        auto it = m_completionProviders.find(name);
        if (it != m_completionProviders.end()) {
            m_completionProviders.erase(it);
            return true;
        }
        return false;
    }
    
    std::vector<std::string> getCompletions(const std::string& input) {
        std::vector<std::string> result;
        
        // Get completions from all providers
        for (const auto& pair : m_completionProviders) {
            auto completions = pair.second(input);
            result.insert(result.end(), completions.begin(), completions.end());
        }
        
        // Sort and remove duplicates
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        
        return result;
    }
    
    void enableSyntaxHighlighting(bool enable) {
        m_syntaxHighlightingEnabled = enable;
    }
    
    bool isSyntaxHighlightingEnabled() const {
        return m_syntaxHighlightingEnabled;
    }
    
    void enableAutocompletion(bool enable) {
        m_autocompletionEnabled = enable;
    }
    
    bool isAutocompletionEnabled() const {
        return m_autocompletionEnabled;
    }
    
    void enableMultiLineInput(bool enable) {
        m_multiLineInputEnabled = enable;
    }
    
    bool isMultiLineInputEnabled() const {
        return m_multiLineInputEnabled;
    }
    
    void setIndentation(const std::string& indentation) {
        m_indentation = indentation;
    }
    
    const std::string& getIndentation() const {
        return m_indentation;
    }
    
    void addSyntaxHighlightingRule(const std::string& pattern, const std::string& colorType) {
        m_syntaxHighlightingRules.push_back({pattern, colorType});
    }
    
    void clearSyntaxHighlightingRules() {
        m_syntaxHighlightingRules.clear();
    }
    
    std::string highlightSyntax(const std::string& text) {
        std::string result = text;
        
        // Apply syntax highlighting rules
        for (const auto& rule : m_syntaxHighlightingRules) {
            try {
                std::regex regex(rule.pattern);
                std::string replacement = "\033[" + getColorCode(rule.colorType) + "m$&\033[0m";
                result = std::regex_replace(result, regex, replacement);
            } catch (const std::regex_error& e) {
                std::cerr << "Regex error: " << e.what() << " in pattern: " << rule.pattern << std::endl;
            }
        }
        
        return result;
    }
    
private:
    std::vector<std::string> m_history;
    size_t m_maxHistorySize;
    bool m_syntaxHighlightingEnabled;
    bool m_autocompletionEnabled;
    bool m_multiLineInputEnabled;
    std::string m_indentation;
    
    struct SyntaxHighlightingRule {
        std::string pattern;
        std::string colorType;
    };
    
    std::vector<SyntaxHighlightingRule> m_syntaxHighlightingRules;
    std::unordered_map<std::string, std::function<std::vector<std::string>(const std::string&)>> m_completionProviders;
    
    std::string getColorCode(const std::string& colorType) {
        // ANSI color codes
        if (colorType == "foreground") return "39";
        if (colorType == "background") return "49";
        if (colorType == "prompt") return "36";
        if (colorType == "result") return "32";
        if (colorType == "error") return "31";
        if (colorType == "warning") return "33";
        if (colorType == "info") return "34";
        
        return "39"; // Default foreground color
    }
};

// InputProcessor implementation
InputProcessor::InputProcessor()
    : m_impl(std::make_unique<Impl>()) {
}

InputProcessor::~InputProcessor() {
}

bool InputProcessor::initialize() {
    return m_impl->initialize();
}

void InputProcessor::shutdown() {
    m_impl->shutdown();
}

std::string InputProcessor::readLine(const std::string& prompt) {
    return m_impl->readLine(prompt);
}

void InputProcessor::addToHistory(const std::string& command) {
    m_impl->addToHistory(command);
}

void InputProcessor::clearHistory() {
    m_impl->clearHistory();
}

bool InputProcessor::saveHistory(const std::string& filename) {
    return m_impl->saveHistory(filename);
}

bool InputProcessor::loadHistory(const std::string& filename) {
    return m_impl->loadHistory(filename);
}

const std::vector<std::string>& InputProcessor::getHistory() const {
    return m_impl->getHistory();
}

void InputProcessor::setMaxHistorySize(size_t size) {
    m_impl->setMaxHistorySize(size);
}

size_t InputProcessor::getMaxHistorySize() const {
    return m_impl->getMaxHistorySize();
}

void InputProcessor::addCompletionProvider(const std::string& name, 
                                         std::function<std::vector<std::string>(const std::string&)> provider) {
    m_impl->addCompletionProvider(name, provider);
}

bool InputProcessor::removeCompletionProvider(const std::string& name) {
    return m_impl->removeCompletionProvider(name);
}

std::vector<std::string> InputProcessor::getCompletions(const std::string& input) {
    return m_impl->getCompletions(input);
}

void InputProcessor::enableSyntaxHighlighting(bool enable) {
    m_impl->enableSyntaxHighlighting(enable);
}

bool InputProcessor::isSyntaxHighlightingEnabled() const {
    return m_impl->isSyntaxHighlightingEnabled();
}

void InputProcessor::enableAutocompletion(bool enable) {
    m_impl->enableAutocompletion(enable);
}

bool InputProcessor::isAutocompletionEnabled() const {
    return m_impl->isAutocompletionEnabled();
}

void InputProcessor::enableMultiLineInput(bool enable) {
    m_impl->enableMultiLineInput(enable);
}

bool InputProcessor::isMultiLineInputEnabled() const {
    return m_impl->isMultiLineInputEnabled();
}

void InputProcessor::setIndentation(const std::string& indentation) {
    m_impl->setIndentation(indentation);
}

const std::string& InputProcessor::getIndentation() const {
    return m_impl->getIndentation();
}

void InputProcessor::addSyntaxHighlightingRule(const std::string& pattern, const std::string& colorType) {
    m_impl->addSyntaxHighlightingRule(pattern, colorType);
}

void InputProcessor::clearSyntaxHighlightingRules() {
    m_impl->clearSyntaxHighlightingRules();
}

std::string InputProcessor::highlightSyntax(const std::string& text) {
    return m_impl->highlightSyntax(text);
}

} // namespace rebelcalc
