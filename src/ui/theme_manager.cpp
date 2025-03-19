#include "theme_manager.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <regex>
#include <yaml-cpp/yaml.h>

namespace rebelcalc {

// Color implementation

// Predefined colors
const Color Color::Black(0, 0, 0);
const Color Color::Red(255, 0, 0);
const Color Color::Green(0, 255, 0);
const Color Color::Yellow(255, 255, 0);
const Color Color::Blue(0, 0, 255);
const Color Color::Magenta(255, 0, 255);
const Color Color::Cyan(0, 255, 255);
const Color Color::White(255, 255, 255);
const Color Color::BrightBlack(128, 128, 128);
const Color Color::BrightRed(255, 128, 128);
const Color Color::BrightGreen(128, 255, 128);
const Color Color::BrightYellow(255, 255, 128);
const Color Color::BrightBlue(128, 128, 255);
const Color Color::BrightMagenta(255, 128, 255);
const Color Color::BrightCyan(128, 255, 255);
const Color Color::BrightWhite(255, 255, 255);

std::optional<Color> Color::fromHex(const std::string& hex) {
    // Remove leading # if present
    std::string hexValue = hex;
    if (hexValue.size() > 0 && hexValue[0] == '#') {
        hexValue = hexValue.substr(1);
    }
    
    // Check if the hex string is valid
    if (hexValue.size() != 6) {
        return std::nullopt;
    }
    
    // Check if all characters are valid hex digits
    for (char c : hexValue) {
        if (!std::isxdigit(c)) {
            return std::nullopt;
        }
    }
    
    // Parse the hex string
    uint8_t r = std::stoi(hexValue.substr(0, 2), nullptr, 16);
    uint8_t g = std::stoi(hexValue.substr(2, 2), nullptr, 16);
    uint8_t b = std::stoi(hexValue.substr(4, 2), nullptr, 16);
    
    return Color(r, g, b);
}

std::string Color::toHex() const {
    std::stringstream ss;
    ss << "#" << std::hex << std::setfill('0')
       << std::setw(2) << static_cast<int>(r)
       << std::setw(2) << static_cast<int>(g)
       << std::setw(2) << static_cast<int>(b);
    return ss.str();
}

std::string Color::toAnsi(bool foreground) const {
    // Convert RGB to ANSI escape sequence
    if (foreground) {
        return "\033[38;2;" + std::to_string(r) + ";" + std::to_string(g) + ";" + std::to_string(b) + "m";
    } else {
        return "\033[48;2;" + std::to_string(r) + ";" + std::to_string(g) + ";" + std::to_string(b) + "m";
    }
}

// Theme implementation

std::optional<Theme> Theme::fromConfig(const std::string& name,
                                     const std::unordered_map<std::string, std::string>& config) {
    // Check if all required colors are present
    std::vector<std::string> requiredColors = {
        "background", "foreground", "prompt", "result", "error", "warning", "info"
    };
    
    for (const auto& color : requiredColors) {
        if (config.find(color) == config.end()) {
            return std::nullopt;
        }
    }
    
    // Parse the colors
    auto background = Color::fromHex(config.at("background"));
    auto foreground = Color::fromHex(config.at("foreground"));
    auto prompt = Color::fromHex(config.at("prompt"));
    auto result = Color::fromHex(config.at("result"));
    auto error = Color::fromHex(config.at("error"));
    auto warning = Color::fromHex(config.at("warning"));
    auto info = Color::fromHex(config.at("info"));
    
    // Check if all colors were parsed successfully
    if (!background || !foreground || !prompt || !result || !error || !warning || !info) {
        return std::nullopt;
    }
    
    // Create the theme
    return Theme(name, *background, *foreground, *prompt, *result, *error, *warning, *info);
}

// ThemeManager implementation

class ThemeManager::Impl {
public:
    Impl() : m_currentTheme("default") {}
    
    bool initialize() {
        // Add default themes
        addDefaultThemes();
        return true;
    }
    
    void shutdown() {
        // Nothing to do
    }
    
    const Theme& getCurrentTheme() const {
        return m_themes.at(m_currentTheme);
    }
    
    bool setCurrentTheme(const std::string& themeName) {
        if (m_themes.find(themeName) == m_themes.end()) {
            return false;
        }
        
        m_currentTheme = themeName;
        return true;
    }
    
    std::vector<std::string> getAvailableThemes() const {
        std::vector<std::string> themes;
        for (const auto& [name, _] : m_themes) {
            themes.push_back(name);
        }
        return themes;
    }
    
    bool addTheme(const Theme& theme) {
        // Check if the theme already exists
        if (m_themes.find(theme.name) != m_themes.end()) {
            return false;
        }
        
        // Add the theme
        m_themes[theme.name] = theme;
        return true;
    }
    
    bool removeTheme(const std::string& themeName) {
        // Check if the theme exists
        if (m_themes.find(themeName) == m_themes.end()) {
            return false;
        }
        
        // Check if the theme is the current theme
        if (m_currentTheme == themeName) {
            return false;
        }
        
        // Remove the theme
        m_themes.erase(themeName);
        return true;
    }
    
    bool loadThemes(const std::string& filename) {
        try {
            // Load the YAML file
            YAML::Node config = YAML::LoadFile(filename);
            
            // Check if the file contains a themes section
            if (!config["themes"]) {
                return false;
            }
            
            // Parse the themes
            const YAML::Node& themes = config["themes"];
            for (const auto& theme : themes) {
                // Get the theme name
                std::string name = theme.first.as<std::string>();
                
                // Get the theme colors
                std::unordered_map<std::string, std::string> colors;
                for (const auto& color : theme.second) {
                    std::string colorName = color.first.as<std::string>();
                    std::string colorValue = color.second.as<std::string>();
                    colors[colorName] = colorValue;
                }
                
                // Create the theme
                auto themeObj = Theme::fromConfig(name, colors);
                if (themeObj) {
                    addTheme(*themeObj);
                }
            }
            
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error loading themes: " << e.what() << std::endl;
            return false;
        }
    }
    
    bool saveThemes(const std::string& filename) {
        try {
            // Create the YAML node
            YAML::Node config;
            
            // Add the themes
            for (const auto& [name, theme] : m_themes) {
                config["themes"][name]["background"] = theme.background.toHex();
                config["themes"][name]["foreground"] = theme.foreground.toHex();
                config["themes"][name]["prompt"] = theme.prompt.toHex();
                config["themes"][name]["result"] = theme.result.toHex();
                config["themes"][name]["error"] = theme.error.toHex();
                config["themes"][name]["warning"] = theme.warning.toHex();
                config["themes"][name]["info"] = theme.info.toHex();
            }
            
            // Save the YAML file
            std::ofstream fout(filename);
            fout << config;
            
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error saving themes: " << e.what() << std::endl;
            return false;
        }
    }
    
    std::string formatText(const std::string& text, const std::string& colorType) const {
        // Get the current theme
        const Theme& theme = getCurrentTheme();
        
        // Get the color
        const Color* color = nullptr;
        if (colorType == "foreground") {
            color = &theme.foreground;
        } else if (colorType == "background") {
            color = &theme.background;
        } else if (colorType == "prompt") {
            color = &theme.prompt;
        } else if (colorType == "result") {
            color = &theme.result;
        } else if (colorType == "error") {
            color = &theme.error;
        } else if (colorType == "warning") {
            color = &theme.warning;
        } else if (colorType == "info") {
            color = &theme.info;
        } else {
            // Unknown color type, use foreground
            color = &theme.foreground;
        }
        
        // Format the text
        return color->toAnsi() + text + "\033[0m";
    }

private:
    void addDefaultThemes() {
        // Default theme
        addTheme(Theme("default",
                      Color::Black,
                      Color::White,
                      Color::Cyan,
                      Color::Green,
                      Color::Red,
                      Color::Yellow,
                      Color::Blue));
        
        // Dark theme
        addTheme(Theme("dark",
                      Color::Black,
                      Color::White,
                      Color::Cyan,
                      Color::Green,
                      Color::Red,
                      Color::Yellow,
                      Color::Blue));
        
        // Light theme
        addTheme(Theme("light",
                      Color::White,
                      Color::Black,
                      Color::Blue,
                      Color::Green,
                      Color::Red,
                      Color(255, 165, 0), // Orange
                      Color::Blue));
    }
    
    std::unordered_map<std::string, Theme> m_themes;
    std::string m_currentTheme;
};

// ThemeManager implementation

ThemeManager::ThemeManager() : m_impl(std::make_unique<Impl>()) {
}

ThemeManager::~ThemeManager() {
    shutdown();
}

bool ThemeManager::initialize() {
    return m_impl->initialize();
}

void ThemeManager::shutdown() {
    if (m_impl) {
        m_impl->shutdown();
    }
}

const Theme& ThemeManager::getCurrentTheme() const {
    return m_impl->getCurrentTheme();
}

bool ThemeManager::setCurrentTheme(const std::string& themeName) {
    return m_impl->setCurrentTheme(themeName);
}

std::vector<std::string> ThemeManager::getAvailableThemes() const {
    return m_impl->getAvailableThemes();
}

bool ThemeManager::addTheme(const Theme& theme) {
    return m_impl->addTheme(theme);
}

bool ThemeManager::removeTheme(const std::string& themeName) {
    return m_impl->removeTheme(themeName);
}

bool ThemeManager::loadThemes(const std::string& filename) {
    return m_impl->loadThemes(filename);
}

bool ThemeManager::saveThemes(const std::string& filename) {
    return m_impl->saveThemes(filename);
}

std::string ThemeManager::formatText(const std::string& text, const std::string& colorType) const {
    return m_impl->formatText(text, colorType);
}

} // namespace rebelcalc
