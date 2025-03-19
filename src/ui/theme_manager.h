#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <optional>

namespace rebelcalc {

/**
 * Represents a color in the terminal
 */
struct Color {
    uint8_t r;
    uint8_t g;
    uint8_t b;
    
    /**
     * Create a color from RGB values
     * @param r Red component (0-255)
     * @param g Green component (0-255)
     * @param b Blue component (0-255)
     */
    Color(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
    
    /**
     * Create a color from a hex string
     * @param hex Hex string in the format "#RRGGBB" or "RRGGBB"
     */
    static std::optional<Color> fromHex(const std::string& hex);
    
    /**
     * Convert the color to a hex string
     * @return Hex string in the format "#RRGGBB"
     */
    std::string toHex() const;
    
    /**
     * Convert the color to an ANSI escape sequence
     * @param foreground Whether this is a foreground color (true) or background color (false)
     * @return ANSI escape sequence
     */
    std::string toAnsi(bool foreground = true) const;
    
    // Predefined colors
    static const Color Black;
    static const Color Red;
    static const Color Green;
    static const Color Yellow;
    static const Color Blue;
    static const Color Magenta;
    static const Color Cyan;
    static const Color White;
    static const Color BrightBlack;
    static const Color BrightRed;
    static const Color BrightGreen;
    static const Color BrightYellow;
    static const Color BrightBlue;
    static const Color BrightMagenta;
    static const Color BrightCyan;
    static const Color BrightWhite;
};

/**
 * Represents a theme for the terminal UI
 */
struct Theme {
    std::string name;
    Color background;
    Color foreground;
    Color prompt;
    Color result;
    Color error;
    Color warning;
    Color info;
    
    /**
     * Create a theme with the given colors
     * @param name The name of the theme
     * @param background The background color
     * @param foreground The foreground color
     * @param prompt The prompt color
     * @param result The result color
     * @param error The error color
     * @param warning The warning color
     * @param info The info color
     */
    Theme(const std::string& name,
          const Color& background,
          const Color& foreground,
          const Color& prompt,
          const Color& result,
          const Color& error,
          const Color& warning,
          const Color& info)
        : name(name),
          background(background),
          foreground(foreground),
          prompt(prompt),
          result(result),
          error(error),
          warning(warning),
          info(info) {}
    
    /**
     * Create a theme from a configuration
     * @param name The name of the theme
     * @param config The configuration map
     * @return The theme, or nullopt if the configuration is invalid
     */
    static std::optional<Theme> fromConfig(const std::string& name,
                                          const std::unordered_map<std::string, std::string>& config);
};

/**
 * Class for managing themes for the terminal UI
 */
class ThemeManager {
public:
    /**
     * Constructor
     */
    ThemeManager();
    
    /**
     * Destructor
     */
    ~ThemeManager();
    
    /**
     * Initialize the theme manager
     * @return true if initialization was successful, false otherwise
     */
    bool initialize();
    
    /**
     * Shutdown the theme manager
     */
    void shutdown();
    
    /**
     * Get the current theme
     * @return The current theme
     */
    const Theme& getCurrentTheme() const;
    
    /**
     * Set the current theme
     * @param themeName The name of the theme to set
     * @return true if the theme was set successfully, false otherwise
     */
    bool setCurrentTheme(const std::string& themeName);
    
    /**
     * Get a list of available themes
     * @return A vector of theme names
     */
    std::vector<std::string> getAvailableThemes() const;
    
    /**
     * Add a theme
     * @param theme The theme to add
     * @return true if the theme was added successfully, false otherwise
     */
    bool addTheme(const Theme& theme);
    
    /**
     * Remove a theme
     * @param themeName The name of the theme to remove
     * @return true if the theme was removed successfully, false otherwise
     */
    bool removeTheme(const std::string& themeName);
    
    /**
     * Load themes from a configuration file
     * @param filename The name of the configuration file
     * @return true if the themes were loaded successfully, false otherwise
     */
    bool loadThemes(const std::string& filename);
    
    /**
     * Save themes to a configuration file
     * @param filename The name of the configuration file
     * @return true if the themes were saved successfully, false otherwise
     */
    bool saveThemes(const std::string& filename);
    
    /**
     * Format text with the current theme
     * @param text The text to format
     * @param colorType The type of color to use (foreground, background, prompt, result, error, warning, info)
     * @return The formatted text
     */
    std::string formatText(const std::string& text, const std::string& colorType) const;

private:
    // Implementation details
    class Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace rebelcalc
