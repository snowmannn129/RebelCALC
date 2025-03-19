#include <gtest/gtest.h>
#include "../ui/theme_manager.h"

namespace rebelcalc {
namespace testing {

class ThemeManagerTest : public ::testing::Test {
protected:
    void SetUp() override {
        themeManager = std::make_shared<ThemeManager>();
        ASSERT_TRUE(themeManager->initialize());
    }
    
    void TearDown() override {
        themeManager->shutdown();
        themeManager.reset();
    }
    
    std::shared_ptr<ThemeManager> themeManager;
};

TEST_F(ThemeManagerTest, GetCurrentTheme) {
    const Theme& theme = themeManager->getCurrentTheme();
    EXPECT_EQ("default", theme.name);
}

TEST_F(ThemeManagerTest, SetCurrentTheme) {
    // Set to an existing theme
    EXPECT_TRUE(themeManager->setCurrentTheme("dark"));
    EXPECT_EQ("dark", themeManager->getCurrentTheme().name);
    
    // Set to a non-existing theme
    EXPECT_FALSE(themeManager->setCurrentTheme("nonexistent"));
    EXPECT_EQ("dark", themeManager->getCurrentTheme().name);
}

TEST_F(ThemeManagerTest, GetAvailableThemes) {
    auto themes = themeManager->getAvailableThemes();
    EXPECT_FALSE(themes.empty());
    
    // Check for default themes
    EXPECT_TRUE(std::find(themes.begin(), themes.end(), "default") != themes.end());
    EXPECT_TRUE(std::find(themes.begin(), themes.end(), "dark") != themes.end());
    EXPECT_TRUE(std::find(themes.begin(), themes.end(), "light") != themes.end());
}

TEST_F(ThemeManagerTest, AddTheme) {
    // Create a new theme
    Theme theme("custom",
               Color(0, 0, 0),
               Color(255, 255, 255),
               Color(0, 255, 255),
               Color(0, 255, 0),
               Color(255, 0, 0),
               Color(255, 255, 0),
               Color(0, 0, 255));
    
    // Add the theme
    EXPECT_TRUE(themeManager->addTheme(theme));
    
    // Check if the theme was added
    auto themes = themeManager->getAvailableThemes();
    EXPECT_TRUE(std::find(themes.begin(), themes.end(), "custom") != themes.end());
    
    // Try to add the same theme again
    EXPECT_FALSE(themeManager->addTheme(theme));
}

TEST_F(ThemeManagerTest, RemoveTheme) {
    // Create and add a new theme
    Theme theme("custom",
               Color(0, 0, 0),
               Color(255, 255, 255),
               Color(0, 255, 255),
               Color(0, 255, 0),
               Color(255, 0, 0),
               Color(255, 255, 0),
               Color(0, 0, 255));
    
    EXPECT_TRUE(themeManager->addTheme(theme));
    
    // Remove the theme
    EXPECT_TRUE(themeManager->removeTheme("custom"));
    
    // Check if the theme was removed
    auto themes = themeManager->getAvailableThemes();
    EXPECT_TRUE(std::find(themes.begin(), themes.end(), "custom") == themes.end());
    
    // Try to remove a non-existing theme
    EXPECT_FALSE(themeManager->removeTheme("nonexistent"));
    
    // Try to remove the current theme
    EXPECT_FALSE(themeManager->removeTheme("default"));
}

TEST_F(ThemeManagerTest, FormatText) {
    // Format text with different color types
    std::string text = "Hello, world!";
    
    std::string formattedText = themeManager->formatText(text, "foreground");
    EXPECT_NE(text, formattedText);
    
    formattedText = themeManager->formatText(text, "background");
    EXPECT_NE(text, formattedText);
    
    formattedText = themeManager->formatText(text, "prompt");
    EXPECT_NE(text, formattedText);
    
    formattedText = themeManager->formatText(text, "result");
    EXPECT_NE(text, formattedText);
    
    formattedText = themeManager->formatText(text, "error");
    EXPECT_NE(text, formattedText);
    
    formattedText = themeManager->formatText(text, "warning");
    EXPECT_NE(text, formattedText);
    
    formattedText = themeManager->formatText(text, "info");
    EXPECT_NE(text, formattedText);
    
    // Format text with an unknown color type
    formattedText = themeManager->formatText(text, "unknown");
    EXPECT_NE(text, formattedText);
}

TEST_F(ThemeManagerTest, ColorFromHex) {
    // Valid hex strings
    auto color = Color::fromHex("#FF0000");
    ASSERT_TRUE(color.has_value());
    EXPECT_EQ(255, color->r);
    EXPECT_EQ(0, color->g);
    EXPECT_EQ(0, color->b);
    
    color = Color::fromHex("00FF00");
    ASSERT_TRUE(color.has_value());
    EXPECT_EQ(0, color->r);
    EXPECT_EQ(255, color->g);
    EXPECT_EQ(0, color->b);
    
    // Invalid hex strings
    color = Color::fromHex("#FF00");
    EXPECT_FALSE(color.has_value());
    
    color = Color::fromHex("FFGG00");
    EXPECT_FALSE(color.has_value());
}

TEST_F(ThemeManagerTest, ColorToHex) {
    Color color(255, 0, 0);
    EXPECT_EQ("#ff0000", color.toHex());
    
    color = Color(0, 255, 0);
    EXPECT_EQ("#00ff00", color.toHex());
    
    color = Color(0, 0, 255);
    EXPECT_EQ("#0000ff", color.toHex());
}

TEST_F(ThemeManagerTest, ColorToAnsi) {
    Color color(255, 0, 0);
    std::string ansi = color.toAnsi();
    EXPECT_FALSE(ansi.empty());
    
    ansi = color.toAnsi(false);
    EXPECT_FALSE(ansi.empty());
}

TEST_F(ThemeManagerTest, ThemeFromConfig) {
    // Valid configuration
    std::unordered_map<std::string, std::string> config = {
        {"background", "#000000"},
        {"foreground", "#FFFFFF"},
        {"prompt", "#00FFFF"},
        {"result", "#00FF00"},
        {"error", "#FF0000"},
        {"warning", "#FFFF00"},
        {"info", "#0000FF"}
    };
    
    auto theme = Theme::fromConfig("custom", config);
    ASSERT_TRUE(theme.has_value());
    EXPECT_EQ("custom", theme->name);
    
    // Invalid configuration (missing color)
    config.erase("error");
    theme = Theme::fromConfig("custom", config);
    EXPECT_FALSE(theme.has_value());
    
    // Invalid configuration (invalid color)
    config["error"] = "#FF00";
    theme = Theme::fromConfig("custom", config);
    EXPECT_FALSE(theme.has_value());
}

} // namespace testing
} // namespace rebelcalc
