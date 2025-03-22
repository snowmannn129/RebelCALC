#ifndef REBELCALC_VISUALIZATION_PLOT2D_H
#define REBELCALC_VISUALIZATION_PLOT2D_H

#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <map>
#include <array>
#include <optional>

namespace RebelCalc {
namespace Visualization {

/**
 * @brief Enum for line styles
 */
enum class LineStyle {
    SOLID,
    DASHED,
    DOTTED,
    DASH_DOT
};

/**
 * @brief Enum for marker styles
 */
enum class MarkerStyle {
    NONE,
    CIRCLE,
    SQUARE,
    TRIANGLE,
    DIAMOND,
    CROSS,
    PLUS
};

/**
 * @brief Struct for RGB color
 */
struct Color {
    unsigned char r, g, b;
    unsigned char a = 255; // Alpha (transparency)

    Color() : r(0), g(0), b(0), a(255) {}
    Color(unsigned char r, unsigned char g, unsigned char b, unsigned char a = 255)
        : r(r), g(g), b(b), a(a) {}

    // Predefined colors
    static Color Black() { return Color(0, 0, 0); }
    static Color White() { return Color(255, 255, 255); }
    static Color Red() { return Color(255, 0, 0); }
    static Color Green() { return Color(0, 255, 0); }
    static Color Blue() { return Color(0, 0, 255); }
    static Color Yellow() { return Color(255, 255, 0); }
    static Color Cyan() { return Color(0, 255, 255); }
    static Color Magenta() { return Color(255, 0, 255); }
    static Color Gray() { return Color(128, 128, 128); }
    static Color Orange() { return Color(255, 165, 0); }
    static Color Purple() { return Color(128, 0, 128); }
    static Color Brown() { return Color(165, 42, 42); }
    static Color Pink() { return Color(255, 192, 203); }
};

/**
 * @brief Class for a data series in a 2D plot
 */
class DataSeries {
public:
    /**
     * @brief Constructor for a data series with x and y values
     * @param x X values
     * @param y Y values
     * @param name Name of the series
     */
    DataSeries(const std::vector<double>& x, const std::vector<double>& y, const std::string& name = "");

    /**
     * @brief Constructor for a data series with a function
     * @param func Function to plot
     * @param xMin Minimum x value
     * @param xMax Maximum x value
     * @param numPoints Number of points to evaluate
     * @param name Name of the series
     */
    DataSeries(std::function<double(double)> func, double xMin, double xMax, size_t numPoints = 100, const std::string& name = "");

    /**
     * @brief Get the x values
     * @return X values
     */
    const std::vector<double>& getX() const;

    /**
     * @brief Get the y values
     * @return Y values
     */
    const std::vector<double>& getY() const;

    /**
     * @brief Get the name of the series
     * @return Name of the series
     */
    const std::string& getName() const;

    /**
     * @brief Set the line style
     * @param style Line style
     */
    void setLineStyle(LineStyle style);

    /**
     * @brief Get the line style
     * @return Line style
     */
    LineStyle getLineStyle() const;

    /**
     * @brief Set the line width
     * @param width Line width
     */
    void setLineWidth(double width);

    /**
     * @brief Get the line width
     * @return Line width
     */
    double getLineWidth() const;

    /**
     * @brief Set the line color
     * @param color Line color
     */
    void setLineColor(const Color& color);

    /**
     * @brief Get the line color
     * @return Line color
     */
    const Color& getLineColor() const;

    /**
     * @brief Set the marker style
     * @param style Marker style
     */
    void setMarkerStyle(MarkerStyle style);

    /**
     * @brief Get the marker style
     * @return Marker style
     */
    MarkerStyle getMarkerStyle() const;

    /**
     * @brief Set the marker size
     * @param size Marker size
     */
    void setMarkerSize(double size);

    /**
     * @brief Get the marker size
     * @return Marker size
     */
    double getMarkerSize() const;

    /**
     * @brief Set the marker color
     * @param color Marker color
     */
    void setMarkerColor(const Color& color);

    /**
     * @brief Get the marker color
     * @return Marker color
     */
    const Color& getMarkerColor() const;

private:
    std::vector<double> m_x;
    std::vector<double> m_y;
    std::string m_name;
    LineStyle m_lineStyle = LineStyle::SOLID;
    double m_lineWidth = 1.0;
    Color m_lineColor = Color::Blue();
    MarkerStyle m_markerStyle = MarkerStyle::NONE;
    double m_markerSize = 5.0;
    Color m_markerColor = Color::Red();
};

/**
 * @brief Class for a 2D plot
 */
class Plot2D {
public:
    /**
     * @brief Constructor for a 2D plot
     * @param title Title of the plot
     */
    Plot2D(const std::string& title = "");

    /**
     * @brief Add a data series to the plot
     * @param series Data series to add
     * @return Index of the added series
     */
    size_t addSeries(const DataSeries& series);

    /**
     * @brief Add a data series to the plot with x and y values
     * @param x X values
     * @param y Y values
     * @param name Name of the series
     * @return Index of the added series
     */
    size_t addSeries(const std::vector<double>& x, const std::vector<double>& y, const std::string& name = "");

    /**
     * @brief Add a data series to the plot with a function
     * @param func Function to plot
     * @param xMin Minimum x value
     * @param xMax Maximum x value
     * @param numPoints Number of points to evaluate
     * @param name Name of the series
     * @return Index of the added series
     */
    size_t addSeries(std::function<double(double)> func, double xMin, double xMax, size_t numPoints = 100, const std::string& name = "");

    /**
     * @brief Get a data series
     * @param index Index of the series
     * @return Data series
     */
    const DataSeries& getSeries(size_t index) const;

    /**
     * @brief Get a data series by name
     * @param name Name of the series
     * @return Data series
     * @throws std::out_of_range if the series doesn't exist
     */
    const DataSeries& getSeriesByName(const std::string& name) const;

    /**
     * @brief Get the number of data series
     * @return Number of data series
     */
    size_t getSeriesCount() const;

    /**
     * @brief Set the title of the plot
     * @param title Title of the plot
     */
    void setTitle(const std::string& title);

    /**
     * @brief Get the title of the plot
     * @return Title of the plot
     */
    const std::string& getTitle() const;

    /**
     * @brief Set the x-axis label
     * @param label X-axis label
     */
    void setXLabel(const std::string& label);

    /**
     * @brief Get the x-axis label
     * @return X-axis label
     */
    const std::string& getXLabel() const;

    /**
     * @brief Set the y-axis label
     * @param label Y-axis label
     */
    void setYLabel(const std::string& label);

    /**
     * @brief Get the y-axis label
     * @return Y-axis label
     */
    const std::string& getYLabel() const;

    /**
     * @brief Set the x-axis range
     * @param min Minimum x value
     * @param max Maximum x value
     */
    void setXRange(double min, double max);

    /**
     * @brief Get the x-axis range
     * @return Pair of minimum and maximum x values
     */
    std::pair<double, double> getXRange() const;

    /**
     * @brief Set the y-axis range
     * @param min Minimum y value
     * @param max Maximum y value
     */
    void setYRange(double min, double max);

    /**
     * @brief Get the y-axis range
     * @return Pair of minimum and maximum y values
     */
    std::pair<double, double> getYRange() const;

    /**
     * @brief Set whether to show the grid
     * @param show Whether to show the grid
     */
    void setShowGrid(bool show);

    /**
     * @brief Get whether to show the grid
     * @return Whether to show the grid
     */
    bool getShowGrid() const;

    /**
     * @brief Set whether to show the legend
     * @param show Whether to show the legend
     */
    void setShowLegend(bool show);

    /**
     * @brief Get whether to show the legend
     * @return Whether to show the legend
     */
    bool getShowLegend() const;

    /**
     * @brief Set the background color
     * @param color Background color
     */
    void setBackgroundColor(const Color& color);

    /**
     * @brief Get the background color
     * @return Background color
     */
    const Color& getBackgroundColor() const;

    /**
     * @brief Set the grid color
     * @param color Grid color
     */
    void setGridColor(const Color& color);

    /**
     * @brief Get the grid color
     * @return Grid color
     */
    const Color& getGridColor() const;

    /**
     * @brief Set the text color
     * @param color Text color
     */
    void setTextColor(const Color& color);

    /**
     * @brief Get the text color
     * @return Text color
     */
    const Color& getTextColor() const;

    /**
     * @brief Save the plot to a file
     * @param filename Filename to save to
     * @param width Width of the image in pixels
     * @param height Height of the image in pixels
     * @return true if the plot was saved successfully, false otherwise
     */
    bool save(const std::string& filename, int width = 800, int height = 600) const;

    /**
     * @brief Show the plot in a window
     * @param width Width of the window in pixels
     * @param height Height of the window in pixels
     */
    void show(int width = 800, int height = 600) const;

    /**
     * @brief Generate SVG representation of the plot
     * @param width Width of the SVG in pixels
     * @param height Height of the SVG in pixels
     * @return SVG string
     */
    std::string toSVG(int width = 800, int height = 600) const;

    /**
     * @brief Generate HTML representation of the plot
     * @param width Width of the plot in pixels
     * @param height Height of the plot in pixels
     * @return HTML string
     */
    std::string toHTML(int width = 800, int height = 600) const;

private:
    std::string m_title;
    std::string m_xLabel;
    std::string m_yLabel;
    std::vector<DataSeries> m_series;
    std::optional<std::pair<double, double>> m_xRange;
    std::optional<std::pair<double, double>> m_yRange;
    bool m_showGrid = true;
    bool m_showLegend = true;
    Color m_backgroundColor = Color::White();
    Color m_gridColor = Color(200, 200, 200); // Light gray
    Color m_textColor = Color::Black();

    // Helper methods
    std::pair<double, double> calculateXRange() const;
    std::pair<double, double> calculateYRange() const;
};

} // namespace Visualization
} // namespace RebelCalc

#endif // REBELCALC_VISUALIZATION_PLOT2D_H
