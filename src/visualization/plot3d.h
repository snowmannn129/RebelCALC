#ifndef REBELCALC_VISUALIZATION_PLOT3D_H
#define REBELCALC_VISUALIZATION_PLOT3D_H

#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <map>
#include <array>
#include <optional>
#include "plot2d.h" // For Color and other shared types

namespace RebelCalc {
namespace Visualization {

/**
 * @brief Enum for 3D surface styles
 */
enum class SurfaceStyle {
    SOLID,
    WIREFRAME,
    POINTS,
    CONTOUR
};

/**
 * @brief Class for a 3D data series
 */
class DataSeries3D {
public:
    /**
     * @brief Constructor for a 3D data series with x, y, and z values
     * @param x X values
     * @param y Y values
     * @param z Z values (must be a 2D grid of size x.size() * y.size())
     * @param name Name of the series
     */
    DataSeries3D(const std::vector<double>& x, const std::vector<double>& y, 
                const std::vector<std::vector<double>>& z, const std::string& name = "");

    /**
     * @brief Constructor for a 3D data series with a function
     * @param func Function to plot (takes x and y, returns z)
     * @param xMin Minimum x value
     * @param xMax Maximum x value
     * @param yMin Minimum y value
     * @param yMax Maximum y value
     * @param xPoints Number of points in x direction
     * @param yPoints Number of points in y direction
     * @param name Name of the series
     */
    DataSeries3D(std::function<double(double, double)> func, 
                double xMin, double xMax, double yMin, double yMax,
                size_t xPoints = 50, size_t yPoints = 50, 
                const std::string& name = "");

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
     * @brief Get the z values
     * @return Z values as a 2D grid
     */
    const std::vector<std::vector<double>>& getZ() const;

    /**
     * @brief Get the name of the series
     * @return Name of the series
     */
    const std::string& getName() const;

    /**
     * @brief Set the surface style
     * @param style Surface style
     */
    void setSurfaceStyle(SurfaceStyle style);

    /**
     * @brief Get the surface style
     * @return Surface style
     */
    SurfaceStyle getSurfaceStyle() const;

    /**
     * @brief Set the color map
     * @param colorMap Color map as a vector of colors
     */
    void setColorMap(const std::vector<Color>& colorMap);

    /**
     * @brief Get the color map
     * @return Color map
     */
    const std::vector<Color>& getColorMap() const;

    /**
     * @brief Set the wireframe color
     * @param color Wireframe color
     */
    void setWireframeColor(const Color& color);

    /**
     * @brief Get the wireframe color
     * @return Wireframe color
     */
    const Color& getWireframeColor() const;

    /**
     * @brief Set the wireframe width
     * @param width Wireframe width
     */
    void setWireframeWidth(double width);

    /**
     * @brief Get the wireframe width
     * @return Wireframe width
     */
    double getWireframeWidth() const;

    /**
     * @brief Set the point size
     * @param size Point size
     */
    void setPointSize(double size);

    /**
     * @brief Get the point size
     * @return Point size
     */
    double getPointSize() const;

    /**
     * @brief Set the contour levels
     * @param levels Contour levels
     */
    void setContourLevels(const std::vector<double>& levels);

    /**
     * @brief Get the contour levels
     * @return Contour levels
     */
    const std::vector<double>& getContourLevels() const;

private:
    std::vector<double> m_x;
    std::vector<double> m_y;
    std::vector<std::vector<double>> m_z;
    std::string m_name;
    SurfaceStyle m_surfaceStyle = SurfaceStyle::SOLID;
    std::vector<Color> m_colorMap;
    Color m_wireframeColor = Color::Black();
    double m_wireframeWidth = 1.0;
    double m_pointSize = 3.0;
    std::vector<double> m_contourLevels;
};

/**
 * @brief Class for a 3D scatter data series
 */
class ScatterSeries3D {
public:
    /**
     * @brief Constructor for a 3D scatter data series
     * @param x X values
     * @param y Y values
     * @param z Z values
     * @param name Name of the series
     */
    ScatterSeries3D(const std::vector<double>& x, const std::vector<double>& y, 
                   const std::vector<double>& z, const std::string& name = "");

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
     * @brief Get the z values
     * @return Z values
     */
    const std::vector<double>& getZ() const;

    /**
     * @brief Get the name of the series
     * @return Name of the series
     */
    const std::string& getName() const;

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

    /**
     * @brief Set the color map for coloring points by z value
     * @param colorMap Color map as a vector of colors
     */
    void setColorMap(const std::vector<Color>& colorMap);

    /**
     * @brief Get the color map
     * @return Color map
     */
    const std::vector<Color>& getColorMap() const;

    /**
     * @brief Set whether to color points by z value
     * @param colorByZ Whether to color points by z value
     */
    void setColorByZ(bool colorByZ);

    /**
     * @brief Get whether to color points by z value
     * @return Whether to color points by z value
     */
    bool getColorByZ() const;

private:
    std::vector<double> m_x;
    std::vector<double> m_y;
    std::vector<double> m_z;
    std::string m_name;
    MarkerStyle m_markerStyle = MarkerStyle::CIRCLE;
    double m_markerSize = 5.0;
    Color m_markerColor = Color::Blue();
    std::vector<Color> m_colorMap;
    bool m_colorByZ = false;
};

/**
 * @brief Class for a 3D plot
 */
class Plot3D {
public:
    /**
     * @brief Constructor for a 3D plot
     * @param title Title of the plot
     */
    Plot3D(const std::string& title = "");

    /**
     * @brief Add a surface data series to the plot
     * @param series Surface data series to add
     * @return Index of the added series
     */
    size_t addSurface(const DataSeries3D& series);

    /**
     * @brief Add a surface data series to the plot with a function
     * @param func Function to plot (takes x and y, returns z)
     * @param xMin Minimum x value
     * @param xMax Maximum x value
     * @param yMin Minimum y value
     * @param yMax Maximum y value
     * @param xPoints Number of points in x direction
     * @param yPoints Number of points in y direction
     * @param name Name of the series
     * @return Index of the added series
     */
    size_t addSurface(std::function<double(double, double)> func, 
                     double xMin, double xMax, double yMin, double yMax,
                     size_t xPoints = 50, size_t yPoints = 50, 
                     const std::string& name = "");

    /**
     * @brief Add a scatter data series to the plot
     * @param series Scatter data series to add
     * @return Index of the added series
     */
    size_t addScatter(const ScatterSeries3D& series);

    /**
     * @brief Add a scatter data series to the plot
     * @param x X values
     * @param y Y values
     * @param z Z values
     * @param name Name of the series
     * @return Index of the added series
     */
    size_t addScatter(const std::vector<double>& x, const std::vector<double>& y, 
                     const std::vector<double>& z, const std::string& name = "");

    /**
     * @brief Get a surface data series
     * @param index Index of the series
     * @return Surface data series
     */
    const DataSeries3D& getSurface(size_t index) const;

    /**
     * @brief Get a scatter data series
     * @param index Index of the series
     * @return Scatter data series
     */
    const ScatterSeries3D& getScatter(size_t index) const;

    /**
     * @brief Get the number of surface data series
     * @return Number of surface data series
     */
    size_t getSurfaceCount() const;

    /**
     * @brief Get the number of scatter data series
     * @return Number of scatter data series
     */
    size_t getScatterCount() const;

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
     * @brief Set the z-axis label
     * @param label Z-axis label
     */
    void setZLabel(const std::string& label);

    /**
     * @brief Get the z-axis label
     * @return Z-axis label
     */
    const std::string& getZLabel() const;

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
     * @brief Set the z-axis range
     * @param min Minimum z value
     * @param max Maximum z value
     */
    void setZRange(double min, double max);

    /**
     * @brief Get the z-axis range
     * @return Pair of minimum and maximum z values
     */
    std::pair<double, double> getZRange() const;

    /**
     * @brief Set the view angles
     * @param azimuth Azimuth angle in degrees
     * @param elevation Elevation angle in degrees
     */
    void setViewAngles(double azimuth, double elevation);

    /**
     * @brief Get the view angles
     * @return Pair of azimuth and elevation angles in degrees
     */
    std::pair<double, double> getViewAngles() const;

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
     * @brief Set whether to show the color bar
     * @param show Whether to show the color bar
     */
    void setShowColorBar(bool show);

    /**
     * @brief Get whether to show the color bar
     * @return Whether to show the color bar
     */
    bool getShowColorBar() const;

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
    std::string m_zLabel;
    std::vector<DataSeries3D> m_surfaces;
    std::vector<ScatterSeries3D> m_scatters;
    std::optional<std::pair<double, double>> m_xRange;
    std::optional<std::pair<double, double>> m_yRange;
    std::optional<std::pair<double, double>> m_zRange;
    double m_azimuth = 30.0;
    double m_elevation = 30.0;
    bool m_showGrid = true;
    bool m_showLegend = true;
    bool m_showColorBar = true;
    Color m_backgroundColor = Color::White();
    Color m_gridColor = Color(200, 200, 200); // Light gray
    Color m_textColor = Color::Black();

    // Helper methods
    std::pair<double, double> calculateXRange() const;
    std::pair<double, double> calculateYRange() const;
    std::pair<double, double> calculateZRange() const;
};

} // namespace Visualization
} // namespace RebelCalc

#endif // REBELCALC_VISUALIZATION_PLOT3D_H
