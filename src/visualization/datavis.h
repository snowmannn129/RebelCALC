#ifndef REBELCALC_VISUALIZATION_DATAVIS_H
#define REBELCALC_VISUALIZATION_DATAVIS_H

#include <vector>
#include <string>
#include <map>
#include <functional>
#include <optional>
#include "plot2d.h" // For Color and other shared types

namespace RebelCalc {
namespace Visualization {

/**
 * @brief Enum for chart types
 */
enum class ChartType {
    BAR,
    STACKED_BAR,
    GROUPED_BAR,
    PIE,
    DONUT,
    HISTOGRAM,
    BOX_PLOT,
    VIOLIN_PLOT,
    HEATMAP,
    TREE_MAP,
    RADAR,
    BUBBLE
};

/**
 * @brief Class for a data point in a chart
 */
class DataPoint {
public:
    /**
     * @brief Constructor for a data point with a single value
     * @param label Label for the data point
     * @param value Value of the data point
     */
    DataPoint(const std::string& label, double value);

    /**
     * @brief Constructor for a data point with multiple values
     * @param label Label for the data point
     * @param values Values of the data point
     */
    DataPoint(const std::string& label, const std::vector<double>& values);

    /**
     * @brief Get the label of the data point
     * @return Label of the data point
     */
    const std::string& getLabel() const;

    /**
     * @brief Get the value of the data point
     * @return Value of the data point
     * @throws std::runtime_error if the data point has multiple values
     */
    double getValue() const;

    /**
     * @brief Get the values of the data point
     * @return Values of the data point
     */
    const std::vector<double>& getValues() const;

    /**
     * @brief Check if the data point has multiple values
     * @return true if the data point has multiple values, false otherwise
     */
    bool hasMultipleValues() const;

    /**
     * @brief Set the color of the data point
     * @param color Color of the data point
     */
    void setColor(const Color& color);

    /**
     * @brief Get the color of the data point
     * @return Color of the data point
     */
    const Color& getColor() const;

    /**
     * @brief Set the highlight state of the data point
     * @param highlight Whether the data point is highlighted
     */
    void setHighlight(bool highlight);

    /**
     * @brief Check if the data point is highlighted
     * @return true if the data point is highlighted, false otherwise
     */
    bool isHighlighted() const;

private:
    std::string m_label;
    std::vector<double> m_values;
    Color m_color = Color::Blue();
    bool m_highlight = false;
};

/**
 * @brief Class for a data series in a chart
 */
class ChartSeries {
public:
    /**
     * @brief Constructor for a chart series
     * @param name Name of the series
     */
    ChartSeries(const std::string& name = "");

    /**
     * @brief Add a data point to the series
     * @param point Data point to add
     */
    void addDataPoint(const DataPoint& point);

    /**
     * @brief Add a data point to the series
     * @param label Label for the data point
     * @param value Value of the data point
     */
    void addDataPoint(const std::string& label, double value);

    /**
     * @brief Add a data point to the series
     * @param label Label for the data point
     * @param values Values of the data point
     */
    void addDataPoint(const std::string& label, const std::vector<double>& values);

    /**
     * @brief Get the data points in the series
     * @return Data points in the series
     */
    const std::vector<DataPoint>& getDataPoints() const;

    /**
     * @brief Get the name of the series
     * @return Name of the series
     */
    const std::string& getName() const;

    /**
     * @brief Set the color of the series
     * @param color Color of the series
     */
    void setColor(const Color& color);

    /**
     * @brief Get the color of the series
     * @return Color of the series
     */
    const Color& getColor() const;

private:
    std::string m_name;
    std::vector<DataPoint> m_dataPoints;
    Color m_color = Color::Blue();
};

/**
 * @brief Class for a chart
 */
class Chart {
public:
    /**
     * @brief Constructor for a chart
     * @param type Type of the chart
     * @param title Title of the chart
     */
    Chart(ChartType type, const std::string& title = "");

    /**
     * @brief Add a data series to the chart
     * @param series Data series to add
     * @return Index of the added series
     */
    size_t addSeries(const ChartSeries& series);

    /**
     * @brief Add a data series to the chart
     * @param name Name of the series
     * @return Index of the added series
     */
    size_t addSeries(const std::string& name = "");

    /**
     * @brief Get a data series
     * @param index Index of the series
     * @return Data series
     */
    const ChartSeries& getSeries(size_t index) const;

    /**
     * @brief Get a data series by name
     * @param name Name of the series
     * @return Data series
     * @throws std::out_of_range if the series doesn't exist
     */
    const ChartSeries& getSeriesByName(const std::string& name) const;

    /**
     * @brief Get a mutable reference to a data series
     * @param index Index of the series
     * @return Mutable reference to the data series
     */
    ChartSeries& getSeries(size_t index);

    /**
     * @brief Get a mutable reference to a data series by name
     * @param name Name of the series
     * @return Mutable reference to the data series
     * @throws std::out_of_range if the series doesn't exist
     */
    ChartSeries& getSeriesByName(const std::string& name);

    /**
     * @brief Get the number of data series
     * @return Number of data series
     */
    size_t getSeriesCount() const;

    /**
     * @brief Set the title of the chart
     * @param title Title of the chart
     */
    void setTitle(const std::string& title);

    /**
     * @brief Get the title of the chart
     * @return Title of the chart
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
     * @brief Set the chart type
     * @param type Chart type
     */
    void setChartType(ChartType type);

    /**
     * @brief Get the chart type
     * @return Chart type
     */
    ChartType getChartType() const;

    /**
     * @brief Set the color palette
     * @param palette Color palette as a vector of colors
     */
    void setColorPalette(const std::vector<Color>& palette);

    /**
     * @brief Get the color palette
     * @return Color palette
     */
    const std::vector<Color>& getColorPalette() const;

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
     * @brief Set whether to show data labels
     * @param show Whether to show data labels
     */
    void setShowDataLabels(bool show);

    /**
     * @brief Get whether to show data labels
     * @return Whether to show data labels
     */
    bool getShowDataLabels() const;

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
     * @brief Save the chart to a file
     * @param filename Filename to save to
     * @param width Width of the image in pixels
     * @param height Height of the image in pixels
     * @return true if the chart was saved successfully, false otherwise
     */
    bool save(const std::string& filename, int width = 800, int height = 600) const;

    /**
     * @brief Show the chart in a window
     * @param width Width of the window in pixels
     * @param height Height of the window in pixels
     */
    void show(int width = 800, int height = 600) const;

    /**
     * @brief Generate SVG representation of the chart
     * @param width Width of the SVG in pixels
     * @param height Height of the SVG in pixels
     * @return SVG string
     */
    std::string toSVG(int width = 800, int height = 600) const;

    /**
     * @brief Generate HTML representation of the chart
     * @param width Width of the chart in pixels
     * @param height Height of the chart in pixels
     * @return HTML string
     */
    std::string toHTML(int width = 800, int height = 600) const;

private:
    ChartType m_type;
    std::string m_title;
    std::string m_xLabel;
    std::string m_yLabel;
    std::vector<ChartSeries> m_series;
    std::vector<Color> m_colorPalette = {
        Color::Blue(),
        Color::Red(),
        Color::Green(),
        Color::Yellow(),
        Color::Cyan(),
        Color::Magenta(),
        Color::Orange(),
        Color::Purple(),
        Color::Brown(),
        Color::Pink()
    };
    bool m_showLegend = true;
    bool m_showDataLabels = false;
    Color m_backgroundColor = Color::White();
    Color m_textColor = Color::Black();

    // Helper methods for different chart types
    std::string generateBarChartSVG(int width, int height) const;
    std::string generateStackedBarChartSVG(int width, int height) const;
    std::string generateGroupedBarChartSVG(int width, int height) const;
    std::string generatePieChartSVG(int width, int height) const;
    std::string generateDonutChartSVG(int width, int height) const;
    std::string generateHistogramSVG(int width, int height) const;
    std::string generateBoxPlotSVG(int width, int height) const;
    std::string generateViolinPlotSVG(int width, int height) const;
    std::string generateHeatmapSVG(int width, int height) const;
    std::string generateTreeMapSVG(int width, int height) const;
    std::string generateRadarChartSVG(int width, int height) const;
    std::string generateBubbleChartSVG(int width, int height) const;
};

/**
 * @brief Class for a dashboard containing multiple charts
 */
class Dashboard {
public:
    /**
     * @brief Constructor for a dashboard
     * @param title Title of the dashboard
     */
    Dashboard(const std::string& title = "");

    /**
     * @brief Add a chart to the dashboard
     * @param chart Chart to add
     * @param row Row position (0-based)
     * @param col Column position (0-based)
     * @param rowSpan Number of rows the chart spans
     * @param colSpan Number of columns the chart spans
     */
    void addChart(const Chart& chart, int row, int col, int rowSpan = 1, int colSpan = 1);

    /**
     * @brief Set the title of the dashboard
     * @param title Title of the dashboard
     */
    void setTitle(const std::string& title);

    /**
     * @brief Get the title of the dashboard
     * @return Title of the dashboard
     */
    const std::string& getTitle() const;

    /**
     * @brief Set the number of rows in the dashboard
     * @param rows Number of rows
     */
    void setRows(int rows);

    /**
     * @brief Get the number of rows in the dashboard
     * @return Number of rows
     */
    int getRows() const;

    /**
     * @brief Set the number of columns in the dashboard
     * @param cols Number of columns
     */
    void setColumns(int cols);

    /**
     * @brief Get the number of columns in the dashboard
     * @return Number of columns
     */
    int getColumns() const;

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
     * @brief Save the dashboard to a file
     * @param filename Filename to save to
     * @param width Width of the image in pixels
     * @param height Height of the image in pixels
     * @return true if the dashboard was saved successfully, false otherwise
     */
    bool save(const std::string& filename, int width = 1200, int height = 800) const;

    /**
     * @brief Show the dashboard in a window
     * @param width Width of the window in pixels
     * @param height Height of the window in pixels
     */
    void show(int width = 1200, int height = 800) const;

    /**
     * @brief Generate HTML representation of the dashboard
     * @param width Width of the dashboard in pixels
     * @param height Height of the dashboard in pixels
     * @return HTML string
     */
    std::string toHTML(int width = 1200, int height = 800) const;

private:
    struct ChartInfo {
        Chart chart;
        int row;
        int col;
        int rowSpan;
        int colSpan;
    };

    std::string m_title;
    std::vector<ChartInfo> m_charts;
    int m_rows = 2;
    int m_cols = 2;
    Color m_backgroundColor = Color::White();
    Color m_textColor = Color::Black();
};

/**
 * @brief Class for a data table
 */
class DataTable {
public:
    /**
     * @brief Constructor for a data table
     * @param title Title of the table
     */
    DataTable(const std::string& title = "");

    /**
     * @brief Set the column names
     * @param columnNames Column names
     */
    void setColumnNames(const std::vector<std::string>& columnNames);

    /**
     * @brief Get the column names
     * @return Column names
     */
    const std::vector<std::string>& getColumnNames() const;

    /**
     * @brief Add a row to the table
     * @param row Row data
     */
    void addRow(const std::vector<std::string>& row);

    /**
     * @brief Get the rows in the table
     * @return Rows in the table
     */
    const std::vector<std::vector<std::string>>& getRows() const;

    /**
     * @brief Get the number of columns in the table
     * @return Number of columns
     */
    size_t getColumnCount() const;

    /**
     * @brief Get the number of rows in the table
     * @return Number of rows
     */
    size_t getRowCount() const;

    /**
     * @brief Set the title of the table
     * @param title Title of the table
     */
    void setTitle(const std::string& title);

    /**
     * @brief Get the title of the table
     * @return Title of the table
     */
    const std::string& getTitle() const;

    /**
     * @brief Set whether to show the header
     * @param show Whether to show the header
     */
    void setShowHeader(bool show);

    /**
     * @brief Get whether to show the header
     * @return Whether to show the header
     */
    bool getShowHeader() const;

    /**
     * @brief Set whether to show alternating row colors
     * @param alternate Whether to show alternating row colors
     */
    void setAlternateRowColors(bool alternate);

    /**
     * @brief Get whether to show alternating row colors
     * @return Whether to show alternating row colors
     */
    bool getAlternateRowColors() const;

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
     * @brief Set the header background color
     * @param color Header background color
     */
    void setHeaderBackgroundColor(const Color& color);

    /**
     * @brief Get the header background color
     * @return Header background color
     */
    const Color& getHeaderBackgroundColor() const;

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
     * @brief Set the header text color
     * @param color Header text color
     */
    void setHeaderTextColor(const Color& color);

    /**
     * @brief Get the header text color
     * @return Header text color
     */
    const Color& getHeaderTextColor() const;

    /**
     * @brief Set the alternate row color
     * @param color Alternate row color
     */
    void setAlternateRowColor(const Color& color);

    /**
     * @brief Get the alternate row color
     * @return Alternate row color
     */
    const Color& getAlternateRowColor() const;

    /**
     * @brief Generate HTML representation of the table
     * @return HTML string
     */
    std::string toHTML() const;

private:
    std::string m_title;
    std::vector<std::string> m_columnNames;
    std::vector<std::vector<std::string>> m_rows;
    bool m_showHeader = true;
    bool m_alternateRowColors = true;
    Color m_backgroundColor = Color::White();
    Color m_headerBackgroundColor = Color(240, 240, 240);
    Color m_textColor = Color::Black();
    Color m_headerTextColor = Color::Black();
    Color m_alternateRowColor = Color(245, 245, 245);
};

} // namespace Visualization
} // namespace RebelCalc

#endif // REBELCALC_VISUALIZATION_DATAVIS_H
