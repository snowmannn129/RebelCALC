#include "plot2d.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <fstream>

namespace RebelCalc {
namespace Visualization {

// DataSeries implementation

DataSeries::DataSeries(const std::vector<double>& x, const std::vector<double>& y, const std::string& name)
    : m_x(x), m_y(y), m_name(name) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("X and Y vectors must have the same size");
    }
}

DataSeries::DataSeries(std::function<double(double)> func, double xMin, double xMax, size_t numPoints, const std::string& name)
    : m_name(name) {
    if (xMin >= xMax) {
        throw std::invalid_argument("xMin must be less than xMax");
    }
    if (numPoints < 2) {
        throw std::invalid_argument("numPoints must be at least 2");
    }

    m_x.resize(numPoints);
    m_y.resize(numPoints);

    double step = (xMax - xMin) / (numPoints - 1);
    for (size_t i = 0; i < numPoints; ++i) {
        double x = xMin + i * step;
        m_x[i] = x;
        m_y[i] = func(x);
    }
}

const std::vector<double>& DataSeries::getX() const {
    return m_x;
}

const std::vector<double>& DataSeries::getY() const {
    return m_y;
}

const std::string& DataSeries::getName() const {
    return m_name;
}

void DataSeries::setLineStyle(LineStyle style) {
    m_lineStyle = style;
}

LineStyle DataSeries::getLineStyle() const {
    return m_lineStyle;
}

void DataSeries::setLineWidth(double width) {
    if (width <= 0) {
        throw std::invalid_argument("Line width must be positive");
    }
    m_lineWidth = width;
}

double DataSeries::getLineWidth() const {
    return m_lineWidth;
}

void DataSeries::setLineColor(const Color& color) {
    m_lineColor = color;
}

const Color& DataSeries::getLineColor() const {
    return m_lineColor;
}

void DataSeries::setMarkerStyle(MarkerStyle style) {
    m_markerStyle = style;
}

MarkerStyle DataSeries::getMarkerStyle() const {
    return m_markerStyle;
}

void DataSeries::setMarkerSize(double size) {
    if (size <= 0) {
        throw std::invalid_argument("Marker size must be positive");
    }
    m_markerSize = size;
}

double DataSeries::getMarkerSize() const {
    return m_markerSize;
}

void DataSeries::setMarkerColor(const Color& color) {
    m_markerColor = color;
}

const Color& DataSeries::getMarkerColor() const {
    return m_markerColor;
}

// Plot2D implementation

Plot2D::Plot2D(const std::string& title)
    : m_title(title) {
}

size_t Plot2D::addSeries(const DataSeries& series) {
    m_series.push_back(series);
    return m_series.size() - 1;
}

size_t Plot2D::addSeries(const std::vector<double>& x, const std::vector<double>& y, const std::string& name) {
    return addSeries(DataSeries(x, y, name));
}

size_t Plot2D::addSeries(std::function<double(double)> func, double xMin, double xMax, size_t numPoints, const std::string& name) {
    return addSeries(DataSeries(func, xMin, xMax, numPoints, name));
}

const DataSeries& Plot2D::getSeries(size_t index) const {
    if (index >= m_series.size()) {
        throw std::out_of_range("Series index out of range");
    }
    return m_series[index];
}

const DataSeries& Plot2D::getSeriesByName(const std::string& name) const {
    for (const auto& series : m_series) {
        if (series.getName() == name) {
            return series;
        }
    }
    throw std::out_of_range("Series not found: " + name);
}

size_t Plot2D::getSeriesCount() const {
    return m_series.size();
}

void Plot2D::setTitle(const std::string& title) {
    m_title = title;
}

const std::string& Plot2D::getTitle() const {
    return m_title;
}

void Plot2D::setXLabel(const std::string& label) {
    m_xLabel = label;
}

const std::string& Plot2D::getXLabel() const {
    return m_xLabel;
}

void Plot2D::setYLabel(const std::string& label) {
    m_yLabel = label;
}

const std::string& Plot2D::getYLabel() const {
    return m_yLabel;
}

void Plot2D::setXRange(double min, double max) {
    if (min >= max) {
        throw std::invalid_argument("xMin must be less than xMax");
    }
    m_xRange = std::make_pair(min, max);
}

std::pair<double, double> Plot2D::getXRange() const {
    if (m_xRange) {
        return *m_xRange;
    }
    return calculateXRange();
}

void Plot2D::setYRange(double min, double max) {
    if (min >= max) {
        throw std::invalid_argument("yMin must be less than yMax");
    }
    m_yRange = std::make_pair(min, max);
}

std::pair<double, double> Plot2D::getYRange() const {
    if (m_yRange) {
        return *m_yRange;
    }
    return calculateYRange();
}

void Plot2D::setShowGrid(bool show) {
    m_showGrid = show;
}

bool Plot2D::getShowGrid() const {
    return m_showGrid;
}

void Plot2D::setShowLegend(bool show) {
    m_showLegend = show;
}

bool Plot2D::getShowLegend() const {
    return m_showLegend;
}

void Plot2D::setBackgroundColor(const Color& color) {
    m_backgroundColor = color;
}

const Color& Plot2D::getBackgroundColor() const {
    return m_backgroundColor;
}

void Plot2D::setGridColor(const Color& color) {
    m_gridColor = color;
}

const Color& Plot2D::getGridColor() const {
    return m_gridColor;
}

void Plot2D::setTextColor(const Color& color) {
    m_textColor = color;
}

const Color& Plot2D::getTextColor() const {
    return m_textColor;
}

bool Plot2D::save(const std::string& filename, int width, int height) const {
    // Extract file extension
    std::string extension;
    size_t dotPos = filename.find_last_of('.');
    if (dotPos != std::string::npos) {
        extension = filename.substr(dotPos + 1);
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    }

    // Save based on file extension
    if (extension == "svg") {
        // Save as SVG
        std::ofstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        file << toSVG(width, height);
        return true;
    } else if (extension == "html") {
        // Save as HTML
        std::ofstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        file << toHTML(width, height);
        return true;
    } else {
        // Unsupported format
        return false;
    }
}

void Plot2D::show(int width, int height) const {
    // Generate a temporary HTML file and open it in the default browser
    std::string tempFilename = "temp_plot.html";
    if (save(tempFilename, width, height)) {
        // Open the file in the default browser (platform-dependent)
#ifdef _WIN32
        std::system(("start " + tempFilename).c_str());
#elif __APPLE__
        std::system(("open " + tempFilename).c_str());
#else
        std::system(("xdg-open " + tempFilename).c_str());
#endif
    }
}

// Helper functions for HTML/Plotly generation
namespace {

std::string colorToRGB(const Color& color) {
    std::stringstream ss;
    ss << "rgb(" << static_cast<int>(color.r) << "," 
       << static_cast<int>(color.g) << "," 
       << static_cast<int>(color.b) << ")";
    return ss.str();
}

std::string lineStyleToPlotly(LineStyle style) {
    switch (style) {
        case LineStyle::SOLID:
            return "solid";
        case LineStyle::DASHED:
            return "dash";
        case LineStyle::DOTTED:
            return "dot";
        case LineStyle::DASH_DOT:
            return "dashdot";
        default:
            return "solid";
    }
}

std::string markerStyleToPlotly(MarkerStyle style) {
    switch (style) {
        case MarkerStyle::CIRCLE:
            return "circle";
        case MarkerStyle::SQUARE:
            return "square";
        case MarkerStyle::TRIANGLE:
            return "triangle-up";
        case MarkerStyle::DIAMOND:
            return "diamond";
        case MarkerStyle::CROSS:
            return "cross";
        case MarkerStyle::PLUS:
            return "x";
        default:
            return "circle";
    }
}

std::string vectorToJSON(const std::vector<double>& data) {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < data.size(); ++i) {
        ss << data[i];
        if (i < data.size() - 1) {
            ss << ",";
        }
    }
    ss << "]";
    return ss.str();
}

} // anonymous namespace

std::string Plot2D::toHTML(int width, int height) const {
    std::ostringstream html;
    
    // HTML header
    html << "<!DOCTYPE html>\n";
    html << "<html>\n";
    html << "<head>\n";
    html << "  <title>" << (m_title.empty() ? "Plot" : m_title) << "</title>\n";
    html << "  <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n";
    html << "  <style>\n";
    html << "    body { font-family: Arial, sans-serif; margin: 20px; }\n";
    html << "    .container { max-width: " << width << "px; margin: 0 auto; }\n";
    html << "    h1 { text-align: center; }\n";
    html << "  </style>\n";
    html << "</head>\n";
    html << "<body>\n";
    html << "  <div class=\"container\">\n";
    
    if (!m_title.empty()) {
        html << "    <h1>" << m_title << "</h1>\n";
    }
    
    // Plot container
    html << "    <div id=\"plot\" style=\"width:" << width << "px;height:" << height << "px;\"></div>\n";
    
    // JavaScript for creating the plot
    html << "    <script>\n";
    html << "        var data = [];\n";
    
    // Add data series
    for (size_t i = 0; i < m_series.size(); ++i) {
        const auto& series = m_series[i];
        
        html << "        var trace" << i << " = {\n";
        html << "            x: " << vectorToJSON(series.getX()) << ",\n";
        html << "            y: " << vectorToJSON(series.getY()) << ",\n";
        html << "            name: '" << series.getName() << "',\n";
        html << "            type: 'scatter',\n";
        
        // Line properties
        html << "            line: {\n";
        html << "                color: '" << colorToRGB(series.getLineColor()) << "',\n";
        html << "                width: " << series.getLineWidth() << ",\n";
        html << "                dash: '" << lineStyleToPlotly(series.getLineStyle()) << "'\n";
        html << "            },\n";
        
        // Marker properties
        if (series.getMarkerStyle() != MarkerStyle::NONE) {
            html << "            mode: 'lines+markers',\n";
            html << "            marker: {\n";
            html << "                color: '" << colorToRGB(series.getMarkerColor()) << "',\n";
            html << "                size: " << series.getMarkerSize() << ",\n";
            html << "                symbol: '" << markerStyleToPlotly(series.getMarkerStyle()) << "'\n";
            html << "            }\n";
        } else {
            html << "            mode: 'lines'\n";
        }
        
        html << "        };\n";
        html << "        data.push(trace" << i << ");\n";
    }
    
    // Layout configuration
    html << "        var layout = {\n";
    html << "            title: '" << m_title << "',\n";
    html << "            autosize: true,\n";
    html << "            width: " << width << ",\n";
    html << "            height: " << height << ",\n";
    html << "            paper_bgcolor: '" << colorToRGB(m_backgroundColor) << "',\n";
    html << "            plot_bgcolor: 'white',\n";
    
    // X-axis configuration
    html << "            xaxis: {\n";
    html << "                title: '" << m_xLabel << "',\n";
    
    // X range
    auto xRange = getXRange();
    html << "                range: [" << xRange.first << ", " << xRange.second << "],\n";
    
    html << "                showgrid: " << (m_showGrid ? "true" : "false") << ",\n";
    html << "                gridcolor: '" << colorToRGB(m_gridColor) << "',\n";
    html << "                color: '" << colorToRGB(m_textColor) << "'\n";
    html << "            },\n";
    
    // Y-axis configuration
    html << "            yaxis: {\n";
    html << "                title: '" << m_yLabel << "',\n";
    
    // Y range
    auto yRange = getYRange();
    html << "                range: [" << yRange.first << ", " << yRange.second << "],\n";
    
    html << "                showgrid: " << (m_showGrid ? "true" : "false") << ",\n";
    html << "                gridcolor: '" << colorToRGB(m_gridColor) << "',\n";
    html << "                color: '" << colorToRGB(m_textColor) << "'\n";
    html << "            },\n";
    
    // Legend
    html << "            showlegend: " << (m_showLegend ? "true" : "false") << ",\n";
    html << "            legend: {\n";
    html << "                x: 1,\n";
    html << "                y: 1,\n";
    html << "                xanchor: 'right',\n";
    html << "                yanchor: 'top',\n";
    html << "                font: {\n";
    html << "                    color: '" << colorToRGB(m_textColor) << "'\n";
    html << "                }\n";
    html << "            },\n";
    
    // Interactive features
    html << "            hovermode: 'closest',\n";
    html << "            dragmode: 'zoom',\n";
    html << "            modebar: {\n";
    html << "                orientation: 'v',\n";
    html << "                activecolor: '#1f77b4'\n";
    html << "            }\n";
    
    html << "        };\n";
    
    // Create the plot
    html << "        Plotly.newPlot('plot', data, layout, {responsive: true});\n";
    
    // Add event listeners for interactivity
    html << "        document.getElementById('plot').on('plotly_click', function(data) {\n";
    html << "            var point = data.points[0];\n";
    html << "            alert('Clicked on point: (' + point.x + ', ' + point.y + ')');\n";
    html << "        });\n";
    
    html << "    </script>\n";
    
    // HTML footer
    html << "  </div>\n";
    html << "</body>\n";
    html << "</html>\n";
    
    return html.str();
}

std::string Plot2D::toSVG(int width, int height) const {
    std::ostringstream svg;
    
    // SVG header
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    svg << "<svg width=\"" << width << "\" height=\"" << height << "\" "
        << "xmlns=\"http://www.w3.org/2000/svg\">\n";
    
    // Background
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"rgb("
        << static_cast<int>(m_backgroundColor.r) << ","
        << static_cast<int>(m_backgroundColor.g) << ","
        << static_cast<int>(m_backgroundColor.b) << ")\" />\n";
    
    // Define plot area dimensions
    int margin = 50;
    int plotWidth = width - 2 * margin;
    int plotHeight = height - 2 * margin;
    int plotX = margin;
    int plotY = margin;
    
    // Draw plot area background
    svg << "<rect x=\"" << plotX << "\" y=\"" << plotY << "\" "
        << "width=\"" << plotWidth << "\" height=\"" << plotHeight << "\" "
        << "fill=\"white\" stroke=\"black\" stroke-width=\"1\" />\n";
    
    // Get data ranges
    auto xRange = getXRange();
    auto yRange = getYRange();
    double xMin = xRange.first;
    double xMax = xRange.second;
    double yMin = yRange.first;
    double yMax = yRange.second;
    
    // Draw grid if enabled
    if (m_showGrid) {
        svg << "<g stroke=\"rgb("
            << static_cast<int>(m_gridColor.r) << ","
            << static_cast<int>(m_gridColor.g) << ","
            << static_cast<int>(m_gridColor.b) << ")\" "
            << "stroke-width=\"0.5\" stroke-dasharray=\"5,5\">\n";
        
        // Vertical grid lines (x-axis)
        int numXGridLines = 10;
        for (int i = 0; i <= numXGridLines; ++i) {
            double x = xMin + (xMax - xMin) * i / numXGridLines;
            int xPos = plotX + static_cast<int>(plotWidth * (x - xMin) / (xMax - xMin));
            svg << "<line x1=\"" << xPos << "\" y1=\"" << plotY << "\" "
                << "x2=\"" << xPos << "\" y2=\"" << plotY + plotHeight << "\" />\n";
        }
        
        // Horizontal grid lines (y-axis)
        int numYGridLines = 10;
        for (int i = 0; i <= numYGridLines; ++i) {
            double y = yMin + (yMax - yMin) * i / numYGridLines;
            int yPos = plotY + plotHeight - static_cast<int>(plotHeight * (y - yMin) / (yMax - yMin));
            svg << "<line x1=\"" << plotX << "\" y1=\"" << yPos << "\" "
                << "x2=\"" << plotX + plotWidth << "\" y2=\"" << yPos << "\" />\n";
        }
        
        svg << "</g>\n";
    }
    
    // Draw data series
    for (const auto& series : m_series) {
        // Skip empty series
        if (series.getX().empty()) {
            continue;
        }
        
        // Set line style
        std::string dashArray;
        switch (series.getLineStyle()) {
            case LineStyle::SOLID:
                dashArray = "";
                break;
            case LineStyle::DASHED:
                dashArray = "10,5";
                break;
            case LineStyle::DOTTED:
                dashArray = "2,2";
                break;
            case LineStyle::DASH_DOT:
                dashArray = "10,5,2,5";
                break;
        }
        
        // Draw line
        svg << "<polyline points=\"";
        for (size_t i = 0; i < series.getX().size(); ++i) {
            double x = series.getX()[i];
            double y = series.getY()[i];
            
            // Skip points outside the plot range
            if (x < xMin || x > xMax || y < yMin || y > yMax) {
                continue;
            }
            
            // Convert data coordinates to SVG coordinates
            int xPos = plotX + static_cast<int>(plotWidth * (x - xMin) / (xMax - xMin));
            int yPos = plotY + plotHeight - static_cast<int>(plotHeight * (y - yMin) / (yMax - yMin));
            
            svg << xPos << "," << yPos << " ";
        }
        svg << "\" fill=\"none\" stroke=\"rgb("
            << static_cast<int>(series.getLineColor().r) << ","
            << static_cast<int>(series.getLineColor().g) << ","
            << static_cast<int>(series.getLineColor().b) << ")\" "
            << "stroke-width=\"" << series.getLineWidth() << "\" ";
        
        if (!dashArray.empty()) {
            svg << "stroke-dasharray=\"" << dashArray << "\" ";
        }
        
        svg << "/>\n";
        
        // Draw markers if enabled
        if (series.getMarkerStyle() != MarkerStyle::NONE) {
            svg << "<g fill=\"rgb("
                << static_cast<int>(series.getMarkerColor().r) << ","
                << static_cast<int>(series.getMarkerColor().g) << ","
                << static_cast<int>(series.getMarkerColor().b) << ")\">\n";
            
            for (size_t i = 0; i < series.getX().size(); ++i) {
                double x = series.getX()[i];
                double y = series.getY()[i];
                
                // Skip points outside the plot range
                if (x < xMin || x > xMax || y < yMin || y > yMax) {
                    continue;
                }
                
                // Convert data coordinates to SVG coordinates
                int xPos = plotX + static_cast<int>(plotWidth * (x - xMin) / (xMax - xMin));
                int yPos = plotY + plotHeight - static_cast<int>(plotHeight * (y - yMin) / (yMax - yMin));
                
                // Draw marker based on style
                double size = series.getMarkerSize();
                switch (series.getMarkerStyle()) {
                    case MarkerStyle::CIRCLE:
                        svg << "<circle cx=\"" << xPos << "\" cy=\"" << yPos << "\" "
                            << "r=\"" << size << "\" />\n";
                        break;
                    case MarkerStyle::SQUARE:
                        svg << "<rect x=\"" << xPos - size << "\" y=\"" << yPos - size << "\" "
                            << "width=\"" << 2 * size << "\" height=\"" << 2 * size << "\" />\n";
                        break;
                    case MarkerStyle::TRIANGLE:
                        svg << "<polygon points=\""
                            << xPos << "," << (yPos - size) << " "
                            << (xPos - size) << "," << (yPos + size) << " "
                            << (xPos + size) << "," << (yPos + size) << "\" />\n";
                        break;
                    case MarkerStyle::DIAMOND:
                        svg << "<polygon points=\""
                            << xPos << "," << (yPos - size) << " "
                            << (xPos + size) << "," << yPos << " "
                            << xPos << "," << (yPos + size) << " "
                            << (xPos - size) << "," << yPos << "\" />\n";
                        break;
                    case MarkerStyle::CROSS:
                        svg << "<line x1=\"" << (xPos - size) << "\" y1=\"" << (yPos - size) << "\" "
                            << "x2=\"" << (xPos + size) << "\" y2=\"" << (yPos + size) << "\" "
                            << "stroke=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" "
                            << "stroke-width=\"" << series.getLineWidth() << "\" />\n";
                        svg << "<line x1=\"" << (xPos - size) << "\" y1=\"" << (yPos + size) << "\" "
                            << "x2=\"" << (xPos + size) << "\" y2=\"" << (yPos - size) << "\" "
                            << "stroke=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" "
                            << "stroke-width=\"" << series.getLineWidth() << "\" />\n";
                        break;
                    case MarkerStyle::PLUS:
                        svg << "<line x1=\"" << xPos << "\" y1=\"" << (yPos - size) << "\" "
                            << "x2=\"" << xPos << "\" y2=\"" << (yPos + size) << "\" "
                            << "stroke=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" "
                            << "stroke-width=\"" << series.getLineWidth() << "\" />\n";
                        svg << "<line x1=\"" << (xPos - size) << "\" y1=\"" << yPos << "\" "
                            << "x2=\"" << (xPos + size) << "\" y2=\"" << yPos << "\" "
                            << "stroke=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" "
                            << "stroke-width=\"" << series.getLineWidth() << "\" />\n";
                        break;
                    default:
                        break;
                }
            }
            
            svg << "</g>\n";
        }
    }
    
    // Draw axes
    svg << "<g stroke=\"black\" stroke-width=\"1\">\n";
    svg << "<line x1=\"" << plotX << "\" y1=\"" << plotY + plotHeight << "\" "
        << "x2=\"" << plotX + plotWidth << "\" y2=\"" << plotY + plotHeight << "\" />\n";
    svg << "<line x1=\"" << plotX << "\" y1=\"" << plotY << "\" "
        << "x2=\"" << plotX << "\" y2=\"" << plotY + plotHeight << "\" />\n";
    svg << "</g>\n";
    
    // Draw axis labels
    svg << "<g fill=\"rgb("
        << static_cast<int>(m_textColor.r) << ","
        << static_cast<int>(m_textColor.g) << ","
        << static_cast<int>(m_textColor.b) << ")\" "
        << "font-family=\"Arial\" font-size=\"12\">\n";
    
    // X-axis labels
    int numXLabels = 10;
    for (int i = 0; i <= numXLabels; ++i) {
        double x = xMin + (xMax - xMin) * i / numXLabels;
        int xPos = plotX + static_cast<int>(plotWidth * i / numXLabels);
        
        svg << "<text x=\"" << xPos << "\" y=\"" << plotY + plotHeight + 20 << "\" "
            << "text-anchor=\"middle\">" << std::fixed << std::setprecision(2) << x << "</text>\n";
    }
    
    // Y-axis labels
    int numYLabels = 10;
    for (int i = 0; i <= numYLabels; ++i) {
        double y = yMin + (yMax - yMin) * i / numYLabels;
        int yPos = plotY + plotHeight - static_cast<int>(plotHeight * i / numYLabels);
        
        svg << "<text x=\"" << plotX - 10 << "\" y=\"" << yPos + 5 << "\" "
            << "text-anchor=\"end\">" << std::fixed << std::setprecision(2) << y << "</text>\n";
    }
    
    // Axis titles
    if (!m_xLabel.empty()) {
        svg << "<text x=\"" << plotX + plotWidth / 2 << "\" y=\"" << plotY + plotHeight + 40 << "\" "
            << "text-anchor=\"middle\" font-size=\"14\">" << m_xLabel << "</text>\n";
    }
    
    if (!m_yLabel.empty()) {
        svg << "<text x=\"" << plotX - 40 << "\" y=\"" << plotY + plotHeight / 2 << "\" "
            << "text-anchor=\"middle\" font-size=\"14\" "
            << "transform=\"rotate(-90 " << plotX - 40 << "," << plotY + plotHeight / 2 << ")\">"
            << m_yLabel << "</text>\n";
    }
    
    // Plot title
    if (!m_title.empty()) {
        svg << "<text x=\"" << plotX + plotWidth / 2 << "\" y=\"" << plotY - 20 << "\" "
            << "text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">"
            << m_title << "</text>\n";
    }
    
    svg << "</g>\n";
    
    // Draw legend if enabled
    if (m_showLegend && !m_series.empty()) {
        int legendX = plotX + plotWidth - 150;
        int legendY = plotY + 20;
        int legendWidth = 140;
        int legendHeight = static_cast<int>(m_series.size()) * 20 + 10;
        
        // Legend background
        svg << "<rect x=\"" << legendX << "\" y=\"" << legendY << "\" "
            << "width=\"" << legendWidth << "\" height=\"" << legendHeight << "\" "
            << "fill=\"white\" fill-opacity=\"0.8\" stroke=\"black\" stroke-width=\"1\" />\n";
        
        // Legend items
        for (size_t i = 0; i < m_series.size(); ++i) {
            const auto& series = m_series[i];
            int itemY = legendY + 15 + static_cast<int>(i) * 20;
            
            // Line sample
            svg << "<line x1=\"" << legendX + 10 << "\" y1=\"" << itemY << "\" "
                << "x2=\"" << legendX + 30 << "\" y2=\"" << itemY << "\" "
                << "stroke=\"rgb("
                << static_cast<int>(series.getLineColor().r) << ","
                << static_cast<int>(series.getLineColor().g) << ","
                << static_cast<int>(series.getLineColor().b) << ")\" "
                << "stroke-width=\"" << series.getLineWidth() << "\" ";
            
            // Set line style
            std::string dashArray;
            switch (series.getLineStyle()) {
                case LineStyle::SOLID:
                    dashArray = "";
                    break;
                case LineStyle::DASHED:
                    dashArray = "5,3";
                    break;
                case LineStyle::DOTTED:
                    dashArray = "1,1";
                    break;
                case LineStyle::DASH_DOT:
                    dashArray = "5,3,1,3";
                    break;
            }
            
            if (!dashArray.empty()) {
                svg << "stroke-dasharray=\"" << dashArray << "\" ";
            }
            
            svg << "/>\n";
            
            // Marker sample if applicable
            if (series.getMarkerStyle() != MarkerStyle::NONE) {
                int markerX = legendX + 20;
                double size = series.getMarkerSize();
                
                switch (series.getMarkerStyle()) {
                    case MarkerStyle::CIRCLE:
                        svg << "<circle cx=\"" << markerX << "\" cy=\"" << itemY << "\" "
                            << "r=\"" << size << "\" fill=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" />\n";
                        break;
                    case MarkerStyle::SQUARE:
                        svg << "<rect x=\"" << markerX - size << "\" y=\"" << itemY - size << "\" "
                            << "width=\"" << 2 * size << "\" height=\"" << 2 * size << "\" fill=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" />\n";
                        break;
                    case MarkerStyle::TRIANGLE:
                        svg << "<polygon points=\""
                            << markerX << "," << (itemY - size) << " "
                            << (markerX - size) << "," << (itemY + size) << " "
                            << (markerX + size) << "," << (itemY + size) << "\" fill=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" />\n";
                        break;
                    case MarkerStyle::DIAMOND:
                        svg << "<polygon points=\""
                            << markerX << "," << (itemY - size) << " "
                            << (markerX + size) << "," << itemY << " "
                            << markerX << "," << (itemY + size) << " "
                            << (markerX - size) << "," << itemY << "\" fill=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" />\n";
                        break;
                    case MarkerStyle::CROSS:
                        svg << "<line x1=\"" << (markerX - size) << "\" y1=\"" << (itemY - size) << "\" "
                            << "x2=\"" << (markerX + size) << "\" y2=\"" << (itemY + size) << "\" "
                            << "stroke=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" "
                            << "stroke-width=\"" << series.getLineWidth() << "\" />\n";
                        svg << "<line x1=\"" << (markerX - size) << "\" y1=\"" << (itemY + size) << "\" "
                            << "x2=\"" << (markerX + size) << "\" y2=\"" << (itemY - size) << "\" "
                            << "stroke=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" "
                            << "stroke-width=\"" << series.getLineWidth() << "\" />\n";
                        break;
                    case MarkerStyle::PLUS:
                        svg << "<line x1=\"" << markerX << "\" y1=\"" << (itemY - size) << "\" "
                            << "x2=\"" << markerX << "\" y2=\"" << (itemY + size) << "\" "
                            << "stroke=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" "
                            << "stroke-width=\"" << series.getLineWidth() << "\" />\n";
                        svg << "<line x1=\"" << (markerX - size) << "\" y1=\"" << itemY << "\" "
                            << "x2=\"" << (markerX + size) << "\" y2=\"" << itemY << "\" "
                            << "stroke=\"rgb("
                            << static_cast<int>(series.getMarkerColor().r) << ","
                            << static_cast<int>(series.getMarkerColor().g) << ","
                            << static_cast<int>(series.getMarkerColor().b) << ")\" "
                            << "stroke-width=\"" << series.getLineWidth() << "\" />\n";
                        break;
                    default:
                        break;
                }
            }
            
            // Series name
            svg << "<text x=\"" << legendX + 40 << "\" y=\"" << itemY + 5 << "\" "
                << "font-family=\"Arial\" font-size=\"12\" fill=\"rgb("
                << static_cast<int>(m_textColor.r) << ","
                << static_cast<int>(m_textColor.g) << ","
                << static_cast<int>(m_textColor.b) << ")\">"
                << series.getName() << "</text>\n";
        }
    }
    
    // SVG footer
    svg << "</svg>\n";
    
    return svg.str();
}

// Helper methods
std::pair<double, double> Plot2D::calculateXRange() const {
    if (m_series.empty()) {
        return std::make_pair(0.0, 1.0);
    }
    
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    
    for (const auto& series : m_series) {
        const auto& x = series.getX();
        if (!x.empty()) {
            min = std::min(min, *std::min_element(x.begin(), x.end()));
            max = std::max(max, *std::max_element(x.begin(), x.end()));
        }
    }
    
    if (min >= max) {
        min = 0.0;
        max = 1.0;
    }
    
    // Add a small margin
    double margin = (max - min) * 0.05;
    return std::make_pair(min - margin, max + margin);
}

std::pair<double, double> Plot2D::calculateYRange() const {
    if (m_series.empty()) {
        return std::make_pair(0.0, 1.0);
    }
    
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    
    for (const auto& series : m_series) {
        const auto& y = series.getY();
        if (!y.empty()) {
            min = std::min(min, *std::min_element(y.begin(), y.end()));
            max = std::max(max, *std::max_element(y.begin(), y.end()));
        }
    }
    
    if (min >= max) {
        min = 0.0;
        max = 1.0;
    }
    
    // Add a small margin
    double margin = (max - min) * 0.05;
    return std::make_pair(min - margin, max + margin);
}

} // namespace Visualization
} // namespace RebelCalc
