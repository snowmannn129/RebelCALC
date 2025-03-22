#include "plot3d.h"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <fstream>

namespace RebelCalc {
namespace Visualization {

// DataSeries3D implementation
DataSeries3D::DataSeries3D(const std::vector<double>& x, const std::vector<double>& y, 
                         const std::vector<std::vector<double>>& z, const std::string& name)
    : m_x(x), m_y(y), m_z(z), m_name(name) {
}

DataSeries3D::DataSeries3D(std::function<double(double, double)> func, 
                         double xMin, double xMax, double yMin, double yMax,
                         size_t xPoints, size_t yPoints, 
                         const std::string& name)
    : m_name(name) {
    // Generate x and y values
    m_x.resize(xPoints);
    m_y.resize(yPoints);
    
    double xStep = (xMax - xMin) / (xPoints - 1);
    double yStep = (yMax - yMin) / (yPoints - 1);
    
    for (size_t i = 0; i < xPoints; ++i) {
        m_x[i] = xMin + i * xStep;
    }
    
    for (size_t i = 0; i < yPoints; ++i) {
        m_y[i] = yMin + i * yStep;
    }
    
    // Generate z values
    m_z.resize(xPoints);
    for (size_t i = 0; i < xPoints; ++i) {
        m_z[i].resize(yPoints);
        for (size_t j = 0; j < yPoints; ++j) {
            m_z[i][j] = func(m_x[i], m_y[j]);
        }
    }
}

const std::vector<double>& DataSeries3D::getX() const {
    return m_x;
}

const std::vector<double>& DataSeries3D::getY() const {
    return m_y;
}

const std::vector<std::vector<double>>& DataSeries3D::getZ() const {
    return m_z;
}

const std::string& DataSeries3D::getName() const {
    return m_name;
}

void DataSeries3D::setSurfaceStyle(SurfaceStyle style) {
    m_surfaceStyle = style;
}

SurfaceStyle DataSeries3D::getSurfaceStyle() const {
    return m_surfaceStyle;
}

void DataSeries3D::setColorMap(const std::vector<Color>& colorMap) {
    m_colorMap = colorMap;
}

const std::vector<Color>& DataSeries3D::getColorMap() const {
    return m_colorMap;
}

void DataSeries3D::setWireframeColor(const Color& color) {
    m_wireframeColor = color;
}

const Color& DataSeries3D::getWireframeColor() const {
    return m_wireframeColor;
}

void DataSeries3D::setWireframeWidth(double width) {
    m_wireframeWidth = width;
}

double DataSeries3D::getWireframeWidth() const {
    return m_wireframeWidth;
}

void DataSeries3D::setPointSize(double size) {
    m_pointSize = size;
}

double DataSeries3D::getPointSize() const {
    return m_pointSize;
}

void DataSeries3D::setContourLevels(const std::vector<double>& levels) {
    m_contourLevels = levels;
}

const std::vector<double>& DataSeries3D::getContourLevels() const {
    return m_contourLevels;
}

// ScatterSeries3D implementation
ScatterSeries3D::ScatterSeries3D(const std::vector<double>& x, const std::vector<double>& y, 
                               const std::vector<double>& z, const std::string& name)
    : m_x(x), m_y(y), m_z(z), m_name(name) {
}

const std::vector<double>& ScatterSeries3D::getX() const {
    return m_x;
}

const std::vector<double>& ScatterSeries3D::getY() const {
    return m_y;
}

const std::vector<double>& ScatterSeries3D::getZ() const {
    return m_z;
}

const std::string& ScatterSeries3D::getName() const {
    return m_name;
}

void ScatterSeries3D::setMarkerStyle(MarkerStyle style) {
    m_markerStyle = style;
}

MarkerStyle ScatterSeries3D::getMarkerStyle() const {
    return m_markerStyle;
}

void ScatterSeries3D::setMarkerSize(double size) {
    m_markerSize = size;
}

double ScatterSeries3D::getMarkerSize() const {
    return m_markerSize;
}

void ScatterSeries3D::setMarkerColor(const Color& color) {
    m_markerColor = color;
}

const Color& ScatterSeries3D::getMarkerColor() const {
    return m_markerColor;
}

void ScatterSeries3D::setColorMap(const std::vector<Color>& colorMap) {
    m_colorMap = colorMap;
}

const std::vector<Color>& ScatterSeries3D::getColorMap() const {
    return m_colorMap;
}

void ScatterSeries3D::setColorByZ(bool colorByZ) {
    m_colorByZ = colorByZ;
}

bool ScatterSeries3D::getColorByZ() const {
    return m_colorByZ;
}

// Plot3D implementation
Plot3D::Plot3D(const std::string& title)
    : m_title(title),
      m_xLabel("X"),
      m_yLabel("Y"),
      m_zLabel("Z") {
}

size_t Plot3D::addSurface(const DataSeries3D& series) {
    m_surfaces.push_back(series);
    return m_surfaces.size() - 1;
}

size_t Plot3D::addSurface(std::function<double(double, double)> func, 
                         double xMin, double xMax, double yMin, double yMax,
                         size_t xPoints, size_t yPoints, 
                         const std::string& name) {
    DataSeries3D series(func, xMin, xMax, yMin, yMax, xPoints, yPoints, name);
    return addSurface(series);
}

size_t Plot3D::addScatter(const ScatterSeries3D& series) {
    m_scatters.push_back(series);
    return m_scatters.size() - 1;
}

size_t Plot3D::addScatter(const std::vector<double>& x, const std::vector<double>& y, 
                         const std::vector<double>& z, const std::string& name) {
    ScatterSeries3D series(x, y, z, name);
    return addScatter(series);
}

const DataSeries3D& Plot3D::getSurface(size_t index) const {
    if (index >= m_surfaces.size()) {
        throw std::out_of_range("Surface index out of range");
    }
    return m_surfaces[index];
}

const ScatterSeries3D& Plot3D::getScatter(size_t index) const {
    if (index >= m_scatters.size()) {
        throw std::out_of_range("Scatter index out of range");
    }
    return m_scatters[index];
}

size_t Plot3D::getSurfaceCount() const {
    return m_surfaces.size();
}

size_t Plot3D::getScatterCount() const {
    return m_scatters.size();
}

void Plot3D::setTitle(const std::string& title) {
    m_title = title;
}

const std::string& Plot3D::getTitle() const {
    return m_title;
}

void Plot3D::setXLabel(const std::string& label) {
    m_xLabel = label;
}

const std::string& Plot3D::getXLabel() const {
    return m_xLabel;
}

void Plot3D::setYLabel(const std::string& label) {
    m_yLabel = label;
}

const std::string& Plot3D::getYLabel() const {
    return m_yLabel;
}

void Plot3D::setZLabel(const std::string& label) {
    m_zLabel = label;
}

const std::string& Plot3D::getZLabel() const {
    return m_zLabel;
}

void Plot3D::setXRange(double min, double max) {
    if (min < max) {
        m_xRange = std::make_pair(min, max);
    }
}

std::pair<double, double> Plot3D::getXRange() const {
    if (m_xRange) {
        return *m_xRange;
    }
    return calculateXRange();
}

void Plot3D::setYRange(double min, double max) {
    if (min < max) {
        m_yRange = std::make_pair(min, max);
    }
}

std::pair<double, double> Plot3D::getYRange() const {
    if (m_yRange) {
        return *m_yRange;
    }
    return calculateYRange();
}

void Plot3D::setZRange(double min, double max) {
    if (min < max) {
        m_zRange = std::make_pair(min, max);
    }
}

std::pair<double, double> Plot3D::getZRange() const {
    if (m_zRange) {
        return *m_zRange;
    }
    return calculateZRange();
}

void Plot3D::setViewAngles(double azimuth, double elevation) {
    m_azimuth = azimuth;
    m_elevation = elevation;
}

std::pair<double, double> Plot3D::getViewAngles() const {
    return std::make_pair(m_azimuth, m_elevation);
}

void Plot3D::setShowGrid(bool show) {
    m_showGrid = show;
}

bool Plot3D::getShowGrid() const {
    return m_showGrid;
}

void Plot3D::setShowLegend(bool show) {
    m_showLegend = show;
}

bool Plot3D::getShowLegend() const {
    return m_showLegend;
}

void Plot3D::setShowColorBar(bool show) {
    m_showColorBar = show;
}

bool Plot3D::getShowColorBar() const {
    return m_showColorBar;
}

void Plot3D::setBackgroundColor(const Color& color) {
    m_backgroundColor = color;
}

const Color& Plot3D::getBackgroundColor() const {
    return m_backgroundColor;
}

void Plot3D::setGridColor(const Color& color) {
    m_gridColor = color;
}

const Color& Plot3D::getGridColor() const {
    return m_gridColor;
}

void Plot3D::setTextColor(const Color& color) {
    m_textColor = color;
}

const Color& Plot3D::getTextColor() const {
    return m_textColor;
}

// Helper methods for HTML generation
namespace {

std::string colorToString(const Color& color) {
    std::stringstream ss;
    ss << "rgb(" << color.r << ", " << color.g << ", " << color.b << ")";
    return ss.str();
}

std::string markerStyleToString(MarkerStyle style) {
    switch (style) {
        case MarkerStyle::CIRCLE:
            return "circle";
        case MarkerStyle::SQUARE:
            return "square";
        case MarkerStyle::DIAMOND:
            return "diamond";
        case MarkerStyle::TRIANGLE:
            return "triangle-up";
        case MarkerStyle::CROSS:
            return "cross";
        default:
            return "circle";
    }
}

std::string colorMapToJSON(const std::vector<Color>& colorMap) {
    std::stringstream ss;
    ss << "[";
    
    for (size_t i = 0; i < colorMap.size(); ++i) {
        ss << "[" << static_cast<double>(i) / (colorMap.size() - 1) << ", '" 
           << colorToString(colorMap[i]) << "']";
        
        if (i < colorMap.size() - 1) {
            ss << ", ";
        }
    }
    
    ss << "]";
    return ss.str();
}

std::string scatterDataToJSON(const std::vector<double>& data) {
    std::stringstream ss;
    ss << "[";
    
    for (size_t i = 0; i < data.size(); ++i) {
        ss << data[i];
        if (i < data.size() - 1) {
            ss << ", ";
        }
    }
    
    ss << "]";
    return ss.str();
}

std::string surfaceXToJSON(const DataSeries3D& surface) {
    std::stringstream ss;
    ss << "[";
    
    const auto& x = surface.getX();
    
    // Create a 2D array of x values
    for (size_t i = 0; i < x.size(); ++i) {
        ss << "[";
        for (size_t j = 0; j < surface.getY().size(); ++j) {
            ss << x[i];
            if (j < surface.getY().size() - 1) {
                ss << ", ";
            }
        }
        ss << "]";
        if (i < x.size() - 1) {
            ss << ", ";
        }
    }
    
    ss << "]";
    return ss.str();
}

std::string surfaceYToJSON(const DataSeries3D& surface) {
    std::stringstream ss;
    ss << "[";
    
    const auto& y = surface.getY();
    
    // Create a 2D array of y values
    for (size_t i = 0; i < surface.getX().size(); ++i) {
        ss << "[";
        for (size_t j = 0; j < y.size(); ++j) {
            ss << y[j];
            if (j < y.size() - 1) {
                ss << ", ";
            }
        }
        ss << "]";
        if (i < surface.getX().size() - 1) {
            ss << ", ";
        }
    }
    
    ss << "]";
    return ss.str();
}

std::string surfaceZToJSON(const DataSeries3D& surface) {
    std::stringstream ss;
    ss << "[";
    
    const auto& z = surface.getZ();
    
    for (size_t i = 0; i < z.size(); ++i) {
        ss << "[";
        for (size_t j = 0; j < z[i].size(); ++j) {
            ss << z[i][j];
            if (j < z[i].size() - 1) {
                ss << ", ";
            }
        }
        ss << "]";
        if (i < z.size() - 1) {
            ss << ", ";
        }
    }
    
    ss << "]";
    return ss.str();
}

} // anonymous namespace

std::string Plot3D::toHTML(int width, int height) const {
    std::stringstream ss;
    
    // HTML header
    ss << "<!DOCTYPE html>\n";
    ss << "<html>\n";
    ss << "<head>\n";
    ss << "    <title>" << m_title << "</title>\n";
    ss << "    <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n";
    ss << "</head>\n";
    ss << "<body>\n";
    
    // Plot container
    ss << "    <div id=\"plot\" style=\"width:" << width << "px;height:" << height << "px;\"></div>\n";
    
    // JavaScript for creating the plot
    ss << "    <script>\n";
    ss << "        var data = [];\n";
    
    // Add surfaces
    for (size_t i = 0; i < m_surfaces.size(); ++i) {
        const DataSeries3D& surface = m_surfaces[i];
        
        ss << "        var surface" << i << " = {\n";
        ss << "            type: 'surface',\n";
        ss << "            name: '" << surface.getName() << "',\n";
        
        // X, Y, Z data
        ss << "            x: " << surfaceXToJSON(surface) << ",\n";
        ss << "            y: " << surfaceYToJSON(surface) << ",\n";
        ss << "            z: " << surfaceZToJSON(surface) << ",\n";
        
        // Appearance
        ss << "            showscale: " << (m_showColorBar ? "true" : "false") << ",\n";
        
        // Surface style
        switch (surface.getSurfaceStyle()) {
            case SurfaceStyle::WIREFRAME:
                ss << "            hidesurface: true,\n";
                ss << "            contours: {\n";
                ss << "                x: { show: true, color: '" << colorToString(surface.getWireframeColor()) << "', width: " << surface.getWireframeWidth() << " },\n";
                ss << "                y: { show: true, color: '" << colorToString(surface.getWireframeColor()) << "', width: " << surface.getWireframeWidth() << " },\n";
                ss << "                z: { show: true, color: '" << colorToString(surface.getWireframeColor()) << "', width: " << surface.getWireframeWidth() << " }\n";
                ss << "            },\n";
                break;
            case SurfaceStyle::POINTS:
                ss << "            hidesurface: true,\n";
                ss << "            surfacecolor: 'rgb(0,0,0)',\n"; // Dummy color, not used
                ss << "            mode: 'markers',\n";
                ss << "            marker: {\n";
                ss << "                size: " << surface.getPointSize() << ",\n";
                ss << "                color: '" << colorToString(surface.getWireframeColor()) << "'\n";
                ss << "            },\n";
                break;
            case SurfaceStyle::CONTOUR:
                ss << "            contours: {\n";
                ss << "                z: { show: true, usecolormap: true, project: { z: true } }\n";
                ss << "            },\n";
                
                // Add contour levels if specified
                if (!surface.getContourLevels().empty()) {
                    ss << "            contours: {\n";
                    ss << "                z: {\n";
                    ss << "                    show: true,\n";
                    ss << "                    usecolormap: true,\n";
                    ss << "                    project: { z: true },\n";
                    ss << "                    start: " << surface.getContourLevels().front() << ",\n";
                    ss << "                    end: " << surface.getContourLevels().back() << ",\n";
                    ss << "                    size: " << (surface.getContourLevels().back() - surface.getContourLevels().front()) / surface.getContourLevels().size() << "\n";
                    ss << "                }\n";
                    ss << "            },\n";
                }
                break;
            case SurfaceStyle::SOLID:
            default:
                // Default is solid surface
                break;
        }
        
        // Add custom color map if specified
        if (!surface.getColorMap().empty()) {
            ss << "            colorscale: " << colorMapToJSON(surface.getColorMap()) << ",\n";
        }
        
        ss << "        };\n";
        ss << "        data.push(surface" << i << ");\n";
    }
    
    // Add scatter data
    for (size_t i = 0; i < m_scatters.size(); ++i) {
        const ScatterSeries3D& scatter = m_scatters[i];
        
        ss << "        var scatter" << i << " = {\n";
        ss << "            type: 'scatter3d',\n";
        ss << "            name: '" << scatter.getName() << "',\n";
        
        // X, Y, Z data
        ss << "            x: " << scatterDataToJSON(scatter.getX()) << ",\n";
        ss << "            y: " << scatterDataToJSON(scatter.getY()) << ",\n";
        ss << "            z: " << scatterDataToJSON(scatter.getZ()) << ",\n";
        
        // Appearance
        ss << "            mode: 'markers',\n";
        ss << "            marker: {\n";
        ss << "                size: " << scatter.getMarkerSize() << ",\n";
        
        // Color by Z value or use a single color
        if (scatter.getColorByZ()) {
            ss << "                color: " << scatterDataToJSON(scatter.getZ()) << ",\n";
            
            // Add custom color map if specified
            if (!scatter.getColorMap().empty()) {
                ss << "                colorscale: " << colorMapToJSON(scatter.getColorMap()) << ",\n";
            }
            
            ss << "                showscale: " << (m_showColorBar ? "true" : "false") << ",\n";
        } else {
            ss << "                color: '" << colorToString(scatter.getMarkerColor()) << "',\n";
        }
        
        // Marker style
        ss << "                symbol: '" << markerStyleToString(scatter.getMarkerStyle()) << "'\n";
        
        ss << "            }\n";
        ss << "        };\n";
        ss << "        data.push(scatter" << i << ");\n";
    }
    
    // Layout configuration
    ss << "        var layout = {\n";
    ss << "            title: '" << m_title << "',\n";
    ss << "            autosize: true,\n";
    ss << "            width: " << width << ",\n";
    ss << "            height: " << height << ",\n";
    ss << "            paper_bgcolor: '" << colorToString(m_backgroundColor) << "',\n";
    ss << "            plot_bgcolor: '" << colorToString(m_backgroundColor) << "',\n";
    
    // Axes
    ss << "            scene: {\n";
    ss << "                xaxis: {\n";
    ss << "                    title: '" << m_xLabel << "',\n";
    
    // X range
    auto xRange = getXRange();
    ss << "                    range: [" << xRange.first << ", " << xRange.second << "],\n";
    
    ss << "                    showgrid: " << (m_showGrid ? "true" : "false") << ",\n";
    ss << "                    gridcolor: '" << colorToString(m_gridColor) << "',\n";
    ss << "                    color: '" << colorToString(m_textColor) << "'\n";
    ss << "                },\n";
    
    ss << "                yaxis: {\n";
    ss << "                    title: '" << m_yLabel << "',\n";
    
    // Y range
    auto yRange = getYRange();
    ss << "                    range: [" << yRange.first << ", " << yRange.second << "],\n";
    
    ss << "                    showgrid: " << (m_showGrid ? "true" : "false") << ",\n";
    ss << "                    gridcolor: '" << colorToString(m_gridColor) << "',\n";
    ss << "                    color: '" << colorToString(m_textColor) << "'\n";
    ss << "                },\n";
    
    ss << "                zaxis: {\n";
    ss << "                    title: '" << m_zLabel << "',\n";
    
    // Z range
    auto zRange = getZRange();
    ss << "                    range: [" << zRange.first << ", " << zRange.second << "],\n";
    
    ss << "                    showgrid: " << (m_showGrid ? "true" : "false") << ",\n";
    ss << "                    gridcolor: '" << colorToString(m_gridColor) << "',\n";
    ss << "                    color: '" << colorToString(m_textColor) << "'\n";
    ss << "                },\n";
    
    // Camera position based on view angles
    ss << "                camera: {\n";
    ss << "                    eye: {\n";
    ss << "                        x: " << std::cos(m_azimuth * 3.14159265358979323846 / 180.0) * std::cos(m_elevation * 3.14159265358979323846 / 180.0) << ",\n";
    ss << "                        y: " << std::sin(m_azimuth * 3.14159265358979323846 / 180.0) * std::cos(m_elevation * 3.14159265358979323846 / 180.0) << ",\n";
    ss << "                        z: " << std::sin(m_elevation * 3.14159265358979323846 / 180.0) << "\n";
    ss << "                    }\n";
    ss << "                },\n";
    
    ss << "                aspectmode: 'cube'\n";
    ss << "            },\n";
    
    // Legend
    ss << "            showlegend: " << (m_showLegend ? "true" : "false") << ",\n";
    ss << "            legend: {\n";
    ss << "                x: 1,\n";
    ss << "                y: 1,\n";
    ss << "                xanchor: 'right',\n";
    ss << "                yanchor: 'top',\n";
    ss << "                font: {\n";
    ss << "                    color: '" << colorToString(m_textColor) << "'\n";
    ss << "                }\n";
    ss << "            }\n";
    
    ss << "        };\n";
    
    // Create the plot
    ss << "        Plotly.newPlot('plot', data, layout);\n";
    ss << "    </script>\n";
    
    // HTML footer
    ss << "</body>\n";
    ss << "</html>\n";
    
    return ss.str();
}

namespace {
// Private helper function to save HTML to a file
bool saveHTMLToFile(const std::string& filename, const std::string& html) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    file << html;
    file.close();
    
    return true;
}
} // anonymous namespace

bool Plot3D::save(const std::string& filename, int width, int height) const {
    // Check file extension
    std::string extension;
    size_t dotPos = filename.find_last_of('.');
    if (dotPos != std::string::npos) {
        extension = filename.substr(dotPos + 1);
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    }

    std::string html = toHTML(width, height);
    
    if (extension == "html") {
        return saveHTMLToFile(filename, html);
    } else {
        // Default to HTML if extension is not recognized
        return saveHTMLToFile(filename + ".html", html);
    }
}

void Plot3D::show(int width, int height) const {
    // Save to a temporary file and open it in the default browser
    std::string tempFile = "temp_plot3d.html";
    std::string html = toHTML(width, height);
    
    if (saveHTMLToFile(tempFile, html)) {
        // Open the file in the default browser
        #ifdef _WIN32
        std::string command = "start " + tempFile;
        #elif __APPLE__
        std::string command = "open " + tempFile;
        #else
        std::string command = "xdg-open " + tempFile;
        #endif
        system(command.c_str());
    }
}


std::pair<double, double> Plot3D::calculateXRange() const {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    bool initialized = false;
    
    // Process surfaces
    for (const auto& surface : m_surfaces) {
        for (const auto& x : surface.getX()) {
            if (!initialized) {
                min = max = x;
                initialized = true;
            } else {
                min = std::min(min, x);
                max = std::max(max, x);
            }
        }
    }
    
    // Process scatter data
    for (const auto& scatter : m_scatters) {
        for (const auto& x : scatter.getX()) {
            if (!initialized) {
                min = max = x;
                initialized = true;
            } else {
                min = std::min(min, x);
                max = std::max(max, x);
            }
        }
    }
    
    // If no data, use default range
    if (!initialized) {
        min = -10.0;
        max = 10.0;
    }
    
    // Add a small margin
    double margin = (max - min) * 0.05;
    min -= margin;
    max += margin;
    
    return std::make_pair(min, max);
}

std::pair<double, double> Plot3D::calculateYRange() const {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    bool initialized = false;
    
    // Process surfaces
    for (const auto& surface : m_surfaces) {
        for (const auto& y : surface.getY()) {
            if (!initialized) {
                min = max = y;
                initialized = true;
            } else {
                min = std::min(min, y);
                max = std::max(max, y);
            }
        }
    }
    
    // Process scatter data
    for (const auto& scatter : m_scatters) {
        for (const auto& y : scatter.getY()) {
            if (!initialized) {
                min = max = y;
                initialized = true;
            } else {
                min = std::min(min, y);
                max = std::max(max, y);
            }
        }
    }
    
    // If no data, use default range
    if (!initialized) {
        min = -10.0;
        max = 10.0;
    }
    
    // Add a small margin
    double margin = (max - min) * 0.05;
    min -= margin;
    max += margin;
    
    return std::make_pair(min, max);
}

std::pair<double, double> Plot3D::calculateZRange() const {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    bool initialized = false;
    
    // Process surfaces
    for (const auto& surface : m_surfaces) {
        for (const auto& zRow : surface.getZ()) {
            for (const auto& z : zRow) {
                if (!initialized) {
                    min = max = z;
                    initialized = true;
                } else {
                    min = std::min(min, z);
                    max = std::max(max, z);
                }
            }
        }
    }
    
    // Process scatter data
    for (const auto& scatter : m_scatters) {
        for (const auto& z : scatter.getZ()) {
            if (!initialized) {
                min = max = z;
                initialized = true;
            } else {
                min = std::min(min, z);
                max = std::max(max, z);
            }
        }
    }
    
    // If no data, use default range
    if (!initialized) {
        min = -10.0;
        max = 10.0;
    }
    
    // Add a small margin
    double margin = (max - min) * 0.05;
    min -= margin;
    max += margin;
    
    return std::make_pair(min, max);
}
