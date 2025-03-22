#include <iostream>
#include <vector>
#include <cmath>
#include "../src/visualization/plot3d.h"
#include "../src/engineering/cad.h"

using namespace RebelCalc::Visualization;
using namespace RebelCalc::Engineering;

// Function to generate a 3D surface
double sinc(double x, double y) {
    double r = std::sqrt(x*x + y*y);
    if (r < 1e-10) {
        return 1.0;
    }
    return std::sin(r) / r;
}

// Function to generate a 3D spiral
void generateSpiral(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, int points) {
    x.resize(points);
    y.resize(points);
    z.resize(points);
    
    for (int i = 0; i < points; ++i) {
        double t = 10.0 * static_cast<double>(i) / points;
        x[i] = std::sin(t) * t;
        y[i] = std::cos(t) * t;
        z[i] = t;
    }
}

int main() {
    std::cout << "RebelCALC Plot3D Demo" << std::endl;
    std::cout << "=====================" << std::endl;
    
    // Create a 3D plot
    Plot3D plot("3D Surface Plot Demo");
    
    // Set labels and other properties
    plot.setXLabel("X Axis");
    plot.setYLabel("Y Axis");
    plot.setZLabel("Z Axis");
    plot.setBackgroundColor(Color(240, 240, 240)); // Light gray
    plot.setGridColor(Color(200, 200, 200));
    plot.setTextColor(Color(50, 50, 50));
    plot.setViewAngles(30.0, 30.0);
    
    // Add a surface using a function
    plot.addSurface(sinc, -8.0, 8.0, -8.0, 8.0, 50, 50, "Sinc Function");
    
    // Get the surface and set its style to wireframe
    DataSeries3D& surface = const_cast<DataSeries3D&>(plot.getSurface(0));
    surface.setSurfaceStyle(SurfaceStyle::WIREFRAME);
    surface.setWireframeColor(Color(0, 0, 255)); // Blue
    surface.setWireframeWidth(2.0);
    
    // Create a color map for the surface
    std::vector<Color> colorMap = {
        Color(0, 0, 255),   // Blue
        Color(0, 255, 255), // Cyan
        Color(0, 255, 0),   // Green
        Color(255, 255, 0), // Yellow
        Color(255, 0, 0)    // Red
    };
    surface.setColorMap(colorMap);
    
    // Add a scatter plot (spiral)
    std::vector<double> spiralX, spiralY, spiralZ;
    generateSpiral(spiralX, spiralY, spiralZ, 100);
    plot.addScatter(spiralX, spiralY, spiralZ, "Spiral");
    
    // Get the scatter and set its properties
    ScatterSeries3D& scatter = const_cast<ScatterSeries3D&>(plot.getScatter(0));
    scatter.setMarkerStyle(MarkerStyle::CIRCLE);
    scatter.setMarkerSize(5.0);
    scatter.setMarkerColor(Color(255, 0, 0)); // Red
    scatter.setColorByZ(true); // Color points by Z value
    scatter.setColorMap(colorMap);
    
    // Save the plot to an HTML file
    std::string filename = "plot3d_demo.html";
    if (plot.save(filename)) {
        std::cout << "Plot saved to " << filename << std::endl;
    } else {
        std::cerr << "Failed to save plot to " << filename << std::endl;
    }
    
    // Show the plot in a browser
    plot.show();
    
    // Create another plot with different surface styles
    Plot3D stylePlot("Surface Style Demo");
    
    // Add surfaces with different styles
    size_t solidIdx = stylePlot.addSurface(sinc, -8.0, 8.0, -8.0, 8.0, 30, 30, "Solid");
    size_t wireframeIdx = stylePlot.addSurface(sinc, -8.0, 8.0, -8.0, 8.0, 30, 30, "Wireframe");
    size_t pointsIdx = stylePlot.addSurface(sinc, -8.0, 8.0, -8.0, 8.0, 30, 30, "Points");
    size_t contourIdx = stylePlot.addSurface(sinc, -8.0, 8.0, -8.0, 8.0, 30, 30, "Contour");
    
    // Set different styles
    const_cast<DataSeries3D&>(stylePlot.getSurface(wireframeIdx)).setSurfaceStyle(SurfaceStyle::WIREFRAME);
    const_cast<DataSeries3D&>(stylePlot.getSurface(pointsIdx)).setSurfaceStyle(SurfaceStyle::POINTS);
    const_cast<DataSeries3D&>(stylePlot.getSurface(contourIdx)).setSurfaceStyle(SurfaceStyle::CONTOUR);
    
    // Set contour levels
    std::vector<double> contourLevels;
    for (double level = -0.2; level <= 1.0; level += 0.1) {
        contourLevels.push_back(level);
    }
    const_cast<DataSeries3D&>(stylePlot.getSurface(contourIdx)).setContourLevels(contourLevels);
    
    // Save the style demo plot
    std::string styleFilename = "surface_style_demo.html";
    if (stylePlot.save(styleFilename)) {
        std::cout << "Style demo plot saved to " << styleFilename << std::endl;
    } else {
        std::cerr << "Failed to save style demo plot to " << styleFilename << std::endl;
    }
    
    return 0;
}
