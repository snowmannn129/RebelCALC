#include <iostream>
#include <vector>
#include <cmath>
#include "../src/visualization/plot2d.h"

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace RebelCalc::Visualization;

// Function to generate a sine wave
double sine(double x) {
    return std::sin(x);
}

// Function to generate a cosine wave
double cosine(double x) {
    return std::cos(x);
}

// Function to generate a damped sine wave
double dampedSine(double x) {
    return std::exp(-0.2 * x) * std::sin(x);
}

int main() {
    std::cout << "RebelCALC Plot2D Demo" << std::endl;
    std::cout << "=====================" << std::endl;
    
    // Create a 2D plot
    Plot2D plot("Interactive 2D Plot Demo");
    
    // Set labels and other properties
    plot.setXLabel("X Axis");
    plot.setYLabel("Y Axis");
    plot.setBackgroundColor(Color(240, 240, 240)); // Light gray
    plot.setGridColor(Color(200, 200, 200));
    plot.setTextColor(Color(50, 50, 50));
    
    // Add a sine wave function
    size_t sineIdx = plot.addSeries(sine, 0.0, 4.0 * M_PI, 100, "Sine Wave");
    
    // Get the series and set its properties
    DataSeries& sineSeries = const_cast<DataSeries&>(plot.getSeries(sineIdx));
    sineSeries.setLineColor(Color::Blue());
    sineSeries.setLineWidth(2.0);
    sineSeries.setMarkerStyle(MarkerStyle::CIRCLE);
    sineSeries.setMarkerSize(5.0);
    sineSeries.setMarkerColor(Color::Red());
    
    // Add a cosine wave function
    size_t cosineIdx = plot.addSeries(cosine, 0.0, 4.0 * M_PI, 100, "Cosine Wave");
    
    // Get the series and set its properties
    DataSeries& cosineSeries = const_cast<DataSeries&>(plot.getSeries(cosineIdx));
    cosineSeries.setLineColor(Color::Green());
    cosineSeries.setLineWidth(2.0);
    cosineSeries.setLineStyle(LineStyle::DASHED);
    
    // Add a damped sine wave function
    size_t dampedIdx = plot.addSeries(dampedSine, 0.0, 4.0 * M_PI, 100, "Damped Sine Wave");
    
    // Get the series and set its properties
    DataSeries& dampedSeries = const_cast<DataSeries&>(plot.getSeries(dampedIdx));
    dampedSeries.setLineColor(Color::Purple());
    dampedSeries.setLineWidth(2.0);
    dampedSeries.setLineStyle(LineStyle::DASH_DOT);
    dampedSeries.setMarkerStyle(MarkerStyle::DIAMOND);
    dampedSeries.setMarkerSize(6.0);
    dampedSeries.setMarkerColor(Color::Orange());
    
    // Save the plot to an HTML file
    std::string filename = "plot2d_demo.html";
    if (plot.save(filename)) {
        std::cout << "Plot saved to " << filename << std::endl;
    } else {
        std::cerr << "Failed to save plot to " << filename << std::endl;
    }
    
    // Show the plot in a browser
    plot.show();
    
    // Create another plot with different line styles
    Plot2D stylePlot("Line Style Demo");
    stylePlot.setXLabel("X Axis");
    stylePlot.setYLabel("Y Axis");
    
    // Add series with different line styles
    std::vector<double> x;
    for (int i = 0; i < 100; ++i) {
        x.push_back(i * 0.1);
    }
    
    std::vector<double> y1, y2, y3, y4;
    for (double val : x) {
        y1.push_back(std::sin(val));
        y2.push_back(std::sin(val) + 1.0);
        y3.push_back(std::sin(val) + 2.0);
        y4.push_back(std::sin(val) + 3.0);
    }
    
    size_t solidIdx = stylePlot.addSeries(x, y1, "Solid Line");
    size_t dashedIdx = stylePlot.addSeries(x, y2, "Dashed Line");
    size_t dottedIdx = stylePlot.addSeries(x, y3, "Dotted Line");
    size_t dashDotIdx = stylePlot.addSeries(x, y4, "Dash-Dot Line");
    
    // Set different line styles
    const_cast<DataSeries&>(stylePlot.getSeries(solidIdx)).setLineStyle(LineStyle::SOLID);
    const_cast<DataSeries&>(stylePlot.getSeries(dashedIdx)).setLineStyle(LineStyle::DASHED);
    const_cast<DataSeries&>(stylePlot.getSeries(dottedIdx)).setLineStyle(LineStyle::DOTTED);
    const_cast<DataSeries&>(stylePlot.getSeries(dashDotIdx)).setLineStyle(LineStyle::DASH_DOT);
    
    // Set different colors
    const_cast<DataSeries&>(stylePlot.getSeries(solidIdx)).setLineColor(Color::Blue());
    const_cast<DataSeries&>(stylePlot.getSeries(dashedIdx)).setLineColor(Color::Green());
    const_cast<DataSeries&>(stylePlot.getSeries(dottedIdx)).setLineColor(Color::Red());
    const_cast<DataSeries&>(stylePlot.getSeries(dashDotIdx)).setLineColor(Color::Purple());
    
    // Save the style demo plot
    std::string styleFilename = "line_style_demo.html";
    if (stylePlot.save(styleFilename)) {
        std::cout << "Style demo plot saved to " << styleFilename << std::endl;
    } else {
        std::cerr << "Failed to save style demo plot to " << styleFilename << std::endl;
    }
    
    // Create a marker style demo
    Plot2D markerPlot("Marker Style Demo");
    markerPlot.setXLabel("X Axis");
    markerPlot.setYLabel("Y Axis");
    
    // Add series with different marker styles
    std::vector<double> markerX;
    for (int i = 0; i < 10; ++i) {
        markerX.push_back(i);
    }
    
    std::vector<double> markerY1, markerY2, markerY3, markerY4, markerY5;
    for (double val : markerX) {
        markerY1.push_back(std::sin(val * 0.5));
        markerY2.push_back(std::sin(val * 0.5) + 1.0);
        markerY3.push_back(std::sin(val * 0.5) + 2.0);
        markerY4.push_back(std::sin(val * 0.5) + 3.0);
        markerY5.push_back(std::sin(val * 0.5) + 4.0);
    }
    
    size_t circleIdx = markerPlot.addSeries(markerX, markerY1, "Circle Markers");
    size_t squareIdx = markerPlot.addSeries(markerX, markerY2, "Square Markers");
    size_t triangleIdx = markerPlot.addSeries(markerX, markerY3, "Triangle Markers");
    size_t diamondIdx = markerPlot.addSeries(markerX, markerY4, "Diamond Markers");
    size_t crossIdx = markerPlot.addSeries(markerX, markerY5, "Cross Markers");
    
    // Set marker styles
    const_cast<DataSeries&>(markerPlot.getSeries(circleIdx)).setMarkerStyle(MarkerStyle::CIRCLE);
    const_cast<DataSeries&>(markerPlot.getSeries(squareIdx)).setMarkerStyle(MarkerStyle::SQUARE);
    const_cast<DataSeries&>(markerPlot.getSeries(triangleIdx)).setMarkerStyle(MarkerStyle::TRIANGLE);
    const_cast<DataSeries&>(markerPlot.getSeries(diamondIdx)).setMarkerStyle(MarkerStyle::DIAMOND);
    const_cast<DataSeries&>(markerPlot.getSeries(crossIdx)).setMarkerStyle(MarkerStyle::CROSS);
    
    // Set marker colors
    const_cast<DataSeries&>(markerPlot.getSeries(circleIdx)).setMarkerColor(Color::Blue());
    const_cast<DataSeries&>(markerPlot.getSeries(squareIdx)).setMarkerColor(Color::Green());
    const_cast<DataSeries&>(markerPlot.getSeries(triangleIdx)).setMarkerColor(Color::Red());
    const_cast<DataSeries&>(markerPlot.getSeries(diamondIdx)).setMarkerColor(Color::Purple());
    const_cast<DataSeries&>(markerPlot.getSeries(crossIdx)).setMarkerColor(Color::Orange());
    
    // Set marker sizes
    const_cast<DataSeries&>(markerPlot.getSeries(circleIdx)).setMarkerSize(8.0);
    const_cast<DataSeries&>(markerPlot.getSeries(squareIdx)).setMarkerSize(8.0);
    const_cast<DataSeries&>(markerPlot.getSeries(triangleIdx)).setMarkerSize(8.0);
    const_cast<DataSeries&>(markerPlot.getSeries(diamondIdx)).setMarkerSize(8.0);
    const_cast<DataSeries&>(markerPlot.getSeries(crossIdx)).setMarkerSize(8.0);
    
    // Set line styles to NONE to show only markers
    const_cast<DataSeries&>(markerPlot.getSeries(circleIdx)).setLineStyle(LineStyle::DOTTED);
    const_cast<DataSeries&>(markerPlot.getSeries(squareIdx)).setLineStyle(LineStyle::DOTTED);
    const_cast<DataSeries&>(markerPlot.getSeries(triangleIdx)).setLineStyle(LineStyle::DOTTED);
    const_cast<DataSeries&>(markerPlot.getSeries(diamondIdx)).setLineStyle(LineStyle::DOTTED);
    const_cast<DataSeries&>(markerPlot.getSeries(crossIdx)).setLineStyle(LineStyle::DOTTED);
    
    // Save the marker demo plot
    std::string markerFilename = "marker_style_demo.html";
    if (markerPlot.save(markerFilename)) {
        std::cout << "Marker demo plot saved to " << markerFilename << std::endl;
    } else {
        std::cerr << "Failed to save marker demo plot to " << markerFilename << std::endl;
    }
    
    return 0;
}
