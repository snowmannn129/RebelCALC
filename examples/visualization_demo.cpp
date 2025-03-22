#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <memory>

#include "../src/visualization/plot2d.h"

using namespace RebelCalc::Visualization;

// Helper function to print a section header
void printHeader(const std::string& header) {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "  " << header << std::endl;
    std::cout << std::string(80, '=') << std::endl;
}

// Demo for 2D plotting
void demo2DPlotting() {
    printHeader("2D Plotting Demonstration");
    
    // Create a 2D plot
    Plot2D plot("Example Plot");
    plot.setXLabel("X Axis");
    plot.setYLabel("Y Axis");
    
    // Create and add a sine wave series
    std::vector<double> x1, y1;
    for (double x = 0; x <= 10; x += 0.1) {
        x1.push_back(x);
        y1.push_back(std::sin(x));
    }
    
    DataSeries sineSeries(x1, y1, "Sine Wave");
    sineSeries.setLineColor(Color::Blue());
    sineSeries.setMarkerStyle(MarkerStyle::CIRCLE);
    sineSeries.setMarkerSize(3.0);
    plot.addSeries(sineSeries);
    
    // Create and add a cosine wave series
    std::vector<double> x2, y2;
    for (double x = 0; x <= 10; x += 0.1) {
        x2.push_back(x);
        y2.push_back(std::cos(x));
    }
    
    DataSeries cosineSeries(x2, y2, "Cosine Wave");
    cosineSeries.setLineColor(Color::Red());
    cosineSeries.setLineStyle(LineStyle::DASHED);
    cosineSeries.setMarkerStyle(MarkerStyle::SQUARE);
    cosineSeries.setMarkerSize(3.0);
    plot.addSeries(cosineSeries);
    
    // Create and add a tangent wave series
    std::vector<double> x3, y3;
    for (double x = 0; x <= 10; x += 0.1) {
        // Skip points where tangent is undefined or very large
        if (std::abs(std::cos(x)) < 0.1) continue;
        
        x3.push_back(x);
        y3.push_back(std::tan(x));
    }
    
    DataSeries tanSeries(x3, y3, "Tangent Wave");
    tanSeries.setLineColor(Color::Green());
    tanSeries.setLineStyle(LineStyle::DOTTED);
    tanSeries.setMarkerStyle(MarkerStyle::TRIANGLE);
    tanSeries.setMarkerSize(3.0);
    plot.addSeries(tanSeries);
    
    // Set the y-axis range to better visualize the data
    plot.setYRange(-2.0, 2.0);
    
    // Save the plot to an SVG file
    std::string filename = "trigonometric_functions.svg";
    if (plot.save(filename)) {
        std::cout << "Plot saved to " << filename << std::endl;
    } else {
        std::cout << "Failed to save plot to " << filename << std::endl;
    }
    
    // Save the plot to an HTML file
    filename = "trigonometric_functions.html";
    if (plot.save(filename)) {
        std::cout << "Plot saved to " << filename << std::endl;
    } else {
        std::cout << "Failed to save plot to " << filename << std::endl;
    }
    
    // Generate and print the SVG code
    std::cout << "\nSVG representation of the plot:" << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << plot.toSVG(400, 300).substr(0, 500) << "..." << std::endl;
    
    // Create a plot with a function
    Plot2D functionPlot("Function Plot");
    functionPlot.setXLabel("X");
    functionPlot.setYLabel("Y");
    
    // Create a function series
    // First, generate the data points
    std::vector<double> xFunc, yFunc;
    for (double x = -5.0; x <= 5.0; x += 0.1) {
        xFunc.push_back(x);
        yFunc.push_back(x * x); // Quadratic function
    }
    
    // Create and configure the series
    DataSeries quadraticSeries(xFunc, yFunc, "Quadratic Function");
    quadraticSeries.setLineColor(Color::Purple());
    quadraticSeries.setLineWidth(2.0);
    
    // Add the series to the plot
    functionPlot.addSeries(quadraticSeries);
    
    // Save the function plot to an SVG file
    filename = "quadratic_function.svg";
    if (functionPlot.save(filename)) {
        std::cout << "\nFunction plot saved to " << filename << std::endl;
    } else {
        std::cout << "\nFailed to save function plot to " << filename << std::endl;
    }
}

int main() {
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "                      RebelCALC Visualization Demo" << std::endl;
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "\nThis program demonstrates the visualization capabilities of RebelCALC." << std::endl;
    
    demo2DPlotting();
    
    std::cout << "\n" << std::string(80, '*') << std::endl;
    std::cout << "                      End of Visualization Demo" << std::endl;
    std::cout << std::string(80, '*') << std::endl;
    
    return 0;
}
