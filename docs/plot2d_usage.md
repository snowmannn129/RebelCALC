# Plot2D Usage Guide

The Plot2D class in RebelCALC provides powerful 2D plotting capabilities with interactive features. This guide explains how to use the Plot2D class to create interactive 2D plots.

## Basic Usage

```cpp
#include "../src/visualization/plot2d.h"
using namespace RebelCalc::Visualization;

// Create a plot
Plot2D plot("My Plot");

// Set axis labels
plot.setXLabel("X Axis");
plot.setYLabel("Y Axis");

// Add a data series from a function
plot.addSeries([](double x) { return std::sin(x); }, 0.0, 2.0 * M_PI, 100, "Sine Wave");

// Add a data series from data points
std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
std::vector<double> y = {2.0, 4.0, 6.0, 8.0, 10.0};
plot.addSeries(x, y, "Linear Data");

// Show the plot in a browser
plot.show();

// Save the plot to an HTML file
plot.save("my_plot.html");
```

## Customizing Plot Appearance

The Plot2D class provides various methods to customize the appearance of the plot:

```cpp
// Set background color
plot.setBackgroundColor(Color(240, 240, 240)); // Light gray

// Set grid color
plot.setGridColor(Color(200, 200, 200));

// Set text color
plot.setTextColor(Color(50, 50, 50));

// Show/hide grid
plot.setShowGrid(true);

// Show/hide legend
plot.setShowLegend(true);

// Set axis ranges
plot.setXRange(0.0, 10.0);
plot.setYRange(-1.0, 1.0);
```

## Customizing Data Series

You can customize the appearance of individual data series:

```cpp
// Add a data series
size_t seriesIndex = plot.addSeries(x, y, "My Series");

// Get the series and customize it
DataSeries& series = const_cast<DataSeries&>(plot.getSeries(seriesIndex));

// Set line properties
series.setLineColor(Color::Blue());
series.setLineWidth(2.0);
series.setLineStyle(LineStyle::DASHED); // SOLID, DASHED, DOTTED, DASH_DOT

// Set marker properties
series.setMarkerStyle(MarkerStyle::CIRCLE); // NONE, CIRCLE, SQUARE, TRIANGLE, DIAMOND, CROSS, PLUS
series.setMarkerSize(5.0);
series.setMarkerColor(Color::Red());
```

## Interactive Features

The Plot2D class generates interactive HTML plots using Plotly.js. These plots provide the following interactive features:

### Zooming and Panning

- **Zoom**: Click and drag to create a zoom box, or use the mouse wheel
- **Pan**: Click and drag in pan mode (toggle using the mode bar)
- **Reset**: Double-click to reset the view

### Data Exploration

- **Tooltips**: Hover over data points to see their values
- **Click Events**: Click on data points to trigger events (customizable)

### Plot Controls

The interactive plot includes a mode bar with the following controls:

- **Zoom**: Zoom in on a selected region
- **Pan**: Pan the plot
- **Box Select**: Select data points within a box
- **Lasso Select**: Select data points with a lasso tool
- **Zoom In/Out**: Zoom in or out on the plot
- **Reset Axes**: Reset the plot to its original view
- **Toggle Spike Lines**: Show spike lines for data points
- **Download Plot as PNG**: Save the plot as a PNG image

## Example

Here's a complete example that demonstrates various features of the Plot2D class:

```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include "../src/visualization/plot2d.h"

using namespace RebelCalc::Visualization;

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    // Create a 2D plot
    Plot2D plot("Interactive 2D Plot Demo");
    
    // Set labels and other properties
    plot.setXLabel("X Axis");
    plot.setYLabel("Y Axis");
    plot.setBackgroundColor(Color(240, 240, 240)); // Light gray
    plot.setGridColor(Color(200, 200, 200));
    plot.setTextColor(Color(50, 50, 50));
    
    // Add a sine wave function
    size_t sineIdx = plot.addSeries([](double x) { return std::sin(x); }, 
                                   0.0, 4.0 * M_PI, 100, "Sine Wave");
    
    // Get the series and set its properties
    DataSeries& sineSeries = const_cast<DataSeries&>(plot.getSeries(sineIdx));
    sineSeries.setLineColor(Color::Blue());
    sineSeries.setLineWidth(2.0);
    sineSeries.setMarkerStyle(MarkerStyle::CIRCLE);
    sineSeries.setMarkerSize(5.0);
    sineSeries.setMarkerColor(Color::Red());
    
    // Add a cosine wave function
    size_t cosineIdx = plot.addSeries([](double x) { return std::cos(x); }, 
                                     0.0, 4.0 * M_PI, 100, "Cosine Wave");
    
    // Get the series and set its properties
    DataSeries& cosineSeries = const_cast<DataSeries&>(plot.getSeries(cosineIdx));
    cosineSeries.setLineColor(Color::Green());
    cosineSeries.setLineWidth(2.0);
    cosineSeries.setLineStyle(LineStyle::DASHED);
    
    // Save the plot to an HTML file
    plot.save("interactive_plot_demo.html");
    
    // Show the plot in a browser
    plot.show();
    
    return 0;
}
```

## Output Formats

The Plot2D class supports the following output formats:

- **HTML**: Interactive plots using Plotly.js
- **SVG**: Static vector graphics

To save a plot in a specific format, use the `save` method with the appropriate file extension:

```cpp
// Save as HTML (interactive)
plot.save("my_plot.html");

// Save as SVG (static)
plot.save("my_plot.svg");
```

## Advanced Usage

### Custom Colors

You can create custom colors using the Color class:

```cpp
// Create a custom color (RGB)
Color customColor(128, 64, 192); // Purple-ish

// Create a custom color with transparency (RGBA)
Color transparentColor(255, 0, 0, 128); // Semi-transparent red
```

### Multiple Plots

You can create multiple plots in the same application:

```cpp
// Create two plots
Plot2D plot1("First Plot");
Plot2D plot2("Second Plot");

// Add data to each plot
plot1.addSeries(x1, y1, "Series 1");
plot2.addSeries(x2, y2, "Series 2");

// Show both plots
plot1.show();
plot2.show();
```

## Conclusion

The Plot2D class provides a powerful and flexible way to create interactive 2D plots in RebelCALC. With its extensive customization options and interactive features, it's suitable for a wide range of visualization needs.
