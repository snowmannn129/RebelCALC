# RebelCALC Engineering Modules

This document provides an overview of the engineering modules in RebelCALC, which extend the calculator's capabilities to handle engineering-specific calculations.

## Table of Contents

1. [Introduction](#introduction)
2. [CAD Module](#cad-module)
3. [Physics Module](#physics-module)
4. [Electrical Module](#electrical-module)
5. [Visualization Module](#visualization-module)
6. [Solvers Module](#solvers-module)
7. [Integration with RebelCALC](#integration-with-rebelcalc)
8. [Examples](#examples)

## Introduction

The RebelCALC engineering modules provide specialized functionality for various engineering disciplines, including:

- **CAD Module**: 3D geometry primitives and operations for computer-aided design
- **Physics Module**: Mechanics, dynamics, and physical simulations
- **Electrical Module**: Circuit analysis and electrical engineering calculations

These modules are designed to work seamlessly with the core RebelCALC functionality, providing engineers with powerful tools for their calculations.

## CAD Module

The CAD module provides 3D geometry primitives and operations for computer-aided design.

### Key Features

- **3D Geometry Primitives**:
  - Point3D: 3D point/vector representation with vector operations
  - LineSegment: Line segment between two points
  - Plane: Infinite plane defined by a point and normal
  - Circle: Circle defined by center, normal, and radius
  - Triangle: Triangle defined by three points
  - Polygon: Polygon defined by a set of points

- **Geometric Calculations**:
  - Distance calculations (point-to-point, point-to-line, point-to-plane)
  - Area calculations (triangle, circle, polygon)
  - Volume calculations (various 3D shapes)
  - Perimeter calculations
  - Normal vector calculations

- **Intersection Detection**:
  - Line-plane intersection
  - Line-triangle intersection
  - Ray-triangle intersection
  - Line-line intersection

- **Spatial Operations**:
  - Bounding box calculation
  - Bounding sphere calculation
  - Convex hull generation
  - Point containment tests
  - Projection operations

- **Advanced Operations**:
  - Plane fitting to points
  - Centroid calculation
  - Inertia tensor calculation
  - Principal axes determination

### Usage Example

```cpp
#include "engineering/cad.h"

using namespace RebelCalc::Engineering;

// Create points
CAD::Point3D p1(1.0, 2.0, 3.0);
CAD::Point3D p2(4.0, 5.0, 6.0);
CAD::Point3D p3(7.0, 8.0, 9.0);

// Create a triangle
CAD::Triangle triangle(p1, p2, p3);

// Calculate area and normal
double area = triangle.area();
CAD::Point3D normal = triangle.normal();

// Test if a point is inside the triangle
CAD::Point3D testPoint = (p1 + p2 + p3) / 3.0; // Centroid
bool inside = triangle.containsPoint(testPoint);
```

## Physics Module

The Physics module provides mechanics, dynamics, and physical simulations.

### Key Features

- **Rigid Body Dynamics**:
  - Position, orientation, velocity tracking
  - Force and torque application
  - Moment of inertia calculations
  - Angular momentum calculations

- **Particle Systems**:
  - Particle creation and management
  - Collision detection and resolution
  - Force application (gravity, springs, etc.)
  - Time-step simulation

- **Spring Systems**:
  - Spring force calculations
  - Damping and constraints
  - Spring-mass networks
  - Oscillation analysis

- **Fluid Dynamics**:
  - Smoothed Particle Hydrodynamics (SPH)
  - Pressure and density calculations
  - Viscosity effects
  - Surface tension

- **Force Calculations**:
  - Gravitational forces
  - Electric forces
  - Magnetic forces
  - Drag forces
  - Spring forces
  - Damping forces

- **Motion Analysis**:
  - Projectile motion
  - Orbital mechanics
  - Pendulum motion
  - Harmonic oscillation

- **Energy and Momentum**:
  - Kinetic energy calculation
  - Potential energy calculation
  - Linear momentum calculation
  - Angular momentum calculation
  - Conservation of energy and momentum

- **Collision Response**:
  - Impulse-based physics
  - Coefficient of restitution
  - Friction modeling
  - Contact point determination

### Usage Example

```cpp
#include "engineering/physics.h"

using namespace RebelCalc::Engineering;

// Projectile motion
CAD::Point3D initialPosition(0.0, 0.0, 0.0);
CAD::Point3D initialVelocity(10.0, 15.0, 0.0);
CAD::Point3D gravity(0.0, -9.81, 0.0);

// Calculate position at time t
double t = 1.5; // seconds
CAD::Point3D position = Physics::projectileTrajectory(initialPosition, initialVelocity, gravity, t);

// Calculate key parameters
double timeOfFlight = Physics::projectileTimeOfFlight(initialPosition.y, initialVelocity.y);
double range = Physics::projectileRange(initialPosition.y, initialVelocity);
double maxHeight = Physics::projectileMaxHeight(initialPosition.y, initialVelocity.y);

// Energy calculations
double mass = 1.0; // kg
double kineticEnergy = Physics::kineticEnergy(mass, initialVelocity);
double potentialEnergy = Physics::potentialEnergy(mass, maxHeight);
```

## Electrical Module

The Electrical module provides circuit analysis and electrical engineering calculations.

### Key Features

- **Circuit Analysis**:
  - DC circuit analysis
  - AC circuit analysis
  - Nodal analysis
  - Mesh analysis
  - Thevenin and Norton equivalents

- **Digital Circuit Simulation**:
  - Logic gates (AND, OR, NOT, etc.)
  - Flip-flops and latches
  - Counters and registers
  - Combinational logic
  - Sequential logic

- **Power System Analysis**:
  - Power calculations (real, reactive, apparent)
  - Power factor correction
  - Three-phase systems
  - Transformer analysis
  - Transmission line modeling

- **Component Calculations**:
  - Resistance calculations
  - Capacitance calculations
  - Inductance calculations
  - Impedance calculations
  - Reactance calculations

- **Filter Design and Analysis**:
  - Low-pass filters
  - High-pass filters
  - Band-pass filters
  - Band-stop filters
  - Filter response analysis

- **Signal Analysis**:
  - Gain calculations
  - Decibel conversions
  - Frequency response
  - Transfer functions
  - Bode plots

- **Circuit Theorems**:
  - Thevenin's theorem
  - Norton's theorem
  - Maximum power transfer theorem
  - Superposition theorem
  - Reciprocity theorem

### Usage Example

```cpp
#include "engineering/electrical.h"

using namespace RebelCalc::Engineering;

// Component values
double resistance = 100.0; // 100 ohms
double capacitance = 1.0e-6; // 1 uF
double inductance = 0.1; // 0.1 H
double frequency = 1000.0; // 1 kHz

// Calculate impedances
Electrical::Complex zR = Electrical::impedanceResistor(resistance);
Electrical::Complex zC = Electrical::impedanceCapacitor(capacitance, frequency);
Electrical::Complex zL = Electrical::impedanceInductor(inductance, frequency);

// Calculate series and parallel RLC impedances
Electrical::Complex zSeries = Electrical::impedanceSeriesRLC(resistance, inductance, capacitance, frequency);
Electrical::Complex zParallel = Electrical::impedanceParallelRLC(resistance, inductance, capacitance, frequency);

// Calculate resonant frequency
double resonantFreq = Electrical::resonantFrequency(inductance, capacitance);
```

## Visualization Module

The Visualization module provides tools for creating 2D and 3D plots and visualizations of data.

### Key Features

- **2D Plotting**:
  - Line plots with customizable appearance
  - Scatter plots with various marker styles
  - Bar charts and histograms
  - Pie charts and polar plots
  - Error bars and confidence intervals
  - Multiple data series in a single plot
  - Customizable line styles, colors, and markers
  - Axis labels, titles, and legends
  - Grid lines and tick marks
  - Logarithmic and linear scales
  - Export to SVG and HTML formats

- **3D Plotting**:
  - Surface plots for 3D functions
  - Scatter plots in 3D space
  - Wireframe and mesh plots
  - Contour plots and heat maps
  - Vector field visualization
  - Isosurface rendering
  - Interactive 3D rotation and zooming
  - Customizable lighting and shading
  - Export to interactive HTML formats

- **Data Visualization**:
  - Box plots and violin plots
  - Heatmaps and correlation matrices
  - Dendrograms and cluster visualization
  - Network graphs and tree structures
  - Geographic maps and spatial data
  - Time series visualization
  - Interactive tooltips and data exploration
  - Customizable color maps and palettes

- **Plot Customization**:
  - Color management with predefined and custom colors
  - Line style management (solid, dashed, dotted, etc.)
  - Marker style management (circle, square, triangle, etc.)
  - Font selection and text formatting
  - Layout management for multiple plots
  - Annotation and text labels
  - Legend positioning and formatting
  - Export options for various formats

### Usage Example

```cpp
#include "visualization/plot2d.h"

using namespace RebelCalc::Visualization;

// Create a 2D plot
Plot2D plot("Example Plot");
plot.setXLabel("X Axis");
plot.setYLabel("Y Axis");

// Create data for a sine wave
std::vector<double> x, y;
for (double t = 0; t <= 10; t += 0.1) {
    x.push_back(t);
    y.push_back(std::sin(t));
}

// Create and add a data series
DataSeries sineSeries(x, y, "Sine Wave");
sineSeries.setLineColor(Color::Blue());
sineSeries.setMarkerStyle(MarkerStyle::CIRCLE);
sineSeries.setMarkerSize(3.0);
plot.addSeries(sineSeries);

// Set the y-axis range
plot.setYRange(-1.5, 1.5);

// Save the plot to an SVG file
plot.save("sine_wave.svg");

// Create a function plot
Plot2D functionPlot("Function Plot");
functionPlot.setXLabel("X");
functionPlot.setYLabel("Y");

// Create and add a function series
std::vector<double> xFunc, yFunc;
for (double x = -5.0; x <= 5.0; x += 0.1) {
    xFunc.push_back(x);
    yFunc.push_back(x * x); // Quadratic function
}

DataSeries quadraticSeries(xFunc, yFunc, "Quadratic Function");
quadraticSeries.setLineColor(Color::Purple());
quadraticSeries.setLineWidth(2.0);
functionPlot.addSeries(quadraticSeries);

// Save the function plot
functionPlot.save("quadratic_function.svg");
```

## Solvers Module

The Solvers module provides specialized numerical solvers for various engineering problems.

### Key Features

- **Finite Element Method (FEM) Solver**:
  - Linear and nonlinear structural analysis
  - Static and dynamic analysis
  - Thermal analysis
  - Mesh generation and refinement
  - Element types (beam, truss, shell, solid)
  - Material models (linear elastic, nonlinear, anisotropic)
  - Boundary conditions and constraints
  - Load cases and combinations
  - Result visualization and post-processing

- **Computational Fluid Dynamics (CFD) Solver**:
  - Incompressible and compressible flow
  - Laminar and turbulent flow
  - Steady-state and transient analysis
  - Heat transfer in fluids
  - Multiphase flow
  - Boundary layer modeling
  - Navier-Stokes equations solver
  - Finite volume discretization
  - Result visualization and streamlines

- **Circuit Simulation Solver**:
  - DC analysis
  - AC analysis
  - Transient analysis
  - Frequency response analysis
  - Noise analysis
  - Distortion analysis
  - Monte Carlo analysis
  - Worst-case analysis
  - Sensitivity analysis
  - Component models (resistors, capacitors, inductors, diodes, transistors)
  - Nonlinear device modeling
  - Circuit topology analysis

- **Optimization Solver**:
  - Unconstrained optimization
  - Constrained optimization
  - Linear programming
  - Quadratic programming
  - Nonlinear programming
  - Multi-objective optimization
  - Gradient-based methods (Gradient Descent, Newton, BFGS, Conjugate Gradient)
  - Direct search methods (Simplex/Nelder-Mead)
  - Global optimization methods (Simulated Annealing, Genetic Algorithm, Particle Swarm)
  - Constraint handling techniques
  - Sensitivity analysis
  - Convergence analysis

### Usage Example

```cpp
#include "solvers/optimization.h"

using namespace RebelCalc::Solvers;

// Define the Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2
auto rosenbrock = [](const std::vector<double>& x) -> double {
    return std::pow(1.0 - x[0], 2) + 100.0 * std::pow(x[1] - x[0] * x[0], 2);
};

// Define the gradient of the Rosenbrock function
auto rosenbrockGradient = [](const std::vector<double>& x) -> std::vector<double> {
    std::vector<double> grad(2);
    grad[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
    grad[1] = 200.0 * (x[1] - x[0] * x[0]);
    return grad;
};

// Create an optimization problem
OptimizationProblem problem(rosenbrock, 2);
problem.setGradient(rosenbrockGradient);

// Create optimization options
OptimizationOptions options;
options.setMaxIterations(1000);
options.setTolerance(1e-6);
options.setStepSize(0.01);

// Create an optimization solver
OptimizationSolver solver(problem, OptimizationAlgorithm::BFGS, options);

// Set the initial point
std::vector<double> initialPoint = {-1.2, 1.0};
solver.setInitialPoint(initialPoint);

// Solve the problem
OptimizationResult result = solver.solve();

// Check if the optimization converged
if (result.hasConverged()) {
    // Get the solution
    std::vector<double> solution = result.getSolution();
    double objectiveValue = result.getObjectiveValue();
    
    // Print the results
    std::cout << "Solution: [" << solution[0] << ", " << solution[1] << "]" << std::endl;
    std::cout << "Objective value: " << objectiveValue << std::endl;
    std::cout << "Iterations: " << result.getIterations() << std::endl;
    std::cout << "Function evaluations: " << result.getFunctionEvaluations() << std::endl;
} else {
    std::cout << "Optimization failed: " << result.getMessage() << std::endl;
}
```

## Integration with RebelCALC

The engineering modules are fully integrated with the core RebelCALC functionality, allowing you to:

- Use engineering functions directly in calculator expressions
- Store engineering objects as variables
- Apply mathematical operations to engineering objects
- Use engineering functions in scripts
- Export engineering calculations to reports

## Integration with RebelSUITE Components

RebelCALC now provides integration with other RebelSUITE components through a unified API. This allows for seamless data exchange and command execution between RebelCALC and other components.

### Integration Architecture

The integration system is built around the following key components:

- **IntegrationInterface**: A base interface that defines the common operations for all component integrations
- **DataExchange**: A class for exchanging data between components using a flexible variant type
- **IntegrationFactory**: A factory for creating integration interfaces for specific components

### Available Integrations

Currently, the following integrations are available:

- **RebelCAD Integration**: Allows RebelCALC to interact with RebelCAD for CAD model operations
  - Import CAD models from RebelCAD
  - Export CAD models to RebelCAD
  - Perform boolean operations on CAD models
  - Measure properties of CAD models (volume, surface area, etc.)

### Using the Integration API

Here's an example of how to use the integration API to interact with RebelCAD:

```cpp
#include "integration/api.h"

using namespace RebelCalc::Integration;

// Get the integration factory
IntegrationFactory& factory = IntegrationFactory::getInstance();

// Create a CAD integration interface
std::shared_ptr<IntegrationInterface> cadInterface = factory.createInterface(Component::REBELCAD);

// Initialize the integration
if (cadInterface->initialize()) {
    // Create a data exchange object
    DataExchange data;
    
    // Import a model from RebelCAD
    data.setData("model_name", std::string("cube"));
    if (cadInterface->executeCommand("import_model", data)) {
        // Get the imported triangles
        std::vector<Engineering::CAD::Triangle> triangles = 
            std::get<std::vector<Engineering::CAD::Triangle>>(data.getData("triangles"));
        
        // Process the triangles...
    }
    
    // Shutdown the integration when done
    cadInterface->shutdown();
}
```

### Lua Script Integration

The integration API is also available in Lua scripts, allowing you to interact with other RebelSUITE components from your scripts:

```lua
-- Get the CAD integration interface
local cad = integration.get_interface("RebelCAD")

-- Initialize the integration
if cad:initialize() then
    -- Import a model from RebelCAD
    local result = cad:execute_command("import_model", {model_name = "cube"})
    
    -- Process the triangles
    local triangles = result.triangles
    for i, triangle in ipairs(triangles) do
        -- Process each triangle...
    end
    
    -- Shutdown the integration when done
    cad:shutdown()
end
```

### Calculator Integration Example

```
// Define components
R = 100 // ohms
C = 1e-6 // farads
L = 0.1 // henries
f = 1000 // Hz

// Calculate impedances
Z_R = impedance_resistor(R)
Z_C = impedance_capacitor(C, f)
Z_L = impedance_inductor(L, f)

// Calculate total impedance (series)
Z_total = Z_R + Z_C + Z_L

// Print results
print("Total impedance: " + Z_total + " ohms")
print("Magnitude: " + abs(Z_total) + " ohms")
print("Phase: " + arg(Z_total) * 180/pi + " degrees")
```

### Lua Script Integration Example

```lua
-- Define a 3D point
function create_point(x, y, z)
    return cad.Point3D.new(x, y, z)
end

-- Calculate distance between two points
function distance(p1, p2)
    return p1:distance(p2)
end

-- Create some points
p1 = create_point(1, 2, 3)
p2 = create_point(4, 5, 6)

-- Calculate and print distance
print("Distance: " .. distance(p1, p2))
```

## Examples

The RebelCALC distribution includes several examples demonstrating the use of the engineering modules:

- **engineering_demo.cpp**: Comprehensive demonstration of all engineering modules
- **integration_demo.cpp**: Demonstration of integration with other RebelSUITE components
- **visualization_demo.cpp**: Demonstration of 2D and 3D plotting capabilities
- **optimization_demo.cpp**: Demonstration of optimization solvers for various problems
- **cad_examples.lua**: Lua script examples for the CAD module
- **physics_simulation.lua**: Particle system simulation using the Physics module
- **circuit_analyzer.lua**: AC circuit analysis using the Electrical module

To run the examples:

```bash
cd build/bin/examples
./engineering_demo
./integration_demo
./visualization_demo
./optimization_demo
```

For more detailed examples and tutorials, see the [RebelCALC User Guide](user_guide.md).
