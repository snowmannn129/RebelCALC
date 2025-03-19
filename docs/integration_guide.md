# RebelCALC Integration Guide

This document outlines how to integrate RebelCALC with other components of RebelSUITE.

## Overview

RebelCALC is designed to be both a standalone application and a component that can be integrated with other parts of RebelSUITE. This guide explains how to integrate RebelCALC with RebelCAD, RebelENGINE, RebelSIM, and other components.

## Integration Methods

There are several ways to integrate RebelCALC with other components:

1. **Library Integration**: Link RebelCALC as a library in your application
2. **IPC Communication**: Use inter-process communication to interact with a running RebelCALC instance
3. **Script Integration**: Use RebelCALC's scripting capabilities to extend your application

## Library Integration

### Prerequisites

- CMake 3.15+
- C++17 compatible compiler
- Lua 5.4+

### Including RebelCALC in Your CMake Project

To include RebelCALC as a library in your CMake project, add the following to your `CMakeLists.txt`:

```cmake
# Add RebelCALC as a subdirectory
add_subdirectory(path/to/RebelCALC)

# Link against RebelCALC library
target_link_libraries(YourTarget PRIVATE RebelCALC)
```

### Using RebelCALC in Your Code

Once you've linked against RebelCALC, you can use it in your code:

```cpp
#include "rebelcalc/calculator.h"

// Create a calculator instance
auto calculator = std::make_shared<rebelcalc::Calculator>();
calculator->initialize();

// Evaluate an expression
auto result = calculator->evaluate("2 + 2");
if (result) {
    // Handle the result
    std::visit([](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, double>) {
            std::cout << "Result: " << arg << std::endl;
        } else if constexpr (std::is_same_v<T, std::string>) {
            std::cout << "Result: " << arg << std::endl;
        } else if constexpr (std::is_same_v<T, std::vector<double>>) {
            std::cout << "Result: ";
            for (size_t i = 0; i < arg.size(); ++i) {
                std::cout << arg[i];
                if (i < arg.size() - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << std::endl;
        }
    }, *result);
}

// Clean up
calculator->shutdown();
```

## IPC Communication

RebelCALC can be run as a separate process and communicated with using inter-process communication (IPC). This approach is useful when you want to keep RebelCALC separate from your main application.

### Starting RebelCALC in IPC Mode

To start RebelCALC in IPC mode, run it with the `--ipc` flag:

```bash
RebelCALC --ipc
```

This will start RebelCALC as a server that listens for IPC connections.

### Connecting to RebelCALC from Your Application

You can connect to RebelCALC from your application using the RebelCALC IPC client library:

```cpp
#include "rebelcalc/ipc_client.h"

// Create an IPC client
auto client = std::make_shared<rebelcalc::IPCClient>();
client->connect("localhost", 8080);

// Evaluate an expression
auto result = client->evaluate("2 + 2");
if (result) {
    // Handle the result
    std::cout << "Result: " << *result << std::endl;
}

// Clean up
client->disconnect();
```

## Script Integration

RebelCALC includes a Lua scripting engine that can be used to extend your application. This approach is useful when you want to add scripting capabilities to your application.

### Embedding the Lua Engine

You can embed the RebelCALC Lua engine in your application:

```cpp
#include "rebelcalc/calculator.h"
#include "rebelcalc/lua_engine.h"

// Create a calculator instance
auto calculator = std::make_shared<rebelcalc::Calculator>();
calculator->initialize();

// Create a Lua engine instance
auto luaEngine = std::make_shared<rebelcalc::LuaEngine>(calculator);
luaEngine->initialize();

// Execute a Lua script
auto result = luaEngine->executeScript("return 2 + 2");
if (result) {
    std::cout << "Result: " << *result << std::endl;
}

// Clean up
luaEngine->shutdown();
calculator->shutdown();
```

### Using RebelCALC Scripts in Your Application

You can also load and execute RebelCALC scripts in your application:

```cpp
#include "rebelcalc/calculator.h"
#include "rebelcalc/lua_engine.h"

// Create a calculator instance
auto calculator = std::make_shared<rebelcalc::Calculator>();
calculator->initialize();

// Create a Lua engine instance
auto luaEngine = std::make_shared<rebelcalc::LuaEngine>(calculator);
luaEngine->initialize();

// Load and execute a script
if (luaEngine->loadScript("path/to/script.lua")) {
    std::cout << "Script loaded successfully" << std::endl;
}

// Clean up
luaEngine->shutdown();
calculator->shutdown();
```

## Integration with RebelCAD

RebelCALC can be integrated with RebelCAD to provide advanced mathematical capabilities for CAD operations.

### Example: Calculating Properties of a CAD Model

```cpp
#include "rebelcad/model.h"
#include "rebelcalc/calculator.h"

// Load a CAD model
auto model = rebelcad::Model::load("path/to/model.rcad");

// Create a calculator instance
auto calculator = std::make_shared<rebelcalc::Calculator>();
calculator->initialize();

// Calculate the volume of the model
double volume = model->getVolume();
calculator->setVariable("volume", volume);

// Calculate the surface area of the model
double surfaceArea = model->getSurfaceArea();
calculator->setVariable("surface_area", surfaceArea);

// Calculate the surface-to-volume ratio
auto result = calculator->evaluate("surface_area / volume");
if (result) {
    std::visit([](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, double>) {
            std::cout << "Surface-to-volume ratio: " << arg << std::endl;
        }
    }, *result);
}

// Clean up
calculator->shutdown();
```

## Integration with RebelENGINE

RebelCALC can be integrated with RebelENGINE to provide advanced mathematical capabilities for game and simulation operations.

### Example: Calculating Physics Properties

```cpp
#include "rebelengine/physics.h"
#include "rebelcalc/calculator.h"

// Create a physics object
auto physics = std::make_shared<rebelengine::Physics>();

// Create a calculator instance
auto calculator = std::make_shared<rebelcalc::Calculator>();
calculator->initialize();

// Calculate the trajectory of a projectile
double initialVelocity = 10.0;
double angle = 45.0;
double gravity = 9.81;

calculator->setVariable("v0", initialVelocity);
calculator->setVariable("angle", angle);
calculator->setVariable("g", gravity);

// Calculate the maximum height
auto maxHeight = calculator->evaluate("(v0^2 * sin(angle * pi/180)^2) / (2 * g)");
if (maxHeight) {
    std::visit([](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, double>) {
            std::cout << "Maximum height: " << arg << " m" << std::endl;
        }
    }, *maxHeight);
}

// Calculate the range
auto range = calculator->evaluate("(v0^2 * sin(2 * angle * pi/180)) / g");
if (range) {
    std::visit([](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, double>) {
            std::cout << "Range: " << arg << " m" << std::endl;
        }
    }, *range);
}

// Clean up
calculator->shutdown();
```

## Integration with RebelSIM

RebelCALC can be integrated with RebelSIM to provide advanced mathematical capabilities for simulation operations.

### Example: Analyzing Simulation Results

```cpp
#include "rebelsim/simulation.h"
#include "rebelcalc/calculator.h"
#include "rebelcalc/lua_engine.h"

// Run a simulation
auto simulation = std::make_shared<rebelsim::Simulation>();
simulation->run();

// Get the simulation results
auto results = simulation->getResults();

// Create a calculator instance
auto calculator = std::make_shared<rebelcalc::Calculator>();
calculator->initialize();

// Create a Lua engine instance
auto luaEngine = std::make_shared<rebelcalc::LuaEngine>(calculator);
luaEngine->initialize();

// Analyze the results using Lua
std::string script = R"(
    -- Get the simulation results
    local results = ...
    
    -- Calculate the mean
    local sum = 0
    for i, value in ipairs(results) do
        sum = sum + value
    end
    local mean = sum / #results
    
    -- Calculate the standard deviation
    local sum_squared_diff = 0
    for i, value in ipairs(results) do
        sum_squared_diff = sum_squared_diff + (value - mean)^2
    end
    local std_dev = math.sqrt(sum_squared_diff / #results)
    
    -- Return the results
    return {mean = mean, std_dev = std_dev}
)";

// Execute the script with the simulation results
auto result = luaEngine->executeScript(script, results);
if (result) {
    std::cout << "Analysis results: " << *result << std::endl;
}

// Clean up
luaEngine->shutdown();
calculator->shutdown();
```

## Conclusion

RebelCALC is designed to be easily integrated with other components of RebelSUITE. Whether you're using it as a library, communicating with it via IPC, or using its scripting capabilities, RebelCALC provides powerful mathematical capabilities for your applications.

For more information, please refer to the RebelCALC API documentation.
