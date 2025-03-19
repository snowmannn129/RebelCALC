# RebelCALC User Guide

## Introduction

RebelCALC is an advanced computational engine and calculator, developed as a standalone yet integrable component of RebelSUITE. It is designed to handle complex engineering, CAD, and simulation-based calculations, providing powerful mathematical capabilities for a wide range of applications.

This user guide will help you get started with RebelCALC and explore its features.

## Installation

### Prerequisites

- CMake 3.15+
- C++17 compatible compiler
- Lua 5.4+

### Windows

1. Clone or download the RebelCALC repository
2. Run the `download_lua.ps1` script to download and build Lua:
   ```powershell
   .\download_lua.ps1
   ```
3. Run the `download_gtest.ps1` script to download and build GoogleTest (optional, for running tests):
   ```powershell
   .\download_gtest.ps1
   ```
4. Build RebelCALC using the provided build script:
   ```powershell
   .\build.bat
   ```
5. The executable will be located at `build\bin\Release\RebelCALC.exe`

### Unix-like Systems (Linux, macOS)

1. Clone or download the RebelCALC repository
2. Install Lua 5.4+ using your package manager:
   ```bash
   # Ubuntu/Debian
   sudo apt-get install liblua5.4-dev
   
   # macOS with Homebrew
   brew install lua
   ```
3. Install GoogleTest (optional, for running tests):
   ```bash
   # Ubuntu/Debian
   sudo apt-get install libgtest-dev
   
   # macOS with Homebrew
   brew install googletest
   ```
4. Build RebelCALC using the provided build script:
   ```bash
   ./build.sh
   ```
5. The executable will be located at `build/bin/RebelCALC`

## Getting Started

### Running RebelCALC

To run RebelCALC, simply execute the built binary:

```bash
# Windows
build\bin\Release\RebelCALC.exe

# Unix-like systems
./build/bin/RebelCALC
```

### Basic Usage

RebelCALC provides a command-line interface for performing calculations. You can enter mathematical expressions directly, or use commands to access more advanced features.

#### Evaluating Expressions

To evaluate a mathematical expression, simply type it at the prompt:

```
> 2 + 2
Result: 4

> sin(pi/2)
Result: 1

> sqrt(16)
Result: 4
```

#### Using Variables

You can define and use variables in your calculations:

```
> x = 5
Variable 'x' set to 5

> y = 10
Variable 'y' set to 10

> x + y
Result: 15

> x * y
Result: 50
```

#### Using Commands

RebelCALC provides a set of commands for accessing more advanced features. Commands are prefixed with a slash (`/`):

```
> /help
Available commands:
  /help - Display help information
  /exit - Exit the application
  /clear - Clear the screen
  /theme - Set or display the current theme
  /vars - Display all defined variables
  /clear_vars - Clear all user-defined variables
  /solve - Solve an equation for a variable
  /diff - Differentiate an expression with respect to a variable
  /int - Integrate an expression with respect to a variable
  /lua - Execute a Lua script
  /load - Load a Lua script from a file
```

## Advanced Features

### Symbolic Mathematics

RebelCALC includes a symbolic mathematics engine that can perform symbolic operations such as simplification, expansion, factorization, differentiation, and integration.

#### Solving Equations

You can solve equations for a specific variable using the `/solve` command:

```
> /solve x + 5 = 10 x
Solution: 5

> /solve 2*x = 10 x
Solution: 5

> /solve x^2 = 4 x
Solutions: -2, 2
```

#### Differentiation

You can differentiate expressions with respect to a variable using the `/diff` command:

```
> /diff x^2 x
Derivative: 2*x

> /diff sin(x) x
Derivative: cos(x)
```

#### Integration

You can integrate expressions with respect to a variable using the `/int` command:

```
> /int x x
Integral: x^2/2

> /int x^2 x
Integral: x^3/3
```

### Numerical Computation

RebelCALC includes a numerical solver that can perform numerical operations such as solving linear systems, finding roots of polynomials, and numerical integration.

#### Solving Linear Systems

You can solve systems of linear equations using the Lua scripting interface:

```
> /lua return solve_system({{2, 3}, {4, 9}}, {8, 22})
Result: {1, 2}
```

#### Finding Roots of Polynomials

You can find the roots of polynomials using the Lua scripting interface:

```
> /lua return find_roots({1, 0, -4})
Result: {-2, 2}
```

### Lua Scripting

RebelCALC includes a Lua scripting engine that allows you to write and execute Lua scripts. This provides a powerful way to extend the functionality of RebelCALC and automate complex calculations.

#### Executing Scripts

You can execute Lua scripts directly using the `/lua` command:

```
> /lua return 2 + 2
Result: 4

> /lua
local x = 10
local y = 20
return x + y
Result: 30
```

#### Loading Scripts from Files

You can load and execute Lua scripts from files using the `/load` command:

```
> /load scripts/examples/basic_math.lua
Script loaded successfully
```

#### Built-in Functions

RebelCALC provides a set of built-in functions that can be called from Lua scripts:

- `evaluate(expression)` - Evaluate a mathematical expression
- `solve(equation, variable)` - Solve an equation for a variable
- `differentiate(expression, variable)` - Differentiate an expression with respect to a variable
- `integrate(expression, variable)` - Integrate an expression with respect to a variable

#### Utility Functions

RebelCALC includes a set of utility functions that provide additional mathematical capabilities:

- Statistical functions: `mean`, `median`, `mode`, `std_dev`, `variance`, `correlation`, `linear_regression`
- Numerical methods: `newton_raphson`, `bisection`, `simpson_integrate`

## Customization

### Themes

RebelCALC supports multiple themes for the user interface. You can set the theme using the `/theme` command:

```
> /theme
Current theme: default

> /theme dark
Theme set to dark

> /theme light
Theme set to light
```

## Integration with RebelSUITE

RebelCALC is designed to integrate seamlessly with other components of RebelSUITE, such as RebelCAD, RebelENGINE, and RebelSIM. This integration allows you to use RebelCALC's powerful mathematical capabilities within these applications.

## Troubleshooting

### Common Issues

- **Error: "Failed to initialize Lua engine"**: Make sure Lua is installed correctly and the Lua libraries are in the correct location.
- **Error: "Failed to load script from file"**: Make sure the script file exists and is readable.
- **Error: "Failed to execute Lua script"**: Check the script for syntax errors or other issues.

### Getting Help

If you encounter any issues or have questions about RebelCALC, please refer to the documentation or contact support.

## Conclusion

RebelCALC is a powerful computational engine that provides a wide range of mathematical capabilities. Whether you're performing simple calculations or complex mathematical operations, RebelCALC has the tools you need to get the job done.

Happy calculating!
