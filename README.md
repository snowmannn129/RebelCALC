# RebelCALC

## Advanced Computational Engine for RebelSUITE

RebelCALC is a next-generation computational engine and calculator, developed as a standalone yet integrable component of RebelSUITE. Designed to handle complex engineering, CAD, and simulation-based calculations, RebelCALC powers advanced workflows in RebelCAD, RebelENGINE, and other components.

## Core Features

- **Advanced Computation** - Symbolic algebra, numerical solvers, and real-time graphing
- **Engineering & Scientific Functions** - CAD-specific calculations, physics engine support, statistical analysis
- **Graphing & Visualization** - 2D/3D function plotting, CAD & simulation visualization
- **Lua Scripting** - Automate calculations, define reusable formulas, and process large datasets
- **Custom ANSI-Inspired UI** - Fully themeable interface with keyboard-first navigation
- **Modular Integration** - Seamlessly integrates with other RebelSUITE components

## Current Capabilities

RebelCALC currently supports:

- **Basic Arithmetic Operations** - Addition, subtraction, multiplication, division, exponentiation
- **Mathematical Functions** - Trigonometric, logarithmic, exponential, and more
- **Variable Management** - Define, use, and manipulate variables in expressions
- **Expression Parsing** - Evaluate complex mathematical expressions with proper precedence
- **Command-Line Interface** - User-friendly terminal interface with comprehensive help
- **Matrix Operations** - Comprehensive matrix algebra with various decompositions and operations
- **Complex Numbers** - Full support for complex arithmetic and functions
- **Statistical Functions** - Descriptive statistics, regression analysis, probability distributions, and statistical tests
- **Unit Conversion** - Comprehensive unit conversion system with support for SI, imperial, and other common unit systems

## Usage Examples

```
> 2 + 3 * 4             # Basic arithmetic with operator precedence
14

> sin(pi/2)             # Trigonometric functions with constants
1

> x = 5                 # Variable assignment
x = 5

> y = 2*x + 3           # Assignment with expressions
y = 13

> sqrt(x^2 + y^2)       # Pythagorean theorem
13.9284

> /vars                 # Display all variables
Variables:
  Built-in constants:
    pi = 3.14159
    e = 2.71828
  User-defined variables:
    x = 5
    y = 13

> /solve x^2 - 4 = 0 x  # Solve equation for x
Solution: x = 2 or x = -2

> /diff x^2 x           # Differentiate with respect to x
Derivative: 2*x

> A = [1, 2; 3, 4]      # Matrix operations
A = [1, 2; 3, 4]

> det(A)                # Matrix determinant
-2

> 3+4i                  # Complex number literal
3 + 4i

> (2+3i) * (1-2i)       # Complex number arithmetic
8 + -1i

> data = [1, 2, 3, 4, 5] # Statistical functions
> mean(data)            # Calculate mean
3
> stddev(data)          # Calculate standard deviation
1.5811
> corr([1,2,3], [2,4,6]) # Calculate correlation
1.0

> convert(5, "km", "mi") # Unit conversion
3.10686
> convert(98.6, "F", "C") # Temperature conversion
37
> convert(1, "L", "gal") # Volume conversion
0.26417
```

## Development Status

RebelCALC is currently in Phase 1 of development. See the [progress.md](progress.md) file for detailed information on the current status and planned features.

## Building and Running

### Prerequisites

- CMake 3.15+
- C++17 compatible compiler
- Lua 5.4+ (for scripting features)

### Build Instructions

#### Windows
```powershell
# Download dependencies
.\download_lua.ps1
.\download_yaml_cpp.ps1
.\download_gtest.ps1

# Build the project
.\build.bat
```

#### Unix-like Systems
```bash
# Download dependencies (if not already installed)
./download_lua.sh
./download_yaml_cpp.sh
./download_gtest.sh

# Build the project
./build.sh
```

### Running

```bash
./build/bin/RebelCALC
```

## Testing

```bash
./run_tests.bat  # Windows
./run_tests.sh   # Unix-like systems
```

## License

Proprietary - All rights reserved.
