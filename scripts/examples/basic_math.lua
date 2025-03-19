-- Basic math operations example for RebelCALC
-- This script demonstrates how to use the Lua scripting capabilities of RebelCALC

-- Print a welcome message
print("RebelCALC Lua Script Example: Basic Math Operations")
print("===================================================")

-- Define some variables
local x = 10
local y = 5

-- Perform basic arithmetic operations
print("Basic Arithmetic:")
print("x = " .. x)
print("y = " .. y)
print("x + y = " .. (x + y))
print("x - y = " .. (x - y))
print("x * y = " .. (x * y))
print("x / y = " .. (x / y))
print("x % y = " .. (x % y))
print("x ^ y = " .. (x ^ y))
print()

-- Use RebelCALC's evaluate function
print("Using RebelCALC's evaluate function:")
print("evaluate('1+1') = " .. evaluate("1+1"))
print("evaluate('sin(pi/2)') = " .. evaluate("sin(pi/2)"))
print("evaluate('sqrt(16)') = " .. evaluate("sqrt(16)"))
print()

-- Solve an equation
print("Solving equations:")
print("solve('x + 5 = 10', 'x') = " .. solve("x + 5 = 10", "x"))
print("solve('2*x = 10', 'x') = " .. solve("2*x = 10", "x"))
print()

-- Differentiate an expression
print("Differentiation:")
print("differentiate('x^2', 'x') = " .. differentiate("x^2", "x"))
print("differentiate('sin(x)', 'x') = " .. differentiate("sin(x)", "x"))
print()

-- Integrate an expression
print("Integration:")
print("integrate('x', 'x') = " .. integrate("x", "x"))
print("integrate('x^2', 'x') = " .. integrate("x^2", "x"))
print()

-- Define a function to calculate the area of a circle
function area_of_circle(radius)
    return evaluate("pi * " .. radius .. "^2")
end

-- Use the function
print("Custom functions:")
print("Area of circle with radius 5 = " .. area_of_circle(5))
print()

-- Define a function to calculate the volume of a sphere
function volume_of_sphere(radius)
    return evaluate("(4/3) * pi * " .. radius .. "^3")
end

-- Use the function
print("Volume of sphere with radius 5 = " .. volume_of_sphere(5))
print()

-- Define a function to solve a quadratic equation
function solve_quadratic(a, b, c)
    local equation = a .. "*x^2 + " .. b .. "*x + " .. c .. " = 0"
    return solve(equation, "x")
end

-- Use the function
print("Solving quadratic equation x^2 + 2x + 1 = 0:")
print(solve_quadratic(1, 2, 1))
print()

print("Script execution completed successfully!")
