-- RebelCALC Initialization Script
-- This script is loaded automatically when RebelCALC starts

-- Print a welcome message
print("Loading RebelCALC initialization script...")

-- Define some useful constants
PI = evaluate("pi")
E = evaluate("e")
PHI = (1 + math.sqrt(5)) / 2
GAMMA = 0.57721566490153286061

-- Define some useful functions

-- Convert degrees to radians
function deg_to_rad(degrees)
    return degrees * PI / 180
end

-- Convert radians to degrees
function rad_to_deg(radians)
    return radians * 180 / PI
end

-- Factorial function
function factorial(n)
    if n == 0 then
        return 1
    else
        return n * factorial(n - 1)
    end
end

-- Combination function (n choose k)
function combination(n, k)
    return factorial(n) / (factorial(k) * factorial(n - k))
end

-- Permutation function
function permutation(n, k)
    return factorial(n) / factorial(n - k)
end

-- Fibonacci function
function fibonacci(n)
    if n <= 0 then
        return 0
    elseif n == 1 then
        return 1
    else
        return fibonacci(n - 1) + fibonacci(n - 2)
    end
end

-- GCD function (greatest common divisor)
function gcd(a, b)
    if b == 0 then
        return a
    else
        return gcd(b, a % b)
    end
end

-- LCM function (least common multiple)
function lcm(a, b)
    return (a * b) / gcd(a, b)
end

-- Function to solve a system of linear equations
-- Example: solve_system({{2, 3}, {4, 9}}, {8, 22}) returns {1, 2}
function solve_system(coefficients, constants)
    local result = evaluate("solve_linear_system(" .. serialize(coefficients) .. ", " .. serialize(constants) .. ")")
    return result
end

-- Function to find roots of a polynomial
-- Example: find_roots({1, 0, -4}) returns {-2, 2}
function find_roots(coefficients)
    local result = evaluate("find_roots(" .. serialize(coefficients) .. ")")
    return result
end

-- Helper function to serialize a table to a string
function serialize(obj)
    local lua_type = type(obj)
    if lua_type == "number" or lua_type == "boolean" then
        return tostring(obj)
    elseif lua_type == "string" then
        return string.format("%q", obj)
    elseif lua_type == "table" then
        local result = "{"
        for k, v in pairs(obj) do
            if type(k) ~= "number" then
                result = result .. "[" .. serialize(k) .. "]="
            end
            result = result .. serialize(v) .. ","
        end
        return result .. "}"
    else
        return "nil"
    end
end

-- Print available functions
print("Available functions:")
print("  deg_to_rad(degrees) - Convert degrees to radians")
print("  rad_to_deg(radians) - Convert radians to degrees")
print("  factorial(n) - Calculate factorial of n")
print("  combination(n, k) - Calculate n choose k")
print("  permutation(n, k) - Calculate permutation of n things taken k at a time")
print("  fibonacci(n) - Calculate the nth Fibonacci number")
print("  gcd(a, b) - Calculate greatest common divisor")
print("  lcm(a, b) - Calculate least common multiple")
print("  solve_system(coefficients, constants) - Solve a system of linear equations")
print("  find_roots(coefficients) - Find roots of a polynomial")

print("RebelCALC initialization completed successfully!")
