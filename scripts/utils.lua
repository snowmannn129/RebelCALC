-- RebelCALC Utility Functions
-- This script provides additional utility functions for RebelCALC

-- Print a welcome message
print("Loading RebelCALC utility functions...")

-- Statistical functions

-- Calculate the mean of a list of numbers
function mean(numbers)
    local sum = 0
    for _, value in ipairs(numbers) do
        sum = sum + value
    end
    return sum / #numbers
end

-- Calculate the median of a list of numbers
function median(numbers)
    local sorted = {}
    for i, value in ipairs(numbers) do
        sorted[i] = value
    end
    table.sort(sorted)
    
    local n = #sorted
    if n % 2 == 0 then
        return (sorted[n/2] + sorted[n/2 + 1]) / 2
    else
        return sorted[math.ceil(n/2)]
    end
end

-- Calculate the mode of a list of numbers
function mode(numbers)
    local counts = {}
    for _, value in ipairs(numbers) do
        counts[value] = (counts[value] or 0) + 1
    end
    
    local max_count = 0
    local modes = {}
    for value, count in pairs(counts) do
        if count > max_count then
            max_count = count
            modes = {value}
        elseif count == max_count then
            table.insert(modes, value)
        end
    end
    
    return modes
end

-- Calculate the standard deviation of a list of numbers
function std_dev(numbers)
    local m = mean(numbers)
    local sum_squared_diff = 0
    for _, value in ipairs(numbers) do
        sum_squared_diff = sum_squared_diff + (value - m)^2
    end
    return math.sqrt(sum_squared_diff / #numbers)
end

-- Calculate the variance of a list of numbers
function variance(numbers)
    local sd = std_dev(numbers)
    return sd^2
end

-- Calculate the correlation coefficient between two lists of numbers
function correlation(x, y)
    if #x ~= #y then
        error("Lists must have the same length")
    end
    
    local n = #x
    local sum_x = 0
    local sum_y = 0
    local sum_xy = 0
    local sum_x2 = 0
    local sum_y2 = 0
    
    for i = 1, n do
        sum_x = sum_x + x[i]
        sum_y = sum_y + y[i]
        sum_xy = sum_xy + x[i] * y[i]
        sum_x2 = sum_x2 + x[i]^2
        sum_y2 = sum_y2 + y[i]^2
    end
    
    local numerator = n * sum_xy - sum_x * sum_y
    local denominator = math.sqrt((n * sum_x2 - sum_x^2) * (n * sum_y2 - sum_y^2))
    
    return numerator / denominator
end

-- Linear regression
function linear_regression(x, y)
    if #x ~= #y then
        error("Lists must have the same length")
    end
    
    local n = #x
    local sum_x = 0
    local sum_y = 0
    local sum_xy = 0
    local sum_x2 = 0
    
    for i = 1, n do
        sum_x = sum_x + x[i]
        sum_y = sum_y + y[i]
        sum_xy = sum_xy + x[i] * y[i]
        sum_x2 = sum_x2 + x[i]^2
    end
    
    local slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x^2)
    local intercept = (sum_y - slope * sum_x) / n
    
    return {slope = slope, intercept = intercept}
end

-- Numerical methods

-- Newton-Raphson method for finding roots
function newton_raphson(f, df, x0, tolerance, max_iterations)
    tolerance = tolerance or 1e-10
    max_iterations = max_iterations or 100
    
    local x = x0
    for i = 1, max_iterations do
        local fx = f(x)
        if math.abs(fx) < tolerance then
            return x
        end
        
        local dfx = df(x)
        if dfx == 0 then
            error("Derivative is zero, cannot continue")
        end
        
        x = x - fx / dfx
    end
    
    error("Failed to converge after " .. max_iterations .. " iterations")
end

-- Bisection method for finding roots
function bisection(f, a, b, tolerance, max_iterations)
    tolerance = tolerance or 1e-10
    max_iterations = max_iterations or 100
    
    local fa = f(a)
    local fb = f(b)
    
    if fa * fb > 0 then
        error("Function must have opposite signs at the interval endpoints")
    end
    
    for i = 1, max_iterations do
        local c = (a + b) / 2
        local fc = f(c)
        
        if math.abs(fc) < tolerance then
            return c
        end
        
        if fa * fc < 0 then
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
    end
    
    error("Failed to converge after " .. max_iterations .. " iterations")
end

-- Simpson's rule for numerical integration
function simpson_integrate(f, a, b, n)
    n = n or 1000
    if n % 2 ~= 0 then
        n = n + 1  -- Ensure n is even
    end
    
    local h = (b - a) / n
    local sum = f(a) + f(b)
    
    for i = 1, n-1 do
        local x = a + i * h
        if i % 2 == 0 then
            sum = sum + 2 * f(x)
        else
            sum = sum + 4 * f(x)
        end
    end
    
    return (h / 3) * sum
end

-- Print available functions
print("Available utility functions:")
print("  Statistical functions:")
print("    mean(numbers) - Calculate the mean of a list of numbers")
print("    median(numbers) - Calculate the median of a list of numbers")
print("    mode(numbers) - Calculate the mode of a list of numbers")
print("    std_dev(numbers) - Calculate the standard deviation of a list of numbers")
print("    variance(numbers) - Calculate the variance of a list of numbers")
print("    correlation(x, y) - Calculate the correlation coefficient between two lists of numbers")
print("    linear_regression(x, y) - Perform linear regression on two lists of numbers")
print("  Numerical methods:")
print("    newton_raphson(f, df, x0, tolerance, max_iterations) - Find a root using the Newton-Raphson method")
print("    bisection(f, a, b, tolerance, max_iterations) - Find a root using the bisection method")
print("    simpson_integrate(f, a, b, n) - Integrate a function using Simpson's rule")

print("RebelCALC utility functions loaded successfully!")
