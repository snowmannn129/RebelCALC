#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <functional>
#include <variant>
#include <optional>
#include "matrix.h"
#include "complex.h"
#include "statistics.h"
#include "units.h"

namespace rebelcalc {

// Forward declarations
class SymbolicEngine;
class NumericSolver;
class Statistics;

/**
 * Result type for calculator operations
 */
using CalculatorValue = std::variant<double, std::string, std::vector<double>, Matrix, Complex>;

/**
 * Main calculator class that handles all computation operations
 */
class Calculator {
public:
    /**
     * Constructor
     */
    Calculator();
    
    /**
     * Destructor
     */
    ~Calculator();
    
    /**
     * Initialize the calculator
     * @return true if initialization was successful, false otherwise
     */
    bool initialize();
    
    /**
     * Shutdown the calculator
     */
    void shutdown();
    
    /**
     * Evaluate a mathematical expression
     * @param expression The expression to evaluate
     * @return The result of the evaluation, or nullopt if evaluation failed
     */
    std::optional<CalculatorValue> evaluate(const std::string& expression);
    
    /**
     * Solve an equation for a variable
     * @param equation The equation to solve
     * @param variable The variable to solve for
     * @return The solution, or nullopt if solving failed
     */
    std::optional<CalculatorValue> solve(const std::string& equation, const std::string& variable);
    
    /**
     * Differentiate an expression with respect to a variable
     * @param expression The expression to differentiate
     * @param variable The variable to differentiate with respect to
     * @return The derivative, or nullopt if differentiation failed
     */
    std::optional<std::string> differentiate(const std::string& expression, const std::string& variable);
    
    /**
     * Integrate an expression with respect to a variable
     * @param expression The expression to integrate
     * @param variable The variable to integrate with respect to
     * @return The integral, or nullopt if integration failed
     */
    std::optional<std::string> integrate(const std::string& expression, const std::string& variable);
    
    /**
     * Set a variable value
     * @param name The variable name
     * @param value The variable value
     * @return true if the variable was set successfully, false otherwise
     */
    bool setVariable(const std::string& name, double value);
    
    /**
     * Set a complex variable value
     * @param name The variable name
     * @param value The complex variable value
     * @return true if the variable was set successfully, false otherwise
     */
    bool setComplexVariable(const std::string& name, const Complex& value);
    
    /**
     * Get a variable value
     * @param name The variable name
     * @return The variable value, or nullopt if the variable doesn't exist
     */
    std::optional<double> getVariable(const std::string& name) const;
    
    /**
     * Get a complex variable value
     * @param name The variable name
     * @return The complex variable value, or nullopt if the variable doesn't exist
     */
    std::optional<Complex> getComplexVariable(const std::string& name) const;
    
    /**
     * Check if a variable exists
     * @param name The variable name
     * @return true if the variable exists, false otherwise
     */
    bool hasVariable(const std::string& name) const;
    
    /**
     * Get all variables
     * @return A map of variable names to values
     */
    std::unordered_map<std::string, double> getAllVariables() const;
    
    /**
     * Clear all variables
     */
    void clearVariables();
    
    /**
     * Process a variable assignment expression
     * @param expression The assignment expression (e.g., "x = 5", "y += 3")
     * @return The result of the assignment, or nullopt if the assignment failed
     */
    std::optional<double> processVariableAssignment(const std::string& expression);
    
    /**
     * Register a custom function
     * @param name The function name
     * @param function The function implementation
     * @param argCount The number of arguments the function takes
     */
    void registerFunction(const std::string& name, 
                         std::function<double(const std::vector<double>&)> function,
                         size_t argCount);
    
    /**
     * Get the symbolic engine
     * @return The symbolic engine
     */
    std::shared_ptr<SymbolicEngine> getSymbolicEngine() const;
    
    /**
     * Get the numeric solver
     * @return The numeric solver
     */
    std::shared_ptr<NumericSolver> getNumericSolver() const;
    
    /**
     * Get the statistics engine
     * @return The statistics engine
     */
    std::shared_ptr<Statistics> getStatistics() const;
    
    /**
     * Get the unit converter
     * @return The unit converter
     */
    std::shared_ptr<UnitConverter> getUnitConverter() const;
    
    /**
     * Convert a value from one unit to another
     * @param value The value to convert
     * @param fromUnit The source unit
     * @param toUnit The target unit
     * @return The converted value, or nullopt if conversion failed
     */
    std::optional<double> convertUnit(double value, const std::string& fromUnit, const std::string& toUnit);
    
    /**
     * Get all available units for a specific quantity
     * @param quantity The physical quantity (e.g., "length", "mass", "time")
     * @return A vector of unit names, or empty vector if quantity is not supported
     */
    std::vector<std::string> getUnitsForQuantity(const std::string& quantity) const;
    
    /**
     * Get all supported physical quantities
     * @return A vector of physical quantity names
     */
    std::vector<std::string> getAllQuantities() const;
    
    /**
     * Information about a registered function
     */
    struct FunctionInfo {
        std::function<double(const std::vector<double>&)> function;
        size_t argCount;
    };

private:
    // Variables
    std::unordered_map<std::string, double> m_variables;
    std::unordered_map<std::string, Complex> m_complexVariables;
    
    // Custom functions
    std::unordered_map<std::string, FunctionInfo> m_functions;
    
    // Components
    std::shared_ptr<SymbolicEngine> m_symbolicEngine;
    std::shared_ptr<NumericSolver> m_numericSolver;
    std::shared_ptr<Statistics> m_statistics;
    std::shared_ptr<UnitConverter> m_unitConverter;
    
    // Helper methods
    bool isNumeric(const std::string& expression) const;
    double evaluateNumeric(const std::string& expression) const;
    std::string evaluateSymbolic(const std::string& expression) const;
    
    // Complex number helper methods
    bool isComplex(const std::string& expression) const;
    Complex evaluateComplex(const std::string& expression) const;
    std::optional<Complex> parseComplex(const std::string& expression) const;
};

} // namespace rebelcalc
