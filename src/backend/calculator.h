#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <functional>
#include <variant>
#include <optional>

namespace rebelcalc {

// Forward declarations
class SymbolicEngine;
class NumericSolver;

/**
 * Result type for calculator operations
 */
using CalculatorValue = std::variant<double, std::string, std::vector<double>>;

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
     */
    void setVariable(const std::string& name, double value);
    
    /**
     * Get a variable value
     * @param name The variable name
     * @return The variable value, or nullopt if the variable doesn't exist
     */
    std::optional<double> getVariable(const std::string& name) const;
    
    /**
     * Clear all variables
     */
    void clearVariables();
    
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

private:
    // Variables
    std::unordered_map<std::string, double> m_variables;
    
    // Custom functions
    struct FunctionInfo {
        std::function<double(const std::vector<double>&)> function;
        size_t argCount;
    };
    std::unordered_map<std::string, FunctionInfo> m_functions;
    
    // Components
    std::shared_ptr<SymbolicEngine> m_symbolicEngine;
    std::shared_ptr<NumericSolver> m_numericSolver;
    
    // Helper methods
    bool isNumeric(const std::string& expression) const;
    double evaluateNumeric(const std::string& expression) const;
    std::string evaluateSymbolic(const std::string& expression) const;
};

} // namespace rebelcalc
