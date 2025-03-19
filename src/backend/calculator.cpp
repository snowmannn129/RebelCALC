#include "calculator.h"
#include "symbolic_engine.h"
#include "numeric_solver.h"

#include <iostream>
#include <cmath>
#include <regex>
#include <stdexcept>

// Define M_PI and M_E if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

namespace rebelcalc {

Calculator::Calculator() 
    : m_symbolicEngine(nullptr), m_numericSolver(nullptr) {
    // Initialize built-in variables
    m_variables["pi"] = M_PI;
    m_variables["e"] = M_E;
    
    // Initialize built-in functions
    registerFunction("sin", [](const std::vector<double>& args) { return std::sin(args[0]); }, 1);
    registerFunction("cos", [](const std::vector<double>& args) { return std::cos(args[0]); }, 1);
    registerFunction("tan", [](const std::vector<double>& args) { return std::tan(args[0]); }, 1);
    registerFunction("asin", [](const std::vector<double>& args) { return std::asin(args[0]); }, 1);
    registerFunction("acos", [](const std::vector<double>& args) { return std::acos(args[0]); }, 1);
    registerFunction("atan", [](const std::vector<double>& args) { return std::atan(args[0]); }, 1);
    registerFunction("sqrt", [](const std::vector<double>& args) { return std::sqrt(args[0]); }, 1);
    registerFunction("log", [](const std::vector<double>& args) { return std::log(args[0]); }, 1);
    registerFunction("log10", [](const std::vector<double>& args) { return std::log10(args[0]); }, 1);
    registerFunction("exp", [](const std::vector<double>& args) { return std::exp(args[0]); }, 1);
    registerFunction("abs", [](const std::vector<double>& args) { return std::abs(args[0]); }, 1);
    registerFunction("floor", [](const std::vector<double>& args) { return std::floor(args[0]); }, 1);
    registerFunction("ceil", [](const std::vector<double>& args) { return std::ceil(args[0]); }, 1);
    registerFunction("round", [](const std::vector<double>& args) { return std::round(args[0]); }, 1);
    registerFunction("pow", [](const std::vector<double>& args) { return std::pow(args[0], args[1]); }, 2);
    registerFunction("max", [](const std::vector<double>& args) { return std::max(args[0], args[1]); }, 2);
    registerFunction("min", [](const std::vector<double>& args) { return std::min(args[0], args[1]); }, 2);
}

Calculator::~Calculator() {
    shutdown();
}

bool Calculator::initialize() {
    try {
        // Create components
        m_symbolicEngine = std::make_shared<SymbolicEngine>();
        m_numericSolver = std::make_shared<NumericSolver>();
        
        // Initialize components
        if (!m_symbolicEngine->initialize()) {
            std::cerr << "Failed to initialize symbolic engine" << std::endl;
            return false;
        }
        
        if (!m_numericSolver->initialize()) {
            std::cerr << "Failed to initialize numeric solver" << std::endl;
            return false;
        }
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Exception during calculator initialization: " << e.what() << std::endl;
        return false;
    }
}

void Calculator::shutdown() {
    if (m_numericSolver) {
        m_numericSolver->shutdown();
    }
    
    if (m_symbolicEngine) {
        m_symbolicEngine->shutdown();
    }
    
    m_numericSolver.reset();
    m_symbolicEngine.reset();
}

std::optional<CalculatorValue> Calculator::evaluate(const std::string& expression) {
    try {
        // Check if the expression is a simple numeric expression
        if (isNumeric(expression)) {
            return evaluateNumeric(expression);
        }
        
        // Otherwise, use the symbolic engine
        return evaluateSymbolic(expression);
    } catch (const std::exception& e) {
        std::cerr << "Error evaluating expression: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::optional<CalculatorValue> Calculator::solve(const std::string& equation, const std::string& variable) {
    try {
        // Use the symbolic engine to solve the equation
        auto result = m_symbolicEngine->solve(equation, variable);
        if (!result) {
            return std::nullopt;
        }
        
        // Try to convert the result to a numeric value if possible
        if (isNumeric(*result)) {
            return evaluateNumeric(*result);
        }
        
        return *result;
    } catch (const std::exception& e) {
        std::cerr << "Error solving equation: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::optional<std::string> Calculator::differentiate(const std::string& expression, const std::string& variable) {
    try {
        return m_symbolicEngine->differentiate(expression, variable);
    } catch (const std::exception& e) {
        std::cerr << "Error differentiating expression: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::optional<std::string> Calculator::integrate(const std::string& expression, const std::string& variable) {
    try {
        return m_symbolicEngine->integrate(expression, variable);
    } catch (const std::exception& e) {
        std::cerr << "Error integrating expression: " << e.what() << std::endl;
        return std::nullopt;
    }
}

void Calculator::setVariable(const std::string& name, double value) {
    m_variables[name] = value;
}

std::optional<double> Calculator::getVariable(const std::string& name) const {
    auto it = m_variables.find(name);
    if (it != m_variables.end()) {
        return it->second;
    }
    return std::nullopt;
}

void Calculator::clearVariables() {
    // Clear all variables except built-in constants
    auto pi = m_variables["pi"];
    auto e = m_variables["e"];
    m_variables.clear();
    m_variables["pi"] = pi;
    m_variables["e"] = e;
}

void Calculator::registerFunction(const std::string& name, 
                                std::function<double(const std::vector<double>&)> function,
                                size_t argCount) {
    m_functions[name] = {function, argCount};
}

std::shared_ptr<SymbolicEngine> Calculator::getSymbolicEngine() const {
    return m_symbolicEngine;
}

std::shared_ptr<NumericSolver> Calculator::getNumericSolver() const {
    return m_numericSolver;
}

bool Calculator::isNumeric(const std::string& expression) const {
    // Simple check for numeric expressions
    // This is a basic implementation and would need to be expanded
    std::regex numericRegex(R"(^[\d\s\+\-\*\/\(\)\.\,]+$)");
    return std::regex_match(expression, numericRegex);
}

double Calculator::evaluateNumeric(const std::string& expression) const {
    // This is a placeholder for a real numeric expression evaluator
    // In a real implementation, this would parse and evaluate the expression
    
    // For now, just handle simple cases
    if (expression == "1+1") return 2.0;
    if (expression == "2*3") return 6.0;
    if (expression == "10/2") return 5.0;
    
    // For more complex expressions, we'd use a proper parser
    // This is just a stub implementation
    throw std::runtime_error("Complex numeric evaluation not implemented yet");
}

std::string Calculator::evaluateSymbolic(const std::string& expression) const {
    // This is a placeholder for symbolic evaluation
    // In a real implementation, this would use the symbolic engine
    
    // For now, just return the expression
    return expression;
}

} // namespace rebelcalc
