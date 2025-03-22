#include "calculator.h"
#include "symbolic_engine.h"
#include "numeric_solver.h"
#include "complex.h"
#include "statistics.h"

#include <iostream>
#include <cmath>
#include <regex>
#include <stdexcept>
#include <complex>

// Define M_PI and M_E if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

namespace rebelcalc {

Calculator::Calculator() 
    : m_symbolicEngine(nullptr), m_numericSolver(nullptr), m_statistics(nullptr), m_unitConverter(nullptr) {
    // Initialize built-in variables
    m_variables["pi"] = M_PI;
    m_variables["e"] = M_E;
    
    // Initialize built-in complex variables
    m_complexVariables["i"] = Complex(0.0, 1.0);
    m_complexVariables["j"] = Complex(0.0, 1.0); // Alternative notation for i
    
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
        m_statistics = std::make_shared<Statistics>();
        m_unitConverter = std::make_shared<UnitConverter>();
        
        // Initialize components
        if (!m_symbolicEngine->initialize()) {
            std::cerr << "Failed to initialize symbolic engine" << std::endl;
            return false;
        }
        
        if (!m_numericSolver->initialize()) {
            std::cerr << "Failed to initialize numeric solver" << std::endl;
            return false;
        }
        
        if (!m_statistics->initialize()) {
            std::cerr << "Failed to initialize statistics engine" << std::endl;
            return false;
        }
        
        if (!m_unitConverter->initialize()) {
            std::cerr << "Failed to initialize unit converter" << std::endl;
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
    
    if (m_statistics) {
        m_statistics->shutdown();
    }
    
    if (m_unitConverter) {
        m_unitConverter->shutdown();
    }
    
    m_numericSolver.reset();
    m_symbolicEngine.reset();
    m_statistics.reset();
    m_unitConverter.reset();
}

std::optional<CalculatorValue> Calculator::evaluate(const std::string& expression) {
    try {
        // Check if the expression is a variable assignment
        std::regex assignmentRegex(R"(^\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*([+\-*/]?=)\s*.+$)");
        if (std::regex_match(expression, assignmentRegex)) {
            auto result = processVariableAssignment(expression);
            if (result) {
                return *result;
            }
            return std::nullopt;
        }
        
        // Check if the expression is a complex number literal
        auto complexResult = parseComplex(expression);
        if (complexResult) {
            return *complexResult;
        }
        
        // Check if the expression involves complex numbers
        if (isComplex(expression)) {
            return evaluateComplex(expression);
        }
        
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

bool Calculator::setVariable(const std::string& name, double value) {
    // Check if the variable name is valid
    if (name.empty()) {
        return false;
    }
    
    // Check if the first character is a letter or underscore
    if (!std::isalpha(name[0]) && name[0] != '_') {
        return false;
    }
    
    // Check if the rest of the characters are letters, digits, or underscores
    for (size_t i = 1; i < name.length(); ++i) {
        if (!std::isalnum(name[i]) && name[i] != '_') {
            return false;
        }
    }
    
    // Check if the variable name is a reserved keyword
    if (name == "pi" || name == "e") {
        // Allow overriding built-in constants, but warn about it
        std::cerr << "Warning: Overriding built-in constant '" << name << "'" << std::endl;
    }
    
    // Check if the variable name is a reserved complex constant
    if (name == "i" || name == "j") {
        std::cerr << "Warning: Overriding built-in complex constant '" << name << "'" << std::endl;
        // Remove from complex variables if it exists
        m_complexVariables.erase(name);
    }
    
    // Set the variable
    m_variables[name] = value;
    return true;
}

bool Calculator::setComplexVariable(const std::string& name, const Complex& value) {
    // Check if the variable name is valid
    if (name.empty()) {
        return false;
    }
    
    // Check if the first character is a letter or underscore
    if (!std::isalpha(name[0]) && name[0] != '_') {
        return false;
    }
    
    // Check if the rest of the characters are letters, digits, or underscores
    for (size_t i = 1; i < name.length(); ++i) {
        if (!std::isalnum(name[i]) && name[i] != '_') {
            return false;
        }
    }
    
    // Check if the variable name is a reserved complex constant
    if (name == "i" || name == "j") {
        // Allow overriding built-in constants, but warn about it
        std::cerr << "Warning: Overriding built-in complex constant '" << name << "'" << std::endl;
    }
    
    // Remove from regular variables if it exists
    m_variables.erase(name);
    
    // Set the complex variable
    m_complexVariables[name] = value;
    return true;
}

std::optional<double> Calculator::getVariable(const std::string& name) const {
    auto it = m_variables.find(name);
    if (it != m_variables.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<Complex> Calculator::getComplexVariable(const std::string& name) const {
    auto it = m_complexVariables.find(name);
    if (it != m_complexVariables.end()) {
        return it->second;
    }
    return std::nullopt;
}

bool Calculator::hasVariable(const std::string& name) const {
    return m_variables.find(name) != m_variables.end() || 
           m_complexVariables.find(name) != m_complexVariables.end();
}

std::unordered_map<std::string, double> Calculator::getAllVariables() const {
    return m_variables;
}

void Calculator::clearVariables() {
    // Clear all variables except built-in constants
    auto pi = m_variables["pi"];
    auto e = m_variables["e"];
    auto i = m_complexVariables["i"];
    auto j = m_complexVariables["j"];
    
    m_variables.clear();
    m_complexVariables.clear();
    
    m_variables["pi"] = pi;
    m_variables["e"] = e;
    m_complexVariables["i"] = i;
    m_complexVariables["j"] = j;
}

std::optional<double> Calculator::processVariableAssignment(const std::string& expression) {
    // Parse the expression to extract the variable name, operator, and value
    std::regex assignmentRegex(R"(^\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*([+\-*/]?=)\s*(.+)$)");
    std::smatch matches;
    
    if (!std::regex_match(expression, matches, assignmentRegex)) {
        return std::nullopt;
    }
    
    std::string variableName = matches[1].str();
    std::string assignmentOp = matches[2].str();
    std::string valueExpr = matches[3].str();
    
    // Evaluate the right-hand side expression
    auto valueResult = evaluate(valueExpr);
    if (!valueResult || !std::holds_alternative<double>(*valueResult)) {
        return std::nullopt;
    }
    
    double value = std::get<double>(*valueResult);
    
    // Process the assignment based on the operator
    if (assignmentOp == "=") {
        // Simple assignment
        if (!setVariable(variableName, value)) {
            return std::nullopt;
        }
    } else {
        // Compound assignment (+=, -=, *=, /=)
        auto currentValue = getVariable(variableName);
        if (!currentValue) {
            // Variable doesn't exist, create it with a default value of 0
            if (!setVariable(variableName, 0.0)) {
                return std::nullopt;
            }
            currentValue = 0.0;
        }
        
        double newValue = 0.0;
        
        if (assignmentOp == "+=") {
            newValue = *currentValue + value;
        } else if (assignmentOp == "-=") {
            newValue = *currentValue - value;
        } else if (assignmentOp == "*=") {
            newValue = *currentValue * value;
        } else if (assignmentOp == "/=") {
            if (value == 0.0) {
                return std::nullopt; // Division by zero
            }
            newValue = *currentValue / value;
        } else {
            return std::nullopt; // Unknown operator
        }
        
        if (!setVariable(variableName, newValue)) {
            return std::nullopt;
        }
        
        value = newValue; // Return the new value
    }
    
    return value;
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

std::shared_ptr<Statistics> Calculator::getStatistics() const {
    return m_statistics;
}

std::shared_ptr<UnitConverter> Calculator::getUnitConverter() const {
    return m_unitConverter;
}

std::optional<double> Calculator::convertUnit(double value, const std::string& fromUnit, const std::string& toUnit) {
    try {
        if (!m_unitConverter) {
            std::cerr << "Unit converter not initialized" << std::endl;
            return std::nullopt;
        }
        
        return m_unitConverter->convert(value, fromUnit, toUnit);
    } catch (const std::exception& e) {
        std::cerr << "Error converting units: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::vector<std::string> Calculator::getUnitsForQuantity(const std::string& quantity) const {
    if (!m_unitConverter) {
        std::cerr << "Unit converter not initialized" << std::endl;
        return {};
    }
    
    return m_unitConverter->getUnitsForQuantity(quantity);
}

std::vector<std::string> Calculator::getAllQuantities() const {
    if (!m_unitConverter) {
        std::cerr << "Unit converter not initialized" << std::endl;
        return {};
    }
    
    return m_unitConverter->getAllQuantities();
}

// Helper class for parsing and evaluating numeric expressions
class ExpressionParser {
public:
    ExpressionParser(const std::unordered_map<std::string, double>& variables,
                    const std::unordered_map<std::string, Calculator::FunctionInfo>& functions)
        : m_variables(variables), m_functions(functions), m_pos(0) {}
    
    double parse(const std::string& expression) {
        m_expression = expression;
        m_pos = 0;
        
        // Skip whitespace
        skipWhitespace();
        
        // Parse the expression
        double result = parseExpression();
        
        // Check if we've consumed the entire expression
        if (m_pos < m_expression.length()) {
            throw std::runtime_error("Unexpected character: " + std::string(1, m_expression[m_pos]));
        }
        
        return result;
    }
    
private:
    std::string m_expression;
    size_t m_pos;
    const std::unordered_map<std::string, double>& m_variables;
    const std::unordered_map<std::string, Calculator::FunctionInfo>& m_functions;
    
    // Skip whitespace characters
    void skipWhitespace() {
        while (m_pos < m_expression.length() && std::isspace(m_expression[m_pos])) {
            ++m_pos;
        }
    }
    
    // Check if the current character matches the expected character
    bool match(char expected) {
        if (m_pos < m_expression.length() && m_expression[m_pos] == expected) {
            ++m_pos;
            skipWhitespace();
            return true;
        }
        return false;
    }
    
    // Parse an expression (addition and subtraction)
    double parseExpression() {
        double left = parseTerm();
        
        while (true) {
            if (match('+')) {
                left += parseTerm();
            } else if (match('-')) {
                left -= parseTerm();
            } else {
                break;
            }
        }
        
        return left;
    }
    
    // Parse a term (multiplication and division)
    double parseTerm() {
        double left = parseFactor();
        
        while (true) {
            if (match('*')) {
                left *= parseFactor();
            } else if (match('/')) {
                double right = parseFactor();
                if (right == 0.0) {
                    throw std::runtime_error("Division by zero");
                }
                left /= right;
            } else if (match('^') || match('p')) { // ^ or p for power
                double right = parseFactor();
                left = std::pow(left, right);
            } else {
                break;
            }
        }
        
        return left;
    }
    
    // Parse a factor (numbers, variables, functions, parentheses, unary operators)
    double parseFactor() {
        // Unary plus
        if (match('+')) {
            return parseFactor();
        }
        
        // Unary minus
        if (match('-')) {
            return -parseFactor();
        }
        
        // Parentheses
        if (match('(')) {
            double result = parseExpression();
            if (!match(')')) {
                throw std::runtime_error("Expected closing parenthesis");
            }
            return result;
        }
        
        // Number
        if (m_pos < m_expression.length() && (std::isdigit(m_expression[m_pos]) || m_expression[m_pos] == '.')) {
            return parseNumber();
        }
        
        // Variable or function
        if (m_pos < m_expression.length() && (std::isalpha(m_expression[m_pos]) || m_expression[m_pos] == '_')) {
            std::string identifier = parseIdentifier();
            
            // Check if it's a function call
            if (m_pos < m_expression.length() && m_expression[m_pos] == '(') {
                return parseFunction(identifier);
            }
            
            // Otherwise, it's a variable
            auto it = m_variables.find(identifier);
            if (it != m_variables.end()) {
                return it->second;
            }
            
            throw std::runtime_error("Unknown variable: " + identifier);
        }
        
        throw std::runtime_error("Unexpected character: " + std::string(1, m_expression[m_pos]));
    }
    
    // Parse a number
    double parseNumber() {
        size_t start = m_pos;
        
        // Parse integer part
        while (m_pos < m_expression.length() && std::isdigit(m_expression[m_pos])) {
            ++m_pos;
        }
        
        // Parse decimal part
        if (m_pos < m_expression.length() && m_expression[m_pos] == '.') {
            ++m_pos;
            
            while (m_pos < m_expression.length() && std::isdigit(m_expression[m_pos])) {
                ++m_pos;
            }
        }
        
        // Parse scientific notation
        if (m_pos < m_expression.length() && (m_expression[m_pos] == 'e' || m_expression[m_pos] == 'E')) {
            ++m_pos;
            
            if (m_pos < m_expression.length() && (m_expression[m_pos] == '+' || m_expression[m_pos] == '-')) {
                ++m_pos;
            }
            
            if (m_pos >= m_expression.length() || !std::isdigit(m_expression[m_pos])) {
                throw std::runtime_error("Invalid scientific notation");
            }
            
            while (m_pos < m_expression.length() && std::isdigit(m_expression[m_pos])) {
                ++m_pos;
            }
        }
        
        skipWhitespace();
        
        // Convert the substring to a double
        std::string numberStr = m_expression.substr(start, m_pos - start);
        try {
            return std::stod(numberStr);
        } catch (const std::exception& e) {
            throw std::runtime_error("Invalid number: " + numberStr);
        }
    }
    
    // Parse an identifier (variable or function name)
    std::string parseIdentifier() {
        size_t start = m_pos;
        
        // First character must be a letter or underscore
        if (m_pos < m_expression.length() && (std::isalpha(m_expression[m_pos]) || m_expression[m_pos] == '_')) {
            ++m_pos;
        } else {
            throw std::runtime_error("Invalid identifier");
        }
        
        // Subsequent characters can be letters, digits, or underscores
        while (m_pos < m_expression.length() && 
               (std::isalnum(m_expression[m_pos]) || m_expression[m_pos] == '_')) {
            ++m_pos;
        }
        
        skipWhitespace();
        
        return m_expression.substr(start, m_pos - start);
    }
    
    // Parse a function call
    double parseFunction(const std::string& name) {
        // Check if the function exists
        auto it = m_functions.find(name);
        if (it == m_functions.end()) {
            throw std::runtime_error("Unknown function: " + name);
        }
        
        const auto& functionInfo = it->second;
        
        // Parse the opening parenthesis
        if (!match('(')) {
            throw std::runtime_error("Expected opening parenthesis after function name");
        }
        
        // Parse the arguments
        std::vector<double> args;
        
        // Handle empty argument list
        if (match(')')) {
            if (functionInfo.argCount != 0) {
                throw std::runtime_error("Function " + name + " expects " + 
                                        std::to_string(functionInfo.argCount) + " arguments");
            }
            return functionInfo.function(args);
        }
        
        // Parse the first argument
        args.push_back(parseExpression());
        
        // Parse additional arguments
        while (match(',')) {
            args.push_back(parseExpression());
        }
        
        // Parse the closing parenthesis
        if (!match(')')) {
            throw std::runtime_error("Expected closing parenthesis after function arguments");
        }
        
        // Check if the number of arguments matches the expected count
        if (args.size() != functionInfo.argCount) {
            throw std::runtime_error("Function " + name + " expects " + 
                                    std::to_string(functionInfo.argCount) + " arguments, but got " + 
                                    std::to_string(args.size()));
        }
        
        // Call the function
        return functionInfo.function(args);
    }
};

bool Calculator::isNumeric(const std::string& expression) const {
    // Check if the expression contains only numeric characters, operators, and parentheses
    // This is a more comprehensive check than the previous regex
    for (char c : expression) {
        if (!std::isdigit(c) && !std::isspace(c) && 
            c != '+' && c != '-' && c != '*' && c != '/' && c != '^' && 
            c != '(' && c != ')' && c != '.' && c != ',') {
            
            // Check if it's a variable or function name
            if (std::isalpha(c) || c == '_') {
                // Check if it's a known variable
                bool isKnownVariable = false;
                for (const auto& [name, _] : m_variables) {
                    if (expression.find(name) != std::string::npos) {
                        isKnownVariable = true;
                        break;
                    }
                }
                
                // Check if it's a known function
                bool isKnownFunction = false;
                for (const auto& [name, _] : m_functions) {
                    if (expression.find(name) != std::string::npos) {
                        isKnownFunction = true;
                        break;
                    }
                }
                
                // Check if it's a complex variable (i or j)
                bool isComplexConstant = false;
                if (expression.find("i") != std::string::npos || 
                    expression.find("j") != std::string::npos) {
                    isComplexConstant = true;
                }
                
                if (!isKnownVariable && !isKnownFunction && !isComplexConstant) {
                    return false;
                }
            } else {
                return false;
            }
        }
    }
    
    // Check if the expression contains complex variables
    for (const auto& [name, _] : m_complexVariables) {
        if (expression.find(name) != std::string::npos) {
            return false; // Contains complex variables, not purely numeric
        }
    }
    
    return true;
}

bool Calculator::isComplex(const std::string& expression) const {
    // Check if the expression contains complex variables or the imaginary unit
    for (const auto& [name, _] : m_complexVariables) {
        if (expression.find(name) != std::string::npos) {
            return true;
        }
    }
    
    // Check for explicit complex number notation (e.g., 3+4i, 2.5i)
    std::regex complexRegex(R"(.*[0-9.][ij].*|.*[+-][^+-]*[ij].*)");
    if (std::regex_match(expression, complexRegex)) {
        return true;
    }
    
    return false;
}

std::optional<Complex> Calculator::parseComplex(const std::string& expression) const {
    // Try to parse a complex number literal
    std::string trimmed = expression;
    
    // Remove whitespace
    trimmed.erase(std::remove_if(trimmed.begin(), trimmed.end(), ::isspace), trimmed.end());
    
    // Check for simple imaginary number (e.g., "i", "j")
    if (trimmed == "i" || trimmed == "j") {
        return Complex(0.0, 1.0);
    }
    
    // Check for simple imaginary number with coefficient (e.g., "3i", "4.5j")
    std::regex imagRegex(R"(([+-]?[0-9]*\.?[0-9]+)([ij]))");
    std::smatch imagMatch;
    if (std::regex_match(trimmed, imagMatch, imagRegex)) {
        double imag = std::stod(imagMatch[1].str());
        return Complex(0.0, imag);
    }
    
    // Check for complex number in form a+bi or a-bi
    std::regex complexRegex(R"(([+-]?[0-9]*\.?[0-9]+)([+-])([0-9]*\.?[0-9]+)([ij]))");
    std::smatch complexMatch;
    if (std::regex_match(trimmed, complexMatch, complexRegex)) {
        double real = std::stod(complexMatch[1].str());
        double imag = std::stod(complexMatch[3].str());
        if (complexMatch[2].str() == "-") {
            imag = -imag;
        }
        return Complex(real, imag);
    }
    
    return std::nullopt;
}

Complex Calculator::evaluateComplex(const std::string& expression) const {
    // This is a placeholder for complex expression evaluation
    // In a real implementation, this would use a more sophisticated parser
    
    // For now, just try to parse it as a complex literal
    auto result = parseComplex(expression);
    if (result) {
        return *result;
    }
    
    // If it's not a complex literal, try to evaluate it as a numeric expression
    // and convert to complex
    try {
        double value = evaluateNumeric(expression);
        return Complex(value);
    } catch (const std::exception&) {
        // If that fails, throw an error
        throw std::runtime_error("Failed to evaluate complex expression: " + expression);
    }
}

double Calculator::evaluateNumeric(const std::string& expression) const {
    try {
        // Create a parser with the current variables and functions
        ExpressionParser parser(m_variables, m_functions);
        
        // Parse and evaluate the expression
        return parser.parse(expression);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error evaluating expression: " + std::string(e.what()));
    }
}

std::string Calculator::evaluateSymbolic(const std::string& expression) const {
    // This is a placeholder for symbolic evaluation
    // In a real implementation, this would use the symbolic engine
    
    // For now, just return the expression
    return expression;
}

} // namespace rebelcalc
