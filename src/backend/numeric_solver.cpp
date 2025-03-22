#include "numeric_solver.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <complex>

// Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace rebelcalc {

// Private implementation class
class NumericSolver::Impl {
public:
    Impl() {}
    ~Impl() {}
    
    bool initialize() {
        // Initialize the numeric solver implementation
        return true;
    }
    
    void shutdown() {
        // Shutdown the numeric solver implementation
    }
    
    std::optional<std::vector<double>> solveLinearSystem(
        const std::vector<std::vector<double>>& coefficients,
        const std::vector<double>& constants) {
        
        // Check if the system is well-formed
        if (coefficients.empty() || constants.empty()) {
            return std::nullopt;
        }
        
        const size_t numEquations = coefficients.size();
        const size_t numVariables = coefficients[0].size();
        
        if (numEquations != constants.size()) {
            return std::nullopt;
        }
        
        // Check that all rows have the same number of columns
        for (const auto& row : coefficients) {
            if (row.size() != numVariables) {
                return std::nullopt;
            }
        }
        
        // Handle simple cases first
        if (numEquations == 1 && numVariables == 1) {
            // Simple one-variable equation: ax = b
            const double a = coefficients[0][0];
            const double b = constants[0];
            
            if (std::abs(a) < 1e-10) {
                // Singular matrix
                return std::nullopt;
            }
            
            return std::vector<double>{b / a};
        }
        
        if (numEquations == 2 && numVariables == 2) {
            // Simple 2x2 system:
            // a1*x + b1*y = c1
            // a2*x + b2*y = c2
            const double a1 = coefficients[0][0];
            const double b1 = coefficients[0][1];
            const double c1 = constants[0];
            
            const double a2 = coefficients[1][0];
            const double b2 = coefficients[1][1];
            const double c2 = constants[1];
            
            const double det = a1 * b2 - a2 * b1;
            if (std::abs(det) < 1e-10) {
                // Singular matrix
                return std::nullopt;
            }
            
            const double x = (c1 * b2 - c2 * b1) / det;
            const double y = (a1 * c2 - a2 * c1) / det;
            
            return std::vector<double>{x, y};
        }
        
        // For larger systems, use Gaussian elimination with partial pivoting
        // Create an augmented matrix [A|b]
        std::vector<std::vector<double>> augmentedMatrix(numEquations, std::vector<double>(numVariables + 1));
        for (size_t i = 0; i < numEquations; ++i) {
            for (size_t j = 0; j < numVariables; ++j) {
                augmentedMatrix[i][j] = coefficients[i][j];
            }
            augmentedMatrix[i][numVariables] = constants[i];
        }
        
        // Forward elimination with partial pivoting
        for (size_t k = 0; k < std::min(numEquations, numVariables); ++k) {
            // Find pivot
            size_t pivotRow = k;
            double pivotValue = std::abs(augmentedMatrix[k][k]);
            
            for (size_t i = k + 1; i < numEquations; ++i) {
                if (std::abs(augmentedMatrix[i][k]) > pivotValue) {
                    pivotRow = i;
                    pivotValue = std::abs(augmentedMatrix[i][k]);
                }
            }
            
            // Check for singularity
            if (pivotValue < 1e-10) {
                continue; // Skip this column
            }
            
            // Swap rows if necessary
            if (pivotRow != k) {
                augmentedMatrix[k].swap(augmentedMatrix[pivotRow]);
            }
            
            // Eliminate below
            for (size_t i = k + 1; i < numEquations; ++i) {
                double factor = augmentedMatrix[i][k] / augmentedMatrix[k][k];
                
                for (size_t j = k; j <= numVariables; ++j) {
                    augmentedMatrix[i][j] -= factor * augmentedMatrix[k][j];
                }
            }
        }
        
        // Check for inconsistent system
        for (size_t i = numVariables; i < numEquations; ++i) {
            bool allZeros = true;
            for (size_t j = 0; j < numVariables; ++j) {
                if (std::abs(augmentedMatrix[i][j]) >= 1e-10) {
                    allZeros = false;
                    break;
                }
            }
            
            if (allZeros && std::abs(augmentedMatrix[i][numVariables]) >= 1e-10) {
                // Inconsistent system (0 = non-zero)
                return std::nullopt;
            }
        }
        
        // Back substitution
        std::vector<double> solution(numVariables, 0.0);
        
        for (int i = static_cast<int>(numVariables) - 1; i >= 0; --i) {
            double sum = 0.0;
            
            for (size_t j = i + 1; j < numVariables; ++j) {
                sum += augmentedMatrix[i][j] * solution[j];
            }
            
            // Check for zero pivot
            if (std::abs(augmentedMatrix[i][i]) < 1e-10) {
                // Underdetermined system, set free variable to 0
                solution[i] = 0.0;
            } else {
                solution[i] = (augmentedMatrix[i][numVariables] - sum) / augmentedMatrix[i][i];
            }
        }
        
        return solution;
    }
    
    std::optional<std::vector<double>> findRoots(const std::vector<double>& coefficients) {
        // Check if the polynomial is well-formed
        if (coefficients.empty()) {
            return std::nullopt;
        }
        
        // Find the degree of the polynomial (ignoring leading zeros)
        size_t degree = coefficients.size() - 1;
        while (degree > 0 && std::abs(coefficients[coefficients.size() - 1 - degree]) < 1e-10) {
            --degree;
        }
        
        // Handle constant polynomial
        if (degree == 0) {
            // No roots for non-zero constant
            if (std::abs(coefficients[coefficients.size() - 1]) >= 1e-10) {
                return std::vector<double>{};
            }
            // All values are roots for zero polynomial
            return std::nullopt;
        }
        
        // Handle linear polynomial: ax + b = 0
        if (degree == 1) {
            const double a = coefficients[coefficients.size() - 2];
            const double b = coefficients[coefficients.size() - 1];
            
            if (std::abs(a) < 1e-10) {
                // Not a proper linear polynomial
                return std::nullopt;
            }
            
            return std::vector<double>{-b / a};
        }
        
        // Handle quadratic polynomial: ax^2 + bx + c = 0
        if (degree == 2) {
            const double a = coefficients[coefficients.size() - 3];
            const double b = coefficients[coefficients.size() - 2];
            const double c = coefficients[coefficients.size() - 1];
            
            if (std::abs(a) < 1e-10) {
                // Not a proper quadratic polynomial
                return findRoots(std::vector<double>{b, c});
            }
            
            const double discriminant = b * b - 4 * a * c;
            if (discriminant < 0) {
                // Complex roots (not supported in this implementation)
                return std::vector<double>{};
            }
            
            const double sqrtDiscriminant = std::sqrt(discriminant);
            const double x1 = (-b + sqrtDiscriminant) / (2 * a);
            const double x2 = (-b - sqrtDiscriminant) / (2 * a);
            
            // Return roots in ascending order
            if (x1 <= x2) {
                return std::vector<double>{x1, x2};
            } else {
                return std::vector<double>{x2, x1};
            }
        }
        
        // For cubic polynomials: ax^3 + bx^2 + cx + d = 0
        if (degree == 3) {
            const double a = coefficients[coefficients.size() - 4];
            const double b = coefficients[coefficients.size() - 3];
            const double c = coefficients[coefficients.size() - 2];
            const double d = coefficients[coefficients.size() - 1];
            
            if (std::abs(a) < 1e-10) {
                // Not a proper cubic polynomial
                return findRoots(std::vector<double>{b, c, d});
            }
            
            // Normalize to x^3 + px^2 + qx + r = 0
            const double p = b / a;
            const double q = c / a;
            const double r = d / a;
            
            // Cardano's method
            const double p_over_3 = p / 3.0;
            const double q_over_3 = q / 3.0;
            const double p_over_3_cubed = p_over_3 * p_over_3 * p_over_3;
            const double p_squared_over_3 = p_over_3 * p_over_3;
            
            const double Q = (3.0 * q - p * p) / 9.0;
            const double R = (9.0 * p * q - 27.0 * r - 2.0 * p * p * p) / 54.0;
            
            const double Q_cubed = Q * Q * Q;
            const double R_squared = R * R;
            
            const double discriminant = R_squared + Q_cubed;
            
            std::vector<double> roots;
            
            if (discriminant > 0) {
                // One real root and two complex conjugate roots
                const double S = std::cbrt(R + std::sqrt(discriminant));
                const double T = std::cbrt(R - std::sqrt(discriminant));
                
                const double x = S + T - p_over_3;
                roots.push_back(x);
            } else if (std::abs(discriminant) < 1e-10) {
                // All roots are real, at least two are equal
                const double S = std::cbrt(R);
                
                const double x1 = 2.0 * S - p_over_3;
                const double x2 = -S - p_over_3;
                
                roots.push_back(x1);
                roots.push_back(x2);
                roots.push_back(x2); // Repeated root
            } else {
                // All roots are real and distinct
                const double theta = std::acos(R / std::sqrt(-Q_cubed));
                const double sqrt_minus_Q = std::sqrt(-Q);
                
                const double x1 = 2.0 * sqrt_minus_Q * std::cos(theta / 3.0) - p_over_3;
                const double x2 = 2.0 * sqrt_minus_Q * std::cos((theta + 2.0 * M_PI) / 3.0) - p_over_3;
                const double x3 = 2.0 * sqrt_minus_Q * std::cos((theta + 4.0 * M_PI) / 3.0) - p_over_3;
                
                roots.push_back(x1);
                roots.push_back(x2);
                roots.push_back(x3);
            }
            
            // Sort the roots
            std::sort(roots.begin(), roots.end());
            
            return roots;
        }
        
        // For higher-degree polynomials, we'll use a simpler approach
        // The Durand-Kerner method is complex and requires complex number support
        // Instead, we'll use a bisection method to find real roots in intervals
        
        // For polynomials of degree > 3, we'll just return an empty vector for now
        // In a real implementation, we would use a more sophisticated method
        
        // This is a placeholder for a more sophisticated root-finding algorithm
        std::vector<double> realRoots;
        
        // Sort the roots
        std::sort(realRoots.begin(), realRoots.end());
        
        return realRoots;
    }
    
    std::optional<double> findMinimum(
        const std::function<double(double)>& function,
        double initialGuess,
        double lowerBound,
        double upperBound) {
        
        // Implementation of the Golden Section Search method
        // This is a more efficient method than the simple grid search
        
        // Golden ratio
        const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
        const double resphi = 2.0 - phi; // 1/phi
        
        // Tolerance for convergence
        const double tolerance = 1e-10;
        
        // Initialize the interval
        double a = lowerBound;
        double b = upperBound;
        
        // Initial points
        double c = b - resphi * (b - a);
        double d = a + resphi * (b - a);
        
        // Function values at the points
        double fc = function(c);
        double fd = function(d);
        
        // Iterate until the interval is small enough
        while (std::abs(b - a) > tolerance) {
            if (fc < fd) {
                // Minimum is in [a, d]
                b = d;
                d = c;
                fd = fc;
                c = b - resphi * (b - a);
                fc = function(c);
            } else {
                // Minimum is in [c, b]
                a = c;
                c = d;
                fc = fd;
                d = a + resphi * (b - a);
                fd = function(d);
            }
        }
        
        // Return the midpoint of the final interval
        return (a + b) / 2.0;
    }
    
    std::optional<double> findMaximum(
        const std::function<double(double)>& function,
        double initialGuess,
        double lowerBound,
        double upperBound) {
        
        // Convert maximization to minimization
        auto negatedFunction = [&function](double x) {
            return -function(x);
        };
        
        auto result = findMinimum(negatedFunction, initialGuess, lowerBound, upperBound);
        return result;
    }
    
    std::optional<double> integrate(
        const std::function<double(double)>& function,
        double lowerBound,
        double upperBound) {
        
        // This is a placeholder for a real numerical integrator
        // In a real implementation, this would use a numerical library
        
        // For now, just use a simple trapezoidal rule
        const int numPoints = 1000;
        const double step = (upperBound - lowerBound) / numPoints;
        
        double sum = 0.5 * (function(lowerBound) + function(upperBound));
        
        for (int i = 1; i < numPoints; ++i) {
            const double x = lowerBound + i * step;
            sum += function(x);
        }
        
        return sum * step;
    }
    
    std::optional<double> differentiate(
        const std::function<double(double)>& function,
        double point) {
        
        // This is a placeholder for a real numerical differentiator
        // In a real implementation, this would use a numerical library
        
        // For now, just use a simple central difference
        const double h = 1e-6;
        return (function(point + h) - function(point - h)) / (2 * h);
    }
};

// NumericSolver implementation

NumericSolver::NumericSolver() : m_impl(std::make_unique<Impl>()) {
}

NumericSolver::~NumericSolver() {
    shutdown();
}

bool NumericSolver::initialize() {
    try {
        return m_impl->initialize();
    } catch (const std::exception& e) {
        std::cerr << "Exception during numeric solver initialization: " << e.what() << std::endl;
        return false;
    }
}

void NumericSolver::shutdown() {
    if (m_impl) {
        m_impl->shutdown();
    }
}

std::optional<std::vector<double>> NumericSolver::solveLinearSystem(
    const std::vector<std::vector<double>>& coefficients,
    const std::vector<double>& constants) {
    
    try {
        return m_impl->solveLinearSystem(coefficients, constants);
    } catch (const std::exception& e) {
        std::cerr << "Error solving linear system: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::optional<std::vector<double>> NumericSolver::findRoots(const std::vector<double>& coefficients) {
    try {
        return m_impl->findRoots(coefficients);
    } catch (const std::exception& e) {
        std::cerr << "Error finding roots: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::optional<double> NumericSolver::findMinimum(
    const std::function<double(double)>& function,
    double initialGuess,
    double lowerBound,
    double upperBound) {
    
    try {
        return m_impl->findMinimum(function, initialGuess, lowerBound, upperBound);
    } catch (const std::exception& e) {
        std::cerr << "Error finding minimum: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::optional<double> NumericSolver::findMaximum(
    const std::function<double(double)>& function,
    double initialGuess,
    double lowerBound,
    double upperBound) {
    
    try {
        return m_impl->findMaximum(function, initialGuess, lowerBound, upperBound);
    } catch (const std::exception& e) {
        std::cerr << "Error finding maximum: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::optional<double> NumericSolver::integrate(
    const std::function<double(double)>& function,
    double lowerBound,
    double upperBound) {
    
    try {
        return m_impl->integrate(function, lowerBound, upperBound);
    } catch (const std::exception& e) {
        std::cerr << "Error integrating function: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::optional<double> NumericSolver::differentiate(
    const std::function<double(double)>& function,
    double point) {
    
    try {
        return m_impl->differentiate(function, point);
    } catch (const std::exception& e) {
        std::cerr << "Error differentiating function: " << e.what() << std::endl;
        return std::nullopt;
    }
}

} // namespace rebelcalc
