#pragma once

#include <string>
#include <vector>
#include <optional>
#include <memory>
#include <functional>

namespace rebelcalc {

/**
 * Class for handling numerical computations and solving
 */
class NumericSolver {
public:
    /**
     * Constructor
     */
    NumericSolver();
    
    /**
     * Destructor
     */
    ~NumericSolver();
    
    /**
     * Initialize the numeric solver
     * @return true if initialization was successful, false otherwise
     */
    bool initialize();
    
    /**
     * Shutdown the numeric solver
     */
    void shutdown();
    
    /**
     * Solve a system of linear equations
     * @param coefficients The coefficient matrix
     * @param constants The constant vector
     * @return The solution vector, or nullopt if solving failed
     */
    std::optional<std::vector<double>> solveLinearSystem(
        const std::vector<std::vector<double>>& coefficients,
        const std::vector<double>& constants);
    
    /**
     * Find roots of a polynomial
     * @param coefficients The polynomial coefficients (highest degree first)
     * @return The roots of the polynomial, or nullopt if root finding failed
     */
    std::optional<std::vector<double>> findRoots(const std::vector<double>& coefficients);
    
    /**
     * Find the minimum of a function
     * @param function The function to minimize
     * @param initialGuess The initial guess
     * @param lowerBound The lower bound for the search
     * @param upperBound The upper bound for the search
     * @return The minimum point, or nullopt if minimization failed
     */
    std::optional<double> findMinimum(
        const std::function<double(double)>& function,
        double initialGuess,
        double lowerBound,
        double upperBound);
    
    /**
     * Find the maximum of a function
     * @param function The function to maximize
     * @param initialGuess The initial guess
     * @param lowerBound The lower bound for the search
     * @param upperBound The upper bound for the search
     * @return The maximum point, or nullopt if maximization failed
     */
    std::optional<double> findMaximum(
        const std::function<double(double)>& function,
        double initialGuess,
        double lowerBound,
        double upperBound);
    
    /**
     * Numerically integrate a function
     * @param function The function to integrate
     * @param lowerBound The lower bound of integration
     * @param upperBound The upper bound of integration
     * @return The integral value, or nullopt if integration failed
     */
    std::optional<double> integrate(
        const std::function<double(double)>& function,
        double lowerBound,
        double upperBound);
    
    /**
     * Numerically differentiate a function at a point
     * @param function The function to differentiate
     * @param point The point at which to differentiate
     * @return The derivative value, or nullopt if differentiation failed
     */
    std::optional<double> differentiate(
        const std::function<double(double)>& function,
        double point);

private:
    // Implementation details
    class Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace rebelcalc
