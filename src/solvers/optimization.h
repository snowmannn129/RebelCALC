#ifndef REBELCALC_SOLVERS_OPTIMIZATION_H
#define REBELCALC_SOLVERS_OPTIMIZATION_H

#include <vector>
#include <string>
#include <map>
#include <functional>
#include <memory>
#include <stdexcept>
#include <limits>
#include <random>
#include "../backend/matrix.h"

namespace RebelCalc {
namespace Solvers {

/**
 * @brief Enum for optimization algorithm types
 */
enum class OptimizationAlgorithm {
    GRADIENT_DESCENT,       // Gradient descent algorithm
    NEWTON,                 // Newton's method
    BFGS,                   // Broyden-Fletcher-Goldfarb-Shanno algorithm
    CONJUGATE_GRADIENT,     // Conjugate gradient method
    SIMPLEX,                // Simplex method (Nelder-Mead)
    SIMULATED_ANNEALING,    // Simulated annealing
    GENETIC_ALGORITHM,      // Genetic algorithm
    PARTICLE_SWARM,         // Particle swarm optimization
    DIFFERENTIAL_EVOLUTION, // Differential evolution
    BAYESIAN_OPTIMIZATION   // Bayesian optimization
};

/**
 * @brief Enum for constraint types
 */
enum class ConstraintType {
    EQUALITY,    // Equality constraint: g(x) = 0
    INEQUALITY   // Inequality constraint: g(x) <= 0
};

/**
 * @brief Class for a constraint in optimization
 */
class Constraint {
public:
    /**
     * @brief Constructor for a constraint
     * @param type Constraint type
     * @param function Constraint function
     * @param gradient Gradient of the constraint function (optional)
     */
    Constraint(ConstraintType type, std::function<double(const std::vector<double>&)> function,
              std::function<std::vector<double>(const std::vector<double>&)> gradient = nullptr);

    /**
     * @brief Get the constraint type
     * @return Constraint type
     */
    ConstraintType getType() const;

    /**
     * @brief Evaluate the constraint function
     * @param x Point at which to evaluate the constraint
     * @return Value of the constraint function
     */
    double evaluate(const std::vector<double>& x) const;

    /**
     * @brief Check if the constraint is satisfied
     * @param x Point at which to check the constraint
     * @param tolerance Tolerance for constraint satisfaction
     * @return true if the constraint is satisfied, false otherwise
     */
    bool isSatisfied(const std::vector<double>& x, double tolerance = 1e-6) const;

    /**
     * @brief Evaluate the gradient of the constraint function
     * @param x Point at which to evaluate the gradient
     * @return Gradient of the constraint function
     */
    std::vector<double> evaluateGradient(const std::vector<double>& x) const;

    /**
     * @brief Check if the gradient is available
     * @return true if the gradient is available, false otherwise
     */
    bool hasGradient() const;

private:
    ConstraintType m_type;
    std::function<double(const std::vector<double>&)> m_function;
    std::function<std::vector<double>(const std::vector<double>&)> m_gradient;
};

/**
 * @brief Class for a bound constraint in optimization
 */
class BoundConstraint {
public:
    /**
     * @brief Constructor for a bound constraint
     * @param lowerBounds Lower bounds for each variable
     * @param upperBounds Upper bounds for each variable
     */
    BoundConstraint(const std::vector<double>& lowerBounds, const std::vector<double>& upperBounds);

    /**
     * @brief Get the lower bounds
     * @return Lower bounds
     */
    const std::vector<double>& getLowerBounds() const;

    /**
     * @brief Get the upper bounds
     * @return Upper bounds
     */
    const std::vector<double>& getUpperBounds() const;

    /**
     * @brief Check if a point is within the bounds
     * @param x Point to check
     * @return true if the point is within the bounds, false otherwise
     */
    bool isWithinBounds(const std::vector<double>& x) const;

    /**
     * @brief Project a point onto the bounds
     * @param x Point to project
     * @return Projected point
     */
    std::vector<double> project(const std::vector<double>& x) const;

private:
    std::vector<double> m_lowerBounds;
    std::vector<double> m_upperBounds;
};

/**
 * @brief Class for an optimization problem
 */
class OptimizationProblem {
public:
    /**
     * @brief Constructor for an optimization problem
     * @param objective Objective function to minimize
     * @param dimension Dimension of the problem (number of variables)
     */
    OptimizationProblem(std::function<double(const std::vector<double>&)> objective, int dimension);

    /**
     * @brief Set the gradient of the objective function
     * @param gradient Gradient of the objective function
     */
    void setGradient(std::function<std::vector<double>(const std::vector<double>&)> gradient);

    /**
     * @brief Set the Hessian of the objective function
     * @param hessian Hessian of the objective function
     */
    void setHessian(std::function<rebelcalc::Matrix(const std::vector<double>&)> hessian);

    /**
     * @brief Add a constraint to the problem
     * @param constraint Constraint to add
     */
    void addConstraint(const Constraint& constraint);

    /**
     * @brief Set the bound constraints
     * @param bounds Bound constraints
     */
    void setBounds(const BoundConstraint& bounds);

    /**
     * @brief Get the objective function
     * @return Objective function
     */
    std::function<double(const std::vector<double>&)> getObjective() const;

    /**
     * @brief Get the gradient of the objective function
     * @return Gradient of the objective function
     */
    std::function<std::vector<double>(const std::vector<double>&)> getGradient() const;

    /**
     * @brief Get the Hessian of the objective function
     * @return Hessian of the objective function
     */
    std::function<rebelcalc::Matrix(const std::vector<double>&)> getHessian() const;

    /**
     * @brief Get the constraints
     * @return Constraints
     */
    const std::vector<Constraint>& getConstraints() const;

    /**
     * @brief Get the bound constraints
     * @return Bound constraints
     */
    const BoundConstraint* getBounds() const;

    /**
     * @brief Get the dimension of the problem
     * @return Dimension of the problem
     */
    int getDimension() const;

    /**
     * @brief Check if the gradient is available
     * @return true if the gradient is available, false otherwise
     */
    bool hasGradient() const;

    /**
     * @brief Check if the Hessian is available
     * @return true if the Hessian is available, false otherwise
     */
    bool hasHessian() const;

    /**
     * @brief Check if the problem has constraints
     * @return true if the problem has constraints, false otherwise
     */
    bool hasConstraints() const;

    /**
     * @brief Check if the problem has bound constraints
     * @return true if the problem has bound constraints, false otherwise
     */
    bool hasBounds() const;

    /**
     * @brief Evaluate the objective function
     * @param x Point at which to evaluate the objective
     * @return Value of the objective function
     */
    double evaluateObjective(const std::vector<double>& x) const;

    /**
     * @brief Evaluate the gradient of the objective function
     * @param x Point at which to evaluate the gradient
     * @return Gradient of the objective function
     */
    std::vector<double> evaluateGradient(const std::vector<double>& x) const;

    /**
     * @brief Evaluate the Hessian of the objective function
     * @param x Point at which to evaluate the Hessian
     * @return Hessian of the objective function
     */
    rebelcalc::Matrix evaluateHessian(const std::vector<double>& x) const;

    /**
     * @brief Check if a point is feasible
     * @param x Point to check
     * @param tolerance Tolerance for constraint satisfaction
     * @return true if the point is feasible, false otherwise
     */
    bool isFeasible(const std::vector<double>& x, double tolerance = 1e-6) const;

private:
    std::function<double(const std::vector<double>&)> m_objective;
    std::function<std::vector<double>(const std::vector<double>&)> m_gradient;
    std::function<rebelcalc::Matrix(const std::vector<double>&)> m_hessian;
    std::vector<Constraint> m_constraints;
    std::unique_ptr<BoundConstraint> m_bounds;
    int m_dimension;
};

/**
 * @brief Class for optimization solver options
 */
class OptimizationOptions {
public:
    /**
     * @brief Constructor for optimization options
     */
    OptimizationOptions();

    /**
     * @brief Set the maximum number of iterations
     * @param maxIterations Maximum number of iterations
     */
    void setMaxIterations(int maxIterations);

    /**
     * @brief Get the maximum number of iterations
     * @return Maximum number of iterations
     */
    int getMaxIterations() const;

    /**
     * @brief Set the convergence tolerance
     * @param tolerance Convergence tolerance
     */
    void setTolerance(double tolerance);

    /**
     * @brief Get the convergence tolerance
     * @return Convergence tolerance
     */
    double getTolerance() const;

    /**
     * @brief Set the step size (for gradient-based methods)
     * @param stepSize Step size
     */
    void setStepSize(double stepSize);

    /**
     * @brief Get the step size
     * @return Step size
     */
    double getStepSize() const;

    /**
     * @brief Set the population size (for population-based methods)
     * @param populationSize Population size
     */
    void setPopulationSize(int populationSize);

    /**
     * @brief Get the population size
     * @return Population size
     */
    int getPopulationSize() const;

    /**
     * @brief Set the mutation rate (for genetic algorithm)
     * @param mutationRate Mutation rate
     */
    void setMutationRate(double mutationRate);

    /**
     * @brief Get the mutation rate
     * @return Mutation rate
     */
    double getMutationRate() const;

    /**
     * @brief Set the crossover rate (for genetic algorithm)
     * @param crossoverRate Crossover rate
     */
    void setCrossoverRate(double crossoverRate);

    /**
     * @brief Get the crossover rate
     * @return Crossover rate
     */
    double getCrossoverRate() const;

    /**
     * @brief Set the initial temperature (for simulated annealing)
     * @param temperature Initial temperature
     */
    void setInitialTemperature(double temperature);

    /**
     * @brief Get the initial temperature
     * @return Initial temperature
     */
    double getInitialTemperature() const;

    /**
     * @brief Set the cooling rate (for simulated annealing)
     * @param coolingRate Cooling rate
     */
    void setCoolingRate(double coolingRate);

    /**
     * @brief Get the cooling rate
     * @return Cooling rate
     */
    double getCoolingRate() const;

    /**
     * @brief Set the inertia weight (for particle swarm optimization)
     * @param inertiaWeight Inertia weight
     */
    void setInertiaWeight(double inertiaWeight);

    /**
     * @brief Get the inertia weight
     * @return Inertia weight
     */
    double getInertiaWeight() const;

    /**
     * @brief Set the cognitive parameter (for particle swarm optimization)
     * @param cognitiveParameter Cognitive parameter
     */
    void setCognitiveParameter(double cognitiveParameter);

    /**
     * @brief Get the cognitive parameter
     * @return Cognitive parameter
     */
    double getCognitiveParameter() const;

    /**
     * @brief Set the social parameter (for particle swarm optimization)
     * @param socialParameter Social parameter
     */
    void setSocialParameter(double socialParameter);

    /**
     * @brief Get the social parameter
     * @return Social parameter
     */
    double getSocialParameter() const;

    /**
     * @brief Set the scaling factor (for differential evolution)
     * @param scalingFactor Scaling factor
     */
    void setScalingFactor(double scalingFactor);

    /**
     * @brief Get the scaling factor
     * @return Scaling factor
     */
    double getScalingFactor() const;

    /**
     * @brief Set whether to use adaptive parameters
     * @param useAdaptiveParameters Whether to use adaptive parameters
     */
    void setUseAdaptiveParameters(bool useAdaptiveParameters);

    /**
     * @brief Get whether to use adaptive parameters
     * @return Whether to use adaptive parameters
     */
    bool getUseAdaptiveParameters() const;

    /**
     * @brief Set the random seed
     * @param seed Random seed
     */
    void setRandomSeed(unsigned int seed);

    /**
     * @brief Get the random seed
     * @return Random seed
     */
    unsigned int getRandomSeed() const;

    /**
     * @brief Set whether to use parallel execution
     * @param useParallel Whether to use parallel execution
     */
    void setUseParallel(bool useParallel);

    /**
     * @brief Get whether to use parallel execution
     * @return Whether to use parallel execution
     */
    bool getUseParallel() const;

    /**
     * @brief Set the number of threads for parallel execution
     * @param numThreads Number of threads
     */
    void setNumThreads(int numThreads);

    /**
     * @brief Get the number of threads for parallel execution
     * @return Number of threads
     */
    int getNumThreads() const;

private:
    int m_maxIterations = 1000;
    double m_tolerance = 1e-6;
    double m_stepSize = 0.01;
    int m_populationSize = 50;
    double m_mutationRate = 0.1;
    double m_crossoverRate = 0.8;
    double m_initialTemperature = 100.0;
    double m_coolingRate = 0.95;
    double m_inertiaWeight = 0.7;
    double m_cognitiveParameter = 1.5;
    double m_socialParameter = 1.5;
    double m_scalingFactor = 0.5;
    bool m_useAdaptiveParameters = false;
    unsigned int m_randomSeed = 0;
    bool m_useParallel = false;
    int m_numThreads = 4;
};

/**
 * @brief Class for optimization results
 */
class OptimizationResult {
public:
    /**
     * @brief Constructor for optimization results
     * @param solution Optimal solution
     * @param objectiveValue Value of the objective function at the solution
     * @param iterations Number of iterations performed
     * @param functionEvaluations Number of function evaluations
     * @param gradientEvaluations Number of gradient evaluations
     * @param hessianEvaluations Number of Hessian evaluations
     * @param converged Whether the optimization converged
     * @param message Message describing the result
     */
    OptimizationResult(const std::vector<double>& solution, double objectiveValue, int iterations,
                      int functionEvaluations, int gradientEvaluations, int hessianEvaluations,
                      bool converged, const std::string& message);

    /**
     * @brief Get the optimal solution
     * @return Optimal solution
     */
    const std::vector<double>& getSolution() const;

    /**
     * @brief Get the value of the objective function at the solution
     * @return Value of the objective function
     */
    double getObjectiveValue() const;

    /**
     * @brief Get the number of iterations performed
     * @return Number of iterations
     */
    int getIterations() const;

    /**
     * @brief Get the number of function evaluations
     * @return Number of function evaluations
     */
    int getFunctionEvaluations() const;

    /**
     * @brief Get the number of gradient evaluations
     * @return Number of gradient evaluations
     */
    int getGradientEvaluations() const;

    /**
     * @brief Get the number of Hessian evaluations
     * @return Number of Hessian evaluations
     */
    int getHessianEvaluations() const;

    /**
     * @brief Check if the optimization converged
     * @return true if the optimization converged, false otherwise
     */
    bool hasConverged() const;

    /**
     * @brief Get the message describing the result
     * @return Message describing the result
     */
    const std::string& getMessage() const;

private:
    std::vector<double> m_solution;
    double m_objectiveValue;
    int m_iterations;
    int m_functionEvaluations;
    int m_gradientEvaluations;
    int m_hessianEvaluations;
    bool m_converged;
    std::string m_message;
};

/**
 * @brief Class for an optimization solver
 */
class OptimizationSolver {
public:
    /**
     * @brief Constructor for an optimization solver
     * @param problem Optimization problem
     * @param algorithm Optimization algorithm
     * @param options Optimization options
     */
    OptimizationSolver(const OptimizationProblem& problem, OptimizationAlgorithm algorithm,
                      const OptimizationOptions& options = OptimizationOptions());

    /**
     * @brief Set the optimization algorithm
     * @param algorithm Optimization algorithm
     */
    void setAlgorithm(OptimizationAlgorithm algorithm);

    /**
     * @brief Get the optimization algorithm
     * @return Optimization algorithm
     */
    OptimizationAlgorithm getAlgorithm() const;

    /**
     * @brief Set the optimization options
     * @param options Optimization options
     */
    void setOptions(const OptimizationOptions& options);

    /**
     * @brief Get the optimization options
     * @return Optimization options
     */
    const OptimizationOptions& getOptions() const;

    /**
     * @brief Set the initial point
     * @param initialPoint Initial point
     */
    void setInitialPoint(const std::vector<double>& initialPoint);

    /**
     * @brief Get the initial point
     * @return Initial point
     */
    const std::vector<double>& getInitialPoint() const;

    /**
     * @brief Solve the optimization problem
     * @return Optimization result
     */
    OptimizationResult solve();

    /**
     * @brief Get the optimization problem
     * @return Optimization problem
     */
    const OptimizationProblem& getProblem() const;

private:
    const OptimizationProblem& m_problem;
    OptimizationAlgorithm m_algorithm;
    OptimizationOptions m_options;
    std::vector<double> m_initialPoint;
    std::mt19937 m_rng;

    // Helper methods for different algorithms
    OptimizationResult solveGradientDescent();
    OptimizationResult solveNewton();
    OptimizationResult solveBFGS();
    OptimizationResult solveConjugateGradient();
    OptimizationResult solveSimplex();
    OptimizationResult solveSimulatedAnnealing();
    OptimizationResult solveGeneticAlgorithm();
    OptimizationResult solveParticleSwarm();
    OptimizationResult solveDifferentialEvolution();
    OptimizationResult solveBayesianOptimization();

    // Helper methods for numerical approximation
    std::vector<double> approximateGradient(const std::vector<double>& x, double h = 1e-6) const;
    rebelcalc::Matrix approximateHessian(const std::vector<double>& x, double h = 1e-4) const;
};

/**
 * @brief Create a linear programming problem
 * @param objectiveCoefficients Coefficients of the objective function
 * @param constraintMatrix Constraint matrix
 * @param constraintRHS Right-hand side of the constraints
 * @param constraintTypes Types of the constraints
 * @param lowerBounds Lower bounds for the variables
 * @param upperBounds Upper bounds for the variables
 * @return Optimization problem
 */
OptimizationProblem createLinearProgrammingProblem(
    const std::vector<double>& objectiveCoefficients,
    const rebelcalc::Matrix& constraintMatrix,
    const std::vector<double>& constraintRHS,
    const std::vector<ConstraintType>& constraintTypes,
    const std::vector<double>& lowerBounds,
    const std::vector<double>& upperBounds);

/**
 * @brief Create a quadratic programming problem
 * @param objectiveMatrix Matrix in the quadratic term of the objective function
 * @param objectiveVector Vector in the linear term of the objective function
 * @param constraintMatrix Constraint matrix
 * @param constraintRHS Right-hand side of the constraints
 * @param constraintTypes Types of the constraints
 * @param lowerBounds Lower bounds for the variables
 * @param upperBounds Upper bounds for the variables
 * @return Optimization problem
 */
OptimizationProblem createQuadraticProgrammingProblem(
    const rebelcalc::Matrix& objectiveMatrix,
    const std::vector<double>& objectiveVector,
    const rebelcalc::Matrix& constraintMatrix,
    const std::vector<double>& constraintRHS,
    const std::vector<ConstraintType>& constraintTypes,
    const std::vector<double>& lowerBounds,
    const std::vector<double>& upperBounds);

/**
 * @brief Create a nonlinear programming problem
 * @param objective Objective function
 * @param dimension Dimension of the problem
 * @param gradient Gradient of the objective function (optional)
 * @param hessian Hessian of the objective function (optional)
 * @param constraints Constraints (optional)
 * @param lowerBounds Lower bounds for the variables (optional)
 * @param upperBounds Upper bounds for the variables (optional)
 * @return Optimization problem
 */
OptimizationProblem createNonlinearProgrammingProblem(
    std::function<double(const std::vector<double>&)> objective,
    int dimension,
    std::function<std::vector<double>(const std::vector<double>&)> gradient = nullptr,
    std::function<rebelcalc::Matrix(const std::vector<double>&)> hessian = nullptr,
    const std::vector<Constraint>& constraints = {},
    const std::vector<double>& lowerBounds = {},
    const std::vector<double>& upperBounds = {});

} // namespace Solvers
} // namespace RebelCalc

#endif // REBELCALC_SOLVERS_OPTIMIZATION_H
