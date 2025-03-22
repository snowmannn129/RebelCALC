#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <functional>

#include "../src/solvers/optimization.h"

using namespace RebelCalc::Solvers;

// Helper function to print a section header
void printHeader(const std::string& header) {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "  " << header << std::endl;
    std::cout << std::string(80, '=') << std::endl;
}

// Helper function to print a vector
void printVector(const std::string& name, const std::vector<double>& vec) {
    std::cout << name << " = [";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::setw(10) << std::setprecision(6) << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}

// Demo for unconstrained optimization
void demoUnconstrainedOptimization() {
    printHeader("Unconstrained Optimization Demonstration");
    
    // Define the Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2
    auto rosenbrock = [](const std::vector<double>& x) -> double {
        return std::pow(1.0 - x[0], 2) + 100.0 * std::pow(x[1] - x[0] * x[0], 2);
    };
    
    // Define the gradient of the Rosenbrock function
    auto rosenbrockGradient = [](const std::vector<double>& x) -> std::vector<double> {
        std::vector<double> grad(2);
        grad[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
        grad[1] = 200.0 * (x[1] - x[0] * x[0]);
        return grad;
    };
    
    // Create an optimization problem
    OptimizationProblem problem(rosenbrock, 2);
    problem.setGradient(rosenbrockGradient);
    
    // Create optimization options
    OptimizationOptions options;
    options.setMaxIterations(1000);
    options.setTolerance(1e-6);
    options.setStepSize(0.01);
    
    // Create an optimization solver
    OptimizationSolver solver(problem, OptimizationAlgorithm::BFGS, options);
    
    // Set the initial point
    std::vector<double> initialPoint = {-1.2, 1.0};
    solver.setInitialPoint(initialPoint);
    
    std::cout << "Minimizing the Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2" << std::endl;
    std::cout << "Initial point: ";
    printVector("x0", initialPoint);
    
    // Solve the problem
    OptimizationResult result = solver.solve();
    
    std::cout << "\nOptimization result:" << std::endl;
    std::cout << "Converged: " << (result.hasConverged() ? "Yes" : "No") << std::endl;
    std::cout << "Message: " << result.getMessage() << std::endl;
    std::cout << "Iterations: " << result.getIterations() << std::endl;
    std::cout << "Function evaluations: " << result.getFunctionEvaluations() << std::endl;
    std::cout << "Gradient evaluations: " << result.getGradientEvaluations() << std::endl;
    std::cout << "Objective value: " << result.getObjectiveValue() << std::endl;
    printVector("Solution", result.getSolution());
    
    std::cout << "\nExpected solution: [1.0, 1.0]" << std::endl;
}

// Demo for constrained optimization
void demoConstrainedOptimization() {
    printHeader("Constrained Optimization Demonstration");
    
    // Define the objective function: f(x,y) = (x-2)^2 + (y-1)^2
    auto objective = [](const std::vector<double>& x) -> double {
        return std::pow(x[0] - 2.0, 2) + std::pow(x[1] - 1.0, 2);
    };
    
    // Define the gradient of the objective function
    auto gradient = [](const std::vector<double>& x) -> std::vector<double> {
        std::vector<double> grad(2);
        grad[0] = 2.0 * (x[0] - 2.0);
        grad[1] = 2.0 * (x[1] - 1.0);
        return grad;
    };
    
    // Define a constraint: g(x,y) = x + y - 2 <= 0
    auto constraint1 = [](const std::vector<double>& x) -> double {
        return x[0] + x[1] - 2.0;
    };
    
    // Define the gradient of the constraint
    auto constraintGradient1 = [](const std::vector<double>& x) -> std::vector<double> {
        return {1.0, 1.0};
    };
    
    // Define another constraint: h(x,y) = x^2 + y^2 - 1 = 0
    auto constraint2 = [](const std::vector<double>& x) -> double {
        return x[0] * x[0] + x[1] * x[1] - 1.0;
    };
    
    // Define the gradient of the second constraint
    auto constraintGradient2 = [](const std::vector<double>& x) -> std::vector<double> {
        return {2.0 * x[0], 2.0 * x[1]};
    };
    
    // Create an optimization problem
    OptimizationProblem problem(objective, 2);
    problem.setGradient(gradient);
    
    // Add constraints
    Constraint ineqConstraint(ConstraintType::INEQUALITY, constraint1, constraintGradient1);
    Constraint eqConstraint(ConstraintType::EQUALITY, constraint2, constraintGradient2);
    
    problem.addConstraint(ineqConstraint);
    problem.addConstraint(eqConstraint);
    
    // Set bound constraints
    BoundConstraint bounds({-10.0, -10.0}, {10.0, 10.0});
    problem.setBounds(bounds);
    
    // Create optimization options
    OptimizationOptions options;
    options.setMaxIterations(1000);
    options.setTolerance(1e-6);
    
    // Create an optimization solver
    OptimizationSolver solver(problem, OptimizationAlgorithm::SIMPLEX, options);
    
    // Set the initial point
    std::vector<double> initialPoint = {0.5, 0.5};
    solver.setInitialPoint(initialPoint);
    
    std::cout << "Minimizing f(x,y) = (x-2)^2 + (y-1)^2" << std::endl;
    std::cout << "Subject to:" << std::endl;
    std::cout << "  x + y - 2 <= 0" << std::endl;
    std::cout << "  x^2 + y^2 - 1 = 0" << std::endl;
    std::cout << "  -10 <= x,y <= 10" << std::endl;
    std::cout << "Initial point: ";
    printVector("x0", initialPoint);
    
    // Solve the problem
    OptimizationResult result = solver.solve();
    
    std::cout << "\nOptimization result:" << std::endl;
    std::cout << "Converged: " << (result.hasConverged() ? "Yes" : "No") << std::endl;
    std::cout << "Message: " << result.getMessage() << std::endl;
    std::cout << "Iterations: " << result.getIterations() << std::endl;
    std::cout << "Function evaluations: " << result.getFunctionEvaluations() << std::endl;
    std::cout << "Gradient evaluations: " << result.getGradientEvaluations() << std::endl;
    std::cout << "Objective value: " << result.getObjectiveValue() << std::endl;
    printVector("Solution", result.getSolution());
}

// Demo for linear programming
void demoLinearProgramming() {
    printHeader("Linear Programming Demonstration");
    
    // Define the objective coefficients: maximize 3x + 4y
    // For minimization, we negate the coefficients
    std::vector<double> objectiveCoefficients = {-3.0, -4.0};
    
    // Define the constraint matrix
    // 2x + 3y <= 6
    // -x + y <= 1
    // x >= 0, y >= 0
    rebelcalc::Matrix constraintMatrix(2, 2);
    constraintMatrix(0, 0) = 2.0;
    constraintMatrix(0, 1) = 3.0;
    constraintMatrix(1, 0) = -1.0;
    constraintMatrix(1, 1) = 1.0;
    
    // Define the right-hand side of the constraints
    std::vector<double> constraintRHS = {6.0, 1.0};
    
    // Define the constraint types
    std::vector<ConstraintType> constraintTypes = {
        ConstraintType::INEQUALITY,
        ConstraintType::INEQUALITY
    };
    
    // Define the bounds
    std::vector<double> lowerBounds = {0.0, 0.0};
    std::vector<double> upperBounds = {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()
    };
    
    // Create a linear programming problem
    OptimizationProblem problem = createLinearProgrammingProblem(
        objectiveCoefficients,
        constraintMatrix,
        constraintRHS,
        constraintTypes,
        lowerBounds,
        upperBounds
    );
    
    // Create optimization options
    OptimizationOptions options;
    options.setMaxIterations(100);
    options.setTolerance(1e-6);
    
    // Create an optimization solver
    OptimizationSolver solver(problem, OptimizationAlgorithm::SIMPLEX, options);
    
    // Set the initial point
    std::vector<double> initialPoint = {0.0, 0.0};
    solver.setInitialPoint(initialPoint);
    
    std::cout << "Maximizing 3x + 4y" << std::endl;
    std::cout << "Subject to:" << std::endl;
    std::cout << "  2x + 3y <= 6" << std::endl;
    std::cout << "  -x + y <= 1" << std::endl;
    std::cout << "  x >= 0, y >= 0" << std::endl;
    
    // Solve the problem
    OptimizationResult result = solver.solve();
    
    std::cout << "\nOptimization result:" << std::endl;
    std::cout << "Converged: " << (result.hasConverged() ? "Yes" : "No") << std::endl;
    std::cout << "Message: " << result.getMessage() << std::endl;
    std::cout << "Iterations: " << result.getIterations() << std::endl;
    std::cout << "Function evaluations: " << result.getFunctionEvaluations() << std::endl;
    std::cout << "Objective value: " << -result.getObjectiveValue() << " (maximization)" << std::endl;
    printVector("Solution", result.getSolution());
    
    std::cout << "\nExpected solution: [0.0, 2.0]" << std::endl;
    std::cout << "Expected objective value: 8.0" << std::endl;
}

int main() {
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "                      RebelCALC Optimization Demo" << std::endl;
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "\nThis program demonstrates the optimization capabilities of RebelCALC." << std::endl;
    
    demoUnconstrainedOptimization();
    demoConstrainedOptimization();
    demoLinearProgramming();
    
    std::cout << "\n" << std::string(80, '*') << std::endl;
    std::cout << "                      End of Optimization Demo" << std::endl;
    std::cout << std::string(80, '*') << std::endl;
    
    return 0;
}
