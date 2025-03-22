#include "matrix.h"
#include <cmath>
#include <random>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <limits>

namespace rebelcalc {

// Default constructor
Matrix::Matrix() {
    // Create an empty matrix
}

// Constructor with dimensions
Matrix::Matrix(size_t rows, size_t cols, double initialValue) {
    m_data.resize(rows, std::vector<double>(cols, initialValue));
}

// Constructor from 2D vector
Matrix::Matrix(const std::vector<std::vector<double>>& data) {
    if (data.empty()) {
        return;
    }
    
    // Check that all rows have the same number of columns
    size_t cols = data[0].size();
    for (const auto& row : data) {
        if (row.size() != cols) {
            throw std::invalid_argument("All rows must have the same number of columns");
        }
    }
    
    m_data = data;
}

// Constructor from 1D vector (creates a column vector)
Matrix::Matrix(const std::vector<double>& data) {
    m_data.resize(data.size(), std::vector<double>(1));
    for (size_t i = 0; i < data.size(); ++i) {
        m_data[i][0] = data[i];
    }
}

// Constructor from initializer list of initializer lists
Matrix::Matrix(std::initializer_list<std::initializer_list<double>> data) {
    if (data.size() == 0) {
        return;
    }
    
    // Check that all rows have the same number of columns
    size_t cols = data.begin()->size();
    for (const auto& row : data) {
        if (row.size() != cols) {
            throw std::invalid_argument("All rows must have the same number of columns");
        }
    }
    
    m_data.resize(data.size(), std::vector<double>(cols));
    
    size_t i = 0;
    for (const auto& row : data) {
        size_t j = 0;
        for (const auto& val : row) {
            m_data[i][j] = val;
            ++j;
        }
        ++i;
    }
}

// Copy constructor
Matrix::Matrix(const Matrix& other) : m_data(other.m_data) {
}

// Move constructor
Matrix::Matrix(Matrix&& other) noexcept : m_data(std::move(other.m_data)) {
}

// Destructor
Matrix::~Matrix() {
}

// Copy assignment operator
Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        m_data = other.m_data;
    }
    return *this;
}

// Move assignment operator
Matrix& Matrix::operator=(Matrix&& other) noexcept {
    if (this != &other) {
        m_data = std::move(other.m_data);
    }
    return *this;
}

// Get the number of rows
size_t Matrix::rows() const {
    return m_data.size();
}

// Get the number of columns
size_t Matrix::cols() const {
    return m_data.empty() ? 0 : m_data[0].size();
}

// Check if the matrix is empty
bool Matrix::empty() const {
    return m_data.empty();
}

// Get the element at the specified position
double& Matrix::at(size_t row, size_t col) {
    validateIndices(row, col);
    return m_data[row][col];
}

// Get the element at the specified position (const version)
const double& Matrix::at(size_t row, size_t col) const {
    validateIndices(row, col);
    return m_data[row][col];
}

// Get the element at the specified position
double& Matrix::operator()(size_t row, size_t col) {
    return m_data[row][col];
}

// Get the element at the specified position (const version)
const double& Matrix::operator()(size_t row, size_t col) const {
    return m_data[row][col];
}

// Get a row of the matrix
std::vector<double> Matrix::getRow(size_t row) const {
    if (row >= rows()) {
        throw std::out_of_range("Row index out of range");
    }
    return m_data[row];
}

// Get a column of the matrix
std::vector<double> Matrix::getColumn(size_t col) const {
    if (col >= cols()) {
        throw std::out_of_range("Column index out of range");
    }
    
    std::vector<double> column(rows());
    for (size_t i = 0; i < rows(); ++i) {
        column[i] = m_data[i][col];
    }
    
    return column;
}

// Set a row of the matrix
void Matrix::setRow(size_t row, const std::vector<double>& data) {
    if (row >= rows()) {
        throw std::out_of_range("Row index out of range");
    }
    
    if (data.size() != cols()) {
        throw std::invalid_argument("Data size doesn't match the number of columns");
    }
    
    m_data[row] = data;
}

// Set a column of the matrix
void Matrix::setColumn(size_t col, const std::vector<double>& data) {
    if (col >= cols()) {
        throw std::out_of_range("Column index out of range");
    }
    
    if (data.size() != rows()) {
        throw std::invalid_argument("Data size doesn't match the number of rows");
    }
    
    for (size_t i = 0; i < rows(); ++i) {
        m_data[i][col] = data[i];
    }
}

// Resize the matrix
void Matrix::resize(size_t rows, size_t cols, double initialValue) {
    m_data.resize(rows);
    for (auto& row : m_data) {
        row.resize(cols, initialValue);
    }
}

// Clear the matrix (set to empty)
void Matrix::clear() {
    m_data.clear();
}

// Fill the matrix with a value
void Matrix::fill(double value) {
    for (auto& row : m_data) {
        std::fill(row.begin(), row.end(), value);
    }
}

// Create an identity matrix
Matrix Matrix::identity(size_t size) {
    Matrix result(size, size, 0.0);
    for (size_t i = 0; i < size; ++i) {
        result(i, i) = 1.0;
    }
    return result;
}

// Create a diagonal matrix
Matrix Matrix::diagonal(const std::vector<double>& diagonal) {
    size_t size = diagonal.size();
    Matrix result(size, size, 0.0);
    for (size_t i = 0; i < size; ++i) {
        result(i, i) = diagonal[i];
    }
    return result;
}

// Create a matrix filled with random values
Matrix Matrix::random(size_t rows, size_t cols, double min, double max) {
    Matrix result(rows, cols);
    
    // Set up random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(min, max);
    
    // Fill the matrix with random values
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result(i, j) = dist(gen);
        }
    }
    
    return result;
}

// Create a matrix from a function
Matrix Matrix::fromFunction(size_t rows, size_t cols, std::function<double(size_t, size_t)> func) {
    Matrix result(rows, cols);
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result(i, j) = func(i, j);
        }
    }
    
    return result;
}

// Transpose the matrix
Matrix Matrix::transpose() const {
    if (empty()) {
        return Matrix();
    }
    
    Matrix result(cols(), rows());
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            result(j, i) = m_data[i][j];
        }
    }
    
    return result;
}

// Calculate the determinant of the matrix
double Matrix::determinant() const {
    validateSquare();
    
    size_t n = rows();
    
    // Handle special cases
    if (n == 0) {
        return 1.0; // Determinant of empty matrix is 1
    } else if (n == 1) {
        return m_data[0][0];
    } else if (n == 2) {
        return m_data[0][0] * m_data[1][1] - m_data[0][1] * m_data[1][0];
    } else if (n == 3) {
        // For 3x3 matrices, use the rule of Sarrus
        return m_data[0][0] * m_data[1][1] * m_data[2][2] +
               m_data[0][1] * m_data[1][2] * m_data[2][0] +
               m_data[0][2] * m_data[1][0] * m_data[2][1] -
               m_data[0][2] * m_data[1][1] * m_data[2][0] -
               m_data[0][1] * m_data[1][0] * m_data[2][2] -
               m_data[0][0] * m_data[1][2] * m_data[2][1];
    }
    
    // For larger matrices, use LU decomposition
    auto [L, U] = luDecomposition();
    
    // The determinant is the product of the diagonal elements of U
    double det = 1.0;
    for (size_t i = 0; i < n; ++i) {
        det *= U(i, i);
    }
    
    return det;
}

// Calculate the inverse of the matrix
Matrix Matrix::inverse() const {
    validateSquare();
    
    size_t n = rows();
    
    // Handle special cases
    if (n == 0) {
        throw std::invalid_argument("Cannot invert an empty matrix");
    } else if (n == 1) {
        if (std::abs(m_data[0][0]) < 1e-10) {
            throw std::runtime_error("Matrix is singular");
        }
        return Matrix({{1.0 / m_data[0][0]}});
    } else if (n == 2) {
        double det = determinant();
        if (std::abs(det) < 1e-10) {
            throw std::runtime_error("Matrix is singular");
        }
        
        Matrix result(2, 2);
        result(0, 0) = m_data[1][1] / det;
        result(0, 1) = -m_data[0][1] / det;
        result(1, 0) = -m_data[1][0] / det;
        result(1, 1) = m_data[0][0] / det;
        
        return result;
    }
    
    // For larger matrices, use Gauss-Jordan elimination
    Matrix augmented(n, 2 * n);
    
    // Set up the augmented matrix [A|I]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented(i, j) = m_data[i][j];
        }
        augmented(i, i + n) = 1.0;
    }
    
    // Perform Gauss-Jordan elimination
    for (size_t i = 0; i < n; ++i) {
        // Find pivot
        size_t pivotRow = i;
        double pivotValue = std::abs(augmented(i, i));
        
        for (size_t j = i + 1; j < n; ++j) {
            if (std::abs(augmented(j, i)) > pivotValue) {
                pivotRow = j;
                pivotValue = std::abs(augmented(j, i));
            }
        }
        
        // Check for singularity
        if (pivotValue < 1e-10) {
            throw std::runtime_error("Matrix is singular");
        }
        
        // Swap rows if necessary
        if (pivotRow != i) {
            for (size_t j = 0; j < 2 * n; ++j) {
                std::swap(augmented(i, j), augmented(pivotRow, j));
            }
        }
        
        // Scale the pivot row
        double pivot = augmented(i, i);
        for (size_t j = 0; j < 2 * n; ++j) {
            augmented(i, j) /= pivot;
        }
        
        // Eliminate other rows
        for (size_t j = 0; j < n; ++j) {
            if (j != i) {
                double factor = augmented(j, i);
                for (size_t k = 0; k < 2 * n; ++k) {
                    augmented(j, k) -= factor * augmented(i, k);
                }
            }
        }
    }
    
    // Extract the inverse from the augmented matrix
    Matrix result(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result(i, j) = augmented(i, j + n);
        }
    }
    
    return result;
}

// Calculate the trace of the matrix (sum of diagonal elements)
double Matrix::trace() const {
    validateSquare();
    
    double result = 0.0;
    for (size_t i = 0; i < rows(); ++i) {
        result += m_data[i][i];
    }
    
    return result;
}

// Calculate the rank of the matrix
size_t Matrix::rank() const {
    if (empty()) {
        return 0;
    }
    
    // Create a copy of the matrix for row reduction
    Matrix temp(*this);
    
    // Perform Gaussian elimination
    size_t r = 0; // Rank
    size_t m = rows();
    size_t n = cols();
    
    for (size_t i = 0; i < m && i < n; ++i) {
        // Find pivot
        size_t pivotRow = i;
        double pivotValue = std::abs(temp(i, i));
        
        for (size_t j = i + 1; j < m; ++j) {
            if (std::abs(temp(j, i)) > pivotValue) {
                pivotRow = j;
                pivotValue = std::abs(temp(j, i));
            }
        }
        
        // Check if the pivot is zero
        if (pivotValue < 1e-10) {
            continue;
        }
        
        // Swap rows if necessary
        if (pivotRow != i) {
            for (size_t j = 0; j < n; ++j) {
                std::swap(temp(i, j), temp(pivotRow, j));
            }
        }
        
        // Eliminate below
        for (size_t j = i + 1; j < m; ++j) {
            double factor = temp(j, i) / temp(i, i);
            
            for (size_t k = i; k < n; ++k) {
                temp(j, k) -= factor * temp(i, k);
            }
        }
        
        ++r; // Increment rank
    }
    
    return r;
}

// Calculate the norm of the matrix
double Matrix::norm(int p) const {
    if (empty()) {
        return 0.0;
    }
    
    if (p == 1) {
        // 1-norm (maximum column sum)
        double maxSum = 0.0;
        
        for (size_t j = 0; j < cols(); ++j) {
            double sum = 0.0;
            for (size_t i = 0; i < rows(); ++i) {
                sum += std::abs(m_data[i][j]);
            }
            maxSum = std::max(maxSum, sum);
        }
        
        return maxSum;
    } else if (p == 2) {
        // 2-norm (spectral norm)
        // This is a placeholder for a real 2-norm calculation
        // In a real implementation, this would use the singular value decomposition
        
        // For now, just return the Frobenius norm
        double sum = 0.0;
        
        for (size_t i = 0; i < rows(); ++i) {
            for (size_t j = 0; j < cols(); ++j) {
                sum += m_data[i][j] * m_data[i][j];
            }
        }
        
        return std::sqrt(sum);
    } else if (p == std::numeric_limits<int>::max()) {
        // Infinity norm (maximum row sum)
        double maxSum = 0.0;
        
        for (size_t i = 0; i < rows(); ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < cols(); ++j) {
                sum += std::abs(m_data[i][j]);
            }
            maxSum = std::max(maxSum, sum);
        }
        
        return maxSum;
    } else {
        throw std::invalid_argument("Invalid norm type");
    }
}

// Calculate the condition number of the matrix
double Matrix::conditionNumber(int p) const {
    validateSquare();
    
    // Calculate the norm of the matrix
    double normA = norm(p);
    
    // Calculate the norm of the inverse
    double normAInv = inverse().norm(p);
    
    return normA * normAInv;
}

// Calculate the eigenvalues of the matrix
std::vector<double> Matrix::eigenvalues() const {
    validateSquare();
    
    // This is a placeholder for a real eigenvalue calculation
    // In a real implementation, this would use a numerical library
    
    // For now, just return a placeholder result
    std::vector<double> result;
    
    // For 2x2 matrices, we can calculate the eigenvalues analytically
    if (rows() == 2) {
        double a = m_data[0][0];
        double b = m_data[0][1];
        double c = m_data[1][0];
        double d = m_data[1][1];
        
        double trace = a + d;
        double det = a * d - b * c;
        
        double discriminant = trace * trace - 4 * det;
        
        if (discriminant >= 0) {
            double lambda1 = (trace + std::sqrt(discriminant)) / 2;
            double lambda2 = (trace - std::sqrt(discriminant)) / 2;
            
            result.push_back(lambda1);
            result.push_back(lambda2);
        }
    }
    
    return result;
}

// Calculate the eigenvectors of the matrix
Matrix Matrix::eigenvectors() const {
    validateSquare();
    
    // This is a placeholder for a real eigenvector calculation
    // In a real implementation, this would use a numerical library
    
    // For now, just return a placeholder result
    return Matrix();
}

// Solve the linear system Ax = b
std::vector<double> Matrix::solve(const std::vector<double>& b) const {
    validateSquare();
    
    if (b.size() != rows()) {
        throw std::invalid_argument("Right-hand side vector has incompatible dimensions");
    }
    
    // For small matrices, use Cramer's rule
    if (rows() <= 3) {
        double det = determinant();
        if (std::abs(det) < 1e-10) {
            throw std::runtime_error("Matrix is singular");
        }
        
        std::vector<double> x(rows());
        
        for (size_t i = 0; i < rows(); ++i) {
            // Create a copy of the matrix
            Matrix A(*this);
            
            // Replace the i-th column with the right-hand side
            for (size_t j = 0; j < rows(); ++j) {
                A(j, i) = b[j];
            }
            
            // Calculate the determinant of the modified matrix
            x[i] = A.determinant() / det;
        }
        
        return x;
    }
    
    // For larger matrices, use LU decomposition
    auto [L, U] = luDecomposition();
    
    // Solve Ly = b
    std::vector<double> y(rows());
    for (size_t i = 0; i < rows(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < i; ++j) {
            sum += L(i, j) * y[j];
        }
        y[i] = (b[i] - sum) / L(i, i);
    }
    
    // Solve Ux = y
    std::vector<double> x(rows());
    for (int i = static_cast<int>(rows()) - 1; i >= 0; --i) {
        double sum = 0.0;
        for (size_t j = i + 1; j < rows(); ++j) {
            sum += U(i, j) * x[j];
        }
        x[i] = (y[i] - sum) / U(i, i);
    }
    
    return x;
}

// Solve the linear system Ax = B
Matrix Matrix::solve(const Matrix& B) const {
    validateSquare();
    
    if (B.rows() != rows()) {
        throw std::invalid_argument("Right-hand side matrix has incompatible dimensions");
    }
    
    // Solve each column of B separately
    Matrix X(rows(), B.cols());
    
    for (size_t j = 0; j < B.cols(); ++j) {
        std::vector<double> b = B.getColumn(j);
        std::vector<double> x = solve(b);
        
        for (size_t i = 0; i < rows(); ++i) {
            X(i, j) = x[i];
        }
    }
    
    return X;
}

// Calculate the LU decomposition of the matrix
std::pair<Matrix, Matrix> Matrix::luDecomposition() const {
    validateSquare();
    
    size_t n = rows();
    
    Matrix L(n, n, 0.0);
    Matrix U(n, n, 0.0);
    
    // Doolittle algorithm
    for (size_t i = 0; i < n; ++i) {
        // Upper triangular matrix
        for (size_t j = i; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L(i, k) * U(k, j);
            }
            U(i, j) = m_data[i][j] - sum;
        }
        
        // Lower triangular matrix
        L(i, i) = 1.0; // Diagonal elements of L are 1
        for (size_t j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L(j, k) * U(k, i);
            }
            
            if (std::abs(U(i, i)) < 1e-10) {
                throw std::runtime_error("Matrix is singular");
            }
            
            L(j, i) = (m_data[j][i] - sum) / U(i, i);
        }
    }
    
    return {L, U};
}

// Calculate the QR decomposition of the matrix
std::pair<Matrix, Matrix> Matrix::qrDecomposition() const {
    if (empty()) {
        return {Matrix(), Matrix()};
    }
    
    size_t m = rows();
    size_t n = cols();
    
    Matrix Q(m, m, 0.0);
    Matrix R(m, n, 0.0);
    
    // Modified Gram-Schmidt process
    std::vector<std::vector<double>> a(n, std::vector<double>(m));
    
    // Initialize a with the columns of the matrix
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < m; ++i) {
            a[j][i] = m_data[i][j];
        }
    }
    
    for (size_t j = 0; j < n; ++j) {
        double r_jj = 0.0;
        for (size_t i = 0; i < m; ++i) {
            r_jj += a[j][i] * a[j][i];
        }
        r_jj = std::sqrt(r_jj);
        
        R(j, j) = r_jj;
        
        if (std::abs(r_jj) < 1e-10) {
            // Linearly dependent column
            continue;
        }
        
        // Normalize the j-th column of a
        for (size_t i = 0; i < m; ++i) {
            Q(i, j) = a[j][i] / r_jj;
        }
        
        // Orthogonalize the remaining columns
        for (size_t k = j + 1; k < n; ++k) {
            double r_jk = 0.0;
            for (size_t i = 0; i < m; ++i) {
                r_jk += Q(i, j) * a[k][i];
            }
            
            R(j, k) = r_jk;
            
            for (size_t i = 0; i < m; ++i) {
                a[k][i] -= r_jk * Q(i, j);
            }
        }
    }
    
    return {Q, R};
}

// Calculate the Cholesky decomposition of the matrix
Matrix Matrix::choleskyDecomposition() const {
    validateSquare();
    
    if (!isSymmetric()) {
        throw std::invalid_argument("Matrix is not symmetric");
    }
    
    size_t n = rows();
    Matrix L(n, n, 0.0);
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double sum = 0.0;
            
            if (j == i) {
                // Diagonal element
                for (size_t k = 0; k < j; ++k) {
                    sum += L(j, k) * L(j, k);
                }
                
                double value = m_data[j][j] - sum;
                if (value <= 0) {
                    throw std::runtime_error("Matrix is not positive definite");
                }
                
                L(j, j) = std::sqrt(value);
            } else {
                // Off-diagonal element
                for (size_t k = 0; k < j; ++k) {
                    sum += L(i, k) * L(j, k);
                }
                
                L(i, j) = (m_data[i][j] - sum) / L(j, j);
            }
        }
    }
    
    return L;
}

// Calculate the singular value decomposition of the matrix
std::tuple<Matrix, Matrix, Matrix> Matrix::singularValueDecomposition() const {
    // This is a placeholder for a real SVD calculation
    // In a real implementation, this would use a numerical library
    
    // For now, just return placeholder matrices
    return {Matrix(), Matrix(), Matrix()};
}

// Apply a function to each element of the matrix
Matrix Matrix::apply(std::function<double(double)> func) const {
    Matrix result(*this);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            result(i, j) = func(m_data[i][j]);
        }
    }
    
    return result;
}

// Check if the matrix is square
bool Matrix::isSquare() const {
    return rows() == cols();
}

// Check if the matrix is symmetric
bool Matrix::isSymmetric() const {
    if (!isSquare()) {
        return false;
    }
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = i + 1; j < cols(); ++j) {
            if (std::abs(m_data[i][j] - m_data[j][i]) > 1e-10) {
                return false;
            }
        }
    }
    
    return true;
}

// Check if the matrix is diagonal
bool Matrix::isDiagonal() const {
    if (!isSquare()) {
        return false;
    }
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            if (i != j && std::abs(m_data[i][j]) > 1e-10) {
                return false;
            }
        }
    }
    
    return true;
}

// Check if the matrix is upper triangular
bool Matrix::isUpperTriangular() const {
    if (!isSquare()) {
        return false;
    }
    
    for (size_t i = 1; i < rows(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (std::abs(m_data[i][j]) > 1e-10) {
                return false;
            }
        }
    }
    
    return true;
}

// Check if the matrix is lower triangular
bool Matrix::isLowerTriangular() const {
    if (!isSquare()) {
        return false;
    }
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = i + 1; j < cols(); ++j) {
            if (std::abs(m_data[i][j]) > 1e-10) {
                return false;
            }
        }
    }
    
    return true;
}

// Check if the matrix is orthogonal
bool Matrix::isOrthogonal() const {
    if (!isSquare()) {
        return false;
    }
    
    // Check if A^T * A = I
    Matrix product = transpose() * (*this);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(product(i, j) - expected) > 1e-10) {
                return false;
            }
        }
    }
    
    return true;
}

// Check if the matrix is positive definite
bool Matrix::isPositiveDefinite() const {
    if (!isSquare() || !isSymmetric()) {
        return false;
    }
    
    try {
        // A matrix is positive definite if and only if its Cholesky decomposition exists
        choleskyDecomposition();
        return true;
    } catch (const std::runtime_error&) {
        return false;
    }
}

// Convert the matrix to a string
std::string Matrix::toString() const {
    std::ostringstream oss;
    
    if (empty()) {
        oss << "[]";
        return oss.str();
    }
    
    oss << "[";
    
    for (size_t i = 0; i < rows(); ++i) {
        if (i > 0) {
            oss << ",\n ";
        }
        
        oss << "[";
        
        for (size_t j = 0; j < cols(); ++j) {
            if (j > 0) {
                oss << ", ";
            }
            
            oss << std::fixed << std::setprecision(6) << m_data[i][j];
        }
        
        oss << "]";
    }
    
    oss << "]";
    
    return oss.str();
}

// Addition operator
Matrix Matrix::operator+(const Matrix& other) const {
    validateDimensions(other);
    
    Matrix result(*this);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            result(i, j) += other(i, j);
        }
    }
    
    return result;
}

// Subtraction operator
Matrix Matrix::operator-(const Matrix& other) const {
    validateDimensions(other);
    
    Matrix result(*this);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            result(i, j) -= other(i, j);
        }
    }
    
    return result;
}

// Multiplication operator (matrix-matrix)
Matrix Matrix::operator*(const Matrix& other) const {
    validateMultiplicationDimensions(other);
    
    size_t m = rows();
    size_t n = other.cols();
    size_t p = cols();
    
    Matrix result(m, n, 0.0);
    
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < p; ++k) {
                result(i, j) += m_data[i][k] * other(k, j);
            }
        }
    }
    
    return result;
}

// Multiplication operator (matrix-scalar)
Matrix Matrix::operator*(double scalar) const {
    Matrix result(*this);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            result(i, j) *= scalar;
        }
    }
    
    return result;
}

// Division operator (matrix-scalar)
Matrix Matrix::operator/(double scalar) const {
    if (std::abs(scalar) < 1e-10) {
        throw std::invalid_argument("Division by zero");
    }
    
    Matrix result(*this);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            result(i, j) /= scalar;
        }
    }
    
    return result;
}

// Unary plus operator
Matrix Matrix::operator+() const {
    return *this;
}

// Unary minus operator
Matrix Matrix::operator-() const {
    Matrix result(*this);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            result(i, j) = -result(i, j);
        }
    }
    
    return result;
}

// Addition assignment operator
Matrix& Matrix::operator+=(const Matrix& other) {
    validateDimensions(other);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            m_data[i][j] += other(i, j);
        }
    }
    
    return *this;
}

// Subtraction assignment operator
Matrix& Matrix::operator-=(const Matrix& other) {
    validateDimensions(other);
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            m_data[i][j] -= other(i, j);
        }
    }
    
    return *this;
}

// Multiplication assignment operator (matrix-scalar)
Matrix& Matrix::operator*=(double scalar) {
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            m_data[i][j] *= scalar;
        }
    }
    
    return *this;
}

// Division assignment operator (matrix-scalar)
Matrix& Matrix::operator/=(double scalar) {
    if (std::abs(scalar) < 1e-10) {
        throw std::invalid_argument("Division by zero");
    }
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            m_data[i][j] /= scalar;
        }
    }
    
    return *this;
}

// Equality operator
bool Matrix::operator==(const Matrix& other) const {
    if (rows() != other.rows() || cols() != other.cols()) {
        return false;
    }
    
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            if (std::abs(m_data[i][j] - other(i, j)) > 1e-10) {
                return false;
            }
        }
    }
    
    return true;
}

// Inequality operator
bool Matrix::operator!=(const Matrix& other) const {
    return !(*this == other);
}

// Helper methods
void Matrix::validateIndices(size_t row, size_t col) const {
    if (row >= rows() || col >= cols()) {
        throw std::out_of_range("Matrix indices out of range");
    }
}

void Matrix::validateDimensions(const Matrix& other) const {
    if (rows() != other.rows() || cols() != other.cols()) {
        throw std::invalid_argument("Matrix dimensions mismatch");
    }
}

void Matrix::validateMultiplicationDimensions(const Matrix& other) const {
    if (cols() != other.rows()) {
        throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
    }
}

void Matrix::validateSquare() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix must be square");
    }
}

// Non-member functions

// Multiplication operator (scalar-matrix)
Matrix operator*(double scalar, const Matrix& matrix) {
    return matrix * scalar;
}

// Output stream operator
std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    os << matrix.toString();
    return os;
}

} // namespace rebelcalc
