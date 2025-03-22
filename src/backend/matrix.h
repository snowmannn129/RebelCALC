#pragma once

#include <vector>
#include <string>
#include <optional>
#include <iostream>
#include <stdexcept>
#include <functional>
#include <initializer_list>

namespace rebelcalc {

/**
 * Class for handling matrix and vector operations
 */
class Matrix {
public:
    /**
     * Default constructor - creates an empty matrix
     */
    Matrix();
    
    /**
     * Constructor with dimensions
     * @param rows Number of rows
     * @param cols Number of columns
     * @param initialValue Initial value for all elements (default: 0.0)
     */
    Matrix(size_t rows, size_t cols, double initialValue = 0.0);
    
    /**
     * Constructor from 2D vector
     * @param data 2D vector of data
     */
    Matrix(const std::vector<std::vector<double>>& data);
    
    /**
     * Constructor from 1D vector (creates a column vector)
     * @param data 1D vector of data
     */
    Matrix(const std::vector<double>& data);
    
    /**
     * Constructor from initializer list of initializer lists
     * @param data Initializer list of initializer lists
     */
    Matrix(std::initializer_list<std::initializer_list<double>> data);
    
    /**
     * Copy constructor
     * @param other Matrix to copy
     */
    Matrix(const Matrix& other);
    
    /**
     * Move constructor
     * @param other Matrix to move
     */
    Matrix(Matrix&& other) noexcept;
    
    /**
     * Destructor
     */
    ~Matrix();
    
    /**
     * Copy assignment operator
     * @param other Matrix to copy
     * @return Reference to this matrix
     */
    Matrix& operator=(const Matrix& other);
    
    /**
     * Move assignment operator
     * @param other Matrix to move
     * @return Reference to this matrix
     */
    Matrix& operator=(Matrix&& other) noexcept;
    
    /**
     * Get the number of rows
     * @return Number of rows
     */
    size_t rows() const;
    
    /**
     * Get the number of columns
     * @return Number of columns
     */
    size_t cols() const;
    
    /**
     * Check if the matrix is empty
     * @return true if the matrix is empty, false otherwise
     */
    bool empty() const;
    
    /**
     * Get the element at the specified position
     * @param row Row index
     * @param col Column index
     * @return Reference to the element
     * @throws std::out_of_range if the indices are out of range
     */
    double& at(size_t row, size_t col);
    
    /**
     * Get the element at the specified position (const version)
     * @param row Row index
     * @param col Column index
     * @return Const reference to the element
     * @throws std::out_of_range if the indices are out of range
     */
    const double& at(size_t row, size_t col) const;
    
    /**
     * Get the element at the specified position
     * @param row Row index
     * @param col Column index
     * @return Reference to the element
     * @note No bounds checking is performed
     */
    double& operator()(size_t row, size_t col);
    
    /**
     * Get the element at the specified position (const version)
     * @param row Row index
     * @param col Column index
     * @return Const reference to the element
     * @note No bounds checking is performed
     */
    const double& operator()(size_t row, size_t col) const;
    
    /**
     * Get a row of the matrix
     * @param row Row index
     * @return Vector containing the row
     * @throws std::out_of_range if the index is out of range
     */
    std::vector<double> getRow(size_t row) const;
    
    /**
     * Get a column of the matrix
     * @param col Column index
     * @return Vector containing the column
     * @throws std::out_of_range if the index is out of range
     */
    std::vector<double> getColumn(size_t col) const;
    
    /**
     * Set a row of the matrix
     * @param row Row index
     * @param data Vector of data to set
     * @throws std::out_of_range if the index is out of range
     * @throws std::invalid_argument if the data size doesn't match the number of columns
     */
    void setRow(size_t row, const std::vector<double>& data);
    
    /**
     * Set a column of the matrix
     * @param col Column index
     * @param data Vector of data to set
     * @throws std::out_of_range if the index is out of range
     * @throws std::invalid_argument if the data size doesn't match the number of rows
     */
    void setColumn(size_t col, const std::vector<double>& data);
    
    /**
     * Resize the matrix
     * @param rows New number of rows
     * @param cols New number of columns
     * @param initialValue Initial value for new elements (default: 0.0)
     */
    void resize(size_t rows, size_t cols, double initialValue = 0.0);
    
    /**
     * Clear the matrix (set to empty)
     */
    void clear();
    
    /**
     * Fill the matrix with a value
     * @param value Value to fill the matrix with
     */
    void fill(double value);
    
    /**
     * Create an identity matrix
     * @param size Size of the identity matrix
     * @return Identity matrix
     */
    static Matrix identity(size_t size);
    
    /**
     * Create a diagonal matrix
     * @param diagonal Vector of diagonal elements
     * @return Diagonal matrix
     */
    static Matrix diagonal(const std::vector<double>& diagonal);
    
    /**
     * Create a matrix filled with random values
     * @param rows Number of rows
     * @param cols Number of columns
     * @param min Minimum value (default: 0.0)
     * @param max Maximum value (default: 1.0)
     * @return Matrix filled with random values
     */
    static Matrix random(size_t rows, size_t cols, double min = 0.0, double max = 1.0);
    
    /**
     * Create a matrix from a function
     * @param rows Number of rows
     * @param cols Number of columns
     * @param func Function that takes row and column indices and returns a value
     * @return Matrix filled with values from the function
     */
    static Matrix fromFunction(size_t rows, size_t cols, 
                              std::function<double(size_t, size_t)> func);
    
    /**
     * Transpose the matrix
     * @return Transposed matrix
     */
    Matrix transpose() const;
    
    /**
     * Calculate the determinant of the matrix
     * @return Determinant of the matrix
     * @throws std::invalid_argument if the matrix is not square
     */
    double determinant() const;
    
    /**
     * Calculate the inverse of the matrix
     * @return Inverse of the matrix
     * @throws std::invalid_argument if the matrix is not square
     * @throws std::runtime_error if the matrix is singular
     */
    Matrix inverse() const;
    
    /**
     * Calculate the trace of the matrix (sum of diagonal elements)
     * @return Trace of the matrix
     * @throws std::invalid_argument if the matrix is not square
     */
    double trace() const;
    
    /**
     * Calculate the rank of the matrix
     * @return Rank of the matrix
     */
    size_t rank() const;
    
    /**
     * Calculate the norm of the matrix
     * @param p Norm type (1, 2, or inf)
     * @return Norm of the matrix
     * @throws std::invalid_argument if p is not 1, 2, or inf
     */
    double norm(int p = 2) const;
    
    /**
     * Calculate the condition number of the matrix
     * @param p Norm type (1, 2, or inf)
     * @return Condition number of the matrix
     * @throws std::invalid_argument if the matrix is not square
     * @throws std::runtime_error if the matrix is singular
     */
    double conditionNumber(int p = 2) const;
    
    /**
     * Calculate the eigenvalues of the matrix
     * @return Vector of eigenvalues
     * @throws std::invalid_argument if the matrix is not square
     * @note This is a placeholder for a real eigenvalue calculation
     */
    std::vector<double> eigenvalues() const;
    
    /**
     * Calculate the eigenvectors of the matrix
     * @return Matrix of eigenvectors (each column is an eigenvector)
     * @throws std::invalid_argument if the matrix is not square
     * @note This is a placeholder for a real eigenvector calculation
     */
    Matrix eigenvectors() const;
    
    /**
     * Solve the linear system Ax = b
     * @param b Right-hand side vector
     * @return Solution vector
     * @throws std::invalid_argument if the matrix is not square
     * @throws std::runtime_error if the matrix is singular
     */
    std::vector<double> solve(const std::vector<double>& b) const;
    
    /**
     * Solve the linear system Ax = B
     * @param B Right-hand side matrix
     * @return Solution matrix
     * @throws std::invalid_argument if the matrix is not square
     * @throws std::runtime_error if the matrix is singular
     */
    Matrix solve(const Matrix& B) const;
    
    /**
     * Calculate the LU decomposition of the matrix
     * @return Pair of matrices (L, U)
     * @throws std::invalid_argument if the matrix is not square
     */
    std::pair<Matrix, Matrix> luDecomposition() const;
    
    /**
     * Calculate the QR decomposition of the matrix
     * @return Pair of matrices (Q, R)
     */
    std::pair<Matrix, Matrix> qrDecomposition() const;
    
    /**
     * Calculate the Cholesky decomposition of the matrix
     * @return Cholesky factor (lower triangular matrix)
     * @throws std::invalid_argument if the matrix is not square
     * @throws std::runtime_error if the matrix is not positive definite
     */
    Matrix choleskyDecomposition() const;
    
    /**
     * Calculate the singular value decomposition of the matrix
     * @return Tuple of matrices (U, S, V)
     * @note This is a placeholder for a real SVD calculation
     */
    std::tuple<Matrix, Matrix, Matrix> singularValueDecomposition() const;
    
    /**
     * Apply a function to each element of the matrix
     * @param func Function to apply
     * @return Matrix with the function applied to each element
     */
    Matrix apply(std::function<double(double)> func) const;
    
    /**
     * Check if the matrix is square
     * @return true if the matrix is square, false otherwise
     */
    bool isSquare() const;
    
    /**
     * Check if the matrix is symmetric
     * @return true if the matrix is symmetric, false otherwise
     */
    bool isSymmetric() const;
    
    /**
     * Check if the matrix is diagonal
     * @return true if the matrix is diagonal, false otherwise
     */
    bool isDiagonal() const;
    
    /**
     * Check if the matrix is upper triangular
     * @return true if the matrix is upper triangular, false otherwise
     */
    bool isUpperTriangular() const;
    
    /**
     * Check if the matrix is lower triangular
     * @return true if the matrix is lower triangular, false otherwise
     */
    bool isLowerTriangular() const;
    
    /**
     * Check if the matrix is orthogonal
     * @return true if the matrix is orthogonal, false otherwise
     */
    bool isOrthogonal() const;
    
    /**
     * Check if the matrix is positive definite
     * @return true if the matrix is positive definite, false otherwise
     */
    bool isPositiveDefinite() const;
    
    /**
     * Convert the matrix to a string
     * @return String representation of the matrix
     */
    std::string toString() const;
    
    /**
     * Addition operator
     * @param other Matrix to add
     * @return Result of addition
     * @throws std::invalid_argument if the matrices have different dimensions
     */
    Matrix operator+(const Matrix& other) const;
    
    /**
     * Subtraction operator
     * @param other Matrix to subtract
     * @return Result of subtraction
     * @throws std::invalid_argument if the matrices have different dimensions
     */
    Matrix operator-(const Matrix& other) const;
    
    /**
     * Multiplication operator (matrix-matrix)
     * @param other Matrix to multiply
     * @return Result of multiplication
     * @throws std::invalid_argument if the matrices have incompatible dimensions
     */
    Matrix operator*(const Matrix& other) const;
    
    /**
     * Multiplication operator (matrix-scalar)
     * @param scalar Scalar to multiply
     * @return Result of multiplication
     */
    Matrix operator*(double scalar) const;
    
    /**
     * Division operator (matrix-scalar)
     * @param scalar Scalar to divide by
     * @return Result of division
     * @throws std::invalid_argument if scalar is zero
     */
    Matrix operator/(double scalar) const;
    
    /**
     * Unary plus operator
     * @return Copy of the matrix
     */
    Matrix operator+() const;
    
    /**
     * Unary minus operator
     * @return Negated matrix
     */
    Matrix operator-() const;
    
    /**
     * Addition assignment operator
     * @param other Matrix to add
     * @return Reference to this matrix
     * @throws std::invalid_argument if the matrices have different dimensions
     */
    Matrix& operator+=(const Matrix& other);
    
    /**
     * Subtraction assignment operator
     * @param other Matrix to subtract
     * @return Reference to this matrix
     * @throws std::invalid_argument if the matrices have different dimensions
     */
    Matrix& operator-=(const Matrix& other);
    
    /**
     * Multiplication assignment operator (matrix-scalar)
     * @param scalar Scalar to multiply
     * @return Reference to this matrix
     */
    Matrix& operator*=(double scalar);
    
    /**
     * Division assignment operator (matrix-scalar)
     * @param scalar Scalar to divide by
     * @return Reference to this matrix
     * @throws std::invalid_argument if scalar is zero
     */
    Matrix& operator/=(double scalar);
    
    /**
     * Equality operator
     * @param other Matrix to compare
     * @return true if the matrices are equal, false otherwise
     */
    bool operator==(const Matrix& other) const;
    
    /**
     * Inequality operator
     * @param other Matrix to compare
     * @return true if the matrices are not equal, false otherwise
     */
    bool operator!=(const Matrix& other) const;

private:
    std::vector<std::vector<double>> m_data; // Matrix data
    
    // Helper methods
    void validateIndices(size_t row, size_t col) const;
    void validateDimensions(const Matrix& other) const;
    void validateMultiplicationDimensions(const Matrix& other) const;
    void validateSquare() const;
};

/**
 * Multiplication operator (scalar-matrix)
 * @param scalar Scalar to multiply
 * @param matrix Matrix to multiply
 * @return Result of multiplication
 */
Matrix operator*(double scalar, const Matrix& matrix);

/**
 * Output stream operator
 * @param os Output stream
 * @param matrix Matrix to output
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const Matrix& matrix);

} // namespace rebelcalc
