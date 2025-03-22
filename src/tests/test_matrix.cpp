#include "../backend/matrix.h"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

namespace rebelcalc {
namespace test {

// Test fixture for Matrix tests
class MatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up test matrices
        m1 = Matrix(2, 2, 1.0);
        m2 = Matrix({{1.0, 2.0}, {3.0, 4.0}});
        m3 = Matrix({{5.0, 6.0}, {7.0, 8.0}});
        
        // Set up test vectors
        v1 = std::vector<double>{1.0, 2.0};
        v2 = std::vector<double>{3.0, 4.0};
    }
    
    Matrix m1; // 2x2 matrix filled with 1.0
    Matrix m2; // 2x2 matrix with values [[1, 2], [3, 4]]
    Matrix m3; // 2x2 matrix with values [[5, 6], [7, 8]]
    std::vector<double> v1; // Vector [1, 2]
    std::vector<double> v2; // Vector [3, 4]
};

// Test constructors
TEST_F(MatrixTest, Constructors) {
    // Default constructor
    Matrix m;
    EXPECT_TRUE(m.empty());
    EXPECT_EQ(m.rows(), 0);
    EXPECT_EQ(m.cols(), 0);
    
    // Constructor with dimensions
    EXPECT_EQ(m1.rows(), 2);
    EXPECT_EQ(m1.cols(), 2);
    EXPECT_EQ(m1(0, 0), 1.0);
    EXPECT_EQ(m1(0, 1), 1.0);
    EXPECT_EQ(m1(1, 0), 1.0);
    EXPECT_EQ(m1(1, 1), 1.0);
    
    // Constructor from 2D vector
    std::vector<std::vector<double>> data = {{1.0, 2.0}, {3.0, 4.0}};
    Matrix m4(data);
    EXPECT_EQ(m4.rows(), 2);
    EXPECT_EQ(m4.cols(), 2);
    EXPECT_EQ(m4(0, 0), 1.0);
    EXPECT_EQ(m4(0, 1), 2.0);
    EXPECT_EQ(m4(1, 0), 3.0);
    EXPECT_EQ(m4(1, 1), 4.0);
    
    // Constructor from 1D vector (creates a column vector)
    Matrix m5(v1);
    EXPECT_EQ(m5.rows(), 2);
    EXPECT_EQ(m5.cols(), 1);
    EXPECT_EQ(m5(0, 0), 1.0);
    EXPECT_EQ(m5(1, 0), 2.0);
    
    // Constructor from initializer list
    Matrix m6({{1.0, 2.0}, {3.0, 4.0}});
    EXPECT_EQ(m6.rows(), 2);
    EXPECT_EQ(m6.cols(), 2);
    EXPECT_EQ(m6(0, 0), 1.0);
    EXPECT_EQ(m6(0, 1), 2.0);
    EXPECT_EQ(m6(1, 0), 3.0);
    EXPECT_EQ(m6(1, 1), 4.0);
    
    // Copy constructor
    Matrix m7(m2);
    EXPECT_EQ(m7.rows(), 2);
    EXPECT_EQ(m7.cols(), 2);
    EXPECT_EQ(m7(0, 0), 1.0);
    EXPECT_EQ(m7(0, 1), 2.0);
    EXPECT_EQ(m7(1, 0), 3.0);
    EXPECT_EQ(m7(1, 1), 4.0);
    
    // Move constructor
    Matrix m8(std::move(m7));
    EXPECT_EQ(m8.rows(), 2);
    EXPECT_EQ(m8.cols(), 2);
    EXPECT_EQ(m8(0, 0), 1.0);
    EXPECT_EQ(m8(0, 1), 2.0);
    EXPECT_EQ(m8(1, 0), 3.0);
    EXPECT_EQ(m8(1, 1), 4.0);
    EXPECT_TRUE(m7.empty()); // m7 should be moved from
}

// Test assignment operators
TEST_F(MatrixTest, AssignmentOperators) {
    // Copy assignment
    Matrix m4;
    m4 = m2;
    EXPECT_EQ(m4.rows(), 2);
    EXPECT_EQ(m4.cols(), 2);
    EXPECT_EQ(m4(0, 0), 1.0);
    EXPECT_EQ(m4(0, 1), 2.0);
    EXPECT_EQ(m4(1, 0), 3.0);
    EXPECT_EQ(m4(1, 1), 4.0);
    
    // Move assignment
    Matrix m5;
    m5 = std::move(m4);
    EXPECT_EQ(m5.rows(), 2);
    EXPECT_EQ(m5.cols(), 2);
    EXPECT_EQ(m5(0, 0), 1.0);
    EXPECT_EQ(m5(0, 1), 2.0);
    EXPECT_EQ(m5(1, 0), 3.0);
    EXPECT_EQ(m5(1, 1), 4.0);
    EXPECT_TRUE(m4.empty()); // m4 should be moved from
}

// Test accessors
TEST_F(MatrixTest, Accessors) {
    // at() with bounds checking
    EXPECT_EQ(m2.at(0, 0), 1.0);
    EXPECT_EQ(m2.at(0, 1), 2.0);
    EXPECT_EQ(m2.at(1, 0), 3.0);
    EXPECT_EQ(m2.at(1, 1), 4.0);
    EXPECT_THROW(m2.at(2, 0), std::out_of_range);
    EXPECT_THROW(m2.at(0, 2), std::out_of_range);
    
    // operator() without bounds checking
    EXPECT_EQ(m2(0, 0), 1.0);
    EXPECT_EQ(m2(0, 1), 2.0);
    EXPECT_EQ(m2(1, 0), 3.0);
    EXPECT_EQ(m2(1, 1), 4.0);
    
    // getRow()
    std::vector<double> row0 = m2.getRow(0);
    EXPECT_EQ(row0.size(), 2);
    EXPECT_EQ(row0[0], 1.0);
    EXPECT_EQ(row0[1], 2.0);
    
    std::vector<double> row1 = m2.getRow(1);
    EXPECT_EQ(row1.size(), 2);
    EXPECT_EQ(row1[0], 3.0);
    EXPECT_EQ(row1[1], 4.0);
    
    EXPECT_THROW(m2.getRow(2), std::out_of_range);
    
    // getColumn()
    std::vector<double> col0 = m2.getColumn(0);
    EXPECT_EQ(col0.size(), 2);
    EXPECT_EQ(col0[0], 1.0);
    EXPECT_EQ(col0[1], 3.0);
    
    std::vector<double> col1 = m2.getColumn(1);
    EXPECT_EQ(col1.size(), 2);
    EXPECT_EQ(col1[0], 2.0);
    EXPECT_EQ(col1[1], 4.0);
    
    EXPECT_THROW(m2.getColumn(2), std::out_of_range);
}

// Test modifiers
TEST_F(MatrixTest, Modifiers) {
    // setRow()
    Matrix m4(m2);
    m4.setRow(0, std::vector<double>{5.0, 6.0});
    EXPECT_EQ(m4(0, 0), 5.0);
    EXPECT_EQ(m4(0, 1), 6.0);
    EXPECT_EQ(m4(1, 0), 3.0);
    EXPECT_EQ(m4(1, 1), 4.0);
    
    EXPECT_THROW(m4.setRow(2, std::vector<double>{7.0, 8.0}), std::out_of_range);
    EXPECT_THROW(m4.setRow(0, std::vector<double>{7.0}), std::invalid_argument);
    
    // setColumn()
    Matrix m5(m2);
    m5.setColumn(0, std::vector<double>{5.0, 6.0});
    EXPECT_EQ(m5(0, 0), 5.0);
    EXPECT_EQ(m5(0, 1), 2.0);
    EXPECT_EQ(m5(1, 0), 6.0);
    EXPECT_EQ(m5(1, 1), 4.0);
    
    EXPECT_THROW(m5.setColumn(2, std::vector<double>{7.0, 8.0}), std::out_of_range);
    EXPECT_THROW(m5.setColumn(0, std::vector<double>{7.0}), std::invalid_argument);
    
    // resize()
    Matrix m6(m2);
    m6.resize(3, 3, 0.0);
    EXPECT_EQ(m6.rows(), 3);
    EXPECT_EQ(m6.cols(), 3);
    EXPECT_EQ(m6(0, 0), 1.0);
    EXPECT_EQ(m6(0, 1), 2.0);
    EXPECT_EQ(m6(0, 2), 0.0);
    EXPECT_EQ(m6(1, 0), 3.0);
    EXPECT_EQ(m6(1, 1), 4.0);
    EXPECT_EQ(m6(1, 2), 0.0);
    EXPECT_EQ(m6(2, 0), 0.0);
    EXPECT_EQ(m6(2, 1), 0.0);
    EXPECT_EQ(m6(2, 2), 0.0);
    
    // clear()
    Matrix m7(m2);
    m7.clear();
    EXPECT_TRUE(m7.empty());
    EXPECT_EQ(m7.rows(), 0);
    EXPECT_EQ(m7.cols(), 0);
    
    // fill()
    Matrix m8(2, 2);
    m8.fill(3.0);
    EXPECT_EQ(m8(0, 0), 3.0);
    EXPECT_EQ(m8(0, 1), 3.0);
    EXPECT_EQ(m8(1, 0), 3.0);
    EXPECT_EQ(m8(1, 1), 3.0);
}

// Test static factory methods
TEST_F(MatrixTest, StaticFactoryMethods) {
    // identity()
    Matrix m4 = Matrix::identity(3);
    EXPECT_EQ(m4.rows(), 3);
    EXPECT_EQ(m4.cols(), 3);
    EXPECT_EQ(m4(0, 0), 1.0);
    EXPECT_EQ(m4(0, 1), 0.0);
    EXPECT_EQ(m4(0, 2), 0.0);
    EXPECT_EQ(m4(1, 0), 0.0);
    EXPECT_EQ(m4(1, 1), 1.0);
    EXPECT_EQ(m4(1, 2), 0.0);
    EXPECT_EQ(m4(2, 0), 0.0);
    EXPECT_EQ(m4(2, 1), 0.0);
    EXPECT_EQ(m4(2, 2), 1.0);
    
    // diagonal()
    Matrix m5 = Matrix::diagonal(std::vector<double>{1.0, 2.0, 3.0});
    EXPECT_EQ(m5.rows(), 3);
    EXPECT_EQ(m5.cols(), 3);
    EXPECT_EQ(m5(0, 0), 1.0);
    EXPECT_EQ(m5(0, 1), 0.0);
    EXPECT_EQ(m5(0, 2), 0.0);
    EXPECT_EQ(m5(1, 0), 0.0);
    EXPECT_EQ(m5(1, 1), 2.0);
    EXPECT_EQ(m5(1, 2), 0.0);
    EXPECT_EQ(m5(2, 0), 0.0);
    EXPECT_EQ(m5(2, 1), 0.0);
    EXPECT_EQ(m5(2, 2), 3.0);
    
    // random()
    Matrix m6 = Matrix::random(2, 2, 0.0, 1.0);
    EXPECT_EQ(m6.rows(), 2);
    EXPECT_EQ(m6.cols(), 2);
    EXPECT_GE(m6(0, 0), 0.0);
    EXPECT_LE(m6(0, 0), 1.0);
    EXPECT_GE(m6(0, 1), 0.0);
    EXPECT_LE(m6(0, 1), 1.0);
    EXPECT_GE(m6(1, 0), 0.0);
    EXPECT_LE(m6(1, 0), 1.0);
    EXPECT_GE(m6(1, 1), 0.0);
    EXPECT_LE(m6(1, 1), 1.0);
    
    // fromFunction()
    Matrix m7 = Matrix::fromFunction(2, 2, [](size_t i, size_t j) { return static_cast<double>(i + j); });
    EXPECT_EQ(m7.rows(), 2);
    EXPECT_EQ(m7.cols(), 2);
    EXPECT_EQ(m7(0, 0), 0.0);
    EXPECT_EQ(m7(0, 1), 1.0);
    EXPECT_EQ(m7(1, 0), 1.0);
    EXPECT_EQ(m7(1, 1), 2.0);
}

// Test matrix operations
TEST_F(MatrixTest, MatrixOperations) {
    // transpose()
    Matrix m4 = m2.transpose();
    EXPECT_EQ(m4.rows(), 2);
    EXPECT_EQ(m4.cols(), 2);
    EXPECT_EQ(m4(0, 0), 1.0);
    EXPECT_EQ(m4(0, 1), 3.0);
    EXPECT_EQ(m4(1, 0), 2.0);
    EXPECT_EQ(m4(1, 1), 4.0);
    
    // determinant()
    EXPECT_EQ(m2.determinant(), -2.0);
    
    // inverse()
    Matrix m5 = m2.inverse();
    EXPECT_EQ(m5.rows(), 2);
    EXPECT_EQ(m5.cols(), 2);
    EXPECT_NEAR(m5(0, 0), -2.0, 1e-10);
    EXPECT_NEAR(m5(0, 1), 1.0, 1e-10);
    EXPECT_NEAR(m5(1, 0), 1.5, 1e-10);
    EXPECT_NEAR(m5(1, 1), -0.5, 1e-10);
    
    // trace()
    EXPECT_EQ(m2.trace(), 5.0);
    
    // rank()
    EXPECT_EQ(m2.rank(), 2);
    
    // norm()
    EXPECT_NEAR(m2.norm(1), 6.0, 1e-10); // 1-norm (maximum column sum)
    EXPECT_NEAR(m2.norm(2), 5.4772, 1e-4); // 2-norm (Frobenius norm)
    EXPECT_NEAR(m2.norm(std::numeric_limits<int>::max()), 7.0, 1e-10); // Infinity norm (maximum row sum)
    
    // apply()
    Matrix m6 = m2.apply([](double x) { return x * 2.0; });
    EXPECT_EQ(m6.rows(), 2);
    EXPECT_EQ(m6.cols(), 2);
    EXPECT_EQ(m6(0, 0), 2.0);
    EXPECT_EQ(m6(0, 1), 4.0);
    EXPECT_EQ(m6(1, 0), 6.0);
    EXPECT_EQ(m6(1, 1), 8.0);
}

// Test matrix properties
TEST_F(MatrixTest, MatrixProperties) {
    // isSquare()
    EXPECT_TRUE(m2.isSquare());
    
    Matrix m4(2, 3);
    EXPECT_FALSE(m4.isSquare());
    
    // isSymmetric()
    Matrix m5({{1.0, 2.0}, {2.0, 3.0}});
    EXPECT_TRUE(m5.isSymmetric());
    EXPECT_FALSE(m2.isSymmetric());
    
    // isDiagonal()
    Matrix m6({{1.0, 0.0}, {0.0, 2.0}});
    EXPECT_TRUE(m6.isDiagonal());
    EXPECT_FALSE(m2.isDiagonal());
    
    // isUpperTriangular()
    Matrix m7({{1.0, 2.0}, {0.0, 3.0}});
    EXPECT_TRUE(m7.isUpperTriangular());
    EXPECT_FALSE(m2.isUpperTriangular());
    
    // isLowerTriangular()
    Matrix m8({{1.0, 0.0}, {2.0, 3.0}});
    EXPECT_TRUE(m8.isLowerTriangular());
    EXPECT_FALSE(m2.isLowerTriangular());
    
    // isOrthogonal()
    Matrix m9({{0.0, 1.0}, {-1.0, 0.0}});
    EXPECT_TRUE(m9.isOrthogonal());
    EXPECT_FALSE(m2.isOrthogonal());
    
    // isPositiveDefinite()
    Matrix m10({{2.0, 1.0}, {1.0, 2.0}});
    EXPECT_TRUE(m10.isPositiveDefinite());
    EXPECT_FALSE(m2.isPositiveDefinite());
}

// Test matrix decompositions
TEST_F(MatrixTest, MatrixDecompositions) {
    // LU decomposition
    auto [L, U] = m2.luDecomposition();
    
    EXPECT_EQ(L.rows(), 2);
    EXPECT_EQ(L.cols(), 2);
    EXPECT_NEAR(L(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(L(0, 1), 0.0, 1e-10);
    EXPECT_NEAR(L(1, 0), 3.0, 1e-10);
    EXPECT_NEAR(L(1, 1), 1.0, 1e-10);
    
    EXPECT_EQ(U.rows(), 2);
    EXPECT_EQ(U.cols(), 2);
    EXPECT_NEAR(U(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(U(0, 1), 2.0, 1e-10);
    EXPECT_NEAR(U(1, 0), 0.0, 1e-10);
    EXPECT_NEAR(U(1, 1), -2.0, 1e-10);
    
    // Verify L * U = A
    Matrix LU = L * U;
    EXPECT_NEAR(LU(0, 0), m2(0, 0), 1e-10);
    EXPECT_NEAR(LU(0, 1), m2(0, 1), 1e-10);
    EXPECT_NEAR(LU(1, 0), m2(1, 0), 1e-10);
    EXPECT_NEAR(LU(1, 1), m2(1, 1), 1e-10);
    
    // Cholesky decomposition
    Matrix m4({{4.0, 2.0}, {2.0, 5.0}});
    Matrix L2 = m4.choleskyDecomposition();
    
    EXPECT_EQ(L2.rows(), 2);
    EXPECT_EQ(L2.cols(), 2);
    EXPECT_NEAR(L2(0, 0), 2.0, 1e-10);
    EXPECT_NEAR(L2(0, 1), 0.0, 1e-10);
    EXPECT_NEAR(L2(1, 0), 1.0, 1e-10);
    EXPECT_NEAR(L2(1, 1), 2.0, 1e-10);
    
    // Verify L2 * L2^T = A
    Matrix L2T = L2.transpose();
    Matrix L2L2T = L2 * L2T;
    EXPECT_NEAR(L2L2T(0, 0), m4(0, 0), 1e-10);
    EXPECT_NEAR(L2L2T(0, 1), m4(0, 1), 1e-10);
    EXPECT_NEAR(L2L2T(1, 0), m4(1, 0), 1e-10);
    EXPECT_NEAR(L2L2T(1, 1), m4(1, 1), 1e-10);
}

// Test linear system solving
TEST_F(MatrixTest, LinearSystemSolving) {
    // solve(vector)
    std::vector<double> b = {5.0, 11.0};
    std::vector<double> x = m2.solve(b);
    
    EXPECT_EQ(x.size(), 2);
    EXPECT_NEAR(x[0], 1.0, 1e-10);
    EXPECT_NEAR(x[1], 2.0, 1e-10);
    
    // Verify A * x = b
    EXPECT_NEAR(m2(0, 0) * x[0] + m2(0, 1) * x[1], b[0], 1e-10);
    EXPECT_NEAR(m2(1, 0) * x[0] + m2(1, 1) * x[1], b[1], 1e-10);
    
    // solve(matrix)
    Matrix B({{5.0, 7.0}, {11.0, 15.0}});
    Matrix X = m2.solve(B);
    
    EXPECT_EQ(X.rows(), 2);
    EXPECT_EQ(X.cols(), 2);
    EXPECT_NEAR(X(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(X(0, 1), 1.0, 1e-10);
    EXPECT_NEAR(X(1, 0), 2.0, 1e-10);
    EXPECT_NEAR(X(1, 1), 3.0, 1e-10);
    
    // Verify A * X = B
    Matrix AX = m2 * X;
    EXPECT_NEAR(AX(0, 0), B(0, 0), 1e-10);
    EXPECT_NEAR(AX(0, 1), B(0, 1), 1e-10);
    EXPECT_NEAR(AX(1, 0), B(1, 0), 1e-10);
    EXPECT_NEAR(AX(1, 1), B(1, 1), 1e-10);
}

// Test arithmetic operators
TEST_F(MatrixTest, ArithmeticOperators) {
    // Addition
    Matrix m4 = m2 + m3;
    EXPECT_EQ(m4.rows(), 2);
    EXPECT_EQ(m4.cols(), 2);
    EXPECT_EQ(m4(0, 0), 6.0);
    EXPECT_EQ(m4(0, 1), 8.0);
    EXPECT_EQ(m4(1, 0), 10.0);
    EXPECT_EQ(m4(1, 1), 12.0);
    
    // Subtraction
    Matrix m5 = m3 - m2;
    EXPECT_EQ(m5.rows(), 2);
    EXPECT_EQ(m5.cols(), 2);
    EXPECT_EQ(m5(0, 0), 4.0);
    EXPECT_EQ(m5(0, 1), 4.0);
    EXPECT_EQ(m5(1, 0), 4.0);
    EXPECT_EQ(m5(1, 1), 4.0);
    
    // Multiplication (matrix-matrix)
    Matrix m6 = m2 * m3;
    EXPECT_EQ(m6.rows(), 2);
    EXPECT_EQ(m6.cols(), 2);
    EXPECT_EQ(m6(0, 0), 19.0);
    EXPECT_EQ(m6(0, 1), 22.0);
    EXPECT_EQ(m6(1, 0), 43.0);
    EXPECT_EQ(m6(1, 1), 50.0);
    
    // Multiplication (matrix-scalar)
    Matrix m7 = m2 * 2.0;
    EXPECT_EQ(m7.rows(), 2);
    EXPECT_EQ(m7.cols(), 2);
    EXPECT_EQ(m7(0, 0), 2.0);
    EXPECT_EQ(m7(0, 1), 4.0);
    EXPECT_EQ(m7(1, 0), 6.0);
    EXPECT_EQ(m7(1, 1), 8.0);
    
    // Multiplication (scalar-matrix)
    Matrix m8 = 2.0 * m2;
    EXPECT_EQ(m8.rows(), 2);
    EXPECT_EQ(m8.cols(), 2);
    EXPECT_EQ(m8(0, 0), 2.0);
    EXPECT_EQ(m8(0, 1), 4.0);
    EXPECT_EQ(m8(1, 0), 6.0);
    EXPECT_EQ(m8(1, 1), 8.0);
    
    // Division (matrix-scalar)
    Matrix m9 = m2 / 2.0;
    EXPECT_EQ(m9.rows(), 2);
    EXPECT_EQ(m9.cols(), 2);
    EXPECT_EQ(m9(0, 0), 0.5);
    EXPECT_EQ(m9(0, 1), 1.0);
    EXPECT_EQ(m9(1, 0), 1.5);
    EXPECT_EQ(m9(1, 1), 2.0);
    
    // Unary plus
    Matrix m10 = +m2;
    EXPECT_EQ(m10.rows(), 2);
    EXPECT_EQ(m10.cols(), 2);
    EXPECT_EQ(m10(0, 0), 1.0);
    EXPECT_EQ(m10(0, 1), 2.0);
    EXPECT_EQ(m10(1, 0), 3.0);
    EXPECT_EQ(m10(1, 1), 4.0);
    
    // Unary minus
    Matrix m11 = -m2;
    EXPECT_EQ(m11.rows(), 2);
    EXPECT_EQ(m11.cols(), 2);
    EXPECT_EQ(m11(0, 0), -1.0);
    EXPECT_EQ(m11(0, 1), -2.0);
    EXPECT_EQ(m11(1, 0), -3.0);
    EXPECT_EQ(m11(1, 1), -4.0);
}

// Test compound assignment operators
TEST_F(MatrixTest, CompoundAssignmentOperators) {
    // Addition assignment
    Matrix m4(m2);
    m4 += m3;
    EXPECT_EQ(m4.rows(), 2);
    EXPECT_EQ(m4.cols(), 2);
    EXPECT_EQ(m4(0, 0), 6.0);
    EXPECT_EQ(m4(0, 1), 8.0);
    EXPECT_EQ(m4(1, 0), 10.0);
    EXPECT_EQ(m4(1, 1), 12.0);
    
    // Subtraction assignment
    Matrix m5(m3);
    m5 -= m2;
    EXPECT_EQ(m5.rows(), 2);
    EXPECT_EQ(m5.cols(), 2);
    EXPECT_EQ(m5(0, 0), 4.0);
    EXPECT_EQ(m5(0, 1), 4.0);
    EXPECT_EQ(m5(1, 0), 4.0);
    EXPECT_EQ(m5(1, 1), 4.0);
    
    // Multiplication assignment (matrix-scalar)
    Matrix m6(m2);
    m6 *= 2.0;
    EXPECT_EQ(m6.rows(), 2);
    EXPECT_EQ(m6.cols(), 2);
    EXPECT_EQ(m6(0, 0), 2.0);
    EXPECT_EQ(m6(0, 1), 4.0);
    EXPECT_EQ(m6(1, 0), 6.0);
    EXPECT_EQ(m6(1, 1), 8.0);
    
    // Division assignment (matrix-scalar)
    Matrix m7(m2);
    m7 /= 2.0;
    EXPECT_EQ(m7.rows(), 2);
    EXPECT_EQ(m7.cols(), 2);
    EXPECT_EQ(m7(0, 0), 0.5);
    EXPECT_EQ(m7(0, 1), 1.0);
    EXPECT_EQ(m7(1, 0), 1.5);
    EXPECT_EQ(m7(1, 1), 2.0);
}

// Test comparison operators
TEST_F(MatrixTest, ComparisonOperators) {
    // Equality
    Matrix m4(m2);
    EXPECT_TRUE(m2 == m4);
    EXPECT_FALSE(m2 == m3);
    
    // Inequality
    EXPECT_FALSE(m2 != m4);
    EXPECT_TRUE(m2 != m3);
}

// Test toString
TEST_F(MatrixTest, ToString) {
    Matrix m4({{1.0, 2.0}, {3.0, 4.0}});
    std::string str = m4.toString();
    EXPECT_FALSE(str.empty());
    // The exact format may vary, so we just check that it contains the values
    EXPECT_TRUE(str.find("1") != std::string::npos);
    EXPECT_TRUE(str.find("2") != std::string::npos);
    EXPECT_TRUE(str.find("3") != std::string::npos);
    EXPECT_TRUE(str.find("4") != std::string::npos);
}

} // namespace test
} // namespace rebelcalc
