#include <gtest/gtest.h>
#include "../backend/complex.h"

namespace rebelcalc {
namespace testing {

class ComplexTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up code here
    }

    void TearDown() override {
        // Tear down code here
    }
};

// Test constructors
TEST_F(ComplexTest, Constructors) {
    // Default constructor
    Complex c1;
    EXPECT_DOUBLE_EQ(0.0, c1.real());
    EXPECT_DOUBLE_EQ(0.0, c1.imag());
    
    // Constructor with real part only
    Complex c2(3.5);
    EXPECT_DOUBLE_EQ(3.5, c2.real());
    EXPECT_DOUBLE_EQ(0.0, c2.imag());
    
    // Constructor with real and imaginary parts
    Complex c3(2.5, 4.5);
    EXPECT_DOUBLE_EQ(2.5, c3.real());
    EXPECT_DOUBLE_EQ(4.5, c3.imag());
    
    // Copy constructor
    Complex c4(c3);
    EXPECT_DOUBLE_EQ(2.5, c4.real());
    EXPECT_DOUBLE_EQ(4.5, c4.imag());
}

// Test getters and setters
TEST_F(ComplexTest, GettersAndSetters) {
    Complex c(1.0, 2.0);
    
    // Test getters
    EXPECT_DOUBLE_EQ(1.0, c.real());
    EXPECT_DOUBLE_EQ(2.0, c.imag());
    
    // Test setters
    c.setReal(3.0);
    c.setImag(4.0);
    EXPECT_DOUBLE_EQ(3.0, c.real());
    EXPECT_DOUBLE_EQ(4.0, c.imag());
}

// Test magnitude and argument
TEST_F(ComplexTest, MagnitudeAndArgument) {
    Complex c(3.0, 4.0);
    
    // Magnitude should be 5.0 (Pythagorean theorem)
    EXPECT_DOUBLE_EQ(5.0, c.magnitude());
    
    // Argument should be atan2(4.0, 3.0)
    EXPECT_DOUBLE_EQ(std::atan2(4.0, 3.0), c.argument());
    
    // Test with negative values
    Complex c2(-3.0, -4.0);
    EXPECT_DOUBLE_EQ(5.0, c2.magnitude());
    EXPECT_DOUBLE_EQ(std::atan2(-4.0, -3.0), c2.argument());
}

// Test conjugate
TEST_F(ComplexTest, Conjugate) {
    Complex c(3.0, 4.0);
    Complex conj = c.conjugate();
    
    EXPECT_DOUBLE_EQ(3.0, conj.real());
    EXPECT_DOUBLE_EQ(-4.0, conj.imag());
}

// Test toString
TEST_F(ComplexTest, ToString) {
    // Real part only
    Complex c1(3.0, 0.0);
    EXPECT_EQ("3.000000", c1.toString());
    
    // Imaginary part only
    Complex c2(0.0, 4.0);
    EXPECT_EQ("4.000000i", c2.toString());
    
    // Both parts, positive imaginary
    Complex c3(3.0, 4.0);
    EXPECT_EQ("3.000000 + 4.000000i", c3.toString());
    
    // Both parts, negative imaginary
    Complex c4(3.0, -4.0);
    EXPECT_EQ("3.000000 - 4.000000i", c4.toString());
}

// Test arithmetic operators
TEST_F(ComplexTest, ArithmeticOperators) {
    Complex c1(3.0, 4.0);
    Complex c2(1.0, 2.0);
    
    // Addition
    Complex sum = c1 + c2;
    EXPECT_DOUBLE_EQ(4.0, sum.real());
    EXPECT_DOUBLE_EQ(6.0, sum.imag());
    
    // Subtraction
    Complex diff = c1 - c2;
    EXPECT_DOUBLE_EQ(2.0, diff.real());
    EXPECT_DOUBLE_EQ(2.0, diff.imag());
    
    // Multiplication
    Complex prod = c1 * c2;
    EXPECT_DOUBLE_EQ(3.0 * 1.0 - 4.0 * 2.0, prod.real());
    EXPECT_DOUBLE_EQ(3.0 * 2.0 + 4.0 * 1.0, prod.imag());
    
    // Division
    Complex quot = c1 / c2;
    double denom = 1.0 * 1.0 + 2.0 * 2.0;
    EXPECT_DOUBLE_EQ((3.0 * 1.0 + 4.0 * 2.0) / denom, quot.real());
    EXPECT_DOUBLE_EQ((4.0 * 1.0 - 3.0 * 2.0) / denom, quot.imag());
    
    // Unary minus
    Complex neg = -c1;
    EXPECT_DOUBLE_EQ(-3.0, neg.real());
    EXPECT_DOUBLE_EQ(-4.0, neg.imag());
}

// Test compound assignment operators
TEST_F(ComplexTest, CompoundAssignmentOperators) {
    Complex c1(3.0, 4.0);
    Complex c2(1.0, 2.0);
    
    // Addition assignment
    Complex c3 = c1;
    c3 += c2;
    EXPECT_DOUBLE_EQ(4.0, c3.real());
    EXPECT_DOUBLE_EQ(6.0, c3.imag());
    
    // Subtraction assignment
    Complex c4 = c1;
    c4 -= c2;
    EXPECT_DOUBLE_EQ(2.0, c4.real());
    EXPECT_DOUBLE_EQ(2.0, c4.imag());
    
    // Multiplication assignment
    Complex c5 = c1;
    c5 *= c2;
    EXPECT_DOUBLE_EQ(3.0 * 1.0 - 4.0 * 2.0, c5.real());
    EXPECT_DOUBLE_EQ(3.0 * 2.0 + 4.0 * 1.0, c5.imag());
    
    // Division assignment
    Complex c6 = c1;
    c6 /= c2;
    double denom = 1.0 * 1.0 + 2.0 * 2.0;
    EXPECT_DOUBLE_EQ((3.0 * 1.0 + 4.0 * 2.0) / denom, c6.real());
    EXPECT_DOUBLE_EQ((4.0 * 1.0 - 3.0 * 2.0) / denom, c6.imag());
}

// Test comparison operators
TEST_F(ComplexTest, ComparisonOperators) {
    Complex c1(3.0, 4.0);
    Complex c2(3.0, 4.0);
    Complex c3(1.0, 2.0);
    
    // Equality
    EXPECT_TRUE(c1 == c2);
    EXPECT_FALSE(c1 == c3);
    
    // Inequality
    EXPECT_FALSE(c1 != c2);
    EXPECT_TRUE(c1 != c3);
}

// Test non-member arithmetic operators
TEST_F(ComplexTest, NonMemberArithmeticOperators) {
    Complex c(3.0, 4.0);
    double d = 2.0;
    
    // Addition
    Complex sum1 = d + c;
    EXPECT_DOUBLE_EQ(5.0, sum1.real());
    EXPECT_DOUBLE_EQ(4.0, sum1.imag());
    
    // Subtraction
    Complex diff1 = d - c;
    EXPECT_DOUBLE_EQ(-1.0, diff1.real());
    EXPECT_DOUBLE_EQ(-4.0, diff1.imag());
    
    // Multiplication
    Complex prod1 = d * c;
    EXPECT_DOUBLE_EQ(6.0, prod1.real());
    EXPECT_DOUBLE_EQ(8.0, prod1.imag());
    
    // Division
    Complex quot1 = d / c;
    double denom = 3.0 * 3.0 + 4.0 * 4.0;
    EXPECT_DOUBLE_EQ(2.0 * 3.0 / denom, quot1.real());
    EXPECT_DOUBLE_EQ(-2.0 * 4.0 / denom, quot1.imag());
}

// Test complex math functions
TEST_F(ComplexTest, MathFunctions) {
    // Square root
    Complex c1(3.0, 4.0);
    Complex sqrt_c1 = sqrt(c1);
    EXPECT_NEAR(2.0, sqrt_c1.real(), 1e-10);
    EXPECT_NEAR(1.0, sqrt_c1.imag(), 1e-10);
    
    // Exponential
    Complex c2(0.0, M_PI);
    Complex exp_c2 = exp(c2);
    EXPECT_NEAR(-1.0, exp_c2.real(), 1e-10);
    EXPECT_NEAR(0.0, exp_c2.imag(), 1e-10);
    
    // Natural logarithm
    Complex c3(1.0, 0.0);
    Complex log_c3 = log(c3);
    EXPECT_NEAR(0.0, log_c3.real(), 1e-10);
    EXPECT_NEAR(0.0, log_c3.imag(), 1e-10);
    
    // Sine
    Complex c4(0.0, 0.0);
    Complex sin_c4 = sin(c4);
    EXPECT_NEAR(0.0, sin_c4.real(), 1e-10);
    EXPECT_NEAR(0.0, sin_c4.imag(), 1e-10);
    
    // Cosine
    Complex cos_c4 = cos(c4);
    EXPECT_NEAR(1.0, cos_c4.real(), 1e-10);
    EXPECT_NEAR(0.0, cos_c4.imag(), 1e-10);
    
    // Tangent
    Complex tan_c4 = tan(c4);
    EXPECT_NEAR(0.0, tan_c4.real(), 1e-10);
    EXPECT_NEAR(0.0, tan_c4.imag(), 1e-10);
    
    // Power
    Complex c5(2.0, 0.0);
    Complex c6(3.0, 0.0);
    Complex pow_c5_c6 = pow(c5, c6);
    EXPECT_NEAR(8.0, pow_c5_c6.real(), 1e-10);
    EXPECT_NEAR(0.0, pow_c5_c6.imag(), 1e-10);
    
    // Power with double exponent
    Complex pow_c5_d = pow(c5, 3.0);
    EXPECT_NEAR(8.0, pow_c5_d.real(), 1e-10);
    EXPECT_NEAR(0.0, pow_c5_d.imag(), 1e-10);
    
    // Power with double base
    Complex pow_d_c6 = pow(2.0, c6);
    EXPECT_NEAR(8.0, pow_d_c6.real(), 1e-10);
    EXPECT_NEAR(0.0, pow_d_c6.imag(), 1e-10);
}

} // namespace testing
} // namespace rebelcalc
