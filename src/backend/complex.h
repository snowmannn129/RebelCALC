#pragma once

#include <iostream>
#include <cmath>
#include <string>

namespace rebelcalc {

/**
 * Class representing a complex number
 */
class Complex {
public:
    /**
     * Default constructor (creates a complex number with real and imaginary parts set to 0)
     */
    Complex() : m_real(0.0), m_imag(0.0) {}
    
    /**
     * Constructor with real and imaginary parts
     * @param real The real part
     * @param imag The imaginary part
     */
    Complex(double real, double imag = 0.0) : m_real(real), m_imag(imag) {}
    
    /**
     * Copy constructor
     * @param other The complex number to copy
     */
    Complex(const Complex& other) = default;
    
    /**
     * Assignment operator
     * @param other The complex number to assign from
     * @return Reference to this complex number
     */
    Complex& operator=(const Complex& other) = default;
    
    /**
     * Get the real part
     * @return The real part
     */
    double real() const { return m_real; }
    
    /**
     * Get the imaginary part
     * @return The imaginary part
     */
    double imag() const { return m_imag; }
    
    /**
     * Set the real part
     * @param real The new real part
     */
    void setReal(double real) { m_real = real; }
    
    /**
     * Set the imaginary part
     * @param imag The new imaginary part
     */
    void setImag(double imag) { m_imag = imag; }
    
    /**
     * Get the magnitude (absolute value) of the complex number
     * @return The magnitude
     */
    double magnitude() const {
        return std::sqrt(m_real * m_real + m_imag * m_imag);
    }
    
    /**
     * Get the argument (phase) of the complex number in radians
     * @return The argument
     */
    double argument() const {
        return std::atan2(m_imag, m_real);
    }
    
    /**
     * Get the complex conjugate
     * @return The complex conjugate
     */
    Complex conjugate() const {
        return Complex(m_real, -m_imag);
    }
    
    /**
     * Convert to string representation
     * @return String representation of the complex number
     */
    std::string toString() const {
        if (m_imag == 0.0) {
            return std::to_string(m_real);
        } else if (m_real == 0.0) {
            return std::to_string(m_imag) + "i";
        } else if (m_imag > 0.0) {
            return std::to_string(m_real) + " + " + std::to_string(m_imag) + "i";
        } else {
            return std::to_string(m_real) + " - " + std::to_string(-m_imag) + "i";
        }
    }
    
    // Arithmetic operators
    
    /**
     * Addition operator
     * @param other The complex number to add
     * @return The sum
     */
    Complex operator+(const Complex& other) const {
        return Complex(m_real + other.m_real, m_imag + other.m_imag);
    }
    
    /**
     * Subtraction operator
     * @param other The complex number to subtract
     * @return The difference
     */
    Complex operator-(const Complex& other) const {
        return Complex(m_real - other.m_real, m_imag - other.m_imag);
    }
    
    /**
     * Multiplication operator
     * @param other The complex number to multiply by
     * @return The product
     */
    Complex operator*(const Complex& other) const {
        return Complex(
            m_real * other.m_real - m_imag * other.m_imag,
            m_real * other.m_imag + m_imag * other.m_real
        );
    }
    
    /**
     * Division operator
     * @param other The complex number to divide by
     * @return The quotient
     */
    Complex operator/(const Complex& other) const {
        double denominator = other.m_real * other.m_real + other.m_imag * other.m_imag;
        if (denominator == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        
        return Complex(
            (m_real * other.m_real + m_imag * other.m_imag) / denominator,
            (m_imag * other.m_real - m_real * other.m_imag) / denominator
        );
    }
    
    /**
     * Unary minus operator
     * @return The negation of this complex number
     */
    Complex operator-() const {
        return Complex(-m_real, -m_imag);
    }
    
    // Compound assignment operators
    
    /**
     * Addition assignment operator
     * @param other The complex number to add
     * @return Reference to this complex number
     */
    Complex& operator+=(const Complex& other) {
        m_real += other.m_real;
        m_imag += other.m_imag;
        return *this;
    }
    
    /**
     * Subtraction assignment operator
     * @param other The complex number to subtract
     * @return Reference to this complex number
     */
    Complex& operator-=(const Complex& other) {
        m_real -= other.m_real;
        m_imag -= other.m_imag;
        return *this;
    }
    
    /**
     * Multiplication assignment operator
     * @param other The complex number to multiply by
     * @return Reference to this complex number
     */
    Complex& operator*=(const Complex& other) {
        double real = m_real * other.m_real - m_imag * other.m_imag;
        double imag = m_real * other.m_imag + m_imag * other.m_real;
        m_real = real;
        m_imag = imag;
        return *this;
    }
    
    /**
     * Division assignment operator
     * @param other The complex number to divide by
     * @return Reference to this complex number
     */
    Complex& operator/=(const Complex& other) {
        double denominator = other.m_real * other.m_real + other.m_imag * other.m_imag;
        if (denominator == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        
        double real = (m_real * other.m_real + m_imag * other.m_imag) / denominator;
        double imag = (m_imag * other.m_real - m_real * other.m_imag) / denominator;
        m_real = real;
        m_imag = imag;
        return *this;
    }
    
    // Comparison operators
    
    /**
     * Equality operator
     * @param other The complex number to compare with
     * @return True if the complex numbers are equal, false otherwise
     */
    bool operator==(const Complex& other) const {
        return m_real == other.m_real && m_imag == other.m_imag;
    }
    
    /**
     * Inequality operator
     * @param other The complex number to compare with
     * @return True if the complex numbers are not equal, false otherwise
     */
    bool operator!=(const Complex& other) const {
        return !(*this == other);
    }
    
private:
    double m_real; // Real part
    double m_imag; // Imaginary part
};

// Non-member arithmetic operators for mixed operations with double

/**
 * Addition operator for double + Complex
 * @param lhs The double value
 * @param rhs The complex number
 * @return The sum
 */
inline Complex operator+(double lhs, const Complex& rhs) {
    return Complex(lhs + rhs.real(), rhs.imag());
}

/**
 * Subtraction operator for double - Complex
 * @param lhs The double value
 * @param rhs The complex number
 * @return The difference
 */
inline Complex operator-(double lhs, const Complex& rhs) {
    return Complex(lhs - rhs.real(), -rhs.imag());
}

/**
 * Multiplication operator for double * Complex
 * @param lhs The double value
 * @param rhs The complex number
 * @return The product
 */
inline Complex operator*(double lhs, const Complex& rhs) {
    return Complex(lhs * rhs.real(), lhs * rhs.imag());
}

/**
 * Division operator for double / Complex
 * @param lhs The double value
 * @param rhs The complex number
 * @return The quotient
 */
inline Complex operator/(double lhs, const Complex& rhs) {
    double denominator = rhs.real() * rhs.real() + rhs.imag() * rhs.imag();
    if (denominator == 0.0) {
        throw std::runtime_error("Division by zero");
    }
    
    return Complex(
        lhs * rhs.real() / denominator,
        -lhs * rhs.imag() / denominator
    );
}

// Complex math functions

/**
 * Calculate the square root of a complex number
 * @param z The complex number
 * @return The square root
 */
inline Complex sqrt(const Complex& z) {
    double r = z.magnitude();
    double theta = z.argument() / 2.0;
    return Complex(std::sqrt(r) * std::cos(theta), std::sqrt(r) * std::sin(theta));
}

/**
 * Calculate the exponential of a complex number
 * @param z The complex number
 * @return The exponential
 */
inline Complex exp(const Complex& z) {
    double expReal = std::exp(z.real());
    return Complex(expReal * std::cos(z.imag()), expReal * std::sin(z.imag()));
}

/**
 * Calculate the natural logarithm of a complex number
 * @param z The complex number
 * @return The natural logarithm
 */
inline Complex log(const Complex& z) {
    return Complex(std::log(z.magnitude()), z.argument());
}

/**
 * Calculate the sine of a complex number
 * @param z The complex number
 * @return The sine
 */
inline Complex sin(const Complex& z) {
    return Complex(
        std::sin(z.real()) * std::cosh(z.imag()),
        std::cos(z.real()) * std::sinh(z.imag())
    );
}

/**
 * Calculate the cosine of a complex number
 * @param z The complex number
 * @return The cosine
 */
inline Complex cos(const Complex& z) {
    return Complex(
        std::cos(z.real()) * std::cosh(z.imag()),
        -std::sin(z.real()) * std::sinh(z.imag())
    );
}

/**
 * Calculate the tangent of a complex number
 * @param z The complex number
 * @return The tangent
 */
inline Complex tan(const Complex& z) {
    return sin(z) / cos(z);
}

/**
 * Calculate the power of a complex number
 * @param base The base
 * @param exponent The exponent
 * @return The power
 */
inline Complex pow(const Complex& base, const Complex& exponent) {
    // z^w = exp(w * log(z))
    return exp(exponent * log(base));
}

/**
 * Calculate the power of a complex number with a double exponent
 * @param base The base
 * @param exponent The exponent
 * @return The power
 */
inline Complex pow(const Complex& base, double exponent) {
    return pow(base, Complex(exponent));
}

/**
 * Calculate the power of a double base with a complex exponent
 * @param base The base
 * @param exponent The exponent
 * @return The power
 */
inline Complex pow(double base, const Complex& exponent) {
    return pow(Complex(base), exponent);
}

} // namespace rebelcalc
