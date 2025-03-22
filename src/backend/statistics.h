#pragma once

#include <vector>
#include <optional>
#include <string>
#include <unordered_map>
#include <functional>

namespace rebelcalc {

/**
 * Statistics class for statistical calculations
 */
class Statistics {
public:
    /**
     * Constructor
     */
    Statistics();
    
    /**
     * Destructor
     */
    ~Statistics();
    
    /**
     * Initialize the statistics engine
     * @return true if initialization was successful, false otherwise
     */
    bool initialize();
    
    /**
     * Shutdown the statistics engine
     */
    void shutdown();
    
    /**
     * Calculate the mean (average) of a dataset
     * @param data The dataset
     * @return The mean value
     */
    double mean(const std::vector<double>& data) const;
    
    /**
     * Calculate the median of a dataset
     * @param data The dataset
     * @return The median value
     */
    double median(const std::vector<double>& data) const;
    
    /**
     * Calculate the mode of a dataset
     * @param data The dataset
     * @return The mode value(s)
     */
    std::vector<double> mode(const std::vector<double>& data) const;
    
    /**
     * Calculate the range of a dataset
     * @param data The dataset
     * @return The range (max - min)
     */
    double range(const std::vector<double>& data) const;
    
    /**
     * Calculate the variance of a dataset
     * @param data The dataset
     * @param sample Whether to calculate sample variance (n-1) or population variance (n)
     * @return The variance
     */
    double variance(const std::vector<double>& data, bool sample = true) const;
    
    /**
     * Calculate the standard deviation of a dataset
     * @param data The dataset
     * @param sample Whether to calculate sample standard deviation (n-1) or population standard deviation (n)
     * @return The standard deviation
     */
    double standardDeviation(const std::vector<double>& data, bool sample = true) const;
    
    /**
     * Calculate the covariance of two datasets
     * @param dataX The first dataset
     * @param dataY The second dataset
     * @param sample Whether to calculate sample covariance (n-1) or population covariance (n)
     * @return The covariance, or nullopt if the datasets have different sizes
     */
    std::optional<double> covariance(const std::vector<double>& dataX, 
                                    const std::vector<double>& dataY, 
                                    bool sample = true) const;
    
    /**
     * Calculate the correlation coefficient of two datasets
     * @param dataX The first dataset
     * @param dataY The second dataset
     * @return The correlation coefficient, or nullopt if the datasets have different sizes
     */
    std::optional<double> correlation(const std::vector<double>& dataX, 
                                     const std::vector<double>& dataY) const;
    
    /**
     * Perform linear regression on two datasets
     * @param dataX The independent variable dataset
     * @param dataY The dependent variable dataset
     * @return A pair containing the slope and y-intercept, or nullopt if the datasets have different sizes
     */
    std::optional<std::pair<double, double>> linearRegression(const std::vector<double>& dataX, 
                                                             const std::vector<double>& dataY) const;
    
    /**
     * Calculate the percentile of a dataset
     * @param data The dataset
     * @param percentile The percentile to calculate (0-100)
     * @return The percentile value
     */
    double percentile(const std::vector<double>& data, double percentile) const;
    
    /**
     * Calculate the quartiles of a dataset
     * @param data The dataset
     * @return A vector containing the first quartile (25%), median (50%), and third quartile (75%)
     */
    std::vector<double> quartiles(const std::vector<double>& data) const;
    
    /**
     * Calculate the interquartile range (IQR) of a dataset
     * @param data The dataset
     * @return The interquartile range
     */
    double interquartileRange(const std::vector<double>& data) const;
    
    /**
     * Calculate the z-score of a value in a dataset
     * @param data The dataset
     * @param value The value to calculate the z-score for
     * @return The z-score
     */
    double zScore(const std::vector<double>& data, double value) const;
    
    /**
     * Calculate the skewness of a dataset
     * @param data The dataset
     * @return The skewness
     */
    double skewness(const std::vector<double>& data) const;
    
    /**
     * Calculate the kurtosis of a dataset
     * @param data The dataset
     * @return The kurtosis
     */
    double kurtosis(const std::vector<double>& data) const;
    
    /**
     * Calculate the geometric mean of a dataset
     * @param data The dataset
     * @return The geometric mean, or nullopt if any value is negative
     */
    std::optional<double> geometricMean(const std::vector<double>& data) const;
    
    /**
     * Calculate the harmonic mean of a dataset
     * @param data The dataset
     * @return The harmonic mean, or nullopt if any value is zero or negative
     */
    std::optional<double> harmonicMean(const std::vector<double>& data) const;
    
    /**
     * Calculate the root mean square (RMS) of a dataset
     * @param data The dataset
     * @return The root mean square
     */
    double rootMeanSquare(const std::vector<double>& data) const;
    
    /**
     * Calculate the mean absolute deviation of a dataset
     * @param data The dataset
     * @return The mean absolute deviation
     */
    double meanAbsoluteDeviation(const std::vector<double>& data) const;
    
    /**
     * Calculate the median absolute deviation of a dataset
     * @param data The dataset
     * @return The median absolute deviation
     */
    double medianAbsoluteDeviation(const std::vector<double>& data) const;
    
    /**
     * Calculate the factorial of a non-negative integer
     * @param n The non-negative integer
     * @return The factorial, or nullopt if n is negative
     */
    std::optional<double> factorial(int n) const;
    
    /**
     * Calculate the binomial coefficient (n choose k)
     * @param n The total number of items
     * @param k The number of items to choose
     * @return The binomial coefficient, or nullopt if n or k is negative or k > n
     */
    std::optional<double> binomialCoefficient(int n, int k) const;
    
    /**
     * Calculate the probability mass function (PMF) of a binomial distribution
     * @param n The number of trials
     * @param p The probability of success
     * @param k The number of successes
     * @return The probability mass function value, or nullopt if parameters are invalid
     */
    std::optional<double> binomialPMF(int n, double p, int k) const;
    
    /**
     * Calculate the cumulative distribution function (CDF) of a binomial distribution
     * @param n The number of trials
     * @param p The probability of success
     * @param k The number of successes
     * @return The cumulative distribution function value, or nullopt if parameters are invalid
     */
    std::optional<double> binomialCDF(int n, double p, int k) const;
    
    /**
     * Calculate the probability density function (PDF) of a normal distribution
     * @param x The value
     * @param mean The mean of the distribution
     * @param stdDev The standard deviation of the distribution
     * @return The probability density function value, or nullopt if stdDev is non-positive
     */
    std::optional<double> normalPDF(double x, double mean, double stdDev) const;
    
    /**
     * Calculate the cumulative distribution function (CDF) of a normal distribution
     * @param x The value
     * @param mean The mean of the distribution
     * @param stdDev The standard deviation of the distribution
     * @return The cumulative distribution function value, or nullopt if stdDev is non-positive
     */
    std::optional<double> normalCDF(double x, double mean, double stdDev) const;
    
    /**
     * Calculate the inverse cumulative distribution function (quantile function) of a normal distribution
     * @param p The probability (0-1)
     * @param mean The mean of the distribution
     * @param stdDev The standard deviation of the distribution
     * @return The inverse cumulative distribution function value, or nullopt if parameters are invalid
     */
    std::optional<double> normalInvCDF(double p, double mean, double stdDev) const;
    
    /**
     * Calculate the probability density function (PDF) of a t-distribution
     * @param x The value
     * @param df The degrees of freedom
     * @return The probability density function value, or nullopt if df is non-positive
     */
    std::optional<double> tPDF(double x, double df) const;
    
    /**
     * Calculate the cumulative distribution function (CDF) of a t-distribution
     * @param x The value
     * @param df The degrees of freedom
     * @return The cumulative distribution function value, or nullopt if df is non-positive
     */
    std::optional<double> tCDF(double x, double df) const;
    
    /**
     * Calculate the inverse cumulative distribution function (quantile function) of a t-distribution
     * @param p The probability (0-1)
     * @param df The degrees of freedom
     * @return The inverse cumulative distribution function value, or nullopt if parameters are invalid
     */
    std::optional<double> tInvCDF(double p, double df) const;
    
    /**
     * Calculate the probability density function (PDF) of a chi-squared distribution
     * @param x The value
     * @param df The degrees of freedom
     * @return The probability density function value, or nullopt if parameters are invalid
     */
    std::optional<double> chiSquaredPDF(double x, double df) const;
    
    /**
     * Calculate the cumulative distribution function (CDF) of a chi-squared distribution
     * @param x The value
     * @param df The degrees of freedom
     * @return The cumulative distribution function value, or nullopt if parameters are invalid
     */
    std::optional<double> chiSquaredCDF(double x, double df) const;
    
    /**
     * Calculate the inverse cumulative distribution function (quantile function) of a chi-squared distribution
     * @param p The probability (0-1)
     * @param df The degrees of freedom
     * @return The inverse cumulative distribution function value, or nullopt if parameters are invalid
     */
    std::optional<double> chiSquaredInvCDF(double p, double df) const;
    
    /**
     * Calculate the probability density function (PDF) of an F-distribution
     * @param x The value
     * @param df1 The numerator degrees of freedom
     * @param df2 The denominator degrees of freedom
     * @return The probability density function value, or nullopt if parameters are invalid
     */
    std::optional<double> fPDF(double x, double df1, double df2) const;
    
    /**
     * Calculate the cumulative distribution function (CDF) of an F-distribution
     * @param x The value
     * @param df1 The numerator degrees of freedom
     * @param df2 The denominator degrees of freedom
     * @return The cumulative distribution function value, or nullopt if parameters are invalid
     */
    std::optional<double> fCDF(double x, double df1, double df2) const;
    
    /**
     * Calculate the inverse cumulative distribution function (quantile function) of an F-distribution
     * @param p The probability (0-1)
     * @param df1 The numerator degrees of freedom
     * @param df2 The denominator degrees of freedom
     * @return The inverse cumulative distribution function value, or nullopt if parameters are invalid
     */
    std::optional<double> fInvCDF(double p, double df1, double df2) const;
    
    /**
     * Perform a one-sample t-test
     * @param data The dataset
     * @param mu The population mean to test against
     * @return A tuple containing the t-statistic, p-value, and degrees of freedom
     */
    std::tuple<double, double, double> tTest(const std::vector<double>& data, double mu) const;
    
    /**
     * Perform a two-sample t-test
     * @param data1 The first dataset
     * @param data2 The second dataset
     * @param equalVariance Whether to assume equal variance
     * @return A tuple containing the t-statistic, p-value, and degrees of freedom, or nullopt if datasets are empty
     */
    std::optional<std::tuple<double, double, double>> tTest(const std::vector<double>& data1, 
                                                           const std::vector<double>& data2, 
                                                           bool equalVariance = false) const;
    
    /**
     * Perform a chi-squared test for independence
     * @param observed The observed frequencies (matrix as flattened vector, row-major order)
     * @param rows The number of rows in the matrix
     * @param cols The number of columns in the matrix
     * @return A pair containing the chi-squared statistic and p-value, or nullopt if parameters are invalid
     */
    std::optional<std::pair<double, double>> chiSquaredTest(const std::vector<double>& observed, 
                                                           size_t rows, 
                                                           size_t cols) const;
    
    /**
     * Perform a one-way ANOVA test
     * @param groups A vector of datasets, one for each group
     * @return A tuple containing the F-statistic, p-value, between-group degrees of freedom, and within-group degrees of freedom,
     *         or nullopt if parameters are invalid
     */
    std::optional<std::tuple<double, double, double, double>> anovaTest(const std::vector<std::vector<double>>& groups) const;
    
    /**
     * Get all available statistical functions
     * @return A map of function names to descriptions
     */
    std::unordered_map<std::string, std::string> getAllFunctions() const;

private:
    // Helper methods
    double calculateMean(const std::vector<double>& data) const;
    double calculateVariance(const std::vector<double>& data, bool sample) const;
    std::vector<double> sortData(const std::vector<double>& data) const;
    
    // Error function (erf) and related functions for normal distribution calculations
    double erf(double x) const;
    double erfc(double x) const;
    double erfInv(double x) const;
};

} // namespace rebelcalc
