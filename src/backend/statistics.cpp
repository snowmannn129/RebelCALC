#include "statistics.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <map>
#include <limits>
#include <iostream>

// Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace rebelcalc {

// Helper method to calculate the mean of a dataset
double Statistics::calculateMean(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate mean of empty dataset");
    }
    
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}

// Helper method to calculate the variance of a dataset
double Statistics::calculateVariance(const std::vector<double>& data, bool sample) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate variance of empty dataset");
    }
    
    if (data.size() == 1 && sample) {
        throw std::invalid_argument("Cannot calculate sample variance with only one data point");
    }
    
    double meanValue = calculateMean(data);
    
    double sum = 0.0;
    for (double value : data) {
        double diff = value - meanValue;
        sum += diff * diff;
    }
    
    double divisor = sample ? data.size() - 1 : data.size();
    return sum / divisor;
}

// Helper method to sort a dataset
std::vector<double> Statistics::sortData(const std::vector<double>& data) const {
    std::vector<double> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());
    return sortedData;
}

// Error function implementation
double Statistics::erf(double x) const {
    // Constants for the approximation
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;
    
    // Save the sign of x
    int sign = (x < 0) ? -1 : 1;
    x = std::abs(x);
    
    // Abramowitz and Stegun approximation
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);
    
    return sign * y;
}

// Complementary error function implementation
double Statistics::erfc(double x) const {
    return 1.0 - erf(x);
}

// Inverse error function implementation
double Statistics::erfInv(double x) const {
    if (x < -1.0 || x > 1.0) {
        throw std::invalid_argument("Argument to inverse error function must be between -1 and 1");
    }
    
    if (x == 0.0) {
        return 0.0;
    }
    
    // Approximation for |x| <= 0.7
    if (std::abs(x) <= 0.7) {
        const double c1 = 1.0;
        const double c2 = 1.0 / 3.0;
        const double c3 = 7.0 / 30.0;
        const double c4 = 127.0 / 630.0;
        const double c5 = 4369.0 / 22680.0;
        const double c6 = 34807.0 / 178200.0;
        
        double z = x * x;
        return x * (c1 + z * (c2 + z * (c3 + z * (c4 + z * (c5 + z * c6)))));
    }
    
    // Approximation for |x| > 0.7
    double z;
    if (x > 0.0) {
        z = std::sqrt(-std::log((1.0 - x) / 2.0));
    } else {
        z = -std::sqrt(-std::log((1.0 + x) / 2.0));
    }
    
    // Rational approximation
    const double p0 = 1.570796288;
    const double p1 = 0.03706987906;
    const double p2 = -0.8364353589e-3;
    const double p3 = -0.2250947176e-3;
    const double p4 = 0.6841218299e-5;
    const double p5 = 0.5824238515e-5;
    const double p6 = -0.1045274970e-5;
    const double p7 = 0.8360937017e-7;
    const double p8 = -0.3231081277e-8;
    const double p9 = 0.3657763036e-10;
    const double p10 = 0.6936233982e-12;
    
    double pz = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * (p8 + z * (p9 + z * p10)))))))));
    
    return z / pz;
}

Statistics::Statistics() {
    // Constructor implementation
}

Statistics::~Statistics() {
    shutdown();
}

bool Statistics::initialize() {
    try {
        // No specific initialization needed for now
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Exception during statistics initialization: " << e.what() << std::endl;
        return false;
    }
}

void Statistics::shutdown() {
    // No specific shutdown needed for now
}

double Statistics::mean(const std::vector<double>& data) const {
    return calculateMean(data);
}

double Statistics::median(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate median of empty dataset");
    }
    
    std::vector<double> sortedData = sortData(data);
    
    size_t size = sortedData.size();
    if (size % 2 == 0) {
        // Even number of elements, average the middle two
        return (sortedData[size / 2 - 1] + sortedData[size / 2]) / 2.0;
    } else {
        // Odd number of elements, return the middle one
        return sortedData[size / 2];
    }
}

std::vector<double> Statistics::mode(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate mode of empty dataset");
    }
    
    // Count occurrences of each value
    std::map<double, int> counts;
    for (double value : data) {
        counts[value]++;
    }
    
    // Find the maximum frequency
    int maxCount = 0;
    for (const auto& pair : counts) {
        maxCount = std::max(maxCount, pair.second);
    }
    
    // If all values appear only once, there is no mode
    if (maxCount == 1) {
        return {};
    }
    
    // Collect all values with the maximum frequency
    std::vector<double> modes;
    for (const auto& pair : counts) {
        if (pair.second == maxCount) {
            modes.push_back(pair.first);
        }
    }
    
    return modes;
}

double Statistics::range(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate range of empty dataset");
    }
    
    auto [minIt, maxIt] = std::minmax_element(data.begin(), data.end());
    return *maxIt - *minIt;
}

double Statistics::variance(const std::vector<double>& data, bool sample) const {
    return calculateVariance(data, sample);
}

double Statistics::standardDeviation(const std::vector<double>& data, bool sample) const {
    return std::sqrt(calculateVariance(data, sample));
}

std::optional<double> Statistics::covariance(const std::vector<double>& dataX, 
                                           const std::vector<double>& dataY, 
                                           bool sample) const {
    if (dataX.size() != dataY.size() || dataX.empty()) {
        return std::nullopt;
    }
    
    double meanX = calculateMean(dataX);
    double meanY = calculateMean(dataY);
    
    double sum = 0.0;
    for (size_t i = 0; i < dataX.size(); ++i) {
        sum += (dataX[i] - meanX) * (dataY[i] - meanY);
    }
    
    double divisor = sample ? dataX.size() - 1 : dataX.size();
    return sum / divisor;
}

std::optional<double> Statistics::correlation(const std::vector<double>& dataX, 
                                            const std::vector<double>& dataY) const {
    if (dataX.size() != dataY.size() || dataX.empty()) {
        return std::nullopt;
    }
    
    double covXY = *covariance(dataX, dataY, true);
    double stdDevX = standardDeviation(dataX, true);
    double stdDevY = standardDeviation(dataY, true);
    
    if (stdDevX == 0.0 || stdDevY == 0.0) {
        return std::nullopt; // Cannot calculate correlation if either standard deviation is zero
    }
    
    return covXY / (stdDevX * stdDevY);
}

std::optional<std::pair<double, double>> Statistics::linearRegression(const std::vector<double>& dataX, 
                                                                    const std::vector<double>& dataY) const {
    if (dataX.size() != dataY.size() || dataX.empty()) {
        return std::nullopt;
    }
    
    double meanX = calculateMean(dataX);
    double meanY = calculateMean(dataY);
    
    double numerator = 0.0;
    double denominator = 0.0;
    
    for (size_t i = 0; i < dataX.size(); ++i) {
        double xDiff = dataX[i] - meanX;
        numerator += xDiff * (dataY[i] - meanY);
        denominator += xDiff * xDiff;
    }
    
    if (denominator == 0.0) {
        return std::nullopt; // Cannot calculate slope if all x values are the same
    }
    
    double slope = numerator / denominator;
    double intercept = meanY - slope * meanX;
    
    return std::make_pair(slope, intercept);
}

double Statistics::percentile(const std::vector<double>& data, double percentile) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate percentile of empty dataset");
    }
    
    if (percentile < 0.0 || percentile > 100.0) {
        throw std::invalid_argument("Percentile must be between 0 and 100");
    }
    
    std::vector<double> sortedData = sortData(data);
    
    if (percentile == 0.0) {
        return sortedData.front();
    }
    
    if (percentile == 100.0) {
        return sortedData.back();
    }
    
    double index = percentile / 100.0 * (sortedData.size() - 1);
    size_t lowerIndex = static_cast<size_t>(std::floor(index));
    size_t upperIndex = static_cast<size_t>(std::ceil(index));
    
    if (lowerIndex == upperIndex) {
        return sortedData[lowerIndex];
    }
    
    double weight = index - lowerIndex;
    return sortedData[lowerIndex] * (1.0 - weight) + sortedData[upperIndex] * weight;
}

std::vector<double> Statistics::quartiles(const std::vector<double>& data) const {
    return {
        percentile(data, 25.0),
        percentile(data, 50.0),
        percentile(data, 75.0)
    };
}

double Statistics::interquartileRange(const std::vector<double>& data) const {
    auto q = quartiles(data);
    return q[2] - q[0];
}

double Statistics::zScore(const std::vector<double>& data, double value) const {
    double meanValue = calculateMean(data);
    double stdDev = standardDeviation(data, true);
    
    if (stdDev == 0.0) {
        throw std::invalid_argument("Cannot calculate z-score when standard deviation is zero");
    }
    
    return (value - meanValue) / stdDev;
}

double Statistics::skewness(const std::vector<double>& data) const {
    if (data.size() < 3) {
        throw std::invalid_argument("Cannot calculate skewness with less than 3 data points");
    }
    
    double meanValue = calculateMean(data);
    double stdDev = standardDeviation(data, false);
    
    if (stdDev == 0.0) {
        throw std::invalid_argument("Cannot calculate skewness when standard deviation is zero");
    }
    
    double sum = 0.0;
    for (double value : data) {
        double diff = value - meanValue;
        sum += std::pow(diff, 3);
    }
    
    return sum / (data.size() * std::pow(stdDev, 3));
}

double Statistics::kurtosis(const std::vector<double>& data) const {
    if (data.size() < 4) {
        throw std::invalid_argument("Cannot calculate kurtosis with less than 4 data points");
    }
    
    double meanValue = calculateMean(data);
    double variance = calculateVariance(data, false);
    
    if (variance == 0.0) {
        throw std::invalid_argument("Cannot calculate kurtosis when variance is zero");
    }
    
    double sum = 0.0;
    for (double value : data) {
        double diff = value - meanValue;
        sum += std::pow(diff, 4);
    }
    
    return sum / (data.size() * std::pow(variance, 2)) - 3.0;
}

std::optional<double> Statistics::geometricMean(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate geometric mean of empty dataset");
    }
    
    // Check if any value is negative
    for (double value : data) {
        if (value < 0.0) {
            return std::nullopt;
        }
    }
    
    double logSum = 0.0;
    for (double value : data) {
        if (value == 0.0) {
            return 0.0; // Geometric mean is 0 if any value is 0
        }
        logSum += std::log(value);
    }
    
    return std::exp(logSum / data.size());
}

std::optional<double> Statistics::harmonicMean(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate harmonic mean of empty dataset");
    }
    
    double sum = 0.0;
    for (double value : data) {
        if (value <= 0.0) {
            return std::nullopt; // Cannot calculate harmonic mean if any value is zero or negative
        }
        sum += 1.0 / value;
    }
    
    return data.size() / sum;
}

double Statistics::rootMeanSquare(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate root mean square of empty dataset");
    }
    
    double sum = 0.0;
    for (double value : data) {
        sum += value * value;
    }
    
    return std::sqrt(sum / data.size());
}

double Statistics::meanAbsoluteDeviation(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate mean absolute deviation of empty dataset");
    }
    
    double meanValue = calculateMean(data);
    
    double sum = 0.0;
    for (double value : data) {
        sum += std::abs(value - meanValue);
    }
    
    return sum / data.size();
}

double Statistics::medianAbsoluteDeviation(const std::vector<double>& data) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot calculate median absolute deviation of empty dataset");
    }
    
    double medianValue = median(data);
    
    std::vector<double> deviations;
    deviations.reserve(data.size());
    
    for (double value : data) {
        deviations.push_back(std::abs(value - medianValue));
    }
    
    return median(deviations);
}

std::optional<double> Statistics::factorial(int n) const {
    if (n < 0) {
        return std::nullopt;
    }
    
    if (n > 170) {
        throw std::invalid_argument("Factorial too large to represent as double");
    }
    
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    
    return result;
}

std::optional<double> Statistics::binomialCoefficient(int n, int k) const {
    if (n < 0 || k < 0 || k > n) {
        return std::nullopt;
    }
    
    // Optimize by using the smaller of k and n-k
    k = std::min(k, n - k);
    
    double result = 1.0;
    for (int i = 1; i <= k; ++i) {
        result *= (n - (i - 1));
        result /= i;
    }
    
    return result;
}

std::optional<double> Statistics::binomialPMF(int n, double p, int k) const {
    if (n < 0 || k < 0 || k > n || p < 0.0 || p > 1.0) {
        return std::nullopt;
    }
    
    auto coef = binomialCoefficient(n, k);
    if (!coef) {
        return std::nullopt;
    }
    
    return (*coef) * std::pow(p, k) * std::pow(1.0 - p, n - k);
}

std::optional<double> Statistics::binomialCDF(int n, double p, int k) const {
    if (n < 0 || k < 0 || k > n || p < 0.0 || p > 1.0) {
        return std::nullopt;
    }
    
    double sum = 0.0;
    for (int i = 0; i <= k; ++i) {
        auto pmf = binomialPMF(n, p, i);
        if (!pmf) {
            return std::nullopt;
        }
        sum += *pmf;
    }
    
    return sum;
}

std::optional<double> Statistics::normalPDF(double x, double mean, double stdDev) const {
    if (stdDev <= 0.0) {
        return std::nullopt;
    }
    
    static const double SQRT_2PI = std::sqrt(2.0 * M_PI);
    double z = (x - mean) / stdDev;
    
    return std::exp(-0.5 * z * z) / (stdDev * SQRT_2PI);
}

std::optional<double> Statistics::normalCDF(double x, double mean, double stdDev) const {
    if (stdDev <= 0.0) {
        return std::nullopt;
    }
    
    double z = (x - mean) / stdDev;
    return 0.5 * (1.0 + erf(z / std::sqrt(2.0)));
}

std::optional<double> Statistics::normalInvCDF(double p, double mean, double stdDev) const {
    if (p < 0.0 || p > 1.0 || stdDev <= 0.0) {
        return std::nullopt;
    }
    
    // Approximation of the inverse error function
    double x = erfInv(2.0 * p - 1.0);
    
    return mean + stdDev * std::sqrt(2.0) * x;
}

std::optional<double> Statistics::tPDF(double x, double df) const {
    if (df <= 0.0) {
        return std::nullopt;
    }
    
    static const double PI = M_PI;
    double numerator = std::tgamma((df + 1.0) / 2.0);
    double denominator = std::sqrt(df * PI) * std::tgamma(df / 2.0);
    
    return numerator / denominator * std::pow(1.0 + x * x / df, -(df + 1.0) / 2.0);
}

std::optional<double> Statistics::tCDF(double x, double df) const {
    if (df <= 0.0) {
        return std::nullopt;
    }
    
    // This is a complex calculation that requires numerical integration
    // For simplicity, we'll use an approximation
    
    // For large df, the t-distribution approaches the normal distribution
    if (df > 100.0) {
        return normalCDF(x, 0.0, 1.0);
    }
    
    // For small df, we'll use a simple approximation
    // This is not accurate and should be replaced with a proper implementation
    double p = 0.5 + 0.5 * x / std::sqrt(x * x + df);
    return p;
}

std::optional<double> Statistics::tInvCDF(double p, double df) const {
    if (p < 0.0 || p > 1.0 || df <= 0.0) {
        return std::nullopt;
    }
    
    // This is a complex calculation that requires numerical methods
    // For simplicity, we'll use an approximation
    
    // For large df, the t-distribution approaches the normal distribution
    if (df > 100.0) {
        auto z = normalInvCDF(p, 0.0, 1.0);
        return z;
    }
    
    // For small df, we'll use a simple approximation
    // This is not accurate and should be replaced with a proper implementation
    double t = std::tan(M_PI * (p - 0.5)) * std::sqrt(df / (2.0 - df));
    return t;
}

std::optional<double> Statistics::chiSquaredPDF(double x, double df) const {
    if (df <= 0.0 || x < 0.0) {
        return std::nullopt;
    }
    
    if (x == 0.0 && df == 2.0) {
        return 0.5;
    }
    
    if (x == 0.0) {
        return (df == 2.0) ? 0.5 : 0.0;
    }
    
    double halfDf = df / 2.0;
    return std::pow(x, halfDf - 1.0) * std::exp(-x / 2.0) / (std::pow(2.0, halfDf) * std::tgamma(halfDf));
}

std::optional<double> Statistics::chiSquaredCDF(double x, double df) const {
    if (df <= 0.0 || x < 0.0) {
        return std::nullopt;
    }
    
    if (x == 0.0) {
        return 0.0;
    }
    
    // This is the regularized lower incomplete gamma function
    // For simplicity, we'll use an approximation
    
    // For large df, the chi-squared distribution approaches the normal distribution
    if (df > 100.0) {
        double z = std::sqrt(2.0 * x) - std::sqrt(2.0 * df - 1.0);
        return normalCDF(z, 0.0, 1.0);
    }
    
    // For small df, we'll use a simple approximation
    // This is not accurate and should be replaced with a proper implementation
    double p = 1.0 - std::exp(-x / 2.0);
    return p;
}

std::optional<double> Statistics::chiSquaredInvCDF(double p, double df) const {
    if (p < 0.0 || p > 1.0 || df <= 0.0) {
        return std::nullopt;
    }
    
    // This is a complex calculation that requires numerical methods
    // For simplicity, we'll use an approximation
    
    // For large df, the chi-squared distribution approaches the normal distribution
    if (df > 100.0) {
        auto z = normalInvCDF(p, 0.0, 1.0);
        if (!z) {
            return std::nullopt;
        }
        return ((*z) + std::sqrt(2.0 * df - 1.0)) * ((*z) + std::sqrt(2.0 * df - 1.0)) / 2.0;
    }
    
    // For small df, we'll use a simple approximation
    // This is not accurate and should be replaced with a proper implementation
    double x = -2.0 * std::log(1.0 - p);
    return x;
}

// Beta function implementation
double beta(double a, double b) {
    return std::tgamma(a) * std::tgamma(b) / std::tgamma(a + b);
}

std::optional<double> Statistics::fPDF(double x, double df1, double df2) const {
    if (df1 <= 0.0 || df2 <= 0.0 || x < 0.0) {
        return std::nullopt;
    }
    
    double numerator = std::pow(df1 / df2, df1 / 2.0) * std::pow(x, df1 / 2.0 - 1.0);
    double denominator = beta(df1 / 2.0, df2 / 2.0) * std::pow(1.0 + df1 * x / df2, (df1 + df2) / 2.0);
    
    return numerator / denominator;
}

std::optional<double> Statistics::fCDF(double x, double df1, double df2) const {
    if (df1 <= 0.0 || df2 <= 0.0 || x < 0.0) {
        return std::nullopt;
    }
    
    // This is a complex calculation that requires the regularized incomplete beta function
    // For simplicity, we'll use an approximation
    
    // For large df, the F-distribution approaches the chi-squared distribution
    if (df1 > 100.0 && df2 > 100.0) {
        double z = std::sqrt(2.0 * std::log(x * df1 / df2)) - std::sqrt(2.0 * std::log(df1 / df2));
        return normalCDF(z, 0.0, 1.0);
    }
    
    // For small df, we'll use a simple approximation
    // This is not accurate and should be replaced with a proper implementation
    double p = x / (x + df2 / df1);
    return p;
}

std::optional<double> Statistics::fInvCDF(double p, double df1, double df2) const {
    if (p < 0.0 || p > 1.0 || df1 <= 0.0 || df2 <= 0.0) {
        return std::nullopt;
    }
    
    // This is a complex calculation that requires numerical methods
    // For simplicity, we'll use an approximation
    
    // For large df, the F-distribution approaches the chi-squared distribution
    if (df1 > 100.0 && df2 > 100.0) {
        auto z = normalInvCDF(p, 0.0, 1.0);
        if (!z) {
            return std::nullopt;
        }
        double logF = ((*z) + std::sqrt(2.0 * std::log(df1 / df2))) / std::sqrt(2.0);
        return std::exp(logF) * df2 / df1;
    }
    
    // For small df, we'll use a simple approximation
    // This is not accurate and should be replaced with a proper implementation
    double x = p / (1.0 - p) * df2 / df1;
    return x;
}

std::tuple<double, double, double> Statistics::tTest(const std::vector<double>& data, double mu) const {
    if (data.empty()) {
        throw std::invalid_argument("Cannot perform t-test on empty dataset");
    }
    
    double meanValue = calculateMean(data);
    double stdDev = standardDeviation(data, true);
    double se = stdDev / std::sqrt(data.size());
    
    double tStat = (meanValue - mu) / se;
    double df = data.size() - 1;
    
    // Calculate p-value (two-tailed)
    auto pValue = tCDF(-std::abs(tStat), df);
    if (!pValue) {
        throw std::runtime_error("Failed to calculate p-value for t-test");
    }
    
    return std::make_tuple(tStat, 2.0 * (*pValue), df);
}

std::optional<std::tuple<double, double, double>> Statistics::tTest(const std::vector<double>& data1, 
                                                                  const std::vector<double>& data2, 
                                                                  bool equalVariance) const {
    if (data1.empty() || data2.empty()) {
        return std::nullopt;
    }
    
    double mean1 = calculateMean(data1);
    double mean2 = calculateMean(data2);
    double var1 = calculateVariance(data1, true);
    double var2 = calculateVariance(data2, true);
    
    double tStat, df;
    
    if (equalVariance) {
        // Pooled variance
        double pooledVar = ((data1.size() - 1) * var1 + (data2.size() - 1) * var2) / 
                          (data1.size() + data2.size() - 2);
        
        double se = std::sqrt(pooledVar * (1.0 / data1.size() + 1.0 / data2.size()));
        tStat = (mean1 - mean2) / se;
        df = data1.size() + data2.size() - 2;
    } else {
        // Welch's t-test
        double se = std::sqrt(var1 / data1.size() + var2 / data2.size());
        tStat = (mean1 - mean2) / se;
        
        // Welch-Satterthwaite equation for degrees of freedom
        double term1 = var1 / data1.size();
        double term2 = var2 / data2.size();
        double numerator = std::pow(term1 + term2, 2);
        double denominator = std::pow(term1, 2) / (data1.size() - 1) + 
                            std::pow(term2, 2) / (data2.size() - 1);
        
        df = numerator / denominator;
    }
    
    // Calculate p-value (two-tailed)
    auto pValue = tCDF(-std::abs(tStat), df);
    if (!pValue) {
        return std::nullopt;
    }
    
    return std::make_tuple(tStat, 2.0 * (*pValue), df);
}

std::optional<std::pair<double, double>> Statistics::chiSquaredTest(const std::vector<double>& observed, 
                                                                  size_t rows, 
                                                                  size_t cols) const {
    if (observed.size() != rows * cols || rows < 2 || cols < 2) {
        return std::nullopt;
    }
    
    // Calculate row and column sums
    std::vector<double> rowSums(rows, 0.0);
    std::vector<double> colSums(cols, 0.0);
    double total = 0.0;
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            double value = observed[i * cols + j];
            rowSums[i] += value;
            colSums[j] += value;
            total += value;
        }
    }
    
    // Calculate expected frequencies
    std::vector<double> expected(rows * cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            expected[i * cols + j] = rowSums[i] * colSums[j] / total;
        }
    }
    
    // Calculate chi-squared statistic
    double chiSquared = 0.0;
    for (size_t i = 0; i < rows * cols; ++i) {
        if (expected[i] > 0.0) {
            chiSquared += std::pow(observed[i] - expected[i], 2) / expected[i];
        }
    }
    
    // Degrees of freedom
    double df = (rows - 1) * (cols - 1);
    
    // Calculate p-value
    auto pValue = chiSquaredCDF(chiSquared, df);
    if (!pValue) {
        return std::nullopt;
    }
    
    return std::make_pair(chiSquared, 1.0 - (*pValue));
}

std::optional<std::tuple<double, double, double, double>> Statistics::anovaTest(const std::vector<std::vector<double>>& groups) const {
    if (groups.size() < 2) {
        return std::nullopt;
    }
    
    // Calculate group means and sizes
    std::vector<double> groupMeans;
    std::vector<size_t> groupSizes;
    double grandTotal = 0.0;
    size_t totalSize = 0;
    
    for (const auto& group : groups) {
        if (group.empty()) {
            return std::nullopt;
        }
        
        double mean = calculateMean(group);
        groupMeans.push_back(mean);
        groupSizes.push_back(group.size());
        
        grandTotal += std::accumulate(group.begin(), group.end(), 0.0);
        totalSize += group.size();
    }
    
    double grandMean = grandTotal / totalSize;
    
    // Calculate between-group sum of squares
    double ssb = 0.0;
    for (size_t i = 0; i < groups.size(); ++i) {
        ssb += groupSizes[i] * std::pow(groupMeans[i] - grandMean, 2);
    }
    
    // Calculate within-group sum of squares
    double ssw = 0.0;
    for (size_t i = 0; i < groups.size(); ++i) {
        for (double value : groups[i]) {
            ssw += std::pow(value - groupMeans[i], 2);
        }
    }
    
    // Degrees of freedom
    double dfb = groups.size() - 1;
    double dfw = totalSize - groups.size();
    
    // Mean squares
    double msb = ssb / dfb;
    double msw = ssw / dfw;
    
    // F-statistic
    double fStat = msb / msw;
    
    // Calculate p-value
    auto pValue = fCDF(fStat, dfb, dfw);
    if (!pValue) {
        return std::nullopt;
    }
    
    return std::make_tuple(fStat, 1.0 - (*pValue), dfb, dfw);
}

std::unordered_map<std::string, std::string> Statistics::getAllFunctions() const {
    return {
        {"mean", "Calculate the mean (average) of a dataset"},
        {"median", "Calculate the median of a dataset"},
        {"mode", "Calculate the mode of a dataset"},
        {"range", "Calculate the range of a dataset"},
        {"variance", "Calculate the variance of a dataset"},
        {"standardDeviation", "Calculate the standard deviation of a dataset"},
        {"covariance", "Calculate the covariance of two datasets"},
        {"correlation", "Calculate the correlation coefficient of two datasets"},
        {"linearRegression", "Perform linear regression on two datasets"},
        {"percentile", "Calculate the percentile of a dataset"},
        {"quartiles", "Calculate the quartiles of a dataset"},
        {"interquartileRange", "Calculate the interquartile range (IQR) of a dataset"},
        {"zScore", "Calculate the z-score of a value in a dataset"},
        {"skewness", "Calculate the skewness of a dataset"},
        {"kurtosis", "Calculate the kurtosis of a dataset"},
        {"geometricMean", "Calculate the geometric mean of a dataset"},
        {"harmonicMean", "Calculate the harmonic mean of a dataset"},
        {"rootMeanSquare", "Calculate the root mean square (RMS) of a dataset"},
        {"meanAbsoluteDeviation", "Calculate the mean absolute deviation of a dataset"},
        {"medianAbsoluteDeviation", "Calculate the median absolute deviation of a dataset"},
        {"factorial", "Calculate the factorial of a non-negative integer"},
        {"binomialCoefficient", "Calculate the binomial coefficient (n choose k)"},
        {"binomialPMF", "Calculate the probability mass function (PMF) of a binomial distribution"},
        {"binomialCDF", "Calculate the cumulative distribution function (CDF) of a binomial distribution"},
        {"normalPDF", "Calculate the probability density function (PDF) of a normal distribution"},
        {"normalCDF", "Calculate the cumulative distribution function (CDF) of a normal distribution"},
        {"normalInvCDF", "Calculate the inverse cumulative distribution function (quantile function) of a normal distribution"},
        {"tPDF", "Calculate the probability density function (PDF) of a t-distribution"},
        {"tCDF", "Calculate the cumulative distribution function (CDF) of a t-distribution"},
        {"tInvCDF", "Calculate the inverse cumulative distribution function (quantile function) of a t-distribution"},
        {"chiSquaredPDF", "Calculate the probability density function (PDF) of a chi-squared distribution"},
        {"chiSquaredCDF", "Calculate the cumulative distribution function (CDF) of a chi-squared distribution"},
        {"chiSquaredInvCDF", "Calculate the inverse cumulative distribution function (quantile function) of a chi-squared distribution"},
        {"fPDF", "Calculate the probability density function (PDF) of an F-distribution"},
        {"fCDF", "Calculate the cumulative distribution function (CDF) of an F-distribution"},
        {"fInvCDF", "Calculate the inverse cumulative distribution function (quantile function) of an F-distribution"},
        {"tTest", "Perform a t-test"},
        {"chiSquaredTest", "Perform a chi-squared test for independence"},
        {"anovaTest", "Perform a one-way ANOVA test"}
    };
}
