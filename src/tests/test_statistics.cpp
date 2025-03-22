#include <gtest/gtest.h>
#include "../backend/statistics.h"
#include <vector>
#include <cmath>
#include <algorithm>

using namespace rebelcalc;

class StatisticsTest : public ::testing::Test {
protected:
    Statistics stats;
    
    void SetUp() override {
        stats.initialize();
    }
    
    void TearDown() override {
        stats.shutdown();
    }
};

// Helper function to check if two doubles are approximately equal
bool approxEqual(double a, double b, double epsilon = 1e-6) {
    return std::abs(a - b) < epsilon;
}

// Helper function to check if two vectors of doubles are approximately equal
bool approxEqualVectors(const std::vector<double>& a, const std::vector<double>& b, double epsilon = 1e-6) {
    if (a.size() != b.size()) {
        return false;
    }
    
    for (size_t i = 0; i < a.size(); ++i) {
        if (!approxEqual(a[i], b[i], epsilon)) {
            return false;
        }
    }
    
    return true;
}

TEST_F(StatisticsTest, MeanTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.mean(data), 3.0));
    
    data = {-1.0, 0.0, 1.0};
    EXPECT_TRUE(approxEqual(stats.mean(data), 0.0));
    
    data = {10.0};
    EXPECT_TRUE(approxEqual(stats.mean(data), 10.0));
    
    EXPECT_THROW(stats.mean({}), std::invalid_argument);
}

TEST_F(StatisticsTest, MedianTest) {
    std::vector<double> data = {1.0, 3.0, 5.0, 7.0, 9.0};
    EXPECT_TRUE(approxEqual(stats.median(data), 5.0));
    
    data = {1.0, 3.0, 5.0, 7.0};
    EXPECT_TRUE(approxEqual(stats.median(data), 4.0));
    
    data = {10.0};
    EXPECT_TRUE(approxEqual(stats.median(data), 10.0));
    
    EXPECT_THROW(stats.median({}), std::invalid_argument);
}

TEST_F(StatisticsTest, ModeTest) {
    std::vector<double> data = {1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0};
    auto modes = stats.mode(data);
    ASSERT_EQ(modes.size(), 1);
    EXPECT_TRUE(approxEqual(modes[0], 4.0));
    
    data = {1.0, 1.0, 2.0, 2.0, 3.0, 3.0};
    modes = stats.mode(data);
    ASSERT_EQ(modes.size(), 3);
    std::sort(modes.begin(), modes.end());
    EXPECT_TRUE(approxEqual(modes[0], 1.0));
    EXPECT_TRUE(approxEqual(modes[1], 2.0));
    EXPECT_TRUE(approxEqual(modes[2], 3.0));
    
    data = {1.0, 2.0, 3.0, 4.0, 5.0};
    modes = stats.mode(data);
    EXPECT_TRUE(modes.empty());
    
    EXPECT_THROW(stats.mode({}), std::invalid_argument);
}

TEST_F(StatisticsTest, RangeTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.range(data), 4.0));
    
    data = {-10.0, 0.0, 10.0};
    EXPECT_TRUE(approxEqual(stats.range(data), 20.0));
    
    data = {5.0};
    EXPECT_TRUE(approxEqual(stats.range(data), 0.0));
    
    EXPECT_THROW(stats.range({}), std::invalid_argument);
}

TEST_F(StatisticsTest, VarianceTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.variance(data, true), 2.5));
    EXPECT_TRUE(approxEqual(stats.variance(data, false), 2.0));
    
    data = {2.0, 2.0, 2.0, 2.0};
    EXPECT_TRUE(approxEqual(stats.variance(data, true), 0.0));
    EXPECT_TRUE(approxEqual(stats.variance(data, false), 0.0));
    
    data = {10.0};
    EXPECT_THROW(stats.variance(data, true), std::invalid_argument);
    EXPECT_TRUE(approxEqual(stats.variance(data, false), 0.0));
    
    EXPECT_THROW(stats.variance({}), std::invalid_argument);
}

TEST_F(StatisticsTest, StandardDeviationTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.standardDeviation(data, true), std::sqrt(2.5)));
    EXPECT_TRUE(approxEqual(stats.standardDeviation(data, false), std::sqrt(2.0)));
    
    data = {2.0, 2.0, 2.0, 2.0};
    EXPECT_TRUE(approxEqual(stats.standardDeviation(data, true), 0.0));
    EXPECT_TRUE(approxEqual(stats.standardDeviation(data, false), 0.0));
    
    data = {10.0};
    EXPECT_THROW(stats.standardDeviation(data, true), std::invalid_argument);
    EXPECT_TRUE(approxEqual(stats.standardDeviation(data, false), 0.0));
    
    EXPECT_THROW(stats.standardDeviation({}), std::invalid_argument);
}

TEST_F(StatisticsTest, CovarianceTest) {
    std::vector<double> dataX = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> dataY = {5.0, 4.0, 3.0, 2.0, 1.0};
    
    auto result = stats.covariance(dataX, dataY, true);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, -2.5));
    
    result = stats.covariance(dataX, dataY, false);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, -2.0));
    
    dataY = {2.0, 4.0, 6.0, 8.0, 10.0};
    result = stats.covariance(dataX, dataY, true);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 2.5));
    
    std::vector<double> emptyData;
    EXPECT_FALSE(stats.covariance(dataX, emptyData, true).has_value());
    EXPECT_FALSE(stats.covariance(emptyData, dataY, true).has_value());
    
    std::vector<double> differentSizeData = {1.0, 2.0};
    EXPECT_FALSE(stats.covariance(dataX, differentSizeData, true).has_value());
}

TEST_F(StatisticsTest, CorrelationTest) {
    std::vector<double> dataX = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> dataY = {5.0, 4.0, 3.0, 2.0, 1.0};
    
    auto result = stats.correlation(dataX, dataY);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, -1.0));
    
    dataY = {2.0, 4.0, 6.0, 8.0, 10.0};
    result = stats.correlation(dataX, dataY);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    dataY = {5.0, 5.0, 5.0, 5.0, 5.0};
    result = stats.correlation(dataX, dataY);
    EXPECT_FALSE(result.has_value()); // Zero standard deviation in dataY
    
    std::vector<double> emptyData;
    EXPECT_FALSE(stats.correlation(dataX, emptyData).has_value());
    EXPECT_FALSE(stats.correlation(emptyData, dataY).has_value());
    
    std::vector<double> differentSizeData = {1.0, 2.0};
    EXPECT_FALSE(stats.correlation(dataX, differentSizeData).has_value());
}

TEST_F(StatisticsTest, LinearRegressionTest) {
    std::vector<double> dataX = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> dataY = {2.0, 4.0, 6.0, 8.0, 10.0};
    
    auto result = stats.linearRegression(dataX, dataY);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(result->first, 2.0)); // Slope
    EXPECT_TRUE(approxEqual(result->second, 0.0)); // Intercept
    
    dataY = {3.0, 5.0, 7.0, 9.0, 11.0};
    result = stats.linearRegression(dataX, dataY);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(result->first, 2.0)); // Slope
    EXPECT_TRUE(approxEqual(result->second, 1.0)); // Intercept
    
    dataX = {1.0, 1.0, 1.0, 1.0, 1.0};
    result = stats.linearRegression(dataX, dataY);
    EXPECT_FALSE(result.has_value()); // All x values are the same
    
    std::vector<double> emptyData;
    EXPECT_FALSE(stats.linearRegression(dataX, emptyData).has_value());
    EXPECT_FALSE(stats.linearRegression(emptyData, dataY).has_value());
    
    std::vector<double> differentSizeData = {1.0, 2.0};
    EXPECT_FALSE(stats.linearRegression(dataX, differentSizeData).has_value());
}

TEST_F(StatisticsTest, PercentileTest) {
    std::vector<double> data = {15.0, 20.0, 35.0, 40.0, 50.0};
    
    EXPECT_TRUE(approxEqual(stats.percentile(data, 0.0), 15.0));
    EXPECT_TRUE(approxEqual(stats.percentile(data, 25.0), 20.0));
    EXPECT_TRUE(approxEqual(stats.percentile(data, 50.0), 35.0));
    EXPECT_TRUE(approxEqual(stats.percentile(data, 75.0), 40.0));
    EXPECT_TRUE(approxEqual(stats.percentile(data, 100.0), 50.0));
    
    // Test interpolation
    EXPECT_TRUE(approxEqual(stats.percentile(data, 10.0), 17.0));
    
    EXPECT_THROW(stats.percentile({}, 50.0), std::invalid_argument);
    EXPECT_THROW(stats.percentile(data, -10.0), std::invalid_argument);
    EXPECT_THROW(stats.percentile(data, 110.0), std::invalid_argument);
}

TEST_F(StatisticsTest, QuartilesTest) {
    std::vector<double> data = {15.0, 20.0, 35.0, 40.0, 50.0};
    
    auto result = stats.quartiles(data);
    ASSERT_EQ(result.size(), 3);
    EXPECT_TRUE(approxEqual(result[0], 20.0)); // Q1
    EXPECT_TRUE(approxEqual(result[1], 35.0)); // Q2 (median)
    EXPECT_TRUE(approxEqual(result[2], 40.0)); // Q3
    
    EXPECT_THROW(stats.quartiles({}), std::invalid_argument);
}

TEST_F(StatisticsTest, InterquartileRangeTest) {
    std::vector<double> data = {15.0, 20.0, 35.0, 40.0, 50.0};
    
    EXPECT_TRUE(approxEqual(stats.interquartileRange(data), 20.0));
    
    EXPECT_THROW(stats.interquartileRange({}), std::invalid_argument);
}

TEST_F(StatisticsTest, ZScoreTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    
    EXPECT_TRUE(approxEqual(stats.zScore(data, 3.0), 0.0));
    EXPECT_TRUE(approxEqual(stats.zScore(data, 5.0), 1.264911064));
    EXPECT_TRUE(approxEqual(stats.zScore(data, 1.0), -1.264911064));
    
    data = {2.0, 2.0, 2.0, 2.0};
    EXPECT_THROW(stats.zScore(data, 3.0), std::invalid_argument); // Zero standard deviation
    
    EXPECT_THROW(stats.zScore({}, 3.0), std::invalid_argument);
}

TEST_F(StatisticsTest, SkewnessTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.skewness(data), 0.0));
    
    data = {1.0, 1.0, 1.0, 1.0, 10.0};
    EXPECT_GT(stats.skewness(data), 0.0); // Positive skewness
    
    data = {1.0, 10.0, 10.0, 10.0, 10.0};
    EXPECT_LT(stats.skewness(data), 0.0); // Negative skewness
    
    data = {2.0, 2.0, 2.0, 2.0};
    EXPECT_THROW(stats.skewness(data), std::invalid_argument); // Zero standard deviation
    
    EXPECT_THROW(stats.skewness({}), std::invalid_argument);
    EXPECT_THROW(stats.skewness({1.0, 2.0}), std::invalid_argument); // Less than 3 data points
}

TEST_F(StatisticsTest, KurtosisTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.kurtosis(data), -1.3)); // Platykurtic
    
    data = {1.0, 3.0, 3.0, 3.0, 5.0};
    EXPECT_LT(stats.kurtosis(data), 0.0); // Platykurtic
    
    data = {1.0, 1.0, 3.0, 5.0, 5.0};
    EXPECT_GT(stats.kurtosis(data), 0.0); // Leptokurtic
    
    data = {2.0, 2.0, 2.0, 2.0};
    EXPECT_THROW(stats.kurtosis(data), std::invalid_argument); // Zero variance
    
    EXPECT_THROW(stats.kurtosis({}), std::invalid_argument);
    EXPECT_THROW(stats.kurtosis({1.0, 2.0, 3.0}), std::invalid_argument); // Less than 4 data points
}

TEST_F(StatisticsTest, GeometricMeanTest) {
    std::vector<double> data = {1.0, 2.0, 4.0, 8.0};
    
    auto result = stats.geometricMean(data);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 2.8284271247));
    
    data = {0.0, 1.0, 2.0, 3.0};
    result = stats.geometricMean(data);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.0)); // Zero in the dataset
    
    data = {-1.0, 2.0, 3.0, 4.0};
    result = stats.geometricMean(data);
    EXPECT_FALSE(result.has_value()); // Negative value in the dataset
    
    EXPECT_THROW(stats.geometricMean({}), std::invalid_argument);
}

TEST_F(StatisticsTest, HarmonicMeanTest) {
    std::vector<double> data = {1.0, 2.0, 4.0, 8.0};
    
    auto result = stats.harmonicMean(data);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 2.0));
    
    data = {0.0, 1.0, 2.0, 3.0};
    result = stats.harmonicMean(data);
    EXPECT_FALSE(result.has_value()); // Zero in the dataset
    
    data = {-1.0, 2.0, 3.0, 4.0};
    result = stats.harmonicMean(data);
    EXPECT_FALSE(result.has_value()); // Negative value in the dataset
    
    EXPECT_THROW(stats.harmonicMean({}), std::invalid_argument);
}

TEST_F(StatisticsTest, RootMeanSquareTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.rootMeanSquare(data), 3.3166247904));
    
    data = {-1.0, -2.0, -3.0, -4.0, -5.0};
    EXPECT_TRUE(approxEqual(stats.rootMeanSquare(data), 3.3166247904));
    
    data = {0.0, 0.0, 0.0, 0.0};
    EXPECT_TRUE(approxEqual(stats.rootMeanSquare(data), 0.0));
    
    EXPECT_THROW(stats.rootMeanSquare({}), std::invalid_argument);
}

TEST_F(StatisticsTest, MeanAbsoluteDeviationTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.meanAbsoluteDeviation(data), 1.2));
    
    data = {3.0, 3.0, 3.0, 3.0, 3.0};
    EXPECT_TRUE(approxEqual(stats.meanAbsoluteDeviation(data), 0.0));
    
    EXPECT_THROW(stats.meanAbsoluteDeviation({}), std::invalid_argument);
}

TEST_F(StatisticsTest, MedianAbsoluteDeviationTest) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_TRUE(approxEqual(stats.medianAbsoluteDeviation(data), 1.0));
    
    data = {3.0, 3.0, 3.0, 3.0, 3.0};
    EXPECT_TRUE(approxEqual(stats.medianAbsoluteDeviation(data), 0.0));
    
    EXPECT_THROW(stats.medianAbsoluteDeviation({}), std::invalid_argument);
}

TEST_F(StatisticsTest, FactorialTest) {
    auto result = stats.factorial(0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    result = stats.factorial(5);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 120.0));
    
    result = stats.factorial(10);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 3628800.0));
    
    result = stats.factorial(-1);
    EXPECT_FALSE(result.has_value());
    
    EXPECT_THROW(stats.factorial(171), std::invalid_argument); // Too large
}

TEST_F(StatisticsTest, BinomialCoefficientTest) {
    auto result = stats.binomialCoefficient(5, 2);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 10.0));
    
    result = stats.binomialCoefficient(10, 0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    result = stats.binomialCoefficient(10, 10);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    result = stats.binomialCoefficient(-1, 2);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialCoefficient(5, -1);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialCoefficient(5, 6);
    EXPECT_FALSE(result.has_value());
}

TEST_F(StatisticsTest, BinomialPMFTest) {
    auto result = stats.binomialPMF(10, 0.5, 5);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.24609375));
    
    result = stats.binomialPMF(10, 0.5, 0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.0009765625));
    
    result = stats.binomialPMF(10, 0.5, 10);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.0009765625));
    
    result = stats.binomialPMF(-1, 0.5, 5);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialPMF(10, -0.1, 5);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialPMF(10, 1.1, 5);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialPMF(10, 0.5, -1);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialPMF(10, 0.5, 11);
    EXPECT_FALSE(result.has_value());
}

TEST_F(StatisticsTest, BinomialCDFTest) {
    auto result = stats.binomialCDF(10, 0.5, 5);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.623046875));
    
    result = stats.binomialCDF(10, 0.5, 0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.0009765625));
    
    result = stats.binomialCDF(10, 0.5, 10);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    result = stats.binomialCDF(-1, 0.5, 5);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialCDF(10, -0.1, 5);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialCDF(10, 1.1, 5);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialCDF(10, 0.5, -1);
    EXPECT_FALSE(result.has_value());
    
    result = stats.binomialCDF(10, 0.5, 11);
    EXPECT_FALSE(result.has_value());
}

TEST_F(StatisticsTest, NormalPDFTest) {
    auto result = stats.normalPDF(0.0, 0.0, 1.0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.3989422804));
    
    result = stats.normalPDF(1.0, 0.0, 1.0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.2419707245));
    
    result = stats.normalPDF(0.0, 1.0, 2.0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.1760326633));
    
    result = stats.normalPDF(0.0, 0.0, 0.0);
    EXPECT_FALSE(result.has_value());
    
    result = stats.normalPDF(0.0, 0.0, -1.0);
    EXPECT_FALSE(result.has_value());
}

TEST_F(StatisticsTest, NormalCDFTest) {
    auto result = stats.normalCDF(0.0, 0.0, 1.0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.5));
    
    result = stats.normalCDF(1.0, 0.0, 1.0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.8413447461));
    
    result = stats.normalCDF(-1.0, 0.0, 1.0);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.1586552539));
    
    result = stats.normalCDF(0.0, 0.0, 0.0);
    EXPECT_FALSE(result.has_value());
    
    result = stats.normalCDF(0.0, 0.0, -1.0);
    EXPECT_FALSE(result.has_value());
}

TEST_F(StatisticsTest, GetAllFunctionsTest) {
    auto functions = stats.getAllFunctions();
    EXPECT_FALSE(functions.empty());
    EXPECT_TRUE(functions.find("mean") != functions.end());
    EXPECT_TRUE(functions.find("median") != functions.end());
    EXPECT_TRUE(functions.find("mode") != functions.end());
    EXPECT_TRUE(functions.find("variance") != functions.end());
    EXPECT_TRUE(functions.find("standardDeviation") != functions.end());
}
