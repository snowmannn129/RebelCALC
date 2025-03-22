#include <gtest/gtest.h>
#include "../backend/units.h"
#include <cmath>
#include <limits>

namespace rebelcalc {
namespace tests {

class UnitConverterTest : public ::testing::Test {
protected:
    void SetUp() override {
        converter.initialize();
    }

    void TearDown() override {
        converter.shutdown();
    }

    UnitConverter converter;
    
    // Helper function to check if two doubles are approximately equal
    bool approxEqual(double a, double b, double epsilon = 1e-9) {
        return std::fabs(a - b) < epsilon;
    }
};

// Test initialization
TEST_F(UnitConverterTest, Initialization) {
    UnitConverter converter;
    EXPECT_TRUE(converter.initialize());
}

// Test length conversions
TEST_F(UnitConverterTest, LengthConversions) {
    // Meter to kilometer
    auto result = converter.convert(1000.0, "m", "km");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    // Kilometer to meter
    result = converter.convert(1.0, "km", "m");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1000.0));
    
    // Meter to centimeter
    result = converter.convert(1.0, "m", "cm");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 100.0));
    
    // Inch to centimeter
    result = converter.convert(1.0, "in", "cm");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 2.54));
    
    // Mile to kilometer
    result = converter.convert(1.0, "mi", "km");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.609344));
    
    // Foot to meter
    result = converter.convert(3.0, "ft", "m");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.9144));
}

// Test mass conversions
TEST_F(UnitConverterTest, MassConversions) {
    // Kilogram to gram
    auto result = converter.convert(1.0, "kg", "g");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1000.0));
    
    // Pound to kilogram
    result = converter.convert(1.0, "lb", "kg");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.45359237));
    
    // Metric ton to kilogram
    result = converter.convert(1.0, "t", "kg");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1000.0));
}

// Test temperature conversions
TEST_F(UnitConverterTest, TemperatureConversions) {
    // Celsius to Kelvin
    auto result = converter.convert(0.0, "C", "K");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 273.15));
    
    // Kelvin to Celsius
    result = converter.convert(273.15, "K", "C");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.0));
    
    // Fahrenheit to Celsius
    result = converter.convert(32.0, "F", "C");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.0));
    
    // Celsius to Fahrenheit
    result = converter.convert(100.0, "C", "F");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 212.0));
}

// Test area conversions
TEST_F(UnitConverterTest, AreaConversions) {
    // Square meter to square centimeter
    auto result = converter.convert(1.0, "m2", "cm2");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 10000.0));
    
    // Square foot to square meter
    result = converter.convert(10.0, "ft2", "m2");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.9290304));
    
    // Acre to square meter
    result = converter.convert(1.0, "acre", "m2");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 4046.8564224));
}

// Test volume conversions
TEST_F(UnitConverterTest, VolumeConversions) {
    // Cubic meter to liter
    auto result = converter.convert(1.0, "m3", "L");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1000.0));
    
    // Gallon (US) to liter
    result = converter.convert(1.0, "gal", "L");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 3.785411784));
    
    // Cubic foot to cubic meter
    result = converter.convert(1.0, "ft3", "m3");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.028316846592));
}

// Test velocity conversions
TEST_F(UnitConverterTest, VelocityConversions) {
    // Meter per second to kilometer per hour
    auto result = converter.convert(1.0, "m/s", "km/h");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 3.6));
    
    // Mile per hour to kilometer per hour
    result = converter.convert(60.0, "mph", "km/h");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 96.56064));
    
    // Knot to meter per second
    result = converter.convert(1.0, "kn", "m/s");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 0.514444444));
}

// Test energy conversions
TEST_F(UnitConverterTest, EnergyConversions) {
    // Joule to kilojoule
    auto result = converter.convert(1000.0, "J", "kJ");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    // Kilowatt-hour to joule
    result = converter.convert(1.0, "kWh", "J");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 3600000.0));
    
    // Calorie to joule
    result = converter.convert(1.0, "cal", "J");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 4.184));
}

// Test power conversions
TEST_F(UnitConverterTest, PowerConversions) {
    // Watt to kilowatt
    auto result = converter.convert(1000.0, "W", "kW");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    // Horsepower to watt
    result = converter.convert(1.0, "hp", "W");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 745.699872));
}

// Test pressure conversions
TEST_F(UnitConverterTest, PressureConversions) {
    // Pascal to kilopascal
    auto result = converter.convert(1000.0, "Pa", "kPa");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 1.0));
    
    // Atmosphere to pascal
    result = converter.convert(1.0, "atm", "Pa");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 101325.0));
    
    // PSI to kilopascal
    result = converter.convert(1.0, "psi", "kPa");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 6.89475729317));
}

// Test angle conversions
TEST_F(UnitConverterTest, AngleConversions) {
    // Radian to degree
    auto result = converter.convert(M_PI, "rad", "deg");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 180.0));
    
    // Degree to radian
    result = converter.convert(180.0, "deg", "rad");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, M_PI));
    
    // Degree to arcminute
    result = converter.convert(1.0, "deg", "arcmin");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(approxEqual(*result, 60.0));
}

// Test invalid conversions
TEST_F(UnitConverterTest, InvalidConversions) {
    // Invalid unit
    auto result = converter.convert(1.0, "invalid", "m");
    EXPECT_FALSE(result.has_value());
    
    // Invalid conversion (different quantities)
    result = converter.convert(1.0, "m", "kg");
    EXPECT_FALSE(result.has_value());
}

// Test unit validation
TEST_F(UnitConverterTest, UnitValidation) {
    EXPECT_TRUE(converter.isValidUnit("m"));
    EXPECT_TRUE(converter.isValidUnit("kg"));
    EXPECT_TRUE(converter.isValidUnit("s"));
    EXPECT_FALSE(converter.isValidUnit("invalid"));
}

// Test quantity validation
TEST_F(UnitConverterTest, QuantityValidation) {
    EXPECT_TRUE(converter.isValidQuantity("length"));
    EXPECT_TRUE(converter.isValidQuantity("mass"));
    EXPECT_TRUE(converter.isValidQuantity("time"));
    EXPECT_FALSE(converter.isValidQuantity("invalid"));
}

// Test getting units for a quantity
TEST_F(UnitConverterTest, GetUnitsForQuantity) {
    auto units = converter.getUnitsForQuantity("length");
    EXPECT_FALSE(units.empty());
    EXPECT_TRUE(std::find(units.begin(), units.end(), "m") != units.end());
    EXPECT_TRUE(std::find(units.begin(), units.end(), "km") != units.end());
    EXPECT_TRUE(std::find(units.begin(), units.end(), "cm") != units.end());
    
    // Invalid quantity
    units = converter.getUnitsForQuantity("invalid");
    EXPECT_TRUE(units.empty());
}

// Test getting all quantities
TEST_F(UnitConverterTest, GetAllQuantities) {
    auto quantities = converter.getAllQuantities();
    EXPECT_FALSE(quantities.empty());
    EXPECT_TRUE(std::find(quantities.begin(), quantities.end(), "length") != quantities.end());
    EXPECT_TRUE(std::find(quantities.begin(), quantities.end(), "mass") != quantities.end());
    EXPECT_TRUE(std::find(quantities.begin(), quantities.end(), "time") != quantities.end());
}

// Test getting base unit for a quantity
TEST_F(UnitConverterTest, GetBaseUnit) {
    EXPECT_EQ(converter.getBaseUnit("length"), "m");
    EXPECT_EQ(converter.getBaseUnit("mass"), "kg");
    EXPECT_EQ(converter.getBaseUnit("time"), "s");
    EXPECT_EQ(converter.getBaseUnit("invalid"), "");
}

// Test getting quantity for a unit
TEST_F(UnitConverterTest, GetQuantityForUnit) {
    EXPECT_EQ(converter.getQuantityForUnit("m"), "length");
    EXPECT_EQ(converter.getQuantityForUnit("kg"), "mass");
    EXPECT_EQ(converter.getQuantityForUnit("s"), "time");
    EXPECT_EQ(converter.getQuantityForUnit("invalid"), "");
}

// Test formatting with unit
TEST_F(UnitConverterTest, FormatWithUnit) {
    EXPECT_EQ(converter.formatWithUnit(1.0, "m", 0), "1 m");
    EXPECT_EQ(converter.formatWithUnit(1.5, "kg", 1), "1.5 kg");
    EXPECT_EQ(converter.formatWithUnit(1.234, "s", 2), "1.23 s");
    EXPECT_EQ(converter.formatWithUnit(1.0, "invalid", 0), "");
}

// Test registering a custom unit
TEST_F(UnitConverterTest, RegisterCustomUnit) {
    // Register a custom unit (light-second)
    bool result = converter.registerUnit("ls", "length",
                                       [](double v) { return v * 299792458.0; }, // to base (m)
                                       [](double v) { return v / 299792458.0; }, // from base (m)
                                       "ls", "light-second");
    EXPECT_TRUE(result);
    
    // Test conversion
    auto convResult = converter.convert(1.0, "ls", "m");
    ASSERT_TRUE(convResult.has_value());
    EXPECT_TRUE(approxEqual(*convResult, 299792458.0));
    
    // Test reverse conversion
    convResult = converter.convert(299792458.0, "m", "ls");
    ASSERT_TRUE(convResult.has_value());
    EXPECT_TRUE(approxEqual(*convResult, 1.0));
    
    // Try to register the same unit again (should fail)
    result = converter.registerUnit("ls", "length",
                                  [](double v) { return v; },
                                  [](double v) { return v; },
                                  "ls", "light-second");
    EXPECT_FALSE(result);
    
    // Try to register a unit for an invalid quantity (should fail)
    result = converter.registerUnit("invalid", "invalid_quantity",
                                  [](double v) { return v; },
                                  [](double v) { return v; },
                                  "invalid", "invalid");
    EXPECT_FALSE(result);
}

} // namespace tests
} // namespace rebelcalc
