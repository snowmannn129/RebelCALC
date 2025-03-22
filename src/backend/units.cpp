#include "units.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <regex>

// Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace rebelcalc {

UnitConverter::UnitConverter() {
    // Constructor
}

UnitConverter::~UnitConverter() {
    shutdown();
}

bool UnitConverter::initialize() {
    try {
        // Register all units and quantities
        registerLengthUnits();
        registerMassUnits();
        registerTimeUnits();
        registerTemperatureUnits();
        registerAreaUnits();
        registerVolumeUnits();
        registerVelocityUnits();
        registerAccelerationUnits();
        registerForceUnits();
        registerEnergyUnits();
        registerPowerUnits();
        registerPressureUnits();
        registerElectricCurrentUnits();
        registerElectricVoltageUnits();
        registerElectricResistanceUnits();
        registerFrequencyUnits();
        registerDataUnits();
        registerAngleUnits();
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Exception during unit converter initialization: " << e.what() << std::endl;
        return false;
    }
}

void UnitConverter::shutdown() {
    // Clear all registered units and quantities
    m_units.clear();
    m_quantities.clear();
}

std::optional<double> UnitConverter::convert(double value, const std::string& fromUnit, const std::string& toUnit) {
    // Check if units are valid
    if (!isValidUnit(fromUnit) || !isValidUnit(toUnit)) {
        return std::nullopt;
    }
    
    // Get unit information
    const auto& fromUnitInfo = m_units[fromUnit];
    const auto& toUnitInfo = m_units[toUnit];
    
    // Check if units are of the same quantity
    if (fromUnitInfo.quantity != toUnitInfo.quantity) {
        return std::nullopt;
    }
    
    // Convert from source unit to base unit, then from base unit to target unit
    double baseValue = fromUnitInfo.toBase(value);
    double result = toUnitInfo.fromBase(baseValue);
    
    return result;
}

std::vector<std::string> UnitConverter::getUnitsForQuantity(const std::string& quantity) const {
    if (!isValidQuantity(quantity)) {
        return {};
    }
    
    return m_quantities.at(quantity).units;
}

std::vector<std::string> UnitConverter::getAllQuantities() const {
    std::vector<std::string> quantities;
    quantities.reserve(m_quantities.size());
    
    for (const auto& [quantity, _] : m_quantities) {
        quantities.push_back(quantity);
    }
    
    return quantities;
}

std::string UnitConverter::getBaseUnit(const std::string& quantity) const {
    if (!isValidQuantity(quantity)) {
        return "";
    }
    
    return m_quantities.at(quantity).baseUnit;
}

std::string UnitConverter::getQuantityForUnit(const std::string& unit) const {
    if (!isValidUnit(unit)) {
        return "";
    }
    
    return m_units.at(unit).quantity;
}

bool UnitConverter::isValidUnit(const std::string& unit) const {
    return m_units.find(unit) != m_units.end();
}

bool UnitConverter::isValidQuantity(const std::string& quantity) const {
    return m_quantities.find(quantity) != m_quantities.end();
}

bool UnitConverter::registerUnit(const std::string& unit, 
                               const std::string& quantity,
                               std::function<double(double)> toBaseFunc,
                               std::function<double(double)> fromBaseFunc) {
    // Check if quantity is valid
    if (!isValidQuantity(quantity)) {
        return false;
    }
    
    // Check if unit already exists
    if (isValidUnit(unit)) {
        return false;
    }
    
    // Add the unit
    UnitInfo unitInfo;
    unitInfo.quantity = quantity;
    unitInfo.toBase = toBaseFunc;
    unitInfo.fromBase = fromBaseFunc;
    unitInfo.symbol = unit;
    unitInfo.name = unit; // Default name is the same as the unit symbol
    
    m_units[unit] = unitInfo;
    
    // Add the unit to the quantity's unit list
    m_quantities[quantity].units.push_back(unit);
    
    return true;
}

std::string UnitConverter::parseUnitExpression(const std::string& expression) const {
    // This is a placeholder for unit expression parsing
    // In a real implementation, this would parse expressions like "m/s", "kg*m^2", etc.
    // For now, just return the expression if it's a valid unit
    if (isValidUnit(expression)) {
        return expression;
    }
    
    return "";
}

std::string UnitConverter::formatWithUnit(double value, const std::string& unit, int precision) const {
    if (!isValidUnit(unit)) {
        return "";
    }
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value << " " << m_units.at(unit).symbol;
    return oss.str();
}

void UnitConverter::addUnit(const std::string& unit, 
                          const std::string& quantity,
                          std::function<double(double)> toBaseFunc,
                          std::function<double(double)> fromBaseFunc,
                          const std::string& symbol,
                          const std::string& name) {
    // Add the unit
    UnitInfo unitInfo;
    unitInfo.quantity = quantity;
    unitInfo.toBase = toBaseFunc;
    unitInfo.fromBase = fromBaseFunc;
    unitInfo.symbol = symbol;
    unitInfo.name = name;
    
    m_units[unit] = unitInfo;
    
    // Add the unit to the quantity's unit list
    m_quantities[quantity].units.push_back(unit);
}

void UnitConverter::addQuantity(const std::string& quantity,
                              const std::string& baseUnit,
                              const std::string& name) {
    // Add the quantity
    QuantityInfo quantityInfo;
    quantityInfo.baseUnit = baseUnit;
    quantityInfo.name = name;
    
    m_quantities[quantity] = quantityInfo;
}

// Register length units (meter is the base unit)
void UnitConverter::registerLengthUnits() {
    // Add the length quantity
    addQuantity("length", "m", "Length");
    
    // Add length units
    // SI units
    addUnit("m", "length", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "m", "meter");
    
    addUnit("km", "length", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "km", "kilometer");
    
    addUnit("cm", "length", 
           [](double v) { return v / 100.0; }, // to base
           [](double v) { return v * 100.0; }, // from base
           "cm", "centimeter");
    
    addUnit("mm", "length", 
           [](double v) { return v / 1000.0; }, // to base
           [](double v) { return v * 1000.0; }, // from base
           "mm", "millimeter");
    
    addUnit("um", "length", 
           [](double v) { return v / 1000000.0; }, // to base
           [](double v) { return v * 1000000.0; }, // from base
           "μm", "micrometer");
    
    addUnit("nm", "length", 
           [](double v) { return v / 1000000000.0; }, // to base
           [](double v) { return v * 1000000000.0; }, // from base
           "nm", "nanometer");
    
    // Imperial units
    addUnit("in", "length", 
           [](double v) { return v * 0.0254; }, // to base
           [](double v) { return v / 0.0254; }, // from base
           "in", "inch");
    
    addUnit("ft", "length", 
           [](double v) { return v * 0.3048; }, // to base
           [](double v) { return v / 0.3048; }, // from base
           "ft", "foot");
    
    addUnit("yd", "length", 
           [](double v) { return v * 0.9144; }, // to base
           [](double v) { return v / 0.9144; }, // from base
           "yd", "yard");
    
    addUnit("mi", "length", 
           [](double v) { return v * 1609.344; }, // to base
           [](double v) { return v / 1609.344; }, // from base
           "mi", "mile");
    
    // Nautical units
    addUnit("nmi", "length", 
           [](double v) { return v * 1852.0; }, // to base
           [](double v) { return v / 1852.0; }, // from base
           "nmi", "nautical mile");
    
    // Astronomical units
    addUnit("au", "length", 
           [](double v) { return v * 149597870700.0; }, // to base
           [](double v) { return v / 149597870700.0; }, // from base
           "AU", "astronomical unit");
    
    addUnit("ly", "length", 
           [](double v) { return v * 9460730472580800.0; }, // to base
           [](double v) { return v / 9460730472580800.0; }, // from base
           "ly", "light-year");
    
    addUnit("pc", "length", 
           [](double v) { return v * 30856775814671900.0; }, // to base
           [](double v) { return v / 30856775814671900.0; }, // from base
           "pc", "parsec");
}

// Register mass units (kilogram is the base unit)
void UnitConverter::registerMassUnits() {
    // Add the mass quantity
    addQuantity("mass", "kg", "Mass");
    
    // Add mass units
    // SI units
    addUnit("kg", "mass", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "kg", "kilogram");
    
    addUnit("g", "mass", 
           [](double v) { return v / 1000.0; }, // to base
           [](double v) { return v * 1000.0; }, // from base
           "g", "gram");
    
    addUnit("mg", "mass", 
           [](double v) { return v / 1000000.0; }, // to base
           [](double v) { return v * 1000000.0; }, // from base
           "mg", "milligram");
    
    addUnit("ug", "mass", 
           [](double v) { return v / 1000000000.0; }, // to base
           [](double v) { return v * 1000000000.0; }, // from base
           "μg", "microgram");
    
    addUnit("t", "mass", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "t", "metric ton");
    
    // Imperial units
    addUnit("oz", "mass", 
           [](double v) { return v * 0.028349523125; }, // to base
           [](double v) { return v / 0.028349523125; }, // from base
           "oz", "ounce");
    
    addUnit("lb", "mass", 
           [](double v) { return v * 0.45359237; }, // to base
           [](double v) { return v / 0.45359237; }, // from base
           "lb", "pound");
    
    addUnit("st", "mass", 
           [](double v) { return v * 6.35029318; }, // to base
           [](double v) { return v / 6.35029318; }, // from base
           "st", "stone");
    
    addUnit("ton", "mass", 
           [](double v) { return v * 907.18474; }, // to base
           [](double v) { return v / 907.18474; }, // from base
           "ton", "short ton");
    
    addUnit("lton", "mass", 
           [](double v) { return v * 1016.0469088; }, // to base
           [](double v) { return v / 1016.0469088; }, // from base
           "lton", "long ton");
    
    // Atomic units
    addUnit("u", "mass", 
           [](double v) { return v * 1.66053906660e-27; }, // to base
           [](double v) { return v / 1.66053906660e-27; }, // from base
           "u", "atomic mass unit");
}

// Register time units (second is the base unit)
void UnitConverter::registerTimeUnits() {
    // Add the time quantity
    addQuantity("time", "s", "Time");
    
    // Add time units
    // SI units
    addUnit("s", "time", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "s", "second");
    
    addUnit("ms", "time", 
           [](double v) { return v / 1000.0; }, // to base
           [](double v) { return v * 1000.0; }, // from base
           "ms", "millisecond");
    
    addUnit("us", "time", 
           [](double v) { return v / 1000000.0; }, // to base
           [](double v) { return v * 1000000.0; }, // from base
           "μs", "microsecond");
    
    addUnit("ns", "time", 
           [](double v) { return v / 1000000000.0; }, // to base
           [](double v) { return v * 1000000000.0; }, // from base
           "ns", "nanosecond");
    
    // Common time units
    addUnit("min", "time", 
           [](double v) { return v * 60.0; }, // to base
           [](double v) { return v / 60.0; }, // from base
           "min", "minute");
    
    addUnit("h", "time", 
           [](double v) { return v * 3600.0; }, // to base
           [](double v) { return v / 3600.0; }, // from base
           "h", "hour");
    
    addUnit("day", "time", 
           [](double v) { return v * 86400.0; }, // to base
           [](double v) { return v / 86400.0; }, // from base
           "d", "day");
    
    addUnit("week", "time", 
           [](double v) { return v * 604800.0; }, // to base
           [](double v) { return v / 604800.0; }, // from base
           "wk", "week");
    
    addUnit("month", "time", 
           [](double v) { return v * 2629746.0; }, // to base (average month)
           [](double v) { return v / 2629746.0; }, // from base
           "mo", "month");
    
    addUnit("year", "time", 
           [](double v) { return v * 31556952.0; }, // to base (average year)
           [](double v) { return v / 31556952.0; }, // from base
           "yr", "year");
}

// Register temperature units (kelvin is the base unit)
void UnitConverter::registerTemperatureUnits() {
    // Add the temperature quantity
    addQuantity("temperature", "K", "Temperature");
    
    // Add temperature units
    // SI unit
    addUnit("K", "temperature", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "K", "kelvin");
    
    // Celsius
    addUnit("C", "temperature", 
           [](double v) { return v + 273.15; }, // to base
           [](double v) { return v - 273.15; }, // from base
           "°C", "Celsius");
    
    // Fahrenheit
    addUnit("F", "temperature", 
           [](double v) { return (v - 32.0) * 5.0 / 9.0 + 273.15; }, // to base
           [](double v) { return (v - 273.15) * 9.0 / 5.0 + 32.0; }, // from base
           "°F", "Fahrenheit");
    
    // Rankine
    addUnit("R", "temperature", 
           [](double v) { return v * 5.0 / 9.0; }, // to base
           [](double v) { return v * 9.0 / 5.0; }, // from base
           "°R", "Rankine");
}

// Register area units (square meter is the base unit)
void UnitConverter::registerAreaUnits() {
    // Add the area quantity
    addQuantity("area", "m2", "Area");
    
    // Add area units
    // SI units
    addUnit("m2", "area", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "m²", "square meter");
    
    addUnit("cm2", "area", 
           [](double v) { return v / 10000.0; }, // to base
           [](double v) { return v * 10000.0; }, // from base
           "cm²", "square centimeter");
    
    addUnit("mm2", "area", 
           [](double v) { return v / 1000000.0; }, // to base
           [](double v) { return v * 1000000.0; }, // from base
           "mm²", "square millimeter");
    
    addUnit("km2", "area", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "km²", "square kilometer");
    
    // Imperial units
    addUnit("in2", "area", 
           [](double v) { return v * 0.00064516; }, // to base
           [](double v) { return v / 0.00064516; }, // from base
           "in²", "square inch");
    
    addUnit("ft2", "area", 
           [](double v) { return v * 0.09290304; }, // to base
           [](double v) { return v / 0.09290304; }, // from base
           "ft²", "square foot");
    
    addUnit("yd2", "area", 
           [](double v) { return v * 0.83612736; }, // to base
           [](double v) { return v / 0.83612736; }, // from base
           "yd²", "square yard");
    
    addUnit("mi2", "area", 
           [](double v) { return v * 2589988.110336; }, // to base
           [](double v) { return v / 2589988.110336; }, // from base
           "mi²", "square mile");
    
    // Other units
    addUnit("ha", "area", 
           [](double v) { return v * 10000.0; }, // to base
           [](double v) { return v / 10000.0; }, // from base
           "ha", "hectare");
    
    addUnit("acre", "area", 
           [](double v) { return v * 4046.8564224; }, // to base
           [](double v) { return v / 4046.8564224; }, // from base
           "ac", "acre");
}

// Register volume units (cubic meter is the base unit)
void UnitConverter::registerVolumeUnits() {
    // Add the volume quantity
    addQuantity("volume", "m3", "Volume");
    
    // Add volume units
    // SI units
    addUnit("m3", "volume", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "m³", "cubic meter");
    
    addUnit("cm3", "volume", 
           [](double v) { return v / 1000000.0; }, // to base
           [](double v) { return v * 1000000.0; }, // from base
           "cm³", "cubic centimeter");
    
    addUnit("mm3", "volume", 
           [](double v) { return v / 1000000000.0; }, // to base
           [](double v) { return v * 1000000000.0; }, // from base
           "mm³", "cubic millimeter");
    
    addUnit("L", "volume", 
           [](double v) { return v / 1000.0; }, // to base
           [](double v) { return v * 1000.0; }, // from base
           "L", "liter");
    
    addUnit("mL", "volume", 
           [](double v) { return v / 1000000.0; }, // to base
           [](double v) { return v * 1000000.0; }, // from base
           "mL", "milliliter");
    
    // Imperial units
    addUnit("in3", "volume", 
           [](double v) { return v * 0.000016387064; }, // to base
           [](double v) { return v / 0.000016387064; }, // from base
           "in³", "cubic inch");
    
    addUnit("ft3", "volume", 
           [](double v) { return v * 0.028316846592; }, // to base
           [](double v) { return v / 0.028316846592; }, // from base
           "ft³", "cubic foot");
    
    addUnit("yd3", "volume", 
           [](double v) { return v * 0.764554857984; }, // to base
           [](double v) { return v / 0.764554857984; }, // from base
           "yd³", "cubic yard");
    
    // US fluid units
    addUnit("floz", "volume", 
           [](double v) { return v * 0.0000295735295625; }, // to base
           [](double v) { return v / 0.0000295735295625; }, // from base
           "fl oz", "fluid ounce (US)");
    
    addUnit("cup", "volume", 
           [](double v) { return v * 0.0002365882365; }, // to base
           [](double v) { return v / 0.0002365882365; }, // from base
           "cup", "cup (US)");
    
    addUnit("pt", "volume", 
           [](double v) { return v * 0.000473176473; }, // to base
           [](double v) { return v / 0.000473176473; }, // from base
           "pt", "pint (US)");
    
    addUnit("qt", "volume", 
           [](double v) { return v * 0.000946352946; }, // to base
           [](double v) { return v / 0.000946352946; }, // from base
           "qt", "quart (US)");
    
    addUnit("gal", "volume", 
           [](double v) { return v * 0.003785411784; }, // to base
           [](double v) { return v / 0.003785411784; }, // from base
           "gal", "gallon (US)");
    
    // UK fluid units
    addUnit("ukfloz", "volume", 
           [](double v) { return v * 0.0000284130625; }, // to base
           [](double v) { return v / 0.0000284130625; }, // from base
           "fl oz (UK)", "fluid ounce (UK)");
    
    addUnit("ukpt", "volume", 
           [](double v) { return v * 0.00056826125; }, // to base
           [](double v) { return v / 0.00056826125; }, // from base
           "pt (UK)", "pint (UK)");
    
    addUnit("ukqt", "volume", 
           [](double v) { return v * 0.0011365225; }, // to base
           [](double v) { return v / 0.0011365225; }, // from base
           "qt (UK)", "quart (UK)");
    
    addUnit("ukgal", "volume", 
           [](double v) { return v * 0.00454609; }, // to base
           [](double v) { return v / 0.00454609; }, // from base
           "gal (UK)", "gallon (UK)");
}

// Register velocity units (meter per second is the base unit)
void UnitConverter::registerVelocityUnits() {
    // Add the velocity quantity
    addQuantity("velocity", "m/s", "Velocity");
    
    // Add velocity units
    // SI units
    addUnit("m/s", "velocity", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "m/s", "meter per second");
    
    addUnit("km/h", "velocity", 
           [](double v) { return v * 0.277777778; }, // to base
           [](double v) { return v / 0.277777778; }, // from base
           "km/h", "kilometer per hour");
    
    // Imperial units
    addUnit("ft/s", "velocity", 
           [](double v) { return v * 0.3048; }, // to base
           [](double v) { return v / 0.3048; }, // from base
           "ft/s", "foot per second");
    
    addUnit("mph", "velocity", 
           [](double v) { return v * 0.44704; }, // to base
           [](double v) { return v / 0.44704; }, // from base
           "mph", "mile per hour");
    
    // Nautical units
    addUnit("kn", "velocity", 
           [](double v) { return v * 0.514444444; }, // to base
           [](double v) { return v / 0.514444444; }, // from base
           "kn", "knot");
    
    // Other units
    addUnit("mach", "velocity", 
           [](double v) { return v * 340.29; }, // to base (at sea level, 15°C)
           [](double v) { return v / 340.29; }, // from base
           "Mach", "Mach number");
    
    addUnit("c", "velocity", 
           [](double v) { return v * 299792458.0; }, // to base
           [](double v) { return v / 299792458.0; }, // from base
           "c", "speed of light");
}

// Register acceleration units (meter per second squared is the base unit)
void UnitConverter::registerAccelerationUnits() {
    // Add the acceleration quantity
    addQuantity("acceleration", "m/s2", "Acceleration");
    
    // Add acceleration units
    // SI units
    addUnit("m/s2", "acceleration", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "m/s²", "meter per second squared");
    
    // Other units
    addUnit("g", "acceleration", 
           [](double v) { return v * 9.80665; }, // to base
           [](double v) { return v / 9.80665; }, // from base
           "g", "standard gravity");
    
    addUnit("ft/s2", "acceleration", 
           [](double v) { return v * 0.3048; }, // to base
           [](double v) { return v / 0.3048; }, // from base
           "ft/s²", "foot per second squared");
}

// Register force units (newton is the base unit)
void UnitConverter::registerForceUnits() {
    // Add the force quantity
    addQuantity("force", "N", "Force");
    
    // Add force units
    // SI units
    addUnit("N", "force", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "N", "newton");
    
    addUnit("kN", "force", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "kN", "kilonewton");
    
    addUnit("MN", "force", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "MN", "meganewton");
    
    // Imperial units
    addUnit("lbf", "force", 
           [](double v) { return v * 4.4482216153; }, // to base
           [](double v) { return v / 4.4482216153; }, // from base
           "lbf", "pound-force");
    
    addUnit("kgf", "force", 
           [](double v) { return v * 9.80665; }, // to base
           [](double v) { return v / 9.80665; }, // from base
           "kgf", "kilogram-force");
    
    addUnit("dyn", "force", 
           [](double v) { return v * 0.00001; }, // to base
           [](double v) { return v / 0.00001; }, // from base
           "dyn", "dyne");
}

// Register energy units (joule is the base unit)
void UnitConverter::registerEnergyUnits() {
    // Add the energy quantity
    addQuantity("energy", "J", "Energy");
    
    // Add energy units
    // SI units
    addUnit("J", "energy", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "J", "joule");
    
    addUnit("kJ", "energy", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "kJ", "kilojoule");
    
    addUnit("MJ", "energy", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "MJ", "megajoule");
    
    // Other energy units
    addUnit("cal", "energy", 
           [](double v) { return v * 4.184; }, // to base
           [](double v) { return v / 4.184; }, // from base
           "cal", "calorie");
    
    addUnit("kcal", "energy", 
           [](double v) { return v * 4184.0; }, // to base
           [](double v) { return v / 4184.0; }, // from base
           "kcal", "kilocalorie");
    
    addUnit("Wh", "energy", 
           [](double v) { return v * 3600.0; }, // to base
           [](double v) { return v / 3600.0; }, // from base
           "Wh", "watt-hour");
    
    addUnit("kWh", "energy", 
           [](double v) { return v * 3600000.0; }, // to base
           [](double v) { return v / 3600000.0; }, // from base
           "kWh", "kilowatt-hour");
    
    addUnit("eV", "energy", 
           [](double v) { return v * 1.602176634e-19; }, // to base
           [](double v) { return v / 1.602176634e-19; }, // from base
           "eV", "electronvolt");
    
    addUnit("BTU", "energy", 
           [](double v) { return v * 1055.05585262; }, // to base
           [](double v) { return v / 1055.05585262; }, // from base
           "BTU", "British thermal unit");
}

// Register power units (watt is the base unit)
void UnitConverter::registerPowerUnits() {
    // Add the power quantity
    addQuantity("power", "W", "Power");
    
    // Add power units
    // SI units
    addUnit("W", "power", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "W", "watt");
    
    addUnit("kW", "power", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "kW", "kilowatt");
    
    addUnit("MW", "power", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "MW", "megawatt");
    
    addUnit("GW", "power", 
           [](double v) { return v * 1000000000.0; }, // to base
           [](double v) { return v / 1000000000.0; }, // from base
           "GW", "gigawatt");
    
    // Other power units
    addUnit("hp", "power", 
           [](double v) { return v * 745.699872; }, // to base
           [](double v) { return v / 745.699872; }, // from base
           "hp", "horsepower");
    
    addUnit("BTU/h", "power", 
           [](double v) { return v * 0.29307107; }, // to base
           [](double v) { return v / 0.29307107; }, // from base
           "BTU/h", "BTU per hour");
}

// Register pressure units (pascal is the base unit)
void UnitConverter::registerPressureUnits() {
    // Add the pressure quantity
    addQuantity("pressure", "Pa", "Pressure");
    
    // Add pressure units
    // SI units
    addUnit("Pa", "pressure", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "Pa", "pascal");
    
    addUnit("kPa", "pressure", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "kPa", "kilopascal");
    
    addUnit("MPa", "pressure", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "MPa", "megapascal");
    
    // Other pressure units
    addUnit("bar", "pressure", 
           [](double v) { return v * 100000.0; }, // to base
           [](double v) { return v / 100000.0; }, // from base
           "bar", "bar");
    
    addUnit("mbar", "pressure", 
           [](double v) { return v * 100.0; }, // to base
           [](double v) { return v / 100.0; }, // from base
           "mbar", "millibar");
    
    addUnit("atm", "pressure", 
           [](double v) { return v * 101325.0; }, // to base
           [](double v) { return v / 101325.0; }, // from base
           "atm", "atmosphere");
    
    addUnit("torr", "pressure", 
           [](double v) { return v * 133.322368421; }, // to base
           [](double v) { return v / 133.322368421; }, // from base
           "Torr", "torr");
    
    addUnit("mmHg", "pressure", 
           [](double v) { return v * 133.322368421; }, // to base
           [](double v) { return v / 133.322368421; }, // from base
           "mmHg", "millimeter of mercury");
    
    addUnit("inHg", "pressure", 
           [](double v) { return v * 3386.38815789; }, // to base
           [](double v) { return v / 3386.38815789; }, // from base
           "inHg", "inch of mercury");
    
    addUnit("psi", "pressure", 
           [](double v) { return v * 6894.75729317; }, // to base
           [](double v) { return v / 6894.75729317; }, // from base
           "psi", "pound per square inch");
}

// Register electric current units (ampere is the base unit)
void UnitConverter::registerElectricCurrentUnits() {
    // Add the electric current quantity
    addQuantity("electric_current", "A", "Electric Current");
    
    // Add electric current units
    // SI units
    addUnit("A", "electric_current", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "A", "ampere");
    
    addUnit("mA", "electric_current", 
           [](double v) { return v / 1000.0; }, // to base
           [](double v) { return v * 1000.0; }, // from base
           "mA", "milliampere");
    
    addUnit("uA", "electric_current", 
           [](double v) { return v / 1000000.0; }, // to base
           [](double v) { return v * 1000000.0; }, // from base
           "μA", "microampere");
    
    addUnit("kA", "electric_current", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "kA", "kiloampere");
}

// Register electric voltage units (volt is the base unit)
void UnitConverter::registerElectricVoltageUnits() {
    // Add the electric voltage quantity
    addQuantity("electric_voltage", "V", "Electric Voltage");
    
    // Add electric voltage units
    // SI units
    addUnit("V", "electric_voltage", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "V", "volt");
    
    addUnit("mV", "electric_voltage", 
           [](double v) { return v / 1000.0; }, // to base
           [](double v) { return v * 1000.0; }, // from base
           "mV", "millivolt");
    
    addUnit("uV", "electric_voltage", 
           [](double v) { return v / 1000000.0; }, // to base
           [](double v) { return v * 1000000.0; }, // from base
           "μV", "microvolt");
    
    addUnit("kV", "electric_voltage", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "kV", "kilovolt");
    
    addUnit("MV", "electric_voltage", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "MV", "megavolt");
}

// Register electric resistance units (ohm is the base unit)
void UnitConverter::registerElectricResistanceUnits() {
    // Add the electric resistance quantity
    addQuantity("electric_resistance", "ohm", "Electric Resistance");
    
    // Add electric resistance units
    // SI units
    addUnit("ohm", "electric_resistance", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "Ω", "ohm");
    
    addUnit("kohm", "electric_resistance", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "kΩ", "kiloohm");
    
    addUnit("Mohm", "electric_resistance", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "MΩ", "megaohm");
    
    addUnit("mOhm", "electric_resistance", 
           [](double v) { return v / 1000.0; }, // to base
           [](double v) { return v * 1000.0; }, // from base
           "mΩ", "milliohm");
}

// Register frequency units (hertz is the base unit)
void UnitConverter::registerFrequencyUnits() {
    // Add the frequency quantity
    addQuantity("frequency", "Hz", "Frequency");
    
    // Add frequency units
    // SI units
    addUnit("Hz", "frequency", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "Hz", "hertz");
    
    addUnit("kHz", "frequency", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "kHz", "kilohertz");
    
    addUnit("MHz", "frequency", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "MHz", "megahertz");
    
    addUnit("GHz", "frequency", 
           [](double v) { return v * 1000000000.0; }, // to base
           [](double v) { return v / 1000000000.0; }, // from base
           "GHz", "gigahertz");
    
    // Other frequency units
    addUnit("rpm", "frequency", 
           [](double v) { return v / 60.0; }, // to base
           [](double v) { return v * 60.0; }, // from base
           "rpm", "revolutions per minute");
    
    addUnit("rad/s", "frequency", 
           [](double v) { return v / (2.0 * M_PI); }, // to base
           [](double v) { return v * (2.0 * M_PI); }, // from base
           "rad/s", "radians per second");
}

// Register data units (byte is the base unit)
void UnitConverter::registerDataUnits() {
    // Add the data quantity
    addQuantity("data", "B", "Data");
    
    // Add data units
    // Base units
    addUnit("B", "data", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "B", "byte");
    
    addUnit("bit", "data", 
           [](double v) { return v / 8.0; }, // to base
           [](double v) { return v * 8.0; }, // from base
           "bit", "bit");
    
    // Binary prefixes
    addUnit("KiB", "data", 
           [](double v) { return v * 1024.0; }, // to base
           [](double v) { return v / 1024.0; }, // from base
           "KiB", "kibibyte");
    
    addUnit("MiB", "data", 
           [](double v) { return v * 1048576.0; }, // to base
           [](double v) { return v / 1048576.0; }, // from base
           "MiB", "mebibyte");
    
    addUnit("GiB", "data", 
           [](double v) { return v * 1073741824.0; }, // to base
           [](double v) { return v / 1073741824.0; }, // from base
           "GiB", "gibibyte");
    
    addUnit("TiB", "data", 
           [](double v) { return v * 1099511627776.0; }, // to base
           [](double v) { return v / 1099511627776.0; }, // from base
           "TiB", "tebibyte");
    
    // Decimal prefixes
    addUnit("KB", "data", 
           [](double v) { return v * 1000.0; }, // to base
           [](double v) { return v / 1000.0; }, // from base
           "KB", "kilobyte");
    
    addUnit("MB", "data", 
           [](double v) { return v * 1000000.0; }, // to base
           [](double v) { return v / 1000000.0; }, // from base
           "MB", "megabyte");
    
    addUnit("GB", "data", 
           [](double v) { return v * 1000000000.0; }, // to base
           [](double v) { return v / 1000000000.0; }, // from base
           "GB", "gigabyte");
    
    addUnit("TB", "data", 
           [](double v) { return v * 1000000000000.0; }, // to base
           [](double v) { return v / 1000000000000.0; }, // from base
           "TB", "terabyte");
    
    // Bit versions
    addUnit("Kbit", "data", 
           [](double v) { return v * 1000.0 / 8.0; }, // to base
           [](double v) { return v * 8.0 / 1000.0; }, // from base
           "Kbit", "kilobit");
    
    addUnit("Mbit", "data", 
           [](double v) { return v * 1000000.0 / 8.0; }, // to base
           [](double v) { return v * 8.0 / 1000000.0; }, // from base
           "Mbit", "megabit");
    
    addUnit("Gbit", "data", 
           [](double v) { return v * 1000000000.0 / 8.0; }, // to base
           [](double v) { return v * 8.0 / 1000000000.0; }, // from base
           "Gbit", "gigabit");
    
    addUnit("Tbit", "data", 
           [](double v) { return v * 1000000000000.0 / 8.0; }, // to base
           [](double v) { return v * 8.0 / 1000000000000.0; }, // from base
           "Tbit", "terabit");
}

// Register angle units (radian is the base unit)
void UnitConverter::registerAngleUnits() {
    // Add the angle quantity
    addQuantity("angle", "rad", "Angle");
    
    // Add angle units
    // SI unit
    addUnit("rad", "angle", 
           [](double v) { return v; }, // to base (identity)
           [](double v) { return v; }, // from base (identity)
           "rad", "radian");
    
    // Other angle units
    addUnit("deg", "angle", 
           [](double v) { return v * M_PI / 180.0; }, // to base
           [](double v) { return v * 180.0 / M_PI; }, // from base
           "°", "degree");
    
    addUnit("grad", "angle", 
           [](double v) { return v * M_PI / 200.0; }, // to base
           [](double v) { return v * 200.0 / M_PI; }, // from base
           "grad", "gradian");
    
    addUnit("arcmin", "angle", 
           [](double v) { return v * M_PI / 10800.0; }, // to base
           [](double v) { return v * 10800.0 / M_PI; }, // from base
           "′", "arcminute");
    
    addUnit("arcsec", "angle", 
           [](double v) { return v * M_PI / 648000.0; }, // to base
           [](double v) { return v * 648000.0 / M_PI; }, // from base
           "″", "arcsecond");
    
    addUnit("turn", "angle", 
           [](double v) { return v * 2.0 * M_PI; }, // to base
           [](double v) { return v / (2.0 * M_PI); }, // from base
           "turn", "turn");
}

} // namespace rebelcalc
