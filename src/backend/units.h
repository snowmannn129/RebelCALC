#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <optional>
#include <functional>

namespace rebelcalc {

/**
 * @class UnitConverter
 * @brief Handles unit conversions between different measurement systems
 * 
 * The UnitConverter class provides functionality for converting between different
 * units of measurement across various physical quantities such as length, mass,
 * time, temperature, etc. It supports both SI and imperial units, as well as
 * common derived units used in engineering and scientific calculations.
 */
class UnitConverter {
public:
    /**
     * @brief Constructor
     */
    UnitConverter();
    
    /**
     * @brief Destructor
     */
    ~UnitConverter();
    
    /**
     * @brief Initialize the unit converter
     * @return true if initialization was successful, false otherwise
     */
    bool initialize();
    
    /**
     * @brief Shutdown the unit converter
     */
    void shutdown();
    
    /**
     * @brief Convert a value from one unit to another
     * @param value The value to convert
     * @param fromUnit The source unit
     * @param toUnit The target unit
     * @return The converted value, or nullopt if conversion failed
     */
    std::optional<double> convert(double value, const std::string& fromUnit, const std::string& toUnit);
    
    /**
     * @brief Get all available units for a specific quantity
     * @param quantity The physical quantity (e.g., "length", "mass", "time")
     * @return A vector of unit names, or empty vector if quantity is not supported
     */
    std::vector<std::string> getUnitsForQuantity(const std::string& quantity) const;
    
    /**
     * @brief Get all supported physical quantities
     * @return A vector of physical quantity names
     */
    std::vector<std::string> getAllQuantities() const;
    
    /**
     * @brief Get the base unit for a specific quantity
     * @param quantity The physical quantity
     * @return The base unit name, or empty string if quantity is not supported
     */
    std::string getBaseUnit(const std::string& quantity) const;
    
    /**
     * @brief Get the quantity for a specific unit
     * @param unit The unit name
     * @return The physical quantity name, or empty string if unit is not supported
     */
    std::string getQuantityForUnit(const std::string& unit) const;
    
    /**
     * @brief Check if a unit is valid
     * @param unit The unit name
     * @return true if the unit is valid, false otherwise
     */
    bool isValidUnit(const std::string& unit) const;
    
    /**
     * @brief Check if a quantity is valid
     * @param quantity The physical quantity name
     * @return true if the quantity is valid, false otherwise
     */
    bool isValidQuantity(const std::string& quantity) const;
    
    /**
     * @brief Register a custom unit
     * @param unit The unit name
     * @param quantity The physical quantity
     * @param toBaseFunc Function to convert from this unit to the base unit
     * @param fromBaseFunc Function to convert from the base unit to this unit
     * @return true if registration was successful, false otherwise
     */
    bool registerUnit(const std::string& unit, 
                     const std::string& quantity,
                     std::function<double(double)> toBaseFunc,
                     std::function<double(double)> fromBaseFunc);
    
    /**
     * @brief Parse a unit expression (e.g., "m/s", "kg*m^2")
     * @param expression The unit expression
     * @return The simplified unit, or empty string if parsing failed
     */
    std::string parseUnitExpression(const std::string& expression) const;
    
    /**
     * @brief Format a value with its unit
     * @param value The value
     * @param unit The unit
     * @param precision The number of decimal places
     * @return The formatted string
     */
    std::string formatWithUnit(double value, const std::string& unit, int precision = 2) const;

private:
    /**
     * @brief Unit information structure
     */
    struct UnitInfo {
        std::string quantity;                  // Physical quantity (e.g., "length", "mass")
        std::function<double(double)> toBase;  // Convert from this unit to base unit
        std::function<double(double)> fromBase; // Convert from base unit to this unit
        std::string symbol;                    // Unit symbol (e.g., "m", "kg")
        std::string name;                      // Full name (e.g., "meter", "kilogram")
    };
    
    /**
     * @brief Quantity information structure
     */
    struct QuantityInfo {
        std::string baseUnit;                  // Base unit for this quantity
        std::vector<std::string> units;        // All units for this quantity
        std::string name;                      // Full name (e.g., "length", "mass")
    };
    
    // Maps for unit and quantity information
    std::unordered_map<std::string, UnitInfo> m_units;
    std::unordered_map<std::string, QuantityInfo> m_quantities;
    
    // Helper methods
    void registerLengthUnits();
    void registerMassUnits();
    void registerTimeUnits();
    void registerTemperatureUnits();
    void registerAreaUnits();
    void registerVolumeUnits();
    void registerVelocityUnits();
    void registerAccelerationUnits();
    void registerForceUnits();
    void registerEnergyUnits();
    void registerPowerUnits();
    void registerPressureUnits();
    void registerElectricCurrentUnits();
    void registerElectricVoltageUnits();
    void registerElectricResistanceUnits();
    void registerFrequencyUnits();
    void registerDataUnits();
    void registerAngleUnits();
    
    // Register a unit with its information
    void addUnit(const std::string& unit, 
                const std::string& quantity,
                std::function<double(double)> toBaseFunc,
                std::function<double(double)> fromBaseFunc,
                const std::string& symbol,
                const std::string& name);
    
    // Register a quantity with its information
    void addQuantity(const std::string& quantity,
                    const std::string& baseUnit,
                    const std::string& name);
};

} // namespace rebelcalc
