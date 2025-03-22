#ifndef REBELCALC_INTEGRATION_API_H
#define REBELCALC_INTEGRATION_API_H

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <memory>
#include <variant>
#include "../backend/matrix.h"
#include "../engineering/cad.h"

namespace RebelCalc {
namespace Integration {

/**
 * @brief Enum for the different RebelSUITE components
 */
enum class Component {
    REBELCAD,
    REBELENGINE,
    REBELFLOW,
    REBELCODE,
    REBELDESK,
    REBELSCRIBE
};

/**
 * @brief Data type for exchanging data between components
 */
using DataVariant = std::variant<
    int,
    double,
    bool,
    std::string,
    std::vector<double>,
    std::vector<std::string>,
    std::vector<Engineering::CAD::Triangle>,
    rebelcalc::Matrix,
    Engineering::CAD::Point3D,
    Engineering::CAD::LineSegment,
    Engineering::CAD::Plane,
    Engineering::CAD::Triangle,
    Engineering::CAD::Circle,
    Engineering::CAD::Polygon
>;

/**
 * @brief Class for handling data exchange between components
 */
class DataExchange {
public:
    /**
     * @brief Constructor
     */
    DataExchange() = default;

    /**
     * @brief Set data with a key
     * 
     * @param key Key for the data
     * @param data Data to store
     */
    void setData(const std::string& key, const DataVariant& data);

    /**
     * @brief Get data by key
     * 
     * @param key Key for the data
     * @return DataVariant Data stored with the key
     * @throws std::out_of_range if the key doesn't exist
     */
    DataVariant getData(const std::string& key) const;

    /**
     * @brief Check if a key exists
     * 
     * @param key Key to check
     * @return true if the key exists, false otherwise
     */
    bool hasKey(const std::string& key) const;

    /**
     * @brief Get all keys
     * 
     * @return std::vector<std::string> All keys
     */
    std::vector<std::string> getKeys() const;

    /**
     * @brief Clear all data
     */
    void clear();

private:
    std::map<std::string, DataVariant> m_data;
};

/**
 * @brief Base interface for component integration
 */
class IntegrationInterface {
public:
    /**
     * @brief Virtual destructor
     */
    virtual ~IntegrationInterface() = default;

    /**
     * @brief Get the component type
     * 
     * @return Component Type of the component
     */
    virtual Component getComponentType() const = 0;

    /**
     * @brief Get the component name
     * 
     * @return std::string Name of the component
     */
    virtual std::string getComponentName() const = 0;

    /**
     * @brief Check if the component is available
     * 
     * @return true if the component is available, false otherwise
     */
    virtual bool isAvailable() const = 0;

    /**
     * @brief Initialize the integration
     * 
     * @return true if initialization was successful, false otherwise
     */
    virtual bool initialize() = 0;

    /**
     * @brief Shutdown the integration
     */
    virtual void shutdown() = 0;

    /**
     * @brief Execute a command on the component
     * 
     * @param command Command to execute
     * @param data Data exchange object for input/output
     * @return true if the command was executed successfully, false otherwise
     */
    virtual bool executeCommand(const std::string& command, DataExchange& data) = 0;

    /**
     * @brief Get the available commands
     * 
     * @return std::vector<std::string> Available commands
     */
    virtual std::vector<std::string> getAvailableCommands() const = 0;

    /**
     * @brief Register a callback for events from the component
     * 
     * @param eventName Name of the event
     * @param callback Callback function
     * @return true if the callback was registered successfully, false otherwise
     */
    virtual bool registerCallback(const std::string& eventName, std::function<void(const DataExchange&)> callback) = 0;

    /**
     * @brief Unregister a callback
     * 
     * @param eventName Name of the event
     * @return true if the callback was unregistered successfully, false otherwise
     */
    virtual bool unregisterCallback(const std::string& eventName) = 0;
};

/**
 * @brief Factory for creating integration interfaces
 */
class IntegrationFactory {
public:
    /**
     * @brief Get the singleton instance
     * 
     * @return IntegrationFactory& Singleton instance
     */
    static IntegrationFactory& getInstance();

    /**
     * @brief Create an integration interface for a component
     * 
     * @param component Component to create an interface for
     * @return std::shared_ptr<IntegrationInterface> Interface for the component
     */
    std::shared_ptr<IntegrationInterface> createInterface(Component component);

    /**
     * @brief Get all available components
     * 
     * @return std::vector<Component> Available components
     */
    std::vector<Component> getAvailableComponents() const;

private:
    /**
     * @brief Private constructor for singleton
     */
    IntegrationFactory() = default;

    /**
     * @brief Private copy constructor
     */
    IntegrationFactory(const IntegrationFactory&) = delete;

    /**
     * @brief Private assignment operator
     */
    IntegrationFactory& operator=(const IntegrationFactory&) = delete;
};

} // namespace Integration
} // namespace RebelCalc

#endif // REBELCALC_INTEGRATION_API_H
