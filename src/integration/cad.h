#ifndef REBELCALC_INTEGRATION_CAD_H
#define REBELCALC_INTEGRATION_CAD_H

#include "api.h"
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace RebelCalc {
namespace Integration {

/**
 * @brief Integration interface for RebelCAD
 */
class CADIntegration : public IntegrationInterface {
public:
    /**
     * @brief Constructor
     */
    CADIntegration();

    /**
     * @brief Destructor
     */
    ~CADIntegration() override;

    /**
     * @brief Get the component type
     * 
     * @return Component Type of the component (REBELCAD)
     */
    Component getComponentType() const override;

    /**
     * @brief Get the component name
     * 
     * @return std::string Name of the component ("RebelCAD")
     */
    std::string getComponentName() const override;

    /**
     * @brief Check if RebelCAD is available
     * 
     * @return true if RebelCAD is available, false otherwise
     */
    bool isAvailable() const override;

    /**
     * @brief Initialize the integration with RebelCAD
     * 
     * @return true if initialization was successful, false otherwise
     */
    bool initialize() override;

    /**
     * @brief Shutdown the integration with RebelCAD
     */
    void shutdown() override;

    /**
     * @brief Execute a command on RebelCAD
     * 
     * @param command Command to execute
     * @param data Data exchange object for input/output
     * @return true if the command was executed successfully, false otherwise
     */
    bool executeCommand(const std::string& command, DataExchange& data) override;

    /**
     * @brief Get the available commands for RebelCAD
     * 
     * @return std::vector<std::string> Available commands
     */
    std::vector<std::string> getAvailableCommands() const override;

    /**
     * @brief Register a callback for events from RebelCAD
     * 
     * @param eventName Name of the event
     * @param callback Callback function
     * @return true if the callback was registered successfully, false otherwise
     */
    bool registerCallback(const std::string& eventName, std::function<void(const DataExchange&)> callback) override;

    /**
     * @brief Unregister a callback
     * 
     * @param eventName Name of the event
     * @return true if the callback was unregistered successfully, false otherwise
     */
    bool unregisterCallback(const std::string& eventName) override;

    // RebelCAD-specific methods

    /**
     * @brief Import a CAD model from RebelCAD
     * 
     * @param modelName Name of the model to import
     * @return std::vector<Engineering::CAD::Triangle> Triangles representing the model
     */
    std::vector<Engineering::CAD::Triangle> importModel(const std::string& modelName);

    /**
     * @brief Export a CAD model to RebelCAD
     * 
     * @param modelName Name to give the model
     * @param triangles Triangles representing the model
     * @return true if the export was successful, false otherwise
     */
    bool exportModel(const std::string& modelName, const std::vector<Engineering::CAD::Triangle>& triangles);

    /**
     * @brief Get the available models in RebelCAD
     * 
     * @return std::vector<std::string> Names of available models
     */
    std::vector<std::string> getAvailableModels() const;

    /**
     * @brief Perform a boolean operation on two models
     * 
     * @param operation Operation to perform ("union", "intersection", "difference")
     * @param model1Name Name of the first model
     * @param model2Name Name of the second model
     * @param resultModelName Name to give the result model
     * @return true if the operation was successful, false otherwise
     */
    bool performBooleanOperation(const std::string& operation, 
                                const std::string& model1Name, 
                                const std::string& model2Name, 
                                const std::string& resultModelName);

    /**
     * @brief Perform a measurement on a model
     * 
     * @param modelName Name of the model
     * @param measurementType Type of measurement ("volume", "surface_area", "bounding_box")
     * @return DataVariant Result of the measurement
     */
    DataVariant performMeasurement(const std::string& modelName, const std::string& measurementType);

private:
    bool m_initialized;
    bool m_available;
    std::map<std::string, std::function<void(const DataExchange&)>> m_callbacks;
    std::vector<std::string> m_availableCommands;
    std::vector<std::string> m_availableModels;
};

} // namespace Integration
} // namespace RebelCalc

#endif // REBELCALC_INTEGRATION_CAD_H
