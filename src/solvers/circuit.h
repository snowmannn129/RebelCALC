#ifndef REBELCALC_SOLVERS_CIRCUIT_H
#define REBELCALC_SOLVERS_CIRCUIT_H

#include <vector>
#include <string>
#include <map>
#include <complex>
#include <functional>
#include <memory>
#include <stdexcept>
#include "../backend/matrix.h"
#include "../engineering/electrical.h"

namespace RebelCalc {
namespace Solvers {

/**
 * @brief Enum for analysis types in circuit simulation
 */
enum class CircuitAnalysisType {
    DC,             // DC analysis
    AC,             // AC analysis
    TRANSIENT,      // Transient analysis
    FREQUENCY,      // Frequency response analysis
    NOISE,          // Noise analysis
    DISTORTION,     // Distortion analysis
    MONTE_CARLO,    // Monte Carlo analysis
    WORST_CASE,     // Worst-case analysis
    SENSITIVITY     // Sensitivity analysis
};

/**
 * @brief Enum for component types in circuit simulation
 */
enum class ComponentType {
    RESISTOR,       // Resistor
    CAPACITOR,      // Capacitor
    INDUCTOR,       // Inductor
    VOLTAGE_SOURCE, // Voltage source
    CURRENT_SOURCE, // Current source
    DIODE,          // Diode
    BJT,            // Bipolar Junction Transistor
    MOSFET,         // MOSFET
    JFET,           // JFET
    OPAMP,          // Operational Amplifier
    TRANSFORMER,    // Transformer
    TRANSMISSION_LINE, // Transmission Line
    SUBCIRCUIT      // Subcircuit
};

/**
 * @brief Class for a node in circuit simulation
 */
class CircuitNode {
public:
    /**
     * @brief Constructor for a node
     * @param id Node ID
     * @param name Node name
     */
    CircuitNode(int id, const std::string& name = "");

    /**
     * @brief Get the node ID
     * @return Node ID
     */
    int getId() const;

    /**
     * @brief Get the node name
     * @return Node name
     */
    const std::string& getName() const;

    /**
     * @brief Set the node name
     * @param name Node name
     */
    void setName(const std::string& name);

    /**
     * @brief Get the node voltage (DC analysis)
     * @return Node voltage
     */
    double getVoltage() const;

    /**
     * @brief Set the node voltage (DC analysis)
     * @param voltage Node voltage
     */
    void setVoltage(double voltage);

    /**
     * @brief Get the node voltage (AC analysis)
     * @return Node voltage as a complex number
     */
    std::complex<double> getComplexVoltage() const;

    /**
     * @brief Set the node voltage (AC analysis)
     * @param voltage Node voltage as a complex number
     */
    void setComplexVoltage(const std::complex<double>& voltage);

    /**
     * @brief Get the node voltage (transient analysis)
     * @param time Time
     * @return Node voltage at the specified time
     */
    double getVoltageAtTime(double time) const;

    /**
     * @brief Set the node voltage (transient analysis)
     * @param time Time
     * @param voltage Node voltage at the specified time
     */
    void setVoltageAtTime(double time, double voltage);

    /**
     * @brief Get the node voltage (frequency response analysis)
     * @param frequency Frequency
     * @return Node voltage at the specified frequency
     */
    std::complex<double> getVoltageAtFrequency(double frequency) const;

    /**
     * @brief Set the node voltage (frequency response analysis)
     * @param frequency Frequency
     * @param voltage Node voltage at the specified frequency
     */
    void setVoltageAtFrequency(double frequency, const std::complex<double>& voltage);

private:
    int m_id;
    std::string m_name;
    double m_voltage = 0.0;
    std::complex<double> m_complexVoltage = std::complex<double>(0.0, 0.0);
    std::map<double, double> m_voltageVsTime;
    std::map<double, std::complex<double>> m_voltageVsFrequency;
};

/**
 * @brief Base class for a component in circuit simulation
 */
class Component {
public:
    /**
     * @brief Constructor for a component
     * @param id Component ID
     * @param type Component type
     * @param name Component name
     * @param nodeIds Node IDs
     */
    Component(int id, ComponentType type, const std::string& name, const std::vector<int>& nodeIds);

    /**
     * @brief Virtual destructor
     */
    virtual ~Component() = default;

    /**
     * @brief Get the component ID
     * @return Component ID
     */
    int getId() const;

    /**
     * @brief Get the component type
     * @return Component type
     */
    ComponentType getType() const;

    /**
     * @brief Get the component name
     * @return Component name
     */
    const std::string& getName() const;

    /**
     * @brief Set the component name
     * @param name Component name
     */
    void setName(const std::string& name);

    /**
     * @brief Get the node IDs
     * @return Node IDs
     */
    const std::vector<int>& getNodeIds() const;

    /**
     * @brief Get the number of nodes
     * @return Number of nodes
     */
    size_t getNodeCount() const;

    /**
     * @brief Get the current through the component (DC analysis)
     * @return Current through the component
     */
    double getCurrent() const;

    /**
     * @brief Set the current through the component (DC analysis)
     * @param current Current through the component
     */
    void setCurrent(double current);

    /**
     * @brief Get the current through the component (AC analysis)
     * @return Current through the component as a complex number
     */
    std::complex<double> getComplexCurrent() const;

    /**
     * @brief Set the current through the component (AC analysis)
     * @param current Current through the component as a complex number
     */
    void setComplexCurrent(const std::complex<double>& current);

    /**
     * @brief Get the current through the component (transient analysis)
     * @param time Time
     * @return Current through the component at the specified time
     */
    double getCurrentAtTime(double time) const;

    /**
     * @brief Set the current through the component (transient analysis)
     * @param time Time
     * @param current Current through the component at the specified time
     */
    void setCurrentAtTime(double time, double current);

    /**
     * @brief Get the current through the component (frequency response analysis)
     * @param frequency Frequency
     * @return Current through the component at the specified frequency
     */
    std::complex<double> getCurrentAtFrequency(double frequency) const;

    /**
     * @brief Set the current through the component (frequency response analysis)
     * @param frequency Frequency
     * @param current Current through the component at the specified frequency
     */
    void setCurrentAtFrequency(double frequency, const std::complex<double>& current);

    /**
     * @brief Get the power dissipated by the component (DC analysis)
     * @return Power dissipated by the component
     */
    double getPower() const;

    /**
     * @brief Get the power dissipated by the component (AC analysis)
     * @return Power dissipated by the component as a complex number
     */
    std::complex<double> getComplexPower() const;

    /**
     * @brief Get the power dissipated by the component (transient analysis)
     * @param time Time
     * @return Power dissipated by the component at the specified time
     */
    double getPowerAtTime(double time) const;

    /**
     * @brief Get the power dissipated by the component (frequency response analysis)
     * @param frequency Frequency
     * @return Power dissipated by the component at the specified frequency
     */
    std::complex<double> getPowerAtFrequency(double frequency) const;

    /**
     * @brief Calculate the component's contribution to the system matrix (DC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     */
    virtual void contributeToSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<double>& rightHandSide,
                                         const std::map<int, CircuitNode>& nodes) const = 0;

    /**
     * @brief Calculate the component's contribution to the system matrix (AC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     * @param frequency Frequency
     */
    virtual void contributeToComplexSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<std::complex<double>>& rightHandSide,
                                               const std::map<int, CircuitNode>& nodes, double frequency) const = 0;

    /**
     * @brief Update the component's state (transient analysis)
     * @param time Time
     * @param timeStep Time step
     * @param nodes Map of node IDs to nodes
     */
    virtual void updateState(double time, double timeStep, const std::map<int, CircuitNode>& nodes) = 0;

protected:
    int m_id;
    ComponentType m_type;
    std::string m_name;
    std::vector<int> m_nodeIds;
    double m_current = 0.0;
    std::complex<double> m_complexCurrent = std::complex<double>(0.0, 0.0);
    std::map<double, double> m_currentVsTime;
    std::map<double, std::complex<double>> m_currentVsFrequency;
};

/**
 * @brief Class for a resistor in circuit simulation
 */
class Resistor : public Component {
public:
    /**
     * @brief Constructor for a resistor
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param resistance Resistance
     */
    Resistor(int id, const std::string& name, const std::vector<int>& nodeIds, double resistance);

    /**
     * @brief Get the resistance
     * @return Resistance
     */
    double getResistance() const;

    /**
     * @brief Set the resistance
     * @param resistance Resistance
     */
    void setResistance(double resistance);

    /**
     * @brief Calculate the component's contribution to the system matrix (DC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     */
    void contributeToSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<double>& rightHandSide,
                                 const std::map<int, CircuitNode>& nodes) const override;

    /**
     * @brief Calculate the component's contribution to the system matrix (AC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     * @param frequency Frequency
     */
    void contributeToComplexSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<std::complex<double>>& rightHandSide,
                                       const std::map<int, CircuitNode>& nodes, double frequency) const override;

    /**
     * @brief Update the component's state (transient analysis)
     * @param time Time
     * @param timeStep Time step
     * @param nodes Map of node IDs to nodes
     */
    void updateState(double time, double timeStep, const std::map<int, CircuitNode>& nodes) override;

private:
    double m_resistance;
};

/**
 * @brief Class for a capacitor in circuit simulation
 */
class Capacitor : public Component {
public:
    /**
     * @brief Constructor for a capacitor
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param capacitance Capacitance
     */
    Capacitor(int id, const std::string& name, const std::vector<int>& nodeIds, double capacitance);

    /**
     * @brief Get the capacitance
     * @return Capacitance
     */
    double getCapacitance() const;

    /**
     * @brief Set the capacitance
     * @param capacitance Capacitance
     */
    void setCapacitance(double capacitance);

    /**
     * @brief Calculate the component's contribution to the system matrix (DC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     */
    void contributeToSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<double>& rightHandSide,
                                 const std::map<int, CircuitNode>& nodes) const override;

    /**
     * @brief Calculate the component's contribution to the system matrix (AC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     * @param frequency Frequency
     */
    void contributeToComplexSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<std::complex<double>>& rightHandSide,
                                       const std::map<int, CircuitNode>& nodes, double frequency) const override;

    /**
     * @brief Update the component's state (transient analysis)
     * @param time Time
     * @param timeStep Time step
     * @param nodes Map of node IDs to nodes
     */
    void updateState(double time, double timeStep, const std::map<int, CircuitNode>& nodes) override;

private:
    double m_capacitance;
    double m_charge = 0.0;
};

/**
 * @brief Class for an inductor in circuit simulation
 */
class Inductor : public Component {
public:
    /**
     * @brief Constructor for an inductor
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param inductance Inductance
     */
    Inductor(int id, const std::string& name, const std::vector<int>& nodeIds, double inductance);

    /**
     * @brief Get the inductance
     * @return Inductance
     */
    double getInductance() const;

    /**
     * @brief Set the inductance
     * @param inductance Inductance
     */
    void setInductance(double inductance);

    /**
     * @brief Calculate the component's contribution to the system matrix (DC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     */
    void contributeToSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<double>& rightHandSide,
                                 const std::map<int, CircuitNode>& nodes) const override;

    /**
     * @brief Calculate the component's contribution to the system matrix (AC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     * @param frequency Frequency
     */
    void contributeToComplexSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<std::complex<double>>& rightHandSide,
                                       const std::map<int, CircuitNode>& nodes, double frequency) const override;

    /**
     * @brief Update the component's state (transient analysis)
     * @param time Time
     * @param timeStep Time step
     * @param nodes Map of node IDs to nodes
     */
    void updateState(double time, double timeStep, const std::map<int, CircuitNode>& nodes) override;

private:
    double m_inductance;
    double m_flux = 0.0;
};

/**
 * @brief Class for a voltage source in circuit simulation
 */
class VoltageSource : public Component {
public:
    /**
     * @brief Constructor for a voltage source
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param voltage Voltage
     */
    VoltageSource(int id, const std::string& name, const std::vector<int>& nodeIds, double voltage);

    /**
     * @brief Get the voltage
     * @return Voltage
     */
    double getVoltage() const;

    /**
     * @brief Set the voltage
     * @param voltage Voltage
     */
    void setVoltage(double voltage);

    /**
     * @brief Set the voltage function (for time-varying sources)
     * @param voltageFunc Voltage function
     */
    void setVoltageFunction(std::function<double(double)> voltageFunc);

    /**
     * @brief Get the voltage at a specific time
     * @param time Time
     * @return Voltage at the specified time
     */
    double getVoltageAtTime(double time) const;

    /**
     * @brief Calculate the component's contribution to the system matrix (DC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     */
    void contributeToSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<double>& rightHandSide,
                                 const std::map<int, CircuitNode>& nodes) const override;

    /**
     * @brief Calculate the component's contribution to the system matrix (AC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     * @param frequency Frequency
     */
    void contributeToComplexSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<std::complex<double>>& rightHandSide,
                                       const std::map<int, CircuitNode>& nodes, double frequency) const override;

    /**
     * @brief Update the component's state (transient analysis)
     * @param time Time
     * @param timeStep Time step
     * @param nodes Map of node IDs to nodes
     */
    void updateState(double time, double timeStep, const std::map<int, CircuitNode>& nodes) override;

private:
    double m_voltage;
    std::function<double(double)> m_voltageFunc;
};

/**
 * @brief Class for a current source in circuit simulation
 */
class CurrentSource : public Component {
public:
    /**
     * @brief Constructor for a current source
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param current Current
     */
    CurrentSource(int id, const std::string& name, const std::vector<int>& nodeIds, double current);

    /**
     * @brief Get the current
     * @return Current
     */
    double getCurrent() const;

    /**
     * @brief Set the current
     * @param current Current
     */
    void setCurrent(double current);

    /**
     * @brief Set the current function (for time-varying sources)
     * @param currentFunc Current function
     */
    void setCurrentFunction(std::function<double(double)> currentFunc);

    /**
     * @brief Get the current at a specific time
     * @param time Time
     * @return Current at the specified time
     */
    double getCurrentAtTime(double time) const;

    /**
     * @brief Calculate the component's contribution to the system matrix (DC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     */
    void contributeToSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<double>& rightHandSide,
                                 const std::map<int, CircuitNode>& nodes) const override;

    /**
     * @brief Calculate the component's contribution to the system matrix (AC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     * @param frequency Frequency
     */
    void contributeToComplexSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<std::complex<double>>& rightHandSide,
                                       const std::map<int, CircuitNode>& nodes, double frequency) const override;

    /**
     * @brief Update the component's state (transient analysis)
     * @param time Time
     * @param timeStep Time step
     * @param nodes Map of node IDs to nodes
     */
    void updateState(double time, double timeStep, const std::map<int, CircuitNode>& nodes) override;

private:
    double m_current;
    std::function<double(double)> m_currentFunc;
};

/**
 * @brief Class for a diode in circuit simulation
 */
class Diode : public Component {
public:
    /**
     * @brief Constructor for a diode
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param saturationCurrent Saturation current
     * @param emissionCoefficient Emission coefficient
     */
    Diode(int id, const std::string& name, const std::vector<int>& nodeIds, double saturationCurrent, double emissionCoefficient);

    /**
     * @brief Get the saturation current
     * @return Saturation current
     */
    double getSaturationCurrent() const;

    /**
     * @brief Set the saturation current
     * @param saturationCurrent Saturation current
     */
    void setSaturationCurrent(double saturationCurrent);

    /**
     * @brief Get the emission coefficient
     * @return Emission coefficient
     */
    double getEmissionCoefficient() const;

    /**
     * @brief Set the emission coefficient
     * @param emissionCoefficient Emission coefficient
     */
    void setEmissionCoefficient(double emissionCoefficient);

    /**
     * @brief Calculate the component's contribution to the system matrix (DC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     */
    void contributeToSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<double>& rightHandSide,
                                 const std::map<int, CircuitNode>& nodes) const override;

    /**
     * @brief Calculate the component's contribution to the system matrix (AC analysis)
     * @param systemMatrix System matrix
     * @param rightHandSide Right-hand side vector
     * @param nodes Map of node IDs to nodes
     * @param frequency Frequency
     */
    void contributeToComplexSystemMatrix(rebelcalc::Matrix& systemMatrix, std::vector<std::complex<double>>& rightHandSide,
                                       const std::map<int, CircuitNode>& nodes, double frequency) const override;

    /**
     * @brief Update the component's state (transient analysis)
     * @param time Time
     * @param timeStep Time step
     * @param nodes Map of node IDs to nodes
     */
    void updateState(double time, double timeStep, const std::map<int, CircuitNode>& nodes) override;

private:
    double m_saturationCurrent;
    double m_emissionCoefficient;
    double m_thermalVoltage = 0.026; // Default: 26 mV at room temperature
};

/**
 * @brief Class for a circuit model
 */
class CircuitModel {
public:
    /**
     * @brief Constructor for a circuit model
     * @param name Model name
     */
    CircuitModel(const std::string& name = "");

    /**
     * @brief Add a node to the model
     * @param id Node ID
     * @param name Node name
     * @return Reference to the added node
     */
    CircuitNode& addNode(int id, const std::string& name = "");

    /**
     * @brief Add a resistor to the model
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param resistance Resistance
     * @return Reference to the added resistor
     */
    Resistor& addResistor(int id, const std::string& name, const std::vector<int>& nodeIds, double resistance);

    /**
     * @brief Add a capacitor to the model
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param capacitance Capacitance
     * @return Reference to the added capacitor
     */
    Capacitor& addCapacitor(int id, const std::string& name, const std::vector<int>& nodeIds, double capacitance);

    /**
     * @brief Add an inductor to the model
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param inductance Inductance
     * @return Reference to the added inductor
     */
    Inductor& addInductor(int id, const std::string& name, const std::vector<int>& nodeIds, double inductance);

    /**
     * @brief Add a voltage source to the model
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param voltage Voltage
     * @return Reference to the added voltage source
     */
    VoltageSource& addVoltageSource(int id, const std::string& name, const std::vector<int>& nodeIds, double voltage);

    /**
     * @brief Add a current source to the model
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param current Current
     * @return Reference to the added current source
     */
    CurrentSource& addCurrentSource(int id, const std::string& name, const std::vector<int>& nodeIds, double current);

    /**
     * @brief Add a diode to the model
     * @param id Component ID
     * @param name Component name
     * @param nodeIds Node IDs (2 nodes)
     * @param saturationCurrent Saturation current
     * @param emissionCoefficient Emission coefficient
     * @return Reference to the added diode
     */
    Diode& addDiode(int id, const std::string& name, const std::vector<int>& nodeIds, double saturationCurrent, double emissionCoefficient);

    /**
     * @brief Get a node by ID
     * @param id Node ID
     * @return Reference to the node
     * @throws std::out_of_range if the node doesn't exist
     */
    CircuitNode& getNode(int id);

    /**
     * @brief Get a node by ID
     * @param id Node ID
     * @return Const reference to the node
     * @throws std::out_of_range if the node doesn't exist
     */
    const CircuitNode& getNode(int id) const;

    /**
     * @brief Get a component by ID
     * @param id Component ID
     * @return Reference to the component
     * @throws std::out_of_range if the component doesn't exist
     */
    Component& getComponent(int id);

    /**
     * @brief Get a component by ID
     * @param id Component ID
     * @return Const reference to the component
     * @throws std::out_of_range if the component doesn't exist
     */
    const Component& getComponent(int id) const;

    /**
     * @brief Get all nodes
     * @return Map of node IDs to nodes
     */
    const std::map<int, CircuitNode>& getNodes() const;

    /**
     * @brief Get all components
     * @return Map of component IDs to components
     */
    const std::map<int, std::unique_ptr<Component>>& getComponents() const;

    /**
     * @brief Get the model name
     * @return Model name
     */
    const std::string& getName() const;

    /**
     * @brief Set the model name
     * @param name Model name
     */
    void setName(const std::string& name);

    /**
     * @brief Import a model from a file
     * @param filename Filename
     * @return true if the model was imported successfully, false otherwise
     */
    bool importFromFile(const std::string& filename);

    /**
     * @brief Export the model to a file
     * @param filename Filename
     * @return true if the model was exported successfully, false otherwise
     */
    bool exportToFile(const std::string& filename) const;

private:
    std::string m_name;
    std::map<int, CircuitNode> m_nodes;
    std::map<int, std::unique_ptr<Component>> m_components;
};

/**
 * @brief Class for a circuit solver
 */
class CircuitSolver {
public:
    /**
     * @brief Constructor for a circuit solver
     * @param model Circuit model
     */
    CircuitSolver(const CircuitModel& model);

    /**
     * @brief Set the analysis type
     * @param type Analysis type
     */
    void setAnalysisType(CircuitAnalysisType type);

    /**
     * @brief Get the analysis type
     * @return Analysis type
     */
    CircuitAnalysisType getAnalysisType() const;

    /**
     * @brief Set the start frequency (for AC and frequency response analysis)
     * @param frequency Start frequency
     */
    void setStartFrequency(double frequency);

    /**
     * @brief Get the start frequency
     * @return Start frequency
     */
    double getStartFrequency() const;

    /**
     * @brief Set the stop frequency (for frequency response analysis)
     * @param frequency Stop frequency
     */
    void setStopFrequency(double frequency);

    /**
     * @brief Get the stop frequency
     * @return Stop frequency
     */
    double getStopFrequency() const;

    /**
     * @brief Set the number of frequency points (for frequency response analysis)
     * @param points Number of frequency points
     */
    void setFrequencyPoints(int points);

    /**
     * @brief Get the number of frequency points
     * @return Number of frequency points
     */
    int getFrequencyPoints() const;

    /**
     * @brief Set the start time (for transient analysis)
     * @param time Start time
     */
    void setStartTime(double time);

    /**
     * @brief Get the start time
     * @return Start time
     */
    double getStartTime() const;

    /**
     * @brief Set the stop time (for transient analysis)
     * @param time Stop time
     */
    void setStopTime(double time);

    /**
     * @brief Get the stop time
     * @return Stop time
     */
    double getStopTime() const;

    /**
     * @brief Set the time step (for transient analysis)
     * @param timeStep Time step
     */
    void setTimeStep(double timeStep);

    /**
     * @brief Get the time step
     * @return Time step
     */
    double getTimeStep() const;

    /**
     * @brief Set the maximum number of iterations (for nonlinear analysis)
     * @param iterations Maximum number of iterations
     */
    void setMaxIterations(int iterations);

    /**
     * @brief Get the maximum number of iterations
     * @return Maximum number of iterations
     */
    int getMaxIterations() const;

    /**
     * @brief Set the convergence tolerance (for nonlinear analysis)
     * @param tolerance Convergence tolerance
     */
    void setConvergenceTolerance(double tolerance);

    /**
     * @brief Get the convergence tolerance
     * @return Convergence tolerance
     */
    double getConvergenceTolerance() const;

    /**
     * @brief Solve the circuit model
     * @return true if the model was solved successfully, false otherwise
     */
    bool solve();

    /**
     * @brief Get the voltage at a node (DC analysis)
     * @param nodeId Node ID
     * @return Voltage at the node
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the node doesn't exist
     */
    double getVoltage(int nodeId) const;

    /**
     * @brief Get the voltage at a node (AC analysis)
     * @param nodeId Node ID
     * @return Voltage at the node as a complex number
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not AC
     * @throws std::out_of_range if the node doesn't exist
     */
    std::complex<double> getComplexVoltage(int nodeId) const;

    /**
     * @brief Get the voltage at a node (transient analysis)
     * @param nodeId Node ID
     * @param time Time
     * @return Voltage at the node at the specified time
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not TRANSIENT
     * @throws std::out_of_range if the node doesn't exist
     */
    double getVoltageAtTime(int nodeId, double time) const;

    /**
     * @brief Get the voltage at a node (frequency response analysis)
     * @param nodeId Node ID
     * @param frequency Frequency
     * @return Voltage at the node at the specified frequency
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not FREQUENCY
     * @throws std::out_of_range if the node doesn't exist
     */
    std::complex<double> getVoltageAtFrequency(int nodeId, double frequency) const;

    /**
     * @brief Get the current through a component (DC analysis)
     * @param componentId Component ID
     * @return Current through the component
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the component doesn't exist
     */
    double getCurrent(int componentId) const;

    /**
     * @brief Get the current through a component (AC analysis)
     * @param componentId Component ID
     * @return Current through the component as a complex number
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not AC
     * @throws std::out_of_range if the component doesn't exist
     */
    std::complex<double> getComplexCurrent(int componentId) const;

    /**
     * @brief Get the current through a component (transient analysis)
     * @param componentId Component ID
     * @param time Time
     * @return Current through the component at the specified time
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not TRANSIENT
     * @throws std::out_of_range if the component doesn't exist
     */
    double getCurrentAtTime(int componentId, double time) const;

    /**
     * @brief Get the current through a component (frequency response analysis)
     * @param componentId Component ID
     * @param frequency Frequency
     * @return Current through the component at the specified frequency
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not FREQUENCY
     * @throws std::out_of_range if the component doesn't exist
     */
    std::complex<double> getCurrentAtFrequency(int componentId, double frequency) const;

    /**
     * @brief Get the power dissipated by a component (DC analysis)
     * @param componentId Component ID
     * @return Power dissipated by the component
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the component doesn't exist
     */
    double getPower(int componentId) const;

    /**
     * @brief Get the power dissipated by a component (AC analysis)
     * @param componentId Component ID
     * @return Power dissipated by the component as a complex number
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not AC
     * @throws std::out_of_range if the component doesn't exist
     */
    std::complex<double> getComplexPower(int componentId) const;

    /**
     * @brief Get the power dissipated by a component (transient analysis)
     * @param componentId Component ID
     * @param time Time
     * @return Power dissipated by the component at the specified time
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not TRANSIENT
     * @throws std::out_of_range if the component doesn't exist
     */
    double getPowerAtTime(int componentId, double time) const;

    /**
     * @brief Get the power dissipated by a component (frequency response analysis)
     * @param componentId Component ID
     * @param frequency Frequency
     * @return Power dissipated by the component at the specified frequency
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not FREQUENCY
     * @throws std::out_of_range if the component doesn't exist
     */
    std::complex<double> getPowerAtFrequency(int componentId, double frequency) const;

    /**
     * @brief Get the total power dissipated by the circuit (DC analysis)
     * @return Total power dissipated by the circuit
     * @throws std::runtime_error if the model hasn't been solved
     */
    double getTotalPower() const;

    /**
     * @brief Get the total power dissipated by the circuit (AC analysis)
     * @return Total power dissipated by the circuit as a complex number
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not AC
     */
    std::complex<double> getTotalComplexPower() const;

    /**
     * @brief Get the total power dissipated by the circuit (transient analysis)
     * @param time Time
     * @return Total power dissipated by the circuit at the specified time
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not TRANSIENT
     */
    double getTotalPowerAtTime(double time) const;

    /**
     * @brief Get the total power dissipated by the circuit (frequency response analysis)
     * @param frequency Frequency
     * @return Total power dissipated by the circuit at the specified frequency
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not FREQUENCY
     */
    std::complex<double> getTotalPowerAtFrequency(double frequency) const;

    /**
     * @brief Export the results to a file
     * @param filename Filename
     * @return true if the results were exported successfully, false otherwise
     */
    bool exportResults(const std::string& filename) const;

private:
    const CircuitModel& m_model;
    CircuitAnalysisType m_analysisType = CircuitAnalysisType::DC;
    double m_startFrequency = 1.0;
    double m_stopFrequency = 1000.0;
    int m_frequencyPoints = 100;
    double m_startTime = 0.0;
    double m_stopTime = 1.0;
    double m_timeStep = 0.001;
    int m_maxIterations = 100;
    double m_convergenceTolerance = 1e-6;
    bool m_solved = false;

    // Helper methods
    void solveDC();
    void solveAC();
    void solveTransient();
    void solveFrequency();
    void solveNoise();
    void solveDistortion();
    void solveMonteCarlo();
    void solveWorstCase();
    void solveSensitivity();
};

} // namespace Solvers
} // namespace RebelCalc

#endif // REBELCALC_SOLVERS_CIRCUIT_H
