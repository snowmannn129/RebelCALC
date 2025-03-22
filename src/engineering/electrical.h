#ifndef REBELCALC_ENGINEERING_ELECTRICAL_H
#define REBELCALC_ENGINEERING_ELECTRICAL_H

#include <vector>
#include <complex>
#include <string>
#include <map>
#include <stdexcept>
#include <functional>
#include "../backend/matrix.h"

namespace RebelCalc {
namespace Engineering {
namespace Electrical {

/**
 * @brief Complex number type for electrical calculations
 */
using Complex = std::complex<double>;

/**
 * @brief Enumeration for circuit element types
 */
enum class ElementType {
    RESISTOR,
    CAPACITOR,
    INDUCTOR,
    VOLTAGE_SOURCE,
    CURRENT_SOURCE,
    DIODE,
    TRANSISTOR
};

/**
 * @brief A class representing a circuit element
 */
class CircuitElement {
public:
    /**
     * @brief Constructor for a circuit element
     * @param type Type of the element
     * @param value Value of the element (resistance, capacitance, inductance, etc.)
     * @param node1 First node
     * @param node2 Second node
     */
    CircuitElement(ElementType type, double value, int node1, int node2);
    
    /**
     * @brief Get the type of the element
     * @return Type of the element
     */
    ElementType getType() const;
    
    /**
     * @brief Get the value of the element
     * @return Value of the element
     */
    double getValue() const;
    
    /**
     * @brief Get the nodes of the element
     * @return Pair of nodes
     */
    std::pair<int, int> getNodes() const;
    
    /**
     * @brief Get the impedance of the element at a given frequency
     * @param frequency Frequency in Hz
     * @return Complex impedance
     */
    Complex getImpedance(double frequency) const;
    
private:
    ElementType m_type;
    double m_value;
    int m_node1;
    int m_node2;
};

/**
 * @brief A class for analyzing DC circuits
 */
class DCCircuit {
public:
    /**
     * @brief Constructor for a DC circuit
     */
    DCCircuit();
    
    /**
     * @brief Add a resistor to the circuit
     * @param resistance Resistance in ohms
     * @param node1 First node
     * @param node2 Second node
     * @return Index of the added element
     */
    size_t addResistor(double resistance, int node1, int node2);
    
    /**
     * @brief Add a voltage source to the circuit
     * @param voltage Voltage in volts
     * @param node1 First node (positive terminal)
     * @param node2 Second node (negative terminal)
     * @return Index of the added element
     */
    size_t addVoltageSource(double voltage, int node1, int node2);
    
    /**
     * @brief Add a current source to the circuit
     * @param current Current in amperes
     * @param node1 First node (current flows from node1 to node2)
     * @param node2 Second node
     * @return Index of the added element
     */
    size_t addCurrentSource(double current, int node1, int node2);
    
    /**
     * @brief Solve the circuit using nodal analysis
     * @return true if the circuit was solved successfully, false otherwise
     */
    bool solve();
    
    /**
     * @brief Get the voltage at a node
     * @param node Node index
     * @return Voltage at the node
     */
    double getNodeVoltage(int node) const;
    
    /**
     * @brief Get the current through an element
     * @param elementIndex Index of the element
     * @return Current through the element
     */
    double getElementCurrent(size_t elementIndex) const;
    
    /**
     * @brief Get the power dissipated by an element
     * @param elementIndex Index of the element
     * @return Power dissipated by the element
     */
    double getElementPower(size_t elementIndex) const;
    
    /**
     * @brief Get the total power dissipated by the circuit
     * @return Total power dissipated by the circuit
     */
    double getTotalPower() const;
    
    /**
     * @brief Get the number of nodes in the circuit
     * @return Number of nodes
     */
    int getNodeCount() const;
    
    /**
     * @brief Get the number of elements in the circuit
     * @return Number of elements
     */
    size_t getElementCount() const;
    
private:
    std::vector<CircuitElement> m_elements;
    std::vector<double> m_nodeVoltages;
    std::vector<double> m_elementCurrents;
    int m_nodeCount;
    bool m_solved;
};

/**
 * @brief A class for analyzing AC circuits
 */
class ACCircuit {
public:
    /**
     * @brief Constructor for an AC circuit
     * @param frequency Frequency in Hz
     */
    ACCircuit(double frequency);
    
    /**
     * @brief Add a resistor to the circuit
     * @param resistance Resistance in ohms
     * @param node1 First node
     * @param node2 Second node
     * @return Index of the added element
     */
    size_t addResistor(double resistance, int node1, int node2);
    
    /**
     * @brief Add a capacitor to the circuit
     * @param capacitance Capacitance in farads
     * @param node1 First node
     * @param node2 Second node
     * @return Index of the added element
     */
    size_t addCapacitor(double capacitance, int node1, int node2);
    
    /**
     * @brief Add an inductor to the circuit
     * @param inductance Inductance in henries
     * @param node1 First node
     * @param node2 Second node
     * @return Index of the added element
     */
    size_t addInductor(double inductance, int node1, int node2);
    
    /**
     * @brief Add a voltage source to the circuit
     * @param voltage Complex voltage in volts
     * @param node1 First node (positive terminal)
     * @param node2 Second node (negative terminal)
     * @return Index of the added element
     */
    size_t addVoltageSource(const Complex& voltage, int node1, int node2);
    
    /**
     * @brief Add a current source to the circuit
     * @param current Complex current in amperes
     * @param node1 First node (current flows from node1 to node2)
     * @param node2 Second node
     * @return Index of the added element
     */
    size_t addCurrentSource(const Complex& current, int node1, int node2);
    
    /**
     * @brief Solve the circuit using nodal analysis
     * @return true if the circuit was solved successfully, false otherwise
     */
    bool solve();
    
    /**
     * @brief Get the complex voltage at a node
     * @param node Node index
     * @return Complex voltage at the node
     */
    Complex getNodeVoltage(int node) const;
    
    /**
     * @brief Get the complex current through an element
     * @param elementIndex Index of the element
     * @return Complex current through the element
     */
    Complex getElementCurrent(size_t elementIndex) const;
    
    /**
     * @brief Get the complex power of an element
     * @param elementIndex Index of the element
     * @return Complex power of the element
     */
    Complex getElementPower(size_t elementIndex) const;
    
    /**
     * @brief Get the total complex power of the circuit
     * @return Total complex power of the circuit
     */
    Complex getTotalPower() const;
    
    /**
     * @brief Get the impedance of the circuit between two nodes
     * @param node1 First node
     * @param node2 Second node
     * @return Complex impedance between the nodes
     */
    Complex getImpedance(int node1, int node2) const;
    
    /**
     * @brief Get the frequency of the circuit
     * @return Frequency in Hz
     */
    double getFrequency() const;
    
    /**
     * @brief Set the frequency of the circuit
     * @param frequency Frequency in Hz
     */
    void setFrequency(double frequency);
    
    /**
     * @brief Get the number of nodes in the circuit
     * @return Number of nodes
     */
    int getNodeCount() const;
    
    /**
     * @brief Get the number of elements in the circuit
     * @return Number of elements
     */
    size_t getElementCount() const;
    
private:
    std::vector<CircuitElement> m_elements;
    std::vector<Complex> m_nodeVoltages;
    std::vector<Complex> m_elementCurrents;
    double m_frequency;
    int m_nodeCount;
    bool m_solved;
};

/**
 * @brief A class for analyzing digital circuits
 */
class DigitalCircuit {
public:
    /**
     * @brief Enumeration for logic gate types
     */
    enum class GateType {
        AND,
        OR,
        NOT,
        NAND,
        NOR,
        XOR,
        XNOR,
        BUFFER
    };
    
    /**
     * @brief Constructor for a digital circuit
     */
    DigitalCircuit();
    
    /**
     * @brief Add an input to the circuit
     * @param name Name of the input
     * @param value Initial value of the input
     * @return Index of the added input
     */
    size_t addInput(const std::string& name, bool value = false);
    
    /**
     * @brief Add a gate to the circuit
     * @param type Type of the gate
     * @param inputs Indices of the inputs to the gate
     * @return Index of the added gate
     */
    size_t addGate(GateType type, const std::vector<size_t>& inputs);
    
    /**
     * @brief Add an output to the circuit
     * @param name Name of the output
     * @param input Index of the input to the output
     * @return Index of the added output
     */
    size_t addOutput(const std::string& name, size_t input);
    
    /**
     * @brief Set the value of an input
     * @param inputIndex Index of the input
     * @param value Value to set
     */
    void setInput(size_t inputIndex, bool value);
    
    /**
     * @brief Get the value of an output
     * @param outputIndex Index of the output
     * @return Value of the output
     */
    bool getOutput(size_t outputIndex) const;
    
    /**
     * @brief Get the value of an output by name
     * @param name Name of the output
     * @return Value of the output
     */
    bool getOutput(const std::string& name) const;
    
    /**
     * @brief Simulate the circuit
     */
    void simulate();
    
    /**
     * @brief Get the truth table of the circuit
     * @return Truth table as a vector of input-output pairs
     */
    std::vector<std::pair<std::vector<bool>, std::vector<bool>>> getTruthTable() const;
    
    /**
     * @brief Get the number of inputs in the circuit
     * @return Number of inputs
     */
    size_t getInputCount() const;
    
    /**
     * @brief Get the number of gates in the circuit
     * @return Number of gates
     */
    size_t getGateCount() const;
    
    /**
     * @brief Get the number of outputs in the circuit
     * @return Number of outputs
     */
    size_t getOutputCount() const;
    
private:
    struct Input {
        std::string name;
        bool value;
    };
    
    struct Gate {
        GateType type;
        std::vector<size_t> inputs;
        bool output;
    };
    
    struct Output {
        std::string name;
        size_t input;
    };
    
    std::vector<Input> m_inputs;
    std::vector<Gate> m_gates;
    std::vector<Output> m_outputs;
    
    bool evaluateGate(const Gate& gate) const;
};

/**
 * @brief A class for analyzing power systems
 */
class PowerSystem {
public:
    /**
     * @brief Constructor for a power system
     * @param baseVoltage Base voltage in volts
     * @param basePower Base power in watts
     */
    PowerSystem(double baseVoltage, double basePower);
    
    /**
     * @brief Add a bus to the system
     * @param name Name of the bus
     * @param type Type of the bus (0 = slack, 1 = PV, 2 = PQ)
     * @param voltage Voltage magnitude in per unit
     * @param angle Voltage angle in radians
     * @param activePower Active power in per unit
     * @param reactivePower Reactive power in per unit
     * @return Index of the added bus
     */
    size_t addBus(const std::string& name, int type, double voltage, double angle, double activePower, double reactivePower);
    
    /**
     * @brief Add a line to the system
     * @param fromBus Index of the from bus
     * @param toBus Index of the to bus
     * @param resistance Resistance in per unit
     * @param reactance Reactance in per unit
     * @param susceptance Susceptance in per unit
     * @return Index of the added line
     */
    size_t addLine(size_t fromBus, size_t toBus, double resistance, double reactance, double susceptance);
    
    /**
     * @brief Add a transformer to the system
     * @param fromBus Index of the from bus
     * @param toBus Index of the to bus
     * @param resistance Resistance in per unit
     * @param reactance Reactance in per unit
     * @param ratio Transformer ratio
     * @return Index of the added transformer
     */
    size_t addTransformer(size_t fromBus, size_t toBus, double resistance, double reactance, double ratio);
    
    /**
     * @brief Solve the power flow using the Newton-Raphson method
     * @param maxIterations Maximum number of iterations
     * @param tolerance Convergence tolerance
     * @return true if the power flow was solved successfully, false otherwise
     */
    bool solvePowerFlow(int maxIterations = 10, double tolerance = 1e-6);
    
    /**
     * @brief Get the voltage magnitude at a bus
     * @param busIndex Index of the bus
     * @return Voltage magnitude in per unit
     */
    double getBusVoltageMagnitude(size_t busIndex) const;
    
    /**
     * @brief Get the voltage angle at a bus
     * @param busIndex Index of the bus
     * @return Voltage angle in radians
     */
    double getBusVoltageAngle(size_t busIndex) const;
    
    /**
     * @brief Get the active power at a bus
     * @param busIndex Index of the bus
     * @return Active power in per unit
     */
    double getBusActivePower(size_t busIndex) const;
    
    /**
     * @brief Get the reactive power at a bus
     * @param busIndex Index of the bus
     * @return Reactive power in per unit
     */
    double getBusReactivePower(size_t busIndex) const;
    
    /**
     * @brief Get the active power flow in a line or transformer
     * @param elementIndex Index of the line or transformer
     * @param fromBus true for from bus to to bus, false for to bus to from bus
     * @return Active power flow in per unit
     */
    double getElementActivePower(size_t elementIndex, bool fromBus = true) const;
    
    /**
     * @brief Get the reactive power flow in a line or transformer
     * @param elementIndex Index of the line or transformer
     * @param fromBus true for from bus to to bus, false for to bus to from bus
     * @return Reactive power flow in per unit
     */
    double getElementReactivePower(size_t elementIndex, bool fromBus = true) const;
    
    /**
     * @brief Get the current magnitude in a line or transformer
     * @param elementIndex Index of the line or transformer
     * @param fromBus true for from bus to to bus, false for to bus to from bus
     * @return Current magnitude in per unit
     */
    double getElementCurrentMagnitude(size_t elementIndex, bool fromBus = true) const;
    
    /**
     * @brief Get the power loss in a line or transformer
     * @param elementIndex Index of the line or transformer
     * @return Power loss in per unit
     */
    double getElementPowerLoss(size_t elementIndex) const;
    
    /**
     * @brief Get the total power loss in the system
     * @return Total power loss in per unit
     */
    double getTotalPowerLoss() const;
    
    /**
     * @brief Get the number of buses in the system
     * @return Number of buses
     */
    size_t getBusCount() const;
    
    /**
     * @brief Get the number of lines and transformers in the system
     * @return Number of lines and transformers
     */
    size_t getElementCount() const;
    
private:
    struct Bus {
        std::string name;
        int type;
        double voltage;
        double angle;
        double activePower;
        double reactivePower;
    };
    
    struct Element {
        size_t fromBus;
        size_t toBus;
        double resistance;
        double reactance;
        double susceptance;
        double ratio;
        bool isTransformer;
    };
    
    std::vector<Bus> m_buses;
    std::vector<Element> m_elements;
    double m_baseVoltage;
    double m_basePower;
    bool m_solved;
};

/**
 * @brief Calculate the resistance of a conductor
 * 
 * @param resistivity Resistivity of the material in ohm-meters
 * @param length Length of the conductor in meters
 * @param area Cross-sectional area of the conductor in square meters
 * @return Resistance in ohms
 */
double resistance(double resistivity, double length, double area);

/**
 * @brief Calculate the resistivity of a material at a given temperature
 * 
 * @param resistivity0 Resistivity at reference temperature in ohm-meters
 * @param alpha Temperature coefficient of resistivity in 1/°C
 * @param temperature Temperature in °C
 * @param temperature0 Reference temperature in °C
 * @return Resistivity at the given temperature in ohm-meters
 */
double resistivityTemperature(double resistivity0, double alpha, double temperature, double temperature0 = 20.0);

/**
 * @brief Calculate the capacitance of a parallel plate capacitor
 * 
 * @param area Area of the plates in square meters
 * @param distance Distance between the plates in meters
 * @param permittivity Permittivity of the dielectric in F/m
 * @return Capacitance in farads
 */
double capacitanceParallelPlate(double area, double distance, double permittivity);

/**
 * @brief Calculate the capacitance of a cylindrical capacitor
 * 
 * @param length Length of the capacitor in meters
 * @param innerRadius Inner radius in meters
 * @param outerRadius Outer radius in meters
 * @param permittivity Permittivity of the dielectric in F/m
 * @return Capacitance in farads
 */
double capacitanceCylindrical(double length, double innerRadius, double outerRadius, double permittivity);

/**
 * @brief Calculate the inductance of a solenoid
 * 
 * @param turns Number of turns
 * @param area Cross-sectional area in square meters
 * @param length Length of the solenoid in meters
 * @param permeability Permeability of the core in H/m
 * @return Inductance in henries
 */
double inductanceSolenoid(int turns, double area, double length, double permeability);

/**
 * @brief Calculate the inductance of a toroid
 * 
 * @param turns Number of turns
 * @param innerRadius Inner radius in meters
 * @param outerRadius Outer radius in meters
 * @param permeability Permeability of the core in H/m
 * @return Inductance in henries
 */
double inductanceToroid(int turns, double innerRadius, double outerRadius, double permeability);

/**
 * @brief Calculate the impedance of a resistor
 * 
 * @param resistance Resistance in ohms
 * @return Complex impedance
 */
Complex impedanceResistor(double resistance);

/**
 * @brief Calculate the impedance of a capacitor
 * 
 * @param capacitance Capacitance in farads
 * @param frequency Frequency in Hz
 * @return Complex impedance
 */
Complex impedanceCapacitor(double capacitance, double frequency);

/**
 * @brief Calculate the impedance of an inductor
 * 
 * @param inductance Inductance in henries
 * @param frequency Frequency in Hz
 * @return Complex impedance
 */
Complex impedanceInductor(double inductance, double frequency);

/**
 * @brief Calculate the impedance of a series RLC circuit
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @param capacitance Capacitance in farads
 * @param frequency Frequency in Hz
 * @return Complex impedance
 */
Complex impedanceSeriesRLC(double resistance, double inductance, double capacitance, double frequency);

/**
 * @brief Calculate the impedance of a parallel RLC circuit
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @param capacitance Capacitance in farads
 * @param frequency Frequency in Hz
 * @return Complex impedance
 */
Complex impedanceParallelRLC(double resistance, double inductance, double capacitance, double frequency);

/**
 * @brief Calculate the resonant frequency of an LC circuit
 * 
 * @param inductance Inductance in henries
 * @param capacitance Capacitance in farads
 * @return Resonant frequency in Hz
 */
double resonantFrequency(double inductance, double capacitance);

/**
 * @brief Calculate the quality factor of a series RLC circuit
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @param capacitance Capacitance in farads
 * @return Quality factor
 */
double qualityFactorSeries(double resistance, double inductance, double capacitance);

/**
 * @brief Calculate the quality factor of a parallel RLC circuit
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @param capacitance Capacitance in farads
 * @return Quality factor
 */
double qualityFactorParallel(double resistance, double inductance, double capacitance);

/**
 * @brief Calculate the bandwidth of an RLC circuit
 * 
 * @param resonantFrequency Resonant frequency in Hz
 * @param qualityFactor Quality factor
 * @return Bandwidth in Hz
 */
double bandwidth(double resonantFrequency, double qualityFactor);

/**
 * @brief Calculate the cutoff frequency of an RC low-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param capacitance Capacitance in farads
 * @return Cutoff frequency in Hz
 */
double cutoffFrequencyLowPassRC(double resistance, double capacitance);

/**
 * @brief Calculate the cutoff frequency of an RC high-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param capacitance Capacitance in farads
 * @return Cutoff frequency in Hz
 */
double cutoffFrequencyHighPassRC(double resistance, double capacitance);

/**
 * @brief Calculate the cutoff frequency of an RL low-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @return Cutoff frequency in Hz
 */
double cutoffFrequencyLowPassRL(double resistance, double inductance);

/**
 * @brief Calculate the cutoff frequency of an RL high-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @return Cutoff frequency in Hz
 */
double cutoffFrequencyHighPassRL(double resistance, double inductance);

/**
 * @brief Calculate the transfer function of an RC low-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param capacitance Capacitance in farads
 * @param frequency Frequency in Hz
 * @return Complex transfer function
 */
Complex transferFunctionLowPassRC(double resistance, double capacitance, double frequency);

/**
 * @brief Calculate the transfer function of an RC high-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param capacitance Capacitance in farads
 * @param frequency Frequency in Hz
 * @return Complex transfer function
 */
Complex transferFunctionHighPassRC(double resistance, double capacitance, double frequency);

/**
 * @brief Calculate the transfer function of an RL low-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @param frequency Frequency in Hz
 * @return Complex transfer function
 */
Complex transferFunctionLowPassRL(double resistance, double inductance, double frequency);

/**
 * @brief Calculate the transfer function of an RL high-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @param frequency Frequency in Hz
 * @return Complex transfer function
 */
Complex transferFunctionHighPassRL(double resistance, double inductance, double frequency);

/**
 * @brief Calculate the transfer function of an RLC band-pass filter
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @param capacitance Capacitance in farads
 * @param frequency Frequency in Hz
 * @return Complex transfer function
 */
Complex transferFunctionBandPassRLC(double resistance, double inductance, double capacitance, double frequency);

/**
 * @brief Calculate the transfer function of an RLC band-stop filter
 * 
 * @param resistance Resistance in ohms
 * @param inductance Inductance in henries
 * @param capacitance Capacitance in farads
 * @param frequency Frequency in Hz
 * @return Complex transfer function
 */
Complex transferFunctionBandStopRLC(double resistance, double inductance, double capacitance, double frequency);

/**
 * @brief Calculate the gain in decibels
 * 
 * @param inputPower Input power in watts
 * @param outputPower Output power in watts
 * @return Gain in dB
 */
double gainDecibels(double inputPower, double outputPower);

/**
 * @brief Calculate the gain in decibels from a voltage ratio
 * 
 * @param inputVoltage Input voltage in volts
 * @param outputVoltage Output voltage in volts
 * @return Gain in dB
 */
double gainDecibelsVoltage(double inputVoltage, double outputVoltage);

/**
 * @brief Calculate the gain in decibels from a current ratio
 * 
 * @param inputCurrent Input current in amperes
 * @param outputCurrent Output current in amperes
 * @return Gain in dB
 */
double gainDecibelsCurrent(double inputCurrent, double outputCurrent);

/**
 * @brief Calculate the power in a DC circuit
 * 
 * @param voltage Voltage in volts
 * @param current Current in amperes
 * @return Power in watts
 */
double powerDC(double voltage, double current);

/**
 * @brief Calculate the power in an AC circuit
 * 
 * @param voltage Voltage magnitude in volts
 * @param current Current magnitude in amperes
 * @param powerFactor Power factor
 * @return Power in watts
 */
double powerAC(double voltage, double current, double powerFactor);

/**
 * @brief Calculate the complex power in an AC circuit
 * 
 * @param voltage Complex voltage in volts
 * @param current Complex current in amperes
 * @return Complex power in VA
 */
Complex powerComplex(const Complex& voltage, const Complex& current);

/**
 * @brief Calculate the active power from complex power
 * 
 * @param power Complex power in VA
 * @return Active power in watts
 */
double activePower(const Complex& power);

/**
 * @brief Calculate the reactive power from complex power
 * 
 * @param power Complex power in VA
 * @return Reactive power in VAR
 */
double reactivePower(const Complex& power);

/**
 * @brief Calculate the apparent power from complex power
 * 
 * @param power Complex power in VA
 * @return Apparent power in VA
 */
double apparentPower(const Complex& power);

/**
 * @brief Calculate the power factor from complex power
 * 
 * @param power Complex power in VA
 * @return Power factor
 */
double powerFactor(const Complex& power);

/**
 * @brief Calculate the power factor angle from complex power
 * 
 * @param power Complex power in VA
 * @return Power factor angle in radians
 */
double powerFactorAngle(const Complex& power);

/**
 * @brief Calculate the energy in a DC circuit
 * 
 * @param power Power in watts
 * @param time Time in seconds
 * @return Energy in joules
 */
double energyDC(double power, double time);

/**
 * @brief Calculate the energy in an AC circuit
 * 
 * @param power Power in watts
 * @param time Time in seconds
 * @return Energy in joules
 */
double energyAC(double power, double time);

/**
 * @brief Calculate the efficiency of a system
 * 
 * @param outputPower Output power in watts
 * @param inputPower Input power in watts
 * @return Efficiency as a fraction
 */
double efficiency(double outputPower, double inputPower);

/**
 * @brief Calculate the voltage divider output
 * 
 * @param inputVoltage Input voltage in volts
 * @param resistance1 First resistance in ohms
 * @param resistance2 Second resistance in ohms
 * @return Output voltage in volts
 */
double voltageDivider(double inputVoltage, double resistance1, double resistance2);

/**
 * @brief Calculate the current divider output
 * 
 * @param inputCurrent Input current in amperes
 * @param resistance1 First resistance in ohms
 * @param resistance2 Second resistance in ohms
 * @return Output current in amperes
 */
double currentDivider(double inputCurrent, double resistance1, double resistance2);

/**
 * @brief Calculate the Thevenin equivalent voltage
 * 
 * @param openCircuitVoltage Open-circuit voltage in volts
 * @return Thevenin equivalent voltage in volts
 */
double theveninVoltage(double openCircuitVoltage);

/**
 * @brief Calculate the Thevenin equivalent resistance
 * 
 * @param openCircuitVoltage Open-circuit voltage in volts
 * @param shortCircuitCurrent Short-circuit current in amperes
 * @return Thevenin equivalent resistance in ohms
 */
double theveninResistance(double openCircuitVoltage, double shortCircuitCurrent);

/**
 * @brief Calculate the Norton equivalent current
 * 
 * @param shortCircuitCurrent Short-circuit current in amperes
 * @return Norton equivalent current in amperes
 */
double nortonCurrent(double shortCircuitCurrent);

/**
 * @brief Calculate the Norton equivalent resistance
 * 
 * @param openCircuitVoltage Open-circuit voltage in volts
 * @param shortCircuitCurrent Short-circuit current in amperes
 * @return Norton equivalent resistance in ohms
 */
double nortonResistance(double openCircuitVoltage, double shortCircuitCurrent);

/**
 * @brief Calculate the maximum power transfer
 * 
 * @param sourceVoltage Source voltage in volts
 * @param sourceResistance Source resistance in ohms
 * @return Maximum power in watts
 */
double maximumPowerTransfer(double sourceVoltage, double sourceResistance);

/**
 * @brief Calculate the load resistance for maximum power transfer
 * 
 * @param sourceResistance Source resistance in ohms
 * @return Load resistance in ohms
 */
double loadResistanceMaximumPowerTransfer(double sourceResistance);

} // namespace Electrical
} // namespace Engineering
} // namespace RebelCalc

#endif // REBELCALC_ENGINEERING_ELECTRICAL_H
