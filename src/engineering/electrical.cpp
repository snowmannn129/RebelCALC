#include "electrical.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace RebelCalc {
namespace Engineering {
namespace Electrical {

//
// CircuitElement implementation
//

CircuitElement::CircuitElement(ElementType type, double value, int node1, int node2)
    : m_type(type), m_value(value), m_node1(node1), m_node2(node2) {
}

ElementType CircuitElement::getType() const {
    return m_type;
}

double CircuitElement::getValue() const {
    return m_value;
}

std::pair<int, int> CircuitElement::getNodes() const {
    return {m_node1, m_node2};
}

Complex CircuitElement::getImpedance(double frequency) const {
    switch (m_type) {
        case ElementType::RESISTOR:
            return impedanceResistor(m_value);
        case ElementType::CAPACITOR:
            return impedanceCapacitor(m_value, frequency);
        case ElementType::INDUCTOR:
            return impedanceInductor(m_value, frequency);
        default:
            return Complex(0.0, 0.0); // Ideal voltage/current sources have zero impedance
    }
}

//
// DCCircuit implementation
//

DCCircuit::DCCircuit() : m_nodeCount(0), m_solved(false) {
}

size_t DCCircuit::addResistor(double resistance, int node1, int node2) {
    m_elements.emplace_back(ElementType::RESISTOR, resistance, node1, node2);
    m_nodeCount = std::max(m_nodeCount, std::max(node1, node2) + 1);
    m_solved = false;
    return m_elements.size() - 1;
}

size_t DCCircuit::addVoltageSource(double voltage, int node1, int node2) {
    m_elements.emplace_back(ElementType::VOLTAGE_SOURCE, voltage, node1, node2);
    m_nodeCount = std::max(m_nodeCount, std::max(node1, node2) + 1);
    m_solved = false;
    return m_elements.size() - 1;
}

size_t DCCircuit::addCurrentSource(double current, int node1, int node2) {
    m_elements.emplace_back(ElementType::CURRENT_SOURCE, current, node1, node2);
    m_nodeCount = std::max(m_nodeCount, std::max(node1, node2) + 1);
    m_solved = false;
    return m_elements.size() - 1;
}

bool DCCircuit::solve() {
    if (m_nodeCount == 0 || m_elements.empty()) {
        return false;
    }
    
    // Count voltage sources
    int voltageSourceCount = 0;
    for (const auto& element : m_elements) {
        if (element.getType() == ElementType::VOLTAGE_SOURCE) {
            voltageSourceCount++;
        }
    }
    
    // Set up the system of equations
    int matrixSize = m_nodeCount - 1 + voltageSourceCount; // -1 because node 0 is ground
    rebelcalc::Matrix A(matrixSize, matrixSize);
    std::vector<double> b(matrixSize, 0.0);
    
    // Fill the conductance matrix and the right-hand side vector
    for (const auto& element : m_elements) {
        auto [node1, node2] = element.getNodes();
        
        switch (element.getType()) {
            case ElementType::RESISTOR: {
                double conductance = 1.0 / element.getValue();
                
                // Add conductance to the matrix
                if (node1 > 0) {
                    A(node1 - 1, node1 - 1) += conductance;
                }
                if (node2 > 0) {
                    A(node2 - 1, node2 - 1) += conductance;
                }
                if (node1 > 0 && node2 > 0) {
                    A(node1 - 1, node2 - 1) -= conductance;
                    A(node2 - 1, node1 - 1) -= conductance;
                }
                break;
            }
            case ElementType::CURRENT_SOURCE: {
                double current = element.getValue();
                
                // Add current to the right-hand side vector
                if (node1 > 0) {
                    b[node1 - 1] -= current;
                }
                if (node2 > 0) {
                    b[node2 - 1] += current;
                }
                break;
            }
            case ElementType::VOLTAGE_SOURCE: {
                // Voltage sources require additional equations
                int voltageIndex = m_nodeCount - 1 + voltageSourceCount--;
                double voltage = element.getValue();
                
                // Add voltage source to the matrix
                if (node1 > 0) {
                    A(node1 - 1, voltageIndex) = 1.0;
                    A(voltageIndex, node1 - 1) = 1.0;
                }
                if (node2 > 0) {
                    A(node2 - 1, voltageIndex) = -1.0;
                    A(voltageIndex, node2 - 1) = -1.0;
                }
                
                // Add voltage to the right-hand side vector
                b[voltageIndex] = voltage;
                break;
            }
            default:
                break;
        }
    }
    
    // Solve the system of equations
    try {
        std::vector<double> x = A.solve(b);
        
        // Extract node voltages
        m_nodeVoltages.resize(m_nodeCount);
        m_nodeVoltages[0] = 0.0; // Ground node
        for (int i = 1; i < m_nodeCount; ++i) {
            m_nodeVoltages[i] = x[i - 1];
        }
        
        // Calculate element currents
        m_elementCurrents.resize(m_elements.size());
        for (size_t i = 0; i < m_elements.size(); ++i) {
            const auto& element = m_elements[i];
            auto [node1, node2] = element.getNodes();
            double voltage1 = (node1 >= 0 && node1 < m_nodeCount) ? m_nodeVoltages[node1] : 0.0;
            double voltage2 = (node2 >= 0 && node2 < m_nodeCount) ? m_nodeVoltages[node2] : 0.0;
            double voltageDiff = voltage1 - voltage2;
            
            switch (element.getType()) {
                case ElementType::RESISTOR:
                    m_elementCurrents[i] = voltageDiff / element.getValue();
                    break;
                case ElementType::VOLTAGE_SOURCE:
                    // Current through voltage source is calculated from the solution
                    m_elementCurrents[i] = x[m_nodeCount - 1 + i];
                    break;
                case ElementType::CURRENT_SOURCE:
                    m_elementCurrents[i] = element.getValue();
                    break;
                default:
                    m_elementCurrents[i] = 0.0;
                    break;
            }
        }
        
        m_solved = true;
        return true;
    } catch (const std::exception&) {
        m_solved = false;
        return false;
    }
}

double DCCircuit::getNodeVoltage(int node) const {
    if (!m_solved || node < 0 || node >= m_nodeCount) {
        throw std::out_of_range("Invalid node index or circuit not solved");
    }
    
    return m_nodeVoltages[node];
}

double DCCircuit::getElementCurrent(size_t elementIndex) const {
    if (!m_solved || elementIndex >= m_elements.size()) {
        throw std::out_of_range("Invalid element index or circuit not solved");
    }
    
    return m_elementCurrents[elementIndex];
}

double DCCircuit::getElementPower(size_t elementIndex) const {
    if (!m_solved || elementIndex >= m_elements.size()) {
        throw std::out_of_range("Invalid element index or circuit not solved");
    }
    
    const auto& element = m_elements[elementIndex];
    double current = m_elementCurrents[elementIndex];
    
    switch (element.getType()) {
        case ElementType::RESISTOR:
            return current * current * element.getValue();
        case ElementType::VOLTAGE_SOURCE: {
            auto [node1, node2] = element.getNodes();
            double voltage1 = (node1 >= 0 && node1 < m_nodeCount) ? m_nodeVoltages[node1] : 0.0;
            double voltage2 = (node2 >= 0 && node2 < m_nodeCount) ? m_nodeVoltages[node2] : 0.0;
            return current * (voltage1 - voltage2);
        }
        case ElementType::CURRENT_SOURCE: {
            auto [node1, node2] = element.getNodes();
            double voltage1 = (node1 >= 0 && node1 < m_nodeCount) ? m_nodeVoltages[node1] : 0.0;
            double voltage2 = (node2 >= 0 && node2 < m_nodeCount) ? m_nodeVoltages[node2] : 0.0;
            return element.getValue() * (voltage1 - voltage2);
        }
        default:
            return 0.0;
    }
}

double DCCircuit::getTotalPower() const {
    if (!m_solved) {
        throw std::runtime_error("Circuit not solved");
    }
    
    double totalPower = 0.0;
    for (size_t i = 0; i < m_elements.size(); ++i) {
        totalPower += getElementPower(i);
    }
    
    return totalPower;
}

int DCCircuit::getNodeCount() const {
    return m_nodeCount;
}

size_t DCCircuit::getElementCount() const {
    return m_elements.size();
}

//
// ACCircuit implementation
//

ACCircuit::ACCircuit(double frequency) : m_frequency(frequency), m_nodeCount(0), m_solved(false) {
}

size_t ACCircuit::addResistor(double resistance, int node1, int node2) {
    m_elements.emplace_back(ElementType::RESISTOR, resistance, node1, node2);
    m_nodeCount = std::max(m_nodeCount, std::max(node1, node2) + 1);
    m_solved = false;
    return m_elements.size() - 1;
}

size_t ACCircuit::addCapacitor(double capacitance, int node1, int node2) {
    m_elements.emplace_back(ElementType::CAPACITOR, capacitance, node1, node2);
    m_nodeCount = std::max(m_nodeCount, std::max(node1, node2) + 1);
    m_solved = false;
    return m_elements.size() - 1;
}

size_t ACCircuit::addInductor(double inductance, int node1, int node2) {
    m_elements.emplace_back(ElementType::INDUCTOR, inductance, node1, node2);
    m_nodeCount = std::max(m_nodeCount, std::max(node1, node2) + 1);
    m_solved = false;
    return m_elements.size() - 1;
}

size_t ACCircuit::addVoltageSource(const Complex& voltage, int node1, int node2) {
    // Store the complex voltage as the real part of the value
    // This is a simplification for this implementation
    m_elements.emplace_back(ElementType::VOLTAGE_SOURCE, voltage.real(), node1, node2);
    m_nodeCount = std::max(m_nodeCount, std::max(node1, node2) + 1);
    m_solved = false;
    return m_elements.size() - 1;
}

size_t ACCircuit::addCurrentSource(const Complex& current, int node1, int node2) {
    // Store the complex current as the real part of the value
    // This is a simplification for this implementation
    m_elements.emplace_back(ElementType::CURRENT_SOURCE, current.real(), node1, node2);
    m_nodeCount = std::max(m_nodeCount, std::max(node1, node2) + 1);
    m_solved = false;
    return m_elements.size() - 1;
}

bool ACCircuit::solve() {
    // This is a simplified implementation of AC circuit analysis
    // A complete implementation would use complex numbers throughout
    
    // For now, we'll just return false to indicate that the circuit was not solved
    m_solved = false;
    return false;
}

Complex ACCircuit::getNodeVoltage(int node) const {
    if (!m_solved || node < 0 || node >= m_nodeCount) {
        throw std::out_of_range("Invalid node index or circuit not solved");
    }
    
    return m_nodeVoltages[node];
}

Complex ACCircuit::getElementCurrent(size_t elementIndex) const {
    if (!m_solved || elementIndex >= m_elements.size()) {
        throw std::out_of_range("Invalid element index or circuit not solved");
    }
    
    return m_elementCurrents[elementIndex];
}

Complex ACCircuit::getElementPower(size_t elementIndex) const {
    if (!m_solved || elementIndex >= m_elements.size()) {
        throw std::out_of_range("Invalid element index or circuit not solved");
    }
    
    const auto& element = m_elements[elementIndex];
    Complex current = m_elementCurrents[elementIndex];
    
    switch (element.getType()) {
        case ElementType::RESISTOR: {
            // P = I^2 * R
            double currentMagnitude = std::abs(current);
            return Complex(currentMagnitude * currentMagnitude * element.getValue(), 0.0);
        }
        case ElementType::CAPACITOR: {
            // P = I^2 * Xc, where Xc = 1/(2*pi*f*C)
            double currentMagnitude = std::abs(current);
            double reactance = 1.0 / (2.0 * M_PI * m_frequency * element.getValue());
            return Complex(0.0, -currentMagnitude * currentMagnitude * reactance);
        }
        case ElementType::INDUCTOR: {
            // P = I^2 * Xl, where Xl = 2*pi*f*L
            double currentMagnitude = std::abs(current);
            double reactance = 2.0 * M_PI * m_frequency * element.getValue();
            return Complex(0.0, currentMagnitude * currentMagnitude * reactance);
        }
        case ElementType::VOLTAGE_SOURCE: {
            auto [node1, node2] = element.getNodes();
            Complex voltage1 = (node1 >= 0 && node1 < m_nodeCount) ? m_nodeVoltages[node1] : Complex(0.0, 0.0);
            Complex voltage2 = (node2 >= 0 && node2 < m_nodeCount) ? m_nodeVoltages[node2] : Complex(0.0, 0.0);
            return current * (voltage1 - voltage2);
        }
        case ElementType::CURRENT_SOURCE: {
            auto [node1, node2] = element.getNodes();
            Complex voltage1 = (node1 >= 0 && node1 < m_nodeCount) ? m_nodeVoltages[node1] : Complex(0.0, 0.0);
            Complex voltage2 = (node2 >= 0 && node2 < m_nodeCount) ? m_nodeVoltages[node2] : Complex(0.0, 0.0);
            return Complex(element.getValue(), 0.0) * (voltage1 - voltage2);
        }
        default:
            return Complex(0.0, 0.0);
    }
}

Complex ACCircuit::getTotalPower() const {
    if (!m_solved) {
        throw std::runtime_error("Circuit not solved");
    }
    
    Complex totalPower(0.0, 0.0);
    for (size_t i = 0; i < m_elements.size(); ++i) {
        totalPower += getElementPower(i);
    }
    
    return totalPower;
}

Complex ACCircuit::getImpedance(int node1, int node2) const {
    if (!m_solved || node1 < 0 || node1 >= m_nodeCount || node2 < 0 || node2 >= m_nodeCount) {
        throw std::out_of_range("Invalid node indices or circuit not solved");
    }
    
    // Apply a test voltage and calculate the resulting current
    Complex testVoltage(1.0, 0.0);
    Complex voltage1 = m_nodeVoltages[node1];
    Complex voltage2 = m_nodeVoltages[node2];
    Complex voltageDiff = voltage1 - voltage2;
    
    // Find a path from node1 to node2
    // This is a simplified approach and may not work for all circuits
    for (size_t i = 0; i < m_elements.size(); ++i) {
        const auto& element = m_elements[i];
        auto [elementNode1, elementNode2] = element.getNodes();
        
        if ((elementNode1 == node1 && elementNode2 == node2) || 
            (elementNode1 == node2 && elementNode2 == node1)) {
            // Direct path found
            return element.getImpedance(m_frequency);
        }
    }
    
    // No direct path found, return a high impedance
    return Complex(1e9, 0.0);
}

double ACCircuit::getFrequency() const {
    return m_frequency;
}

void ACCircuit::setFrequency(double frequency) {
    m_frequency = frequency;
    m_solved = false;
}

int ACCircuit::getNodeCount() const {
    return m_nodeCount;
}

size_t ACCircuit::getElementCount() const {
    return m_elements.size();
}

//
// DigitalCircuit implementation
//

DigitalCircuit::DigitalCircuit() {
}

size_t DigitalCircuit::addInput(const std::string& name, bool value) {
    m_inputs.push_back({name, value});
    return m_inputs.size() - 1;
}

size_t DigitalCircuit::addGate(GateType type, const std::vector<size_t>& inputs) {
    m_gates.push_back({type, inputs, false});
    return m_gates.size() - 1;
}

size_t DigitalCircuit::addOutput(const std::string& name, size_t input) {
    m_outputs.push_back({name, input});
    return m_outputs.size() - 1;
}

void DigitalCircuit::setInput(size_t inputIndex, bool value) {
    if (inputIndex >= m_inputs.size()) {
        throw std::out_of_range("Invalid input index");
    }
    
    m_inputs[inputIndex].value = value;
}

bool DigitalCircuit::getOutput(size_t outputIndex) const {
    if (outputIndex >= m_outputs.size()) {
        throw std::out_of_range("Invalid output index");
    }
    
    size_t input = m_outputs[outputIndex].input;
    
    if (input < m_inputs.size()) {
        return m_inputs[input].value;
    } else if (input - m_inputs.size() < m_gates.size()) {
        return m_gates[input - m_inputs.size()].output;
    } else {
        throw std::out_of_range("Invalid output input index");
    }
}

bool DigitalCircuit::getOutput(const std::string& name) const {
    for (size_t i = 0; i < m_outputs.size(); ++i) {
        if (m_outputs[i].name == name) {
            return getOutput(i);
        }
    }
    
    throw std::out_of_range("Output not found");
}

void DigitalCircuit::simulate() {
    // Evaluate all gates
    for (auto& gate : m_gates) {
        gate.output = evaluateGate(gate);
    }
}

std::vector<std::pair<std::vector<bool>, std::vector<bool>>> DigitalCircuit::getTruthTable() const {
    std::vector<std::pair<std::vector<bool>, std::vector<bool>>> truthTable;
    
    // Generate all possible input combinations
    size_t numCombinations = 1 << m_inputs.size();
    
    for (size_t i = 0; i < numCombinations; ++i) {
        // Set inputs
        std::vector<bool> inputs(m_inputs.size());
        for (size_t j = 0; j < m_inputs.size(); ++j) {
            inputs[j] = (i >> j) & 1;
        }
        
        // Create a copy of the circuit and set the inputs
        DigitalCircuit circuit = *this;
        for (size_t j = 0; j < m_inputs.size(); ++j) {
            circuit.setInput(j, inputs[j]);
        }
        
        // Simulate the circuit
        circuit.simulate();
        
        // Get outputs
        std::vector<bool> outputs(m_outputs.size());
        for (size_t j = 0; j < m_outputs.size(); ++j) {
            outputs[j] = circuit.getOutput(j);
        }
        
        truthTable.emplace_back(inputs, outputs);
    }
    
    return truthTable;
}

size_t DigitalCircuit::getInputCount() const {
    return m_inputs.size();
}

size_t DigitalCircuit::getGateCount() const {
    return m_gates.size();
}

size_t DigitalCircuit::getOutputCount() const {
    return m_outputs.size();
}

bool DigitalCircuit::evaluateGate(const Gate& gate) const {
    std::vector<bool> inputValues;
    
    for (size_t input : gate.inputs) {
        if (input < m_inputs.size()) {
            inputValues.push_back(m_inputs[input].value);
        } else if (input - m_inputs.size() < m_gates.size()) {
            inputValues.push_back(m_gates[input - m_inputs.size()].output);
        } else {
            throw std::out_of_range("Invalid gate input index");
        }
    }
    
    switch (gate.type) {
        case GateType::AND: {
            if (inputValues.empty()) {
                return false;
            }
            return std::all_of(inputValues.begin(), inputValues.end(), [](bool v) { return v; });
        }
        case GateType::OR: {
            if (inputValues.empty()) {
                return false;
            }
            return std::any_of(inputValues.begin(), inputValues.end(), [](bool v) { return v; });
        }
        case GateType::NOT: {
            if (inputValues.empty()) {
                return true;
            }
            return !inputValues[0];
        }
        case GateType::NAND: {
            if (inputValues.empty()) {
                return true;
            }
            return !std::all_of(inputValues.begin(), inputValues.end(), [](bool v) { return v; });
        }
        case GateType::NOR: {
            if (inputValues.empty()) {
                return true;
            }
            return !std::any_of(inputValues.begin(), inputValues.end(), [](bool v) { return v; });
        }
        case GateType::XOR: {
            if (inputValues.empty()) {
                return false;
            }
            return std::count(inputValues.begin(), inputValues.end(), true) % 2 == 1;
        }
        case GateType::XNOR: {
            if (inputValues.empty()) {
                return true;
            }
            return std::count(inputValues.begin(), inputValues.end(), true) % 2 == 0;
        }
        case GateType::BUFFER: {
            if (inputValues.empty()) {
                return false;
            }
            return inputValues[0];
        }
        default:
            return false;
    }
}

//
// PowerSystem implementation
//

PowerSystem::PowerSystem(double baseVoltage, double basePower)
    : m_baseVoltage(baseVoltage), m_basePower(basePower), m_solved(false) {
}

size_t PowerSystem::addBus(const std::string& name, int type, double voltage, double angle, double activePower, double reactivePower) {
    m_buses.push_back({name, type, voltage, angle, activePower, reactivePower});
    m_solved = false;
    return m_buses.size() - 1;
}

size_t PowerSystem::addLine(size_t fromBus, size_t toBus, double resistance, double reactance, double susceptance) {
    if (fromBus >= m_buses.size() || toBus >= m_buses.size()) {
        throw std::out_of_range("Invalid bus index");
    }
    
    m_elements.push_back({fromBus, toBus, resistance, reactance, susceptance, 1.0, false});
    m_solved = false;
    return m_elements.size() - 1;
}

size_t PowerSystem::addTransformer(size_t fromBus, size_t toBus, double resistance, double reactance, double ratio) {
    if (fromBus >= m_buses.size() || toBus >= m_buses.size()) {
        throw std::out_of_range("Invalid bus index");
    }
    
    m_elements.push_back({fromBus, toBus, resistance, reactance, 0.0, ratio, true});
    m_solved = false;
    return m_elements.size() - 1;
}

bool PowerSystem::solvePowerFlow(int maxIterations, double tolerance) {
    // This is a simplified implementation of power flow analysis
    // A complete implementation would use the Newton-Raphson method
    
    // For now, we'll just return false to indicate that the power flow was not solved
    m_solved = false;
    return false;
}

double PowerSystem::getBusVoltageMagnitude(size_t busIndex) const {
    if (!m_solved || busIndex >= m_buses.size()) {
        throw std::out_of_range("Invalid bus index or power flow not solved");
    }
    
    return m_buses[busIndex].voltage;
}

double PowerSystem::getBusVoltageAngle(size_t busIndex) const {
    if (!m_solved || busIndex >= m_buses.size()) {
        throw std::out_of_range("Invalid bus index or power flow not solved");
    }
    
    return m_buses[busIndex].angle;
}

double PowerSystem::getBusActivePower(size_t busIndex) const {
    if (!m_solved || busIndex >= m_buses.size()) {
        throw std::out_of_range("Invalid bus index or power flow not solved");
    }
    
    return m_buses[busIndex].activePower;
}

double PowerSystem::getBusReactivePower(size_t busIndex) const {
    if (!m_solved || busIndex >= m_buses.size()) {
        throw std::out_of_range("Invalid bus index or power flow not solved");
    }
    
    return m_buses[busIndex].reactivePower;
}

double PowerSystem::getElementActivePower(size_t elementIndex, bool fromBus) const {
    if (!m_solved || elementIndex >= m_elements.size()) {
        throw std::out_of_range("Invalid element index or power flow not solved");
    }
    
    // This is a simplified implementation
    // A complete implementation would calculate the power flow based on the solved voltages
    
    return 0.0;
}

double PowerSystem::getElementReactivePower(size_t elementIndex, bool fromBus) const {
    if (!m_solved || elementIndex >= m_elements.size()) {
        throw std::out_of_range("Invalid element index or power flow not solved");
    }
    
    // This is a simplified implementation
    // A complete implementation would calculate the power flow based on the solved voltages
    
    return 0.0;
}

double PowerSystem::getElementCurrentMagnitude(size_t elementIndex, bool fromBus) const {
    if (!m_solved || elementIndex >= m_elements.size()) {
        throw std::out_of_range("Invalid element index or power flow not solved");
    }
    
    // This is a simplified implementation
    // A complete implementation would calculate the current based on the solved voltages
    
    return 0.0;
}

double PowerSystem::getElementPowerLoss(size_t elementIndex) const {
    if (!m_solved || elementIndex >= m_elements.size()) {
        throw std::out_of_range("Invalid element index or power flow not solved");
    }
    
    // This is a simplified implementation
    // A complete implementation would calculate the power loss based on the solved voltages
    
    return 0.0;
}

double PowerSystem::getTotalPowerLoss() const {
    if (!m_solved) {
        throw std::runtime_error("Power flow not solved");
    }
    
    double totalLoss = 0.0;
    for (size_t i = 0; i < m_elements.size(); ++i) {
        totalLoss += getElementPowerLoss(i);
    }
    
    return totalLoss;
}

size_t PowerSystem::getBusCount() const {
    return m_buses.size();
}

size_t PowerSystem::getElementCount() const {
    return m_elements.size();
}

//
// Utility functions
//

double resistance(double resistivity, double length, double area) {
    return resistivity * length / area;
}

double resistivityTemperature(double resistivity0, double alpha, double temperature, double temperature0) {
    return resistivity0 * (1.0 + alpha * (temperature - temperature0));
}

double capacitanceParallelPlate(double area, double distance, double permittivity) {
    return permittivity * area / distance;
}

double capacitanceCylindrical(double length, double innerRadius, double outerRadius, double permittivity) {
    return 2.0 * M_PI * permittivity * length / std::log(outerRadius / innerRadius);
}

double inductanceSolenoid(int turns, double area, double length, double permeability) {
    return permeability * turns * turns * area / length;
}

double inductanceToroid(int turns, double innerRadius, double outerRadius, double permeability) {
    return permeability * turns * turns * std::log(outerRadius / innerRadius) / (2.0 * M_PI);
}

Complex impedanceResistor(double resistance) {
    return Complex(resistance, 0.0);
}

Complex impedanceCapacitor(double capacitance, double frequency) {
    if (frequency < 1e-10) {
        return Complex(0.0, -1e10); // Very high impedance at DC
    }
    
    double reactance = 1.0 / (2.0 * M_PI * frequency * capacitance);
    return Complex(0.0, -reactance);
}

Complex impedanceInductor(double inductance, double frequency) {
    double reactance = 2.0 * M_PI * frequency * inductance;
    return Complex(0.0, reactance);
}

Complex impedanceSeriesRLC(double resistance, double inductance, double capacitance, double frequency) {
    return impedanceResistor(resistance) + impedanceInductor(inductance, frequency) + impedanceCapacitor(capacitance, frequency);
}

Complex impedanceParallelRLC(double resistance, double inductance, double capacitance, double frequency) {
    Complex zR = impedanceResistor(resistance);
    Complex zL = impedanceInductor(inductance, frequency);
    Complex zC = impedanceCapacitor(capacitance, frequency);
    
    return 1.0 / (1.0 / zR + 1.0 / zL + 1.0 / zC);
}

double resonantFrequency(double inductance, double capacitance) {
    return 1.0 / (2.0 * M_PI * std::sqrt(inductance * capacitance));
}

double qualityFactorSeries(double resistance, double inductance, double capacitance) {
    double omega0 = 1.0 / std::sqrt(inductance * capacitance);
    return omega0 * inductance / resistance;
}

double qualityFactorParallel(double resistance, double inductance, double capacitance) {
    double omega0 = 1.0 / std::sqrt(inductance * capacitance);
    return resistance / (omega0 * inductance);
}

double bandwidth(double resonantFrequency, double qualityFactor) {
    return resonantFrequency / qualityFactor;
}

double cutoffFrequencyLowPassRC(double resistance, double capacitance) {
    return 1.0 / (2.0 * M_PI * resistance * capacitance);
}

double cutoffFrequencyHighPassRC(double resistance, double capacitance) {
    return 1.0 / (2.0 * M_PI * resistance * capacitance);
}

double cutoffFrequencyLowPassRL(double resistance, double inductance) {
    return resistance / (2.0 * M_PI * inductance);
}

double cutoffFrequencyHighPassRL(double resistance, double inductance) {
    return resistance / (2.0 * M_PI * inductance);
}

Complex transferFunctionLowPassRC(double resistance, double capacitance, double frequency) {
    Complex impedance = impedanceResistor(resistance) + impedanceCapacitor(capacitance, frequency);
    return impedanceCapacitor(capacitance, frequency) / impedance;
}

Complex transferFunctionHighPassRC(double resistance, double capacitance, double frequency) {
    Complex impedance = impedanceResistor(resistance) + impedanceCapacitor(capacitance, frequency);
    return impedanceResistor(resistance) / impedance;
}

Complex transferFunctionLowPassRL(double resistance, double inductance, double frequency) {
    Complex impedance = impedanceResistor(resistance) + impedanceInductor(inductance, frequency);
    return impedanceResistor(resistance) / impedance;
}

Complex transferFunctionHighPassRL(double resistance, double inductance, double frequency) {
    Complex impedance = impedanceResistor(resistance) + impedanceInductor(inductance, frequency);
    return impedanceInductor(inductance, frequency) / impedance;
}

Complex transferFunctionBandPassRLC(double resistance, double inductance, double capacitance, double frequency) {
    Complex impedance = impedanceResistor(resistance) + impedanceInductor(inductance, frequency) + impedanceCapacitor(capacitance, frequency);
    return impedanceResistor(resistance) / impedance;
}

Complex transferFunctionBandStopRLC(double resistance, double inductance, double capacitance, double frequency) {
    Complex impedanceParallel = 1.0 / (1.0 / impedanceInductor(inductance, frequency) + 1.0 / impedanceCapacitor(capacitance, frequency));
    Complex impedance = impedanceResistor(resistance) + impedanceParallel;
    return impedanceParallel / impedance;
}

double gainDecibels(double inputPower, double outputPower) {
    return 10.0 * std::log10(outputPower / inputPower);
}

double gainDecibelsVoltage(double inputVoltage, double outputVoltage) {
    return 20.0 * std::log10(outputVoltage / inputVoltage);
}

double gainDecibelsCurrent(double inputCurrent, double outputCurrent) {
    return 20.0 * std::log10(outputCurrent / inputCurrent);
}

double powerDC(double voltage, double current) {
    return voltage * current;
}

double powerAC(double voltage, double current, double powerFactor) {
    return voltage * current * powerFactor;
}

Complex powerComplex(const Complex& voltage, const Complex& current) {
    return voltage * std::conj(current);
}

double activePower(const Complex& power) {
    return power.real();
}

double reactivePower(const Complex& power) {
    return power.imag();
}

double apparentPower(const Complex& power) {
    return std::abs(power);
}

double powerFactor(const Complex& power) {
    return power.real() / std::abs(power);
}

double powerFactorAngle(const Complex& power) {
    return std::atan2(power.imag(), power.real());
}

double energyDC(double power, double time) {
    return power * time;
}

double energyAC(double power, double time) {
    return power * time;
}

double efficiency(double outputPower, double inputPower) {
    return outputPower / inputPower;
}

double voltageDivider(double inputVoltage, double resistance1, double resistance2) {
    return inputVoltage * resistance2 / (resistance1 + resistance2);
}

double currentDivider(double inputCurrent, double resistance1, double resistance2) {
    return inputCurrent * resistance1 / (resistance1 + resistance2);
}

double theveninVoltage(double openCircuitVoltage) {
    return openCircuitVoltage;
}

double theveninResistance(double openCircuitVoltage, double shortCircuitCurrent) {
    return openCircuitVoltage / shortCircuitCurrent;
}

double nortonCurrent(double shortCircuitCurrent) {
    return shortCircuitCurrent;
}

double nortonResistance(double openCircuitVoltage, double shortCircuitCurrent) {
    return openCircuitVoltage / shortCircuitCurrent;
}

double maximumPowerTransfer(double sourceVoltage, double sourceResistance) {
    return sourceVoltage * sourceVoltage / (4.0 * sourceResistance);
}

double loadResistanceMaximumPowerTransfer(double sourceResistance) {
    return sourceResistance;
}

} // namespace Electrical
} // namespace Engineering
} // namespace RebelCalc
