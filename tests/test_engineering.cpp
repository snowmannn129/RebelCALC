#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "../src/engineering/cad.h"
#include "../src/engineering/physics.h"
#include "../src/engineering/electrical.h"

using namespace RebelCalc::Engineering;

// Helper function to print a Point3D
void printPoint(const CAD::Point3D& point, const std::string& name) {
    std::cout << name << ": (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
}

// Helper function to print a Complex number
void printComplex(const Electrical::Complex& complex, const std::string& name) {
    std::cout << name << ": " << complex.real() << " + " << complex.imag() << "i" << std::endl;
}

// Test CAD module
void testCAD() {
    std::cout << "\n=== Testing CAD Module ===\n" << std::endl;
    
    // Create some points
    CAD::Point3D p1(1.0, 2.0, 3.0);
    CAD::Point3D p2(4.0, 5.0, 6.0);
    CAD::Point3D p3(7.0, 8.0, 9.0);
    
    printPoint(p1, "Point 1");
    printPoint(p2, "Point 2");
    printPoint(p3, "Point 3");
    
    // Test point operations
    CAD::Point3D sum = p1 + p2;
    CAD::Point3D diff = p2 - p1;
    CAD::Point3D scaled = p1 * 2.0;
    
    printPoint(sum, "Point 1 + Point 2");
    printPoint(diff, "Point 2 - Point 1");
    printPoint(scaled, "Point 1 * 2.0");
    
    // Test dot and cross product
    double dot = p1.dot(p2);
    CAD::Point3D cross = p1.cross(p2);
    
    std::cout << "Dot product: " << dot << std::endl;
    printPoint(cross, "Cross product");
    
    // Test distance
    double distance = p1.distance(p2);
    std::cout << "Distance between Point 1 and Point 2: " << distance << std::endl;
    
    // Create a triangle
    CAD::Triangle triangle(p1, p2, p3);
    
    // Calculate area and normal
    double area = triangle.area();
    CAD::Point3D normal = triangle.normal();
    
    std::cout << "Triangle area: " << area << std::endl;
    printPoint(normal, "Triangle normal");
    
    // Test if a point is inside the triangle
    CAD::Point3D testPoint = (p1 + p2 + p3) / 3.0; // Centroid should be inside
    bool inside = triangle.containsPoint(testPoint);
    
    printPoint(testPoint, "Test point");
    std::cout << "Is test point inside triangle? " << (inside ? "Yes" : "No") << std::endl;
    
    // Create a circle
    CAD::Point3D center(0.0, 0.0, 0.0);
    CAD::Point3D circleNormal(0.0, 0.0, 1.0);
    double radius = 5.0;
    CAD::Circle circle(center, circleNormal, radius);
    
    // Calculate circle properties
    double circleArea = circle.area();
    double circumference = circle.circumference();
    
    std::cout << "Circle area: " << circleArea << std::endl;
    std::cout << "Circle circumference: " << circumference << std::endl;
    
    // Test if a point is on the circle
    CAD::Point3D circlePoint = circle.pointAtAngle(M_PI / 4.0);
    bool onCircle = circle.containsPoint(circlePoint);
    
    printPoint(circlePoint, "Point on circle");
    std::cout << "Is point on circle? " << (onCircle ? "Yes" : "No") << std::endl;
}

// Test Physics module
void testPhysics() {
    std::cout << "\n=== Testing Physics Module ===\n" << std::endl;
    
    // Test projectile motion
    CAD::Point3D initialPosition(0.0, 0.0, 0.0);
    CAD::Point3D initialVelocity(10.0, 15.0, 0.0);
    CAD::Point3D gravity(0.0, -9.81, 0.0);
    
    std::cout << "Projectile Motion:" << std::endl;
    std::cout << "Initial position: (" << initialPosition.x << ", " << initialPosition.y << ", " << initialPosition.z << ")" << std::endl;
    std::cout << "Initial velocity: (" << initialVelocity.x << ", " << initialVelocity.y << ", " << initialVelocity.z << ")" << std::endl;
    
    // Calculate trajectory at different times
    for (double t = 0.0; t <= 3.0; t += 0.5) {
        CAD::Point3D position = Physics::projectileTrajectory(initialPosition, initialVelocity, gravity, t);
        std::cout << "Position at t = " << t << "s: (" << position.x << ", " << position.y << ", " << position.z << ")" << std::endl;
    }
    
    // Calculate time of flight, range, and maximum height
    double timeOfFlight = Physics::projectileTimeOfFlight(initialPosition.y, initialVelocity.y);
    double range = Physics::projectileRange(initialPosition.y, initialVelocity);
    double maxHeight = Physics::projectileMaxHeight(initialPosition.y, initialVelocity.y);
    
    std::cout << "Time of flight: " << timeOfFlight << " s" << std::endl;
    std::cout << "Range: " << range << " m" << std::endl;
    std::cout << "Maximum height: " << maxHeight << " m" << std::endl;
    
    // Test energy calculations
    double mass = 1.0; // 1 kg
    double kineticEnergy = Physics::kineticEnergy(mass, initialVelocity);
    double potentialEnergy = Physics::potentialEnergy(mass, maxHeight);
    
    std::cout << "\nEnergy Calculations:" << std::endl;
    std::cout << "Initial kinetic energy: " << kineticEnergy << " J" << std::endl;
    std::cout << "Maximum potential energy: " << potentialEnergy << " J" << std::endl;
    
    // Test spring force
    CAD::Point3D springPos1(0.0, 0.0, 0.0);
    CAD::Point3D springPos2(1.0, 0.0, 0.0);
    double springConstant = 10.0; // 10 N/m
    double restLength = 0.5; // 0.5 m
    
    CAD::Point3D springForce = Physics::springForce(springConstant, restLength, springPos1, springPos2);
    
    std::cout << "\nSpring Force:" << std::endl;
    std::cout << "Spring constant: " << springConstant << " N/m" << std::endl;
    std::cout << "Rest length: " << restLength << " m" << std::endl;
    std::cout << "Current length: " << springPos1.distance(springPos2) << " m" << std::endl;
    std::cout << "Spring force: (" << springForce.x << ", " << springForce.y << ", " << springForce.z << ") N" << std::endl;
    
    // Test gravitational force
    CAD::Point3D pos1(0.0, 0.0, 0.0);
    CAD::Point3D pos2(1.0, 0.0, 0.0);
    double mass1 = 1.0e6; // 1,000,000 kg
    double mass2 = 1.0e6; // 1,000,000 kg
    
    CAD::Point3D gravForce = Physics::gravitationalForce(mass1, mass2, pos1, pos2);
    
    std::cout << "\nGravitational Force:" << std::endl;
    std::cout << "Mass 1: " << mass1 << " kg" << std::endl;
    std::cout << "Mass 2: " << mass2 << " kg" << std::endl;
    std::cout << "Distance: " << pos1.distance(pos2) << " m" << std::endl;
    std::cout << "Gravitational force: (" << gravForce.x << ", " << gravForce.y << ", " << gravForce.z << ") N" << std::endl;
}

// Test Electrical module
void testElectrical() {
    std::cout << "\n=== Testing Electrical Module ===\n" << std::endl;
    
    // Test impedance calculations
    double resistance = 100.0; // 100 ohms
    double capacitance = 1.0e-6; // 1 uF
    double inductance = 0.1; // 0.1 H
    double frequency = 1000.0; // 1 kHz
    
    Electrical::Complex zR = Electrical::impedanceResistor(resistance);
    Electrical::Complex zC = Electrical::impedanceCapacitor(capacitance, frequency);
    Electrical::Complex zL = Electrical::impedanceInductor(inductance, frequency);
    
    std::cout << "Impedance Calculations at " << frequency << " Hz:" << std::endl;
    printComplex(zR, "Resistor impedance");
    printComplex(zC, "Capacitor impedance");
    printComplex(zL, "Inductor impedance");
    
    // Test series and parallel RLC circuits
    Electrical::Complex zSeries = Electrical::impedanceSeriesRLC(resistance, inductance, capacitance, frequency);
    Electrical::Complex zParallel = Electrical::impedanceParallelRLC(resistance, inductance, capacitance, frequency);
    
    printComplex(zSeries, "Series RLC impedance");
    printComplex(zParallel, "Parallel RLC impedance");
    
    // Test resonant frequency
    double resonantFreq = Electrical::resonantFrequency(inductance, capacitance);
    double qFactorSeries = Electrical::qualityFactorSeries(resistance, inductance, capacitance);
    double bandwidth = Electrical::bandwidth(resonantFreq, qFactorSeries);
    
    std::cout << "\nResonant Circuit Properties:" << std::endl;
    std::cout << "Resonant frequency: " << resonantFreq << " Hz" << std::endl;
    std::cout << "Quality factor (series): " << qFactorSeries << std::endl;
    std::cout << "Bandwidth: " << bandwidth << " Hz" << std::endl;
    
    // Test filter calculations
    double cutoffFreqLP = Electrical::cutoffFrequencyLowPassRC(resistance, capacitance);
    double cutoffFreqHP = Electrical::cutoffFrequencyHighPassRC(resistance, capacitance);
    
    std::cout << "\nFilter Calculations:" << std::endl;
    std::cout << "Low-pass RC cutoff frequency: " << cutoffFreqLP << " Hz" << std::endl;
    std::cout << "High-pass RC cutoff frequency: " << cutoffFreqHP << " Hz" << std::endl;
    
    // Test transfer functions
    Electrical::Complex tfLP = Electrical::transferFunctionLowPassRC(resistance, capacitance, frequency);
    Electrical::Complex tfHP = Electrical::transferFunctionHighPassRC(resistance, capacitance, frequency);
    
    printComplex(tfLP, "Low-pass RC transfer function at " + std::to_string(frequency) + " Hz");
    printComplex(tfHP, "High-pass RC transfer function at " + std::to_string(frequency) + " Hz");
    
    // Test power calculations
    double voltage = 10.0; // 10 V
    double current = 0.1; // 0.1 A
    double powerFactor = 0.8; // 0.8 lagging
    
    double powerDC = Electrical::powerDC(voltage, current);
    double powerAC = Electrical::powerAC(voltage, current, powerFactor);
    
    std::cout << "\nPower Calculations:" << std::endl;
    std::cout << "DC power: " << powerDC << " W" << std::endl;
    std::cout << "AC power: " << powerAC << " W" << std::endl;
    
    // Test voltage divider
    double r1 = 1000.0; // 1 kOhm
    double r2 = 2000.0; // 2 kOhm
    double inputVoltage = 12.0; // 12 V
    
    double outputVoltage = Electrical::voltageDivider(inputVoltage, r1, r2);
    
    std::cout << "\nVoltage Divider:" << std::endl;
    std::cout << "Input voltage: " << inputVoltage << " V" << std::endl;
    std::cout << "R1: " << r1 << " Ohm" << std::endl;
    std::cout << "R2: " << r2 << " Ohm" << std::endl;
    std::cout << "Output voltage: " << outputVoltage << " V" << std::endl;
}

int main() {
    std::cout << "RebelCALC Engineering Modules Test" << std::endl;
    std::cout << "=================================" << std::endl;
    
    testCAD();
    testPhysics();
    testElectrical();
    
    return 0;
}
