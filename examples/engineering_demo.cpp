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

// Helper function to print a section header
void printHeader(const std::string& header) {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "  " << header << std::endl;
    std::cout << std::string(80, '=') << std::endl;
}

// Demo for CAD module
void demoCAD() {
    printHeader("CAD Module Demonstration");
    
    // Create some points
    CAD::Point3D p1(1.0, 2.0, 3.0);
    CAD::Point3D p2(4.0, 5.0, 6.0);
    CAD::Point3D p3(7.0, 8.0, 9.0);
    
    std::cout << "Created three points:" << std::endl;
    printPoint(p1, "Point 1");
    printPoint(p2, "Point 2");
    printPoint(p3, "Point 3");
    
    // Demonstrate point operations
    std::cout << "\nPoint operations:" << std::endl;
    printPoint(p1 + p2, "Point 1 + Point 2");
    printPoint(p2 - p1, "Point 2 - Point 1");
    printPoint(p1 * 2.0, "Point 1 * 2.0");
    printPoint(p1 / 2.0, "Point 1 / 2.0");
    
    // Demonstrate vector operations
    std::cout << "\nVector operations:" << std::endl;
    std::cout << "Dot product of Point 1 and Point 2: " << p1.dot(p2) << std::endl;
    printPoint(p1.cross(p2), "Cross product of Point 1 and Point 2");
    std::cout << "Magnitude of Point 1: " << p1.magnitude() << std::endl;
    printPoint(p1.normalize(), "Normalized Point 1");
    
    // Demonstrate distance calculation
    std::cout << "\nDistance between Point 1 and Point 2: " << p1.distance(p2) << std::endl;
    
    // Create and demonstrate a line segment
    CAD::LineSegment line(p1, p2);
    std::cout << "\nLine segment operations:" << std::endl;
    std::cout << "Line length: " << line.length() << std::endl;
    printPoint(line.midpoint(), "Line midpoint");
    printPoint(line.direction(), "Line direction (normalized)");
    printPoint(line.pointAt(0.25), "Point at 25% along the line");
    
    // Create and demonstrate a plane
    CAD::Plane plane(p1, p1.cross(p2).normalize());
    std::cout << "\nPlane operations:" << std::endl;
    std::cout << "Distance from Point 3 to plane: " << plane.distanceToPoint(p3) << std::endl;
    printPoint(plane.projectPoint(p3), "Projection of Point 3 onto plane");
    
    // Create and demonstrate a circle
    CAD::Point3D center(0.0, 0.0, 0.0);
    CAD::Point3D normal(0.0, 0.0, 1.0);
    double radius = 5.0;
    CAD::Circle circle(center, normal, radius);
    
    std::cout << "\nCircle operations:" << std::endl;
    std::cout << "Circle area: " << circle.area() << std::endl;
    std::cout << "Circle circumference: " << circle.circumference() << std::endl;
    
    // Create points on the circle at different angles
    std::cout << "\nPoints on the circle:" << std::endl;
    for (int i = 0; i < 8; ++i) {
        double angle = i * M_PI / 4.0;
        CAD::Point3D point = circle.pointAtAngle(angle);
        printPoint(point, "Point at angle " + std::to_string(angle) + " radians");
    }
    
    // Create and demonstrate a triangle
    CAD::Triangle triangle(p1, p2, p3);
    std::cout << "\nTriangle operations:" << std::endl;
    std::cout << "Triangle area: " << triangle.area() << std::endl;
    std::cout << "Triangle perimeter: " << triangle.perimeter() << std::endl;
    printPoint(triangle.normal(), "Triangle normal");
    
    // Test if a point is inside the triangle
    CAD::Point3D testPoint = (p1 + p2 + p3) / 3.0; // Centroid should be inside
    std::cout << "\nTesting if a point is inside the triangle:" << std::endl;
    printPoint(testPoint, "Test point (centroid)");
    std::cout << "Is the test point inside the triangle? " << (triangle.containsPoint(testPoint) ? "Yes" : "No") << std::endl;
    
    // Demonstrate intersection detection
    CAD::LineSegment testLine(center, testPoint);
    CAD::Point3D intersection;
    bool intersects = CAD::lineTriangleIntersection(testLine, triangle, intersection);
    
    std::cout << "\nIntersection detection:" << std::endl;
    std::cout << "Does the line intersect the triangle? " << (intersects ? "Yes" : "No") << std::endl;
    if (intersects) {
        printPoint(intersection, "Intersection point");
    }
    
    // Demonstrate bounding box calculation
    std::vector<CAD::Point3D> points = {p1, p2, p3, center};
    std::array<CAD::Point3D, 2> bbox = CAD::boundingBox(points);
    
    std::cout << "\nBounding box of points:" << std::endl;
    printPoint(bbox[0], "Minimum corner");
    printPoint(bbox[1], "Maximum corner");
}

// Demo for Physics module
void demoPhysics() {
    printHeader("Physics Module Demonstration");
    
    // Demonstrate projectile motion
    std::cout << "Projectile Motion Simulation:" << std::endl;
    CAD::Point3D initialPosition(0.0, 0.0, 0.0);
    CAD::Point3D initialVelocity(20.0, 15.0, 0.0);
    CAD::Point3D gravity(0.0, -9.81, 0.0);
    
    std::cout << "Initial position: (" << initialPosition.x << ", " << initialPosition.y << ", " << initialPosition.z << ")" << std::endl;
    std::cout << "Initial velocity: (" << initialVelocity.x << ", " << initialVelocity.y << ", " << initialVelocity.z << ")" << std::endl;
    std::cout << "Gravity: (" << gravity.x << ", " << gravity.y << ", " << gravity.z << ")" << std::endl;
    
    // Calculate and display trajectory at different times
    std::cout << "\nTrajectory:" << std::endl;
    std::cout << std::setw(10) << "Time (s)" << std::setw(15) << "X Position (m)" << std::setw(15) << "Y Position (m)" << std::setw(15) << "Z Position (m)" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    
    double timeStep = 0.5;
    double maxTime = 3.5;
    
    for (double t = 0.0; t <= maxTime; t += timeStep) {
        CAD::Point3D position = Physics::projectileTrajectory(initialPosition, initialVelocity, gravity, t);
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::setw(10) << t << std::setw(15) << position.x << std::setw(15) << position.y << std::setw(15) << position.z << std::endl;
    }
    
    // Calculate and display key projectile parameters
    double timeOfFlight = Physics::projectileTimeOfFlight(initialPosition.y, initialVelocity.y);
    double range = Physics::projectileRange(initialPosition.y, initialVelocity);
    double maxHeight = Physics::projectileMaxHeight(initialPosition.y, initialVelocity.y);
    
    std::cout << "\nProjectile Parameters:" << std::endl;
    std::cout << "Time of flight: " << timeOfFlight << " s" << std::endl;
    std::cout << "Range: " << range << " m" << std::endl;
    std::cout << "Maximum height: " << maxHeight << " m" << std::endl;
    
    // Demonstrate energy calculations
    double mass = 1.0; // 1 kg
    double kineticEnergy = Physics::kineticEnergy(mass, initialVelocity);
    double potentialEnergy = Physics::potentialEnergy(mass, maxHeight);
    
    std::cout << "\nEnergy Calculations:" << std::endl;
    std::cout << "Initial kinetic energy: " << kineticEnergy << " J" << std::endl;
    std::cout << "Maximum potential energy: " << potentialEnergy << " J" << std::endl;
    std::cout << "Total energy: " << kineticEnergy + potentialEnergy << " J" << std::endl;
    
    // Demonstrate spring force calculation
    CAD::Point3D springPos1(0.0, 0.0, 0.0);
    CAD::Point3D springPos2(1.5, 0.0, 0.0);
    double springConstant = 10.0; // 10 N/m
    double restLength = 1.0; // 1.0 m
    
    CAD::Point3D springForce = Physics::springForce(springConstant, restLength, springPos1, springPos2);
    
    std::cout << "\nSpring Force Calculation:" << std::endl;
    std::cout << "Spring constant: " << springConstant << " N/m" << std::endl;
    std::cout << "Rest length: " << restLength << " m" << std::endl;
    std::cout << "Current length: " << springPos1.distance(springPos2) << " m" << std::endl;
    std::cout << "Spring force: (" << springForce.x << ", " << springForce.y << ", " << springForce.z << ") N" << std::endl;
    std::cout << "Spring force magnitude: " << springForce.magnitude() << " N" << std::endl;
    
    // Demonstrate gravitational force calculation
    CAD::Point3D pos1(0.0, 0.0, 0.0);
    CAD::Point3D pos2(10.0, 0.0, 0.0);
    double mass1 = 5.97e24; // Earth's mass in kg
    double mass2 = 7.35e22; // Moon's mass in kg
    
    CAD::Point3D gravForce = Physics::gravitationalForce(mass1, mass2, pos1, pos2);
    
    std::cout << "\nGravitational Force Calculation:" << std::endl;
    std::cout << "Mass 1: " << mass1 << " kg" << std::endl;
    std::cout << "Mass 2: " << mass2 << " kg" << std::endl;
    std::cout << "Distance: " << pos1.distance(pos2) << " m" << std::endl;
    std::cout << "Gravitational force: (" << gravForce.x << ", " << gravForce.y << ", " << gravForce.z << ") N" << std::endl;
    std::cout << "Gravitational force magnitude: " << gravForce.magnitude() << " N" << std::endl;
    
    // Demonstrate particle system
    std::cout << "\nParticle System Simulation:" << std::endl;
    Physics::ParticleSystem particleSystem;
    
    // Add particles
    particleSystem.addParticle(CAD::Point3D(0.0, 0.0, 0.0), CAD::Point3D(1.0, 0.0, 0.0), 1.0, 0.5);
    particleSystem.addParticle(CAD::Point3D(5.0, 0.0, 0.0), CAD::Point3D(-1.0, 0.0, 0.0), 1.0, 0.5);
    
    std::cout << "Created a particle system with " << particleSystem.getParticleCount() << " particles" << std::endl;
    
    // Apply gravity to all particles
    particleSystem.applyForceToAll(CAD::Point3D(0.0, -9.81, 0.0));
    
    // Simulate for a few time steps
    std::cout << "\nParticle positions over time:" << std::endl;
    std::cout << std::setw(10) << "Time (s)" << std::setw(20) << "Particle 1 Position" << std::setw(20) << "Particle 2 Position" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    double dt = 0.1;
    for (int i = 0; i <= 5; ++i) {
        double t = i * dt;
        
        // Print current positions
        CAD::Point3D p1 = particleSystem.getParticlePosition(0);
        CAD::Point3D p2 = particleSystem.getParticlePosition(1);
        
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::setw(10) << t 
                  << std::setw(20) << "(" << p1.x << ", " << p1.y << ", " << p1.z << ")"
                  << std::setw(20) << "(" << p2.x << ", " << p2.y << ", " << p2.z << ")" << std::endl;
        
        // Update the system
        particleSystem.update(dt);
    }
}

// Demo for Electrical module
void demoElectrical() {
    printHeader("Electrical Module Demonstration");
    
    // Demonstrate impedance calculations
    std::cout << "Impedance Calculations:" << std::endl;
    
    double resistance = 100.0; // 100 ohms
    double capacitance = 1.0e-6; // 1 uF
    double inductance = 0.1; // 0.1 H
    
    std::cout << "Component values:" << std::endl;
    std::cout << "Resistance: " << resistance << " Ω" << std::endl;
    std::cout << "Capacitance: " << capacitance * 1e6 << " μF" << std::endl;
    std::cout << "Inductance: " << inductance * 1e3 << " mH" << std::endl;
    
    // Calculate impedance at different frequencies
    std::cout << "\nImpedance vs. Frequency:" << std::endl;
    std::cout << std::setw(15) << "Frequency (Hz)" 
              << std::setw(20) << "Resistor (Ω)" 
              << std::setw(20) << "Capacitor (Ω)" 
              << std::setw(20) << "Inductor (Ω)" << std::endl;
    std::cout << std::string(75, '-') << std::endl;
    
    for (int i = 0; i <= 5; ++i) {
        double frequency = std::pow(10.0, i); // 1, 10, 100, 1000, 10000, 100000 Hz
        
        Electrical::Complex zR = Electrical::impedanceResistor(resistance);
        Electrical::Complex zC = Electrical::impedanceCapacitor(capacitance, frequency);
        Electrical::Complex zL = Electrical::impedanceInductor(inductance, frequency);
        
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::setw(15) << frequency 
                  << std::setw(20) << std::abs(zR) 
                  << std::setw(20) << std::abs(zC) 
                  << std::setw(20) << std::abs(zL) << std::endl;
    }
    
    // Demonstrate series and parallel RLC circuits
    std::cout << "\nSeries and Parallel RLC Circuits:" << std::endl;
    std::cout << std::setw(15) << "Frequency (Hz)" 
              << std::setw(20) << "Series RLC (Ω)" 
              << std::setw(20) << "Parallel RLC (Ω)" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    
    for (int i = 0; i <= 5; ++i) {
        double frequency = std::pow(10.0, i); // 1, 10, 100, 1000, 10000, 100000 Hz
        
        Electrical::Complex zSeries = Electrical::impedanceSeriesRLC(resistance, inductance, capacitance, frequency);
        Electrical::Complex zParallel = Electrical::impedanceParallelRLC(resistance, inductance, capacitance, frequency);
        
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::setw(15) << frequency 
                  << std::setw(20) << std::abs(zSeries) 
                  << std::setw(20) << std::abs(zParallel) << std::endl;
    }
    
    // Demonstrate resonant frequency calculation
    double resonantFreq = Electrical::resonantFrequency(inductance, capacitance);
    double qFactorSeries = Electrical::qualityFactorSeries(resistance, inductance, capacitance);
    double qFactorParallel = Electrical::qualityFactorParallel(resistance, inductance, capacitance);
    double bandwidth = Electrical::bandwidth(resonantFreq, qFactorSeries);
    
    std::cout << "\nResonant Circuit Properties:" << std::endl;
    std::cout << "Resonant frequency: " << resonantFreq << " Hz" << std::endl;
    std::cout << "Quality factor (series): " << qFactorSeries << std::endl;
    std::cout << "Quality factor (parallel): " << qFactorParallel << std::endl;
    std::cout << "Bandwidth: " << bandwidth << " Hz" << std::endl;
    
    // Demonstrate filter calculations
    std::cout << "\nFilter Calculations:" << std::endl;
    
    double cutoffFreqLP = Electrical::cutoffFrequencyLowPassRC(resistance, capacitance);
    double cutoffFreqHP = Electrical::cutoffFrequencyHighPassRC(resistance, capacitance);
    
    std::cout << "Low-pass RC cutoff frequency: " << cutoffFreqLP << " Hz" << std::endl;
    std::cout << "High-pass RC cutoff frequency: " << cutoffFreqHP << " Hz" << std::endl;
    
    // Demonstrate transfer function calculations
    std::cout << "\nTransfer Function vs. Frequency:" << std::endl;
    std::cout << std::setw(15) << "Frequency (Hz)" 
              << std::setw(20) << "Low-pass RC" 
              << std::setw(20) << "High-pass RC" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    
    for (int i = -1; i <= 4; ++i) {
        double frequency = cutoffFreqLP * std::pow(10.0, i); // 0.1fc, fc, 10fc, 100fc, 1000fc, 10000fc
        
        Electrical::Complex tfLP = Electrical::transferFunctionLowPassRC(resistance, capacitance, frequency);
        Electrical::Complex tfHP = Electrical::transferFunctionHighPassRC(resistance, capacitance, frequency);
        
        // Convert to dB
        double gainLP = 20.0 * std::log10(std::abs(tfLP));
        double gainHP = 20.0 * std::log10(std::abs(tfHP));
        
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::setw(15) << frequency 
                  << std::setw(20) << gainLP << " dB" 
                  << std::setw(20) << gainHP << " dB" << std::endl;
    }
    
    // Demonstrate voltage divider
    std::cout << "\nVoltage Divider:" << std::endl;
    
    double r1 = 1000.0; // 1 kOhm
    double r2 = 2000.0; // 2 kOhm
    double inputVoltage = 12.0; // 12 V
    
    double outputVoltage = Electrical::voltageDivider(inputVoltage, r1, r2);
    
    std::cout << "Input voltage: " << inputVoltage << " V" << std::endl;
    std::cout << "R1: " << r1 << " Ω" << std::endl;
    std::cout << "R2: " << r2 << " Ω" << std::endl;
    std::cout << "Output voltage: " << outputVoltage << " V" << std::endl;
    
    // Demonstrate power calculations
    std::cout << "\nPower Calculations:" << std::endl;
    
    double voltage = 120.0; // 120 V
    double current = 2.0; // 2 A
    double powerFactor = 0.8; // 0.8 lagging
    
    double powerDC = Electrical::powerDC(voltage, current);
    double powerAC = Electrical::powerAC(voltage, current, powerFactor);
    
    std::cout << "Voltage: " << voltage << " V" << std::endl;
    std::cout << "Current: " << current << " A" << std::endl;
    std::cout << "Power factor: " << powerFactor << std::endl;
    std::cout << "DC power: " << powerDC << " W" << std::endl;
    std::cout << "AC power: " << powerAC << " W" << std::endl;
    
    // Demonstrate Thevenin and Norton equivalents
    std::cout << "\nThevenin and Norton Equivalents:" << std::endl;
    
    double openCircuitVoltage = 10.0; // 10 V
    double shortCircuitCurrent = 0.1; // 0.1 A
    
    double theveninVoltage = Electrical::theveninVoltage(openCircuitVoltage);
    double theveninResistance = Electrical::theveninResistance(openCircuitVoltage, shortCircuitCurrent);
    double nortonCurrent = Electrical::nortonCurrent(shortCircuitCurrent);
    double nortonResistance = Electrical::nortonResistance(openCircuitVoltage, shortCircuitCurrent);
    
    std::cout << "Open-circuit voltage: " << openCircuitVoltage << " V" << std::endl;
    std::cout << "Short-circuit current: " << shortCircuitCurrent << " A" << std::endl;
    std::cout << "Thevenin voltage: " << theveninVoltage << " V" << std::endl;
    std::cout << "Thevenin resistance: " << theveninResistance << " Ω" << std::endl;
    std::cout << "Norton current: " << nortonCurrent << " A" << std::endl;
    std::cout << "Norton resistance: " << nortonResistance << " Ω" << std::endl;
    
    // Demonstrate maximum power transfer
    double maxPower = Electrical::maximumPowerTransfer(theveninVoltage, theveninResistance);
    double loadResistance = Electrical::loadResistanceMaximumPowerTransfer(theveninResistance);
    
    std::cout << "\nMaximum Power Transfer:" << std::endl;
    std::cout << "Maximum power: " << maxPower << " W" << std::endl;
    std::cout << "Load resistance for maximum power: " << loadResistance << " Ω" << std::endl;
}

int main() {
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "                      RebelCALC Engineering Modules Demo" << std::endl;
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "\nThis program demonstrates the capabilities of the RebelCALC engineering modules:" << std::endl;
    std::cout << "1. CAD Module - 3D geometry and spatial operations" << std::endl;
    std::cout << "2. Physics Module - Mechanics, dynamics, and physical simulations" << std::endl;
    std::cout << "3. Electrical Module - Circuit analysis and electrical engineering calculations" << std::endl;
    
    demoCAD();
    demoPhysics();
    demoElectrical();
    
    std::cout << "\n" << std::string(80, '*') << std::endl;
    std::cout << "                      End of Engineering Modules Demo" << std::endl;
    std::cout << std::string(80, '*') << std::endl;
    
    return 0;
}
