#ifndef REBELCALC_ENGINEERING_PHYSICS_H
#define REBELCALC_ENGINEERING_PHYSICS_H

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <stdexcept>
#include <functional>
#include "../backend/matrix.h"
#include "cad.h" // For Point3D

namespace RebelCalc {
namespace Engineering {
namespace Physics {

/**
 * @brief A class representing a rigid body for physics calculations
 */
class RigidBody {
public:
    /**
     * @brief Constructor for a rigid body
     * @param mass Mass of the body
     * @param position Initial position
     * @param velocity Initial velocity
     * @param orientation Initial orientation as a quaternion (w, x, y, z)
     * @param angularVelocity Initial angular velocity
     * @param inertiaTensor Inertia tensor of the body
     */
    RigidBody(double mass, 
              const CAD::Point3D& position, 
              const CAD::Point3D& velocity,
              const std::array<double, 4>& orientation,
              const CAD::Point3D& angularVelocity,
              const rebelcalc::Matrix& inertiaTensor);
    
    /**
     * @brief Apply a force to the body at a specific point
     * @param force Force vector
     * @param point Point of application (in world coordinates)
     */
    void applyForce(const CAD::Point3D& force, const CAD::Point3D& point);
    
    /**
     * @brief Apply a torque to the body
     * @param torque Torque vector
     */
    void applyTorque(const CAD::Point3D& torque);
    
    /**
     * @brief Update the body's state for a time step
     * @param dt Time step
     */
    void update(double dt);
    
    /**
     * @brief Get the body's position
     * @return Position vector
     */
    const CAD::Point3D& getPosition() const;
    
    /**
     * @brief Get the body's velocity
     * @return Velocity vector
     */
    const CAD::Point3D& getVelocity() const;
    
    /**
     * @brief Get the body's orientation as a quaternion
     * @return Orientation quaternion (w, x, y, z)
     */
    const std::array<double, 4>& getOrientation() const;
    
    /**
     * @brief Get the body's angular velocity
     * @return Angular velocity vector
     */
    const CAD::Point3D& getAngularVelocity() const;
    
    /**
     * @brief Get the body's mass
     * @return Mass
     */
    double getMass() const;
    
    /**
     * @brief Get the body's inertia tensor
     * @return Inertia tensor
     */
    const rebelcalc::Matrix& getInertiaTensor() const;
    
    /**
     * @brief Get the body's transformation matrix
     * @return 4x4 transformation matrix
     */
    rebelcalc::Matrix getTransformationMatrix() const;
    
private:
    double m_mass;                    // Mass of the body
    CAD::Point3D m_position;          // Position of the center of mass
    CAD::Point3D m_velocity;          // Linear velocity
    std::array<double, 4> m_orientation; // Orientation as a quaternion (w, x, y, z)
    CAD::Point3D m_angularVelocity;   // Angular velocity
    rebelcalc::Matrix m_inertiaTensor; // Inertia tensor (3x3 matrix)
    
    CAD::Point3D m_force;             // Accumulated force
    CAD::Point3D m_torque;            // Accumulated torque
    
    // Helper methods
    void normalizeQuaternion();
    rebelcalc::Matrix quaternionToRotationMatrix() const;
    std::array<double, 4> quaternionMultiply(const std::array<double, 4>& q1, const std::array<double, 4>& q2) const;
};

/**
 * @brief A class for simulating particle systems
 */
class ParticleSystem {
public:
    /**
     * @brief Constructor for a particle system
     */
    ParticleSystem();
    
    /**
     * @brief Add a particle to the system
     * @param position Initial position
     * @param velocity Initial velocity
     * @param mass Mass of the particle
     * @param radius Radius of the particle
     * @return Index of the added particle
     */
    size_t addParticle(const CAD::Point3D& position, const CAD::Point3D& velocity, double mass, double radius);
    
    /**
     * @brief Apply a force to a specific particle
     * @param index Index of the particle
     * @param force Force vector
     */
    void applyForce(size_t index, const CAD::Point3D& force);
    
    /**
     * @brief Apply a force to all particles
     * @param force Force vector
     */
    void applyForceToAll(const CAD::Point3D& force);
    
    /**
     * @brief Apply a force field to all particles
     * @param forceField Function that takes a position and returns a force vector
     */
    void applyForceField(std::function<CAD::Point3D(const CAD::Point3D&)> forceField);
    
    /**
     * @brief Update the particle system for a time step
     * @param dt Time step
     */
    void update(double dt);
    
    /**
     * @brief Get the number of particles in the system
     * @return Number of particles
     */
    size_t getParticleCount() const;
    
    /**
     * @brief Get the position of a particle
     * @param index Index of the particle
     * @return Position vector
     */
    const CAD::Point3D& getParticlePosition(size_t index) const;
    
    /**
     * @brief Get the velocity of a particle
     * @param index Index of the particle
     * @return Velocity vector
     */
    const CAD::Point3D& getParticleVelocity(size_t index) const;
    
    /**
     * @brief Get the mass of a particle
     * @param index Index of the particle
     * @return Mass
     */
    double getParticleMass(size_t index) const;
    
    /**
     * @brief Get the radius of a particle
     * @param index Index of the particle
     * @return Radius
     */
    double getParticleRadius(size_t index) const;
    
    /**
     * @brief Set the collision detection and response function
     * @param collisionFunc Function that handles collisions between particles
     */
    void setCollisionFunction(std::function<void(size_t, size_t, ParticleSystem&)> collisionFunc);
    
    /**
     * @brief Set the boundary condition function
     * @param boundaryFunc Function that handles boundary conditions for particles
     */
    void setBoundaryFunction(std::function<void(size_t, ParticleSystem&)> boundaryFunc);
    
private:
    struct Particle {
        CAD::Point3D position;
        CAD::Point3D velocity;
        CAD::Point3D force;
        double mass;
        double radius;
    };
    
    std::vector<Particle> m_particles;
    std::function<void(size_t, size_t, ParticleSystem&)> m_collisionFunc;
    std::function<void(size_t, ParticleSystem&)> m_boundaryFunc;
    
    // Helper methods
    void detectAndResolveCollisions();
    void applyBoundaryConditions();
};

/**
 * @brief A class for simulating springs and spring systems
 */
class SpringSystem {
public:
    /**
     * @brief Constructor for a spring system
     */
    SpringSystem();
    
    /**
     * @brief Add a particle to the system
     * @param position Initial position
     * @param velocity Initial velocity
     * @param mass Mass of the particle
     * @param fixed Whether the particle is fixed in place
     * @return Index of the added particle
     */
    size_t addParticle(const CAD::Point3D& position, const CAD::Point3D& velocity, double mass, bool fixed = false);
    
    /**
     * @brief Add a spring between two particles
     * @param index1 Index of the first particle
     * @param index2 Index of the second particle
     * @param springConstant Spring constant (stiffness)
     * @param dampingFactor Damping factor
     * @param restLength Rest length of the spring (if 0, uses current distance)
     * @return Index of the added spring
     */
    size_t addSpring(size_t index1, size_t index2, double springConstant, double dampingFactor, double restLength = 0.0);
    
    /**
     * @brief Apply a force to a specific particle
     * @param index Index of the particle
     * @param force Force vector
     */
    void applyForce(size_t index, const CAD::Point3D& force);
    
    /**
     * @brief Apply a force to all particles
     * @param force Force vector
     */
    void applyForceToAll(const CAD::Point3D& force);
    
    /**
     * @brief Update the spring system for a time step
     * @param dt Time step
     */
    void update(double dt);
    
    /**
     * @brief Get the number of particles in the system
     * @return Number of particles
     */
    size_t getParticleCount() const;
    
    /**
     * @brief Get the number of springs in the system
     * @return Number of springs
     */
    size_t getSpringCount() const;
    
    /**
     * @brief Get the position of a particle
     * @param index Index of the particle
     * @return Position vector
     */
    const CAD::Point3D& getParticlePosition(size_t index) const;
    
    /**
     * @brief Get the velocity of a particle
     * @param index Index of the particle
     * @return Velocity vector
     */
    const CAD::Point3D& getParticleVelocity(size_t index) const;
    
    /**
     * @brief Get the indices of the particles connected by a spring
     * @param index Index of the spring
     * @return Pair of particle indices
     */
    std::pair<size_t, size_t> getSpringParticles(size_t index) const;
    
    /**
     * @brief Get the current length of a spring
     * @param index Index of the spring
     * @return Current length
     */
    double getSpringLength(size_t index) const;
    
    /**
     * @brief Get the rest length of a spring
     * @param index Index of the spring
     * @return Rest length
     */
    double getSpringRestLength(size_t index) const;
    
    /**
     * @brief Get the spring constant of a spring
     * @param index Index of the spring
     * @return Spring constant
     */
    double getSpringConstant(size_t index) const;
    
private:
    struct Particle {
        CAD::Point3D position;
        CAD::Point3D velocity;
        CAD::Point3D force;
        double mass;
        bool fixed;
    };
    
    struct Spring {
        size_t index1;
        size_t index2;
        double springConstant;
        double dampingFactor;
        double restLength;
    };
    
    std::vector<Particle> m_particles;
    std::vector<Spring> m_springs;
    
    // Helper methods
    void calculateSpringForces();
};

/**
 * @brief A class for simulating fluid dynamics using Smoothed Particle Hydrodynamics (SPH)
 */
class FluidSystem {
public:
    /**
     * @brief Constructor for a fluid system
     * @param smoothingLength Smoothing length for SPH
     * @param particleMass Mass of each fluid particle
     * @param restDensity Rest density of the fluid
     * @param gasConstant Gas constant for pressure calculation
     * @param viscosity Viscosity coefficient
     * @param surfaceTension Surface tension coefficient
     */
    FluidSystem(double smoothingLength, double particleMass, double restDensity,
                double gasConstant, double viscosity, double surfaceTension);
    
    /**
     * @brief Add a fluid particle to the system
     * @param position Initial position
     * @param velocity Initial velocity
     * @return Index of the added particle
     */
    size_t addParticle(const CAD::Point3D& position, const CAD::Point3D& velocity);
    
    /**
     * @brief Update the fluid system for a time step
     * @param dt Time step
     */
    void update(double dt);
    
    /**
     * @brief Get the number of particles in the system
     * @return Number of particles
     */
    size_t getParticleCount() const;
    
    /**
     * @brief Get the position of a particle
     * @param index Index of the particle
     * @return Position vector
     */
    const CAD::Point3D& getParticlePosition(size_t index) const;
    
    /**
     * @brief Get the velocity of a particle
     * @param index Index of the particle
     * @return Velocity vector
     */
    const CAD::Point3D& getParticleVelocity(size_t index) const;
    
    /**
     * @brief Get the density of a particle
     * @param index Index of the particle
     * @return Density
     */
    double getParticleDensity(size_t index) const;
    
    /**
     * @brief Get the pressure of a particle
     * @param index Index of the particle
     * @return Pressure
     */
    double getParticlePressure(size_t index) const;
    
    /**
     * @brief Set the boundary condition function
     * @param boundaryFunc Function that handles boundary conditions for particles
     */
    void setBoundaryFunction(std::function<void(size_t, FluidSystem&)> boundaryFunc);
    
    /**
     * @brief Set the external force function
     * @param forceFunc Function that applies external forces to particles
     */
    void setExternalForceFunction(std::function<CAD::Point3D(const CAD::Point3D&)> forceFunc);
    
private:
    struct Particle {
        CAD::Point3D position;
        CAD::Point3D velocity;
        CAD::Point3D force;
        double density;
        double pressure;
    };
    
    double m_smoothingLength;
    double m_particleMass;
    double m_restDensity;
    double m_gasConstant;
    double m_viscosity;
    double m_surfaceTension;
    
    std::vector<Particle> m_particles;
    std::function<void(size_t, FluidSystem&)> m_boundaryFunc;
    std::function<CAD::Point3D(const CAD::Point3D&)> m_externalForceFunc;
    
    // Helper methods
    void findNeighbors(std::vector<std::vector<size_t>>& neighbors);
    void calculateDensityAndPressure(const std::vector<std::vector<size_t>>& neighbors);
    void calculateForces(const std::vector<std::vector<size_t>>& neighbors);
    void applyBoundaryConditions();
    
    // SPH kernel functions
    double kernelPoly6(double r, double h);
    CAD::Point3D kernelPoly6Gradient(const CAD::Point3D& r, double h);
    double kernelPoly6Laplacian(double r, double h);
    CAD::Point3D kernelSpikyGradient(const CAD::Point3D& r, double h);
    double kernelViscosityLaplacian(double r, double h);
};

/**
 * @brief Calculate the gravitational force between two masses
 * 
 * @param mass1 First mass
 * @param mass2 Second mass
 * @param position1 Position of the first mass
 * @param position2 Position of the second mass
 * @param G Gravitational constant (default: 6.67430e-11)
 * @return Gravitational force vector (acting on mass1)
 */
CAD::Point3D gravitationalForce(double mass1, double mass2, const CAD::Point3D& position1, 
                               const CAD::Point3D& position2, double G = 6.67430e-11);

/**
 * @brief Calculate the electric force between two charges
 * 
 * @param charge1 First charge
 * @param charge2 Second charge
 * @param position1 Position of the first charge
 * @param position2 Position of the second charge
 * @param k Coulomb constant (default: 8.9875517923e9)
 * @return Electric force vector (acting on charge1)
 */
CAD::Point3D electricForce(double charge1, double charge2, const CAD::Point3D& position1, 
                          const CAD::Point3D& position2, double k = 8.9875517923e9);

/**
 * @brief Calculate the magnetic force on a moving charge
 * 
 * @param charge Charge
 * @param velocity Velocity of the charge
 * @param magneticField Magnetic field vector
 * @return Magnetic force vector
 */
CAD::Point3D magneticForce(double charge, const CAD::Point3D& velocity, const CAD::Point3D& magneticField);

/**
 * @brief Calculate the drag force on an object
 * 
 * @param density Fluid density
 * @param velocity Velocity of the object
 * @param area Cross-sectional area
 * @param dragCoefficient Drag coefficient
 * @return Drag force vector
 */
CAD::Point3D dragForce(double density, const CAD::Point3D& velocity, double area, double dragCoefficient);

/**
 * @brief Calculate the spring force
 * 
 * @param springConstant Spring constant
 * @param restLength Rest length of the spring
 * @param position1 Position of the first end of the spring
 * @param position2 Position of the second end of the spring
 * @return Spring force vector (acting on position1)
 */
CAD::Point3D springForce(double springConstant, double restLength, const CAD::Point3D& position1, 
                        const CAD::Point3D& position2);

/**
 * @brief Calculate the damping force
 * 
 * @param dampingCoefficient Damping coefficient
 * @param velocity Velocity
 * @return Damping force vector
 */
CAD::Point3D dampingForce(double dampingCoefficient, const CAD::Point3D& velocity);

/**
 * @brief Calculate the projectile trajectory
 * 
 * @param initialPosition Initial position
 * @param initialVelocity Initial velocity
 * @param acceleration Acceleration (e.g., gravity)
 * @param time Time
 * @return Position at the given time
 */
CAD::Point3D projectileTrajectory(const CAD::Point3D& initialPosition, const CAD::Point3D& initialVelocity, 
                                 const CAD::Point3D& acceleration, double time);

/**
 * @brief Calculate the time of flight for a projectile
 * 
 * @param initialHeight Initial height
 * @param initialVelocity Initial velocity
 * @param gravity Gravity (default: 9.81)
 * @return Time of flight
 */
double projectileTimeOfFlight(double initialHeight, double initialVelocityY, double gravity = 9.81);

/**
 * @brief Calculate the range of a projectile
 * 
 * @param initialHeight Initial height
 * @param initialVelocity Initial velocity
 * @param gravity Gravity (default: 9.81)
 * @return Range of the projectile
 */
double projectileRange(double initialHeight, const CAD::Point3D& initialVelocity, double gravity = 9.81);

/**
 * @brief Calculate the maximum height of a projectile
 * 
 * @param initialHeight Initial height
 * @param initialVelocityY Initial vertical velocity
 * @param gravity Gravity (default: 9.81)
 * @return Maximum height
 */
double projectileMaxHeight(double initialHeight, double initialVelocityY, double gravity = 9.81);

/**
 * @brief Calculate the kinetic energy of an object
 * 
 * @param mass Mass
 * @param velocity Velocity
 * @return Kinetic energy
 */
double kineticEnergy(double mass, const CAD::Point3D& velocity);

/**
 * @brief Calculate the potential energy of an object in a gravitational field
 * 
 * @param mass Mass
 * @param height Height
 * @param gravity Gravity (default: 9.81)
 * @return Potential energy
 */
double potentialEnergy(double mass, double height, double gravity = 9.81);

/**
 * @brief Calculate the elastic potential energy of a spring
 * 
 * @param springConstant Spring constant
 * @param displacement Displacement from rest length
 * @return Elastic potential energy
 */
double elasticPotentialEnergy(double springConstant, double displacement);

/**
 * @brief Calculate the momentum of an object
 * 
 * @param mass Mass
 * @param velocity Velocity
 * @return Momentum vector
 */
CAD::Point3D momentum(double mass, const CAD::Point3D& velocity);

/**
 * @brief Calculate the angular momentum of an object
 * 
 * @param position Position vector (relative to the axis of rotation)
 * @param momentum Linear momentum vector
 * @return Angular momentum vector
 */
CAD::Point3D angularMomentum(const CAD::Point3D& position, const CAD::Point3D& momentum);

/**
 * @brief Calculate the moment of inertia of a point mass
 * 
 * @param mass Mass
 * @param radius Distance from the axis of rotation
 * @return Moment of inertia
 */
double momentOfInertia(double mass, double radius);

/**
 * @brief Calculate the rotational kinetic energy
 * 
 * @param momentOfInertia Moment of inertia
 * @param angularVelocity Angular velocity
 * @return Rotational kinetic energy
 */
double rotationalKineticEnergy(double momentOfInertia, double angularVelocity);

/**
 * @brief Calculate the centripetal force
 * 
 * @param mass Mass
 * @param velocity Tangential velocity
 * @param radius Radius of the circular path
 * @return Centripetal force
 */
double centripetalForce(double mass, double velocity, double radius);

/**
 * @brief Calculate the centrifugal force
 * 
 * @param mass Mass
 * @param velocity Tangential velocity
 * @param radius Radius of the circular path
 * @return Centrifugal force
 */
double centrifugalForce(double mass, double velocity, double radius);

/**
 * @brief Calculate the Coriolis force
 * 
 * @param mass Mass
 * @param velocity Velocity in the rotating frame
 * @param angularVelocity Angular velocity of the rotating frame
 * @return Coriolis force vector
 */
CAD::Point3D coriolisForce(double mass, const CAD::Point3D& velocity, const CAD::Point3D& angularVelocity);

/**
 * @brief Calculate the work done by a force
 * 
 * @param force Force vector
 * @param displacement Displacement vector
 * @return Work done
 */
double work(const CAD::Point3D& force, const CAD::Point3D& displacement);

/**
 * @brief Calculate the power
 * 
 * @param force Force vector
 * @param velocity Velocity vector
 * @return Power
 */
double power(const CAD::Point3D& force, const CAD::Point3D& velocity);

/**
 * @brief Calculate the impulse
 * 
 * @param force Force vector
 * @param time Time
 * @return Impulse vector
 */
CAD::Point3D impulse(const CAD::Point3D& force, double time);

/**
 * @brief Calculate the change in momentum from an impulse
 * 
 * @param impulse Impulse vector
 * @return Change in momentum vector
 */
CAD::Point3D momentumChange(const CAD::Point3D& impulse);

/**
 * @brief Calculate the coefficient of restitution from initial and final velocities
 * 
 * @param initialVelocity1 Initial velocity of the first object
 * @param initialVelocity2 Initial velocity of the second object
 * @param finalVelocity1 Final velocity of the first object
 * @param finalVelocity2 Final velocity of the second object
 * @return Coefficient of restitution
 */
double coefficientOfRestitution(double initialVelocity1, double initialVelocity2, 
                               double finalVelocity1, double finalVelocity2);

/**
 * @brief Calculate the final velocities after a collision
 * 
 * @param mass1 Mass of the first object
 * @param mass2 Mass of the second object
 * @param initialVelocity1 Initial velocity of the first object
 * @param initialVelocity2 Initial velocity of the second object
 * @param coefficientOfRestitution Coefficient of restitution
 * @return Pair of final velocities (first, second)
 */
std::pair<double, double> collisionVelocities(double mass1, double mass2, 
                                             double initialVelocity1, double initialVelocity2, 
                                             double coefficientOfRestitution);

/**
 * @brief Calculate the pressure in a fluid
 * 
 * @param density Fluid density
 * @param height Height
 * @param gravity Gravity (default: 9.81)
 * @return Pressure
 */
double fluidPressure(double density, double height, double gravity = 9.81);

/**
 * @brief Calculate the buoyant force
 * 
 * @param fluidDensity Fluid density
 * @param volume Submerged volume
 * @param gravity Gravity (default: 9.81)
 * @return Buoyant force
 */
double buoyantForce(double fluidDensity, double volume, double gravity = 9.81);

/**
 * @brief Calculate the Reynolds number
 * 
 * @param density Fluid density
 * @param velocity Flow velocity
 * @param length Characteristic length
 * @param viscosity Dynamic viscosity
 * @return Reynolds number
 */
double reynoldsNumber(double density, double velocity, double length, double viscosity);

/**
 * @brief Calculate the Bernoulli's equation
 * 
 * @param pressure1 Pressure at point 1
 * @param density Fluid density
 * @param velocity1 Velocity at point 1
 * @param height1 Height at point 1
 * @param pressure2 Pressure at point 2
 * @param velocity2 Velocity at point 2
 * @param height2 Height at point 2
 * @param gravity Gravity (default: 9.81)
 * @return true if Bernoulli's equation is satisfied, false otherwise
 */
bool bernoulliEquation(double pressure1, double density, double velocity1, double height1,
                      double pressure2, double velocity2, double height2, double gravity = 9.81);

/**
 * @brief Calculate the flow rate
 * 
 * @param area Cross-sectional area
 * @param velocity Flow velocity
 * @return Flow rate
 */
double flowRate(double area, double velocity);

/**
 * @brief Calculate the terminal velocity
 * 
 * @param mass Mass
 * @param area Cross-sectional area
 * @param dragCoefficient Drag coefficient
 * @param fluidDensity Fluid density
 * @param gravity Gravity (default: 9.81)
 * @return Terminal velocity
 */
double terminalVelocity(double mass, double area, double dragCoefficient, 
                       double fluidDensity, double gravity = 9.81);

} // namespace Physics
} // namespace Engineering
} // namespace RebelCalc

#endif // REBELCALC_ENGINEERING_PHYSICS_H
