#include "physics.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include <cmath>

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace RebelCalc {
namespace Engineering {
namespace Physics {

//
// RigidBody implementation
//

RigidBody::RigidBody(double mass, 
                     const CAD::Point3D& position, 
                     const CAD::Point3D& velocity,
                     const std::array<double, 4>& orientation,
                     const CAD::Point3D& angularVelocity,
                     const rebelcalc::Matrix& inertiaTensor)
    : m_mass(mass), 
      m_position(position), 
      m_velocity(velocity),
      m_orientation(orientation),
      m_angularVelocity(angularVelocity),
      m_inertiaTensor(inertiaTensor),
      m_force(0, 0, 0),
      m_torque(0, 0, 0) {
    
    normalizeQuaternion();
}

void RigidBody::applyForce(const CAD::Point3D& force, const CAD::Point3D& point) {
    // Add the force to the accumulated force
    m_force = m_force + force;
    
    // Calculate the torque
    CAD::Point3D r = point - m_position;
    CAD::Point3D torque = r.cross(force);
    m_torque = m_torque + torque;
}

void RigidBody::applyTorque(const CAD::Point3D& torque) {
    m_torque = m_torque + torque;
}

void RigidBody::update(double dt) {
    // Update linear motion
    CAD::Point3D acceleration = m_force * (1.0 / m_mass);
    m_velocity = m_velocity + acceleration * dt;
    m_position = m_position + m_velocity * dt;
    
    // Update angular motion
    // Convert inertia tensor from local to world coordinates
    rebelcalc::Matrix R = quaternionToRotationMatrix();
    rebelcalc::Matrix worldInertiaTensor = R * m_inertiaTensor * R.transpose();
    
    // Calculate angular acceleration
    rebelcalc::Matrix inverseWorldInertiaTensor = worldInertiaTensor.inverse();
    rebelcalc::Matrix torqueMatrix(3, 1);
    torqueMatrix(0, 0) = m_torque.x;
    torqueMatrix(1, 0) = m_torque.y;
    torqueMatrix(2, 0) = m_torque.z;
    
    rebelcalc::Matrix angularAccelerationMatrix = inverseWorldInertiaTensor * torqueMatrix;
    CAD::Point3D angularAcceleration(angularAccelerationMatrix(0, 0), 
                                    angularAccelerationMatrix(1, 0), 
                                    angularAccelerationMatrix(2, 0));
    
    // Update angular velocity
    m_angularVelocity = m_angularVelocity + angularAcceleration * dt;
    
    // Update orientation
    double angle = m_angularVelocity.magnitude() * dt;
    if (angle > 1e-10) {
        CAD::Point3D axis = m_angularVelocity.normalize();
        double s = std::sin(angle / 2.0);
        std::array<double, 4> q = {
            std::cos(angle / 2.0),
            axis.x * s,
            axis.y * s,
            axis.z * s
        };
        m_orientation = quaternionMultiply(m_orientation, q);
        normalizeQuaternion();
    }
    
    // Reset accumulated forces and torques
    m_force = CAD::Point3D(0, 0, 0);
    m_torque = CAD::Point3D(0, 0, 0);
}

const CAD::Point3D& RigidBody::getPosition() const {
    return m_position;
}

const CAD::Point3D& RigidBody::getVelocity() const {
    return m_velocity;
}

const std::array<double, 4>& RigidBody::getOrientation() const {
    return m_orientation;
}

const CAD::Point3D& RigidBody::getAngularVelocity() const {
    return m_angularVelocity;
}

double RigidBody::getMass() const {
    return m_mass;
}

const rebelcalc::Matrix& RigidBody::getInertiaTensor() const {
    return m_inertiaTensor;
}

rebelcalc::Matrix RigidBody::getTransformationMatrix() const {
    rebelcalc::Matrix R = quaternionToRotationMatrix();
    rebelcalc::Matrix T(4, 4);
    
    // Set the rotation part
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T(i, j) = R(i, j);
        }
    }
    
    // Set the translation part
    T(0, 3) = m_position.x;
    T(1, 3) = m_position.y;
    T(2, 3) = m_position.z;
    
    // Set the bottom row
    T(3, 0) = 0.0;
    T(3, 1) = 0.0;
    T(3, 2) = 0.0;
    T(3, 3) = 1.0;
    
    return T;
}

void RigidBody::normalizeQuaternion() {
    double norm = std::sqrt(m_orientation[0] * m_orientation[0] +
                           m_orientation[1] * m_orientation[1] +
                           m_orientation[2] * m_orientation[2] +
                           m_orientation[3] * m_orientation[3]);
    
    if (norm < 1e-10) {
        m_orientation = {1.0, 0.0, 0.0, 0.0}; // Default to identity quaternion
    } else {
        m_orientation[0] /= norm;
        m_orientation[1] /= norm;
        m_orientation[2] /= norm;
        m_orientation[3] /= norm;
    }
}

rebelcalc::Matrix RigidBody::quaternionToRotationMatrix() const {
    double w = m_orientation[0];
    double x = m_orientation[1];
    double y = m_orientation[2];
    double z = m_orientation[3];
    
    rebelcalc::Matrix R(3, 3);
    
    R(0, 0) = 1.0 - 2.0 * (y * y + z * z);
    R(0, 1) = 2.0 * (x * y - w * z);
    R(0, 2) = 2.0 * (x * z + w * y);
    
    R(1, 0) = 2.0 * (x * y + w * z);
    R(1, 1) = 1.0 - 2.0 * (x * x + z * z);
    R(1, 2) = 2.0 * (y * z - w * x);
    
    R(2, 0) = 2.0 * (x * z - w * y);
    R(2, 1) = 2.0 * (y * z + w * x);
    R(2, 2) = 1.0 - 2.0 * (x * x + y * y);
    
    return R;
}

std::array<double, 4> RigidBody::quaternionMultiply(const std::array<double, 4>& q1, const std::array<double, 4>& q2) const {
    double w1 = q1[0], x1 = q1[1], y1 = q1[2], z1 = q1[3];
    double w2 = q2[0], x2 = q2[1], y2 = q2[2], z2 = q2[3];
    
    return {
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    };
}

//
// ParticleSystem implementation
//

ParticleSystem::ParticleSystem() {
    // Default collision function (elastic collision)
    m_collisionFunc = [](size_t i, size_t j, ParticleSystem& system) {
        Particle& p1 = system.m_particles[i];
        Particle& p2 = system.m_particles[j];
        
        CAD::Point3D r = p1.position - p2.position;
        double distance = r.magnitude();
        double minDistance = p1.radius + p2.radius;
        
        if (distance < minDistance) {
            // Collision detected
            CAD::Point3D normal = r.normalize();
            
            // Relative velocity
            CAD::Point3D relativeVelocity = p1.velocity - p2.velocity;
            
            // Relative velocity along the normal
            double vn = relativeVelocity.dot(normal);
            
            // Do not resolve if objects are moving away from each other
            if (vn > 0) {
                return;
            }
            
            // Coefficient of restitution (1.0 for elastic collision)
            double e = 1.0;
            
            // Impulse
            double j = -(1.0 + e) * vn / (1.0 / p1.mass + 1.0 / p2.mass);
            CAD::Point3D impulse = normal * j;
            
            // Apply impulse
            p1.velocity = p1.velocity + impulse * (1.0 / p1.mass);
            p2.velocity = p2.velocity - impulse * (1.0 / p2.mass);
            
            // Resolve penetration
            double penetration = minDistance - distance;
            CAD::Point3D correction = normal * (penetration * 0.5);
            p1.position = p1.position + correction;
            p2.position = p2.position - correction;
        }
    };
    
    // Default boundary function (no boundaries)
    m_boundaryFunc = [](size_t, ParticleSystem&) {};
}

size_t ParticleSystem::addParticle(const CAD::Point3D& position, const CAD::Point3D& velocity, double mass, double radius) {
    Particle particle;
    particle.position = position;
    particle.velocity = velocity;
    particle.force = CAD::Point3D(0, 0, 0);
    particle.mass = mass;
    particle.radius = radius;
    
    m_particles.push_back(particle);
    return m_particles.size() - 1;
}

void ParticleSystem::applyForce(size_t index, const CAD::Point3D& force) {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    m_particles[index].force = m_particles[index].force + force;
}

void ParticleSystem::applyForceToAll(const CAD::Point3D& force) {
    for (auto& particle : m_particles) {
        particle.force = particle.force + force;
    }
}

void ParticleSystem::applyForceField(std::function<CAD::Point3D(const CAD::Point3D&)> forceField) {
    for (auto& particle : m_particles) {
        particle.force = particle.force + forceField(particle.position);
    }
}

void ParticleSystem::update(double dt) {
    // Update positions and velocities
    for (auto& particle : m_particles) {
        CAD::Point3D acceleration = particle.force * (1.0 / particle.mass);
        particle.velocity = particle.velocity + acceleration * dt;
        particle.position = particle.position + particle.velocity * dt;
        particle.force = CAD::Point3D(0, 0, 0); // Reset force
    }
    
    // Detect and resolve collisions
    detectAndResolveCollisions();
    
    // Apply boundary conditions
    applyBoundaryConditions();
}

size_t ParticleSystem::getParticleCount() const {
    return m_particles.size();
}

const CAD::Point3D& ParticleSystem::getParticlePosition(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].position;
}

const CAD::Point3D& ParticleSystem::getParticleVelocity(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].velocity;
}

double ParticleSystem::getParticleMass(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].mass;
}

double ParticleSystem::getParticleRadius(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].radius;
}

void ParticleSystem::setCollisionFunction(std::function<void(size_t, size_t, ParticleSystem&)> collisionFunc) {
    m_collisionFunc = collisionFunc;
}

void ParticleSystem::setBoundaryFunction(std::function<void(size_t, ParticleSystem&)> boundaryFunc) {
    m_boundaryFunc = boundaryFunc;
}

void ParticleSystem::detectAndResolveCollisions() {
    for (size_t i = 0; i < m_particles.size(); ++i) {
        for (size_t j = i + 1; j < m_particles.size(); ++j) {
            m_collisionFunc(i, j, *this);
        }
    }
}

void ParticleSystem::applyBoundaryConditions() {
    for (size_t i = 0; i < m_particles.size(); ++i) {
        m_boundaryFunc(i, *this);
    }
}

//
// SpringSystem implementation
//

SpringSystem::SpringSystem() {
}

size_t SpringSystem::addParticle(const CAD::Point3D& position, const CAD::Point3D& velocity, double mass, bool fixed) {
    Particle particle;
    particle.position = position;
    particle.velocity = velocity;
    particle.force = CAD::Point3D(0, 0, 0);
    particle.mass = mass;
    particle.fixed = fixed;
    
    m_particles.push_back(particle);
    return m_particles.size() - 1;
}

size_t SpringSystem::addSpring(size_t index1, size_t index2, double springConstant, double dampingFactor, double restLength) {
    if (index1 >= m_particles.size() || index2 >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    Spring spring;
    spring.index1 = index1;
    spring.index2 = index2;
    spring.springConstant = springConstant;
    spring.dampingFactor = dampingFactor;
    
    if (restLength <= 0.0) {
        // Calculate rest length from current positions
        spring.restLength = m_particles[index1].position.distance(m_particles[index2].position);
    } else {
        spring.restLength = restLength;
    }
    
    m_springs.push_back(spring);
    return m_springs.size() - 1;
}

void SpringSystem::applyForce(size_t index, const CAD::Point3D& force) {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    m_particles[index].force = m_particles[index].force + force;
}

void SpringSystem::applyForceToAll(const CAD::Point3D& force) {
    for (auto& particle : m_particles) {
        if (!particle.fixed) {
            particle.force = particle.force + force;
        }
    }
}

void SpringSystem::update(double dt) {
    // Calculate spring forces
    calculateSpringForces();
    
    // Update positions and velocities
    for (auto& particle : m_particles) {
        if (!particle.fixed) {
            CAD::Point3D acceleration = particle.force * (1.0 / particle.mass);
            particle.velocity = particle.velocity + acceleration * dt;
            particle.position = particle.position + particle.velocity * dt;
        }
        
        particle.force = CAD::Point3D(0, 0, 0); // Reset force
    }
}

size_t SpringSystem::getParticleCount() const {
    return m_particles.size();
}

size_t SpringSystem::getSpringCount() const {
    return m_springs.size();
}

const CAD::Point3D& SpringSystem::getParticlePosition(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].position;
}

const CAD::Point3D& SpringSystem::getParticleVelocity(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].velocity;
}

std::pair<size_t, size_t> SpringSystem::getSpringParticles(size_t index) const {
    if (index >= m_springs.size()) {
        throw std::out_of_range("Spring index out of range");
    }
    
    return {m_springs[index].index1, m_springs[index].index2};
}

double SpringSystem::getSpringLength(size_t index) const {
    if (index >= m_springs.size()) {
        throw std::out_of_range("Spring index out of range");
    }
    
    const Spring& spring = m_springs[index];
    return m_particles[spring.index1].position.distance(m_particles[spring.index2].position);
}

double SpringSystem::getSpringRestLength(size_t index) const {
    if (index >= m_springs.size()) {
        throw std::out_of_range("Spring index out of range");
    }
    
    return m_springs[index].restLength;
}

double SpringSystem::getSpringConstant(size_t index) const {
    if (index >= m_springs.size()) {
        throw std::out_of_range("Spring index out of range");
    }
    
    return m_springs[index].springConstant;
}

void SpringSystem::calculateSpringForces() {
    for (const auto& spring : m_springs) {
        Particle& p1 = m_particles[spring.index1];
        Particle& p2 = m_particles[spring.index2];
        
        // Calculate spring force
        CAD::Point3D direction = p2.position - p1.position;
        double distance = direction.magnitude();
        
        if (distance < 1e-10) {
            continue; // Avoid division by zero
        }
        
        direction = direction * (1.0 / distance); // Normalize
        
        double displacement = distance - spring.restLength;
        double springForce = spring.springConstant * displacement;
        
        // Calculate damping force
        CAD::Point3D relativeVelocity = p2.velocity - p1.velocity;
        double dampingForce = spring.dampingFactor * relativeVelocity.dot(direction);
        
        // Total force
        CAD::Point3D force = direction * (springForce + dampingForce);
        
        // Apply forces
        if (!p1.fixed) {
            p1.force = p1.force + force;
        }
        
        if (!p2.fixed) {
            p2.force = p2.force - force;
        }
    }
}

//
// FluidSystem implementation
//

FluidSystem::FluidSystem(double smoothingLength, double particleMass, double restDensity,
                         double gasConstant, double viscosity, double surfaceTension)
    : m_smoothingLength(smoothingLength),
      m_particleMass(particleMass),
      m_restDensity(restDensity),
      m_gasConstant(gasConstant),
      m_viscosity(viscosity),
      m_surfaceTension(surfaceTension) {
    
    // Default boundary function (no boundaries)
    m_boundaryFunc = [](size_t, FluidSystem&) {};
    
    // Default external force function (no external forces)
    m_externalForceFunc = [](const CAD::Point3D&) { return CAD::Point3D(0, 0, 0); };
}

size_t FluidSystem::addParticle(const CAD::Point3D& position, const CAD::Point3D& velocity) {
    Particle particle;
    particle.position = position;
    particle.velocity = velocity;
    particle.force = CAD::Point3D(0, 0, 0);
    particle.density = m_restDensity;
    particle.pressure = 0.0;
    
    m_particles.push_back(particle);
    return m_particles.size() - 1;
}

void FluidSystem::update(double dt) {
    // Find neighbors for each particle
    std::vector<std::vector<size_t>> neighbors(m_particles.size());
    findNeighbors(neighbors);
    
    // Calculate density and pressure
    calculateDensityAndPressure(neighbors);
    
    // Calculate forces
    calculateForces(neighbors);
    
    // Apply boundary conditions
    applyBoundaryConditions();
    
    // Update positions and velocities
    for (auto& particle : m_particles) {
        CAD::Point3D acceleration = particle.force * (1.0 / m_particleMass);
        particle.velocity = particle.velocity + acceleration * dt;
        particle.position = particle.position + particle.velocity * dt;
        particle.force = CAD::Point3D(0, 0, 0); // Reset force
    }
}

size_t FluidSystem::getParticleCount() const {
    return m_particles.size();
}

const CAD::Point3D& FluidSystem::getParticlePosition(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].position;
}

const CAD::Point3D& FluidSystem::getParticleVelocity(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].velocity;
}

double FluidSystem::getParticleDensity(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].density;
}

double FluidSystem::getParticlePressure(size_t index) const {
    if (index >= m_particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    return m_particles[index].pressure;
}

void FluidSystem::setBoundaryFunction(std::function<void(size_t, FluidSystem&)> boundaryFunc) {
    m_boundaryFunc = boundaryFunc;
}

void FluidSystem::setExternalForceFunction(std::function<CAD::Point3D(const CAD::Point3D&)> forceFunc) {
    m_externalForceFunc = forceFunc;
}

void FluidSystem::findNeighbors(std::vector<std::vector<size_t>>& neighbors) {
    double h2 = m_smoothingLength * m_smoothingLength;
    
    for (size_t i = 0; i < m_particles.size(); ++i) {
        neighbors[i].clear();
        
        for (size_t j = 0; j < m_particles.size(); ++j) {
            if (i != j) {
                double distance2 = (m_particles[i].position - m_particles[j].position).dot(
                                   m_particles[i].position - m_particles[j].position);
                
                if (distance2 < h2) {
                    neighbors[i].push_back(j);
                }
            }
        }
    }
}

void FluidSystem::calculateDensityAndPressure(const std::vector<std::vector<size_t>>& neighbors) {
    for (size_t i = 0; i < m_particles.size(); ++i) {
        // Calculate density
        double density = 0.0;
        
        for (size_t j : neighbors[i]) {
            CAD::Point3D r = m_particles[i].position - m_particles[j].position;
            double r2 = r.dot(r);
            
            if (r2 < m_smoothingLength * m_smoothingLength) {
                density += m_particleMass * kernelPoly6(std::sqrt(r2), m_smoothingLength);
            }
        }
        
        // Add self-contribution
        density += m_particleMass * kernelPoly6(0.0, m_smoothingLength);
        
        m_particles[i].density = density;
        
        // Calculate pressure using equation of state
        m_particles[i].pressure = m_gasConstant * (density - m_restDensity);
    }
}

void FluidSystem::calculateForces(const std::vector<std::vector<size_t>>& neighbors) {
    for (size_t i = 0; i < m_particles.size(); ++i) {
        CAD::Point3D pressureForce(0, 0, 0);
        CAD::Point3D viscosityForce(0, 0, 0);
        CAD::Point3D surfaceTensionForce(0, 0, 0);
        
        for (size_t j : neighbors[i]) {
            CAD::Point3D r = m_particles[i].position - m_particles[j].position;
            double distance = r.magnitude();
            
            if (distance > 0.0 && distance < m_smoothingLength) {
                // Pressure force
                double pressureTerm = -m_particleMass * 
                                     (m_particles[i].pressure / (m_particles[i].density * m_particles[i].density) +
                                      m_particles[j].pressure / (m_particles[j].density * m_particles[j].density));
                
                pressureForce = pressureForce + kernelSpikyGradient(r, m_smoothingLength) * pressureTerm;
                
                // Viscosity force
                CAD::Point3D velocityDiff = m_particles[j].velocity - m_particles[i].velocity;
                viscosityForce = viscosityForce + velocityDiff * 
                                (m_viscosity * m_particleMass / m_particles[j].density * 
                                 kernelViscosityLaplacian(distance, m_smoothingLength));
                
                // Surface tension force (simplified)
                surfaceTensionForce = surfaceTensionForce + 
                                     kernelPoly6Gradient(r, m_smoothingLength) * 
                                     (m_surfaceTension * m_particleMass / m_particles[j].density);
            }
        }
        
        // External forces (e.g., gravity)
        CAD::Point3D externalForce = m_externalForceFunc(m_particles[i].position);
        
        // Total force
        m_particles[i].force = pressureForce + viscosityForce + surfaceTensionForce + externalForce;
    }
}

void FluidSystem::applyBoundaryConditions() {
    for (size_t i = 0; i < m_particles.size(); ++i) {
        m_boundaryFunc(i, *this);
    }
}

double FluidSystem::kernelPoly6(double r, double h) {
    if (r < 0.0 || r > h) {
        return 0.0;
    }
    
    double h2 = h * h;
    double h9 = h2 * h2 * h2 * h2 * h;
    double r2 = r * r;
    
    return 315.0 / (64.0 * M_PI * h9) * std::pow(h2 - r2, 3);
}

CAD::Point3D FluidSystem::kernelPoly6Gradient(const CAD::Point3D& r, double h) {
    double distance = r.magnitude();
    
    if (distance < 1e-10 || distance > h) {
        return CAD::Point3D(0, 0, 0);
    }
    
    double h2 = h * h;
    double h9 = h2 * h2 * h2 * h2 * h;
    double r2 = distance * distance;
    
    double factor = -945.0 / (32.0 * M_PI * h9) * std::pow(h2 - r2, 2);
    return r * factor;
}

double FluidSystem::kernelPoly6Laplacian(double r, double h) {
    if (r < 0.0 || r > h) {
        return 0.0;
    }
    
    double h2 = h * h;
    double h9 = h2 * h2 * h2 * h2 * h;
    double r2 = r * r;
    
    return -945.0 / (32.0 * M_PI * h9) * (h2 - r2) * (3.0 * h2 - 7.0 * r2);
}

CAD::Point3D FluidSystem::kernelSpikyGradient(const CAD::Point3D& r, double h) {
    double distance = r.magnitude();
    
    if (distance < 1e-10 || distance > h) {
        return CAD::Point3D(0, 0, 0);
    }
    
    double h6 = h * h * h * h * h * h;
    
    double factor = -45.0 / (M_PI * h6) * std::pow(h - distance, 2) / distance;
    return r * factor;
}

double FluidSystem::kernelViscosityLaplacian(double r, double h) {
    if (r < 0.0 || r > h) {
        return 0.0;
    }
    
    return 45.0 / (M_PI * std::pow(h, 6)) * (h - r);
}

// Utility functions

CAD::Point3D gravitationalForce(double mass1, double mass2, const CAD::Point3D& position1, 
                               const CAD::Point3D& position2, double G) {
    CAD::Point3D r = position2 - position1;
    double distance = r.magnitude();
    
    if (distance < 1e-10) {
        return CAD::Point3D(0, 0, 0); // Avoid division by zero
    }
    
    double forceMagnitude = G * mass1 * mass2 / (distance * distance);
    return r.normalize() * forceMagnitude;
}

CAD::Point3D electricForce(double charge1, double charge2, const CAD::Point3D& position1, 
                          const CAD::Point3D& position2, double k) {
    CAD::Point3D r = position2 - position1;
    double distance = r.magnitude();
    
    if (distance < 1e-10) {
        return CAD::Point3D(0, 0, 0); // Avoid division by zero
    }
    
    double forceMagnitude = k * charge1 * charge2 / (distance * distance);
    return r.normalize() * forceMagnitude;
}

CAD::Point3D magneticForce(double charge, const CAD::Point3D& velocity, const CAD::Point3D& magneticField) {
    return velocity.cross(magneticField) * charge;
}

CAD::Point3D dragForce(double density, const CAD::Point3D& velocity, double area, double dragCoefficient) {
    double velocityMagnitude = velocity.magnitude();
    
    if (velocityMagnitude < 1e-10) {
        return CAD::Point3D(0, 0, 0); // Avoid division by zero
    }
    
    double forceMagnitude = 0.5 * density * velocityMagnitude * velocityMagnitude * area * dragCoefficient;
    return velocity.normalize() * -forceMagnitude; // Negative because drag opposes motion
}

CAD::Point3D springForce(double springConstant, double restLength, const CAD::Point3D& position1, 
                        const CAD::Point3D& position2) {
    CAD::Point3D r = position2 - position1;
    double distance = r.magnitude();
    
    if (distance < 1e-10) {
        return CAD::Point3D(0, 0, 0); // Avoid division by zero
    }
    
    double displacement = distance - restLength;
    double forceMagnitude = springConstant * displacement;
    return r.normalize() * forceMagnitude;
}

CAD::Point3D dampingForce(double dampingCoefficient, const CAD::Point3D& velocity) {
    return velocity * -dampingCoefficient; // Negative because damping opposes motion
}

CAD::Point3D projectileTrajectory(const CAD::Point3D& initialPosition, const CAD::Point3D& initialVelocity, 
                                 const CAD::Point3D& acceleration, double time) {
    return initialPosition + initialVelocity * time + acceleration * (0.5 * time * time);
}

double projectileTimeOfFlight(double initialHeight, double initialVelocityY, double gravity) {
    // Solve for time when y = 0: initialHeight + initialVelocityY * t - 0.5 * gravity * t^2 = 0
    // Using quadratic formula: t = (-b ± sqrt(b^2 - 4ac)) / 2a
    // where a = -0.5 * gravity, b = initialVelocityY, c = initialHeight
    
    double a = -0.5 * gravity;
    double b = initialVelocityY;
    double c = initialHeight;
    
    double discriminant = b * b - 4 * a * c;
    
    if (discriminant < 0) {
        // No real solutions (projectile never reaches the ground)
        return 0.0;
    }
    
    double t1 = (-b + std::sqrt(discriminant)) / (2 * a);
    double t2 = (-b - std::sqrt(discriminant)) / (2 * a);
    
    // Return the positive time
    return (t1 > 0) ? t1 : t2;
}

double projectileRange(double initialHeight, const CAD::Point3D& initialVelocity, double gravity) {
    double timeOfFlight = projectileTimeOfFlight(initialHeight, initialVelocity.y, gravity);
    return initialVelocity.x * timeOfFlight;
}

double projectileMaxHeight(double initialHeight, double initialVelocityY, double gravity) {
    // Maximum height is reached when vertical velocity is zero
    // Time to reach maximum height: t = initialVelocityY / gravity
    // Maximum height: initialHeight + initialVelocityY * t - 0.5 * gravity * t^2
    
    double timeToMaxHeight = initialVelocityY / gravity;
    return initialHeight + initialVelocityY * timeToMaxHeight - 0.5 * gravity * timeToMaxHeight * timeToMaxHeight;
}

double kineticEnergy(double mass, const CAD::Point3D& velocity) {
    double velocitySquared = velocity.dot(velocity);
    return 0.5 * mass * velocitySquared;
}

double potentialEnergy(double mass, double height, double gravity) {
    return mass * gravity * height;
}

double elasticPotentialEnergy(double springConstant, double displacement) {
    return 0.5 * springConstant * displacement * displacement;
}

CAD::Point3D momentum(double mass, const CAD::Point3D& velocity) {
    return velocity * mass;
}

CAD::Point3D angularMomentum(const CAD::Point3D& position, const CAD::Point3D& momentum) {
    return position.cross(momentum);
}

double momentOfInertia(double mass, double radius) {
    return mass * radius * radius;
}

double rotationalKineticEnergy(double momentOfInertia, double angularVelocity) {
    return 0.5 * momentOfInertia * angularVelocity * angularVelocity;
}

double centripetalForce(double mass, double velocity, double radius) {
    return mass * velocity * velocity / radius;
}

double centrifugalForce(double mass, double velocity, double radius) {
    return mass * velocity * velocity / radius;
}

CAD::Point3D coriolisForce(double mass, const CAD::Point3D& velocity, const CAD::Point3D& angularVelocity) {
    return velocity.cross(angularVelocity) * (2.0 * mass);
}

double work(const CAD::Point3D& force, const CAD::Point3D& displacement) {
    return force.dot(displacement);
}

double power(const CAD::Point3D& force, const CAD::Point3D& velocity) {
    return force.dot(velocity);
}

CAD::Point3D impulse(const CAD::Point3D& force, double time) {
    return force * time;
}

CAD::Point3D momentumChange(const CAD::Point3D& impulse) {
    return impulse; // By definition, impulse equals change in momentum
}

double coefficientOfRestitution(double initialVelocity1, double initialVelocity2, 
                               double finalVelocity1, double finalVelocity2) {
    return -(finalVelocity2 - finalVelocity1) / (initialVelocity2 - initialVelocity1);
}

std::pair<double, double> collisionVelocities(double mass1, double mass2, 
                                             double initialVelocity1, double initialVelocity2, 
                                             double coefficientOfRestitution) {
    double finalVelocity1 = ((mass1 - coefficientOfRestitution * mass2) * initialVelocity1 + 
                            (1.0 + coefficientOfRestitution) * mass2 * initialVelocity2) / 
                           (mass1 + mass2);
    
    double finalVelocity2 = ((1.0 + coefficientOfRestitution) * mass1 * initialVelocity1 + 
                            (mass2 - coefficientOfRestitution * mass1) * initialVelocity2) / 
                           (mass1 + mass2);
    
    return {finalVelocity1, finalVelocity2};
}

double fluidPressure(double density, double height, double gravity) {
    return density * gravity * height;
}

double buoyantForce(double fluidDensity, double volume, double gravity) {
    return fluidDensity * volume * gravity;
}

double reynoldsNumber(double density, double velocity, double length, double viscosity) {
    return density * velocity * length / viscosity;
}

bool bernoulliEquation(double pressure1, double density, double velocity1, double height1,
                      double pressure2, double velocity2, double height2, double gravity) {
    double leftSide = pressure1 + 0.5 * density * velocity1 * velocity1 + density * gravity * height1;
    double rightSide = pressure2 + 0.5 * density * velocity2 * velocity2 + density * gravity * height2;
    
    // Check if the equation is satisfied within a small tolerance
    return std::abs(leftSide - rightSide) < 1e-10;
}

double flowRate(double area, double velocity) {
    return area * velocity;
}

double terminalVelocity(double mass, double area, double dragCoefficient, 
                       double fluidDensity, double gravity) {
    return std::sqrt(2.0 * mass * gravity / (fluidDensity * area * dragCoefficient));
}

} // namespace Physics
} // namespace Engineering
} // namespace RebelCalc
