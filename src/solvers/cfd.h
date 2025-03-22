#ifndef REBELCALC_SOLVERS_CFD_H
#define REBELCALC_SOLVERS_CFD_H

#include <vector>
#include <string>
#include <map>
#include <functional>
#include <memory>
#include <stdexcept>
#include "../backend/matrix.h"
#include "../engineering/cad.h"

namespace RebelCalc {
namespace Solvers {

/**
 * @brief Enum for fluid models in CFD
 */
enum class FluidModel {
    INCOMPRESSIBLE,         // Incompressible fluid
    COMPRESSIBLE,           // Compressible fluid
    NON_NEWTONIAN,          // Non-Newtonian fluid
    MULTIPHASE              // Multiphase fluid
};

/**
 * @brief Enum for turbulence models in CFD
 */
enum class TurbulenceModel {
    LAMINAR,                // Laminar flow (no turbulence model)
    RANS_K_EPSILON,         // Reynolds-Averaged Navier-Stokes with k-epsilon model
    RANS_K_OMEGA,           // Reynolds-Averaged Navier-Stokes with k-omega model
    RANS_SPALART_ALLMARAS,  // Reynolds-Averaged Navier-Stokes with Spalart-Allmaras model
    LES_SMAGORINSKY,        // Large Eddy Simulation with Smagorinsky model
    LES_DYNAMIC,            // Large Eddy Simulation with dynamic model
    DNS                     // Direct Numerical Simulation
};

/**
 * @brief Enum for boundary condition types in CFD
 */
enum class BoundaryConditionType {
    INLET_VELOCITY,         // Inlet with specified velocity
    INLET_PRESSURE,         // Inlet with specified pressure
    OUTLET_PRESSURE,        // Outlet with specified pressure
    WALL_NO_SLIP,           // Wall with no-slip condition
    WALL_SLIP,              // Wall with slip condition
    SYMMETRY,               // Symmetry boundary
    PERIODIC,               // Periodic boundary
    OPEN                    // Open boundary
};

/**
 * @brief Class for a fluid in CFD
 */
class Fluid {
public:
    /**
     * @brief Constructor for a fluid
     * @param density Density
     * @param viscosity Dynamic viscosity
     * @param model Fluid model
     */
    Fluid(double density, double viscosity, FluidModel model = FluidModel::INCOMPRESSIBLE);

    /**
     * @brief Get the density
     * @return Density
     */
    double getDensity() const;

    /**
     * @brief Set the density
     * @param value Density
     */
    void setDensity(double value);

    /**
     * @brief Get the dynamic viscosity
     * @return Dynamic viscosity
     */
    double getViscosity() const;

    /**
     * @brief Set the dynamic viscosity
     * @param value Dynamic viscosity
     */
    void setViscosity(double value);

    /**
     * @brief Get the kinematic viscosity
     * @return Kinematic viscosity
     */
    double getKinematicViscosity() const;

    /**
     * @brief Get the fluid model
     * @return Fluid model
     */
    FluidModel getModel() const;

    /**
     * @brief Set the fluid model
     * @param model Fluid model
     */
    void setModel(FluidModel model);

    /**
     * @brief Set the specific heat ratio (for compressible fluids)
     * @param value Specific heat ratio
     */
    void setSpecificHeatRatio(double value);

    /**
     * @brief Get the specific heat ratio (for compressible fluids)
     * @return Specific heat ratio
     */
    double getSpecificHeatRatio() const;

    /**
     * @brief Set the speed of sound (for compressible fluids)
     * @param value Speed of sound
     */
    void setSpeedOfSound(double value);

    /**
     * @brief Get the speed of sound (for compressible fluids)
     * @return Speed of sound
     */
    double getSpeedOfSound() const;

    /**
     * @brief Set the power law index (for non-Newtonian fluids)
     * @param value Power law index
     */
    void setPowerLawIndex(double value);

    /**
     * @brief Get the power law index (for non-Newtonian fluids)
     * @return Power law index
     */
    double getPowerLawIndex() const;

    /**
     * @brief Set the consistency index (for non-Newtonian fluids)
     * @param value Consistency index
     */
    void setConsistencyIndex(double value);

    /**
     * @brief Get the consistency index (for non-Newtonian fluids)
     * @return Consistency index
     */
    double getConsistencyIndex() const;

private:
    double m_density;
    double m_viscosity;
    FluidModel m_model;
    double m_specificHeatRatio = 1.4;  // Default for air
    double m_speedOfSound = 343.0;     // Default for air at 20°C
    double m_powerLawIndex = 1.0;      // Default for Newtonian fluid
    double m_consistencyIndex = 0.0;   // Default for Newtonian fluid
};

/**
 * @brief Class for a mesh in CFD
 */
class CFDMesh {
public:
    /**
     * @brief Constructor for a mesh
     */
    CFDMesh();

    /**
     * @brief Add a node to the mesh
     * @param id Node ID
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     */
    void addNode(int id, double x, double y, double z);

    /**
     * @brief Add a cell to the mesh
     * @param id Cell ID
     * @param nodeIds Node IDs
     */
    void addCell(int id, const std::vector<int>& nodeIds);

    /**
     * @brief Add a boundary face to the mesh
     * @param id Face ID
     * @param nodeIds Node IDs
     * @param boundaryType Boundary condition type
     */
    void addBoundaryFace(int id, const std::vector<int>& nodeIds, BoundaryConditionType boundaryType);

    /**
     * @brief Get the number of nodes in the mesh
     * @return Number of nodes
     */
    size_t getNodeCount() const;

    /**
     * @brief Get the number of cells in the mesh
     * @return Number of cells
     */
    size_t getCellCount() const;

    /**
     * @brief Get the number of boundary faces in the mesh
     * @return Number of boundary faces
     */
    size_t getBoundaryFaceCount() const;

    /**
     * @brief Get a node by ID
     * @param id Node ID
     * @return Node coordinates
     * @throws std::out_of_range if the node doesn't exist
     */
    Engineering::CAD::Point3D getNode(int id) const;

    /**
     * @brief Get a cell by ID
     * @param id Cell ID
     * @return Node IDs of the cell
     * @throws std::out_of_range if the cell doesn't exist
     */
    const std::vector<int>& getCell(int id) const;

    /**
     * @brief Get a boundary face by ID
     * @param id Face ID
     * @return Node IDs of the face
     * @throws std::out_of_range if the face doesn't exist
     */
    const std::vector<int>& getBoundaryFace(int id) const;

    /**
     * @brief Get the boundary condition type of a face
     * @param id Face ID
     * @return Boundary condition type
     * @throws std::out_of_range if the face doesn't exist
     */
    BoundaryConditionType getBoundaryConditionType(int id) const;

    /**
     * @brief Set the boundary condition type of a face
     * @param id Face ID
     * @param type Boundary condition type
     * @throws std::out_of_range if the face doesn't exist
     */
    void setBoundaryConditionType(int id, BoundaryConditionType type);

    /**
     * @brief Generate a structured mesh for a rectangular domain
     * @param xMin Minimum x coordinate
     * @param xMax Maximum x coordinate
     * @param yMin Minimum y coordinate
     * @param yMax Maximum y coordinate
     * @param zMin Minimum z coordinate
     * @param zMax Maximum z coordinate
     * @param nx Number of cells in x direction
     * @param ny Number of cells in y direction
     * @param nz Number of cells in z direction
     */
    void generateStructuredMesh(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax,
                               int nx, int ny, int nz);

    /**
     * @brief Import a mesh from a file
     * @param filename Filename
     * @return true if the mesh was imported successfully, false otherwise
     */
    bool importFromFile(const std::string& filename);

    /**
     * @brief Export the mesh to a file
     * @param filename Filename
     * @return true if the mesh was exported successfully, false otherwise
     */
    bool exportToFile(const std::string& filename) const;

private:
    std::map<int, Engineering::CAD::Point3D> m_nodes;
    std::map<int, std::vector<int>> m_cells;
    std::map<int, std::vector<int>> m_boundaryFaces;
    std::map<int, BoundaryConditionType> m_boundaryConditionTypes;
};

/**
 * @brief Class for a boundary condition in CFD
 */
class BoundaryCondition {
public:
    /**
     * @brief Constructor for a boundary condition
     * @param type Boundary condition type
     * @param faceIds Face IDs
     */
    BoundaryCondition(BoundaryConditionType type, const std::vector<int>& faceIds);

    /**
     * @brief Get the boundary condition type
     * @return Boundary condition type
     */
    BoundaryConditionType getType() const;

    /**
     * @brief Get the face IDs
     * @return Face IDs
     */
    const std::vector<int>& getFaceIds() const;

    /**
     * @brief Set the velocity (for INLET_VELOCITY)
     * @param velocity Velocity vector
     */
    void setVelocity(const Engineering::CAD::Point3D& velocity);

    /**
     * @brief Get the velocity (for INLET_VELOCITY)
     * @return Velocity vector
     * @throws std::runtime_error if the boundary condition type is not INLET_VELOCITY
     */
    Engineering::CAD::Point3D getVelocity() const;

    /**
     * @brief Set the pressure (for INLET_PRESSURE or OUTLET_PRESSURE)
     * @param pressure Pressure
     */
    void setPressure(double pressure);

    /**
     * @brief Get the pressure (for INLET_PRESSURE or OUTLET_PRESSURE)
     * @return Pressure
     * @throws std::runtime_error if the boundary condition type is not INLET_PRESSURE or OUTLET_PRESSURE
     */
    double getPressure() const;

    /**
     * @brief Set the wall temperature (for WALL_NO_SLIP or WALL_SLIP)
     * @param temperature Wall temperature
     */
    void setWallTemperature(double temperature);

    /**
     * @brief Get the wall temperature (for WALL_NO_SLIP or WALL_SLIP)
     * @return Wall temperature
     * @throws std::runtime_error if the boundary condition type is not WALL_NO_SLIP or WALL_SLIP
     */
    double getWallTemperature() const;

    /**
     * @brief Set the wall heat flux (for WALL_NO_SLIP or WALL_SLIP)
     * @param heatFlux Wall heat flux
     */
    void setWallHeatFlux(double heatFlux);

    /**
     * @brief Get the wall heat flux (for WALL_NO_SLIP or WALL_SLIP)
     * @return Wall heat flux
     * @throws std::runtime_error if the boundary condition type is not WALL_NO_SLIP or WALL_SLIP
     */
    double getWallHeatFlux() const;

private:
    BoundaryConditionType m_type;
    std::vector<int> m_faceIds;
    Engineering::CAD::Point3D m_velocity;
    double m_pressure = 0.0;
    double m_wallTemperature = 0.0;
    double m_wallHeatFlux = 0.0;
};

/**
 * @brief Class for a CFD model
 */
class CFDModel {
public:
    /**
     * @brief Constructor for a CFD model
     * @param name Model name
     */
    CFDModel(const std::string& name = "");

    /**
     * @brief Set the fluid
     * @param fluid Fluid
     */
    void setFluid(const Fluid& fluid);

    /**
     * @brief Get the fluid
     * @return Fluid
     */
    const Fluid& getFluid() const;

    /**
     * @brief Set the mesh
     * @param mesh Mesh
     */
    void setMesh(const CFDMesh& mesh);

    /**
     * @brief Get the mesh
     * @return Mesh
     */
    const CFDMesh& getMesh() const;

    /**
     * @brief Add a boundary condition
     * @param condition Boundary condition
     */
    void addBoundaryCondition(const BoundaryCondition& condition);

    /**
     * @brief Get the boundary conditions
     * @return Boundary conditions
     */
    const std::vector<BoundaryCondition>& getBoundaryConditions() const;

    /**
     * @brief Set the turbulence model
     * @param model Turbulence model
     */
    void setTurbulenceModel(TurbulenceModel model);

    /**
     * @brief Get the turbulence model
     * @return Turbulence model
     */
    TurbulenceModel getTurbulenceModel() const;

    /**
     * @brief Set the reference pressure
     * @param pressure Reference pressure
     */
    void setReferencePressure(double pressure);

    /**
     * @brief Get the reference pressure
     * @return Reference pressure
     */
    double getReferencePressure() const;

    /**
     * @brief Set the reference temperature
     * @param temperature Reference temperature
     */
    void setReferenceTemperature(double temperature);

    /**
     * @brief Get the reference temperature
     * @return Reference temperature
     */
    double getReferenceTemperature() const;

    /**
     * @brief Set the reference length
     * @param length Reference length
     */
    void setReferenceLength(double length);

    /**
     * @brief Get the reference length
     * @return Reference length
     */
    double getReferenceLength() const;

    /**
     * @brief Set the reference velocity
     * @param velocity Reference velocity
     */
    void setReferenceVelocity(double velocity);

    /**
     * @brief Get the reference velocity
     * @return Reference velocity
     */
    double getReferenceVelocity() const;

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
    Fluid m_fluid;
    CFDMesh m_mesh;
    std::vector<BoundaryCondition> m_boundaryConditions;
    TurbulenceModel m_turbulenceModel = TurbulenceModel::LAMINAR;
    double m_referencePressure = 101325.0;  // Default: 1 atm
    double m_referenceTemperature = 293.15; // Default: 20°C
    double m_referenceLength = 1.0;
    double m_referenceVelocity = 1.0;
};

/**
 * @brief Class for a CFD solver
 */
class CFDSolver {
public:
    /**
     * @brief Constructor for a CFD solver
     * @param model CFD model
     */
    CFDSolver(const CFDModel& model);

    /**
     * @brief Set the maximum number of iterations
     * @param iterations Maximum number of iterations
     */
    void setMaxIterations(int iterations);

    /**
     * @brief Get the maximum number of iterations
     * @return Maximum number of iterations
     */
    int getMaxIterations() const;

    /**
     * @brief Set the convergence tolerance
     * @param tolerance Convergence tolerance
     */
    void setConvergenceTolerance(double tolerance);

    /**
     * @brief Get the convergence tolerance
     * @return Convergence tolerance
     */
    double getConvergenceTolerance() const;

    /**
     * @brief Set the time step
     * @param timeStep Time step
     */
    void setTimeStep(double timeStep);

    /**
     * @brief Get the time step
     * @return Time step
     */
    double getTimeStep() const;

    /**
     * @brief Set the maximum simulation time
     * @param time Maximum simulation time
     */
    void setMaxTime(double time);

    /**
     * @brief Get the maximum simulation time
     * @return Maximum simulation time
     */
    double getMaxTime() const;

    /**
     * @brief Set the CFL number
     * @param cfl CFL number
     */
    void setCFL(double cfl);

    /**
     * @brief Get the CFL number
     * @return CFL number
     */
    double getCFL() const;

    /**
     * @brief Set whether to use adaptive time stepping
     * @param adaptive Whether to use adaptive time stepping
     */
    void setAdaptiveTimeStep(bool adaptive);

    /**
     * @brief Get whether to use adaptive time stepping
     * @return Whether to use adaptive time stepping
     */
    bool getAdaptiveTimeStep() const;

    /**
     * @brief Set the under-relaxation factor for pressure
     * @param factor Under-relaxation factor
     */
    void setPressureUnderRelaxation(double factor);

    /**
     * @brief Get the under-relaxation factor for pressure
     * @return Under-relaxation factor
     */
    double getPressureUnderRelaxation() const;

    /**
     * @brief Set the under-relaxation factor for velocity
     * @param factor Under-relaxation factor
     */
    void setVelocityUnderRelaxation(double factor);

    /**
     * @brief Get the under-relaxation factor for velocity
     * @return Under-relaxation factor
     */
    double getVelocityUnderRelaxation() const;

    /**
     * @brief Set the discretization scheme for convection terms
     * @param scheme Discretization scheme (1=First-order upwind, 2=Second-order upwind, 3=QUICK, 4=Central)
     */
    void setConvectionScheme(int scheme);

    /**
     * @brief Get the discretization scheme for convection terms
     * @return Discretization scheme
     */
    int getConvectionScheme() const;

    /**
     * @brief Solve the CFD model
     * @return true if the model was solved successfully, false otherwise
     */
    bool solve();

    /**
     * @brief Get the velocity at a cell
     * @param cellId Cell ID
     * @return Velocity vector
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the cell doesn't exist
     */
    Engineering::CAD::Point3D getVelocity(int cellId) const;

    /**
     * @brief Get the pressure at a cell
     * @param cellId Cell ID
     * @return Pressure
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the cell doesn't exist
     */
    double getPressure(int cellId) const;

    /**
     * @brief Get the temperature at a cell
     * @param cellId Cell ID
     * @return Temperature
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the cell doesn't exist
     */
    double getTemperature(int cellId) const;

    /**
     * @brief Get the density at a cell
     * @param cellId Cell ID
     * @return Density
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the cell doesn't exist
     */
    double getDensity(int cellId) const;

    /**
     * @brief Get the turbulent kinetic energy at a cell
     * @param cellId Cell ID
     * @return Turbulent kinetic energy
     * @throws std::runtime_error if the model hasn't been solved or if the turbulence model is not RANS_K_EPSILON or RANS_K_OMEGA
     * @throws std::out_of_range if the cell doesn't exist
     */
    double getTurbulentKineticEnergy(int cellId) const;

    /**
     * @brief Get the turbulent dissipation rate at a cell
     * @param cellId Cell ID
     * @return Turbulent dissipation rate
     * @throws std::runtime_error if the model hasn't been solved or if the turbulence model is not RANS_K_EPSILON
     * @throws std::out_of_range if the cell doesn't exist
     */
    double getTurbulentDissipationRate(int cellId) const;

    /**
     * @brief Get the specific dissipation rate at a cell
     * @param cellId Cell ID
     * @return Specific dissipation rate
     * @throws std::runtime_error if the model hasn't been solved or if the turbulence model is not RANS_K_OMEGA
     * @throws std::out_of_range if the cell doesn't exist
     */
    double getSpecificDissipationRate(int cellId) const;

    /**
     * @brief Get the wall shear stress at a boundary face
     * @param faceId Face ID
     * @return Wall shear stress
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the face doesn't exist
     */
    double getWallShearStress(int faceId) const;

    /**
     * @brief Get the heat flux at a boundary face
     * @param faceId Face ID
     * @return Heat flux
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the face doesn't exist
     */
    double getHeatFlux(int faceId) const;

    /**
     * @brief Get the drag coefficient
     * @return Drag coefficient
     * @throws std::runtime_error if the model hasn't been solved
     */
    double getDragCoefficient() const;

    /**
     * @brief Get the lift coefficient
     * @return Lift coefficient
     * @throws std::runtime_error if the model hasn't been solved
     */
    double getLiftCoefficient() const;

    /**
     * @brief Get the moment coefficient
     * @return Moment coefficient
     * @throws std::runtime_error if the model hasn't been solved
     */
    double getMomentCoefficient() const;

    /**
     * @brief Get the Nusselt number
     * @return Nusselt number
     * @throws std::runtime_error if the model hasn't been solved
     */
    double getNusseltNumber() const;

    /**
     * @brief Get the Reynolds number
     * @return Reynolds number
     */
    double getReynoldsNumber() const;

    /**
     * @brief Get the Mach number
     * @return Mach number
     */
    double getMachNumber() const;

    /**
     * @brief Get the Prandtl number
     * @return Prandtl number
     */
    double getPrandtlNumber() const;

    /**
     * @brief Export the results to a file
     * @param filename Filename
     * @return true if the results were exported successfully, false otherwise
     */
    bool exportResults(const std::string& filename) const;

private:
    const CFDModel& m_model;
    int m_maxIterations = 1000;
    double m_convergenceTolerance = 1e-6;
    double m_timeStep = 0.001;
    double m_maxTime = 1.0;
    double m_cfl = 0.9;
    bool m_adaptiveTimeStep = true;
    double m_pressureUnderRelaxation = 0.3;
    double m_velocityUnderRelaxation = 0.7;
    int m_convectionScheme = 2;  // Default: Second-order upwind
    bool m_solved = false;
    std::map<int, Engineering::CAD::Point3D> m_velocities;
    std::map<int, double> m_pressures;
    std::map<int, double> m_temperatures;
    std::map<int, double> m_densities;
    std::map<int, double> m_turbulentKineticEnergies;
    std::map<int, double> m_turbulentDissipationRates;
    std::map<int, double> m_specificDissipationRates;
    std::map<int, double> m_wallShearStresses;
    std::map<int, double> m_heatFluxes;
    double m_dragCoefficient = 0.0;
    double m_liftCoefficient = 0.0;
    double m_momentCoefficient = 0.0;
    double m_nusseltNumber = 0.0;

    // Helper methods
    void initializeFields();
    void solveNavierStokes();
    void solveEnergy();
    void solveTurbulence();
    void calculateForces();
    void calculateHeatTransfer();
};

} // namespace Solvers
} // namespace RebelCalc

#endif // REBELCALC_SOLVERS_CFD_H
