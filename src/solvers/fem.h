#ifndef REBELCALC_SOLVERS_FEM_H
#define REBELCALC_SOLVERS_FEM_H

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
 * @brief Enum for element types in finite element analysis
 */
enum class ElementType {
    LINEAR_1D,      // 1D linear element (2 nodes)
    QUADRATIC_1D,   // 1D quadratic element (3 nodes)
    LINEAR_2D_TRI,  // 2D linear triangular element (3 nodes)
    LINEAR_2D_QUAD, // 2D linear quadrilateral element (4 nodes)
    QUADRATIC_2D_TRI, // 2D quadratic triangular element (6 nodes)
    QUADRATIC_2D_QUAD, // 2D quadratic quadrilateral element (8 nodes)
    LINEAR_3D_TET,  // 3D linear tetrahedral element (4 nodes)
    LINEAR_3D_HEX,  // 3D linear hexahedral element (8 nodes)
    QUADRATIC_3D_TET, // 3D quadratic tetrahedral element (10 nodes)
    QUADRATIC_3D_HEX  // 3D quadratic hexahedral element (20 nodes)
};

/**
 * @brief Enum for analysis types in finite element analysis
 */
enum class AnalysisType {
    STATIC,         // Static analysis
    MODAL,          // Modal analysis (natural frequencies and mode shapes)
    HARMONIC,       // Harmonic analysis (frequency response)
    TRANSIENT,      // Transient analysis (time response)
    BUCKLING,       // Buckling analysis
    THERMAL,        // Thermal analysis
    COUPLED_THERMAL_STRUCTURAL // Coupled thermal-structural analysis
};

/**
 * @brief Enum for material models in finite element analysis
 */
enum class MaterialModel {
    LINEAR_ELASTIC,         // Linear elastic material
    NONLINEAR_ELASTIC,      // Nonlinear elastic material
    PLASTIC,                // Plastic material
    VISCOELASTIC,           // Viscoelastic material
    HYPERELASTIC,           // Hyperelastic material
    ORTHOTROPIC,            // Orthotropic material
    ANISOTROPIC             // Anisotropic material
};

/**
 * @brief Class for a node in finite element analysis
 */
class Node {
public:
    /**
     * @brief Constructor for a node
     * @param id Node ID
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     */
    Node(int id, double x, double y, double z);

    /**
     * @brief Get the node ID
     * @return Node ID
     */
    int getId() const;

    /**
     * @brief Get the node coordinates
     * @return Node coordinates as a 3D point
     */
    Engineering::CAD::Point3D getCoordinates() const;

    /**
     * @brief Set a degree of freedom constraint
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @param value Constraint value
     */
    void setConstraint(int dof, double value);

    /**
     * @brief Check if a degree of freedom is constrained
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @return true if the degree of freedom is constrained, false otherwise
     */
    bool isConstrained(int dof) const;

    /**
     * @brief Get the constraint value for a degree of freedom
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @return Constraint value
     * @throws std::out_of_range if the degree of freedom is not constrained
     */
    double getConstraintValue(int dof) const;

    /**
     * @brief Apply a force to the node
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @param value Force value
     */
    void applyForce(int dof, double value);

    /**
     * @brief Get the force applied to a degree of freedom
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @return Force value
     */
    double getForce(int dof) const;

    /**
     * @brief Get the displacement for a degree of freedom
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @return Displacement value
     */
    double getDisplacement(int dof) const;

    /**
     * @brief Set the displacement for a degree of freedom
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @param value Displacement value
     */
    void setDisplacement(int dof, double value);

private:
    int m_id;
    Engineering::CAD::Point3D m_coordinates;
    std::map<int, double> m_constraints;
    std::map<int, double> m_forces;
    std::map<int, double> m_displacements;
};

/**
 * @brief Class for a material in finite element analysis
 */
class Material {
public:
    /**
     * @brief Constructor for a linear elastic material
     * @param id Material ID
     * @param name Material name
     * @param youngModulus Young's modulus
     * @param poissonRatio Poisson's ratio
     * @param density Density
     */
    Material(int id, const std::string& name, double youngModulus, double poissonRatio, double density);

    /**
     * @brief Get the material ID
     * @return Material ID
     */
    int getId() const;

    /**
     * @brief Get the material name
     * @return Material name
     */
    const std::string& getName() const;

    /**
     * @brief Get the material model
     * @return Material model
     */
    MaterialModel getModel() const;

    /**
     * @brief Set the material model
     * @param model Material model
     */
    void setModel(MaterialModel model);

    /**
     * @brief Get Young's modulus
     * @return Young's modulus
     */
    double getYoungModulus() const;

    /**
     * @brief Set Young's modulus
     * @param value Young's modulus
     */
    void setYoungModulus(double value);

    /**
     * @brief Get Poisson's ratio
     * @return Poisson's ratio
     */
    double getPoissonRatio() const;

    /**
     * @brief Set Poisson's ratio
     * @param value Poisson's ratio
     */
    void setPoissonRatio(double value);

    /**
     * @brief Get density
     * @return Density
     */
    double getDensity() const;

    /**
     * @brief Set density
     * @param value Density
     */
    void setDensity(double value);

    /**
     * @brief Get shear modulus
     * @return Shear modulus
     */
    double getShearModulus() const;

    /**
     * @brief Get bulk modulus
     * @return Bulk modulus
     */
    double getBulkModulus() const;

    /**
     * @brief Set thermal expansion coefficient
     * @param value Thermal expansion coefficient
     */
    void setThermalExpansionCoefficient(double value);

    /**
     * @brief Get thermal expansion coefficient
     * @return Thermal expansion coefficient
     */
    double getThermalExpansionCoefficient() const;

    /**
     * @brief Set thermal conductivity
     * @param value Thermal conductivity
     */
    void setThermalConductivity(double value);

    /**
     * @brief Get thermal conductivity
     * @return Thermal conductivity
     */
    double getThermalConductivity() const;

    /**
     * @brief Set specific heat
     * @param value Specific heat
     */
    void setSpecificHeat(double value);

    /**
     * @brief Get specific heat
     * @return Specific heat
     */
    double getSpecificHeat() const;

    /**
     * @brief Set yield strength
     * @param value Yield strength
     */
    void setYieldStrength(double value);

    /**
     * @brief Get yield strength
     * @return Yield strength
     */
    double getYieldStrength() const;

    /**
     * @brief Set ultimate strength
     * @param value Ultimate strength
     */
    void setUltimateStrength(double value);

    /**
     * @brief Get ultimate strength
     * @return Ultimate strength
     */
    double getUltimateStrength() const;

private:
    int m_id;
    std::string m_name;
    MaterialModel m_model = MaterialModel::LINEAR_ELASTIC;
    double m_youngModulus;
    double m_poissonRatio;
    double m_density;
    double m_thermalExpansionCoefficient = 0.0;
    double m_thermalConductivity = 0.0;
    double m_specificHeat = 0.0;
    double m_yieldStrength = 0.0;
    double m_ultimateStrength = 0.0;
};

/**
 * @brief Base class for an element in finite element analysis
 */
class Element {
public:
    /**
     * @brief Constructor for an element
     * @param id Element ID
     * @param type Element type
     * @param materialId Material ID
     * @param nodeIds Node IDs
     */
    Element(int id, ElementType type, int materialId, const std::vector<int>& nodeIds);

    /**
     * @brief Virtual destructor
     */
    virtual ~Element() = default;

    /**
     * @brief Get the element ID
     * @return Element ID
     */
    int getId() const;

    /**
     * @brief Get the element type
     * @return Element type
     */
    ElementType getType() const;

    /**
     * @brief Get the material ID
     * @return Material ID
     */
    int getMaterialId() const;

    /**
     * @brief Get the node IDs
     * @return Node IDs
     */
    const std::vector<int>& getNodeIds() const;

    /**
     * @brief Calculate the element stiffness matrix
     * @param nodes Map of node IDs to nodes
     * @param material Material
     * @return Element stiffness matrix
     */
    virtual rebelcalc::Matrix calculateStiffnessMatrix(
        const std::map<int, Node>& nodes, const Material& material) const = 0;

    /**
     * @brief Calculate the element mass matrix
     * @param nodes Map of node IDs to nodes
     * @param material Material
     * @return Element mass matrix
     */
    virtual rebelcalc::Matrix calculateMassMatrix(
        const std::map<int, Node>& nodes, const Material& material) const = 0;

    /**
     * @brief Calculate the element load vector
     * @param nodes Map of node IDs to nodes
     * @param material Material
     * @return Element load vector
     */
    virtual std::vector<double> calculateLoadVector(
        const std::map<int, Node>& nodes, const Material& material) const = 0;

    /**
     * @brief Calculate the element stress
     * @param nodes Map of node IDs to nodes
     * @param material Material
     * @param displacements Element displacements
     * @return Element stress
     */
    virtual std::vector<double> calculateStress(
        const std::map<int, Node>& nodes, const Material& material,
        const std::vector<double>& displacements) const = 0;

    /**
     * @brief Calculate the element strain
     * @param nodes Map of node IDs to nodes
     * @param material Material
     * @param displacements Element displacements
     * @return Element strain
     */
    virtual std::vector<double> calculateStrain(
        const std::map<int, Node>& nodes, const Material& material,
        const std::vector<double>& displacements) const = 0;

protected:
    int m_id;
    ElementType m_type;
    int m_materialId;
    std::vector<int> m_nodeIds;
};

/**
 * @brief Class for a finite element model
 */
class FEModel {
public:
    /**
     * @brief Constructor for a finite element model
     * @param name Model name
     */
    FEModel(const std::string& name = "");

    /**
     * @brief Add a node to the model
     * @param id Node ID
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     * @return Reference to the added node
     */
    Node& addNode(int id, double x, double y, double z);

    /**
     * @brief Add a material to the model
     * @param id Material ID
     * @param name Material name
     * @param youngModulus Young's modulus
     * @param poissonRatio Poisson's ratio
     * @param density Density
     * @return Reference to the added material
     */
    Material& addMaterial(int id, const std::string& name, double youngModulus, double poissonRatio, double density);

    /**
     * @brief Add an element to the model
     * @param id Element ID
     * @param type Element type
     * @param materialId Material ID
     * @param nodeIds Node IDs
     * @return Reference to the added element
     */
    Element& addElement(int id, ElementType type, int materialId, const std::vector<int>& nodeIds);

    /**
     * @brief Get a node by ID
     * @param id Node ID
     * @return Reference to the node
     * @throws std::out_of_range if the node doesn't exist
     */
    Node& getNode(int id);

    /**
     * @brief Get a node by ID
     * @param id Node ID
     * @return Const reference to the node
     * @throws std::out_of_range if the node doesn't exist
     */
    const Node& getNode(int id) const;

    /**
     * @brief Get a material by ID
     * @param id Material ID
     * @return Reference to the material
     * @throws std::out_of_range if the material doesn't exist
     */
    Material& getMaterial(int id);

    /**
     * @brief Get a material by ID
     * @param id Material ID
     * @return Const reference to the material
     * @throws std::out_of_range if the material doesn't exist
     */
    const Material& getMaterial(int id) const;

    /**
     * @brief Get an element by ID
     * @param id Element ID
     * @return Reference to the element
     * @throws std::out_of_range if the element doesn't exist
     */
    Element& getElement(int id);

    /**
     * @brief Get an element by ID
     * @param id Element ID
     * @return Const reference to the element
     * @throws std::out_of_range if the element doesn't exist
     */
    const Element& getElement(int id) const;

    /**
     * @brief Get all nodes
     * @return Map of node IDs to nodes
     */
    const std::map<int, Node>& getNodes() const;

    /**
     * @brief Get all materials
     * @return Map of material IDs to materials
     */
    const std::map<int, Material>& getMaterials() const;

    /**
     * @brief Get all elements
     * @return Map of element IDs to elements
     */
    const std::map<int, std::unique_ptr<Element>>& getElements() const;

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
     * @brief Apply a constraint to a node
     * @param nodeId Node ID
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @param value Constraint value
     */
    void applyConstraint(int nodeId, int dof, double value);

    /**
     * @brief Apply a force to a node
     * @param nodeId Node ID
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @param value Force value
     */
    void applyForce(int nodeId, int dof, double value);

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
    std::map<int, Node> m_nodes;
    std::map<int, Material> m_materials;
    std::map<int, std::unique_ptr<Element>> m_elements;
};

/**
 * @brief Class for a finite element solver
 */
class FESolver {
public:
    /**
     * @brief Constructor for a finite element solver
     * @param model Finite element model
     */
    FESolver(const FEModel& model);

    /**
     * @brief Set the analysis type
     * @param type Analysis type
     */
    void setAnalysisType(AnalysisType type);

    /**
     * @brief Get the analysis type
     * @return Analysis type
     */
    AnalysisType getAnalysisType() const;

    /**
     * @brief Solve the finite element model
     * @return true if the model was solved successfully, false otherwise
     */
    bool solve();

    /**
     * @brief Get the displacement at a node
     * @param nodeId Node ID
     * @param dof Degree of freedom (0=x, 1=y, 2=z, 3=rx, 4=ry, 5=rz)
     * @return Displacement value
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the node doesn't exist
     */
    double getDisplacement(int nodeId, int dof) const;

    /**
     * @brief Get the stress at an element
     * @param elementId Element ID
     * @return Element stress
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the element doesn't exist
     */
    std::vector<double> getStress(int elementId) const;

    /**
     * @brief Get the strain at an element
     * @param elementId Element ID
     * @return Element strain
     * @throws std::runtime_error if the model hasn't been solved
     * @throws std::out_of_range if the element doesn't exist
     */
    std::vector<double> getStrain(int elementId) const;

    /**
     * @brief Get the natural frequencies
     * @return Natural frequencies
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not MODAL
     */
    std::vector<double> getNaturalFrequencies() const;

    /**
     * @brief Get the mode shapes
     * @return Mode shapes
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not MODAL
     */
    std::vector<std::vector<double>> getModeShapes() const;

    /**
     * @brief Get the frequency response
     * @return Frequency response
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not HARMONIC
     */
    std::vector<std::vector<double>> getFrequencyResponse() const;

    /**
     * @brief Get the time response
     * @return Time response
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not TRANSIENT
     */
    std::vector<std::vector<double>> getTimeResponse() const;

    /**
     * @brief Get the buckling factors
     * @return Buckling factors
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not BUCKLING
     */
    std::vector<double> getBucklingFactors() const;

    /**
     * @brief Get the buckling modes
     * @return Buckling modes
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not BUCKLING
     */
    std::vector<std::vector<double>> getBucklingModes() const;

    /**
     * @brief Get the temperature at a node
     * @param nodeId Node ID
     * @return Temperature value
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not THERMAL or COUPLED_THERMAL_STRUCTURAL
     * @throws std::out_of_range if the node doesn't exist
     */
    double getTemperature(int nodeId) const;

    /**
     * @brief Get the heat flux at an element
     * @param elementId Element ID
     * @return Element heat flux
     * @throws std::runtime_error if the model hasn't been solved or if the analysis type is not THERMAL or COUPLED_THERMAL_STRUCTURAL
     * @throws std::out_of_range if the element doesn't exist
     */
    std::vector<double> getHeatFlux(int elementId) const;

private:
    const FEModel& m_model;
    AnalysisType m_analysisType = AnalysisType::STATIC;
    bool m_solved = false;
    std::map<int, std::map<int, double>> m_displacements;
    std::map<int, std::vector<double>> m_stresses;
    std::map<int, std::vector<double>> m_strains;
    std::vector<double> m_naturalFrequencies;
    std::vector<std::vector<double>> m_modeShapes;
    std::vector<std::vector<double>> m_frequencyResponse;
    std::vector<std::vector<double>> m_timeResponse;
    std::vector<double> m_bucklingFactors;
    std::vector<std::vector<double>> m_bucklingModes;
    std::map<int, double> m_temperatures;
    std::map<int, std::vector<double>> m_heatFluxes;

    // Helper methods
    void assembleGlobalMatrices(rebelcalc::Matrix& stiffnessMatrix, rebelcalc::Matrix& massMatrix, std::vector<double>& loadVector);
    void applyBoundaryConditions(rebelcalc::Matrix& stiffnessMatrix, rebelcalc::Matrix& massMatrix, std::vector<double>& loadVector);
    void solveStatic();
    void solveModal();
    void solveHarmonic();
    void solveTransient();
    void solveBuckling();
    void solveThermal();
    void solveCoupledThermalStructural();
    void calculateStressAndStrain();
};

} // namespace Solvers
} // namespace RebelCalc

#endif // REBELCALC_SOLVERS_FEM_H
