#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <memory>

#include "../src/integration/api.h"
#include "../src/engineering/cad.h"

using namespace RebelCalc::Integration;
using namespace RebelCalc::Engineering;

// Helper function to print a Point3D
void printPoint(const CAD::Point3D& point, const std::string& name) {
    std::cout << name << ": (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
}

// Helper function to print a section header
void printHeader(const std::string& header) {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "  " << header << std::endl;
    std::cout << std::string(80, '=') << std::endl;
}

// Demo for integration with RebelCAD
void demoCADIntegration() {
    printHeader("RebelCAD Integration Demonstration");
    
    // Get the integration factory
    IntegrationFactory& factory = IntegrationFactory::getInstance();
    
    // Get the available components
    std::vector<Component> components = factory.getAvailableComponents();
    
    std::cout << "Available components:" << std::endl;
    for (const auto& component : components) {
        switch (component) {
            case Component::REBELCAD:
                std::cout << "- RebelCAD" << std::endl;
                break;
            case Component::REBELENGINE:
                std::cout << "- RebelENGINE" << std::endl;
                break;
            case Component::REBELFLOW:
                std::cout << "- RebelFLOW" << std::endl;
                break;
            case Component::REBELCODE:
                std::cout << "- RebelCODE" << std::endl;
                break;
            case Component::REBELDESK:
                std::cout << "- RebelDESK" << std::endl;
                break;
            case Component::REBELSCRIBE:
                std::cout << "- RebelSCRIBE" << std::endl;
                break;
        }
    }
    
    // Create a CAD integration interface
    std::shared_ptr<IntegrationInterface> cadInterface = factory.createInterface(Component::REBELCAD);
    
    if (!cadInterface) {
        std::cout << "RebelCAD integration is not available." << std::endl;
        return;
    }
    
    std::cout << "\nCreated integration interface for " << cadInterface->getComponentName() << std::endl;
    
    // Check if RebelCAD is available
    if (!cadInterface->isAvailable()) {
        std::cout << "RebelCAD is not available on this system." << std::endl;
        return;
    }
    
    std::cout << "RebelCAD is available on this system." << std::endl;
    
    // Initialize the integration
    if (!cadInterface->initialize()) {
        std::cout << "Failed to initialize RebelCAD integration." << std::endl;
        return;
    }
    
    std::cout << "RebelCAD integration initialized successfully." << std::endl;
    
    // Get the available commands
    std::vector<std::string> commands = cadInterface->getAvailableCommands();
    
    std::cout << "\nAvailable commands:" << std::endl;
    for (const auto& command : commands) {
        std::cout << "- " << command << std::endl;
    }
    
    // List available models
    DataExchange data;
    if (cadInterface->executeCommand("list_models", data)) {
        std::vector<std::string> models = std::get<std::vector<std::string>>(data.getData("models"));
        
        std::cout << "\nAvailable models in RebelCAD:" << std::endl;
        for (const auto& model : models) {
            std::cout << "- " << model << std::endl;
        }
    } else {
        std::cout << "Failed to list models." << std::endl;
    }
    
    // Import a model
    data.clear();
    data.setData("model_name", std::string("cube"));
    
    if (cadInterface->executeCommand("import_model", data)) {
        std::vector<CAD::Triangle> triangles = std::get<std::vector<CAD::Triangle>>(data.getData("triangles"));
        
        std::cout << "\nImported model 'cube' with " << triangles.size() << " triangles." << std::endl;
        
        // Print the first triangle
        if (!triangles.empty()) {
            std::cout << "First triangle:" << std::endl;
            printPoint(triangles[0].p1, "Point 1");
            printPoint(triangles[0].p2, "Point 2");
            printPoint(triangles[0].p3, "Point 3");
        }
    } else {
        std::cout << "Failed to import model 'cube'." << std::endl;
    }
    
    // Perform a measurement
    data.clear();
    data.setData("model_name", std::string("cube"));
    data.setData("measurement_type", std::string("volume"));
    
    if (cadInterface->executeCommand("measure", data)) {
        double volume = std::get<double>(data.getData("result"));
        
        std::cout << "\nVolume of 'cube': " << volume << " cubic units." << std::endl;
    } else {
        std::cout << "Failed to measure volume of 'cube'." << std::endl;
    }
    
    // Perform a boolean operation
    data.clear();
    data.setData("operation", std::string("union"));
    data.setData("model1_name", std::string("cube"));
    data.setData("model2_name", std::string("sphere"));
    data.setData("result_model_name", std::string("combined"));
    
    if (cadInterface->executeCommand("boolean_operation", data)) {
        bool success = std::get<bool>(data.getData("success"));
        
        if (success) {
            std::cout << "\nPerformed boolean union of 'cube' and 'sphere' to create 'combined'." << std::endl;
        } else {
            std::cout << "\nFailed to perform boolean union." << std::endl;
        }
    } else {
        std::cout << "Failed to execute boolean operation command." << std::endl;
    }
    
    // List models again to see the new model
    data.clear();
    if (cadInterface->executeCommand("list_models", data)) {
        std::vector<std::string> models = std::get<std::vector<std::string>>(data.getData("models"));
        
        std::cout << "\nUpdated list of models in RebelCAD:" << std::endl;
        for (const auto& model : models) {
            std::cout << "- " << model << std::endl;
        }
    } else {
        std::cout << "Failed to list models." << std::endl;
    }
    
    // Shutdown the integration
    cadInterface->shutdown();
    
    std::cout << "\nRebelCAD integration shut down successfully." << std::endl;
}

int main() {
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "                      RebelCALC Integration Demo" << std::endl;
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "\nThis program demonstrates the integration capabilities of RebelCALC with other RebelSUITE components." << std::endl;
    
    demoCADIntegration();
    
    std::cout << "\n" << std::string(80, '*') << std::endl;
    std::cout << "                      End of Integration Demo" << std::endl;
    std::cout << std::string(80, '*') << std::endl;
    
    return 0;
}
