#include "api.h"
#include "cad.h"
#include <stdexcept>
#include <algorithm>

namespace RebelCalc {
namespace Integration {

// DataExchange implementation
void DataExchange::setData(const std::string& key, const DataVariant& data) {
    m_data[key] = data;
}

DataVariant DataExchange::getData(const std::string& key) const {
    auto it = m_data.find(key);
    if (it == m_data.end()) {
        throw std::out_of_range("Key not found in DataExchange: " + key);
    }
    return it->second;
}

bool DataExchange::hasKey(const std::string& key) const {
    return m_data.find(key) != m_data.end();
}

std::vector<std::string> DataExchange::getKeys() const {
    std::vector<std::string> keys;
    keys.reserve(m_data.size());
    for (const auto& pair : m_data) {
        keys.push_back(pair.first);
    }
    return keys;
}

void DataExchange::clear() {
    m_data.clear();
}

// Forward declarations of other integration interfaces
class EngineIntegration;
class FlowIntegration;
class CodeIntegration;
class DeskIntegration;
class ScribeIntegration;

// IntegrationFactory implementation
IntegrationFactory& IntegrationFactory::getInstance() {
    static IntegrationFactory instance;
    return instance;
}

std::shared_ptr<IntegrationInterface> IntegrationFactory::createInterface(Component component) {
    // Create and return the appropriate interface based on the component type
    
    switch (component) {
        case Component::REBELCAD:
            return std::make_shared<CADIntegration>();
        case Component::REBELENGINE:
            // return std::make_shared<EngineIntegration>();
            break;
        case Component::REBELFLOW:
            // return std::make_shared<FlowIntegration>();
            break;
        case Component::REBELCODE:
            // return std::make_shared<CodeIntegration>();
            break;
        case Component::REBELDESK:
            // return std::make_shared<DeskIntegration>();
            break;
        case Component::REBELSCRIBE:
            // return std::make_shared<ScribeIntegration>();
            break;
    }
    
    // Return nullptr for components that are not yet implemented
    return nullptr;
}

std::vector<Component> IntegrationFactory::getAvailableComponents() const {
    // This is a placeholder implementation
    // In a real implementation, we would check which components are actually available
    
    return {
        Component::REBELCAD,
        Component::REBELENGINE,
        Component::REBELFLOW,
        Component::REBELCODE,
        Component::REBELDESK,
        Component::REBELSCRIBE
    };
}

} // namespace Integration
} // namespace RebelCalc
