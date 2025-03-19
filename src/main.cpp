#include <iostream>
#include <string>
#include <memory>

#include "backend/calculator.h"
#include "ui/terminal_ui.h"
#include "scripts/lua_engine.h"

void printBanner() {
    std::cout << "\n";
    std::cout << "  ██████╗ ███████╗██████╗ ███████╗██╗      ██████╗ █████╗ ██╗      ██████╗\n";
    std::cout << "  ██╔══██╗██╔════╝██╔══██╗██╔════╝██║     ██╔════╝██╔══██╗██║     ██╔════╝\n";
    std::cout << "  ██████╔╝█████╗  ██████╔╝█████╗  ██║     ██║     ███████║██║     ██║     \n";
    std::cout << "  ██╔══██╗██╔══╝  ██╔══██╗██╔══╝  ██║     ██║     ██╔══██║██║     ██║     \n";
    std::cout << "  ██║  ██║███████╗██████╔╝███████╗███████╗╚██████╗██║  ██║███████╗╚██████╗\n";
    std::cout << "  ╚═╝  ╚═╝╚══════╝╚═════╝ ╚══════╝╚══════╝ ╚═════╝╚═╝  ╚═╝╚══════╝ ╚═════╝\n";
    std::cout << "                                                                           \n";
    std::cout << "  Advanced Computational Engine for RebelSUITE - v" << REBELCALC_VERSION << "\n";
    std::cout << "  Type 'help' for a list of commands or 'exit' to quit\n";
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    // Create components
    auto calculator = std::make_shared<rebelcalc::Calculator>();
    auto luaEngine = std::make_shared<rebelcalc::LuaEngine>(calculator);
    auto ui = std::make_shared<rebelcalc::TerminalUI>(calculator, luaEngine);
    
    // Initialize components
    if (!calculator->initialize()) {
        std::cerr << "Failed to initialize calculator engine" << std::endl;
        return 1;
    }
    
    if (!luaEngine->initialize()) {
        std::cerr << "Failed to initialize Lua scripting engine" << std::endl;
        return 1;
    }
    
    if (!ui->initialize()) {
        std::cerr << "Failed to initialize user interface" << std::endl;
        return 1;
    }
    
    // Print welcome banner
    printBanner();
    
    // Run the main application loop
    ui->run();
    
    // Cleanup
    ui->shutdown();
    luaEngine->shutdown();
    calculator->shutdown();
    
    return 0;
}
