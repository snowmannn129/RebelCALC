# Source files
set(SOURCES
    main.cpp
    backend/calculator.cpp
    backend/symbolic_engine.cpp
    backend/numeric_solver.cpp
    ui/terminal_ui.cpp
    ui/theme_manager.cpp
    scripts/lua_engine.cpp
)

# Header files
set(HEADERS
    backend/calculator.h
    backend/symbolic_engine.h
    backend/numeric_solver.h
    ui/terminal_ui.h
    ui/theme_manager.h
    scripts/lua_engine.h
)

# Create executable
add_executable(RebelCALC ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(RebelCALC PRIVATE 
    ${LUA_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
)

# Define preprocessor macros
target_compile_definitions(RebelCALC PRIVATE
    REBELCALC_VERSION="${PROJECT_VERSION}"
    REBELCALC_CONFIG_DIR="${CMAKE_INSTALL_PREFIX}/etc/RebelCALC"
)

# Add include directories
target_include_directories(RebelCALC PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${YAML_CPP_INCLUDE_DIR}
)

# Set properties
set_target_properties(RebelCALC PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Install
install(TARGETS RebelCALC
    RUNTIME DESTINATION bin
)
