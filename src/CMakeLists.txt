# Source files
set(SOURCES
    main.cpp
    backend/calculator.cpp
    backend/symbolic_engine.cpp
    backend/numeric_solver.cpp
    backend/complex.cpp
    backend/matrix.cpp
    backend/statistics.cpp
    backend/units.cpp
    ui/terminal_ui.cpp
    ui/theme_manager.cpp
    ui/input_processor.cpp
    scripts/lua_engine.cpp
    engineering/cad.cpp
    engineering/physics.cpp
    engineering/electrical.cpp
    integration/api.cpp
    integration/cad.cpp
    visualization/plot2d.cpp
    visualization/plot3d.cpp
)

# Header files
set(HEADERS
    backend/calculator.h
    backend/symbolic_engine.h
    backend/numeric_solver.h
    backend/complex.h
    backend/matrix.h
    backend/statistics.h
    backend/units.h
    ui/terminal_ui.h
    ui/theme_manager.h
    ui/input_processor.h
    scripts/lua_engine.h
    engineering/cad.h
    engineering/physics.h
    engineering/electrical.h
    integration/api.h
    integration/cad.h
    visualization/plot2d.h
    visualization/plot3d.h
    visualization/datavis.h
    solvers/fem.h
    solvers/cfd.h
    solvers/circuit.h
    solvers/optimization.h
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
