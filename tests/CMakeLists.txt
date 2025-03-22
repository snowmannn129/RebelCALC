# Test source files
set(TEST_SOURCES
    test_calculator.cpp
    test_symbolic_engine.cpp
    test_numeric_solver.cpp
    test_lua_engine.cpp
    test_engineering.cpp
)

# Create test executable
add_executable(RebelCALC_Tests ${TEST_SOURCES})

# Link libraries
target_link_libraries(RebelCALC_Tests PRIVATE 
    ${GTEST_LIBRARIES}
    ${LUA_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
)

# Add include directories
target_include_directories(RebelCALC_Tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${GTEST_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
)

# Set properties
set_target_properties(RebelCALC_Tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Add test
add_test(NAME RebelCALC_Tests COMMAND RebelCALC_Tests)
