# Find GTest package
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Test source files
set(TEST_SOURCES
    test_main.cpp
    test_calculator.cpp
    test_symbolic_engine.cpp
    test_numeric_solver.cpp
    test_lua_engine.cpp
    test_theme_manager.cpp
)

# Create test executable
add_executable(RebelCALC_tests ${TEST_SOURCES})

# Link libraries
target_link_libraries(RebelCALC_tests PRIVATE
    ${GTEST_LIBRARIES}
    ${LUA_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
    pthread
)

# Add include directories
target_include_directories(RebelCALC_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${GTEST_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
)

# Add tests
add_test(NAME RebelCALC_tests COMMAND RebelCALC_tests)

# Set properties
set_target_properties(RebelCALC_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
