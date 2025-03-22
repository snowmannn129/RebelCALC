# Example source files
set(EXAMPLE_SOURCES
    engineering_demo.cpp
    integration_demo.cpp
    visualization_demo.cpp
    optimization_demo.cpp
    plot3d_demo.cpp
    plot2d_demo.cpp
    repl_demo.cpp
)

# Create example executables
foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_SOURCE} NAME_WE)
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_SOURCE})
    
    # Link libraries
    target_link_libraries(${EXAMPLE_NAME} PRIVATE 
        ${LUA_LIBRARIES}
        ${YAML_CPP_LIBRARIES}
    )
    
    # Add include directories
    target_include_directories(${EXAMPLE_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${YAML_CPP_INCLUDE_DIR}
    )
    
    # Set properties
    set_target_properties(${EXAMPLE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/examples
    )
    
    # Install
    install(TARGETS ${EXAMPLE_NAME}
        RUNTIME DESTINATION bin/examples
    )
endforeach()
