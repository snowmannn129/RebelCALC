cmake_minimum_required(VERSION 3.15)
project(RebelCALC VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external/include
)

# Link directories
link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/external/lib
)

# Find Lua package
find_package(Lua 5.4 REQUIRED)
include_directories(${LUA_INCLUDE_DIR})

# Find yaml-cpp package
find_package(yaml-cpp QUIET)
if(NOT yaml-cpp_FOUND)
    # If yaml-cpp is not found, use the one we built
    set(YAML_CPP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/include)
    set(YAML_CPP_LIBRARIES yaml-cpp)
    message(STATUS "Using bundled yaml-cpp")
else()
    message(STATUS "Using system yaml-cpp")
endif()

# Add subdirectories
add_subdirectory(src)

# Enable testing
enable_testing()
add_subdirectory(tests)

# Add examples
add_subdirectory(examples)

# Installation
install(TARGETS RebelCALC
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Install configuration files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/config/
    DESTINATION etc/RebelCALC
)

# Print configuration summary
message(STATUS "")
message(STATUS "RebelCALC Configuration Summary")
message(STATUS "------------------------------")
message(STATUS "Version:              ${PROJECT_VERSION}")
message(STATUS "C++ Compiler:         ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Build type:           ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix:       ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Lua version:          ${LUA_VERSION_STRING}")
message(STATUS "yaml-cpp:             ${yaml-cpp_FOUND}")
message(STATUS "")
