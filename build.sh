#!/bin/bash
# Build script for RebelCALC on Unix-like systems

echo "Building RebelCALC..."

# Check for required tools
command -v cmake >/dev/null 2>&1 || { echo "CMake is required but not installed. Aborting."; exit 1; }

# Check for Lua
if [ ! -f "/usr/include/lua5.4/lua.h" ] && [ ! -f "/usr/local/include/lua5.4/lua.h" ] && [ ! -f "external/include/lua/lua.h" ]; then
    echo "Lua 5.4 is required but not found. Please install it using your package manager."
    echo "For example:"
    echo "  Ubuntu/Debian: sudo apt-get install liblua5.4-dev"
    echo "  macOS with Homebrew: brew install lua"
    exit 1
fi

# Check for yaml-cpp
if [ ! -f "/usr/include/yaml-cpp/yaml.h" ] && [ ! -f "/usr/local/include/yaml-cpp/yaml.h" ] && [ ! -f "external/include/yaml-cpp/yaml.h" ]; then
    echo "yaml-cpp is required but not found. Please install it using your package manager."
    echo "For example:"
    echo "  Ubuntu/Debian: sudo apt-get install libyaml-cpp-dev"
    echo "  macOS with Homebrew: brew install yaml-cpp"
    exit 1
fi

# Check for GoogleTest (optional)
if [ ! -f "/usr/include/gtest/gtest.h" ] && [ ! -f "/usr/local/include/gtest/gtest.h" ] && [ ! -f "external/include/gtest/gtest.h" ]; then
    echo "GoogleTest is not found. Tests will not be built."
    echo "To install GoogleTest:"
    echo "  Ubuntu/Debian: sudo apt-get install libgtest-dev"
    echo "  macOS with Homebrew: brew install googletest"
fi

# Create build directory if it doesn't exist
mkdir -p build

# Navigate to build directory
cd build

# Configure with CMake
echo "Configuring with CMake..."
cmake ..

# Build the project
echo "Building..."
if [ "$(uname)" == "Darwin" ]; then
    # macOS
    cmake --build . -- -j$(sysctl -n hw.ncpu)
else
    # Linux and other Unix-like systems
    cmake --build . -- -j$(nproc)
fi

# Check if build was successful
if [ $? -ne 0 ]; then
    echo "Build failed with error code $?"
    cd ..
    exit 1
fi

echo "Build completed successfully!"
echo "Executable is located at: build/bin/RebelCALC"

# Return to original directory
cd ..

exit 0
