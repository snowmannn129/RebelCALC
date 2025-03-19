#!/bin/bash
# Script to run RebelCALC tests on Unix-like systems

echo "Running RebelCALC tests..."

# Check if build directory exists
if [ ! -d "build" ]; then
    echo "Build directory not found. Building RebelCALC first..."
    ./build.sh
    if [ $? -ne 0 ]; then
        echo "Build failed with error code $?"
        exit 1
    fi
fi

# Navigate to build directory
cd build

# Run the tests
echo "Running tests..."
ctest --output-on-failure

# Check if tests were successful
if [ $? -ne 0 ]; then
    echo "Tests failed with error code $?"
    cd ..
    exit 1
fi

echo "All tests passed successfully!"

# Return to original directory
cd ..

exit 0
