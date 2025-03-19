@echo off
REM Script to run RebelCALC tests on Windows

echo Running RebelCALC tests...

REM Check if build directory exists
if not exist build (
    echo Build directory not found. Building RebelCALC first...
    call build.bat
    if %ERRORLEVEL% NEQ 0 (
        echo Build failed with error code %ERRORLEVEL%
        exit /b %ERRORLEVEL%
    )
)

REM Navigate to build directory
cd build

REM Run the tests
echo Running tests...
ctest -C Release --output-on-failure

REM Check if tests were successful
if %ERRORLEVEL% NEQ 0 (
    echo Tests failed with error code %ERRORLEVEL%
    cd ..
    exit /b %ERRORLEVEL%
)

echo All tests passed successfully!

REM Return to original directory
cd ..

exit /b 0
