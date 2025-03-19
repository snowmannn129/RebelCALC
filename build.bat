@echo off
REM Build script for RebelCALC on Windows

echo Building RebelCALC...

REM Download and build dependencies if they don't exist
if not exist external\lib\lua.lib (
    echo Downloading and building Lua...
    powershell -ExecutionPolicy Bypass -File download_lua.ps1
    if %ERRORLEVEL% NEQ 0 (
        echo Failed to download and build Lua with error code %ERRORLEVEL%
        exit /b %ERRORLEVEL%
    )
)

if not exist external\lib\yaml-cpp.lib (
    echo Downloading and building yaml-cpp...
    powershell -ExecutionPolicy Bypass -File download_yaml_cpp.ps1
    if %ERRORLEVEL% NEQ 0 (
        echo Failed to download and build yaml-cpp with error code %ERRORLEVEL%
        exit /b %ERRORLEVEL%
    )
)

if not exist external\lib\gtest.lib (
    echo Downloading and building GoogleTest...
    powershell -ExecutionPolicy Bypass -File download_gtest.ps1
    if %ERRORLEVEL% NEQ 0 (
        echo Failed to download and build GoogleTest with error code %ERRORLEVEL%
        exit /b %ERRORLEVEL%
    )
)

REM Create build directory if it doesn't exist
if not exist build mkdir build

REM Navigate to build directory
cd build

REM Configure with CMake
echo Configuring with CMake...
cmake .. -G "Visual Studio 17 2022" -A x64

REM Build the project
echo Building...
cmake --build . --config Release

REM Check if build was successful
if %ERRORLEVEL% NEQ 0 (
    echo Build failed with error code %ERRORLEVEL%
    cd ..
    exit /b %ERRORLEVEL%
)

echo Build completed successfully!
echo Executable is located at: build\bin\Release\RebelCALC.exe

REM Return to original directory
cd ..

exit /b 0
