# PowerShell script to download and extract GoogleTest for RebelCALC

# Configuration
$GTEST_VERSION = "1.14.0"
$DOWNLOAD_URL = "https://github.com/google/googletest/archive/refs/tags/v$GTEST_VERSION.zip"
$DOWNLOAD_PATH = "external/downloads"
$EXTRACT_PATH = "external/googletest-$GTEST_VERSION"

# Create directories if they don't exist
Write-Host "Creating directories..."
New-Item -ItemType Directory -Force -Path $DOWNLOAD_PATH | Out-Null
New-Item -ItemType Directory -Force -Path "external" | Out-Null

# Download GoogleTest
Write-Host "Downloading GoogleTest $GTEST_VERSION..."
$OUTPUT_FILE = "$DOWNLOAD_PATH/googletest-$GTEST_VERSION.zip"

try {
    New-Item -ItemType Directory -Force -Path $DOWNLOAD_PATH | Out-Null
    Invoke-WebRequest -Uri $DOWNLOAD_URL -OutFile $OUTPUT_FILE
    Write-Host "Download completed successfully."
} catch {
    Write-Host "Error downloading GoogleTest: $_"
    exit 1
}

# Extract the archive
Write-Host "Extracting GoogleTest..."
try {
    # Check if 7-Zip is available
    if (Get-Command "7z" -ErrorAction SilentlyContinue) {
        # Use 7-Zip to extract
        7z x $OUTPUT_FILE -o"external" -y
    } else {
        # Use built-in Expand-Archive if available (PowerShell 5.0+)
        Expand-Archive -Path $OUTPUT_FILE -DestinationPath "external" -Force
    }
    Write-Host "Extraction completed successfully."
} catch {
    Write-Host "Error extracting GoogleTest: $_"
    exit 1
}

# Build GoogleTest
Write-Host "Building GoogleTest..."
try {
    # Create build directory
    New-Item -ItemType Directory -Force -Path "$EXTRACT_PATH/build" | Out-Null
    
    # Navigate to build directory
    $currentDir = Get-Location
    Set-Location "$EXTRACT_PATH/build"
    
    # Configure with CMake
    cmake .. -G "Visual Studio 17 2022" -A x64 -DBUILD_GMOCK=ON
    
    # Build
    cmake --build . --config Release
    
    # Return to original directory
    Set-Location $currentDir
    
    Write-Host "GoogleTest build completed successfully."
} catch {
    # Return to original directory
    Set-Location $currentDir
    Write-Host "Error building GoogleTest: $_"
    exit 1
}

# Copy headers and libraries to the appropriate locations
Write-Host "Copying GoogleTest headers and libraries..."
try {
    # Create include and lib directories
    New-Item -ItemType Directory -Force -Path "external/include/gtest" | Out-Null
    New-Item -ItemType Directory -Force -Path "external/include/gmock" | Out-Null
    New-Item -ItemType Directory -Force -Path "external/lib" | Out-Null
    
    # Copy headers
    Copy-Item "$EXTRACT_PATH/googletest/include/gtest/*" -Destination "external/include/gtest" -Recurse -Force
    Copy-Item "$EXTRACT_PATH/googlemock/include/gmock/*" -Destination "external/include/gmock" -Recurse -Force
    
    # Copy libraries
    Copy-Item "$EXTRACT_PATH/build/lib/Release/*" -Destination "external/lib" -Force
    
    Write-Host "GoogleTest headers and libraries copied successfully."
} catch {
    Write-Host "Error copying GoogleTest headers and libraries: $_"
    exit 1
}

Write-Host "GoogleTest $GTEST_VERSION has been downloaded, built, and installed successfully."
Write-Host "Headers: external/include/gtest and external/include/gmock"
Write-Host "Libraries: external/lib"
