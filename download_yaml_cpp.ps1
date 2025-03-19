# PowerShell script to download and extract yaml-cpp for RebelCALC

# Configuration
$YAML_CPP_VERSION = "0.8.0"
$DOWNLOAD_URL = "https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-$YAML_CPP_VERSION.zip"
$DOWNLOAD_PATH = "external/downloads"
$EXTRACT_PATH = "external/yaml-cpp-$YAML_CPP_VERSION"

# Create directories if they don't exist
Write-Host "Creating directories..."
New-Item -ItemType Directory -Force -Path $DOWNLOAD_PATH | Out-Null
New-Item -ItemType Directory -Force -Path "external" | Out-Null

# Download yaml-cpp
Write-Host "Downloading yaml-cpp $YAML_CPP_VERSION..."
$OUTPUT_FILE = "$DOWNLOAD_PATH/yaml-cpp-$YAML_CPP_VERSION.zip"

try {
    New-Item -ItemType Directory -Force -Path $DOWNLOAD_PATH | Out-Null
    Invoke-WebRequest -Uri $DOWNLOAD_URL -OutFile $OUTPUT_FILE
    Write-Host "Download completed successfully."
} catch {
    Write-Host "Error downloading yaml-cpp: $_"
    exit 1
}

# Extract the archive
Write-Host "Extracting yaml-cpp..."
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
    Write-Host "Error extracting yaml-cpp: $_"
    exit 1
}

# Build yaml-cpp
Write-Host "Building yaml-cpp..."
try {
    # Create build directory
    New-Item -ItemType Directory -Force -Path "$EXTRACT_PATH/build" | Out-Null
    
    # Navigate to build directory
    $currentDir = Get-Location
    Set-Location "$EXTRACT_PATH/build"
    
    # Configure with CMake
    cmake .. -G "Visual Studio 17 2022" -A x64 -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
    
    # Build
    cmake --build . --config Release
    
    # Return to original directory
    Set-Location $currentDir
    
    Write-Host "yaml-cpp build completed successfully."
} catch {
    # Return to original directory
    Set-Location $currentDir
    Write-Host "Error building yaml-cpp: $_"
    exit 1
}

# Copy headers and libraries to the appropriate locations
Write-Host "Copying yaml-cpp headers and libraries..."
try {
    # Create include and lib directories
    New-Item -ItemType Directory -Force -Path "external/include/yaml-cpp" | Out-Null
    New-Item -ItemType Directory -Force -Path "external/lib" | Out-Null
    
    # Copy headers
    Copy-Item "$EXTRACT_PATH/include/yaml-cpp/*" -Destination "external/include/yaml-cpp" -Recurse -Force
    
    # Copy libraries
    Copy-Item "$EXTRACT_PATH/build/Release/*" -Destination "external/lib" -Force
    
    Write-Host "yaml-cpp headers and libraries copied successfully."
} catch {
    Write-Host "Error copying yaml-cpp headers and libraries: $_"
    exit 1
}

Write-Host "yaml-cpp $YAML_CPP_VERSION has been downloaded, built, and installed successfully."
Write-Host "Headers: external/include/yaml-cpp"
Write-Host "Libraries: external/lib"
