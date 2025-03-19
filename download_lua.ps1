# PowerShell script to download and extract Lua for RebelCALC

# Configuration
$LUA_VERSION = "5.4.6"
$DOWNLOAD_URL = "https://www.lua.org/ftp/lua-$LUA_VERSION.tar.gz"
$DOWNLOAD_PATH = "external/downloads"
$EXTRACT_PATH = "external/lua-$LUA_VERSION"

# Create directories if they don't exist
Write-Host "Creating directories..."
New-Item -ItemType Directory -Force -Path $DOWNLOAD_PATH | Out-Null
New-Item -ItemType Directory -Force -Path $EXTRACT_PATH | Out-Null

# Download Lua
Write-Host "Downloading Lua $LUA_VERSION..."
$OUTPUT_FILE = "$DOWNLOAD_PATH/lua-$LUA_VERSION.tar.gz"

try {
    Invoke-WebRequest -Uri $DOWNLOAD_URL -OutFile $OUTPUT_FILE
    Write-Host "Download completed successfully."
} catch {
    Write-Host "Error downloading Lua: $_"
    exit 1
}

# Extract the archive
Write-Host "Extracting Lua..."
try {
    # Check if 7-Zip is available
    if (Get-Command "7z" -ErrorAction SilentlyContinue) {
        # Use 7-Zip to extract
        7z x $OUTPUT_FILE -o"$DOWNLOAD_PATH" -y
        7z x "$DOWNLOAD_PATH/lua-$LUA_VERSION.tar" -o"external" -y
    } else {
        # Use built-in tar command if available (Windows 10 1803+)
        if (Get-Command "tar" -ErrorAction SilentlyContinue) {
            tar -xzf $OUTPUT_FILE -C external
        } else {
            Write-Host "Error: Neither 7-Zip nor tar command is available. Please install 7-Zip or use Windows 10 1803+."
            exit 1
        }
    }
    Write-Host "Extraction completed successfully."
} catch {
    Write-Host "Error extracting Lua: $_"
    exit 1
}

# Build Lua (Windows)
Write-Host "Building Lua..."
try {
    $currentDir = Get-Location
    Set-Location "$EXTRACT_PATH/src"
    
    # Create a simple makefile for Windows with MSVC
    @"
# Makefile for building Lua with MSVC
CC= cl
CFLAGS= /O2 /W3 /MD /DLUA_BUILD_AS_DLL
LIBS= 
MYLIBS= 
MYOBJS= 

PLATS= 

LUA_A= lua.lib
LUA_T= lua.exe
LUAC_T= luac.exe

ALL_O= lapi.obj lcode.obj lctype.obj ldebug.obj ldo.obj ldump.obj lfunc.obj lgc.obj llex.obj \
	lmem.obj lobject.obj lopcodes.obj lparser.obj lstate.obj lstring.obj ltable.obj \
	ltm.obj lundump.obj lvm.obj lzio.obj \
	lauxlib.obj lbaselib.obj lcorolib.obj ldblib.obj liolib.obj \
	lmathlib.obj loslib.obj lstrlib.obj ltablib.obj lutf8lib.obj loadlib.obj linit.obj
CORE_O= lapi.obj lcode.obj lctype.obj ldebug.obj ldo.obj ldump.obj lfunc.obj lgc.obj llex.obj \
	lmem.obj lobject.obj lopcodes.obj lparser.obj lstate.obj lstring.obj ltable.obj \
	ltm.obj lundump.obj lvm.obj lzio.obj
LIB_O= lauxlib.obj lbaselib.obj lcorolib.obj ldblib.obj liolib.obj \
	lmathlib.obj loslib.obj lstrlib.obj ltablib.obj lutf8lib.obj loadlib.obj linit.obj
BASE_O= $(CORE_O) $(LIB_O) $(MYOBJS)

LUA_O= lua.obj
LUAC_O= luac.obj

ALL_T= $(LUA_A) $(LUA_T) $(LUAC_T)
ALL_A= $(LUA_A)

default: $(LUA_A)

$(LUA_A): $(BASE_O)
	$(AR) /OUT:$@ $(BASE_O) $(LIBS)

clean:
	del *.obj *.lib *.exp *.dll *.exe
"@ | Out-File -FilePath "Makefile.win" -Encoding ASCII

    # Build Lua
    if (Get-Command "nmake" -ErrorAction SilentlyContinue) {
        nmake -f Makefile.win
    } else {
        Write-Host "Error: nmake command is not available. Please install Visual Studio or Visual C++ Build Tools."
        exit 1
    }
    
    Set-Location $currentDir
    Write-Host "Lua build completed successfully."
} catch {
    Set-Location $currentDir
    Write-Host "Error building Lua: $_"
    exit 1
}

# Copy headers and libraries to the appropriate locations
Write-Host "Copying Lua headers and libraries..."
try {
    # Create include and lib directories
    New-Item -ItemType Directory -Force -Path "external/include/lua" | Out-Null
    New-Item -ItemType Directory -Force -Path "external/lib" | Out-Null
    
    # Copy headers
    Copy-Item "$EXTRACT_PATH/src/*.h" -Destination "external/include/lua" -Force
    
    # Copy library
    Copy-Item "$EXTRACT_PATH/src/lua.lib" -Destination "external/lib" -Force
    
    Write-Host "Lua headers and libraries copied successfully."
} catch {
    Write-Host "Error copying Lua headers and libraries: $_"
    exit 1
}

Write-Host "Lua $LUA_VERSION has been downloaded, built, and installed successfully."
Write-Host "Headers: external/include/lua"
Write-Host "Library: external/lib/lua.lib"
