@echo off
setlocal enabledelayedexpansion

echo ============================================
echo   Install & Check Project Dependencies
echo   MOSEK, MATLAB, WSL (Ubuntu) - DACE, CSPICE, Astrotools, nlhomann-Json
echo ============================================
echo.


REM ============================================
REM  Locate MOSEK automatically
REM ============================================
:MOSEK
echo Searching for MOSEK installation...

set FOUND_MOSEK=
set SEARCH_DIRS="C:\Program Files\Mosek" "C:\Program Files (x86)\Mosek" "%LOCALAPPDATA%\Mosek"

for %%D in (%SEARCH_DIRS%) do (
    if exist %%D (
        for /d %%V in (%%D\*) do (
            if exist "%%V\tools\platform\win64x86\bin\mosek.exe" (
                set FOUND_MOSEK=%%V\tools\platform\win64x86\bin\mosek.exe
            )
        )
    )
)

if "%FOUND_MOSEK%"=="" (
    echo [WARNING] MOSEK not found
    echo Please install MOSEK manually from https://www.mosek.com/
    echo.
    pause
    goto :MATLAB
)   
    echo MOSEK found at: %FOUND_MOSEK%
    echo Please, ensure you have a valid MOSEK license before proceeding.
    echo.
    
REM ============================================
REM  Locate MATLAB automatically
REM ============================================

:MATLAB
echo Searching for MATLAB installation...

set FOUND_MATLAB=
set MATLAB_SEARCH_DIRS="C:\Program Files\MATLAB" "C:\Program Files (x86)\MATLAB"

for %%D in (%MATLAB_SEARCH_DIRS%) do (
    if exist %%D (
        for /d %%V in (%%D\R*) do (
            if exist "%%V\bin\matlab.exe" (
                set FOUND_MATLAB=%%V\bin\matlab.exe
            )
        )
    )
)

if "%FOUND_MATLAB%"=="" (
    echo [WARNING] MATLAB not found
    echo Please install MATLAB manually from https://www.mathworks.com/
    echo.
    pause
    goto :WSL
) 
    echo MATLAB found at: %FOUND_MATLAB%
    echo Please, ensure you have a valid MATLAB license.
    echo.
    goto :WSL


REM ==========================================================
REM  Check WSL Ubuntu installation and CMAKE installation
REM ==========================================================
:WSL
set "ROOT=%~dp0"
set "LOG=%ROOT%install_log.txt"
del "%LOG%" >nul 2>&1

REM Detect WSL user/home
for /f "delims=" %%u in ('wsl -d Ubuntu echo $USER') do set "WSLUSER=%%u"
if "%WSLUSER%"=="" (
  echo [ERROR] Could not detect WSL user. >> "%LOG%"
  echo [ERROR] Could not detect WSL user. Ensure Ubuntu WSL is installed.
  pause
  exit /b 1
)
set "WSLHOME=\\wsl$\Ubuntu\home\%WSLUSER%"
echo Detected WSL user: %WSLUSER%
echo WSL home: %WSLHOME%
echo Log file: %LOG%
echo.

REM Helper macro for WSL execution
set "WSL=cmd /c wsl -d Ubuntu bash -lc"
goto :CMAKE


:CMAKE
%WSL% "sudo apt-get update" >> "%LOG%" 2>&1
%WSL% "sudo apt install cmake" >> "%LOG%" 2>&1
if "%ERRORLEVEL%"=="0" (
    echo   CMAKE installed
    echo.
) else (
    echo Manually install CMake before proceeding
    echo.
    pause
    exit /b
)

goto :DACE


REM ==========================================================
REM ==    DACE detection & installation (Linux procedure)
REM ==    Requires -DWITH_ALGEBRAICMATRIX=ON for Astrotools
REM ==========================================================

:DACE
echo ==== Checking DACE (/usr/local) ====

REM Header
%WSL% "sh -c 'test -f /usr/local/include/dace/dace.h || test -f /usr/include/dace/dace.h'" >> "%LOG%" 2>&1
set DACE_HDR=%ERRORLEVEL%

REM Library
%WSL% "sh -c '[ -f /usr/local/lib/libdace_s.a ] || [ -f /usr/local/lib/libdace.so ] || [ -f /usr/local/lib/libdace.a ]'" >> "%LOG%" 2>&1
set DACE_LIB=%ERRORLEVEL%

if "%DACE_HDR%"=="0" if "%DACE_LIB%"=="0" (
    echo   DACE already installed
    echo.
    goto :CSPICE
)  
    echo   DACE not found 
    choice /M "Install DACE now?"
    if errorlevel 2 (
        echo Skipped DACE installation by user.
        echo.
        goto :CSPICE
    )

    REM Remove and clone fresh
    %WSL% "cd /home/%WSLUSER% && git clone https://github.com/dacelib/dace.git dace" >> "%LOG%" 2>&1

    REM Configure (with algebraic matrix support)
    %WSL% "cmake -S /home/%WSLUSER%/dace -B /home/%WSLUSER%/dace-build -DWITH_ALGEBRAICMATRIX=ON" >> "%LOG%" 2>&1
    if errorlevel 1 goto :DACE_FAIL

    REM Build
    %WSL% "cmake --build /home/%WSLUSER%/dace-build -j" >> "%LOG%" 2>&1
    if errorlevel 1 goto :DACE_FAIL

    REM Install into /usr/local
    %WSL% "sudo cmake --install /home/%WSLUSER%/dace-build" >> "%LOG%" 2>&1
    if errorlevel 1 goto :DACE_FAIL

    %WSL% "sudo ldconfig" >> "%LOG%"

    echo   DACE installed
    echo.
    goto :CSPICE

:DACE_FAIL
echo [ERROR] DACE installation failed. See %LOG%
pause
exit /b 10

REM ==========================================================
REM == 2) CSPICE detection & installation (Linux toolkit)
REM ==    Download cspice.tar.Z (Linux GCC 64‑bit)
REM ==    Install headers/libs into /usr/local
REM ==========================================================
:CSPICE
echo ==== Checking CSPICE (/usr/local) ====

REM Headers (allow /usr/local/include or /usr/local/include/cspice)
%WSL% "sh -lc 'test -f /usr/include/SpiceUsr.h || test -f /usr/include/cspice/SpiceUsr.h || test -f /usr/local/include/SpiceUsr.h || test -f /usr/local/include/cspice/SpiceUsr.h'" >> "%LOG%" 2>&1
set CSPICE_HDR=%ERRORLEVEL%

REM Libraries
%WSL% "sh -lc '[ -f /usr/local/lib/cspice.a ] || [ -f /usr/local/lib/libcspice.so ] || [ -f /usr/local/lib/csupport.a ]'" >> "%LOG%" 2>&1
set CSPICE_LIB=%ERRORLEVEL%

if "%CSPICE_HDR%"=="0" if "%CSPICE_LIB%"=="0" (
    echo   CSPICE already installed
    echo.
    goto :ASTROTOOLS
) 
    echo   CSPICE not found 
    choice /M "Download and install CSPICE now?"
    if errorlevel 2 (
        echo Skipped CSPICE installation by user.
        echo.
        goto :ASTROTOOLS
    )

    REM Download Linux GCC64 package
    echo   Downloading Linux CSPICE package...
    %WSL% "mkdir -p /home/%WSLUSER%/naif && cd /home/%WSLUSER%/naif && rm -f cspice.tar* && curl -L -O https://naif.jpl.nasa.gov/pub/naif/toolkit/C/PC_Linux_GCC_64bit/packages/cspice.tar.Z" >> "%LOG%" 2>&1

    REM Extract (.Z requires `uncompress`)
    %WSL% "cd /home/%WSLUSER%/naif && uncompress -f cspice.tar.Z && tar xf cspice.tar" >> "%LOG%" 2>&1

    REM Copy include + lib to system
    echo   Installing CSPICE headers/libs into /usr/local...
    %WSL% "sudo mkdir -p /usr/local/include/cspice" >> "%LOG%"
    %WSL% "sudo cp -r /home/%WSLUSER%/naif/cspice/include/* /usr/local/include/cspice/" >> "%LOG%"
    %WSL% "if [ -d /home/%WSLUSER%/naif/cspice/lib ]; then sudo cp -r /home/%WSLUSER%/naif/cspice/lib/* /usr/local/lib/; fi" >> "%LOG%"
    %WSL% "if [ -d /home/%WSLUSER%/naif/cspice/lib64 ]; then sudo cp -r /home/%WSLUSER%/naif/cspice/lib64/* /usr/local/lib/; fi" >> "%LOG%"

    REM Symlink SpiceUsr.h for convenience
    %WSL% "sudo ln -sf /usr/local/include/cspice/SpiceUsr.h /usr/local/include/SpiceUsr.h" >> "%LOG%"

    %WSL% "sudo ldconfig" >> "%LOG%"
    echo   CSPICE installed 
    echo.

goto :ASTROTOOLS

REM ==========================================================
REM == 3) ASTROTOOLS detection, extraction, dependencies, build
REM ==    Extract from astrotools.zip in parent directory
REM ==    Install headers/libs into /usr/local
REM ==========================================================
:ASTROTOOLS

echo ==== Checking Astrotools installation ====

REM Header check (e.g., dynorb/AIDA.h)
%WSL% "sh -lc 'test -f /usr/local/include/dynorb/AIDA.h || test -f /usr/include/dynorb/AIDA.h'" >> "%LOG%" 2>&1
set AT_HDR=%ERRORLEVEL%

REM Library check (astro, geco)
%WSL% "sh -lc '[ -f /usr/local/lib/libastro.so ] || [ -f /usr/local/lib/libgeco.so ]'" >> "%LOG%" 2>&1
set AT_LIB=%ERRORLEVEL%

if "%AT_HDR%"=="0" if "%AT_LIB%"=="0" (
    echo   Astrotools is already installed
    echo.
    goto :JSON
)

echo   Astrotools not fully installed 
choice /M "Extract and install Astrotools now?"
if errorlevel 2 (
        echo Skipped Astrotools installation by user.
        echo.
        goto :JSON
    )
    REM Check if astrotools.zip exists in parent directory
    set "ASTROTOOLS_ZIP=%ROOT%astrotools.zip"
    if not exist "%ASTROTOOLS_ZIP%" (
        echo [ERROR] astrotools.zip not found at %ASTROTOOLS_ZIP%
        pause
        exit /b 20
    )

    REM Extract astrotools.zip to WSL home
    echo   Extracting Astrotools from zip...
    %WSL% "unzip -q -o /mnt/c/$(echo '%ASTROTOOLS_ZIP:C:\=\%' | sed 's|\\|/|g') -d /home/%WSLUSER%" >> "%LOG%" 2>&1
    if errorlevel 1 (
        echo [ERROR] Failed to extract astrotools.zip
        pause
        exit /b 21
    )

    REM Install prerequisites via apt
    echo   Installing Astrotools build prerequisites (Eigen3, JSONcpp, DLib)...
    %WSL% "sudo apt-get install -y libeigen3-dev libjsoncpp-dev libdlib-dev" >> "%LOG%" 2>&1

    REM Build + install
    echo   Building Astrotools...
    %WSL% "cd /home/%WSLUSER%/astrotools && cmake -S . -B _build" >> "%LOG%" 2>&1
    if errorlevel 1 goto :ASTROTOOLS_FAIL

    %WSL% "cmake --build /home/%WSLUSER%/astrotools/_build -j" >> "%LOG%" 2>&1
    if errorlevel 1 goto :ASTROTOOLS_FAIL

    echo   Installing Astrotools into /usr/local...
    %WSL% "sudo cmake --install /home/%WSLUSER%/astrotools/_build" >> "%LOG%" 2>&1
    if errorlevel 1 goto :ASTROTOOLS_FAIL

    %WSL% "sudo ldconfig" >> "%LOG%"
    echo   Astrotools installed 
    goto :JSON

:ASTROTOOLS_FAIL
echo [ERROR] Astrotools installation failed. See %LOG%
pause
exit /b 22

REM ==========================================================
REM == 4) nlohmann-Json detection and installation
REM ==========================================================
:JSON
REM 2026: Code Improved using Claude (Sonnet 4.6)
echo ==== Checking nlohmann-Json ====

REM Header check (allow /usr/include or /usr/local/include)
%WSL% "sh -lc 'test -f /usr/include/nlohmann/json.hpp || test -f /usr/local/include/nlohmann/json.hpp'" >> "%LOG%" 2>&1
set JSON_HDR=%ERRORLEVEL%

if "%JSON_HDR%"=="0" (
    echo   nlohmann-Json is already installed
    echo.
    goto :DONE
) 
    echo   nlohmann-Json not found 
    choice /M "Install nlohmann-Json now?"
    if errorlevel 2 (
        echo Skipped nlohmann-Json installation by user.
        echo.
        goto :DONE
    )
    
    echo   Installing nlohmann-Json...
    %WSL% "sudo apt-get install -y nlohmann-json3-dev" >> "%LOG%" 2>&1
    if errorlevel 1 (
        echo [ERROR] Failed to install nlohmann-Json
        echo.
        pause
        exit /b 30
    )
    echo   nlohmann-Json installed 
    echo.

:DONE
echo ============================================
echo   All dependencies checked and installed 
echo   See log: %LOG%
echo ============================================
pause
exit /b 0

