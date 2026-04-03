@echo off
setlocal enabledelayedexpansion

echo ============================================
echo       BUILD WITH WSL AND CMAKE
echo ============================================
echo.

echo Starting WSL build...

set WIN_PATH=%cd%
for /f "delims=" %%a in ('wsl wslpath "%WIN_PATH%"') do set WSL_PATH=%%a

wsl bash -c " cd '%WSL_PATH%' && rm -rf build && cmake -S . -B build && cmake --build build --parallel"

if errorlevel 1 (
    echo.
    echo ERROR: Build failed inside WSL.
    pause
    exit /b 1
)

pause

