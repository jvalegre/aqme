@echo off
:: ============================================================
:: labcode environment setup (Windows CMD)
:: Creates virtualenv, installs CUDA/CPU PyTorch + requirements
:: ============================================================

setlocal enabledelayedexpansion
title LabCode Environment Setup

:: --- CONFIGURATION ---
set ENV_NAME=lab_env
set PYTHON_VERSION=3.11
set REQUIREMENTS=requirements.txt

echo ============================================================
echo ?? Setting up environment: %ENV_NAME%
echo ============================================================

:: --- Step 1: Verify Python installation ---
where python >nul 2>nul
if %errorlevel% neq 0 (
    echo ? Python not found. Please install Python %PYTHON_VERSION% first.
    echo https://www.python.org/downloads/release/python-%PYTHON_VERSION:.=%
    pause
    exit /b 1
)
echo ? Python found: %PYTHON_VERSION% or later assumed.

:: --- Step 2: Create virtual environment (if missing) ---
if not exist "%ENV_NAME%" (
    echo Creating new virtual environment "%ENV_NAME%"...
    python -m venv "%ENV_NAME%"
    if %errorlevel% neq 0 (
        echo ? Failed to create virtual environment.
        exit /b 1
    )
) else (
    echo Environment "%ENV_NAME%" already exists.
)

:: --- Step 3: Activate environment ---
echo Activating "%ENV_NAME%"...
call "%ENV_NAME%\Scripts\activate"
if %errorlevel% neq 0 (
    echo ? Failed to activate virtual environment.
    exit /b 1
)
echo ? Virtual environment activated.

:: --- Step 4: Ensure pip/setuptools/wheel ---
echo Upgrading pip, setuptools, wheel...
python -m ensurepip
python -m pip install --upgrade pip setuptools wheel

:: --- Step 5: Check CUDA availability ---
echo Checking for NVIDIA GPU / CUDA...
nvidia-smi >nul 2>nul
if %errorlevel%==0 (
    echo ? NVIDIA GPU detected. Installing CUDA-enabled PyTorch...
    python -m pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
) else (
    echo ?? No GPU detected or CUDA missing. Installing CPU-only PyTorch...
    python -m pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
)

:: --- Step 6: Install requirements.txt packages ---
if exist "%REQUIREMENTS%" (
    echo Installing dependencies from "%REQUIREMENTS%"...
    pip install -r "%REQUIREMENTS%"
) else (
    echo ?? File "%REQUIREMENTS%" not found. Skipping package installation.
)

:: --- Step 7: Verify installation ---
echo.
echo ============================================================
echo ?? Verifying core packages...
python - <<EOF
import sys
print(f"Python {sys.version.split()[0]}")
try:
    import torch, rdkit, pandas, numpy
    print("? Core imports successful.")
except Exception as e:
    print(f"? Import error: {e}")
EOF

echo.
echo ============================================================
echo ? Setup complete!
echo To activate your environment later, run:
echo     call %ENV_NAME%\Scripts\activate
echo ============================================================
pause
