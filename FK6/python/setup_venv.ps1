$ErrorActionPreference = "Stop"

$pythonRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$venvPath = Join-Path $pythonRoot ".venv"
$requirements = Join-Path $pythonRoot "requirements.txt"

if (-not (Test-Path (Join-Path $venvPath "Scripts\python.exe"))) {
    Write-Host "Creating project virtual environment..."
    python -m venv $venvPath
}

$venvPython = Join-Path $venvPath "Scripts\python.exe"
& $venvPython -m pip install --upgrade pip
& $venvPython -m pip install -r $requirements
Write-Host "Project environment is ready: $venvPython"
