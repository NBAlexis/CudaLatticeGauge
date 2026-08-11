@echo off
rem runtest.bat TestName GpuCount -- multi-GPU-improve1.md 3.10 test driver.
rem Same contract as runtest.sh: rankCount/GpuGrid come from
rem `CLGTest --mg-config TestName` (no test-name table here); mpiexec exit
rem code is propagated unchanged. Overrides: CLGTEST_BIN, MPIEXEC.
setlocal EnableDelayedExpansion

if "%~2"=="" (
    echo usage: %~nx0 TestName GpuCount 1>&2
    exit /b 2
)
set "TEST_NAME=%~1"
set "GPU_COUNT=%~2"

echo(%GPU_COUNT%|findstr /r "^[1-9][0-9]*$" >nul
if errorlevel 1 (
    echo runtest: GpuCount must be a positive integer, got '%GPU_COUNT%' 1>&2
    exit /b 2
)

set "SCRIPT_DIR=%~dp0"
if "%SCRIPT_DIR:~-1%"=="\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
if not defined CLGTEST_BIN set "CLGTEST_BIN=%SCRIPT_DIR%\Bin\Ubuntu\CLGTest.exe"
if not defined MPIEXEC set "MPIEXEC=mpiexec"
if not exist "%CLGTEST_BIN%" (
    echo runtest: CLGTest binary not found: %CLGTEST_BIN% 1>&2
    exit /b 2
)

rem CLGTest loads the test YAML blocks via paths relative to Bin\Ubuntu.
pushd "%SCRIPT_DIR%\Bin\Ubuntu"
if errorlevel 1 (
    echo runtest: cannot cd %SCRIPT_DIR%\Bin\Ubuntu 1>&2
    exit /b 2
)

set "CFG_FILE=%TEMP%\clgtest_mgconfig_%RANDOM%%RANDOM%.tmp"
"%CLGTEST_BIN%" --mg-config "%TEST_NAME%" > "%CFG_FILE%"
if errorlevel 1 (
    del "%CFG_FILE%" >nul 2>&1
    echo runtest: --mg-config failed for %TEST_NAME% ^(unknown test, not _TEST_MULTIGPU, or invalid metadata^) 1>&2
    popd
    exit /b 1
)
set /p CONFIG=<"%CFG_FILE%"
del "%CFG_FILE%" >nul 2>&1

for /f "tokens=1-5" %%a in ("%CONFIG%") do (
    set "RANK_COUNT=%%a"
    set "GX=%%b"
    set "GY=%%c"
    set "GZ=%%d"
    set "GT=%%e"
)
if not defined GT (
    echo runtest: malformed --mg-config output for %TEST_NAME%: '%CONFIG%' 1>&2
    popd
    exit /b 1
)
echo(%RANK_COUNT%|findstr /r "^[1-9][0-9]*$" >nul
if errorlevel 1 (
    echo runtest: malformed rankCount in --mg-config output: '%CONFIG%' 1>&2
    popd
    exit /b 1
)
if %GPU_COUNT% GTR %RANK_COUNT% (
    echo runtest: GpuCount %GPU_COUNT% ^> rankCount %RANK_COUNT% for %TEST_NAME% 1>&2
    popd
    exit /b 2
)

echo runtest: %TEST_NAME% ranks=%RANK_COUNT% grid=[%GX%,%GY%,%GZ%,%GT%] devicePerNode=%GPU_COUNT% 1>&2
"%MPIEXEC%" -n %RANK_COUNT% "%CLGTEST_BIN%" "%TEST_NAME%" --mg-worker --gpu-grid "%GX%,%GY%,%GZ%,%GT%" --device-per-node %GPU_COUNT%
set "RC=%ERRORLEVEL%"
popd
exit /b %RC%
