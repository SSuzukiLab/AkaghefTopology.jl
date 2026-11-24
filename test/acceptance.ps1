# Daily acceptance run of Topology.jl: every tier up to acceptance, all modules.
# Registered as the Windows scheduled task "Topology acceptance (daily)"; policy in
# test/TESTING.md. Writes test/acceptance-logs/<date>.log, appends per-file times to
# test/acceptance-logs/timings.tsv, and writes the last result line to
# test/acceptance-logs/latest.txt. Skips the run if one is already in progress.

$ErrorActionPreference = "Stop"
$test = $PSScriptRoot
$pkg = Split-Path $test -Parent
$logs = Join-Path $test "acceptance-logs"
New-Item -ItemType Directory -Force $logs | Out-Null
$lock = Join-Path $logs "running.lock"
if (Test-Path $lock) {
    $age = (Get-Date) - (Get-Item $lock).LastWriteTime
    if ($age.TotalHours -lt 12) { exit 0 }
}
Set-Content $lock (Get-Date -Format o)
$stamp = Get-Date -Format "yyyy-MM-dd"
$log = Join-Path $logs "$stamp.log"
$env:TOPOLOGY_TEST_TIMINGS = Join-Path $logs "timings.tsv"
$julia = (Get-Command julia -ErrorAction Stop).Source
try {
    Push-Location $pkg
    $start = Get-Date
    & $julia --project=. test/run.jl --tier acceptance *> $log
    $code = $LASTEXITCODE
    Pop-Location
    $last = (Get-Content $log | Where-Object { $_ -match '^(ALL_PASS|FAILED)' } | Select-Object -Last 1)
    $mins = [math]::Round(((Get-Date) - $start).TotalMinutes, 1)
    Set-Content (Join-Path $logs "latest.txt") "$stamp exit=$code minutes=$mins $last"
} finally {
    Remove-Item $lock -Force -ErrorAction SilentlyContinue
}
exit $code
