# usage (from lib/topology.jl): pwsh scripts/build_conventions.ps1
# Renders every module convention page docs/conventions/<Module>.md to <Module>.html beside it with pandoc.
# The Markdown page is the canonical text; the HTML page is regenerated from it and never edited by hand.
$ErrorActionPreference = 'Stop'
$dir = Join-Path $PSScriptRoot '..\docs\conventions'
# The `pandoc` on PATH may be a bash wrapper that PowerShell cannot run; prefer the Windows executable.
$pandoc = @("$env:LOCALAPPDATA\Pandoc\pandoc.exe", "$env:ProgramFiles\Pandoc\pandoc.exe") | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $pandoc) { $pandoc = (Get-Command pandoc.exe -ErrorAction Stop).Source }
Push-Location $dir
try {
    foreach ($md in Get-ChildItem -Filter *.md) {
        $html = [IO.Path]::ChangeExtension($md.Name, '.html')
        & $pandoc -f gfm -t html5 -s --embed-resources --css=assets/page.css --metadata pagetitle="$($md.BaseName) の規約" -o $html $md.Name
        if ($LASTEXITCODE -ne 0) { throw "pandoc failed on $($md.Name)" }
        "wrote $html"
    }
} finally { Pop-Location }
