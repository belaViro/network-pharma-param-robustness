param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ReproduceArgs
)

$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

$python = Get-Command python -ErrorAction SilentlyContinue
if (-not $python) {
    $python = Get-Command py -ErrorAction SilentlyContinue
}
if (-not $python) {
    throw "Python >= 3.10 is required but was not found on PATH."
}

& $python.Source .\reproduce.py @ReproduceArgs
exit $LASTEXITCODE
