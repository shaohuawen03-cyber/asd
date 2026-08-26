# sync_results_to_github.ps1
# Push ONLY the small analysis/result files (figures + xvg/csv/dat/txt/ndx/log/mdp/top)
# from the four MD result folders to GitHub branch "results-sync".
# Big MD binaries (md.xtc / md.tpr / *.edr / *.cpt / *.gro / *.trr) are NEVER uploaded.
# Does NOT touch the running mdrun jobs. Safe to re-run any time (e.g. after
# ylsllqr / ache finish their 100 ns MD, just run this script again).

$ErrorActionPreference = "Stop"

# 1. Go to the repo root (where this script lives)
Set-Location (Split-Path -Parent $MyInvocation.MyCommand.Path)

$smallExts = @(".png", ".pdf", ".svg", ".xvg", ".csv", ".tsv", ".dat", ".txt", ".ndx", ".log", ".mdp", ".top", ".itp")
$mdDirs    = @("gromacs_md\md_alllhrc", "gromacs_md\md_fllhttr", "gromacs_md\md_ylsllqr", "gromacs_md\md_ache")

$toAdd = New-Object System.Collections.Generic.List[string]

foreach ($d in $mdDirs) {
    if (Test-Path $d) {
        Get-ChildItem -Path $d -Recurse -File -ErrorAction SilentlyContinue |
            Where-Object { ($smallExts -contains $_.Extension.ToLowerInvariant()) -and -not ($_.Name -like "#*") } |
            ForEach-Object { $toAdd.Add($_.FullName) }
        Write-Host ("[SCAN] {0}: found {1} files" -f $d, (($toAdd | Where-Object { $_ -like "*$d*" }) | Measure-Object).Count)
    } else {
        Write-Host "[SKIP] folder not found: $d"
    }
}

# input PDBs (structure sources, small)
if (Test-Path "input") {
    Get-ChildItem -Path "input" -File -Filter *.pdb -ErrorAction SilentlyContinue |
        ForEach-Object { $toAdd.Add($_.FullName) }
}

# mdp templates (small text files)
if (Test-Path "mdp_templates") {
    Get-ChildItem -Path "mdp_templates" -Recurse -File -Include *.mdp -ErrorAction SilentlyContinue |
        ForEach-Object { $toAdd.Add($_.FullName) }
}

if ($toAdd.Count -eq 0) {
    Write-Host "[ABORT] no small result files found to sync."
    exit 1
}

Write-Host ("[1/3] staging {0} small result files (NO xtc/tpr/edr/cpt/gro/trr) ..." -f $toAdd.Count)
& git add -- $toAdd.ToArray()
if ($LASTEXITCODE -ne 0) { Write-Host "[FAIL] git add failed"; exit 1 }

# 2. Commit (only if something is staged)
& git diff --cached --quiet
if ($LASTEXITCODE -eq 0) {
    Write-Host "[2/3] nothing new to commit"
} else {
    $msg = "Sync real MD results (figures + small analysis outputs) " + (Get-Date -Format "yyyy-MM-dd HH:mm")
    & git commit -m $msg
    if ($LASTEXITCODE -ne 0) { Write-Host "[FAIL] git commit failed"; exit 1 }
    Write-Host "[2/3] committed: $msg"
}

# 3. Push current HEAD to GitHub branch "results-sync" (created on first run)
Write-Host "[3/3] pushing to origin results-sync ..."
& git push origin HEAD:refs/heads/results-sync
if ($LASTEXITCODE -ne 0) { Write-Host "[FAIL] push failed - check GitHub connection"; exit 1 }

Write-Host ""
Write-Host "[DONE] results are on GitHub branch: results-sync"
Write-Host "       Tell the analysis side to fetch, then it will verify and merge into arena/019ff90e-asd."
