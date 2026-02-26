# PRS Execution Runbook (WSL)

## Runtime policy

1. Preferred: container runtime (`docker`, then `apptainer`/`singularity`).
2. Fallback: `conda` only when containers are unavailable.
3. Never run `wsl` from inside WSL.

`workflow/prs/run_pgsc_calc.sh` auto-selects runtime by default.

## Open the correct shell

From Windows PowerShell:

```powershell
wsl -d Ubuntu
```

Then in Ubuntu:

```bash
cd /mnt/c/path/to/Genomics
```

## One-time checks

```bash
which nextflow
which docker || true
which apptainer || true
which singularity || true
which conda || true
```

## Run reliably (tmux)

```bash
tmux new -s prs
cd /mnt/c/path/to/Genomics
bash workflow/prs/run_pgsc_calc.sh
```

Detach without stopping run:

```bash
Ctrl+b d
```

Re-attach later:

```bash
tmux attach -t prs
```

## Monitor in a second shell

```bash
cd /mnt/c/path/to/Genomics
tail -f .nextflow.log
nextflow log
```

## If interrupted

The script supports resume; rerun the same command:

```bash
cd /mnt/c/path/to/Genomics
bash workflow/prs/run_pgsc_calc.sh
```

## Known failure patterns and fixes

1. `Command 'wsl' not found`
   - Cause: running `wsl` inside Linux shell.
   - Fix: use `wsl` only from Windows PowerShell/CMD.

2. Session lock error (`Unable to acquire lock on session...`)
   - Cause: stale/stopped Nextflow process.
   - Fix: stop old process, then rerun with resume.

3. `Process requirement exceeds available memory`
   - Fix: lower memory, for example:
   - `PGSC_MAX_MEMORY=12.GB bash workflow/prs/run_pgsc_calc.sh`

4. `chrX` sex / PAR errors
   - Built-in script options:
   - `--split-par b38`
   - Optional `config/prs_update_sex.tsv` with strict `FID IID SEX` format.

5. Conda solve mismatch for `pgsc_calc v2.2.0`
   - Affects only `-profile conda`.
   - Follow script guidance to patch cached env pin if needed.