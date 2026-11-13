# ORFscape — project instructions

This repository is a proteogenomics pipeline (Python 3.9+, no build step). Run
everything from the repository root with `python -m orfscape <command>`.

## Ground rules

- **Numbers, not narration.** After every stage, run the checker and report the
  numbers it printed. A stage is not done until its check passes.
- **Stop on the first silent failure** — empty database, missing layer, zero PSM,
  unlocated event. Do not "try anyway" and do not continue to see what happens.
- **Never re-run a completed stage** to make sure it worked. Use `--force` only when
  an upstream input changed or the user asked.
- **Never touch** `orfscape_projects/<proj>/{P13_results,02_db,02_tx,00_logs}` of an
  existing project unless the user explicitly asks: those are frozen results.
- Ask before deleting a stage directory or overwriting `novel_validated.tsv`; move the
  old run to `_staging/<name>_<date>` instead.

## Commands

```bash
python -m orfscape doctor --workdir orfscape_projects/envcheck    # environment
python -m orfscape gui                                           # desktop GUI
python -m orfscape run --workdir <proj> --raw <raw> [--from P0 --to P8c] [--skip P6b]
python -m orfscape eukbuild --workdir <proj> --preset HSAP       # eukaryotic layers
python integrations/scripts/preflight.py --workdir <proj> --raw <raw>
python integrations/scripts/stage_check.py --workdir <proj> [--stage P6]
python integrations/scripts/run_pipeline.py --workdir <proj> --raw <raw>
```

Tests (run before claiming anything works):

```bash
python tests/test_gui_core.py && python tests/test_gui_state.py && python tests/gui_theme.py
python tests/gui_visual.py && python tests/gui_run.py && python tests/test_gui_lint.py
python tests/test_eukbuild.py
```

## Conventions that bite

- Layer whitelist covers `L1_`, `L2_`, `L3|`, `L6|`. A non-empty layer with zero novel
  peptides is a bug, not a result.
- pFind `Positions` is 0-based; aa→nt is `pos*3` (never `(pos-1)*3`).
- Enzyme tags match `TL` / `LL` exactly — never `"L" in "TL"`.
- Event classification uses the in-frame overlapping CDS (`rel % 3 == 0`), and the
  event classes are mutually exclusive (N-terminal > frame > antisense > intergenic).
- Card thresholds are recomputed from Tier1 each round; KEEP counts are not
  comparable across rounds.
- Deliverable tables and figures are English; lab records are Chinese.
- External tool paths come from `ORFSCAPE_*` environment variables or
  `orfscape/paths.local.json` — never hard-code a drive letter.

## Layout

`orfscape/` package (cli, orchestrator, progress/runlog/estimator, checks, core/,
gui/), `tests/`, `scripts/`, `docs/`, `integrations/`, `GUI.md`, `README.md`.
Full operating contract: `integrations/AGENT_RUNBOOK.md`.
