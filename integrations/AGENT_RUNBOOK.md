# ORFscape agent runbook

This is the shell-agnostic operating contract for an agent that runs ORFscape.
Adapters in `zcode/`, `claude-code/` and `generic/` all reduce to this document.

## Role

You run ORFscape, a proteogenomics pipeline. Your contract is **numbers, not
narration**: every stage ends with a self-check, and you stop at the first sign of a
silent failure (empty database, missing layer, zero PSM, unlocated event).

## First moves, always

1. Read this runbook and `../orfscape/tools/README.md` (what is bundled, what must be
   installed).
2. Run the preflight:

   ```bash
   python integrations/scripts/preflight.py --workdir <proj> --raw <raw folder>
   ```

   Any FAIL → report which dependency is missing and **stop**. ORFscape's failure
   mode is "it ran and produced nothing", which costs hours; do not "try anyway".
3. Read `<proj>/run_state.json` and any `00_logs/manifest_*.json` to see what already
   ran. Do not repeat completed stages.

## Driving the pipeline

Prefer the driver — it runs one stage, checks it, and stops on failure:

```bash
python integrations/scripts/run_pipeline.py --workdir <proj> --raw <raw> \
    [--from P0] [--to P8c] [--skip P6b,QC] [--thread 32] [--workers 4] [--dry-run]
```

It delegates to `orfscape.orchestrator`, so the run also writes
`00_logs/progress.jsonl` and a human can attach the GUI (`python -m orfscape gui
--workdir <proj>`) without disturbing it.

If you run a stage by hand, use `--dry-run` first when available, then execute, then
**always**:

```bash
python integrations/scripts/stage_check.py --workdir <proj> --stage <S>
```

A stage is not done until its check passes.

## Stage table

| Stage | Command | Gate |
|---|---|---|
| P0 | `python -m orfscape scan --raw <raw> --workdir <proj>` | raw > 0, scans > 0 |
| P1 | `python -m orfscape genome --workdir <proj>` (or drop fna/gff/faa into `01_genome/<sp>/`) | three files present, non-empty |
| D0 | `python -m orfscape split --workdir <proj>` | MGF count == raw count |
| C | `python -m orfscape dbbuild --workdir <proj>` | L0/L2 non-empty; `L1_coords.tsv` if L1 exists |
| E1 | `python -m orfscape e1 --workdir <proj> --thread T` | tasks == MGFs, PSM > 0 |
| QC | `python -m orfscape qc --workdir <proj>` (+ `qcfigs`, `e1figs`) | figures and tables exist |
| E2 | `python -m orfscape e2 --workdir <proj> --dbs aug,6f --thread T` | tasks == MGFs × DBs, PSM > 0 |
| P6 | `python -m orfscape novel --workdir <proj>` | rows > 0; layer hits match non-empty layers |
| P6b | `python -m orfscape research --workdir <proj> --thread T --workers W` | merged ≥ P6; small DB |
| P7 | `python -m orfscape validate --workdir <proj> --from-research --min-cos auto` | KEEP > 0; every MISSING has a reason |
| P8a | `python -m orfscape novgene --workdir <proj> --from 06b_validate_research` | events > 0, none unlocated |
| P8b | `python -m orfscape novannot --workdir <proj> --from 06b_validate_research [--no-esm]` | attribute/homolog/feature tables exist |
| P8c | `python -m orfscape novreport --workdir <proj> --from 06b_validate_research` | report has non-zero statistics |
| P9 | `python -m orfscape novelfigs --workdir <proj> --from 06b_validate_research` | funnel + parameter figures |

Eukaryotic layer build (before any search):

```bash
python -m orfscape eukbuild --workdir <proj> --sp CREI --fna genome.fna --gff genomic.gff
python -m orfscape.checks --workdir <proj> --stage EUK
```

## Conventions you must not break

1. The novel-peptide layer whitelist covers `L1_`, `L2_`, `L3|`, `L6|`. A layer whose
   FASTA is non-empty but contributes zero novel peptides is a **bug**, not a result.
2. pFind `Positions` is a **0-based** amino-acid offset; aa→nt is always `pos*3`.
   A peptide drawn one codon off means someone used `(pos-1)*3`.
3. Enzyme calls match `TL` / `LL` exactly — never `"L" in "TL"`.
4. Tier1 = same core across enzymes (K/R trimmed); Tier2 = same locus, two enzymes;
   Tier3 = single enzyme.
5. Card thresholds are recomputed from Tier1 **every round**. Never compare KEEP
   counts between rounds without saying so explicitly.
6. Event classification: compare same-strand CDS only; among several overlapping CDS
   take the in-frame one (`rel % 3 == 0`). Include `05b_research` in the locus scan.
7. ORF counting: same contig + strand + frame and >80% overlap = one locus.
8. Events are mutually exclusive: N-terminal correction > different frame >
   antisense > intergenic.
9. Report **unique peptides** externally; when attributing GSSPs to ORFs, the four
   event classes must sum to the total.
10. All output paths are normalised with `realpath`; never write through `..`.
11. isoPSV has a single-instance lock: clear it before a rerun, and keep exactly one
    `.0.dta` per `File_Name` or it silently produces nothing.
12. Eukaryotic codon tables come from each transcript's `transl_table` (1/2/3), never
    a global constant.

## What to report after each stage

A short block: stage, command, wall time, and the concrete numbers the checker
printed (e.g. `E1: tasks 36/36, PSM 59,369`).

- WARN → quote the line and say what it implies; it is usually the early signal of a
  silently empty layer.
- FAIL → stop. Give symptom → likely cause → fix. Do not continue "to see what
  happens".

## Stop and escalation rules

- Stop on: preflight FAIL, non-zero exit, failed gate, missing layer, zero PSM.
- Do not re-run a completed stage to "make sure". Use `--force` only when the user
  asks, or when an upstream input changed (genome, DB, MGF, layer whitelist).
- Before any destructive action (deleting a stage directory, overwriting
  `novel_validated.tsv`), confirm with the user; move the old run under
  `_staging/<name>_<date>` instead of deleting it.
- If two runs disagree on a number, report both and say which configuration produced
  each — never silently pick one.

## Deliverables

Point the user at:

- `06b_validate_research/novel_validated.tsv` + `validate_summary.md` (cards, DROP/MISSING reasons)
- `07_novgene/NG_table.tsv` (events)
- `08_annotation/` (structure/function/homology) and `08_report/*.docx`
- `P13_results/` if it exists (statistics, QC figures, event tables, peptide booklets)

Summarise with the numbers the checker produced and list anything still WARN with a
one-line next step.
