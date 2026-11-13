---
name: orfscape-proteogenomics
description: End-to-end proteogenomics with ORFscape (prokaryotic P0–P9 plus eukaryotic layer building). Use when a genome plus DDA mass spectra must be turned into novel peptides / novel ORFs / N-terminal corrections: data confirmation, six-layer database (L0–L6), E1/E2 searches, P6 novel peptides, P6b candidate re-search, P7 per-spectrum arbitration, P8 event classification and annotation, P9 figures and report. Also use to resume, audit or debug an existing ORFscape project.
---

# ORFscape proteogenomics

This skill turns the validated ORFscape pipeline into explicit steps with a gate
after every stage. Full operating contract: [`../AGENT_RUNBOOK.md`](../AGENT_RUNBOOK.md).
Read it once; the table below is the short form.

## When to use / not use

**Use** — prokaryotic genome + DDA MS data, looking for novel peptides, novel ORFs,
N-terminal start corrections, antisense or different-frame expression; or resuming,
re-running and auditing an existing ORFscape project.

**Do not use** — quantitative proteomics that needs no custom database, or
eukaryotic *search* (only `eukbuild` layer construction is supported in v0.1).

## 0. Preflight

```bash
python integrations/scripts/preflight.py --workdir orfscape_projects/<proj> --raw <raw folder>
```

Required: Python ≥3.9, `orfscape` importable, pFind, isoPSV, writable project
directory, raw files. Any FAIL — stop. Every later stage would produce empty results
silently.

## 1. Directory contract

```
orfscape_projects/<proj>/
├── 00_input/          scan outputs (groups.csv / raw_info.tsv / species.tsv)
├── 01_genome/<sp>/    genomic.fna / genomic.gff / protein.faa
├── 02_db/<sp>/        L0_annot, L1_*, L2_6frame, L3_altstart, L4, L5, L6_tx_orfs,
│                      db_annot, db_aug, db_6f{,1..3}, L1_coords.tsv
├── 03_mgf/            one MGF per raw (file name = task name)
├── 04_search_annot/   E1 tasks      04_search_disc/  E2 tasks
├── 04a_qc/            QC figures
├── 05_novel/<sp>/     novel_peptides.tsv (P6)
├── 05b_research/<sp>/ novel_after_research.tsv (P6b)
├── 06b_validate_research/novel_validated.tsv (P7)
├── 07_novgene/        NG_table.tsv (events)
├── 08_annotation/     structure / function / homology + report
└── P13_results/       result pack (statistics, QC, events, booklets)
```

## 2. Stages

| Stage | Command | Gate |
|---|---|---|
| P0 | `python -m orfscape scan --raw <raw> --workdir <proj>` | raw > 0, scans > 0 |
| P1 | `python -m orfscape genome --workdir <proj>` (or place the triplet by hand) | fna/gff/faa non-empty |
| D0 | `python -m orfscape split --workdir <proj>` | MGF count == raw count |
| C | `python -m orfscape dbbuild --workdir <proj>` | L0/L2 non-empty; `L1_coords.tsv` if L1 |
| E1 | `python -m orfscape e1 --workdir <proj> --thread T` | tasks == MGFs, PSM > 0 |
| QC | `python -m orfscape qc --workdir <proj>` (+ `qcfigs`, `e1figs`) | figures and tables exist |
| E2 | `python -m orfscape e2 --workdir <proj> --dbs aug,6f --thread T` | tasks == MGFs × DBs, PSM > 0 |
| P6 | `python -m orfscape novel --workdir <proj>` | rows > 0; layer hits match non-empty layers |
| P6b | `python -m orfscape research --workdir <proj> --thread T --workers W` | merged ≥ P6; small DB |
| P7 | `python -m orfscape validate --workdir <proj> --from-research --min-cos auto` | KEEP > 0; MISSING has reasons |
| P8 | `python -m orfscape novgene|novannot|novreport --workdir <proj> --from 06b_validate_research` | events > 0, annotations and report present |
| P9 | `python -m orfscape novelfigs --workdir <proj> --from 06b_validate_research` | funnel + parameter figures |

Eukaryotic layers (before any search):

```bash
python -m orfscape eukbuild --workdir <proj> --sp CREI \
    --fna genome.fna --gff genomic.gff [--transcripts tx.fa] [--ref-protein protein.faa]
python -m orfscape.checks --workdir <proj> --stage EUK
```

## 3. Conventions (do not change casually)

1. Layer whitelist covers `L1_`, `L2_`, `L3|`, `L6|` — a non-empty layer with zero
   novel peptides is a bug.
2. pFind `Positions` is 0-based; aa→nt is always `pos*3`.
3. Enzyme calls match `TL`/`LL` exactly, never `"L" in "TL"`.
4. Tier1 = same core across enzymes; Tier2 = same locus, two enzymes; Tier3 = one enzyme.
5. Card thresholds are recomputed from Tier1 every round; never compare KEEP counts
   across rounds without saying so.
6. Event classification compares same-strand CDS only, and takes the in-frame CDS
   (`rel % 3 == 0`) when several overlap.
7. Same contig + strand + frame and >80% overlap = one ORF locus.
8. Events are mutually exclusive: N-terminal correction > different frame > antisense > intergenic.
9. Report unique peptides; event classes must sum to the total.
10. isoPSV has a single-instance lock; keep exactly one `.0.dta` per `File_Name`.

## 4. Stage check

```bash
python integrations/scripts/stage_check.py --workdir <proj>            # all stages
python integrations/scripts/stage_check.py --workdir <proj> --stage P6 # one stage
```

FAIL means stop. WARN goes into the lab notebook — most WARNs are the early signal
of a silently lost layer.

## 5. One-shot run (resumable)

```bash
python integrations/scripts/run_pipeline.py \
  --workdir orfscape_projects/<proj> --raw <raw folder> \
  [--from P0] [--to P8c] [--thread 32] [--workers 4] [--skip P6b] [--dry-run]
```

Each step: run → `stage_check` → stop on failure with the log tail. Logs and manifest
go to `<proj>/00_logs/`, together with `progress.jsonl`, which the GUI can attach to
(`python -m orfscape gui --workdir <proj>`).

## 6. Deliverables

- `06b_validate_research/novel_validated.tsv` + `validate_summary.md`
- `07_novgene/NG_table.tsv`
- `08_annotation/` + `08_report/*.docx`
- `P13_results/`
- Lab record and handover note (numbers, cards, conventions, open issues)
