# Bundled tools

`orfscape/tools/` holds the helpers that are **not installable as a package** and
therefore ship with the source. Everything else is detected from your system.

| Path | What | Where it comes from |
|---|---|---|
| `rawdump/` | `rawdump.exe` + Thermo CommonCore DLLs — reads a Thermo `.raw` (header, scan table) without the Thermo runtime | locally built from [AI4MSOmics/eval_ms_pipeline_repo](https://github.com/AI4MSOmics/eval_ms_pipeline_repo) (Apache-2.0); upstream ships no prebuilt binary |
| `rawmethod/` | `rawmethod.exe` + OpenMcdf — extracts the instrument method embedded in a `.raw` | same upstream (Apache-2.0) |
| `pFind_modification.ini`, `restrict_pFind.cfg` | modification table and search configuration for pFind | derived from pFind 3.1; edit if your pFind uses a different table |

See `../NOTICE` for licenses. The Thermo DLLs are Thermo Fisher Scientific
property, redistributed only as the runtime `rawdump` requires.

## Not here: Prodigal

`prodigal.exe` is **not** bundled — it is GPL-3.0 and every user can install it:

```bat
conda install -c bioconda prodigal        :: or grab prodigal.exe from the release page
```

`orfscape doctor` reports it as `Prodigal` and `L1` degradation is explicit: without
it the ab initio layer is skipped and the run log says so, rather than producing an
empty layer silently.

## Detected, not bundled

pFind 3.1, pParse, msconvert, BLAST+, isoPSV, MaxQuant, Comet, MSFragger, Sage,
MetaMorpheus, AlphaPeptDeep, Augustus, GeneMarkS-2, rnaSPAdes, ThermoRawFileParser.

Point ORFscape at them in either way:

```bat
:: 1) environment variables (highest priority)
set ORFSCAPE_PFIND=D:\Bioinformatics\pFind3.1
set ORFSCAPE_BLAST=D:\Bioinformatics\BLAST\ncbi-blast-2.15.0+\bin
set ORFSCAPE_ISOPSV=D:\isoPSV_multiengine
set ORFSCAPE_TOOL_ROOTS=D:\Bioinformatics;C:\Program Files\ProteoWizard
```

```json
// 2) orfscape/paths.local.json  (git-ignored; copy paths.local.example.json)
{
  "pfind": "D:/Bioinformatics/pFind3.1",
  "blast": "D:/Bioinformatics/BLAST/ncbi-blast-2.15.0+/bin",
  "isopsv": "D:/isoPSV_multiengine"
}
```

Then check the whole environment:

```bat
python -m orfscape doctor --workdir orfscape_projects/envcheck
```
