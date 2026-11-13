"""Per-peptide AlphaPeptDeep theoretical similarity table.

One row per novel peptide: the measured spectrum used for the mirror plot, the
AlphaPeptDeep cosine, how many predicted ions matched inside 20 ppm, plus the
Tier / db / event / keep context so the column can be sorted and filtered.

    python sjh_simtable.py --table1 A.tsv --spectra B C --cache D.pkl \
        --predict E.txt F.txt [G.txt] --tsv OUT.tsv --xlsx OUT.xlsx

The peak cache is read through a whitelist unpickler (plain containers and
numbers only); nothing else is de-serialised.
"""
import argparse
import pickle
import sys

import numpy as np
import pandas as pd

import sjh_pdfs as S

_SAFE_TYPES = ("list", "dict", "tuple", "set", "frozenset", "float", "int",
               "str", "bytes", "bool")


class SafeUnpickler(pickle.Unpickler):
    """Only plain containers / scalars may be reconstructed."""

    def find_class(self, module, name):
        if module == "builtins" and name in _SAFE_TYPES:
            return super().find_class(module, name)
        raise pickle.UnpicklingError("blocked global: %s.%s" % (module, name))


def load_cache(path):
    with open(path, "rb") as fh:
        return SafeUnpickler(fh).load()


def quant(v, q):
    return float(np.percentile(v, q)) if len(v) else float("nan")


def build_rows(table1, spectra, cache_path, pred_paths):
    t1 = S.load_table(table1)
    cache = load_cache(cache_path)
    wanted = {S.iln(r["peptide"]) for r in t1}
    psms = S.best_psms(spectra, wanted)

    pred = {}
    for p in pred_paths:
        pred.update(S.read_predictions(p))

    rows = []
    for r in t1:
        psm = psms.get(S.iln(r["peptide"]))
        if psm is None:
            continue
        stem, scan = S.split_title(psm.get("File_Name", ""))
        peaks = cache.get((stem, scan)) or []
        mz = np.array([x[0] for x in peaks], float) if peaks else np.array([])
        it = np.array([x[1] for x in peaks], float) if peaks else np.array([])
        if len(mz):
            o = np.argsort(mz)
            mz, it = mz[o], it[o]
        try:
            chg = int(float(psm.get("Charge") or 0))
        except ValueError:
            chg = 0
        sticks = S.pred_sticks(r["peptide"], r["modification"],
                               pred.get((r["peptide"], chg)))
        cos = S.theory_cosine(mz, it, sticks)
        n_hit = 0
        if sticks and len(mz):
            for (m, _v, _i) in sticks:
                lo, hi = m * (1 - 20e-6), m * (1 + 20e-6)
                if np.searchsorted(mz, hi, side="right") > \
                        np.searchsorted(mz, lo, side="left"):
                    n_hit += 1
        tail = stem.split("_")
        tag = tail[-2].upper() if len(tail) >= 2 else ""
        enzyme = ("Trypsin" if tag.startswith("T")
                  else ("LysargiNase" if tag.startswith("L") else "-"))
        rows.append({
            "peptide": r["peptide"],
            "modification": r["modification"],
            "charge": chg,
            "enzyme": enzyme,
            "psm_score": psm.get("Final_Score", ""),
            "n_pred_ions": len(sticks),
            "n_matched_ions": n_hit,
            "sim_alphapeptdeep": np.nan if cos is None else round(cos, 4),
            "tier": r.get("tier", ""),
            "db_class": r.get("db_class", ""),
            "event_type": r.get("event_type", ""),
            "event_sub": r.get("event_sub", ""),
            "final_keep": r.get("final_keep", ""),
            "keep_reason": r.get("keep_reason", ""),
            "orf": r.get("orf", ""),
            "pep_nt_start": r.get("pep_nt_start", ""),
            "pep_nt_end": r.get("pep_nt_end", ""),
            "spectrum": psm.get("File_Name", "").split("/")[-1],
        })
    df = pd.DataFrame(rows)
    return df.sort_values("sim_alphapeptdeep", ascending=False, na_position="last")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--table1", required=True)
    ap.add_argument("--spectra", nargs="+", required=True)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--predict", nargs="+", required=True)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--xlsx", default=None)
    args = ap.parse_args(argv)

    df = build_rows(args.table1, args.spectra, args.cache, args.predict)
    df.to_csv(args.tsv, sep="\t", index=False, encoding="utf-8-sig")
    if args.xlsx:
        df.to_excel(args.xlsx, index=False, sheet_name="similarity")

    v = df["sim_alphapeptdeep"].dropna().to_numpy(float)
    print("rows %d (with similarity %d)" % (len(df), len(v)))
    print("median %.3f  p25 %.3f  p75 %.3f  min %.3f  max %.3f"
          % (quant(v, 50), quant(v, 25), quant(v, 75), v.min(), v.max()))
    for grp in ("tier", "final_keep", "db_class"):
        print("--", grp)
        for k, g in df.groupby(grp):
            gv = g["sim_alphapeptdeep"].dropna().to_numpy(float)
            print("   %-26s n=%-4d median %.3f" % (k, len(gv), quant(gv, 50)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
