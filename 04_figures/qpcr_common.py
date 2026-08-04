"""Shared loading / relative-expression maths for the qPCR figures.

All values derive from ``results/qpcr/supplementary_qpcr_data.csv``, the same file the
log2FC tables in the manuscript are built from.  Keeping one implementation here
guarantees the figures and the tables cannot drift apart.
"""

import csv
import os
from statistics import mean

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
DATA = os.path.join(ROOT, "results", "qpcr", "supplementary_qpcr_data.csv")

MATURE = "Mature neuron"
MODEL = "H2O2-induced AD model"
VORINO = "H2O2-induced AD model + Vorinostat"

MARKERS = ["APP", "PSEN1", "PSEN2"]
HUBS = ["HDAC1", "CREBBP", "NFKB1", "KAT2B", "SP1", "UBB",
        "ATP5F1A", "PRKACA", "PRKACB", "PSMD14"]


def load():
    """gene -> condition -> list of per-well records (in replicate order)."""
    out = {}
    with open(DATA, newline="") as fh:
        for row in csv.DictReader(fh):
            rec = {
                "replicate": int(row["well_replicate"]),
                "delta_Ct": float(row["delta_Ct"]) if row["delta_Ct"] else None,
                "flag": row["Cp_flag"],
            }
            out.setdefault(row["gene"], {}).setdefault(row["condition"], []).append(rec)
    for conds in out.values():
        for wells in conds.values():
            wells.sort(key=lambda r: r["replicate"])
    return out


def mean_dct(wells):
    vals = [w["delta_Ct"] for w in wells if w["delta_Ct"] is not None]
    return mean(vals) if vals else None


def log2fc(data, gene, test, calibrator):
    """log2FC = -ddCt, positive meaning higher expression in `test`."""
    t = mean_dct(data[gene][test])
    c = mean_dct(data[gene][calibrator])
    if t is None or c is None:
        return None
    return -(t - c)


def relative(data, gene, condition, calibrator):
    """2^-ddCt for one condition against the calibrator (calibrator = 1.0)."""
    v = log2fc(data, gene, condition, calibrator)
    return None if v is None else 2.0 ** v


def paired_log2fc(data, gene, test, calibrator):
    """Per-replicate log2FC, pairing each test well with the calibrator well
    that shares its replicate index.

    Both conditions of a replicate sit in the same well group of the plate, so
    pairing them cancels the systematic offset between the two well groups.
    The mean of these paired values is algebraically identical to
    ``log2fc``, which is what the manuscript tables report.

    Returns a list of ``(value, censored)`` tuples.
    """
    cal = {w["replicate"]: w for w in data[gene][calibrator]}
    out = []
    for w in data[gene][test]:
        c = cal.get(w["replicate"])
        if c is None or w["delta_Ct"] is None or c["delta_Ct"] is None:
            continue
        censored = w["flag"] != "detected" or c["flag"] != "detected"
        out.append((-(w["delta_Ct"] - c["delta_Ct"]), censored))
    return out
