#!/usr/bin/env python3
"""
Area compressibility modulus (K_A) from box-area fluctuations.

    K_A = k_B * T * <A> / Var(A)

where A(t) = Box_X(t) * Box_Y(t) is the total projected area of the bilayer
(Feller & Pastor, J. Chem. Phys. 1999). The number of lipids does not enter
K_A; it is only used to report the area per lipid (APL).

The error is the standard error of K_A over consecutive blocks of the
trajectory (default 5 blocks).

Input: one or more GROMACS .xvg files with columns  time  Box-X  Box-Y,
for example from

    gmx energy -f md.edr -o box.xvg      # select Box-X and Box-Y

Usage:

    python ka_from_box_fluctuations.py box.xvg --n-leaflet 144
    python ka_from_box_fluctuations.py box_*.xvg --n-leaflet 144 --csv ka.csv
"""
import argparse
import csv
import sys

import numpy as np

K_B = 1.38065e-23  # Boltzmann constant (J/K)


def read_box_xvg(path):
    """Return time (ps), Box-X (nm) and Box-Y (nm) arrays from a GROMACS .xvg file."""
    time, box_x, box_y = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "@")) or not line.strip():
                continue
            vals = line.split()
            if len(vals) >= 3:
                time.append(float(vals[0]))
                box_x.append(float(vals[1]))
                box_y.append(float(vals[2]))
    if not time:
        raise ValueError(f"no data rows found in {path}")
    return np.array(time), np.array(box_x), np.array(box_y)


def ka_from_area(area_nm2, temperature=298.15):
    """K_A in mN/m from an area time series in nm^2."""
    mean_a = np.mean(area_nm2) * 1e-18          # m^2
    var_a = np.var(area_nm2, ddof=1) * 1e-36    # m^4
    return K_B * temperature * mean_a / var_a * 1000.0


def block_error(area_nm2, temperature=298.15, n_blocks=5):
    """Standard error of K_A over consecutive blocks, plus the block values."""
    size = len(area_nm2) // n_blocks
    blocks = []
    for i in range(n_blocks):
        block = area_nm2[i * size:(i + 1) * size]
        if len(block) > 1 and np.var(block, ddof=1) > 0:
            blocks.append(ka_from_area(block, temperature))
    blocks = np.array(blocks)
    if len(blocks) < 2:
        return float("nan"), blocks
    return np.std(blocks, ddof=1) / np.sqrt(len(blocks)), blocks


def analyse(path, n_leaflet, temperature=298.15, n_blocks=5):
    time, box_x, box_y = read_box_xvg(path)
    area = box_x * box_y                       # nm^2
    ka = ka_from_area(area, temperature)
    err, blocks = block_error(area, temperature, n_blocks)
    result = {
        "file": path,
        "frames": len(time),
        "time_ns": (time[-1] - time[0]) / 1000.0,
        "mean_area_nm2": float(np.mean(area)),
        "var_area_nm4": float(np.var(area, ddof=1)),
        "KA_mN_per_m": float(ka),
        "KA_sem_mN_per_m": float(err),
        "block_KA_mN_per_m": [float(b) for b in blocks],
    }
    if n_leaflet:
        result["APL_A2"] = float(np.mean(area) / n_leaflet * 100.0)
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("xvg", nargs="+", help=".xvg file(s) with time, Box-X, Box-Y columns")
    p.add_argument("--n-leaflet", type=int, default=None,
                   help="lipids per leaflet (only used for the area per lipid)")
    p.add_argument("--temperature", type=float, default=298.15, help="temperature in K (default 298.15)")
    p.add_argument("--blocks", type=int, default=5, help="number of blocks for the error (default 5)")
    p.add_argument("--csv", help="also write the results to this CSV file")
    args = p.parse_args(argv)

    rows = []
    for path in args.xvg:
        r = analyse(path, args.n_leaflet, args.temperature, args.blocks)
        rows.append(r)
        apl = f"  APL = {r['APL_A2']:.1f} A^2" if "APL_A2" in r else ""
        print(f"{path}: {r['frames']} frames, {r['time_ns']:.1f} ns  "
              f"K_A = {r['KA_mN_per_m']:.1f} +/- {r['KA_sem_mN_per_m']:.1f} mN/m{apl}")
        print("    blocks: " + ", ".join(f"{b:.1f}" for b in r["block_KA_mN_per_m"]))

    if args.csv:
        fields = ["file", "frames", "time_ns", "mean_area_nm2", "var_area_nm4",
                  "APL_A2", "KA_mN_per_m", "KA_sem_mN_per_m"]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
