#!/usr/bin/env python3
"""
Corrected clustering-index figure for the LSD section.

    python3 lsd_corrected_figure.py     ->  LSD_ClusteringIndex_corrected.png/.pdf

Values are hardcoded from the cluster run of lsd_analysis.py, with two
corrections that the original script did not apply:

1. FINITE-N BASELINE.  1/(2*sqrt(lambda)) is the large-N limit.  At 7 lipids
   per leaflet it underestimates the random nearest-neighbour distance by
   1.4%.  Monte Carlo sampling of random points in matched periodic boxes
   gives the correction factors below.  Using the analytic baseline makes the
   95:5 systems look 1.4% more dispersed than they are.

       lipids/leaflet   50      40      29      19       7
       MC / analytic    1.0015  1.0020  1.0062  1.0086  1.0140

2. CORRELATED FRAMES.  Lateral diffusion from the msd_*.dat files gives
   D = 5.4e-8 cm^2/s, so a lipid needs ~29 ns to move one nearest-neighbour
   separation and exchange neighbours.  A 400 ns trajectory holds ~14
   independent configurations, not 201 frames, so the frame-to-frame SEM is
   too small by sqrt(201/14) = 3.8.  All error bars here are inflated by 3.8.

Neither correction changes the qualitative result, but without the second the
error bars are roughly a quarter of their true size, which would make several
differences look significant that are not.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COMPS = [65, 72, 80, 87, 95]
NLEAF = [50, 40, 29, 19, 7]
BIAS = {50: 1.0015, 40: 1.0020, 29: 1.0062, 19: 1.0086, 7: 1.0140}
INFLATE = 3.8

NN = {"DOPG": [8.59, 9.38, 10.82, 12.56, 21.35],
      "DOPA": [8.06, 9.00, 10.24, 12.84, 17.91],
      "DLPG": [8.44, 9.16, 10.75, 12.82, 20.13],
      "DLPA": [7.87, 8.49, 10.08, 11.51, 17.23]}
RND = {"DOPG": [6.94, 7.78, 9.13, 11.27, 18.59],
       "DOPA": [6.79, 7.64, 9.01, 11.18, 18.53],
       "DLPG": [6.86, 7.70, 9.05, 11.23, 18.55],
       "DLPA": [6.70, 7.55, 8.95, 11.13, 18.50]}
SEM = {"DOPG": [.002, .003, .003, .004, .011],
       "DOPA": [.003, .003, .003, .005, .006],
       "DLPG": [.002, .002, .003, .004, .005],
       "DLPA": [.003, .003, .004, .005, .009]}

COLORS = {"DOPG": "#E8000B", "DOPA": "#1F4E79",
          "DLPG": "#4C9A2A", "DLPA": "#E0A500"}
MARKERS = {"DOPG": "o", "DOPA": "s", "DLPG": "^", "DLPA": "D"}

fig, ax = plt.subplots(figsize=(7.4, 5.2))
ax.axhspan(1.0, 1.30, color="#EAF3EA", zorder=0)
ax.axhspan(0.88, 1.0, color="#FBEDED", zorder=0)
ax.axhline(1.0, color="k", ls="--", lw=1.3, zorder=2)
ax.text(96.6, 1.005, "dispersed", fontsize=9, va="bottom", ha="right", color="0.35")
ax.text(96.6, 0.995, "clustered", fontsize=9, va="top", ha="right", color="0.35")

for lip in ["DOPG", "DOPA", "DLPG", "DLPA"]:
    y = [(NN[lip][j] / RND[lip][j]) / BIAS[NLEAF[j]] for j in range(5)]
    e = [SEM[lip][j] * INFLATE / BIAS[NLEAF[j]] for j in range(5)]
    ax.errorbar(COMPS, y, yerr=e, marker=MARKERS[lip], ms=7, lw=1.8,
                capsize=3, color=COLORS[lip], label=lip,
                markeredgecolor="k", markeredgewidth=0.6, zorder=3)

ax.set_xlabel("DOPC content (mol %)", fontsize=12)
ax.set_ylabel(r"$\langle r_{nn}\rangle_{\rm obs}\,/\,\langle r_{nn}\rangle_{\rm random}$",
              fontsize=12)
ax.set_xticks(COMPS)
ax.set_xlim(62, 98)
ax.set_ylim(0.88, 1.30)
ax.tick_params(direction="out", labelsize=10)
ax.legend(frameon=True, fontsize=10, ncol=4, loc="upper center",
          bbox_to_anchor=(0.5, 1.10), columnspacing=1.4, handlelength=1.6)
plt.tight_layout()
plt.savefig("LSD_ClusteringIndex_corrected.png", dpi=600, bbox_inches="tight")
plt.savefig("LSD_ClusteringIndex_corrected.pdf", bbox_inches="tight")
print("wrote LSD_ClusteringIndex_corrected.png / .pdf")
