# Tuning the mechanical properties of DOPC-based bilayers via the sn-1 and sn-2 acyl chains

Figure scripts for *Tuning the Mechanical Properties of DOPC-Based Lipid
Bilayers through Multiscale Simulations and Targeted Modification of the sn-1
and sn-2 Acyl Chains* (S. Adepoju and F. R. Hung; under review).

The study covers 27 all-atom bilayers in three groups: PC/PA at 80:20 and
PC/PG at 65:35 and 95:5, with DOPC, POPC or SOPC as the PC lipid and the PA or
PG lipid carrying DO, PO or SO tails. A coarse-grained screen of 17 systems
was used to select the all-atom compositions.

## Data

`data/KA_PC_PA_PG_Mixtures.csv` is a copy of the published dataset:
Adepoju & Hung, Zenodo, [10.5281/zenodo.23850654](https://doi.org/10.5281/zenodo.23850654).
`plot_ka_groups.py` and `plot_ka_vs_interaction_energy.py` read K<sub>A</sub>,
area per lipid and errors from this file, so the figures always match the
published values. The other scripts carry their values in the script itself.

## Scripts

Run from the `scripts/` folder; each writes PNG and PDF files to the current
directory.

| Script | Figures |
|---|---|
| `plot_ka_groups.py` | K<sub>A</sub> per group and for all 27 systems, with pure DOPC as reference |
| `plot_ka_vs_interaction_energy.py` | K<sub>A</sub> against total lipid–lipid interaction energy, per group and combined; prints the regression statistics |
| `plot_order_parameters.py` | Order parameters for all systems and by chain type (18:1, 16:0, 18:0) |
| `plot_order_sn1.py`, `plot_order_sn2.py` | sn-1 and sn-2 order parameters for all 27 systems |
| `plot_partial_interaction_energies.py` | Partial interaction energies: secondary–secondary, PC–PC, PC–secondary, total |
| `plot_decomposed_interaction_energies.py` | Decomposed total interaction energies, 4-panel and stacked |
| `plot_ion_lipid_energies.py` | Na⁺ and Cl⁻ interactions with the PA or PG lipid |
| `plot_lipid_spacing_distribution.py` | Lipid spacing distribution: profiles, nearest-neighbour distances, clustering index |
| `plot_lsd_clustering_corrected.py` | Clustering index with the finite-N random baseline correction (supersedes the clustering-index panel above) |
| `plot_cg_ka.py` | K<sub>A</sub> of the 17 coarse-grained systems |
| `voronoi_diagrams.py` | Representative Voronoi diagrams per leaflet (needs `Prod_0.tpr` and `Prod_merged.xtc`) |
