# DOPC composition series with anionic lipids

Analysis of all-atom CHARMM36 bilayers of DOPC mixed with DOPA, DOPG, DLPA or
DLPG at 65, 72, 80, 87 and 95 mol% DOPC.

| Notebook | Contents |
|---|---|
| `composition_series_analysis.ipynb` | K<sub>A</sub> from box fluctuations, then RDF, tail splay, water depth, order parameters, diffusion, Voronoi and spatial distribution for each of the four series, and the K<sub>A</sub> bar chart |
| `rdf_all_systems.ipynb` | Lateral RDFs for all systems and compositions |
| `rdf_DOPC_DOPA_composition.ipynb` | RDF for the DOPC/DOPA series |
| `tail_splay_DOPC_DOPA_composition.ipynb` | Tail splay for the DOPC/DOPA series |
| `interaction_energies.ipynb` | Decomposed partial and total interaction energies, ion–lipid energies, K<sub>A</sub> against interaction energy |
| `interaction_energy_figures.ipynb` | Publication figures for the interaction energies |

## Expected input files

Most cells expect, per composition (65, 72, 80, 87, 95):
`Prod_<c>.tpr`, `Prod_merged_<c>.edr`, `Prod_merged_center_<c>.xtc` and
`Prod_merged_nojump_<c>.xtc`, and `box_<lipid>_<c>.xvg` for K<sub>A</sub>.
Edit the file names in each cell's config block to match your runs.

## Notes

- The K<sub>A</sub> bar chart uses corrected DOPA values written into the cell (explained in its header), not the original spreadsheet.
- The "load your simulation results" cells in `interaction_energy_figures.ipynb` import a helper module, `load_energies.py`, which is not yet in this repository. The figure cells run on the values written into the notebook.
