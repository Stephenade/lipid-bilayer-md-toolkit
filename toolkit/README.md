# Toolkit

General-purpose analyses for all-atom CHARMM36 binary lipid bilayers simulated
in GROMACS. Every notebook has a config cell at the top: set the topology
(`.tpr`), the structure or trajectory, and the two lipid residue names, then
run. With a single `.gro` frame the results are a snapshot; with a trajectory
they are averaged over frames.

| File | Property |
|---|---|
| `ka_from_box_fluctuations.py` | Area compressibility modulus K<sub>A</sub> = k<sub>B</sub>T⟨A⟩/Var(A) from Box-X and Box-Y, with block-averaged standard error and area per lipid. Command-line tool. |
| `lipid_analysis.ipynb` | All of the analyses below in one notebook; the trajectory is loaded once. |
| `01_interaction_energies.ipynb` | Lipid–lipid short-range Coulomb and Lennard-Jones energies from the `.edr` (needs `energygrps` set in the run). |
| `02_order_parameter.ipynb` | Deuterium order parameter −S<sub>CD</sub> per carbon for the sn-1 (`C3x`) and sn-2 (`C2x`) chains. |
| `03_spatial_distribution.ipynb` | Lateral density maps per leaflet. |
| `04_voronoi.ipynb` | Periodic 2D Voronoi tessellation of phosphates: area per lipid, neighbour count, like-neighbour fraction. |
| `05_water_depth.ipynb` | Water penetration into the hydrophobic core. |
| `06_tail_splay.ipynb` | Splay angle between the two acyl chains. |
| `07_rdf.ipynb` | 2D lateral P–P radial distribution functions. |
| `08_diffusion.ipynb` | Lateral diffusion coefficient from the centre-of-mass MSD (needs a `-pbc nojump` trajectory). |

## Preparing trajectories

```bash
gmx trjconv -s md.tpr -f md.xtc        -o md_whole.xtc  -pbc whole
gmx trjconv -s md.tpr -f md_whole.xtc  -o md_center.xtc -pbc mol -center   # structural analyses
gmx trjconv -s md.tpr -f md_whole.xtc  -o md_nojump.xtc -pbc nojump        # diffusion
gmx energy  -f md.edr -o box.xvg                                            # Box-X, Box-Y for K_A
```

## Conventions and caveats

- CHARMM36 atom names: sn-1 chain `C3x`, sn-2 chain `C2x` (carbonyls `C21`/`C31` excluded), phosphate `P`, TIP3 water oxygen `OH2`.
- Interaction energies are the short-range Coul-SR + LJ-SR terms; the PME reciprocal part cannot be split by group, so keep cut-offs identical across systems you compare.
- Diffusion coefficients depend on box size (Yeh–Hummer correction).
- Check the Voronoi area per lipid against the box estimate the notebook prints.
