# Lipid bilayer MD analysis toolkit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23851716.svg)](https://doi.org/10.5281/zenodo.23851716)

Analysis code for all-atom molecular dynamics (MD) simulations of mixed lipid
bilayers, used to design mechanically tuned liposomal drug carriers. It
computes the area compressibility modulus (K<sub>A</sub>), area per lipid,
acyl-chain order parameters, lipid–lipid and ion–lipid interaction energies,
lateral organisation (Voronoi, spatial distribution, RDF, lipid spacing),
water penetration, tail splay and lateral diffusion. It also contains the
figure scripts for the associated publications.

Developed in the group of Prof. Francisco R. Hung, Department of Chemical
Engineering, Northeastern University.

## Repository layout

| Folder | What it contains |
|---|---|
| [`toolkit/`](toolkit/) | General-purpose analyses for any binary bilayer (CHARMM36, GROMACS). One notebook per property, an all-in-one notebook, and a command-line K<sub>A</sub> tool. |
| [`papers/sn1-sn2-acyl-chains/`](papers/sn1-sn2-acyl-chains/) | Figure scripts for *Tuning the Mechanical Properties of DOPC-Based Lipid Bilayers through Multiscale Simulations and Targeted Modification of the sn-1 and sn-2 Acyl Chains* (27 PC/PA and PC/PG systems), with the published K<sub>A</sub> dataset. |
| [`papers/dopc-composition-series/`](papers/dopc-composition-series/) | Analysis of DOPC mixed with DOPA, DOPG, DLPA or DLPG at 65–95 mol% DOPC. |

## Installation

```bash
git clone https://github.com/Stephenade/lipid-bilayer-md-toolkit.git
cd lipid-bilayer-md-toolkit
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Python 3.9 or newer. Trajectory analyses need GROMACS output files (`.tpr`,
`.xtc`, `.edr`); these are not included because of their size.

## Quick start

**K<sub>A</sub> from box-area fluctuations** (K<sub>A</sub> = k<sub>B</sub>T⟨A⟩/Var(A), with a block-averaged standard error):

```bash
gmx energy -f md.edr -o box.xvg            # select Box-X and Box-Y
python toolkit/ka_from_box_fluctuations.py box.xvg --n-leaflet 144
```

**Structural properties:** open `toolkit/lipid_analysis.ipynb`, set the
topology, trajectory and lipid names in the config cell, and run the sections
you need. Each property is also available as a standalone notebook
(`toolkit/01_…` to `08_…`).

**Paper figures:**

```bash
cd papers/sn1-sn2-acyl-chains/scripts
python plot_ka_groups.py
```

## Data and citation

| Record | DOI |
|---|---|
| K<sub>A</sub> and APL dataset for the sn-1/sn-2 study (Zenodo) | [10.5281/zenodo.23850654](https://doi.org/10.5281/zenodo.23850654) |
| This software (Zenodo, all versions) | [10.5281/zenodo.23851716](https://doi.org/10.5281/zenodo.23851716) |

If you use this code, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub shows a *Cite this repository* button),
together with the relevant dataset and article.

## Licence

MIT. See [`LICENSE`](LICENSE).
