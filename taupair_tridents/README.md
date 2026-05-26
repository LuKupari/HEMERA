# Tau Pair Tridents Cross-Sections

This code builds the Fortran integrations with `f2py` and plots the outgoing
muon energy and invariant mass differential cross sections. These cross-sections are calculated for
SM (Bethe-Heitler-like) and BSM (scalar mediated) channels.

## Setup

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate taupair
```

## Build

Build the compiled modules:

```bash
python build_modules.py
```

This writes the compiled extension modules into `modules/`.

## Run
Python scripts are included to generate two separate differential cross sections as examples of how to use the modules. 
They can be run with:
```bash
python invariant_mass_distribution.py
python outgoing_muon_energy_distribution.py
```
These distribution files are also saved in the `data/` directory.


