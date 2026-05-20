# Tau Pair Tridents Cross-Sections

This code builds the Fortran integrations with `f2py` and plots the outgoing
muon energy differential cross sections against the reference curves.

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

Python scripts are included to generate two seperate differential cross sections,
as an example for using the modules; they can be run with.
```bash
python invariant_mass_distribution.py
python outgoing_muon_energy_distribution.py
```
These distribution files are also saved in the `data/` directory.


