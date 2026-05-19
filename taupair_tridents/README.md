# Tau Pair Supplemental Code

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

Run the comparison plot from this directory:

```bash
python differential_cross_section_check.py
```

The script computes the three cross sections, prints the integrated totals, and
opens a plot comparing the `f2py` results with the CSV reference curves in
`data/reference/`.

