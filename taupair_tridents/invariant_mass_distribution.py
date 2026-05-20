"""Example: plot d sigma / dm_tau tau for tau-pair trident production."""

import numpy as np
import matplotlib.pyplot as plt

import cross_sections as xs


emu = 1000.0
mtautau_min = 1.5
bin_width = 1.0
nbins = 13
ncall = 1_000_000
itmx = 5

bsm_points = [
    (4.0, 5.0),
    (10.0, 20.0),
]

a_sm = xs.a_sm_cross_section(
    emu=emu,
    e1_min=mtautau_min,
    bin_width=bin_width,
    nbins=nbins,
    ncall=ncall,
    itmx=itmx,
    invariant_mass=True,
)
b_sm = xs.b_sm_cross_section(
    emu=emu,
    e1_min=mtautau_min,
    bin_width=bin_width,
    nbins=nbins,
    ncall=ncall,
    itmx=itmx,
    invariant_mass=True,
)
b_bsm_results = [
    (
        mx,
        xi,
        xs.b_bsm_cross_section(
            emu=emu,
            mx=mx,
            xi=xi,
            e1_min=mtautau_min,
            bin_width=bin_width,
            nbins=nbins,
            ncall=ncall,
            itmx=itmx,
            invariant_mass=True,
        ),
    )
    for mx, xi in bsm_points
]

fig, ax = plt.subplots()

for label, result in [
    ("SM a diagram", a_sm),
    ("SM b diagram", b_sm),
]:
    bin_edges = np.append(result["bin_lo"], result["bin_hi"][-1])
    differential_sigma = result["sigma"] / bin_width
    ax.step(
        bin_edges,
        np.append(differential_sigma, differential_sigma[-1]),
        where="post",
        label=label,
    )

for mx, xi, result in b_bsm_results:
    bin_edges = np.append(result["bin_lo"], result["bin_hi"][-1])
    differential_sigma = result["sigma"] / bin_width
    ax.step(
        bin_edges,
        np.append(differential_sigma, differential_sigma[-1]),
        where="post",
        label=rf"BSM b diagram, $m_S={mx:g}$ GeV, $\xi={xi:g}$",
    )

ax.set_yscale("log")
ax.set_ylim(1e-42, 1e-32)
ax.set_xlim(3,14.5)
ax.set_xlabel(r"Invariant mass $m_{\tau\tau}$ [GeV]")
ax.set_ylabel(r"$d\sigma / dm_{\tau\tau}$ [cm$^2$ GeV$^{-1}$]")
ax.set_title(rf"Si, $E_\mu = {emu:.0f}$ GeV")
ax.grid(True, which="both", alpha=0.3)
ax.legend()

plt.show()
