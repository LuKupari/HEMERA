"""Example: plot d sigma / dE1 for tau-pair trident production."""

import numpy as np
import matplotlib.pyplot as plt

import cross_sections as xs


emu = 1000.0
e1_min = 0.0
bin_width = 50.0
nbins = 20
ncall = 1_000_000
itmx = 5

mx = 4.0
xi = 5.0

a_sm = xs.a_sm_cross_section(
    emu=emu,
    e1_min=e1_min,
    bin_width=bin_width,
    nbins=nbins,
    ncall=ncall,
    itmx=itmx,
)
b_sm = xs.b_sm_cross_section(
    emu=emu,
    e1_min=e1_min,
    bin_width=bin_width,
    nbins=nbins,
    ncall=ncall,
    itmx=itmx,
)
b_bsm = xs.b_bsm_cross_section(
    emu=emu,
    mx=mx,
    xi=xi,
    e1_min=e1_min,
    bin_width=bin_width,
    nbins=nbins,
    ncall=ncall,
    itmx=itmx,
)

fig, ax = plt.subplots()

for label, result in [
    ("SM a diagram", a_sm),
    ("SM b diagram", b_sm),
    (rf"BSM b diagram, $m_S={mx:g}$ GeV, $\xi={xi:g}$", b_bsm),
]:
    bin_edges = np.append(result["bin_lo"], result["bin_hi"][-1])
    differential_sigma = result["sigma"] / bin_width
    ax.step(
        bin_edges,
        np.append(differential_sigma, differential_sigma[-1]),
        where="post",
        label=label,
    )

total = (a_sm["sigma"] + b_sm["sigma"] + b_bsm["sigma"]) / bin_width
bin_edges = np.append(a_sm["bin_lo"], a_sm["bin_hi"][-1])
ax.step(
    bin_edges,
    np.append(total, total[-1]),
    where="post",
    color="black",
    linewidth=2.0,
    label="Total",
)

ax.set_yscale("log")
ax.set_xlabel(r"Outgoing muon energy $E_1$ [GeV]")
ax.set_ylabel(r"$d\sigma / dE_1$ [cm$^2$ GeV$^{-1}$]")
ax.set_ylim(1e-44,1e-34)
ax.set_title(rf"Si, $E_\mu = {emu:.0f}$ GeV")
ax.grid(True, which="both", alpha=0.3)
ax.legend()

plt.show()
