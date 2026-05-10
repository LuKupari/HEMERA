"""Plot outgoing-muon-energy cross sections against reference curves."""

import numpy as np
import matplotlib.pyplot as plt

import cross_sections as xs


PARAMS = {
    "emu": 1000.0,
    "e1_min": 0.0,
    "bin_width": 50.0,
    "nbins": 20,
    "ncall": 1000000,
    "itmx": 5,
    "verbose": False,
}

MX = 4.0
XI = 5.0


def main() -> None:
    a_sm = xs.a_sm_cross_section(**PARAMS)
    b_sm = xs.b_sm_cross_section(**PARAMS)
    b_bsm = xs.b_bsm_cross_section(mx=MX, xi=XI, **PARAMS)

    total, err = xs.total_cross_section(a_sm)
    print(f"SM (a): {total:.6e} +/- {err:.2e} cm^2")
    total, err = xs.total_cross_section(b_sm)
    print(f"SM (b): {total:.6e} +/- {err:.2e} cm^2")
    total, err = xs.total_cross_section(b_bsm)
    print(f"BSM (b), m_S={MX:g} GeV, xi={XI:g}: {total:.6e} +/- {err:.2e} cm^2")

    plots = [
        ("SM (a)", a_sm, "data/reference/a-sm.csv"),
        ("SM (b)", b_sm, "data/reference/b-sm.csv"),
        (rf"BSM (b), $m_S={MX:g}$ GeV, $\xi={XI:g}$", b_bsm, "data/reference/b-bsm.csv"),
    ]

    for title, result, reference_file in plots:
        reference = np.loadtxt(reference_file, delimiter=",")
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.semilogy(result["e1_lo"], result["sigma"], linewidth=1.8, label="f2py")
        ax.semilogy(reference[:, 0], reference[:, 1], color="red", linestyle=":", linewidth=1.8, label="reference")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r"Outgoing muon energy $E_1$ [GeV]")
        ax.set_ylabel(r"$d\sigma / dE_1$ [cm$^2$ GeV$^{-1}$]")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(frameon=False, loc="best")
        fig.suptitle(r"Si, $E_\mu = %.0f$ GeV" % PARAMS["emu"])

    plt.show()


if __name__ == "__main__":
    main()
