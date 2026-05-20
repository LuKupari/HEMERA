import numpy as np
from modules import a_diagram_ext
from modules import b_diagram_sm_ext
from modules import b_diagram_bsm_ext


def a_sm_cross_section(
    emu: float = 1000.0,
    z: float = 14.0,
    a: float = 28.0,
    mp: float = 0.94,
    mmu: float = 0.106,
    me: float = 1.78,
    e1_min: float = 0.0,
    bin_width: float = 50.0,
    nbins: int = 20,
    ncall: int = 1_000_000,
    itmx: int = 10,
    verbose: bool = False,
    invariant_mass: bool = False,
) -> dict:
    """Compute the SM a-diagram cross section in selected distribution bins."""

    if nbins <= 0:
        raise ValueError("nbins must be positive")
    if bin_width <= 0.0:
        raise ValueError("bin_width must be positive")

    bin_lo = e1_min + np.arange(nbins, dtype=np.float64) * bin_width
    bin_hi = bin_lo + bin_width
    distribution = "invariant_mass" if invariant_mass else "outgoing_muon_energy"
    sigma, avgi, sd = a_diagram_ext.compute_a_sm(
        z, a, mp, mmu, me, emu, bin_lo, bin_hi, ncall, itmx,
        1 if verbose else 0,
        int(invariant_mass),
    )
    return {
        "bin_lo": bin_lo,
        "bin_hi": bin_hi,
        "e1_lo": bin_lo,
        "e1_hi": bin_hi,
        "distribution": distribution,
        "sigma": sigma,
        "avgi": avgi,
        "sd": sd,
        "covered_range": (float(bin_lo[0]), float(bin_hi[-1])),
    }


def b_sm_cross_section(
    emu: float = 1000.0,
    z: float = 14.0,
    a: float = 28.0,
    mp: float = 0.94,
    mmu: float = 0.106,
    me: float = 1.78,
    e1_min: float = 0.0,
    bin_width: float = 50.0,
    nbins: int = 20,
    ncall: int = 1_000_000,
    itmx: int = 10,
    verbose: bool = False,
    invariant_mass: bool = False,
) -> dict:
    """Compute the SM b-diagram cross section in selected distribution bins."""

    if nbins <= 0:
        raise ValueError("nbins must be positive")
    if bin_width <= 0.0:
        raise ValueError("bin_width must be positive")

    bin_lo = e1_min + np.arange(nbins, dtype=np.float64) * bin_width
    bin_hi = bin_lo + bin_width
    distribution = "invariant_mass" if invariant_mass else "outgoing_muon_energy"
    sigma, avgi, sd = b_diagram_sm_ext.compute_b_sm(
        z, a, mp, mmu, me, emu, bin_lo, bin_hi, ncall, itmx,
        1 if verbose else 0,
        int(invariant_mass),
    )
    return {
        "bin_lo": bin_lo,
        "bin_hi": bin_hi,
        "e1_lo": bin_lo,
        "e1_hi": bin_hi,
        "distribution": distribution,
        "sigma": sigma,
        "avgi": avgi,
        "sd": sd,
        "covered_range": (float(bin_lo[0]), float(bin_hi[-1])),
    }


def b_bsm_cross_section(
    emu: float = 1000.0,
    mx: float = 4.0,
    xi: float = 5.0,
    z: float = 14.0,
    a: float = 28.0,
    mp: float = 0.94,
    mmu: float = 0.106,
    me: float = 1.78,
    e1_min: float = 0.0,
    bin_width: float = 50.0,
    nbins: int = 20,
    ncall: int = 1_000_000,
    itmx: int = 10,
    verbose: bool = False,
    invariant_mass: bool = False,
) -> dict:
    """Compute the BSM b-diagram cross section in selected distribution bins."""

    if nbins <= 0:
        raise ValueError("nbins must be positive")
    if bin_width <= 0.0:
        raise ValueError("bin_width must be positive")

    bin_lo = e1_min + np.arange(nbins, dtype=np.float64) * bin_width
    bin_hi = bin_lo + bin_width
    distribution = "invariant_mass" if invariant_mass else "outgoing_muon_energy"
    sigma, avgi, sd = b_diagram_bsm_ext.compute_b_bsm(
        z, a, mp, mmu, me, emu, mx, xi, bin_lo, bin_hi, ncall, itmx,
        1 if verbose else 0,
        int(invariant_mass),
    )
    return {
        "bin_lo": bin_lo,
        "bin_hi": bin_hi,
        "e1_lo": bin_lo,
        "e1_hi": bin_hi,
        "distribution": distribution,
        "sigma": sigma,
        "avgi": avgi,
        "sd": sd,
        "covered_range": (float(bin_lo[0]), float(bin_hi[-1])),
    }
