import numpy as np
from modules import a_diagram_ext
from modules import b_diagram_sm_ext
from modules import b_diagram_bsm_ext

def result_dict(e1_lo, e1_hi, sigma, avgi, sd):
    return {
        "e1_lo": e1_lo,
        "e1_hi": e1_hi,
        "sigma": sigma,
        "avgi": avgi,
        "sd": sd,
        "covered_range": (float(e1_lo[0]), float(e1_hi[-1])),
    }
    
def energy_edges(e1_min: float, bin_width: float, nbins: int):
    if nbins <= 0:
        raise ValueError("nbins must be positive")
    if bin_width <= 0.0:
        raise ValueError("bin_width must be positive")
    e1_lo = e1_min + np.arange(nbins, dtype=np.float64) * bin_width
    return e1_lo, e1_lo + bin_width

def total_cross_section(result: dict) -> tuple[float, float]:
    sigma = result["sigma"]
    sd = result["sd"]
    return float(sigma.sum()), float(np.sqrt((sd**2).sum()))

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
) -> dict:
    """Compute the SM a-diagram cross section in outgoing-energy bins."""

    e1_lo, e1_hi = energy_edges(e1_min, bin_width, nbins)
    sigma, avgi, sd = a_diagram_ext.compute_a_sm(
        z, a, mp, mmu, me, emu, e1_lo, e1_hi, ncall, itmx,
        1 if verbose else 0,
    )
    return result_dict(e1_lo, e1_hi, sigma, avgi, sd)


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
) -> dict:
    """Compute the SM b-diagram cross section in outgoing-energy bins."""

    e1_lo, e1_hi = energy_edges(e1_min, bin_width, nbins)
    sigma, avgi, sd = b_diagram_sm_ext.compute_b_sm(
        z, a, mp, mmu, me, emu, e1_lo, e1_hi, ncall, itmx,
        1 if verbose else 0,
    )
    return result_dict(e1_lo, e1_hi, sigma, avgi, sd)


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
) -> dict:
    """Compute the BSM b-diagram cross section in outgoing-energy bins."""

    e1_lo, e1_hi = energy_edges(e1_min, bin_width, nbins)
    sigma, avgi, sd = b_diagram_bsm_ext.compute_b_bsm(
        z, a, mp, mmu, me, emu, mx, xi, e1_lo, e1_hi, ncall, itmx,
        1 if verbose else 0,
    )
    return result_dict(e1_lo, e1_hi, sigma, avgi, sd)
