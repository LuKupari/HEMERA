"""
Python frontend for the a/b diagram Fortran integrations.

Build all extensions once, or after editing the Fortran:
    python3 build_modules.py

Then import and use:
    import cross_sections as xs
    res = xs.a_sm_cross_section(ncall=10000, itmx=2)
"""

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent
FORTRAN_DIR = ROOT_DIR / "fortran"
MODULES_DIR = ROOT_DIR / "modules"

FORTRAN_SOURCE = {
    "a_sm": {
        "module": "a_diagram_ext",
        "source": FORTRAN_DIR / "a_sm" / "a_diagram_f2py.f",
        "workdir": FORTRAN_DIR / "a_sm",
        "routine": "compute_a_sm",
    },
    "b_sm": {
        "module": "b_diagram_sm_ext",
        "source": FORTRAN_DIR / "b_sm" / "b_diagram_sm_f2py.f",
        "workdir": FORTRAN_DIR / "b_sm",
        "routine": "compute_b_sm",
    },
    "b_bsm": {
        "module": "b_diagram_bsm_ext",
        "source": FORTRAN_DIR / "b_bsm" / "b_diagram_bsm_f2py.f",
        "workdir": FORTRAN_DIR / "b_bsm",
        "routine": "compute_b_bsm",
    },
}


def clean_build_dir(cfg: dict, verbose: bool) -> None:
    """Remove stale f2py outputs from the Fortran source directory."""
    for so in cfg["workdir"].glob(f"{cfg['module']}*.so"):
        if verbose:
            print(f"  Removing stale {so.name}")
        so.unlink()


def build(which: str = "all", verbose: bool = True) -> None:
    MODULES_DIR.mkdir(exist_ok=True)

    if which == "all":
        diagrams_to_build = list(FORTRAN_SOURCE)
    else:
        diagrams_to_build = [which]
    
    for diagram in diagrams_to_build:
        if diagram not in FORTRAN_SOURCE:
            raise ValueError(f"unkown extension {diagram!r}")
        
        files = FORTRAN_SOURCE[diagram]

        clean_build_dir(files, verbose=verbose)
        
        cmd = [
            sys.executable,
            "-m",
            "numpy.f2py",
            "-c",
            files["source"].name,
            "-m",
            files["module"],
            "only:",
            files["routine"],
            ":",
            "--f77flags=-std=legacy",
        ]

        if verbose:
            print(f"Building {diagram}:  " + " ".join(cmd) )
            
            
        result = subprocess.run(cmd,
                                cwd=files["workdir"],
                                capture_output=True,
                                text=True,
        )
        
        importlib.invalidate_caches()
        sys.modules.pop(files["module"], None)
        sys.modules.pop(f"modules.{files['module']}", None)
        
        if result.returncode != 0:
            output = ""
            if result.stdout:
                output += result.stdout
            if result.stderr:
                output += result.stderr
            raise RuntimeError(f"f2py build failed for {files}.\n{output}")

        built_modules = list(files["workdir"].glob(f"{files['module']}*.so"))
        if not built_modules:
            raise RuntimeError(
                f"Build finished, but no {files['module']}*.so file was found"
            )

        for built_module in built_modules:
            for stale_module in MODULES_DIR.glob(f"{files['module']}*.so"):
                stale_module.unlink()
            destination = MODULES_DIR / built_module.name
            shutil.move(str(built_module), destination)

if __name__ == "__main__":
    build(which="all")
