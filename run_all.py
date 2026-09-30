#!/usr/bin/env python3
"""Recompute finite checks and/or build the article and proof companion."""
from __future__ import annotations
import argparse, shutil, subprocess, sys, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
DOCS=["Proof_Companion","Seven_Coincidences_and_Sixty_Symmetries_Intelligencer"]

def run(args,cwd=ROOT):
    print("+"," ".join(map(str,args)),flush=True)
    subprocess.run(args,cwd=cwd,check=True)

def main():
    ap=argparse.ArgumentParser()
    g=ap.add_mutually_exclusive_group()
    g.add_argument("--build-only",action="store_true")
    g.add_argument("--math-only",action="store_true")
    ap.add_argument("--all-math",action="store_true",
                    help="also run the specialist spin/conformal/native-spectrum calculations")
    ns=ap.parse_args()
    if not ns.build_only:
        programs=["verification/verify_revision.py",
                  "verification/src/verify_finite_core.py",
                  "verification/src/verify_geometry.py",
                  "verification/constructions/verify_faces_1540.py",
                  "verification/constructions/verify_vertex_triples.py",
                  "verification/constructions/verify_character_lines.py",
                  "verification/induction/verify_wave_i.py",
                  "verification/verify_integration.py",
                  "verification/reconstruction/reconstruct.py",
                  "verification/reconstruction/mathieu.py",
                  "verification/reconstruction/klein.py",
                  "verification/reconstruction/compare.py"]
        if ns.all_math:
            programs += ["verification/src/verify_spin_arithmetic.py",
                         "verification/src/verify_marked_basis.py",
                         "verification/src/verify_modular_conformal.py",
                         "verification/src/verify_native_spectrum.py"]
        for x in programs: run([sys.executable,str(ROOT/x)])
    if not ns.math_only:
        run([sys.executable,str(ROOT/"figures/build_figures.py")])
        pdflatex=shutil.which("pdflatex")
        if not pdflatex: raise RuntimeError("pdflatex not found")
        with tempfile.TemporaryDirectory(prefix="seven-coincidences-") as td:
            build=Path(td)
            for stem in DOCS:
                for _ in range(3):
                    run([pdflatex,"-interaction=nonstopmode","-halt-on-error",
                         f"-output-directory={build}",f"{stem}.tex"])
                shutil.copy2(build/f"{stem}.pdf",ROOT/f"{stem}.pdf")
        run([sys.executable,str(ROOT/"verification/check_consistency.py")])
    print("Completed successfully.")
if __name__=="__main__":
    try: main()
    except Exception as e:
        print(e,file=sys.stderr); sys.exit(1)
