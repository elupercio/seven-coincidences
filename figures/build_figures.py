#!/usr/bin/env python3
"""Rebuild the four figure PDFs from the supplied raster art and TeX labels."""
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
SOURCES = {
    'fig01_pascal': 'fig01_pascal_labels',
    'fig02_vertices_axes': 'fig02_vertices_axes',
    'fig03_spin_intersection': 'fig03_spin_intersection_labels',
    'fig04_divisor_cusps': 'fig04_divisor_cusps_labels',
}

def main():
    latex = shutil.which('pdflatex')
    if latex is None:
        raise RuntimeError('pdflatex not found')
    with tempfile.TemporaryDirectory(prefix='seven-figures-') as directory:
        for target, source in SOURCES.items():
            subprocess.run([latex, '-interaction=nonstopmode', '-halt-on-error',
                            f'-output-directory={directory}', f'{source}.tex'],
                           cwd=ROOT/'source', check=True)
            shutil.copy2(Path(directory)/f'{source}.pdf', ROOT/f'{target}.pdf')
    print('Rebuilt all four figure PDFs.')

if __name__ == '__main__':
    main()
