# Seven Coincidences and Sixty Symmetries

**Ernesto Lupercio and Evgeniya Shaydurova**

Manuscripts, illustrations and exact calculations for seven exterior-power
identities of the icosahedral rotation group A5. The constructions use simplex
edges, vertices, axes, faces, vertex triples and stabilizer character lines.
The proof companion gives the finite enumeration, rigidity statements,
spinorial and arithmetic calculations, Klein residue construction, and modular
and Conway interpretations.

## Contents

| Material | Location |
|---|---|
| Main article, PDF and LaTeX | `Seven_Coincidences_and_Sixty_Symmetries_Intelligencer.{pdf,tex}` |
| Proof companion, PDF and LaTeX | `Proof_Companion.{pdf,tex}` |
| Four figures and their sources | `figures/`, `figures/source/` |
| Exact programs, inputs and certificates | `verification/` |
| Package and citation metadata | `MANIFEST.json`, `CITATION.cff` |

## Reproduction

Use Python 3.10 or later and the packages pinned in `requirements.txt`. A TeX
installation must provide `pdflatex`, `newpxtext`, `newpxmath`, `standalone`,
`tikz` and the other packages loaded by the two manuscript sources. No network
access is used by the verification programs or the manuscript build.

```bash
python -m pip install -r requirements.txt
python run_all.py --all-math
```

This runs all sixteen mathematical programs, rebuilds the four figure PDFs,
builds both manuscripts with three LaTeX passes each, and checks agreement of
the manuscript tables, certificates and references. The mathematical programs
are unchanged by the editorial revision. Both manuscripts have been rebuilt
and their tables checked against the certificates.

```bash
python run_all.py                 # twelve finite/reconstruction programs and both PDFs
python run_all.py --math-only     # the twelve programs without PDF builds
python run_all.py --all-math --math-only
python run_all.py --build-only    # figure and manuscript builds with consistency check
python figures/build_figures.py  # figures alone
```

## Calculations

All paths in this table are relative to `verification/`.

| Program | Calculation |
|---|---|
| `verify_revision.py` | Complete complex realization lists, exclusions at 3003, harmonics and Mathieu subgroup |
| `src/verify_finite_core.py` | Finite characters, signed orbits and semisimple obstruction |
| `src/verify_geometry.py` | Klein residue and geometric calculations |
| `constructions/verify_faces_1540.py` | Face construction and permutation uniqueness at 1540 |
| `constructions/verify_vertex_triples.py` | Triple orbits, transported orientations and 24310 |
| `constructions/verify_character_lines.py` | Cubic character lines, 7140 and its specified exterior-square root |
| `induction/verify_wave_i.py` | All 59 subgroups, 159 linear characters and restricted realization counts |
| `verify_integration.py` | Manuscript rows and certificate agreement |
| `reconstruction/reconstruct.py` | Independent subgroup, induction, inverse-character and explicit-action reconstruction |
| `reconstruction/mathieu.py` | Independent right-action word paths, orbit pairing and Golay preservation |
| `reconstruction/klein.py` | Independent weight enumeration and positive harmonic reductions |
| `reconstruction/compare.py` | Comparison of completed reconstruction, earlier results and manuscript tables |
| `src/verify_spin_arithmetic.py` | Spin restrictions and arithmetic identities |
| `src/verify_marked_basis.py` | Marked basis and Klein currents |
| `src/verify_modular_conformal.py` | Exact modular data and finite character-series checks |
| `src/verify_native_spectrum.py` | Exact sparse-spinor and Krylov spectral calculation |

The last four programs require `--all-math`. The four reconstruction programs
use the Python standard library; only their final comparison reads earlier
certificates. The independent search returns the same 102 nontrivial ordered
complex module pairs and all five classes of realization counts. The article's
rows at 1540, 24310 and 7140 are checked against explicit actions; other valid
alternatives remain in the companion and complete solution lists.

Input tensors and matrices are in `verification/prior_inputs/`. Generated JSON
certificates are in `verification/certificates/`, `constructions/`, `induction/`
and `reconstruction/`, with the complete original enumeration in
`verification/revision_results.json`. Its `duration_seconds` field records
performance only and is excluded from mathematical result comparisons.

The programs prove the stated finite computations. Classical arithmetic,
spinorial, modular and conformal results imported from the literature retain
their citations; finite series checks do not replace those theorems.

## Figures

Figure 1 adds a vector diagram of the normalized half-triangle above the
original Pascal-addition artwork. The original artwork is supplied as PNG
images. Revised labels are typeset as vector text over that artwork;
Figure 2 is unchanged raster art.
The figure source README gives the build instructions and source locations.

## Citation

Use the article title, authors and package version in `CITATION.cff`. No public
repository URL, release DOI or license has been assigned in this package.
