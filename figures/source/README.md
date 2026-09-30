# Figure sources

Run `python figures/build_figures.py` from the repository root to rebuild the
four PDF assets. The manuscript build runs this command automatically.

Figure 1(a) is a TikZ diagram of the normalized half-triangle, with binomial
entries computed from its row and position indices. Figure 1(b) retains the
original Pascal-addition artwork. The original PNGs for Figures 1, 3 and 4 are
in this directory; their `*_labels.tex` files overlay the revised labels as
vector text. Figure 2 uses `../fig02_vertices_axes.png` without alteration,
wrapped by `fig02_vertices_axes.tex`. No editable vector source for the
original artwork is included.

The PNGs in the parent directory are viewing copies. The article includes the
rebuilt PDFs. Rebuilding requires `pdflatex`, `standalone`, `graphicx`, `tikz`,
`newpxtext` and `newpxmath` and their dependencies.
