#!/usr/bin/env python3
"""Cross-check the article and proof companion after a repository build."""
from pathlib import Path
import re, sys
ROOT=Path(__file__).resolve().parents[1]
A=(ROOT/"Seven_Coincidences_and_Sixty_Symmetries_Intelligencer.tex").read_text()
C=(ROOT/"Proof_Companion.tex").read_text()

required_article=[
r"10,\quad2,\quad11,\quad19,\quad3,\quad26,\quad31",
r"4,\quad2,\quad5,\quad3,\quad1,\quad6,\quad11",
r"46H_0\oplus103H_1\oplus96H_2\oplus103H_3\oplus96H_4\oplus53H_5",
r"M_2(\Gamma(5))\cong B\oplus C",
r"\mathfrak a|_{A_5}\cong Q\oplus Q",
]
missing=[x for x in required_article if x not in A]
if missing: raise SystemExit("Missing article invariants: "+repr(missing))

labels=set(re.findall(r"\\label\{([^}]+)\}",C))
needed={"sec:finite","subsec:A1","subsec:A4","sec:orbits","prop:C1","subsec:D1","subsec:D2","subsec:D4","subsec:E2","subsec:F1","subsec:F2","subsec:F3","subsec:F6","subsec:F7","subsec:G1","subsec:H1"}
if not needed<=labels: raise SystemExit("Missing companion labels: "+repr(sorted(needed-labels)))
assert {f'subsec:A{i}' for i in range(5,9)} <= labels
from verify_integration import check
check()

for doc,text in [("article",A),("companion",C)]:
    bib=set(re.findall(r"\\bibitem\{([^}]+)\}",text))
    cites=set()
    for m in re.finditer(r"\\cite(?:p|t|alp|author|year)?\*?(?:\[[^]]*\])*\{([^}]+)\}",text):
        cites.update(x.strip() for x in m.group(1).split(","))
    if bib!=cites:
        raise SystemExit(f"{doc} bibliography mismatch: uncited={sorted(bib-cites)}, missing={sorted(cites-bib)}")
print("Article/companion consistency check: PASS")
