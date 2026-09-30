#!/usr/bin/env python3
"""Independently enumerate both 3003 torus-weight bases and their reductions."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
q=[0]+list(range(-5,6));p=[0]+list(range(-2,3))
symq=[q[i]+q[j] for i,j in combinations_with_replacement(range(12),2)]
extp=[p[i]+p[j] for i,j in combinations(range(6),2)]
left=Counter(sum(symq[i] for i in pair) for pair in combinations(range(78),2))
right=Counter(sum(extp[i] for i in subset) for subset in combinations(range(15),5))
assert sum(left.values())==sum(right.values())==3003
def harmonics(weights):
    assert all(weights[w]==weights[-w] for w in weights)
    result={j:weights[j]-weights[j+1] for j in range(max(weights)+1)}
    assert min(result.values())>=0
    return result
def reduce_terms(harmonic):
    counts=Counter(harmonic)
    for j in range(max(counts),5,-1):
        multiplicity=counts.pop(j,0);r=j%5
        counts[j-6]+=multiplicity;counts[r]+=multiplicity;counts[5-r]+=multiplicity
    return [counts[i] for i in range(6)]
hl,hr=harmonics(left),harmonics(right)
rl,rr=reduce_terms(hl),reduce_terms(hr)
assert rl==rr
def roots(c):
    if len(c)==2:return (c[0]-c[1],0)
    if len(c)==3:
        assert c[1]==c[2];return (c[0]-c[1],0)
    assert c[1]==c[4] and c[2]==c[3]
    return (c[0]-c[1],c[1]-c[2])
def hchar(j):
    answer=[(2*j+1,0)]
    for n,multiplier in [(2,1),(3,1),(5,1),(5,2)]:
        c=[0]*n
        for w in range(-j,j+1):c[w*multiplier%n]+=1
        answer.append(roots(c))
    return answer
def multiply(x,y):return (x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
independent=json.loads((ROOT/'independent_results.json').read_text())
irreps=independent['irreducible_characters'];sizes=independent['class_sizes']
matrix=[]
for ir in irreps:
    row=[]
    for j in range(6):
        terms=[multiply(a,b) for a,b in zip(ir,hchar(j))]
        value=tuple(sum(n*x[k] for n,x in zip(sizes,terms)) for k in range(2))
        assert value[1]==0 and value[0]%60==0
        row.append(value[0]//60)
    matrix.append(row)
mat=[[Fraction(a) for a in row] for row in matrix];pivots=[];row=0
for col in range(6):
    pivot=next((r for r in range(row,5) if mat[r][col]),None)
    if pivot is None:continue
    mat[row],mat[pivot]=mat[pivot],mat[row];d=mat[row][col]
    mat[row]=[a/d for a in mat[row]]
    for r in range(5):
        if r!=row:
            d=mat[r][col];mat[r]=[a-d*b for a,b in zip(mat[r],mat[row])]
    pivots.append(col);row+=1
free=[i for i in range(6) if i not in pivots]
assert len(free)==1
kernel=[Fraction(0)]*6;kernel[free[0]]=Fraction(-1)
for r,col in enumerate(pivots):kernel[col]=mat[r][free[0]]
assert all(x.denominator==1 for x in kernel)
kernel=list(map(int,kernel))
common=[sum(a*b for a,b in zip(row,rl)) for row in matrix]
hom_dimension=sum(a*a for a in common)
fs=[]
for ir in irreps:
    vals=[ir[c[2]] for c in independent['power_classes']]
    result=tuple(sum(n*x[k] for n,x in zip(sizes,vals)) for k in range(2))
    assert result==(60,0);fs.append(1)
record={'left_weight_multiplicities':dict(sorted(left.items())),
    'right_weight_multiplicities':dict(sorted(right.items())),
    'left_harmonics':hl,'right_harmonics':hr,'left_reduced':rl,'right_reduced':rr,
    'restriction_matrix':matrix,'restriction_rank':len(pivots),'primitive_integral_kernel':kernel,
    'common_complex_irrep_multiplicities':common,'equivariant_hom_dimension':hom_dimension,
    'Frobenius_Schur_indicators':fs,
    'choices_required':['invariant complements','Clebsch-Gordan decompositions','ordering repeated copies',
                        'identifications of multiplicity spaces']}
(ROOT/'klein_results.json').write_text(json.dumps(record,indent=2)+'\n')
print('3003 basis weights enumerated on each side; reductions:',rl,rr)
print('Restriction rank:',len(pivots),'; primitive kernel:',kernel)
print('Common irreducible multiplicities:',common,'; Hom dimension:',hom_dimension)
print('Frobenius-Schur indicators:',fs)
