#!/usr/bin/env python3
"""Exact MI-S4 reconstruction. Python 3 standard library; no stored answers imported."""
from itertools import permutations, combinations, product
from collections import Counter
from pathlib import Path
from math import comb
import json

OUT = Path(__file__).resolve().parent
ID = tuple(range(5))
def compose(p, q): return tuple(p[q[i]] for i in range(len(p)))
def inverse(p): return tuple(p.index(i) for i in range(len(p)))
def parity(p): return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))
G = sorted(p for p in permutations(range(5)) if parity(p) == 1)
def perm(*cycles):
    p = list(ID)
    for c in cycles:
        c = tuple(i-1 for i in c)
        for a, b in zip(c, c[1:]+c[:1]): p[a] = b
    return tuple(p)
def closure(gens):
    seen = {ID}; todo = [ID]
    for p in todo:
        for q in gens:
            r = compose(p, q)
            if r not in seen: seen.add(r); todo.append(r)
    return frozenset(seen)
def power(g, k):
    r = ID
    for _ in range(k): r = compose(r, g)
    return r
def order(g): return next(k for k in (1,2,3,5) if power(g,k) == ID)
def cycles(p):
    left = set(range(len(p))); ans = []
    while left:
        i = min(left); c = []
        while i in left: left.remove(i); c.append(i); i = p[i]
        ans.append(c)
    return ans
def cosets(H):
    left = set(G); sets = []
    while left:
        g = min(left); c = frozenset(compose(g,h) for h in H)
        sets.append(c); left -= c
    index = {g:i for i,c in enumerate(sets) for g in c}
    acts = {g:tuple(index[compose(g,min(c))] for c in sets) for g in G}
    assert all(compose(acts[g],acts[h]) == acts[compose(g,h)] for g in G for h in G)
    return sets, acts
subgroups = {
    5:closure((perm((1,2,3)),perm((1,2,4)))),
    6:closure((perm((1,2,3,4,5)),perm((2,5),(3,4)))),
    12:closure((perm((1,2,3,4,5)),)),
    20:closure((perm((1,2,3)),))}
actions = {}; coset_data = {}
for d,H in subgroups.items():
    cs, actions[d] = cosets(H); assert len(cs) == d
    coset_data[d] = {'stabilizer':sorted(H), 'representatives':[min(c) for c in cs]}
left = set(G); classes = []
while left:
    g = min(left); c = {compose(compose(x,g),inverse(x)) for x in G}
    classes.append(c); left -= c
classes.sort(key=lambda c:(order(min(c)),min(c)))
reps = [min(c) for c in classes]
assert list(map(len,classes)) == [1,15,20,12,12]
classof = {g:i for i,c in enumerate(classes) for g in c}
def fixed(p): return sum(i == x for i,x in enumerate(p))
def cycle_coeff(lengths, k):
    a = [1]+[0]*k
    for r in lengths:
        for j in range(k,r-1,-1): a[j] += (-1)**(r+1)*a[j-r]
    return a[k]
def newton(trace, g, k):
    a = [1]
    for j in range(1,k+1):
        n = sum((-1)**(r-1)*trace(power(g,r))*a[j-r] for r in range(1,j+1))
        assert n % j == 0; a.append(n//j)
    return a[k]

# Construct every orbit and an actual sign change making its wedge basis permuted.
triples = set(combinations(range(12),3)); orbits = []; gauge = {}
while triples:
    base = min(triples); transported = {}
    for g in G:
        image = tuple(actions[12][g][i] for i in base)
        t = tuple(sorted(image)); s = parity(image)
        if t in transported: assert transported[t] == s
        transported[t] = s
    stab = [g for g in G if tuple(sorted(actions[12][g][i] for i in base)) == base]
    signs = [parity(tuple(actions[12][g][i] for i in base)) for g in stab]
    assert all(s == 1 for s in signs) and len(stab) in (1,3)
    assert len(transported)*len(stab) == 60
    gauge.update(transported)
    orbits.append({'representative':base, 'size':len(transported),
                   'stabilizer':stab, 'stabilizer_signs':signs,
                   'oriented_basis':[{'triple':t,'sign':s} for t,s in sorted(transported.items())]})
    triples -= set(transported)
assert sorted(o['size'] for o in orbits) == [20,20,60,60,60]
for g in G:
    for t,s in gauge.items():
        image = tuple(actions[12][g][i] for i in t)
        assert s*parity(image) == gauge[tuple(sorted(image))]
sub3 = {closure((g,)) for g in G if order(g) == 3}
assert len(sub3) == 10
assert all(sum(all(tuple(sorted(actions[12][g][i] for i in t)) == t for g in H)
                   for t in gauge) == 4 for H in sub3)

# Full element-by-element characters, two independent exterior-power calculations.
qtrace = lambda g: fixed(actions[12][g])
triple_trace = lambda g: sum(parity(tuple(actions[12][g][i] for i in t)) for t in gauge
                            if tuple(sorted(actions[12][g][i] for i in t)) == t)
utrace = lambda g: 1+triple_trace(g)
vtrace = lambda g: fixed(actions[5][g])+2*fixed(actions[6][g])
alttrace = lambda g: 3*fixed(actions[6][g])-1  # 2*1+3B, B=P-1
records = []
for g in G:
    qcycles = [len(c) for c in cycles(actions[12][g])]
    vcycles = [len(c) for d in (5,6,6) for c in cycles(actions[d][g])]
    tc = triple_trace(g)
    assert tc == cycle_coeff(qcycles,3) == newton(qtrace,g,3)
    assert tc == 2*fixed(actions[20][g])+3*(60 if g == ID else 0)
    lhs = newton(utrace,g,2); rhs = cycle_coeff(vcycles,8)
    assert lhs == rhs == newton(vtrace,g,8) == newton(alttrace,g,8)
    records.append({'g':g,'class':classof[g], 'Q_cycles':qcycles,'V_cycles':vcycles,
                    'triple_trace':tc, 'U_trace':utrace(g), 'V_trace':vtrace(g),
                    'alternative_V_trace':alttrace(g), 'common_exterior_trace':lhs})
common = [newton(utrace,g,2) for g in reps]
assert common == [24310,-110,10,0,0]
assert [utrace(g) for g in reps] == [221,1,5,1,1]
assert [vtrace(g) for g in reps] == [17,5,2,2,2]
assert [alttrace(g) for g in reps] == [17,5,-1,2,2]
assert sum(len(c)*utrace(min(c)) for c in classes)//60 == 6
assert sum(len(c)*vtrace(min(c)) for c in classes)//60 == 3
assert comb(221,2) == comb(17,8) == 24310
def article_character(v):
    a,b,c,d = v
    return [a+4*b+5*c+6*d, a+c-2*d, a+b-c, a-b+d, a-b+d]
assert article_character((5,16,17,11)) == [triple_trace(g) for g in reps]
assert article_character((6,16,17,11)) == [utrace(g) for g in reps]
assert article_character((3,1,2,0)) == [vtrace(g) for g in reps]
assert article_character((2,0,3,0)) == [alttrace(g) for g in reps]
old_u_character = article_character((3,13,20,11))
old_u_trace = lambda g: old_u_character[classof[g]]
assert all(newton(old_u_trace,g,2) == newton(alttrace,g,8) for g in G)

# Independent Euclidean geometry in Z[phi], phi^2=phi+1.
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def sub(x,y): return add(x,neg(y))
def mul(x,y): return (x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def dot(x,y):
    z = (0,0)
    for a,b in zip(x,y): z = add(z,mul(a,b))
    return z
def det(a,b,c):
    return add(sub(mul(a[0],sub(mul(b[1],c[2]),mul(b[2],c[1]))),
                   mul(a[1],sub(mul(b[0],c[2]),mul(b[2],c[0])))),
               mul(a[2],sub(mul(b[0],c[1]),mul(b[1],c[0]))))
vertices = sorted({v for s,t in product((-1,1),repeat=2) for v in
                   (((0,0),(s,0),(0,t)), ((s,0),(0,t),(0,0)), ((0,t),(0,0),(s,0)))})
def distance2(i,j):
    diff = tuple(sub(a,b) for a,b in zip(vertices[i],vertices[j])); return dot(diff,diff)
types = {(4,0):'short', (4,4):'long', (8,4):'diameter'}
assert len(vertices) == 12
pairtype = {ij:types[distance2(*ij)] for ij in combinations(range(12),2)}
for kind in ('short','long'):
    for i in range(12):
        neighbors = [j for j in range(12) if j != i and pairtype[tuple(sorted((i,j)))] == kind]
        assert len(neighbors) == 5
        assert all(sum(k != j and pairtype[tuple(sorted((j,k)))] == kind for k in neighbors) == 2
                   for j in neighbors)
shape = lambda t: tuple(sorted(pairtype[tuple(sorted(ij))] for ij in combinations(t,2)))
shape_counts = Counter(shape(t) for t in combinations(range(12),3))
assert sorted(shape_counts.values()) == [20,20,60,60,60]

# Build the actual 60 rotations by transporting an ordered triangular face.
source = next(t for t in combinations(range(12),3) if shape(t) == ('short',)*3)
src = [vertices[i] for i in source]; srcdet = det(*src); assert srcdet != (0,0)
target_gram = [[dot(a,b) for b in src] for a in src]
src_signatures = [tuple(dot(v,b) for b in src) for v in vertices]
rotations = set()
for dest in permutations(range(12),3):
    bs = [vertices[i] for i in dest]
    if [[dot(a,b) for b in bs] for a in bs] != target_gram or det(*bs) != srcdet: continue
    signatures = {tuple(dot(v,b) for b in bs):i for i,v in enumerate(vertices)}
    assert len(signatures) == 12 and all(k in signatures for k in src_signatures)
    r = tuple(signatures[k] for k in src_signatures)
    assert all(distance2(i,j) == distance2(r[i],r[j]) for i,j in combinations(range(12),2))
    rotations.add(r)
assert len(rotations) == 60
assert all(compose(g,h) in rotations for g in rotations for h in rotations)
remaining = set(combinations(range(12),3)); geo_orbits = []
while remaining:
    t = min(remaining); orb = {tuple(sorted(g[i] for i in t)) for g in rotations}
    assert len({shape(s) for s in orb}) == 1
    geo_orbits.append({'representative':t,'size':len(orb),'side_types':shape(t),'triples':sorted(orb)})
    remaining -= orb
assert len(geo_orbits) == 5
assert Counter({tuple(o['side_types']):o['size'] for o in geo_orbits}) == shape_counts

result = {'conventions':{'permutations':'zero-based; left composition',
                         'classes':'1,2A,3A,5A,5B', 'quadratic_coordinates':'[a,b] means a+b*phi'},
          'coset_models':coset_data, 'class_representatives':reps,
          'class_sizes':list(map(len,classes)), 'triple_orbits':orbits,
          'all_60_element_checks':records, 'common_character':common,
          'article_coordinates':{'U':[6,16,17,11],'V':[3,1,2,0],'alternative_V':[2,0,3,0]},
          'old_selected_pair':{'U_article':[3,13,20,11], 'V_article':[2,0,3,0],
                               'U_character':old_u_character,'verified_on_all_60_elements':True},
          'geometric_vertices':vertices, 'geometric_rotations':sorted(rotations),
          'geometric_triple_orbits':geo_orbits, 'invariant_dimensions':{'U':6,'V':3}}
(OUT/'MI_S4_exact_results.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: triple orbit sizes 20,20,60,60,60; all stabilizer signs positive.')
print('PASS: explicit orientation transport checked on all 60 x 220 actions.')
print('PASS: independent exact icosahedral rotations and all five geometric shapes.')
print('PASS: exterior identity and alternative input on all 60 elements.')
print('PASS: old selected 24310 pair remains valid, with the same exterior character.')
print('Common exterior character:', common)
for o in sorted(geo_orbits,key=lambda o:(o['size'],o['side_types'])):
    print('Geometric orbit:',o['size'],', '.join(o['side_types']))
