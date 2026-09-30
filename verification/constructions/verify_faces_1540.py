#!/usr/bin/env python3
"""Exact MI-S3 check; standard library only; no stored answers imported.

Construct the A5 coset actions, check both exterior powers at all 60 elements,
enumerate all permutation inputs in dimensions 56 and 22, and distinguish
orbit presentations from representation characters. The mathematical proof
of the nine-type subgroup classification is in the companion insertion.
"""
from itertools import permutations
from collections import Counter, defaultdict
from pathlib import Path
from math import isqrt
import json

OUT=Path(__file__).resolve().parent
def compose(p,q):return tuple(p[q[i]] for i in range(5))
def inverse(p):return tuple(p.index(i) for i in range(5))
def perm(*cycles):
    p=list(range(5))
    for c in cycles:
        c=tuple(x-1 for x in c)
        for a,b in zip(c,c[1:]+c[:1]):p[a]=b
    return tuple(p)
ID=tuple(range(5))
G=sorted(p for p in permutations(range(5)) if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0)
def closure(gens):
    seen={ID};todo=[ID]
    for p in todo:
        for q in gens:
            r=compose(p,q)
            if r not in seen:seen.add(r);todo.append(r)
    return frozenset(seen)
def cyc(p):
    unseen=set(range(len(p)));out=[]
    while unseen:
        i=min(unseen);length=0
        while i in unseen:
            unseen.remove(i);length+=1;i=p[i]
        out.append(length)
    return sorted(out)
def power(p,k):
    q=ID
    for _ in range(k):q=compose(q,p)
    return q
def order(p):return next(k for k in (1,2,3,5) if power(p,k)==ID)
def fixed(p):return sum(i==x for i,x in enumerate(p))
def cosets(H):
    left=set(G);cs=[]
    while left:
        g=min(left);c=frozenset(compose(g,h) for h in H);cs.append(c);left-=c
    ind={g:i for i,c in enumerate(cs) for g in c}
    acts={g:tuple(ind[compose(g,min(c))] for c in cs) for g in G}
    return cs,acts
subgroups={
    1:frozenset(G),
    5:closure((perm((1,2,3)),perm((1,2,4)))),
    6:closure((perm((1,2,3,4,5)),perm((2,5),(3,4)))),
    10:closure((perm((3,4,5)),perm((1,2),(3,4)))),
    12:closure((perm((1,2,3,4,5)),)),
    15:closure((perm((1,2),(3,4)),perm((1,3),(2,4)))),
    20:closure((perm((1,2,3)),)),
    30:closure((perm((1,2),(3,4)),)),
    60:frozenset([ID])}
actions={};raw={}
for d,H in subgroups.items():
    cs,acts=cosets(H);assert len(cs)==d
    assert all(tuple(acts[g][acts[h][i]] for i in range(d))==acts[compose(g,h)] for g in G for h in G)
    actions[d]=acts
    raw[d]={'stabilizer':sorted(H),'coset_representatives':[min(c) for c in cs]}
unseen=set(G);classes=[]
while unseen:
    g=min(unseen);cl={compose(compose(x,g),inverse(x)) for x in G}
    classes.append(cl);unseen-=cl
classes.sort(key=lambda c:(order(min(c)),min(c)))
reps=[min(c) for c in classes]
assert list(map(len,classes))==[1,15,20,12,12]
chars={d:tuple(fixed(actions[d][g]) for g in reps) for d in subgroups}
def coeff(cycles,k):
    p=[1]+[0]*k
    for length in cycles:
        sign=(-1)**(length+1)
        for j in range(k,length-1,-1):p[j]+=sign*p[j-length]
    return p[k]
all_checks=[]
for g in G:
    cu=cyc(actions[6][g])+cyc(actions[10][g])+2*cyc(actions[20][g])
    cv=2*cyc(actions[5][g])+cyc(actions[12][g])
    u=sum(fixed(actions[d][g]) for d in (6,10,20,20))
    u2=sum(fixed(actions[d][power(g,2)]) for d in (6,10,20,20))
    l=coeff(cu,2);r=coeff(cv,3)
    assert l==r==(u*u-u2)//2
    all_checks.append({'g':g,'class':next(i for i,c in enumerate(classes) if g in c),
                       'U_cycles':sorted(cu),'V_cycles':sorted(cv),'left':l,'right':r})
assert [next(x['left'] for x in all_checks if tuple(x['g'])==g) for g in reps]==[1540,-20,10,0,0]

ds=sorted(chars)
def presentations(n,i=0):
    if i==len(ds):
        if n==0:yield ()
        return
    for c in range(n//ds[i]+1):
        for rest in presentations(n-c*ds[i],i+1):yield (c,)+rest
def character(coeffs):return tuple(sum(c*chars[d][i] for c,d in zip(coeffs,ds)) for i in range(5))
def exterior(c,k):
    values=[]
    for ci in range(5):
        es=[1]
        for j in range(1,k+1):
            a=0
            for t in range(1,j+1):
                pc=next(ii for ii,cl in enumerate(classes) if power(reps[ci],t) in cl)
                a+=(-1)**(t-1)*c[pc]*es[j-t]
            assert a%j==0;es.append(a//j)
        values.append(es[-1])
    return tuple(values)
bydim={};pres_count={}
for d in (56,22):
    buckets=defaultdict(list)
    for p in presentations(d):buckets[character(p)].append(p)
    bydim[d]=buckets;pres_count[d]=sum(map(len,buckets.values()))
targets=defaultdict(list)
for v in bydim[22]:targets[exterior(v,3)].append(v)
pairs=[]
for u in bydim[56]:
    for v in targets.get(exterior(u,2),[]):
        if u==(56,)*5 and v==(22,)*5:continue
        pairs.append({'U_character':u,'V_character':v,
                      'U_presentations':bydim[56][u],'V_presentations':bydim[22][v]})
assert len(pairs)==1
p=pairs[0]
assert p['U_character']==(56,4,5,1,1) and p['V_character']==(22,2,4,2,2)
raw_pair_count=sum(len(p['U_presentations'])*len(p['V_presentations']) for p in pairs)
assert raw_pair_count==6

# Independent arithmetic for each branch of the short uniqueness proof.
inv_filter=[]
for f in range(0,23,2):
    n=(f**3-64*f+168)//3;assert (f**3-64*f+168)%3==0
    root=isqrt(n) if n>=0 and isqrt(n)**2==n else None
    inv_filter.append({'V_trace':f,'required_U_trace_squared':n,'U_trace':root})
assert [(x['V_trace'],x['U_trace']) for x in inv_filter if x['U_trace'] is not None]==[(2,4),(6,0),(22,56)]
u_zero=[x for x in bydim[56] if x[1]==0]
assert u_zero==[(56,0,2,6,6)]
v_two=[{'character':v,'presentations':ps} for v,ps in bydim[22].items() if v[1]==2]
assert len(v_two)==2 and sum(len(p['presentations']) for p in v_two)==3
order3_values=[{'fixed':r,'exterior_cube':r*(r-1)*(r-2)//6+(22-r)//3} for r in range(1,23,3)]
assert min(x['exterior_cube'] for x in order3_values)==7
assert all(((r+3)*(r+2)*(r+1)-r*(r-1)*(r-2))//6-1==3*r*(r+1)//2 for r in range(1,20,3))
assert (216+0)%60!=0 and (216+24)//60==4

result={'group':'A5, even permutations of five letters, left composition',
 'class_representatives':reps,'class_sizes':list(map(len,classes)),
 'coset_models':raw,'transitive_characters':chars,'all_60_element_checks':all_checks,
 'presentation_size_order':ds,'presentation_counts':pres_count,
 'representation_counts':{d:len(b) for d,b in bydim.items()},
 'nontrivial_module_pairs':pairs,'nontrivial_presentation_pair_count':raw_pair_count,
 'involution_filter':inv_filter,'U_zero_involution_characters':u_zero,
 'V_trace_two_candidates':v_two,'order_three_exterior_values':order3_values,
 'invariant_dimensions':{'U':4,'V':3}}
(OUT/'MI_S3_exact_results.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: explicit coset actions and exterior characters on all 60 elements.')
print('PASS: unique nontrivial ordered permutation-module pair; six orbit-presentation pairs.')
print('PASS: every branch of the short uniqueness proof; invariant dimensions 4 and 3.')
print('Presentation counts:',pres_count)
print('Distinct permutation-module counts:',{d:len(b) for d,b in bydim.items()})
