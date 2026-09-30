#!/usr/bin/env python3
"""MI-S5 exact reconstruction, using Python's standard library only.

Matrices are monomial column actions: e_i -> zeta6^a_i e_p(i), where
zeta6=-omega^2=1+omega and omega^2+omega+1=0. Cubic actions have even a_i.
No mathematical result or source code from previous stages is imported.
"""
from itertools import permutations, combinations, product
from pathlib import Path
from math import comb
import json

OUT = Path(__file__).resolve().parent
def compose(p,q): return tuple(p[q[i]] for i in range(len(p)))
def inverse(p): return tuple(p.index(i) for i in range(len(p)))
def parity(p): return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
ID = tuple(range(5))
G = sorted(p for p in permutations(range(5)) if parity(p)==1)
def perm(*cycles):
    p = list(ID)
    for c in cycles:
        c = tuple(i-1 for i in c)
        for a,b in zip(c,c[1:]+c[:1]): p[a]=b
    return tuple(p)
def closure(gens):
    seen={ID}; todo=[ID]
    for p in todo:
        for q in gens:
            r=compose(p,q)
            if r not in seen: seen.add(r); todo.append(r)
    return frozenset(seen)
def power(g,k):
    r=tuple(range(len(g)))
    for _ in range(k): r=compose(r,g)
    return r
def order(g): return next(k for k in (1,2,3,5) if power(g,k)==tuple(range(len(g))))
left=set(G); classes=[]
while left:
    g=min(left); cl={compose(compose(x,g),inverse(x)) for x in G}
    classes.append(cl); left-=cl
classes.sort(key=lambda c:(order(min(c)),min(c)))
reps=[min(c) for c in classes]
assert list(map(len,classes))==[1,15,20,12,12]

# Cubic character of A4 from its action on the three pair partitions.
H4=frozenset(g for g in G if g[4]==4)
pairings=sorted({tuple(sorted((tuple(c),tuple(i for i in range(4) if i not in c))))
                 for c in combinations(range(4),2)})
def pairing_action(g):
    return tuple(pairings.index(tuple(sorted(tuple(sorted(g[i] for i in c)) for c in p)))
                 for p in pairings)
c=perm((1,2,3)); p0=pairing_action(c)
eta={g:next(a for a in range(3) if pairing_action(g)==power(p0,a)) for g in H4}
H3=closure((c,)); omega={power(c,a):a for a in range(3)}
assert len(H4)==12 and len(H3)==3
assert sorted(list(eta.values()).count(a) for a in range(3))==[4,4,4]

def induced(group,H,exponents):
    identity=tuple(range(len(group[0])))
    assert all((exponents[g]+exponents[h])%3==exponents[compose(g,h)] for g in H for h in H)
    remaining=set(group); cosets=[]
    while remaining:
        r=min(remaining); cs=frozenset(compose(r,h) for h in H)
        cosets.append(cs); remaining-=cs
    rs=[min(cs) for cs in cosets]; ix={g:i for i,cs in enumerate(cosets) for g in cs}
    result={}
    for g in group:
        ps=[]; es=[]
        for r in rs:
            gr=compose(g,r); j=ix[gr]; h=compose(inverse(rs[j]),gr)
            assert h in H
            ps.append(j); es.append(2*exponents[h]%6)
        result[g]=(tuple(ps),tuple(es))
    return result,rs
def direct_sum(*acts):
    result={}
    for g in acts[0]:
        ps=[]; es=[]; offset=0
        for act in acts:
            p,e=act[g]; ps.extend(x+offset for x in p); es.extend(e); offset+=len(p)
        result[g]=(tuple(ps),tuple(es))
    return result
def action_product(a,b):
    p,e=a; q,f=b
    return compose(p,q),tuple((f[i]+e[q[i]])%6 for i in range(len(p)))
def verify_law(act,group):
    assert all(action_product(act[g],act[h])==act[compose(g,h)] for g in group for h in group)
def wedge2(act):
    pairs=list(combinations(range(len(next(iter(act.values()))[0])),2)); ix={p:i for i,p in enumerate(pairs)}
    out={}
    for g,(p,e) in act.items():
        ps=[]; es=[]
        for i,j in pairs:
            a,b=p[i],p[j]; ps.append(ix[tuple(sorted((a,b)))])
            es.append((e[i]+e[j]+(3 if a>b else 0))%6)
        out[g]=(tuple(ps),tuple(es))
    return out
B,b_reps=induced(G,H4,eta)
Bbar,_=induced(G,H4,{g:-a%3 for g,a in eta.items()})
T,t_reps=induced(G,H3,omega)
Tbar,_=induced(G,H3,{g:-a%3 for g,a in omega.items()})
R20,_=induced(G,H3,{g:0 for g in H3})
One={g:((0,),(0,)) for g in G}
W=direct_sum(One,B,B,B)
U_wedge=wedge2(W)
U_cubic=direct_sum(R20,R20,R20,T,T,T)
V=direct_sum(W,T)
for act in (B,T,R20,W,U_wedge,U_cubic,V): verify_law(act,G)
assert all(a%2==0 for act in (B,T,W,U_cubic,V) for _,es in act.values() for a in es)

# Arithmetic in Z[omega] for traces and weighted-cycle determinant products.
ZERO=(0,0); ONE=(1,0)
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def scale(n,x): return (n*x[0],n*x[1])
def mul(x,y): return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
ROOTS=[(1,0),(1,1),(0,1),(-1,0),(-1,-1),(0,-1)]
def trace(a):
    p,es=a; x=ZERO
    for i,j in enumerate(p):
        if i==j: x=add(x,ROOTS[es[i]])
    return x
def determinant_coefficients(a,degree=None):
    p,es=a; degree=len(p) if degree is None else degree
    left=set(range(len(p))); out=[ONE]+[ZERO]*degree
    while left:
        i=min(left); n=0; phase=0
        while i in left: left.remove(i); n+=1; phase+=es[i]; i=p[i]
        z=scale((-1)**(n+1),ROOTS[phase%6])
        for j in range(degree,n-1,-1): out[j]=add(out[j],mul(z,out[j-n]))
    return out
def newton(act,g,k):
    es=[ONE]
    for j in range(1,k+1):
        z=ZERO
        for r in range(1,j+1): z=add(z,scale((-1)**(r-1),mul(trace(act[power(g,r)]),es[j-r])))
        assert z[0]%j==0 and z[1]%j==0
        es.append((z[0]//j,z[1]//j))
    return es[k]
def integral(x): assert x[1]==0; return x[0]
table={name:[integral(trace(act[g])) for g in reps]
       for name,act in [('B',B),('T',T),('W',W),('U',U_cubic),('V',V)]}
assert table=={'B':[5,1,-1,0,0],'T':[20,0,-1,0,0],'W':[16,4,-2,1,1],
               'U':[120,0,3,0,0],'V':[36,4,-3,1,1]}
assert all(trace(B[g])==trace(Bbar[g]) and trace(T[g])==trace(Tbar[g]) for g in G)
assert all(trace(T[g])==scale(sum(g[i]==i for i in range(5))-1,trace(B[g])) for g in G)
checks=[]
for g in G:
    assert trace(U_wedge[g])==trace(U_cubic[g])==newton(W,g,2)
    l=determinant_coefficients(U_wedge[g],2)[2]
    r=determinant_coefficients(V[g],3)[3]
    assert l==r==newton(U_cubic,g,2)==newton(V,g,3)
    checks.append({'g':g,'order':order(g),'W':integral(trace(W[g])),
                   'U':integral(trace(U_cubic[g])),'V':integral(trace(V[g])),
                   'common_exterior_trace':integral(l)})
common=[integral(newton(V,g,3)) for g in reps]
assert common==[7140,-60,3,0,0]
polynomials={name:[[integral(x) for x in determinant_coefficients(act[g])] for g in reps]
             for name,act in [('W',W),('U',U_cubic),('V',V)]}

# Independent exterior-square-root search, allowing unequal 3-dimensional multiplicities.
# Here [a,b] means a+b*phi, with phi^2=phi+1, not a+b*omega.
def phimul(x,y): return (x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def irreducible_character(v):
    t,x,y,a,b=v
    return [(t+3*x+3*y+4*a+5*b,0),(t-x-y+b,0),(t+a-b,0),
            (t-a+y,x-y),(t-a+x,y-x)]
def square_char(ch):
    ans=[]
    for ci,sq in enumerate((0,0,2,4,3)):
        z=add(phimul(ch[ci],ch[ci]),scale(-1,ch[sq]))
        assert all(a%2==0 for a in z); ans.append(tuple(a//2 for a in z))
    return ans
candidates=[]; roots=[]
for x in range(6):
    for y in range(6):
        for a in range(5):
            for b in range(4):
                t=16-3*x-3*y-4*a-5*b
                if t<0: continue
                v=(t,x,y,a,b); candidates.append(v)
                if square_char(irreducible_character(v))==[(a,0) for a in table['U']]: roots.append(v)
assert len(candidates)==73 and roots==[(1,0,0,0,3)]
article_coordinates={'W':(1,0,3,0),'T':(0,1,2,1),'U':(3,9,9,6),'V':(1,1,5,1)}
for name,(a,b,c,d) in article_coordinates.items():
    assert irreducible_character((a,d,d,b,c))==[(n,0) for n in table[name]]
root_filter=[{'r':r,'q':q,'A_multiplicity_numerator':64+20*r-24*q} for r in (3,-2) for q in (0,1)]
assert [x for x in root_filter if x['A_multiplicity_numerator']%60==0]==[root_filter[-1]]

# Geometric tangent characters from exact icosahedral rotations in Z[phi].
def sub(x,y): return add(x,scale(-1,y))
def dot(x,y):
    z=ZERO
    for a,b in zip(x,y): z=add(z,phimul(a,b))
    return z
def cross(a,b):
    return (sub(phimul(a[1],b[2]),phimul(a[2],b[1])),
            sub(phimul(a[2],b[0]),phimul(a[0],b[2])),
            sub(phimul(a[0],b[1]),phimul(a[1],b[0])))
def det(a,b,c): return dot(a,cross(b,c))
def sign_phi(z):
    # 2(a+b*phi)=(2a+b)+b*sqrt(5); exact sign by comparison of squares.
    a,b=z; u=2*a+b
    if b==0: return (u>0)-(u<0)
    if u==0: return (b>0)-(b<0)
    if (u>0)==(b>0): return (u>0)-(u<0)
    s=u*u-5*b*b
    return ((u>0)-(u<0))*((s>0)-(s<0))
vertices=sorted({v for s,t in product((-1,1),repeat=2) for v in
                 (((0,0),(s,0),(0,t)),((s,0),(0,t),(0,0)),((0,t),(0,0),(s,0)))})
def dist(i,j): return dot(tuple(sub(a,b) for a,b in zip(vertices[i],vertices[j])),
                         tuple(sub(a,b) for a,b in zip(vertices[i],vertices[j])))
faces=[t for t in combinations(range(12),3) if all(dist(i,j)==(4,0) for i,j in combinations(t,2))]
assert len(faces)==20
base=faces[0]; src=[vertices[i] for i in base]; srcdet=det(*src)
gram=[[dot(a,b) for b in src] for a in src]
srcsig=[tuple(dot(v,b) for b in src) for v in vertices]
rotations=set()
for dst in permutations(range(12),3):
    bs=[vertices[i] for i in dst]
    if det(*bs)!=srcdet or [[dot(a,b) for b in bs] for a in bs]!=gram: continue
    ix={tuple(dot(v,b) for b in bs):i for i,v in enumerate(vertices)}
    assert len(ix)==12 and all(x in ix for x in srcsig)
    rotations.add(tuple(ix[x] for x in srcsig))
rotations=sorted(rotations)
assert len(rotations)==60 and all(compose(g,h) in rotations for g in rotations for h in rotations)
oriented=[]
for f in faces:
    a,b,c=f
    if sign_phi(det(vertices[a],vertices[b],vertices[c]))<0: b,c=c,b
    oriented.append((a,b,c))
    u=tuple(sub(x,y) for x,y in zip(vertices[b],vertices[a]))
    v=tuple(sub(x,y) for x,y in zip(vertices[c],vertices[b]))
    n=tuple(add(add(x,y),z) for x,y,z in zip(vertices[a],vertices[b],vertices[c]))
    assert dot(n,u)==dot(n,v)==ZERO
    assert dot(u,u)==dot(v,v)==(4,0) and dot(u,v)==(-2,0)
    assert sign_phi(dot(n,cross(u,v)))>0
tangent_checks=[]
for g in rotations:
    phases=[]
    for f in oriented:
        image=tuple(g[i] for i in f)
        if set(image)!=set(f): continue
        shift=next(k for k in range(3) if image==f[k:]+f[:k]); phases.append(shift)
    tr=ZERO
    for a in phases: tr=add(tr,ROOTS[2*a])
    assert integral(tr)=={1:20,2:0,3:-1,5:0}[order(g)]
    if order(g)==3: assert sorted(phases)==[1,2]
    tangent_checks.append({'rotation':g,'order':order(g),'fixed_face_phases_mod_3':phases,'trace':integral(tr)})

def stored_action(act):
    return [{'g':g,'permutation':p,'exponents_mod_6':e} for g,(p,e) in sorted(act.items())]
result={'conventions':{'group':'even permutations of five letters, left composition, zero-based labels',
                       'class_order':['1','2A','3A','5A','5B'],
                       'monomial_exponents':'zeta6=-omega^2; cubic phases have even exponents',
                       'geometric_scalars':'[a,b] = a+b*phi, phi^2=phi+1'},
        'class_representatives':reps,'pair_partitions':pairings,
        'A4_character':[{'g':g,'exponent_mod_3':a} for g,a in sorted(eta.items())],
        'C3_character':[{'g':g,'exponent_mod_3':a} for g,a in sorted(omega.items())],
        'coset_representatives':{'B':b_reps,'T':t_reps},
        'induced_actions':{'B':stored_action(B),'T':stored_action(T),'R20':stored_action(R20)},
        'characters':table,'article_coordinates':article_coordinates,
        'determinant_polynomials_by_class':polynomials,
        'all_60_element_checks':checks,'common_exterior_character':common,
        'root_candidates_in_code_order':candidates,'exterior_square_roots':roots,'root_filter':root_filter,
        'geometric_vertices':vertices,'outward_oriented_faces':oriented,'tangent_line_checks':tangent_checks}
(OUT/'MI_S5_exact_results.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: both cubic stabilizer characters and their induced actions, with full group laws.')
print('PASS: U=Lambda^2 W=3R20+3T and V=W+T on all 60 elements.')
print('PASS: weighted-cycle and Newton calculations; common character',common)
print('PASS: cube-root bases for W,U,V; V has order-three trace -3, excluding signed bases.')
print('PASS: all 73 sixteen-dimensional modules checked; unique exterior-square root 1+3B.')
print('PASS: exact tangent-line character from 60 rotations; fixed phases omega and omega^2.')
