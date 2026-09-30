#!/usr/bin/env python3
"""MI-S0--S2, independent exact reconstruction. Python standard library only.

No received script or stored mathematical answer is imported. Construct every
subgroup and every linear character; form exterior characters by multiplying
eigenvalue factors in Z[z]/(z^r-1), independently of Newton recurrence.
All coordinate tuples use (1,D_plus,D_minus,A,B); a+b*phi is stored as [a,b].
"""
from itertools import permutations, combinations
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from math import gcd, comb
import json

OUT = Path(__file__).resolve().parent
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def mul(x,y): return (x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def scale(n,x): return (n*x[0],n*x[1])
def divide(x,n):
    assert all(a%n==0 for a in x)
    return tuple(a//n for a in x)
def ch(v):
    t,x,y,a,b=v
    return ((t+3*x+3*y+4*a+5*b,0),(t-x-y+b,0),(t+a-b,0),
            (t-a+y,x-y),(t-a+x,y-x))
IRR=[ch(tuple(int(i==j) for i in range(5))) for j in range(5)]
SIZES=(1,15,20,12,12)
def decomp(c):
    result=[]
    for r in IRR:
        z=(0,0)
        for a,b,n in zip(c,r,SIZES): z=add(z,scale(n,mul(a,b)))
        z=divide(z,60);assert z[1]==0 and z[0]>=0
        result.append(z[0])
    assert ch(result)==tuple(c)
    return tuple(result)
assert all(decomp(r)==tuple(int(i==j) for i in range(5)) for j,r in enumerate(IRR))

# Group construction uses left composition p*q, q applied first.
G=sorted(p for p in permutations(range(5)) if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0)
ix={p:i for i,p in enumerate(G)}; ID=ix[tuple(range(5))]
M=[[ix[tuple(p[q[i]] for i in range(5))] for q in G] for p in G]
INV=[next(j for j in range(60) if M[i][j]==ID) for i in range(60)]
def pw(g,k):
    r=ID
    for _ in range(k):r=M[r][g]
    return r
def conjugate(g,h): return M[M[g][h]][INV[g]]
def closure(gens):
    seen={ID};queue=[ID]
    for x in queue:
        for g in gens:
            y=M[x][g]
            if y not in seen:seen.add(y);queue.append(y)
    return frozenset(seen)
orders=[next(k for k in range(1,6) if pw(g,k)==ID) for g in range(60)]
remaining=set(range(60));classes=[]
while remaining:
    g=min(remaining);c=frozenset(conjugate(h,g) for h in range(60));classes.append(c);remaining-=c
classes.sort(key=lambda c:(orders[min(c)],min(c)))
assert tuple(map(len,classes))==SIZES
reps=[min(c) for c in classes]
classof={g:i for i,c in enumerate(classes) for g in c}
subgens={frozenset([ID]):()};queue=[frozenset([ID])]
for H in queue:
    for g in range(60):
        if g not in H:
            K=closure(subgens[H]+(g,))
            if K not in subgens:subgens[K]=subgens[H]+(g,);queue.append(K)
subs=sorted(subgens,key=lambda H:(len(H),tuple(sorted(H))))
assert len(subs)==59

def root_sum(a):
    """Reduce a symmetric sum of r-th roots to Z[phi]."""
    if len(a)==1:return (a[0],0)
    if len(a)==2:return (a[0]-a[1],0)
    if len(a)==3:
        assert a[1]==a[2]
        return (a[0]-a[1],0)
    assert len(a)==5 and a[1]==a[4] and a[2]==a[3]
    return (a[0]-a[1],a[1]-a[2])

# Every linear character has a kernel K with cyclic quotient H/K. Inspect all
# normal K and all possible quotient orders, not a preselected list of primes.
catalogue=[];generators={n:set() for n in ['permutation','signed','monomial']}
subdata=[]
for hi,H in enumerate(subs):
    linears=[]
    for K in subs:
        if not K<=H:continue
        if any(conjugate(h,k) not in K for h in H for k in K):continue
        L=len(H)//len(K)
        generator=next((g for g in sorted(H) if closure(tuple(K)+(g,))==H),None)
        if generator is None:continue
        vals={M[pw(generator,j)][k]:j for j in range(L) for k in K}
        assert len(vals)==len(H)
        for a in range(L):
            if gcd(a,L)!=1:continue
            exps={h:(a*vals[h])%L for h in H}
            assert all(exps[M[h][k]]==(exps[h]+exps[k])%L for h in H for k in H)
            linears.append((L,exps))
    assert len({tuple(v[h] for h in sorted(H))+(L,) for L,v in linears})==len(linears)
    comm=closure(tuple(M[M[M[h][k]][INV[h]]][INV[k]] for h in H for k in H))
    assert len(linears)==len(H)//len(comm)
    for L,vals in linears:
        cc=[]
        for c in classes:
            a=[0]*L
            for h in H&c:a[vals[h]]+=1
            cc.append(divide(scale(60//len(c),root_sum(a)),len(H)))
        v=decomp(cc);generators['monomial'].add(v)
        if L<=2:generators['signed'].add(v)
        if L==1:generators['permutation'].add(v)
        catalogue.append({'subgroup':hi,'character_order':L,'exponents':[[h,vals[h]] for h in sorted(H)],'multiplicities':v,'character':cc})
    subdata.append({'id':hi,'elements':sorted(H),'generators':subgens[H],
                    'order':len(H),'commutator_order':len(comm),'linear_characters':len(linears)})
generators={n:sorted(v) for n,v in generators.items()}
assert [len(generators[n]) for n in generators]==[9,13,17]
assert len(catalogue)==159

def permutation_from_cycles(*cycles):
    p=list(range(5))
    for c in cycles:
        c=tuple(x-1 for x in c)
        for a,b in zip(c,c[1:]+c[:1]):p[a]=b
    return ix[tuple(p)]
c5=permutation_from_cycles((1,2,3,4,5));c3=permutation_from_cycles((1,2,3))
H5=closure((c5,)); H3=closure((c3,))
H10=frozenset(g for g in range(60) if frozenset(conjugate(g,h) for h in H5)==H5)
H12=frozenset(g for g,p in enumerate(G) if p[4]==4)
H6=frozenset(g for g,p in enumerate(G) if set(p[:2])=={0,1})
actions={};cosetdata={}
for d,H in [(5,H12),(6,H10),(10,H6),(12,H5),(20,H3),(60,frozenset([ID]))]:
    left=set(range(60));cosets=[]
    while left:
        g=min(left);c=frozenset(M[g][h] for h in H);cosets.append(c);left-=c
    index={h:i for i,c in enumerate(cosets) for h in c}
    act=[tuple(index[M[g][min(c)]] for c in cosets) for g in range(60)]
    assert len(cosets)==d
    assert all(tuple(act[g][act[h][i]] for i in range(d))==act[M[g][h]] for g in range(60) for h in range(60))
    cc=tuple((sum(i==v for i,v in enumerate(act[g])),0) for g in reps)
    actions[d]=act;cosetdata[d]={'stabilizer_elements':sorted(H),'cosets':[sorted(c) for c in cosets],
                               'actions':act,'character':cc,'multiplicities':decomp(cc)}

# Oriented triples in Q: stabilizer signs and all five orbits explicitly.
triples=set(combinations(range(12),3));triple_orbits=[]
while triples:
    s=min(triples);orbit={tuple(sorted(actions[12][g][i] for i in s)) for g in range(60)}
    stab=[g for g in range(60) if tuple(sorted(actions[12][g][i] for i in s))==s]
    signs=[(-1)**sum(actions[12][g][s[i]]>actions[12][g][s[j]] for i in range(3) for j in range(i+1,3)) for g in stab]
    assert set(signs)=={1}
    triple_orbits.append({'representative':s,'size':len(orbit),'stabilizer':stab,'signs':signs,'elements':sorted(orbit)})
    triples-=orbit
assert sorted(x['size'] for x in triple_orbits)==[20,20,60,60,60]

# Eigenvalue exponents on (1,D+,D-,A,B); class 5B doubles 5A exponents.
EIG={1:[[0],[0,0,0],[0,0,0],[0,0,0,0],[0,0,0,0,0]],
     2:[[0],[0,1,1],[0,1,1],[0,0,1,1],[0,0,0,1,1]],
     3:[[0],[0,1,2],[0,1,2],[0,0,1,2],[0,1,1,2,2]],
     5:[[0],[0,1,4],[0,2,3],[1,2,3,4],[0,1,2,3,4]]}
def exterior(v,k):
    cc=[]
    for ci,r in enumerate((1,2,3,5,5)):
        eig=[(e*(2 if ci==4 else 1))%r for n,ee in zip(v,EIG[r]) for _ in range(n) for e in ee]
        dp=[[0]*r for _ in range(k+1)];dp[0][0]=1
        for i,e in enumerate(eig):
            for j in range(min(k,i+1),0,-1):
                for a in range(r):dp[j][(a+e)%r]+=dp[j-1][a]
        cc.append(root_sum(dp[k]))
    return tuple(cc)
def square(v):
    c=ch(v)
    return tuple(divide(add(mul(x,x),scale(-1,c[j])),2) for x,j in zip(c,(0,0,2,4,3)))
def mods(n):
    for x in range(n//3+1):
        for y in range((n-3*x)//3+1):
            for a in range((n-3*x-3*y)//4+1):
                for b in range((n-3*x-3*y-4*a)//5+1):yield (n-3*x-3*y-4*a-5*b,x,y,a,b)
def witness_factory(gens):
    @lru_cache(None)
    def f(v,i=0):
        if not any(v):return (0,)*(len(gens)-i)
        if i==len(gens):return None
        g=gens[i];bound=min(v[j]//g[j] for j in range(5) if g[j])
        for a in range(bound,-1,-1):
            z=f(tuple(v[j]-a*g[j] for j in range(5)),i+1)
            if z is not None:return (a,)+z
        return None
    return f
members={n:witness_factory(g) for n,g in generators.items()}
ROWS=[(120,16,10,3),(210,21,10,4),(1540,56,22,3),(3003,78,15,5),(7140,120,36,3),(11628,153,19,5),(24310,221,17,8)]
results={}
for val,m,n,k in ROWS:
    targets=defaultdict(list)
    for v in mods(n):targets[exterior(v,k)].append(v)
    early={(c[1][0],c[2][0]) for c in targets};pairs=[];nu=0
    for u in mods(m):
        nu+=1;c=ch(u);s=c[1][0];t=c[2][0]
        if ((s*s-m)//2,(t*t-t)//2) not in early:continue
        for v in targets.get(square(u),[]):
            if u==(m,0,0,0,0) and v==(n,0,0,0,0):continue
            assert exterior(u,2)==exterior(v,k)
            witnesses={name:{'U':f(u),'V':f(v)} for name,f in members.items()}
            pairs.append({'U':u,'V':v,'exterior_character':square(u),'rational':u[1]==u[2] and v[1]==v[2],'witnesses':witnesses})
    counts={'complex':len(pairs),'rational':sum(p['rational'] for p in pairs)}
    counts.update({name:sum(all(x is not None for x in p['witnesses'][name].values()) for p in pairs) for name in members})
    results[val]={'dimensions':[m,n],'degree':k,'U_candidates':nu,'V_candidates':sum(1 for _ in mods(n)),'counts':counts,'pairs':pairs}
    print(val,counts,flush=True)

def full(v):a,b,c,d=v;return(a,d,d,b,c)
newrows={1540:((4,5,4,2),(3,2,1,1)),7140:((3,9,9,6),(1,1,5,1)),24310:((6,16,17,11),(3,1,2,0))}
oldrows={120:((0,0,2,1),(1,1,1,0)),210:((3,2,2,0),(1,1,1,0)),1540:((0,2,6,3),(4,2,2,0)),3003:((4,4,8,3),(0,1,1,1)),7140:((3,9,9,6),(1,1,5,1)),11628:((6,9,15,6),(1,0,0,3)),24310:((3,13,20,11),(2,0,3,0))}
for rows in [newrows,oldrows]:
    for val,(u,v) in rows.items():assert any(p['U']==full(u) and p['V']==full(v) for p in results[val]['pairs'])
roots=[v for v in mods(16) if square(v)==ch(full(newrows[7140][0]))]
assert roots==[(1,0,0,0,3)]
assert all(ch(p['V'])[2]==(-3,0) for p in results[7140]['pairs'])
aux={}
for name,v in [('W',(1,0,0,0,3)),('T',(0,1,1,1,2)),('Lambda3Q',decomp(exterior((1,1,1,0,1),3)))]:
    aux[name]={'multiplicities':v,'character':ch(v)}
assert aux['Lambda3Q']['multiplicities']==(5,11,11,16,17)
def vsum(*vs):return tuple(sum(v[i] for v in vs) for i in range(5))
R={d:tuple(data['multiplicities']) for d,data in cosetdata.items()}
one=(1,0,0,0,0);W=aux['W']['multiplicities'];T=aux['T']['multiplicities']
assert tuple(mul(a,b) for a,b in zip(IRR[3],IRR[4]))==ch(T)
assert vsum(R[6],R[10],R[20],R[20])==full(newrows[1540][0])
assert vsum(R[5],R[5],R[12])==full(newrows[1540][1])
assert vsum(one,aux['Lambda3Q']['multiplicities'])==full(newrows[24310][0])
assert vsum(R[5],R[6],R[6])==full(newrows[24310][1])
assert aux['Lambda3Q']['multiplicities']==vsum(R[20],R[20],R[60],R[60],R[60])
assert vsum(W,T)==full(newrows[7140][1])
assert tuple(3*x for x in vsum(R[20],T))==full(newrows[7140][0])
assert exterior(vsum(one,W),8)==exterior(full(newrows[24310][1]),8)
assert vsum(one,W)!=full(newrows[24310][1])
for H,target in [(H12,(0,0,0,0,1)),(H3,T)]:
    cs=[x for x in catalogue if subs[x['subgroup']]==H and x['character_order']==3]
    assert len(cs)==2 and all(x['multiplicities']==target for x in cs)
construction_witnesses={val:{'U':full(u),'V':full(v),'U_character':ch(full(u)),
                            'V_character':ch(full(v)),'exterior_character':square(full(u))}
                        for val,(u,v) in newrows.items()}
record={'tuple_order':['1','D_plus','D_minus','A','B'],'class_order':['1','2A','3A','5A','5B'],
        'group_permutations':G,'class_representatives':reps,'subgroups':subdata,
        'subgroups_by_order':dict(sorted(Counter(map(len,subs)).items())),
        'linear_character_total':len(catalogue),'induction_catalogue':catalogue,
        'generators':generators,'coset_actions':cosetdata,'vertex_triple_orbits':triple_orbits,
        'rows':results,'old_selected_rows':oldrows,'new_selected_rows':newrows,
        'root_candidates':sum(1 for _ in mods(16)),'exterior_square_roots':roots,'auxiliary':aux,
        'construction_witnesses':construction_witnesses}
(OUT/'wave_i_results.json').write_text(json.dumps(record,indent=2)+'\n')
print('PASS: exact subgroup, induction, coset, exterior, counting, and preservation checks.')
