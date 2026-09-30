#!/usr/bin/env python3
"""Exact finite checks through MI-S6, using integers in Z[phi] (phi^2=phi+1).
No numerical character evaluation is used. Run from any directory.
The finite search ranges are dimension equations, not empirical cutoffs.
"""
from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import comb
from pathlib import Path
import json, time

ROOT=Path(__file__).resolve().parent
# Algebraic integers a+b*phi.
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def mul(x,y): return (x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def smul(c,x): return (c*x[0],c*x[1])
def div(x,n):
    assert x[0]%n==x[1]%n==0, (x,n)
    return (x[0]//n,x[1]//n)
Z=(0,0); ONE=(1,0)
# 1,D,D',A,B, classes 1,2A,3A,5A,5B
IRR=[[(1,0)]*5,[(3,0),(-1,0),Z,(0,1),(1,-1)],
     [(3,0),(-1,0),Z,(1,-1),(0,1)],
     [(4,0),Z,(1,0),(-1,0),(-1,0)],[(5,0),(1,0),(-1,0),Z,Z]]
SIZES=[1,15,20,12,12]
def pc(c,i):
    if c==0:return 0
    if c==1:return 0 if i%2==0 else 1
    if c==2:return 0 if i%3==0 else 2
    return (0,c,7-c,7-c,c)[i%5]
def character(v):
    t,x,y,a,b=v
    return ((t+3*x+3*y+4*a+5*b,0),(t-x-y+b,0),(t+a-b,0),
            (t-a+y,x-y),(t-a+x,y-x))
def exterior_ch(ch,k):
    ans=[]
    for c in range(5):
        e=[ONE]
        for r in range(1,k+1):
            z=Z
            for i in range(1,r+1):z=add(z,smul((-1)**(i-1),mul(e[r-i],ch[pc(c,i)])))
            e.append(div(z,r))
        ans.append(e[k])
    return tuple(ans)
def modules(n):
    for x in range(n//3+1):
        for y in range((n-3*x)//3+1):
            rem=n-3*x-3*y
            for a in range(rem//4+1):
                for b in range((rem-4*a)//5+1):
                    yield (rem-4*a-5*b,x,y,a,b)
def enum_pairs(n,m,k):
    buckets=defaultdict(list)
    for v in modules(m):buckets[exterior_ch(character(v),k)].append(v)
    early={(ch[1][0],ch[2][0]) for ch in buckets}
    pairs=[]; nU=0
    for u in modules(n):
        nU+=1;t,x,y,a,b=u
        c2=t-x-y+b;c3=t+a-b
        if ((c2*c2-n)//2,(c3*c3-c3)//2) not in early:continue
        ch=character(u)
        e2=tuple(div(add(mul(ch[c],ch[c]),smul(-1,ch[pc(c,2)])),2) for c in range(5))
        pairs.extend((u,v) for v in buckets.get(e2,[]))
    trivial=((n,0,0,0,0),(m,0,0,0,0))
    assert trivial in pairs
    return [p for p in pairs if p!=trivial], nU, sum(1 for _ in modules(m))
def scalar(ch,irr):
    z=Z
    for n,a,b in zip(SIZES,ch,irr):z=add(z,smul(n,mul(a,b)))
    q=div(z,60); assert q[1]==0
    return q[0]
ROWS=[(120,16,10,3),(210,21,10,4),(1540,56,22,3),(3003,78,15,5),
      (7140,120,36,3),(11628,153,19,5),(24310,221,17,8)]
PRINTED=[((0,0,2,1),(1,1,1,0)),((3,2,2,0),(1,1,1,0)),
         ((4,5,4,2),(3,2,1,1)),((4,4,8,3),(0,1,1,1)),
         ((3,9,9,6),(1,1,5,1)),((6,9,15,6),(1,0,0,3)),
         ((6,16,17,11),(3,1,2,0))]
EARLIER={1540:((0,2,6,3),(4,2,2,0)),24310:((3,13,20,11),(2,0,3,0))}
def full(v):t,a,b,c=v;return(t,c,c,a,b)

def weights_ext(ws,k):
    dp=[Counter() for _ in range(k+1)];dp[0][0]=1
    for i,w in enumerate(ws):
        for r in range(min(k,i+1),0,-1):
            for v,m in list(dp[r-1].items()):dp[r][v+w]+=m
    return dp[k]
def wt_list(counter):return [w for w,m in counter.items() for _ in range(m)]
def harmonics(wt):
    assert all(wt[-j]==wt[j] for j in wt)
    out={j:wt[j]-wt[j+1] for j in range(max(wt)+1)}
    assert all(n>=0 for n in out.values())
    return out
def reduce_h(j):
    if j<6:return [int(r==j) for r in range(6)]
    a=reduce_h(j-6);r=j%5;a[r]+=1;a[5-r]+=1;return a
def reduce_sum(h):return [sum(m*reduce_h(j)[i] for j,m in h.items()) for i in range(6)]
def hchar(j):
    ns=Counter(w%5 for w in range(-j,j+1));assert ns[1]==ns[4] and ns[2]==ns[3]
    th=Counter(w%3 for w in range(-j,j+1));assert th[1]==th[2]
    return ((2*j+1,0),(sum(1 if w % 2 == 0 else -1 for w in range(-j,j+1)),0),
       (th[0]-th[1],0),(ns[0]-ns[1],ns[1]-ns[2]),(ns[0]-ns[2],ns[2]-ns[1]))
def soh_models(dim,target,j=0):
    # exhaustive sums of H_j in a fixed dimension; only j <= (dim-1)/2 can occur.
    ans=[]
    def rec(remaining,start,counts):
        if remaining==0:
            ch=tuple((sum(counts.get(r,0)*hchar(r)[c][0] for r in counts),
                      sum(counts.get(r,0)*hchar(r)[c][1] for r in counts)) for c in range(5))
            if ch==target:ans.append(counts.copy())
            return
        for r in range(start,(remaining-1)//2+1):
            counts[r]=counts.get(r,0)+1;rec(remaining-2*r-1,r,counts)
            counts[r]-=1
            if not counts[r]:del counts[r]
    rec(dim,0,{})
    return ans

def cycles_perm(n,cs):
    p=list(range(n))
    for cyc in cs:
        for x,y in zip(cyc,cyc[1:]+cyc[:1]):p[x-1]=y-1
    return tuple(p)
def rprod(p,q):return tuple(q[p[i]] for i in range(len(p)))
def inv(p):return tuple(p.index(i) for i in range(len(p)))
def word(w,a,b):
    p=tuple(range(len(a)))
    for x in w:p=rprod(p,{'a':a,'b':b}[x])
    return p
def perm_order(p):
    from math import lcm
    left=set(range(len(p)));ls=[]
    while left:
        a=next(iter(left));n=0
        while a in left:left.remove(a);n+=1;a=p[a]
        ls.append(n)
    return lcm(*ls)
def group_check():
    a=cycles_perm(24,[[1,4],[2,7],[3,17],[5,13],[6,9],[8,15],[10,19],[11,18],[12,21],[14,16],[20,24],[22,23]])
    b=cycles_perm(24,[[1,4,6],[2,21,14],[3,9,15],[5,18,10],[13,17,16],[19,24,23]])
    y=word('ababbababbababbababb',a,b);w=word('bbbbabbbbabbabbbbababbbaabbbbabababababb',a,b)
    z=rprod(rprod(inv(w),y),w)
    e=tuple(range(24));G={e};todo=[e]
    for g in todo:
        for h in [a,z]:
            gh=rprod(g,h)
            if gh not in G:G.add(gh);todo.append(gh)
    assert len(G)==60
    dist=Counter(perm_order(g) for g in G);assert dist=={1:1,2:15,3:20,5:24}
    remaining=set(range(24));orbits=[]
    while remaining:
        x=min(remaining);orb={g[x] for g in G};orbits.append(sorted(i+1 for i in orb));remaining-=orb
    assert [len(o) for o in orbits]==[12,12]
    class_sizes=[];left=set(G)
    while left:
        g=next(iter(left));cl={rprod(rprod(inv(h),g),h) for h in G};class_sizes.append(len(cl));left-=cl
    assert sorted(class_sizes)==[1,12,12,15,20]
    return {'order':len(G),'order_distribution':dict(dist),'class_sizes':sorted(class_sizes),'orbits':orbits}

def main():
    start=time.time();res={'arithmetic':'Z[phi], phi^2=phi+1','tuple_order':['1','D','Dprime','A','B'],'rows':{}}
    for i,(value,n,m,k) in enumerate(ROWS):
        pp,nu,nv=enum_pairs(n,m,k);rr=[(u,v) for u,v in pp if u[1]==u[2] and v[1]==v[2]]
        assert len(pp)==[10,2,11,19,3,26,31][i] and len(rr)==[4,2,5,3,1,6,11][i]
        assert tuple(map(full,PRINTED[i])) in pp
        res['rows'][str(value)]={'complex':len(pp),'rational':len(rr),'U_candidates':nu,'V_candidates':nv,'all_pairs':pp,'rational_pairs':rr}
        def describe(pair):
            u,v=map(full,pair)
            cu,cv=character(u),character(v)
            ext=exterior_ch(cu,2)
            assert ext==exterior_ch(cv,k)
            return {'article_coordinates':pair,'U_character':cu,'V_character':cv,'common_character':ext}
        res['rows'][str(value)]['selected']=describe(PRINTED[i])
        if value in EARLIER:
            assert tuple(map(full,EARLIER[value])) in pp
            res['rows'][str(value)]['earlier_valid_alternative']=describe(EARLIER[value])
        print(value,len(pp),len(rr),nu,nv,flush=True)
    pp,nu,nv=enum_pairs(78,14,6);assert pp==[]
    by=defaultdict(list)
    for v in modules(14):by[exterior_ch(character(v),6)].append(v)
    pv=[(u,v) for u in modules(15) for v in by.get(exterior_ch(character(u),5),[])]
    assert pv==[((15,0,0,0,0),(14,0,0,0,0))]
    res['exclusions']={'78_14_nontrivial':pp,'15_14_all':pv}
    wsp=[0]+list(range(-2,3));wsq=[0]+list(range(-5,6))
    wsU=[wsq[i]+wsq[j] for i,j in combinations_with_replacement(range(12),2)]
    wsV=[wsp[i]+wsp[j] for i,j in combinations(range(6),2)]
    hL=harmonics(weights_ext(wsU,2));hR=harmonics(weights_ext(wsV,5))
    rL=reduce_sum(hL);rR=reduce_sum(hR)
    assert rL==rR==[46,103,96,103,96,53]
    ht=[[scalar(hchar(j),ir) for j in range(6)] for ir in IRR]
    relation=[0,1,2,1,-1,-1]
    assert all(sum(a*b for a,b in zip(row,relation))==0 for row in ht)
    # rank over Q by exact row elimination
    mat=[[Fraction(x) for x in row] for row in ht];rank=0
    for c in range(6):
        p=next((r for r in range(rank,5) if mat[r][c]),None)
        if p is None:continue
        mat[rank],mat[p]=mat[p],mat[rank];d=mat[rank][c];mat[rank]=[x/d for x in mat[rank]]
        for r in range(5):
            if r!=rank:
                d=mat[r][c];mat[r]=[x-d*y for x,y in zip(mat[r],mat[rank])]
        rank+=1
        if rank==5:break
    assert rank==5
    # Kernel integral saturation follows since last coefficient of relation is -1.
    target=tuple([(3003,0),(-21,0),(0,0),(3,0),(3,0)])
    mult=[scalar(target,ir) for ir in IRR];assert mult==[46,156,156,199,245]
    assert sum(x*x for x in mult)==150414
    Pch=character((1,0,0,0,1));Qch=character((1,1,1,0,1))
    modelsP=soh_models(6,Pch);modelsQ=soh_models(12,Qch)
    assert modelsP==[{0:1,2:1}] and modelsQ==[{0:1,5:1}]
    res['harmonics']={'left':hL,'right':hR,'reduced':rL,'restriction_matrix':ht,'rank':rank,'kernel_generator':relation,'P_models':modelsP,'Q_models':modelsQ,'common_irrep_multiplicities':mult,'hom_dimension':150414}
    res['conway_subgroup']=group_check()
    res['duration_seconds']=round(time.time()-start,3)
    (ROOT/'revision_results.json').write_text(json.dumps(res,indent=2))
    print('all checks passed',res['duration_seconds'],flush=True)
if __name__=='__main__':main()
