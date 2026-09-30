#!/usr/bin/env python3
"""Fresh MI-N1 computation. No archived code or character answers are imported.
Construct A5 as even permutations; construct its 5,6,12,10 point actions;
compare exterior characters by Newton recursion AND exact characteristic polynomials.
Also enumerate signed basis orbits and all sl2 partition obstructions.
"""
from itertools import permutations,combinations
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path
import json
OUT=Path(__file__).resolve().parents[1]/'certificates'

def mul(p,q): return tuple(p[q[i]] for i in range(len(p)))
def inverse(p): return tuple(p.index(i) for i in range(len(p)))
def power(p,k):
 r=tuple(range(len(p)))
 for _ in range(k):r=mul(r,p)
 return r
def parity(p):return sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2
def cycles(p):
 left=set(range(len(p)));res=[]
 while left:
  x=min(left);k=0
  while x in left:left.remove(x); k+=1;x=p[x]
  res.append(k)
 return sorted(res)
def order(p):
 from math import lcm
 return lcm(*cycles(p))
G=sorted(p for p in permutations(range(5)) if not parity(p));eid=tuple(range(5))
assert len(G)==60
h=(1,2,3,4,0); H={power(h,i) for i in range(5)}
N={g for g in G if {mul(mul(g,k),inverse(g)) for k in H}==H}
assert len(N)==10

def coset_action(H):
 left=set(G);cos=[]
 while left:
  g=min(left);s=frozenset(mul(g,k) for k in H);cos.append(s);left-=s
 ix={s:i for i,s in enumerate(cos)}
 actions={g:tuple(ix[frozenset(mul(g,k) for k in c)] for c in cos) for g in G}
 return cos,actions
_,a6=coset_action(N);_,a12=coset_action(H)
edges=list(combinations(range(5),2));ei={e:i for i,e in enumerate(edges)}
a10={g:tuple(ei[tuple(sorted((g[a],g[b])))] for a,b in edges) for g in G}

def pmul(a,b,cap=None):
 n=len(a)+len(b)-1
 if cap is not None:n=min(n,cap+1)
 out=[0]*n
 for i,x in enumerate(a):
  for j,y in enumerate(b[:max(0,n-i)]):out[i+j]+=x*y
 return out

def pdiv(a,b):
 assert b[0]==1
 q=[0]*(len(a)-len(b)+1)
 for k in range(len(q)):
  q[k]=a[k]-sum(b[i]*q[k-i] for i in range(1,min(k,len(b)-1)+1))
 assert pmul(q,b)==a
 return q

def detperm(p):
 a=[1]
 for l in cycles(p):
  b=[0]*(l+1);b[0]=1;b[l]=(-1)**(l+1);a=pmul(a,b)
 return a

def fixed(p):return sum(i==j for i,j in enumerate(p))
def basics(g):
 t=[1,1]
 return [t,pdiv(detperm(g),t),pdiv(detperm(a6[g]),t),pdiv(detperm(a12[g]),detperm(a6[g]))]
Bpol={g:basics(g) for g in G}
Bchar={g:[p[1] for p in Bpol[g]] for g in G}

def char_direct(t,g):return sum(a*b for a,b in zip(t,Bchar[g]))
def exterior_newton(t,g,k):
 e=[1]
 for j in range(1,k+1):
  a=sum((-1)**(i-1)*char_direct(t,power(g,i))*e[j-i] for i in range(1,j+1))
  assert a%j==0;e.append(a//j)
 return e[k]
def exterior_product(t,g,k):
 out=[1]
 for n,p in zip(t,Bpol[g]):
  for _ in range(n):out=pmul(out,p,k)
 return out[k] if k<len(out) else 0
rows=[(120,16,10,3,(0,0,2,1),(1,1,1,0)),(210,21,10,4,(3,2,2,0),(1,1,1,0)),(1540,56,22,3,(4,5,4,2),(3,2,1,1)),(3003,78,15,5,(4,4,8,3),(0,1,1,1)),(7140,120,36,3,(3,9,9,6),(1,1,5,1)),(11628,153,19,5,(6,9,15,6),(1,0,0,3)),(24310,221,17,8,(6,16,17,11),(3,1,2,0))]
remaining=set(G);classes=[]
while remaining:
 g=min(remaining);cls={mul(mul(x,g),inverse(x)) for x in G};remaining-=cls
 classes.append((order(g),g,cls))
classes.sort(key=lambda c:(c[0],c[1])); reps=[c[1] for c in classes]
assert [(o,len(c)) for o,g,c in classes]==[(1,1),(2,15),(3,20),(5,12),(5,12)]
cert={'construction':'even permutations, independently built cosets and edge action','class_sizes':[len(c) for o,g,c in classes],'module_characters':{name:[Bchar[g][i] for g in reps] for i,name in enumerate(['1','A','B','C'])},'rows':[]}
for value,m,n,k,u,v in rows:
 assert char_direct(u,eid)==m and char_direct(v,eid)==n and comb(m,2)==comb(n,k)==value
 allvals={}
 for g in G:
  a=exterior_newton(u,g,2); b=exterior_newton(v,g,k)
  assert a==b==exterior_product(u,g,2)==exterior_product(v,g,k)
  allvals[g]=a
 chi=[allvals[g] for g in reps]
 cert['rows'].append(dict(value=value,m=m,n=n,k=k,U=u,V=v,character=chi,group_elements_checked=60))
h3=next(g for g in G if order(g)==3)
_,a20=coset_action({power(h3,i) for i in range(3)})
for g in G:
 q=fixed(a12[g]);q2=fixed(a12[power(g,2)]); p=fixed(a6[g]);p2=fixed(a6[power(g,2)])
 assert (q*q+q2)//2==char_direct(rows[3][4],g)
 assert (p*p-p2)//2==char_direct(rows[3][5],g)
 w=q+Bchar[g][2]; w2=fixed(a12[power(g,2)])+Bchar[power(g,2)][2]
 assert (w*w+w2)//2==char_direct(rows[5][4],g)
 assert 1+2*fixed(a10[g])==char_direct(rows[1][4],g)
 assert fixed(a10[g])==char_direct((1,1,1,0),g)
 ca=Bchar[g][1]; ca2=Bchar[power(g,2)][1]
 assert (ca*ca+ca2)//2==fixed(a10[g])
 assert (ca*ca-ca2)//2==Bchar[g][3]
 assert ca*(Bchar[g][3]-2)==char_direct(rows[0][4],g)
 # New inputs are checked from their constructions, independently of row data.
 assert fixed(a6[g])+fixed(a10[g])+2*fixed(a20[g])==char_direct(rows[2][4],g)
 assert 2*fixed(g)+fixed(a12[g])==char_direct(rows[2][5],g)
 assert 1+detperm(a12[g])[3]==char_direct(rows[6][4],g)
 assert fixed(g)+2*fixed(a6[g])==char_direct(rows[6][5],g)
 wc=1+3*Bchar[g][2];wc2=1+3*Bchar[power(g,2)][2]
 assert (wc*wc-wc2)//2==char_direct(rows[4][4],g)
 assert wc+ca*Bchar[g][2]==char_direct(rows[4][5],g)
cert['new_constructions']={
 '1540':{'U':'P+E+2R20','V':'2R5+Q'},
 '7140':{'W':'1+3B','T':'A tensor B','U':'exterior square W','V':'W+T'},
 '24310':{'U':'1+exterior cube Q','V':'R5+2P'},
 'group_elements_checked':60}
cert['earlier_valid_alternatives']=[]
for value,k,u,v in [(1540,3,(0,2,6,3),(4,2,2,0)),(24310,8,(3,13,20,11),(2,0,3,0))]:
 for g in G:
  assert exterior_newton(u,g,2)==exterior_newton(v,g,k)==exterior_product(u,g,2)==exterior_product(v,g,k)
 cert['earlier_valid_alternatives'].append({'value':value,'U':u,'V':v,'character':[exterior_newton(u,g,2) for g in reps],'group_elements_checked':60})
cert['geometric_identifications']='all 60 group elements checked'
q_pairs=[(i,j) for i in range(12) for j in range(i,12)];qi={e:i for i,e in enumerate(q_pairs)}
q_sym={g:tuple(qi[tuple(sorted((a12[g][i],a12[g][j])))] for i,j in q_pairs) for g in G}
p_pairs=list(combinations(range(6),2));pi={e:i for i,e in enumerate(p_pairs)}
p_wedge={};p_sign={}
for g in G:
 pp=[];ss=[]
 for i,j in p_pairs:
  a,b=a6[g][i],a6[g][j];pp.append(pi[tuple(sorted((a,b)))]);ss.append(1 if a<b else -1)
 p_wedge[g]=pp;p_sign[g]=ss

def signed_orbits(labels,kind):
 unseen=set(labels);counts=Counter();sizes=Counter()
 while unseen:
  x=min(unseen);images=set();stab=[]
  for g in G:
   if kind=='L':
    a,b=(q_sym[g][i] for i in x);y=tuple(sorted((a,b)));sgn=1 if a<b else -1
   else:
    ys=[p_wedge[g][i] for i in x];sgn=(-1)**parity(ys)
    for i in x:sgn*=p_sign[g][i]
    y=tuple(sorted(ys))
   images.add(y)
   if y==x:stab.append((g,sgn))
  unseen-=images
  typ={1:'1',2:'C2',4:'V4',5:'C5',10:'D5'}[len(stab)]
  sig='trivial' if all(s==1 for g,s in stab) else 'nontrivial'
  counts[(typ,sig)]+=1;sizes[len(images)]+=1
 return {'types':{'/'.join(k):v for k,v in sorted(counts.items())},'unsigned_sizes':dict(sorted(sizes.items()))}
cert['signed_orbits']={'L':signed_orbits(combinations(range(78),2),'L'),'R':signed_orbits(combinations(range(15),5),'R')}
OUT.mkdir(exist_ok=True)
(OUT/'finite_core.json').write_text(json.dumps(cert,indent=2))
print(json.dumps(cert,indent=2),flush=True)

def parts(n,maxpart=None):
 if n==0:yield ();return
 for d in range(min(n,maxpart or n),0,-1):
  for rest in parts(n-d,d):yield (d,)+rest

def moments_exterior(weights,k,rmax=6):
 dp=[[0]*(rmax+1) for _ in range(k+1)];dp[0][0]=1
 used=0
 for w in weights:
  for l in range(min(k,used+1),0,-1):
   for r in range(rmax+1):
    dp[l][r]+=sum(comb(r,s)*w**(r-s)*dp[l-1][s] for s in range(r+1))
  used+=1
 return dp[k]
obs=[]
triples=[(21,10,4),(56,22,3),(78,14,6),(78,15,5),(120,36,3),(153,19,5),(221,17,8)]
for m,n,k in triples:
 total=0;pass2=[];pass4=[];pass6=[]
 for partition in parts(n):
  if max(partition)==1:continue
  total+=1
  I=sum(comb(d+1,3) for d in partition)
  a2=Fraction(2*comb(n-2,k-1)*I,m-2)
  if a2.denominator!=1 or a2<=0 or a2%2:continue
  pass2.append(partition)
  ww=[d-1-2*i for d in partition for i in range(d)]
  M=moments_exterior(ww,k)
  assert M[0]==comb(n,k) and M[2]==(m-2)*a2
  assert M[1]==M[3]==M[5]==0
  a4=Fraction(M[4]-3*a2*a2,m-8)
  if a4.denominator!=1 or a4<=0 or (a4-a2)%24:continue
  a6=Fraction(M[6]-15*a2*a4,m-32)
  pass4.append({'partition':partition,'S2':str(a2),'S4':str(a4),'S6':str(a6)})
  if a6.denominator!=1 or a6<=0 or (a6-a2)%120:continue
  pass6.append(partition)
 obs.append(dict(m=m,n=n,k=k,nontrivial_partitions=total,pass2=len(pass2),pass4=len(pass4),pass6=len(pass6),fourth_survivors=pass4))
 assert not pass6
 print('obstruction',m,n,k,total,len(pass2),len(pass4),len(pass6),flush=True)
(OUT/'semisimple_obstruction.json').write_text(json.dumps(obs,indent=2))
