#!/usr/bin/env python3
"""Fresh exact checks of congruence data, affine branches, M24 subgroup and
low-energy protected sector. Uses no archived script or expected character series.
The character identity proofs additionally use the external conformal embedding.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
from math import comb,lcm
from pathlib import Path
import json
OUT=Path(__file__).resolve().parents[1]/'certificates'
# SL(2,F5), projectivization and cusp action on nonzero vectors modulo sign.
SL=[(a,b,c,d) for a,b,c,d in product(range(5),repeat=4) if (a*d-b*c)%5==1]
def canmat(m):return min(m,tuple(-x%5 for x in m))
PSL=sorted({canmat(m) for m in SL});assert len(SL)==120 and len(PSL)==60
def canvec(v):return min(v,tuple(-x%5 for x in v))
C=sorted({canvec(v) for v in product(range(5),repeat=2) if v!=(0,0)})
assert len(C)==12
cusp_action=[]
for a,b,c,d in PSL:
 act=[C.index(canvec(((a*x+b*y)%5,(c*x+d*y)%5))) for x,y in C]
 cusp_action.append(act)
assert len({a[0] for a in cusp_action})==12
assert sum(p[0]==0 for p in cusp_action)==5
cert={'Gamma5':{'SL2_order':120,'PSL2_order':60,'cusps_mod_sign':12,'cusp_stabilizer':5,'genus_formula':str(F(1)+F(60,12)-F(12,2))}}
# Direct coefficient recurrences.
N=24

def conv(a,b):
 c=[0]*(N+1)
 for i,x in enumerate(a[:N+1]):
  for j,y in enumerate(b[:N+1-i]):c[i+j]+=x*y
 return c
# Three independent bosonic oscillator species, built with the coin recurrence.
denom=[1]+[0]*N
for m in range(1,N+1):
 for species in range(3):
  for n in range(m,N+1):denom[n]+=denom[n-m]
chars=[]
for ell in range(7):
 numerator=[0]*(N+1)
 for k in range(-10,11):
  e=8*k*k+(ell+1)*k
  if 0<=e<=N:numerator[e]+=16*k+ell+1
 chars.append(conv(denom,numerator))
h=[F(ell*(ell+2),32) for ell in range(7)];hb=[F(0),F(1,2),F(9,16)]
branches=[[(a,b) for a,b in product(range(7),repeat=2) if h[a]+h[b]>=h0 and (h[a]+h[b]-h0).denominator==1] for h0 in hb]
# Independent fermion Fock occupancy DP, with t=q^(1/2).
ns=[1]+[0]*(2*N+1)
for mode in range(1,2*N+2,2):
 for sp in range(9):
  for energy in range(2*N+1,mode-1,-1):ns[energy]+=ns[energy-mode]
r=[16]+[0]*N
for mode in range(1,N+1):
 for sp in range(9):
  for e in range(N,mode-1,-1):r[e]+=r[e-mode]
ferm=[ns[0::2],ns[1::2],r]
for idx,pairs in enumerate(branches):
 out=[0]*(N+1)
 for a,b in pairs:
  k=int(h[a]+h[b]-hb[idx]);ser=conv(chars[a],chars[b])
  for j in range(N+1-k):out[j+k]+=ser[j]
 assert out==ferm[idx][:N+1]
cert['affine']={'level':6,'grade_checked':N,'branches':branches,'spin_first_five':r[:5]}
# Exact S-matrix arithmetic in Q[z]/(z^8+1).
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,c):return tuple(c*x for x in a)
def times(a,b):
 out=[F(0)]*8
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[(i+j)%8]+=x*y*(-1 if i+j>=8 else 1)
 return tuple(out)
def root(n):
 n%=16;a=[F(0)]*8;a[n%8]=F(1 if n<8 else -1);return tuple(a)
z0=root(0);zero=scale(z0,0)
def summ(xs):
 out=zero
 for x in xs:out=add(out,x)
 return out
S6=[[scale(times(root(-4),add(root((a+1)*(b+1)),scale(root(-(a+1)*(b+1)),-1))),F(1,4)) for b in range(7)] for a in range(7)]
sq2=add(root(2),root(-2));S9=[[scale(z0,F(1,2)),scale(z0,F(1,2)),scale(sq2,F(1,2))],[scale(z0,F(1,2)),scale(z0,F(1,2)),scale(sq2,F(-1,2))],[scale(sq2,F(1,2)),scale(sq2,F(-1,2)),zero]]
for i in range(3):
 for a,b in product(range(7),repeat=2):
  assert summ(times(S6[x][a],S6[y][b]) for x,y in branches[i])==summ(S9[i][j] for j in range(3) if (a,b) in branches[j])
cert['affine']['modular_S_entries_checked']=147
cert['affine']['T_checked']='exact conformal-weight congruences; central charges equal'
# ATLAS generators, transcribed from primary GAP files verified in this stage.
def perm(cycles):
 p=list(range(24))
 for cyc in cycles:
  c=[a-1 for a in cyc]
  for a,b in zip(c,c[1:]+c[:1]):p[a]=b
 return tuple(p)
def compose(a,b):return tuple(b[a[i]] for i in range(len(a))) # right action
A=perm([(1,4),(2,7),(3,17),(5,13),(6,9),(8,15),(10,19),(11,18),(12,21),(14,16),(20,24),(22,23)])
B=perm([(1,4,6),(2,21,14),(3,9,15),(5,18,10),(13,17,16),(19,24,23)])
identity=tuple(range(24))
def word(w):
 g=identity
 for l in w:g=compose(g,A if l=='a' else B)
 return g
y=word('ababbababbababbababb');w=word('bbbbabbbbabbabbbbababbbaabbbbabababababb')
wi=tuple(w.index(i) for i in range(24));Z=compose(compose(wi,y),w)
G={identity};todo=[identity]
while todo:
 g=todo.pop()
 for a in [A,Z]:
  q=compose(g,a)
  if q not in G:G.add(q);todo.append(q)
assert len(G)==60
orbits=[];left=set(range(24))
while left:
 x=min(left);orbit={g[x] for g in G};orbits.append(sorted(orbit));left-=orbit
assert list(map(len,orbits))==[12,12]
def cycletype(p):
 todo=set(range(24));lens=[]
 while todo:
  x=min(todo);k=0
  while x in todo:todo.remove(x);k+=1;x=p[x]
  lens.append(k)
 return tuple(sorted(lens))
ct=Counter(cycletype(g) for g in G)
for g in G:
 o=lcm(*cycletype(g))
 assert all(sum(g[x]==x for x in orbit)=={1:12,2:0,3:0,5:2}[o] for orbit in orbits)
cert['M24']={'group_order':60,'orbits_0_based':orbits,'cycle_inventory':{str(k):v for k,v in sorted(ct.items())},'generators_source':'ATLAS M24G1-p24B0.g1 and .g2'}
# Topological zero-energy ground states in D12 plus spinor lattice.
# h_top=(p-e1)^2/2 -1/2+n. At 0 it forces n=0 and p-e1 = +-e_i.
p0=[]
for i in range(12):
 for sg in [-1,1]:
  p=[0]*12;p[0]=1;p[i]+=sg
  assert sum(p)%2==0 and F(sum(x*x for x in p),2)-p[0]==0
  p0.append(p)
charges=Counter(2*p[0] for p in p0);assert charges=={0:1,2:22,4:1}
cert['topological_cohomology_ground_charges']=dict(charges)
# Complete first positive-grade complex by neutral-mode exterior charges.
chargesV=[3]*5+[10]*3+[-15]*3+[-3]*5+[-10]*3+[15]*3
from itertools import combinations
c0=Counter(sum(xs) for xs in combinations(chargesV,2));c0[0]+=1
c2=Counter(sum(xs) for xs in combinations(chargesV,3));c2.update(Counter({k:2*v for k,v in Counter(chargesV).items()}))
spin=[Counter(),Counter()]
for a,b,c in product(range(6),range(4),range(4)):
 spin[(a+b+c)%2][3*a+10*b-15*c]+=comb(5,a)*comb(3,b)*comb(3,c)
complexes=[c0,spin[0],c2,spin[1],c0]
assert [sum(x.values()) for x in complexes]==[232,1024,1584,1024,232]
block=[co[-5+2*r] for r,co in enumerate(complexes)];assert block==[9,5,0,0,0]
cert['native_small_sector']={'total_dimensions':[sum(x.values()) for x in complexes],'C_equals_minus5_complete_dimensions':block,'spin_K2_occupation':[a for a in product(range(6),range(4),range(4)) if sum(a)%2==0 and 3*a[0]+10*a[1]-15*a[2]==2]}
cert['levels_and_central_charges']={'SU2_tensor_level':6,'principal_Spin9_level':sum(i*i for i in range(-4,5)),'harmonic_Conway_level':2*sum(i*i for i in range(-5,6)),'Spin9_coset_c':str(F(9,2)-F(180,62)),'Conway_coset_c':str(12-F(660,222))}
(OUT/'modular_conformal.json').write_text(json.dumps(cert,indent=2));print(json.dumps(cert,indent=2))
