#!/usr/bin/env python3
"""Independent reconstruction of the saved orbit pairing and Klein class-sum
operator, plus finite checks accompanying the Conway chirality correction."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from math import lcm
import sympy as s
import json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'prior_inputs/Full_Conway_Klein_Comparison_Package';OUT=ROOT/'certificates'
def perm(cycles):
 p=list(range(24))
 for cc in cycles:
  c=[a-1 for a in cc]
  for a,b in zip(c,c[1:]+c[:1]):p[a]=b
 return tuple(p)
def mult(a,b):return tuple(b[a[i]] for i in range(24))
def inverse(a):return tuple(a.index(i) for i in range(24))
def order(p):
 left=set(range(24));nn=[]
 while left:
  i=min(left);n=0
  while i in left:left.remove(i);n+=1;i=p[i]
  nn.append(n)
 return lcm(*nn)
a=perm([(1,4),(2,7),(3,17),(5,13),(6,9),(8,15),(10,19),(11,18),(12,21),(14,16),(20,24),(22,23)])
b=perm([(1,4,6),(2,21,14),(3,9,15),(5,18,10),(13,17,16),(19,24,23)])
e=tuple(range(24))
def word(w):
 r=e
 for c in w:r=mult(r,a if c=='a' else b)
 return r
y=word('ababbababbababbababb');w=word('bbbbabbbbabbabbbbababbbaabbbbabababababb');z=mult(mult(inverse(w),y),w)
G={e};todo=[e]
while todo:
 h=todo.pop()
 for g in (a,z):
  k=mult(h,g)
  if k not in G:G.add(k);todo.append(k)
cl=[];left=set(G)
while left:
 g=min(left);cc={mult(mult(inverse(h),g),h) for h in G};left-=cc;cl.append(sorted(cc))
cl.sort(key=lambda c:(order(c[0]),c[0]))
d=json.loads((B/'basis_data.json').read_text());X,Y=d['orbits'];pair=dict(zip(X,Y))
assert all(pair[g[x]]==g[pair[x]] for g in G for x in X)
ss=[]
for cc in cl:
 M=s.zeros(12)
 for g in cc:
  for j,x in enumerate(X):M[X.index(g[x]),j]+=1
 ss.append(M)
M0=2*ss[0]+2*ss[1]-ss[2]-(ss[3]+ss[4])/2
M1=s.Rational(5,2)*(ss[3]-ss[4])
assert M0==s.Matrix(d['A0_times4']) and 2*M1==s.Matrix(d['A1_times8'])
def prod(A,B):
 a,b=A;c,d=B
 return a*c+5*b*d,a*d+b*c
AA=(M0/4,M1/4);eye=(s.eye(12),s.zeros(12));traces={}
for l in [0,3,10,-15]:
 P=eye;den=1
 for k in [0,3,10,-15]:
  if k!=l:P=prod(P,(AA[0]-k*s.eye(12),AA[1]));den*=l-k
 assert s.trace(P[1])==0
 traces[str(l)]=int(s.trace(P[0])/den)
assert traces=={'0':1,'3':5,'10':3,'-15':3}
# Fixed vectors and a reflection needed to compare the two Conway conventions.
fixed=s.zeros(24,1)
for i in X:fixed[i]=1
for g in G:assert all(fixed[g[i]]==fixed[i] for i in range(24))
R=s.eye(24)-s.Rational(1,6)*fixed*fixed.T
assert R*R==s.eye(24) and R.det()==-1
for gen in [a,z]:
 P=s.zeros(24)
 for i,j in enumerate(gen):P[j,i]=1
 assert R*P==P*R
# Two chiral restrictions have equal characters because of a fixed vector;
# the Pin intertwiner is invertible and exchanges the chirality eigenspaces.
tau=json.loads((B/'tau_exact_0.json').read_text())['nonzero']
assert all(n.bit_count()%2==0 for n,r,i in tau)
ordered=[i for pair_xy in zip(X,Y) for i in pair_xy]
orientation=(-1)**sum(ordered[i]>ordered[j] for i in range(24) for j in range(i+1,24))
assert orientation==1
cert={'marked_pairing_equivariant':True,'paired_Majorana_orientation':orientation,'class_sum_Klein_operator_rebuilt':True,'Klein_eigenspace_dimensions':traces,'fixed_vector_support':X,'reflection_det':-1,'reflection_commutes_with_A5':True,'sparse_supercurrent_Clifford_chirality':'even','convention_correction':'V^f and V^s are not the same Co0 module; an A5-fixed-vector reflection provides an A5-equivariant transport after choosing its extension phase.'}
(OUT/'marked_basis.json').write_text(json.dumps(cert,indent=2));print(json.dumps(cert,indent=2))
