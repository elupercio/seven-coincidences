#!/usr/bin/env python3
"""Fresh exact weight and arithmetic verification, independent of archived scripts."""
from collections import Counter
from itertools import combinations,combinations_with_replacement,product
from fractions import Fraction as R
from pathlib import Path
from math import comb
import sympy as s
import json
OUT=Path(__file__).resolve().parents[1]/'certificates'
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def weights_ext(ww,k):
 dim=len(ww[0]);z=(0,)*dim;dp=[Counter() for _ in range(k+1)];dp[0][z]=1
 for w in ww:
  for j in range(k,0,-1):
   for v,m in list(dp[j-1].items()):dp[j][add(v,w)]+=m
 return +dp[k]
def sym2(ww):return Counter(add(ww[i],ww[j]) for i in range(len(ww)) for j in range(i,len(ww)))
def expanded(c):return [w for w,m in c.items() for _ in range(m)]
def vweights(n):
 return [tuple(sg*2 if j==i else 0 for j in range(n)) for i in range(n) for sg in [-1,1]]
def spweights(n,parity=0):return [ss for ss in product([-1,1],repeat=n) if sum(x<0 for x in ss)%2==parity]
spin=[]
for n in [5,6]:
 rhs=weights_ext(vweights(n),n-2)
 if n==6:rhs[(0,)*n]+=1
 for chirality in [0,1]:
  lhs=weights_ext(spweights(n,chirality),2);assert lhs==rhs
  spin.append({'rank':n,'chirality':chirality,'dimension':sum(lhs.values()),'distinct_weights':len(lhs)})
Aw=[(1,0),(-1,0),(0,1),(0,-1)]
planes=[(2,0),(0,2),(1,1),(1,-1),(0,0)]
restring=[]
for parity in [0,1]:
 out=Counter()
 for signs in spweights(5,parity):
  v=tuple(sum(signs[i]*planes[i][j] for i in range(5)) for j in range(2))
  assert all(x%2==0 for x in v);out[tuple(x//2 for x in v)]+=1
 tensor=Counter(add(a,b) for a in Aw for b in expanded(weights_ext(Aw,2)))
 tensor.subtract(Counter({w:2 for w in Aw}));tensor=+tensor
 assert out==tensor
 assert weights_ext(expanded(out),2)==weights_ext(expanded(sym2(Aw)),3)
 restring.append({'chirality':parity,'tensor_dimension':sum(out.values())})
principal=Counter((sum(a*b for a,b in zip(signs,[1,2,3,4]))//2,) for signs in product([-1,1],repeat=4))
assert principal==Counter((x,) for j in [2,5] for x in range(-j,j+1))
cert={'full_spin_identities':spin,'SO4_tensor_branching':restring,'principal_spinor_weights':{str(k[0]):v for k,v in sorted(principal.items())},'lift_winding_parities':{'SO4_generator':(2+1+1)%2,'SO4_minus_identity_lift':(1+1+1)%2,'principal_SO3_generator':sum(range(1,5))%2}}
# Elliptic curve: group law for y^2+y=x^3-9x+20.
def neg(P):return None if P is None else (P[0],-1-P[1])
def plus(P,Q):
 if P is None:return Q
 if Q is None:return P
 x,y=P;xx,yy=Q
 if x==xx and yy==-1-y:return None
 lam=R(3*x*x-9,2*y+1) if P==Q else R(yy-y,xx-x)
 xx3=lam*lam-x-xx
 return xx3,-1-y+lam*(x-xx3)
def times(k,P):
 if k<0:return times(-k,neg(P))
 out=None
 for _ in range(k):out=plus(out,P)
 return out
P=(R(0),R(4));Q=(R(3),R(4))
pts=[plus(times(2,Q),neg(P)),times(3,P),neg(plus(P,times(3,Q)))]
checks=[]
for x,y in pts:
 assert y*y+y==x**3-9*x+20
 n=x/3+1;m=(y+5)/9;assert n.denominator==m.denominator==1
 n=int(n);m=int(m);assert comb(m,2)==comb(n,3)
 checks.append({'point':[str(x),str(y)],'n':n,'m':m,'value':comb(m,2)})
cert['elliptic_points']=checks
m,n,t=s.symbols('m n t')
expr=(9*m-5)**2+(9*m-5)-(3*(n-1))**3+9*(3*(n-1))-20
assert s.expand(expr-162*(m*(m-1)/2-n*(n-1)*(n-2)/6))==0
families=[]
for shift in [0,4]:
 nn=15*t*t+shift
 mm=(t*(nn**2-5*nn+5)+1)/2 if shift==0 else (t*(nn**2-3*nn+1)+1)/2
 rr=(t+1)/2
 c5=s.prod(nn-i for i in range(5))/120
 defect=s.factor(mm*(mm-1)/2-c5-rr*(rr-1)/2)
 assert defect==0
 families.append({'shift':shift,'symbolic_defect':str(defect),'n_at_1':int(nn.subs(t,1)),'m_at_1':int(mm.subs(t,1))})
cert['fifth_order_families']=families
F=[0,1];L=[2,1]
for i in range(2,80):F.append(F[-1]+F[-2]);L.append(L[-1]+L[-2])
for i in range(80):assert L[i]**2-5*F[i]**2==4*(-1)**i
ls=[]
for i in range(1,7):
 nn=F[2*i+2]*F[2*i+3];kk=F[2*i]*F[2*i+3]
 assert nn*(kk+1)==(nn-kk)*(nn-kk-1)
 assert L[4*i+3]-4==5*kk
 ls.append([nn,kk])
cert['Lind_Singmaster_first_six']=ls
(OUT/'spin_arithmetic.json').write_text(json.dumps(cert,indent=2));print(json.dumps(cert,indent=2))
