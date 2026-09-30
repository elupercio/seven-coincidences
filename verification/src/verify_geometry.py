#!/usr/bin/env python3
"""Fresh exact plethysms, transvectants, and an external generator test of
archived residue operators. All group-equivariance tests use exact Q(sqrt(5))
arithmetic, not sampling or modular equality tests.
"""
from collections import Counter
from itertools import combinations_with_replacement
from math import comb,factorial,lcm
from pathlib import Path
import sympy as s
import json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'certificates'

def ext(C,k):
 dp=[Counter() for _ in range(k+1)];dp[0][0]=1
 for w,m in C.items():
  for _ in range(m):
   for j in range(k,0,-1):
    for v,a in list(dp[j-1].items()):dp[j][v+w]+=a
 return +dp[k]
def sym2(C):
 ww=[w for w,m in C.items() for _ in range(m)]
 return Counter(ww[i]+ww[j] for i in range(len(ww)) for j in range(i,len(ww)))
def harms(C):
 assert all(C[w]==C[-w] for w in C)
 ans={j:C[j]-C[j+1] for j in range(max(C)+1)}
 assert all(a>=0 for a in ans.values())
 return {j:a for j,a in ans.items() if a}
def H(j):return Counter(range(-j,j+1))

def reduction(j):
 if j<6:return Counter({j:1})
 r=j%5;return reduction(j-6)+Counter([r,5-r])
def redall(hs):
 out=Counter()
 for j,m in hs.items():
  for a,b in reduction(j).items():out[a]+=b*m
 return out
L=harms(ext(sym2(H(0)+H(5)),2));R=harms(ext(ext(H(0)+H(2),2),5))
assert redall(L)==redall(R)
L210=harms(ext(Counter({0:3})+Counter({w:2*m for w,m in H(4).items()}),2))
R210=harms(ext(H(0)+H(4),4));assert redall(L210)==redall(R210)
cert={'L_3003':L,'R_3003':R,'common_reduction_3003':dict(redall(L)),'common_reduction_210':dict(redall(L210))}
# Polynomial arithmetic and exact residue sequence ranks.
z=s.symbols('z');f=z**11+11*z**6-z
assert s.gcd(f,s.diff(f,z))==1
# Formal CRT representative [polynomial degree<=10, value at infinity].
def remainder(expr):return s.Poly(s.rem(expr,f,z),z)
def ev(j):
 return s.Matrix(12,2*j+1,lambda r,c: remainder(z**c).nth(r) if r<11 else int(c==2*j))
D=s.zeros(12)
for col in range(11):
 p=remainder(s.diff(f,z)*z**col)
 for row in range(11):D[row,col]=p.nth(row)
D[11,11]=-1
cert['delta_determinant']=str(D.det());assert D.det()==5**25
low=[]
for r in range(6):
 E=ev(r); beta=s.zeros(11-2*r,12)
 for b in range(11-2*r):
  for c in range(11):beta[b,c]=remainder(z**(b+c)).nth(10)
  beta[b,11]=-int(b==10-2*r)
 assert beta*E==s.zeros(11-2*r,2*r+1)
 assert E.rank()==2*r+1 and beta.rank()==11-2*r
 low.append({'r':r,'restriction_rank':E.rank(),'residue_rank':beta.rank()})
cert['low_divisor_sequences']=low
for j in range(6,20):
 E=ev(j); Fmul=s.zeros(2*j+1,2*(j-6)+1)
 for col in range(Fmul.cols):
  for a,b in [(11,1),(6,11),(1,-1)]:Fmul[a+col,col]=b
 assert E*Fmul==s.zeros(12,Fmul.cols)
 assert E.rank()==12 and Fmul.rank()==2*j-11
# Integral pairs p+q sqrt(5).
def pairmul(a,b):return a[0]*b[0]+5*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def pairpow(a,n):
 r=(1,0)
 for _ in range(n):r=pairmul(r,a)
 return r
_cache={}
def involution_numerator(j):
 # Sym^(2j) of i/sqrt(1+phi^2) [[1,phi],[phi,-1]], ascending X powers.
 # Common denominator is 40^j. Entries are exact integer pairs.
 if j in _cache:return _cache[j]
 n=2*j;A=[[0]*(n+1) for _ in range(n+1)];B=[[0]*(n+1) for _ in range(n+1)]
 u=pairpow((5,-1),j)
 pp=[pairpow((1,1),a) for a in range(n+1)]
 for c in range(n+1):
  for r in range(n+1):
   v=(0,0)
   for a in range(max(0,r-(n-c)),min(c,r)+1):
    b=r-a;e=c-a+b
    coef=(-1)**(j+n-c-b)*comb(c,a)*comb(n-c,b)*2**(n-e)
    t=pp[e];v=(v[0]+coef*t[0],v[1]+coef*t[1])
   A[r][c],B[r][c]=pairmul(u,v)
 _cache[j]=(A,B);return A,B

def matprod(A,B):
 nr=len(A);nc=len(B[0]);nk=len(B);out=[[0]*nc for _ in range(nr)]
 for i in range(nr):
  for k in range(nk):
   a=A[i][k]
   if a:
    for j,b in enumerate(B[k]):
     if b:out[i][j]+=a*b
 return out
# The generator fixes the homogeneous groundform, exact over Q(sqrt5).
A6,B6=involution_numerator(6);Fv=[0]*13
for i,c in [(11,1),(6,11),(1,-1)]:Fv[i]=c
assert [sum(a*b for a,b in zip(row,Fv)) for row in A6]==[40**6*c for c in Fv]
assert all(sum(a*b for a,b in zip(row,Fv))==0 for row in B6)
cert['Klein_F_generator_invariance']='exact over Q(sqrt5); diagonal 5-fold invariance follows from monomial weights'
# Audit the archived maps against both independent generators.
archive=ROOT/'prior_inputs/Klein_Residue_Intertwiner/residue_reduction_operators.json'
data=json.loads(archive.read_text());checks=[]
for jj,rows in data.items():
 j=int(jj);r=j%5;ks=[j-6,r,5-r]
 T=s.Matrix([[s.Rational(x) for x in row] for row in rows])
 denom=lcm(*(int(x.q) for x in T));Ti=[[int(x*denom) for x in row] for row in T.tolist()]
 # diagonal root-of-unity generator
 rowlabels=[(k,i) for k in ks for i in range(2*k+1)]
 for row,(k,i) in enumerate(rowlabels):
  for c,t in enumerate(Ti[row]):
   assert not t or (c-j-(i-k))%5==0
 for q in [0,1]:
  SJ=involution_numerator(j)[q];lhs=matprod(Ti,SJ);rhs=[];offset=0
  for k in ks:
   nk=2*k+1;K=involution_numerator(k)[q];p=matprod(K,Ti[offset:offset+nk]);sc=40**(j-k)
   rhs.extend([[sc*a for a in row] for row in p]);offset+=nk
  assert lhs==rhs,('equivariance failed',j,q)
 # Exact rational rank, independent of old modular checks.
 rank=T.rank();assert rank==2*j+1
 checks.append({'j':j,'rank':rank,'two_generator_equivariance':'exact','field':'Q(sqrt5), fifth-root weight congruences'})
 print('residue',j,'exact equivariance and rank',rank,flush=True)
cert['archived_residue_maps_independent_checks']=checks
# Linear transvectants built directly from monomial derivatives.
def ff(a,b):return factorial(a)//factorial(a-b) if 0<=b<=a else 0
def trans_fixed(d,r):
 out=s.zeros(12+d-2*r+1,d+1)
 for c in range(d+1):
  for i,coef in [(11,1),(6,11),(1,-1)]:
   row=i+c-r
   if not 0<=row<out.rows:continue
   out[row,c]+=sum(coef*(-1)**a*comb(r,a)*ff(i,r-a)*ff(12-i,a)*ff(c,a)*ff(d-c,r-a) for a in range(r+1))
 return out
K=trans_fixed(10,6)/399168000
metric=s.diag(*[s.Rational(1,comb(10,i)) for i in range(11)])
assert metric*K==K.T*metric
cp=s.factor(K.charpoly().as_expr());cert['Klein_operator_charpoly']=str(cp)
assert cp==(s.Symbol('lambda')-10)**3*(s.Symbol('lambda')-3)**5*(s.Symbol('lambda')+15)**3 or s.expand(cp-(s.Symbol('lambda')-10)**3*(s.Symbol('lambda')-3)**5*(s.Symbol('lambda')+15)**3)==0
T20=(s.Matrix.vstack(trans_fixed(10,10),trans_fixed(10,8),trans_fixed(10,9),s.zeros(5,11))).row_join(s.Matrix.vstack(s.zeros(3,9),trans_fixed(8,7),s.zeros(5,9),trans_fixed(8,8)))
det20=T20.det();assert det20!=0;cert['transvectant20_det']=str(det20);cert['transvectant20_factorization']={str(k):v for k,v in s.factorint(det20).items()}
# Previous note uses descending X powers in every source and target summand.
def reversal_sign(n):return (-1)**(n*(n-1)//2)
basis_sign=reversal_sign(11)*reversal_sign(9)*reversal_sign(3)*reversal_sign(7)*reversal_sign(5)**2
assert basis_sign==-1
cert['transvectant20_basis_convention']={'current':'ascending X powers','previous':'descending X powers','sign_of_change':basis_sign,'previous_convention_determinant':str(basis_sign*det20)}
Dp=s.Matrix.hstack(*(K-10*s.eye(11)).nullspace());Dm=s.Matrix.hstack(*(K+15*s.eye(11)).nullspace())
def trans_pair(u,v,r):
 out=s.zeros(21-2*r,1)
 for i,a in enumerate(u):
  if not a:continue
  for j,b in enumerate(v):
   if not b:continue
   row=i+j-r
   if 0<=row<out.rows:
    out[row]+=a*b*sum((-1)**k*comb(r,k)*ff(i,r-k)*ff(10-i,k)*ff(j,k)*ff(10-j,r-k) for k in range(r+1))
 return out
T8=s.Matrix.hstack(*(trans_pair(Dp[:,a],Dm[:,b],8) for a in range(3) for b in range(3)))
T6=s.Matrix.hstack(*(trans_pair(Dp[:,a],Dm[:,b],6) for a in range(3) for b in range(3)))
Amap=trans_fixed(8,8)
assert T8.rank()==5 and T6.rank()==4 and Amap*T6==s.zeros(5,9)
null=s.Matrix.hstack(*T8.nullspace());assert (T6*null).rank()==4
cert['cubic_geometry']={'T8_rank':T8.rank(),'T6_rank':T6.rank(),'F8_T6_zero':True,'T6_on_kernel_T8_rank':4}
# The monic restriction polynomial, with independent character recurrence.
t=s.symbols('t');Hs=[s.Integer(1),t]
for j in range(1,5):Hs.append(s.expand((t-1)*Hs[-1]-Hs[-2]))
p=s.factor(Hs[5]+Hs[4]-Hs[3]-2*Hs[2]-Hs[1]);assert s.expand(p-t*(t+1)*(t-3)*(t*t-t-1))==0
cert['restriction_polynomial']=str(p)
(OUT/'geometry.json').write_text(json.dumps(cert,indent=2))
print('all geometry checks passed')
