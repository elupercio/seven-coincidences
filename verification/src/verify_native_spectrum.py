#!/usr/bin/env python3
"""Independent audit of the retained sparse Conway spinor.
Input tensors only are reused. Majorana action, the Klein operator, and a Krylov
Gram calculation are rebuilt. The minimal polynomial is DISCOVERED from the
Gram kernel, not imported from the expected spectral list.
Coefficients are integers in Q(i,sqrt(5)); no floating point spectral calculations.
"""
from pathlib import Path
from itertools import combinations
from fractions import Fraction
from collections import Counter
import numpy as np
import sympy as sp
import json
ROOT=Path(__file__).resolve().parents[1]
B=ROOT/'prior_inputs/Full_Conway_Klein_Comparison_Package';OUT=ROOT/'certificates'
N=4096
# Row coordinates (1,i,sqrt5,i sqrt5).
def vi(v):return np.stack((-v[:,1],v[:,0],-v[:,3],v[:,2]),axis=1)
def vs(v):return np.stack((5*v[:,2],5*v[:,3],v[:,0],v[:,1]),axis=1)
def blank():return np.zeros((N,4),dtype=object)
# Direct source -> destination Jordan-Wigner transformations.
idx=np.arange(N,dtype=int);flips=[];signs=[]
for k in range(12):
 flips.append(idx^(1<<k))
 signs.append(np.array([(-1)**((j&((1<<k)-1)).bit_count()) for j in range(N)],dtype=object))

def gamma(a,v):
 k=a//2;u=v*signs[k][:,None]
 if a%2:
  bit=np.array([1 if not (j>>k)&1 else -1 for j in range(N)],dtype=object)
  u=vi(u*bit[:,None])
 out=blank();out[flips[k],:]=u
 return out

def word(seq,v):
 for a in reversed(seq):v=gamma(a,v)
 return v

def inn(v,w):
 a,b,c,d=[v[:,i] for i in range(4)];e,f,g,h=[w[:,i] for i in range(4)]
 return tuple(int(x) for x in (sum(a*e+b*f+5*c*g+5*d*h),sum(a*f-b*e+5*c*h-5*d*g),sum(a*g+b*h+c*e+d*f),sum(a*h-b*g+c*f-d*e)))
tau=blank();data=json.loads((B/'tau_exact_0.json').read_text())
for n,a,b in data['nonzero']:tau[n,0]=a;tau[n,1]=b
assert inn(tau,tau)==(32,0,0,0)
# Sanity checks of Clifford conventions on a dense deterministic test vector.
v=blank()
for j in range(N):v[j,0]=(j*7+3)%11-5
for a in range(24):assert np.array_equal(gamma(a,gamma(a,v)),v)
for a,b in [(0,1),(1,2),(6,21),(15,23)]:assert np.array_equal(gamma(a,gamma(b,v)),-gamma(b,gamma(a,v)))
# Standard permutations as products of reflection lifts: twelve factors => scale64.
facA=[(0,3),(1,6),(2,16),(4,12),(5,8),(7,14),(9,18),(10,17),(11,20),(13,15),(19,23),(21,22)]
facB=[(0,3),(3,5),(1,20),(20,13),(2,8),(8,14),(4,17),(17,9),(12,16),(16,15),(18,23),(23,22)]
for fac in [facA,facB]:
 t=tau.copy()
 for a,b in reversed(fac):t=gamma(a,t)-gamma(b,t)
 assert np.array_equal(t,64*tau)
# Candidate stabilizer masks are inputs; all properties and eigenactions are tested.
candidates=json.loads((B/'golay_stabilizer_masks.json').read_text())
codegen=candidates['Golay_generator_bitmasks'];span={0};eigen=[]
for g in codegen:
 span|={a^g for a in list(span)}
 q=word([a for a in range(24) if (g>>a)&1],tau)
 assert np.array_equal(q,tau) or np.array_equal(q,-tau)
 eigen.append(1 if np.array_equal(q,tau) else -1)
assert len(span)==4096
assert all((a&b).bit_count()%2==0 for a in codegen for b in codegen)
wd=Counter(a.bit_count() for a in span)
assert min(n for n in wd if n)>4 and all(n%4==0 for n in wd)
for k in [2,4]:
 for cc in combinations(range(24),k):
  mask=sum(1<<a for a in cc)
  assert any((mask&g).bit_count()%2 for g in codegen)
# Reality phase.
con=tau.copy();con[:,1]*=-1
assert np.array_equal(word(list(range(1,24,2)),con),vi(tau))
print('Clifford relations, generator lifts, code and N1 constraints checked',flush=True)
# Retained marked orbit pairing and Klein operator, independently tested below.
bd=json.loads((B/'basis_data.json').read_text());X,Y=bd['orbits']
M0=np.array(bd['A0_times4'],dtype=object);M1=np.array(bd['A1_times8'],dtype=object)
assert all(x%2==0 for x in M1.flat);M1=M1//2
assert np.array_equal(M0,M0.T) and np.array_equal(M1,M1.T)
# Fresh coefficient pair action of 8K; precompute source permutation/sign
# by tracking a symbolic unit coefficient rather than reusing archived transforms.
terms=[]
for a in range(12):
 for b in range(12):
  c0=int(M0[a,b]);c1=int(M1[a,b])
  if not(c0 or c1):continue
  x,y=X[a],Y[b];dest=[];sg=[];iphase=(x%2+y%2+1)%4
  for j in range(N):
   mid=j^(1<<(y//2));target=mid^(1<<(x//2))
   e=((j&((1<<(y//2))-1)).bit_count()+(mid&((1<<(x//2))-1)).bit_count())
   if y%2:e+=(j>>(y//2))&1
   if x%2:e+=(mid>>(x//2))&1
   sign=(-1)**e
   if iphase>=2:sign=-sign
   dest.append(target);sg.append(sign)
  terms.append((np.array(dest),np.array(sg,dtype=object),iphase%2,c0,c1))

def D(v):
 out=blank()
 for dest,sg,ip,c0,c1 in terms:
  u=v*sg[:,None]
  if ip:u=vi(u)
  out[dest,:]+=c0*u+c1*vs(u)
 return out
raw=blank()
for x in X:
 for y in Y:raw+=gamma(x,gamma(y,tau))
v=12*tau-vi(raw)
assert inn(v,v)==(9216,0,0,0)
powers=[v]
for k in range(10):
 powers.append(D(powers[-1]));print('Krylov vector',k+1,flush=True)
# Hermitian Gram matrix determines the minimal polynomial on this cyclic span.
GG=[]
for i in range(11):
 row=[]
 for j in range(11):
  a=inn(powers[i],powers[j]);assert a[1:]==(0,0,0);row.append(a[0])
 GG.append(row)
Gram=sp.Matrix(GG);null=Gram.nullspace();assert len(null)==1
coef=null[0]/null[0][-1];assert all(a.q==1 for a in coef)
res=blank()
for a,p in zip(coef,powers):res+=int(a)*p
assert not any(res.flat)
x=sp.symbols('x');poly=sum(int(a)*x**i for i,a in enumerate(coef));fac=sp.factor_list(poly)
eigs=[]
for f,m in fac[1]:
 assert sp.degree(f)==1 and m==1
 eigs.append(-f.subs(x,0)/f.coeff(x))
eigs=sorted(eigs);assert all(a.q==1 and a%8==0 for a in eigs)
charges=[int(a)//8 for a in eigs];weights=[]
for a in eigs:
 pp=sp.Poly(sp.cancel(poly/(x-a)),x);den=pp.eval(a)
 num=sum(pp.nth(k)*Gram[0,k] for k in range(10))
 wt=sp.cancel(num/(den*Gram[0,0]));assert wt>0;weights.append(str(wt))
assert sum(sp.Rational(a) for a in weights)==1
support_pairs=[(a,b) for a in charges for b in charges if a+b==4]
assert support_pairs==[(2,2)]
cert={'spinor_input':'tau_exact_0.json','norm_before_normalization':32,'M24_lift_invariance':True,'Majorana_reality_phase':'i','code_dimension':12,'weight_enumerator':dict(wd),'degree2_degree4_constraints':276+10626,'marked_orbit_pairing':'basis_data.json','Krylov_Gram_rank':11-len(null),'minimal_polynomial_8K':str(sp.factor(poly)),'K_charge_spectrum_discovered':charges,'normalized_projection_norms':weights,'charge4_pairs':support_pairs,'charge2_physical_norm_squared':str(8*sp.Rational(weights[charges.index(2)])),'checks_are_exact':True,'native_operator_scope':'projective small-sector theorem after stated topological twist; not full 3003'}
(OUT/'native_spectrum.json').write_text(json.dumps(cert,indent=2));print(json.dumps(cert,indent=2))
