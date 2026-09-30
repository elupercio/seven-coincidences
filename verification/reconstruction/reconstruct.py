#!/usr/bin/env python3
"""Independent MI-S9 reconstruction; standard library, exact arithmetic.

No earlier verifier or result is imported. Subgroup characters are propagated
from generator assignments into the exponent of the abelianization. Induction
is computed on right cosets. Large inputs are recovered from exterior-square
characters, rather than enumerated. Comparison with prior files is separate.
"""
from collections import Counter, defaultdict
from functools import cache
from itertools import permutations, product, combinations
from math import comb, isqrt, lcm
from pathlib import Path
import json

OUT=Path(__file__).resolve().parent
ZERO=(0,0); ONE=(1,0)
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[0],-a[1])
def mul(a,b):return (a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def scale(n,a):return (n*a[0],n*a[1])
def div(a,n):
    assert a[0]%n==a[1]%n==0
    return (a[0]//n,a[1]//n)
def conjugate(a):return (a[0]+a[1],-a[1])
def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def compose(p,q):return tuple(q[x] for x in p)  # p acts first
def inverse(p):return tuple(p.index(i) for i in range(len(p)))

G=sorted(p for p in permutations(range(5)) if sign(p)==1)
idx={p:i for i,p in enumerate(G)}; identity=idx[tuple(range(5))]
table=[[idx[compose(p,q)] for q in G] for p in G]
inverses=[idx[inverse(p)] for p in G]
def closure(gs):
    found={identity}; queue=[identity]
    for p in queue:
        for q in gs:
            r=table[p][q]
            if r not in found:found.add(r);queue.append(r)
    return frozenset(found)
def power(g,n):
    answer=identity
    for _ in range(n):answer=table[answer][g]
    return answer
orders=[next(n for n in range(1,61) if power(g,n)==identity) for g in range(60)]
def conj(g,h):return table[table[inverses[h]][g]][h]
remaining=set(range(60)); classes=[]
while remaining:
    g=min(remaining); cl=frozenset(conj(g,h) for h in range(60))
    classes.append(cl);remaining-=cl
classes.sort(key=lambda c:(orders[min(c)],min(c)))
classof={g:i for i,c in enumerate(classes) for g in c}
reps=list(map(min,classes)); sizes=list(map(len,classes))
power_classes=[[classof[power(g,j)] for j in range(9)] for g in reps]

# Subgroup lattice closure, with its completeness criterion checked explicitly.
subgroup_gens={frozenset([identity]):()}; queue=list(subgroup_gens)
for H in queue:
    for g in range(60):
        if g in H:continue
        K=closure(subgroup_gens[H]+(g,))
        if K not in subgroup_gens:subgroup_gens[K]=subgroup_gens[H]+(g,);queue.append(K)
subgroups=sorted(subgroup_gens,key=lambda H:(len(H),tuple(sorted(H))))

def coset_action(H):
    left=set(range(60)); cosets=[]; representatives=[]
    while left:
        r=min(left); C=frozenset(table[h][r] for h in H)
        cosets.append(C);representatives.append(r);left-=C
    locate={g:i for i,C in enumerate(cosets) for g in C}
    maps=[]; multipliers=[]
    for g in range(60):
        image=[]; elements=[]
        for r in representatives:
            rg=table[r][g]; j=locate[rg]
            h=table[rg][inverses[representatives[j]]]
            assert h in H
            image.append(j);elements.append(h)
        maps.append(image);multipliers.append(elements)
    return cosets,representatives,maps,multipliers

def rootsum(coefficients):
    n=len(coefficients);a=coefficients
    if n==1:return (a[0],0)
    if n==2:return (a[0]-a[1],0)
    if n==3:
        assert a[1]==a[2]
        return (a[0]-a[1],0)
    if n==5:
        assert a[1]==a[4] and a[2]==a[3]
        return (a[0]-a[1],a[1]-a[2])
    if n==6:
        assert a[1]+a[2]-a[4]-a[5]==0
        return (a[0]+a[1]-a[3]-a[4],0)
    raise ValueError(n)

c5=next(g for g in range(60) if orders[g]==5)
c3=next(g for g in range(60) if orders[g]==3)
H5=closure((c5,)); H3=closure((c3,))
H10=frozenset(g for g in range(60) if frozenset(conj(h,g) for h in H5)==H5)
H12=frozenset(g for g,p in enumerate(G) if p[4]==4)
H6=frozenset(g for g,p in enumerate(G) if set(p[:2])=={0,1})
models={}
for name,H in [('R5',H12),('P',H10),('E',H6),('Q',H5),('R20',H3),('R60',frozenset([identity]))]:
    cs,rs,acts,hs=coset_action(H)
    chars=tuple((sum(i==j for i,j in enumerate(acts[g])),0) for g in reps)
    models[name]={'H':H,'cosets':cs,'representatives':rs,'actions':acts,'multipliers':hs,'character':chars}

# Rational models A=R5-1, B=P-1, C=Q-P; D+/D- are the two rotation characters.
trivial=(ONE,)*5
A=tuple(add(x,neg(y)) for x,y in zip(models['R5']['character'],trivial))
B=tuple(add(x,neg(y)) for x,y in zip(models['P']['character'],trivial))
C=tuple(add(x,neg(y)) for x,y in zip(models['Q']['character'],models['P']['character']))
D=((3,0),(-1,0),ZERO,(0,1),(1,-1))
Dp=tuple(conjugate(x) for x in D)
irreducibles=(trivial,D,Dp,A,B)
dimensions=tuple(c[0][0] for c in irreducibles)
def inner(c,d):
    v=ZERO
    for n,x,y in zip(sizes,c,d):v=add(v,scale(n,mul(x,y)))
    return div(v,60)
assert sum(d*d for d in dimensions)==60
assert all(inner(c,d)==(int(i==j),0) for i,c in enumerate(irreducibles) for j,d in enumerate(irreducibles))
assert tuple(add(x,y) for x,y in zip(D,Dp))==C
def character(v):
    return tuple(tuple(sum(v[j]*irreducibles[j][i][k] for j in range(5)) for k in range(2)) for i in range(5))
def decompose(c):
    nums=[]
    for d in irreducibles:
        v=ZERO
        for n,x,y in zip(sizes,c,d):v=add(v,scale(n,mul(x,y)))
        if v[1] or v[0]%60 or v[0]<0:return None
        nums.append(v[0]//60)
    nums=tuple(nums)
    return nums if character(nums)==tuple(c) else None
def exterior(c,k):
    result=[]
    for ci in range(5):
        e=[ONE]
        for degree in range(1,k+1):
            total=ZERO
            for j in range(1,degree+1):
                total=add(total,scale((-1)**(j+1),mul(e[degree-j],c[power_classes[ci][j]])))
            e.append(div(total,degree))
        result.append(e[-1])
    return tuple(result)
def modules(d):
    def visit(i,left,prefix):
        if i==len(dimensions)-1:
            if left%dimensions[i]==0:yield tuple(prefix+[left//dimensions[i]])
            return
        for n in range(left//dimensions[i]+1):yield from visit(i+1,left-n*dimensions[i],prefix+[n])
    return visit(0,d,[])
def module_count(d):
    dp=[1]+[0]*d
    for n in dimensions:
        for k in range(n,d+1):dp[k]+=dp[k-n]
    return dp[d]

# All homomorphisms H -> mu_L, L the exponent of H/[H,H].
catalogue=[]; subgroup_records=[]
generators={name:set() for name in ['permutation','signed','monomial']}
for hi,H in enumerate(subgroups):
    gens=subgroup_gens[H]
    comms=[table[table[table[inverses[h]][inverses[k]]][h]][k] for h in H for k in H]
    K=closure(comms)
    L=lcm(*(next(j for j in range(1,61) if power(h,j) in K) for h in H))
    linears=[]
    for assignment in product(range(L),repeat=len(gens)):
        vals={identity:0};queue=[identity];valid=True
        for h in queue:
            for g,value in zip(gens,assignment):
                hg=table[h][g];v=(vals[h]+value)%L
                if hg in vals:
                    if vals[hg]!=v:valid=False;break
                else:vals[hg]=v;queue.append(hg)
            if not valid:break
        if valid:
            assert set(vals)==set(H)
            assert all(vals[table[h][k]]==(vals[h]+vals[k])%L for h in H for k in H)
            linears.append(vals)
    assert len(linears)==len(H)//len(K)
    cs,rs,acts,hs=coset_action(H)
    for vals in linears:
        actual_order=lcm(*(L//__import__('math').gcd(L,v) for v in vals.values()))
        phases=[[vals[h] for h in row] for row in hs]
        # Right-action law against generators of the whole group, at every element.
        whole_gens=subgroup_gens[frozenset(range(60))]
        for g in range(60):
            for h in whole_gens:
                gh=table[g][h]
                assert all(acts[h][acts[g][i]]==acts[gh][i] and
                    (phases[g][i]+phases[h][acts[g][i]])%L==phases[gh][i] for i in range(len(cs)))
        cc=[]
        for g in reps:
            coefficients=[0]*L
            for i,j in enumerate(acts[g]):
                if i==j:coefficients[phases[g][i]]+=1
            cc.append(rootsum(coefficients))
        v=decompose(cc);assert v is not None
        generators['monomial'].add(v)
        if actual_order<=2:generators['signed'].add(v)
        if actual_order==1:generators['permutation'].add(v)
        catalogue.append({'subgroup_index':hi,'order':len(H),'abelian_exponent':L,
            'linear_order':actual_order,'values':[[h,vals[h]] for h in sorted(H)],
            'character':cc,'module':v})
    subgroup_records.append({'elements':sorted(H),'generators':gens,'order':len(H),
        'commutator':sorted(K),'abelian_exponent':L,'linear_count':len(linears)})
generators={name:sorted(values) for name,values in generators.items()}
print('Subgroups:',len(subgroups),'linear characters:',len(catalogue),'generators:',{k:len(v) for k,v in generators.items()},flush=True)

# Membership by subtraction. Every successful step lowers total multiplicity;
# all eligible generators are visited, making failures exhaustive as well.
def membership_factory(basis):
    ordered=sorted(enumerate(basis),key=lambda p:(-sum(p[1]),p[1]))
    @cache
    def solve(v):
        if not any(v):return (0,)*len(basis)
        for index,g in ordered:
            if all(a>=b for a,b in zip(v,g)):
                answer=solve(tuple(a-b for a,b in zip(v,g)))
                if answer is not None:
                    out=list(answer);out[index]+=1;return tuple(out)
        return None
    return solve
members={name:membership_factory(gs) for name,gs in generators.items()}

def integer_square_roots(n):
    if n<0:return ()
    r=isqrt(n)
    return tuple(sorted({r,-r})) if r*r==n else ()
def triangle_roots(n):
    return tuple(sorted({(1+r)//2 for r in integer_square_roots(1+8*n) if (1+r)%2==0}))
def lift_exterior_square(target,m):
    """Recover all characters of dimension m whose exterior square is target.

    At order 5 write x=a+b*phi. Then 2c=a²+b²-a-b and
    2d=b(2a+b+1). Necessarily |b|=|mult(D+)-mult(D-)| <= m/3.
    """
    if target[0]!=(comb(m,2),0):return []
    assert target[1][1]==target[2][1]==0 and target[4]==conjugate(target[3])
    c,d=target[3];fifths=set()
    for b in range(-(m//3),m//3+1):
        if b==0:
            possible=triangle_roots(c) if d==0 else ()
        elif (2*d)%b==0:
            numerator=2*d//b-b-1
            possible=(numerator//2,) if numerator%2==0 else ()
        else:possible=()
        for a in possible:
            if a*a+b*b-a-b==2*c and b*(2*a+b+1)==2*d:fifths.add((a,b))
    answer=set()
    for s,t,x in product(integer_square_roots(m+2*target[1][0]),triangle_roots(target[2][0]),sorted(fifths)):
        cc=((m,0),(s,0),(t,0),x,conjugate(x));v=decompose(cc)
        if v is not None:
            assert exterior(cc,2)==target
            answer.add(v)
    return sorted(answer)

ROWS=[(120,16,10,3),(210,21,10,4),(1540,56,22,3),(3003,78,15,5),
      (7140,120,36,3),(11628,153,19,5),(24310,221,17,8)]
rows={}
for value,m,n,k in ROWS:
    found=[]
    for v in modules(n):
        target=exterior(character(v),k)
        for u in lift_exterior_square(target,m):
            if u==(m,0,0,0,0) and v==(n,0,0,0,0):continue
            witnesses={name:{'U':check(u),'V':check(v)} for name,check in members.items()}
            for name,uv in witnesses.items():
                for label,target_tuple in [('U',u),('V',v)]:
                    if uv[label] is not None:
                        assert tuple(sum(a*g[j] for a,g in zip(uv[label],generators[name])) for j in range(5))==target_tuple
            found.append({'U':u,'V':v,'exterior_character':target,'rational':u[1]==u[2] and v[1]==v[2],
                'witnesses':witnesses})
    found.sort(key=lambda p:(p['U'],p['V']))
    counts={'complex':len(found),'rational':sum(p['rational'] for p in found)}
    counts.update({name:sum(p['witnesses'][name]['U'] is not None and p['witnesses'][name]['V'] is not None for p in found) for name in members})
    rows[str(value)]={'dimensions':[m,n],'degree':k,'counts':counts,'pairs':found,
        'U_count_by_generating_function':module_count(m),'V_enumerated':module_count(n)}
    print(value,counts,flush=True)

excluded_78_14=[(u,v) for v in modules(14) for u in lift_exterior_square(exterior(character(v),6),78)]
v14=defaultdict(list)
for v in modules(14):v14[exterior(character(v),6)].append(v)
excluded_15_14=[(u,v) for u in modules(15) for v in v14[exterior(character(u),5)]]
root=lift_exterior_square(character((3,6,6,9,9)),16)

# Explicit monomial matrices, encoded as permutation and sixth-root exponents.
def permutation_model(name):
    acts=models[name]['actions']
    return [(p,[0]*len(p)) for p in acts]
def trivial_model(n=1):return [(list(range(n)),[0]*n) for _ in G]
def direct_sum(*parts):
    answer=[]
    for g in range(60):
        perm=[];phase=[]
        for part in parts:
            p,e=part[g];offset=len(perm)
            perm.extend(j+offset for j in p);phase.extend(e)
        answer.append((perm,phase))
    return answer
def exterior_model(model,k):
    basis=list(combinations(range(len(model[0][0])),k));indices={b:i for i,b in enumerate(basis)}
    answer=[]
    for p,e in model:
        perm=[];phase=[]
        for b in basis:
            image=tuple(p[i] for i in b)
            perm.append(indices[tuple(sorted(image))])
            phase.append((sum(e[i] for i in b)+(3 if sign(image)<0 else 0))%6)
        answer.append((perm,phase))
    return basis,answer
def induced_model(H):
    entry=next(c for c in catalogue if frozenset(subgroups[c['subgroup_index']])==H and c['linear_order']==3)
    L=entry['abelian_exponent'];values=dict(entry['values'])
    cs,rs,acts,hs=coset_action(H)
    return [(p,[values[h]*(6//L) for h in hh]) for p,hh in zip(acts,hs)],entry
def monomial_trace(action):
    p,e=action;coeff=[0]*6
    for i,j in enumerate(p):
        if i==j:coeff[e[i]]+=1
    return rootsum(coeff)
def determinant_coefficient(action,k):
    """[t^k] det(1+tM), by weighted cycles, independent of Newton's recurrence."""
    p,e=action;unseen=set(range(len(p)));cycles=[]
    while unseen:
        i=min(unseen);length=0;phase=0
        while i in unseen:
            unseen.remove(i);length+=1;phase+=e[i];i=p[i]
        cycles.append((length,phase%6))
    poly=[[0]*6 for _ in range(k+1)];poly[0][0]=1
    for length,phase in cycles:
        coefficient=(-1)**(length+1)
        for degree in range(k,length-1,-1):
            for a in range(6):poly[degree][(a+phase)%6]+=coefficient*poly[degree-length][a]
    return rootsum(poly[k])
def check_action(model):
    for g in range(60):
        p,e=model[g]
        for h in subgroup_gens[frozenset(range(60))]:
            q,f=model[h];r,z=model[table[g][h]]
            assert all(q[p[i]]==r[i] and (e[i]+f[p[i]])%6==z[i] for i in range(len(p)))
    return True

pm={name:permutation_model(name) for name in models}
Bmodel,Bentry=induced_model(H12);Tmodel,Tentry=induced_model(H3)
Wmodel=direct_sum(trivial_model(),Bmodel,Bmodel,Bmodel)
wedges,W2model=exterior_model(Wmodel,2)
triples,Q3model=exterior_model(pm['Q'],3)
explicit={
    '1540':(direct_sum(pm['P'],pm['E'],pm['R20'],pm['R20']),direct_sum(pm['R5'],pm['R5'],pm['Q']),3),
    '24310':(direct_sum(trivial_model(),Q3model),direct_sum(pm['R5'],pm['P'],pm['P']),8),
    '7140':(W2model,direct_sum(Wmodel,Tmodel),3)}
action_records={}
for value,(u,v,k) in explicit.items():
    assert check_action(u) and check_action(v)
    chars=[];uch=[];vch=[]
    for g in range(60):
        a=determinant_coefficient(u[g],2);b=determinant_coefficient(v[g],k)
        assert a==b,(value,g,a,b)
        chars.append(a);uch.append(monomial_trace(u[g]));vch.append(monomial_trace(v[g]))
    cm=tuple(chars[g] for g in reps);cu=tuple(uch[g] for g in reps);cv=tuple(vch[g] for g in reps)
    ut=decompose(cu);vt=decompose(cv)
    assert any(tuple(row['U'])==ut and tuple(row['V'])==vt and tuple(row['exterior_character'])==cm for row in rows[value]['pairs'])
    action_records[value]={'U_module':ut,'V_module':vt,'U_character':cu,'V_character':cv,'exterior_character':cm,
        'all_60_U_traces':uch,'all_60_V_traces':vch,'all_60_exterior_traces':chars,
        'U_action':u,'V_action':v,'group_law_checked':True}
    print('Weighted-cycle check at',value,cm,'on all 60 elements.',flush=True)

# Transport orientations on the actual 220 wedge labels, rather than replace
# the action by a previously known permutation decomposition.
unseen=set(range(len(triples)));triple_orbits=[]
while unseen:
    i=min(unseen);orb={Q3model[g][0][i] for g in range(60)}
    stabilizer=[g for g in range(60) if Q3model[g][0][i]==i]
    assert all(Q3model[g][1][i]==0 for g in stabilizer)
    orientation={j:Q3model[next(g for g in range(60) if Q3model[g][0][i]==j)][1][i] for j in orb}
    assert all((orientation[j]+Q3model[g][1][j]-orientation[Q3model[g][0][j]])%6==0 for g in range(60) for j in orb)
    triple_orbits.append({'representative':triples[i],'size':len(orb),'stabilizer':stabilizer,
        'elements':[triples[j] for j in sorted(orb)],'orientation_exponents':[[j,orientation[j]] for j in sorted(orb)]})
    unseen-=orb
assert sorted(o['size'] for o in triple_orbits)==[20,20,60,60,60]
for g in range(60):
    assert monomial_trace(Q3model[g])==add(scale(2,monomial_trace(pm['R20'][g])),scale(3,monomial_trace(pm['R60'][g])))
    assert monomial_trace(W2model[g])==add(scale(3,monomial_trace(pm['R20'][g])),scale(3,monomial_trace(Tmodel[g])))
vtilde=direct_sum(trivial_model(),Wmodel)
assert all(determinant_coefficient(vtilde[g],8)==determinant_coefficient(explicit['24310'][1][g],8) for g in range(60))
assert any(monomial_trace(vtilde[g])!=monomial_trace(explicit['24310'][1][g]) for g in range(60))

# The rational odd vertex module: Q -> P has fibers of size two, equivariantly.
axis_of={g:i for i,coset in enumerate(models['P']['cosets']) for g in coset}
projection=[axis_of[r] for r in models['Q']['representatives']]
assert sorted(Counter(projection).values())==[2]*6
assert all(projection[pm['Q'][g][0][i]]==pm['P'][g][0][projection[i]] for g in range(60) for i in range(12))

result={'method':'commutator quotients, right-coset induction, inverse exterior-square reconstruction',
    'tuple_order':['1','D_plus','D_minus','A','B'],'group':G,'classes':[sorted(c) for c in classes],
    'class_representatives':[G[r] for r in reps],'class_sizes':sizes,'power_classes':power_classes,
    'irreducible_characters':irreducibles,'subgroups':subgroup_records,'linear_characters':catalogue,
    'generators':generators,'rows':rows,'78_14_all_pairs':excluded_78_14,'15_14_all_pairs':excluded_15_14,
    'roots_of_selected_120_input':root,'explicit_actions':action_records,'triple_orbits':triple_orbits,
    'B_cubic_induction':Bentry,'T_cubic_induction':Tentry,
    'alternative_17_character':[monomial_trace(vtilde[g]) for g in reps],
    'alternative_17_equal_eighth_exteriors':True,'vertex_to_axis_projection':projection}
(OUT/'independent_results.json').write_text(json.dumps(result,indent=2)+'\n')
print('Independent enumeration written; prior certificates have not been read.',flush=True)
