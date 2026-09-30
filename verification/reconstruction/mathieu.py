#!/usr/bin/env python3
"""Reconstruct the printed right-action words point by point and certify them."""
from collections import Counter
from itertools import combinations
from math import lcm
from pathlib import Path
import hashlib,json

ROOT=Path(__file__).resolve().parent
def cycles(*cs):
    p=list(range(24))
    for c in cs:
        for x,y in zip(c,c[1:]+c[:1]):p[x-1]=y-1
    return tuple(p)
a=cycles((1,4),(2,7),(3,17),(5,13),(6,9),(8,15),(10,19),(11,18),(12,21),(14,16),(20,24),(22,23))
b=cycles((1,4,6),(2,21,14),(3,9,15),(5,18,10),(13,17,16),(19,24,23))
word_y='ababbababbababbababb'
word_w='bbbbabbbbabbabbbbababbbaabbbbabababababb'
def evaluate(text):
    image=[];paths=[]
    for start in range(24):
        x=start;path=[x]
        for letter in text:x={'a':a,'b':b}[letter][x];path.append(x)
        image.append(x);paths.append(path)
    return tuple(image),paths
y,y_paths=evaluate(word_y);w,w_paths=evaluate(word_w)
wi=tuple(w.index(i) for i in range(24))
z=tuple(w[y[wi[i]]] for i in range(24))
def compose(p,q):return tuple(q[x] for x in p)
def inv(p):return tuple(p.index(i) for i in range(24))
e=tuple(range(24));group={e};queue=[e]
for p in queue:
    for q in [a,z]:
        pq=compose(p,q)
        if pq not in group:group.add(pq);queue.append(pq)
def cycle_shape(p):
    unseen=set(range(24));lengths=[]
    while unseen:
        i=min(unseen);n=0
        while i in unseen:unseen.remove(i);n+=1;i=p[i]
        lengths.append(n)
    return tuple(sorted(lengths))
def order(p):return lcm(*cycle_shape(p))
classes=[];left=set(group)
while left:
    p=min(left);cl={compose(compose(inv(q),p),q) for q in group}
    classes.append(cl);left-=cl
classes.sort(key=lambda cl:(order(min(cl)),min(cl)))
orbits=[];left=set(range(24))
while left:
    i=min(left);o={g[i] for g in group};orbits.append(sorted(o));left-=o
assert len(group)==60 and sorted(map(len,classes))==[1,12,12,15,20]
assert Counter(map(order,group))=={1:1,2:15,3:20,5:24}
assert [len(o) for o in orbits]==[12,12]
class_sizes=[len(cl) for cl in classes]
normal_candidates=[1+sum(c) for n in range(1,4) for c in combinations(class_sizes[1:],n)]
assert not any(60%n==0 for n in normal_candidates)
orbit_characters=[[sum(g[i]==i for i in o) for g in map(min,classes)] for o in orbits]
assert orbit_characters==[[12,0,0,2,2],[12,0,0,2,2]]

data_path=ROOT.parent/'prior_inputs/Full_Conway_Klein_Comparison_Package/basis_data.json'
marked=json.loads(data_path.read_text())['orbits']
pairing=dict(zip(*marked))
assert set(pairing)==set(orbits[0]) and set(pairing.values())==set(orbits[1])
assert all(pairing[g[i]]==g[pairing[i]] for g in group for i in orbits[0])
golay_path=data_path.with_name('golay_stabilizer_masks.json')
gens=json.loads(golay_path.read_text())['Golay_generator_bitmasks']
code={0}
for mask in gens:code|={x^mask for x in list(code)}
def move(mask,p):return sum(1<<p[i] for i in range(24) if mask>>i&1)
preservation={name:all(move(mask,p) in code for mask in gens) for name,p in [('a',a),('b',b),('y',y),('w',w),('z',z)]}
assert len(code)==4096
assert Counter(x.bit_count() for x in code)=={0:1,8:759,12:2576,16:759,24:1}
assert preservation['a'] and preservation['z']
result={'convention':'right action; a acts before b in ab','generating_pair':['a','z'],
    'words':{'y':word_y,'w':word_w,'z':'w^(-1) y w'},'permutations':{'a':a,'b':b,'y':y,'w':w,'z':z},
    'point_paths':{'y':y_paths,'w':w_paths},'order':len(group),'elements':sorted(group),
    'element_order_distribution':dict(Counter(map(order,group))),'class_sizes':class_sizes,
    'class_representatives':list(map(min,classes)),'orbits':orbits,'orbit_characters':orbit_characters,
    'stabilizer_orders':[sum(g[o[0]]==o[0] for g in group) for o in orbits],
    'marked_orbit_pairing':sorted(pairing.items()),'marked_pairing_equivariant':True,
    'y_in_subgroup':y in group,'w_in_subgroup':w in group,
    'golay_dimension':12,'golay_weight_distribution':dict(Counter(x.bit_count() for x in code)),
    'golay_preservation':preservation,'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [data_path,golay_path]}}
(ROOT/'mathieu_results.json').write_text(json.dumps(result,indent=2)+'\n')
print('Right-action <a,z>: order',len(group),'; classes',class_sizes,'; orbits',[len(o) for o in orbits])
print('Orbit characters:',orbit_characters,'; marked pairing equivariant.')
print('y in <a,z>:',y in group,'; w in <a,z>:',w in group)
print('Golay code preservation:',preservation)
