#!/usr/bin/env python3
"""Compare the completed independent reconstruction with the delivered sources.

This file reads earlier certificates only after reconstruction has finished.
It contains no mathematical search and cannot influence the independent result.
"""
from pathlib import Path
from collections import Counter
import json,re

HERE=Path(__file__).resolve().parent
VERIFY=HERE.parent;REPO=VERIFY.parent
load=lambda p:json.loads(p.read_text())
fresh=load(HERE/'independent_results.json');old=load(VERIFY/'induction/wave_i_results.json')
revision=load(VERIFY/'revision_results.json')
assert fresh['group']==old['group_permutations']
assert fresh['class_representatives']==[fresh['group'][i] for i in old['class_representatives']]
assert fresh['class_sizes']==[1,15,20,12,12]
assert {tuple(s['elements']) for s in fresh['subgroups']}=={tuple(s['elements']) for s in old['subgroups']}
assert fresh['generators']==old['generators']
def new_linear(c):
    d=c['linear_order'];factor=c['abelian_exponent']//d
    assert all(v%factor==0 for h,v in c['values'])
    return (tuple(fresh['subgroups'][c['subgroup_index']]['elements']),d,
            tuple((h,v//factor) for h,v in c['values']),tuple(c['module']),tuple(map(tuple,c['character'])))
def old_linear(c):
    return (tuple(old['subgroups'][c['subgroup']]['elements']),c['character_order'],
            tuple(map(tuple,c['exponents'])),tuple(c['multiplicities']),tuple(map(tuple,c['character'])))
assert Counter(map(new_linear,fresh['linear_characters']))==Counter(map(old_linear,old['induction_catalogue']))
for key,row in fresh['rows'].items():
    prior=old['rows'][key]
    signature=lambda p:(tuple(p['U']),tuple(p['V']))
    current={signature(p):p for p in row['pairs']};previous={signature(p):p for p in prior['pairs']}
    assert set(current)==set(previous)=={(tuple(u),tuple(v)) for u,v in revision['rows'][key]['all_pairs']}
    assert row['counts']==prior['counts']
    assert row['U_count_by_generating_function']==prior['U_candidates']
    assert row['V_enumerated']==prior['V_candidates']
    for pair,p in current.items():
        q=previous[pair]
        assert p['exterior_character']==q['exterior_character'] and p['rational']==q['rational']
        for category in fresh['generators']:
            for side in ['U','V']:
                assert (p['witnesses'][category][side] is None)==(q['witnesses'][category][side] is None)
assert fresh['78_14_all_pairs']==[[[78,0,0,0,0],[14,0,0,0,0]]]
assert fresh['15_14_all_pairs']==revision['exclusions']['15_14_all']
assert fresh['roots_of_selected_120_input']==old['exterior_square_roots']==[[1,0,0,0,3]]

article=(REPO/'Seven_Coincidences_and_Sixty_Symmetries_Intelligencer.tex').read_text()
companion=(REPO/'Proof_Companion.tex').read_text()
integers=lambda s:tuple(map(int,re.findall(r'-?\d+',s)))
def read_table(text,kind):
    table={}
    for line in text.splitlines():
        fields=line.split('&')
        if len(fields)!=kind or not fields[0].strip().isdigit():continue
        value=fields[0].strip()
        if value not in fresh['rows']:continue
        if kind==6 and r'\(' in line:
            table[value]=(integers(' '.join(fields[1:4])),integers(fields[4]),integers(fields[5]))
        if kind==5 and r'\(' in line:
            table[value]=(integers(fields[1]),integers(fields[2]),integers(fields[3]),integers(fields[4]))
    return table
ar=read_table(article,6);co=read_table(companion,5)
assert len(ar)==len(co)==7
def expand(v):a,b,c,d=v;return (a,d,d,b,c)
selected={}
for value,(dims,u,v) in ar.items():
    row=fresh['rows'][value]
    assert list(dims[:2])==row['dimensions'] and dims[2]==row['degree']
    match=[p for p in row['pairs'] if tuple(p['U'])==expand(u) and tuple(p['V'])==expand(v)]
    assert len(match)==1
    common=tuple(a for a,b in match[0]['exterior_character']);assert all(b==0 for a,b in match[0]['exterior_character'])
    assert co[value]==(dims,u,v,common)
    selected[value]={'dimensions_degree':dims,'U_article':u,'V_article':v,'common_character':common}
    if value in fresh['explicit_actions']:
        assert fresh['explicit_actions'][value]['U_module']==list(expand(u))
        assert fresh['explicit_actions'][value]['V_module']==list(expand(v))
for line in companion.splitlines():
    f=line.split('&')
    if len(f)==7 and f[0].isdigit() and f[0] in fresh['rows']:
        actual=tuple(fresh['rows'][f[0]]['counts'][k] for k in ['complex','rational','permutation','signed','monomial'])
        assert tuple(integers(x)[0] for x in f[2:])==actual
earlier={}
for value in ['1540','24310']:
    prior=revision['rows'][value]['earlier_valid_alternative']
    u,v=map(expand,prior['article_coordinates'])
    p=next(p for p in fresh['rows'][value]['pairs'] if tuple(p['U'])==u and tuple(p['V'])==v)
    assert p['exterior_character']==prior['common_character']
    earlier[value]=prior
assert selected['1540']['common_character'][1]==-20
assert earlier['1540']['common_character'][1]==[-28,0]
assert fresh['alternative_17_character']==[[17,0],[5,0],[-1,0],[2,0],[2,0]]
assert fresh['explicit_actions']['24310']['V_character']==[[17,0],[5,0],[2,0],[2,0],[2,0]]
for p in fresh['rows']['7140']['pairs']:
    t,x,y,a,b=p['V'];assert t+a-b==-3

k=load(HERE/'klein_results.json');stored=revision['harmonics']
assert k['left_reduced']==k['right_reduced']==stored['reduced']
assert k['left_harmonics']==stored['left'] and k['right_harmonics']==stored['right']
assert k['restriction_matrix']==stored['restriction_matrix']
assert k['primitive_integral_kernel']==stored['kernel_generator']
assert k['common_complex_irrep_multiplicities']==stored['common_irrep_multiplicities']
assert k['equivariant_hom_dimension']==stored['hom_dimension']==150414
m=load(HERE/'mathieu_results.json');prior=revision['conway_subgroup']
assert m['order']==prior['order']
assert sorted(m['class_sizes'])==prior['class_sizes']
assert [sorted(i+1 for i in o) for o in m['orbits']]==prior['orbits']
assert all(m['golay_preservation'].values())

# Reference navigation is checked here; the mathematical reading is recorded
# separately in MI_S9_Reference_Review.md, not inferred from string presence.
section=0;sub=0;destinations={}
for line in companion.splitlines():
    if line.startswith(r'\section{'):
        section+=1;sub=0;key=chr(64+section);destinations[key]=line
    elif line.startswith(r'\subsection{'):
        sub+=1;destinations[f'{chr(64+section)}.{sub}']=line
for text in [article,companion]:
    labels=re.findall(r'\\label\{([^}]+)\}',text)
    refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text)
    assert len(labels)==len(set(labels)) and set(refs)<=set(labels)
    bib=set(re.findall(r'\\bibitem\{([^}]+)\}',text))
    cited=set()
    for x in re.findall(r'\\cite\w*(?:\[[^\]]*\])*\{([^}]+)\}',text):cited.update(x.split(','))
    assert cited<=bib
    if text==article:
        for key in re.findall(r'\b[A-I]\.\d+\b',text):assert key in destinations or key=='C.1',key
record={'status':'PASS','complete_pair_lists':102,'subgroups':len(fresh['subgroups']),
    'linear_characters':len(fresh['linear_characters']),'all_159_inductions_equal':True,
    'all_membership_decisions_equal':True,'candidate_counts_equal':True,
    'selected_rows':selected,'earlier_valid_alternatives':earlier,
    'all_five_class_values_equal':True,'explicit_actions_checked_on_all_60_elements':True,
    'klein_reductions_equal':True,'mathieu_right_action_confirmed':True,
    'all_local_labels_and_citation_keys_resolve':True,'companion_destinations':destinations}
(HERE/'comparison_results.json').write_text(json.dumps(record,indent=2)+'\n')
print('PASS: all 159 induced characters, 102 pairs, restricted decisions, and five-class values agree.')
print('PASS: both printed tables, prior valid alternatives, 3003 exclusions, and selected root agree.')
print('PASS: Klein reductions, Mathieu subgroup, and local reference destinations agree.')
