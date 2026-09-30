#!/usr/bin/env python3
"""Compare printed rows with independently regenerated mathematical certificates."""
from pathlib import Path
import json, re
from verify_revision import ROWS, PRINTED, full, character, exterior_ch

ROOT = Path(__file__).resolve().parents[1]

def read_json(path):
    return json.loads((ROOT/path).read_text())

def printed_rows(text, companion=False):
    result={}
    for line in text.splitlines():
        if not re.match(r'^(120|210|1540|3003|7140|11628|24310) & ', line):
            continue
        fields=line.split('&')
        def tup(field):
            return tuple(map(int,re.search(r'\((-?\d+(?:,-?\d+)+)\)',field)[1].split(',')))
        if companion and len(fields)==5:
            value=int(fields[0]); dims=tuple(map(int,re.search(r'\((\d+),(\d+);(\d+)\)',fields[1]).groups()))
            result[value]=(dims,tup(fields[2]),tup(fields[3]),tup(fields[4]))
        elif not companion and len(fields)==6:
            result[int(fields[0])]=(tuple(map(int,fields[1:4])),tup(fields[4]),tup(fields[5]))
    return result

def check():
    article=(ROOT/'Seven_Coincidences_and_Sixty_Symmetries_Intelligencer.tex').read_text()
    companion=(ROOT/'Proof_Companion.tex').read_text()
    ar=printed_rows(article); cr=printed_rows(companion,True)
    revision=read_json('verification/revision_results.json')
    finite=read_json('verification/certificates/finite_core.json')
    induction=read_json('verification/induction/wave_i_results.json')
    expected_counts=[(10,4,0,1,3),(2,2,1,2,2),(11,5,1,3,11),
                     (19,3,0,2,3),(3,1,0,0,1),(26,6,0,2,8),(31,11,2,4,15)]
    summaries=[]
    for i,(value,m,n,k) in enumerate(ROWS):
        u,v=PRINTED[i]; key=str(value)
        left=exterior_ch(character(full(u)),2)
        right=exterior_ch(character(full(v)),k)
        assert left==right and all(b==0 for a,b in left)
        common=tuple(a for a,b in left)
        assert ar[value]==((m,n,k),u,v)
        assert cr[value]==((m,n,k),u,v,common)
        f=next(x for x in finite['rows'] if x['value']==value)
        assert tuple(f['U'])==u and tuple(f['V'])==v and tuple(f['character'])==common
        r=revision['rows'][key]; w=induction['rows'][key]
        rp={(tuple(a),tuple(b)) for a,b in r['all_pairs']}
        wp={(tuple(p['U']),tuple(p['V'])) for p in w['pairs']}
        assert rp==wp and (full(u),full(v)) in rp
        counts=tuple(w['counts'][s] for s in ['complex','rational','permutation','signed','monomial'])
        assert counts==expected_counts[i]
        # Check the actual restricted-count table in the companion source.
        expected=f'{value}&$({m},{n};{k})$&'+'&'.join(map(str,counts))
        assert expected in companion
        assert (r['complex'],r['rational'])==counts[:2]
        assert (r['U_candidates'],r['V_candidates'])==(w['U_candidates'],w['V_candidates'])
        if 'earlier_valid_alternative' in r:
            old=r['earlier_valid_alternative']
            assert tuple(full(x) for x in old['article_coordinates']) in rp
            alt=next(x for x in finite['earlier_valid_alternatives'] if x['value']==value)
            assert alt['character']==[x[0] for x in old['common_character']]
        summaries.append({'value':value,'U':u,'V':v,'common_character':common,'counts':counts})
    # Cross-document references and the construction equations must survive integration.
    for label in ['eq:faces-1540','eq:triples-24310','eq:lines-7140']:
        assert r'\label{'+label+'}' in article
    for text in (article,companion):
        labels=re.findall(r'\\label\{([^}]+)\}',text)
        assert len(labels)==len(set(labels)), 'duplicate labels'
        refs=set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text))
        assert refs<=set(labels), sorted(refs-set(labels))
    out={'status':'PASS','arithmetic':'exact Z[phi] and independent eigenvalue factors',
         'complete_pair_lists_equal':True,'total_nontrivial_complex_pairs':sum(x[0] for x in expected_counts),
         'selected_rows':summaries}
    (ROOT/'verification/certificates/integration.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: manuscript tables, selected rows, earlier alternatives, and all 102 pairs agree.')
    return out

if __name__=='__main__':
    check()
