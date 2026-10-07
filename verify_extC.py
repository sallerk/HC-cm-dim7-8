# verify_extC.py -- INDEPENDENT VERIFIER (pure Python) for ext_weil.py Stage C mode (g=7, D up to Dmax):
# per-level (D, rank, saturation index) from verify_generalC.chars_by_D (own subgroups / cosets / lattice code).
import json, sys
from verify_enum import hnf, sat_index, rank
from verify_generalC import chars_by_D
Dmax = int(sys.argv[1]); recs = json.load(open(sys.argv[2]))
E = json.load(open('enum_g7.json')); groups = {G['gid']: G for G in E['groups']}
agree = 0
for r in recs:
    byD = chars_by_D(groups[r['gid']], 7, r['Phi'], Dmax)
    acc = []; lv = []
    for D in sorted(byD):
        for w in byD[D]:
            assert rank(r['LU'] + [list(w)]) == r['d']
        acc += [list(w) for w in byD[D]]
        L = hnf(acc)
        lv.append((D, len(L), sat_index(L) if len(L) == r['d'] else None))
    prod_ = [(l['D'], l['rank'], l['satindex']) for l in r['levels']]
    ok = lv == prod_; agree += ok
    print('gid', r['gid'], r['Phi'], 'AGREE' if ok else 'DISAGREE', lv, '' if ok else prod_, flush=True)
print('Stage C extended verification (Dmax=%d): %d/%d agree' % (Dmax, agree, len(recs)))
