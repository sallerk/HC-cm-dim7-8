# compare_blocks_g7.py -- case-by-case comparison of the two blocks checks at g = 7:
#   blocks_check/used_blocks_check_g7_cases.json  (used_blocks_check_g7.py, the CM-sixfold program run at g = 7)
#   blocks_check/verify_blocks_g7.json            (verify_blocks.py, written separately)
# Checks that each file has one record per (gid, Phi), that both key sets equal the set of reduced cases with d >= 1
# in enum_g7.json (stored flags), and that the flags strict / all / easy agree record by record.
#   usage:  python compare_blocks_g7.py > blocks_check/compare_blocks_g7_out.txt
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    recs = json.load(open(os.path.join(HERE, 'blocks_check', name)))
    keys = [(r['gid'], tuple(r['Phi'])) for r in recs]
    assert len(set(keys)) == len(keys), 'duplicate keys in ' + name
    return dict(zip(keys, recs))


ubc = load('used_blocks_check_g7_cases.json')
vb = load('verify_blocks_g7.json')
enum = json.load(open(os.path.join(HERE, 'enum_g7.json')))
expected = {(c['gid'], tuple(c['Phi'])) for c in enum['cases'] if c['reduced'] and c['d'] >= 1}
print('records: used_blocks_check %d, verify_blocks %d, reduced cases with d >= 1 in enum_g7.json %d'
      % (len(ubc), len(vb), len(expected)))
print('key sets equal to the expected set: used_blocks_check %s, verify_blocks %s'
      % (set(ubc) == expected, set(vb) == expected))
bad = 0
for key in sorted(expected):
    for f in ('strict', 'all', 'easy'):
        if ubc[key][f] != vb[key][f]:
            bad += 1
            print('DISAGREE', key, f, ubc[key][f], vb[key][f])
print('disagreements (strict, all, easy over %d cases): %d' % (len(expected), bad))
for f in ('strict', 'all', 'easy'):
    print('  %-6s true in %d cases (both programs)' % (f, sum(1 for k in expected if vb[k][f] and ubc[k][f])))
ok = set(ubc) == expected and set(vb) == expected and bad == 0
print('RESULT:', 'agree' if ok else 'DIFFER')
sys.exit(0 if ok else 1)
