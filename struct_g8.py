# struct_g8.py -- GAP StructureDescription / IsAbelian for the (G,rho) pairs carrying primitive degenerate types.
import json, os
from sage.all import libgap
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(HERE + '/trans_g8.json'))
out = {}
for G in D['groups']:
    if any(c['primitive'] and c['d'] >= 1 for c in G['cases']):
        H = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in G['gens']])
        out['%s/%d' % (G['TI'], G['j'])] = dict(struct=str(libgap.StructureDescription(H)), abelian=bool(libgap.IsAbelian(H)))
json.dump(out, open(HERE + '/struct_g8.json', 'w'))
print(len(out))
