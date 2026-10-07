# export_g8.py -- one JSON record per W(B8)-class (N_W(G)-class) of primitive degenerate degree-16 CM types.
import json, os
from summarize_g8 import classify, deg0
T = json.load(open('trans_g8.json'))
ST = json.load(open('struct_g8.json'))
ext = {}
for fn in ('ext_weil_D32.json', 'ext_weil_D48_sel.json'):
    for r in json.load(open(fn)):
        ext[(r['TI'], r['j'], tuple(r['Phi']))] = r['levels']
out = []
for G in T['groups']:
    seen = {}
    for c in G['cases']:
        if c['primitive'] and c['d'] >= 1:
            seen.setdefault(c['Nclass'], []).append(c)
    for ncl, lst in seen.items():
        c = lst[0]
        out.append(dict(TI=G['TI'], j=G['j'], order=G['order'], structure=ST['%s/%d' % (G['TI'], G['j'])]['struct'],
                        gens=G['gens'], rho=[(i + 8) % 16 for i in range(16)], Phi=c['Phi'],
                        Phi_all_G_orbit_reps=[x['Phi'] for x in lst], d=c['d'], dimMT=c['dimMT'], LU=c['LU'],
                        IQ=c['IQ'], CMsub=c['CMsub'], weil_D16=c['weil'], n0=c['n0'], n0_vec=c['n0_vec'],
                        deg0=deg0(c['LU']), ext_levels=ext.get((G['TI'], G['j'], tuple(c['Phi']))),
                        classification=classify(c['weil'])))
json.dump(out, open('deg_classes_g8.json', 'w'), indent=None)
print(len(out))
