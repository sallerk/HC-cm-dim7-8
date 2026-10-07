# compact_g8.py -- compact one-line-per-class table of the 257 primitive degenerate degree-16 W-classes -> table_g8_compact.txt.
import json, collections
R = json.load(open('deg_classes_g8.json'))
f = lambda x: '-' if x is None else str(x)
out = open('table_g8_compact.txt', 'w')
def P(s):
    print(s); out.write(s + '\n')
P("TI/j      |G|   orb Phi      d IQ-t     K'(deg:#,*=bal)  Dlat/Dq IQ|4|8|all n0 dg cls   LambdaU (HNF rows, signed digits)")
R.sort(key=lambda r: (r['d'], r['classification'], r['order'], int(r['TI'][3:]), r['j'], r['Phi']))
for r in R:
    w = r['weil_D16']
    degs = collections.Counter(s['deg'] for s in r['CMsub'])
    bal = collections.Counter(s['deg'] for s in r['CMsub'] if s['balanced'])
    ks = ' '.join('%d:%d%s' % (dg, degs[dg], '*' * bal.get(dg, 0)) for dg in sorted(degs)) or '-'
    tiers = '%s/%s|%s/%s|%s/%s|%s/%s' % (f(w['Dlat_IQ']), f(w['Dq_IQ']), f(w['Dlat_le4']), f(w['Dq_le4']), f(w['Dlat_le8']), f(w['Dq_le8']), f(w['Dlat_all']), f(w['Dq_all']))
    lu = ';'.join(''.join('%+d' % x for x in row) for row in r['LU'])
    iq = ','.join(str(x['t']) for x in r['IQ']) or '-'
    P("%-9s %-5d %-3d %-8s %d %-8s %-16s %-19s %s  %s  %-5s %s" % (
        '%s/%d' % (r['TI'], r['j']), r['order'], len(r['Phi_all_G_orbit_reps']), ''.join('%x' % x for x in r['Phi']), r['d'], iq, ks, tiers, r['n0'], r['deg0'], r['classification'], lu))
