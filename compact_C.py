# compact_C.py -- compact list of Stage C reduced configurations with minimal Weil dimension > 6 (one line per N-class).
import json, collections
D = json.load(open('enum_g7.json')); groups = {G['gid']: G for G in D['groups']}
ext = {(r['gid'], tuple(r['Phi'])): r['levels'] for r in json.load(open('ext_weil_g7_D24.json'))}
part = lambda c: tuple(sorted((len(o) // 2 for o in groups[c['gid']]['orbits']), reverse=True))
big = [c for c in D['cases'] if 'general' in c and (c['general']['Dmin'] is None or c['general']['Dmin'] > 6)]
cls = collections.OrderedDict()
for c in sorted(big, key=lambda c: (part(c), c['d'], c['general']['Dmin'] or 99, groups[c['gid']]['order'], c['gid'], c['Phi'])):
    cls.setdefault((c['gid'], c['Nclass']), []).append(c)
out = open('list_C_compact.txt', 'w')
def P(s):
    print(s); out.write(s + '\n')
P("dims        gid   |G|   orbit TIs              Phi(hex)  d Dmin Dq  IQmult LambdaU (signed digits)")
for (gid, ncl), lst in cls.items():
    c = lst[0]; G = groups[gid]
    Dmin = c['general']['Dmin']; Dq = c['general']['Dq']
    if Dmin is None and (gid, tuple(c['Phi'])) in ext:
        lv = ext[(gid, tuple(c['Phi']))]
        Dmin = next((l['D'] for l in lv if l['satindex'] == 1), None); Dq = next((l['D'] for l in lv if l['rank'] == c['d']), None)
        Dmin = '%s*' % Dmin; Dq = '%s*' % Dq
    P("%-11s %-5d %-5d %-22s %-9s %d %-4s %-3s %-6s %s" % (','.join(map(str, part(c))), gid, G['order'], '+'.join(G['orbit_TI']),
      ''.join('%x' % x for x in c['Phi']), c['d'], Dmin, Dq, c['mult']['Dmin'] if c['mult']['Dmin'] is not None else '-',
      ';'.join(''.join('%+d' % x for x in r) for r in c['LU'])))
P("(* = from the extended search ext_weil.py Dmax=24; the stage-C tier search itself stopped at D=16)")
