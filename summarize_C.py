# summarize_C.py -- Stage C summary (g=7, all subgroups of W(B7) containing rho) from enum_g7.json.  pure Python.
import json, collections
D = json.load(open('enum_g7.json'))
out = open('summary_C.txt', 'w')
def P(*a):
    print(*a); print(*a, file=out)
groups = {G['gid']: G for G in D['groups']}
cs = D['cases']
P("g=7: W(B7)-classes of subgroups containing rho: %d ; cases (G-orbits of types) %d ; N-classes %d" % (
    len(groups), len(cs), len(set((c['gid'], c['Nclass']) for c in cs))))
red = [c for c in cs if c['reduced']]
P("reduced cases %d (N-classes %d) ; d-dist %s" % (len(red), len(set((c['gid'], c['Nclass']) for c in red)),
  dict(sorted(collections.Counter(c['d'] for c in red).items()))))
P("transitive reduced: %d, d-dist %s" % (sum(1 for c in red if groups[c['gid']]['transitive']),
  dict(collections.Counter(c['d'] for c in red if groups[c['gid']]['transitive']))))
P("strict KEY TEST on reduced: pass %d / fail %d ; non-reduced strict fails %d" % (
    sum(c['passed'] for c in red), sum(not c['passed'] for c in red), sum(1 for c in cs if not c['reduced'] and not c['passed'])))
def part(c):
    return tuple(sorted((len(o) // 2 for o in groups[c['gid']]['orbits']), reverse=True))
tier = collections.Counter()
tierN = collections.defaultdict(set)
fails = []
for c in red:
    if c['d'] == 0:
        continue
    if c['passed']:
        t = 'strict'
    else:
        gm = c['general']['Dmin']
        im = c['mult']['Dmin']
        t = 'IQmult D=%s' % im if (im is not None and im <= 6) else 'general D=%s (IQmult %s, Dq %s)' % (gm, im, c['general']['Dq'])
        fails.append(c)
    key = (part(c), c['d'], t)
    tier[key] += 1
    tierN[key].add((c['gid'], c['Nclass']))
P("")
P("reduced d>=1 by (factor dims, d, tier): #cases / #N-classes")
for k in sorted(tier, key=str):
    P("   %-22s d=%d  %-40s %4d / %d" % (k[0], k[1], k[2], tier[k], len(tierN[k])))
P("")
P("consistency: general Dmin <= IQmult Dmin whenever both exist: %s" % all(
    c['general']['Dmin'] is None or c['mult']['Dmin'] is None or c['general']['Dmin'] <= c['mult']['Dmin'] for c in fails))
big = [c for c in fails if c['general']['Dmin'] is None or c['general']['Dmin'] > 6]
P("REDUCED CONFIGURATIONS WITH MINIMAL WEIL DIMENSION > 6 (one line per N-class; generators in enum_g7.json):")
seen = set()
for c in sorted(big, key=lambda c: (part(c), c['d'], c['general']['Dmin'] or 99, c['gid'])):
    if (c['gid'], c['Nclass']) in seen:
        continue
    seen.add((c['gid'], c['Nclass']))
    G = groups[c['gid']]
    norb = sum(1 for x in big if (x['gid'], x['Nclass']) == (c['gid'], c['Nclass']))
    P("  gid %-5d %-22s |G|=%-6d TI=%-28s Phi=%-22s d=%d LU=%s IQmultDmin=%s IQmultQrank=%s generalDmin=%s Dq=%s Dlat_IQ/4/8=%s/%s/%s  (#orb %d)" % (
        c['gid'], part(c), G['order'], '+'.join(G['orbit_TI']), c['Phi'], c['d'], c['LU'], c['mult']['Dmin'], c['mult']['qrank'],
        c['general']['Dmin'], c['general']['Dq'], c['general']['Dlat_IQ'], c['general']['Dlat_le4'], c['general']['Dlat_le8'], norb))
P("total: %d cases / %d N-classes with minimal Weil dimension > 6" % (len(big), len(seen)))
P("by (factor dims, d, general Dmin):", dict(collections.Counter((part(c), c['d'], c['general']['Dmin']) for c in big)))
