# summarize_g8.py -- tables from trans_g8.json (pure Python).  usage: python summarize_g8.py [trans_g8.json] [out.txt]
import json, sys, collections, itertools, os
from verify_enum import rank
fn = sys.argv[1] if len(sys.argv) > 1 else 'trans_g8.json'
outfn = sys.argv[2] if len(sys.argv) > 2 else 'table_g8.txt'
D = json.load(open(fn))
out = open(outfn, 'w')
def P(*a):
    print(*a); print(*a, file=out)


def classify(w):
    if w['Dlat_IQ'] == 8:
        return '(i)'
    if w['Dlat_all'] == 8:
        return '(ii)'
    if w['Dlat_all'] == 16:
        return '(iii)'
    return '(iv)'


def mindeg(w, D):
    for cname, deg in (('IQ', 2), ('le4', 4), ('le8', 8), ('all', 16)):
        if w['Dlat_' + cname] is not None and w['Dlat_' + cname] <= D:
            return deg
    return None


def subsig(c):
    parts = []
    for s in sorted(c['CMsub'], key=lambda s: (s['deg'], s['sig'])):
        rs = sorted(r for B, r in s['sig'])
        parts.append("%d:%s%s" % (s['deg'], ''.join(str(r) for r in rs), '*' if s['balanced'] else ''))
    return ' '.join(parts)


def deg0(LU):
    # minimal cohomological degree of a non-Lefschetz Hodge class on A itself: min sum|v_i| over 0 != v in Lambda_U, |v_i|<=1
    d = len(LU)
    best = None
    for v in itertools.product((-1, 0, 1), repeat=len(LU[0])):
        if any(v) and rank(LU + [list(v)]) == d:
            w = sum(map(abs, v))
            if best is None or w < best:
                best = w
    return best


ST = json.load(open('struct_g8.json')) if os.path.exists('struct_g8.json') else {}
def st(G):
    x = ST.get('%s/%d' % (G['TI'], G['j']))
    return '-' if x is None else x['struct'].replace(' ', '')
tot = collections.Counter()
rows = []
for G in D['groups']:
    cs = G['cases']
    tot['(G,rho) pairs'] += 1
    tot['type G-orbits'] += len(cs)
    tot['type N-classes'] += len(set(c['Nclass'] for c in cs))
    tot['primitive G-orbits'] += sum(c['primitive'] for c in cs)
    tot['primitive N-classes'] += len(set(c['Nclass'] for c in cs if c['primitive']))
    deg = [c for c in cs if c['primitive'] and c['d'] >= 1]
    tot['prim. degenerate G-orbits'] += len(deg)
    seen = {}
    for c in deg:
        seen.setdefault(c['Nclass'], []).append(c)
    tot['prim. degenerate N-classes'] += len(seen)
    for ncl, lst in seen.items():
        c = lst[0]
        rows.append((G, c, len(lst)))
P("COUNTS:", dict(tot))
byd = collections.Counter(c['d'] for G, c, m in rows)
byd_orb = collections.Counter()
for G, c, m in rows:
    byd_orb[c['d']] += m
P("prim. degenerate by d (N-classes):", dict(sorted(byd.items())), " (G-orbits):", dict(sorted(byd_orb.items())))
cl = collections.Counter((c['d'], classify(c['weil']), mindeg(c['weil'], 8), mindeg(c['weil'], 16), c['weil']['Dq_IQ'], c['weil']['Dq_all']) for G, c, m in rows)
P("classification (d, class, min field degree for lattice at D=8, at D<=16, Dq_IQ, Dq_all): #N-classes")
for k, v in sorted(cl.items(), key=str):
    P("   ", k, v)
P("n0 (index of degeneracy) distribution:", dict(collections.Counter(c['n0'] for G, c, m in rows)))
P("")
P("LEGEND: one row per N_W(G)-class (W(B8)-class of the pair (G,Phi)); #orb = number of G-orbits of types in the class.")
P("  Phi: 0-based points, rho(i)=i+8.  LU: HNF basis of Lambda_U in coordinates chi -> (chi_i - chi_{i+8})_{i<8}.")
P("  IQ t: t = 2|Phi cap B_K| - 8 for each imaginary quadratic subfield K (0 = balanced (4,4)).")
P("  CM subfields: 'deg:r1r2..' = degree 2k of K' and, for each of the k rho-pairs of blocks, r = |Phi cap B| (K'-signature (r, b-r));")
P("     '*' = balanced at every embedding (a Weil structure for K' on A itself).")
P("  tiers: Dlat_X / Dq_X = smallest Weil-variety dimension whose Weil characters (fields of index category X) generate")
P("     Lambda_U as a lattice / span it over Q; X = IQ (quadratic), 4 (deg<=4), 8 (deg<=8), all (deg<=16); '-' = none up to 16.")
P("  n0: smallest n such that A^n carries a Hodge class outside the Lefschetz algebra.")
P("")
hdr = "%-8s %-2s %-6s %-4s %-26s %-2s %-40s %-10s %-34s %-27s %-3s %-4s %-5s" % (
    "TI", "j", "|G|", "#orb", "Phi", "d", "LU", "IQ t", "CM subfields", "Dlat/Dq IQ|4|8|all", "n0", "deg0", "class") + " structure(G)"
P(hdr)
rows.sort(key=lambda r: (r[1]['d'], classify(r[1]['weil']), r[0]['order'], r[0]['k'], r[0]['j'], r[1]['Phi']))
f = lambda x: '-' if x is None else str(x)
for G, c, m in rows:
    w = c['weil']
    tiers = "%s/%s|%s/%s|%s/%s|%s/%s" % (f(w['Dlat_IQ']), f(w['Dq_IQ']), f(w['Dlat_le4']), f(w['Dq_le4']),
                                         f(w['Dlat_le8']), f(w['Dq_le8']), f(w['Dlat_all']), f(w['Dq_all']))
    iq = ','.join(str(x['t']) for x in c['IQ']) or '-'
    P("%-8s %-2d %-6d %-4d %-26s %-2d %-40s %-10s %-34s %-27s %-3s %-4s %-5s" % (
        G['TI'], G['j'], G['order'], m, ','.join(map(str, c['Phi'])), c['d'],
        ';'.join(''.join('%+d' % x for x in r) for r in c['LU']), iq, subsig(c), tiers, f(c['n0']), f(deg0(c['LU'])), classify(w)) + " " + st(G))
P("")
if ST:
    ab = [(G, c) for G, c, m in rows if ST['%s/%d' % (G['TI'], G['j'])]['abelian']]
    P("abelian Galois groups: %d N-classes ; classes %s ; n0 %s ; groups %s" % (len(ab),
      dict(collections.Counter(classify(c['weil']) for G, c in ab)), dict(collections.Counter(c['n0'] for G, c in ab)),
      sorted(set('%s=%s' % (G['TI'], st(G)) for G, c in ab))))
P("deg0 distribution by d:", dict(collections.Counter((c['d'], deg0(c['LU'])) for G, c, m in rows)))
# extended Weil search for the d>=3 classes (ext_weil.py), if available
ext = {}
for fn in ('ext_weil_D32.json', 'ext_weil_D48_sel.json'):
    if os.path.exists(fn):
        for r in json.load(open(fn)):
            ext[(r['TI'], r['j'], tuple(r['Phi']))] = (fn, r['levels'])
if ext:
    P("")
    P("EXTENDED WEIL SEARCH beyond D=16 (ext_weil.py; per level: D, rank, saturation index of the lattice generated by all")
    P("Weil characters of dimension <= D, field indices [G:M] = [K':Q] used so far):")
    for G, c, m in rows:
        if c['d'] < 3:
            continue
        key = (G['TI'], G['j'], tuple(c['Phi']))
        if key in ext:
            fn, lv = ext[key]
            P("  %-8s j=%d |G|=%-4d Phi=%-26s %s  [%s]" % (G['TI'], G['j'], G['order'], ','.join(map(str, c['Phi'])),
              ' ; '.join('D=%d rank %d sat %s idx %s' % (l['D'], l['rank'], l['satindex'], l['indices']) for l in lv), fn))
