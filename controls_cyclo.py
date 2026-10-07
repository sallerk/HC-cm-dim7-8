# controls_cyclo.py -- controls (Sage+GAP) for degree-16 cyclotomic CM fields Q(zeta_m), phi(m)=16.
#  G = (Z/m)^x acting regularly (sigma_a -> sigma_{ha}); points: u_0..u_7 = the units < m/2, point i+8 = -u_i.
#  (a) all CM types of Q(zeta_m): d, primitivity, tiers (analyse_trans);
#  (b) named types: J(y^2 = x^17 + 1) (Phi = {1..8} in (Z/17)^x), and every Fermat-curve factor type
#      H_{a,b,c} = { t : <ta/m> + <tb/m> + <tc/m> = 1 },  a+b+c = 0 mod m, abc != 0, gcd(a,b,c,m) = 1;
#  (c) location of each cyclotomic (G,Phi) inside the transitive enumeration trans_g8.json (W(B8)-conjugacy),
#      comparing invariants computed twice (direct vs stored).
import sys, os, json, itertools, collections
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sage.all import libgap
from cmenum import rho_perm
from trans_enum import analyse_trans

g = 8
n = 16
out = open(HERE + '/controls_cyclo_out.txt', 'w')
def P(*a):
    print(*a, flush=True); print(*a, file=out); out.flush()

St = json.load(open(HERE + '/trans_g8.json'))
rho = rho_perm(g)
W = libgap.Centralizer(libgap.SymmetricGroup(n), rho)
stored = {}
for Gr in St['groups']:
    stored.setdefault(Gr['TI'], []).append(Gr)

def gapgroup(gens):
    return libgap.Group([libgap.PermList([x + 1 for x in p]) for p in gens])

def inv_of(c):
    w = c.get('weil')
    return dict(d=c['d'], dimMT=c['dimMT'], primitive=c['primitive'], n0=c.get('n0'),
                IQ=sorted(abs(x['t']) for x in c.get('IQ', [])),
                CMsub=sorted((s['deg'], tuple(sorted(min(r, s['b'] - r) for B, r in s['sig']))) for s in c.get('CMsub', [])),
                tiers=None if w is None else tuple(w[k] for k in ('Dlat_IQ', 'Dq_IQ', 'Dlat_le4', 'Dq_le4', 'Dlat_le8', 'Dq_le8', 'Dlat_all', 'Dq_all')))

summary = []
for m in [17, 32, 34, 40, 48, 60]:
    units = [a for a in range(1, m) if gcd(a, m) == 1]
    half = [a for a in units if a < m / 2]
    assert len(half) == 8
    lab = {}
    for i, a in enumerate(half):
        lab[a] = i
        lab[(m - a) % m] = i + 8
    gens = []
    for h in units:
        gens.append([lab[(h * a) % m] for a in sorted(lab, key=lambda a: lab[a])])
    G = gapgroup(gens)
    assert int(libgap.Size(G)) == 16 and bool(libgap.IsSubset(G, [rho]))
    TI = "16T%d" % int(libgap.TransitiveIdentification(G))
    gi, cases = analyse_trans(G, g)
    P("=== m=%d  Gal = %s  %s ; type G-orbits %d (N-classes %d) ; primitive %d ; primitive degenerate %d (N-classes %d) ; d-dist(primitive) %s ; d-dist(all) %s" % (
        m, libgap.StructureDescription(G), TI, len(cases), len(set(c['Nclass'] for c in cases)),
        sum(c['primitive'] for c in cases), sum(1 for c in cases if c['primitive'] and c['d'] >= 1),
        len(set(c['Nclass'] for c in cases if c['primitive'] and c['d'] >= 1)),
        dict(collections.Counter(c['d'] for c in cases if c['primitive'])), dict(collections.Counter(c['d'] for c in cases))))
    # location in the stored enumeration
    cand = stored.get(TI, [])
    hit = None
    for Gr in cand:
        Gs = gapgroup(Gr['gens'])
        x = libgap.RepresentativeAction(W, G, Gs)
        if x != libgap.eval('fail'):
            assert hit is None
            hit = (Gr, Gs, x)
    if hit is None:
        P("   NOT LOCATED in trans_g8.json (candidates %d)" % len(cand))
        continue
    Gr, Gs, x = hit
    P("   located: %s j=%d (stored |G|=%d)" % (Gr['TI'], Gr['j'], Gr['order']))
    sc = {tuple(c['Phi']): c for c in Gr['cases']}
    nmis = 0
    for c in cases:
        img = libgap.OnSets(libgap.Set([p + 1 for p in c['Phi']]), x)
        orb = libgap.Orbit(Gs, img, libgap.OnSets)
        mn = min(tuple(sorted(int(y) - 1 for y in s)) for s in orb)
        if inv_of(c) != inv_of(sc[mn]):
            nmis += 1
            P("   MISMATCH", c['Phi'], inv_of(c), inv_of(sc[mn]))
    P("   direct-vs-located invariant mismatches: %d of %d type orbits" % (nmis, len(cases)))
    def find(Phi_units):
        S = tuple(sorted(lab[a % m] for a in Phi_units))
        GS =libgap.Orbit(G, libgap.Set([p + 1 for p in S]), libgap.OnSets)
        mn = min(tuple(sorted(int(y) - 1 for y in s)) for s in GS)
        return next(c for c in cases if tuple(c['Phi']) == mn)
    def show(tag, c):
        w = c.get('weil')
        P("   %-34s Phi(min rep)=%s d=%d dimMT=%d primitive=%s n0=%s IQ t=%s CMsub=%s tiers=%s" % (
            tag, c['Phi'], c['d'], c['dimMT'], c['primitive'], c.get('n0'), [x['t'] for x in c.get('IQ', [])],
            inv_of(c)['CMsub'], inv_of(c)['tiers']))
    if m in (17, 34):
        base = 17
        Phi = [a for a in range(1, 9)]
        if m == 34:
            Phi = [a if a % 2 else a + 17 for a in Phi]   # zeta_34 = -zeta_17^9: use the odd representatives
        show("y^2=x^17+1, Phi={1..8}", find(Phi))
    # Fermat factor types
    seen = set()
    ferm = []
    for a in range(1, m):
        for b in range(1, m):
            cc = (-a - b) % m
            if cc == 0 or gcd(gcd(gcd(a, b), cc), m) != 1:
                continue
            Phi = [t for t in units if ((t * a) % m + (t * b) % m + (t * cc) % m) == m]
            assert len(Phi) == 8, (a, b, cc, Phi)
            c = find(Phi)
            key = (tuple(c['Phi']),)
            if key in seen:
                continue
            seen.add(key)
            ferm.append(((a, b, cc), c))
    P("   Fermat factor types H_{a,b,c} (one per G-orbit of types): %d ; d-dist %s ; primitive %d" % (
        len(ferm), dict(collections.Counter(c['d'] for _, c in ferm)), sum(c['primitive'] for _, c in ferm)))
    for abc, c in ferm:
        if c['d'] >= 1:
            show("Fermat (a,b,c)=%s" % (abc,), c)
    summary.append((m, TI, len(cases), sum(1 for c in cases if c['primitive'] and c['d'] >= 1)))
P("SUMMARY (m, TI, type orbits, primitive degenerate orbits):", summary)
