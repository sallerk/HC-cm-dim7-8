# controls_lit.py -- literature controls (Sage+GAP); the sources are cited below and in README.txt.
#  (1) Ribet 1980 (3.12) (Lenstra): Q(zeta32), S={1,7,13,15,21,23,27,29} degenerate; S'={1,7,9,11,13,15,27,29}: both odd
#      characters of order 2 vanish  -> expect d=1 and d=2, both primitive.
#  (2) Yanai JTNB 27 (2015) sec. 5, p=17: G = Z/2 x (Z/17)^x/{+-1}, S={(0,1),(1,2),...,(1,8)}: Weil type wrt Q(i), index 1.
#  (3) Kida, Moscow J. Comb. Number Theory 8 (2019), Ex. 6.4 (k=3): primitive degenerate (rank 5) CM types:
#      SD16: 64, M16: 32, D16 (dihedral of order 16): none.   Dodson J. Algebra 1987 (via zbMATH): min rank B(8)=5.
#  (4) Kida-Yanai IJNT 16 (2020) (abstract): index of degeneracy 2 occurs for C3 x D4 and SmallGroup(24,8) (degree 24).
import sys, os, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sage.all import libgap
from cmenum import rho_perm
from trans_enum import analyse_trans

out = open(HERE + '/controls_lit_out.txt', 'w')
def P(*a):
    print(*a, flush=True); print(*a, file=out); out.flush()


def regular(G, g):
    """regular action of the abstract group G on 2g points; returns (GAP perm group conjugated so that the unique
    central fpf involution is the standard rho, element list, map element->point)."""
    n = 2 * g
    els = list(libgap.Elements(G))
    act = libgap.Action(G, els, libgap.OnLeftInverse)   # x -> g x  (left multiplication)
    hom = libgap.ActionHomomorphism(G, els, libgap.OnLeftInverse)
    Z = libgap.Centre(act)
    invs = [z for z in libgap.Elements(Z) if int(libgap.Order(z)) == 2 and int(libgap.NrMovedPoints(z)) == n]
    return act, els, invs


def conj_to_rho(act, z, g):
    n = 2 * g
    pi = libgap.RepresentativeAction(libgap.SymmetricGroup(n), z, rho_perm(g))
    return libgap.ConjugateGroup(act, pi), pi


def find_case(cases, Gc, Phi):
    orb = libgap.Orbit(Gc, libgap.Set([p + 1 for p in Phi]), libgap.OnSets)
    mn = min(tuple(sorted(int(y) - 1 for y in s)) for s in orb)
    return next(c for c in cases if tuple(c['Phi']) == mn)


def show(tag, c):
    w = c.get('weil')
    P("   %-40s d=%d dimMT=%d primitive=%s n0=%s IQ |t|=%s tiers(Dlat/Dq IQ,4,8,all)=%s" % (
        tag, c['d'], c['dimMT'], c['primitive'], c.get('n0'), sorted(abs(x['t']) for x in c.get('IQ', [])),
        None if w is None else (w['Dlat_IQ'], w['Dq_IQ'], w['Dlat_le4'], w['Dq_le4'], w['Dlat_le8'], w['Dq_le8'], w['Dlat_all'], w['Dq_all'])))


# (1) Q(zeta32): G=(Z/32)^x acting regularly, points labelled as in controls_cyclo.py
m, g = 32, 8
units = [a for a in range(1, m) if a % 2]
half = [a for a in units if a < m / 2]
lab = {}
for i, a in enumerate(half):
    lab[a] = i; lab[m - a] = i + 8
gens = [[lab[(h * a) % m] for a in sorted(lab, key=lambda a: lab[a])] for h in units]
G = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in gens])
gi, cases = analyse_trans(G, g)
P("(1) Ribet 1980 (3.12), Q(zeta32), %s:" % libgap.StructureDescription(G))
for tag, S in [("S  = {1,7,13,15,21,23,27,29}", [1, 7, 13, 15, 21, 23, 27, 29]), ("S' = {1,7,9,11,13,15,27,29}", [1, 7, 9, 11, 13, 15, 27, 29])]:
    show(tag, find_case(cases, G, [lab[a] for a in S]))

# (2) Yanai sec. 5, p = 17: G = Z/2 x (Z/17)^x/{+-1}; element (e, j), j in 1..8 representing +-j
P("(2) Yanai JTNB 2015 sec. 5 (Tautz-Top-Verberkmoes curve, p=17):")
elsY = [(e, j) for e in (0, 1) for j in range(1, 9)]
def red(x):
    x %= 17
    return min(x, 17 - x)
def mulY(a, b):
    return ((a[0] + b[0]) % 2, red(a[1] * b[1]))
# labelling: point i = (0, i+1) for i<8, point i+8 = rho*(0,i+1) = (1, i+1); rho = (1,1)
labY = {(0, j): j - 1 for j in range(1, 9)}
labY.update({(1, j): j + 7 for j in range(1, 9)})
gensY = [[labY[mulY(h, x)] for x in sorted(labY, key=lambda x: labY[x])] for h in elsY]
GY = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in gensY])
assert int(libgap.Size(GY)) == 16 and bool(libgap.IsSubset(GY, [rho_perm(8)]))
giY, casesY = analyse_trans(GY, 8)
SY = [(0, 1), (1, 2), (0, 3), (1, 4), (0, 5), (1, 6), (0, 7), (1, 8)]
P("   G = %s (16T%d)" % (libgap.StructureDescription(GY), int(libgap.TransitiveIdentification(GY))))
show("S = {(0,1),(1,2),(0,3),...,(1,8)}", find_case(casesY, GY, [labY[x] for x in SY]))

# (3) Kida Ex. 6.4 (k=3): SD16, M16, D16 regular; count CM TYPES (not orbits) by (primitive, d)
P("(3) Kida 2019 Ex. 6.4 / Dodson B(8)=5: regular Galois CM fields of degree 16:")
for name, sid in [("SD16", libgap.SmallGroup(16, 8)), ("M16", libgap.SmallGroup(16, 6)), ("D16", libgap.SmallGroup(16, 7))]:
    act, els, invs = regular(sid, 8)
    assert len(invs) == 1
    Gc, pi = conj_to_rho(act, invs[0], 8)
    gic, cs = analyse_trans(Gc, 8, weil=False)
    cnt = collections.Counter()
    for c in cs:
        cnt[(c['primitive'], c['d'])] += c['norb']
    P("   %-5s = %-12s 16T%-3d CM types by (primitive, d): %s ; primitive degenerate types: %d" % (
        name, libgap.StructureDescription(sid), int(libgap.TransitiveIdentification(Gc)), dict(sorted(cnt.items())),
        sum(v for (p, d), v in cnt.items() if p and d >= 1)))

# (4) Kida-Yanai: index of degeneracy 2 in degree 24
P("(4) Kida-Yanai IJNT 2020: degree-24 Galois CM fields, index of degeneracy (n0) of primitive degenerate types:")
for name, sid in [("C3 x D4", libgap.DirectProduct(libgap.CyclicGroup(3), libgap.DihedralGroup(8))), ("SmallGroup(24,8)", libgap.SmallGroup(24, 8))]:
    act, els, invs = regular(sid, 12)
    for z in invs:
        Gc, pi = conj_to_rho(act, z, 12)
        gic, cs = analyse_trans(Gc, 12, want_N=False, weil=False)
        cnt = collections.Counter((c['d'], c['n0']) for c in cs if c['primitive'] and c['d'] >= 1)
        P("   %-16s (%s) central fpf involutions %d: primitive degenerate G-orbits by (d, n0): %s" % (
            name, libgap.StructureDescription(sid), len(invs), dict(sorted(cnt.items()))))
