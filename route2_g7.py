# route2_g7.py -- second GAP route for the Stage C class list (g = 7).
# Stage C ("stageC.py prep") took ConjugacyClassesSubgroups of the quotient W(B7)/<rho> and lifted the classes.
# Here GAP computes ConjugacyClassesSubgroups of W(B7) itself and keeps the classes containing rho (rho is central,
# so this is a property of the class).  The two lists are compared through the multiset of the invariants
# (order, orbit lengths on the 14 points, order of the centre, order of the derived subgroup).
#   usage (sagemath container):  sage -python route2_g7.py > log_route2_g7.txt
import json, time, collections
from sage.all import libgap

g = 7
n = 2 * g
rho = libgap.PermList([((i + g) % n) + 1 for i in range(n)])
W = libgap.Centralizer(libgap.SymmetricGroup(n), rho)
R = libgap.Group([rho])
pts = libgap.eval('[1..%d]' % n)


def inv(H):
    orbs = tuple(sorted(int(libgap.Size(o)) for o in libgap.Orbits(H, pts)))
    return (int(libgap.Size(H)), orbs, int(libgap.Size(libgap.Centre(H))),
            int(libgap.Size(libgap.DerivedSubgroup(H))))


print("g", g, "|W|", libgap.Size(W), flush=True)
t = time.time()
cc = libgap.ConjugacyClassesSubgroups(W)
print("W(B7): classes of subgroups", len(cc), "time %.1fs" % (time.time() - t), flush=True)
keep = [H for H in (libgap.Representative(c) for c in cc) if libgap.IsSubgroup(H, R)]
print("classes containing rho", len(keep), flush=True)
A = collections.Counter(inv(H) for H in keep)

stored = json.load(open('stageC_classes.json'))
B = collections.Counter()
bad_order = 0
for r in stored:
    H = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in r['gens']])
    v = inv(H)
    bad_order += (v[0] != r['order'])
    B[v] += 1
print("stored classes", len(stored), "| stored orders wrong:", bad_order, flush=True)
print("distinct invariants: route 2", len(A), "| stored", len(B))
print("invariant multisets equal:", A == B)
d = list(((A - B) + (B - A)).items())
print("differences (first 20):", d[:20])
print("total time %.1fs" % (time.time() - t))
