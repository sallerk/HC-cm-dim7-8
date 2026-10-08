# route2_match_g7.py -- class-by-class match of the subgroup list used by the g = 7 enumeration with a second GAP
# route.  Route 1 (stageC.py prep): ConjugacyClassesSubgroups of W(B7)/<rho>, lifted; the groups the enumeration
# uses are those of enum_g7.json.  Route 2: ConjugacyClassesSubgroups of W(B7) itself, keeping the classes that
# contain rho (rho is central, so this is a property of the class).  Each group of enum_g7.json is compared, with
# IsConjugate in W(B7), against the route-2 representatives with the same invariants (order, orbit lengths on the
# 14 points, numbers of elements of each cycle type, abelian invariants, orders of the centre and the derived
# subgroup) that are not yet matched.  The two lists are the same list of classes iff every group of enum_g7.json
# is matched and the numbers of classes are equal.
#   usage (sagemath container):  sage -python route2_match_g7.py > log_route2_match_g7.txt
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
    ct = collections.Counter()
    for c in libgap.ConjugacyClasses(H):
        ct[tuple(sorted(int(x) for x in libgap.CycleLengths(libgap.Representative(c), pts)))] += int(libgap.Size(c))
    return (int(libgap.Size(H)), orbs, tuple(sorted(ct.items())),
            tuple(int(x) for x in libgap.AbelianInvariants(H)),
            int(libgap.Size(libgap.Centre(H))), int(libgap.Size(libgap.DerivedSubgroup(H))))


t0 = time.time()
print("g", g, "|W|", libgap.Size(W), flush=True)
cc = libgap.ConjugacyClassesSubgroups(W)
reps = [H for H in (libgap.Representative(c) for c in cc) if libgap.IsSubgroup(H, R)]
print("route 2: classes of subgroups of W(B7) %d, containing rho %d, time %.1fs"
      % (len(cc), len(reps), time.time() - t0), flush=True)
bucket = collections.defaultdict(list)
for j, H in enumerate(reps):
    bucket[inv(H)].append(j)
print("route 2 invariants: %d distinct values, largest bucket %d, time %.1fs"
      % (len(bucket), max(len(v) for v in bucket.values()), time.time() - t0), flush=True)

stored = json.load(open('enum_g7.json'))['groups']
matched = {}
problems = []
tests = 0
for i, r in enumerate(stored):
    H = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in r['gens']])
    ok_in = bool(libgap.IsSubgroup(W, H)) and bool(libgap.IsSubgroup(H, R)) and int(libgap.Size(H)) == r['order']
    if not ok_in:
        problems.append(('not in W, no rho, or wrong order', r['gid']))
    cand = bucket.get(inv(H), [])
    hit = None
    for j in cand:
        tests += 1
        if libgap.IsConjugate(W, H, reps[j]):
            hit = j
            break
    if hit is None:
        problems.append(('no unmatched conjugate class', r['gid']))
    else:
        cand.remove(hit)
        matched[r['gid']] = hit
    if (i + 1) % 500 == 0:
        print("  %d/%d stored groups done, %d IsConjugate tests, time %.1fs"
              % (i + 1, len(stored), tests, time.time() - t0), flush=True)
print("stored groups (enum_g7.json): %d; matched to distinct route-2 classes: %d; IsConjugate tests: %d"
      % (len(stored), len(matched), tests))
print("problems:", len(problems), problems[:20])
same = len(problems) == 0 and len(matched) == len(stored) == len(reps)
print("RESULT: the two class lists %s" % ("are the same (bijection by conjugacy in W(B7))" if same else "DIFFER"))
print("total time %.1fs" % (time.time() - t0))
