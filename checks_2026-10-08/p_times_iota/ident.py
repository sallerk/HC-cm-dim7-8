# Part 3: identify P = image of G on the conjugation pairs (GAP TransitiveIdentification).
import json, collections, time
from sage.all import libgap
T0 = time.time()
import os
CL = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")) + "/"


def tid(gens, n):
    G = libgap.Group([libgap.PermList([x + 1 for x in g]) for g in gens] or [libgap.PermList([])])
    return int(libgap.TransitiveIdentification(G)), int(libgap.Size(G)), str(libgap.StructureDescription(G))


def closure(gens, n):
    idt = tuple(range(n)); seen = {idt}; fr = [idt]
    gens = [tuple(g) for g in gens]
    while fr:
        nx = []
        for a in fr:
            for g in gens:
                b = tuple(g[a[i]] for i in range(n))
                if b not in seen:
                    seen.add(b); nx.append(b)
        fr = nx
    return seen


def orbits(gens, n):
    seen = set(); res = []
    for s in range(n):
        if s in seen:
            continue
        o = {s}; fr = [s]
        while fr:
            x = fr.pop()
            for g in gens:
                if g[x] not in o:
                    o.add(g[x]); fr.append(g[x])
        seen |= o; res.append(sorted(o))
    return res


out = {}
# ---------------- g = 8 ----------------
D = json.load(open(CL + "deg_classes_g8.json"))
g8 = []
for n, x in enumerate(D):
    gens = x["gens"]
    P8 = [[g[j] % 8 for j in range(8)] for g in gens]
    k, o, sd = tid(P8, 8)
    g8.append(dict(n=n, TI=x["TI"], j=x["j"], d=x["d"], order=x["order"], Pk=k, Porder=o, ratio=x["order"] // o,
                   norb=len(x["Phi_all_G_orbit_reps"])))
out["g8"] = g8
for d in (1, 2, 4):
    c = collections.Counter(r["Pk"] for r in g8 if r["d"] == d)
    r2 = collections.Counter(r["ratio"] for r in g8 if r["d"] == d)
    print("g8 d=%d: classes %d; P=8Tk counts %s; |G|/|P| counts %s" % (d, sum(c.values()), dict(sorted(c.items())), dict(r2)))
print("g8 d=1 per class:", [(r["n"], r["TI"], r["j"], r["order"], "8T%d" % r["Pk"], r["norb"]) for r in g8 if r["d"] == 1])

# ---------------- g = 7 ----------------
E = json.load(open(CL + "enum_g7.json"))["groups"]
rows = []
for line in open(CL + "list_C_compact.txt"):
    f = line.split()
    if f and f[0] in ("6,1", "5,1,1"):
        rows.append(dict(dims=f[0], gid=int(f[1]), G=int(f[2]), TIs=f[3], Phi=[int(c, 16) for c in f[4]],
                         d=int(f[5]), Dmin=f[6], Dq=f[7], IQmult=f[8], LU=f[9]))
g7 = []
for r in rows:
    grp = E[r["gid"]]
    assert grp["gid"] == r["gid"] and grp["order"] == r["G"]
    gens = grp["gens"]
    els = closure(gens, 14)
    assert len(els) == r["G"]
    orbs = sorted(orbits(gens, 14), key=len, reverse=True)
    big = orbs[0]; small = orbs[1:]
    assert len(big) == {"6,1": 12, "5,1,1": 10}[r["dims"]] and all(len(s) == 2 for s in small)
    m = len(big) // 2
    reps = sorted(i for i in big if i < 7)
    pidx = {}
    for a, i in enumerate(reps):
        pidx[i] = a; pidx[i + 7] = a
    Pg = [[pidx[g[i]] for i in reps] for g in gens]
    k, Po, Psd = tid(Pg, m)
    bidx = {p: a for a, p in enumerate(big)}
    k12, o12, _ = tid([[bidx[g[p]] for p in big] for g in gens], len(big))
    # characters of the size-2 orbits; image on pairs per element
    def pairimg(h):
        return tuple(pidx[h[i]] for i in reps)
    def chi(h, s):
        return int(h[s[0]] != s[0])
    rec = dict(dims=r["dims"], gid=r["gid"], G=r["G"], TIs=r["TIs"], d=r["d"], Dmin=r["Dmin"],
               Phi="".join("%x" % c for c in r["Phi"]), P="%dT%d" % (m, k), Porder=Po, Pstruct=Psd,
               big_TI="%dT%d" % (len(big), k12), big_order=o12, G_over_P=r["G"] // Po)
    # rho-swapped 2-block systems of the big orbit = index-2 subgroups of the image avoiding rho; check
    # whether some character of the big-orbit image equals chi of a small orbit
    imgs = collections.defaultdict(set)
    for h in els:
        key = tuple(h[p] for p in big)
        imgs[key].add(tuple(chi(h, s) for s in small))
    rec["small_chars_determined_by_big_orbit"] = all(len(v) == 1 for v in imgs.values())
    if r["dims"] == "5,1,1":
        sa, sb = small
        # chi_b determined by (big orbit, chi_a)?
        dd = collections.defaultdict(set)
        for h in els:
            dd[(tuple(h[p] for p in big), chi(h, sa))].add(chi(h, sb))
        dd2 = collections.defaultdict(set)
        for h in els:
            dd2[(tuple(h[p] for p in big), chi(h, sb))].add(chi(h, sa))
        rec["chi_b_det_by_big_and_chi_a"] = all(len(v) == 1 for v in dd.values())
        rec["chi_a_det_by_big_and_chi_b"] = all(len(v) == 1 for v in dd2.values())
        rec["order_img_big_plus_a"] = len({(tuple(h[p] for p in big), chi(h, sa)) for h in els})
        rec["order_img_big_plus_b"] = len({(tuple(h[p] for p in big), chi(h, sb)) for h in els})
        rec["order_img_big"] = len({tuple(h[p] for p in big) for h in els})
        # psi = chi_a * chi_b : does it factor through P (action on the 5 pairs)?
        fp = collections.defaultdict(set)
        for h in els:
            fp[pairimg(h)].add(chi(h, sa) ^ chi(h, sb))
        fac = all(len(v) == 1 for v in fp.values())
        rec["chi_a_chi_b_factors_through_P"] = fac
        if fac:
            H = [p for p, v in fp.items() if v == {0}]
            rec["index_of_ker_in_P"] = Po // len(H)
            if len(H) < Po:
                Hg = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in H])
                tr = bool(libgap.IsTransitive(Hg, libgap(list(range(1, m + 1)))))
                rec["ker_in_P"] = ("%dT%d" % (m, int(libgap.TransitiveIdentification(Hg))) if tr else "intransitive") + \
                    " order %d %s" % (len(H), libgap.StructureDescription(Hg))
        # chi_a, chi_b alone factor through P?
        for nm, s in (("a", sa), ("b", sb)):
            fq = collections.defaultdict(set)
            for h in els:
                fq[pairimg(h)].add(chi(h, s))
            rec["chi_%s_factors_through_P" % nm] = all(len(v) == 1 for v in fq.values())
    g7.append(rec)
out["g7"] = g7
for rec in g7:
    print(json.dumps(rec))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ident_out.json"), "w"), indent=0)
print("time", time.time() - T0)
