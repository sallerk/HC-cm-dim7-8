# Part 2: independent control "G = P x C2" for g = 8 degenerate CM types.
# Pure Python 3 + numpy. Input: transgroups_5_6_8.json (dumped from GAP by dump.py).
# Stored data (read-only) only used for the final comparison.
import json, sys, time, itertools, collections
from fractions import Fraction
import numpy as np

T0 = time.time()
import os
CLONE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")) + "/"
TG = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "transgroups_5_6_8.json")))["8"]
SH16 = (np.arange(16, dtype=np.uint64) * np.uint64(4))
SH8 = (np.arange(8, dtype=np.uint64) * np.uint64(4))


def closure(gens, n):
    idt = tuple(range(n))
    seen = {idt}
    frontier = [idt]
    gens = [tuple(g) for g in gens]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = tuple(g[a[i]] for i in range(n))  # g o a
                if b not in seen:
                    seen.add(b); nxt.append(b)
        frontier = nxt
    return np.array(sorted(seen), dtype=np.int64)


def codes(arr, sh):
    return (arr.astype(np.uint64) << sh).sum(axis=1)


def member(sorted_codes, c):
    idx = np.searchsorted(sorted_codes, c)
    idx2 = np.minimum(idx, len(sorted_codes) - 1)
    return (idx < len(sorted_codes)) & (sorted_codes[idx2] == c)


def rank_null(M):
    """exact rank and rational nullspace basis of an integer matrix (list of rows)"""
    A = [[Fraction(x) for x in r] for r in M]
    m, n = len(A), len(A[0])
    piv = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == m:
            break
    free = [c for c in range(n) if c not in piv]
    null = []
    for fc in free:
        v = [Fraction(0)] * n
        v[fc] = Fraction(1)
        for i, pc in enumerate(piv):
            v[pc] = -A[i][fc]
        null.append(v)
    return r, null


def prim_int(v):
    from math import gcd
    den = 1
    for x in v:
        den = den * x.denominator // gcd(den, x.denominator)
    w = [int(x * den) for x in v]
    g = 0
    for x in w:
        g = gcd(g, abs(x))
    w = [x // g for x in w]
    if next(x for x in w if x != 0) < 0:
        w = [-x for x in w]
    return w


# ---------- block systems on 16 points (bitmasks) ----------
def block_systems(gens16):
    tabs = []
    for g in gens16:
        lo = np.zeros(256, dtype=np.int64); hi = np.zeros(256, dtype=np.int64)
        for m in range(256):
            a = b = 0
            for i in range(8):
                if m >> i & 1:
                    a |= 1 << g[i]; b |= 1 << g[i + 8]
            lo[m] = a; hi[m] = b
        tabs.append((lo, hi))
    systems = set()
    for size in (2, 4, 8):
        for rest in itertools.combinations(range(1, 16), size - 1):
            B = 1
            for x in rest:
                B |= 1 << x
            orb = {B}; frontier = [B]; ok = True
            while frontier and ok:
                nf = []
                for X in frontier:
                    for lo, hi in tabs:
                        Y = int(lo[X & 255] | hi[X >> 8])
                        if Y in orb:
                            continue
                        if any(Y & Z for Z in orb):
                            ok = False; break
                        orb.add(Y); nf.append(Y)
                    if not ok:
                        break
                frontier = nf
            if ok and len(orb) * size == 16:
                systems.add(frozenset(orb))
    return [sorted(s) for s in systems]


def primitive_mask(phimask, systems):
    for S in systems:
        if all((phimask & b) == 0 or (phimask & b) == b for b in S):
            return False
    return True


# ---------- W(B8) ----------
SIG = np.array(list(itertools.permutations(range(8))), dtype=np.int64)   # 40320 x 8
EPS = ((np.arange(256)[:, None] >> np.arange(8)[None, :]) & 1).astype(np.int64)  # 256 x 8


def w_perm16(sig, eps):
    """rows: j+8e -> sig(j) + 8(e xor eps_j)"""
    W = np.empty((sig.shape[0], 16), dtype=np.int64)
    W[:, :8] = sig + 8 * eps
    W[:, 8:] = sig + 8 * (1 - eps)
    return W


def normalizer(gens16, Gcodes, batch=32):
    surv_sig, surv_eps = [], []
    nsig = SIG.shape[0]
    for e0 in range(0, 256, batch):
        eb = EPS[e0:e0 + batch]
        sig = np.tile(SIG, (eb.shape[0], 1))
        eps = np.repeat(eb, nsig, axis=0)
        W = w_perm16(sig, eps)
        keep = np.ones(W.shape[0], dtype=bool)
        for g in gens16:
            idx = np.nonzero(keep)[0]
            if idx.size == 0:
                break
            Wk = W[idx]
            R = np.empty_like(Wk)
            np.put_along_axis(R, Wk, Wk[:, g], axis=1)   # R[w(i)] = w(g(i))
            keep[idx] = member(Gcodes, codes(R, SH16))
        idx = np.nonzero(keep)[0]
        surv_sig.append(sig[idx]); surv_eps.append(eps[idx])
    return np.concatenate(surv_sig), np.concatenate(surv_eps)


E256 = EPS.copy()   # type t <-> e-vector bits (e_j = bit j of t); Phi = {j + 8 e_j}
PHIMASK = np.array([sum(1 << (j + 8 * int(E256[t, j])) for j in range(8)) for t in range(256)])


def min_image(sig, eps, chunk=2048):
    """for each type t: min over elements (sig,eps) of image code; image e'_{sig(j)} = e_j xor eps_j"""
    best = np.full(256, 1 << 30, dtype=np.int64)
    for a in range(0, sig.shape[0], chunk):
        s = sig[a:a + chunk]; e = eps[a:a + chunk]
        X = E256[None, :, :] ^ e[:, None, :]
        wts = (1 << s)[:, None, :]
        c = (X * wts).sum(-1)
        best = np.minimum(best, c.min(axis=0))
    return best


def to_e16(phi_list):
    """CM type on 16 points (pairs {j,j+8}) -> type index t"""
    s = set(phi_list)
    return sum(1 << j for j in range(8) if (j + 8) in s)


# ================= main loop over P = 8Tk =================
RES = {}
CREP = {}
CLASSES = {}   # (k, class_rep_t) -> info
Pcodes_all = {}
neg_control = None
lam_ok = 0; lam_bad = []
for P in TG:
    k = P["k"]; t1 = time.time()
    Pel = closure(P["gens"], 8)
    assert len(Pel) == P["order"]
    Pcodes_all[k] = np.sort(codes(Pel, SH8))
    Pinv = np.argsort(Pel, axis=1)
    # G = P x C2
    Gsig = np.concatenate([Pel, Pel]); Geps = np.concatenate([np.zeros_like(Pel), np.ones_like(Pel)])
    Gel = w_perm16(Gsig, Geps)
    Gcodes = np.sort(codes(Gel, SH16))
    assert len(np.unique(Gcodes)) == 2 * len(Pel)
    gens16 = [list(w_perm16(np.array([g]), np.zeros((1, 8), dtype=np.int64))[0]) for g in P["gens"]]
    rho = [(i + 8) % 16 for i in range(16)]
    gensG = gens16 + [rho]
    systems = block_systems(gensG)
    Nsig, Neps = normalizer(gensG, Gcodes)
    NWo = len(Nsig)
    assert NWo % len(Gcodes) == 0
    # sanity: G inside N
    Ncodes = np.sort(codes(w_perm16(Nsig, Neps), SH16))
    assert member(Ncodes, Gcodes).all()
    orep = min_image(Gsig, Geps)
    crep = min_image(Nsig, Neps)
    CREP[k] = crep
    # per orbit: d, primitivity, signature
    info = {}
    for t in sorted(set(orep.tolist())):
        s = 1 - 2 * E256[t]                 # mu_Phi: +1 iff e_j = 0
        V = s[Pinv]                          # V[p, i] = s[p^-1(i)] = mu_{p Phi}
        Gram = (V.T @ V).tolist()
        rk, null = rank_null(Gram)
        d = 8 - rk
        prim = primitive_mask(int(PHIMASK[t]), systems)
        a = int((E256[t] == 0).sum()); sig = (max(a, 8 - a), min(a, 8 - a))
        info[t] = dict(d=d, prim=prim, sig=sig)
        if d == 1 and prim:
            v = prim_int(null[0])
            mus = np.unique(np.concatenate([V, -V]), axis=0)
            wK = np.ones(8, dtype=np.int64)
            if v == [1] * 8 and (mus @ wK == 0).all() and rk == 7:
                lam_ok += 1
            else:
                lam_bad.append((k, t, v))
        if neg_control is None and sig == (5, 3):
            mus = np.concatenate([V, -V])
            neg_control = dict(k=k, t=t, Phi=[j + 8 * int(E256[t, j]) for j in range(8)],
                               wK_dot_mu_values=sorted(set((mus @ np.ones(8, dtype=np.int64)).tolist())))
    for t in range(256):
        assert info[int(orep[t])]["sig"] == info[int(orep[t])]["sig"]
    # classes
    cls = collections.defaultdict(set)
    for t in range(256):
        cls[int(crep[t])].add(int(orep[t]))
    # consistency: d/prim/sig constant on classes
    tab = collections.defaultdict(lambda: [0, 0])
    sigs = collections.defaultdict(set)
    for c, orbs in cls.items():
        vals = {(info[o]["d"], info[o]["prim"]) for o in orbs}
        assert len(vals) == 1, (k, c, vals)
        d, prim = vals.pop()
        sig = tuple(sorted({info[o]["sig"] for o in orbs}))   # signatures w.r.t. the K of the construction
        CLASSES[(k, c)] = dict(k=k, d=d, prim=prim, sig=sig, norb=len(orbs), G=len(Gcodes), NW=NWo,
                               Phi=[j + 8 * int(E256[c, j]) for j in range(8)])
        if prim:
            tab[d][0] += 1; tab[d][1] += len(orbs); sigs[d].update(sig)
    RES[k] = dict(G=len(Gcodes), NW=NWo, nsys=len(systems),
                  prim={d: dict(classes=v[0], orbits=v[1], sigs=sorted(sigs[d])) for d, v in sorted(tab.items())},
                  time=round(time.time() - t1, 1))
    print(k, RES[k], flush=True)

print("loop time", round(time.time() - T0, 1), flush=True)

# ================= totals =================
tot = collections.Counter(); totorb = collections.Counter(); sigtot = collections.defaultdict(collections.Counter)
for (k, c), x in CLASSES.items():
    if x["prim"]:
        tot[x["d"]] += 1; totorb[x["d"]] += x["norb"]; sigtot[x["d"]][x["sig"]] += 1
print("TOTAL primitive classes by d:", dict(sorted(tot.items())))
print("TOTAL primitive orbits  by d:", dict(sorted(totorb.items())))
print("signatures by d:", {d: dict(v) for d, v in sorted(sigtot.items())})
print("Lambda_U = Z wK checks ok (orbits, d=1 prim):", lam_ok, "bad:", lam_bad)
print("negative control:", neg_control)

# ================= comparison with stored data =================
D = json.load(open(CLONE + "deg_classes_g8.json"))
TR = json.load(open(CLONE + "trans_g8.json"))["groups"]
NWst = {(g["TI"], g["j"]): g["NW_order"] for g in TR}
mine = collections.defaultdict(collections.Counter)
for x in CLASSES.values():
    if x["prim"] and x["d"] >= 1:
        mine[x["d"]][(x["G"], x["norb"], x["NW"])] += 1
st = collections.defaultdict(collections.Counter)
for x in D:
    st[x["d"]][(x["order"], len(x["Phi_all_G_orbit_reps"]), NWst[(x["TI"], x["j"])])] += 1
for d in sorted(set(mine) | set(st)):
    a, b = mine[d], st[d]
    print("d=%d: mine %d classes, stored %d; mine<=stored: %s; equal: %s; mine-stored: %s" %
          (d, sum(a.values()), sum(b.values()), not (a - b), a == b, dict(a - b)))
    print("    stored-mine:", dict(sorted((b - a).items())))

# exact matching of stored classes into the P x C2 table
def ident8(gens8, order):
    cands = [P["k"] for P in TG if P["order"] == order]
    hits = []
    for k in cands:
        keep = np.ones(SIG.shape[0], dtype=bool)
        for g in gens8:
            R = np.empty_like(SIG)
            np.put_along_axis(R, SIG, SIG[:, g], axis=1)
            keep &= member(Pcodes_all[k], codes(R, SH8))
        idx = np.nonzero(keep)[0]
        if idx.size:
            hits.append((k, SIG[idx[0]]))
    return hits

match = collections.Counter(); unmatched = []; notPxC2 = collections.Counter(); hitcount = collections.Counter()
per_k_stored = collections.defaultdict(collections.Counter)
mism = []
for n, x in enumerate(D):
    gens = x["gens"]
    systems = block_systems(gens + [x["rho"]])
    prim = primitive_mask(sum(1 << i for i in x["Phi"]), systems)
    iq = [S for S in systems if len(S) == 2 and all(((b >> i) & 1) != ((b >> ((i + 8) % 16)) & 1)
                                                    for b in S for i in range(16))]
    if not prim:
        mism.append(("stored class not primitive by my check", n))
    if not iq:
        notPxC2[x["d"]] += 1
        continue
    keys = set()
    for S in iq:
        B0 = S[0]
        f = [(i % 8) + 8 * (0 if (B0 >> i) & 1 else 1) for i in range(16)]
        finv = [0] * 16
        for i in range(16):
            finv[f[i]] = i
        gp = [[f[g[finv[y]]] for y in range(16)] for g in gens]
        P8 = [[gg[j] % 8 for j in range(8)] for gg in gp]
        assert x["order"] % 2 == 0
        hits = ident8(P8, x["order"] // 2)
        if len(hits) != 1:
            mism.append(("ident8 not unique/none", n, [h[0] for h in hits])); continue
        k, sg = hits[0]
        w = [int(sg[f[i] % 8]) + 8 * (f[i] // 8) for i in range(16)]
        winv = [0] * 16
        for i in range(16):
            winv[w[i]] = i
        # check w G w^-1 inside P_k x C2 (equal orders => equal)
        Gk = np.sort(codes(w_perm16(np.concatenate([closure(TG[k - 1]["gens"], 8)] * 2),
                                    np.concatenate([np.zeros((TG[k - 1]["order"], 8), dtype=np.int64),
                                                    np.ones((TG[k - 1]["order"], 8), dtype=np.int64)])), SH16))
        cg = np.array([[w[g[winv[y]]] for y in range(16)] for g in gens + [x["rho"]]])
        assert member(Gk, codes(cg, SH16)).all() and len(Gk) == x["order"]
        phi2 = [w[i] for i in x["Phi"]]
        t = to_e16(phi2)
        keys.add((k, t))
    # resolve class rep for each (k, t) found
    reps = {(k, int(CREP[k][t])) for (k, t) in keys}
    hitcount[len(reps)] += 1
    if len(reps) != 1:
        mism.append(("different IQ systems give different classes", n, reps)); continue
    (k, c), = reps
    y = CLASSES[(k, c)]
    if (y["d"], y["norb"], y["G"], y["NW"]) != (x["d"], len(x["Phi_all_G_orbit_reps"]), x["order"], NWst[(x["TI"], x["j"])]):
        mism.append(("invariant mismatch", n, x["TI"], x["d"], y))
    match[(k, c)] += 1
    per_k_stored[x["d"]][k] += 1
print("stored classes with no rho-swapped 2-block system (not P x C2), by d:", dict(notPxC2))
print("number of distinct classes obtained from the different IQ systems of one stored class:", dict(hitcount))
dup = [kc for kc, v in match.items() if v > 1]
print("my classes hit by >1 stored class:", dup)
for d in (1, 2, 3, 4):
    mineset = {kc for kc, y in CLASSES.items() if y["prim"] and y["d"] == d}
    hit = {kc for kc in match if CLASSES[kc]["d"] == d}
    print("d=%d: my primitive classes %d, matched by stored %d, my classes not in stored: %s" %
          (d, len(mineset), len(hit), sorted((kc[0], CLASSES[kc]["Phi"]) for kc in mineset - hit)))
    print("    stored per k:", dict(sorted(per_k_stored[d].items())))
print("mismatches:", mism)
json.dump(dict(RES={str(k): v for k, v in RES.items()},
               classes=[dict(v, key=[k, c]) for (k, c), v in CLASSES.items() if v["prim"] and v["d"] >= 1]),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pxc_results.json"), "w"), default=str, indent=0)
print("total time", round(time.time() - T0, 1))
