# verify_trans.py -- INDEPENDENT VERIFIER (pure Python; no Sage/GAP) for trans_g{g}.json (simple CM fields).
#   usage:  python verify_trans.py g [nsample] [maxS]
# From the stored generators alone it recomputes:
#   |G| (own deterministic Schreier-Sims), rho central, transitivity;
#   completeness: the stored Phi's represent ALL G-orbits on the 2^g CM types, each exactly once;
#   per case: orbit size, dim MT, d, Lambda_U (integer kernel, HNF, saturation), primitivity (ALL block systems,
#   enumerated by closure from the base point), and for primitive degenerate cases: imaginary-quadratic subfields and
#   t-values, all CM subfields K' with signatures, the index of degeneracy (free-variable enumeration), and the
#   general Weil tiers (Dlat / Dq for IQ, index<=4, index<=8, all) by an enumeration of the relevant transitive
#   G-sets Z = G/M that is independent of GAP's LowIndexSubgroups:
#     a unit of weight <= 16 is a G-orbit Omega in O x Z of size 16 or 32; size 16 => Omega ~ O and Z is a quotient
#     of O; size 32 => Omega is a double cover Y = G/H of O (H of index 2 in S = Stab(0)) and Z is a quotient of Y.
#     So Z runs over block-system quotients of O and of all double covers Y (index-2 subgroups of S from S/<squares>).
import json, sys, itertools, time, random
from verify_enum import hnf, rank, int_kernel, sat_index, set_orbit


# ---------------- permutations (tuples; p[x] = image of x; rmul(a,b) = first a then b) ----------------
def rmul(a, b):
    return tuple(b[x] for x in a)


def pinv(a):
    r = [0] * len(a)
    for i, x in enumerate(a):
        r[x] = i
    return tuple(r)


def is_id(a):
    return all(i == x for i, x in enumerate(a))


class BSGS:
    """Deterministic Schreier-Sims (Holt, Handbook of CGT, alg. SCHREIERSIMS), right actions."""
    def __init__(self, gens, n, base0=()):
        self.n = n
        gens = [tuple(g) for g in gens if not is_id(g)]
        self.B = list(base0)
        for s in gens:
            if all(s[b] == b for b in self.B):
                self.B.append(next(x for x in range(n) if s[x] != x))
        k = len(self.B)
        self.S = [[s for s in gens if all(s[self.B[j]] == self.B[j] for j in range(i))] for i in range(k)]
        self.U = [None] * k
        for i in range(k):
            self._orbit(i)
        i = k - 1
        while i >= 0:
            jump = None
            for beta in list(self.U[i].keys()):
                ub = self.U[i][beta]
                for x in self.S[i]:
                    bx = x[beta]
                    h = rmul(rmul(ub, x), pinv(self.U[i][bx]))
                    if is_id(h):
                        continue
                    h2, j = self.strip(h, i + 1)
                    if j < len(self.B) or not is_id(h2):
                        if j == len(self.B):
                            self.B.append(next(y for y in range(n) if h2[y] != y))
                            self.S.append([])
                            self.U.append(None)
                        for l in range(i + 1, j + 1):
                            self.S[l].append(h2)
                            self._orbit(l)
                        jump = j
                        break
                if jump is not None:
                    break
            if jump is not None:
                i = jump
            else:
                i -= 1

    def _orbit(self, i):
        b = self.B[i]
        e = tuple(range(self.n))
        T = {b: e}
        st = [b]
        while st:
            x = st.pop()
            for s in self.S[i]:
                y = s[x]
                if y not in T:
                    T[y] = rmul(T[x], s)
                    st.append(y)
        self.U[i] = T

    def strip(self, g, start=0):
        for l in range(start, len(self.B)):
            beta = g[self.B[l]]
            if beta not in self.U[l]:
                return g, l
            g = rmul(g, pinv(self.U[l][beta]))
        return g, len(self.B)

    def order(self):
        o = 1
        for T in self.U:
            o *= len(T)
        return o

    def elements_from(self, level):
        """all elements of the stabilizer of B[:level] (as products of transversals)."""
        els = [tuple(range(self.n))]
        for l in range(len(self.B) - 1, level - 1, -1):
            els = [rmul(e, u) for u in self.U[l].values() for e in els]
        return els


# ---------------- block systems ----------------
def finest_partition(gens, pts, ident):
    """finest G-invariant partition of pts in which all of `ident` lie in one class (union-find + propagation)."""
    par = {x: x for x in pts}
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    def un(x, y):
        x, y = f(x), f(y)
        if x != y:
            par[x] = y
            return True
        return False
    ident = list(ident)
    for y in ident[1:]:
        un(ident[0], y)
    changed = True
    while changed:
        changed = False
        cls = {}
        for x in pts:
            cls.setdefault(f(x), []).append(x)
        for p in gens:
            for cl in cls.values():
                im = [p[x] for x in cl]
                for z in im[1:]:
                    if un(im[0], z):
                        changed = True
    cls = {}
    for x in pts:
        cls.setdefault(f(x), []).append(x)
    return tuple(sorted(tuple(sorted(c)) for c in cls.values()))


def all_block_systems(gens, pts, base):
    """all block systems (as partitions), including the trivial ones, via closure from {base}."""
    start = finest_partition(gens, pts, [base])
    seen = {start}
    fr = [start]
    while fr:
        nf = []
        for part in fr:
            B = next(b for b in part if base in b)
            for y in pts:
                if y in B:
                    continue
                q = finest_partition(gens, pts, list(B) + [y])
                if q not in seen:
                    seen.add(q)
                    nf.append(q)
        fr = nf
    return sorted(seen, key=lambda p: (len(p[0]), p))


# ---------------- G-sets ----------------
def canon_gset(gens_act, npts):
    best = None
    for z0 in range(npts):
        lab = {z0: 0}
        order = [z0]
        i = 0
        while i < len(order):
            z = order[i]
            for p in gens_act:
                y = p[z]
                if y not in lab:
                    lab[y] = len(order)
                    order.append(y)
            i += 1
        if len(order) != npts:
            return None
        key = tuple(tuple(lab[p[order[t]]] for t in range(npts)) for p in gens_act)
        if best is None or key < best:
            best = key
    return best


def wvec(B, g):
    w = [0] * g
    for x in B:
        if x < g:
            w[x] += 1
        else:
            w[x - g] -= 1
    return w


def lat_ok(vs, d):
    if not vs:
        return False, 0
    L = hnf(vs)
    if len(L) < d:
        return False, len(L)
    return sat_index(L) == 1, len(L)


def index_of_degeneracy_v(mu, g, d, nmax=8):
    """free-variable enumeration: v in Z^g with mu v = 0, solved over Q from the free coordinates."""
    from fractions import Fraction
    # RREF over Q
    A = [[Fraction(x) for x in r] for r in mu]
    piv = []
    r = 0
    for c in range(g):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(g) if c not in piv]
    assert len(free) == d
    for nn in range(1, nmax + 1):
        for vals in itertools.product(range(-nn, nn + 1), repeat=d):
            if not any(vals):
                continue
            v = [Fraction(0)] * g
            for c, x in zip(free, vals):
                v[c] = Fraction(x)
            ok = True
            for i, c in enumerate(piv):
                val = -sum(A[i][fc] * Fraction(x) for fc, x in zip(free, vals))
                if val.denominator != 1 or abs(val) > nn:
                    ok = False
                    break
                v[c] = val
            if ok:
                return nn
    return None


# ---------------- general Weil tiers (independent enumeration of Z) ----------------
def gsets_Z(gens, g, n, bs, Dmax, maxS):
    """returns list of (Zgens, Zrho, nZ) for transitive G-sets Z (|Z|<=Dmax, rho free) through which a unit of
    weight <= 16 can factor; None if |S| > maxS."""
    rho = tuple((i + g) % n for i in range(n))
    pts = list(range(n))
    Zs = {}
    def add(Zg, Zr, nz):
        if nz < 2 or nz > Dmax:
            return
        if any(Zr[z] == z for z in range(nz)):
            return
        key = canon_gset(Zg + [Zr], nz)
        if key not in Zs:
            Zs[key] = (Zg, Zr, nz)
    def quotient(gens_act, r_act, part):
        bl = {}
        for i, B in enumerate(part):
            for x in B:
                bl[x] = i
        Zg = [tuple(bl[p[B[0]]] for B in part) for p in gens_act]
        Zr = tuple(bl[r_act[B[0]]] for B in part)
        return Zg, Zr, len(part)
    # quotients of O (including O itself)
    for part in all_block_systems(gens, pts, 0):
        add(*quotient(gens, rho, part))
    # double covers
    sizeS = bs.order() // n
    if sizeS > maxS:
        return None, sizeS
    Sel = bs.elements_from(1)
    assert len(Sel) == sizeS and all(s[0] == 0 for s in Sel)
    Sset = set(Sel)
    # Q = <squares>
    e = tuple(range(n))
    Q = {e}
    Qgens = []
    for s in Sel:
        s2 = rmul(s, s)
        if s2 in Q:
            continue
        Qgens.append(s2)
        fr = list(Q)
        Q = set(Q)
        # closure of <Q, s2>: multiply until stable
        while fr:
            nf = []
            for a in fr:
                for q in Qgens:
                    c = rmul(a, q)
                    if c not in Q:
                        Q.add(c)
                        nf.append(c)
            fr = nf
    # S/Q elementary abelian: coordinates
    lab = {q: 0 for q in Q}
    coset_gens = []
    cur = set(Q)
    for s in Sel:
        if s in lab:
            continue
        bit = 1 << len(coset_gens)
        coset_gens.append(s)
        new = {}
        for x, v in lab.items():
            new[rmul(x, s)] = v | bit
        for y, v in new.items():
            assert y not in lab
        lab.update(new)
    assert len(lab) == sizeS
    r = len(coset_gens)
    # transversal t_x: 0 -> x
    T = {0: e}
    st = [0]
    while st:
        x = st.pop()
        for p in gens:
            y = p[x]
            if y not in T:
                T[y] = rmul(T[x], p)
                st.append(y)
    Tinv = {x: pinv(t) for x, t in T.items()}
    def label_of(gel, x):
        # S-part of gel at x:  t_x * gel * t_{x^gel}^{-1}  (right actions) lies in S
        y = gel[x]
        s = rmul(rmul(T[x], gel), Tinv[y])
        assert s[0] == 0
        return lab[s], y
    for f in range(1, 1 << r):
        chi = lambda v: bin(v & f).count('1') & 1
        def lift(gel):
            out = [None] * (2 * n)
            for x in range(n):
                v, y = label_of(gel, x)
                c = chi(v)
                for eps in (0, 1):
                    out[2 * x + eps] = 2 * y + (eps ^ c)
            return tuple(out)
        Yg = [lift(p) for p in gens]
        Yr = lift(rho)
        ypts = list(range(2 * n))
        for part in all_block_systems(Yg, ypts, 0):
            if len(part) > Dmax:
                continue
            add(*quotient(Yg, Yr, part))
    return list(Zs.values()), sizeS


def weil_tiers_v(gens, g, n, Phi, d, LU, Zlist, Dmax=16):
    Ps = set(Phi)
    allv = []
    O = list(range(n))
    for Zg, Zr, nz in Zlist:
        seen = set()
        units = []
        done = set()
        for x in O:
            for t in range(nz):
                if (x, t) in seen:
                    continue
                om = {(x, t)}
                stck = [(x, t)]
                while stck:
                    a, b = stck.pop()
                    for p, q in zip(gens, Zg):
                        y = (p[a], q[b])
                        if y not in om:
                            om.add(y)
                            stck.append(y)
                seen |= om
                om = frozenset(om)
                if om in done:
                    continue
                done.add(om)
                done.add(frozenset((a, Zr[b]) for a, b in om))
                if len(om) // 2 > Dmax:
                    continue
                Fs = [[a for a, b in om if b == tau] for tau in range(nz)]
                beta = [wvec(F, g) for F in Fs]
                if not any(any(b) for b in beta):
                    continue
                s = [2 * len(Ps & set(F)) - len(F) for F in Fs]
                units.append((len(om) // 2, beta, s))
        if not units:
            continue
        ws = [u[0] for u in units]
        for c in itertools.product(*[range(-(Dmax // w), Dmax // w + 1) for w in ws]):
            D = sum(abs(ci) * w for ci, w in zip(c, ws))
            if D == 0 or D > Dmax:
                continue
            if any(sum(ci * u[2][tau] for ci, u in zip(c, units)) for tau in range(nz)):
                continue
            for tau in range(nz):
                w = [sum(ci * u[1][tau][j] for ci, u in zip(c, units)) for j in range(g)]
                if any(w):
                    if rank(LU + [w]) != d:
                        raise RuntimeError("Weil character outside Lambda_U")
                    allv.append((D, nz, w))
    out = {}
    Ds = sorted(set(z[0] for z in allv))
    for cname, f in dict(IQ=lambda i: i == 2, le4=lambda i: i <= 4, le8=lambda i: i <= 8, all=lambda i: True).items():
        Dlat = Dq = None
        for D in Ds:
            ok, r = lat_ok([z[2] for z in allv if z[0] <= D and f(z[1])], d)
            if Dq is None and r == d:
                Dq = D
            if ok:
                Dlat = D
                break
        out['Dlat_' + cname] = Dlat
        out['Dq_' + cname] = Dq
    return out


# ---------------- main ----------------
def verify_group(G, g, deg_weil=True, sample=None, maxS=400000, Dmax=16):
    n = 2 * g
    gens = [tuple(p) for p in G['gens']]
    rho = tuple((i + g) % n for i in range(n))
    errs = []
    bs = BSGS(gens, n, base0=(0,))
    if bs.order() != G['order']:
        errs.append('order %d vs %d' % (bs.order(), G['order']))
    if bs.strip(rho)[1] != len(bs.B) or not is_id(bs.strip(rho)[0]):
        errs.append('rho not in G')
    for p in gens:
        if any(p[rho[i]] != rho[p[i]] for i in range(n)):
            errs.append('rho not central')
    if len(bs.U[0]) != n:
        errs.append('not transitive')
    systems = [p for p in all_block_systems(gens, list(range(n)), 0) if 1 < len(p[0]) < n]
    if sorted(map(lambda p: tuple(map(tuple, p)), G['systems'])) != sorted(systems):
        errs.append('block systems')
    cm = [p for p in systems if rho[p[0][0]] not in p[0]]
    # completeness of type orbits
    alltypes = set(frozenset(i + g * ((m >> i) & 1) for i in range(g)) for m in range(2 ** g))
    cov = []
    orbs = {}
    for c in G['cases']:
        o = set_orbit(gens, c['Phi'])
        orbs[tuple(c['Phi'])] = o
        cov += list(o)
    if len(cov) != len(set(cov)) or set(cov) != alltypes:
        errs.append('type completeness')
    Zlist = None
    sizeS = None
    res = []
    for c in G['cases']:
        if sample is not None and not (c['primitive'] and c['d'] >= 1) and tuple(c['Phi']) not in sample:
            continue
        ce = []
        Phi = c['Phi']; Ps = set(Phi)
        orb = orbs[tuple(Phi)]
        if len(orb) != c['norb']:
            ce.append('norb')
        if min(tuple(sorted(T)) for T in orb) != tuple(Phi):
            ce.append('Phi not orbit min')
        otypes = [sorted(T) for T in orb]
        mu = [[1 if i in T else -1 for i in range(g)] for T in otypes]
        rmu = rank(mu)
        d = g - rmu
        if d != c['d'] or rmu + 1 != c['dimMT']:
            ce.append('d/dimMT')
        LU = int_kernel(mu, g)
        if hnf(c['LU']) != LU or (LU and sat_index(LU) != 1):
            ce.append('LU')
        prim = not any(all(set(b) <= Ps or not (set(b) & Ps) for b in part) for part in systems)
        if prim != c['primitive']:
            ce.append('primitive')
        weil_checked = False
        if prim and d >= 1:
            iq = sorted((list(p[0]), 2 * len(Ps & set(p[0])) - len(p[0])) for p in cm if len(p) == 2)
            if iq != sorted((x['B'], x['t']) for x in c['IQ']):
                ce.append('IQ')
            sub = []
            for p in cm:
                sig = sorted((list(B), len(Ps & set(B))) for B in p if min(B) < min((x + g) % n for x in B))
                sub.append((len(p), len(p[0]), sig))
            if sorted(sub) != sorted((x['deg'], x['b'], sorted((B, r) for B, r in x['sig'])) for x in c['CMsub']):
                ce.append('CMsub')
            n0 = index_of_degeneracy_v(mu, g, d)
            if n0 != c['n0']:
                ce.append('n0 %s vs %s' % (n0, c['n0']))
            if deg_weil and 'weil' in c:
                if Zlist is None and sizeS is None:
                    Zlist, sizeS = gsets_Z(gens, g, n, bs, Dmax, maxS)
                if Zlist is not None:
                    wt = weil_tiers_v(gens, g, n, Phi, d, LU, Zlist, Dmax)
                    for key, v in wt.items():
                        if c['weil'].get(key) != v:
                            ce.append('weil %s: %s vs %s' % (key, v, c['weil'].get(key)))
                    weil_checked = True
        res.append((tuple(Phi), prim and d >= 1, ce, weil_checked))
    return errs, res, (len(Zlist) if Zlist else None), sizeS


if __name__ == '__main__':
    g = int(sys.argv[1])
    nsample = int(sys.argv[2]) if len(sys.argv) > 2 else -1     # -1: all nondegenerate cases
    maxS = int(sys.argv[3]) if len(sys.argv) > 3 else 400000
    fn = sys.argv[4] if len(sys.argv) > 4 else 'trans_g%d.json' % g
    D = json.load(open(fn))
    out = open('verify_trans_g%d.txt' % g, 'w')
    def P(*a):
        print(*a, flush=True); print(*a, file=out); out.flush()
    t0 = time.time()
    random.seed(20261005)
    nondeg = [(gi, tuple(c['Phi'])) for gi, G in enumerate(D['groups']) for c in G['cases'] if not (c['primitive'] and c['d'] >= 1)]
    if nsample >= 0 and nsample < len(nondeg):
        samp = set(random.sample(nondeg, nsample))
    else:
        samp = set(nondeg)
    stats = dict(groups=0, group_errs=0, deg_cases=0, deg_agree=0, deg_weil_checked=0, nondeg_checked=0, nondeg_agree=0)
    bad = []
    sampd = {}
    for gj, Phi in samp:
        sampd.setdefault(gj, set()).add(Phi)
    for gi, G in enumerate(D['groups']):
        sample = sampd.get(gi, set())
        errs, res, nZ, sizeS = verify_group(G, g, True, sample, maxS)
        stats['groups'] += 1
        if errs:
            stats['group_errs'] += 1
            bad.append((G['TI'], G['j'], errs))
        for Phi, isdeg, ce, wc in res:
            if isdeg:
                stats['deg_cases'] += 1
                stats['deg_agree'] += (not ce)
                stats['deg_weil_checked'] += wc
                if not wc:
                    bad.append((G['TI'], G['j'], list(Phi), 'weil NOT checked |S|=%s' % sizeS))
            else:
                stats['nondeg_checked'] += 1
                stats['nondeg_agree'] += (not ce)
            if ce:
                bad.append((G['TI'], G['j'], list(Phi), ce))
        if gi % 100 == 0:
            P("  %d/%d groups, %.1fs, %s" % (gi, len(D['groups']), time.time() - t0, stats))
    P("FINAL g=%d: %s ; %.1fs" % (g, stats, time.time() - t0))
    for b in bad[:60]:
        P("   ", b)
    P("number of issue lines: %d" % len(bad))
