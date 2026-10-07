# trans_enum.py -- main computation (Sage + libgap).  Simple CM fields of degree n = 2g (transitive G <= W(B_g)).
# usage:  sage -python trans_enum.py g [kmin kmax] [outtag]
# For each TransitiveGroup(n,k) with a central fixed-point-free involution rho' (each choice up to N_{S_n}(T)):
#   conjugate rho' -> standard rho (i -> i+g mod n), G = T^pi <= W(B_g);
#   all 2^g CM types up to G; d, Lambda_U (HNF over Z), primitivity (AllBlocks), N_W(G)-class;
#   for primitive degenerate types: imaginary quadratic subfields (t-values), all CM subfields K' (block systems with
#   rho(B) != B) with signature data, general Weil tiers (weil_tiers.py), index of degeneracy.
# Conventions as in the CM-sixfold enumeration (cmenum.py): X={0..n-1}, rho(i)=i+g mod n, chi -> (chi_i - chi_{i+g})_{i<g}.
import sys, os, json, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sage.all import libgap, matrix, ZZ, QQ, vector, prod
from cmenum import rho_perm, perm_list, all_types, hnf_rows
from weil_tiers import weil_tiers


def wvec(B, g):
    w = [0] * g
    for x in B:
        if x < g:
            w[x] += 1
        else:
            w[x - g] -= 1
    return w


def index_of_degeneracy(LU, g, nmax=8):
    """min over nonzero v in Lambda_U of max_i |v_i|  (= smallest n with a non-Lefschetz Hodge class on A^n)."""
    d = len(LU)
    if d == 0:
        return None, None
    B = matrix(ZZ, LU)
    piv = B.pivots()
    Bp = B.matrix_from_columns(piv)
    Bpi = Bp.inverse()
    for nn in range(1, nmax + 1):
        best = None
        for vals in itertools.product(range(-nn, nn + 1), repeat=d):
            if not any(vals):
                continue
            c = vector(QQ, vals) * Bpi
            if any(x not in ZZ for x in c):
                continue
            v = c * B
            if max(abs(x) for x in v) <= nn:
                v = [int(x) for x in v]
                if best is None or v > best:
                    best = v
        if best is not None:
            return nn, best
    return None, None


def group_info(G, g):
    n = 2 * g
    rho = rho_perm(g)
    gens = list(libgap.SmallGeneratingSet(G))
    gl = [perm_list(x, n) for x in gens]
    systems = []          # list of partitions (tuple of sorted blocks), 1 < block size < n
    for blk in libgap.AllBlocks(G):
        sysorb = libgap.Orbit(G, blk, libgap.OnSets)
        part = tuple(sorted(tuple(sorted(int(x) - 1 for x in b)) for b in sysorb))
        if 1 < len(part[0]) < n:
            systems.append(part)
    systems = sorted(set(systems), key=lambda p: (len(p[0]), p))
    cm = []               # CM subfields K': block systems with rho(B) != B
    for part in systems:
        b0 = part[0]
        if ((b0[0] + g) % n) in b0:
            continue      # rho-stable blocks: totally real subfield
        cm.append(part)
    return gl, systems, cm


def analyse_trans(G, g, want_N=True, weil=True, Dmax=16):
    n = 2 * g
    rho = rho_perm(g)
    gl, systems, cm = group_info(G, g)
    order = int(libgap.Size(G))
    types = all_types(g)
    gtypes = [libgap.Set([x + 1 for x in t]) for t in types]
    torbs = libgap.Orbits(G, gtypes, libgap.OnSets)
    torbs = [sorted(tuple(sorted(int(x) - 1 for x in s)) for s in o) for o in torbs]
    torbs.sort(key=lambda o: o[0])
    nmap = None
    NW_order = None
    if want_N:
        NS = libgap.Normalizer(libgap.SymmetricGroup(n), G)
        NW = libgap.Centralizer(NS, rho)
        NW_order = int(libgap.Size(NW))
        norbs = libgap.Orbits(NW, gtypes, libgap.OnSets)
        nmap = {}
        for k, o in enumerate(norbs):
            for s in o:
                nmap[tuple(sorted(int(x) - 1 for x in s))] = k
    cases = []
    for ot in torbs:
        Phi = ot[0]
        Ps = set(Phi)
        mu = [[1 if i in T else -1 for i in range(g)] for T in ot]
        Mmu = matrix(ZZ, mu)
        rmu = Mmu.rank()
        d = g - rmu
        KU = Mmu.right_kernel_matrix()
        LU = hnf_rows(KU) if KU.nrows() else []
        assert len(LU) == d
        prim = not any(all(set(b) <= Ps or not (set(b) & Ps) for b in part) for part in systems)
        c = dict(Phi=list(Phi), norb=len(ot), dimMT=rmu + 1, d=d, LU=LU, primitive=prim)
        if nmap is not None:
            c['Nclass'] = nmap[tuple(Phi)]
        if prim and d >= 1:
            # imaginary quadratic subfields (2-block CM systems): block containing 0, t = 2|Phi cap B| - |B|
            c['IQ'] = [dict(B=list(p[0]), t=2 * len(Ps & set(p[0])) - len(p[0])) for p in cm if len(p) == 2]
            # all CM subfields K' (block size b, 2k = n/b blocks): for each rho-pair of blocks (B, rho B) with
            # min(B) < min(rho B): r = |Phi cap B| (signature (r, b-r) at the embedding B)
            sub = []
            for p in cm:
                b = len(p[0])
                sig = []
                for B in p:
                    rB = tuple(sorted((x + g) % n for x in B))
                    if min(B) < min(rB):
                        sig.append([list(B), len(Ps & set(B))])
                sub.append(dict(deg=len(p), b=b, sig=sig, balanced=all(2 * r == b for _, r in sig)))
            c['CMsub'] = sub
            nd, vmin = index_of_degeneracy(LU, g)
            c['n0'] = nd
            c['n0_vec'] = vmin
            if weil:
                c['weil'] = weil_tiers(G, gl, g, Phi, d, LU, Dmax=Dmax, maxindex=Dmax)
        cases.append(c)
    return dict(gens=gl, order=order, systems=[[list(b) for b in p] for p in systems],
                cm_systems=[[list(b) for b in p] for p in cm], NW_order=NW_order,
                ntypeorbits=len(torbs)), cases


if __name__ == '__main__':
    g = int(sys.argv[1])
    n = 2 * g
    NT = int(libgap.NrTransitiveGroups(n))
    kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else NT
    tag = sys.argv[4] if len(sys.argv) > 4 else 'all'
    t0 = time.time()
    Sn = libgap.SymmetricGroup(n)
    rho0 = rho_perm(g)
    out = []
    log = open(HERE + '/log_trans_g%d_%s.txt' % (g, tag), 'w')
    def P(*a):
        print(*a, flush=True); print(*a, file=log); log.flush()
    P("g=%d n=%d NrTransitiveGroups=%d range %d..%d" % (g, n, NT, kmin, kmax))
    for k in range(kmin, kmax + 1):
        tk = time.time()
        T = libgap.TransitiveGroup(n, k)
        Z = libgap.Centre(T)
        invs = [z for z in libgap.Elements(Z) if int(libgap.Order(z)) == 2 and int(libgap.NrMovedPoints(z)) == n]
        if not invs:
            continue
        cls = [invs[0]]
        if len(invs) > 1:
            NT_ = libgap.Normalizer(Sn, T)
            cls = []
            for z in invs:
                if not any(bool(libgap.IsConjugate(NT_, z, y)) for y in cls):
                    cls.append(z)
        for j, z in enumerate(cls):
            pi = libgap.RepresentativeAction(Sn, z, rho0)
            G = libgap.ConjugateGroup(T, pi)
            assert bool(libgap.IsSubset(G, [rho0]))
            assert bool(libgap.IsTransitive(G, list(range(1, n + 1))))
            gi, cases = analyse_trans(G, g)
            rec = dict(k=k, TI="%dT%d" % (n, k), j=j, n_central_fpf_inv=len(invs), n_rho_classes=len(cls), **gi,
                       cases=cases, time=time.time() - tk)
            out.append(rec)
            ndeg = sum(1 for c in cases if c['primitive'] and c['d'] >= 1)
            if ndeg or k % 50 == 0:
                P("  %dT%d j=%d |G|=%d typeorbits=%d prim=%d primdeg=%d (%.1fs, total %.1fs)" %
                  (n, k, j, gi['order'], len(cases), sum(c['primitive'] for c in cases), ndeg, time.time() - tk, time.time() - t0))
    json.dump(dict(g=g, kmin=kmin, kmax=kmax, groups=out), open(HERE + '/trans_g%d_%s.json' % (g, tag), 'w'), default=int)
    P("done g=%d: (G,rho) pairs %d, time %.1fs" % (g, len(out), time.time() - t0))
