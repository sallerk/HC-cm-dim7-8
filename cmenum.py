# cmenum.py -- main computation library (Sage + GAP via libgap).  Run inside the sagemath docker image with `sage -python`.
#
# Conventions: X = {0..2g-1}, rho(i) = i+g mod 2g.  W = W(B_g) = centralizer of rho in S_2g.
# A CM type Phi = one point of each rho-pair.  Coordinates of X*(U_E) = Z^X/even  ~=  Z^g:
#   chi |-> (chi_i - chi_{i+g})_{i<g};  [i] -> e_i, [i+g] -> -e_i.
# Hodge cocharacters mu_h = 1_{hPhi} - 1_{rho hPhi}  -> in Z^g: +1 at i if i in hPhi, -1 if i+g in hPhi.
# Lambda_U = {chi in Z^g : chi . mu_h = 0 for all h}.
import itertools
from sage.all import libgap, matrix, ZZ, QQ, prod


def rho_perm(g):
    return libgap.PermList([((i + g) % (2 * g)) + 1 for i in range(2 * g)])


def W_group(g):
    rho = rho_perm(g)
    return libgap.Centralizer(libgap.SymmetricGroup(2 * g), rho), rho


def perm_list(p, n):
    return [int(x) - 1 for x in libgap.ListPerm(p, n)]


def all_types(g):
    return [tuple(sorted(i + g * ((m >> i) & 1) for i in range(g))) for m in range(2 ** g)]


def group_data(G, g, W=None, want_normalizer=True):
    """Group-level data for G <= W(B_g) containing rho (GAP group on 1..2g)."""
    n = 2 * g
    rho = rho_perm(g)
    assert bool(libgap.IsSubset(G, [rho])), "rho not in G"
    gens = list(libgap.SmallGeneratingSet(G))
    if len(gens) == 0:
        gens = [rho]
    gl = [perm_list(x, n) for x in gens]
    orbs = sorted([sorted(int(x) - 1 for x in o) for o in libgap.Orbits(G, list(range(1, n + 1)))], key=min)
    oinfo = []
    for oi, O in enumerate(orbs):
        O1 = [x + 1 for x in O]
        act = libgap.Action(G, O1)
        ti = int(libgap.TransitiveIdentification(act))
        img = int(libgap.Size(act))
        systems = set()
        for blk in libgap.AllBlocks(act):
            sysorb = libgap.Orbit(act, blk, libgap.OnSets)
            part = tuple(sorted(tuple(sorted(O[int(k) - 1] for k in b)) for b in sysorb))
            bs = len(part[0])
            if 1 < bs < len(O):
                systems.add(part)
        systems = sorted(systems)
        quad = []
        cand = list(systems)
        if len(O) == 2:   # the 2-block system of a size-2 orbit is the (trivial) singleton partition
            cand.append(((O[0],), (O[1],)))
        for part in cand:
            if len(part) == 2:
                b0, b1 = set(part[0]), set(part[1])
                if all(((x + g) % n) in b1 for x in b0):
                    quad.append(part)
        oinfo.append(dict(orbit=O, size=len(O), TI="%dT%d" % (len(O), ti), img_order=img,
                          systems=[list(map(list, p)) for p in systems], quad=quad))
    # imaginary quadratic characters: sign vectors on generators
    Kdict = {}
    for oi, inf in enumerate(oinfo):
        O = inf['orbit']
        for part in inf['quad']:
            Bc = part[0] if min(O) in part[0] else part[1]
            Bs = set(Bc)
            other = set(part[0]) | set(part[1])
            other -= Bs
            sv = []
            for p in gl:
                im = set(p[x] for x in Bc)
                if im == Bs:
                    sv.append(1)
                elif im == other:
                    sv.append(-1)
                else:
                    raise RuntimeError("block system not invariant")
            sv = tuple(sv)
            Kdict.setdefault(sv, {})
            assert oi not in Kdict[sv]
            Kdict[sv][oi] = sorted(Bc)
    Ks = [dict(signs=list(sv), blocks=Kdict[sv]) for sv in sorted(Kdict)]
    # G-set isomorphisms between distinct orbits (via point stabilizers)
    isos = {}
    for i in range(len(orbs)):
        for j in range(i + 1, len(orbs)):
            if len(orbs[i]) != len(orbs[j]):
                continue
            a = orbs[i][0] + 1
            Sa = libgap.Stabilizer(G, a)
            trans = {x: libgap.RepresentativeAction(G, a, x + 1) for x in orbs[i]}
            fl = []
            for b in orbs[j]:
                Sb = libgap.Stabilizer(G, b + 1)
                if bool(Sa == Sb):
                    f = {x: int(libgap.OnPoints(b + 1, trans[x])) - 1 for x in orbs[i]}
                    assert sorted(f.values()) == orbs[j]
                    fl.append(f)
            if fl:
                isos[(i, j)] = fl
    # orbits of G on the 2^g CM types
    types = all_types(g)
    gtypes = [libgap.Set([x + 1 for x in t]) for t in types]
    torbs = libgap.Orbits(G, gtypes, libgap.OnSets)
    torbs = [sorted(tuple(sorted(int(x) - 1 for x in s)) for s in o) for o in torbs]
    torbs.sort(key=lambda o: o[0])
    nmap = None
    if want_normalizer and W is not None:
        N = libgap.Normalizer(W, G)
        norbs = libgap.Orbits(N, gtypes, libgap.OnSets)
        nmap = {}
        for k, o in enumerate(norbs):
            for s in o:
                nmap[tuple(sorted(int(x) - 1 for x in s))] = k
    return dict(g=g, gens=gl, order=int(libgap.Size(G)), orbits=orbs, oinfo=oinfo, Ks=Ks,
                isos=isos, torbs=torbs, nmap=nmap)


def hnf_rows(M):
    """Nonzero rows of the HNF (Sage echelon form over ZZ)."""
    if M.nrows() == 0:
        return []
    E = M.echelon_form()
    return [list(map(int, r)) for r in E.rows() if not r.is_zero()]


def type_analysis(gd, otypes):
    g = gd['g']
    n = 2 * g
    Phi = otypes[0]
    Ps = set(Phi)
    mu = [[1 if i in T else -1 for i in range(g)] for T in otypes]
    S = [[1 if j in T else 0 for j in range(n)] for T in otypes]
    Mmu = matrix(ZZ, mu)
    rmu = Mmu.rank()
    rS = matrix(ZZ, S).rank()
    assert rS == rmu + 1
    d = g - rmu
    KU = Mmu.right_kernel_matrix()
    LU = hnf_rows(KU) if KU.nrows() else []
    assert len(LU) == d
    # admissible Weil characters
    adm = []
    for kidx, K in enumerate(gd['Ks']):
        oids = sorted(K['blocks'].keys())
        t = {o: 2 * len(Ps & set(K['blocks'][o])) - len(K['blocks'][o]) for o in oids}
        for r in range(1, len(oids) + 1):
            for J in itertools.combinations(oids, r):
                for signs in itertools.product([1, -1], repeat=r - 1):
                    eps = (1,) + signs
                    if sum(e * t[o] for e, o in zip(eps, J)) != 0:
                        continue
                    B = []
                    for e, o in zip(eps, J):
                        blk = K['blocks'][o]
                        B += blk if e == 1 else [(x + g) % n for x in blk]
                    XJ = sorted(x for o in J for x in gd['orbits'][o])
                    if XJ[0] not in B:
                        B = [(x + g) % n for x in B]
                    B = sorted(B)
                    w = [0] * g
                    for x in B:
                        if x < g:
                            w[x] += 1
                        else:
                            w[x - g] -= 1
                    assert all(sum(a * b for a, b in zip(w, m)) == 0 for m in mu)
                    assert 2 * len(Ps & set(B)) == len(B)
                    adm.append(dict(K=kidx, J=list(J), B=B, size=len(B), w=w))
    if adm:
        MW = matrix(ZZ, [a['w'] for a in adm])
        rW = MW.rank()
        LW = hnf_rows(MW)
        ed = [int(x) for x in MW.elementary_divisors() if x != 0]
        satidx = int(prod(ed)) if ed else 1
    else:
        rW, LW, satidx = 0, [], 1
    index = None
    if rW == d:
        if d == 0:
            index = 1
        else:
            C = matrix(QQ, LU).solve_left(matrix(QQ, LW))
            assert all(x in ZZ for x in C.list())
            index = abs(int(C.det()))
            assert index == satidx
    passed = (rW == d and index == 1)
    # factor types
    prim, coarsest = [], []
    for oi, inf in enumerate(gd['oinfo']):
        Pi = Ps & set(inf['orbit'])
        quals = [p for p in inf['systems'] if all(set(b) <= Pi or not (set(b) & Pi) for b in p)]
        prim.append(len(quals) == 0)
        if quals:
            c = max(quals, key=lambda p: len(p[0]))
            # theory check: every qualifying system refines the coarsest one
            for p in quals:
                for b in p:
                    assert any(set(b) <= set(cb) for cb in c)
            coarsest.append(c)
        else:
            coarsest.append(None)
    isog = []
    for (i, j), fl in gd['isos'].items():
        Pi = Ps & set(gd['orbits'][i])
        Pj = Ps & set(gd['orbits'][j])
        if any(set(f[x] for x in Pi) == Pj for f in fl):
            isog.append([i, j])
    isog.sort()
    reduced = all(prim) and not isog
    return dict(Phi=list(Phi), norb=len(otypes), dimMT=rS, d=d, LU=LU, adm=adm, rankW=rW, LW=LW,
                index=index, satidx=satidx, passed=passed, primitive=prim, isog=isog, reduced=reduced,
                coarsest=coarsest)


def analyse_group(G, g, W=None, want_normalizer=True):
    gd = group_data(G, g, W, want_normalizer)
    res = []
    for ot in gd['torbs']:
        ta = type_analysis(gd, ot)
        if gd['nmap'] is not None:
            ta['Nclass'] = gd['nmap'][tuple(ot[0])]
        res.append(ta)
    return gd, res


# ---------------- reduction of a non-reduced (G,Phi) to its reduced model ----------------
def _gset_iso_exists(gl, Oa, Ob, Pa, Pb, act):
    """act(p, x) = image of x under generator p (on quotient points). BFS-based G-set iso Oa->Ob with f(Pa)=Pb."""
    a0 = Oa[0]
    for b0 in Ob:
        f = {a0: b0}
        stack = [a0]
        ok = True
        while stack and ok:
            x = stack.pop()
            for p in gl:
                y, fy = act(p, x), act(p, f[x])
                if y in f:
                    if f[y] != fy:
                        ok = False
                        break
                else:
                    f[y] = fy
                    stack.append(y)
        if ok and len(f) == len(Oa) and len(set(f.values())) == len(Oa):
            if set(f[x] for x in Pa) == set(Pb):
                return True
    return False


def reduced_model(gd, ta):
    """Build reduced model (G', Phi') : quotient each orbit by its coarsest Phi-compatible block system,
    keep one orbit per isogeny class.  Returns (gens', g', Phi') with standard labelling rho(i)=i+g'."""
    g = gd['g']
    n = 2 * g
    gl = gd['gens']
    Ps = set(ta['Phi'])
    # quotient points: (orbit index, block tuple)
    qorbs = []
    for oi, inf in enumerate(gd['oinfo']):
        c = ta['coarsest'][oi]
        if c is None:
            c = [[x] for x in inf['orbit']]
        qorbs.append([tuple(b) for b in c])
    pt2blk = {}
    for oi, blks in enumerate(qorbs):
        for b in blks:
            for x in b:
                pt2blk[x] = b

    def act(p, b):
        return pt2blk[p[b[0]]]

    qPhi = [[b for b in blks if set(b) <= Ps] for blks in qorbs]
    # isogeny classes among quotient factors
    reps = []
    for oi in range(len(qorbs)):
        iso = False
        for rj in reps:
            if len(qorbs[rj]) == len(qorbs[oi]) and _gset_iso_exists(gl, qorbs[rj], qorbs[oi], qPhi[rj], qPhi[oi], act):
                iso = True
                break
        if not iso:
            reps.append(oi)
    pts = [b for oi in reps for b in qorbs[oi]]
    gp = len(pts) // 2
    rho_b = lambda b: pt2blk[(b[0] + g) % n]
    # standard labels: choose first of each rho-pair as i, partner as i+gp
    lab = {}
    k = 0
    for b in pts:
        if b in lab:
            continue
        lab[b] = k
        lab[rho_b(b)] = k + gp
        k += 1
    assert k == gp
    gens2 = []
    for p in gl:
        q = [None] * (2 * gp)
        for b in pts:
            q[lab[b]] = lab[act(p, b)]
        gens2.append(q)
    Phi2 = sorted(lab[b] for oi in reps for b in qPhi[oi])
    return gens2, gp, Phi2


def analyse_reduced_model(gd, ta):
    gens2, gp, Phi2 = reduced_model(gd, ta)
    G2 = libgap.Group([libgap.PermList([x + 1 for x in q]) for q in gens2])
    gd2 = group_data(G2, gp, None, False)
    for ot in gd2['torbs']:
        if tuple(Phi2) in ot:
            ta2 = type_analysis(gd2, ot)
            return gd2, ta2, Phi2
    raise RuntimeError("Phi2 not found")


def reduced_dimMT(gd, ta, otypes):
    """Cheap check: rank of S for the reduced model computed from the image of the G-orbit of Phi."""
    gens2, gp, Phi2 = reduced_model(gd, ta)
    # rebuild the point->label map exactly as reduced_model does
    g = gd['g']; n = 2 * g
    Ps = set(ta['Phi'])
    qorbs = []
    for oi, inf in enumerate(gd['oinfo']):
        c = ta['coarsest'][oi]
        if c is None:
            c = [[x] for x in inf['orbit']]
        qorbs.append([tuple(b) for b in c])
    # S' rows: for each h Phi in the orbit, indicator of blocks inside hPhi, over all quotient blocks (all orbits);
    # dropping isogenous copies does not change the rank (isogenous rows are permuted copies), so use all blocks.
    blocks = [b for blks in qorbs for b in blks]
    rows = [[1 if set(b) <= set(T) else 0 for b in blocks] for T in otypes]
    return matrix(ZZ, rows).rank(), gp


# ---------------- extension: imaginary-quadratic Weil characters WITH MULTIPLICITIES ----------------
# A Weil class for K on a sub-product of a POWER of A, prod_i A_i^{m_i} (i in J_K), with K acting on m_i^+ copies
# via iota and m_i^- copies via conj(iota):  character  w_c = sum_i c_i b_i  (c_i = m_i^+ - m_i^-),
# balanced iff sum_i c_i t_i = 0, dimension D(c) = sum_i |c_i| g_i  (g_i = |O_i|/2).
# The strict definition (each simple factor used at most once) is the special case c_i in {-1,0,1}.
def _cvecs(sizes, Dmax):
    k = len(sizes)
    out = []
    def rec(i, cur, used):
        if i == k:
            if any(cur):
                out.append(tuple(cur))
            return
        m = (Dmax - used) // sizes[i]
        for c in range(-m, m + 1):
            cur.append(c)
            rec(i + 1, cur, used + abs(c) * sizes[i])
            cur.pop()
    rec(0, [], 0)
    return out


def mult_weil(gd, ta, Dmax=12):
    g = gd['g']
    n = 2 * g
    Ps = set(ta['Phi'])
    d = ta['d']
    LU = matrix(ZZ, ta['LU']) if d else None
    vecs = []   # (D, w, K, c, J)
    for kidx, K in enumerate(gd['Ks']):
        oids = sorted(K['blocks'].keys())
        b, t, gs = [], [], []
        for o in oids:
            blk = K['blocks'][o]
            w = [0] * g
            for x in blk:
                if x < g:
                    w[x] += 1
                else:
                    w[x - g] -= 1
            b.append(w)
            t.append(2 * len(Ps & set(blk)) - len(blk))
            gs.append(len(blk))
        for c in _cvecs(gs, Dmax):
            if sum(ci * ti for ci, ti in zip(c, t)) != 0:
                continue
            w = [sum(ci * bi[j] for ci, bi in zip(c, b)) for j in range(g)]
            if not any(w):
                continue
            D = sum(abs(ci) * gi for ci, gi in zip(c, gs))
            vecs.append((D, w, kidx, c, oids))
    vecs.sort(key=lambda v: v[0])
    Dmin = None
    qrank = matrix(ZZ, [v[1] for v in vecs]).rank() if vecs else 0
    if d == 0:
        Dmin = 0
    else:
        for D in sorted(set(v[0] for v in vecs)):
            M = matrix(ZZ, [v[1] for v in vecs if v[0] <= D])
            if M.rank() == d:
                ed = [int(x) for x in M.elementary_divisors() if x != 0]
                if prod(ed) == 1:
                    Dmin = D
                    break
    # smallest-D generators that realize Dmin (for reporting)
    ex = [dict(D=v[0], w=v[1], K=v[2], c=list(v[3]), J=v[4]) for v in vecs if Dmin is not None and v[0] <= Dmin][:6]
    return dict(Dmin=Dmin, qrank=qrank, examples=ex, Dmax=Dmax)
