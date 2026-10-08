# verify_blocks.py -- separately written check of the "blocks" step: for which reduced CM pairs (G, Phi) do the Weil
# characters of sub-products of powers of A of dimension at most 6, of the kinds whose Weil classes are known to be
# algebraic (assuming Markman's theorem), generate Lambda_U?
#
# Pure Python, standard library only; it imports nothing from cmenum.py, verify_enum.py, used_blocks_check*.py or the
# Sage/GAP programs.  From enum_g<g>.json it reads the generators of G, Phi, and the stored fields it compares with
# (order, d, LU, reduced, passed, mult.Dmin); the case list itself (one Phi per G-orbit of CM types) is taken from the
# file.  Everything else is recomputed:
#   - for EVERY case: the G-orbit of Phi and Lambda_U = {w in Z^g : <w, mu_Psi> = 0 for all Psi in G.Phi}, as an
#     integer kernel (unimodular column reduction), so d = 0 is re-checked too; and whether the pair is reduced:
#     the CM type on each G-orbit O of X is primitive (Phi & O is not a union of blocks of a non-trivial block system
#     of G on O; it suffices to test the minimal block systems in which x0 and y lie in one block, y in O), and no
#     two orbits are isomorphic as G-sets with CM types (the G-set isomorphisms O -> O' are found by extending x0 -> y);
#   - for each group used below: a base and strong generating set (Schreier-Sims), the order of G (compared with the
#     stored order), and the hypotheses of the F_2 argument below: rho lies in G and commutes with every generator;
#   - the imaginary quadratic subfields K of the Galois closure, as characters eps: G -> {+1,-1} with
#     eps(rho) = -1.  On each G-orbit O of X, K lies in the factor E_O if and only if the F_2-linear system
#         p(s x) = p(x) + e(s)  (s a generator, x in O),   p(rho x) = p(x) + 1,   p(x0) = 0
#     has a solution; e = eps on the generators, and {x : p(x) = 0} is the K-block of O containing x0.
#     Orbits with the same e belong to the same K.  Each orbit's block is the one containing its smallest point, so the
#     blocks of different orbits are not oriented by one embedding of K; this is harmless, since c_i runs over a
#     symmetric range and w_i and t_i both change sign under B_i <-> rho(B_i);
#   - for each K, the balanced Weil characters w = sum_i c_i w_i (c_i integers, D = sum |c_i| dim A_i <= 6, found by
#     a recursive search), where w_i is the character of the K-block B_i and balance means sum_i c_i t_i = 0,
#     t_i = 2|Phi & B_i| - |B_i|;
#   - their kinds (sixfold note, Sec. 5): (a) D = 2; (b) D = 4; (d) D = 6 and a single simple factor used once;
#     (c) D = 6 and a factor of odd dimension among those with c_i != 0; every other B is "hard".  For reduced pairs
#     with g <= 7 no hard B can occur: it would need a factor of dimension 2 whose CM field contains K; that field is
#     then biquadratic, and no CM type of a biquadratic field is primitive;
#   - lattice generation: the characters generate Lambda_U iff their echelon basis has d rows and spans a saturated
#     lattice, i.e. the unimodular column reduction of the basis to [H | 0] has |det H| = 1 (Lambda_U is saturated,
#     and every w is checked to lie in it).
# The program exits with status 1 if a recomputed d, reduced flag or group order differs from the stored one.
#   usage:  python verify_blocks.py <enum_g<g>.json> <out.json>
import itertools, json, sys, collections

DMAX = 6


# ---------------------------------------------------------------- permutations (tuples p with p[x] = image of x)
def mul(p, q):
    """p first, then q."""
    return tuple(q[x] for x in p)


def inv(p):
    r = [0] * len(p)
    for x, y in enumerate(p):
        r[y] = x
    return tuple(r)


def orbit_trans(b, S, ident):
    """{x: u} with u in <S> mapping b to x."""
    T = {b: ident}
    queue = [b]
    for x in queue:
        for s in S:
            y = s[x]
            if y not in T:
                T[y] = mul(T[x], s)
                queue.append(y)
    return T


def sift(h, base, trans, start=0):
    """Strip h through the levels start, start+1, ...; returns (residue, level of drop-out or len(base))."""
    for i in range(start, len(base)):
        x = h[base[i]]
        if x not in trans[i]:
            return h, i
        h = mul(h, inv(trans[i][x]))
    return h, len(base)


def schreier_sims(gens, n):
    """Deterministic Schreier-Sims: base and transversals of the stabilizer chain of <gens>."""
    ident = tuple(range(n))
    gens = [g0 for g0 in (tuple(s) for s in gens) if g0 != ident]
    base = []
    for s in gens:
        if all(s[b] == b for b in base):
            base.append(next(x for x in range(n) if s[x] != x))
    S = [[s for s in gens if all(s[base[j]] == base[j] for j in range(i))] for i in range(len(base))]
    T = [orbit_trans(base[i], S[i], ident) for i in range(len(base))]
    i = len(base) - 1
    while i >= 0:
        new = None
        for x in list(T[i]):
            for s in S[i]:
                sg = mul(mul(T[i][x], s), inv(T[i][s[x]]))
                if sg == ident:
                    continue
                h, j = sift(sg, base, T, i + 1)
                if j < len(base) or h != ident:
                    new = (h, j)
                    break
            if new:
                break
        if new is None:
            i -= 1
            continue
        h, j = new
        if j == len(base):
            base.append(next(x for x in range(n) if h[x] != x))
            S.append([])
            T.append({base[-1]: ident})
        for lev in range(i + 1, j + 1):
            S[lev].append(h)
            T[lev] = orbit_trans(base[lev], S[lev], ident)
        i = j
    return base, T


# ---------------------------------------------------------------- linear algebra over F_2 and Z
def gf2_solutions(eqs, nvars):
    """All solutions of a GF(2) system; eqs = list of (mask, rhs).  Returns a list of bit lists."""
    basis = {}
    for m, r in eqs:
        while m:
            lb = m.bit_length() - 1
            if lb in basis:
                bm, br = basis[lb]
                m ^= bm
                r ^= br
            else:
                basis[lb] = (m, r)
                break
        else:
            if r:
                return []
    free = [v for v in range(nvars) if v not in basis]
    sols = []
    for bits in itertools.product((0, 1), repeat=len(free)):
        val = [0] * nvars
        for v, b in zip(free, bits):
            val[v] = b
        for lb in sorted(basis):
            m, r = basis[lb]
            rest = m ^ (1 << lb)
            par = r
            v = 0
            while rest:
                if rest & 1:
                    par ^= val[v]
                rest >>= 1
                v += 1
            val[lb] = par
        sols.append(val)
    return sols


def col_reduce(rows, n):
    """Unimodular column reduction of an integer matrix.  Returns (M, U, pivots): M = rows * U is in column echelon
    form, U is unimodular, pivots = the pivot entries (one per row with a pivot), and the columns of M after the
    last pivot column are zero."""
    M = [list(r) for r in rows]
    U = [[int(i == j) for j in range(n)] for i in range(n)]   # U[i][j]: column j of the transform

    def colop(j, k, q):          # column j -= q * column k
        for row in M:
            row[j] -= q * row[k]
        for row in U:
            row[j] -= q * row[k]

    def colswap(j, k):
        for row in M:
            row[j], row[k] = row[k], row[j]
        for row in U:
            row[j], row[k] = row[k], row[j]

    pc = 0
    pivots = []
    for r in range(len(M)):
        if pc >= n:
            break
        while True:
            nz = [j for j in range(pc, n) if M[r][j] != 0]
            if not nz:
                break
            k = min(nz, key=lambda j: abs(M[r][j]))
            if k != pc:
                colswap(k, pc)
            for j in range(pc + 1, n):
                if M[r][j]:
                    colop(j, pc, M[r][j] // M[r][pc])
            if all(M[r][j] == 0 for j in range(pc + 1, n)):
                break
        if any(M[r][j] for j in range(pc, n)):
            pivots.append(M[r][pc])
            pc += 1
    assert all(M[r][j] == 0 for r in range(len(M)) for j in range(pc, n))
    return M, U, pivots


def int_kernel(rows, n):
    """Z-basis of {x in Z^n : r.x = 0 for all r in rows}: the last n - rank columns of U."""
    M, U, pivots = col_reduce(rows, n)
    pc = len(pivots)
    return [[U[i][j] for i in range(n)] for j in range(pc, n)]


def egcd(a, b):
    if b == 0:
        return (abs(a), (1 if a > 0 else -1), 0)
    g0, x, y = egcd(b, a % b)
    return (g0, y, x - (a // b) * y)


def echelon(vecs, n):
    """Echelon Z-basis of the lattice spanned by vecs (unimodular row operations)."""
    A = [list(v) for v in vecs if any(v)]
    out = []
    for c in range(n):
        piv, rest = None, []
        for r in A:
            if r[c] == 0:
                rest.append(r)
            elif piv is None:
                piv = r
            else:
                a, b = piv[c], r[c]
                g0, x, y = egcd(a, b)
                newp = [x * u + y * v for u, v in zip(piv, r)]
                other = [(b // g0) * u - (a // g0) * v for u, v in zip(piv, r)]
                assert newp[c] == g0 and other[c] == 0
                piv = newp
                if any(other):
                    rest.append(other)
        if piv is not None:
            out.append(piv)
        A = rest
    assert not any(any(r) for r in A)
    return out


def generates(vecs, d, n):
    """True iff vecs (all in the saturated lattice Lambda_U of rank d) generate Lambda_U: the echelon basis B has d
    rows, and B * U = [H | 0] with |det H| = product of the pivots = 1 (the row lattice of B is saturated)."""
    B = echelon(vecs, n)
    if len(B) != d:
        return False
    M, U, pivots = col_reduce(B, n)
    assert len(pivots) == d
    return all(abs(p) == 1 for p in pivots)


# ---------------------------------------------------------------- G-sets
def orbits_of(gens, n):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for s in gens:
        for x in range(n):
            a, b = find(x), find(s[x])
            if a != b:
                parent[a] = b
    orbs = collections.defaultdict(list)
    for x in range(n):
        orbs[find(x)].append(x)
    return sorted(orbs.values())


def min_block_system(O, gens, a, b):
    """The finest G-invariant partition of the orbit O in which a and b lie in one block (Atkinson); returns the
    list of blocks."""
    parent = {x: x for x in O}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx == ry:
            return False
        parent[rx] = ry
        return True
    union(a, b)
    queue = [(a, b)]
    while queue:
        u, v = queue.pop()
        for s in gens:
            if union(s[u], s[v]):
                queue.append((s[u], s[v]))
    blocks = collections.defaultdict(list)
    for x in O:
        blocks[find(x)].append(x)
    return sorted(blocks.values())


def gset_isos(O, O2, gens):
    """All G-set isomorphisms O -> O2, as dicts."""
    out = []
    x0 = O[0]
    for y in O2:
        f = {x0: y}
        queue = [x0]
        ok = True
        for x in queue:
            for s in gens:
                sx, sfx = s[x], s[f[x]]
                if sx in f:
                    if f[sx] != sfx:
                        ok = False
                        break
                else:
                    f[sx] = sfx
                    queue.append(sx)
            if not ok:
                break
        if ok and len(f) == len(O) and set(f.values()) == set(O2):
            out.append(f)
    return out


def group_data(gens, n):
    orbs = orbits_of(gens, n)
    systems = []
    for O in orbs:
        sy = {tuple(map(tuple, min_block_system(O, gens, O[0], y))) for y in O[1:]}
        systems.append([[set(b) for b in s] for s in sy])
    isos = []
    for i, j in itertools.combinations(range(len(orbs)), 2):
        if len(orbs[i]) == len(orbs[j]):
            isos.extend((i, j, f) for f in gset_isos(orbs[i], orbs[j], gens))
    return dict(orbs=orbs, systems=systems, isos=isos)


def is_reduced(P, gd):
    for O, sys_O in zip(gd['orbs'], gd['systems']):
        for blocks in sys_O:
            if all(b <= P or not (b & P) for b in blocks):
                return False          # Phi & O is a union of blocks: not primitive
    for i, j, f in gd['isos']:
        if all((x in P) == (f[x] in P) for x in gd['orbs'][i]):
            return False              # two isomorphic orbits with CM types: isogenous factors
    return True


# ---------------------------------------------------------------- Weil characters
def balanced_vectors(sz, tv, dmax):
    """All integer vectors c != 0 with D = sum |c_i| sz_i <= dmax and sum c_i t_i = 0 (recursive search)."""
    out = []
    cur = []

    def rec(i, D, bal):
        if i == len(sz):
            if D > 0 and bal == 0:
                out.append(tuple(cur))
            return
        m = (dmax - D) // sz[i]
        for c in range(-m, m + 1):
            cur.append(c)
            rec(i + 1, D + abs(c) * sz[i], bal + c * tv[i])
            cur.pop()
    rec(0, 0, 0)
    return out


def kind(cv, sz):
    D = sum(abs(c) * s for c, s in zip(cv, sz))
    if D == 2:
        return 'a'
    if D == 4:
        return 'b'
    assert D == 6
    if sum(abs(c) for c in cv) == 1:
        return 'd'
    if any(c != 0 and s % 2 == 1 for c, s in zip(cv, sz)):
        return 'c'
    return 'hard'


def main():
    D = json.load(open(sys.argv[1]))
    g = D['g']
    n = 2 * g
    ident = tuple(range(n))
    rho = tuple((x + g) % n for x in range(n))
    groups = {G['gid']: G for G in D['groups']}
    gcache = {}
    cnt = collections.Counter()
    tab = collections.defaultdict(collections.Counter)
    records = []
    bad = []
    for c in D['cases']:
        G = groups[c['gid']]
        gens = [tuple(s) for s in G['gens']]
        Phi = c['Phi']
        P = set(Phi)
        assert len(Phi) == g and all((j in P) != (j + g in P) for j in range(g))
        # G-orbit of Phi and Lambda_U
        start = frozenset(Phi)
        seen, todo = {start}, [start]
        while todo:
            T = todo.pop()
            for s in gens:
                U2 = frozenset(s[x] for x in T)
                if U2 not in seen:
                    seen.add(U2)
                    todo.append(U2)
        mus = [[1 if j in T else -1 for j in range(g)] for T in seen]
        LU = int_kernel(mus, g)
        d = len(LU)
        cnt['cases'] += 1
        if d == c['d']:
            cnt['d agrees with stored'] += 1
        else:
            bad.append(('d', c['gid'], Phi, d, c['d']))
        # reduced?
        if c['gid'] not in gcache:
            gcache[c['gid']] = group_data(gens, n)
        gd = gcache[c['gid']]
        red = is_reduced(P, gd)
        cnt['reduced (recomputed)'] += red
        if red == c['reduced']:
            cnt['reduced agrees with stored'] += 1
        else:
            bad.append(('reduced', c['gid'], Phi, red, c['reduced']))
        if not (red and d >= 1):
            continue
        # stored LU spans the same lattice
        sLU = c['LU']
        assert all(sum(a * b for a, b in zip(w, mu)) == 0 for w in sLU for mu in mus)
        cnt['stored LU is a basis of Lambda_U'] += generates(sLU, d, g) and len(sLU) == d
        # hypotheses of the F_2 argument, and the imaginary quadratic characters, once per group
        if 'Ks' not in gd:
            orbs = gd['orbs']
            assert all(set(rho[x] for x in O) == set(O) for O in orbs)
            assert all(mul(s, rho) == mul(rho, s) for s in gens)
            base, trans = schreier_sims(gens, n)
            order = 1
            for t in trans:
                order *= len(t)
            h, lev = sift(rho, base, trans)
            assert lev == len(base) and h == ident, 'rho not in G'
            cnt['groups used'] += 1
            if order == G['order']:
                cnt['group order agrees with stored'] += 1
            else:
                bad.append(('order', c['gid'], order, G['order']))
            k = len(gens)
            Ks = collections.defaultdict(dict)        # e-vector -> {orbit index: block}
            for oi, O in enumerate(orbs):
                idx = {x: i for i, x in enumerate(O)}
                m = len(O)
                eqs = []
                for si, s in enumerate(gens):
                    for x in O:
                        eqs.append(((1 << idx[s[x]]) ^ (1 << idx[x]) ^ (1 << (m + si)), 0))
                for x in O:
                    eqs.append(((1 << idx[rho[x]]) ^ (1 << idx[x]), 1))
                eqs.append((1 << idx[O[0]], 0))
                for sol in gf2_solutions(eqs, m + k):
                    e = tuple(sol[m:])
                    Ks[e][oi] = frozenset(x for x in O if sol[idx[x]] == 0)
            gd['Ks'] = dict(Ks)
            cnt['K-blocks on orbits of size 4'] += sum(1 for bl in Ks.values() for o in bl if len(orbs[o]) == 4)
        orbs, Ks = gd['orbs'], gd['Ks']
        part = tuple(sorted((len(O) // 2 for O in orbs), reverse=True))
        ws = []          # (w, strict, kind)
        for e, blocks in Ks.items():
            oids = sorted(blocks)
            sz = [len(orbs[o]) // 2 for o in oids]
            wv = []
            tv = []
            for o in oids:
                B = blocks[o]
                wv.append([1 if j in B else (-1 if j + g in B else 0) for j in range(g)])
                tv.append(2 * len(P & B) - len(B))
            for cv in balanced_vectors(sz, tv, DMAX):
                w = [sum(x * wi[j] for x, wi in zip(cv, wv)) for j in range(g)]
                assert any(w) and all(sum(a * b for a, b in zip(w, mu)) == 0 for mu in mus)
                ws.append((w, all(abs(x) <= 1 for x in cv), kind(cv, sz)))
        strict_ok = generates([w for w, st, kd in ws if st], d, g)
        all_ok = generates([w for w, st, kd in ws], d, g)
        ab_ok = generates([w for w, st, kd in ws if kd in ('a', 'b')], d, g)
        abd_ok = generates([w for w, st, kd in ws if kd in ('a', 'b', 'd')], d, g)
        easy_ok = generates([w for w, st, kd in ws if kd != 'hard'], d, g)
        hard_any = any(kd == 'hard' for w, st, kd in ws)
        cnt['reduced d>=1'] += 1
        cnt['strict flag agrees'] += (strict_ok == c['passed'])
        if 'mult' in c:
            md = c['mult']['Dmin']
            cnt['mult Dmin<=6 agrees'] += (all_ok == (md is not None and md <= DMAX))
        cnt['pass (D<=6)'] += all_ok
        cnt['pass strict'] += strict_ok
        cnt['pass by kinds (a),(b)'] += ab_ok
        cnt['pass by kinds (a),(b),(d)'] += abd_ok
        cnt['pass by kinds (a)-(d)'] += easy_ok
        cnt['pass but not by kinds (a)-(d)'] += (all_ok and not easy_ok)
        cnt['hard B present'] += hard_any
        key = (part, d, 'strict' if strict_ok else ('mult<=6' if all_ok else 'fail'))
        tab[key]['cases'] += 1
        tab[key]['easy'] += easy_ok
        records.append(dict(gid=c['gid'], Phi=Phi, d=d, part=part, strict=strict_ok, all=all_ok, easy=easy_ok,
                            hard_present=hard_any, nK=len(Ks)))
    print('g=%d' % g)
    for k2 in ('cases', 'd agrees with stored', 'reduced (recomputed)', 'reduced agrees with stored', 'reduced d>=1',
               'groups used', 'group order agrees with stored', 'K-blocks on orbits of size 4',
               'stored LU is a basis of Lambda_U', 'strict flag agrees', 'mult Dmin<=6 agrees', 'pass (D<=6)',
               'pass strict', 'pass by kinds (a),(b)', 'pass by kinds (a),(b),(d)', 'pass by kinds (a)-(d)',
               'pass but not by kinds (a)-(d)', 'hard B present'):
        print('  %-36s %d' % (k2, cnt[k2]))
    print('  by (factor dims, d, tier): cases ; generated by kinds (a)-(d)')
    for key in sorted(tab, key=lambda t: (-t[0][0], t)):
        print('    %-16s d=%d %-8s %4d ; %4d' % (key[0], key[1], key[2], tab[key]['cases'], tab[key]['easy']))
    print('  disagreements with stored data: %d' % len(bad))
    for b in bad[:20]:
        print('   ', b)
    json.dump(records, open(sys.argv[2], 'w'), indent=0)
    if bad:
        sys.exit(1)


if __name__ == '__main__':
    main()
