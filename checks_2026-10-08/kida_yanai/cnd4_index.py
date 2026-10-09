# Independent check: defect and index of degeneracy for CM types of G = C_n x D_4 (Galois CM field, [K:Q]=|G|).
# Criterion (Pohlmann/Kubota): a in Z[G] with a(rho g) = -a(g) is a character of U_K trivial on Hg(A)
# iff sum_{g in tau S} a(g) = 0 for all tau in G (left translates).  Index of degeneracy n0 = min ||a||_inf over
# nonzero such a (a Hodge monomial on A^k with multiplicity m, m - m o rho = a, is exceptional iff a != 0).
# Simplicity: right stabilizer {u : S u = S} trivial.
import itertools, sys
from fractions import Fraction
import sympy

def group(n):
    els = [(i, j, e) for i in range(n) for j in range(4) for e in range(2)]
    def mul(x, y):
        return ((x[0]+y[0]) % n, (x[1] + (y[1] if x[2] == 0 else -y[1])) % 4, (x[2]+y[2]) % 2)
    return els, mul

def analyse(n, kmax=4):
    G, mul = group(n)
    idx = {g: t for t, g in enumerate(G)}
    one = (0, 0, 0)
    central_inv = [z for z in G if z != one and mul(z, z) == one and all(mul(z, g) == mul(g, z) for g in G)]
    # automorphisms (for counting classes up to Aut(G) and right translation)
    c, r, s = (1 % n, 0, 0), (0, 1, 0), (0, 0, 1)
    def order(x):
        k, y = 1, x
        while y != one:
            y = mul(y, x); k += 1
        return k
    def word(x, i, j, e, a, b, cc):  # image of c^i r^j s^e
        y = one
        for _ in range(i): y = mul(y, a)
        for _ in range(j): y = mul(y, b)
        for _ in range(e): y = mul(y, cc)
        return y
    autos = []
    for a in G:
        if order(a) != n: continue
        for b in G:
            if order(b) != 4: continue
            for cc in G:
                if order(cc) != 2: continue
                if not (mul(a, b) == mul(b, a) and mul(a, cc) == mul(cc, a)): continue
                if mul(mul(cc, b), cc) != mul(mul(b, b), b): continue
                img = {g: word(g, g[0], g[1], g[2], a, b, cc) for g in G}
                if len(set(img.values())) == len(G):
                    autos.append(img)
    out = {}
    for rho in central_inv:
        R = sorted({min(g, mul(rho, g)) for g in G})
        rows_tau = G
        inv = {g: next(h for h in G if mul(g, h) == one) for g in G}
        results = []
        for bits in itertools.product((0, 1), repeat=len(R)):
            S = frozenset(r_ if b == 0 else mul(rho, r_) for r_, b in zip(R, bits))
            if any(u != one and frozenset(mul(x, u) for x in S) == S for u in G):
                continue  # not simple
            # row entry: 1[r in tau S] - 1[rho r in tau S]; tau S contains r  <=>  tau^{-1} r in S
            M = [[(1 if mul(inv[tau], x) in S else -1) for x in R] for tau in G]
            A = sympy.Matrix(M)
            rk = A.rank()
            d = len(R) - rk
            if d == 0:
                results.append((S, 0, None)); continue
            ns = A.nullspace()
            # rref basis: pick pivot (free) coordinates
            B = sympy.Matrix.hstack(*ns).T  # d x |R|
            Bq, piv = B.rref()
            n0 = None
            for k in range(1, kmax+1):
                found = False
                for coeffs in itertools.product(range(-k, k+1), repeat=d):
                    if all(c_ == 0 for c_ in coeffs): continue
                    v = sum((coeffs[t] * Bq.row(t) for t in range(d)), sympy.zeros(1, len(R)))
                    if all(x.is_integer and abs(x) <= k for x in v):
                        found = True; break
                if found:
                    n0 = k; break
            results.append((S, d, n0))
        # classes up to Aut(G) fixing rho and right translation
        def canon(S):
            best = None
            for au in autos:
                if au[rho] != rho: continue
                T = [au[x] for x in S]
                for u in G:
                    key = tuple(sorted(idx[mul(x, u)] for x in T))
                    if best is None or key < best: best = key
            return best
        cls = {}
        for S, d, n0 in results:
            cls.setdefault(canon(S), (d, n0))
        from collections import Counter
        out[rho] = (Counter((d, n0) for _, d, n0 in results), Counter(cls.values()), len(autos))
    return out

if __name__ == '__main__':
    for n in map(int, sys.argv[1:]):
        res = analyse(n)
        for rho, (cnt_types, cnt_cls, na) in res.items():
            print(f"n={n} |G|={8*n} rho={rho} |Aut G|={na}")
            print("   simple CM types by (defect, n0):", dict(sorted(cnt_types.items(), key=str)))
            print("   classes up to Aut(G,rho) x right translation:", dict(sorted(cnt_cls.items(), key=str)))
