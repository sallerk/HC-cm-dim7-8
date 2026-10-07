# verify_ext.py -- INDEPENDENT VERIFIER (pure Python) for ext_weil.py (Weil search beyond D=16) on the d>=3 classes.
# Own element enumeration, ALL subgroups (joins of cyclic subgroups, from verify_general.py), left-coset actions,
# own orbit / lattice code; budget DFS over multiplicities.  Compares the per-level data (D, rank, saturation index).
#   usage: python verify_ext.py Dmax ext_weil_D{Dmax}[...].json
import json, sys, time
from verify_enum import hnf, sat_index, rank, wvec
from verify_general import elements, all_subgroups


def levels_for(G, g, Phi, d, LU, Dmax):
    n = 2 * g
    gens = G['gens']
    el, idx = elements(gens, n)
    N = len(el)
    rho = tuple((i + g) % n for i in range(n))
    r = idx[rho]
    subs, mul = all_subgroups(el, idx, n)
    invs = [next(y for y in range(N) if mul[x][y] == 0) for x in range(N)]
    Ps = set(Phi)
    gi = [idx[tuple(p)] for p in gens]
    O = list(range(n))
    byD = {}
    seenclass = set()
    for M in subs:
        if r in M or N // len(M) > Dmax:
            continue
        conj = min(tuple(sorted(mul[mul[x][m]][invs[x]] for m in M)) for x in range(N))
        if conj in seenclass:
            continue
        seenclass.add(conj)
        cos = {}
        reps = []
        for x in range(N):
            if x in cos:
                continue
            c = len(reps)
            reps.append(x)
            for m in M:
                cos[mul[x][m]] = c
        k2 = len(reps)
        act = [[cos[mul[gg][reps[c]]] for c in range(k2)] for gg in gi]
        ract = [cos[mul[r][reps[c]]] for c in range(k2)]
        units = []
        seen = set()
        done = set()
        for x in O:
            for t in range(k2):
                if (x, t) in seen:
                    continue
                om = {(x, t)}
                st = [(x, t)]
                while st:
                    a, b = st.pop()
                    for p, q in zip(gens, act):
                        y = (p[a], q[b])
                        if y not in om:
                            om.add(y); st.append(y)
                seen |= om
                om = frozenset(om)
                if om in done:
                    continue
                done.add(om); done.add(frozenset((a, ract[b]) for a, b in om))
                wgt = len(om) // 2
                if wgt > Dmax:
                    continue
                Fs = [[a for a, b in om if b == tau] for tau in range(k2)]
                beta = [wvec(F, g) for F in Fs]
                if not any(any(b) for b in beta):
                    continue
                s = [2 * len(Ps & set(F)) - len(F) for F in Fs]
                units.append((wgt, beta, s))
        if not units:
            continue
        nu = len(units)
        def rec(j, budget, ssum, chis, D):
            if j == nu:
                if D and not any(ssum):
                    for w in chis:
                        if any(w):
                            byD.setdefault(D, {}).setdefault(k2, set()).add(tuple(w))
                return
            wgt, beta, s = units[j]
            mm = budget // wgt
            for c in range(-mm, mm + 1):
                if c == 0:
                    rec(j + 1, budget, ssum, chis, D)
                else:
                    rec(j + 1, budget - abs(c) * wgt, [a + c * b for a, b in zip(ssum, s)],
                        [[a + c * b for a, b in zip(ch, be)] for ch, be in zip(chis, beta)], D + abs(c) * wgt)
        rec(0, Dmax, [0] * k2, [[0] * g for _ in range(k2)], 0)
    levels = []
    acc = []
    used = set()
    for D in sorted(byD):
        for k2, ws in byD[D].items():
            for w in ws:
                if rank(LU + [list(w)]) != d:
                    raise RuntimeError("Weil character outside Lambda_U")
            acc += [list(w) for w in ws]
            used.add(k2)
        L = hnf(acc)
        sat = sat_index(L) if len(L) == d else None
        levels.append(dict(D=D, rank=len(L), satindex=sat, indices=sorted(used)))
    return levels


if __name__ == '__main__':
    Dmax = int(sys.argv[1])
    recs = json.load(open(sys.argv[2]))
    T = json.load(open('trans_g8.json'))
    byTI = {(G['TI'], G['j']): G for G in T['groups']}
    t0 = time.time()
    agree = 0
    for rcd in recs:
        G = byTI[(rcd['TI'], rcd['j'])]
        lv = levels_for(G, 8, rcd['Phi'], rcd['d'], rcd['LU'], Dmax)
        mine = [(l['D'], l['rank'], l['satindex'], tuple(l['indices'])) for l in lv]
        prod_ = [(l['D'], l['rank'], l['satindex'], tuple(l['indices'])) for l in rcd['levels']]
        ok = (mine == prod_)
        agree += ok
        print(rcd['TI'], rcd['j'], rcd['Phi'], 'AGREE' if ok else 'DISAGREE', mine, '' if ok else prod_, '%.1fs' % (time.time() - t0), flush=True)
    print("extended-Weil verification (Dmax=%d): %d/%d agree ; %.1fs" % (Dmax, agree, len(recs), time.time() - t0))
