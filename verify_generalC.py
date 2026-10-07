# verify_generalC.py -- INDEPENDENT VERIFIER (pure Python) for the Stage C general CM-field tier (stageC.py 'general'):
# for every reduced strict failure at g=7, recompute the minimal dimension D at which Weil characters (any CM field
# K' <-> M <= G with rho not in M, [G:M] <= D, any multiplicities) generate Lambda_U as a lattice, by: own element
# enumeration, ALL subgroups (joins of cyclic subgroups), left-coset actions, budget DFS over multiplicities.
#   usage: python verify_generalC.py [maxorder]
import json, sys, time
from verify_enum import hnf, sat_index, rank, wvec
from verify_general import elements, all_subgroups


def chars_by_D(G, g, Phi, Dmax):
    n = 2 * g
    gens = G['gens']
    el, idx = elements(gens, n)
    N = len(el)
    assert N == G['order']
    rho = tuple((i + g) % n for i in range(n))
    r = idx[rho]
    subs, mul = all_subgroups(el, idx, n)
    invs = [None] * N
    for x in range(N):
        for y in range(N):
            if mul[x][y] == 0:
                invs[x] = y
                break
    Ps = set(Phi)
    gi = [idx[tuple(p)] for p in gens]
    orbs = G['orbits']
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
        for O in orbs:
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
                            byD.setdefault(D, set()).add(tuple(w))
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
    return byD


def dmin(byD, LU, d):
    acc = []
    Dlat = Dq = None
    for D in sorted(byD):
        for w in byD[D]:
            if rank(LU + [list(w)]) != d:
                raise RuntimeError("Weil character outside Lambda_U")
        acc += [list(w) for w in byD[D]]
        L = hnf(acc)
        if Dq is None and len(L) == d:
            Dq = D
        if len(L) == d and sat_index(L) == 1:
            Dlat = D
            break
    return Dlat, Dq


if __name__ == '__main__':
    maxorder = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    D = json.load(open('enum_g7.json'))
    groups = {G['gid']: G for G in D['groups']}
    t0 = time.time()
    tot = agree = skipped = 0
    out = open('log_verify_generalC.txt', 'w')
    def P(*a):
        print(*a, flush=True); print(*a, file=out); out.flush()
    for c in D['cases']:
        if 'general' not in c:
            continue
        G = groups[c['gid']]
        tot += 1
        if G['order'] > maxorder:
            skipped += 1
            P(7, c['gid'], c['Phi'], 'SKIPPED |G|=%d' % G['order'])
            continue
        Dm = c['general']['Dmin'] if c['general']['Dmin'] is not None else c['general']['Dmax']
        byD = chars_by_D(G, 7, c['Phi'], Dm)
        Dlat, Dq = dmin(byD, c['LU'], c['d'])
        ok = (Dlat == c['general']['Dmin'])
        agree += ok
        P(7, c['gid'], c['Phi'], '|G|=%d verifier Dmin=%s producer Dmin=%s %s' % (G['order'], Dlat, c['general']['Dmin'], 'AGREE' if ok else 'DISAGREE'))
    P("general-Weil verification g=7: %d cases, %d agree, %d skipped (|G|>%d) ; %.1fs" % (tot, agree, skipped, maxorder, time.time() - t0))
