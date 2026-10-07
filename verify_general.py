# verify_general.py -- INDEPENDENT VERIFIER (pure Python) for the general CM-field Weil analysis (general_weil2.py)
# on all reduced (G,Phi) that FAIL the strict key test.  Own subgroup enumeration (joins of cyclic subgroups),
# own coset actions, own orbit computations, own lattice routines (from verify_enum).
# Checks: minimal dimension D at which Weil characters (any CM field K' with [K':Q] <= 16, any multiplicities on
# sub-products of powers) generate Lambda_U, searched up to Dmax = 8.  (Units have dimension >= k = [K':Q]/2, so
# [K':Q] <= 16 is exhaustive for D <= 8.)
import json, sys, os, itertools, time
from verify_enum import hnf, sat_index, rank, wvec, set_orbit

HERE = os.path.dirname(os.path.abspath(__file__))


def elements(gens, n):
    e = tuple(range(n))
    el = [e]
    seen = {e: 0}
    i = 0
    while i < len(el):
        x = el[i]
        for p in gens:
            y = tuple(p[x[j]] for j in range(n))   # p o x
            if y not in seen:
                seen[y] = len(el)
                el.append(y)
        i += 1
    return el, seen


def all_subgroups(el, idx, n):
    N = len(el)
    mul = [[idx[tuple(a[b[j]] for j in range(n))] for b in el] for a in el]   # mul[a][b] = a o b
    def close(gs):
        S = {0}
        fr = [0]
        S |= set(gs)
        fr = list(S)
        while fr:
            nf = []
            for a in fr:
                for b in list(gs):
                    c = mul[a][b]
                    if c not in S:
                        S.add(c); nf.append(c)
            fr = nf
        return frozenset(S)
    cyc = set(close([a]) for a in range(N))
    subs = set(cyc)
    frontier = set(cyc)
    while frontier:
        new = set()
        for H in frontier:
            for C in cyc:
                if C <= H:
                    continue
                J = close(list(H | C))
                if J not in subs:
                    new.add(J)
        subs |= new
        frontier = new
    return subs, mul


def analyse(G, g, Phi, d, LU, Dmax=8, maxindex=16):
    n = 2 * g
    gens = G['gens']
    el, idx = elements(gens, n)
    N = len(el)
    assert N == G['order']
    rho = tuple((i + g) % n for i in range(n))
    r = idx[rho]
    subs, mul = all_subgroups(el, idx, n)
    orbs = G['orbits']
    Ps = set(Phi)
    gi = [idx[tuple(p)] for p in gens]
    allv = []
    # conjugacy classes not needed: using every subgroup is harmless (conjugate M give the same units up to relabelling)
    seenclass = set()
    for M in subs:
        if r in M or N // len(M) > maxindex:
            continue
        # conjugacy-class dedupe (cheap): canonical = min over conjugates of sorted tuple
        conj = min(tuple(sorted(mul[mul[x][m]][inv] for m in M)) for x in range(N)
                   for inv in [next(y for y in range(N) if mul[x][y] == 0)])
        if conj in seenclass:
            continue
        seenclass.add(conj)
        # left cosets xM, left action
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
        taus = []
        for t in range(k2):
            if t not in taus and ract[t] not in taus:
                taus.append(t)
        units = []
        for O in orbs:
            seen = set()
            oms = []
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
                    oms.append(frozenset(om))
            done = set()
            for om in oms:
                if om in done:
                    continue
                done.add(om)
                done.add(frozenset((a, ract[b]) for a, b in om))
                m = len(om) // len(O)
                Fs = [[a for a, b in om if b == t] for t in taus]
                beta = [wvec(F, g) for F in Fs]
                if not any(any(b) for b in beta):
                    continue
                s = [2 * len(Ps & set(F)) - len(F) for F in Fs]
                units.append((m * len(O) // 2, beta, s))
        if not units:
            continue
        ws = [u[0] for u in units]
        ranges = [range(-(Dmax // w), Dmax // w + 1) for w in ws]
        for c in itertools.product(*ranges):
            D = sum(abs(ci) * w for ci, w in zip(c, ws))
            if D == 0 or D > Dmax:
                continue
            if any(sum(ci * u[2][ti] for ci, u in zip(c, units)) for ti in range(len(taus))):
                continue
            for ti in range(len(taus)):
                w = [sum(ci * u[1][ti][j] for ci, u in zip(c, units)) for j in range(g)]
                if any(w):
                    if rank(LU + [w]) != d:
                        raise RuntimeError("Weil character outside Lambda_U")
                    allv.append((D, N // len(M), w))
    Dmin = None
    for D in sorted(set(z[0] for z in allv)):
        L = hnf([z[2] for z in allv if z[0] <= D])
        if len(L) == d and sat_index(L) == 1:
            Dmin = D
            break
    fields = sorted(set((z[1], z[0]) for z in allv if Dmin is not None and z[0] <= Dmin))
    return Dmin, fields


if __name__ == '__main__':
    t0 = time.time()
    gw = json.load(open(os.path.join(HERE, 'general_weil2_out_D12.json')))
    tot = agree = 0
    for g in [5, 6]:
        D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g)))
        groups = {G['gid']: G for G in D['groups']}
        prod_res = {(r['gid'], tuple(r['Phi'])): r for r in gw[str(g)]}
        for c in D['cases']:
            if not (c['reduced'] and not c['passed']):
                continue
            Dmin, fields = analyse(groups[c['gid']], g, c['Phi'], c['d'], c['LU'])
            pr = prod_res[(c['gid'], tuple(c['Phi']))]
            ok = (Dmin == pr['Dmin'])
            tot += 1
            agree += ok
            print(g, c['gid'], c['Phi'], 'verifier Dmin=%s general_weil2 genDmin=%s %s fields(index,D)=%s' %
                  (Dmin, pr['Dmin'], 'AGREE' if ok else 'DISAGREE', fields), flush=True)
    print("general-Weil verification: %d cases, %d agree ; %.1fs" % (tot, agree, time.time() - t0))
