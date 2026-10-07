# verify_generalC3.py -- INDEPENDENT VERIFIER (pure Python) for the general Weil tier of the (6,1) Stage C cases
# (orbits O1 of size 12 and O2 of size 2), without enumerating subgroups.  For D <= Dmax < 12:
#   a unit on O1 x Z has |Omega| = 12m, weight 6m, so weight <= Dmax forces m = 1: Omega ~ O1 and Z is a quotient of O1;
#   a combination using only O2-units has characters supported on the O2 coordinate, and these must be 0 when
#   e_{O2} is not in Lambda_U (checked).  So Z runs over the rho-free block-system quotients of O1 (|Z| <= Dmax).
# Units on O1 x Z and O2 x Z (all G-orbits), budget DFS over multiplicities, own lattice code.
import json, sys, time
from verify_enum import hnf, sat_index, rank, wvec
from verify_trans import all_block_systems
from verify_generalC import dmin

g = 7; n = 14
D = json.load(open('enum_g7.json'))
groups = {G['gid']: G for G in D['groups']}
tot = agree = 0
t0 = time.time()
for c in D['cases']:
    if 'general' not in c:
        continue
    G = groups[c['gid']]
    if sorted(len(o) for o in G['orbits']) != [2, 12]:
        continue
    O1 = next(o for o in G['orbits'] if len(o) == 12); O2 = next(o for o in G['orbits'] if len(o) == 2)
    Dm = c['general']['Dmin'] if c['general']['Dmin'] is not None else c['general']['Dmax']
    assert Dm < 12
    e2 = [0] * g; e2[min(O2) % g] = 1
    assert rank(c['LU'] + [e2]) > c['d'], "O2-only characters could matter"
    gens = G['gens']; rho = [(i + g) % n for i in range(n)]
    Ps = set(c['Phi'])
    byD = {}
    for part in all_block_systems(gens, O1, O1[0]):
        nz = len(part)
        if nz > Dm or nz < 2:
            continue
        bl = {x: i for i, B in enumerate(part) for x in B}
        Zg = [[bl[p[B[0]]] for B in part] for p in gens]
        Zr = [bl[rho[B[0]]] for B in part]
        if any(Zr[z] == z for z in range(nz)):
            continue
        units = []
        for O in (O1, O2):
            seen = set(); done = set()
            for x in O:
                for t in range(nz):
                    if (x, t) in seen:
                        continue
                    om = {(x, t)}; st = [(x, t)]
                    while st:
                        a, b = st.pop()
                        for p, q in zip(gens, Zg):
                            y = (p[a], q[b])
                            if y not in om:
                                om.add(y); st.append(y)
                    seen |= om; om = frozenset(om)
                    if om in done:
                        continue
                    done.add(om); done.add(frozenset((a, Zr[b]) for a, b in om))
                    if len(om) // 2 > Dm:
                        continue
                    Fs = [[a for a, b in om if b == tau] for tau in range(nz)]
                    beta = [wvec(F, g) for F in Fs]
                    if any(any(b) for b in beta):
                        units.append((len(om) // 2, beta, [2 * len(Ps & set(F)) - len(F) for F in Fs]))
        nu = len(units)
        def rec(j, budget, ssum, chis, Dd):
            if j == nu:
                if Dd and not any(ssum):
                    for w in chis:
                        if any(w):
                            byD.setdefault(Dd, set()).add(tuple(w))
                return
            wgt, beta, s = units[j]; mm = budget // wgt
            for cc in range(-mm, mm + 1):
                if cc == 0:
                    rec(j + 1, budget, ssum, chis, Dd)
                else:
                    rec(j + 1, budget - abs(cc) * wgt, [a + cc * b for a, b in zip(ssum, s)],
                        [[a + cc * b for a, b in zip(ch, be)] for ch, be in zip(chis, beta)], Dd + abs(cc) * wgt)
        rec(0, Dm, [0] * nz, [[0] * g for _ in range(nz)], 0)
    Dlat, Dq = dmin(byD, c['LU'], c['d'])
    ok = Dlat == c['general']['Dmin']; tot += 1; agree += ok
    print(7, c['gid'], c['Phi'], '|G|=%d verifier Dmin=%s producer Dmin=%s %s' % (G['order'], Dlat, c['general']['Dmin'], 'AGREE' if ok else 'DISAGREE'), flush=True)
print('verify_generalC3 ((6,1) cases, block-quotient method): %d cases, %d agree ; %.1fs' % (tot, agree, time.time() - t0))
