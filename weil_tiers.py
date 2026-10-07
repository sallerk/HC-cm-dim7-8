# weil_tiers.py -- main computation (Sage+GAP): general CM-field Weil characters on sub-products of powers of A
# (logic of general_weil2.py, CM-sixfold enumeration).  For K' <-> M <= G (rho not in M, [G:M] = 2k <= maxindex), a "unit" is a
# G-orbit Omega on O x G/M (K' inside End^0(A_O^m) = M_m(E_O), m = |Omega|/|O|).  For tau in G/M:
#   beta_tau = sum_{sigma : (sigma,tau) in Omega} [sigma],  s_tau = 2|Phi cap F_tau| - |F_tau|,  weight m*|O|/2.
# Combination c (multiplicities, sign = conjugate action) is of Weil type iff sum_j c_j s^j_tau = 0 for all tau;
# Weil characters chi_tau = sum_j c_j beta^j_tau;  dimension D = sum_j |c_j| weight_j.
# A balanced K'-action needs D/k even, so [K':Q] = 2k <= D: maxindex = Dmax is exhaustive.
# Subgroups of index <= maxindex: GAP LowIndexSubgroups (conjugacy-class representatives).
import itertools
from sage.all import libgap, matrix, ZZ, prod
from cmenum import rho_perm, perm_list, _cvecs


def units_for(G, gl, ggens, orbs, M, g, rho):
    n = 2 * g
    idx = int(libgap.Index(G, M))
    cos = libgap.RightCosets(G, M)
    hom = libgap.ActionHomomorphism(G, cos, libgap.OnRight)
    imgs = [[int(libgap.OnPoints(p, libgap.Image(hom, x))) - 1 for p in range(1, idx + 1)] for x in ggens]
    rimg = [int(libgap.OnPoints(p, libgap.Image(hom, rho))) - 1 for p in range(1, idx + 1)]
    taus = []
    for t in range(idx):
        if t not in taus and rimg[t] not in taus:
            taus.append(t)
    units = []
    for i, O in enumerate(orbs):
        seen = set()
        oms = []
        for x in O:
            for t in range(idx):
                if (x, t) in seen:
                    continue
                om = {(x, t)}
                st = [(x, t)]
                while st:
                    a, b = st.pop()
                    for q, qc in zip(gl, imgs):
                        y = (q[a], qc[b])
                        if y not in om:
                            om.add(y)
                            st.append(y)
                seen |= om
                oms.append(frozenset(om))
        done = set()
        for om in oms:
            if om in done:
                continue
            conj = frozenset((a, rimg[b]) for a, b in om)
            done.add(om); done.add(conj)
            m = len(om) // len(O)
            beta, sig = [], []
            for t in taus:
                F = [a for a, b in om if b == t]
                w = [0] * g
                for a in F:
                    if a < g:
                        w[a] += 1
                    else:
                        w[a - g] -= 1
                beta.append(w)
                sig.append(F)
            units.append(dict(orbit=i, m=m, weight=m * len(O) // 2, beta=beta, F=sig))
    return idx, taus, units


def _lat_ok(vs, d):
    if not vs:
        return False, 0, None
    Mx = matrix(ZZ, vs)
    r = Mx.rank()
    if r < d:
        return False, r, None
    ed = [int(x) for x in Mx.elementary_divisors() if x != 0]
    s = int(prod(ed)) if ed else 1
    return s == 1, r, s


def weil_tiers(G, gl_unused, g, Phi, d, LU, Dmax=16, maxindex=16):
    n = 2 * g
    rho = rho_perm(g)
    ggens = list(libgap.GeneratorsOfGroup(G))
    gl = [perm_list(x, n) for x in ggens]
    orbs = sorted([sorted(int(x) - 1 for x in o) for o in libgap.Orbits(G, list(range(1, n + 1)))], key=min)
    Ps = set(Phi)
    allv = []           # (D, index, w, Mid)
    per_field = []
    subs = libgap.LowIndexSubgroups(G, maxindex)
    nsub = 0
    for Mid, M in enumerate(subs):
        if bool(libgap.IsSubset(M, [rho])):
            continue
        nsub += 1
        idx, taus, units = units_for(G, gl, ggens, orbs, M, g, rho)
        k = idx // 2
        for u in units:
            u['s'] = [2 * len(Ps & set(F)) - len(F) for F in u['F']]
        units = [u for u in units if u['weight'] <= Dmax and any(any(b) for b in u['beta'])]
        if not units:
            continue
        best = None
        nchar = 0
        for c in _cvecs([u['weight'] for u in units], Dmax):
            if any(sum(ci * units[j]['s'][ti] for j, ci in enumerate(c)) for ti in range(len(taus))):
                continue
            D = sum(abs(ci) * units[j]['weight'] for j, ci in enumerate(c))
            chis = [[sum(ci * units[j]['beta'][ti][x] for j, ci in enumerate(c)) for x in range(g)] for ti in range(len(taus))]
            chis = [w for w in chis if any(w)]
            if not chis:
                continue
            assert d > 0, "nonzero Weil character with d=0"
            assert matrix(ZZ, LU + chis).rank() == d, "Weil character not in Lambda_U"
            assert D % k == 0 and (D // k) % 2 == 0
            for w in chis:
                allv.append((D, idx, w, Mid))
                nchar += 1
            best = D if best is None else min(best, D)
        if best is not None:
            ok, r, s = _lat_ok([z[2] for z in allv if z[3] == Mid], d)
            per_field.append(dict(Mid=Mid, index=idx, Morder=int(libgap.Size(M)), bestD=best, nchar=nchar,
                                  alone_rank=r, alone_generates=ok,
                                  units=[(u['m'], u['weight']) for u in units]))
    cats = dict(IQ=lambda i: i == 2, le4=lambda i: i <= 4, le8=lambda i: i <= 8, all=lambda i: True)
    Ds = sorted(set(z[0] for z in allv))
    res = dict(nsub_rho_free=nsub, Dmax=Dmax, maxindex=maxindex, Dlevels=Ds)
    for cname, f in cats.items():
        Dlat = Dq = None
        for D in Ds:
            vs = [z[2] for z in allv if z[0] <= D and f(z[1])]
            ok, r, s = _lat_ok(vs, d)
            if Dq is None and r == d:
                Dq = D
            if ok:
                Dlat = D
                break
        res['Dlat_' + cname] = Dlat
        res['Dq_' + cname] = Dq
    # saturation index at each level (all fields)
    lev = []
    for D in Ds:
        vs = [z[2] for z in allv if z[0] <= D]
        ok, r, s = _lat_ok(vs, d)
        lev.append(dict(D=D, rank=r, satindex=s, generates=ok,
                        indices=sorted(set(z[1] for z in allv if z[0] <= D))))
    res['levels'] = lev
    res['per_field'] = per_field
    # minimal index set at the minimal D (all fields): smallest max-index among categories
    return res
