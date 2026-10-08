# used_blocks_check.py -- pure Python (uses verify_enum.py; no Sage/GAP).
# For every REDUCED (G,Phi) with d >= 1 and g <= 6 that passes with imaginary-quadratic Weil characters of
# dimension <= 6 (strict KEY TEST, or Weil characters of sub-products of powers of A with multiplicities, D <= 6),
# test whether the Weil characters of the following "easy" balanced sub-products B of powers of A already generate
# Lambda_U as a lattice (index 1):
#   (i)   dim B <= 4;
#   (ii)  B = a single simple factor, used once;
#   (iii) dim B = 6 and B has a factor of odd dimension (factors counted with multiplicity, e.g. A_E x E_K^2).
# All other B ("hard") have dim B = 6, at least two factor copies, and only factors of even dimension
# (e.g. A_4 x A_2, A_2 x A_2' x A_2'', A_2^3).
# B is encoded as in cmenum.mult_weil: for an imaginary quadratic K (a G-character with rho -> -1) and the orbits O_i
# on which K has a block B_i, a vector c in Z^{orbits}; character w = sum_i c_i [B_i], dim B = sum_i |c_i| g_i,
# balanced iff sum_i c_i t_i = 0 with t_i = 2|Phi cap B_i| - |B_i|.  The strict definition is c_i in {-1,0,1}.
# Everything is recomputed from the generators in enum_g{g}.json; the stored flags (passed, mult.Dmin) and the
# general-CM-field minimal dimension of general_weil2_out_D12.json are only compared.
#   usage:  python used_blocks_check.py > used_blocks_check_out.txt
# g = 7 version (t10_cm_dim78/used_blocks_check_g7.py): copied from the CM-sixfold enumeration with five changes:
# the loop runs over g = 7 only; there is no general_weil2 output for g = 7; the stored 'mult' field is used only
# where present; one record per reduced case with d >= 1 is written to blocks_check/used_blocks_check_g7_cases.json;
# and the general_weil2 comparison is printed as n/a.
#   usage:  python used_blocks_check_g7.py > blocks_check/used_blocks_check_g7_out.txt
import json, os, itertools, collections
from verify_enum import group_struct, wvec, hnf, sat_index, set_orbit

HERE = os.path.dirname(os.path.abspath(__file__))
DMAX = 6


def weil_vectors(G, g, Phi):
    gs = group_struct(G, g)
    Ps = set(Phi)
    mus = [[1 if i in T else -1 for i in range(g)] for T in set_orbit(G['gens'], Phi)]
    out = []
    for K in gs['Ks']:
        oids = sorted(K)
        b = [wvec(K[o], g) for o in oids]
        t = [2 * len(Ps & set(K[o])) - len(K[o]) for o in oids]
        sz = [len(K[o]) for o in oids]                      # |B_i| = dim A_i
        for cv in itertools.product(*[range(-(DMAX // s), DMAX // s + 1) for s in sz]):
            D = sum(abs(x) * s for x, s in zip(cv, sz))
            if D == 0 or D > DMAX or sum(x * y for x, y in zip(cv, t)) != 0:
                continue
            w = [sum(x * bb[j] for x, bb in zip(cv, b)) for j in range(g)]
            if not any(w):
                continue
            assert all(sum(a * m for a, m in zip(w, mu)) == 0 for mu in mus), "Weil character not in Lambda_U"
            copies = sum(abs(x) for x in cv)
            odd = any(x != 0 and s % 2 == 1 for x, s in zip(cv, sz))
            kind = ('i' if D <= 4 else 'ii' if copies == 1 else 'iii' if odd else 'hard')
            out.append(dict(w=w, D=D, strict=all(abs(x) <= 1 for x in cv), kind=kind,
                            B=[(o, x, s) for o, x, s in zip(oids, cv, sz) if x]))
    return out


def generates(vecs, LU, sel):
    ws = [v['w'] for v in vecs if sel(v)]
    if not ws:
        return False
    L = hnf(ws)
    ok = (L == LU)
    assert ok == (len(L) == len(LU) and sat_index(L) == 1)   # w in Lambda_U, Lambda_U saturated
    return ok


gw = {}   # no general_weil2 output for g = 7
records = []
hard_needed_all = []
for g in (7,):
    D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g)))
    groups = {G['gid']: G for G in D['groups']}
    gwr = {(r['gid'], tuple(r['Phi'])): r for r in gw.get(str(g), [])}
    cnt = collections.Counter()
    tab = collections.defaultdict(collections.Counter)
    tabN = collections.defaultdict(set)
    for c in D['cases']:
        if not (c['reduced'] and c['d'] >= 1):
            continue
        G = groups[c['gid']]
        part = tuple(sorted((len(o) // 2 for o in G['orbits']), reverse=True))
        vecs = weil_vectors(G, g, c['Phi'])
        LU = hnf(c['LU'])
        assert len(LU) == c['d']
        strict_ok = generates(vecs, LU, lambda v: v['strict'])
        all_ok = generates(vecs, LU, lambda v: True)
        cnt['reduced d>=1'] += 1
        cnt['stored passed flag agrees'] += (strict_ok == c['passed'])
        if 'mult' in c:
            cnt['stored mult.Dmin<=6 agrees'] += (all_ok == (c['mult']['Dmin'] is not None and c['mult']['Dmin'] <= 6))
        rec = dict(gid=c['gid'], Phi=c['Phi'], strict=strict_ok, all=all_ok, easy=False)
        records.append(rec)
        if not all_ok:
            cnt['failing (no generation at dim <= 6)'] += 1
            continue
        tier = 'strict' if strict_ok else 'mult<=6'
        r = gwr.get((c['gid'], tuple(c['Phi'])))
        cnt['general_weil2 Dmin<=6'] += (r is not None and r['Dmin'] is not None and r['Dmin'] <= 6)
        e1 = generates(vecs, LU, lambda v: v['kind'] == 'i')
        e12 = generates(vecs, LU, lambda v: v['kind'] in ('i', 'ii'))
        e123 = generates(vecs, LU, lambda v: v['kind'] != 'hard')
        rec['easy'] = e123
        hard_avail = any(v['kind'] == 'hard' for v in vecs)
        cnt['passing'] += 1
        cnt['passing ' + tier] += 1
        cnt['easy: (i) alone'] += e1
        cnt['easy: (i)+(ii)'] += e12
        cnt['easy: (i)+(ii)+(iii)'] += e123
        cnt['some hard B balanced (available)'] += hard_avail
        key = (part, c['d'], tier)
        tab[key]['cases'] += 1
        tab[key]['(i)'] += e1
        tab[key]['(i)+(ii)'] += e12
        tab[key]['easy'] += e123
        tab[key]['hard avail'] += hard_avail
        tabN[key].add((c['gid'], c['Nclass']))
        if not e123:
            hard = sorted(set(str(v['B']) for v in vecs if v['kind'] == 'hard'))
            hard_needed_all.append((g, c['gid'], c['Phi'], part, c['d'], tier, hard))
    print("=" * 100)
    print("g=%d: reduced (G,Phi) with d>=1: %d ; passing at dim<=6: %d (strict %d, multiplicities only %d) ; failing: %d"
          % (g, cnt['reduced d>=1'], cnt['passing'], cnt['passing strict'], cnt['passing mult<=6'],
             cnt['failing (no generation at dim <= 6)']))
    print("  generated by easy B only: %d / %d   [(i) alone: %d ; (i)+(ii): %d ; (i)+(ii)+(iii): %d] ; "
          "passing cases with some hard B available: %d ; needing a hard B: %d"
          % (cnt['easy: (i)+(ii)+(iii)'], cnt['passing'], cnt['easy: (i) alone'], cnt['easy: (i)+(ii)'],
             cnt['easy: (i)+(ii)+(iii)'], cnt['some hard B balanced (available)'],
             cnt['passing'] - cnt['easy: (i)+(ii)+(iii)']))
    print("  consistency with stored data: passed flag %d/%d, mult.Dmin<=6 %d/%d, general_weil2 Dmin<=6 %s"
          % (cnt['stored passed flag agrees'], cnt['reduced d>=1'], cnt['stored mult.Dmin<=6 agrees'],
             cnt['reduced d>=1'], ('%d/%d (passing)' % (cnt['general_weil2 Dmin<=6'], cnt['passing'])) if gw
             else 'n/a (no general_weil2 data for g = 7)'))   # g = 7 change 5
    if tab:
        print("  by (factor dims, d, tier): G-orbit cases / N-classes ; generated by (i) / (i)+(ii) / easy ; hard B available")
        for k in sorted(tab, key=lambda k: (-k[0][0], k)):
            t = tab[k]
            print("    %-14s d=%d %-8s %3d / %3d ; %3d / %3d / %3d ; %3d"
                  % (k[0], k[1], k[2], t['cases'], len(tabN[k]), t['(i)'], t['(i)+(ii)'], t['easy'], t['hard avail']))
print("=" * 100)
print("passing cases that need a 6-dimensional B with only even-dimensional factors: %d" % len(hard_needed_all))
for x in hard_needed_all:
    print("   g=%d gid=%d Phi=%s factor dims=%s d=%d tier=%s hard B's (orbit, c, dim)=%s" % x)
json.dump(records, open(os.path.join(HERE, 'blocks_check', 'used_blocks_check_g7_cases.json'), 'w'), indent=0)
