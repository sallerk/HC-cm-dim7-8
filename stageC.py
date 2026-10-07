# stageC.py -- main computation (Sage+GAP), g = 7, ALL subgroups G <= W(B7) containing rho (non-simple CM 7-folds included).
#   sage -python stageC.py prep                 -> stageC_classes.json (generators of one G per W-class containing rho,
#                                                  as preimages of ConjugacyClassesSubgroups(W/<rho>))
#   sage -python stageC.py work i0 i1 tag       -> enum_g7_<tag>.json (run_enum.py format, + general Weil tier)
#   python3 stageC.py merge                     -> enum_g7.json
# Per class: all CM types up to G, d, Lambda_U, strict KEY TEST, primitivity/isogeny/reduced (cmenum.analyse_group);
# for reduced d>=1: IQ-multiplicity minimal dimension (cmenum.mult_weil, Dmax 14); for reduced strict failures:
# general CM-field minimal dimension (weil_tiers, LowIndexSubgroups), searched at Dmax = 6, 8, ..., 16.
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
g = 7
mode = sys.argv[1]
if mode == 'merge':
    import glob
    groups, cases = [], []
    for fn in sorted(glob.glob(HERE + '/enum_g7_w*.json')):
        D = json.load(open(fn))
        groups += D['groups']; cases += D['cases']
    groups.sort(key=lambda G: G['gid']); cases.sort(key=lambda c: (c['gid'], c['Phi']))
    json.dump(dict(g=g, groups=groups, cases=cases), open(HERE + '/enum_g7.json', 'w'))
    print("merged groups %d cases %d" % (len(groups), len(cases)))
    sys.exit(0)

from sage.all import libgap
from cmenum import *
from weil_tiers import weil_tiers

t0 = time.time()
W, rho = W_group(g)
if mode == 'prep':
    hom = libgap.NaturalHomomorphismByNormalSubgroup(W, libgap.Subgroup(W, [rho]))
    Q = libgap.ImagesSource(hom)
    cc = libgap.ConjugacyClassesSubgroups(Q)
    out = []
    for c in cc:
        H = libgap.PreImages(hom, libgap.Representative(c))
        assert bool(libgap.IsSubset(H, [rho]))
        gens = list(libgap.SmallGeneratingSet(H)) or [rho]
        out.append(dict(gens=[perm_list(x, 2 * g) for x in gens], order=int(libgap.Size(H))))
    json.dump(out, open(HERE + '/stageC_classes.json', 'w'))
    print("prep: |Q|=%d classes %d time %.1fs" % (int(libgap.Size(Q)), len(out), time.time() - t0))
    sys.exit(0)

i0, i1, tag = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
CL = json.load(open(HERE + '/stageC_classes.json'))
groups_out, cases_out = [], []
for gid in range(i0, min(i1, len(CL))):
    G = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in CL[gid]['gens']])
    assert int(libgap.Size(G)) == CL[gid]['order']
    gd, res = analyse_group(G, g, W, True)
    groups_out.append(dict(gid=gid, gens=gd['gens'], order=gd['order'], orbits=gd['orbits'],
                           orbit_TI=[o['TI'] for o in gd['oinfo']], orbit_img_order=[o['img_order'] for o in gd['oinfo']],
                           transitive=(len(gd['orbits']) == 1),
                           Ks=[dict(signs=K['signs'], blocks={str(k): v for k, v in K['blocks'].items()}) for K in gd['Ks']],
                           ntypeorbits=len(gd['torbs'])))
    for ot, ta in zip(gd['torbs'], res):
        c = dict(gid=gid, Phi=ta['Phi'], norb=ta['norb'], dimMT=ta['dimMT'], d=ta['d'], LU=ta['LU'],
                 adm=[dict(K=a['K'], J=a['J'], B=a['B'], size=a['size']) for a in ta['adm']],
                 rankW=ta['rankW'], LW=ta['LW'], index=ta['index'], satidx=ta['satidx'], passed=ta['passed'],
                 primitive=ta['primitive'], isog=ta['isog'], reduced=ta['reduced'], Nclass=ta.get('Nclass'))
        if ta['reduced'] and ta['d'] >= 1:
            c['mult'] = mult_weil(gd, ta, 14)
            if not ta['passed']:
                gen = None
                for Dm in (6, 8, 10, 12, 14, 16):
                    gen = weil_tiers(G, None, g, ta['Phi'], ta['d'], ta['LU'], Dmax=Dm, maxindex=Dm)
                    if gen['Dlat_all'] is not None:
                        break
                c['general'] = dict(Dmin=gen['Dlat_all'], Dq=gen['Dq_all'], Dmax=Dm, Dlat_IQ=gen['Dlat_IQ'],
                                    Dlat_le4=gen['Dlat_le4'], Dlat_le8=gen['Dlat_le8'], levels=gen['levels'])
        cases_out.append(c)
    if gid % 200 == 0:
        print("  [%s] group %d  cases so far %d  (%.1fs)" % (tag, gid, len(cases_out), time.time() - t0), flush=True)
json.dump(dict(g=g, groups=groups_out, cases=cases_out), open(HERE + '/enum_g7_%s.json' % tag, 'w'), default=int)
print("[%s] done groups %d..%d cases %d time %.1fs" % (tag, i0, i1, len(cases_out), time.time() - t0), flush=True)
