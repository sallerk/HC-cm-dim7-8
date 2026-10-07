# verify_mutation.py -- sensitivity test of verify_trans.py: apply single random corruptions to trans_g8.json records
# and count how many are detected.  usage: python verify_mutation.py [ntrials]
import json, random, copy, sys, collections
from verify_trans import verify_group

D = json.load(open('trans_g8.json'))
random.seed(7)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
degidx = [(gi, ci) for gi, G in enumerate(D['groups']) for ci, c in enumerate(G['cases']) if c['primitive'] and c['d'] >= 1]
allidx = [(gi, ci) for gi, G in enumerate(D['groups']) for ci, c in enumerate(G['cases'])]
kinds = ['d', 'LU', 'primitive', 'order', 'Phi_swap', 'IQ_t', 'CMsub', 'n0', 'weil_Dlat', 'weil_Dq', 'gens', 'drop_case']
det = collections.Counter(); tot = collections.Counter(); missed = []
for trial in range(N):
    kind = kinds[trial % len(kinds)]
    pool = degidx if kind in ('IQ_t', 'CMsub', 'n0', 'weil_Dlat', 'weil_Dq') else allidx
    gi, ci = random.choice(pool)
    G = copy.deepcopy(D['groups'][gi])
    c = G['cases'][ci]
    if kind == 'd':
        c['d'] += random.choice([-1, 1])
    elif kind == 'LU':
        if not c['LU']:
            c['LU'] = [[1] + [0] * 7]
        else:
            r = random.randrange(len(c['LU'])); j = random.randrange(8)
            c['LU'][r][j] += random.choice([-1, 1])
    elif kind == 'primitive':
        c['primitive'] = not c['primitive']
    elif kind == 'order':
        G['order'] *= 2
    elif kind == 'Phi_swap':
        i = random.randrange(8)
        P = set(c['Phi'])
        x = next(x for x in P if x % 8 == i)
        P.remove(x); P.add((x + 8) % 16)
        c['Phi'] = sorted(P)
    elif kind == 'IQ_t':
        if not c['IQ']:
            c['IQ'] = [dict(B=list(range(8)), t=0)]
        else:
            c['IQ'][0]['t'] = -c['IQ'][0]['t'] + (0 if c['IQ'][0]['t'] else 2)
    elif kind == 'CMsub':
        if c['CMsub']:
            s = random.choice(c['CMsub'])
            s['sig'][0][1] = s['b'] - s['sig'][0][1] if 2 * s['sig'][0][1] != s['b'] else s['sig'][0][1] + 1
        else:
            c['CMsub'] = [dict(deg=2, b=8, sig=[[list(range(8)), 4]], balanced=True)]
    elif kind == 'n0':
        c['n0'] += 1
    elif kind == 'weil_Dlat':
        key = random.choice(['Dlat_IQ', 'Dlat_le4', 'Dlat_le8', 'Dlat_all'])
        c['weil'][key] = {None: 8, 8: 16, 16: None}[c['weil'][key]]
    elif kind == 'weil_Dq':
        key = random.choice(['Dq_IQ', 'Dq_le4', 'Dq_le8', 'Dq_all'])
        c['weil'][key] = {None: 16, 8: None, 16: 8}[c['weil'][key]]
    elif kind == 'gens':
        p = G['gens'][0][:]
        a, b = random.sample(range(16), 2)
        p[a], p[b] = p[b], p[a]
        G['gens'][0] = p
    elif kind == 'drop_case':
        if len(G['cases']) > 1:
            del G['cases'][ci]
        else:
            G['cases'][0]['Phi'] = [0, 1, 2, 3, 4, 5, 6, 7] if G['cases'][0]['Phi'] != [0, 1, 2, 3, 4, 5, 6, 7] else [0, 1, 2, 3, 4, 5, 6, 15]
    sample = set(tuple(x['Phi']) for x in G['cases'])
    try:
        errs, res, nZ, sS = verify_group(G, 8, True, sample)
        found = bool(errs) or any(ce for _, _, ce, _ in res)
    except Exception as e:
        found = True
    tot[kind] += 1
    det[kind] += found
    if not found:
        missed.append((kind, G['TI'], G['j'], c.get('Phi')))
print("mutation test: %d trials, %d detected" % (sum(tot.values()), sum(det.values())))
for k in kinds:
    print("  %-10s %d/%d" % (k, det[k], tot[k]))
for m in missed[:20]:
    print("  MISSED", m)
