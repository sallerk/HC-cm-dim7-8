# run_case.py -- run one case of verify_generalC.py (same functions, same output line), for parallel runs.
#   copy to the repository root (it imports verify_generalC.py and reads enum_g7.json), then
#   python run_case.py <index>      index into the 340 cases with a 'general' record, 0-based, in file order
#   (README.txt section 10: indices 337, 338, 339 = cases 338-340, gid 6402)
import sys, json, time
from verify_generalC import chars_by_D, dmin

i = int(sys.argv[1])
D = json.load(open('enum_g7.json'))
groups = {G['gid']: G for G in D['groups']}
cases = [c for c in D['cases'] if 'general' in c]
c = cases[i]
G = groups[c['gid']]
t0 = time.time()
Dm = c['general']['Dmin'] if c['general']['Dmin'] is not None else c['general']['Dmax']
byD = chars_by_D(G, 7, c['Phi'], Dm)
Dlat, Dq = dmin(byD, c['LU'], c['d'])
ok = (Dlat == c['general']['Dmin'])
line = ' '.join(str(a) for a in (7, c['gid'], c['Phi'], '|G|=%d verifier Dmin=%s producer Dmin=%s %s' % (
    G['order'], Dlat, c['general']['Dmin'], 'AGREE' if ok else 'DISAGREE')))
with open('log_case_%d.txt' % i, 'w') as f:
    f.write(line + '\n')
print(line)
print('case %d: %.1f s' % (i, time.time() - t0))
