# summarize.py -- summary of the GAP-variation experiment (runs of "stageC.py prep" and of prep.g):
# for each run, md5 of the list of classes of subgroups, number of entries that differ from the stored
# stageC_classes.json of the repository, and whether the multiset of (order, orbit lengths on 14 points)
# equals the stored one.
#   python summarize.py <stored stageC_classes.json> <run dir> ...
import sys, os, json, hashlib, collections


def orbits(gens, n=14):
    seen = set(); res = []
    for s in range(n):
        if s in seen:
            continue
        orb = {s}; st = [s]
        while st:
            x = st.pop()
            for g in gens:
                y = g[x]
                if y not in orb:
                    orb.add(y); st.append(y)
        seen |= orb; res.append(len(orb))
    return tuple(sorted(res))


def load(path):
    if path.endswith('.json'):
        return json.load(open(path))
    # gap_out.txt of prep.g: one line per class, [order, gens] with gens as lists of 0-based images
    out = []
    for line in open(path):
        line = line.strip()
        if line:
            o, g = json.loads(line)
            out.append(dict(order=o, gens=g))
    return out


stored = json.load(open(sys.argv[1]))
inv0 = collections.Counter((e['order'], orbits(e['gens'])) for e in stored)
print('stored: %d classes, md5 %s' % (len(stored), hashlib.md5(open(sys.argv[1], 'rb').read()).hexdigest()))
for d in sys.argv[2:]:
    f = os.path.join(d, 'stageC_classes.json')
    if not os.path.exists(f):
        f = os.path.join(d, 'gap_out.txt')
    if not os.path.exists(f):
        continue
    L = load(f)
    md5 = hashlib.md5(open(f, 'rb').read()).hexdigest()
    ndiff = sum(1 for a, b in zip(L, stored) if a['gens'] != b['gens']) if len(L) == len(stored) else None
    inv = collections.Counter((e['order'], orbits(e['gens'])) for e in L)
    pidx = sum(1 for a, b in zip(L, stored) if (a['order'], orbits(a['gens'])) != (b['order'], orbits(b['gens'])))
    print('%-6s %-19s classes %d  md5 %s  entries differing from stored: %s  same (order, orbits) multiset: %s'
          '  index-wise (order, orbits) mismatches: %d'
          % (os.path.basename(d.rstrip('/')), os.path.basename(f), len(L), md5, ndiff, inv == inv0, pidx))
