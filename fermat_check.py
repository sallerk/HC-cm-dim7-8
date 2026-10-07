# fermat_check.py -- pure-Python recheck of Fermat-factor types H_{a,b,c} on Q(zeta_m), phi(m)=16 (rank = dim MT,
# stabilizer |W| = #{u : uH = H}), listed per class of triples up to scaling and permutation (the convention of
# section B of cm16_control_data.txt, an independent brute-force table), to reconcile with
# controls_cyclo_out.txt (which lists one line per G-orbit of types).
from math import gcd
import itertools, collections
from verify_enum import rank
for m in [17, 32, 34, 40, 48, 60]:
    U = [a for a in range(1, m) if gcd(a, m) == 1]
    seen = set(); res = collections.Counter(); ex = collections.defaultdict(list)
    for a in range(1, m):
        for b in range(a, m):
            c = (-a - b) % m
            if c == 0 or c < b or gcd(gcd(gcd(a, b), c), m) != 1:
                continue
            key = min(tuple(sorted(((t * a) % m, (t * b) % m, (t * c) % m))) for t in U)
            if key in seen:
                continue
            seen.add(key)
            H = frozenset(t for t in U if (t * a) % m + (t * b) % m + (t * c) % m == m)
            orb = set(frozenset((h * x) % m for x in H) for h in U)
            mu = [[1 if u in T else -1 for u in U if u < m / 2] for T in orb]
            rk = rank(mu) + 1
            stab = sum(1 for u in U if frozenset((u * x) % m for x in H) == H)
            res[(stab, rk)] += 1
            ex[(stab, rk)].append(key)
    print(m, dict(sorted(res.items())), {k: v[:4] for k, v in sorted(ex.items()) if k[0] == 1 and k[1] < 9})
