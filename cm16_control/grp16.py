import itertools, numpy as np
from collections import Counter
# group C8 x| C2 = <t,s | t^8, s^2, s t s = t^u>; elements (i,j) = t^i s^j
def make(u):
    def mul(a,b):
        (i,j),(k,l)=a,b
        # t^i s^j t^k s^l = t^{i + k*u^j} s^{j+l}
        return ((i + k*pow(u,j))%8, (j+l)%2)
    G=[(i,j) for i in range(8) for j in range(2)]
    return G,mul
def inv(G,mul,a):
    for b in G:
        if mul(a,b)==(0,0): return b
def analyze(name,u):
    G,mul=make(u); e=(0,0); eps=(4,0)
    assert all(mul(eps,g)==mul(g,eps) for g in G)
    idx={g:n for n,g in enumerate(G)}
    # pairs {g, eps g}
    pairs=[]; seen=set()
    for g in G:
        if g in seen: continue
        h=mul(eps,g); pairs.append((g,h)); seen|={g,h}
    cnt=Counter()
    for choice in itertools.product((0,1),repeat=8):
        S=frozenset(p[c] for p,c in zip(pairs,choice))
        st=[g for g in G if frozenset(mul(g,s) for s in S)==S]
        rf=[g for g in G if frozenset(mul(s,g) for s in S)==S]
        A=np.array([[1 if mul(inv(G,mul,g),h) in S else 0 for h in G] for g in G])
        r=np.linalg.matrix_rank(A)
        cnt[(len(st),len(rf),r)]+=1
    print(name)
    for k in sorted(cnt): print("   |stab|=%d |reflexsubgp|=%d rank=%d : %d CM types"%(k+(cnt[k],)))
analyze("D8 (dihedral order 16), u=7",7)
analyze("SD16 (semidihedral), u=3",3)
analyze("M16 (modular), u=5",5)
analyze("C8xC2, u=1",1)
