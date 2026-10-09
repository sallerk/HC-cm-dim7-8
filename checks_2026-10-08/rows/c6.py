# G = Z6 x Z2 written as Z2(a) x Z2(b) x Z3(c); iota = (1,0,0) wlog by symmetry.
import itertools, numpy as np
from fractions import Fraction
G=[(a,b,c) for a in range(2) for b in range(2) for c in range(3)]
add=lambda x,y:((x[0]+y[0])%2,(x[1]+y[1])%2,(x[2]+y[2])%3)
iota=(1,0,0)
# index-2 subgroups: kernels of nontrivial chars of Z2xZ2: (a,b)->a, b, a+b
chars={'a':lambda x:x[0],'b':lambda x:x[1],'ab':lambda x:(x[0]+x[1])%2}
Ks=[k for k,f in chars.items() if f(iota)==1]   # imaginary quadratic: iota not in kernel
print('IQ chars',Ks)
reps=[g for g in G if g[0]==0]  # one per {s, iota s}
def rank(M):
    return np.linalg.matrix_rank(np.array(M,dtype=float))
res={}
for bits in itertools.product([0,1],repeat=6):
    Phi6={ (r if s==0 else add(r,iota)) for r,s in zip(reps,bits)}
    # primitive: no nontrivial h with h+Phi6=Phi6
    if any(all(add(h,x) in Phi6 for x in Phi6) for h in G if h!=(0,0,0)): continue
    for K in Ks:
        f=chars[K]
        for c in (1,-1):
            # coordinates: reps (6) + tau (1) where tau = coset f=0
            rows=[]
            for h in G:
                v=[1 if add(h,r) in Phi6 else -1 for r in reps]  # mu_{h^-1 Phi}? sign irrelevant for rank
                # elliptic coordinate: Phi1 = {tau} if c=1; h moves tau to tau-bar if f(h)=1
                v.append(c*(1 if f(h)==0 else -1))
                rows.append(v)
            d=7-rank(rows)
            # signature wrt K: count of Phi6 elements with f=0
            p=sum(1 for x in Phi6 if f(x)==0)
            res.setdefault((K,d,(p,6-p)),0); res[(K,d,(p,6-p))]+=1
for k,v in sorted(res.items()): print(k,v)
