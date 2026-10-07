from cyc import *
from collections import Counter
for m in (17,34,32,40,48,60):
    U=units(m); CH=chars(m)
    seen=set(); res=Counter(); examples={}
    for a in range(1,m):
        for b in range(1,m):
            c=(-a-b)%m
            if c==0: continue
            if gcd(gcd(gcd(a,b),c),m)!=1: continue
            # canonical rep: min over t-scaling and permutations
            key=min(tuple(sorted((t*a%m,t*b%m,t*c%m))) for t in U)
            if key in seen: continue
            seen.add(key)
            H=sorted(t for t in U if (t*a%m+t*b%m+t*c%m)==m)
            assert len(H)==len(U)//2
            st=stab(m,H); r=rank(m,H,CH)
            k=(len(st), r)
            res[k]+=1
            examples.setdefault(k,[]).append(key)
    print("m=%d: #classes(a,b,c) up to scaling+perm = %d"%(m,len(seen)))
    for k in sorted(res):
        print("   |W_alpha|=%d rank=%d : %d classes; e.g. %s"%(k[0],k[1],res[k],examples[k][:6]))
