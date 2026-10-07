import itertools, cmath, math
from math import gcd
def units(m): return [a for a in range(1,m) if gcd(a,m)==1]
def chars(m):
    # brute-force characters of (Z/m)^* via generators: build from discrete logs using a generating set
    U=units(m); n=len(U)
    # find a basis by brute force: represent group as product of cyclic groups
    # simple approach: compute all homomorphisms to roots of unity of order exp
    exp=1
    def order(a):
        k=1;x=a
        while x!=1: x=x*a%m;k+=1
        return k
    for a in U: exp=exp*order(a)//gcd(exp,order(a))
    # pick generators greedily
    gens=[];sub={1}
    for a in sorted(U,key=lambda a:-order(a)):
        if a in sub: continue
        gens.append(a)
        new=set()
        for s in sub:
            x=1
            for _ in range(order(a)):
                new.add(s*x%m); x=x*a%m
        sub=new
        if len(sub)==n: break
    # express each element as exponent vector (not unique if not direct basis) -> use homomorphism enumeration
    res=[]
    ords=[order(g) for g in gens]
    for exps in itertools.product(*[range(exp) for _ in gens]):
        # candidate chi(g_i)=zeta_exp^{e_i}; check well-defined
        val={}
        ok=True
        for vec in itertools.product(*[range(o) for o in ords]):
            x=1; e=0
            for g,k,ei in zip(gens,vec,exps):
                x=x*pow(g,k,m)%m; e+=k*ei
            e%=exp
            if x in val and val[x]!=e: ok=False;break
            val[x]=e
        if ok and len(val)==n: res.append(val)
    # dedupe
    uniq=[]
    seen=set()
    for v in res:
        key=tuple(v[a] for a in U)
        if key not in seen: seen.add(key); uniq.append(v)
    assert len(uniq)==n, (len(uniq),n)
    return uniq,exp
def rank(m,S,CH=None):
    U=units(m)
    if CH is None: CH=chars(m)
    ch,exp=CH
    r=0
    for v in ch:
        s=sum(cmath.exp(2j*math.pi*v[a]/exp) for a in S)
        if abs(s)>1e-9: r+=1
    return r
def stab(m,S):
    S=set(S); return [t for t in units(m) if {t*s%m for s in S}==S]
if __name__=="__main__":
    m=32; CH=chars(m)
    for S in ([1,7,13,15,21,23,27,29],[1,7,9,11,13,15,27,29],[1,7,13,21,23,27,29,17]):
        print(m,S,"rank",rank(m,S,CH),"stab",stab(m,S))
    S17=list(range(1,9)); CH17=chars(17)
    print(17,S17,"rank",rank(17,S17,CH17),"stab",stab(17,S17))
