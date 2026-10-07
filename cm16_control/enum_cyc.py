from cyc import *
from collections import Counter
for m in (17,32,40,48,60):
    U=units(m); CH=chars(m)
    pairs=sorted({tuple(sorted((a,m-a))) for a in U})
    allS=[]
    for choice in itertools.product((0,1),repeat=len(pairs)):
        S=tuple(sorted(p[c] for p,c in zip(pairs,choice)))
        allS.append(S)
    seen=set(); orbits=[]
    for S in allS:
        if S in seen: continue
        orb={tuple(sorted(t*s%m for s in S)) for t in U}
        seen|=orb; orbits.append((S,len(orb)))
    cnt=Counter(); ex={}
    for S,L in orbits:
        st=stab(m,S); r=rank(m,S,CH)
        key=("prim" if len(st)==1 else "imprim|stab=%d"%len(st), r)
        cnt[key]+=1
        if len(st)==1 and r<9 and key not in ex: ex[key]=S
    print("m=",m,"#CMtypes",len(allS),"#orbits",len(orbits))
    for k in sorted(cnt): print("   ",k,"orbits:",cnt[k])
    for k,S in ex.items(): print("   example",k,S)
