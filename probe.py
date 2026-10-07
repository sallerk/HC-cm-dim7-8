# probe.py -- feasibility probe (Sage + libgap): prints the GAP version, whether LowIndexSubgroups is bound, and
# for n = 14, 16 the number of transitive groups of degree n with a central fixed-point-free involution, how many
# of them have more than one, the largest order among them, and the elapsed time.   usage: sage -python probe.py
import time
from sage.all import libgap
print(libgap.eval('GAPInfo.Version'))
print(libgap.eval('IsBound(LowIndexSubgroups)'))
t=time.time()
for n in [14,16]:
    N=int(libgap.NrTransitiveGroups(n)); cnt=0; multi=0; orders=[]
    for k in range(1,N+1):
        T=libgap.TransitiveGroup(n,k)
        Z=libgap.Centre(T)
        invs=[z for z in libgap.Elements(Z) if int(libgap.Order(z))==2 and int(libgap.NrMovedPoints(z))==n]
        if invs:
            cnt+=1; orders.append(int(libgap.Size(T)))
            if len(invs)>1: multi+=1
    print(n,N,'with central fpf inv',cnt,'multi',multi,'max order',max(orders), time.time()-t, flush=True)
