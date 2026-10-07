# probe_norm.py -- timing probe (Sage + libgap): time of Normalizer(S_n, T) for each transitive group T of
# degree n with a central fixed-point-free involution; prints progress and the 15 slowest cases
# (seconds, k, |T|, |N(T)|).   usage: sage -python probe_norm.py n   (stored output for n = 16: probe_norm16.txt)
import time, sys
from sage.all import libgap
n=int(sys.argv[1])
S=libgap.SymmetricGroup(n)
t0=time.time(); worst=[]
N=int(libgap.NrTransitiveGroups(n))
for k in range(1,N+1):
    T=libgap.TransitiveGroup(n,k)
    Z=libgap.Centre(T)
    invs=[z for z in libgap.Elements(Z) if int(libgap.Order(z))==2 and int(libgap.NrMovedPoints(z))==n]
    if not invs: continue
    t=time.time()
    NT=libgap.Normalizer(S,T)
    dt=time.time()-t
    worst.append((dt,k,int(libgap.Size(T)),int(libgap.Size(NT))))
    if k%100==0: print(k, 'elapsed %.1f'%(time.time()-t0), flush=True)
worst.sort(reverse=True)
print('total %.1f'%(time.time()-t0)); print('worst', worst[:15])
