# timing_c.py -- Stage C runtime estimate: ConjugacyClassesSubgroups of W(B_g)/<rho> (GAP), g given.
import time, sys
from sage.all import libgap
g = int(sys.argv[1]); n = 2*g
rho = libgap.PermList([((i+g) % n)+1 for i in range(n)])
W = libgap.Centralizer(libgap.SymmetricGroup(n), rho)
hom = libgap.NaturalHomomorphismByNormalSubgroup(W, libgap.Subgroup(W, [rho]))
Q = libgap.ImagesSource(hom)
print("g", g, "|W|", libgap.Size(W), "|Q|", libgap.Size(Q), "Q degree", libgap.NrMovedPoints(Q), flush=True)
t = time.time()
cc = libgap.ConjugacyClassesSubgroups(Q)
print("Q classes", len(cc), "time %.1fs" % (time.time()-t), flush=True)
