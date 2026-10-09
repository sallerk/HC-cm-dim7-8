# Part 1: dump TransitiveGroup(n,k) generators (0-based), orders, names. Only GAP use of Parts 1-2.
import json, time
from sage.all import libgap
t0 = time.time()
out = {}
for n, K in ((5, 5), (6, 16), (8, 50)):
    assert int(libgap.NrTransitiveGroups(n)) == K
    L = []
    for k in range(1, K + 1):
        G = libgap.TransitiveGroup(n, k)
        gens = [[int(x) - 1 for x in libgap.ListPerm(g, n)] for g in libgap.GeneratorsOfGroup(G)]
        L.append({"k": k, "order": int(libgap.Size(G)), "name": str(libgap.Name(G)),
                  "structure": str(libgap.StructureDescription(G)), "gens": gens})
    out[str(n)] = L
json.dump(out, open("transgroups_5_6_8.json", "w"))
print("done", time.time() - t0)
