# Post-check of the Lambda_U control for d = 1 primitive types, allowing for several imaginary quadratic
# subfields (rho-swapped 2-block systems). Reuses the function definitions of pxc_control.py (not its main loop).
import collections
import os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pxc_control.py")).read().split("# ================= main loop")[0]
exec(src)

rep = collections.Counter(); bad = []; info_k = {}
for P in TG:
    k = P["k"]
    Pel = closure(P["gens"], 8); Pinv = np.argsort(Pel, axis=1)
    gens16 = [list(w_perm16(np.array([g]), np.zeros((1, 8), dtype=np.int64))[0]) for g in P["gens"]]
    systems = block_systems(gens16 + [[(i + 8) % 16 for i in range(16)]])
    iq = [S for S in systems if len(S) == 2 and all(((b >> i) & 1) != ((b >> ((i + 8) % 16)) & 1)
                                                    for b in S for i in range(16))]
    wB = [np.array([1 if (S[0] >> i) & 1 else -1 for i in range(8)]) for S in iq]   # w_{K'} for K' <-> S
    nd1 = 0
    for t in range(256):
        s = 1 - 2 * E256[t]; V = s[Pinv]
        rk, null = rank_null((V.T @ V).tolist())
        if rk != 7 or not primitive_mask(int(PHIMASK[t]), systems):
            continue
        nd1 += 1
        v = np.array(prim_int(null[0]))
        mus = np.concatenate([V, -V])
        assert (mus @ v == 0).all()
        which = [a for a, w in enumerate(wB) if (v == w).all() or (v == -w).all()]
        phi = int(PHIMASK[t])
        bal = [a for a, S in enumerate(iq) if bin(phi & S[0]).count("1") == 4]
        ok = len(which) == 1 and which == bal
        kidx = [a for a, S in enumerate(iq) if 0x00FF in S]
        rep[(ok, len(iq), which == kidx)] += 1
        if not ok:
            bad.append((k, t, v.tolist(), which, bal))
    info_k[k] = (len(iq), nd1)
print("d=1 primitive CM types (all 256 per k, not up to G): (ok, #IQ systems, LU = Z w_K for the constructed K) ->", dict(rep))
print("bad:", bad)
print("per k: (#IQ subfields, #d=1 primitive types):", {k: v for k, v in info_k.items() if v[1]})
