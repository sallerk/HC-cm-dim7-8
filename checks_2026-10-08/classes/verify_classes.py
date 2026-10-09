# verify_classes.py -- independent (no GAP) check of the grouping of stored G-orbits of CM types into
# N_W(G)-classes (field "Nclass") in enum_g7.json (g = 7) and trans_g8.json (g = 8).
#
# Method (pure Python 3 + numpy):
#   * W(B_g) = all 2^g g! signed permutations of {0..2g-1} (centralizer of rho(i) = i+g mod 2g), as a uint8 array.
#   * G = closure of the stored generators (vectorized BFS; permutations encoded as uint64, 4 bits per point).
#   * N_W(G) by brute force over ALL of W: w is kept iff w h w^-1 lies in G for every stored generator h.
#     (w h w^-1)(y) = w(h(w^-1(y))).  Since every element of W is determined by its images of 0..g-1, the test is done
#     point by point y = 0..g-1 on the partial codes (images of 0..y), filtering survivors at each step; the full code
#     of each surviving conjugate is then checked against the sorted codes of G once more.
#   * class invariant of a CM type Phi: min over n in N of bitmask(n(Phi)).  Two stored orbits of the same group are
#     in the same class iff the invariants agree; the induced partition is compared with the stored Nclass partition.
#   * also checked per group: |G| = stored order, G subset of N, |N| = NW_order (g = 8), N closed under products
#     (random sample), every stored Phi a CM type, the stored Phi represent every G-orbit of CM types exactly once,
#     stored orbit size norb, and that d (and reduced / primitive) is constant on each computed class.
#
# Usage:
#   python verify_classes.py controls                      positive/negative controls + timing probes
#   python verify_classes.py run 7 sel  [log]              g = 7, gids with a reduced d >= 1 case
#   python verify_classes.py run 8 sel  [log]              g = 8, groups with a primitive d >= 1 case
#   python verify_classes.py run 7 all  [log]              all 6665 groups of enum_g7.json
#   python verify_classes.py run 8 all  [log]              all 1943 groups of trans_g8.json
#   python verify_classes.py run 8 keys 16T1/0,16T5/0 [log]   (g = 7: gids, e.g. "run 7 keys 66,592")
#   python verify_classes.py summary 7|8 sel|all [log]     re-summarize a (possibly partial) JSONL log
# Each group's result is appended to the log (JSONL, flushed per group) so a long run can be checked/resumed
# (groups already present in the log are skipped).
import sys, os, json, time, itertools, collections, random
import numpy as np

DATA = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- permutations / codes
def encode(P):
    """rows of P (uint8, n <= 16 columns) -> uint64 codes, point x at bits 4x..4x+3"""
    P = np.ascontiguousarray(P, dtype=np.uint8)
    m, n = P.shape
    if n < 16:
        P = np.concatenate([P, np.zeros((m, 16 - n), np.uint8)], axis=1)
    Q = np.ascontiguousarray(P[:, 0::2] | (P[:, 1::2] << 4))
    return Q.view('<u8').ravel().copy()

def isin_sorted(a, b):
    """membership of a in sorted 1-D array b"""
    if len(b) == 0:
        return np.zeros(len(a), bool)
    pos = np.searchsorted(b, a)
    pos[pos == len(b)] = 0
    return b[pos] == a

_W = {}
def get_W(g):
    if g in _W:
        return _W[g]
    n = 2 * g
    perms = np.array(list(itertools.permutations(range(g))), dtype=np.uint8)
    signs = np.array(list(itertools.product((0, 1), repeat=g)), dtype=np.uint8)
    top = (perms[:, None, :] + np.uint8(g) * signs[None, :, :]).reshape(-1, g)
    W = np.empty((top.shape[0], n), np.uint8)
    W[:, :g] = top
    W[:, g:] = (top + g) % n
    del top
    # columns y < g of the inverses: winv[r, W[r, x]] = x
    inv_cols = [np.empty(len(W), np.uint8) for _ in range(g)]
    CH = 1 << 20
    ar = np.arange(n, dtype=np.uint8)
    for s in range(0, len(W), CH):
        blk = W[s:s + CH]
        inv = np.empty_like(blk)
        np.put_along_axis(inv, blk.astype(np.intp), np.broadcast_to(ar, blk.shape), axis=1)
        for y in range(g):
            inv_cols[y][s:s + CH] = inv[:, y]
    assert len(np.unique(encode(W))) == len(W) == 2 ** g * int(np.prod(range(1, g + 1)))
    _W[g] = (W, inv_cols)
    return _W[g]

def closure(gens, n):
    """all elements of <gens> (uint8 rows) and their sorted codes"""
    gens = [np.asarray(s, np.uint8) for s in gens]
    idp = np.arange(n, dtype=np.uint8)[None, :]
    elems = [idp]
    known = encode(idp)
    frontier = idp
    while len(frontier):
        new = np.concatenate([s[frontier] for s in gens])          # s o f
        c = encode(new)
        cu, first = np.unique(c, return_index=True)
        fresh = ~isin_sorted(cu, known)
        frontier = new[first[fresh]]
        known = np.sort(np.concatenate([known, cu[fresh]]))
        elems.append(frontier)
    E = np.concatenate(elems)
    return E, known

def normalizer(g, gens, Gcodes):
    """indices into W of N_W(G) (brute force over all of W)"""
    W, inv_cols = get_W(g)
    n = 2 * g
    Wflat = W.ravel()
    partial = []
    for y in range(g):
        mask = np.uint64((1 << (4 * (y + 1))) - 1)
        pc = np.unique(Gcodes & mask)
        if 4 * (y + 1) <= 20:                                         # small partial keys: boolean lookup table
            tab = np.zeros(1 << (4 * (y + 1)), bool); tab[pc.astype(np.int64)] = True
            partial.append(('tab', tab))
        else:
            partial.append(('srt', pc))
    idx = np.arange(len(W), dtype=np.int64)
    for h in gens:
        h = np.asarray(h, np.uint8)
        base = idx * n
        code = np.zeros(len(idx), np.uint64)
        for y in range(g):
            b = h[inv_cols[y][idx]]                                   # h(w^-1(y))
            c = Wflat[base + b]                                       # w(h(w^-1(y)))
            code |= c.astype(np.uint64) << np.uint64(4 * y)
            kind, P = partial[y]
            keep = P[code.astype(np.int64)] if kind == 'tab' else isin_sorted(code, P)
            idx = idx[keep]; base = base[keep]; code = code[keep]
        # full check of the surviving conjugates (all 2g points)
        Wi = W[idx]
        conj = np.empty_like(Wi)
        np.put_along_axis(conj, Wi.astype(np.intp), Wi[:, h], axis=1)  # conj[w(x)] = w(h(x))
        ok = isin_sorted(encode(conj), Gcodes)
        assert ok.all(), 'partial test inconsistent with full test'
    return idx

def masks(E, Phi):
    """bitmasks of e(Phi) for all rows e of E"""
    cols = E[:, Phi].astype(np.uint32)
    return (np.uint32(1) << cols).sum(axis=1, dtype=np.uint32) if len(Phi) else np.zeros(len(E), np.uint32)

def is_cm_type(Phi, g):
    return (len(Phi) == g and len(set(Phi)) == g and all(0 <= x < 2 * g for x in Phi)
            and sorted(x % g for x in Phi) == list(range(g)))

def same_partition(a, b):
    return len(set(zip(a, b))) == len(set(a)) == len(set(b))

# ---------------------------------------------------------------- data
def load(g):
    if g == 7:
        D = json.load(open(DATA + '/enum_g7.json'))
        byg = collections.defaultdict(list)
        for c in D['cases']:
            byg[c['gid']].append(c)
        groups = collections.OrderedDict()
        for G in D['groups']:
            groups[str(G['gid'])] = dict(G, cases=byg[G['gid']])
        return groups
    D = json.load(open(DATA + '/trans_g8.json'))
    groups = collections.OrderedDict()
    for G in D['groups']:
        groups['%s/%d' % (G['TI'], G['j'])] = G
    return groups

def in_scope(g, G):
    if g == 7:
        return any(c['reduced'] and c['d'] >= 1 for c in G['cases'])
    return any(c['primitive'] and c['d'] >= 1 for c in G['cases'])

# ---------------------------------------------------------------- per group
def analyse(g, key, G, rng=random.Random(1)):
    t0 = time.time()
    n = 2 * g
    E, Gcodes = closure(G['gens'], n)
    t1 = time.time()
    Nidx = normalizer(g, G['gens'], Gcodes)
    W, _ = get_W(g)
    N = W[Nidx]
    t2 = time.time()
    Ncodes = np.sort(encode(N))
    rec = dict(key=key, order=len(E), order_ok=(len(E) == G['order']), N_order=len(N))
    rec['G_in_N'] = bool(isin_sorted(Gcodes, Ncodes).all())
    rec['div_ok'] = (len(N) % len(E) == 0 and len(W) % len(N) == 0)
    if len(N) > 0:                                                    # closure of N under products (sample)
        a = N[[rng.randrange(len(N)) for _ in range(64)]]; b = N[[rng.randrange(len(N)) for _ in range(64)]]
        ab = np.take_along_axis(a, b.astype(np.intp), axis=1)           # a o b
        rec['N_closed_sample'] = bool(isin_sorted(encode(ab), Ncodes).all())
    if 'NW_order' in G:
        rec['NW_order'] = G['NW_order']; rec['NW_ok'] = (len(N) == G['NW_order'])
    cases = G['cases']
    rec['cm_ok'] = all(is_cm_type(c['Phi'], g) for c in cases)
    seen = set(); tot = 0; norb_ok = True; disjoint = True
    inv = []
    for c in cases:
        Phi = np.array(c['Phi'], np.intp)
        om = np.unique(masks(E, Phi))
        if 'norb' in c and len(om) != c['norb']:
            norb_ok = False
        s = set(om.tolist())
        if seen & s:
            disjoint = False
        seen |= s; tot += len(om)
        inv.append(int(masks(N, Phi).min()))
    rec['orbits_complete'] = (disjoint and tot == 2 ** g and len(seen) == 2 ** g)
    rec['norb_ok'] = norb_ok
    rec['inv'] = inv
    rec['Nclass'] = [c['Nclass'] for c in cases]
    rec['partition_ok'] = same_partition(inv, rec['Nclass'])
    rec['n_classes'] = len(set(inv)); rec['n_Nclass'] = len(set(rec['Nclass']))
    flag = 'reduced' if g == 7 else 'primitive'
    by = collections.defaultdict(set)
    for c, v in zip(cases, inv):
        by[v].add((c['d'], json.dumps(c[flag])))
    rec['d_flag_const'] = all(len(s) == 1 for s in by.values())
    rec['t_closure'] = round(t1 - t0, 3); rec['t_normalizer'] = round(t2 - t1, 3); rec['t_total'] = round(time.time() - t0, 3)
    return rec, N

GOOD = ('order_ok', 'G_in_N', 'div_ok', 'N_closed_sample', 'cm_ok', 'orbits_complete', 'norb_ok', 'partition_ok',
        'd_flag_const')

def bad_fields(rec):
    b = [k for k in GOOD if k in rec and not rec[k]]
    if 'NW_ok' in rec and not rec['NW_ok']:
        b.append('NW_ok')
    return b

# ---------------------------------------------------------------- runs
def run(g, scope, log, keys=None):
    groups = load(g)
    if keys is not None:
        sel = [k for k in keys]
    elif scope == 'sel':
        sel = [k for k, G in groups.items() if in_scope(g, G)]
    else:
        sel = list(groups)
    done = set()
    if os.path.exists(log):
        for line in open(log):
            if line.startswith('{'):
                done.add(json.loads(line)['key'])
    t0 = time.time(); get_W(g)
    print('W(B_%d) built: %d elements, %.1f s; %d groups to do (%d already in log)' % (g, len(get_W(g)[0]), time.time() - t0,
          len(sel), len(done & set(sel))), flush=True)
    f = open(log, 'a')
    T = time.time(); k = 0
    for key in sel:
        if key in done:
            continue
        rec, _ = analyse(g, key, groups[key])
        f.write(json.dumps(rec) + '\n'); f.flush()
        k += 1
        b = bad_fields(rec)
        if b or k % 100 == 0:
            print('%5d %-10s |G|=%-8d |N|=%-8d t=%.2fs elapsed %.0fs %s' % (k, key, rec['order'], rec['N_order'], rec['t_total'],
                  time.time() - T, ('BAD ' + ','.join(b)) if b else ''), flush=True)
    f.close()
    print('done %d groups in %.1f s' % (k, time.time() - T), flush=True)
    summary(g, scope if keys is None else 'keys', log, groups, sel)

def read_log(log):
    R = collections.OrderedDict()
    for line in open(log):
        if line.startswith('{'):
            r = json.loads(line); R[r['key']] = r
    return R

def hex2phi(s):
    return [int(ch, 16) for ch in s]

def summary(g, scope, log, groups=None, sel=None):
    if groups is None:
        groups = load(g)
    R = read_log(log)
    if sel is None:
        sel = [k for k, G in groups.items() if in_scope(g, G)] if scope == 'sel' else list(groups)
    out = []
    P = lambda s: (print(s), out.append(s))
    P('=== summary g = %d, scope %s, log %s' % (g, scope, os.path.basename(log)))
    missing = [k for k in sel if k not in R]
    P('groups in scope: %d; analysed: %d; missing: %d' % (len(sel), len(sel) - len(missing), len(missing)))
    bad = {k: bad_fields(R[k]) for k in sel if k in R and bad_fields(R[k])}
    for k in GOOD + ('NW_ok',):
        nf = sum(1 for kk in sel if kk in R and k in R[kk] and not R[kk][k])
        nt = sum(1 for kk in sel if kk in R and k in R[kk])
        P('  %-16s failed in %d of %d groups' % (k, nf, nt))
    for k, b in list(bad.items())[:50]:
        r = R[k]
        P('  BAD %s: %s  |G|=%d stored %s |N|=%d NW_order=%s inv=%s Nclass=%s' % (k, b, r['order'], groups[k]['order'], r['N_order'],
          r.get('NW_order'), r['inv'], r['Nclass']))
    allc = sum(r['n_classes'] for k, r in R.items() if k in sel)
    alln = sum(r['n_Nclass'] for k, r in R.items() if k in sel)
    allo = sum(len(r['inv']) for k, r in R.items() if k in sel)
    P('all stored orbits of the analysed groups: %d orbits, %d classes (brute force), %d distinct (group, Nclass)' % (allo, allc, alln))
    ttot = sum(r['t_total'] for k, r in R.items() if k in sel)
    big = max((R[k] for k in sel if k in R), key=lambda r: r['t_total'], default=None)
    if big:
        P('time: sum of per-group times %.1f s; slowest group %s (|G| = %d, |N| = %d) %.2f s' % (ttot, big['key'], big['order'],
          big['N_order'], big['t_total']))
    # flagged orbits
    rows = collections.OrderedDict()
    if g == 7:
        f1 = lambda c: c['reduced'] and c['d'] >= 1
        P('-- reduced orbits with d >= 1:')
    else:
        f1 = lambda c: c['primitive'] and c['d'] >= 1
        P('-- primitive orbits with d >= 1:')
    S = collections.defaultdict(set); So = collections.Counter(); Sn = collections.defaultdict(set)
    for k in sel:
        if k not in R:
            continue
        for c, v in zip(groups[k]['cases'], R[k]['inv']):
            if f1(c):
                S[c['d']].add((k, v)); So[c['d']] += 1; Sn[c['d']].add((k, c['Nclass']))
    P('   total: %d orbits, %d classes (brute force), %d distinct (group, Nclass)' % (sum(So.values()),
      sum(len(s) for s in S.values()), sum(len(s) for s in Sn.values())))
    for d in sorted(S):
        P('   d = %d: %d classes / %d orbits  (Nclass: %d)' % (d, len(S[d]), So[d], len(Sn[d])))
    if g == 7:
        summary7(groups, R, sel, P)
    else:
        summary8(groups, R, sel, P)
    name = os.path.join(OUT, 'summary_g%d_%s.txt' % (g, scope))
    open(name, 'w').write('\n'.join(out) + '\n')
    print('written', name)

def summary7(groups, R, sel, P):
    fail = lambda c: 'general' in c and not (c.get('mult') and c['mult'].get('Dmin') is not None and c['mult']['Dmin'] <= 6)
    fail2 = lambda c: 'general' in c and (c['general']['Dmin'] is None or c['general']['Dmin'] > 6)   # compact_C.py
    part = lambda G: tuple(sorted((len(o) // 2 for o in G['orbits']), reverse=True))
    F = []; n2 = 0; same_sel = True
    for k in sel:
        if k not in R:
            continue
        for c, v in zip(groups[k]['cases'], R[k]['inv']):
            if fail(c) != fail2(c):
                same_sel = False
            if fail(c):
                F.append((k, c, v))
    P('-- failing orbits (general field, not mult Dmin <= 6): %d orbits; selection identical to compact_C.py\'s: %s' % (len(F), same_sel))
    cls = set((k, v) for k, c, v in F); ncl = set((k, c['Nclass']) for k, c, v in F)
    P('   classes (brute force): %d; distinct (gid, Nclass): %d' % (len(cls), len(ncl)))
    rows = collections.OrderedDict(); rowo = collections.Counter()
    for k, c, v in sorted(F, key=lambda t: (tuple(-x for x in part(groups[t[0]])), t[1]['d'], t[1]['general']['Dmin'] or 99)):
        r = (part(groups[k]), c['d'], c['general']['Dmin'])
        rows.setdefault(r, set()).add((k, v)); rowo[r] += 1
    for r, s in rows.items():
        P('   %-10s d=%d Dmin=%-5s %3d classes (%3d orbits)' % (','.join(map(str, r[0])), r[1], r[2], len(s), rowo[r]))
    # cross-check list_C_compact.txt
    lines = [l.split() for l in open(DATA + '/list_C_compact.txt') if l[:1].isdigit()]
    LC = collections.Counter(); keyset = set(); ok = True; notes = []
    for t in lines:
        dims, gid, phihex, d, dmin = t[0], t[1], t[4], int(t[5]), t[6]
        G = groups[gid]
        Phi = hex2phi(phihex)
        if gid not in R:
            notes.append('gid %s not analysed' % gid); ok = False; continue
        # invariant of the listed Phi: it must be a stored orbit representative
        idx = [i for i, c in enumerate(G['cases']) if c['Phi'] == Phi]
        if len(idx) != 1:
            notes.append('gid %s Phi %s: %d stored matches' % (gid, phihex, len(idx))); ok = False; continue
        c = G['cases'][idx[0]]; v = R[gid]['inv'][idx[0]]
        if (gid, v) in keyset:
            notes.append('gid %s Phi %s: class listed twice' % (gid, phihex)); ok = False
        keyset.add((gid, v))
        if ','.join(map(str, part(G))) != dims or c['d'] != d:
            notes.append('gid %s Phi %s: dims/d mismatch' % (gid, phihex)); ok = False
        dm = c['general']['Dmin']
        if not (str(dm) == dmin or (dm is None and dmin.endswith('*'))):
            notes.append('gid %s Phi %s: Dmin %s vs %s' % (gid, phihex, dm, dmin)); ok = False
        LC[(dims, d, dmin)] += 1
    P('-- list_C_compact.txt: %d class lines; set of (gid, class) equal to the brute-force classes of the failing orbits: %s; '
      'all lines consistent: %s' % (len(lines), keyset == cls, ok))
    for n_ in notes[:20]:
        P('   ' + n_)
    for r, k in LC.items():
        P('   file row %-10s d=%d Dmin=%-4s %3d lines' % (r[0], r[1], r[2], k))

def summary8(groups, R, sel, P):
    nw = [(k, R[k]['N_order'], R[k].get('NW_order')) for k in sel if k in R]
    P('-- |N| vs NW_order: %d groups compared, %d agree' % (len(nw), sum(1 for k, a, b in nw if a == b)))
    lines = [l.split() for l in open(DATA + '/table_g8_compact.txt') if l[:2] == '16']
    ok = True; keyset = set(); notes = []
    for t in lines:
        key, order, orb, phihex, d = t[0], int(t[1]), int(t[2]), t[3], int(t[4])
        if key not in R:
            notes.append('%s not analysed' % key); ok = False; continue
        G = groups[key]; Phi = hex2phi(phihex)
        idx = [i for i, c in enumerate(G['cases']) if c['Phi'] == Phi]
        if len(idx) != 1:
            notes.append('%s Phi %s: %d stored matches' % (key, phihex, len(idx))); ok = False; continue
        v = R[key]['inv'][idx[0]]
        if (key, v) in keyset:
            notes.append('%s: class listed twice' % key); ok = False
        keyset.add((key, v))
        norb = sum(1 for vv in R[key]['inv'] if vv == v)
        if norb != orb or G['cases'][idx[0]]['d'] != d or G['order'] != order:
            notes.append('%s Phi %s: orb %d vs %d, d %d vs %d' % (key, phihex, norb, orb, G['cases'][idx[0]]['d'], d)); ok = False
    cls = set()
    for k in sel:
        if k in R:
            for c, v in zip(groups[k]['cases'], R[k]['inv']):
                if c['primitive'] and c['d'] >= 1:
                    cls.add((k, v))
    P('-- table_g8_compact.txt: %d class lines; set of (group, class) equal to the brute-force classes of the primitive d >= 1 '
      'orbits: %s; orb column and d consistent on all lines: %s' % (len(lines), keyset == cls, ok))
    for n_ in notes[:20]:
        P('   ' + n_)

# ---------------------------------------------------------------- controls
def controls():
    out = []
    P = lambda s: (print(s, flush=True), out.append(s))
    for g in (7, 8):
        n = 2 * g
        t = time.time(); W, _ = get_W(g)
        P('g = %d: W(B_%d) built, %d elements (2^g g! = %d), %.1f s' % (g, g, len(W), 2 ** g * int(np.prod(range(1, g + 1))), time.time() - t))
        rho = [(i + g) % n for i in range(n)]
        # positive control 1: N_W(<rho>) = W
        t = time.time(); E, C = closure([rho], n); Ni = normalizer(g, [rho], C)
        P('  positive control N_W(<rho>): |<rho>| = %d, |N| = %d, = |W|: %s (%.1f s)' % (len(E), len(Ni), len(Ni) == len(W), time.time() - t))
        # positive control 2: N_W(W) = W, W generated by a transposition pair, a g-cycle pair and a sign change
        tr = list(range(n)); tr[0], tr[1], tr[g], tr[g + 1] = 1, 0, g + 1, g
        cy = [(i + 1) % g for i in range(g)] + [g + (i + 1) % g for i in range(g)]
        sg = list(range(n)); sg[0], sg[g] = g, 0
        t = time.time(); E, C = closure([tr, cy, sg], n); Ni = normalizer(g, [tr, cy, sg], C)
        P('  positive control N_W(W): |<gens>| = %d (= |W|: %s), |N| = %d (= |W|: %s) (%.1f s)' % (len(E), len(E) == len(W), len(Ni),
          len(Ni) == len(W), time.time() - t))
        # control 3: trivial group generated by the identity -> N = W; the group <cycle pair> -> N should be smaller
        t = time.time(); E, C = closure([cy], n); Ni = normalizer(g, [cy], C)
        P('  control <g-cycle on both halves>: |G| = %d, |N| = %d (%.1f s)' % (len(E), len(Ni), time.time() - t))
    # negative controls on stored data
    for g in (7, 8):
        groups = load(g)
        key = next(k for k, G in groups.items() if in_scope(g, G) and len(set(c['Nclass'] for c in G['cases'])) >= 2
                   and max(collections.Counter(c['Nclass'] for c in G['cases']).values()) >= 2)
        rec, N = analyse(g, key, groups[key])
        lab = list(rec['Nclass'])
        P('negative controls g = %d on group %s: stored partition agrees: %s (%d orbits, %d classes)' % (g, key, rec['partition_ok'],
          len(lab), rec['n_classes']))
        cnt = collections.Counter(lab)
        i = next(i for i, x in enumerate(lab) if cnt[x] >= 2)
        split = list(lab); split[i] = max(lab) + 1                    # split one orbit off its class
        j = next(j for j, x in enumerate(lab) if x != lab[0])
        merge = list(lab); merge[j] = lab[0]                          # move one orbit into another class
        swap = list(lab); swap[0], swap[j] = lab[j], lab[0]           # same counts, different partition
        for name, L in (('split', split), ('merge', merge), ('swap (same class sizes multiset)', swap)):
            r = same_partition(rec['inv'], L)
            P('  perturbed Nclass (%s): comparison says agree = %s -> %s' % (name, r, 'FLAGGED' if not r else 'NOT FLAGGED'))
    # timing probes
    for g, key in ((7, None), (8, None)):
        groups = load(g)
        sc = [k for k, G in groups.items() if in_scope(g, G)]
        for lab, ks in (('largest |G| in scope', max(sc, key=lambda k: groups[k]['order'])),
                        ('largest |G| overall', max(groups, key=lambda k: groups[k]['order'])),
                        ('most orbits overall', max(groups, key=lambda k: len(groups[k]['cases'])))):
            rec, N = analyse(g, ks, groups[ks])
            P('timing g = %d %s: %s |G| = %d, |N| = %d, orbits %d: closure %.2f s, normalizer %.2f s, total %.2f s; bad: %s' % (
              g, lab, ks, rec['order'], rec['N_order'], len(rec['inv']), rec['t_closure'], rec['t_normalizer'], rec['t_total'],
              bad_fields(rec)))
    open(os.path.join(OUT, 'controls_out.txt'), 'w').write('\n'.join(out) + '\n')

if __name__ == '__main__':
    a = sys.argv[1:]
    if a[0] == 'controls':
        controls()
    elif a[0] == 'run':
        g = int(a[1])
        if a[2] == 'keys':
            keys = a[3].split(','); log = a[4] if len(a) > 4 else os.path.join(OUT, 'log_g%d_keys.jsonl' % g)
            run(g, 'keys', log, keys)
        else:
            log = a[3] if len(a) > 3 else os.path.join(OUT, 'log_g%d_%s.jsonl' % (g, a[2]))
            run(g, a[2], log)
    elif a[0] == 'summary':
        g = int(a[1]); log = a[3] if len(a) > 3 else os.path.join(OUT, 'log_g%d_%s.jsonl' % (g, a[2]))
        summary(g, a[2], log)
