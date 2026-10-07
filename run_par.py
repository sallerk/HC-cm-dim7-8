# run_par.py -- run trans_enum.py g over k-chunks with at most 7 worker processes; then merge the chunk JSONs.
# usage (inside the sage container):  python3 run_par.py g chunk_size [workers]
import sys, subprocess, json, time, os
HERE = os.path.dirname(os.path.abspath(__file__))
g = int(sys.argv[1]); cs = int(sys.argv[2]); W = int(sys.argv[3]) if len(sys.argv) > 3 else 7
NT = {14: 63, 16: 1954}[2 * g]
chunks = [(a, min(a + cs - 1, NT)) for a in range(1, NT + 1, cs)]
t0 = time.time()
running = []
todo = list(chunks)
while todo or running:
    while todo and len(running) < W:
        a, b = todo.pop(0)
        tag = 'c%04d' % a
        lf = open(HERE + '/par_logs/stdout_%s.txt' % tag, 'w')
        p = subprocess.Popen(['sage', '-python', HERE + '/trans_enum.py', str(g), str(a), str(b), tag], stdout=lf, stderr=subprocess.STDOUT)
        running.append((p, a, b, tag, lf))
    time.sleep(2)
    for r in list(running):
        if r[0].poll() is not None:
            running.remove(r)
            r[4].close()
            print("chunk %d-%d done rc=%d at %.1fs" % (r[1], r[2], r[0].returncode, time.time() - t0), flush=True)
groups = []
for a, b in chunks:
    D = json.load(open(HERE + '/trans_g%d_c%04d.json' % (g, a)))
    groups += D['groups']
json.dump(dict(g=g, groups=groups), open(HERE + '/trans_g%d.json' % g, 'w'), default=int)
print("merged %d (G,rho) pairs ; total %.1fs" % (len(groups), time.time() - t0), flush=True)
