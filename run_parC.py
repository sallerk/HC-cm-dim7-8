# run_parC.py -- run stageC.py work over class-index chunks with at most Wk worker processes (default 7), then
# merge.   usage (inside the sage container):  sage -python run_parC.py [chunk_size [Wk]]
import subprocess, time, json, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
N = len(json.load(open(HERE + '/stageC_classes.json'))); cs = int(sys.argv[1]) if len(sys.argv) > 1 else 100; Wk = int(sys.argv[2]) if len(sys.argv) > 2 else 7
chunks = [(a, min(a + cs, N)) for a in range(0, N, cs)]
t0 = time.time(); running = []; todo = list(chunks)
while todo or running:
    while todo and len(running) < Wk:
        a, b = todo.pop(0); tag = 'w%05d' % a
        lf = open(HERE + '/par_logs/stageC_%s.txt' % tag, 'w')
        running.append((subprocess.Popen(['sage', '-python', HERE + '/stageC.py', 'work', str(a), str(b), tag], stdout=lf, stderr=subprocess.STDOUT), a, b, lf))
    time.sleep(2)
    for r in list(running):
        if r[0].poll() is not None:
            running.remove(r); r[3].close()
            print("chunk %d-%d rc=%d at %.1fs" % (r[1], r[2], r[0].returncode, time.time() - t0), flush=True)
subprocess.run(['python3', HERE + '/stageC.py', 'merge'])
print("total %.1fs" % (time.time() - t0), flush=True)
