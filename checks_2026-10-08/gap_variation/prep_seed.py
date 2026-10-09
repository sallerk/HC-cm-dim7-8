# prep_seed.py -- run "stageC.py prep" after optionally fixing Sage's random seed; print the seeds and a hash of
# GAP's random-source states before the run (State() does not consume random numbers).
#   sage -python prep_seed.py none|<seed>
import sys, runpy, hashlib
from sage.all import set_random_seed, initial_seed, libgap
arg = sys.argv[1]
if arg != 'none':
    set_random_seed(int(arg))
def h(x):
    return hashlib.md5(str(x).encode()).hexdigest()[:12]
print('arg', arg, 'sage initial_seed', initial_seed())
print('GAP GlobalMersenneTwister state', h(libgap.eval('State(GlobalMersenneTwister)')))
print('GAP GlobalRandomSource state', h(libgap.eval('State(GlobalRandomSource)')), flush=True)
sys.argv = ['stageC.py', 'prep']
runpy.run_path('stageC.py', run_name='__main__')
