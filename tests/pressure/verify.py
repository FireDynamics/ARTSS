import sys
import numpy as np

try:
    p = np.loadtxt("p.dat")
except:
    print "test failed, could not load data file p.dat"
    sys.exit(1)

try:
    p_ref = np.loadtxt("p_ref.dat")
except:
    print "test failed, could not load data file p_ref.dat"
    sys.exit(1)


l0 = p.shape

if (p_ref.shape != l0):
    print "test failed, data sizes do not match"
    sys.exit(1)

d = np.sum(np.sqrt((p-p_ref)**2))

if d < 1e-6:
    print "test passed"
    sys.exit(0)
else:
    print "test failed, difference: ", d
    sys.exit(1)
