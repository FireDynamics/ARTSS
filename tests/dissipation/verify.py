import sys
import numpy as np

try:
    T = np.loadtxt("T.dat")
except:
    print "test failed, could not load data file T.dat"
    sys.exit(1)

try:
    T_ref = np.loadtxt("T_ref.dat")
except:
    print "test failed, could not load data file T_ref.dat"
    sys.exit(1)


l0 = T.shape

if (T_ref.shape != l0):
    print "test failed, data sizes do not match"
    sys.exit(1)

d = np.sum(np.sqrt((T-T_ref)**2))

if d < 1e-6:
    print "test passed"
    sys.exit(0)
else:
    print "test failed, difference: ", d
    sys.exit(1)
