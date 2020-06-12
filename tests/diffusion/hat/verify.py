import sys
import numpy as np

try:
    u = np.loadtxt("u.dat")
except:
    print("test failed, could not load data file u.dat")
    sys.exit(1)

try:
    v = np.loadtxt("v.dat")
except:
    print("test failed, could not load data file v.dat")
    sys.exit(1)
    
try:
	w = np.loadtxt("w.dat")
except:
    print("test failed, could not load data file w.dat")
    sys.exit(1)

try:
    u_ref = np.loadtxt("u_ref.dat")
except:
    print("test failed, could not load data file u_ref.dat")
    sys.exit(1)

try:
    v_ref = np.loadtxt("v_ref.dat")
except:
    print("test failed, could not load data file v_ref.dat")
    sys.exit(1)

try:
    w_ref = np.loadtxt("w_ref.dat")
except:
    print("test failed, could not load data file w_ref.dat")
    sys.exit(1)
    
l0 = u.shape

if (v.shape != l0 or w.shape != l0 or u_ref.shape != l0 or v_ref.shape != l0 or w_ref.shape != l0):
    print("test failed, data sizes do not match")
    sys.exit(1)

d = np.sum(np.sqrt((u-u_ref)**2 + (v-v_ref)**2 + (w-w_ref)**2))

if d < 1e-6:
    print("test passed")
    sys.exit(0)
else:
    print("test failed, difference:", d)
    sys.exit(1)
