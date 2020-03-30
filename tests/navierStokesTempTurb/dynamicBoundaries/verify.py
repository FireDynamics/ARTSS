import sys
import numpy as np

try:
    u = np.loadtxt("u.dat")
except:
    print "velocity test failed, could not load data file u.dat"
    sys.exit(1)

try:
    v = np.loadtxt("v.dat")
except:
    print "velocity test failed, could not load data file v.dat"
    sys.exit(1)
    
try:
	w = np.loadtxt("w.dat")
except:
    print "velocity test failed, could not load data file w.dat"
    sys.exit(1)

try:
    u_ref = np.loadtxt("u_ref.dat")
except:
    print "velocity test failed, could not load data file u_ref.dat"
    sys.exit(1)

try:
    v_ref = np.loadtxt("v_ref.dat")
except:
    print "velocity test failed, could not load data file v_ref.dat"
    sys.exit(1)

try:
    w_ref = np.loadtxt("w_ref.dat")
except:
    print "velocity test failed, could not load data file w_ref.dat"
    sys.exit(1)
    
l0 = u.shape

if (v.shape != l0 or w.shape != l0 or u_ref.shape != l0 or v_ref.shape != l0 or w_ref.shape != l0):
    print "velocity test failed, data sizes do not match"
    sys.exit(1)

d = np.sum(np.sqrt((u-u_ref)**2 + (v-v_ref)**2 + (w-w_ref)**2))

if d < 1e-6:
    print "velocity test passed"
#sys.exit(0)
else:
    print "velocity test failed, difference: ", d
    sys.exit(1)

try:
    p = np.loadtxt("p.dat")
except:
    print "pressure test failed, could not load data file p.dat"
    sys.exit(1)

try:
    p_ref = np.loadtxt("p_ref.dat")
except:
    print "pressure test failed, could not load data file p_ref.dat"
    sys.exit(1)

l0 = p.shape

if (p_ref.shape != l0):
    print "pressure test failed, data sizes do not match"
    sys.exit(1)

d = np.sum(np.sqrt((p-p_ref)**2))

if d < 1e-6:
    print "pressure test passed"
#sys.exit(0)
else:
    print "pressure test failed, difference: ", d
    sys.exit(1)


try:
    T = np.loadtxt("T.dat")
except:
    print "temperature test failed, could not load data file T.dat"
    sys.exit(1)

try:
    T_ref = np.loadtxt("T_ref.dat")
except:
    print "temperature test failed, could not load data file T_ref.dat"
    sys.exit(1)

l0 = T.shape

if (T_ref.shape != l0):
    print "temperature test failed, data sizes do not match"
    sys.exit(1)

d = np.sum(np.sqrt((T-T_ref)**2))

if d < 1e-6:
    print "temperature test passed"
    sys.exit(0)
else:
    print "temperature test failed, difference: ", d
    sys.exit(1)


