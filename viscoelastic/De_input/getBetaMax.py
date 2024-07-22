# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import sys

nGFS = 100

We = int(sys.argv[1])
#ci = int(sys.argv[1])
#name = "%d_getForce.dat" % ci
name = "%d_We_pforce.dat" % We

if os.path.exists(name):
    print("File %s found! New data will be appended to the file" % name)
for ti in range(nGFS):
    t = 0.010 * ti
    place = "intermediate/snapshot-%5.4f" % t
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        exe = "./getBetaMax %s %s" % (place, name)
        os.system(exe)
    print(("Done %d of %d" % (ti, nGFS)))
