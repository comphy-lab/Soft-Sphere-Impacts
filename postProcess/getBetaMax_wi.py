# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import sys

from snapshot_tables import extract_snapshot_table

extract_snapshot_table("getBetaMax", f"{float(sys.argv[1]):3.2f}_El_pforce.dat")
