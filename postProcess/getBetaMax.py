# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import sys

from snapshot_tables import extract_snapshot_table, weber_label

extract_snapshot_table("getBetaMax", f"{weber_label(sys.argv[1])}_We_pforce.dat")
