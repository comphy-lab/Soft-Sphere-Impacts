#!/bin/bash
 We="20"
 Ec="1.0"

 qcc -Wall -O2 -disable-dimensions getBetaMax.c -o getBetaMax -lm
 python3 getBetaMax.py $We
 python3 plot_force.py $We $Ec
 python3 f_vcm.py $We $Ec
 python3 getRfoot.py
