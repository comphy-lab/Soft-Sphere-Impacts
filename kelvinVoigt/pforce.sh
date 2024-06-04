#!/bin/bash
 We="20"

 qcc -Wall -O2 -disable-dimensions getBetaMax.c -o getBetaMax -lm
 python3 getBetaMax.py $We
 python3 plot_force.py $We
 python3 f_vcm.py $We
