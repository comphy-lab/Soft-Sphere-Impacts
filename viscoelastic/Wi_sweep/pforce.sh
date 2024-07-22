#!/bin/bash
We="100"
El="1.0"

qcc -Wall -O2 -disable-dimensions getBetaMax.c -o getBetaMax -lm
python3 getBetaMax.py $El
python3 plot_force.py $We $El
python3 f_vcm.py $We $El
