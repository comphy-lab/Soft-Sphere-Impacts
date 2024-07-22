#!/bin/bash
We="100.0"
El="0.1"
Wi="0.1"

qcc -fopenmp -Wall -O2 -disable-dimensions bounceVE_WeissenbergSweep.c -o bounceVE_WeissenbergSweep -lm
export OMP_NUM_THREADS=8
./bounceVE_WeissenbergSweep $We $El $Wi
