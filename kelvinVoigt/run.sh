#!/bin/bash
 We="20"
 Ec="0.1"

qcc -fopenmp -Wall -O2 -disable-dimensions bounceElastic.c -o bounceElastic -lm
export OMP_NUM_THREADS=8
./bounceElastic $We $Ec
