#!/bin/bash 

rm director_twister
icc -xHost -O2  parser.c director.c -lm -lgsl -lgslcblas -o director_twister

./director_twister <data.input

gnuplot plot.gp
evince transmitance.pdf &
