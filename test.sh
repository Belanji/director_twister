#!/bin/bash 

rm director_twister
gcc parser.c director.c -lm -lgsl -lgslcblas -o director_twister

./director_twister <data.input
