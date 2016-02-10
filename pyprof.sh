#!/bin/bash
#https://github.com/martinxyz/config/blob/master/scripts/pyprof
set -e # exit on error
rm -f profile.png profile.dat
python -m cProfile -o profile.dat $@
gprof2dot.py -f pstats profile.dat | dot -Tpng -o profile.png
eog profile.png
