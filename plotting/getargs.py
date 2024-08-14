#!/usr/bin/python3

import numpy as np
import sys, getopt

def help():
    print('plot-phosphorus.py -h     (help)')
    print('plot-phosphorus.py n      (n = number of plots in one row)')
    sys.exit()

print(sys.argv)
try:
    opts, args = getopt.getopt(sys.argv[1:],"abhc:d:")
except getopt.GetoptError:
    help()

print(opts)
print(args)

for t in opts:
    print(t)
    if t[0]=='-c':
        print(t[1])
        cols=[int(x) for x in t[1].split(',')]
        print(cols)
