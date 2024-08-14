#!/usr/bin/python3

import numpy as np
import sys, getopt

def help():
    print('plot-conc.py -h     (help)')
    print('plot-conc.py n      (n = number of node for plotting)')
    sys.exit()

prog=sys.argv[0]
argv=sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"h")
except getopt.GetoptError:
    help()

for opt, arg in opts:
    if opt == '-h':
        help()

try:
    node=argv[-1]
except:
    node=input("Please enter node number: ")

name='conc.out'
finp=open(name,'r')
name='conc-node{}.gnudata'.format(node)
fout=open(name,'w')
first=True
s=finp.readline()
s=finp.readline()
s=finp.readline()
s=finp.readline()
n=len(s.split())

while s:
    if s.startswith('#'):
        first=False
        data=s.split()
        time=data[1]
    else:
        data=s.split()
        if data[0]==node:
            fout.write(' '.join(data[2:]) + "\n")
    s=finp.readline()
finp.close()  
fout.close()  
