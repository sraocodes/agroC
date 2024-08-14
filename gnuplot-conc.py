#!/usr/bin/python3

import numpy as np
import sys, getopt

def help():
    print('plot-conc.py -h     (help)')
    print('plot-conc.py n      (n = number of column for plotting)')
    print('plot-conc.py -l n   (-l for log10(data))')
    sys.exit()

prog=sys.argv[0]
argv=sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hl")
except getopt.GetoptError:
    help()

log=False
for opt, arg in opts:
    if opt == '-h':
        help()
    elif opt=="-l":
        log=True

name='conc.out'
f=open(name,'r')
pdata=[]
pdepth=[]
first=True
s=f.readline()
s=f.readline()
s=f.readline()
s=f.readline()
n=len(s.split())
try:
    col=int(argv[-1])
except:
    c=input("Please enter column number from 1 to {}: ".format(n))
    try:
        col=int(c) 
    except:
        print('Missing column number, run: plot-conc.py -h')
        sys.exit(2)
col-=1
while s:
    if s.startswith('#'):
        first=False
    else:
        data=s.split()
        if first:
            pdepth.append(float(data[1]))
        pdata.append(float(data[col]))
    s=f.readline()
f.close()  
depth=np.array(pdepth)
ny=depth.size
nl=len(pdata)
#print(nl)
#print(ny)
nx=nl//ny
data=np.array(np.transpose(np.reshape(pdata,(nx,ny))))
#sh=data.shape
#print(sh)
pdata=[]
pdepth=[]
#print(data)
#print(depth)
time=np.linspace(0,nx-1,nx)
if log:
    z=np.log10(data)
else:
    z=data

name='conc-col{}.gnudata'.format(col+1)
f=open(name,'w')
for iy in range(ny):
    for ix in range(nx):
        f.write('{} {} {}\n'.format(time[ix],depth[iy],z[iy,ix]))
    f.write("\n")
f.close()
