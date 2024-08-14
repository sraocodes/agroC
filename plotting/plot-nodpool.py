#!/usr/bin/python3

import numpy as np
import sys, getopt
import matplotlib.pyplot as p
import matplotlib.colors as c

def help():
    print('plot-nodpool.py -h     (help)')
    print('plot-nodpool.py n      (n = number of column for plotting, starting from 1)')
    print('plot-nodpool.py -l n   (-l for log10(data))')
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

try:
    col=int(argv[-1])-1
except:
    print('Missing column number, run: plot-nodpool.py -h')
    sys.exit(2)

name='nod_pool.out'
f=open(name,'r')
pdata=[]
pdepth=[]
first=True
s=f.readline()
s=f.readline()
s=f.readline()
s=f.readline()
while s:
    if s.startswith('#'):
        first=False
    else:
        data=s.split()
        if first:
            pdepth.append(float(data[1]))
        pdata.append(float(data[col]))
    s=f.readline()

depth=np.array(pdepth)
ny=depth.size
nl=len(pdata)
#print(nl)
#print(ny)
nx=nl//ny
data=np.array(np.transpose(np.reshape(pdata,(nx,ny))))
sh=data.shape
#print(sh)
pdata=[]
pdepth=[]
#print(data)
#print(depth)
time=np.linspace(1,nx,nx)
if log:
    z=np.log10(data)
else:
    z=data
p.contourf(time, depth, z, 200, cmap='coolwarm');
#p.imshow(data, extent=[data.min(), 2*data.max(), data.min(), data.max()], origin='lower', cmap='coolwarm');
p.colorbar();
p.show()
