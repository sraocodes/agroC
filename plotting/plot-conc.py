#!/usr/bin/python3

import numpy as np
import sys, getopt
import matplotlib.pyplot as p
import matplotlib.colors as c

def help():
    print('plot-conc.py n       (n=node to plot)')
    print('plot-conc.py -h      (help)')
    print('plot-conc.py -l n    (-l for log10(data))')
    sys.exit()

prog=sys.argv[0]
argv=sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hl")
except getopt.GetoptError:
    help()

log=False
ncolumns=3
for opt, arg in opts:
    if opt == '-h':
        help()
    elif opt == "-l":
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
    print(
        """
         1 = node
         2 = depth
         3 = Urea in solution
         4 = Urea sorbed
         5 = NH4 in solution
         6 = NH4 sorbed
         7 = NO3 in solution
         8 = NO3 sorbed
         9 = first additional solute in solution
        10 = first additional solute sorbed
        ...
        """
    )
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
#p.imshow(z, extent=[z.min(), 2*z.max(), z.min(), z.max()], origin='lower', cmap='coolwarm');
p.colorbar();
(name1,name2)=name.split('.')
pdfout="{}-col{}.pdf".format(name1,col+1)
p.savefig(pdfout, bbox_inches='tight', papertype='A0', dpi=600)
p.show()
