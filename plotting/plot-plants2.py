#!/usr/bin/python3

import numpy as np
import sys, getopt
import matplotlib.pyplot as p

def help():
    print('plot-plants2.py -h     (help)')
    print('plot-plants2.py -p 3   (number of plots in one row)')
    sys.exit()

prog=sys.argv[0]
argv=sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hp:")
except getopt.GetoptError:
    help()

ncolumns=3
for opt, arg in opts:
    if opt == '-h':
        help()
    if opt == '-p':
        ncolumns=int(arg)

name='plants2.out'
f=open(name,'r')
for n in range(0,15):
    s=f.readline()
c=s.split()
colnames=c[3:]
ncol=len(colnames)

pdata=[]
tdata=[]
s=f.readline()
while s:
    data=s.split()
    x=[]
    for d in data:
        x.append(float(d))
    pdata.append(x[4:])
    tdata.append(x[3])
    s=f.readline()

data=np.array(pdata)
time=np.array(tdata)
pdata=[]
tdata=[]


if ncolumns==1:
    nlines=ncol
else:
    nlines=int( (ncol+(ncolumns-0.1))/ncolumns )

fig, ax = p.subplots(nlines, ncolumns, sharex=True,figsize=(16,9)) #,subplot_kw={ 'yticks': []})
k=0
for i in range(0,nlines):
    for j in range(0,ncolumns):
        if ncolumns==1:
            ax1=ax[k]
        else:
            ax1=ax[i,j]
        ax1.plot(time,data[:,k],label=colnames[k])
        ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
        ax1.grid()
        ax1.legend()
        k+=1
        if k>=ncol:
            break
#ax1.set_xlabel("time")
#title.set_y(0.95)
fig.tight_layout(pad=0.0)
#fig.subplots_adjust(top=0.9)
fig.savefig('plants2.out.pdf', bbox_inches='tight', papertype='A0', dpi=600)
p.show()
