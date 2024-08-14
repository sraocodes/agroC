#!/usr/bin/python3

import numpy as np
import sys, getopt
import matplotlib.pyplot as p

def help():
    print('plot-nitrogen.py -h     (help)')
    print('plot-nitrogen.py -p 3   (number of plots in one row)')
    sys.exit()

prog=sys.argv[0]
argv=sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hp:")
except getopt.GetoptError:
    help()

ncolumns=5
for opt, arg in opts:
    if opt == '-h':
        help()
    if opt == '-p':
        ncolumns=int(arg)

name='nitrogen.out'
f=open(name,'r')
for n in range(0,14):
    s=f.readline()
c=s.split()
colnames1=c[3:]

pdata=[]
s=f.readline()
while s:
    data=s.split()
    x=[]
    for d in data:
        x.append(float(d))
    pdata.append(x[4:])
    s=f.readline()
f.close()
data1=np.array(pdata)

pdata=[]
name='plantupt.out'
f=open(name,'r')
for n in range(0,26):
    s=f.readline()
c=s.split()
colnames2=c[3:]

pdata=[]
tdata=[]
s=f.readline()
while s:
    data=s.split()
    x=[]
    for d in data:
        x.append(float(d))
    tdata.append(x[3])
    pdata.append(x[4:])
    s=f.readline()
f.close()
data2=np.array(pdata)
time=np.array(tdata)
pdata=[]
tdata=[]

data=np.concatenate((data1,data2),axis=1)
colnames=np.concatenate((colnames1,colnames2))
ncol=len(colnames)

if ncolumns==1:
    nlines=ncol
else:
    nlines=int( (ncol+(ncolumns-0.1))/ncolumns )

fig, ax = p.subplots(nlines, ncolumns, sharex=True,figsize=(16,9)) #,subplot_kw={ 'yticks': []})
k=0
for i in range(0,nlines):
    for j in range(0,ncolumns):
        if ncolumns==1 or nlines==1:
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

#title.set_y(0.95)
fig.tight_layout()
#fig.subplots_adjust(top=0.9)
fig.savefig('nitrogen.out.pdf', bbox_inches='tight', papertype='A0', dpi=600)
#fig.savefig('nitrogen.out.png', bbox_inches='tight', dpi=300)
p.show()
