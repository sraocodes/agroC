#!/usr/bin/python3

import numpy as np
import sys, getopt
import matplotlib.pyplot as p

def help():
    print('plot-mass.py -h     (help)')
    print('plot-mass.py -p 3   (number of plots in one row)')
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

name='mass.out'
f=open(name,'r')
s=f.readline()
colnames=s.split()
colnames[0]=colnames[0].split('#')[1] # remove first char
ncol=len(colnames)

pdata=[]
s=f.readline()
s=f.readline()
while s:
    data=s.split()
    x=[]
    for d in data:
        x.append(float(d))
    pdata.append(x)
    s=f.readline()
f.close()
data=np.array(pdata)
time=data[:,0]
pdata=[]

if ncolumns==1:
    nlines=ncol-1
else:
    nlines=int( (ncol+(ncolumns-1.1))/ncolumns )
fig, ax = p.subplots(nlines, ncolumns, sharex=True,figsize=(16,9)) #,subplot_kw={ 'yticks': []})
#title=fig.suptitle(name)
# row 0
k=1
for i in range(0,nlines):
    for j in range(0,ncolumns):
        if ncolumns==1:
            ax1=ax[k-1]
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
fig.savefig('mass.out.pdf', bbox_inches='tight', papertype='A0', dpi=600)
#fig.savefig('mass.out.png', bbox_inches='tight', dpi=300)
p.show()
