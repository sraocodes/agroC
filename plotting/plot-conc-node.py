#!/usr/bin/python3

import numpy as np
import sys, getopt
import matplotlib.pyplot as p
import matplotlib.colors as c

def help():
    print('plot-conc-node.py -h          (help)')
    print('plot-conc-node.py -n 10,20    (nodes, starting from 1)')
    print('plot-conc-node.py -p 3        (number of plots in one row)')
    sys.exit()

prog=sys.argv[0]
argv=sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hn:p:",["help","nodes=","ncolumns="])
except getopt.GetoptError:
    help()

nodes=""
ncolumns=3
for opt, arg in opts:
    if opt == '-h':
        help()
    if opt == '-n':
        nodes=arg
    if opt == '-p':
        ncolumns=int(arg)

if nodes=="":
    help()

name='conc.out'
f=open(name,'r')
with open(name, 'r') as f:
    lines=f.readlines()
names=lines[0].split()
names[0]=names[0][1:]
ncol=len(names)
max_nodes=lines[3].split()[0]

nodes=nodes.split(',')
nodenum=[]
for n in nodes:
    nodenum.append(int(n))
nodenum.sort(reverse=True)
nodenum=np.array(nodenum)
ns=nodenum.size
nodename=[]
for n in nodenum:
    nodename.append('n=' + str(n))

conc=[]
for i in range(0,ns):
    conc.append([])

i=0
t=0
for s in lines[2:]:
    if s[0]=='#':
        i=0
        t+=1
    else:
        data=s.split()
        if i<ns and int(data[0])==nodenum[i]:
            fdata=[]
            for f in data:
                fdata.append(float(f))
            conc[i].append(fdata)
            i+=1

data=np.zeros((ns,t,ncol))
i=0
for d in conc:
    data[i,:,:]=d
    i+=1

if ncolumns==1:
    nlines=ncol-2
else:
    nlines=int( (ncol-2+ncolumns-1.1)/ncolumns )
    
fig, ax = p.subplots(nlines, ncolumns, sharex=True,figsize=(16,9)) #,subplot_kw={ 'yticks': []})

k=2
for i in range(0,nlines):
    for j in range(0,ncolumns):
        if ncolumns==1:
            ax1=ax[k-2]
        else:
            ax1=ax[i,j]
        for c in range(ns):
            ax1.plot(data[c,:,k],'-',label=nodename[c]) # '.-' for symbols
        ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
        ax1.grid()
        ax1.legend()
        ax1.set_title(names[k])
        k+=1
        if k>=ncol:
            break

#title.set_y(0.95)
fig.tight_layout()
#fig.subplots_adjust(top=0.9)
fig.savefig('conc-node.out.pdf', bbox_inches='tight', papertype='A0', dpi=600)
#fig.savefig('mass.out.png', bbox_inches='tight', dpi=300)
p.show()
