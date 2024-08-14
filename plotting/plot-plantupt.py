#!/usr/bin/python3

import numpy as np
import sys
import matplotlib.pyplot as p

name='plantupt.out'
f=open(name,'r')
while True:
    s=f.readline()
    if s[0]!='#':
        break

pdata=[]
while s:
    data=s.split()
    x=[]
    for d in data:
        x.append(float(d))
    pdata.append(x)
    s=f.readline()

data=np.array(pdata)
pdata=[]

time=data[:,3]
fig, ax = p.subplots(5, 4, sharex=True,figsize=(16,9)) #,subplot_kw={ 'yticks': []})
#title=fig.suptitle(name)
# row 0
ax1=ax[0,0]
ax1.plot(time,data[:,4],label='rndemlv')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[0,1]
ax1.plot(time,data[:,5],label='rndemst')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[0,2]
ax1.plot(time,data[:,6],label='rndemrt')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[0,3]
ax1.plot(time,data[:,7],label='rndemso')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
# row 1
ax1=ax[1,0]
ax1.plot(time,data[:,10],label='tunc')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[1,1]
ax1.plot(time,data[:,11],label='tund')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[1,2]
ax1.plot(time,data[:,9],label='totdem')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[1,3]
ax1.plot(time,data[:,12],label='totup')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
# row 2
ax1=ax[2,0]
ax1.plot(time,data[:,13],label='anlv')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[2,1]
ax1.plot(time,data[:,14],label='anst')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[2,2]
ax1.plot(time,data[:,15],label='anrt')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[2,3]
ax1.plot(time,data[:,16],label='anso')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
# row 3
ax1=ax[3,0]
ax1.plot(time,data[:,19],label='anclv')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[3,1]
ax1.plot(time,data[:,20],label='ancst')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[3,2]
ax1.plot(time,data[:,21],label='ancrt')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[3,3]
ax1.plot(time,data[:,22],label='ancso')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
# row 4
ax1=ax[4,0]
ax1.plot(time,data[:, 8],label='rndemcrn')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[4,1]
ax1.plot(time,data[:,17],label='ancrn')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[4,2]
ax1.plot(time,data[:,23],label='anccrn')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()
ax1=ax[4,3]
ax1.plot(time,data[:,18],label='NitReduct')
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useOffset=False)
ax1.grid()
ax1.legend()

#title.set_y(0.95)
fig.tight_layout()
#fig.subplots_adjust(top=0.9)
fig.savefig(name+'.pdf', bbox_inches='tight', papertype='A0', dpi=600)
#fig.savefig(name+'.png', bbox_inches='tight', dpi=300)
p.show()
