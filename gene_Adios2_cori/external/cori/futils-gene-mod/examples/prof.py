#!/usr/bin/env python

from tables import *
from pylab import *

f=openFile("prof.h5", mode="r")

timeid=f.getNode('/profile_1d','time')
xid=f.getNode('/profile_1d','xgrid')
phiid=f.getNode('/profile_1d','phi')

timeid=f.getNode('/profile_2d','time')
xgid=f.getNode('/profile_2d','xg')
ygid=f.getNode('/profile_2d','yg')
phi2id=f.getNode('/profile_2d','phi')


title_ann = f.getNodeAttr('/', 'title')
time_ann  = timeid.getAttr('title')
x_ann  = xid.getAttr('title')
phi_ann  = phiid.getAttr('title')

#
#  Profiles of 1D potential
hold(True)
nsteps=size(timeid)
for i in range(10):
    label="t=%.1f" % (timeid[i],)
    plot(xid[:], phiid[i,:], label=label)

legend(loc='best')
xlabel(x_ann)
ylabel(phi_ann)
title(title_ann)

#
#  Time evolution of 1D phi at some positions x
figure()
hold(True)
n=size(xid)
for i in range(0,n/2,8):
    plot(timeid[:10], phiid[:10,i],label="X=%.2f"%xid[i])

legend(loc='best')
xlabel(time_ann)
ylabel(phi_ann)
title(title_ann)
grid(True)

#
#  Contour plots of phi(x,t)
figure()
hold(True)
#contourf(xid[:],timeid[:],phiid[:,:],10)
cset=contour(xid[:],timeid[:10],phiid[:10,:],10, linewidths=2)
clabel(cset, inline=1, fmt='%1.1f', fontsize=10)
#hot()
#colorbar()
xlabel(time_ann)
ylabel(x_ann)
title(phi_ann)

show()
