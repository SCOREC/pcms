#!/usr/bin/env python

from tables import *
from pylab import *
import time

ion()

f=openFile("pprof.h5", mode="r")

timeid=f.getNode('/profile_2d','time')
xgid=f.getNode('/profile_2d','xg')
ygid=f.getNode('/profile_2d','yg')
phiid=f.getNode('/profile_2d','phi')


title_ann = f.getNodeAttr('/', 'title')
time_ann  = timeid.getAttr('title')
xg_ann  = xgid.getAttr('title')
yg_ann  = ygid.getAttr('title')
phi_ann  = phiid.getAttr('title')

nsteps=size(timeid)
extent=(xgid[0],xgid[-1],ygid[0],ygid[-1])
im=imshow(phiid[0,:,:], extent=extent,
          interpolation='nearest')
#          interpolation='nearest', cmap=cm.hot)
xlabel(xg_ann)
ylabel(yg_ann)
#colorbar(orientation='horizontal')
tstart = time.time()                 # for profiling
msteps=min(nsteps,200)
for i in range(msteps):
    label="t=%.1f" % (timeid[i],)
    title(label)
    im.set_data(phiid[i,:,:])
    draw()

print 'FPS:' , nsteps/(time.time()-tstart)
