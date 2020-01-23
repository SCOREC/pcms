#!/usr/bin/env python

from tables import *
from pylab import *
import os, sys

f=openFile("prof.h5", mode="r")

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
          interpolation='bilinear')
#          interpolation='bilinear', cmap=cm.gray)
xlabel(xg_ann)
ylabel(yg_ann)
files = []
msteps=min(nsteps,500)
for i in range(msteps):
    label="t=%.1f/%.1f" % (timeid[i],timeid[msteps-1])
    title(label)
    im.set_data(phiid[i,:,:])
    fname = '_tmp%03d.png'%i
    print 'Saving frame', fname
    savefig(fname)
    files.append(fname)

print 'Making movie animation.mpg - this make take a while'
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o demo.mpg")
#os.system("convert _tmp*.png animation.mng")

# cleanup
for fname in files: os.remove(fname)
