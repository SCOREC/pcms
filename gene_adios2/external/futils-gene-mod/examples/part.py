#!/usr/bin/env python

from tables import *
from pylab import *
import mpl3d.mplot3d as p3
import time

f=openFile("part.h5", mode="r")

step=1
partid=f.getNode('/part', '%3.3d' % step)

r=partid[:,0]
z=partid[:,1]
phi=partid[:,2]

x=r*cos(phi)
y=r*sin(phi)

fig=figure()
ax = p3.Axes3D(fig)
ax.scatter3D(x, y, z, s=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
#  Try to get an unscaled aspect ratio
vlim=ax.get_w_xlim(); ax.set_w_ylim(vlim); ax.set_w_zlim(vlim)
fig.add_axes(ax)
show()
