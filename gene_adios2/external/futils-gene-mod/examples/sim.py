#!/usr/bin/env python

from tables import *
from pylab import *

f=openFile("sim.h5", mode="r")
timeid=f.getNode('/scalars','time')
epotid=f.getNode('/scalars','epot')
ekinid=f.getNode('/scalars','ekin')

print "Number of steps = %d"% (timeid.nrows,)

time_ann  = timeid.getAttr('title')
epot_ann = epotid.getAttr('title')
ekin_ann = ekinid.getAttr('title')
title_ann = f.getNodeAttr('/', 'title')

time=timeid[:]
epot=epotid[:]
ekin=ekinid[:]

f.close()

hold(True)
plot(time, epot, 'b', label=epot_ann)
plot(time, ekin, 'r', label=ekin_ann)
legend(loc='best')
xlabel(time_ann);
title(title_ann);
grid(True)

show()
