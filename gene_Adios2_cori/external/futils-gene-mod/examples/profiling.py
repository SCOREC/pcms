#!/usr/bin/env python

from tables import *
from pylab import *
import time, os, sys

file='demo.h5'

def mem():
    pid = os.getpid()
    arch = sys.platform
    if arch[:5] == 'linux':
        a2 = open('/proc/%d/statm'%(pid,)).read()
        mb = 4*int(a2.split()[1])/1024.0
    elif arch[:6] == 'darwin':
        a2 = os.popen('ps -p %d -o rss,vsz' % pid).readlines()
        mb = int(a2[1].split()[1])/1024.0
    else:
        mb = 0
    return mb


t0=time.time()

f=openFile(file, mode="r")
timeid=f.getNode('/profile_2d','time')
xgid=f.getNode('/profile_2d','xg')
ygid=f.getNode('/profile_2d','yg')
phiid=f.getNode('/profile_2d','phi')

print "\nelapsed time %.3f s" % (time.time()-t0,)
print "memory used %.3f MB\n" % (mem(),)

t0=time.time()
phi=phiid[0,:,:]
print "elapsed time %.3f s" % (time.time()-t0,)
print "size of phi %.3f MB" % (4*phi.size()/1024.0/1024.0,)
print "memory used %.3f MB\n" % (mem(),)

t0=time.time()
phi=phiid[-1,:,:]
print "elapsed time %.3f s" % (time.time()-t0,)
print "size of phi %.3f MB" % (4*phi.size()/1024.0/1024.0,)
print "memory used %.3f MB\n" % (mem(),)

t0=time.time()
phi=phiid[:,:,:]
print "elapsed time %.3f s" % (time.time()-t0,)
print "size of phi %.3f MB" % (4*phi.size()/1024.0/1024.0,)
print "memory used %.3f MB\n" % (mem(),)

print "Using Matlab/hdf5read"
matlab_com = "tic;phi=hdf5read('%s','/profile_2d/phi');toc;quit;"% (file,)
com="""matlab -nojvm -nosplash -r "%s" """ % matlab_com
stat=os.system(com)
