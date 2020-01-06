#!/usr/bin/env python

import sys, os.path
from tables import *

command = os.path.basename(sys.argv[0])
if len(sys.argv) != 3:
    print "Usage: %s <file> <full path of dataset>" % (command,)
    sys.exit()

file = sys.argv[1]
dset = sys.argv[2]

f=openFile(file, mode='r')
a=f.getNode(dset).read()
s=a.tostring()
sys.stdout.write(s)
f.close()
