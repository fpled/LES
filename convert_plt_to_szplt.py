#!/usr/bin/python

import os

for g in [16,32,64,128,256]:
    print "grid",g
    for l in range(1,41):
        print "sample",l
        for t in range(0,10):
            os.system('tec360 -convert ./Grid'+str(g)+'/Cas-'+str(l)+'/validation/InvPhase3d_0000'+str(t)+'00.plt -o ./Grid'+str(g)+'/Cas-'+str(l)+'/validation/InvPhase3d_0000'+str(t)+'00.szplt')
            print "time",t
        for t in range(10,51):
            os.system('tec360 -convert ./Grid'+str(g)+'/Cas-'+str(l)+'/validation/InvPhase3d_000'+str(t)+'00.plt -o ./Grid'+str(g)+'/Cas-'+str(l)+'/validation/InvPhase3d_000'+str(t)+'00.szplt')
            print "time",t



