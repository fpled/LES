#!/usr/bin/python

import os

for g in [32]:
    print "grid",g
    for s in ["10-1","10-2","510-1"]:
        print "std=",s
        for l in range(1,7):
            print "sample",l
            for t in range(0,10):
                os.system('tec360 -convert ./Grid'+str(g)+'/EC0'+s+'/Set'+str(l)+'/PhaseInv3D_ramdom_i00000'+str(t)+'00.plt -o ./Grid'+str(g)+'/EC0'+s+'/Set'+str(l)+'/PhaseInv3D_ramdom_i00000'+str(t)+'00.szplt')
            for t in range(10,21):
                os.system('tec360 -convert ./Grid'+str(g)+'/EC0'+s+'/Set'+str(l)+'/PhaseInv3D_ramdom_i0000'+str(t)+'00.plt -o ./Grid'+str(g)+'/EC0'+s+'/Set'+str(l)+'/PhaseInv3D_ramdom_i0000'+str(t)+'00.szplt')



