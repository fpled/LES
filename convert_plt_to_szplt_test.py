#!/usr/bin/python

import os

for g in [32]:
    print "grid",g
    for s in [2,3,4]:
        print "std",s
        for l in range(1,7):
            print "sample",l
            for t in range(0,20):
                os.system('tec360 -convert ./Grid'+str(g)+'/EC010-'+str(s)+'/Set'+str(l)+'/PhaseInv3D_ramdom_i0000'+str(t)+'00.plt -o ./Grid'+str(g)+'/EC010-'+str(s)+'/Set'+str(l)+'/PhaseInv3D_random_i0000'+str(t)+'00.szplt')



