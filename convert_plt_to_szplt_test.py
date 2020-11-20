#!/usr/bin/python

import os

for g in [32]:
    print "grid",g
    for l in range(1,6):
        print "sample",l
        os.system('tec360 -convert ./Grid'+str(g)+'/Set'+str(l)+'/PhaseInv3D_random_i00000000.plt -o ./Grid'+str(g)+'/Set-'+str(l)+'/PhaseInv3D_random_i00000000.szplt')
        os.system('tec360 -convert ./Grid'+str(g)+'/Set-'+str(l)+'/PhaseInv3D_random_i00001300.plt -o ./Grid'+str(g)+'/Set-'+str(l)+'/PhaseInv3D_random_i00001300.szplt')



