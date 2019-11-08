
#############################################################################
# RADIA Python Example #1: Magnetic field created by rectangular parallelepiped with constant magnetization over volume
# v 0.02
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import radia as rad
print('RADIA Library Version:', rad.UtiVer(), '\n')

print('RADIA Python Example #1:')
print('This is the simplest example. A magnetized cube is placed at position [0,0,0].')
print('It is 1 mm in size and is magnetized according to the vector [-0.5,1,0.7] in Tesla.')
print('The three components of the field at position [0.52,0.6,0.7] are computed.')
print('Values close to [0.12737,0.028644,0.077505] are expected.')
print('')

m = rad.ObjRecMag([0,0,0], [1,1,1], [-0.5,1,0.7])
B = rad.Fld(m, "b", [0.52,0.6,0.7])

print(B)
