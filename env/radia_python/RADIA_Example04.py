
#############################################################################
# RADIA Python Example #4: This example illustrates the use of a polyhedron shape by means of the "radObjMltExtPgn" function.
# A uniformly magnetized polyhedron is created with a magnetization of 1 Tesla.
# The field produced by this polyhedron is computed and shown to be uniform inside the volume
# of the polyhedron and equal to 2/3 Tesla as expected from an analytical integration.
# v 0.02
#############################################################################

from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility
import radia as rad
from math import *
from uti_plot import *

print('RADIA Python Example #4:')
print('This example illustrates the use of a polyhedron shape by means of the "radObjMltExtPgn" function.')
print('A uniformly magnetized polyhedron is created with a magnetization of 1 Tesla.')
print('The field produced by this polyhedron is computed and shown to be uniform ')
print('inside the volume of the polyhedron and equal to 2/3 Tesla as expected from an analytical integration.')
print('')

#*********************************Build the Geometry
def SphericalVolume(_r, _n_phi, _nz, _M):
    dz = 2.*_r/_nz
    z = -_r + dz
    dPhi = 2.*pi/_n_phi
    allSlicePgns=[[[[0.,0.]], -_r]]
    for i in range(1, _nz):
        #print(i)
        theta = asin(z/_r)
        cosTheta = cos(theta)
        phi = dPhi
        slicePgn = [[_r*cosTheta, 0.]]
        for k in range(1, _n_phi):
            slicePgn.append([_r*cos(phi)*cosTheta, _r*sin(phi)*cosTheta])
            phi += dPhi
        allSlicePgns.append([slicePgn, z])
        z += dz
    allSlicePgns.append([[[0.,0.]], _r])
    return rad.ObjMltExtPgn(allSlicePgns, _M)

#*********************************Entry Point
if __name__=="__main__":

    #Build the Geometry 
    aSpherMag = SphericalVolume(1, 15, 15, [1,0,0])
    #Apply Color to it
    rad.ObjDrwAtr(aSpherMag, [0,0.5,0.8])

    #Display the Geometry in 3D Viewer
    rad.ObjDrwOpenGL(aSpherMag)

    #Calculate Magnetic Field
    print('Field in the Center = ', rad.Fld(aSpherMag, 'b', [0,0,0]))

    #Horizontal Field vs Longitudinal Position
    yMin = -0.99; yMax = 0.99; ny = 301
    yStep = (yMax - yMin)/(ny - 1)
    BxVsY = rad.Fld(aSpherMag, 'bx', [[0,yMin+i*yStep,0] for i in range(ny)])

    #Plot the Results
    uti_plot1d(BxVsY, [yMin, yMax, ny], ['Longitudinal Position [mm]', 'Bx [T]', 'Horizontal Magnetic Field'])
    uti_plot_show() #show all graphs
