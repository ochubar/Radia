
#############################################################################
# RADIA Python Example #2: This example consists in the creation of a set of racetrack
# and circular coils, plotting the coil geometry and the magnetic field produced.
# v 0.04
#############################################################################

from __future__ import print_function #Py 2.*/3.* compatibility
import radia as rad
from uti_plot import *

print('RADIA Python Example #2:')
print('This example consists in the creation of a set of racetrack and circular coils,')
print('plotting the coil geometry and the magnetic field produced. This geometry corresponds ')
print('to a 4T superconducting wiggler that used to be in operation at the ESRF.')
print('')

#*********************************Build the Geometry
def BuildGeometry():

    #Current Densities in A/mm^2
    j1 = 128; j2 = 256

    #Coil Presentation Parameters
    n1 = 3; n2 = 6; c2 = [1,0,0]; c1 = [0,1,1]; thcn = 0.001

    #Create 5 Coils
    Rt1 = rad.ObjRaceTrk([0.,0.,38.], [9.5,24.5], [120.,0.], 36, n1, j1)
    rad.ObjDrwAtr(Rt1, c1, thcn)
    Rt3 = rad.ObjRaceTrk([0.,0.,76.], [10.,25.], [90.,0.], 24, n1, j1)
    rad.ObjDrwAtr(Rt3, c1, thcn)
    Rt2 = rad.ObjRaceTrk([0.,0.,38.], [24.5,55.5], [120.,0.], 36, n1, j2)
    rad.ObjDrwAtr(Rt2, c2, thcn)
    Rt4 = rad.ObjRaceTrk([0.,0.,76.], [25.,55.], [90.,0.], 24, n1, j2)
    rad.ObjDrwAtr(Rt4, c2, thcn)
    Rt5 = rad.ObjRaceTrk([0.,0.,60.], [150.,166.3], [0.,0.], 39, n2, -j2)
    rad.ObjDrwAtr(Rt5, c2, thcn)

    Grp = rad.ObjCnt([Rt1, Rt2, Rt3, Rt4, Rt5])

    #Define Mirror Coils
    rad.TrfZerPara(Grp, [0,0,0], [0,0,1])

    return Grp

#*********************************Calculate Magnetic Field and Field Integrals
def CalcField(g):

    #Vertical Magnetic Field vs Longitudinal Position
    yMin = 0.; yMax = 300.; ny = 301
    yStep = (yMax - yMin)/(ny - 1)
    xc = 0.; zc = 0.
    #y = yMin
    #Points = []
    #for i in range(ny):
    #    Points.append([xc,y,zc])
    #    y += yStep
    #BzVsY = rad.Fld(g, 'bz', Points)
    #More compact method to do the above: 
    BzVsY = rad.Fld(g, 'bz', [[xc,yMin+iy*yStep,zc] for iy in range(ny)])

    #Vertical Magnetic Field Integral (along Longitudinal Position) vs Horizontal Position
    xMin = 0.; xMax = 400.; nx = 201
    xStep = (xMax - xMin)/(nx - 1)
    zc = 0.
    #x = xMin
    #IBzVsX = []
    #for i in range(ny):
    #    IBzVsX.append(rad.FldInt(g, 'inf', 'ibz', [x,-300.,zc], [x,300.,zc]))
    #    x += xStep
    #More compact method to do the above: 
    IBzVsX = [rad.FldInt(g, 'inf', 'ibz', [xMin+ix*xStep,-300.,zc], [xMin+ix*xStep,300.,zc]) for ix in range(nx)]
    
    return BzVsY, [yMin, yMax, ny], IBzVsX, [xMin, xMax, nx]

#*********************************Entry Point
if __name__=="__main__":

    #Build the Geometry
    g = BuildGeometry()
    print('SCW Geometry Index:', g)

    #Display the Geometry in 3D Viewer
    rad.ObjDrwOpenGL(g)

    #Calculate Magnetic Field
    BzVsY, MeshY, IBzVsX, MeshX = CalcField(g)

    print('Field in Center:', BzVsY[0], 'T')
    print('Field Integral in Center:', IBzVsX[0], 'T.mm')

    #Plot the Results
    uti_plot1d(BzVsY, MeshY, ['Longitudinal Position [mm]', 'Bz [T]', 'Vertical Magnetic Field'])
    uti_plot1d(IBzVsX, MeshX, ['Horizontal Position [mm]', 'Integral of Bz [T.mm]', 'Vertical Magnetic Field Integral'])
    uti_plot_show() #show all graphs (and block further execution, if any)


