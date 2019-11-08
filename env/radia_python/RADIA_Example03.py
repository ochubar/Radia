
#############################################################################
# RADIA Python Example #3: Generic Hybrid Undulator
# v 0.03
#############################################################################

from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility
import radia as rad
from uti_plot import *
import time

print('RADIA Python Example #3:')
print('This example creates and solves a simple U46 Hybrid Undulator made with rectangular magnet blocks.')
print('')

#*********************************Geometry
def Und(lp, mp, np, cp, lm, mm, nm, cm, gap, gapOffset, numPer):

    zer = [0,0,0]
    Grp = rad.ObjCnt([])

    #Principal Poles and Magnets
    y = 0.25*lp[1];

    Pole = rad.ObjFullMag([0.25*lp[0],y,-0.5*(lp[2]+gap)], [0.5*lp[0],0.5*lp[1],lp[2]], zer, np, Grp, mp, cp)
    y += 0.25*lp[1];

    mDir = -1
    for i in range(0, numPer):
        initM = [0, mDir, 0]; mDir *= -1
        y += 0.5*lm[1]
        
        Magnet = rad.ObjFullMag([0.25*lm[0],y,-0.5*(lm[2]+gap)-gapOffset], [0.5*lm[0],lm[1],lm[2]], initM, nm, Grp, mm, cm)
        y += 0.5*(lm[1] + lp[1])

        Pole = rad.ObjFullMag([0.25*lp[0],y,-0.5*(lp[2]+gap)], [0.5*lp[0],lp[1],lp[2]], zer, np, Grp, mp, cp)
        y += 0.5*lp[1]

    initM = [0, mDir, 0]
    y += 0.25*lm[1];
    
    Magnet = rad.ObjFullMag([0.25*lm[0],y,-0.5*(lm[2]+gap)-gapOffset], [0.5*lm[0],0.5*lm[1],lm[2]], initM, nm, Grp, mm, cm)

    #Mirrors
    rad.TrfZerPerp(Grp, [0,0,0], [1,0,0])
    rad.TrfZerPara(Grp, zer, [0,0,1])
    rad.TrfZerPerp(Grp, zer, [0,1,0])
    
    return Grp, Pole, Magnet

#*********************************Materials
def Materials():
    #Defines magnetic materials for the Undulators Poles and Magnets
    #Pole (~iron type Va Permendur) material data
    H = [0.8, 1.5, 2.2, 3.6, 5, 6.8, 9.8, 18, 28, 37.5, 42, 55, 71.5, 80, 85, 88, 92, 100, 120, 150, 200, 300, 400, 600, 800, 1000, 2000, 4000, 6000, 10000, 25000, 40000]
    M = [0.000998995, 0.00199812, 0.00299724, 0.00499548,0.00699372, 0.00999145, 0.0149877, 0.0299774, 0.0499648, 0.0799529, 0.0999472, 0.199931, 0.49991, 0.799899, 0.999893, 1.09989, 1.19988, 1.29987, 1.41985, 1.49981, 1.59975, 1.72962, 1.7995, 1.89925, 1.96899, 1.99874,  2.09749, 2.19497, 2.24246, 2.27743, 2.28958, 2.28973]
    convH = 4.*3.141592653589793e-07
    #ma = []
    #for i in range(len(H)): ma.append([H[i]*convH, M[i]])
    #mp = rad.MatSatIsoTab(ma)
    #A more compact way:
    mp = rad.MatSatIsoTab([[H[i]*convH, M[i]] for i in range(len(H))])

    #(Permanent) Magnet material: NdFeB with 1.2 Tesla Remanent Magnetization
    mm = rad.MatStd('NdFeB', 1.2)

    return mp, mm

def GetMagnMaterCompMvsH(MeshH, ind, cmpnH, cmpnM):
    #Extracting Magnetization vs Field Strength magnetic material data
    hMin = MeshH[0]; hMax = MeshH[1]; nh = MeshH[2]
    hStep = (hMax - hMin)/(nh - 1)
    #M = [0]*nh
    #sCmpnM = 'm' + cmpnM
    #h = hMin
    #H = [0,0,0]
    #for i in range(nh):
    #    if(cmpnH == 'x'): H[0] = h
    #    elif(cmpnH == 'y'): H[1] = h
    #    elif(cmpnH == 'z'): H[2] = h
    #    M[i] = rad.MatMvsH(ind, sCmpnM, H)
    #    h += hStep
    #A more compact way:
    if(cmpnH == 'x'): M = [rad.MatMvsH(ind, 'm'+cmpnM, [hMin+i*hStep,0,0]) for i in range(nh)]
    elif(cmpnH == 'y'): M = [rad.MatMvsH(ind, 'm'+cmpnM, [0,hMin+i*hStep,0]) for i in range(nh)]
    elif(cmpnH == 'z'): M = [rad.MatMvsH(ind, 'm'+cmpnM, [0,0,hMin+i*hStep]) for i in range(nh)]
    else: M = None
    return M

#*********************************Calculate Magnetic Field and Field Integrals
def CalcField(g, per, numPer):

    #Vertical Magnetic Field vs Longitudinal Position
    yMax = per*(numPer+1)/2; yMin = -yMax; ny = 301
    yStep = (yMax - yMin)/(ny - 1)
    xc = 0; zc = 0
    #y = yMin
    #Points = []
    #for i in range(ny):
    #    Points.append([xc,y,zc])
    #    y += yStep
    #BzVsY = rad.Fld(g, 'bz', Points)
    #A more compact way:
    BzVsY = rad.Fld(g, 'bz', [[xc,yMin+iy*yStep,zc] for iy in range(ny)])
    return BzVsY, [yMin, yMax, ny]

#*********************************Entry Point
if __name__=="__main__":

    #General Undulator Parameters
    gap = 20; numPer = 2; per = 46; gapOffset = 1

    #Pole Parameters
    lp = [45,5,25]; np = [2,2,5]; cp = [1,0,1]
    #ll = per/2 - lp[1]
    ll = 0.5*per - lp[1]

    #Magnet Parameters
    lm = [65,ll,45]; nm = [1,3,1]; cm = [0,1,1]

    #Magnetic Materials
    mp, mm = Materials()

    #Build the Structure
    und, pole, magnet = Und(lp, mp, np, cp, lm, mm, nm, cm, gap, gapOffset, numPer)
    #Show the Structure in 3D Viewer
    rad.ObjDrwOpenGL(und)

    #Solve the Magnetization Problem
    t0 = time.time()
    res = rad.Solve(und, 0.0003, 1000)
    print('Solved for Magnetization in', round(time.time() - t0, 2), 's')
    print('Relaxation Results:', res)
    print('Peak Magnetic Field:', round(rad.Fld(und, 'bz', [0,0,0]), 5), 'T') 

    #Calculate Magnetic Field
    BzVsY, MeshY = CalcField(und, per, numPer)

    #Extracting Characteristics of Magnetic Materials
    MeshH_Pole = [-0.002, 0.002, 201]
    M_Pole = GetMagnMaterCompMvsH(MeshH_Pole, pole, 'x', 'x')
    MeshH_Mag = [-1, 1, 201]
    Mpar_Mag = GetMagnMaterCompMvsH(MeshH_Mag, magnet, 'y', 'y')
    Mper_Mag = GetMagnMaterCompMvsH(MeshH_Mag, magnet, 'x', 'x')

    #Plot the Results
    uti_plot1d(M_Pole, MeshH_Pole, ['Magnetic Field Strength (mu0*H)', 'Magnetization', 'Pole Material M vs H'], ['T', 'T'])
    uti_plot1d(Mpar_Mag, MeshH_Mag, ['Magnetic Field Strength (mu0*H)', 'Magnetization Parallel to Easy Axis', 'Permanent Magnet M vs H Parallel to Easy Axis'], ['T', 'T'])
    uti_plot1d(Mper_Mag, MeshH_Mag, ['Magnetic Field Strength (mu0*H)', 'Magnetization Perpendicular to Easy Axis', 'Permanent Magnet M vs H Perpendicular to Easy Axis'], ['T', 'T'])
    uti_plot1d(BzVsY, MeshY, ['Longitudinal Position', 'Magnetic Field', 'Magnetic Field in Central Part of Undulator'], ['mm', 'T'])
    uti_plot_show() #show all graphs (and block further execution, if any)
