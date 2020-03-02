
#############################################################################
# RADIA Python Example #6: Simple Quadrupole Magnet
# v 0.05
#############################################################################

from __future__ import absolute_import, division, print_function
import radia as rad

from uti_plot import *
from time import *
from math import *
from array import *

print('RADIA Python Example #6:')
print('')
print('This example considers a simple quadrupole.')
print('It is an iron dominated structure.')
print('The quadrupole presents hyperbolic pole faces with a chamber.')
print('')
print('The field computation with Radia in the case of iron dominated geometries ')
print('presents a few specific difficulties and is usually less accurate than ')
print('in the case of coil or permanent magnet dominated structure. Nevertheless, ')
print('special methods have been developed to reach a reasonable precision within ')
print('a reasonable CPU time and memory.')
print('')
print('This example is in some way similar to the simple dipole magnet of Example#5.nb')
print('and all the remarks made there are also valid here.')
print('The user should at least first read the introduction of Example#5.nb ')
print('before playing with this one.')
print('')
print('Among the most important things to remember if one wants to reach a good ')
print('precision in a reasonable time are: ')
print('- Always segment the corners of the iron circuits as parallel or perpendicular ')
print('to the flux line as possible. For right angle corners, this can be done with ')
print('the circular or ellipsoidal mode of segmentation (see below). ')
print('This is extremely important. ')
print('In this example, we make use of the circular segmentation twice (the other ')
print('two corners are segmented similarly by symmetry) it is made in the section ')
print('entitled "Defining the Model" when calling ')
print('  "rad.ObjDivMag(g3, [nr3,np3,nx], \'cyl\', cy)" and ')
print('  "rad.ObjDivMag(g5, [nr5,np5,nx], \'cyl\', cy)". ')
print('- Use a narrower segmentation on the iron faces close to the region ')
print('of interest. ')
print('- Start with low segmentation numbers and increase them gradually to check that ')
print('the field is stable. Beware that the memory and cpu time tends to grow like ')
print('the square number of elements. ')
print('- Use symmetries as much as possible. It saves both memory and CPU time. ')
print('')

#*********************************Geometry
def Geom():

    #Pole faces
    rap = 0.5
    ct = [0,0,0]
    z0 = gap/2
    y0 = width/2
    amax = hyp*asinh(y0/z0)
    dz = z0*(cosh(amax)-1)
    aStep = amax/np
    na = int(amax*(1+2/np)/aStep) + 1
    qq = [[(z0*sinh(ia*aStep/hyp)), (z0*cosh(ia*aStep))] for ia in range(na)]
    hh = qq[np][1] + height*rap - dz
    qq[np+1] = [qq[np][0],hh]
    qq[np+2] = [0,hh]
    g1 = rad.ObjThckPgn(thick/4, thick/2, qq)   
    rad.ObjDivMag(g1, n1)

    #Vertical segment on top of pole faces
    g2 = rad.ObjRecMag([thick/4,width/4,gap/2+height*(1/2+rap/2)], [thick/2,width/2,height*(1-rap)])
    rad.ObjDivMag(g2, n2)

    #Corner
    gg = rad.ObjCnt([g1, g2])
    gp = rad.ObjCutMag(gg, [thick/2-chamfer-gap/2,0,0],[1,0,-1])[0]
    g3 = rad.ObjRecMag([thick/4,width/4,gap/2+height+depth/2], [thick/2,width/2,depth])
    cy = [[[0,width/2,gap/2+height],[1,0,0]], [0,0,gap/2+height], 2*depth/width]
    rad.ObjDivMag(g3, [nr3,np3,nx], 'cyl', cy)

    #Horizontal segment between the corners
    tan_n = tan(2*pi/2/Nn)
    length = tan_n*(height+gap/2)-width/2
    g4 = rad.ObjRecMag([thick/4,width/2+length/2,gap/2+height+depth/2], [thick/2,length,depth])
    rad.ObjDivMag(g4, n4)

    #The other corner
    posy = width/2 + length
    posz = posy/tan_n
    g5 = rad.ObjThckPgn(thick/4, thick/2, [[posy,posz],[posy,posz+depth],[posy+depth*tan_n,posz+depth]])
    cy = [[[0,posy,posz],[1,0,0]], [0,posy,posz+depth], 1]
    rad.ObjDivMag(g5, [nr5,np5,nx], 'cyl', cy)

    #Generation of the coil
    Rmax = Rmin - width/2 + gap/2 + offset - 2
    coil1 = rad.ObjRaceTrk([0,0,gap/2+height/2+offset/2], [Rmin,Rmax], [thick,width-2*Rmin], height-offset, 3, CurDens)
    rad.ObjDrwAtr(coil1, coilcolor)
    hh = (height - offset)/2
    coil2 = rad.ObjRaceTrk([0,0,gap/2+height-hh/2], [Rmax,Rmax+hh*0.8], [thick,width-2*Rmin], hh, 3, CurDens)
    rad.ObjDrwAtr(coil2, coilcolor)

    #Make container, set the colors and define symmetries
    g = rad.ObjCnt([gp,g3,g4,g5])
    rad.ObjDrwAtr(g, ironcolor)
    gd=rad.ObjCnt([g])
    
    rad.TrfZerPerp(gd, ct, [1,0,0])
    rad.TrfZerPerp(gd, ct, [0,1,0])
    
    t=rad.ObjCnt([gd, coil1, coil2])
    rad.TrfZerPara(t, ct, [0, cos(pi/Nn), sin(pi/Nn)])

    rad.TrfMlt(t, rad.TrfRot(ct, [1,0,0], 4*pi/Nn), int(round(Nn/2)))
    rad.MatApl(g, ironmat)
    rad.TrfOrnt(t, rad.TrfRot([0,0,0],[1,0,0],pi/Nn))
    
    return t

#*********************************For Harmonic Analysis
def harm(obj, y, z, ro, np):
    arHarm = [complex(0,0)]*np
    d_tet = 2*pi/np
    tet = 0
    for i in range(np):
        cosTet = cos(tet); sinTet = sin(tet)
        re = rad.FldInt(obj, 'inf', 'ibz', [-1,y+ro*cosTet,z+ro*sinTet], [1,y+ro*cosTet,z+ro*sinTet])
        im = rad.FldInt(obj, 'inf', 'iby', [-1,y+ro*cosTet,z+ro*sinTet], [1,y+ro*cosTet,z+ro*sinTet])
        arHarm[i] = complex(re, im)
        tet += d_tet
    return [np, ro, arHarm]

def multipole(w, q):
    np = w[0]; ro = w[1]; arHarm = w[2]
    s = 0
    for p in range(np):
        arg = -2*pi/np*q*p
        s += arHarm[p]*complex(cos(arg), sin(arg))
    return s/np/(ro**q)
    
#*********************************Entry Point
if __name__=="__main__":

    #Creating Geometry and Solving
    #General Geometry Params
    np = 4
    hyp = 1 #must be >0.2 !
    gap = 40 #magnetic gap
    width = 30 #pole width
    height = 50 #height of iron yoke
    thick = 60 #length of the quad
    chamfer = 8
    
    Nn = 4
    ironcolor = [0,0.5,1]
    coilcolor = [1,0,0]
    depth = 1.2*width/2

    #Coils
    CurDens = -3 #current density [A/mm^2]
    Rmin = 2
    offset = 10

    #Segmentation Params
    nx = 2
    ny = 2
    n1 = [nx,ny,3]
    n2 = [nx,ny,3]
    np3 = 2
    nr3 = ny
    n4 = [nx,3,ny]
    np5 = ceil(np3/2); print()
    nr5 = ny

    t0 = time()
    rad.UtiDelAll()
    ironmat = rad.MatSatIsoFrm([2000,2],[0.1,2],[0.1,2])
    g = Geom()
    size = rad.ObjDegFre(g)

    t1 = time()
    Nmax = 10000
    res = rad.Solve(g, 0.00001, Nmax)
    t2 = time()
    
    Bz = rad.Fld(g,'Bz',[0,1,0])*1000
    Iz = rad.FldInt(g,'inf','ibz',[-1,1,0],[1,1,0])
    Iz1 = rad.FldInt(g,'inf','ibz',[-1,10,0],[1,10,0])/10

    print('M_Max H_Max N_Iter = ', round(res[1],4), 'T', round(res[2],4), 'T', round(res[3]))
    if(res[3] == Nmax): print('Unstable or Incomplete Relaxation')
    print('Built & Solved in ', round(t1-t0,3),'&', round(t2-t1, 3),' seconds')
    print('Interaction Matrix : ', size, 'X', size, 'or', round(size*size*4/1000000, 3),'MB')
    print('Gradient = ', round(Bz,4),'T/m')
    print('Int. Quad. @ 1 mm = ', round(Iz,5), 'T')
    print('Delta Int. Quad. @ 10 mm = ', round(100*(Iz1/Iz-1),2), '%')
    
    #Display the Geometry
    rad.ObjDrwOpenGL(g)

    #Magnetic Field Plots
    z = 0; x1 = 0; x2 = 30; ymax = 40; np = 20
    Bz1 = rad.FldLst(g, 'bz', [x1,0,z], [x1,ymax,z], np, 'arg', 0)
    Bz2 = rad.FldLst(g, 'bz', [x2,0,z], [x2,ymax,z], np, 'arg',0)
    uti_plot1d_m([Bz1,Bz2],
                 labels=['Y', 'Vertical Magnetic Field', 'Vertical Magnetic Field vs. Vertical Position'], units=['mm', 'T'],
                 styles=['-b.', '--r.'], legend=['X = {} mm'.format(x1), 'X = {} mm'.format(x2)])

    #Magnetic Field Integral Plots
    z = 0; ymin = 0.001; ymax = 10; npy = 20
    yStep = (ymax - ymin)/(npy - 1)
    IBz0 = rad.FldInt(g, 'inf', 'ibz', [-1,1,z], [1,1,z])

    IBzVsY = [(rad.FldInt(g, 'inf', 'ibz', [-1,ymin+iy*yStep,z], [1,ymin+iy*yStep,z])/((ymin+iy*yStep)*IBz0) - 1)*100 for iy in range(npy)]
    uti_plot1d(IBzVsY, [ymin,ymax,npy],
               ['Y', 'dIz', 'Rel. Variation of Vertical Field Integral along X at Z = ' + repr(z) + ' mm'], ['mm', '%'])

    #Harmonic Analysis of the Field Integrals
    nharm = 10; radius = 2; y = 0; z = 0;
    w = harm(g, y, z, radius, nharm);

    mm = [multipole(w, i) for i in range(nharm)];
    round_mm = [complex(round(mm[i].real, 9), round(mm[i].imag, 9)) for i in range(nharm)];
    print(round_mm)
    
    uti_plot_show()
