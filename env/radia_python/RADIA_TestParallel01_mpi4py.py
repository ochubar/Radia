
#############################################################################
# RADIA Python Test Example: Parallel (MPI based) Computation of Magnetic Field
# created by Hybrid Undulator
# v 0.01
#############################################################################

from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility
import radia as rad
from uti_plot import *
from math import *
import time

#*********************************Geometry
def HybridUndCenPart(_gap, _gap_ofst, _nper, _air, _lp, _ch_p, _np, _np_tip, _mp, _cp, _lm, _ch_m_xz, _ch_m_yz, _ch_m_yz_r, _nm, _mm, _cm, _use_ex_sym=False):

    zer = [0,0,0]
    grp = rad.ObjCnt([])

    y = _lp[1]/4
    initM = [0,-1,0]

    pole = rad.ObjFullMag([_lp[0]/4,y,-_lp[2]/2-_gap/2-_ch_p], [_lp[0]/2,_lp[1]/2,_lp[2]], zer, [_np[0],int(_np[1]/2+0.5),_np[2]], grp, _mp, _cp)

    if(_ch_p > 0.): #Pole Tip
        poleTip = rad.ObjThckPgn(_lp[0]/4, _lp[0]/2,
                                 [[y-_lp[1]/4,-_gap/2-_ch_p],
                                  [y-_lp[1]/4,-_gap/2],
                                  [y+_lp[1]/4-_ch_p,-_gap/2],
                                  [y+_lp[1]/4,-_gap/2-_ch_p]], zer)
        rad.ObjDivMag(poleTip, [_np_tip[0],int(_np_tip[1]/2+0.5),_np_tip[2]]); rad.MatApl(poleTip, mp); rad.ObjDrwAtr(poleTip, cp); rad.ObjAddToCnt(grp, [poleTip])

    y += _lp[1]/4+_air+_lm[1]/2

    for i in range(_nper):
        magnet = rad.ObjThckPgn(_lm[0]/4, _lm[0]/2,
                                [[y+_lm[1]/2-_ch_m_yz_r*_ch_m_yz,-_gap/2-_gap_ofst],
                                 [y+_lm[1]/2,-_gap/2-_gap_ofst-_ch_m_yz],
                                 [y+_lm[1]/2,-_gap/2-_gap_ofst-_lm[2]+_ch_m_yz],
                                 [y+_lm[1]/2-_ch_m_yz_r*_ch_m_yz,-_gap/2-_gap_ofst-_lm[2]],
                                 [y-_lm[1]/2+_ch_m_yz_r*_ch_m_yz,-_gap/2-_gap_ofst-_lm[2]],
                                 [y-_lm[1]/2,-_gap/2-_gap_ofst-_lm[2]+_ch_m_yz],
                                 [y-_lm[1]/2,-_gap/2-_gap_ofst-_ch_m_yz],
                                 [y-_lm[1]/2+_ch_m_yz_r*_ch_m_yz,-_gap/2-_gap_ofst]], initM)
        #Cuting Magnet Corners
        magnet = rad.ObjCutMag(magnet, [_lm[0]/2-_ch_m_xz,0,-_gap/2-_gap_ofst], [1,0,1])[0]
        magnet = rad.ObjCutMag(magnet, [_lm[0]/2-_ch_m_xz,0,-_gap/2-_gap_ofst-_lm[2]], [1,0,-1])[0]
        
        rad.ObjDivMag(magnet, _nm); rad.MatApl(magnet, _mm); rad.ObjDrwAtr(magnet, _cm); rad.ObjAddToCnt(grp, [magnet])

        initM[1] *= -1
        y += _lm[1]/2+_lp[1]/2+_air;

        if(i < _nper-1):
            pole = rad.ObjFullMag([_lp[0]/4,y,-_lp[2]/2-_gap/2-_ch_p], [_lp[0]/2,_lp[1],_lp[2]], zer, _np, grp, _mp, _cp)

            if(_ch_p > 0.): #Pole Tip
                poleTip = rad.ObjThckPgn(_lp[0]/4, _lp[0]/2,
                                         [[y-_lp[1]/2,-_gap/2-_ch_p],
                                          [y-_lp[1]/2+_ch_p,-_gap/2],
                                          [y+_lp[1]/2-_ch_p,-_gap/2],
                                          [y+_lp[1]/2,-_gap/2-_ch_p]], zer)
                rad.ObjDivMag(poleTip, _np_tip); rad.MatApl(poleTip, mp); rad.ObjDrwAtr(poleTip, cp); rad.ObjAddToCnt(grp, [poleTip])
	
            y += _lm[1]/2+_lp[1]/2+_air;

    y -= _lp[1]/4
    pole = rad.ObjFullMag([_lp[0]/4,y,-_lp[2]/2-_gap/2-_ch_p], [_lp[0]/2,_lp[1]/2,_lp[2]], zer, [_np[0],int(_np[1]/2+0.5),_np[2]], grp, _mp, _cp)
    if(_ch_p > 0.): #Pole Tip
        poleTip = rad.ObjThckPgn(_lp[0]/4, _lp[0]/2,
                                 [[y-_lp[1]/4,-_gap/2-_ch_p],
                                  [y-_lp[1]/4+_ch_p,-_gap/2],
                                  [y+_lp[1]/4,-_gap/2],
                                  [y+_lp[1]/4,-_gap/2-_ch_p]], zer)
        rad.ObjDivMag(poleTip, [_np_tip[0],int(_np_tip[1]/2+0.5),_np_tip[2]]); rad.MatApl(poleTip, mp); rad.ObjDrwAtr(poleTip, cp); rad.ObjAddToCnt(grp, [poleTip])

    #Symmetries
    if(_use_ex_sym): #Some "non-physical" mirroring (applicable for calculation of central field only)
        y += _lp[1]/4
        rad.TrfZerPerp(grp, [0,y,0], [0,1,0]) #Mirror left-right
        rad.TrfZerPerp(grp, [0,2*y,0], [0,1,0])
        
    #"Physical" symmetries (applicable also for calculation of total structure with terminations)
    rad.TrfZerPerp(grp, zer, [0,1,0]) #Mirror left-right
    #Mirror front-back
    rad.TrfZerPerp(grp, zer, [1,0,0])
    #Mirror top-bottom
    rad.TrfZerPara(grp, zer, [0,0,1])

    return grp

#*********************************All Calculations
if __name__=="__main__":

    #Initialize MPI
    #Initializing of MPI can be done in two ways.
    #The first (standard) method initializes MPI in C++. It doesn't need mp4py installed.
    #Uncomment the following line to use it:
    #rank = rad.UtiMPI('on')

    #The second method uses mpi4py. It mimics the initialization and the usage of MPI from Python script.
    #Uncomment the following 4 lines to use it:
    from mpi4py import MPI
    comMPI = MPI.COMM_WORLD
    rank = comMPI.Get_rank()
    rad.UtiMPI('in')

    if(rank <= 0):
        print('Python Test Example:')
        print('Parallel (MPI based) Computation of Magnetic Field created by a Hybrid Undulator.')
        print('To run this example in parallel mode (on a multi-core server), execute from command line:')
        print('mpiexec -n 11 python RADIA_TestParallel01.py')
        print('(in place of 11, one can specify another number of parallel processes that can be supported by server).')
        print('')

    #General Undulator Parameters
    per = 14.
    nPer = 15 #10 
    gap = 4.
    gapOffset = 0.
    air = 0.05

    poleThickPer20 = 3.
    poleHeight = 19.5
    poleWidth = 30. #40.
    lp = [poleWidth, poleThickPer20*per/20., poleHeight]
    chPole = 0.1 #If cham_pole > 0, it adds pole tip

    npx = 6; npy = 6 #5
    np = [npx,npy,[12,0.18]] #Pole Subdivision Params
    #np = [[6,1.],[5,1.],[12,0.18]]
    #np = [npx,npy,12]
    npTip = [npx,npy,1] #Pole Tip Subdivision Params
    cp = [1.,0.,1.] #Pole Color

    magWidth = 45. #58.
    lmy = per/2-lp[1]-2*air
    lm = [magWidth,lmy,27]

    nm = [3,2,[6,1./3.]] #Magnet Subdivision Params
    #nm = [[3,1.],[2,1.],[6,1./3.]]
    #nm = [3,2,6]
    cm = [0.,1.,1.] #Magnet Color

    chMagXZ = 3. #Magnet Chamfer in the XZ plane
    chMagYZ = 0.05 #Magnet Chamfer in the YZ plane
    chMagYZrat = sqrt(3.) #Magnet Chamfer Ratio: Longitudinal/Vertical

    #Pole Material 
    #B [G] vs H [G] data from NEOMAX
    BvsH_G = [[0.,0],[0.5,5000],[1,10000],[1.5,13000],[2,15000],[3,16500],[4,17400],[6,18500],[8,19250],[10,19800],
              [12,20250],[14,20600],[16,20900],[18,21120],[20,21250],[25,21450],[30,21590],[40,21850],[50,22000],
              [70,22170],[100,22300],[200,22500],[300,22650],[500,23000],[1000,23900],[2000,24900]]
    MvsH_T = [[BvsH_G[i][0]*1.e-04, (BvsH_G[i][1]-BvsH_G[i][0])*1.e-04] for i in range(len(BvsH_G))]
    mp = rad.MatSatIsoTab(MvsH_T)

    #Magnet Material
    magBr = 1.67 #Remanent Magnetization
    mm = rad.MatLin({0.05, 0.15}, magBr)
    
    grp = HybridUndCenPart(_gap=gap, _gap_ofst=gapOffset, _nper=nPer, _air=air,
                           _lp=lp, _ch_p=chPole, _np=np, _np_tip=npTip, _mp=mp, _cp=cp,
                           _lm=lm, _ch_m_xz=chMagXZ, _ch_m_yz=chMagYZ, _ch_m_yz_r=chMagYZrat, _nm=nm, _mm=mm, _cm=cm)

    #Display the Geometry
    if(rank <= 0): rad.ObjDrwOpenGL(grp)

    #Construct Interaction Matrix
    t0 = time.time()
    IM = rad.RlxPre(grp)
    if(rank <= 0): print('Interaction Matrix was set up in:', round(time.time() - t0, 2), 's')

    #Perform the Relaxation
    t0 = time.time()
    res = rad.RlxAuto(IM, 0.001, 5000)
    if(rank <= 0): 
        print('Relaxation took:', round(time.time() - t0, 2), 's')
        print('Relaxation Results:', res)

    #Synchronizing all processes
    rad.UtiMPI('barrier')
    print('   After Relaxation: rank=', rank, ' t=', time.time())

    #Calculate Magnetic Field after the Relaxation
    Bz0 = rad.Fld(grp, 'bz', [0,0,0])
    if(rank <= 0): print('Bz0=', Bz0)

    t0 = time.time()
    nyPer = 200
    yMax = 0.7*nPer*per
    ny = int(1.4*nPer*nyPer) + 1
    yStep = 2*yMax/(ny - 1)
    
    Bz = rad.Fld(grp, 'bz', [[0,-yMax+iy*yStep,0] for iy in range(ny)])

    #Synchronizing all processes
    rad.UtiMPI('barrier')
    print('   After Field Calculation: rank=', rank, ' t=', time.time())

    if(rank <= 0):
        print('Field was calculated in:', round(time.time() - t0, 2), 's')

        #Plot the calculated Magnetic Field
        uti_plot1d(Bz, [-yMax*0.001,yMax*0.001,ny], ['Longitudinal Position [m]', 'Bz [T]', 'Vertical Magnetic Field']) #, ['mm', 'T'])
        uti_plot_show()

    #End MPI
    rad.UtiMPI('off')
