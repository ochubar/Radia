
#############################################################################
# RADIA Python Test Script
# v 0.01
#############################################################################

from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility
import radia as rad
import sys
import time

print('RADIA Version:', rad.UtiVer())

print('RADIA Python Test Script 01')
print('')

prb = (10.0**-15)
prcrd = (10.0**-14)
fld_arg1 = 'PrcB->' + repr(prb)
fld_arg2 = 'PrcCoord->' + repr(prcrd)
fld_arg3 = fld_arg1+','+fld_arg2 
rad.FldCmpPrc(fld_arg1)
rad.FldCmpPrc(fld_arg2)
rad.FldCmpPrc(fld_arg3)


##help(rad.ObjRecMag)

mag00 = rad.ObjRecMag([5,0,0], [5,8,10], [-0.5,1,0.7])

#rad.ObjDrwOpenGL(mag00)

#data = rad.ObjDrwVTK(mag00, 'Axes->False')
##print(data)

#print(data['polygons']['vertices'])

#print(6*4*3)

#print('Number of Polygon Vertex Coords:', len(data['polygons']['vertices']))



mag01 = rad.ObjThckPgn(20, 10., [[-10,-10], [-12,5], [5,0], [7,-15]], 'x', [0,0,1])

vertices = [[-10,0,0], [0,0,0], [0,10,0], [0,0,10]]
faces = [[1,2,3], [1,4,2], [2,4,3], [1,3,4]]
mag02 = rad.ObjPolyhdr(vertices, faces, [0,0,1])

mag03 = rad.ObjArcPgnMag([0,5], 'z', [[2,0],[2,10],[10,10],[10,5]], [0,0.5], 5, 'nosym', [0,0,1])

mag04 = rad.ObjMltExtPgn([[[[-10,-10],[-15,-5],[-5,5],[5,5],[10,-15]], -15], [[[-5,-5],[-7.5,-2.5],[-2.5,2.5],[2.5,2.5],[5,-7.5]], -7]], [0,0,1])

mag05 = rad.ObjMltExtRtg([[[0,0,12],[5,10]], [[5,10,20],[15,5]]], [0,0,1])

mag06 = rad.ObjMltExtTri(25, 8, [[0,-15],[-15,0],[0,15],[15,0]], [[5,1],[5,2],[5,3],[5,1]], 'z', [0,0,1], 'ki->Numb,TriAngMin->20,TriAreaMax->10')

mag07 = rad.ObjCylMag([0,20,0], 5, 10, 21, 'z', [0,0,1])

mag08 = rad.ObjRecCur([-15,0,0], [5,7,15], [0.5,0.5,1.])

mag09 = rad.ObjArcCur([0,0,-5], [10,13], [0,2.5], 10, 15, 1.7, 'man', 'z')

mag10 = rad.ObjRaceTrk([0,0,0], [27,28], [1,2.5], 5, 15, 1.7, 'man', 'z')

mag11 = rad.ObjFlmCur([[-10,-30,-10],[30,-30,-10],[30,25,25],[-30,25,25],[-30,-30,-10]], 10.2)

magBkg = rad.ObjBckg([1,2,3])

mag = rad.ObjCnt([mag00, mag01, mag02, mag03, mag04, mag05, mag06, mag07, mag08, mag09])
print('Container Content:', rad.ObjCntStuf(mag))
print('Container Size:', rad.ObjCntSize(mag))

rad.ObjAddToCnt(mag, [mag10, mag11])
cnt02 = rad.ObjCnt([mag00, mag])

mat = rad.MatLin([1.01, 1.2], [0, 0, 1.3])
#mat = rad.MatStd('NdFeB', 1.2)
rad.MatApl(mag01, mat)
print('Magn. Material index:', mat, ' appled to object:', mag01)

mag00a = rad.ObjFullMag([10,0,40],[12,18,5],[0,0,1],[2,2,2],cnt02,mat,[0.5,0,0])

rad.ObjDrwOpenGL(cnt02)

data_cnt = rad.ObjDrwVTK(cnt02)
print(data_cnt)


objAfterCut = rad.ObjCutMag(mag00a,[10,0,40],[1,1,1]) #,'Frame->Lab')
print('Indexes of objects after cutting:', objAfterCut)
#rad.ObjDrwOpenGL(objAfterCut[0])

print(rad.UtiDmp(mag01, 'asc'))
print(rad.UtiDmp(mat, 'asc'))
#print(rad.UtiDmp(107, 'asc'))

magDpl = rad.ObjDpl(mag, 'FreeSym->False')

print('Number of objects in the container:', rad.ObjCntSize(mag))
print('Number of objects in 2nd container:', rad.ObjCntSize(cnt02))
print('Number of objects in fake container:', rad.ObjCntSize(mag04))

print('Indices of elements in the container:', rad.ObjCntStuf(mag))
print('Indices of elements in the duplicated container:', rad.ObjCntStuf(magDpl))

#rad.ObjDrwOpenGL(mag)
#rad.ObjDrwOpenGL(magDpl)

#mag01sbd = rad.ObjDivMag(mag01, [[2,0.5],[3,0.2],[4,0.1]], 'pln', [[1,0.4,0.1],[0.4,1,0.2],[0,0,1]], 'Frame->Lab')

#mag01sbd = rad.ObjDivMag(mag01, [[2,0.5],[3,0.2],[4,0.1]], 'Frame->Lab')
#rad.ObjDrwOpenGL(mag01sbd)

#mag01sbd = rad.ObjDivMagPln(mag01, [[2,0.5],[3,0.2],[4,0.1]], [1,0.4,0.1], [0.4,1,0.2], [0,0,1], 'Frame->Lab')
mag01sbd = rad.ObjDivMagPln(mag01, [[2,0.5],[3,0.2],[4,0.1]])
#rad.ObjDrwOpenGL(mag01sbd)

#mag00sbd = rad.ObjDivMag(mag00, [[2,0.5],[3,0.2],[4,0.1]], 'cyl', [[2.5,4,0],[0,0,1],[8,0,0],3], 'Frame->Lab')
#mag00sbd = rad.ObjDivMagCyl(mag00, [[2,0.5],[3,0.2],[4,0.1]], [2.5,4,0], [0,0,1], [8,0,0], 3, 'Frame->Lab')

print('Volume of 3D object:', rad.ObjGeoVol(mag01sbd))
print('Geom. Limits of 3D object:', rad.ObjGeoLim(mag01sbd))

#rad.ObjDrwOpenGL(mag01)

trf01 = rad.TrfPlSym([0,10,0], [0,1,0])
trf02 = rad.TrfRot([0,10,0], [0,0,1], 1.)
trf03 = rad.TrfTrsl([30,10,0])
trf04 = rad.TrfInv()
trf05 = rad.TrfCmbL(trf01, trf04)
trf06 = rad.TrfCmbR(trf01, trf04)

#rad.TrfMlt(mag01, trf03, 3)
rad.TrfOrnt(mag01, trf06)

#rad.ObjDrwOpenGL(mag01)

matNdFeB = rad.MatStd('NdFeB')
M = rad.MatMvsH(matNdFeB, 'M', [0,0,0])
print('NdFeB material index:', matNdFeB, ' Magnetization:', M)

matLin01 = rad.MatLin([0.1,0.2],1.1)
matLin02 = rad.MatLin([0.1,0.2],[0,0,1.1])
print('Linear material indexes:', matLin01, matLin02)

dmp = rad.UtiDmp([mag01, trf02], 'bin')
#print(dmp)

elemsRest = rad.UtiDmpPrs(dmp)
print('Indexes of restored elements:', elemsRest)
#rad.ObjDrwOpenGL(elemsRest[0])

print(rad.UtiDmp(elemsRest[0], 'asc'))
print(rad.UtiDmp(elemsRest[1], 'asc'))

#rad.UtiDel(elemsRest[0])
#print(rad.UtiDmp(elemsRest[0], 'asc'))

print(rad.UtiDmp(trf03, 'asc'))
rad.UtiDelAll()
print(rad.UtiDmp(trf03, 'asc'))


#cutMag = rad.ObjCutMag(mag, [0,0,50], [0,0,1], "Frame->Lab")
#print('Obj. indexes after magnet cut:', cutMag)

#rad.ObjDrwOpenGL(cutMag[0])

#time.sleep(0.1)
#rad.ObjDrwOpenGL(cutMag[1])

#rad.ObjSetM(mag00, [1,2,3])
#M = rad.ObjM(mag00)
#print(M)

#H = rad.ObjCenFld(mag, 'A')
#print(H)

#B = rad.Fld(mag, "bha", [[0,0,0],[1,1,1]])
#print(B)

#print(rad.Fld(magDpl, "b", [1,1,1]))
