
#############################################################################
# RADIA Python Test Script
# v 0.01
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import radia as rad

print('RADIA Python Test Script 01')
print('')

mag00 = rad.ObjRecMag([5,0,0], [5,8,10], [-0.5,1,0.7])

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

mag = rad.ObjCnt([mag00, mag01, mag02, mag03, mag04, mag05, mag06, mag07, mag08, mag09])
rad.ObjAddToCnt(mag, [mag10, mag11])

rad.ObjDrwOpenGL(mag)

B = rad.Fld(mag, "b", [1,1,1])

print(B)
