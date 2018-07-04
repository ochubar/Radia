/*-------------------------------------------------------------------------
*
* File name:      radvlpgn.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 polyhedron with constant magnetization
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADVLPGN_H
#define __RADVLPGN_H

#include "radg3d.h"
#include "radplnr.h"
#include "radtrans.h"
#include "radcnvrg.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------

struct radTHandlePgnAndTrans {
	radTHandle<radTPolygon> PgnHndl;
	radTHandle<radTrans> TransHndl;
	bool FaceIsInternalAfterCut;

	radTHandlePgnAndTrans(radTHandle<radTPolygon>& inPgnHndl, radTHandle<radTrans>& inTransHndl) 
	{
		PgnHndl = inPgnHndl; TransHndl = inTransHndl;
	}
	radTHandlePgnAndTrans() { FaceIsInternalAfterCut = false;}
	inline friend int operator <(const radTHandlePgnAndTrans&, const radTHandlePgnAndTrans&);
	inline friend int operator ==(const radTHandlePgnAndTrans&, const radTHandlePgnAndTrans&);
};

//-------------------------------------------------------------------------

inline int operator <(const radTHandlePgnAndTrans&, const radTHandlePgnAndTrans&) { return 1;}

//-------------------------------------------------------------------------

inline int operator ==(const radTHandlePgnAndTrans& h1, const radTHandlePgnAndTrans& h2)
{
	return (h1.PgnHndl.rep == h2.PgnHndl.rep) && (h1.TransHndl.rep == h2.TransHndl.rep);
}

//-------------------------------------------------------------------------

#ifdef __GCC__
typedef vector<radTHandlePgnAndTrans> radTVectHandlePgnAndTrans;
typedef vector<TVector3d*> radTVectOfPtrToVect3d;
typedef list<TVector3d> radTListOfVector3d;
typedef vector<TVector3d> radTVectVect3d;
#else
typedef vector<radTHandlePgnAndTrans, allocator<radTHandlePgnAndTrans> > radTVectHandlePgnAndTrans;
typedef vector<TVector3d*, allocator<TVector3d*> > radTVectOfPtrToVect3d;
typedef list<TVector3d, allocator<TVector3d> > radTListOfVector3d;
typedef vector<TVector3d, allocator<TVector3d> > radTVectVect3d;
#endif

#ifdef __MWERKS__
/*
null_template
struct iterator_traits <TVector3d*> {
     typedef ptrdiff_t difference_type;
     typedef TVector3d value_type;
     typedef TVector3d* pointer;
     typedef TVector3d& reference;
     typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTHandlePgnAndTrans*> {
     typedef ptrdiff_t difference_type;
     typedef radTHandlePgnAndTrans value_type;
     typedef radTHandlePgnAndTrans* pointer;
     typedef radTHandlePgnAndTrans& reference;
     typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <TVector3d**> {
     typedef ptrdiff_t difference_type;
     typedef TVector3d* value_type;
     typedef TVector3d** pointer;
     typedef TVector3d*& reference;
     typedef random_access_iterator_tag iterator_category;
};
*/
#endif

//-------------------------------------------------------------------------

class radTPolyhedron : public radTg3dRelax {
	
	//const TMatrix3d* pJ_LinCoef;
	TMatrix3d* pJ_LinCoef;
	char mLinTreat; //0- treat as relative

public:
	radTVectHandlePgnAndTrans VectHandlePgnAndTrans;
	int AmOfFaces;

	TVector3d J; //to move to base?
	bool J_IsNotZero;

	short SomethingIsWrong;
	radTPairOfDouble AuxPairOfDouble; // Used for cylindrical subdivision

	radTPolyhedron(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, int lenArrayOfFaces, const TVector3d& InMagn) 
		: radTg3dRelax(InMagn)
	{
		//Magn = InMagn; AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		pJ_LinCoef = 0; mLinTreat = 0;
		J_IsNotZero = false;

		//DefineCentrPoint(ArrayOfPoints, lenArrayOfPoints); //OC090908
		ShiftFacesNumeration(ArrayOfFaces, ArrayOfLengths);
		FillInVectHandlePgnAndTrans(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfLengths);
		if(SomethingIsWrong) return;
		DefineCentrPoint(ArrayOfPoints, lenArrayOfPoints);
	}
	radTPolyhedron(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, int lenArrayOfFaces, 
		const TVector3d& InMagn, TMatrix3d& InM_LinCoef, TVector3d& InJ, TMatrix3d& InJ_LinCoef, char LinTreat) 
		: radTg3dRelax(InMagn, InM_LinCoef)
	{
		AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		//DefineCentrPoint(ArrayOfPoints, lenArrayOfPoints); //OC090908
		ShiftFacesNumeration(ArrayOfFaces, ArrayOfLengths);
		FillInVectHandlePgnAndTrans(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfLengths);
		if(SomethingIsWrong) return;
		DefineCentrPoint(ArrayOfPoints, lenArrayOfPoints);

		J = InJ;
		bool J_LinCoefIsNotZero = !InJ_LinCoef.isZero();
		J_IsNotZero = (!InJ.isZero()) || J_LinCoefIsNotZero; //??
		if(J_LinCoefIsNotZero) pJ_LinCoef = new TMatrix3d(InJ_LinCoef);
		else pJ_LinCoef = 0;

		mLinTreat = LinTreat;
	}
	radTPolyhedron(TVector3d** ArrayOfFaces, int* ArrayOfLengths, int lenArrayOfFaces, const TVector3d& InMagn) 
		: radTg3dRelax(InMagn)
	{
		//Magn = InMagn; AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		AmOfFaces = lenArrayOfFaces; SomethingIsWrong = 0;
		pJ_LinCoef = 0; mLinTreat = 0;
		J_IsNotZero = false;
		TVector3d* OutArrayOfPoints;
		int lenArrayOfPoints;
		int** OutArrayOfFaces;
		MakeNormalPresentation(ArrayOfFaces, ArrayOfLengths, OutArrayOfPoints, lenArrayOfPoints, OutArrayOfFaces);
		if(SomethingIsWrong) { DeleteInputArrays(OutArrayOfPoints, OutArrayOfFaces); return;}

		//DefineCentrPoint(OutArrayOfPoints, lenArrayOfPoints); //OC090908
		FillInVectHandlePgnAndTrans(OutArrayOfPoints, lenArrayOfPoints, OutArrayOfFaces, ArrayOfLengths);
		if(SomethingIsWrong) { DeleteInputArrays(OutArrayOfPoints, OutArrayOfFaces); return;}
		DefineCentrPoint(OutArrayOfPoints, lenArrayOfPoints);
		DeleteInputArrays(OutArrayOfPoints, OutArrayOfFaces);
	}
	radTPolyhedron(const radTVectHandlePgnAndTrans& InVectHandlePgnAndTrans, 
		const TVector3d* pInMagn, TMatrix3d* pInM_LinCoef, const radThg& InMatHandle, 
		const TVector3d* pInJ, TMatrix3d* pInJ_LinCoef, char LinTreat, const TVector3d* pPrevLinRefP) //used at cutting / subdivision
	//radTPolyhedron(const radTVectHandlePgnAndTrans& InVectHandlePgnAndTrans, const TVector3d& InMagn, const radThg& InMatHandle)
	//radTPolyhedron(const radTVectHandlePgnAndTrans& InVectHandlePgnAndTrans, const TVector3d& InCentrPoint, const TVector3d& InMagn, const radThg& InMatHandle)
		: radTg3dRelax(pInMagn, pInM_LinCoef, InMatHandle)
	{
		//pJ_LinCoef = 0; mLinTreat = 0;
		//J_IsNotZero = false;

		AmOfFaces = (int)InVectHandlePgnAndTrans.size();
		for(int i=0; i<AmOfFaces; i++) VectHandlePgnAndTrans.push_back(InVectHandlePgnAndTrans[i]);
		SomethingIsWrong = 0;
		DefineCentrPoint();

		J_IsNotZero = false;
		J.Zero();
		pJ_LinCoef = 0;
		if(pInJ != 0)
		{
			J = *pInJ;
			if(!J.isZero()) J_IsNotZero = true;
		}
		bool J_LinCoef_AreNotZero = false;
		if(pInJ_LinCoef != 0)
		{
			if(!pInJ_LinCoef->isZero())
			{
				pJ_LinCoef = new TMatrix3d(*pInJ_LinCoef);
				J_LinCoef_AreNotZero = true;
				J_IsNotZero = true;
			}
		}
		mLinTreat = LinTreat;

		if((LinTreat == 0) && (pPrevLinRefP != 0) && ((pM_LinCoef != 0) || J_LinCoef_AreNotZero)) //Linear dependence is Relative w.r. to object Center
		{//Do this correction only after M and/or J are defined!
		 //Consider moving this to base class(es)
			TVector3d dCenP = CentrPoint - (*pPrevLinRefP);
			if(pM_LinCoef != 0)
			{
				if(!pM_LinCoef->isZero()) Magn += ((*pM_LinCoef)*dCenP);
			}
			if(J_LinCoef_AreNotZero) J += ((*pJ_LinCoef)*dCenP);
		}
	}
	radTPolyhedron(const radTHandlePgnAndTrans& inHandleBasePgnAndTrf1, const radTHandlePgnAndTrans& inHandleBasePgnAndTrf2, double avgCur=0, double* arMagComp=0)
	{//tries to generate a convex polyhedron from two base face polygons
		SomethingIsWrong = 0;
		mLinTreat = 0; J_IsNotZero = false;
		J.x = J.y = J.z = 0.; pJ_LinCoef = 0; 
		Magn.x = Magn.y = Magn.z = 0.; pM_LinCoef = 0;

		AttemptToCreateConvexPolyhedronFromTwoBaseFaces(inHandleBasePgnAndTrf1, inHandleBasePgnAndTrf2);
		if(SomethingIsWrong) return;
		if(avgCur != 0)
		{
			SetCurrentDensityForConstCurrent(avgCur, 0, 1); //may set J and pJ_LinCoef
		}
		if(arMagComp != 0)
		{
			Magn.x = *arMagComp;
			Magn.y = *(arMagComp + 1);
			Magn.z = *(arMagComp + 2);
		}
	}
	radTPolyhedron(radTPolyhedron& aPlhdr) : radTg3dRelax(aPlhdr)
	{
		*this = aPlhdr;
		if(aPlhdr.pJ_LinCoef != 0) pJ_LinCoef = new TMatrix3d(*(aPlhdr.pJ_LinCoef));
	}
	radTPolyhedron(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{
		DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
		DumpBinParse_g3dRelax(inStr, mKeysOldNew, gMapOfHandlers);
		DumpBinParse_Polyhedron(inStr);
	}
	radTPolyhedron() : radTg3dRelax()
	{ 
		pJ_LinCoef = 0; mLinTreat = 0;
		J_IsNotZero = false;
		SomethingIsWrong = 0;
	}
	~radTPolyhedron() 
	{
		if(pJ_LinCoef != 0) delete pJ_LinCoef;
	}

	int Type_g3dRelax() { return 5;}
	int NumberOfDegOfFreedom() { return 3;}

	void FillInVectHandlePgnAndTrans(TVector3d*, int, int**, int*);
	void MakeNormalPresentation(TVector3d**, int*, TVector3d*&, int&, int**&);
	int CheckIfFacePolygonsArePlanar(TVector3d*, int**, int*, TVector3d*);
	int DetermineActualFacesNormals(TVector3d*, int, int**, int*, TVector3d*);
	int FillInTransAndFacesInLocFrames(TVector3d*, int**, int*, TVector3d*);

	//void B_comp(radTField*);
	//void B_intComp(radTField*);

	void B_comp_frM(radTField*);
	void B_comp_frJ(radTField*);
	void B_intComp_frM(radTField*);
	void B_intComp_frJ(radTField*);

	int CutItself(TVector3d*, radThg&, radTPair_int_hg&, radTPair_int_hg&, radTApplication*, radTSubdivOptions*);
	int FindIntersectionWithFace(int, TVector3d*, radTVectOfPtrToVect3d&, radTVectHandlePgnAndTrans&, radTVectHandlePgnAndTrans&, char&, double*);
	int SetUpUpperAndLowerPolygon(TVector2d*, int*, radTHandlePgnAndTrans&, radTVectHandlePgnAndTrans&, radTVectHandlePgnAndTrans&, double*);
	int DetermineNewFaceAndTrans(radTVectOfPtrToVect3d&, TVector3d&, radTHandlePgnAndTrans&, radTHandlePgnAndTrans&, double*);
	int FillInNewHandlePgnAndTransFrom3d(TVector3d**, int, TVector3d&, radTHandlePgnAndTrans&, double*);
	int CheckIfTwoPointAlreadyMapped(TVector3d&, TVector3d&, radTVectOfPtrToVect3d&, double*);
	void DefineRelAndAbsTol(double*);

	int SubdivideItselfByParPlanes(double*, int, radThg&, radTApplication*, radTSubdivOptions*);
	int KsFromSizeToNumb(double*, int, radTSubdivOptions*);

	int SubdivideItselfByOneSetOfParPlanes(TVector3d&, TVector3d*, int, radThg&, radTApplication*, radTSubdivOptions*, radTvhg*);
	int DeterminePointsOnCuttingPlanes(TVector3d&, double*, short, TVector3d*);
	int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions*, TVector3d&, TVector3d&, radTrans&, char&, char&);

	int CheckForSpecialShapes(radTVectHandlePgnAndTrans&, radThg&, double*);
	void IntrsctOfTwoLines(const TVector2d&, const TVector2d&, const TVector2d&, const TVector2d&, TVector2d&, TLinesIntrsctCase&, double*);

	double Volume();
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	void EstimateSize(TVector3d*, double*, int);
	void FindTypicalSize(TVector3d*, int, double&);

	radTg3dGraphPresent* CreateGraphPresent();
	void Dump(std::ostream&, int ShortSign =0);
	void DumpPureObjInfo(std::ostream&, int);
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey);
	void DumpBin_Polyhedron(CAuxBinStrVect& oStr);
	void DumpBinParse_Polyhedron(CAuxBinStrVect& inStr);

	int SubdivideItselfByEllipticCylinder(double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	void FindEdgePointsOverPhiAndAxForCylSubd(radTCylindricSubdivSpec*, TVector3d*, double*);
	int FindEdgePointsOverEllipseSet0(double*, radTCylindricSubdivSpec*, radThg, TVector3d*, double*, radTSubdivOptions*);
	int FindEdgePointsOverEllipseSet(double*, radTCylindricSubdivSpec*, radThg, TVector3d*, double*, radTSubdivOptions*);
	void FindLocalEllipticCoord(double, double, double, double&, double&);

	int SubdivideItselfOverAzimuth(double*, double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideByEllipses(double*, double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*);

	void Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl);
	
	void AttemptToCreateConvexPolyhedronFromTwoBaseFaces(const radTHandlePgnAndTrans& inHandleBasePgnAndTrf1, const radTHandlePgnAndTrans& inHandleBasePgnAndTrf2);
	bool CheckIfAllPolygonVerticesAreOnOneSideOfPlane(const radTHandlePgnAndTrans& hPgnAndTrf, const TVector3d& vPoint, const TVector3d& vNorm, double AbsTol);
	void CollectAndMapUniquePolygonPoints(const radTHandlePgnAndTrans& hPgnAndTrf, vector<TVector3d>& vectPoints, vector<int>& vectInd, double AbsTol);
	void GenerateSideFacesContainingSegmentsOfBaseFace(const vector<TVector3d>& vectVertexPoints, vector<vector<int> >& vectIndAllFaces, int indBaseFace, double* arTol);
	void ReorderPointsToEnsureNonSelfIntersectingPolygon(const vector<TVector3d>& vectPoints, vector<int>& vectIndPgnPoints, double* arTol);
	void SetCurrentDensityForConstCurrent(double avgCur, int indBaseFace1, int indBaseFace2);

	int ConvertToPolyhedron(radThg&, radTApplication*, char) { return 1;} // 1 is essential
	void B_comp(radTField* pField)
	{
		bool M_IsNotZero = !Magn.isZero();
		if(M_IsNotZero || (pField->FieldKey.PreRelax_)) B_comp_frM(pField);
		if(J_IsNotZero) B_comp_frJ(pField);
	}
	void B_intComp(radTField* pField)
	{
		bool M_IsNotZero = !Magn.isZero();
		if(M_IsNotZero) B_intComp_frM(pField);
		if(J_IsNotZero) B_intComp_frJ(pField);
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTPolyhedron(*this), hg);
	}
	int SizeOfThis()
	{
		int GenSize = sizeof(radTPolyhedron);
		int BufSize = sizeof(radTrans);
		GenSize += AmOfFaces*BufSize;
		for(int i=0; i<AmOfFaces; i++) GenSize += (VectHandlePgnAndTrans[i].PgnHndl.rep)->SizeOfThis();
		return GenSize;
	}
	void DefineCentrPoint(TVector3d* ArrayOfPoints, int AmOfPoints)
	{// Modify later if necessary
		TVector3d Sum(0.,0.,0.);
		for(int k=0; k<AmOfPoints; k++) Sum = Sum + ArrayOfPoints[k];
		double Buf = 1./AmOfPoints;
		CentrPoint = TVector3d(radCR.Double(Buf*Sum.x), radCR.Double(Buf*Sum.y), radCR.Double(Buf*Sum.z));
	}
	void DefineCentrPoint()
	{// This algorithm differs from the above (and may give different results)
		TVector3d Sum(0.,0.,0.);
		for(int k=0; k<AmOfFaces; k++)
		{
			radTHandlePgnAndTrans& HandlePgnAndTrans = VectHandlePgnAndTrans[k];
			radTPolygon* PgnPtr = HandlePgnAndTrans.PgnHndl.rep;
			TVector2d& LocFaceCP = PgnPtr->CentrPoint;
			TVector3d FaceCP(LocFaceCP.x, LocFaceCP.y, PgnPtr->CoordZ);
			Sum = Sum + HandlePgnAndTrans.TransHndl.rep->TrBiPoint(FaceCP);
		}
		double Buf = 1./double(AmOfFaces);
		CentrPoint = TVector3d(radCR.Double(Buf*Sum.x), radCR.Double(Buf*Sum.y), radCR.Double(Buf*Sum.z));
	}

	//void ReCalcCentrPointFromPgnAndTrans(const radTVectHandlePgnAndTrans& vHandlePgnAndTrans, const TVector3d& oldCenPoint, TVector3d& newCenPoint) //OC090908
	//{// This algorithm differs from the above (and may give different results)
	//	TVector3d Sum(0.,0.,0.);
	//	int locAmOfFaces = (int)vHandlePgnAndTrans.size();
	//	if(locAmOfFaces <= 0) return;
	//	for(int k=0; k<locAmOfFaces; k++)
	//	{
	//		const radTHandlePgnAndTrans& hPgnAndTrans = vHandlePgnAndTrans[k];
	//		radTPolygon* pPgn = hPgnAndTrans.PgnHndl.rep;
	//		TVector2d& LocFaceCP = pPgn->CentrPoint;
	//		TVector3d FaceCP(LocFaceCP.x, LocFaceCP.y, pPgn->CoordZ);
	//		Sum = Sum + hPgnAndTrans.TransHndl.rep->TrBiPoint(FaceCP);
	//	}
	//	double Buf = 1./double(locAmOfFaces);
	//	//CentrPoint = TVector3d(radCR.Double(Buf*Sum.x), radCR.Double(Buf*Sum.y), radCR.Double(Buf*Sum.z));
	//	newCenPoint.x = radCR.Double(Buf*Sum.x); //??
	//	newCenPoint.y = radCR.Double(Buf*Sum.y); 
	//	newCenPoint.z = radCR.Double(Buf*Sum.z);
	//	newCenPoint += oldCenPoint;
	//}

	//void CorrectFacePolygonsForNewCenPoint(radTVectHandlePgnAndTrans& vHandlePgnAndTrans, const TVector3d& difCenPoints) //OC090908
	//{//difCenPoints = Old - New
	//	int locAmOfFaces = (int)vHandlePgnAndTrans.size();
	//	if(locAmOfFaces <= 0) return;
	//	for(int k=0; k<locAmOfFaces; k++)
	//	{
	//		radTHandlePgnAndTrans& hPgnAndTrans = vHandlePgnAndTrans[k];
	//		radTPolygon* pPgn = hPgnAndTrans.PgnHndl.rep;
	//		radTrans* pTrans = hPgnAndTrans.TransHndl.rep;
	//		TVector3d addVertexPoint = pTrans->TrBiPoint(difCenPoints);
	//		pPgn->CoordZ += addVertexPoint.z;
	//		TVector2d addVertexPoint2d(addVertexPoint.x, addVertexPoint.y);
	//		pPgn->CentrPoint += addVertexPoint2d;
	//		radTVect2dVect& curVectEdgePoints = pPgn->EdgePointsVector;
	//		for(int j=0; j<pPgn->AmOfEdgePoints; j++)
	//		{
	//			curVectEdgePoints[j] += addVertexPoint2d;
	//		}
	//	}
	//}

	void ShiftFacesNumeration(int** ArrayOfFaces, int* ArrayOfLengths)
	{
		for(int i=0; i<AmOfFaces; i++)
		{
			int* CurrentFace = ArrayOfFaces[i];
			for(int j=0; j<ArrayOfLengths[i]; j++) (CurrentFace[j])--;
		}
	}
	void DefineNormalVia3Points(const TVector3d& P1, const TVector3d& P2, const TVector3d& P3, TVector3d& Normal)
	{
		TVector3d R1 = P2 - P1, R2 = P3 - P1;
		Normal.x = R1.y*R2.z - R2.y*R1.z;
		Normal.y = R2.x*R1.z - R1.x*R2.z;
		Normal.z = R1.x*R2.y - R2.x*R1.y;
	}
	double Vect3dNorm(const TVector3d& V)
	{
		double AbsX = Abs(V.x), AbsY = Abs(V.y), AbsZ = Abs(V.z);
		double MaxXY = (AbsX>AbsY)? AbsX : AbsY;
		return (MaxXY>AbsZ)? MaxXY : AbsZ;
	}
	int NextCircularNumber(int CurrentNo, int Total)
	{
		return (CurrentNo == Total-1)? 0 : CurrentNo + 1;
	}
	void ReverseArrayOfInt(int* ArrayOfInt, int lenArrayOfInt)
	{
		int *DirPtr = ArrayOfInt, *RevPtr = &(ArrayOfInt[lenArrayOfInt-1]);
		for(int i=0; i < (lenArrayOfInt >> 1); i++)
		{
			int Buf = *DirPtr; *(DirPtr++) = *RevPtr; *(RevPtr--) = Buf;
		}
	}
	void ReverseArrayOfVect3dPtr(TVector3d** ArrayOfVect3dPtr, int lenArrayOfVect3dPtr)
	{
		TVector3d** DirPtr = ArrayOfVect3dPtr;
		TVector3d** RevPtr = &(ArrayOfVect3dPtr[lenArrayOfVect3dPtr-1]);
		for(int i=0; i < (lenArrayOfVect3dPtr >> 1); i++)
		{
			TVector3d* Buf = *DirPtr; *(DirPtr++) = *RevPtr; *(RevPtr--) = Buf;
		}
	}
	int CheckIfJunctionIsConvex(const TVector3d& N1, const TVector3d& JointSegm, const TVector3d& N2)
	{
		TVector3d N1_vect_by_N2(N1.y*N2.z - N2.y*N1.z, N2.x*N1.z - N1.x*N2.z, N1.x*N2.y - N2.x*N1.y);
		return (N1_vect_by_N2*JointSegm >= 0.)? 1 : 0; // Or ">"?
	}
	void DeleteInputArrays(TVector3d* ArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths =NULL)
	{
		if(ArrayOfFaces != NULL)
		{
			for(int i=0; i<AmOfFaces; i++) delete[] (ArrayOfFaces[i]);
			delete[] ArrayOfFaces;
		}
		if(ArrayOfLengths != NULL) delete[] ArrayOfLengths;
		if(ArrayOfPoints != NULL) delete[] ArrayOfPoints;

		if(pJ_LinCoef != 0) delete pJ_LinCoef;
		pJ_LinCoef = 0;
	}
	void DeleteAuxInputArrays(TVector3d** ArrayOfFaces)
	{
		if(ArrayOfFaces != NULL)
		{
			for(int i=0; i<AmOfFaces; i++) delete[] (ArrayOfFaces[i]);
			delete[] ArrayOfFaces;
		}
	}
	void DeleteAuxInputArrays(short** ArrayOfFaces)
	{
		if(ArrayOfFaces != NULL)
		{
			for(int i=0; i<AmOfFaces; i++) delete[] (ArrayOfFaces[i]);
			delete[] ArrayOfFaces;
		}
	}
	int ItemIsNotFullyInternalAfterCut()
	{
		for(int k=0; k<AmOfFaces; k++)
			if(!VectHandlePgnAndTrans[k].FaceIsInternalAfterCut) return 1;
		return 0;
	}
	int CreateNewEntity(radTVectHandlePgnAndTrans& vHandlePgnAndTrans, radThg& hg, short RecognizeRecMagsInPolyhedrons, double* RelAbsTol)
	{
		short CreateA_Polyhedron = 1;
		if(RecognizeRecMagsInPolyhedrons)
			if(CheckForSpecialShapes(vHandlePgnAndTrans, hg, RelAbsTol)) CreateA_Polyhedron = 0;
		if(CreateA_Polyhedron)
		{
			radTSend Send;
			radTPolyhedron* PolyhedronPtr = new radTPolyhedron(vHandlePgnAndTrans, &Magn, pM_LinCoef, MaterHandle, &J, pJ_LinCoef, mLinTreat, &CentrPoint);
			//radTPolyhedron* PolyhedronPtr = new radTPolyhedron(vHandlePgnAndTrans, Magn, MaterHandle);

			if(PolyhedronPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
			hg = radThg(PolyhedronPtr);
		}
		return 1;
	}
	int FindTwoOrtogonalVectors(TVector3d& InV, TVector3d* TwoVect)
	{
		TVector3d V1(0.,0.,0.);
		if(!((InV.x==0.) && (InV.y==0.))) { V1.x = -InV.y; V1.y = InV.x;}
		else if(!((InV.x==0.) && (InV.z==0.))) { V1.x = -InV.z; V1.z = InV.x;}
		else if(!((InV.y==0.) && (InV.z==0.))) { V1.y = -InV.z; V1.z = InV.y;}
		else return 0;
		TVector3d V2 = InV^V1;
		*TwoVect = V1; *(TwoVect+1) = V2; 
		return 1;
	}
	char CheckIfOnlyNeighbouringEdgePointsTrapped(int* IntersectingBoundsNos, int AmOfEdgePo) //OC291003
	{
		if((AmOfEdgePo <= 0) || (IntersectingBoundsNos == 0)) return 0;
		int FirstNo = IntersectingBoundsNos[0], SecondNo = IntersectingBoundsNos[1];

		if((FirstNo == SecondNo) || (SecondNo == NextCircularNumber(FirstNo, AmOfEdgePo)) || (FirstNo == NextCircularNumber(SecondNo, AmOfEdgePo))) return 1;
		else return 0;
	}
	int ScaleCurrent(double scaleCoef) //OC250713 //virtual in g3d
	{//note: if(scaleCoef == 0) this still doesn't change J_IsNotZero
		if(J_IsNotZero) 
		{
			J *= scaleCoef; 

			if(pJ_LinCoef != 0)
			{
				*pJ_LinCoef *= scaleCoef; //??
			}
			return 1;
		}
		else return 0;
	}
};

//-------------------------------------------------------------------------

#endif
