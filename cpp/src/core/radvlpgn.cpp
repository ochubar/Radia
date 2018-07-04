/*-------------------------------------------------------------------------
*
* File name:      radvlpgn.cpp
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

#include "radsend.h"
#include "radvlpgn.h"
#include "radg3dgr.h"
#include "radsbdvp.h"
#include "radg3da1.h"
#include "radappl.h"
#include "auxparse.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTPolyhedron::FillInVectHandlePgnAndTrans(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths)
{
	radTSend Send;

	TVector3d* ArrayOfFacesNormals = new TVector3d[AmOfFaces];
	if(ArrayOfFacesNormals == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return;}

	if(!CheckIfFacePolygonsArePlanar(ArrayOfPoints, ArrayOfFaces, ArrayOfLengths, ArrayOfFacesNormals)) return;
	if(!DetermineActualFacesNormals(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfLengths, ArrayOfFacesNormals)) return;
	if(!FillInTransAndFacesInLocFrames(ArrayOfPoints, ArrayOfFaces, ArrayOfLengths, ArrayOfFacesNormals)) return;

	delete[] ArrayOfFacesNormals;
}

//-------------------------------------------------------------------------

void radTPolyhedron::MakeNormalPresentation(TVector3d** ArrayOfFaces, int* ArrayOfLengths, TVector3d*& OutArrayOfPoints, int& lenArrayOfPoints, int**& OutArrayOfFaces)
{
	radTVectVect3d VectOfPoints;
	radTSend Send;

	OutArrayOfFaces = 0;
	OutArrayOfFaces = new int*[AmOfFaces];
	if(OutArrayOfFaces == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); DeleteAuxInputArrays(ArrayOfFaces); return;}

	int PointsCount = 0;
	for(int i=0; i<AmOfFaces; i++)
	{
		TVector3d* CurrentFace = ArrayOfFaces[i];

		int* CurrentFaceInt = new int[ArrayOfLengths[i]];
		if(CurrentFaceInt == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); DeleteAuxInputArrays(ArrayOfFaces); return;}
		OutArrayOfFaces[i] = CurrentFaceInt;

		for(int k=0; k<ArrayOfLengths[i]; k++)
		{
			TVector3d& CurrentPoint = CurrentFace[k];
			short PointAlreadyThere = 0;
			int LocCount = -1;
			for(radTVectVect3d::iterator Iter = VectOfPoints.begin(); Iter != VectOfPoints.end(); ++Iter)
			{
				++LocCount;
				TVector3d& PointInCont = *Iter;
				if((PointInCont.x == CurrentPoint.x) && (PointInCont.y == CurrentPoint.y) && (PointInCont.z == CurrentPoint.z))
				{
					PointAlreadyThere = 1; break;
				}
			}
			if(PointAlreadyThere) CurrentFaceInt[k] = LocCount;
			else
			{
				VectOfPoints.push_back(CurrentPoint);
				CurrentFaceInt[k] = PointsCount;
				++PointsCount;
			}
		}
	}
	lenArrayOfPoints = PointsCount;
	OutArrayOfPoints = new TVector3d[lenArrayOfPoints];
	if(OutArrayOfPoints == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); DeleteAuxInputArrays(ArrayOfFaces); return;}
	for(int j=0; j<lenArrayOfPoints; j++)
	{
		OutArrayOfPoints[j] = VectOfPoints[j];
	}
}

//-------------------------------------------------------------------------

int radTPolyhedron::CheckIfFacePolygonsArePlanar(TVector3d* ArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, TVector3d* ArrayOfFacesNormals)
{
	radTSend Send;
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	const double RelLengthTol = RelAbsTol[0];

	for(int i=0; i<AmOfFaces; i++)
	{
		int* CurrentFace = ArrayOfFaces[i];

		TVector3d& P0 = ArrayOfPoints[CurrentFace[0]];
		TVector3d& P1 = ArrayOfPoints[CurrentFace[1]];
		TVector3d Normal;

		double BufAbsLenTol0 = RelLengthTol*Vect3dNorm(P0);
		double BufAbsLenTol1 = RelLengthTol*Vect3dNorm(P1);
		double AbsLenTol = (BufAbsLenTol1>BufAbsLenTol0)? BufAbsLenTol1 : BufAbsLenTol0;

		int CurrentLength = ArrayOfLengths[i];
		if(CurrentLength < 3) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error048"); return 0;}
		
		int j;
		short NormalDetermined = 0;
		for(j=2; j<CurrentLength; j++)
		{
			TVector3d& P2 = ArrayOfPoints[CurrentFace[j]];

			BufAbsLenTol0 = RelLengthTol*Vect3dNorm(P2);
			if(BufAbsLenTol0>AbsLenTol) AbsLenTol = BufAbsLenTol0;

			double AbsLenTol_e2 = AbsLenTol*AbsLenTol;

			DefineNormalVia3Points(P0, P1, P2, Normal);
			if(Abs(Normal.x)>AbsLenTol_e2 || Abs(Normal.y)>AbsLenTol_e2 || Abs(Normal.z)>AbsLenTol_e2)
			{
				NormalDetermined = 1; break;
			}
		}
		if(!NormalDetermined) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error049"); return 0;}

		for(int k=j+1; k<CurrentLength; k++)
		{
			TVector3d& aPoint = ArrayOfPoints[CurrentFace[k]];
			TVector3d anR = aPoint - P0;

			double AbsBufVal = Abs(Normal*anR);
			double CompareValue = Vect3dNorm(Normal)*AbsLenTol;

			if(AbsBufVal > CompareValue) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error047"); return 0;}
		}
		ArrayOfFacesNormals[i] = Normal;
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTPolyhedron::FindTypicalSize(TVector3d* ArrayOfPoints, int AmOfPoints, double& TypicalSize)
{
	double MaxX = 0., MaxY = 0., MaxZ = 0.;
	TVector3d& FirstPo = *ArrayOfPoints;
	for(int i=1; i<AmOfPoints; i++)
	{
		TVector3d& aPoint = ArrayOfPoints[i];
		double rX = fabs(aPoint.x - FirstPo.x), rY = fabs(aPoint.y - FirstPo.y), rZ = fabs(aPoint.z - FirstPo.z);
		if(rX > MaxX) MaxX = rX;
		if(rY > MaxY) MaxY = rY;
		if(rZ > MaxZ) MaxZ = rZ;
	}
	TypicalSize = (MaxX > MaxY)? MaxX : MaxY;
	if(TypicalSize < MaxZ) TypicalSize = MaxZ;
}

//-------------------------------------------------------------------------

int radTPolyhedron::DetermineActualFacesNormals(TVector3d* ArrayOfPoints, int AmOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, TVector3d* ArrayOfFacesNormals)
{
	radTSend Send;
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	const double RelLengthTol = RelAbsTol[0];
	double TypicalSize;
	FindTypicalSize(ArrayOfPoints, AmOfPoints, TypicalSize);
	double AbsLengthTol = RelLengthTol*TypicalSize;

	short GoodFaceFound = 0, NormalJustReversed = 0;

// Looking for the first good face
	int i;
	for(i=0; i<AmOfFaces; i++)
	{
		int* CurrentFace = ArrayOfFaces[i];
		TVector3d& Normal = ArrayOfFacesNormals[i];

		double InvNorm = 1./sqrt(Normal.x*Normal.x + Normal.y*Normal.y + Normal.z*Normal.z);
		Normal = InvNorm*Normal;

		TVector3d& P0 = ArrayOfPoints[CurrentFace[0]];

		short AtLeastOnePointIsOK = 0, SomePointsAreNotOK = 0;
		NormalJustReversed = 0;

		int j;
		for(j=0; j<AmOfPoints; j++)
		{
			TVector3d anR = P0 - ArrayOfPoints[j];
			double ScalProd = Normal*anR;

			if(Abs(ScalProd) > AbsLengthTol)
			{
				if(ScalProd < -AbsLengthTol) 
				{
					if(!AtLeastOnePointIsOK)
					{
						Normal = (-1.)*Normal; NormalJustReversed = 1; 
					}
					else SomePointsAreNotOK = 1;
					break;
				}
				else AtLeastOnePointIsOK = 1;
			}
		}
		if(NormalJustReversed)
			for(j=0; j<AmOfPoints; j++)
			{
				TVector3d anR = P0 - ArrayOfPoints[j];
				double ScalProd = Normal*anR;
				if(Normal*anR < -AbsLengthTol) { SomePointsAreNotOK = 1; break;}
			}
		if(!SomePointsAreNotOK) { GoodFaceFound = 1; break;}
	}
	if(!GoodFaceFound) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error050"); return 0;}
	int NoOfFirstGoodFace = i;

// Determining all the rest normals
	short** SegmentPassed = new short*[AmOfFaces];
	if(SegmentPassed == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return 0;}
	for(i=0; i<AmOfFaces; i++)
	{
		int CurrentLength = ArrayOfLengths[i];
		SegmentPassed[i] = new short[CurrentLength];
		short* LocSegmentPassed = SegmentPassed[i];
		if(LocSegmentPassed == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return 0;}
		for(int k=0; k<CurrentLength; k++) LocSegmentPassed[k] = 0;
	}
	int ii = NoOfFirstGoodFace;

	char* GenFacesPassed = new char[AmOfFaces];
	if(GenFacesPassed == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return 0;}
	char* PossibleNextFaces = new char[AmOfFaces];
	if(PossibleNextFaces == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return 0;}
	for(int p=0; p<AmOfFaces; p++) GenFacesPassed[p] = 0;

	for(i=0; i<AmOfFaces; i++)
	{
		for(int p=0; p<AmOfFaces; p++) PossibleNextFaces[p] = 0;

		int CurrentLength = ArrayOfLengths[ii];
		int* CurrentFace = ArrayOfFaces[ii];

		TVector3d& FirstNormal = ArrayOfFacesNormals[ii];
		if(NormalJustReversed)
		{ 
			ReverseArrayOfInt(CurrentFace, CurrentLength); NormalJustReversed = 0;
		}

		short* SegmentPassedOnCurrentFace = SegmentPassed[ii];
		for(int j=0; j<CurrentLength; j++)
		{
			int NoOfSegmStPo = CurrentFace[j]; 
			int NoOfSegmFiPo = CurrentFace[NextCircularNumber(j, CurrentLength)];
			
			if(!(SegmentPassedOnCurrentFace[j]))
			{
				short ThisSegmentIsSomeWhereElse = 0;
				for(int iii=0; iii<AmOfFaces; iii++)
				{
					if(iii != ii)
					{
						int LocCurrentLength = ArrayOfLengths[iii];
						int* LocCurrentFace = ArrayOfFaces[iii];

						int jjSt = -1, jjFi = -1;
						for(int jj=0; jj<LocCurrentLength; jj++)
						{
							int CurrentNo = LocCurrentFace[jj];
							if(CurrentNo == NoOfSegmStPo)
							{
								jjSt = jj; if(jjFi > -1) break;
							}
							if(CurrentNo == NoOfSegmFiPo)
							{
								jjFi = jj; if(jjSt > -1) break;
							}
						}
						if((jjSt>-1) && (jjFi>-1))
						{
							ThisSegmentIsSomeWhereElse = 1;
							short& ThisSegmentAlreadyPassed = (SegmentPassed[iii])[jjFi];
							if(ThisSegmentAlreadyPassed) break;

							PossibleNextFaces[iii] = 1;

							if(NextCircularNumber(jjFi, LocCurrentLength) != jjSt)
							{
								ReverseArrayOfInt(LocCurrentFace, LocCurrentLength);
								TVector3d& Normal = ArrayOfFacesNormals[iii];
								Normal = (-1.)*Normal;
							}
							TVector3d SegmVect = ArrayOfPoints[NoOfSegmFiPo] - ArrayOfPoints[NoOfSegmStPo];
							if(!CheckIfJunctionIsConvex(FirstNormal, SegmVect, ArrayOfFacesNormals[iii])) 
							{// Modify this if non-convex volumes are supported or treated specially
								SomethingIsWrong = 1; 
								Send.ErrorMessage("Radia::Error106"); 
								DeleteAuxInputArrays(SegmentPassed); return 0;
							}
							ThisSegmentAlreadyPassed = 1; break;
						}
					}
				}
				if(!ThisSegmentIsSomeWhereElse) // Modify this if unclosed volumes are not supported
				{
					Send.WarningMessage("Radia::Warning013"); 
				}
			}
		}
		GenFacesPassed[ii] = 1;

		//int NextFaceNo;
		int NextFaceNo = -1; //OC 100902
		for(int pp=0; pp<AmOfFaces; pp++)
		{
			if((PossibleNextFaces[pp] == 1) && (!GenFacesPassed[pp])) { NextFaceNo = pp; break;}
		}

		if(NextFaceNo >= 0)//OC 100902
		{
            ii = NextFaceNo;
		}
	}
	DeleteAuxInputArrays(SegmentPassed);

	if(GenFacesPassed != 0) delete[] GenFacesPassed;
	if(PossibleNextFaces != 0) delete[] PossibleNextFaces;
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::FillInTransAndFacesInLocFrames(TVector3d* ArrayOfPoints, int** ArrayOfFaces, int* ArrayOfLengths, TVector3d* ArrayOfFacesNormals)
{
	radTSend Send;
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	const double RelLengthTol = RelAbsTol[0];

	for(int i=0; i<AmOfFaces; i++)
	{
		TVector3d& N = ArrayOfFacesNormals[i];
		//double AbsLengthTol = RelLengthTol*Vect3dNorm(N); //OC
		double InvLen = 1./sqrt(N.x*N.x + N.y*N.y + N.z*N.z);
		N = InvLen*N;

		TVector3d St1, St2, St3;
		//if(Abs(N.z+1.) > AbsLengthTol) //OC
		if(Abs(N.z+1.) > RelLengthTol) //OC
		{
			double InvNzp1 = 1./(N.z + 1.);
			St1 = TVector3d(N.y*N.y*InvNzp1 + N.z, -N.x*N.y*InvNzp1, -N.x);
			St2 = TVector3d(St1.y, N.x*N.x*InvNzp1 + N.z, -N.y);
			St3 = TVector3d(-St1.z, -St2.z, N.z);
		}
		else
		{
			St1 = TVector3d(1., 0., 0.);
			St2 = TVector3d(0., -1., 0.);
			St3 = TVector3d(0., 0., -1.);
		}
		TMatrix3d R1(St1, St2, St3);

		int k;
		int* CurrentFace = ArrayOfFaces[i];
		int CurrentLength = ArrayOfLengths[i];

		//TVector3d EdgeVect; //OC
		//for(k=0; k<CurrentLength; k++)
		//{
		//	EdgeVect = ArrayOfPoints[CurrentFace[NextCircularNumber(k, CurrentLength)]] - ArrayOfPoints[CurrentFace[k]];
		//	double EdgeVectLength = sqrt(EdgeVect.x*EdgeVect.x + EdgeVect.y*EdgeVect.y + EdgeVect.z*EdgeVect.z);
		//	if(EdgeVectLength > AbsLengthTol)
		//	{
		//		EdgeVect = (1./EdgeVectLength)*EdgeVect; break;
		//	}
		//}
		//TVector3d V = R1*EdgeVect;

// Uncomment the following to speed-up a bit (however, loss of precision in FieldInt computation was encountered due to this)
		//St1.x = V.y; St1.y = -V.x; St1.z = 0.;
		//St2.x = V.x; St2.y = V.y; St2.z = 0.;
		//St3.x = 0.; St3.y = 0.; St3.z = 1.;
		//TMatrix3d R2(St1, St2, St3);
		//TMatrix3d R = R2*R1;

		TMatrix3d R = R1;
		TVector3d Zero(0.,0.,0.);

		radTrans* RotationPtr = new radTrans(R, Zero, 1., 1., 2);
		if(RotationPtr == 0) { SomethingIsWrong=1; Send.ErrorMessage("Radia::Error900"); return 0;}

		radTVect2dVect Vect2dVect;
		double LocCoordZ;
		for(k=0; k<CurrentLength; k++)
		{
			TVector3d P = ArrayOfPoints[CurrentFace[k]];
			//TVector3d P = ArrayOfPoints[CurrentFace[k]] - CentrPoint; //OC090908 assuming that CentrPoint was already defined
			//ATTENTION: this modification requires updates in all Field computation and Subdivision routines !!! //OC090908

			TVector3d P_loc = RotationPtr->TrPoint(P);
			TVector2d P2d(P_loc.x, P_loc.y);
			Vect2dVect.push_back(P2d);
			LocCoordZ = P_loc.z;
		}

		TVector3d LocMagn = RotationPtr->TrVectField(Magn);
		radTPolygon* FacePgnPtr = new radTPolygon(Vect2dVect, LocCoordZ, LocMagn);
		if(FacePgnPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
		radTHandle<radTPolygon> HandlePgn(FacePgnPtr);

		RotationPtr->Invert();
		radTHandle<radTrans> HandleRotat(RotationPtr);

		radTHandlePgnAndTrans HandlePgnAndTrans;
		HandlePgnAndTrans.PgnHndl = HandlePgn;
		HandlePgnAndTrans.TransHndl = HandleRotat;
		VectHandlePgnAndTrans.push_back(HandlePgnAndTrans);
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_comp_frM(radTField* FieldPtr)
//void radTPolyhedron::B_comp(radTField* FieldPtr)
{
	TVector3d Zero(0.,0.,0.);
	short PointIsInside = 1;

	//TVector3d OrigP = FieldPtr->P; //OC090908
	//FieldPtr->P -= CentrPoint; //OC090908

	radTFieldKey LocFieldKey = FieldPtr->FieldKey;
	if(LocFieldKey.B_) LocFieldKey.H_ = 1;
	radTField SumLocField(LocFieldKey, FieldPtr->CompCriterium, FieldPtr->P, Zero, Zero, Zero, Zero, Zero);

	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans HandlePgnAndTrans = VectHandlePgnAndTrans[i];

		radTPolygon* PgnPtr = HandlePgnAndTrans.PgnHndl.rep;
		radTrans* TransPtr = HandlePgnAndTrans.TransHndl.rep;

		radTField LocField(LocFieldKey, SumLocField.CompCriterium, Zero, Zero, Zero, Zero, Zero, Zero);
		LocField.P = TransPtr->TrPoint_inv(SumLocField.P);
		TVector3d PrevP = LocField.P;

		if(!LocFieldKey.PreRelax_) PgnPtr->Magn = TransPtr->TrVectField_inv(Magn);
		PgnPtr->B_comp(&LocField); //OC040504 test

		if(LocField.P != PrevP)
		{
			SumLocField.P = TransPtr->TrPoint(LocField.P); //OC040504 test
		}

		if(!LocField.PointIsInsideFrame) PointIsInside = 0;

		if(!LocFieldKey.PreRelax_) SumLocField += TransPtr->TrField(LocField);
		else
		{
			TMatrix3d Q(LocField.B, LocField.H, LocField.A);
			TransPtr->TrMatrixLeft_inv(Q);
			TransPtr->TrMatrix(Q);
			SumLocField.B += Q.Str0;
			SumLocField.H += Q.Str1;
			SumLocField.A += Q.Str2;
		}
	}

	//FieldPtr->P = OrigP; //OC090908

	radTFieldKey& FldKey = FieldPtr->FieldKey;
	if(FldKey.PreRelax_)
	{
		FieldPtr->B += SumLocField.B;
		FieldPtr->H += SumLocField.H;
		FieldPtr->A += SumLocField.A;
		return;
	}
	if(FldKey.H_) FieldPtr->H += SumLocField.H;
	if(FldKey.M_) if(PointIsInside) FieldPtr->M += Magn;
	if(FldKey.B_)
	{
		FieldPtr->B += SumLocField.H;
		if(PointIsInside) FieldPtr->B += Magn;
	}
	if(FldKey.A_) FieldPtr->A += SumLocField.A;
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_comp_frJ(radTField* pField)
{
	TVector3d vEx(1.,0.,0.), vEy(0.,1.,0.), vEz(0.,0.,1.);
	const double Pi = 3.14159265358979;
	const double relPrecSwitchRootDecomp = 1E-08;

	TVector3d &PobsLab = pField->P;
	TMatrix3d QT; //for linear terms of J in local frames of faces
	TVector3d vSumFaces0(0,0,0), vSumFaces1(0,0,0), vSumFaces3(0,0,0); 
	double sumFaces2 = 0;
	TVector2d vQTX, vQTY, vI2, vSegm, vSegmUnit, vSegmExtNorm, vPobsProjToVertex, vProjToVertex1, vProjToVertex2;
	TVector2d vShiftCenPointVertex, PobsProj;
	double qtXZ, qtYZ;
	double I1, partI3;
	double R1, R2, s1_p_R1, s2_p_R2, PiMult1;

	bool PointIsInside = true;

	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans hPgnAndTrans = VectHandlePgnAndTrans[i];
		radTPolygon* pPgn = hPgnAndTrans.PgnHndl.rep;
		radTrans* pTrans = hPgnAndTrans.TransHndl.rep;

		radTVect2dVect &vPgnVertices = pPgn->EdgePointsVector;

		TVector2d &pgnFirstP = vPgnVertices[0];
		double AbsRandX = radCR.AbsRandMagnitude(pgnFirstP.x - pPgn->CentrPoint.x);
		double AbsRandY = radCR.AbsRandMagnitude(pgnFirstP.y - pPgn->CentrPoint.y);
		double AbsRandXY = AbsRandX;
		if(AbsRandXY < AbsRandY) AbsRandXY = AbsRandY;
		double AbsRandZ = radCR.AbsRandMagnitude(pPgn->CoordZ);
		if(AbsRandZ < AbsRandXY) AbsRandZ = AbsRandXY;

		TVector3d Pobs = pTrans->TrPoint_inv(PobsLab);
		double hi = pPgn->CoordZ - Pobs.z; //to check

		// Artificial shift of an observation point
		// if the point is exactly on the border (to avoid "divide by zero" error):
		if(hi == 0.) 
		{
			hi = AbsRandZ;
		}

		if(hi < 0.) PointIsInside = false; //to check

		double hiE2 = hi*hi;
		double abs_hi = ::fabs(hi);

		//TVector2d PobsProj(Pobs.x, Pobs.y);
		PobsProj.x = Pobs.x; PobsProj.y = Pobs.y;
		bool PobsProjIsInside = true;

		int pgnAmOfVertices = pPgn->AmOfEdgePoints;
		int pgnAmOfVertices_mi_1 = pgnAmOfVertices - 1;
		int j2 = 1;
		I1 = 0.;
		partI3 = 0.;
		vI2.x = vI2.y = 0.;
		double sum_hij = 0.;

		for(int j=0; j<pgnAmOfVertices; j++)
		{
			if(j == pgnAmOfVertices_mi_1) j2 = 0;

			TVector2d &vPgnVertex1 = vPgnVertices[j];

			vSegm = vPgnVertices[j2++] - vPgnVertex1; //don't use j2 after this!
			vSegmUnit = vSegm;
			vSegmUnit.Normalize();
			vSegmExtNorm.x = vSegmUnit.y; vSegmExtNorm.y = -vSegmUnit.x;

			vPobsProjToVertex = vPgnVertex1 - PobsProj;
			//vPobsProjToVertex = (vPgnVertex1 - vShiftCenPointVertex) - PobsProj;
			double hij = vPobsProjToVertex*vSegmExtNorm;
			// Artificial shift of an observation point
			// if the point is exactly on the border (to avoid "divide by zero" error):
			if(hij == 0.)
			{
				hij = AbsRandXY;
				//slight displacement of the observation point is also required...?
			}
			double abs_hij = ::fabs(hij);
			double sign_hij = Sign(hij);

			if(hij < 0)
			{
				PobsProjIsInside = false;
			}
			sum_hij += hij;

			vProjToVertex1 = vPobsProjToVertex - (hij*vSegmExtNorm);
			vProjToVertex2 = vProjToVertex1 + vSegm;

			double s1 = vProjToVertex1*vSegmUnit;
			double s2 = vProjToVertex2*vSegmUnit;

			//double scalProdProjToVert = vProjToVertex1*vProjToVertex2;
			//if(scalProdProjToVert >= 0)
			//{
			if(s2 < s1) //??
			{
				double sBuf = s2;
				s2 = s1; s1 = sBuf;
			}
			//}
			//else 
			//{//to check
			//	s1 = -s1;
			//}

			double s1e2 = s1*s1;
			double s2e2 = s2*s2;
			double hiE2_p_hijE2 = hiE2 + hij*hij;

			if((hiE2_p_hijE2 < s1e2*relPrecSwitchRootDecomp) && (s1 < 0.))
			{
				s1_p_R1 = -0.5*hiE2_p_hijE2/s1;
				R1 = -s1 + s1_p_R1;
			}
			else
			{
				R1 = sqrt(hiE2_p_hijE2 + s1e2);
				s1_p_R1 = s1 + R1; 
			}

			if((hiE2_p_hijE2 < s2e2*relPrecSwitchRootDecomp) && (s2 < 0.))
			{
				s2_p_R2 = -0.5*hiE2_p_hijE2/s2;
				R2 = -s2 + s2_p_R2;
			}
			else
			{
				R2 = sqrt(hiE2_p_hijE2 + s2e2);
				s2_p_R2 = s2 + R2; 
			}

			//R1 = sqrt(hiE2_p_hijE2 + s1*s1);
			//R2 = sqrt(hiE2_p_hijE2 + s2*s2);
			//if(s1_p_R1 <= 0) 
			//{
			//	s1_p_R1 = AbsRandXY;
			//}
			//if(s2_p_R2 <= 0) 
			//{
			//	s2_p_R2 = AbsRandXY;
			//}

			double LogDif = log(s2_p_R2/(s1_p_R1)); 
			
			double ArgAtan1 = hi*s1/(hij*R1);
			double ArgAtan2 = hi*s2/(hij*R2);

			double SumAtan1 = atan(TransAtans(ArgAtan2, -ArgAtan1, PiMult1));
			SumAtan1 += Pi*PiMult1;
			double addAtanI1 = hi*SumAtan1 + hij*LogDif;

			I1 += addAtanI1;

			double multI2 = 0.5*(s2*R2 - s1*R1 + hiE2_p_hijE2*LogDif);
			vI2 += multI2*vSegmExtNorm;

			partI3 += hij*multI2;
		}

		double I1_0 = I1;
		if(PobsProjIsInside)
		{
			I1 -= 2*Pi*abs_hi; //!!!
		}

		TVector3d vFaceN = pTrans->TrBiPoint(vEz); //face normal in laboratory frame?
		TVector3d vExLab = pTrans->TrBiPoint(vEx);
		TVector3d vEyLab = pTrans->TrBiPoint(vEy);

		vSumFaces0 += I1*vFaceN;
		
		double hi_I1 = hi*I1; //required for A
		sumFaces2 += hi_I1;

		if(pJ_LinCoef != 0)
		{
			QT = *pJ_LinCoef;
			//pTrans->TrMatrixGeom(QT); //linear terms of J in local frame of the current face - to check!!
			//pTrans->TrMatrixGeomLeft_inv(QT); //linear terms of J in local frame of the current face - to check!!
			pTrans->TrMatrixGeomLeft(QT); //linear terms of J in local frame of the current face - to check!!
			pTrans->TrMatrixGeom_inv(QT); //linear terms of J in local frame of the current face - to check!!

			vQTX.x = QT.Str0.x; vQTX.y = QT.Str0.y;
			vQTY.x = QT.Str1.x; vQTY.y = QT.Str1.y;
			qtXZ = QT.Str0.z;
			qtYZ = QT.Str1.z;

			vSumFaces1 += ((vI2*vQTY)*vExLab) - ((vI2*vQTX)*vEyLab);
			vSumFaces1 += hi_I1*((qtYZ*vExLab) - (qtXZ*vEyLab));

			if(pField->FieldKey.A_)
			{
				double I3 = hi*hi_I1 + partI3;
				vSumFaces3 += I3*vFaceN;
			}
		}
	}

	const double ConstForJ = 0.0001;
	TVector3d Jmain = J;
	if((mLinTreat == 0) && (pJ_LinCoef != 0)) //treat as being relative 
	{
		Jmain -= ((*pJ_LinCoef)*CentrPoint);
	}
	if(pJ_LinCoef != 0)
	{
		if(pField->FieldKey.A_ || pField->FieldKey.B_ || pField->FieldKey.H_ || pField->FieldKey.J_)
		{
			Jmain += (*pJ_LinCoef)*PobsLab;
		}
	}

	if(pField->FieldKey.J_) 
	{
		if(PointIsInside) pField->J += Jmain;
	}
	if(pField->FieldKey.A_)
	{
		TVector3d vA(0,0,0);
		if(pJ_LinCoef != 0)
		{
			//Jmain += (*pJ_LinCoef)*PobsLab;
			vA += (*pJ_LinCoef)*((ConstForJ/3.)*vSumFaces3);
		}
		vA += (0.5*ConstForJ*sumFaces2)*Jmain;
		pField->A += vA;
	}
	if(pField->FieldKey.B_ || pField->FieldKey.H_)
	{
		TVector3d vB(0,0,0);
		if(pJ_LinCoef != 0)
		{
			//Jmain += (*pJ_LinCoef)*PobsLab;
			TVector3d vQc(pJ_LinCoef->Str1.z - pJ_LinCoef->Str2.y, 
						  pJ_LinCoef->Str2.x - pJ_LinCoef->Str0.z, 
						  pJ_LinCoef->Str0.y - pJ_LinCoef->Str1.x);
			vB += ConstForJ*(vSumFaces1 - ((0.5*sumFaces2)*vQc));
		}
		vB += ConstForJ*(Jmain^vSumFaces0);
		
		pField->B += vB;
		pField->H += vB;
	}
	//if(pField->FieldKey.J_) 
	//{
	//	if(PointIsInside) 
	//	{
	//		pField->J += J;
	//		if(pJ_LinCoef != 0)
	//		{
	//			pField->J += (mLinTreat == 0)? (*pJ_LinCoef)*(PobsLab - CentrPoint) : (*pJ_LinCoef)*PobsLab;
	//		}
	//	}
	//}
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_intComp_frM(radTField* FieldPtr)
//void radTPolyhedron::B_intComp(radTField* FieldPtr)
{ 
	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans HandlePgnAndTrans = VectHandlePgnAndTrans[i];

		radTPolygon* PgnPtr = HandlePgnAndTrans.PgnHndl.rep;
		radTrans* TransPtr = HandlePgnAndTrans.TransHndl.rep;

		radTField LocField(FieldPtr->FieldKey, FieldPtr->CompCriterium, Zero, Zero, Zero, Zero, Zero, Zero);

		PgnPtr->Magn = TransPtr->TrVectField_inv(Magn);

		LocField.P = TransPtr->TrPoint_inv(FieldPtr->P);
		LocField.NextP = TransPtr->TrPoint_inv(FieldPtr->NextP);
		//LocField.P = TransPtr->TrPoint_inv(FieldPtr->P - CentrPoint); //OC090908
		//LocField.NextP = TransPtr->TrPoint_inv(FieldPtr->NextP - CentrPoint); //OC090908

		PgnPtr->B_intComp(&LocField);
		*FieldPtr += TransPtr->TrField(LocField);
	}
}

//-------------------------------------------------------------------------

void radTPolyhedron::B_intComp_frJ(radTField* pField)
{
	if(pField->FieldKey.FinInt_) { B_intCompFinNum(pField); return;} //to check if this still works

	//- find rotation which transforms integration line/vector to Ez;
	//- for each face polygon:
		//* find coordinates of vertex points in the frame where the integration line is along Ez;
		//* calculate face contribution to the total field integral value (in the above frame);
	//- transform the total field integral vector to the laboratory frame

	const double Pi = 3.14159265358979;
	const double inv12 = 1./12.;
	double relTolPerp = 10.*radCR.RelRand; //1.E-11; //to tune
	TVector3d vIntAx = pField->NextP - pField->P, Ez(0,0,1.);
	TVector3d vIntAxUnit = vIntAx; vIntAxUnit.Normalize();

	radTrans trIntAxis2Ez, trAuxRot, trAuxTransl;
	trAuxRot.SetupRotation(CentrPoint, vIntAxUnit, Ez);

	TVector3d vIntAxCenPoint = 0.5*(pField->NextP + pField->P);
	TVector3d vTestIntAxCenPoint = trAuxRot.TrPoint(vIntAxCenPoint);
	TVector3d vTransl(-vTestIntAxCenPoint.x, -vTestIntAxCenPoint.y, 0.);
	//TVector3d vTransl(-vTestIntAxCenPoint.x, -vTestIntAxCenPoint.y, -vTestIntAxCenPoint.z);
	trAuxTransl.SetupTranslation(vTransl);
	TrProduct(&trAuxTransl, &trAuxRot, trIntAxis2Ez);

	TVector3d Jrot = trIntAxis2Ez.TrVectField(J); //to check this !!!
	double qxx, qxy, qxz;
	double qyx, qyy, qyz;
	double qzx, qzy, qzz;
	double qxxmqyy, qxypqyx, twoqxx, twoqzx, twoqzy;

	if(pJ_LinCoef != 0)
	{
		TMatrix3d Qrot = *pJ_LinCoef;
		trIntAxis2Ez.TrMatrixGeomLeft_inv(Qrot); //{ Qrot = Qrot*M_inv;} //to check
		trIntAxis2Ez.TrMatrix(Qrot); //{ Qrot = s*M*Qrot;}

		TVector3d vRefP(0,0,0);

		if(mLinTreat == 0) //treat as being relative 
		{
			vRefP = CentrPoint;
			//TVector3d vCenPtIntFrame = trIntAxis2Ez.TrPoint(CentrPoint);
			//Jrot -= (Qrot*vCenPtIntFrame); //to check
		}
		TVector3d vRefPIntFrame = trIntAxis2Ez.TrPoint(vRefP);
		Jrot -= (Qrot*vRefPIntFrame); //to check

		TVector3d &st0 = Qrot.Str0, &st1 = Qrot.Str1, &st2 = Qrot.Str2;
		qxx = st0.x; qxy = st0.y; qxz = st0.z;
		qyx = st1.x; qyy = st1.y; qyz = st1.z;
		qzx = st2.x; qzy = st2.y; qzz = st2.z;

		qxxmqyy = qxx - qyy; qxypqyx = qxy + qyx;
		twoqxx = 2.*qxx;
		twoqzx = 2.*qzx; twoqzy = 2.*qzy;
	}
	double axe2, aye2, cze2;
	double trecze2qzz, trecze2qxz, tricze2qyz, axe2qzz, axe2qxz, aye2qyz, axqzz, axqyz, ayqzz, ayqxz, ayqyz;
	double axqyzp2qyxp2qxy, qxxpaxqxzmqyy, qxxpaxqxzmqyymayqyz, qxypayqxzpqyxpaxqyz;
	double qzypayqzz, twoqzypayqzz, aytwoqzypayqzz, twoqxypayqxzptwoqyx, twoqxxmqyymayqyz, twoqxxpaxqxzmqyy;
	double twoax, twocz, tricz, tricze2, triczqzz;
	double czqxz, triczqxz, czqyz, triczqyz;
	double qzxpaxqzz, twoaxqzz;

	double jx0 = Jrot.x, jy0 = Jrot.y, jz0 = Jrot.z;
	double quartJz0 = 0.25*jz0;
	double twojx0 = 2.*jx0, trijx0 = 3.*jx0, twojy0 = 2.*jy0, trijy0 = 3.*jy0, twojz0 = 2.*jz0, trijz0 = 3*jz0;

	radTrans trFace2Int;
	TVector3d vR1loc, vR1int, vR2loc, vR2int, vSegm3d, vSegm3dUnit, vSegm3dExtNorm, vResLocIB(0,0,0);
	TVector2d vSegmExtNorm;
	double auxBuf, piMult1 = 0;
	double IBxLoc, IByLoc, IBzLoc;

	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans hPgnAndTrans = VectHandlePgnAndTrans[i];
		radTPolygon* pPgn = hPgnAndTrans.PgnHndl.rep;
		radTrans* pTrans = hPgnAndTrans.TransHndl.rep;
		TrProduct(&trIntAxis2Ez, pTrans, trFace2Int);

		TVector3d vFaceNormIntFrame = trFace2Int.TrBiPoint(Ez);
		//if this vector is perpendicular to Ez, then don't compute contribution from this face
		if(::fabs(vFaceNormIntFrame.z) < relTolPerp) continue; //to check consistency with the subsequent
		
		radTVect2dVect &vPgnVertices = pPgn->EdgePointsVector;
		TVector2d &vert0_2d = vPgnVertices[0];
		TVector3d vR0loc(vert0_2d.x, vert0_2d.y, pPgn->CoordZ);
		TVector3d vR0int = trFace2Int.TrPoint(vR0loc);

		double signZ_FaceNorm = Sign(vFaceNormIntFrame.z);
		double ax = -vFaceNormIntFrame.x/vFaceNormIntFrame.z;
		double ay = -vFaceNormIntFrame.y/vFaceNormIntFrame.z;
		double cz = vR0int.z - ax*vR0int.x - ay*vR0int.y;

		if(pJ_LinCoef != 0)
		{
			axe2 = ax*ax; aye2 = ay*ay; cze2 = cz*cz;
			twoax = 2.*ax;
			twocz = 2.*cz; tricz = 3.*cz; tricze2 = 3.*cze2;
			axe2qzz = axe2*qzz; axe2qxz = axe2*qxz; aye2qyz = aye2*qyz;
			trecze2qzz = tricze2*qzz; trecze2qxz = tricze2*qxz; tricze2qyz = tricze2*qyz;
			triczqzz = tricz*qzz;
			czqxz = cz*qxz; triczqxz = tricz*qxz;
			czqyz = cz*qyz; triczqyz = tricz*qyz; 
			axqzz = ax*qzz; axqyz = ax*qyz;
			ayqzz = ay*qzz; ayqxz = ay*qxz; ayqyz = ay*qyz;

			axqyzp2qyxp2qxy = axqyz + 2.*qyx + 2.*qxy;
			qxxpaxqxzmqyy = qxxmqyy + ax*qxz;
			qxxpaxqxzmqyymayqyz = qxxpaxqxzmqyy - ayqyz;
			qzxpaxqzz = qzx + axqzz; twoaxqzz = 2.*axqzz;
			qzypayqzz = qzy + ayqzz; twoqzypayqzz = 2.*qzy + ayqzz;
			aytwoqzypayqzz = ay*twoqzypayqzz;
			qxypayqxzpqyxpaxqyz = qxypqyx + ayqxz + axqyz;
			twoqxypayqxzptwoqyx = 2.*qxypqyx + ayqxz;
			twoqxxmqyymayqyz = 2.*(qxxmqyy - ayqyz);
			twoqxxpaxqxzmqyy = 2.*qxxpaxqxzmqyy;

		}

		int pgnAmOfVertices = pPgn->AmOfEdgePoints;
		int pgnAmOfVertices_mi_1 = pgnAmOfVertices - 1;
		int j2 = 1;
		vR1loc.z = pPgn->CoordZ;
		vR2loc.z = pPgn->CoordZ;

		double AbsRandX = radCR.AbsRandMagnitude(vert0_2d.x - pPgn->CentrPoint.x);
		double AbsRandY = radCR.AbsRandMagnitude(vert0_2d.y - pPgn->CentrPoint.y);
		double AbsRandXY = AbsRandX;
		if(AbsRandXY < AbsRandY) AbsRandXY = AbsRandY;

		TVector3d vFaceContribIB(0,0,0);
		//bool IntLineCrossesFace = true;

		for(int j=0; j<pgnAmOfVertices; j++)
		{
			if(j == pgnAmOfVertices_mi_1) j2 = 0;

			TVector2d &vPgnVertex1 = vPgnVertices[j];
			TVector2d &vPgnVertex2 = vPgnVertices[j2++]; //don't use j2 after this!

			vR1loc.x = vPgnVertex1.x; vR1loc.y = vPgnVertex1.y;
			vR2loc.x = vPgnVertex2.x; vR2loc.y = vPgnVertex2.y;

			vR1int = trFace2Int.TrPoint(vR1loc);
			vR2int = trFace2Int.TrPoint(vR2loc);
			vSegm3d = vR2int - vR1int;
			vSegm3dUnit = vSegm3d; vSegm3dUnit.Normalize();
			//vSegm3dExtNorm = vSegm3dUnit^vFaceNormIntFrame;
			vSegm3dExtNorm = vSegm3dUnit^(signZ_FaceNorm*Ez); //!!
			
			vSegmExtNorm.x = vSegm3dExtNorm.x; vSegmExtNorm.y = vSegm3dExtNorm.y;
			//vSegmExtNorm.Normalize();
			if(::fabs(vSegmExtNorm.y) < relTolPerp) continue; //to check consistency with the subsequent
			
			double signY_EdgeNorm = Sign(vSegmExtNorm.y);

			double x1 = vR1int.x, x2 = vR2int.x;
			double y1 = vR1int.y, y2 = vR2int.y;
			if(x1 > x2)
			{
				auxBuf = x1; x1 = x2; x2 = auxBuf;
				auxBuf = y1; y1 = y2; y2 = auxBuf;
			}
			else if(x1 == x2)
			{
				x1 -= AbsRandXY; //Artificial shift
			}
			if(x1 == 0.) x1 = -AbsRandXY; //Artificial shift
			if(x2 == 0.) x2 = AbsRandXY;

			//if((x1*vSegmExtNorm.x + y1*vSegmExtNorm.y) < 0.)
			//{
			//	IntLineCrossesFace = false; //to check whether this is true and whether it is necessary
			//}

			double kx = (y2 - y1)/(x2 - x1);
			double by = y1 - kx*x1;

			if(by == 0.) by = AbsRandXY; //?

			double kxe2 = kx*kx;
			double kxe2p1 = kxe2 + 1., kxe2m1 = kxe2 - 1.;
			double y1dx1 = y1/x1, y2dx2 = y2/x2;
			double invkxe2p1 = 1./kxe2p1;
			double invkxe2p1e2 = invkxe2p1*invkxe2p1;
			double x1e2 = x1*x1, x2e2 = x2*x2, y1e2 = y1*y1, y2e2 = y2*y2;
			double x1e2py1e2 = x1e2 + y1e2, x2e2py2e2 = x2e2 + y2e2, x2e2mx1e2 = x2e2 - x1e2, x2mx1 = x2 - x1;
			double x1e2py1e2dx1e2 = x1e2py1e2/x1e2, x2e2py2e2dx2e2 = x2e2py2e2/x2e2, x2e2py2e2dx1e2py1e2 = x2e2py2e2/x1e2py1e2;
			double kxy1px1 = kx*y1 + x1, kxy2px2 = kx*y2 + x2;
			double kxy1px1dby = kxy1px1/by, kxy2px2dby = kxy2px2/by;
			double axx1p2cz = ax*x1 + twocz, axx2p2cz = ax*x2 + twocz;
			double czkxe2p1 = cz*kxe2p1;

			piMult1 = 0;
			double sumAtan1 = atan(TransAtans(kxy2px2dby, -kxy1px1dby, piMult1)); //-ArcTan(kxy1px1dby) + ArcTan(kxy2px2dby)
			sumAtan1 += Pi*piMult1;
			double atan_y1dx1 = atan(y1dx1), atan_y2dx2 = atan(y2dx2);
			double log_x2e2py2e2dx2e2 = log(x2e2py2e2dx2e2), log_x1e2py1e2dx1e2 = log(x1e2py1e2dx1e2), log_x2e2py2e2dx1e2py1e2 = log(x2e2py2e2dx1e2py1e2);
			//consider storing logs and atans to minimize their re-calculation ?

			if(pJ_LinCoef == 0) //without linear term
			{
				IBxLoc = quartJz0*(2*ay*kx*x2e2mx1e2 + 2*by*invkxe2p1*(ay + ax*kx + 2*ay*kxe2)*x2mx1 - 2*by*invkxe2p1e2*(ay*by*kxe2m1 - 2*(czkxe2p1 - ax*by*kx))*sumAtan1 
					+ 2*ay*(x1e2*atan_y1dx1 - x2e2*atan_y2dx2) + axx2p2cz*x2*log_x2e2py2e2dx2e2 - axx1p2cz*x1*log_x1e2py1e2dx1e2 - by*invkxe2p1e2*(ax*by*kxe2m1 - 2*kx*(ay*by + czkxe2p1))*log_x2e2py2e2dx1e2py1e2);
				IByLoc = quartJz0*(-2*by*invkxe2p1*(ax + ay*kx)*x2mx1 - 2*by*invkxe2p1e2*(ax*by*kxe2m1 - 2*kx*(ay*by + czkxe2p1))*sumAtan1 
					+ 2*(axx1p2cz*x1*atan_y1dx1 - axx2p2cz*x2*atan_y2dx2) + ay*(x1e2*log_x1e2py1e2dx1e2 - x2e2*log_x2e2py2e2dx2e2) + by*invkxe2p1e2*(ay*by*kxe2m1 - 2*(czkxe2p1 - ax*by*kx))*log_x2e2py2e2dx1e2py1e2);
				IBzLoc = 0.25*(-2*ay*jx0*kx*x2e2mx1e2 - 2*by*invkxe2p1*(ax*(-jy0 + jx0*kx) + ay*(jx0 - jy0*kx + twojx0*kxe2))*x2mx1 
					+ 2*by*invkxe2p1e2*(ay*by*(-twojy0*kx + jx0*kxe2m1) + ax*by*(twojx0*kx + jy0*kxe2m1) - twocz*(jx0 + jy0*kx)*kxe2p1)*sumAtan1 
					+ 2*x2*(twocz*jy0 + ay*jx0*x2 + ax*jy0*x2)*atan_y2dx2 - 2*x1*(twocz*jy0 + ay*jx0*x1 + ax*jy0*x1)*atan_y1dx1 
					+ x1*(twocz*jx0 + ax*jx0*x1 - ay*jy0*x1)*log_x1e2py1e2dx1e2 - x2*(twocz*jx0 + ax*jx0*x2 - ay*jy0*x2)*log_x2e2py2e2dx2e2 
					- by*invkxe2p1e2*(ax*by*(jx0 + twojy0*kx - jx0*kxe2) + ay*by*(twojx0*kx + jy0*kxe2m1) - twocz*(jy0 - jx0*kx)*kxe2p1)*log_x2e2py2e2dx1e2py1e2);
			}
			else
			{
				double kxe2p1kx = kxe2p1*kx, kxe4 = kxe2*kxe2, invkxe2p1e3 = invkxe2p1e2*invkxe2p1, kxe2m3 = kxe2 - 3.;
				double kxe3 = kxe2*kx, kxe2p1e2 = kxe2p1*kxe2p1, bye2 = by*by;
				double kxe4m1 = kxe4 - 1, trekxe2m1 = 3.*kxe2 - 1, trekxe2p2 = 3.*kxe2 + 2.;
				double kxqxy = kx*qxy, kxe2qzy = kxe2*qzy;
				double x1e3 = x1e2*x1, x2e3 = x2e2*x2;
				double x2e3mx1e3 = x2e3 - x1e3;
				double czkxkxe2p1 = cz*kxe2p1kx, bykxqyy = by*kx*qyy, bykxqxy = by*kxqxy;
				double sixczkxkxe2p1 = 6.*czkxkxe2p1, triczkxkxe2p1 = 3.*czkxkxe2p1, triczkxe2p1 = 3.*czkxe2p1;

				IBxLoc = inv12*(invkxe2p1*(sixczkxkxe2p1*qzy + twoax*by*(kx*qzx + twoqzy + 3*kxe2qzy) + axe2*by*kx*qzz - aye2*by*kx*qzz 
					+ ay*(6*jz0*kxe2p1kx + sixczkxkxe2p1*qzz + by*(4*qzx + 6*kxe2*qzx - twoqzy*kx + 4*axqzz + 6*ax*kxe2*qzz)))*x2e2mx1e2 + 4*kx*(ay*qzxpaxqzz + ax*qzy)*x2e3mx1e3 
					+ (2*by*invkxe2p1e2*(triczkxe2p1*(kx*qzx + qzy + 2*kxe2qzy) - axe2*by*kxe2m1*qzz + aye2*by*kxe2m1*qzz + ax*(trijz0*kxe2p1kx + by*(twoqzx - twoqzx*kxe2 + 4*qzy*kx) + triczkxkxe2p1*qzz) 
					+ ay*(jz0*(3 + 9*kxe2 + 6*kxe4) + 2*by*(2*kx*qzxpaxqzz - qzy + kxe2qzy) + tricz*(1 + 3*kxe2 + 2*kxe4)*qzz)) + aytwoqzypayqzz*(3*bye2 + 3*by*kx*(x1 + x2) + kxe2*(x1e2 + x1*x2 + x2e2)))*x2mx1 
					+ 2*by*invkxe2p1e3*(6*czkxe2p1*(jz0*kxe2p1 - by*kx*qzxpaxqzz) + tricze2*kxe2p1e2*qzz + by*((twoax*by*kx*kxe2m3 - tricz*kxe4m1)*qzy 
					+ ay*(-trijz0*kxe4m1 + 2*by*kx*kxe2m3*qzxpaxqzz - triczqzz*kxe4m1)) + by*(axe2*by*qzz*trekxe2m1 - by*aytwoqzypayqzz*trekxe2m1 - twoax*(trijz0*kxe2p1kx - by*qzx*trekxe2m1)))*sumAtan1 
					+ 2*(qzy*(tricz + twoax*x1) + ay*(trijz0 + triczqzz + twoqzx*x1 + twoaxqzz*x1))*x1e2*atan_y1dx1 - 2*(qzy*(tricz + twoax*x2) + ay*(trijz0 + triczqzz + twoqzx*x2 + twoaxqzz*x2))*x2e2*atan_y2dx2 
					+ x2*(trecze2qzz + tricz*(twojz0 + qzxpaxqzz*x2) + x2*(axe2qzz*x2 - aytwoqzypayqzz*x2 + ax*(trijz0 + twoqzx*x2)))*log_x2e2py2e2dx2e2
					- x1*(trecze2qzz + tricz*(twojz0 + qzxpaxqzz*x1) + x1*(axe2qzz*x1 - aytwoqzypayqzz*x1 + ax*(trijz0 + twoqzx*x1)))*log_x1e2py1e2dx1e2 
					+ by*invkxe2p1e3*(axe2*bye2*kx*kxe2m3*qzz - aye2*bye2*kx*kxe2m3*qzz + triczkxe2p1*(twojz0*kxe2p1kx + by*(qzx - kxe2*qzx + twoqzy*kx) + czkxkxe2p1*qzz) 
					+ 2*ay*by*(trijz0*kxe2p1kx + by*(qzx - 3*kxe2*qzx + 3*qzy*kx - kxe3*qzy) + triczkxkxe2p1*qzz) 
					+ ax*by*(-trijz0*kxe4m1 - triczqzz*kxe4m1 + 2*by*(-3*kx*qzx + kxe3*qzx + qzypayqzz - 3*kxe2*qzypayqzz)))*log_x2e2py2e2dx1e2py1e2);
				IByLoc = inv12*(-(by*invkxe2p1*(axe2qzz + twoax*(qzx + kx*qzypayqzz) + ay*(2*kx*qzx + 4*qzy + 2*ayqzz + 3*kxe2*twoqzypayqzz))*x2e2mx1e2) - 2*kx*aytwoqzypayqzz*x2e3mx1e3 
					+ 2*by*invkxe2p1e2*(-triczkxe2p1*(qzx + kx*qzy) + 2*axe2*by*kx*qzz - 2*aye2*by*kx*qzz - ay*(trijz0*kxe2p1kx + by*(twoqzx - twoqzx*kxe2 + 4*qzy*kx) + triczkxkxe2p1*qzz) 
					- ax*(trijz0*kxe2p1 + triczkxe2p1*qzz - 2*by*(twoqzx*kx - qzy - ayqzz + kxe2*qzypayqzz)))*x2mx1 
					+ 2*by*invkxe2p1e3*(axe2*bye2*kx*kxe2m3*qzz - aye2*bye2*kx*kxe2m3*qzz + triczkxe2p1*(twojz0*kxe2p1kx + by*(qzx - kxe2*qzx + twoqzy*kx) + czkxkxe2p1*qzz) 
					+ 2*ay*by*(trijz0*kxe2p1kx + by*(qzx - 3*qzx*kxe2 + 3*qzy*kx - kxe3*qzy) + triczkxkxe2p1*qzz) 
					+ ax*by*(-trijz0*kxe4m1 - triczqzz*kxe4m1 + 2*by*(-3*qzx*kx + kxe3*qzx + qzypayqzz - 3*kxe2*qzypayqzz)))*sumAtan1 
					+ 2*x1*(trecze2qzz + tricz*(twojz0 + qzxpaxqzz*x1) + x1*(axe2qzz*x1 - aytwoqzypayqzz*x1 + ax*(trijz0 + twoqzx*x1)))*atan_y1dx1 
					- 2*x2*(trecze2qzz + tricz*(twojz0 + qzxpaxqzz*x2) + x2*(axe2qzz*x2 - aytwoqzypayqzz*x2 + ax*(trijz0 + twoqzx*x2)))*atan_y2dx2 
					+ (qzy*(tricz + twoax*x1) + ay*(trijz0 + triczqzz + twoqzx*x1 + twoaxqzz*x1))*x1e2*log_x1e2py1e2dx1e2 
					- (qzy*(tricz + twoax*x2) + ay*(trijz0 + triczqzz + twoqzx*x2 + twoaxqzz*x2))*x2e2*log_x2e2py2e2dx2e2
					+ by*invkxe2p1e3*(-triczkxe2p1*(twojz0*kxe2p1 + by*(-2*kx*qzxpaxqzz + qzy - kxe2qzy)) - tricze2*kxe2p1e2*qzz 
					+ ay*by*(trijz0*kxe4m1 - 2*by*(-3*kx*qzxpaxqzz + kxe3*qzxpaxqzz + qzy - 3*kxe2qzy) + triczqzz*kxe4m1) + aye2*bye2*qzz*trekxe2m1 
					+ ax*by*(6*jz0*kxe2p1kx + by*(6*qzy*kx - twoqzy*kxe3 - axqzz*trekxe2m1 - twoqzx*trekxe2m1)))*log_x2e2py2e2dx1e2py1e2);
				IBzLoc = inv12*(-(invkxe2p1*(sixczkxkxe2p1*qxy + twoax*by*(kx*qxxmqyy + 2*qxy + 3*kxe2*qxy - qyx) + axe2*by*(kx*qxz - qyz) + aye2*by*(kx*qxz - qyz)*trekxe2p2 
					+ ay*(6*jx0*kxe2p1kx + sixczkxkxe2p1*qxz + by*(6*kxe3*qxy + 4*ax*qxz + 6*kxe2*(ax*qxz - qyy) - 4*qyy - 2*kx*(-2*qxy + qyx + axqyz) + twoqxx*trekxe2p2)))*x2e2mx1e2) 
					- kx*(4*ax*qxy + 2*ay*(2*qxxmqyy + kxqxy + twoax*qxz) + aye2*(kx*qxz - 2*qyz))*x2e3mx1e3 
					- by*invkxe2p1e2*(aye2*by*(qxz + 8*kxe2*qxz + 3*kxe4*qxz - 4*kx*qyz) + 2*ay*(-trijy0*kxe2p1kx + jx0*(3 + 9*kxe2 + 6*kxe4) + 4*qxx*by*kx + by*qxy + 8*by*kxe2*qxy + 3*by*kxe4*qxy 
					+ triczqxz + 4*ax*by*kx*qxz + 9*czqxz*kxe2 + 6*czqxz*kxe4 - 2*by*qyx + 2*by*kxe2*qyx - 4*bykxqyy - twoax*by*qyz - triczqyz*kx + twoax*by*kxe2*qyz - triczqyz*kxe3) 
					+ 2*(triczkxe2p1*(kx*qxxmqyy + qxy + 2*kxe2*qxy - qyx) + axe2*by*(qxz - kxe2*qxz + 2*kx*qyz) 
					+ ax*(-trijy0*kxe2p1 + trijx0*kxe2p1kx + twoqxx*by - twoqxx*by*kxe2 + 4*bykxqxy + triczqxz*kx + triczqxz*kxe3 + 4*by*kx*qyx - 2*by*qyy + 2*by*kxe2*qyy - triczqyz - triczqyz*kxe2)))*x2mx1 
					- 2*by*invkxe2p1e3*(ax*by*(-6*jx0*kxe2p1kx - trijy0*kxe4m1 + by*(-3*axqyzp2qyxp2qxy*kx + axqyzp2qyxp2qxy*kxe3 + (-2 + 6*kxe2)*qxx - ax*qxz + 3*kxe2*(ax*qxz - 2*qyy) + 2*qyy)) 
					+ tricze2*kxe2p1e2*(qxz + kx*qyz) + aye2*bye2*(qxz - 3*kxe2*qxz - kx*kxe2m3*qyz) + ay*by*(6*jy0*kxe2p1kx - trijx0*kxe4m1 - 6*qxx*by*kx + twoqxx*by*kxe3 + 2*by*qxy - 6*by*kxe2*qxy 
					+ triczqxz - 6*ax*by*kx*qxz + twoax*by*kxe3*qxz - triczqxz*kxe4 + 2*by*qyx - 6*by*kxe2*qyx + 6*bykxqyy - 2*by*kxe3*qyy + twoax*by*qyz + 6*czqyz*kx - 6*ax*by*kxe2*qyz + 6*czqyz*kxe3) 
					+ triczkxe2p1*(twojx0*kxe2p1 + twojy0*kxe2p1kx + by*(-2*kx*qxxpaxqxzmqyy + qxypqyx + axqyz - kxe2*(qxypqyx + axqyz))))*sumAtan1  
					+ 2*x1*(-tricze2qyz - tricz*(twojy0 + qxypayqxzpqyxpaxqyz*x1) + x1*(aye2qyz*x1 - ax*(trijy0 + axqyzp2qyxp2qxy*x1) - ay*(trijx0 + twoqxxpaxqxzmqyy*x1)))*atan_y1dx1 
					- 2*x2*(-tricze2qyz - tricz*(twojy0 + qxypayqxzpqyxpaxqyz*x2) + x2*(aye2qyz*x2 - ax*(trijy0 + axqyzp2qyxp2qxy*x2) - ay*(trijx0 + twoqxxpaxqxzmqyy*x2)))*atan_y2dx2 
					+ x1*(trecze2qxz + tricz*(twojx0 + qxxpaxqxzmqyymayqyz*x1) + x1*(axe2qxz*x1 - ay*(trijy0 + twoqxypayqxzptwoqyx*x1) + ax*(trijx0 + twoqxxmqyymayqyz*x1)))*log_x1e2py1e2dx1e2 
					- x2*(trecze2qxz + tricz*(twojx0 + qxxpaxqxzmqyymayqyz*x2) + x2*(axe2qxz*x2 - ay*(trijy0 + twoqxypayqxzptwoqyx*x2) + ax*(trijx0 + twoqxxmqyymayqyz*x2)))*log_x2e2py2e2dx2e2
					- by*invkxe2p1e3*(axe2*bye2*(-3*kx*qxz + kxe3*qxz + qyz - 3*kxe2*qyz) - aye2*bye2*(-3*kx*qxz + kxe3*qxz + qyz - 3*kxe2*qyz) 
					+ triczkxe2p1*(-twojy0*kxe2p1 + twojx0*kxe2p1kx + by*qxx - by*kxe2*qxx + 2*bykxqxy + czqxz*kx + czqxz*kxe3 + 2*by*kx*qyx - by*qyy + by*kxe2*qyy - czqyz - czqyz*kxe2) 
					+ ax*by*(6*jy0*kxe2p1kx - trijx0*kxe4m1 - 6*qxx*by*kx + twoqxx*by*kxe3 + 2*by*qxy - 6*by*kxe2*qxy + 2*ay*by*qxz + triczqxz - 6*ay*by*kxe2*qxz - triczqxz*kxe4 + 2*by*qyx - 6*by*kxe2*qyx 
					+ 6*bykxqyy - 2*by*kxe3*qyy + 6*ay*by*kx*qyz + 6*czqyz*kx - 2*ay*by*kxe3*qyz + 6*czqyz*kxe3) + ay*by*(6*jx0*kxe2p1kx + trijy0*kxe4m1 + twoqxx*by - 6*qxx*by*kxe2 
					+ 6*bykxqxy - 2*by*kxe3*qxy + 6*czqxz*kx + 6*czqxz*kxe3 + 6*by*kx*qyx - 2*by*kxe3*qyx - 2*by*qyy + 6*by*kxe2*qyy - triczqyz + triczqyz*kxe4))*log_x2e2py2e2dx1e2py1e2);
			}
			vFaceContribIB.x += signY_EdgeNorm*IBxLoc;
			vFaceContribIB.y += signY_EdgeNorm*IByLoc;
			vFaceContribIB.z += signY_EdgeNorm*IBzLoc;
		}

		vResLocIB += signZ_FaceNorm*vFaceContribIB;
	}
	const double ConstForJ = 0.0002;
	vResLocIB *= ConstForJ;
	TVector3d vResIB = trIntAxis2Ez.TrVectField_inv(vResLocIB);
	
	if(pField->FieldKey.Ib_) pField->Ib += vResIB;
	if(pField->FieldKey.Ih_) pField->Ih += vResIB;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTPolyhedron::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTPolyhedronGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

void radTPolyhedron::Dump(std::ostream& o, int ShortSign) // Porting
{
	radTg3dRelax::Dump(o);
	DumpPureObjInfo(o, ShortSign);
	if(ShortSign==1) return;

	DumpMaterApplied(o);
	DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTPolyhedron::DumpPureObjInfo(std::ostream& o, int ShortSign)
{
	o << "Relaxable: ";
	o << "Polyhedron";

	if(ShortSign==1) return;

	o << endl;
	o << "   Number of faces: " << AmOfFaces << endl;
	o << "   {x,y,z}= {" << CentrPoint.x << ',' << CentrPoint.y << ',' << CentrPoint.z << '}' << endl;
	o << "   {mx,my,mz}= {" << Magn.x << ',' << Magn.y << ',' << Magn.z << '}';
	
	o << endl;
	o << "   Face Vertices:" << endl;

	int AmOfFaces_mi_1 = AmOfFaces;
	for(int i=0; i<AmOfFaces; i++)
	{
		o << "   {";

		radTHandlePgnAndTrans &FacePgnAndTrans = VectHandlePgnAndTrans[i];
		radTPolygon *FacePgnPtr = FacePgnAndTrans.PgnHndl.rep;
		double zcFace = FacePgnPtr->CoordZ;
		int numVertInCurFace = FacePgnPtr->AmOfEdgePoints;
		radTVect2dVect &vVertexPoints = FacePgnPtr->EdgePointsVector;
		radTrans *FaceTransPtr = FacePgnAndTrans.TransHndl.rep;

		TVector3d vP(0,0,zcFace);

		int numVertInCurFace_mi_1 = numVertInCurFace - 1;
		for(int j=0; j<numVertInCurFace; j++)
		{
			TVector2d &vP2d = vVertexPoints[j];
			vP.x = vP2d.x; vP.y = vP2d.y;
			TVector3d vPlab = FaceTransPtr->TrPoint(vP);
			o << '{' << vPlab.x << ',' << vPlab.y << ',' << vPlab.z << '}';
			if(j < numVertInCurFace_mi_1) o << ',';
		}

		o << '}';
		if(i < AmOfFaces_mi_1) o << "," << endl;
	}
	//o << '}';
}

//-------------------------------------------------------------------------

void radTPolyhedron::DumpBin_Polyhedron(CAuxBinStrVect& oStr)
{
	//TMatrix3d* pJ_LinCoef;
	if(pJ_LinCoef != 0)
	{
		oStr << (char)1;
		oStr << *pJ_LinCoef;
	}
	else oStr << (char)0;

	//char mLinTreat; //0- treat as relative
	oStr << mLinTreat;

	//int AmOfFaces;
	oStr << AmOfFaces;
	//radTVectHandlePgnAndTrans VectHandlePgnAndTrans;
	for(int i=0; i<AmOfFaces; i++)
	{
		radTHandlePgnAndTrans &hPgnAndTrans = VectHandlePgnAndTrans[i];
		
		radTHandle<radTPolygon> &hPgn = hPgnAndTrans.PgnHndl;
		if(hPgn.rep != 0)
		{
			oStr << (char)1;
			hPgn.rep->DumpBin_Polygon(oStr);

		}
		else oStr << (char)0;

		radTHandle<radTrans> &hTrans = hPgnAndTrans.TransHndl;
		if(hTrans.rep != 0)
		{
			oStr << (char)1;
			hTrans.rep->DumpBin_Trans(oStr);
		}
		else oStr << (char)0;

		//bool FaceIsInternalAfterCut;
		oStr << hPgnAndTrans.FaceIsInternalAfterCut;
	}

	//TVector3d J; //to move to base?
	oStr << J;

	//bool J_IsNotZero;
	oStr << J_IsNotZero;

	//short SomethingIsWrong;
	oStr << SomethingIsWrong;
	
	//radTPairOfDouble AuxPairOfDouble; // Used for cylindrical subdivision
	oStr << AuxPairOfDouble.First;
	oStr << AuxPairOfDouble.Second;
}

//-------------------------------------------------------------------------

void radTPolyhedron::DumpBinParse_Polyhedron(CAuxBinStrVect& inStr)
{
	//TMatrix3d* pJ_LinCoef;
	pJ_LinCoef = 0;
	char cJ_LinCoefDefined = 0;
	inStr >> cJ_LinCoefDefined;
	if(cJ_LinCoefDefined)
	{
		if(pJ_LinCoef != 0) delete pJ_LinCoef;
		pJ_LinCoef = new TMatrix3d();
		inStr >> (*pJ_LinCoef);
	}

	//char mLinTreat; //0- treat as relative
	inStr >> mLinTreat;

	//int AmOfFaces;
	inStr >> AmOfFaces;

	//radTVectHandlePgnAndTrans VectHandlePgnAndTrans;
	char cPgnDef=0, cTrfDef=0;
	for(int i=0; i<AmOfFaces; i++)
	{
		inStr >> cPgnDef;
		radTPolygon *pPgn=0;
		if(cPgnDef) pPgn = new radTPolygon(inStr);
		radTHandle<radTPolygon> hPgn(pPgn);

		inStr >> cTrfDef;
		radTrans *pTrf=0;
		if(cTrfDef) pTrf = new radTrans(inStr);
		radTHandle<radTrans> hTrf(pTrf);
		radTHandlePgnAndTrans hPgnAndTrans(hPgn, hTrf);
			
		//bool FaceIsInternalAfterCut;
		inStr >> hPgnAndTrans.FaceIsInternalAfterCut;

		VectHandlePgnAndTrans.push_back(hPgnAndTrans);
	}

	//TVector3d J; //to move to base?
	inStr >> J;

	//bool J_IsNotZero;
	inStr >> J_IsNotZero;

	//short SomethingIsWrong;
	inStr >> SomethingIsWrong;

	//radTPairOfDouble AuxPairOfDouble; // Used for cylindrical subdivision
	inStr >> AuxPairOfDouble.First;
	inStr >> AuxPairOfDouble.Second;
}

//-------------------------------------------------------------------------

void radTPolyhedron::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//Dumping objects that may be used by this object
	vector<pair<int, int> > vTrfKeys;
	DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys);

	int matKey=0;
	DumpBin_g3dRelax_TreatMat(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, matKey);

	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)Type_g3d();
	oStr << (char)Type_g3dRelax();
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTg3d
	DumpBin_g3d(oStr, vTrfKeys);

	//Members of radTg3dRelax
	DumpBin_g3dRelax(oStr, matKey);

	//Members of radTPolyhedron
	DumpBin_Polyhedron(oStr);
}

//-------------------------------------------------------------------------

void radTPolyhedron::DefineRelAndAbsTol(double* RelAbsTol)
{
	double RelZeroToler = 1.E-09;
	//RelZeroToler = 500.*max(RelZeroToler, radCR.RelRand);

	double MaxVal = RelZeroToler;
	//if(radCR.RelRand < radCR.RelRand) MaxVal = radCR.RelRand;
	if(MaxVal < radCR.RelRand) MaxVal = radCR.RelRand;
	//RelZeroToler = 500.*MaxVal;
	RelZeroToler = 100.*MaxVal; //OC291003

	RelAbsTol[0] = RelZeroToler;
	if(VectHandlePgnAndTrans.empty()) return;

	radTPolygon* PgnPtr = VectHandlePgnAndTrans[0].PgnHndl.rep;
	TVector2d& vpBuf = (PgnPtr->EdgePointsVector)[0];
	TVector3d aVertexPoint = TVector3d(vpBuf.x, vpBuf.y, PgnPtr->CoordZ);
	aVertexPoint = VectHandlePgnAndTrans[0].TransHndl.rep->TrPoint(aVertexPoint);
	TVector3d VectToCenter = CentrPoint - aVertexPoint;
	RelAbsTol[1] = RelZeroToler*NormAbs(VectToCenter);
	//RelAbsTol[1] = RelZeroToler*NormAbs(aVertexPoint); //OC090908
}

//-------------------------------------------------------------------------

int radTPolyhedron::CutItself(TVector3d* InCuttingPlane, radThg& In_hg, radTPair_int_hg& LowerNewPair_int_hg, radTPair_int_hg& UpperNewPair_int_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;
	char CutCoils = pSubdivOptions->SubdivideCoils;
	char AddNewElemsToGenCont = pSubdivOptions->PutNewStuffIntoGenCont;
	char SeparatePiecesAtCutting = pSubdivOptions->SeparatePiecesAtCutting;
	char MapInternalFacesAfterCut = pSubdivOptions->MapInternalFacesAfterCut;

	if(J_IsNotZero && (!CutCoils)) return 1; //OC290908

	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double& RelZeroToler = RelAbsTol[0];
	double& AbsZeroToler = RelAbsTol[1];
	TVector3d CuttingPlane[] = { *InCuttingPlane, InCuttingPlane[1] };
	//TVector3d CuttingPlane[] = { *InCuttingPlane - CentrPoint, InCuttingPlane[1] }; //OC090908

	radTSend Send;
	radTVectHandlePgnAndTrans LowerNewVectPgnAndTrans, UpperNewVectPgnAndTrans;
	radTVectOfPtrToVect3d VectOfTwoPoints3d;
	radTHandlePgnAndTrans LowerNewHandlePgnAndTrans, UpperNewHandlePgnAndTrans;

	TVector3d& PointOnCutPlane = *CuttingPlane;
	TVector3d& CutPlaneNormal = CuttingPlane[1];
	double SqLen = CutPlaneNormal.x*CutPlaneNormal.x + CutPlaneNormal.y*CutPlaneNormal.y + CutPlaneNormal.z*CutPlaneNormal.z;
	if(fabs(SqLen - 1.) > RelZeroToler) CutPlaneNormal = (1./sqrt(SqLen))*CutPlaneNormal;

	radTrans ResTransf;
	short SomethingFound = 0;
	if((!SubdInLocFrame) || ConsiderOnlyWithTrans)
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	if((!SubdInLocFrame) || ConsiderOnlyWithTrans)
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	if(SomethingFound) 
	{
		CutPlaneNormal = ResTransf.TrBiPoint_inv(CutPlaneNormal);
		PointOnCutPlane = ResTransf.TrPoint_inv(PointOnCutPlane);
	}

	for(int i=0; i<AmOfFaces; i++)
	{
		char FaceCoinsidenceNoticed = 0;
		int FoundOK = FindIntersectionWithFace(i, CuttingPlane, VectOfTwoPoints3d, LowerNewVectPgnAndTrans, UpperNewVectPgnAndTrans, FaceCoinsidenceNoticed, RelAbsTol);
		if(!FoundOK) { SomethingIsWrong = 1; return 0;}

		if(SeparatePiecesAtCutting && (!MapInternalFacesAfterCut))
		{
			if(FaceCoinsidenceNoticed) VectHandlePgnAndTrans[i].FaceIsInternalAfterCut = false;
		}
	}

	char ThereIsActualCut = ((VectOfTwoPoints3d.size() > 2) && (!LowerNewVectPgnAndTrans.empty()) && (!UpperNewVectPgnAndTrans.empty()));
	if(ThereIsActualCut)
	{
		radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedPolyhedron(this);
		if(GroupInPlaceOfThisPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg NewHandle(GroupInPlaceOfThisPtr);

		if(!DetermineNewFaceAndTrans(VectOfTwoPoints3d, CutPlaneNormal, LowerNewHandlePgnAndTrans, UpperNewHandlePgnAndTrans, RelAbsTol))
		{ SomethingIsWrong = 1; return 0;}

		if(MapInternalFacesAfterCut) LowerNewHandlePgnAndTrans.FaceIsInternalAfterCut = UpperNewHandlePgnAndTrans.FaceIsInternalAfterCut = true;
		else LowerNewHandlePgnAndTrans.FaceIsInternalAfterCut = UpperNewHandlePgnAndTrans.FaceIsInternalAfterCut = false;

		LowerNewVectPgnAndTrans.push_back(LowerNewHandlePgnAndTrans);
		UpperNewVectPgnAndTrans.push_back(UpperNewHandlePgnAndTrans);

		radThg hg1;
		if(!CreateNewEntity(LowerNewVectPgnAndTrans, hg1, radPtr->RecognizeRecMagsInPolyhedrons, RelAbsTol)) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
		int LowerNewElKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hg1) : 1;
		GroupInPlaceOfThisPtr->AddElement(LowerNewElKey, hg1);

		radThg hg2;
		if(!CreateNewEntity(UpperNewVectPgnAndTrans, hg2, radPtr->RecognizeRecMagsInPolyhedrons, RelAbsTol)) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
		int UpperNewElKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hg2) : 2;
		GroupInPlaceOfThisPtr->AddElement(UpperNewElKey, hg2);

		In_hg = NewHandle;

		if(SeparatePiecesAtCutting)
		{
			if(!g3dListOfTransform.empty())
			{
				radThg hgNewDpl = hg1;
				if(!(hg1.rep)->DuplicateItself(hgNewDpl, radPtr, AddNewElemsToGenCont)) return 0;
				((radTg3d*)(hgNewDpl.rep))->g3dListOfTransform = g3dListOfTransform;
				LowerNewElKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hgNewDpl) : 1;
				hg1 = hgNewDpl;

				hgNewDpl = hg2;
				if(!(hg2.rep)->DuplicateItself(hgNewDpl, radPtr, AddNewElemsToGenCont)) return 0;
				((radTg3d*)(hgNewDpl.rep))->g3dListOfTransform = g3dListOfTransform;
				UpperNewElKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hgNewDpl) : 2;
				hg2 = hgNewDpl;
			}
			LowerNewPair_int_hg.m = LowerNewElKey;
			LowerNewPair_int_hg.Handler_g = hg1;
			UpperNewPair_int_hg.m = UpperNewElKey;
			UpperNewPair_int_hg.Handler_g = hg2;

			if(AddNewElemsToGenCont)
			{
				int NewGroupElemKey = radPtr->RetrieveElemKey(GroupInPlaceOfThisPtr);
				radPtr->CopyDrawAttr(NewGroupElemKey, LowerNewElKey);
				radPtr->CopyDrawAttr(NewGroupElemKey, UpperNewElKey);
			}
		}
	}
	else
	{
		if(SeparatePiecesAtCutting)
		{
			char CenPtPositionChar;
			CheckCenPtPositionWithRespectToPlane(CuttingPlane, CenPtPositionChar);

			if(CenPtPositionChar == 'L')
			{
				LowerNewPair_int_hg.m = radPtr->RetrieveElemKey(this);
				LowerNewPair_int_hg.Handler_g = In_hg;
			}
			else
			{
				UpperNewPair_int_hg.m = radPtr->RetrieveElemKey(this);
				UpperNewPair_int_hg.Handler_g = In_hg;
			}
		}
	}

	for(int k=0; k<(int)VectOfTwoPoints3d.size(); k++) delete[] VectOfTwoPoints3d[k];
	VectOfTwoPoints3d.erase(VectOfTwoPoints3d.begin(), VectOfTwoPoints3d.end());
	LowerNewVectPgnAndTrans.erase(LowerNewVectPgnAndTrans.begin(), LowerNewVectPgnAndTrans.end());
	UpperNewVectPgnAndTrans.erase(UpperNewVectPgnAndTrans.begin(), UpperNewVectPgnAndTrans.end());
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::FindIntersectionWithFace(int FaceNo, TVector3d* CuttingPlane, 
	radTVectOfPtrToVect3d& VectOfTwoPoints3d, radTVectHandlePgnAndTrans& LowerNewVectPgnAndTrans, 
	radTVectHandlePgnAndTrans& UpperNewVectPgnAndTrans, char& FaceCoinsidenceNoticed, double* RelAbsTol)
{
	radTHandlePgnAndTrans& FacePgnAndTrans = VectHandlePgnAndTrans[FaceNo];
	radTPolygon* FacePgnPtr = FacePgnAndTrans.PgnHndl.rep;
	radTrans* FaceTransPtr = FacePgnAndTrans.TransHndl.rep;

	FaceCoinsidenceNoticed = 0;

	const double RelZeroToler = RelAbsTol[0];
	const double AbsZeroToler = RelAbsTol[1];

	TVector3d PointOnCutPlane = FaceTransPtr->TrPoint_inv(CuttingPlane[0]);
	TVector3d CutPlaneNormal = FaceTransPtr->TrBiPoint_inv(CuttingPlane[1]);

	double BufZ = PointOnCutPlane.z - FacePgnPtr->CoordZ;

	double AbsNx = Abs(CutPlaneNormal.x), AbsNy = Abs(CutPlaneNormal.y);
	if(Abs(Abs(CutPlaneNormal.z) - 1.) < RelZeroToler)
	{
		if(Abs(BufZ) < AbsZeroToler)
		{
			TVector3d VectorToCenter = CentrPoint - CuttingPlane[0];
			VectorToCenter = FaceTransPtr->TrBiPoint_inv(VectorToCenter);
			if(VectorToCenter*CutPlaneNormal < 0.) LowerNewVectPgnAndTrans.push_back(FacePgnAndTrans);
			else UpperNewVectPgnAndTrans.push_back(FacePgnAndTrans);

			FaceCoinsidenceNoticed = 1;
		}
		else if(BufZ < 0.)
		{
			(Abs(CutPlaneNormal.z - 1.) < RelZeroToler)?
				UpperNewVectPgnAndTrans.push_back(FacePgnAndTrans) :
				LowerNewVectPgnAndTrans.push_back(FacePgnAndTrans);
		}
		else
		{
			(Abs(CutPlaneNormal.z - 1.) < RelZeroToler)?
				LowerNewVectPgnAndTrans.push_back(FacePgnAndTrans) :
				UpperNewVectPgnAndTrans.push_back(FacePgnAndTrans);
		}
		return 1;
	}

	TVector2d VectIntrsctLine;
	VectIntrsctLine.x = CutPlaneNormal.y;
	VectIntrsctLine.y = -CutPlaneNormal.x;

	TVector2d PointOnIntrsctLine;
	if(Abs(CutPlaneNormal.y) >= RelZeroToler)
	{
		PointOnIntrsctLine.x = PointOnCutPlane.x;
		PointOnIntrsctLine.y = BufZ*(CutPlaneNormal.z/CutPlaneNormal.y) + PointOnCutPlane.y;
	}
	else
	{
		PointOnIntrsctLine.x = BufZ*(CutPlaneNormal.z/CutPlaneNormal.x) + PointOnCutPlane.x;
		PointOnIntrsctLine.y = 0.;
	}

	TVector2d Points2dAndNormProj[3];
	int IntersectingBoundsNos[2];
	TVector2d& FirstIntersectPoint2d = Points2dAndNormProj[0];
	TVector2d& SecondIntersectPoint2d = Points2dAndNormProj[1];
	int AmOfEdgePo = FacePgnPtr->AmOfEdgePoints;
	TVector2d R0_Bound = FacePgnPtr->EdgePointsVector[0], R1_Bound;
	short OneEdgePointAlreadyTrapped = 0, OnePointAlreadyFound = 0, TwoGoodPointsFound = 0;

	char AtLeastOnePointIsWithinBounds = 0; //OC291003

	for(int i_Bound = 0; i_Bound < AmOfEdgePo; i_Bound++)
	{
		TVector2d IntrsctPo;
		TLinesIntrsctCase IntrsctCase;

		R1_Bound = FacePgnPtr->EdgePointsVector[NextCircularNumber(i_Bound, AmOfEdgePo)];
		IntrsctOfTwoLines(VectIntrsctLine, PointOnIntrsctLine, R0_Bound, R1_Bound, IntrsctPo, IntrsctCase, RelAbsTol);

		if(IntrsctCase == PointOnBoundEdge)
		{
			if(OneEdgePointAlreadyTrapped)
			{
				if((Abs(FirstIntersectPoint2d.x - IntrsctPo.x) > AbsZeroToler) || (Abs(FirstIntersectPoint2d.y - IntrsctPo.y) > AbsZeroToler))
				{
					SecondIntersectPoint2d = IntrsctPo; TwoGoodPointsFound = 1; 
					IntersectingBoundsNos[1] = ((Abs(R0_Bound.x - IntrsctPo.x) < AbsZeroToler) && (Abs(R0_Bound.y - IntrsctPo.y) < AbsZeroToler))? i_Bound : i_Bound + 1;
				}
			}
			else
			{
				if(OnePointAlreadyFound)
				{
					SecondIntersectPoint2d = IntrsctPo; TwoGoodPointsFound = 1; 
					IntersectingBoundsNos[1] = ((Abs(R0_Bound.x - IntrsctPo.x) < AbsZeroToler) && (Abs(R0_Bound.y - IntrsctPo.y) < AbsZeroToler))? i_Bound : i_Bound + 1; 
					break;
				}
				else
				{
					FirstIntersectPoint2d = IntrsctPo;
					IntersectingBoundsNos[0] = ((Abs(R0_Bound.x - IntrsctPo.x) < AbsZeroToler) && (Abs(R0_Bound.y - IntrsctPo.y) < AbsZeroToler))? i_Bound : i_Bound + 1;
					OneEdgePointAlreadyTrapped = 1; OnePointAlreadyFound = 1;
				}
			}
		}
		else if(IntrsctCase == PointWithinBound)
		{
			AtLeastOnePointIsWithinBounds = 1; //OC291003
			if(OnePointAlreadyFound)
			{
				SecondIntersectPoint2d = IntrsctPo; TwoGoodPointsFound = 1; 
				IntersectingBoundsNos[1] = i_Bound; break;
			}
			else
			{
				FirstIntersectPoint2d = IntrsctPo; OnePointAlreadyFound = 1;
				IntersectingBoundsNos[0] = i_Bound;
			}
		}
		else if(IntrsctCase == LineIsIntrsct) 
		{ 
			TwoGoodPointsFound = 0; 

			TVector3d FirstPo3d(R0_Bound.x, R0_Bound.y, FacePgnPtr->CoordZ);
			TVector3d SecondPo3d(R1_Bound.x, R1_Bound.y, FacePgnPtr->CoordZ);
			FirstPo3d = FaceTransPtr->TrPoint(FirstPo3d);
			SecondPo3d = FaceTransPtr->TrPoint(SecondPo3d);
			if(!CheckIfTwoPointAlreadyMapped(FirstPo3d, SecondPo3d, VectOfTwoPoints3d, RelAbsTol))
			{
				TVector3d* TwoPoints3d = new TVector3d[2];
				TwoPoints3d[0] = FirstPo3d;
				TwoPoints3d[1] = SecondPo3d;
				VectOfTwoPoints3d.push_back(TwoPoints3d);
			}
			break;
		}
		R0_Bound = R1_Bound;
	}

	if((!AtLeastOnePointIsWithinBounds) && (TwoGoodPointsFound)) //OC291003
	{
		if(CheckIfOnlyNeighbouringEdgePointsTrapped(IntersectingBoundsNos, AmOfEdgePo)) TwoGoodPointsFound = 0; //OC291003
	}

	if(TwoGoodPointsFound)
	{
		TVector2d& CutPlaneNormProj = Points2dAndNormProj[2];
		CutPlaneNormProj.x = CutPlaneNormal.x;
		CutPlaneNormProj.y = CutPlaneNormal.y;
		if(!SetUpUpperAndLowerPolygon(Points2dAndNormProj, IntersectingBoundsNos, FacePgnAndTrans, LowerNewVectPgnAndTrans, UpperNewVectPgnAndTrans, RelAbsTol))
		{ SomethingIsWrong = 1; return 0;}
		
		TVector3d FirstPo3d(FirstIntersectPoint2d.x, FirstIntersectPoint2d.y, FacePgnPtr->CoordZ);
		TVector3d SecondPo3d(SecondIntersectPoint2d.x, SecondIntersectPoint2d.y, FacePgnPtr->CoordZ);
		TVector3d* TwoPoints3d = new TVector3d[2];
		TwoPoints3d[0] = FaceTransPtr->TrPoint(FirstPo3d);
		TwoPoints3d[1] = FaceTransPtr->TrPoint(SecondPo3d);
		VectOfTwoPoints3d.push_back(TwoPoints3d);
	}
	else
	{
		TVector2d CenterOfPgn(0.,0.);
		for(int i = 0; i < AmOfEdgePo; i++)
		{
			CenterOfPgn = CenterOfPgn + FacePgnPtr->EdgePointsVector[i];
		}
		CenterOfPgn = (1./double(AmOfEdgePo))*CenterOfPgn;

		TVector2d VectorToCenter = CenterOfPgn - PointOnIntrsctLine;
		if((VectorToCenter.x*CutPlaneNormal.x + VectorToCenter.y*CutPlaneNormal.y) < 0.) 
			LowerNewVectPgnAndTrans.push_back(FacePgnAndTrans);
		else UpperNewVectPgnAndTrans.push_back(FacePgnAndTrans);
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTPolyhedron::IntrsctOfTwoLines(const TVector2d& V1, const TVector2d& R01, const TVector2d& R02, const TVector2d& R12, TVector2d& IntrsctPo, TLinesIntrsctCase& IntrsctCase, double* RelAbsTol)
{// This is used at subdivision
	//const double t_Toler = RelAbsTol[0];
	const double t_Toler = 100*RelAbsTol[0]; //OC291003 - "safety margin" - to avoid problems at cuting by planes, etc.

	double AbsV1x = fabs(V1.x), AbsV1y = fabs(V1.y);
	double V1Norm = ((AbsV1x > AbsV1y)? AbsV1x : AbsV1y);
	//double V1Norm = sqrt(V1.x*V1.x + V1.y*V1.y); //OC081008??
	double V_Toler = V1Norm*t_Toler;

	TVector2d V2 = R12 - R02;

	double AbsV2x = fabs(V2.x), AbsV2y = fabs(V2.y);
	double MaxR = (AbsV2x > AbsV2y)? AbsV2x : AbsV2y;
	//double MaxR = V2.Abs(); //OC081008??
	double D_Toler = MaxR*V_Toler;

	double V1yV2x = V1.y*V2.x;
	double V1xV2y = V1.x*V2.y;
	double D = V1xV2y - V1yV2x;

	if(fabs(D) > D_Toler)
	{
		IntrsctPo.x = -(-R01.y*V1.x*V2.x + R02.y*V1.x*V2.x + R01.x*V1yV2x - R02.x*V1xV2y)/D;
		IntrsctPo.y = -(R02.y*V1yV2x - R01.y*V1xV2y + R01.x*V1.y*V2.y - R02.x*V1.y*V2.y)/D;
		double t_Intrsct = (Abs(V2.x) > V_Toler)? (IntrsctPo.x - R02.x)/V2.x : (IntrsctPo.y - R02.y)/V2.y;
		//IntrsctCase = ((t_Intrsct > t_Toler) && (t_Intrsct + t_Toler < 1.))? PointWithinBound : (((Abs(t_Intrsct) < t_Toler) || (Abs(t_Intrsct-1.) < t_Toler))? PointOnBoundEdge : PointOutsideBound);

			//DEBUG test
			//TVector2d V2u = V2, V1u = V1;
			//V2u.Normalize(); V1u.Normalize();
			//double t_IntrsctTest = ((R02.y - R01.y)*V1u.x - (R02.x - R01.x)*V1u.y)/(-V1u.x*V2u.y + V1u.y*V2u.x);
			//END DEBUG test

		if((t_Intrsct > t_Toler) && (t_Intrsct + t_Toler < 1.))
		{
			IntrsctCase = PointWithinBound; return;
		}
		else
		{
			if(Abs(t_Intrsct) < t_Toler)
			{
				IntrsctCase = PointOnBoundEdge;
				IntrsctPo = R02;
				return;
			}
			else if(Abs(t_Intrsct-1.) < t_Toler)
			{
				IntrsctCase = PointOnBoundEdge;
				IntrsctPo = R12;
				return;
			}
			else 
			{
				IntrsctCase = PointOutsideBound;
				return;
			}
		}
	}
	else 
	{
		IntrsctPo = R02;

		double V3x = R01.x-R02.x, V3y = R01.y-R02.y;
		double AbsV3x = fabs(V3x), AbsV3y = fabs(V3y);
		double V3Norm = ((AbsV3x > AbsV3y)? AbsV3x : AbsV3y);
		//double V3Norm = sqrt(V3x*V3x + V3y*V3y); //OC081008??

		double AbsR01x = fabs(R01.x), AbsR01y = fabs(R01.y);
		double R01Norm = ((AbsR01x > AbsR01y)? AbsR01x : AbsR01y);
		double AbsR02x = fabs(R02.x), AbsR02y = fabs(R02.y);
		double R02Norm = ((AbsR02x > AbsR02y)? AbsR02x : AbsR02y);
		double MaxNormR01R02 = (R01Norm > R02Norm)? R01Norm : R02Norm;

		if(V3Norm < MaxNormR01R02*t_Toler)
		{//is this really important?
			IntrsctCase = LineIsIntrsct; return;
		}

		//double LineCoinsBufToler = V1Norm*t_Toler;
		double LineCoinsToler = V_Toler*V3Norm;
		double CompareVal = fabs(V1.x*V3y - V1.y*V3x);
		IntrsctCase = (CompareVal < LineCoinsToler)? LineIsIntrsct : Zero; // To check !!!
	}
}

//-------------------------------------------------------------------------

int radTPolyhedron::CheckIfTwoPointAlreadyMapped(
	TVector3d& pp0, TVector3d& pp1, 
	radTVectOfPtrToVect3d& VectOfTwoPoints3d, double* RelAbsTol)
{
	const double AbsZeroToler = RelAbsTol[1];
	short CaseFound = 0;
	for(radTVectOfPtrToVect3d::iterator iter = VectOfTwoPoints3d.begin(); iter != VectOfTwoPoints3d.end(); ++iter)
	{
		TVector3d* aPair = *iter;
		TVector3d& p0 = aPair[0];
		TVector3d& p1 = aPair[1];
		if((((Abs(p0.x - pp0.x) < AbsZeroToler) && (Abs(p0.y - pp0.y) < AbsZeroToler) && (Abs(p0.z - pp0.z) < AbsZeroToler)) &&
			((Abs(p1.x - pp1.x) < AbsZeroToler) && (Abs(p1.y - pp1.y) < AbsZeroToler) && (Abs(p1.z - pp1.z) < AbsZeroToler))) ||
			(((Abs(p0.x - pp1.x) < AbsZeroToler) && (Abs(p0.y - pp1.y) < AbsZeroToler) && (Abs(p0.z - pp1.z) < AbsZeroToler)) &&
			((Abs(p1.x - pp0.x) < AbsZeroToler) && (Abs(p1.y - pp0.y) < AbsZeroToler) && (Abs(p1.z - pp0.z) < AbsZeroToler))))
		{ CaseFound = 1; break;}
	}
	return CaseFound;
}

//-------------------------------------------------------------------------

int radTPolyhedron::SetUpUpperAndLowerPolygon(
	TVector2d* Points2dAndNormProj, int* IntersectingBoundsNos, 
	radTHandlePgnAndTrans& FacePgnAndTrans, 
	radTVectHandlePgnAndTrans& LowerNewVectPgnAndTrans, 
	radTVectHandlePgnAndTrans& UpperNewVectPgnAndTrans, double* RelAbsTol)
{
	const double RelZeroToler = RelAbsTol[0];
	const double AbsZeroToler = RelAbsTol[1];

	radTPolygon* FacePgnPtr = FacePgnAndTrans.PgnHndl.rep;
	radTrans* FaceTransPtr = FacePgnAndTrans.TransHndl.rep;

	radTVect2dVect EdgePoVect[2];
	int AmOfEdgePo = FacePgnPtr->AmOfEdgePoints;
	int PointToLookFor = 0;

	TVector2d R0_Bound = FacePgnPtr->EdgePointsVector[0], R1_Bound;
	for(int i_Bound = 0; i_Bound < AmOfEdgePo; i_Bound++)
	{
		R1_Bound = FacePgnPtr->EdgePointsVector[NextCircularNumber(i_Bound, AmOfEdgePo)];
		EdgePoVect[PointToLookFor].push_back(R0_Bound);
		if(i_Bound == IntersectingBoundsNos[PointToLookFor])
		{
			TVector2d& IntersectPo = Points2dAndNormProj[PointToLookFor];

			if((Abs(IntersectPo.x - R0_Bound.x) > AbsZeroToler) || (Abs(IntersectPo.y - R0_Bound.y) > AbsZeroToler))
			{
				EdgePoVect[PointToLookFor].push_back(Points2dAndNormProj[PointToLookFor]);
			}

			if(PointToLookFor == 0)
			{
				EdgePoVect[1].push_back(Points2dAndNormProj[0]);
				PointToLookFor = 1;
			}
			else
			{
				EdgePoVect[0].push_back(Points2dAndNormProj[1]);
				PointToLookFor = 0;
			}
		}
		R0_Bound = R1_Bound;
	}

	int LowerPgnNo = 0;
	TVector2d VectBoundPgn1 = (EdgePoVect[1])[1] - (EdgePoVect[1])[0];
	if(VectBoundPgn1*Points2dAndNormProj[2] < 0.) LowerPgnNo = 1;
	int UpperPgnNo = (LowerPgnNo==0)? 1 : 0;

	radTSend Send;
	radTHandlePgnAndTrans LowerNewHandlePgnAndTrans, UpperNewHandlePgnAndTrans;

	if(EdgePoVect[LowerPgnNo].size() > 2) //OC
	{
		radTPolygon* PolygonPtr = new radTPolygon(EdgePoVect[LowerPgnNo], FacePgnPtr->CoordZ, FacePgnPtr->Magn);
		if(PolygonPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
		radTHandle<radTPolygon> LowerPgnHndl(PolygonPtr);

		//DEBUG
		//radTVect2dVect TestVect = EdgePoVect[LowerPgnNo];
		//TVector2d *Arr = 0;
		//if(TestVect.size() > 0) 
		//{
		//	Arr = new TVector2d[TestVect.size()];
		//	for(int i=0; i<(int)(TestVect.size()); i++) Arr[i] = TestVect[i];
		//}
		LowerNewHandlePgnAndTrans.PgnHndl = LowerPgnHndl;
		LowerNewHandlePgnAndTrans.TransHndl = FacePgnAndTrans.TransHndl;
		LowerNewHandlePgnAndTrans.FaceIsInternalAfterCut = FacePgnAndTrans.FaceIsInternalAfterCut;
		LowerNewVectPgnAndTrans.push_back(LowerNewHandlePgnAndTrans);
	}

	if(EdgePoVect[UpperPgnNo].size() > 2) //OC
	{
		radTPolygon* PolygonPtr = new radTPolygon(EdgePoVect[UpperPgnNo], FacePgnPtr->CoordZ, FacePgnPtr->Magn);
		if(PolygonPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
		radTHandle<radTPolygon> UpperPgnHndl(PolygonPtr);

		UpperNewHandlePgnAndTrans.PgnHndl = UpperPgnHndl;
		UpperNewHandlePgnAndTrans.TransHndl = FacePgnAndTrans.TransHndl;
		UpperNewHandlePgnAndTrans.FaceIsInternalAfterCut = FacePgnAndTrans.FaceIsInternalAfterCut;
		UpperNewVectPgnAndTrans.push_back(UpperNewHandlePgnAndTrans);
	}

	//DEBUG
	//if(Arr != 0) delete[] Arr;

	return 1;
}


//-------------------------------------------------------------------------

int radTPolyhedron::DetermineNewFaceAndTrans(
	radTVectOfPtrToVect3d& InVectOfTwoPoints3d, TVector3d& CutPlaneNormal, 
	radTHandlePgnAndTrans& LowerNewHandlePgnAndTrans, 
	radTHandlePgnAndTrans& UpperNewHandlePgnAndTrans, double* RelAbsTol)
{//I.e. determine and set up the new face which belongs to the cutting plane
	const double RelZeroToler = RelAbsTol[0];
	const double AbsZeroToler = RelAbsTol[1];

	radTVectOfPtrToVect3d VectOfTwoPoints3d = InVectOfTwoPoints3d;
	int AmOfEdgePoints = (int)VectOfTwoPoints3d.size();

	radTVectOfPtrToVect3d::iterator VectOfPtrToVect3dIter = VectOfTwoPoints3d.begin();

	TVector3d* P0Ptr = *VectOfPtrToVect3dIter;
	TVector3d* P1Ptr = P0Ptr + 1;
	++VectOfPtrToVect3dIter;

	radTSend Send;
	TVector3d** NewFaceFirstList = new TVector3d*[AmOfEdgePoints];
	if(NewFaceFirstList == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
	TVector3d** NewFaceSecondList = new TVector3d*[AmOfEdgePoints];
	if(NewFaceSecondList == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
	NewFaceSecondList[0] = P0Ptr;
	NewFaceSecondList[1] = P1Ptr;

		//DEBUG
		//TVector3d &p00 = *(InVectOfTwoPoints3d[0]);
		//TVector3d &p01 = *(InVectOfTwoPoints3d[0] + 1);
		//TVector3d &p10 = *(InVectOfTwoPoints3d[1]);
		//TVector3d &p11 = *(InVectOfTwoPoints3d[1] + 1);
		//TVector3d &p20 = *(InVectOfTwoPoints3d[2]);
		//TVector3d &p21 = *(InVectOfTwoPoints3d[2] + 1);
		//TVector3d &p30 = *(InVectOfTwoPoints3d[3]);
		//TVector3d &p31 = *(InVectOfTwoPoints3d[3] + 1);
		//END DEBUG

	int FirstCount = 0, SecondCount = 2;
	TVector3d BufVect;
	for(int i=2; i<AmOfEdgePoints; i++)
	{
		VectOfPtrToVect3dIter = VectOfTwoPoints3d.begin(); //OC23032007
		++VectOfPtrToVect3dIter; //OC23032007
		//after each pass if the following loop, the size of VectOfTwoPoints3d may be reduced!

		bool neighbSegmFound = false; //OC091008

		for(radTVectOfPtrToVect3d::iterator iter = VectOfPtrToVect3dIter; iter != VectOfTwoPoints3d.end(); ++iter)
		{
			TVector3d* TestPoint0Ptr = *iter;
			TVector3d* TestPoint1Ptr = TestPoint0Ptr + 1;
			
			//TVector3d BufVect00 = *TestPoint0Ptr - *P0Ptr;
			//if(NormAbs(BufVect00) < AbsZeroToler) 
			BufVect = *TestPoint0Ptr - *P0Ptr;
			if(NormAbs(BufVect) < AbsZeroToler) 
			{
				NewFaceFirstList[FirstCount++] = P0Ptr = TestPoint1Ptr;
				neighbSegmFound = true;
				VectOfTwoPoints3d.erase(iter); break;
			}
			//TVector3d BufVect01 = *TestPoint0Ptr - *P1Ptr;
			//if(NormAbs(BufVect01) < AbsZeroToler)
			BufVect = *TestPoint0Ptr - *P1Ptr;
			if(NormAbs(BufVect) < AbsZeroToler)
			{
				NewFaceSecondList[SecondCount++] = P1Ptr = TestPoint1Ptr;
				neighbSegmFound = true;
				VectOfTwoPoints3d.erase(iter); break;
			}
			//TVector3d BufVect10 = *TestPoint1Ptr - *P0Ptr;
			//if(NormAbs(BufVect10) < AbsZeroToler)
			BufVect = *TestPoint1Ptr - *P0Ptr;
			if(NormAbs(BufVect) < AbsZeroToler)
			{
				NewFaceFirstList[FirstCount++] = P0Ptr = TestPoint0Ptr;
				neighbSegmFound = true;
				VectOfTwoPoints3d.erase(iter); break;
			}
			//TVector3d BufVect11 = *TestPoint1Ptr - *P1Ptr;
			//if(NormAbs(BufVect11) < AbsZeroToler)
			//if(NormAbs(*TestPoint1Ptr - *P1Ptr) < AbsZeroToler)
			BufVect = *TestPoint1Ptr - *P1Ptr;
			if(NormAbs(BufVect) < AbsZeroToler)
			{
				NewFaceSecondList[SecondCount++] = P1Ptr = TestPoint0Ptr;
				neighbSegmFound = true;
				VectOfTwoPoints3d.erase(iter); break;
			}
		}

		if((!neighbSegmFound) && (VectOfTwoPoints3d.size() > 1)) //OC091008
		{//processing rare case: no exactly "touching" segment(s) found; trying to find closest among remaining

			TVector3d auxR = (*VectOfPtrToVect3dIter)[0] - *P0Ptr;
			double curMinDistE2 = auxR.AmpE2();
			double testDistE2 = curMinDistE2;
			int choiceCase = 0;
			radTVectOfPtrToVect3d::iterator iterClosest = VectOfPtrToVect3dIter;
			for(radTVectOfPtrToVect3d::iterator iter = VectOfPtrToVect3dIter; iter != VectOfTwoPoints3d.end(); ++iter)
			{
				TVector3d* TestPoint0Ptr = *iter;
				TVector3d* TestPoint1Ptr = TestPoint0Ptr + 1;

				auxR = *TestPoint0Ptr - *P0Ptr; testDistE2 = auxR.AmpE2();
				if(curMinDistE2 > testDistE2) { curMinDistE2 = testDistE2; choiceCase = 0; iterClosest = iter;}

				auxR = *TestPoint0Ptr - *P1Ptr; testDistE2 = auxR.AmpE2();
				if(curMinDistE2 > testDistE2) { curMinDistE2 = testDistE2; choiceCase = 1; iterClosest = iter;}

				auxR = *TestPoint1Ptr - *P0Ptr; testDistE2 = auxR.AmpE2();
				if(curMinDistE2 > testDistE2) { curMinDistE2 = testDistE2; choiceCase = 2; iterClosest = iter;}

				auxR = *TestPoint1Ptr - *P1Ptr; testDistE2 = auxR.AmpE2();
				if(curMinDistE2 > testDistE2) { curMinDistE2 = testDistE2; choiceCase = 3; iterClosest = iter;}
			}
			if(choiceCase == 0)
			{
				NewFaceFirstList[FirstCount++] = P0Ptr = (*iterClosest) + 1; //TestPoint1Ptr;
			}
			else if(choiceCase == 1)
			{
				NewFaceSecondList[SecondCount++] = P1Ptr = (*iterClosest) + 1; //TestPoint1Ptr;
			}
			else if(choiceCase == 2)
			{
				NewFaceFirstList[FirstCount++] = P0Ptr = *iterClosest; //TestPoint0Ptr;
			}
			else //if(choiceCase == 3)
			{
				NewFaceSecondList[SecondCount++] = P1Ptr = *iterClosest; //TestPoint0Ptr;
			}
			VectOfTwoPoints3d.erase(iterClosest);
		}
	}
	TVector3d** ArrayOf3dEdgePo = new TVector3d*[AmOfEdgePoints];
	if(ArrayOf3dEdgePo == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}

	//Filling-in the array of actual edge points of the new polygon
	//TVector3d** TraversMain = &(ArrayOf3dEdgePo[FirstCount-1]);
	//TVector3d** TraversBuf = &(NewFaceFirstList[0]);
	//int j;
	//for(j=0; j<FirstCount; j++) *(TraversMain--) = *(TraversBuf++);
	int j;
	TVector3d **TraversMain, **TraversBuf; //OC081008: modified because of observed crash at (FirstCount == 0)
	if(FirstCount > 0)
	{
		TraversMain = &(ArrayOf3dEdgePo[FirstCount-1]);
		TraversBuf = &(NewFaceFirstList[0]);
		for(j=0; j<FirstCount; j++) *(TraversMain--) = *(TraversBuf++);
	}
	TraversMain = &(ArrayOf3dEdgePo[FirstCount]);
	TraversBuf = NewFaceSecondList;
	for(j=0; j<SecondCount; j++) *(TraversMain++) = *(TraversBuf++);

	//Determining "upper" and "lower" normal vectors to the new polygon
	TVector3d& p0 = *(ArrayOf3dEdgePo[0]);
	TVector3d& p1 = *(ArrayOf3dEdgePo[1]);
	TVector3d v1 = p1-p0, v2;
	double BufAbs1 = NormAbs(v1);
	double RelZeroToler_e2 = RelZeroToler*RelZeroToler;
	TVector3d v1_VectBy_v2;
	for(j=2; j<AmOfEdgePoints; j++)
	{
		TVector3d& TestP = *(ArrayOf3dEdgePo[j]);
		v2 = TestP - p1;
		double AbsTol_e2 = RelZeroToler_e2*BufAbs1*NormAbs(v2);
		v1_VectBy_v2 = TVector3d(v1.y*v2.z - v2.y*v1.z, v2.x*v1.z - v1.x*v2.z, v1.x*v2.y - v2.x*v1.y);
		if(NormAbs(v1_VectBy_v2) > AbsTol_e2) break;
	}
	TVector3d& NormalForLower = CutPlaneNormal;
	short ReorderingNeededForLower = (v1_VectBy_v2*CutPlaneNormal < 0.)? 1 : 0;
	TVector3d NormalForUpper = (-1.)*NormalForLower;

	int FilledOK = 0;
	if(ReorderingNeededForLower) FilledOK = FillInNewHandlePgnAndTransFrom3d(ArrayOf3dEdgePo, AmOfEdgePoints, NormalForUpper, UpperNewHandlePgnAndTrans, RelAbsTol);
	else FilledOK = FillInNewHandlePgnAndTransFrom3d(ArrayOf3dEdgePo, AmOfEdgePoints, NormalForLower, LowerNewHandlePgnAndTrans, RelAbsTol);
	if(!FilledOK) { SomethingIsWrong = 1; return 0;}

	ReverseArrayOfVect3dPtr(ArrayOf3dEdgePo, AmOfEdgePoints);
	if(ReorderingNeededForLower) FilledOK = FillInNewHandlePgnAndTransFrom3d(ArrayOf3dEdgePo, AmOfEdgePoints, NormalForLower, LowerNewHandlePgnAndTrans, RelAbsTol);
	else FilledOK = FillInNewHandlePgnAndTransFrom3d(ArrayOf3dEdgePo, AmOfEdgePoints, NormalForUpper, UpperNewHandlePgnAndTrans, RelAbsTol);
	if(!FilledOK) { SomethingIsWrong = 1; return 0;}

	delete[] ArrayOf3dEdgePo;
	delete[] NewFaceSecondList;
	delete[] NewFaceFirstList;
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::FillInNewHandlePgnAndTransFrom3d(
	TVector3d** ArrayOf3dEdgePo, int AmOfEdgePo, TVector3d& N, 
	radTHandlePgnAndTrans& aHandlePgnAndTrans, double* RelAbsTol)
{
	radTSend Send;
	const double RelLenToler = RelAbsTol[0];

	TVector3d St1, St2, St3;
	if(Abs(N.z+1.) > RelLenToler)
	{
		double InvNzp1 = 1./(N.z + 1.);
		St1 = TVector3d(N.y*N.y*InvNzp1 + N.z, -N.x*N.y*InvNzp1, -N.x);
		St2 = TVector3d(St1.y, N.x*N.x*InvNzp1 + N.z, -N.y);
		St3 = TVector3d(-St1.z, -St2.z, N.z);
	}
	else
	{
		St1 = TVector3d(1., 0., 0.);
		St2 = TVector3d(0., -1., 0.);
		St3 = TVector3d(0., 0., -1.);
	}
	TMatrix3d R(St1, St2, St3);
	TVector3d Zero(0.,0.,0.);
	radTrans* RotationPtr = new radTrans(R, Zero, 1., 1., 2);
	if(RotationPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}

	radTVect2dVect Vect2dVect;

	//TVector2d *TestArr = new TVector2d[AmOfEdgePo]; //DEBUG

	double LocCoordZ;
	//double AbsTol = RelAbsTol[1]; //OC170902
    //TVector2d FirstP2d, PrevP2d; //OC170902

	int AmOfEdgePo_mi_1 = AmOfEdgePo - 1;
	for(int k=0; k<AmOfEdgePo; k++)
	{
		TVector3d& P = *(ArrayOf3dEdgePo[k]);
		TVector3d P_loc = RotationPtr->TrPoint(P);
		TVector2d P2d(P_loc.x, P_loc.y);

		//if(k == 0) 
		//{
		//	FirstP2d = P2d; //OC170902
		//}
		//else 
		//{
		//	//if((Abs(P2d.x - PrevP2d.x) < AbsTol) && (Abs(P2d.y - PrevP2d.y) < AbsTol)) continue;
		//	//if(k == AmOfEdgePo_mi_1)  //OC281003
		//	//{
		//	//	if((Abs(P2d.x - FirstP2d.x) < AbsTol) && (Abs(P2d.y - FirstP2d.y) < AbsTol)) break;
		//	//}
		//}

		Vect2dVect.push_back(P2d);
		//PrevP2d = P2d; //OC281003

		LocCoordZ = P_loc.z;

		//TestArr[k] = P2d; //DEBUG
	}

	//int ActAmOfEdgePts = Vect2dVect.size(); //OC281003
	//if(ActAmOfEdgePts <= 2) return 1; //OC281003

	TVector3d LocMagn = RotationPtr->TrVectField(Magn);
	RotationPtr->Invert();
	radTHandle<radTrans> HandleRotat(RotationPtr);

	radTPolygon* FacePgnPtr = new radTPolygon(Vect2dVect, LocCoordZ, LocMagn);
	if(FacePgnPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
	radTHandle<radTPolygon> HandlePgn(FacePgnPtr);

	aHandlePgnAndTrans.PgnHndl = HandlePgn;
	aHandlePgnAndTrans.TransHndl = HandleRotat;

	//delete[] TestArr; //Debug

	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::KsFromSizeToNumb(double* SubdivArray, int AmOfDir, radTSubdivOptions* pSubdivOptions)
{// Attention: this does not modify frame, the frame is modified (if necessary) deeply in SubdivideItselfByParPlanes
	radTSend Send;
	char TransfShouldBeTreated = (((pSubdivOptions->SubdivisionFrame != 0) || ConsiderOnlyWithTrans) && (!g3dListOfTransform.empty()));
	radTrans ResTransf;
	short SomethingFound = 0;
	if(TransfShouldBeTreated && (pSubdivOptions->SubdivisionFrame != 0)) FindResTransfWithMultOne(ResTransf, SomethingFound);
	else if(TransfShouldBeTreated && ConsiderOnlyWithTrans) FindInnerTransfWithMultOne(ResTransf, SomethingFound);

	TransfShouldBeTreated = (TransfShouldBeTreated && SomethingFound);

	TVector3d* Directions = new TVector3d[AmOfDir];
	if(Directions == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
	TVector3d* tDirections = Directions;
	double* tSubdivArray = SubdivArray;
	for(int i=0; i<AmOfDir; i++)
	{
		tDirections->x = *(tSubdivArray++);
		tDirections->y = *(tSubdivArray++);
		tDirections->z = *(tSubdivArray++);
		if(TransfShouldBeTreated) *tDirections = ResTransf.TrBiPoint_inv(*tDirections);

		tDirections++;
		tSubdivArray += 2;
	}

	double* Sizes = new double[AmOfDir];
	if(Sizes == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}

	EstimateSize(Directions, Sizes, AmOfDir);

	tDirections = Directions;
	tSubdivArray = SubdivArray;
	double* tSizes = Sizes;
	for(int k=0; k<AmOfDir; k++)
	{
		tSubdivArray += 3;
		*tSubdivArray = (*tSubdivArray < *tSizes)? Round((*tSizes)/(*tSubdivArray)) : 1.;
		tSizes++;
		tSubdivArray +=2;
	}

	if(Directions != 0) delete[] Directions;
	if(Sizes != 0) delete[] Sizes;
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::SubdivideItselfByParPlanes(double* InSubdivArray, int AmOfDir, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	double* SubdivArray = InSubdivArray;
	double* SubdivArrayToDelete = 0;
	radTSend Send;
	if(pSubdivOptions->SubdivisionParamCode == 1) 
	{
		SubdivArrayToDelete = new double[AmOfDir*5];
		if(SubdivArrayToDelete == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
		SubdivArray = SubdivArrayToDelete;

		double* tSubdivArray = SubdivArray;
		double* tInSubdivArray = InSubdivArray;
		for(int i=0; i<AmOfDir*5; i++) *(tSubdivArray++) = *(tInSubdivArray++);

		if(!KsFromSizeToNumb(SubdivArray, AmOfDir, pSubdivOptions)) return 0;
	}

	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	short MainGroupNotCreatedYet = 1, OneDirectionAlreadyPassed = 0;
	radTGroup* GroupInPlaceOfThisPtr = 0;
	radThg GroupHandle;
	radThg AuxGroupHandle = In_hg;

	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);

	radTSubdivOptions LocSubdivOptions = *pSubdivOptions;
	LocSubdivOptions.SubdivisionFrame = 0; // Subd in Local frame for SubdivideItselfByOneSetOfParPlanes here !

	radTvhg DummyVectOfHgChanged;

	int LocElemCount = 0;
	for(int kk=0; kk<AmOfDir; kk++)
	{
		double SubdivisionParam[2];

		int kk_mu_5 = kk*5;
		int kk_mu_5_p_3 = kk_mu_5 + 3;
		TVector3d PlanesNormal = TVector3d(SubdivArray[kk_mu_5], SubdivArray[kk_mu_5+1], SubdivArray[kk_mu_5+2]);
		
		double SqLen = PlanesNormal.x*PlanesNormal.x + PlanesNormal.y*PlanesNormal.y + PlanesNormal.z*PlanesNormal.z;
		if(fabs(SqLen - 1.) > RelAbsTol[0]) PlanesNormal = (1./sqrt(SqLen))*PlanesNormal;

		int AmOfPieces = int(SubdivArray[kk_mu_5_p_3] + 1.E-10);
		SubdivisionParam[0] = SubdivArray[kk_mu_5_p_3];
		SubdivisionParam[1] = SubdivArray[kk_mu_5 + 4];

		if(AmOfPieces > 1)
		{
			int AmOfPieces_mi_1 = AmOfPieces - 1;
			TVector3d* PointsOnCuttingPlanes = new TVector3d[AmOfPieces_mi_1];
			if(PointsOnCuttingPlanes == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}

			// The following may change the PlanesNormal
			if(!DeterminePointsOnCuttingPlanes(PlanesNormal, SubdivisionParam, SubdInLocFrame, PointsOnCuttingPlanes))
			{ SomethingIsWrong = 1; return 0;}

			if(MainGroupNotCreatedYet)
			{
				GroupInPlaceOfThisPtr = new radTSubdividedPolyhedron(this);
				if(GroupInPlaceOfThisPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
				GroupHandle = radThg(GroupInPlaceOfThisPtr);

				MainGroupNotCreatedYet = 0;
			}
			if(OneDirectionAlreadyPassed)
			{
				radTGroup aGroup;
				radTmhg& AuxGroupMapOfHandlers = ((radTGroup*)(AuxGroupHandle.rep))->GroupMapOfHandlers;
				for(radTmhg::iterator AuxGroupIter = AuxGroupMapOfHandlers.begin(); AuxGroupIter != AuxGroupMapOfHandlers.end(); ++AuxGroupIter)
				{
					radThg AuxTmpGroupHandle = (*AuxGroupIter).second;
					radTg* Old_gPtr = AuxTmpGroupHandle.rep;

					radThg BufHandle = AuxTmpGroupHandle;
					((radTg3d*)(BufHandle.rep))->SubdivideItselfByOneSetOfParPlanes(PlanesNormal, PointsOnCuttingPlanes, AmOfPieces_mi_1, AuxTmpGroupHandle, radPtr, &LocSubdivOptions, &DummyVectOfHgChanged);

					char AuxTmpGroupIsTreatedOnlyWithTrans = ((radTg3d*)(AuxTmpGroupHandle.rep))->ConsiderOnlyWithTrans;

					radTlphg& AuxTmpGroup_g3dListOfTransform = ((radTg3d*)(AuxTmpGroupHandle.rep))->g3dListOfTransform; // New
					char AuxTmpGroupHas_aTrans = !(AuxTmpGroup_g3dListOfTransform.empty()); // New

					if(AuxTmpGroupHandle.rep != Old_gPtr)
					{
						radTmhg& CurrentNewGroupMapOfHandlers = ((radTGroup*)(AuxTmpGroupHandle.rep))->GroupMapOfHandlers;
						for(radTmhg::iterator Iter = CurrentNewGroupMapOfHandlers.begin(); Iter != CurrentNewGroupMapOfHandlers.end(); ++Iter)
						{
							// Transformations (if appeared during the subdivision process) are copied to elements 
							if(AuxTmpGroupHas_aTrans) 
							{
								radTg3d* g3dLocPtr = (radTg3d*)((*Iter).second.rep);

								g3dLocPtr->g3dListOfTransform = AuxTmpGroup_g3dListOfTransform; // New
								g3dLocPtr->ConsiderOnlyWithTrans = AuxTmpGroupIsTreatedOnlyWithTrans; // New
							}
						
							aGroup.AddElement(++LocElemCount, (*Iter).second);
						}
					}
					else aGroup.AddElement((*AuxGroupIter).first, AuxTmpGroupHandle);
				}
				AuxGroupMapOfHandlers.erase(AuxGroupMapOfHandlers.begin(), AuxGroupMapOfHandlers.end());
				AuxGroupMapOfHandlers = aGroup.GroupMapOfHandlers;
			}
			else
			{
				SubdivideItselfByOneSetOfParPlanes(PlanesNormal, PointsOnCuttingPlanes, AmOfPieces_mi_1, AuxGroupHandle, radPtr, &LocSubdivOptions, &DummyVectOfHgChanged);
				OneDirectionAlreadyPassed = 1;
				LocElemCount = (int)(((radTGroup*)(AuxGroupHandle.rep))->GroupMapOfHandlers.size());
			}
			delete[] PointsOnCuttingPlanes;
		}
	}
	if(GroupInPlaceOfThisPtr != 0) 
	{
		int NewStuffCounter = 0;
		radTmhg& AuxGroupMapOfHandlers = ((radTGroup*)(AuxGroupHandle.rep))->GroupMapOfHandlers;
		for(radTmhg::iterator AuxGroupIter = AuxGroupMapOfHandlers.begin(); AuxGroupIter != AuxGroupMapOfHandlers.end(); ++AuxGroupIter)
		{
			radThg& a_hg = (*AuxGroupIter).second;
			if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(a_hg), a_hg);
			else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, a_hg);
		}
		In_hg = GroupHandle;
	}

	if(SubdivArrayToDelete != 0) delete[] SubdivArrayToDelete;
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::SubdivideItselfByOneSetOfParPlanes(
	TVector3d& InPlanesNormal, TVector3d* InPointsOnCuttingPlanes, int AmOfPieces_mi_1, 
	radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions, radTvhg* pVectOfHgChanged)
{// This does not put new sub-elements to general container !
	TVector3d PlanesNormal, *PointsOnCuttingPlanes;
	if(!TransferSubdivisionStructToLocalFrame(InPlanesNormal, InPointsOnCuttingPlanes, AmOfPieces_mi_1, pSubdivOptions, PlanesNormal, PointsOnCuttingPlanes)) return 0;
	
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	radTSend Send;
	radTGroup* GroupInPlaceOfThisPtr = 0;
	radThg NewHandle;
	short MainGroupNotCreatedYet = 1;

	radTSubdivOptions SubdivOptionsForCut;
	SubdivOptionsForCut.SubdivisionFrame = 0;
	SubdivOptionsForCut.SubdivisionParamCode = 0; // Is not used by Cut
	SubdivOptionsForCut.SubdivideCoils = SubdivideCoils; //OC071008
	SubdivOptionsForCut.PutNewStuffIntoGenCont = 0;
	SubdivOptionsForCut.SeparatePiecesAtCutting = 0;
	SubdivOptionsForCut.MapInternalFacesAfterCut = 1;

	radTPair_int_hg DummyLowerNewPair_int_hg, DummyUpperNewPair_int_hg; // To change sometime...

	radThg hgTmp = In_hg;
	radThg hgTmpBuf;
	int ElemCount = 0;
	for(int k=0; k<AmOfPieces_mi_1; k++)
	{
		TVector3d CuttingPlane[] = { PointsOnCuttingPlanes[k], PlanesNormal };

		radTg* Old_g3dPtr = hgTmp.rep;
		radThg hgTmpLocBuf = hgTmp;

		((radTg3d*)(hgTmpLocBuf.rep))->CutItself(CuttingPlane, hgTmp, DummyLowerNewPair_int_hg, DummyUpperNewPair_int_hg, radPtr, &SubdivOptionsForCut);

		if(hgTmp.rep != Old_g3dPtr)
		{
			if(MainGroupNotCreatedYet)
			{
				GroupInPlaceOfThisPtr = new radTSubdividedPolyhedron(this);
				if(GroupInPlaceOfThisPtr == 0) { SomethingIsWrong = 1; Send.ErrorMessage("Radia::Error900"); return 0;}
				NewHandle = radThg(GroupInPlaceOfThisPtr);
				
				MainGroupNotCreatedYet = 0;
			}
			radTmhg& TmpGroupMapOfHandlers = ((radTGroup*)(hgTmp.rep))->GroupMapOfHandlers;
			radTmhg::iterator GroupMapIter = TmpGroupMapOfHandlers.begin();
			radThg HandleOfLowerPiece = (*GroupMapIter).second;

			GroupInPlaceOfThisPtr->AddElement(++ElemCount, HandleOfLowerPiece);
			((radTg3d*)(HandleOfLowerPiece.rep))->SetMessageChar(1);

			hgTmpBuf = (*(++GroupMapIter)).second;
			hgTmp = hgTmpBuf;
		}
	}
	if(GroupInPlaceOfThisPtr != 0)
	{
		GroupInPlaceOfThisPtr->AddElement(++ElemCount, hgTmp);
		((radTg3d*)(hgTmp.rep))->SetMessageChar(1);

		In_hg = NewHandle;

		if(MessageChar==0) pVectOfHgChanged->push_back(In_hg);
	}

	if(PointsOnCuttingPlanes != InPointsOnCuttingPlanes) delete[] PointsOnCuttingPlanes;
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::DeterminePointsOnCuttingPlanes(
	TVector3d& PlanesNormal, double* SubdivisionParam, short SubdInLocFrame, 
	TVector3d* PointsOnCuttingPlanes)
{
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double& RelZeroTol = RelAbsTol[0];

	TVector3d LowestVertexPoint, UppestVertexPoint;
	radTrans DummyTrans;
	char DummyTransWasSet, DummyIgnore;

	radTSubdivOptions SubdivOptions;
	SubdivOptions.SubdivisionFrame = SubdInLocFrame? 0 : 1;
	SubdivOptions.SubdivisionParamCode = 0;
	SubdivOptions.SubdivideCoils = 0;

	FindLowestAndUppestVertices(PlanesNormal, &SubdivOptions, LowestVertexPoint, UppestVertexPoint, DummyTrans, DummyTransWasSet, DummyIgnore);
	
	TVector3d V = UppestVertexPoint - LowestVertexPoint;

	double& kk = SubdivisionParam[0];
	double& qq = SubdivisionParam[1];
	double q0 = (fabs(kk-1.)>RelZeroTol)? pow(qq, 1./(kk-1.)) : qq;
	double Buf = qq*q0 - 1.;
	double a1 = (fabs(Buf) > RelZeroTol)? (q0 - 1.)/Buf : 1./kk;

	int AmOfPieces_mi_1 = int(kk + 1.E-10 - 1.);
	double dTau = a1;
	double Tau = dTau;

	for(int j=0; j<AmOfPieces_mi_1; j++)
	{
		PointsOnCuttingPlanes[j] = LowestVertexPoint + Tau*V;
		dTau *= q0; Tau += dTau;
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::FindLowestAndUppestVertices(TVector3d& PlanesNormal, radTSubdivOptions* pSubdivOptions, 
	TVector3d& LowestVertexPoint, TVector3d& UppestVertexPoint, radTrans& Trans, char& TransWasSet, char& Ignore)
{
	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;
	Ignore = 0;
	TransWasSet = 0;
	TVector3d& ActualPlanesNormal = PlanesNormal;

	radTrans ResTransf;
	short SomethingFound = 0;

	if(!SubdInLocFrame)
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	else if(ConsiderOnlyWithTrans)
	{
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
	}
	if(SomethingFound) 
	{
		ActualPlanesNormal = ResTransf.TrBiPoint_inv(ActualPlanesNormal);
		Trans = ResTransf;
		TransWasSet = 1;
	}

	short Starting = 1;
	TVector3d LowestPo, UppestPo;
	for(radTVectHandlePgnAndTrans::iterator FaceIter = VectHandlePgnAndTrans.begin(); FaceIter != VectHandlePgnAndTrans.end(); ++FaceIter)
	{
		radTPolygon* PgnPtr = (*FaceIter).PgnHndl.rep;
		radTrans* TransPtr = (*FaceIter).TransHndl.rep;
		double LocZ = PgnPtr->CoordZ;
		for(radTVect2dVect::iterator PointIter = PgnPtr->EdgePointsVector.begin(); PointIter != PgnPtr->EdgePointsVector.end(); ++PointIter)
		{
			TVector2d& p2d = *PointIter;
			TVector3d p3d = TVector3d(p2d.x, p2d.y, LocZ);
			p3d = TransPtr->TrPoint(p3d);
			//p3d = TransPtr->TrPoint(p3d) + CentrPoint; //OC090908 ?

			if(Starting) { LowestPo = p3d; UppestPo = p3d; Starting = 0;}
			else
			{
				TVector3d TestLoV = p3d - LowestPo;
				TVector3d TestUpV = p3d - UppestPo;
				if(TestLoV*ActualPlanesNormal < 0.) LowestPo = p3d;
				if(TestUpV*ActualPlanesNormal > 0.) UppestPo = p3d;
			}
		}
	}
	LowestVertexPoint = LowestPo;
	UppestVertexPoint = UppestPo;
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::CheckForSpecialShapes(radTVectHandlePgnAndTrans& VectHandlePgnAndTrans, radThg& In_hg, double* RelAbsTol)
{// Currently, this distinguishes only parallelepipeds
	int NumberOfFaces = (int)VectHandlePgnAndTrans.size();
	if(NumberOfFaces != 6) return 0;

	double CPx, CPy, CPz, DimX, DimY, DimZ;
	TVector3d aNormal(0.,0.,1.);
	TVector3d AllNormals[6];
	double xSt, xFi, ySt, yFi, zSt, zFi;
	double RelTol = RelAbsTol[0], AbsTol = RelAbsTol[1];
	TVector3d Zero(0.,0.,0.);

	short InternalFacesMap[] = { 0,0,0,0,0,0 };
	short OrientationIsNotGood = 0;
	for(int k=0; k<NumberOfFaces; k++)
	{
		radTHandlePgnAndTrans& HandlePgnAndTrans = VectHandlePgnAndTrans[k];
		radTPolygon* PgnPtr = HandlePgnAndTrans.PgnHndl.rep;
		radTrans* TransPtr = HandlePgnAndTrans.TransHndl.rep;
		if(PgnPtr->AmOfEdgePoints != 4) return 0;

		TVector3d CurrentNormal = TransPtr->TrBiPoint(aNormal);

		double AbsNx = fabs(CurrentNormal.x), AbsNy = fabs(CurrentNormal.y), AbsNz = fabs(CurrentNormal.z);

		if((AbsNy < RelTol) && (AbsNz < RelTol) && (fabs(CurrentNormal.x + 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			xSt = AnEdgePo3d.x; 
			InternalFacesMap[0] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNy < RelTol) && (AbsNz < RelTol) && (fabs(CurrentNormal.x - 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			xFi = AnEdgePo3d.x;
			InternalFacesMap[1] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNx < RelTol) && (AbsNz < RelTol) && (fabs(CurrentNormal.y + 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			ySt = AnEdgePo3d.y;
			InternalFacesMap[2] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNx < RelTol) && (AbsNz < RelTol) && (fabs(CurrentNormal.y - 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			yFi = AnEdgePo3d.y;
			InternalFacesMap[3] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNx < RelTol) && (AbsNy < RelTol) && (fabs(CurrentNormal.z + 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			zSt = AnEdgePo3d.z;
			InternalFacesMap[4] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else if((AbsNx < RelTol) && (AbsNy < RelTol) && (fabs(CurrentNormal.z - 1.) < RelTol))
		{
			TVector2d& AnEdgePo2d = PgnPtr->EdgePointsVector[0];
			TVector3d AnEdgePo3d(AnEdgePo2d.x, AnEdgePo2d.y, PgnPtr->CoordZ);
			AnEdgePo3d = TransPtr->TrPoint(AnEdgePo3d);
			zFi = AnEdgePo3d.z;
			InternalFacesMap[5] = HandlePgnAndTrans.FaceIsInternalAfterCut;
		}
		else OrientationIsNotGood = 1;
		AllNormals[k] = CurrentNormal;
	}

	radTrans TransForRecMag;
	if(!OrientationIsNotGood)
	{
		CPx = 0.5*(xSt+xFi); DimX = xFi-xSt;
		CPy = 0.5*(ySt+yFi); DimY = yFi-ySt;
		CPz = 0.5*(zSt+zFi); DimZ = zFi-zSt;
	}
	else
	{
		int NoOfFacePerpToFirst=-1;
		TVector3d* TestNormalPtr = AllNormals;
		for(int j=0; j<NumberOfFaces; j++)
		{
			TVector3d* LocPtr = TestNormalPtr;
			for(int jj=j+1; jj<NumberOfFaces; jj++)
			{
				double ScalProd = (*TestNormalPtr)*(*(++LocPtr));
				short NormalsArePerp = (fabs(ScalProd) < RelTol);
				if(j==0) if(NoOfFacePerpToFirst < 0) if(NormalsArePerp) NoOfFacePerpToFirst = jj;

				if(!(NormalsArePerp || (fabs(ScalProd+1.) < RelTol))) return 0;
			}
			TestNormalPtr++;
		}
		if(NoOfFacePerpToFirst < 0) return 0;
		radTHandlePgnAndTrans& FirstHandlePgnAndTrans = VectHandlePgnAndTrans[0];
		radTPolygon* FirstPgnPtr = FirstHandlePgnAndTrans.PgnHndl.rep;
		radTrans* FirstTransPtr = FirstHandlePgnAndTrans.TransHndl.rep;

		radTHandlePgnAndTrans& NextHandlePgnAndTrans = VectHandlePgnAndTrans[NoOfFacePerpToFirst];
		radTPolygon* NextPgnPtr = NextHandlePgnAndTrans.PgnHndl.rep;
		radTrans* NextTransPtr = NextHandlePgnAndTrans.TransHndl.rep;

		TVector2d& p0 = FirstPgnPtr->EdgePointsVector[0];
		TVector2d& p1 = FirstPgnPtr->EdgePointsVector[1];
		TVector2d& p2 = FirstPgnPtr->EdgePointsVector[2];
		TVector2d& p3 = FirstPgnPtr->EdgePointsVector[3];

		TVector2d v0 = p1-p0;
		double InvNorm = 1./sqrt(v0.x*v0.x + v0.y*v0.y);
		v0 = InvNorm*v0;
		TVector3d St0(v0.x, v0.y, 0.);
		TVector3d St1(-v0.y, v0.x, 0.);
		TVector3d St2(0., 0., 1.);
		TMatrix3d M_fr_v0_to_i(St0, St1, St2);
		radTrans R(M_fr_v0_to_i, Zero, 1., 1., 2); // From v0 to i

		double FirstZ = FirstPgnPtr->CoordZ;
		TVector3d P0(p0.x, p0.y, FirstZ), P1(p1.x, p1.y, FirstZ), P2(p2.x, p2.y, FirstZ), P3(p3.x, p3.y, FirstZ);
		P0 = R.TrPoint(P0); P1 = R.TrPoint(P1); P2 = R.TrPoint(P2); P3 = R.TrPoint(P3);

		CPx = 0.5*(P0.x + P1.x), CPy = 0.5*(P1.y + P2.y);
		DimX = fabs(P1.x - P0.x), DimY = fabs(P2.y - P1.y);

		R.Invert();
		TrProduct(FirstTransPtr, &R, TransForRecMag);

		TVector3d OneMorePoint;
		double NextZ = NextPgnPtr->CoordZ;
		for(int kk=0; kk<3; kk++)
		{
			TVector2d& p = NextPgnPtr->EdgePointsVector[kk];
			TVector3d Pt(p.x, p.y, NextZ);
			Pt = NextTransPtr->TrPoint(Pt); Pt = TransForRecMag.TrPoint_inv(Pt);

			if((fabs(Pt.z-P0.z)>AbsTol) && (fabs(Pt.z-P1.z)>AbsTol) && (fabs(Pt.z-P2.z)>AbsTol) && (fabs(Pt.z-P1.z)>AbsTol))
			{
				OneMorePoint = Pt; break;
			}
		}
		CPz = 0.5*(P0.z + OneMorePoint.z);
		DimZ = fabs(P1.z - OneMorePoint.z);

		for(int i=0; i<NumberOfFaces; i++)
		{
			TVector3d NormalInLocFrame = TransForRecMag.TrBiPoint_inv(AllNormals[i]);
			short FaceNoForRecMag;
			if(fabs(NormalInLocFrame.x + 1.) < RelTol) FaceNoForRecMag = 0;
			else if(fabs(NormalInLocFrame.x - 1.) < RelTol) FaceNoForRecMag = 1;
			else if(fabs(NormalInLocFrame.y + 1.) < RelTol) FaceNoForRecMag = 2;
			else if(fabs(NormalInLocFrame.y - 1.) < RelTol) FaceNoForRecMag = 3;
			else if(fabs(NormalInLocFrame.z + 1.) < RelTol) FaceNoForRecMag = 4;
			else if(fabs(NormalInLocFrame.z - 1.) < RelTol) FaceNoForRecMag = 5;
			else return 0;
			InternalFacesMap[FaceNoForRecMag] = VectHandlePgnAndTrans[i].FaceIsInternalAfterCut;
		}
	}

	TVector3d CPoiVect(CPx, CPy, CPz), DimsVect(DimX, DimY, DimZ);
	TVector3d MagnForRecMag = OrientationIsNotGood? TransForRecMag.TrVectField_inv(Magn) : Magn;

	TVector3d J_ForRecMag = Zero;
	short J_IsNotZeroLoc = 0;
	if(J_IsNotZero)
	{
		J_ForRecMag = OrientationIsNotGood? TransForRecMag.TrVectField_inv(J) : J; //??
		J_IsNotZeroLoc = 1;
	}

	//short J_IsNotZero = 0;
	//radTRecMag* NewRecMagPtr = new radTRecMag(CPoiVect, DimsVect, MagnForRecMag, Zero, MaterHandle, J_IsNotZero);
	radTRecMag* NewRecMagPtr = new radTRecMag(CPoiVect, DimsVect, MagnForRecMag, J_ForRecMag, MaterHandle, J_IsNotZeroLoc);
	if(NewRecMagPtr==0) return 0;
	NewRecMagPtr->J_IsNotZero = 0;
	NewRecMagPtr->SetFacesInternalAfterCut(InternalFacesMap);

	if(OrientationIsNotGood)
	{
		radThg hTrans(new radTrans(TransForRecMag));
		NewRecMagPtr->AddTransform(1, hTrans);
		NewRecMagPtr->ConsiderOnlyWithTrans = 1;
	}

	radThg hRecMag(NewRecMagPtr);
	In_hg = hRecMag;
	return 1;
}

//-------------------------------------------------------------------------

double radTPolyhedron::Volume()
{
	const double kMax = 1.E+10;
	double VolSum = 0.;

	radTPolygon* PgnPtr;
	radTrans* TransPtr;
	radTHandlePgnAndTrans* HandlePgnAndTransPtr;

	for(int i=0; i<AmOfFaces; i++)
	{
		HandlePgnAndTransPtr = &(VectHandlePgnAndTrans[i]);
		PgnPtr = HandlePgnAndTransPtr->PgnHndl.rep;
		TransPtr = HandlePgnAndTransPtr->TransHndl.rep;

		TVector2d p1 = PgnPtr->EdgePointsVector[0], p2 = PgnPtr->EdgePointsVector[1], p3 = PgnPtr->EdgePointsVector[PgnPtr->IndexOfGoodThirdPoint()];
		TVector3d P1 = TVector3d(p1.x, p1.y, PgnPtr->CoordZ); P1 = TransPtr->TrPoint(P1);
		TVector3d P2 = TVector3d(p2.x, p2.y, PgnPtr->CoordZ); P2 = TransPtr->TrPoint(P2);
		TVector3d P3 = TVector3d(p3.x, p3.y, PgnPtr->CoordZ); P3 = TransPtr->TrPoint(P3);

		double x1 = P1.x, x2 = P2.x, x3 = P3.x;
		double y1 = P1.y, y2 = P2.y, y3 = P3.y;
		double z1 = P1.z, z2 = P2.z, z3 = P3.z;

		double x2mx1 = x2 - x1, y2my1 = y2 - y1;

		double BufDenom = (x2-x3)*y1 + (x3-x1)*y2 - x2mx1*y3;
		double AbsBufDenom_kMax = fabs(BufDenom)*kMax;
		double aSigBuf = y3*(z1-z2)+y1*(z2-z3)+y2*(z3-z1), aSig;
		double bSigBuf = x3*(z2-z1)+x2*(z1-z3)+x1*(z3-z2), bSig;

		if((fabs(aSigBuf) < AbsBufDenom_kMax) && (fabs(bSigBuf) < AbsBufDenom_kMax))
		{
			aSig = aSigBuf/BufDenom; bSig = bSigBuf/BufDenom;
			double cSig = z1 - aSig*x1 - bSig*y1;

			double xk1 = x1, yk1 = y1, xk2 = x2, yk2 = y2;
			for(int k=1; k<=PgnPtr->AmOfEdgePoints; k++)
			{
				if(k>1)
				{
					p3 = PgnPtr->EdgePointsVector[(k==PgnPtr->AmOfEdgePoints)? 0 : k];
					P3 = TVector3d(p3.x, p3.y, PgnPtr->CoordZ); P3 = TransPtr->TrPoint(P3);

					xk2 = P3.x; yk2 = P3.y;
					x2mx1 = xk2 - xk1; y2my1 = yk2 - yk1;
				}
				if(fabs(y2my1) < fabs(x2mx1)*kMax)
				{
					double aSigK = y2my1/x2mx1, bSigK = yk1 - aSigK*xk1;
					double xk1pxk2 = xk1 + xk2;
					double xk1e2pxk1xk2pxk2e2 = xk1*xk1 + xk1*xk2 + xk2*xk2;

					VolSum += (x2mx1/6.)*(3.*bSigK*(2.*cSig + aSig*xk1pxk2) + aSigK*(3.*cSig*xk1pxk2 + 2.*aSig*xk1e2pxk1xk2pxk2e2) 
							+ bSig*(3.*bSigK*(bSigK + aSigK*xk1pxk2) + aSigK*aSigK*xk1e2pxk1xk2pxk2e2));
				}
				xk1 = xk2; yk1 = yk2;
			}
		}
	}
	return fabs(VolSum);
}

//-------------------------------------------------------------------------

void radTPolyhedron::EstimateSize(TVector3d* Directions, double* Sizes, int AmOfDir)
{
	double Min[50], Max[50];
	int j;
	for(j=0; j<AmOfDir; j++)
	{
		Min[j] = 1.E+23; Max[j] = -1.E+23;
	}

	radTPolygon* PgnPtr;
	radTrans* TransPtr;
	radTHandlePgnAndTrans* HandlePgnAndTransPtr;

	for(int i=0; i<AmOfFaces; i++)
	{
		HandlePgnAndTransPtr = &(VectHandlePgnAndTrans[i]);
		PgnPtr = HandlePgnAndTransPtr->PgnHndl.rep;
		TransPtr = HandlePgnAndTransPtr->TransHndl.rep;

		for(int k=0; k<PgnPtr->AmOfEdgePoints; k++)
		{
			TVector2d& p2d = PgnPtr->EdgePointsVector[k];
			TVector3d p3d(p2d.x, p2d.y, PgnPtr->CoordZ); p3d = TransPtr->TrPoint(p3d);

			for(j=0; j<AmOfDir; j++)
			{
				double& MinProd = Min[j];
				double& MaxProd = Max[j];
				TVector3d& CurrentDir = Directions[j];
				double Prod = p3d*CurrentDir;
				if(Prod < MinProd) MinProd = Prod;
				if(Prod > MaxProd) MaxProd = Prod;
			}
		}
	}
	for(j=0; j<AmOfDir; j++)
	{
		TVector3d& Dir = Directions[j];
		double InvLenDirVect = 1./sqrt(Dir.x*Dir.x + Dir.y*Dir.y + Dir.z*Dir.z);
		Sizes[j] = (Max[j] - Min[j])*InvLenDirVect;
	}
}

//-------------------------------------------------------------------------

int radTPolyhedron::SubdivideItselfByEllipticCylinder(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	short OldRecognizeRecMagsInPolyhedrons = radPtr->RecognizeRecMagsInPolyhedrons;
	radPtr->RecognizeRecMagsInPolyhedrons = 0;

	double kPhi = SubdivArray[2], kz = SubdivArray[4];
	double qPhi = SubdivArray[3], qz = SubdivArray[5];

	const double PI = 3.1415926535898;

	radTCylindricSubdivSpec LocSubdivSpec = *pSubdivSpec;

	radTrans ResTransf;
	short SomethingFound = 0;
	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty()))
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	else if(ConsiderOnlyWithTrans && (!g3dListOfTransform.empty()))
	{
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
	}
	if(SomethingFound)
	{
		LocSubdivSpec.PointOnCylAx = ResTransf.TrPoint_inv(pSubdivSpec->PointOnCylAx);
		LocSubdivSpec.CylAxVect = ResTransf.TrAxialVect_inv(pSubdivSpec->CylAxVect);
		LocSubdivSpec.PointOnEllAx = ResTransf.TrPoint_inv(pSubdivSpec->PointOnEllAx);
	}

	TVector3d EdgePointsOverPhiAndAxForCylSubd[6];
	double Limits[6];

	SetMessageChar(0); // This is used to map AxisCrossesVolume for "on-top"
	FindEdgePointsOverPhiAndAxForCylSubd(&LocSubdivSpec, EdgePointsOverPhiAndAxForCylSubd, Limits);
	char AxisCrossesVolume = MessageChar;

	radTSend Send;
	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		double SizePhi;
		if(!AxisCrossesVolume)
		{
			double aCenPt, PhiCenPt;
			FindEllipticCoordOfPoint(&LocSubdivSpec, CentrPoint, aCenPt, PhiCenPt);
			SizePhi = EstimateLengthAlongEllipse(aCenPt, LocSubdivSpec.EllAxRatio, Limits[0], Limits[1]);
		}
		else
		{
			TVector3d MaxOffsetVertPt = EdgePointsOverPhiAndAxForCylSubd[5];
			double APt, PhiPt;
			FindEllipticCoordOfPoint(&LocSubdivSpec, MaxOffsetVertPt, APt, PhiPt);
			double a = 0.5*APt;
			SizePhi = PI*a*(1.5*(LocSubdivSpec.EllAxRatio + 1.) - sqrt(LocSubdivSpec.EllAxRatio));
		}
		double SizeAx = Limits[3] - Limits[2];

		kPhi = (kPhi < SizePhi)? Round(SizePhi/kPhi) : 1.;
		kz = (kz < SizeAx)? Round(SizeAx/kz) : 1.;
	}
	if(AxisCrossesVolume && (int(kPhi)==1)) { Send.ErrorMessage("Radia::Error068"); return 0;}

	double kPhi_qPhi[2];
	kPhi_qPhi[0] = kPhi; kPhi_qPhi[1] = qPhi;

	radTSubdivOptions LocSubdivOptions = *pSubdivOptions;
	LocSubdivOptions.SubdivisionFrame = 0;
	LocSubdivOptions.SubdivisionParamCode = 0;
	LocSubdivOptions.SubdivideCoils = 0;
	LocSubdivOptions.PutNewStuffIntoGenCont = 0;
	LocSubdivOptions.ReplaceOldStuff = 0;
	LocSubdivOptions.SeparatePiecesAtCutting = 1;
	LocSubdivOptions.MapInternalFacesAfterCut = 1;

	radThg AuxGroupHandle = In_hg;
	if(!SubdivideItselfOverAzimuth(kPhi_qPhi, Limits, &LocSubdivSpec, AuxGroupHandle, radPtr, &LocSubdivOptions)) return 0;

	TVector3d EdgePointsOverEllipseSet[2];
	double EllipticCoordOfEdgePoints[2];

	int EdgePointsOverEllipseSetFound = 0;
	if(LocSubdivOptions.MethForRadialSegmAtEllCylSubd == 0)
		EdgePointsOverEllipseSetFound = FindEdgePointsOverEllipseSet0(SubdivArray, &LocSubdivSpec, AuxGroupHandle, EdgePointsOverEllipseSet, EllipticCoordOfEdgePoints, &LocSubdivOptions);
	else if(LocSubdivOptions.MethForRadialSegmAtEllCylSubd == 1)
		EdgePointsOverEllipseSetFound = FindEdgePointsOverEllipseSet(SubdivArray, &LocSubdivSpec, AuxGroupHandle, EdgePointsOverEllipseSet, EllipticCoordOfEdgePoints, &LocSubdivOptions);

	if(!EdgePointsOverEllipseSetFound) { Send.ErrorMessage("Radia::Error999"); return 0;}

	double ka = SubdivArray[0];
	double qa = SubdivArray[1];
	if((int(ka)==1) && (int(kPhi)==1) && (int(kz)==1)) return 1;

	double ka_qa[2];
	ka_qa[0] = ka; ka_qa[1] = qa;
	if(!SubdivideByEllipses(ka_qa, EllipticCoordOfEdgePoints, &LocSubdivSpec, AuxGroupHandle, radPtr, &LocSubdivOptions)) return 0;

	if(kz > 1.)
	{
		radThg OldAuxGroupHandle = AuxGroupHandle;
		radTvhg VectOfHgChanged;
		int AmOfPieces_mi_1 = int(kz) - 1;
		TVector3d* PointsOnCuttingPlanes = new TVector3d[AmOfPieces_mi_1];
		if(PointsOnCuttingPlanes==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		const double AbsZeroTol = 5.E-12;
		double DelZ = Limits[3] - Limits[2];
		double q0z = (fabs(kz-1.) > AbsZeroTol)? pow(qz, 1./(kz-1.)) : qz;
		double BufZ = qz*q0z - 1.;
		double NewDelZ = (fabs(BufZ) > AbsZeroTol)? DelZ*(q0z - 1.)/BufZ : DelZ/kz;

		TVector3d Pz = EdgePointsOverPhiAndAxForCylSubd[2];
		TVector3d* tPointsOnCuttingPlanes = PointsOnCuttingPlanes;
		for(int k=0; k<AmOfPieces_mi_1; k++)
		{
			Pz += NewDelZ*LocSubdivSpec.CylAxVect;
			NewDelZ *= q0z;
			*(tPointsOnCuttingPlanes++) = Pz;
		}

		if(!((radTg3d*)(OldAuxGroupHandle.rep))->SubdivideItselfByOneSetOfParPlanes(LocSubdivSpec.CylAxVect, PointsOnCuttingPlanes, AmOfPieces_mi_1, AuxGroupHandle, radPtr, &LocSubdivOptions, &VectOfHgChanged)) return 0;

		radTg3d* g3dPtr = (radTg3d*)(AuxGroupHandle.rep);
		//radTCast Cast;
		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup != 0)
		{
			char RespectKeys = 0;
			radTApplication* Loc_radPtr = 0; // No operations on Global Container !
			pGroup->FlattenNestedStructure(Loc_radPtr, RespectKeys);
		}
		if(PointsOnCuttingPlanes!=0) delete[] PointsOnCuttingPlanes;
	}

	if(AuxGroupHandle.rep != this) //OC030504
	{
		radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedPolyhedron(this);
		if(GroupInPlaceOfThisPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg NewHandle(GroupInPlaceOfThisPtr);

		int NewStuffCounter = 0;
		radTmhg& AuxGroupMapOfHandlers = ((radTGroup*)(AuxGroupHandle.rep))->GroupMapOfHandlers;
		for(radTmhg::iterator AuxGroupIter = AuxGroupMapOfHandlers.begin(); AuxGroupIter != AuxGroupMapOfHandlers.end(); ++AuxGroupIter)
		{
			radThg& a_hg = (*AuxGroupIter).second;
			if(pSubdivOptions->PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(a_hg), a_hg);
			else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, a_hg);
		}

		In_hg = NewHandle;
	}

	radPtr->RecognizeRecMagsInPolyhedrons = OldRecognizeRecMagsInPolyhedrons;
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::SubdivideByEllipses(double* ka_qa, double* aLimits, radTCylindricSubdivSpec* pCylSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	if(int(*ka_qa)==1) return 1;

	short OldRecognizeRecMagsInPolyhedrons = radPtr->RecognizeRecMagsInPolyhedrons;
	radPtr->RecognizeRecMagsInPolyhedrons = 0;

	radTSend Send;

	radTg3d* g3dPtr = (radTg3d*)(In_hg.rep);
	//radTCast Cast;
	radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
	if(pGroup != 0)
	{
		for(radTmhg::iterator iter = pGroup->GroupMapOfHandlers.begin(); iter != pGroup->GroupMapOfHandlers.end(); ++iter)
		{
			radThg& hgLoc = (*iter).second;
			radThg hgOld = hgLoc;
			radTPolyhedron* pPolyhdr = (radTPolyhedron*)((radTg3dRelax*)((radTg3d*)(hgLoc.rep)));
			if(!pPolyhdr->SubdivideByEllipses(ka_qa, aLimits, pCylSubdivSpec, hgLoc, radPtr, pSubdivOptions)) return 0;
		}

		char RespectKeys = 0;
		radTApplication* Loc_radPtr = 0; // No operations on Global Container !
		pGroup->FlattenNestedStructure(Loc_radPtr, RespectKeys);
		return 1;
	}

	radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast(g3dPtr);
	if(g3dRelaxPtr == 0) return 0;
	radTPolyhedron* pPolyhdr = radTCast::PolyhedronCast(g3dRelaxPtr);
	if(pPolyhdr == 0) return 0;

	TVector3d& CylAxVect = pCylSubdivSpec->CylAxVect;
	TVector3d& PointOnCylAx = pCylSubdivSpec->PointOnCylAx;
	TVector3d& PointOnEllAx = pCylSubdivSpec->PointOnEllAx;
	double Rat = pCylSubdivSpec->EllAxRatio;

	TVector3d PointOnEllAx_mi_PointOnCylAx = PointOnEllAx - PointOnCylAx;
	TVector3d Ex = PointOnEllAx_mi_PointOnCylAx - (PointOnEllAx_mi_PointOnCylAx*CylAxVect)*CylAxVect;
	double InvLenEx = 1./sqrt(Ex.x*Ex.x + Ex.y*Ex.y + Ex.z*Ex.z);
	Ex = InvLenEx*Ex;
	TVector3d Ey = CylAxVect^Ex;

	double ka = *ka_qa, qa = ka_qa[1];
	const double AbsZeroTol = 5.E-12;
	double DelA = aLimits[1] - aLimits[0];
	double q0a = (fabs(ka-1.) > AbsZeroTol)? pow(qa, 1./(ka-1.)) : qa;
	double BufA = qa*q0a - 1.;
	double a1a = (fabs(BufA) > AbsZeroTol)? DelA*(q0a - 1.)/BufA : DelA/ka;
  	double NewDelA = a1a;

	radTSubdivOptions SubdivOptionsForCut = *pSubdivOptions;
	SubdivOptionsForCut.SubdivisionFrame = 0;
	SubdivOptionsForCut.SubdivideCoils = 0;
	SubdivOptionsForCut.PutNewStuffIntoGenCont = 0;
	SubdivOptionsForCut.SeparatePiecesAtCutting = 1;
	SubdivOptionsForCut.MapInternalFacesAfterCut = 1;

	radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedPolyhedron(pPolyhdr);
	if(GroupInPlaceOfThisPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
	radThg NewHandle(GroupInPlaceOfThisPtr);

	radThg hgDupl = In_hg;

	TVector3d CuttingPlane[2];
	TVector3d& PointInCutPlane = *CuttingPlane;
	TVector3d& CutPlaneNormal = CuttingPlane[1];

	double a = aLimits[0];
	double aPrev = a;
	double xLoc, yLoc;

	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double SmallestPossiblePieceSize = 5.*RelAbsTol[1];

	radTPairOfDouble& PhiStFi = AuxPairOfDouble;

	double tSt, tFi, aDummy1, aDummy2;
	FindLocalEllipticCoord(cos(PhiStFi.First), sin(PhiStFi.First), Rat, aDummy1, tSt);
	FindLocalEllipticCoord(cos(PhiStFi.Second), sin(PhiStFi.Second), Rat, aDummy2, tFi);

	char ThereWereNoActualCuts = 1;
	int AmOfCuts = int(ka) - 1;
	int AmOfCuts_mi_1 = AmOfCuts - 1;
	int ElemCount = 0;

	int AmOfExtraCutAttempts = 100*AmOfCuts; // To steer manually !
	int TotAmOfCutAttempts = AmOfCuts + AmOfExtraCutAttempts;

	for(int ia=0; ia<TotAmOfCutAttempts; ia++)
	{
		aPrev = a;
		a += NewDelA;
		NewDelA *= q0a;

		double b = Rat*a;
		xLoc = a*cos(tSt); yLoc = b*sin(tSt);
		TVector3d P1 = (xLoc*Ex + yLoc*Ey) + PointOnCylAx;
		xLoc = a*cos(tFi); yLoc = b*sin(tFi);
		TVector3d P2 = (xLoc*Ex + yLoc*Ey) + PointOnCylAx;

		PointInCutPlane = P1;
		TVector3d P2_mi_P1 = P2 - P1;
		CutPlaneNormal = P2_mi_P1^CylAxVect;

		radThg hgDuplCutted = hgDupl;
		radTPair_int_hg LowerNewPair_int_hg, UpperNewPair_int_hg;
		((radTg3d*)(hgDupl.rep))->CutItself(CuttingPlane, hgDuplCutted, LowerNewPair_int_hg, UpperNewPair_int_hg, radPtr, &SubdivOptionsForCut);

		if(UpperNewPair_int_hg.Handler_g.rep != 0)
		{
			hgDupl = UpperNewPair_int_hg.Handler_g;

			if(LowerNewPair_int_hg.Handler_g.rep != 0)
			{
				((radTg3d*)(LowerNewPair_int_hg.Handler_g.rep))->EraseAllTransformations();

				int NewElKey = (pSubdivOptions->PutNewStuffIntoGenCont)? radPtr->AddElementToContainer(LowerNewPair_int_hg.Handler_g) : ++ElemCount;
				GroupInPlaceOfThisPtr->AddElement(NewElKey, LowerNewPair_int_hg.Handler_g);

				ThereWereNoActualCuts = 0;
			}

			if(NewDelA <= SmallestPossiblePieceSize) break;
		}
		else
		{
			hgDupl = LowerNewPair_int_hg.Handler_g; break;
		}
	}
	//if(hgDupl.rep != 0)
	if((hgDupl.rep != 0) && (!ThereWereNoActualCuts)) //OC030504
	{
		((radTg3d*)(hgDupl.rep))->EraseAllTransformations();

		int NewElKey = (pSubdivOptions->PutNewStuffIntoGenCont)? radPtr->AddElementToContainer(hgDupl) : ++ElemCount;
		GroupInPlaceOfThisPtr->AddElement(NewElKey, hgDupl);

		ThereWereNoActualCuts = 0;
	}

	SetMessageChar(0);
	if(!ThereWereNoActualCuts)
	{
		In_hg = NewHandle;
		GroupInPlaceOfThisPtr->MessageChar = 1;
	}

	radPtr->RecognizeRecMagsInPolyhedrons = OldRecognizeRecMagsInPolyhedrons;
	return 1;
}

//-------------------------------------------------------------------------

void radTPolyhedron::FindLocalEllipticCoord(double xLoc, double yLoc, double Ratio, double& a, double& t)
{
	const double HalfPI = 1.5707963267949;
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double TolForZero = 10.*RelAbsTol[1];

	if(fabs(xLoc) < TolForZero)
	{
		if(yLoc > TolForZero) { t = HalfPI; a = yLoc/Ratio; return;}
		else if(yLoc < -TolForZero) { t = 3.*HalfPI; a = -yLoc/Ratio; return;}
		else { a = 0.; t = 0.; return;}
	}
	else if(xLoc > TolForZero)
	{
		if(yLoc > TolForZero) t = atan(yLoc/(Ratio*xLoc));
		else if(yLoc < -TolForZero) t = 4.*HalfPI + atan(yLoc/(Ratio*xLoc));
		else { a = xLoc; t = 0.; return;}
	}
	else if(xLoc < -TolForZero)
	{
		t = 2.*HalfPI + atan(yLoc/(Ratio*xLoc));
	}
	a = xLoc/cos(t);
}

//-------------------------------------------------------------------------

int radTPolyhedron::SubdivideItselfOverAzimuth(double* kPhi_qPhi, double* Limits, radTCylindricSubdivSpec* pCylSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	char AxisCrossesVolume = MessageChar;

	radTSend Send;
	int kPh = int(*kPhi_qPhi);
	if(AxisCrossesVolume && (kPh==1)) { Send.ErrorMessage("Radia::Error068"); return 0;}

	const double PI = 3.1415926535898;
	const double TwoPI = 2.*PI;

	if((Limits[1] < Limits[0]) && (!AxisCrossesVolume)) 
	{
		Limits[1] += TwoPI;
		if(Limits[1] < Limits[0]) Limits[1] += TwoPI;
	}
	if(AxisCrossesVolume) 
	{
		Limits[0] = 0.;
		Limits[1] = TwoPI;
	}

	if(kPh==1)
	{
		radTPairOfDouble PhiMinMax(Limits[0], Limits[1]);

		AuxPairOfDouble = PhiMinMax;
		return 1;
	}

	short OldRecognizeRecMagsInPolyhedrons = radPtr->RecognizeRecMagsInPolyhedrons;
	radPtr->RecognizeRecMagsInPolyhedrons = 0;

	char CompareWithActualLimits = 0;
	double ActualLimits[6];	
	if(!AxisCrossesVolume)
	{
		ActualLimits[0] = AuxPairOfDouble.First;
		ActualLimits[1] = AuxPairOfDouble.Second;

		if(ActualLimits[0] > TwoPI) ActualLimits[0] -= TwoPI;
		if(ActualLimits[1] > TwoPI) ActualLimits[1] -= TwoPI;
		CompareWithActualLimits = 1;
	}
	char ThereWereNoActualCuts = 1;

	TVector3d CuttingPlane[2];
	TVector3d& PointOnCutPlane = *CuttingPlane;
	PointOnCutPlane = pCylSubdivSpec->PointOnCylAx;
	TVector3d& CutPlaneNormal = CuttingPlane[1];
	TVector3d RadialVectInCutPlane;

	TVector3d& CylAxVect = pCylSubdivSpec->CylAxVect;
	TVector3d& PointOnEllAx = pCylSubdivSpec->PointOnEllAx;

	TVector3d PointOnEllAx_mi_PointOnCylAx = PointOnEllAx - PointOnCutPlane;
	TVector3d Ex = PointOnEllAx_mi_PointOnCylAx - (PointOnEllAx_mi_PointOnCylAx*CylAxVect)*CylAxVect;
	double InvLenEx = 1./sqrt(Ex.x*Ex.x + Ex.y*Ex.y + Ex.z*Ex.z);
	Ex = InvLenEx*Ex;
	TVector3d Ey = CylAxVect^Ex;

	double kPhi = *kPhi_qPhi, qPhi = kPhi_qPhi[1];
	const double AbsZeroTol = 5.E-12;
	double DelPhi = Limits[1] - Limits[0];
	double q0Phi = (fabs(kPhi-1.) > AbsZeroTol)? pow(qPhi, 1./(kPhi-1.)) : qPhi;
	double BufPhi = qPhi*q0Phi - 1.;
	double a1Phi = (fabs(BufPhi) > AbsZeroTol)? DelPhi*(q0Phi - 1.)/BufPhi : DelPhi/kPhi;
  	double NewDelPhi = a1Phi;

	if(fabs(NewDelPhi) >= PI) { Send.ErrorMessage("Radia::Error068"); return 0;}

	double Phi = Limits[0];
	double PrevPhi = Phi;
	double xLoc, yLoc;

	double CheckPhi = Phi;
	for(;;)
	{
		if(CheckPhi >= TwoPI) CheckPhi -= TwoPI;
		else break;
	}
	char CuttingGoesFromMiddle = AngleIsBetween(CheckPhi, ActualLimits[0], ActualLimits[1]);
	double BoundingPhi1 =0., BoundingPhi2 =0.;

	radTSubdivOptions SubdivOptionsForCut = *pSubdivOptions;
	SubdivOptionsForCut.SubdivisionFrame = 0;
	SubdivOptionsForCut.SubdivideCoils = 0;
	SubdivOptionsForCut.PutNewStuffIntoGenCont = 0;
	SubdivOptionsForCut.SeparatePiecesAtCutting = 1;
	SubdivOptionsForCut.MapInternalFacesAfterCut = 1;

	radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedPolyhedron(this);
	if(GroupInPlaceOfThisPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
	radThg NewHandle(GroupInPlaceOfThisPtr);

	radThg hgDupl = In_hg;
	radTPair_int_hg LowerNewPair_int_hg, UpperNewPair_int_hg, FirstLowerPair_int_hg;

	char PreviousCutWasOK = 0;
	int NewElKey;
	int AmOfFurtherPhiCuts = int(kPhi) - 1;
	for(int iPhi=0; iPhi<AmOfFurtherPhiCuts; iPhi++)
	{
		PrevPhi = Phi;
		Phi += NewDelPhi;
		NewDelPhi *= q0Phi;

		if(fabs(NewDelPhi) >= PI) { Send.ErrorMessage("Radia::Error068"); return 0;}

		double CheckPhi = Phi;
		for(;;)
		{
			if(CheckPhi >= TwoPI) CheckPhi -= TwoPI;
			else break;
		}
		char TryToMakeCut = CompareWithActualLimits? AngleIsBetween(CheckPhi, ActualLimits[0], ActualLimits[1]) : 1;
		if(TryToMakeCut)
		{
			xLoc = cos(Phi); yLoc = sin(Phi);
			RadialVectInCutPlane = xLoc*Ex + yLoc*Ey;
			CutPlaneNormal = CylAxVect^RadialVectInCutPlane;

			radThg hgDuplCutted = hgDupl;
			((radTg3d*)(hgDupl.rep))->CutItself(CuttingPlane, hgDuplCutted, LowerNewPair_int_hg, UpperNewPair_int_hg, radPtr, &SubdivOptionsForCut);

			PreviousCutWasOK = (hgDuplCutted.rep != hgDupl.rep);
			if(PreviousCutWasOK) ThereWereNoActualCuts = 0;

			if(AxisCrossesVolume || (CuttingGoesFromMiddle && (iPhi == 0)))
			{
				xLoc = cos(PrevPhi); yLoc = sin(PrevPhi);

				RadialVectInCutPlane = xLoc*Ex + yLoc*Ey;
				CutPlaneNormal = CylAxVect^RadialVectInCutPlane;

				radThg Loc_hgDuplCutted = LowerNewPair_int_hg.Handler_g;
				radTPair_int_hg LocLowerNewPair_int_hg, LocUpperNewPair_int_hg;
				((radTg3d*)(LowerNewPair_int_hg.Handler_g.rep))->CutItself(CuttingPlane, Loc_hgDuplCutted, LocLowerNewPair_int_hg, LocUpperNewPair_int_hg, radPtr, &SubdivOptionsForCut);

				if(Loc_hgDuplCutted.rep != LowerNewPair_int_hg.Handler_g.rep) ThereWereNoActualCuts = 0;

				if(LocUpperNewPair_int_hg.Handler_g.rep != 0)
					LowerNewPair_int_hg.Handler_g = LocUpperNewPair_int_hg.Handler_g;

				if(CuttingGoesFromMiddle && (iPhi == 0))
				{
					if(LocLowerNewPair_int_hg.Handler_g.rep != 0)
					{
						FirstLowerPair_int_hg.Handler_g = LocLowerNewPair_int_hg.Handler_g;
					}
				}
			}
			if(!AxisCrossesVolume)
			{
				if(UpperNewPair_int_hg.Handler_g.rep != 0) hgDupl = UpperNewPair_int_hg.Handler_g;
			}
			if(LowerNewPair_int_hg.Handler_g.rep != 0)
			{
				radTg3d* g3dNewPtr = (radTg3d*)(LowerNewPair_int_hg.Handler_g.rep);
				g3dNewPtr->EraseAllTransformations();
				((radTPolyhedron*)((radTg3dRelax*)g3dNewPtr))->AuxPairOfDouble = radTPairOfDouble(PrevPhi, Phi);

				NewElKey = (pSubdivOptions->PutNewStuffIntoGenCont)? radPtr->AddElementToContainer(LowerNewPair_int_hg.Handler_g) : (iPhi + 1);
				GroupInPlaceOfThisPtr->AddElement(NewElKey, LowerNewPair_int_hg.Handler_g);
			}
		}
		else
		{
			if(ThereWereNoActualCuts)
			{
				double CheckPrevPhi = PrevPhi;
				for(;;)
				{
					if(CheckPrevPhi >= TwoPI) CheckPrevPhi -= TwoPI;
					else break;
				}
				if(AngleIsBetween(ActualLimits[0], CheckPrevPhi, CheckPhi) && AngleIsBetween(ActualLimits[1], CheckPrevPhi, CheckPhi))
				{
					BoundingPhi1 = CheckPrevPhi;
					BoundingPhi2 = CheckPhi;
				}
			}

			if(CuttingGoesFromMiddle && (iPhi == 0))
			{
				xLoc = cos(PrevPhi); yLoc = sin(PrevPhi);

				RadialVectInCutPlane = xLoc*Ex + yLoc*Ey;
				CutPlaneNormal = CylAxVect^RadialVectInCutPlane;

				radThg Loc_hgDuplCutted = hgDupl;
				radTPair_int_hg LocLowerNewPair_int_hg, LocUpperNewPair_int_hg;
				((radTg3d*)(hgDupl.rep))->CutItself(CuttingPlane, Loc_hgDuplCutted, LocLowerNewPair_int_hg, LocUpperNewPair_int_hg, radPtr, &SubdivOptionsForCut);

				if(Loc_hgDuplCutted.rep != hgDupl.rep) ThereWereNoActualCuts = 0;
			
				if(LocUpperNewPair_int_hg.Handler_g.rep != 0)
					UpperNewPair_int_hg.Handler_g = LocUpperNewPair_int_hg.Handler_g;

				if(LocLowerNewPair_int_hg.Handler_g.rep != 0)
				{
					FirstLowerPair_int_hg.Handler_g = LocLowerNewPair_int_hg.Handler_g;
				}
				PreviousCutWasOK = 1;
			}
		}

		if(CompareWithActualLimits && PreviousCutWasOK && (!TryToMakeCut))
		{
			if(UpperNewPair_int_hg.Handler_g.rep != 0)
			{
				radTg3d* g3dNewPtr = (radTg3d*)(UpperNewPair_int_hg.Handler_g.rep);
				g3dNewPtr->EraseAllTransformations();
				((radTPolyhedron*)((radTg3dRelax*)g3dNewPtr))->AuxPairOfDouble = radTPairOfDouble(PrevPhi, Phi);

				NewElKey = (pSubdivOptions->PutNewStuffIntoGenCont)? radPtr->AddElementToContainer(UpperNewPair_int_hg.Handler_g) : (iPhi + 1);
				GroupInPlaceOfThisPtr->AddElement(NewElKey, UpperNewPair_int_hg.Handler_g);
			}

			if(CuttingGoesFromMiddle)
			{
				if(FirstLowerPair_int_hg.Handler_g.rep != 0) hgDupl = FirstLowerPair_int_hg.Handler_g;
			}
			PreviousCutWasOK = 0;

			radThg hgZero;
			UpperNewPair_int_hg.Handler_g = hgZero;
		}
	}

	if(AxisCrossesVolume)
	{
		xLoc = cos(Limits[0]); yLoc = sin(Limits[0]);

		RadialVectInCutPlane = xLoc*Ex + yLoc*Ey;
		CutPlaneNormal = CylAxVect^RadialVectInCutPlane;

		radTPair_int_hg LocLowerNewPair_int_hg, LocUpperNewPair_int_hg;
		radThg Loc_hgDuplCutted;
		if(UpperNewPair_int_hg.Handler_g.rep != 0)
		{
			Loc_hgDuplCutted = UpperNewPair_int_hg.Handler_g;
			((radTg3d*)(UpperNewPair_int_hg.Handler_g.rep))->CutItself(CuttingPlane, Loc_hgDuplCutted, LocLowerNewPair_int_hg, LocUpperNewPair_int_hg, radPtr, &SubdivOptionsForCut);

			if(Loc_hgDuplCutted.rep != UpperNewPair_int_hg.Handler_g.rep) ThereWereNoActualCuts = 0;
		
			if(LocLowerNewPair_int_hg.Handler_g.rep != 0)
				UpperNewPair_int_hg.Handler_g = LocLowerNewPair_int_hg.Handler_g;
		}
	}

	if(hgDupl == FirstLowerPair_int_hg.Handler_g) UpperNewPair_int_hg.Handler_g = hgDupl;

	if(UpperNewPair_int_hg.Handler_g.rep != 0)
	{
		radTg3d* g3dNewPtr = (radTg3d*)(UpperNewPair_int_hg.Handler_g.rep);
		g3dNewPtr->EraseAllTransformations();
		((radTPolyhedron*)((radTg3dRelax*)g3dNewPtr))->AuxPairOfDouble = radTPairOfDouble(Phi, Limits[1]);

		NewElKey = (pSubdivOptions->PutNewStuffIntoGenCont)? radPtr->AddElementToContainer(UpperNewPair_int_hg.Handler_g) : (AmOfFurtherPhiCuts + 1);
		GroupInPlaceOfThisPtr->AddElement(NewElKey, UpperNewPair_int_hg.Handler_g);
	}

	if(ThereWereNoActualCuts)
	{
		double CheckPrevPhi = Phi;
		for(;;)
		{
			if(CheckPrevPhi >= TwoPI) CheckPrevPhi -= TwoPI;
			else break;
		}
		double CheckPhi = Limits[0];
		for(;;)
		{
			if(CheckPhi >= TwoPI) CheckPhi -= TwoPI;
			else break;
		}
		if(AngleIsBetween(ActualLimits[0], CheckPrevPhi, CheckPhi) && AngleIsBetween(ActualLimits[1], CheckPrevPhi, CheckPhi))
		{
			BoundingPhi1 = CheckPrevPhi;
			BoundingPhi2 = CheckPhi;
		}
	}

	SetMessageChar(0); // Releasing this member variable

	if(!ThereWereNoActualCuts)
	{
		In_hg = NewHandle;

		GroupInPlaceOfThisPtr->SetMessageChar(0);
		GroupInPlaceOfThisPtr->MessageChar = 1;
	}
	else
	{
		if((BoundingPhi1 != 0.) || (BoundingPhi2 != 0.))
		{
			AuxPairOfDouble.First = BoundingPhi1;
			AuxPairOfDouble.Second = BoundingPhi2;
		}
	}

	radPtr->RecognizeRecMagsInPolyhedrons = OldRecognizeRecMagsInPolyhedrons;
	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::FindEdgePointsOverEllipseSet0(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg In_hg, TVector3d* EdgePoints, double* Limits, radTSubdivOptions* pSubdivOptions)
{
	double& ka = SubdivArray[0];
	double& qa = SubdivArray[1];
	if((pSubdivOptions->SubdivisionParamCode == 0) && (int(ka) == 1)) return 1;

	TVector3d& CylAxVect = pSubdivSpec->CylAxVect;
	TVector3d& PointOnCylAx = pSubdivSpec->PointOnCylAx;
	TVector3d& PointOnEllAx = pSubdivSpec->PointOnEllAx;

	EdgePoints[0] = pSubdivSpec->PointOnCylAx;
	EdgePoints[1] = pSubdivSpec->PointOnEllAx;
	Limits[0] = 0.;

	TVector3d PointOnEllAx_mi_PointOnCylAx = PointOnEllAx - PointOnCylAx;
	TVector3d Ex = PointOnEllAx_mi_PointOnCylAx - (PointOnEllAx_mi_PointOnCylAx*CylAxVect)*CylAxVect;

	double DelA = sqrt(Ex.x*Ex.x + Ex.y*Ex.y + Ex.z*Ex.z);
	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		ka = (ka < DelA)? Round(DelA/ka) : 1.;
	}
	Limits[1] = DelA;

	return 1;
}

//-------------------------------------------------------------------------

int radTPolyhedron::FindEdgePointsOverEllipseSet(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg In_hg, TVector3d* EdgePoints, double* Limits, radTSubdivOptions* pSubdivOptions)
{
	double& ka = SubdivArray[0];
	double& qa = SubdivArray[1];
	if((pSubdivOptions->SubdivisionParamCode == 0) && (int(ka) == 1)) return 1;

	radTg3d* g3dPtr = (radTg3d*)(In_hg.rep);
	//radTCast Cast;
	radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
	if(pGroup != 0)
	{
		TVector3d LocEdgePoints[2];
		double LocLimits[2];
		double aMin = 1.E+23, aMax = 0.;
		for(radTmhg::const_iterator iter = pGroup->GroupMapOfHandlers.begin(); iter != pGroup->GroupMapOfHandlers.end(); ++iter)
		{
			radThg hgLoc = (*iter).second;
			if(!FindEdgePointsOverEllipseSet(SubdivArray, pSubdivSpec, hgLoc, LocEdgePoints, LocLimits, pSubdivOptions)) return 0;
			if(*LocLimits < aMin)
			{
				aMin = *LocLimits; *EdgePoints = *LocEdgePoints;
			}
			if(LocLimits[1] > aMax)
			{
				aMax = LocLimits[1]; EdgePoints[1] = LocEdgePoints[1];
			}
		}
		*Limits = aMin; Limits[1] = aMax;
		
		if(pSubdivOptions->SubdivisionParamCode == 1)
		{
			double SizeA = Limits[1] - Limits[0];
			double& ka = *SubdivArray;
			ka = (ka < SizeA)? Round(SizeA/ka) : 1.;
		}
		return (Limits[1] >= *Limits)? 1 : 0;
	}

	radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast(g3dPtr);
	if(g3dRelaxPtr == 0) return 0;
	radTPolyhedron* pPolyhdr = radTCast::PolyhedronCast(g3dRelaxPtr);
	if(pPolyhdr == 0) return 0;

	double aMin = 1.E+23, aMax = 0.;
	TVector3d& CylAxVect = pSubdivSpec->CylAxVect;
	TVector3d& PointOnCylAx = pSubdivSpec->PointOnCylAx;
	TVector3d& PointOnEllAx = pSubdivSpec->PointOnEllAx;
	double Ratio = pSubdivSpec->EllAxRatio;

	TVector3d PointOnEllAx_mi_PointOnCylAx = PointOnEllAx - PointOnCylAx;
	TVector3d Ex = PointOnEllAx_mi_PointOnCylAx - (PointOnEllAx_mi_PointOnCylAx*CylAxVect)*CylAxVect;
	double InvLenEx = 1./sqrt(Ex.x*Ex.x + Ex.y*Ex.y + Ex.z*Ex.z);
	Ex = InvLenEx*Ex;
	TVector3d Ey = CylAxVect^Ex;

	radTPolygon* PgnPtr;
	radTrans* TransPtr;
	radTHandlePgnAndTrans* HandlePgnAndTransPtr;

	const double HalfPI = 1.5707963267949;
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double TolForZero = 10.*RelAbsTol[1];

	for(int i=0; i<pPolyhdr->AmOfFaces; i++)
	{
		HandlePgnAndTransPtr = &(pPolyhdr->VectHandlePgnAndTrans[i]);
		PgnPtr = HandlePgnAndTransPtr->PgnHndl.rep;
		TransPtr = HandlePgnAndTransPtr->TransHndl.rep;

		for(int k=0; k<PgnPtr->AmOfEdgePoints; k++)
		{
			TVector2d& p2d = PgnPtr->EdgePointsVector[k];
			TVector3d p3d(p2d.x, p2d.y, PgnPtr->CoordZ); p3d = TransPtr->TrPoint(p3d);
			//TVector3d p3d(p2d.x, p2d.y, PgnPtr->CoordZ); p3d = TransPtr->TrPoint(p3d) + CentrPoint; //OC090908 ??

			TVector3d aVect = p3d - PointOnCylAx;
			TVector3d OrtComp = aVect - (aVect*CylAxVect)*CylAxVect;

			double Phi, a;
			char aDefined = 0;
			double xLoc = OrtComp*Ex, yLoc = OrtComp*Ey;
			if(fabs(xLoc) < TolForZero)
			{
				if(yLoc > TolForZero) { Phi = HalfPI; a = yLoc/Ratio;}
				else if(yLoc < -TolForZero) { Phi = 3.*HalfPI; a = -yLoc/Ratio;}
				else { a = 0.; Phi = 0.;}
				aDefined = 1;
			}
			else if(xLoc > TolForZero)
			{
				if(yLoc > TolForZero) Phi = atan(yLoc/(Ratio*xLoc));
				else if(yLoc < -TolForZero) Phi = 4.*HalfPI + atan(yLoc/(Ratio*xLoc));
				else { a = xLoc; Phi = 0.; aDefined = 1;}
			}
			else if(xLoc < -TolForZero)
			{
				Phi = 2.*HalfPI + atan(yLoc/(Ratio*xLoc));
			}
			if(!aDefined) a = xLoc/cos(Phi);

			if(a < aMin) { aMin = a; *EdgePoints = p3d;}
			if(a > aMax) { aMax = a; EdgePoints[1] = p3d;}
		}
	}
	*Limits = aMin; Limits[1] = aMax;

	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		double SizeA = Limits[1] - Limits[0];
		double& ka = *SubdivArray;
		ka = (ka < SizeA)? Round(SizeA/ka) : 1.;
	}
	return (Limits[1] >= *Limits)? 1 : 0;
}

//-------------------------------------------------------------------------

void radTPolyhedron::FindEdgePointsOverPhiAndAxForCylSubd(radTCylindricSubdivSpec* pSubdivSpec, TVector3d* EdgePointsForCylSubd, double* PhiAndAxLimits)
{
	double MinRe2=1.E+23, MaxRe2=0., MinFunPh=1.E+23, MaxFunPh=-1.E+23, MinAxProd=1.E+23, MaxAxProd=-1.E+23;
	
	TVector3d MinVertPoCylAx, MaxVertPoCylAx;
	TVector3d MinVertPoR, MaxVertPoR;
	TVector3d MinVertPoPh, MaxVertPoPh;

	TVector3d& CylAxVect = pSubdivSpec->CylAxVect;
	TVector3d& PointOnCylAx = pSubdivSpec->PointOnCylAx;
	TVector3d MinOrtComp(0.,0.,0.), MaxOrtComp(0.,0.,0.);

	TVector3d OrtBasis[2];
	FindTwoOrtogonalVectors(CylAxVect, OrtBasis);

	char AxisCrossesVolume = 0;
	SetMessageChar(0);

	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double AbsZeroToler = 50.*RelAbsTol[1];
	double RelToler = 50.*RelAbsTol[0];

	radTPolygon* PgnPtr;
	radTrans* TransPtr;
	radTHandlePgnAndTrans* HandlePgnAndTransPtr;

	for(int i=0; i<AmOfFaces; i++)
	{
		HandlePgnAndTransPtr = &(VectHandlePgnAndTrans[i]);
		PgnPtr = HandlePgnAndTransPtr->PgnHndl.rep;
		TransPtr = HandlePgnAndTransPtr->TransHndl.rep;

		for(int k=0; k<PgnPtr->AmOfEdgePoints; k++)
		{
			TVector2d& p2d = PgnPtr->EdgePointsVector[k];
			TVector3d p3d(p2d.x, p2d.y, PgnPtr->CoordZ); p3d = TransPtr->TrPoint(p3d);
			//TVector3d p3d(p2d.x, p2d.y, PgnPtr->CoordZ); p3d = TransPtr->TrPoint(p3d) + CentrPoint; //OC090908

			TVector3d aVect = p3d - PointOnCylAx;

			double AxProd = aVect*CylAxVect;
			if(AxProd < MinAxProd) { MinAxProd = AxProd; MinVertPoCylAx = p3d;}
			if(AxProd > MaxAxProd) { MaxAxProd = AxProd; MaxVertPoCylAx = p3d;}

			TVector3d OrtComp = aVect - AxProd*CylAxVect;

			double Re2 = OrtComp*OrtComp;
			if(Re2 < MinRe2) { MinRe2 = Re2; MinVertPoR = p3d;}
			if(Re2 > MaxRe2) { MaxRe2 = Re2; MaxVertPoR = p3d;}

			if(!((fabs(OrtComp.x) < AbsZeroToler) && (fabs(OrtComp.y) < AbsZeroToler) && (fabs(OrtComp.z) < AbsZeroToler)))
			{
				if((NormAbs(MinOrtComp) < AbsZeroToler) && (NormAbs(MaxOrtComp) < AbsZeroToler))
				{
					MinOrtComp = OrtComp; MaxOrtComp = OrtComp;
					MinVertPoPh = p3d; MaxVertPoPh = p3d;
				}
				else
				{
					double NormOrtComp = NormAbs(OrtComp);
					
					TVector3d MinOrtComp_x_OrtComp = MinOrtComp^OrtComp;
					if(NormAbs(MinOrtComp_x_OrtComp) > RelToler*NormAbs(MinOrtComp)*NormOrtComp)
					{
						double ProdMixedMin = MinOrtComp_x_OrtComp*CylAxVect;
						if(ProdMixedMin < 0.) 
						{
							MinOrtComp = OrtComp; MinVertPoPh = p3d;
							
							TVector3d MinOrtComp_x_MaxOrtComp = MinOrtComp^MaxOrtComp;
							if(NormAbs(MinOrtComp_x_MaxOrtComp) > RelToler*NormAbs(MinOrtComp)*NormAbs(MaxOrtComp))
							{
								if((!AxisCrossesVolume) && (MinOrtComp_x_MaxOrtComp*CylAxVect < 0.)) 
								{
									AxisCrossesVolume = 1;
								}
							}
						}
					}
					
					TVector3d MaxOrtComp_x_OrtComp = MaxOrtComp^OrtComp;
					if(NormAbs(MaxOrtComp_x_OrtComp) > RelToler*NormAbs(MaxOrtComp)*NormOrtComp)
					{
						double ProdMixedMax = MaxOrtComp_x_OrtComp*CylAxVect;
						if(ProdMixedMax > 0.) 
						{ 
							MaxOrtComp = OrtComp; MaxVertPoPh = p3d;

							TVector3d MinOrtComp_x_MaxOrtComp = MinOrtComp^MaxOrtComp;
							if(NormAbs(MinOrtComp_x_MaxOrtComp) > RelToler*NormAbs(MinOrtComp)*NormAbs(MaxOrtComp))
							{
								if((!AxisCrossesVolume) && (MinOrtComp_x_MaxOrtComp*CylAxVect < 0.)) 
								{
									AxisCrossesVolume = 1;
								}
							}
						}
					}
				}
			}
		}
	}
	if(pSubdivSpec->EllAxNotDefined)
	{
		pSubdivSpec->PointOnEllAx = MinVertPoPh;
		pSubdivSpec->EllAxNotDefined = 0;
	}

	if(AxisCrossesVolume) 
	{
		MaxVertPoPh = MinVertPoPh;
		SetMessageChar(1); // MessageChar in g3d is used to inform "on-top" algorithms that Axis Crosses this Volume
	}

	radTCylindricSubdivSpec LocSpec = *pSubdivSpec;
	LocSpec.EllAxRatio = 1.;

	double aDummy;
	EdgePointsForCylSubd[0] = MinVertPoPh;
	FindEllipticCoordOfPoint(&LocSpec, MinVertPoPh, aDummy, PhiAndAxLimits[0]);
	EdgePointsForCylSubd[1] = MaxVertPoPh;
	FindEllipticCoordOfPoint(&LocSpec, MaxVertPoPh, aDummy, PhiAndAxLimits[1]);

	EdgePointsForCylSubd[2] = MinVertPoCylAx;
	PhiAndAxLimits[2] = MinAxProd;
	EdgePointsForCylSubd[3] = MaxVertPoCylAx;
	PhiAndAxLimits[3] = MaxAxProd;

	EdgePointsForCylSubd[4] = MinVertPoR;
	PhiAndAxLimits[4] = MinRe2;
	EdgePointsForCylSubd[5] = MaxVertPoR;
	PhiAndAxLimits[5] = MaxRe2;

	AuxPairOfDouble.First = PhiAndAxLimits[0]; // To keep limits over phi; values should be between 0 and 2*PI
	AuxPairOfDouble.Second = PhiAndAxLimits[1];

	const double TwoPI = 2.*3.1415926535898;
	if((PhiAndAxLimits[1] < PhiAndAxLimits[0]) && (!AxisCrossesVolume)) 
	{
		PhiAndAxLimits[1] += TwoPI;
	}
	if(AxisCrossesVolume) 
	{
		PhiAndAxLimits[0] = 0.;
		PhiAndAxLimits[1] = TwoPI;
	}
}

//-------------------------------------------------------------------------

int radTPolyhedron::SubdivideItself(double* InSubdArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	double kx = InSubdArray[0], ky = InSubdArray[2], kz = InSubdArray[4];
	double qx = InSubdArray[1], qy = InSubdArray[3], qz = InSubdArray[5];
	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		TVector3d Directions[3];
		*Directions = TVector3d(1.,0.,0.); Directions[1] = TVector3d(0.,1.,0.); Directions[2] = TVector3d(0.,0.,1.);

		radTrans ResTransf;
		short SomethingFound = 0;
		if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty()))
		{
			FindResTransfWithMultOne(ResTransf, SomethingFound);
		}
		if(ConsiderOnlyWithTrans && (!g3dListOfTransform.empty()))
		{
			FindInnerTransfWithMultOne(ResTransf, SomethingFound);
		}
		if(SomethingFound)
		{
			TVector3d* tDirections = Directions;
			*tDirections = ResTransf.TrBiPoint_inv(*tDirections); tDirections++;
			*tDirections = ResTransf.TrBiPoint_inv(*tDirections); tDirections++;
			*tDirections = ResTransf.TrBiPoint_inv(*tDirections);
		}

		double Sizes[3];
		EstimateSize(Directions, Sizes, 3);
		kx = (kx < *Sizes)? Round(*Sizes/kx) : 1.;
		ky = (ky < Sizes[1])? Round(Sizes[1]/ky) : 1.;
		kz = (kz < Sizes[2])? Round(Sizes[2]/kz) : 1.;
	}
	double SubdivArray[] = {1.,0.,0.,kx,qx, 0.,1.,0.,ky,qy, 0.,0.,1.,kz,qz};

	radTSubdivOptions LocSubdivOptions = *pSubdivOptions;
	LocSubdivOptions.SubdivisionParamCode = 0;

	int SubdOK = SubdivideItselfByParPlanes(SubdivArray, 3, In_hg, radPtr, &LocSubdivOptions);
	if(!SubdOK) SomethingIsWrong = 1;
	return SubdOK;
}

//-------------------------------------------------------------------------

void radTPolyhedron::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique)
{
	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double AbsTol = RelAbsTol[1];
	double AbsTolE2 = AbsTol*AbsTol;

	radTVectorOfVector3d LocPts;

	for(radTVectHandlePgnAndTrans::iterator FaceIter = VectHandlePgnAndTrans.begin(); FaceIter != VectHandlePgnAndTrans.end(); ++FaceIter)
	{
		radTPolygon* PgnPtr = (*FaceIter).PgnHndl.rep;
		radTrans* TransPtr = (*FaceIter).TransHndl.rep;
		double LocZ = PgnPtr->CoordZ;
		for(radTVect2dVect::iterator PointIter = PgnPtr->EdgePointsVector.begin(); PointIter != PgnPtr->EdgePointsVector.end(); ++PointIter)
		{
			TVector2d& p2d = *PointIter;
			TVector3d p3d = TVector3d(p2d.x, p2d.y, LocZ);
			p3d = TransPtr->TrPoint(p3d);
			//p3d = TransPtr->TrPoint(p3d) + CentrPoint; //OC090908

			bool ThisPtAlreadyExists = false;
			if(EnsureUnique)
			{
				for(int k=0; k<(int)LocPts.size(); k++)
				{
					TVector3d dp = LocPts[k] - p3d;
					if((dp.x*dp.x + dp.y*dp.y + dp.z*dp.z) <= AbsTolE2)
					{
						ThisPtAlreadyExists = true;
						break;
					}
				}
			}
			if(!ThisPtAlreadyExists) LocPts.push_back(p3d);
		}
	}
	int SizeLocPts = (int)LocPts.size();
	if(SizeLocPts == 0) return;

	for(int j=0; j<SizeLocPts; j++) OutVect.push_back(LocPts[j]);

	LocPts.erase(LocPts.begin(), LocPts.end());
}

//-------------------------------------------------------------------------

void radTPolyhedron::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
 //to move to base class?
	radTrans bufTrans, *pTrans=0;
	TVector3d cenPointInLabFr, vZero(0.,0.,0.);
	GetTrfAndCenPointInLabFrame(pBaseTrans, bufTrans, pTrans, cenPointInLabFr);
	radTPairOfVect3d Pair(cenPointInLabFr, vZero);

	if(pFieldKey->M_)
	{
		TVector3d resM = Magn;
		if((pM_LinCoef != 0) && (mLinTreat != 0)) //treat Lin. as absolute
		{
			resM += ((*pM_LinCoef)*CentrPoint); //CentrPoint in Loc. Frame - to check!
		}
		Pair.V2 = (pTrans == 0)? resM : pTrans->TrVectField(resM);
	}
	else if(pFieldKey->J_)
	{
		TVector3d resJ = J;
		if((pJ_LinCoef != 0) && (mLinTreat != 0)) //treat Lin. as absolute
		{
			resJ += ((*pJ_LinCoef)*CentrPoint); //CentrPoint in Loc. Frame - to check!
		}
		Pair.V2 = (pTrans == 0)? resJ : pTrans->TrVectField(resJ);
	}
	else
	{
		radTCompCriterium CompCriterium;
		radTField Field(*pFieldKey, CompCriterium, cenPointInLabFr, vZero, vZero, vZero, vZero, 0.);
		g3dSrcPtr->B_genComp(&Field);
		Pair.V2 = (pFieldKey->B_)? Field.B : ((pFieldKey->H_)? Field.H : ((pFieldKey->A_)? Field.A : vZero));
	}
	pVectPairOfVect3d->push_back(Pair);
}

//-------------------------------------------------------------------------

void radTPolyhedron::AttemptToCreateConvexPolyhedronFromTwoBaseFaces(const radTHandlePgnAndTrans& inHandleBasePgnAndTrf1, const radTHandlePgnAndTrans& inHandleBasePgnAndTrf2)
{//sets SomethingIsWrong = 1 in case of any problem

	radTPolygon *pPgn1 = inHandleBasePgnAndTrf1.PgnHndl.rep;
	radTrans *pTrf1 = inHandleBasePgnAndTrf1.TransHndl.rep;

	radTPolygon *pPgn2 = inHandleBasePgnAndTrf2.PgnHndl.rep;
	radTrans *pTrf2 = inHandleBasePgnAndTrf2.TransHndl.rep;

	double zc1 = pPgn1->CoordZ, zc2 = pPgn2->CoordZ;
	TVector2d vCenPoint1 = pPgn1->CentrPoint, vCenPoint2 = pPgn2->CentrPoint;

	TVector3d vNormZ(0,0,1);
	TVector3d vNorm1 = pTrf1->TrBiPoint(vNormZ), vNorm2 = pTrf2->TrBiPoint(vNormZ);
	TVector3d vCenPoint1_3d(vCenPoint1.x, vCenPoint1.y, zc1), vCenPoint2_3d(vCenPoint2.x, vCenPoint2.y, zc2);
	vCenPoint1_3d = pTrf1->TrPoint(vCenPoint1_3d); 
	vCenPoint2_3d = pTrf2->TrPoint(vCenPoint2_3d);

	double RelZeroToler = 1.E-09;
	RelZeroToler = 10.*((RelZeroToler>radCR.RelRand)? RelZeroToler : radCR.RelRand);
	double AbsTol = RelZeroToler*0.5*(pPgn1->EstimateTypSize() + pPgn2->EstimateTypSize());
	double arTol[] = {RelZeroToler, AbsTol};

	if((!CheckIfAllPolygonVerticesAreOnOneSideOfPlane(inHandleBasePgnAndTrf2, vCenPoint1_3d, vNorm1, AbsTol)) ||
	   (!CheckIfAllPolygonVerticesAreOnOneSideOfPlane(inHandleBasePgnAndTrf1, vCenPoint2_3d, vNorm2, AbsTol)))
	{
		SomethingIsWrong = 1;
		radTSend::ErrorMessage("Radia::Error110"); 
		return;
	}

	//generating mantle (side faces)
	//generate array of unique vertex points and indexes of these points for two bases faces
	vector<TVector3d> vectVertexPoints;
	vector<vector<int> > vectIndAllFaces;
	vector<int> vectIndFacePgn;

	CollectAndMapUniquePolygonPoints(inHandleBasePgnAndTrf1, vectVertexPoints, vectIndFacePgn, AbsTol);
	vectIndAllFaces.push_back(vectIndFacePgn);
	vectIndFacePgn.erase(vectIndFacePgn.begin(), vectIndFacePgn.end());

	CollectAndMapUniquePolygonPoints(inHandleBasePgnAndTrf2, vectVertexPoints, vectIndFacePgn, AbsTol);
	vectIndAllFaces.push_back(vectIndFacePgn);
	vectIndFacePgn.erase(vectIndFacePgn.begin(), vectIndFacePgn.end());

	GenerateSideFacesContainingSegmentsOfBaseFace(vectVertexPoints, vectIndAllFaces, 0, arTol);
	GenerateSideFacesContainingSegmentsOfBaseFace(vectVertexPoints, vectIndAllFaces, 1, arTol);

	int lenArVertexPoints = (int)vectVertexPoints.size();
	TVector3d *arVertexPoints = CAuxParse::Vect2Ar(vectVertexPoints);

	AmOfFaces = (int)vectIndAllFaces.size();
	int *arIndAllFacesLengths = 0;
	int **arIndAllFaces = CAuxParse::Vect2Ar2D(vectIndAllFaces, arIndAllFacesLengths);

	FillInVectHandlePgnAndTrans(arVertexPoints, lenArVertexPoints, arIndAllFaces, arIndAllFacesLengths);
	if(SomethingIsWrong) return;
	DefineCentrPoint(arVertexPoints, lenArVertexPoints);

	if(arVertexPoints != 0) delete[] arVertexPoints;
	if(arIndAllFacesLengths != 0) delete[] arIndAllFacesLengths;
	if(arIndAllFaces != 0)
	{
		for(int i=0; i<AmOfFaces; i++)
		{
			if(arIndAllFaces[i] != 0) delete[] arIndAllFaces[i];
		}
		delete[] arIndAllFaces;
	}
	vectVertexPoints.erase(vectVertexPoints.begin(), vectVertexPoints.end());
	vectIndAllFaces.erase(vectIndAllFaces.begin(), vectIndAllFaces.end());
	vectIndFacePgn.erase(vectIndFacePgn.begin(), vectIndFacePgn.end());
}

//-------------------------------------------------------------------------

void radTPolyhedron::GenerateSideFacesContainingSegmentsOfBaseFace(const vector<TVector3d>& vectVertexPoints, vector<vector<int> >& vectIndAllFaces, int indBaseFace, double* arTol)
{
	int numExistingFaces = (int)vectIndAllFaces.size();
	if((indBaseFace < 0) || (indBaseFace >= numExistingFaces)) return;

	int indOtherBaseFace = (indBaseFace == 0)? 1 : 0;
	if(indOtherBaseFace >= numExistingFaces) return;

	vector<int> vectBaseFace = vectIndAllFaces[indBaseFace];
	int numPointsInBaseFace = (int)vectBaseFace.size();

	vector<int> vectOtherBaseFace = vectIndAllFaces[indOtherBaseFace];
	int numPointsInOtherBaseFace = (int)vectOtherBaseFace.size();

	int totNumVertexPoints = (int)vectVertexPoints.size();
	double AbsTol = arTol[1];

	//loop over segments of base face
	vector<int> vectNewFace, vectNewFaceAux;
	int indP1 = vectBaseFace[0];
	const TVector3d *pP1 = &(vectVertexPoints[indP1]), *pP2;
	TVector3d vPlaneNorm;
	for(int i=0; i<numPointsInBaseFace; i++)
	{
		int i2 = i + 1;
		if(i2 == numPointsInBaseFace) i2 = 0;

		int indP2 = vectBaseFace[i2];
		pP2 = &(vectVertexPoints[indP2]);

		//look for vertex points in other base face, which can define new face
		for(int j=0; j<numPointsInOtherBaseFace; j++)
		{
			int indP = vectOtherBaseFace[j];
			const TVector3d &P = vectVertexPoints[indP];

			//check whether this point doesn't belong to the line passing through *pP1, *pP2
			if(CheckIfThreePointsAreOnOneLine(*pP1, *pP2, P, AbsTol)) continue;

			DefineNormalVia3Points(*pP1, *pP2, P, vPlaneNorm);
			vPlaneNorm.Normalize();

			vectNewFace.erase(vectNewFace.begin(), vectNewFace.end());
			int signPrevScalProd = 0;
			bool newFaceIsGood = true;
			//check whether all other points of the other base face are located on the same side of the test plane as the points of the main base 
			for(int k=0; k<totNumVertexPoints; k++)
			{
				const TVector3d &testP = vectVertexPoints[k];
				if(PracticallyEqual(testP, *pP1, AbsTol) || PracticallyEqual(testP, *pP2, AbsTol) || PracticallyEqual(testP, P, AbsTol)) continue;

				TVector3d dTestP = testP - P;
				double testScalProd = dTestP*vPlaneNorm;
				int signTestScalProd = 0;
				if(testScalProd < -AbsTol) signTestScalProd = -1;
				else if(testScalProd > AbsTol) signTestScalProd = 1;

				//if(signPrevScalProd*testScalProd < 0) { newFaceIsGood = false; break;}
				if(signPrevScalProd*signTestScalProd < 0) //OC140209
				{ 
					newFaceIsGood = false; break;
				}

				if(signTestScalProd != 0) signPrevScalProd = signTestScalProd;
				else
				{
					vectNewFace.push_back(k); //one more vertex point belongs to this face
				}
			}
			if(newFaceIsGood)
			{
				vectNewFace.push_back(indP1);
				vectNewFace.push_back(indP2);
				vectNewFace.push_back(indP);
				//break;

				//check whether the "new" face is not already present in vectIndAllFaces
				if(!CAuxParse::CheckIfVectElemArePresent(vectNewFace, vectIndAllFaces)) //OC150209
				{
					int numVertexPointsInNewFace = (int)vectNewFace.size();
					if(numVertexPointsInNewFace > 3)
					{
						//make sure that vertex indices are listed continuously, without creating self-intersecting polygon
						ReorderPointsToEnsureNonSelfIntersectingPolygon(vectVertexPoints, vectNewFace, arTol);
					}
					vectIndAllFaces.push_back(vectNewFace); //add only if not already present
					vectNewFace.erase(vectNewFace.begin(), vectNewFace.end());
					break;
				}
				else
				{
					vectNewFace.erase(vectNewFace.begin(), vectNewFace.end());
					//continue looking for points to create new face with given segment (*pP1, *pP2)
				}
			}
		}
		
		//check whether the "new" face is not already present in vectIndAllFaces
		//if(!CAuxParse::CheckIfVectElemArePresent(vectNewFace, vectIndAllFaces)) //OC150209
		//{
		//	int numVertexPointsInNewFace = (int)vectNewFace.size();
		//	if(numVertexPointsInNewFace > 3)
		//	{
		//		//make sure that vertex indices are listed continuously, without creating self-intersecting polygon
		//		ReorderPointsToEnsureNonSelfIntersectingPolygon(vectVertexPoints, vectNewFace, arTol);
		//	}
		//	vectIndAllFaces.push_back(vectNewFace); //add only if not already present
		//}
		//vectNewFace.erase(vectNewFace.begin(), vectNewFace.end());

		pP1 = pP2;
		indP1 = indP2;
	}
}

//-------------------------------------------------------------------------

void radTPolyhedron::ReorderPointsToEnsureNonSelfIntersectingPolygon(const vector<TVector3d>& vectPoints, vector<int>& vectIndPgnPoints, double* arTol)
{//may not work for non-convex polygons
	int numInd = (int)vectIndPgnPoints.size();
	if(numInd <= 0) return;

	double RelTol = arTol[0];
	double RelTolE2 = RelTol*RelTol;

	//find face center point
	TVector3d vCenPoint(0,0,0);
	for(int i=0; i<numInd; i++) vCenPoint += vectPoints[vectIndPgnPoints[i]];
	vCenPoint *= (1./numInd);

	TVector3d r0 = vectPoints[vectIndPgnPoints[0]] - vCenPoint;
	r0.Normalize();

	//calculating reference normal to polygon plane
	TVector3d vRefNorm;
	for(int j=1; j<numInd; j++)
	{
		TVector3d rj = vectPoints[vectIndPgnPoints[j]] - vCenPoint;
		rj.Normalize();
		vRefNorm = r0^rj;
		if(vRefNorm.AmpE2() > RelTolE2) break;
	}
	vRefNorm.Normalize();

	//re-ordering points
	list<radTPairIntDouble> listIndAngle;
	radTPairIntDouble firstPair(vectIndPgnPoints[0], 0.);
	listIndAngle.push_back(firstPair);

	for(int k=1; k<numInd; k++)
	{
		int curInd = vectIndPgnPoints[k];
		TVector3d rk = vectPoints[curInd] - vCenPoint;
		rk.Normalize();
		double curAngle = AngleBwUnitVectors(r0, rk, &vRefNorm);
		radTPairIntDouble curPair(curInd, curAngle);
		listIndAngle.push_back(curPair);
	}
	listIndAngle.sort(radTPairIntDouble::less);

	vectIndPgnPoints.erase(vectIndPgnPoints.begin(), vectIndPgnPoints.end());
	for(list<radTPairIntDouble>::const_iterator iter = listIndAngle.begin(); iter != listIndAngle.end(); ++iter)
	{
		vectIndPgnPoints.push_back(iter->mInt);
	}
	listIndAngle.erase(listIndAngle.begin(), listIndAngle.end());
}

//-------------------------------------------------------------------------

void radTPolyhedron::CollectAndMapUniquePolygonPoints(const radTHandlePgnAndTrans& hPgnAndTrf, vector<TVector3d>& vectPoints, vector<int>& vectInd, double AbsTol)
{
	radTPolygon *pPgn = hPgnAndTrf.PgnHndl.rep;
	radTVect2dVect &vPointsPgn = pPgn->EdgePointsVector;
	double zc = pPgn->CoordZ;
	radTrans *pTrf = hPgnAndTrf.TransHndl.rep;
	int numVertexPoints = (int)vPointsPgn.size();
	double AbsTolE2 = AbsTol*AbsTol;

	for(int i=0; i<numVertexPoints; i++)
	{
		TVector2d &vP2d = vPointsPgn[i];
		TVector3d vP(vP2d.x, vP2d.y, zc);
		vP = pTrf->TrPoint(vP);

		int curSizeVectPoints = (int)vectPoints.size();
		int indPointFound = -1;
		for(int j=0; j<curSizeVectPoints; j++)
		{
			TVector3d dP = vP - vectPoints[j];
			if(dP.AmpE2() <= AbsTolE2)
			{
				indPointFound = j; break;
			}
		}
		if(indPointFound < 0)
		{
			vectPoints.push_back(vP);
			indPointFound = curSizeVectPoints;
		}
		vectInd.push_back(indPointFound);
	}
}

//-------------------------------------------------------------------------

bool radTPolyhedron::CheckIfAllPolygonVerticesAreOnOneSideOfPlane(const radTHandlePgnAndTrans& hPgnAndTrf, const TVector3d& vPoint, const TVector3d& vNorm, double AbsTol)
{//vNorm should be unit vector
	radTPolygon *pPgn = hPgnAndTrf.PgnHndl.rep;
	radTVect2dVect &vPointsPgn = pPgn->EdgePointsVector;
	double zc = pPgn->CoordZ;
	radTrans *pTrf = hPgnAndTrf.TransHndl.rep;
	int numVertexPoints = (int)vPointsPgn.size();

	int signPrevScalProd = 0;
	for(int i=0; i<numVertexPoints; i++)
	{
		TVector2d &vP2d = vPointsPgn[i];
		TVector3d vP(vP2d.x, vP2d.y, zc);
		TVector3d vR = pTrf->TrPoint(vP) - vPoint;
		double testScalProd = vR*vNorm;
		int signTestScalProd = 0;
		if(testScalProd < -AbsTol) signTestScalProd = -1;
		else if(testScalProd > AbsTol) signTestScalProd = 1;

		//if(signPrevScalProd*testScalProd < 0) return false;
		if(signPrevScalProd*signTestScalProd < 0) return false; //OC140209
		if(signTestScalProd != 0) signPrevScalProd = signTestScalProd;
	}
	return true;
}

//-------------------------------------------------------------------------

void radTPolyhedron::SetCurrentDensityForConstCurrent(double avgCur, int indBaseFace1, int indBaseFace2)
{//assumes that the member:
 //radTVectHandlePgnAndTrans VectHandlePgnAndTrans;
 //was already set up.
	if(avgCur == 0) return;
	int numFaces = (int)VectHandlePgnAndTrans.size();
	if(numFaces <= 0) { SomethingIsWrong = 1; return;}
	if((indBaseFace1 < 0) || (indBaseFace1 >= numFaces)) { SomethingIsWrong = 1; return;}
	if((indBaseFace2 < 0) || (indBaseFace2 >= numFaces)) { SomethingIsWrong = 1; return;}

	TVector3d Ez(0,0,1);

	radTHandlePgnAndTrans hPgnAndTrans1 = VectHandlePgnAndTrans[indBaseFace1];
	radTPolygon* pPgn1 = hPgnAndTrans1.PgnHndl.rep;
	radTrans* pTrans1 = hPgnAndTrans1.TransHndl.rep;
	TVector3d vNorm1 = (pTrans1->TrBiPoint(Ez));
	TVector3d cenP1(pPgn1->CentrPoint.x, pPgn1->CentrPoint.y, pPgn1->CoordZ);
	cenP1 = pTrans1->TrPoint(cenP1);
	double S1 = pPgn1->Area();

	radTHandlePgnAndTrans hPgnAndTrans2 = VectHandlePgnAndTrans[indBaseFace2];
	radTPolygon* pPgn2 = hPgnAndTrans2.PgnHndl.rep;
	radTrans* pTrans2 = hPgnAndTrans2.TransHndl.rep;
	TVector3d vNorm2 = (pTrans2->TrBiPoint(Ez));
	TVector3d cenP2(pPgn2->CentrPoint.x, pPgn2->CentrPoint.y, pPgn2->CoordZ);
	cenP2 = pTrans2->TrPoint(cenP2);
	double S2 = pPgn2->Area();

	TVector3d vCurAxis = cenP2 - cenP1;
	double height = vCurAxis.Abs();
	if(height == 0)
	{
		SomethingIsWrong = 1;
		radTSend::ErrorMessage("Radia::Error122"); return;
	}
	double invHeight = 1./height;
	vCurAxis *= invHeight;

	double vNorm1_vCurAxis = vNorm1*vCurAxis, vNorm2_vCurAxis = vNorm2*vCurAxis;

	if((S1 <= 0) || (S2 <= 0) || (vNorm1_vCurAxis == 0) || (vNorm2_vCurAxis == 0))
	{
		SomethingIsWrong = 1;
		radTSend::ErrorMessage("Radia::Error122"); return;
	}

	if(vNorm1_vCurAxis < 0) { vNorm1_vCurAxis = -vNorm1_vCurAxis; vNorm1 *= (-1);}
	if(vNorm2_vCurAxis < 0) { vNorm2_vCurAxis = -vNorm2_vCurAxis; vNorm2 *= (-1);}

	double J1 = avgCur/(vNorm1_vCurAxis*S1);
	double J2 = avgCur/(vNorm2_vCurAxis*S2);
	double Js = (J2 - J1)*invHeight;

	J = J1*vCurAxis;

	if(Js != 0.)
	{//Matrix of linear coefficients for J is necessary
		
		J += ((-Js)*(cenP1*vCurAxis))*vCurAxis;
		
		if(pJ_LinCoef == 0) pJ_LinCoef = new TMatrix3d();
		
		double vxvy = Js*vCurAxis.x*vCurAxis.y, vxvz = Js*vCurAxis.x*vCurAxis.z, vyvz = Js*vCurAxis.y*vCurAxis.z;
		pJ_LinCoef->Str0.x = Js*vCurAxis.x*vCurAxis.x; pJ_LinCoef->Str0.y = vxvy; pJ_LinCoef->Str0.z = vxvz;
		pJ_LinCoef->Str1.x = vxvy; pJ_LinCoef->Str1.y = Js*vCurAxis.y*vCurAxis.y; pJ_LinCoef->Str1.z = vyvz;
		pJ_LinCoef->Str2.x = vxvz; pJ_LinCoef->Str2.y = vyvz; pJ_LinCoef->Str2.z = Js*vCurAxis.z*vCurAxis.z;

		J += (*pJ_LinCoef)*CentrPoint;
		J_IsNotZero = true;
		mLinTreat = 0; //treat J and lin. coef. as Relative
	}

	if(!J.isZero()) J_IsNotZero = true;
}

//-------------------------------------------------------------------------
