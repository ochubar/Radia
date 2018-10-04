/*-------------------------------------------------------------------------
*
* File name:      radapl1.cpp
*
* Project:        RADIA
*
* Description:    Wrapping RADIA application function calls
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radappl.h"
#include "radg3da1.h"
#include "radrec.h"
#include "radarccu.h"
#include "radplnr.h"
#include "radexpgn.h"
#include "radvlpgn.h"
#include "radcnvrg.h"
#include "radopnam.h"

#include <math.h>
#include <string.h>

//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------

int radTApplication::ValidateVector3d(double* ArrayToCheck, long LenArray, TVector3d* VectorPtr)
{
	if(ArrayToCheck == 0) return 1; //?OC20082008
	if(LenArray==3)
	{
		*VectorPtr = TVector3d(ArrayToCheck);
		return 1;
	}
	Send.ErrorMessage("Radia::Error000"); return 0;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateVector2d(double* ArrayToCheck, long LenArray, TVector2d* VectorPtr)
{
	if(LenArray==2)
	{
		TVector2d aVect2d; 
		aVect2d.x = ArrayToCheck[0]; aVect2d.y = ArrayToCheck[1]; 
		*VectorPtr = aVect2d;
		return 1;
	}
	Send.ErrorMessage("Radia::Error000"); return 0;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateMatrix3d(double* arToCheck, long LenAr, TMatrix3d* MatrixPtr)
{
	if(arToCheck == 0) return 1; //?OC20082008
	if(LenAr==9)
	{
		*MatrixPtr = TMatrix3d(arToCheck);
		return 1;
	}
	Send.ErrorMessage("Radia::Error000"); return 0;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateElemKey(long ElemKey, radThg& hg)
{
	radTmhg::const_iterator iter = GlobalMapOfHandlers.find(ElemKey);
	if(iter == GlobalMapOfHandlers.end())
	{
		Send.ErrorMessage("Radia::Error002");
		return 0;
	}
	hg = (*iter).second;
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateFieldChar(char* FieldChar, radTFieldKey* FieldKeyPtr, bool LocSendRequired)
{
	char* BufChar = FieldChar;
	if(*FieldChar == '\0')
	{
		FieldKeyPtr->B_=1; FieldKeyPtr->H_ = 1; FieldKeyPtr->A_ = 1; FieldKeyPtr->M_ = 1; FieldKeyPtr->J_ = 1; FieldKeyPtr->Phi_ = 1;
	}
	else
	{
		while(*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b')) FieldKeyPtr->B_ = 1;
			if((*BufChar == 'H') || (*BufChar == 'h')) FieldKeyPtr->H_ = 1;
			if((*BufChar == 'A') || (*BufChar == 'a')) FieldKeyPtr->A_ = 1;
			if((*BufChar == 'M') || (*BufChar == 'm')) FieldKeyPtr->M_ = 1;
			if((*BufChar == 'J') || (*BufChar == 'j')) FieldKeyPtr->J_ = 1;
			if((*BufChar == 'P') || (*BufChar == 'p')) FieldKeyPtr->Phi_ = 1;
			
			if((*BufChar == 'Q') || (*BufChar == 'q')) //OC161005
			{
				FieldKeyPtr->B_ = FieldKeyPtr->H_ = FieldKeyPtr->PreRelax_ = FieldKeyPtr->Q_ = 1;
			}
			
			BufChar++;
		}
	}
	if((!FieldKeyPtr->B_) && (!FieldKeyPtr->H_) && (!FieldKeyPtr->A_) && (!FieldKeyPtr->M_) && (!FieldKeyPtr->J_) && (!FieldKeyPtr->Phi_))
	{
		if(LocSendRequired) Send.ErrorMessage("Radia::Error009"); 
		return 0;
	}
	else return 1;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateFieldIntChar(char* FieldIntChar, char* FinOrInfChar, radTFieldKey* FieldKeyPtr, bool LocSendRequired)
{
	char* BufChar = FieldIntChar;
	short I_used = 0;
	if(*FieldIntChar == '\0') {	FieldKeyPtr->Ib_=1; FieldKeyPtr->FinInt_=0;}
	else
	{
		while(*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b')) FieldKeyPtr->Ib_ = 1;
			if((*BufChar == 'H') || (*BufChar == 'h')) FieldKeyPtr->Ih_ = 1;
			if((*BufChar == 'I') || (*BufChar == 'i')) I_used = 1;
			BufChar++;
		}
	}
	if((!FieldKeyPtr->Ib_) && (!FieldKeyPtr->Ih_))
	{
		if(I_used) { FieldKeyPtr->Ib_=1; FieldKeyPtr->FinInt_=0;}
		else 
		{ 
			if(LocSendRequired) Send.ErrorMessage("Radia::Error032"); 
			return 0;
		}
	}

	if((!strcmp(FinOrInfChar, "inf")) || (!strcmp(FinOrInfChar, "Inf")) || (!strcmp(FinOrInfChar, "INF"))) 
		FieldKeyPtr->FinInt_=0;
	else if((!strcmp(FinOrInfChar, "fin")) || (!strcmp(FinOrInfChar, "Fin")) || (!strcmp(FinOrInfChar, "FIN"))) 
		FieldKeyPtr->FinInt_=1;
	else 
	{ 
		if(LocSendRequired) Send.ErrorMessage("Radia::Error033"); 
		return 0;
	}

	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateFieldEnergyForceChar(char* ComponIDChar, radTFieldKey* FieldKeyPtr)
{
	char* BufChar = ComponIDChar;
	if(*ComponIDChar == '\0') {	FieldKeyPtr->Energy_=1; FieldKeyPtr->ForceEnr_=1; FieldKeyPtr->Torque_=1;}
	else
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'E') || (*BufChar == 'e')) FieldKeyPtr->Energy_ = 1;
			if((*BufChar == 'F') || (*BufChar == 'f')) FieldKeyPtr->ForceEnr_ = 1;
			if((*BufChar == 'T') || (*BufChar == 't')) FieldKeyPtr->Torque_ = 1;
			BufChar++;
		}
	if((!FieldKeyPtr->Energy_) && (!FieldKeyPtr->ForceEnr_) && (!FieldKeyPtr->Torque_))
	{
		Send.ErrorMessage("Radia::Error056"); return 0;
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateMagnChar(char* MagnChar)
{
	char* BufChar = MagnChar;
	if(*BufChar == '\0') { return 1;}
	
	short MagnCharIsValid = 0;
	while (*BufChar != '\0') 
	{
		if((*BufChar == 'M') || (*BufChar == 'm') ||
		   (*BufChar == 'X') || (*BufChar == 'x') ||
		   (*BufChar == 'Y') || (*BufChar == 'y') ||
		   (*BufChar == 'Z') || (*BufChar == 'z')) MagnCharIsValid = 1;
		BufChar++;
	}

	if(!MagnCharIsValid)
	{
		Send.ErrorMessage("Radia::Error026"); return 0;
	}
	else return 1;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateForceChar(char* ForceChar)
{
	char* BufChar = ForceChar;
	if(*BufChar == '\0') { return 1;}
	
	short CharIsValid = 0;
	while (*BufChar != '\0') 
	{
		//if((*BufChar == 'M') || (*BufChar == 'm') ||
		if((*BufChar == 'F') || (*BufChar == 'f') || //OC061008
		   (*BufChar == 'X') || (*BufChar == 'x') ||
		   (*BufChar == 'Y') || (*BufChar == 'y') ||
		   (*BufChar == 'Z') || (*BufChar == 'z')) CharIsValid = 1;
		BufChar++;
	}
	if(!CharIsValid) { Send.ErrorMessage("Radia::Error058"); return 0;}
	else return 1;
}

//-------------------------------------------------------------------------

int radTApplication::ValidateTorqueChar(char* TorqueChar)
{
	char* BufChar = TorqueChar;
	if(*BufChar == '\0') { return 1;}
	
	short CharIsValid = 0;
	while (*BufChar != '\0') 
	{
		if((*BufChar == 'T') || (*BufChar == 't') ||
		   (*BufChar == 'X') || (*BufChar == 'x') ||
		   (*BufChar == 'Y') || (*BufChar == 'y') ||
		   (*BufChar == 'Z') || (*BufChar == 'z')) CharIsValid = 1;
		BufChar++;
	}
	if(!CharIsValid) { Send.ErrorMessage("Radia::Error059"); return 0;}
	else return 1;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTApplication::SetRecMag(double* CPoi, long lenCPoi, double* Dims, long lenDims, double* Magn, long lenMagn, double* J, long lenJ, short J_IsNotZero)
{
	try
	{
		TVector3d CPoiVect, DimsVect, MagnVect, J_vect;
		if(!ValidateVector3d(CPoi, lenCPoi, &CPoiVect)) return 0;
		if(!ValidateVector3d(Dims, lenDims, &DimsVect)) return 0;
		if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return 0;
		if(!ValidateVector3d(J, lenJ, &J_vect)) return 0;
		// May be not necessary?
		if((Dims[0]<0.) || (Dims[1]<0.) || (Dims[2]<0.))
		{
			Send.ErrorMessage("Radia::Error001"); return 0;
		}

		int ElemKey;
		const double ZeroTolCurrDens = 1.E-10;

		radTRecMag* RecMagPtr = 0;
		radTExtrPolygon* ExtrPolygonPtr = 0;
		radThg hg;

		if((!TreatRecMagsAsExtrPolygons) || ((fabs(J_vect.x)>ZeroTolCurrDens) || (fabs(J_vect.y)>ZeroTolCurrDens) || (fabs(J_vect.z)>ZeroTolCurrDens)))
		{
			RecMagPtr = new radTRecMag(CPoiVect, DimsVect, MagnVect, J_vect, J_IsNotZero);
			if(RecMagPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			hg = radThg(RecMagPtr);
		}
		else
		{
			TVector3d FirstPoiVect = CPoiVect - (0.5*DimsVect);
			TVector2d ArrayOfPoints2d[4];
			ArrayOfPoints2d[0] = TVector2d(FirstPoiVect.y, FirstPoiVect.z);
			ArrayOfPoints2d[1] = TVector2d(FirstPoiVect.y + DimsVect.y, FirstPoiVect.z);
			ArrayOfPoints2d[2] = TVector2d(FirstPoiVect.y + DimsVect.y, FirstPoiVect.z + DimsVect.z);
			ArrayOfPoints2d[3] = TVector2d(FirstPoiVect.y, FirstPoiVect.z + DimsVect.z);

			ExtrPolygonPtr = new radTExtrPolygon(FirstPoiVect, ParallelToX, DimsVect.x, ArrayOfPoints2d, 4, MagnVect);
			if(ExtrPolygonPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			hg = radThg(ExtrPolygonPtr);

			if(TreatExtrPgnsAsPolyhedrons) if(!ExtrPolygonPtr->ConvertToPolyhedron(hg, this, 1)) return 0;
		}

		if(TreatRecMagsAsPolyhedrons) if(!RecMagPtr->ConvertToPolyhedron(hg, this, 1)) return 0;

		ElemKey = AddElementToContainer(hg);

		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetExtrudedPolygon(double* FirstPoi, long lenFirstPoi, double Lx, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* InMagn, long lenMagn, const char* OrientStr)
{
	try
	{
		//TVector3d FirstPoiVect, MagnVect; //OC090106
		//if(!ValidateVector3d(FirstPoi, lenFirstPoi, &FirstPoiVect)) return 0;
		//if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return 0;

		if(lenMagn != 3) { Send.ErrorMessage("Radia::Error000"); return 0;}
		double Magn[] = { InMagn[0], InMagn[1], InMagn[2]}; //OC180210 to prevent InMagn from modification on output

		if(CheckIfExtrudedPolygonIsRecMag(ArrayOfPoints2d, lenArrayOfPoints2d)) //OC190210
		{//code from void RecMag(double xc, double yc, double zc, double Lx, double Ly, double Lz, double Mx, double My, double Mz) in radinter.cpp
			double J[] = {0.,0.,0.};
			double CP[] = {0.,0.,0.}, Dims[] = {0.,0.,0.};
			TVector2d &FirstPoint2d = ArrayOfPoints2d[0];

			char a = *OrientStr;
			if((a == 'x') || (a == 'X'))
			{
				CP[0] = FirstPoi[0] + 0.5*Lx; //FirstPoi[0] = xc - 0.5*Lx;
				//CP[1] = FirstPoint2d.x; //FirstPoi[1] = FirstPoint2d.x;
				//CP[2] = FirstPoint2d.y; //FirstPoi[2] = FirstPoint2d.y;
				Dims[0] = Lx;

				double dim1a = fabs(ArrayOfPoints2d[1].x - ArrayOfPoints2d[0].x);
				double dim1b = fabs(ArrayOfPoints2d[2].x - ArrayOfPoints2d[1].x);
				//Dims[1] = (dim1a > dim1b)? dim1a : dim1b;
				if(dim1a > dim1b)
				{
					CP[1] = 0.5*(ArrayOfPoints2d[1].x + ArrayOfPoints2d[0].x);
					Dims[1] = dim1a;
				}
				else
				{
					CP[1] = 0.5*(ArrayOfPoints2d[2].x + ArrayOfPoints2d[1].x);
					Dims[1] = dim1b;
				}

				double dim2a = fabs(ArrayOfPoints2d[1].y - ArrayOfPoints2d[0].y);
				double dim2b = fabs(ArrayOfPoints2d[2].y - ArrayOfPoints2d[1].y);
				//Dims[2] = (dim2a > dim2b)? dim2a : dim2b;
				if(dim2a > dim2b)
				{
					CP[2] = 0.5*(ArrayOfPoints2d[1].y + ArrayOfPoints2d[0].y);
					Dims[2] = dim2a;
				}
				else
				{
					CP[2] = 0.5*(ArrayOfPoints2d[2].y + ArrayOfPoints2d[1].y);
					Dims[2] = dim2b;
				}
			}
			else if((a == 'y') || (a == 'Y'))
			{
				//CP[0] = FirstPoint2d.y; //FirstPoi[0] = FirstPoint2d.y;
				CP[1] = FirstPoi[1] + 0.5*Lx; //FirstPoi[1] = xc - 0.5*Lx;
				//CP[2] = FirstPoint2d.x; //FirstPoi[2] = FirstPoint2d.x;
				Dims[1] = Lx;

				double dim0a = fabs(ArrayOfPoints2d[1].y - ArrayOfPoints2d[0].y);
				double dim0b = fabs(ArrayOfPoints2d[2].y - ArrayOfPoints2d[1].y);
				//Dims[0] = (dim0a > dim0b)? dim0a : dim0b;
				if(dim0a > dim0b)
				{
					CP[0] = 0.5*(ArrayOfPoints2d[1].y + ArrayOfPoints2d[0].y);
					Dims[0] = dim0a;
				}
				else
				{
					CP[0] = 0.5*(ArrayOfPoints2d[2].y + ArrayOfPoints2d[1].y);
					Dims[0] = dim0b;
				}

				double dim2a = fabs(ArrayOfPoints2d[1].x - ArrayOfPoints2d[0].x);
				double dim2b = fabs(ArrayOfPoints2d[2].x - ArrayOfPoints2d[1].x);
				//Dims[2] = (dim2a > dim2b)? dim2a : dim2b;
				if(dim2a > dim2b)
				{
					CP[2] = 0.5*(ArrayOfPoints2d[1].x + ArrayOfPoints2d[0].x);
					Dims[2] = dim2a;
				}
				else
				{
					CP[2] = 0.5*(ArrayOfPoints2d[2].x + ArrayOfPoints2d[1].x);
					Dims[2] = dim2b;
				}
			}
			else
			{
				//CP[0] = FirstPoint2d.x; //FirstPoi[0] = FirstPoint2d.x;
				//CP[1] = FirstPoint2d.y; //FirstPoi[1] = FirstPoint2d.y;
				CP[2] = FirstPoi[2] + 0.5*Lx; //FirstPoi[2] = xc - 0.5*Lx;
				Dims[2] = Lx;

				double dim0a = fabs(ArrayOfPoints2d[1].x - ArrayOfPoints2d[0].x);
				double dim0b = fabs(ArrayOfPoints2d[2].x - ArrayOfPoints2d[1].x);
				//Dims[0] = (dim0a > dim0b)? dim0a : dim0b;
				if(dim0a > dim0b)
				{
					CP[0] = 0.5*(ArrayOfPoints2d[1].x + ArrayOfPoints2d[0].x);
					Dims[0] = dim0a;
				}
				else
				{
					CP[0] = 0.5*(ArrayOfPoints2d[2].x + ArrayOfPoints2d[1].x);
					Dims[0] = dim0b;
				}
				
				double dim1a = fabs(ArrayOfPoints2d[1].y - ArrayOfPoints2d[0].y);
				double dim1b = fabs(ArrayOfPoints2d[2].y - ArrayOfPoints2d[1].y);
				//Dims[1] = (dim1a > dim1b)? dim1a : dim1b;
				if(dim1a > dim1b)
				{
					CP[1] = 0.5*(ArrayOfPoints2d[1].y + ArrayOfPoints2d[0].y);
					Dims[1] = dim1a;
				}
				else
				{
					CP[1] = 0.5*(ArrayOfPoints2d[2].y + ArrayOfPoints2d[1].y);
					Dims[1] = dim1b;
				}
			}

			short OldActOnDoubles = radCR.ActOnDoubles;
			if(TreatRecMagsAsExtrPolygons) radCR.ActOnDoubles = 0;

			double *tDims = Dims;
			for(int ii=0; ii<3; ii++) { *tDims = radCR.Double(*tDims); tDims++;}

			int outRes = SetRecMag(CP, 3, Dims, 3, Magn, 3, J, 3, 0);
			if(TreatRecMagsAsExtrPolygons) radCR.ActOnDoubles = OldActOnDoubles;
			return outRes;
		}

		short PrevSendingIsRequired = SendingIsRequired; SendingIsRequired = 0; //OC301207
		//double CPoi[] = {FirstPoi[0] + 0.5*Lx, 0, 0};
		double CPoi[] = {0, 0, 0};

		int IndTr = FindSpaceTransToOrientObjAlongMainAxis(CPoi, 'X', *OrientStr);
		if((Magn != 0) && (lenMagn >= 3) && (IndTr != 0)) TransformBackMagnOrCurDensArr(IndTr, Magn, lenMagn);
		
		TransformBackPointArr(IndTr, FirstPoi, lenFirstPoi); //OC090106
		
		TVector3d FirstPoiVect, MagnVect; //OC090106
		if(!ValidateVector3d(FirstPoi, lenFirstPoi, &FirstPoiVect)) return 0;
		if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return 0;

		TAxisOrient AxOrnt = ParallelToX;

		radTExtrPolygon* ExtrPgnPtr = new radTExtrPolygon(FirstPoiVect, AxOrnt, Lx, ArrayOfPoints2d, int(lenArrayOfPoints2d), MagnVect);
		if(ExtrPgnPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		if(((radTPolygon*)(ExtrPgnPtr->BasePolygonHandle.rep))->SomethingIsWrong)
		{
			SendingIsRequired = PrevSendingIsRequired;
			delete ExtrPgnPtr; return 0;
		}
		else
		{
			radThg hg(ExtrPgnPtr);
			if(TreatExtrPgnsAsPolyhedrons) if(!ExtrPgnPtr->ConvertToPolyhedron(hg, this, 1)) return 0;

			int ElemKey = AddElementToContainer(hg);

            //short PrevSendingIsRequired = SendingIsRequired; SendingIsRequired = 0; //OC301207
            if(IndTr != 0) ElemKey = ApplySymmetry(ElemKey, IndTr, 1);
			SendingIsRequired = PrevSendingIsRequired;

			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetPlanarPolygon(double CoordZ, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* Magn, long lenMagn)
{
	TVector3d MagnVect;
	if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return 0;

	try
	{
		radTPolygon* PgnPtr = new radTPolygon(CoordZ, ArrayOfPoints2d, int(lenArrayOfPoints2d), MagnVect);

		if(PgnPtr->SomethingIsWrong)
		{
			delete PgnPtr; return 0;
		}
		else
		{
			radThg hg(PgnPtr);
			int ElemKey = AddElementToContainer(hg);
			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetPolyhedron1(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfNumOfPoInFaces, int lenArrayOfFaces, double* Magn, double* arM_LinCoef, double* arJ, double* arJ_LinCoef, const char** OptionNames, const char** OptionValues, int OptionCount)
//int radTApplication::SetPolyhedron1(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfNumOfPoInFaces, int lenArrayOfFaces, double* Magn, long lenMagn)
{
	TVector3d MagnVect, J_Vect;
	if(!ValidateVector3d(Magn, 3, &MagnVect)) return 0;
	if(!ValidateVector3d(arJ, 3, &J_Vect)) return 0;

	TMatrix3d M_LinCoefMatr, J_LinCoefMatr;
	if(!ValidateMatrix3d(arM_LinCoef, 9, &M_LinCoefMatr)) return 0; //arM_LinCoef is expected to contain matrix elements string by string 
	if(!ValidateMatrix3d(arJ_LinCoef, 9, &J_LinCoefMatr)) return 0; //arJ_LinCoef is expected to contain matrix elements string by string 

	try
	{
		radTOptionNames OptNam;
		const char* OptNamesToFind[] = {OptNam.LinTreat};
		char OptValsFoundParsed[] = {0};
		char &LinTreat = OptValsFoundParsed[0]; // 0- Relative to center of the polyhedron; 1- Absolute;

		if(!OptNam.findParseOptionValues(OptionNames, OptionValues, OptionCount, OptNamesToFind, 1, OptValsFoundParsed, 0, 0))
		{
			Send.ErrorMessage("Radia::Error062"); return 0;
		}

		//const char** BufNameString = OptionNames;
		//const char** BufValString = OptionValues;
		//for(int i=0; i<OptionCount; i++)
		//{
		//	if(!strcmp(*BufNameString, OptNam.LinCoef))
		//	{
		//		if(!strcmp(*BufValString, (OptNam.LinCoefValues)[0])) LinCoefTreat = 0;
		//		else if(!strcmp(*BufValString, (OptNam.LinCoefValues)[1])) LinCoefTreat = 0;
		//		else if(!strcmp(*BufValString, (OptNam.LinCoefValues)[2])) LinCoefTreat = 1;
		//		else if(!strcmp(*BufValString, (OptNam.LinCoefValues)[3])) LinCoefTreat = 1;
		//		else { Send.ErrorMessage("Radia::Error062"); return 0;}
		//	}
		//	else { Send.ErrorMessage("Radia::Error062"); return 0;}
		//	BufNameString++; BufValString++;
		//}

		//radTPolyhedron* VolLimByPgnsPtr = new radTPolyhedron(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfNumOfPoInFaces, lenArrayOfFaces, MagnVect);
		radTPolyhedron* VolLimByPgnsPtr = new radTPolyhedron(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfNumOfPoInFaces, lenArrayOfFaces, MagnVect, M_LinCoefMatr, J_Vect, J_LinCoefMatr, LinTreat);
		if(VolLimByPgnsPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		if(VolLimByPgnsPtr->SomethingIsWrong) { delete VolLimByPgnsPtr; return 0;}
		else
		{
			radThg hg(VolLimByPgnsPtr);
			if(RecognizeRecMagsInPolyhedrons)
			{
				double RelAbsTol[] = { radCR.AbsRand, radCR.RelRand };
				VolLimByPgnsPtr->CheckForSpecialShapes(VolLimByPgnsPtr->VectHandlePgnAndTrans, hg, RelAbsTol);
			}
			int ElemKey = AddElementToContainer(hg);
			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetPolyhedron2(TVector3d** ArrayOfFaces, int* ArrayOfNumOfPoInFaces, long lenArrayOfFaces, double* Magn, long lenMagn)
{
	TVector3d MagnVect;
	if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return 0;

	try
	{
		radTPolyhedron* VolLimByPgnsPtr = new radTPolyhedron(ArrayOfFaces, ArrayOfNumOfPoInFaces, lenArrayOfFaces, MagnVect);
		if(VolLimByPgnsPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		if(VolLimByPgnsPtr->SomethingIsWrong) { delete VolLimByPgnsPtr; return 0;}
		else
		{
			radThg hg(VolLimByPgnsPtr);
			if(RecognizeRecMagsInPolyhedrons)
			{
				double RelAbsTol[] = { radCR.AbsRand, radCR.RelRand };
				VolLimByPgnsPtr->CheckForSpecialShapes(VolLimByPgnsPtr->VectHandlePgnAndTrans, hg, RelAbsTol);
			}
			int ElemKey = AddElementToContainer(hg);
			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::RecMagsAsExtrPolygons(char* OnOrOff)
{
	try
	{
		if((!strcmp(OnOrOff, "on")) || (!strcmp(OnOrOff, "On")) || (!strcmp(OnOrOff, "ON")))
		{
			TreatRecMagsAsExtrPolygons = 1; TreatRecMagsAsPolyhedrons = 0;
		}
		else if((!strcmp(OnOrOff, "off")) || (!strcmp(OnOrOff, "Off")) || (!strcmp(OnOrOff, "OFF")))
			TreatRecMagsAsExtrPolygons = 0;
		else { Send.ErrorMessage("Radia::Error043"); return 0;}
		if(SendingIsRequired) Send.Int(int(TreatRecMagsAsExtrPolygons));
		return 1;
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::RecMagsAsPolyhedrons(char* OnOrOff)
{
	try
	{
		if((!strcmp(OnOrOff, "on")) || (!strcmp(OnOrOff, "On")) || (!strcmp(OnOrOff, "ON")))
		{
			TreatRecMagsAsPolyhedrons = 1; TreatRecMagsAsExtrPolygons = 0;
		}
		else if((!strcmp(OnOrOff, "off")) || (!strcmp(OnOrOff, "Off")) || (!strcmp(OnOrOff, "OFF")))
			TreatRecMagsAsPolyhedrons = 0;
		else { Send.ErrorMessage("Radia::Error043"); return 0;}
		if(SendingIsRequired) Send.Int(int(TreatRecMagsAsPolyhedrons));
		return 1;
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::RecognizeRecMags(char* OnOrOff)
{
	try
	{
		if((!strcmp(OnOrOff, "on")) || (!strcmp(OnOrOff, "On")) || (!strcmp(OnOrOff, "ON")))
			RecognizeRecMagsInPolyhedrons = 1;
		else if((!strcmp(OnOrOff, "off")) || (!strcmp(OnOrOff, "Off")) || (!strcmp(OnOrOff, "OFF")))
			RecognizeRecMagsInPolyhedrons = 0;
		else { Send.ErrorMessage("Radia::Error043"); return 0;}
		if(SendingIsRequired) Send.Int(int(RecognizeRecMagsInPolyhedrons));
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ExtPgnsAsPolyhedrons(char* OnOrOff)
{
	try
	{
		if((!strcmp(OnOrOff, "on")) || (!strcmp(OnOrOff, "On")) || (!strcmp(OnOrOff, "ON")))
			TreatExtrPgnsAsPolyhedrons = 1;
		else if((!strcmp(OnOrOff, "off")) || (!strcmp(OnOrOff, "Off")) || (!strcmp(OnOrOff, "OFF")))
			TreatExtrPgnsAsPolyhedrons = 0;
		else { Send.ErrorMessage("Radia::Error043"); return 0;}
		if(SendingIsRequired) Send.Int(int(TreatExtrPgnsAsPolyhedrons));
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
} 

//-------------------------------------------------------------------------

//int radTApplication::SetMultGenExtrPolygonCur(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, const char** arOptionNames, const char** arOptionValues, int numOptions)
//int radTApplication::SetMultGenExtrPolygonCur(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInSteps, const char** arOptionNames, const char** arOptionValues, int numOptions)
int radTApplication::SetMultGenExtrPolygonCur(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double* arSubdData, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInSteps, const char** arOptionNames, const char** arOptionValues, int numOptions)
{
	try
	{
		if((arPoints2d == 0) || (lenArPoints2d <= 0)) return 0;
		//if((arPtrRotParInExtrSteps == 0) || (arNumRotInExtrSteps == 0)) return 0;
		//if((arTrslParInExtrSteps == 0) || (arHomParInExtrSteps == 0) || (NumSteps <= 0)) return 0;
		if((arPtrTrfParInExtrSteps == 0) || (arStrTrfOrderInExtrSteps == 0) || (arNumTrfInExtrSteps == 0) || (NumSteps <= 0)) return 0;

		//check types of transformations (e.g. one homothety is not allowed)
		for(int i=0; i<NumSteps; i++)
		{
			double **arPtrTrfParInCurExtrStep = arPtrTrfParInExtrSteps[i];
			char *arStrTrfOrderInCurExtrStep = arStrTrfOrderInExtrSteps[i];
			int numTrfInCurExtrStep = arNumTrfInExtrSteps[i];
			if((arPtrTrfParInCurExtrStep == 0) || (arStrTrfOrderInCurExtrStep == 0) || (numTrfInCurExtrStep <= 0))
			{
				Send.ErrorMessage("Radia::Error123"); return 0;
			}

			char *t_arStrTrfOrderInCurExtrStep = arStrTrfOrderInCurExtrStep;
			bool RotWasDefined = false, TranslWasDefined = false, HomWasDefined = false;
			for(int j=0; j<numTrfInCurExtrStep; j++)
			{
				double *pCurTrfPar = arPtrTrfParInCurExtrStep[j];
				if(pCurTrfPar == 0) { Send.ErrorMessage("Radia::Error123"); return 0;}

				if((*t_arStrTrfOrderInCurExtrStep == 'r') || (*t_arStrTrfOrderInCurExtrStep == 'R')) RotWasDefined = true;
				else if((*t_arStrTrfOrderInCurExtrStep == 't') || (*t_arStrTrfOrderInCurExtrStep == 'T')) TranslWasDefined = true;
				else if((*t_arStrTrfOrderInCurExtrStep == 'h') || (*t_arStrTrfOrderInCurExtrStep == 'H')) HomWasDefined = true;
				else { Send.ErrorMessage("Radia::Error123"); return 0;}
				t_arStrTrfOrderInCurExtrStep++;
			}
			if(HomWasDefined && (!RotWasDefined) && (!TranslWasDefined)) 
			{
				Send.ErrorMessage("Radia::Error124"); return 0;
			}
			//add more costraints here, if necessary
		}

		//treat option "Frame->Loc|Lab":
		//radTOptionNames OptNam;
		//const char* OptNamesToFind[] = {OptNam.Frame};
		//char OptValsFoundParsed[] = {0};
		//char &testFrame = OptValsFoundParsed[0]; // 0- Local; 1,2- Laboratory;
		//if(!OptNam.findParseOptionValues(arOptionNames, arOptionValues, numOptions, OptNamesToFind, 1, OptValsFoundParsed, 0, 0))
		//{
		//	Send.ErrorMessage("Radia::Error062"); return 0;
		//}

		//treat options Frame->Loc|Lab,ki->Numb|Size,TriAngMin->...,TriAreaMax->...,ExtOpt->\"...\"
		radTOptionNames OptNam;
		const char* OptNamesToFind[] = {OptNam.Frame, OptNam.SubdParamBorderCode, OptNam.TriAngMin, OptNam.TriAreaMax, OptNam.TriExtOpt};
		char cOptValsFoundParsed[] = {0, 0}; //Frame->Loc|Lab,ki->Numb|Size
		char &frame = cOptValsFoundParsed[0], &triSubdParamBorderCode = cOptValsFoundParsed[1];
		double dArOptValsFoundParsed[] = {0, 0}; //TriAngMin->...,TriAreaMax->...
		double &triAngMin = dArOptValsFoundParsed[0], &triAreaMax = dArOptValsFoundParsed[1];
		char sTriExtOpt[256]; 
		sTriExtOpt[0] = '\0';
		char *sArOptValsFoundParsed[] = {sTriExtOpt}; //TriExtOpt->\"...\"
		if(!OptNam.findParseOptionValues(arOptionNames, arOptionValues, numOptions, OptNamesToFind, 5, cOptValsFoundParsed, dArOptValsFoundParsed, sArOptValsFoundParsed))
		{
			Send.ErrorMessage("Radia::Error062"); return 0;
		}

/**
		char triSubdParamBorderCode=0; // 0- ki are subdiv. numbers; 1- ki are average sizes of pieces;
		double triAngMin=0, triAreaMax=0;
		const char *sTriExtOpt=0;
		char frame;

		radTOptionNames OptNames;
		const char** BufNameString = arOptionNames;
		const char** BufValString = arOptionValues;

		for(int i=0; i<numOptions; i++)
		{
			if(!strcmp(*BufNameString, OptNames.SubdParamBorderCode))
			{
				if(!strcmp(*BufValString, (OptNames.SubdParamCodeValues)[0])) triSubdParamBorderCode = 0;
				else if(!strcmp(*BufValString, (OptNames.SubdParamCodeValues)[1])) triSubdParamBorderCode = 1;
				else { Send.ErrorMessage("Radia::Error062"); throw 0;}
			}
			else if(!strcmp(*BufNameString, OptNames.TriAngMin))
			{
				triAngMin = atof(*BufValString);
			}
			else if(!strcmp(*BufNameString, OptNames.TriAreaMax))
			{
				triAreaMax = atof(*BufValString);
			}
			else if(!strcmp(*BufNameString, OptNames.TriExtOpt))
			{
				sTriExtOpt = *BufValString;
			}
			else if(!strcmp(*BufNameString, OptNames.Frame))
			{
				if(!strcmp(*BufValString, "Loc")) frame = 0;
				else if(!strcmp(*BufValString, "Lab")) frame = 1;
				else if(!strcmp(*BufValString, "LabTot")) frame = 2;
			}
			else { Send.ErrorMessage("Radia::Error062"); throw 0;}

			BufNameString++; BufValString++;
		}
**/
		if(triAngMin > 35.) { Send.ErrorMessage("Radia::Error098"); throw 0;} //keep the constraint coherent with the error message!

		TVector2d *arTriVertPt=0;
		int numTriVertPt=0, *arTriVertInd=0, numTri=0;
		if(arSubdData != 0)
		{
			//if(!TriangulatePolygon(arPoints2d, lenArPoints2d, arSubdData, arOptionNames, arOptionValues, numOptions, arTriVertPt, numTriVertPt, arTriVertInd, numTri)) return 0;
			if(!TriangulatePolygon(arPoints2d, lenArPoints2d, arSubdData, triSubdParamBorderCode, triAngMin, triAreaMax, sTriExtOpt, arTriVertPt, numTriVertPt, arTriVertInd, numTri)) return 0;
		}

		radThg hg;
		int resOK = 0;

		if((numTriVertPt > 2) && (numTri > 1))
			resOK = SetUpPolyhedronsFromBaseFacePolygonsTri(zc, strOrient, arPoints2d, lenArPoints2d, arTriVertPt, numTriVertPt, arTriVertInd, numTri, arPtrTrfParInExtrSteps, arStrTrfOrderInExtrSteps, arNumTrfInExtrSteps, NumSteps, avgCur, arMagnCompInSteps, frame, hg);
		else
			resOK = SetUpPolyhedronsFromBaseFacePolygons(zc, strOrient, arPoints2d, lenArPoints2d, arPtrTrfParInExtrSteps, arStrTrfOrderInExtrSteps, arNumTrfInExtrSteps, NumSteps, avgCur, arMagnCompInSteps, frame, hg);

		if(arTriVertPt != 0) delete[] arTriVertPt;
		if(arTriVertInd != 0) delete[] arTriVertInd;

		//if(!SetUpPolyhedronsFromBaseFacePolygons(zc, strOrient, arPoints2d, lenArPoints2d, arPtrTrfParInExtrSteps, arStrTrfOrderInExtrSteps, arNumTrfInExtrSteps, NumSteps, avgCur, frame, hg)) return 0;
		//if(!SetUpPolyhedronsFromBaseFacePolygons(zc, strOrient, arPoints2d, lenArPoints2d, arPtrTrfParInExtrSteps, arStrTrfOrderInExtrSteps, arNumTrfInExtrSteps, NumSteps, avgCur, arMagnCompInSteps, frame, hg)) return 0;
		if(!resOK) return 0;
		else
		{
			int ElemKey = AddElementToContainer(hg);
			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}

		//if(SendingIsRequired) Send.Int(0); //to remove
		//return 0;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetMultGenExtrPolygon(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons, double* Magn, long lenMagn)
{
	try
	{
		TVector3d MagnVect;
		if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return 0;
		if(!CheckLayerPolygonStructures(LayerPolygons, PtsNumbersInLayerPgns, CoordsZ, AmOfLayerPolygons)) return 0;

		radThg hg;
		if(!SetUpPolyhedronsFromLayerPolygons(LayerPolygons, PtsNumbersInLayerPgns, CoordsZ, AmOfLayerPolygons, MagnVect, hg)) return 0;
		else
		{
			int ElemKey = AddElementToContainer(hg);
			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetMultGenExtrRectangle(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect, double* Magn, long lenMagn)
{
	try
	{
		TVector3d MagnVect;
		if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return 0;
		if(!CheckLayerRectangleStructures(RectCenPoints, RectDims, AmOfLayerRect)) return 0;

		radThg hg;
		if(!SetUpPolyhedronsFromLayerRectangles(RectCenPoints, RectDims, AmOfLayerRect, MagnVect, hg)) return 0;
		else
		{
			int ElemKey = AddElementToContainer(hg);
			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetArcMag(double* CPoi, long lenCPoi, double* Radii, long lenRadii, double* Angles, long lenAngles, double h, int NumberOfSegm, double* Magn, long lenMagn, char* Orient)
{
	try
	{
		TVector3d CPoiVect;
		if(!ValidateVector3d(CPoi, lenCPoi, &CPoiVect)) return 0;
		// May be not necessary?
		if((lenRadii!=2) || (Radii[0]<0.) || (Radii[1]<0.) || (Radii[1]<Radii[0]))
		{
			Send.ErrorMessage("Radia::Error010"); return 0;
		}
		const double TwoPi = 2.*3.141592653589793238;
		if(Angles == 0)
		{
			Send.ErrorMessage("Radia::Error000"); return 0;
		}
		if((Angles[1] < Angles[0]) || ((Angles[1] - Angles[0]) > TwoPi))
		{
			Send.ErrorMessage("Radia::Error011"); return 0;
		}
		//if((lenAngles!=2) || (Angles[0]<0. || Angles[0]>TwoPi) || (Angles[0]<0. || Angles[0]>TwoPi) || (Angles[1]<Angles[0]))
		//{
		//	Send.ErrorMessage("Radia::Error011"); return 0;
		//}
		if(h <= 0)
		{
			Send.ErrorMessage("Radia::Error012"); return 0;
		}
		if(NumberOfSegm <= 0) { Send.ErrorMessage("Radia::Error090"); return 0;}

		double HalfH = 0.5*h;
		double AngSegm = (Angles[1] - Angles[0])/NumberOfSegm;
		double CurAng = Angles[0];

		double r1 = Radii[0], r2 = Radii[1];
		double cos_ang = cos(CurAng), sin_ang = sin(CurAng);
		double x11 = r1*cos_ang, y11 = r1*sin_ang;
		double x12 = r2*cos_ang, y12 = r2*sin_ang;

		TVector2d ArrayOfPoints2d[4];
		ArrayOfPoints2d->x = x11; ArrayOfPoints2d->y = y11;
		ArrayOfPoints2d[1].x = x12; ArrayOfPoints2d[1].y = y12;

		int* ArrayOfKeys = new int[NumberOfSegm];

		short PrevSendingIsRequired = SendingIsRequired; //OC301207
		SendingIsRequired = 0;

        int IndTr = FindSpaceTransToOrientObjAlongMainAxis(CPoi, 'X', *Orient);
        if((Magn != 0) && (lenMagn >= 3) && (IndTr != 0)) TransformBackMagnOrCurDensArr(IndTr, Magn, lenMagn);

		//short PrevSendingIsRequired = SendingIsRequired;
		SendingIsRequired = PrevSendingIsRequired; //OC301207

		for(int i=0; i<NumberOfSegm; i++)
		{
			CurAng += AngSegm;
			cos_ang = cos(CurAng); sin_ang = sin(CurAng);
			double x21 = r1*cos_ang, y21 = r1*sin_ang;
			double x22 = r2*cos_ang, y22 = r2*sin_ang;

			ArrayOfPoints2d[2].x = x22; ArrayOfPoints2d[2].y = y22;
			ArrayOfPoints2d[3].x = x21; ArrayOfPoints2d[3].y = y21;

			double FirstPoi[] = {CPoiVect.x - HalfH, CPoiVect.y + x11, CPoiVect.z + y11}; //assuming *Orient == 'x'
			if((*Orient == 'y') || (*Orient == 'Y')) //OC090106
			{
				FirstPoi[0] = CPoiVect.x + y11;
				FirstPoi[1] = CPoiVect.y - HalfH;
				FirstPoi[2] = CPoiVect.z + x11;
			}
			else if((*Orient == 'z') || (*Orient == 'Z'))
			{
				FirstPoi[0] = CPoiVect.x + x11;
				FirstPoi[1] = CPoiVect.y + y11;
				FirstPoi[2] = CPoiVect.z - HalfH;
			}

			SendingIsRequired = 0;
			ArrayOfKeys[i] = SetExtrudedPolygon(FirstPoi, 3, h, ArrayOfPoints2d, 4, Magn, lenMagn, "x");
			SendingIsRequired = PrevSendingIsRequired;
			if(ArrayOfKeys[i] <= 0) 
			{ 
				if(ArrayOfKeys != 0) { delete[] ArrayOfKeys; ArrayOfKeys = 0;}
				Send.Int(0); return 0;
			}

			x11 = x21; y11 = y21;
			x12 = x22, y12 = y22;
			ArrayOfPoints2d->x = x11; ArrayOfPoints2d->y = y11;
			ArrayOfPoints2d[1].x = x12; ArrayOfPoints2d[1].y = y12;
		}

		SendingIsRequired = 0;
		int GrpInd = SetGroup(ArrayOfKeys, NumberOfSegm);
		if(GrpInd <= 0) 
		{ 
			//if(SendingIsRequired) Send.Int(0); 
			Send.ErrorMessage("Radia::Error000"); return 0;
		}
		if(IndTr != 0) GrpInd = ApplySymmetry(GrpInd, IndTr, 1);
		SendingIsRequired = PrevSendingIsRequired;

		if(ArrayOfKeys != 0) { delete[] ArrayOfKeys; ArrayOfKeys = 0;}

		if(SendingIsRequired) Send.Int(GrpInd); 
		return GrpInd;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetArcPolygon(double* p2d, const char* OrientStr, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* Angles, int NumberOfSegm, const char* SymOrNoSymStr, double* Magn)
{
	TVector3d* VertexPointsArr = 0;
	int** ArrayOfFaces = 0;
	int AmOfFaces = 0, AmOfFacesOrig = 0;
	int* ArrayOfNumOfPoInFaces = 0;
	int* arIndBase2 = 0;
	try
	{
		TVector3d CPoiVect;
		//if(!ValidateVector3d(CenP, 3, &CPoiVect)) return 0;

		if(p2d == 0) { Send.ErrorMessage("Radia::Error000"); return 0;}
		if((OrientStr == 0) || ((*OrientStr != 'x') && (*OrientStr != 'X') && (*OrientStr != 'y') && (*OrientStr != 'Y') && (*OrientStr != 'z') && (*OrientStr != 'Z')))
		{
			Send.ErrorMessage("Radia::Error092"); return 0;
		}

		if((*OrientStr == 'x') || (*OrientStr == 'X'))
		{
			CPoiVect.x = 0.; CPoiVect.y = p2d[0]; CPoiVect.z = p2d[1]; 
		}
		else if((*OrientStr == 'y') || (*OrientStr == 'Y')) //OC081007_BNL
		{
			CPoiVect.x = p2d[1]; CPoiVect.y = 0; CPoiVect.z = p2d[0]; 
		}
		else if((*OrientStr == 'z') || (*OrientStr == 'Z'))
		{
			CPoiVect.x = p2d[0]; CPoiVect.y = p2d[1]; CPoiVect.z = 0; 
		}
	
		if((ArrayOfPoints2d == 0) || (lenArrayOfPoints2d <= 0)) { Send.ErrorMessage("Radia::Error000"); return 0;}

		const double Pi = 3.141592653589793238;
		const double TwoPi = 2.*Pi;
		if(Angles == 0) { Send.ErrorMessage("Radia::Error000"); return 0;}
		if((Angles[1] <= Angles[0]) || ((Angles[1] - Angles[0]) > TwoPi))
		{
			Send.ErrorMessage("Radia::Error011"); return 0;
		}
		if(NumberOfSegm <= 0) { Send.ErrorMessage("Radia::Error090"); return 0;}

		double SegmAngle = (Angles[1] - Angles[0])/NumberOfSegm;
		if((SegmAngle <= 0) || (SegmAngle >= Pi)) { Send.ErrorMessage("Radia::Error094"); return 0;}
		double AzimAngle1 = Angles[0];
		double AzimAngle2 = AzimAngle1 + SegmAngle;

		if((SymOrNoSymStr == 0) || ((strcmp(SymOrNoSymStr, "sym") != 0) && (strcmp(SymOrNoSymStr, "SYM") != 0) && (strcmp(SymOrNoSymStr, "nosym") != 0) && (strcmp(SymOrNoSymStr, "NOSYM") != 0)))
		{
			Send.ErrorMessage("Radia::Error093"); return 0;
		}
		bool FreeSym = true;
		if((strcmp(SymOrNoSymStr, "sym") == 0) || (strcmp(SymOrNoSymStr, "SYM") == 0)) FreeSym = false;

		TVector3d AxisVect(OrientStr);
		//double RotAng = 0;
		//if(NumberOfSegm > 1)
		//{
		//	RotAng = (Angles[1] - Angles[0])/(NumberOfSegm - 1);
		//}

		TVector3d AzAxVect0(0,0,0), RadVect1, RadVect2;
		if((*OrientStr == 'x') || (*OrientStr == 'X')) AzAxVect0.y = 1;
		else if((*OrientStr == 'y') || (*OrientStr == 'Y')) AzAxVect0.z = 1;
		else if((*OrientStr == 'z') || (*OrientStr == 'Z')) AzAxVect0.x = 1;

		gmTrans Rot1, Rot2;
		if(AzimAngle1 != 0) 
		{
			Rot1.SetupRotation(CPoiVect, AxisVect, AzimAngle1);
			RadVect1 = Rot1.TrBiPoint(AzAxVect0);
		}
		else RadVect1 = AzAxVect0;
		if(AzimAngle2 != 0) 
		{
			Rot2.SetupRotation(CPoiVect, AxisVect, AzimAngle2);
			RadVect2 = Rot2.TrBiPoint(AzAxVect0);
		}
		else RadVect2 = AzAxVect0;

		int AmOfVertexPoints = 2*lenArrayOfPoints2d;
		VertexPointsArr = new TVector3d[AmOfVertexPoints];

		bool RadiusOverlaps = false;
		bool RadiusIsZeroAtOnePoint = false, RadiusIsZeroAtTwoPoints = false;

		arIndBase2 = new int[lenArrayOfPoints2d];

		TVector2d *tPoints2d = ArrayOfPoints2d;
		TVector3d *tVertexPoints = VertexPointsArr, *tVertexPoints2 = VertexPointsArr + lenArrayOfPoints2d;
		int BasePointCount1 = 0, BasePointCount2 = 0;
		for(int i=0; i<lenArrayOfPoints2d; i++)
		{
			double r = tPoints2d->x, z = tPoints2d->y;

			//*tVertexPoints = r*RadVect1 + z*AxisVect;
			*tVertexPoints = CPoiVect + (r*RadVect1 + z*AxisVect); //OC081007_BNL
			tVertexPoints++; BasePointCount1++;

			if(r < 0.) RadiusOverlaps = true;
			else if(r == 0.) 
			{
				if(RadiusIsZeroAtTwoPoints) RadiusOverlaps = true;
				else if(RadiusIsZeroAtOnePoint) { RadiusIsZeroAtOnePoint = false; RadiusIsZeroAtTwoPoints = true;}
				else RadiusIsZeroAtOnePoint = true;

				arIndBase2[i] = BasePointCount1;
			}
			else
			{
				//*tVertexPoints2 = r*RadVect2 + z*AxisVect;
				*tVertexPoints2 = CPoiVect + (r*RadVect2 + z*AxisVect); //OC081007_BNL
                tVertexPoints2++; BasePointCount2++;

				arIndBase2[i] = lenArrayOfPoints2d + BasePointCount2;
			}
			tPoints2d++;
		}
		if(RadiusOverlaps) { Send.ErrorMessage("Radia::Error095"); throw(0);}

		AmOfVertexPoints = BasePointCount1 + BasePointCount2;

		bool InvertingBasePointsOrderIsNecessary = false;
		TVector3d vBase1 = VertexPointsArr[1] - VertexPointsArr[0];
		TVector3d vBase2 = VertexPointsArr[2] - VertexPointsArr[1];
		TVector3d ExtNormBase1 = vBase1^vBase2;
		ExtNormBase1 = (1./ExtNormBase1.Abs())*ExtNormBase1;
		if((RadVect1^AxisVect)*ExtNormBase1 < 0.) InvertingBasePointsOrderIsNecessary = true;

		AmOfFaces = lenArrayOfPoints2d + 2;
		AmOfFacesOrig = AmOfFaces;
		if(RadiusIsZeroAtTwoPoints) AmOfFaces--;

		//ArrayOfFaces = new int*[AmOfFaces];
		ArrayOfFaces = new int*[AmOfFacesOrig];
		ArrayOfFaces[0] = new int[lenArrayOfPoints2d];
		ArrayOfFaces[1] = new int[lenArrayOfPoints2d];

		//ArrayOfNumOfPoInFaces = new int[AmOfFaces];
		ArrayOfNumOfPoInFaces = new int[AmOfFacesOrig];
        ArrayOfNumOfPoInFaces[0] = lenArrayOfPoints2d;
        ArrayOfNumOfPoInFaces[1] = lenArrayOfPoints2d;

		int *tBaseInd1 = ArrayOfFaces[0];
		int *tBaseInd2 = ArrayOfFaces[1];
		int **tFaces = ArrayOfFaces + 2;
		int *tMantlePoints = 0;
		int *tNumOfPoInFaces = ArrayOfNumOfPoInFaces + 2;

		int lenArrayOfPoints2d_mi_1 = lenArrayOfPoints2d - 1;
		for(int j=0; j<lenArrayOfPoints2d; j++)
		{
			*tFaces = new int[4];
			tMantlePoints = *tFaces;
            *tNumOfPoInFaces = 4;

			int indInP2d_1=-1, indInP2d_2=-1;

			if(InvertingBasePointsOrderIsNecessary) 
			{
				*tBaseInd1 = lenArrayOfPoints2d - j;
				//*tBaseInd2 = lenArrayOfPoints2d + 1 + j;
				*tBaseInd2 = arIndBase2[j];

				if(j != lenArrayOfPoints2d_mi_1)
				{
					tMantlePoints[0] = lenArrayOfPoints2d - j - 1;
					tMantlePoints[1] = lenArrayOfPoints2d - j;
					//tMantlePoints[2] = 2*lenArrayOfPoints2d - j;
					//tMantlePoints[3] = 2*lenArrayOfPoints2d - j - 1;
					tMantlePoints[2] = arIndBase2[lenArrayOfPoints2d - j - 1];
					tMantlePoints[3] = arIndBase2[lenArrayOfPoints2d - j - 2];

					indInP2d_1 = lenArrayOfPoints2d - j - 2;
					indInP2d_2 = lenArrayOfPoints2d - j - 1;
				}
				else
				{
                    tMantlePoints[0] = lenArrayOfPoints2d;
                    tMantlePoints[1] = 1;
                    //tMantlePoints[2] = lenArrayOfPoints2d + 1;
					//tMantlePoints[3] = 2*lenArrayOfPoints2d;
                    tMantlePoints[2] = arIndBase2[0];
					tMantlePoints[3] = arIndBase2[lenArrayOfPoints2d - 1];

					indInP2d_1 = lenArrayOfPoints2d - 1;
					indInP2d_2 = 0;
                }
			}
			else 
			{
				*tBaseInd1 = j + 1;
                //*tBaseInd2 = 2*lenArrayOfPoints2d - j;
                *tBaseInd2 = arIndBase2[lenArrayOfPoints2d - j - 1];

				if(j != lenArrayOfPoints2d_mi_1)
				{
					//if(ArrayOfPoints2d[j + 1] > 0.)
					//{
					//	tMantlePoints[0] = j + 2;
					//	tMantlePoints[3] = lenArrayOfPoints2d + j + 2;
					//}
					//else
					//{
					//	if(ArrayOfPoints2d[j] > 0.)
					//	{
					//		tMantlePoints[0] = j + 2;
					//		tMantlePoints[1] = j + 1;
					//		tMantlePoints[2] = lenArrayOfPoints2d + j + 1;
					//		*tNumOfPoInFaces = 3;
					//	}
					//}

					tMantlePoints[0] = j + 2;
					tMantlePoints[1] = j + 1;
					//tMantlePoints[2] = lenArrayOfPoints2d + j + 1;
					//tMantlePoints[3] = lenArrayOfPoints2d + j + 2;
					tMantlePoints[2] = arIndBase2[j];
					tMantlePoints[3] = arIndBase2[j + 1];

					indInP2d_1 = j + 1;
					indInP2d_2 = j;
				}
				else
				{
					tMantlePoints[0] = 1;
                    tMantlePoints[1] = lenArrayOfPoints2d;
					//tMantlePoints[2] = 2*lenArrayOfPoints2d;
					//tMantlePoints[3] = lenArrayOfPoints2d + 1;
					tMantlePoints[2] = arIndBase2[lenArrayOfPoints2d - 1];
                    tMantlePoints[3] = arIndBase2[0];

					indInP2d_1 = 0;
					indInP2d_2 = lenArrayOfPoints2d - 1;
				}
			}

			//TVector2d &InP2d_1 = ArrayOfPoints2d[tMantlePoints[0] - 1];
			//TVector2d &InP2d_2 = ArrayOfPoints2d[tMantlePoints[1] - 1];
			TVector2d &InP2d_1 = ArrayOfPoints2d[indInP2d_1];
			TVector2d &InP2d_2 = ArrayOfPoints2d[indInP2d_2];

			if((InP2d_1.x <= 0.) && (InP2d_2.x <= 0.))
			{
				delete[] *tFaces; *tFaces = 0;
			}
			else
			{
				if(InP2d_1.x <= 0.)
				{
					*tNumOfPoInFaces = 3;
				}
				else if(InP2d_2.x <= 0.)
				{
					*tNumOfPoInFaces = 3;
					tMantlePoints[2] = tMantlePoints[3];
				}
				tFaces++; tNumOfPoInFaces++;
			}

			tBaseInd1++; tBaseInd2++; 
		}

		short PrevSendingIsRequired = SendingIsRequired; 
		SendingIsRequired = 0;
        //int IndPolyhedr = SetPolyhedron1(VertexPointsArr, AmOfVertexPoints, ArrayOfFaces, ArrayOfNumOfPoInFaces, AmOfFaces, Magn, 3);
        int IndPolyhedr = SetPolyhedron1(VertexPointsArr, AmOfVertexPoints, ArrayOfFaces, ArrayOfNumOfPoInFaces, AmOfFaces, Magn);

		if(NumberOfSegm > 1)
		{
			radTrans *pRotMlt = new radTrans();
			pRotMlt->SetupRotation(CPoiVect, AxisVect, SegmAngle);
            radThg hg(pRotMlt);
            int IndRotMlt = AddElementToContainer(hg);
			if(ApplySymmetry(IndPolyhedr, IndRotMlt, NumberOfSegm) == 0) throw(0);

			if(FreeSym)
			{
				const char* OptionNames[] = {"FreeSym"};
				const char* OptionValues[] = {"True"};
                IndPolyhedr = DuplicateElement_g3d(IndPolyhedr, OptionNames, OptionValues, 1);

				if((Magn[0] != 0.) || (Magn[1] != 0.) || (Magn[2] != 0.))
				{//OC04082010: we need to "set back" magnetizations in former "symmetry childs"
					TVector3d vM;
					if(!ValidateVector3d(Magn, 3, &vM)) { Send.ErrorMessage("Radia::Error125"); return 0;}

					radThg hgRes;
					if(!ValidateElemKey(IndPolyhedr, hgRes)) { Send.ErrorMessage("Radia::Error125"); return 0;}
					radTg3d* g3dResP = Cast.g3dCast(hgRes.rep); if(g3dResP==0) { Send.ErrorMessage("Radia::Error125"); return 0;}
					radTGroup* GroupResP = Cast.GroupCast(g3dResP); if(GroupResP==0) { Send.ErrorMessage("Radia::Error125"); return 0;}
					for(radTmhg::const_iterator iter = GroupResP->GroupMapOfHandlers.begin(); iter != GroupResP->GroupMapOfHandlers.end(); ++iter)
					{
						radTg3d* pG3D = (radTg3d*)(((*iter).second).rep);
						radTlphg &curListTrf = pG3D->g3dListOfTransform;
						TVector3d vMb(vM);

						if(!curListTrf.empty())
						{
							radTvhg vhFlatTrfs;
							pG3D->FlattenSpaceTransforms(vhFlatTrfs);
							if(vhFlatTrfs.size() >= 1)
							{
								radTrans* pTrf = (radTrans*)((vhFlatTrfs[0]).rep);
								if(pTrf != 0)
								{
									vMb = pTrf->TrBiPoint_inv(vM);
								}
							}
						}
						pG3D->SetM(vMb);
					}
				}
			}
		}
        SendingIsRequired = PrevSendingIsRequired;

		if(VertexPointsArr != 0) delete[] VertexPointsArr;
		if((ArrayOfFaces != 0) && (AmOfFaces > 0)) 
		{
			int** tFaces = ArrayOfFaces;
			//for(int i=0; i<AmOfFaces; i++)
			for(int i=0; i<AmOfFacesOrig; i++)
			{
				if(*tFaces != 0) 
				{ 
					delete[] (*tFaces); *tFaces = 0;
				}
				tFaces++;
			}
			delete[] ArrayOfFaces;
		}
		if(ArrayOfNumOfPoInFaces != 0) delete[] ArrayOfNumOfPoInFaces;
		if(arIndBase2 != 0) delete[] arIndBase2;

		if(SendingIsRequired) Send.Int(IndPolyhedr); 
		return IndPolyhedr;
	}
	catch(...)
	{
		if(VertexPointsArr != 0) delete[] VertexPointsArr;
		if((ArrayOfFaces != 0) && (AmOfFaces > 0)) 
		{
			int** tFaces = ArrayOfFaces;
			for(int i=0; i<AmOfFaces; i++)
			{
				if(*tFaces != 0) { delete[] (*tFaces); *tFaces = 0;}
				tFaces++;
			}
			delete[] ArrayOfFaces;
		}
		if(ArrayOfNumOfPoInFaces != 0) delete[] ArrayOfNumOfPoInFaces;
		if(arIndBase2 != 0) delete[] arIndBase2;

		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::TransformBackMagnOrCurDensArr(int IndTr, double* Magn, long lenMagn)
{//Assumes Magn[3] 
	if((IndTr <= 0) || (Magn == 0)) return;

	radThg hg;
	if(!ValidateElemKey(IndTr, hg)) return;
	radTrans* TransPtr = Cast.TransCast(hg.rep);
	if(TransPtr==0) { Send.ErrorMessage("Radia::Error006"); throw(0);}

	TVector3d MagnVect, NewMagnVect;
	if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return;
	NewMagnVect = TransPtr->TrBiPoint_inv(MagnVect);
	Magn[0] = NewMagnVect.x;
	Magn[1] = NewMagnVect.y;
	Magn[2] = NewMagnVect.z;
}

//-------------------------------------------------------------------------

void radTApplication::TransformBackPointArr(int IndTr, double* arP, long lenP)
{
	if((IndTr <= 0) || (arP == 0)) return;

	radThg hg;
	if(!ValidateElemKey(IndTr, hg)) return;
	radTrans* TransPtr = Cast.TransCast(hg.rep);
	if(TransPtr==0) { Send.ErrorMessage("Radia::Error006"); throw(0);}

	TVector3d vP, new_vP;
	if(!ValidateVector3d(arP, lenP, &vP)) return;
	new_vP = TransPtr->TrPoint_inv(vP);
	arP[0] = new_vP.x;
	arP[1] = new_vP.y;
	arP[2] = new_vP.z;
}

//-------------------------------------------------------------------------

int radTApplication::SetCylMag(double* CPoi, long lenCPoi, double r, double h, int NumberOfSegm, double* Magn, long lenMagn, char* Orient)
{
	try
	{
		TVector3d CPoiVect; //, MagnVect;
		if(!ValidateVector3d(CPoi, lenCPoi, &CPoiVect)) return 0;
		if(r <= 0) { Send.ErrorMessage("Radia::Error089"); return 0;}
		if(h <= 0) { Send.ErrorMessage("Radia::Error012"); return 0;}
		if(NumberOfSegm <= 0) { Send.ErrorMessage("Radia::Error090"); return 0;}
		//if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) return 0;

		double FirstPoi[] = {CPoiVect.x - 0.5*h, CPoiVect.y + r, CPoiVect.z}; //assuming *Orient == 'x'
		if((*Orient == 'y') || (*Orient == 'Y')) //OC090106
		{
			FirstPoi[0] = CPoiVect.x;
			FirstPoi[1] = CPoiVect.y - 0.5*h;
			FirstPoi[2] = CPoiVect.z + r;
		}
		else if((*Orient == 'z') || (*Orient == 'Z'))
		{
			FirstPoi[0] = CPoiVect.x + r;
			FirstPoi[1] = CPoiVect.y;
			FirstPoi[2] = CPoiVect.z - 0.5*h;
		}
		
		for(int i=0; i<3; i++) FirstPoi[i] = radCR.Double(FirstPoi[i]);

		const double Pi = 3.141592653589793238;
		const double TwoPi = 2.*Pi;
		double AngSegm = TwoPi/NumberOfSegm;

		TVector2d* ArrayOfPoints2d = new TVector2d[NumberOfSegm];
		TVector2d *tPoint2D = ArrayOfPoints2d;
		double CurAng = 0.;

		for(int k=0; k<NumberOfSegm; k++)
		{
			tPoint2D->x = r*cos(CurAng);
			tPoint2D->y = r*sin(CurAng);
			CurAng += AngSegm;
			tPoint2D++;
		}

		short PrevSendingIsRequired = SendingIsRequired; 
		SendingIsRequired = 0;
		int CylInd = SetExtrudedPolygon(FirstPoi, 3, h, ArrayOfPoints2d, NumberOfSegm, Magn, lenMagn, Orient);
		SendingIsRequired = PrevSendingIsRequired;

		if(ArrayOfPoints2d != 0) { delete[] ArrayOfPoints2d; ArrayOfPoints2d = 0;}

		if(SendingIsRequired) Send.Int(CylInd);
		return CylInd;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::FindSpaceTransToOrientObjAlongMainAxis(double* CPoi, char DefOrient, char Orient)
{//this does not send anything
	try
	{
		DefOrient = (char)toupper(DefOrient);
		Orient = (char)toupper(Orient);

		//short PrevSendingIsRequired = SendingIsRequired;
		if((CPoi == 0) || (Orient == DefOrient)) return 0;

		const double Pi = 3.141592653589793238;
		const double HalfPi = 0.5*Pi;
		const double QuartPi = 0.25*Pi;

		double AxisVect[] = {1, 1, 1};
		double RotAng = HalfPi;
		double NormVect[] = {0, 0, 0};
		int TrfInd = 0;

		if(DefOrient == 'X')
		{
			if(Orient == 'Y') 
			{
				//NormVect[0] = -1;
				//NormVect[1] = 1;
				//TrfInd = SetPlaneSym(CPoi, 3, NormVect, 3, 1);
				RotAng = 2.*Pi/3.;
				TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
			}
			else if(Orient == 'Z') 
			{
				//AxisVect[0] = AxisVect[1] = AxisVect[2] = -1;
				//RotAng = 2.*Pi/3.;
				RotAng = 4.*Pi/3.;
				TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
			}
		}
		else if(DefOrient == 'Y')
		{
			if(Orient == 'X') 
			{
				//NormVect[0] = -1;
				//NormVect[1] = 1;
				//TrfInd = SetPlaneSym(CPoi, 3, NormVect, 3, 1);
				RotAng = 4.*Pi/3.;
				TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
			}
			else if(Orient == 'Z') 
			{
				//NormVect[1] = -1;
				//NormVect[2] = 1;
				//TrfInd = SetPlaneSym(CPoi, 3, NormVect, 3, 1);
				RotAng = 2.*Pi/3.;
				TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
			}
		}
		else if(DefOrient == 'Z')
		{
			if(Orient == 'X') 
			{
				RotAng = 2.*Pi/3.;
                TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
			}
			else if(Orient == 'Y') 
			{
				//NormVect[1] = -1;
				//NormVect[2] = 1;
				//TrfInd = SetPlaneSym(CPoi, 3, NormVect, 3, 1);
				RotAng = 4.*Pi/3.;
                TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
			}
		}
		else 
		{
			return 0;
		}

		return TrfInd;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

//int radTApplication::OrientObjAlongMainAxis(int ObjInd, double* CPoi, char DefOrient, char Orient)
//{
//	try
//	{
//		DefOrient = (char)toupper(DefOrient);
//		Orient = (char)toupper(Orient);
//
//		//short PrevSendingIsRequired = SendingIsRequired;
//		if((ObjInd <= 0) || (CPoi == 0) || (Orient == DefOrient)) 
//		{
//			//if(PrevSendingIsRequired) Send.Int(ObjInd); 
//			return ObjInd;
//		}
//
//		const double Pi = 3.141592653589793238;
//		const double HalfPi = 0.5*Pi;
//
//		double AxisVect[] = {0, 0, 0};
//		double RotAng = HalfPi;
//		//double TranslVect[] = {0, 0, 0};
//
//		if(DefOrient == 'X')
//		{
//			if(Orient == 'Y') 
//			{
//				AxisVect[2] = -1;
//				//TranslVect[1] = 2*CPoi[0];
//			}
//			else if(Orient == 'Z') 
//			{
//				AxisVect[0] = AxisVect[1] = AxisVect[2] = -1;
//				RotAng = 2.*Pi/3.;
//			}
//		}
//		else if(DefOrient == 'Y')
//		{
//			if(Orient == 'X') 
//			{
//				AxisVect[2] = 1;
//				//TranslVect[0] = 2*CPoi[1];
//			}
//			else if(Orient == 'Z') 
//			{
//				AxisVect[0] = AxisVect[1] = AxisVect[2] = 1;
//				RotAng = 2.*Pi/3.;
//			}
//		}
//		else if(DefOrient == 'Z')
//		{
//			if(Orient == 'X') 
//			{
//				AxisVect[0] = AxisVect[1] = AxisVect[2] = 1;
//				RotAng = 2.*Pi/3.;
//			}
//			else if(Orient == 'Y') 
//			{
//				AxisVect[0] = AxisVect[1] = AxisVect[2] = -1;
//				RotAng = 2.*Pi/3.;
//			}
//		}
//		else 
//		{
//			//if(PrevSendingIsRequired) Send.Int(ObjInd); 
//			return ObjInd;
//		}
//
//		short PrevSendingIsRequired = SendingIsRequired;
//		SendingIsRequired = 0;
//
//		//double CenRot[] = {0, 0, 0};
//		int RotInd = 0;
//		if((AxisVect[0] != 0) || (AxisVect[1] != 0) || (AxisVect[2] != 0))
//		{
//			//RotInd = SetRotation(CenRot, 3, AxisVect, 3, RotAng);
//			RotInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
//		}
//
//		int res = ObjInd;
//		if(RotInd != 0) res = ApplySymmetry(ObjInd, RotInd, 1);	
//
//		SendingIsRequired = PrevSendingIsRequired;
//		return res;
//	}
//	catch(...)
//	{
//		Initialize(); return 0;
//	}
//}

//-------------------------------------------------------------------------

int radTApplication::SetArcCur(double* CPoi, long lenCPoi, double* Radii, long lenRadii, double* Angles, long lenAngles, double InHeight, double InJ_azim, int InNumberOfSegm, char* ManOrAuto, char* Orient)
{
	try
	{
		TVector3d CPoiVect;
		if(!ValidateVector3d(CPoi, lenCPoi, &CPoiVect)) return 0;
		// May be not necessary?
		if((lenRadii!=2) || (Radii[0]<0.) || (Radii[1]<0.) || (Radii[1]<Radii[0]))
		{
			Send.ErrorMessage("Radia::Error010"); return 0;
		}
		const double TwoPi = 2.*3.141592653589793238;
		//if((lenAngles!=2) || (Angles[0]<0. || Angles[0]>TwoPi) || (Angles[0]<0. || Angles[0]>TwoPi) || (Angles[1]<Angles[0]))
		if((lenAngles != 2) || (Angles[1] < Angles[0]) || ((Angles[1] - Angles[0]) > TwoPi)) //OC270308
		{
			Send.ErrorMessage("Radia::Error011"); return 0;
		}
		if(InHeight <= 0)
		{
			Send.ErrorMessage("Radia::Error012");	return 0;
		}

		short InBasedOnPrecFlag;
		if((!strcmp(ManOrAuto, "auto")) || (!strcmp(ManOrAuto, "Auto")) || (!strcmp(ManOrAuto, "AUTO"))) InBasedOnPrecFlag = 1;
		else if((!strcmp(ManOrAuto, "man")) || (!strcmp(ManOrAuto, "Man")) || (!strcmp(ManOrAuto, "MAN"))) InBasedOnPrecFlag = 0;
		else { Send.ErrorMessage("Radia::Error045"); return 0;}

		short PrevSendingIsRequired = SendingIsRequired; SendingIsRequired = 0; //OC301207
		int IndTr = FindSpaceTransToOrientObjAlongMainAxis(CPoi, 'Z', *Orient);
		//if((Magn != 0) && (lenMagn >= 3) && (IndTr != 0)) TransformBackMagnOrCurDensArr(IndTr, Magn, lenMagn);

		radThg hg(new radTArcCur(CPoiVect, Radii, Angles, InHeight, InJ_azim, InNumberOfSegm, InBasedOnPrecFlag));
		int ElemKey = AddElementToContainer(hg);

		if(ElemKey <= 0) { Send.ErrorMessage("Radia::Error000"); return 0;}

		//short PrevSendingIsRequired = SendingIsRequired; SendingIsRequired = 0; //OC301207
        if(IndTr != 0) ElemKey = ApplySymmetry(ElemKey, IndTr, 1);
        SendingIsRequired = PrevSendingIsRequired;

		//ElemKey = OrientObjAlongMainAxis(ElemKey, CPoi, 'Z', *Orient);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetRectangle(double* CPoi, long lenCPoi, double* Dims, long lenDims)
{
	try
	{
		TVector3d CPoiVect;
		TVector2d DimsVect;
		if(!ValidateVector3d(CPoi, lenCPoi, &CPoiVect)) return 0;
		if(!ValidateVector2d(Dims, lenDims, &DimsVect)) return 0;
		// May be not necessary?
		if((Dims[0]<0.) || (Dims[1]<0.))
		{
			Send.ErrorMessage("Radia::Error001"); return 0;
		}
		radThg hg(new radTRectangle(CPoiVect, DimsVect));
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetBackgroundFieldSource(double* B, long lenB)
{
	TVector3d B_Vect;
	if(!ValidateVector3d(B, lenB, &B_Vect)) return 0;

	try
	{
		radThg hg(new radTBackgroundFieldSource(B_Vect));
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetGroup(int* ArrayOfKeys, long lenArrayOfKeys)
{
	try
	{
		radTGroup* GroupPtr = new radTGroup();
		radThg hg;
		for(long i = 0; i < lenArrayOfKeys; i++)
		{
			if(!ValidateElemKey(ArrayOfKeys[i], hg)) { delete GroupPtr; return 0; }
			if(Cast.g3dCast(hg.rep)==0) { Send.ErrorMessage("Radia::Error003"); delete GroupPtr; return 0; }
			GroupPtr->AddElement(ArrayOfKeys[i], hg); 
		}
		radThg hgGroup(GroupPtr);
		int ElemKey = AddElementToContainer(hgGroup);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::AddToGroup(int GroupKey, int* ArrayOfKeys, long lenArrayOfKeys)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(GroupKey, hg)) return 0;
		radTg3d* g3dP = Cast.g3dCast(hg.rep); if(g3dP==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
		radTGroup* GroupP = Cast.GroupCast(g3dP); if(GroupP==0) { Send.ErrorMessage("Radia::Error004"); return 0;}
		radThg hgInner;
		for(long i = 0; i < lenArrayOfKeys; i++)
		{
			if(!ValidateElemKey(ArrayOfKeys[i], hgInner)) return 0;
			if(Cast.g3dCast(hgInner.rep)==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
			if(ArrayOfKeys[i]==GroupKey) { Send.ErrorMessage("Radia::Error005"); return 0;}
			//Maybe more sophisticated testing for nested groups is needed
			GroupP->AddElement(ArrayOfKeys[i], hgInner);
		}

		if(SendingIsRequired) Send.Int(GroupKey);
		return GroupKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::OutGroupSize(int ElemKey)
//int radTApplication::OutGroupSize(int ElemKey, char deep)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		radTGroup* GroupPtr = Cast.GroupCast(g3dPtr); 
		if(GroupPtr==0) { if(SendingIsRequired) Send.Int(0); return 0;}
		else
		{
			int lenGroupSubObjArray = (int)(GroupPtr->GroupMapOfHandlers.size());

			if(SendingIsRequired) Send.Int(lenGroupSubObjArray);
			return lenGroupSubObjArray;
		}
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::OutGroupSubObjectKeys(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		radTGroup* GroupPtr = Cast.GroupCast(g3dPtr); 
		if(GroupPtr==0) { if(SendingIsRequired) Send.IntList(&ElemKey, 1); return ElemKey;}
		else
		{
			int* GroupSubObjArray = NULL;
			int lenGroupSubObjArray = (int)(GroupPtr->GroupMapOfHandlers.size());
			GroupSubObjArray = new int[lenGroupSubObjArray];

			int* ArrayTravers = GroupSubObjArray;
			for(radTmhg::const_iterator iter = GroupPtr->GroupMapOfHandlers.begin();
				iter != GroupPtr->GroupMapOfHandlers.end(); ++iter)
				*(ArrayTravers++) = (*iter).first;

			if(SendingIsRequired) Send.IntList(GroupSubObjArray, lenGroupSubObjArray);
			delete[] GroupSubObjArray;
			return ElemKey;
		}
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetRaceTrack(double* CPoi, long lenCPoi,double* Radii, long lenRadii, double* StrPartDims, long lenStrPartDims, double InHeight, double InJ_azim, int InNumberOfSegm, char* ManOrAuto, char* Orient)
{
	try
	{
		TVector3d CPoiVect;
		if(!ValidateVector3d(CPoi, lenCPoi, &CPoiVect)) return 0;
		// May be not necessary?
		if((lenRadii!=2) || (Radii[0]<0.) || (Radii[1]<0.) || (Radii[1]<Radii[0]))
		{
			Send.ErrorMessage("Radia::Error010"); return 0;
		}
		if((lenStrPartDims!=2) || (StrPartDims[0]<0.) || (StrPartDims[1]<0.))
		{
			Send.ErrorMessage("Radia::Error029"); return 0;
		}
		if(InHeight <= 0)
		{
			Send.ErrorMessage("Radia::Error012"); return 0;
		}

		short PrevSendingFlag = SendingIsRequired;
		SendingIsRequired = 0;

		const double Pi = 3.14159265358979324;

        int IndTr = FindSpaceTransToOrientObjAlongMainAxis(CPoi, 'Z', *Orient);

		double ArcCurCPoi[] = {CPoi[0]+0.5*StrPartDims[0], CPoi[1]+0.5*StrPartDims[1], CPoi[2]};
		double ArcCurAngles[] = {0., 0.5*Pi};
		//int FirstArcCur = SetArcCur(ArcCurCPoi, 3, Radii, 2, ArcCurAngles, 2, InHeight, InJ_azim, InNumberOfSegm, ManOrAuto, "Z");
		int FirstArcCur = SetArcCur(ArcCurCPoi, 3, Radii, 2, ArcCurAngles, 2, InHeight, InJ_azim, InNumberOfSegm, ManOrAuto, (char*)"Z"); //OC01052013 to walk-around of a warning in GCC
		if(FirstArcCur==0) return 0;
		ArcCurCPoi[0] = CPoi[0]-0.5*StrPartDims[0]; ArcCurCPoi[1] = CPoi[1]+0.5*StrPartDims[1];
		ArcCurAngles[0] = 0.5*Pi; ArcCurAngles[1] = Pi;
		int SecondArcCur = SetArcCur(ArcCurCPoi, 3, Radii, 2, ArcCurAngles, 2, InHeight, InJ_azim, InNumberOfSegm, ManOrAuto, (char*)"Z");
		if(SecondArcCur==0) return 0;
		ArcCurCPoi[0] = CPoi[0]-0.5*StrPartDims[0]; ArcCurCPoi[1] = CPoi[1]-0.5*StrPartDims[1];
		ArcCurAngles[0] = Pi; ArcCurAngles[1] = 1.5*Pi;
		int ThirdArcCur = SetArcCur(ArcCurCPoi, 3, Radii, 2, ArcCurAngles, 2, InHeight, InJ_azim, InNumberOfSegm, ManOrAuto, (char*)"Z");
		if(ThirdArcCur==0) return 0;
		ArcCurCPoi[0] = CPoi[0]+0.5*StrPartDims[0]; ArcCurCPoi[1] = CPoi[1]-0.5*StrPartDims[1];
		ArcCurAngles[0] = 1.5*Pi; ArcCurAngles[1] = 1.999999999999*Pi;
		int FourthArcCur = SetArcCur(ArcCurCPoi, 3, Radii, 2, ArcCurAngles, 2, InHeight, InJ_azim, InNumberOfSegm, ManOrAuto, (char*)"Z");
		if(FourthArcCur==0) return 0;

		int ArrayOfArcs[] = {FirstArcCur, SecondArcCur, ThirdArcCur, FourthArcCur};
		int RaceTrack = SetGroup(ArrayOfArcs, 4);

		double ZeroArray[] = {0.,0.,0.};
		if(StrPartDims[1]!=0.)
		{
			double EffRad = 0.5*(StrPartDims[0]+Radii[0]+Radii[1]);
			double RecCurCPoi[] = {CPoi[0]+EffRad, CPoi[1], CPoi[2]};
			double RecCurDims[] = {Radii[1]-Radii[0], StrPartDims[1], InHeight};
			double RecCurCur[] = {0., InJ_azim, 0.};
			int FirstRecCur = SetRecMag(RecCurCPoi, 3, RecCurDims, 3, ZeroArray, 3, RecCurCur, 3, 1);
			if(FirstRecCur==0) return 0;
			RecCurCPoi[0] = CPoi[0]-EffRad; RecCurCur[1] = -InJ_azim;
			int SecondRecCur = SetRecMag(RecCurCPoi, 3, RecCurDims, 3, ZeroArray, 3, RecCurCur, 3, 1);
			if(SecondRecCur==0) return 0;
			int ArrayOfRecs[] = {FirstRecCur, SecondRecCur};
			AddToGroup(RaceTrack, ArrayOfRecs, 2);
		}
		if(StrPartDims[0]!=0.)
		{
			double EffRad = 0.5*(StrPartDims[1]+Radii[0]+Radii[1]);
			double RecCurCPoi[] = {CPoi[0], CPoi[1]+EffRad, CPoi[2]};
			double RecCurDims[] = {StrPartDims[0], Radii[1]-Radii[0], InHeight};
			double RecCurCur[] = {-InJ_azim, 0., 0.};
			int ThirdRecCur = SetRecMag(RecCurCPoi, 3, RecCurDims, 3, ZeroArray, 3, RecCurCur, 3, 1);
			if(ThirdRecCur==0) return 0;
			RecCurCPoi[1] = CPoi[1]-EffRad; RecCurCur[0] = InJ_azim;
			int FourthRecCur = SetRecMag(RecCurCPoi, 3, RecCurDims, 3, ZeroArray, 3, RecCurCur, 3, 1);
			if(FourthRecCur==0) return 0;
			int ArrayOfRecs[] = {ThirdRecCur, FourthRecCur};
			AddToGroup(RaceTrack, ArrayOfRecs, 2);
		}

        //short PrevSendingIsRequired = SendingIsRequired; SendingIsRequired = 0;
        if(IndTr != 0) RaceTrack = ApplySymmetry(RaceTrack, IndTr, 1);

		SendingIsRequired = PrevSendingFlag;

		if(RaceTrack <= 0)
		{ 
			Send.ErrorMessage("Radia::Error000"); return 0;
		}
		//return OrientObjAlongMainAxis(RaceTrack, CPoi, 'Z', *Orient);
		//RaceTrack = OrientObjAlongMainAxis(RaceTrack, CPoi, 'Z', *Orient);

		if(SendingIsRequired) Send.Int(RaceTrack);
		return RaceTrack;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetFlmCur(double I, TVector3d* ArrayOfPoints, int lenArrayOfPoints)
{
	int FlmCurKey;
	try
	{
		if(lenArrayOfPoints<2) { Send.ErrorMessage("Radia::Error035"); return 0;}
		else if(lenArrayOfPoints==2)
		{
			radThg hg(new radTFlmLinCur(ArrayOfPoints[0], ArrayOfPoints[1], I));
			FlmCurKey = AddElementToContainer(hg);
		}
		else
		{
			int AmOfFlmLinCur = lenArrayOfPoints-1;
			int* ArrayOfFlmLinCur = NULL;
			ArrayOfFlmLinCur = new int[AmOfFlmLinCur];
			for(int i=0; i<AmOfFlmLinCur; i++)
			{
				radThg hg(new radTFlmLinCur(ArrayOfPoints[i], ArrayOfPoints[i+1], I));
				ArrayOfFlmLinCur[i] = AddElementToContainer(hg);
			}
			short PrevSendingFlag = SendingIsRequired;
			SendingIsRequired = 0;
			FlmCurKey = SetGroup(ArrayOfFlmLinCur, AmOfFlmLinCur);
			SendingIsRequired = PrevSendingFlag;
		}
		if(SendingIsRequired) Send.Int(FlmCurKey);
		return FlmCurKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ComputeNumberOfDegOfFreedom(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		int NumDegFr = g3dPtr->NumberOfDegOfFreedom();

		if(SendingIsRequired) Send.Int(NumDegFr);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeMagnOrJ_InCenter(int ElemKey, char MorJ)
//void radTApplication::ComputeMagnInCenter(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radTFieldKey FieldKey;
		if((MorJ == 'M') || (MorJ == 'm')) FieldKey.M_ = 1;
		else if((MorJ == 'J') || (MorJ == 'j')) FieldKey.J_ = 1;
		else if((MorJ == 'A') || (MorJ == 'a')) FieldKey.A_ = 1;
		else if((MorJ == 'B') || (MorJ == 'b')) FieldKey.B_ = 1;
		else if((MorJ == 'H') || (MorJ == 'h')) FieldKey.H_ = 1;
		else { Send.ErrorMessage("Radia::Error070"); return;}

		radThg hgDplWithoutSym;
		char PutNewStuffIntoGenCont = 0;
		if(!g3dPtr->CreateFromSym(hgDplWithoutSym, this, PutNewStuffIntoGenCont)) return;
		radTg3d* g3dDplWithoutSymPtr = (radTg3d*)(hgDplWithoutSym.rep);

		radTvhg vhFlatTransforms; //OC061007_BNL
		g3dDplWithoutSymPtr->FlattenSpaceTransforms(vhFlatTransforms);
		if(vhFlatTransforms.size() > 0)
		{
			g3dDplWithoutSymPtr->EraseAllTransformations();
			g3dDplWithoutSymPtr->AddTransform(1, vhFlatTransforms[0]);
		}

		radTVectPairOfVect3d VectPairOfVect3d;
		radTrans* pBaseTrans = 0;
		g3dDplWithoutSymPtr->Push_backCenterPointAndField(&FieldKey, &VectPairOfVect3d, pBaseTrans, g3dPtr, this);
		//g3dPtr->Push_backCenterPointAndField(&FieldKey, &VectPairOfVect3d, pBaseTrans, g3dPtr); //use this to calc. M without symmetry childs

		//if(SendingIsRequired) OutCenFieldCompRes(&VectPairOfVect3d);
		if(SendingIsRequired) Send.ArrayOfPairOfVect3d(&VectPairOfVect3d);

		VectPairOfVect3d.erase(VectPairOfVect3d.begin(), VectPairOfVect3d.end());
	}
	catch(...)
	{
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ScaleCurrent(int ElemKey, double scaleCoef)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		int curWasScaled = g3dPtr->ScaleCurrent(scaleCoef); 
		if(!curWasScaled) Send.WarningMessage("Radia::Warning018");

		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetObjMagn(int ElemKey, double mx, double my, double mz)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d *g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
		
		radTg3dRelax *g3dRelaxPtr = Cast.g3dRelaxCast(g3dPtr);
		if(g3dRelaxPtr==0) 
		{
			radTGroup* GroupPtr = Cast.GroupCast(g3dPtr);
			if(GroupPtr==0) { Send.ErrorMessage("Radia::Error015"); return 0;}
		}
		
		TVector3d M(mx, my, mz);
		g3dPtr->SetM(M);

		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTApplication::SetTranslation(double* TranslArray, long lenTranslArray)
{
	try
	{
		TVector3d TranslVect;
		if(!ValidateVector3d(TranslArray, lenTranslArray, &TranslVect)) return 0;

		TVector3d St0(1.,0.,0.), St1(0.,1.,0.), St2(0.,0.,1.);
		TMatrix3d E(St0, St1, St2);

		radTrans* TransPtr = new radTrans(E, E, TranslVect, 1., 1., 1); // ID_No = 1
		if(TransPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(TransPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetRotation(double* PoiOnAxArray, long lenPoiOnAxArray, double* AxVectArray, long lenAxVectArray, double Angle)
{
	try
	{
		TVector3d PoiOnAxVect, AxVect;
		if(!ValidateVector3d(PoiOnAxArray, lenPoiOnAxArray, &PoiOnAxVect)) return 0;
		if(!ValidateVector3d(AxVectArray, lenAxVectArray, &AxVect)) return 0;

		double NormFact = 1./sqrt(AxVect.x*AxVect.x+AxVect.y*AxVect.y+AxVect.z*AxVect.z);
		AxVect = NormFact*AxVect;
		double VxVx, VyVy, VzVz;
		VxVx=AxVect.x*AxVect.x; VyVy=AxVect.y*AxVect.y; VzVz=AxVect.z*AxVect.z;

		double cosAng, sinAng, One_m_cos;
		cosAng = cos(Angle); sinAng = sin(Angle); One_m_cos = 1. - cosAng;
		double One_m_cosVxVy, One_m_cosVxVz, One_m_cosVyVz, sinVx, sinVy, sinVz;
		One_m_cosVxVy = One_m_cos*AxVect.x*AxVect.y;
		One_m_cosVxVz = One_m_cos*AxVect.x*AxVect.z;
		One_m_cosVyVz = One_m_cos*AxVect.y*AxVect.z;
		sinVx = sinAng*AxVect.x; sinVy = sinAng*AxVect.y; sinVz = sinAng*AxVect.z;

		TVector3d St0(VxVx+cosAng*(VyVy+VzVz), One_m_cosVxVy-sinVz, One_m_cosVxVz+sinVy);
		TVector3d St1(One_m_cosVxVy+sinVz, VyVy+cosAng*(VxVx+VzVz), One_m_cosVyVz-sinVx);
		TVector3d St2(One_m_cosVxVz-sinVy, One_m_cosVyVz+sinVx, VzVz+cosAng*(VxVx+VyVy));
		TMatrix3d M(St0, St1, St2);
		TVector3d St00(1.-St0.x, -St0.y, -St0.z);
		TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
		TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
		TMatrix3d M0(St00, St01, St02);

		radTrans* TransPtr = new radTrans(M, M0*PoiOnAxVect, 1., 1., 2); // ID_No = 2
		if(TransPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(TransPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetPlaneSym(double* PoiOnPlane, long lenPoiOnPlane, double* PlaneNormal, long lenPlaneNormal, int s)
{
	try
	{
		TVector3d PoiOnPlaneVect, N;
		if(!ValidateVector3d(PoiOnPlane, lenPoiOnPlane, &PoiOnPlaneVect)) return 0;
		if(!ValidateVector3d(PlaneNormal, lenPlaneNormal, &N)) return 0;

		double NormFact = 1./sqrt(N.x*N.x+N.y*N.y+N.z*N.z);
		N = NormFact*N;

		TVector3d St0(1.-2.*N.x*N.x, -2.*N.x*N.y, -2.*N.x*N.z);
		TVector3d St1(St0.y, 1.-2.*N.y*N.y, -2.*N.y*N.z);
		TVector3d St2(St0.z, St1.z, 1.-2.*N.z*N.z);
		TMatrix3d M(St0, St1, St2);
		TVector3d St00(1.-St0.x, -St0.y, -St0.z);
		TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
		TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
		TMatrix3d M0(St00, St01, St02);

		radTrans* TransPtr = new radTrans(M, M0*PoiOnPlaneVect, -1., double(s), 3); // ID_No = 3
		if(TransPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(TransPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetFieldInversion()
{
	try
	{
		TVector3d St0(1.,0.,0.), St1(0.,1.,0.), St2(0.,0.,1.), ZeroVect(0.,0.,0.);
		TMatrix3d E(St0, St1, St2);

		radTrans* TransPtr = new radTrans(E, E, ZeroVect, 1., -1., 4); // ID_No = 4
		if(TransPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(TransPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::CombineTransformations(int ThisElemKey, int AnotherElemKey, char L_or_R)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ThisElemKey, hg)) return 0;
		radTrans* ThisTransPtr = Cast.TransCast(hg.rep);
		if(ThisTransPtr==0) { Send.ErrorMessage("Radia::Error006"); return 0;}

		if(!ValidateElemKey(AnotherElemKey, hg)) return 0;
		radTrans* AnotherTransPtr = Cast.TransCast(hg.rep);
		if(AnotherTransPtr==0) { Send.ErrorMessage("Radia::Error006"); return 0;}

		if(L_or_R == 'L') *ThisTransPtr = Product(*AnotherTransPtr, *ThisTransPtr);
		else if(L_or_R == 'R') *ThisTransPtr = Product(*ThisTransPtr, *AnotherTransPtr);
		else return 0;
		if(SendingIsRequired) Send.Int(ThisElemKey);
		return ThisElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ApplySymmetry(int g3dElemKey, int TransElemKey, int Multiplicity)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(g3dElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		if(!ValidateElemKey(TransElemKey, hg)) return 0;
		if(Cast.TransCast(hg.rep)==0) { Send.ErrorMessage("Radia::Error006"); return 0;}
		g3dPtr->AddTransform(Multiplicity, hg);
		if(SendingIsRequired) Send.Int(g3dElemKey);
		return g3dElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------
