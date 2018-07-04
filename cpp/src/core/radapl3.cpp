/*-------------------------------------------------------------------------
*
* File name:      radapl3.cpp
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
#include "radptrj.h"
#include "radg3da1.h"
#include "radopnam.h"

#include <math.h>
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTApplication::ComputeField(int ElemKey, char* FieldChar, double* StObsPoi, long lenStObsPoi, 
								   double* FiObsPoi, long lenFiObsPoi, int Np, char* ShowArgFlag, double StrtArg)
{
	radTField* FieldArray = NULL;
	double* ArgArray = NULL;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		radTFieldKey FieldKey;
		if(!ValidateFieldChar(FieldChar, &FieldKey)) return;
		TVector3d StObsPoiVect, FiObsPoiVect;
		if(!ValidateVector3d(StObsPoi, lenStObsPoi, &StObsPoiVect)) return;
		if(!ValidateVector3d(FiObsPoi, lenFiObsPoi, &FiObsPoiVect)) return;

		short ArgumentNeeded = 0;
		if(!strcmp(ShowArgFlag, "arg")) ArgumentNeeded = 1;
		else if(strcmp(ShowArgFlag, "noarg")) { Send.ErrorMessage("Radia::Error034"); return;}

		if(Np==1 && (FiObsPoiVect.x < 1.E+22) && (FiObsPoiVect.y < 1.E+22) && (FiObsPoiVect.z < 1.E+22)) Np = 101; // New Default 

		TVector3d ZeroVect(0.,0.,0.);
		TVector3d ObsPoiVect = StObsPoiVect;

		radTField Field(FieldKey, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		g3dPtr->B_genComp(&Field);

		if(Np>1)
		{
			FieldArray = new radTField[Np];
			if(ArgumentNeeded) ArgArray = new double[Np];

			FieldArray[0] = Field;
			TVector3d TranslVect = (1./double(Np-1))*(FiObsPoiVect-StObsPoiVect);
			double StepArg;
			if(ArgumentNeeded) 
			{
				ArgArray[0] = StrtArg;
				StepArg	= sqrt(TranslVect.x*TranslVect.x + TranslVect.y*TranslVect.y + TranslVect.z*TranslVect.z);
			}

			for(int i=1; i<Np; i++)
			{
				ObsPoiVect += TranslVect;
				FieldArray[i] = radTField(FieldKey, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
				g3dPtr->B_genComp(&(FieldArray[i]));
				if(ArgumentNeeded) ArgArray[i] = ArgArray[i-1] + StepArg;
			}
		}
		else FieldArray = &Field;

		if(SendingIsRequired) Send.OutFieldCompRes(FieldChar, FieldArray, ArgArray, Np);
		if(Np>1) 
		{
			delete[] FieldArray;
			FieldArray = NULL;
			if(ArgumentNeeded) delete[] ArgArray;
			ArgArray = NULL;
		}
	}
	catch(...) 
	{ 
		if(FieldArray != NULL) delete[] FieldArray;
		if(ArgArray != NULL) delete[] ArgArray;
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeField(int ElemKey, char* FieldChar, radTVectorOfVector3d& VectorOfVector3d, radTVectInputCell& VectInputCell)
{
	radTField* FieldArray = 0;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		radTFieldKey FieldKey;
		if(!ValidateFieldChar(FieldChar, &FieldKey)) return;

		long Np = (long)(VectorOfVector3d.size());
		FieldArray = new radTField[Np];
		if(FieldArray == 0) { Send.ErrorMessage("Radia::Error900"); return;}
		radTField* tField = FieldArray;

		TVector3d ZeroVect(0.,0.,0.);
		for(long i=0; i<Np; i++)
		{
			*tField = radTField(FieldKey, CompCriterium, VectorOfVector3d[i], ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
			g3dPtr->B_genComp(tField); tField++;
		}

		if(SendingIsRequired) OutFieldCompRes(FieldChar, FieldArray, Np, VectInputCell);
		if(FieldArray != 0) 
		{
			delete[] FieldArray;
			FieldArray = 0;
		}
	}
	catch(...) 
	{ 
		if(FieldArray != 0) delete[] FieldArray;
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeField(int ElemKey, char* FieldChar, double** Points, long Np)
{
	radTField* FieldArray = 0;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		radTFieldKey FieldKey;
		if(!ValidateFieldChar(FieldChar, &FieldKey)) return;

		FieldArray = new radTField[Np];
		if(FieldArray == 0) { Send.ErrorMessage("Radia::Error900"); return;}
		radTField* tField = FieldArray;

		TVector3d ZeroVect(0.,0.,0.), v;
		double **tPoints = Points;
		for(long i=0; i<Np; i++)
		{
			double *t = *tPoints;
			v.x = *(t++); v.y = *(t++); v.z = *t;

			*tField = radTField(FieldKey, CompCriterium, v, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
			g3dPtr->B_genComp(tField); 
			tField++; tPoints++;
		}

		if(SendingIsRequired) OutFieldCompRes(FieldChar, FieldArray, Np);
		if(FieldArray != 0) { delete[] FieldArray; FieldArray = 0;}
	}
	catch(...) 
	{ 
		if(FieldArray != 0) delete[] FieldArray;
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldInt(int ElemKey, char* FinOrInfChar, char* FieldIntChar, double* StPoi, long lenStPoi, double* FiPoi, long lenFiPoi)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		TVector3d StPoiVect;
		if(!ValidateVector3d(StPoi, lenStPoi, &StPoiVect)) return;
		TVector3d FiPoiVect;
		if(!ValidateVector3d(FiPoi, lenFiPoi, &FiPoiVect)) return;

		radTFieldKey FieldKey;
		if(!ValidateFieldIntChar(FieldIntChar, FinOrInfChar, &FieldKey)) return;

		TVector3d ZeroVect(0.,0.,0.);
		radTField Field(FieldKey, CompCriterium, StPoiVect, FiPoiVect, ZeroVect, ZeroVect);

		g3dPtr->B_genComp(&Field);

		if(SendingIsRequired) Send.OutFieldIntCompRes(FieldIntChar, &Field);
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldForce(int SourceElemKey, int ShapeElemKey)
{
	try
	{
		radThg hSource;
		if(!ValidateElemKey(SourceElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radThg hShape;
		if(!ValidateElemKey(ShapeElemKey, hShape)) return;
		radTg3d* ShapePtr = Cast.g3dCast(hShape.rep);
		if(ShapePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		radTRectangle* RectanglePtr = Cast.RectangleCast(ShapePtr); 
		if(RectanglePtr==0)
		{
			radTg3dRelax* g3dRelaxPtr = Cast.g3dRelaxCast(ShapePtr);
			if(g3dRelaxPtr!=0)
			{
				radTRecMag* RecMagPtr = Cast.RecMagCast(g3dRelaxPtr);
				if(RecMagPtr==0) { Send.ErrorMessage("Radia::Error036"); return;}
			}
			else { Send.ErrorMessage("Radia::Error036"); return;}
			// Modify this later (incl. Error message), as integration methods for other primitives are ready
			// What about Group Shape?
		}

		radTFieldKey FieldKey;
		FieldKey.Force_= 1;

		radTStructForShapeInt ShapeIntData;
		ShapeIntData.HandleOfSource = hSource;
		ShapeIntData.HandleOfShape = hShape;
		ShapeIntData.IntegrandLength = 1; // Number of elements in TVector3d* to be integrated over a Shape
		ShapeIntData.IntegrandFunPtr = &radTg3d::NormStressTensor;
		ShapeIntData.IntOverLine_= ShapeIntData.IntOverVol_= 0;
		ShapeIntData.IntOverSurf_= 1;
		ShapeIntData.AbsPrecArray = &(CompCriterium.AbsPrecForce);
		TVector3d LocForce(0.,0.,0.);
		ShapeIntData.VectArray = &LocForce;
		char LocForceType = 'r'; // 'a' - axial, 'r' - regular
		ShapeIntData.VectTypeArray = &LocForceType;

		TVector3d ZeroVect(0.,0.,0.);
		radTField Field(FieldKey, CompCriterium, &ShapeIntData /*, ZeroVect */);

		ShapePtr->B_genComp(&Field);

		if(SendingIsRequired) Send.Vector3d(&LocForce);
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldEnergy(int DestElemKey, int SourceElemKey, int* SubdivArray, long lenSubdivArray)
{
	try
	{
		radThg hDest;
		if(!ValidateElemKey(DestElemKey, hDest)) return;
		radTg3d* DestPtr = Cast.g3dCast(hDest.rep);
		if(DestPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radThg hSource;
		if(!ValidateElemKey(SourceElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radTFieldKey FieldKey;
		FieldKey.Energy_ = 1;

		if((lenSubdivArray != 3) || (SubdivArray[0] < 0) || (SubdivArray[1] < 0) || (SubdivArray[2] < 0))
		{
			Send.ErrorMessage("Radia::Error021"); return;
		}
		//double ActualSubdivisionArray[] = {SubdivArray[0], 1., SubdivArray[1], 1., SubdivArray[2], 1.};
		double ActualSubdivisionArray[] = {(double)SubdivArray[0], 1., (double)SubdivArray[1], 1., (double)SubdivArray[2], 1.}; //OC101015

		radTStructForEnergyForceTorqueComp* StructForEnergyForceTorqueCompPtr = new radTStructForEnergyForceTorqueComp();
		StructForEnergyForceTorqueCompPtr->hSource = hSource;
		StructForEnergyForceTorqueCompPtr->hDest = hDest;
		StructForEnergyForceTorqueCompPtr->radPtr = this;
		StructForEnergyForceTorqueCompPtr->DestSubdivArray = ActualSubdivisionArray;
		StructForEnergyForceTorqueCompPtr->AutoDestSubdivision = CheckForAutoDestSubdivision(ActualSubdivisionArray);
		radTHandleStructForEnergyForceTorqueComp HandleStructForEnergyForceTorqueComp(StructForEnergyForceTorqueCompPtr);

		radTField Field(FieldKey, CompCriterium, HandleStructForEnergyForceTorqueComp);
		DestPtr->EnergyForceTorqueComp(&Field);

		if(Field.HandleEnergyForceTorqueCompData.rep->SomethingIsWrong) return;
		if(SendingIsRequired) Send.Double(Field.Energy);
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldForceThroughEnergy(int DestElemKey, int SourceElemKey, char* ForceComponID, int* SubdivArray, long lenSubdivArray)
{
	try
	{
		radThg hDest;
		if(!ValidateElemKey(DestElemKey, hDest)) return;
		radTg3d* DestPtr = Cast.g3dCast(hDest.rep);
		if(DestPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radThg hSource;
		if(!ValidateElemKey(SourceElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		if(!ValidateForceChar(ForceComponID)) return;

		radTFieldKey FieldKey;
		FieldKey.ForceEnr_ = 1;

		if((lenSubdivArray != 3) || (SubdivArray[0] < 0) || (SubdivArray[1] < 0) || (SubdivArray[2] < 0))
		{
			Send.ErrorMessage("Radia::Error021"); return;
		}
		//double ActualSubdivisionArray[] = {SubdivArray[0], 1., SubdivArray[1], 1., SubdivArray[2], 1.};
		double ActualSubdivisionArray[] = {(double)SubdivArray[0], 1., (double)SubdivArray[1], 1., (double)SubdivArray[2], 1.}; //OC101015

		radTStructForEnergyForceTorqueComp* StructForEnergyForceTorqueCompPtr = new radTStructForEnergyForceTorqueComp();
		StructForEnergyForceTorqueCompPtr->hSource = hSource;
		StructForEnergyForceTorqueCompPtr->hDest = hDest;
		StructForEnergyForceTorqueCompPtr->radPtr = this;
		StructForEnergyForceTorqueCompPtr->DestSubdivArray = ActualSubdivisionArray;
		StructForEnergyForceTorqueCompPtr->AutoDestSubdivision = CheckForAutoDestSubdivision(ActualSubdivisionArray);
		radTHandleStructForEnergyForceTorqueComp HandleStructForEnergyForceTorqueComp(StructForEnergyForceTorqueCompPtr);

		radTField Field(FieldKey, CompCriterium, HandleStructForEnergyForceTorqueComp);
		DestPtr->EnergyForceTorqueComp(&Field);

		if(Field.HandleEnergyForceTorqueCompData.rep->SomethingIsWrong) return;
		if(SendingIsRequired) Send.OutFieldForceOrTorqueThroughEnergyCompRes(ForceComponID, Field.Force, 'f');
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldTorqueThroughEnergy(int DestElemKey, int SourceElemKey, char* TorqueComponID, int* SubdivArray, long lenSubdivArray, double* TorqueCenPo, long lenTorqueCenPo)
{
	try
	{
		radThg hDest;
		if(!ValidateElemKey(DestElemKey, hDest)) return;
		radTg3d* DestPtr = Cast.g3dCast(hDest.rep);
		if(DestPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radThg hSource;
		if(!ValidateElemKey(SourceElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		if(!ValidateTorqueChar(TorqueComponID)) return;

		radTFieldKey FieldKey;
		FieldKey.Torque_ = 1;

		if((lenSubdivArray != 3) || (SubdivArray[0] < 0) || (SubdivArray[1] < 0) || (SubdivArray[2] < 0))
		{
			Send.ErrorMessage("Radia::Error021"); return;
		}
		//double ActualSubdivisionArray[] = {SubdivArray[0], 1., SubdivArray[1], 1., SubdivArray[2], 1.};
		double ActualSubdivisionArray[] = {(double)SubdivArray[0], 1., (double)SubdivArray[1], 1., (double)SubdivArray[2], 1.}; //OC101015

		radTStructForEnergyForceTorqueComp* StructForEnergyForceTorqueCompPtr = new radTStructForEnergyForceTorqueComp();
		StructForEnergyForceTorqueCompPtr->hSource = hSource;
		StructForEnergyForceTorqueCompPtr->hDest = hDest;
		StructForEnergyForceTorqueCompPtr->radPtr = this;
		StructForEnergyForceTorqueCompPtr->DestSubdivArray = ActualSubdivisionArray;
		StructForEnergyForceTorqueCompPtr->AutoDestSubdivision = CheckForAutoDestSubdivision(ActualSubdivisionArray);
		radTHandleStructForEnergyForceTorqueComp HandleStructForEnergyForceTorqueComp(StructForEnergyForceTorqueCompPtr);

		radTField Field(FieldKey, CompCriterium, HandleStructForEnergyForceTorqueComp);
		if(!ValidateVector3d(TorqueCenPo, lenTorqueCenPo, &(Field.P))) return;

		DestPtr->EnergyForceTorqueComp(&Field);

		if(Field.HandleEnergyForceTorqueCompData.rep->SomethingIsWrong) return;
		if(SendingIsRequired) Send.OutFieldForceOrTorqueThroughEnergyCompRes(TorqueComponID, Field.Torque, 't');
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------
/**
void radTApplication::OutFieldForceOrTorqueThroughEnergyCompRes(char* ForceComponID, TVector3d& Vect, char ID)
{// This is only for Force and Torque!
	char* BufChar = ForceComponID;
	char* EqEmptyStr = (ID=='f')? "FxFyFz" : "TxTyTz";

	char SmallID = ID;
	char CapitalID = (SmallID=='f')? 'F' : 'T';

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			if((((*BufChar==CapitalID) || (*BufChar==SmallID)) && 
			   (*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) ||
			   (*BufChar == 'X') || (*BufChar == 'x') ||
			   (*BufChar == 'Y') || (*BufChar == 'y') ||
			   (*BufChar == 'Z') || (*BufChar == 'z')) ItemCount++;
			BufChar++;
		}
		BufChar = ForceComponID;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}
	if(ItemCount > 1) Send.InitOutList(ItemCount);

	while (*BufChar != '\0') 
	{
		if((*(BufChar)==CapitalID) || (*(BufChar)==SmallID))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) Send.Vector3d(&Vect);
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) Send.Double(Vect.x);
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) Send.Double(Vect.y);
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) Send.Double(Vect.z);
		BufChar++;
	}
}
**/
//-------------------------------------------------------------------------

void radTApplication::ComputeParticleTrajectory(int ElemKey, double Energy, double x0, double dxdy0, double z0, double dzdy0, double y0, double y1, int Np)
{
	double *TrjData = 0;
	try
	{
		radThg hSource;
		if(!ValidateElemKey(ElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		TrjData = new double[Np*5];
		if(TrjData==0) { Send.ErrorMessage("Radia::Error900"); return;}

		short OnPrc = (CompCriterium.AbsPrecTrjCoord > 0.) || (CompCriterium.AbsPrecTrjAngle > 0.);
		double PrecArray[] = { CompCriterium.AbsPrecTrjCoord, CompCriterium.AbsPrecTrjAngle,
			CompCriterium.AbsPrecTrjCoord, CompCriterium.AbsPrecTrjAngle };

		radTPrtclTrj PrtclTrj(Energy, SourcePtr, CompCriterium, OnPrc, PrecArray);

		double F0[] = { x0, dxdy0, z0, dzdy0 };
		PrtclTrj.Tabulation(F0, y0, y1, Np, TrjData);

		int Depth = 2;
		int Dims[] = { 5, Np };
		if(SendingIsRequired) Send.ArbNestedArrays(TrjData, Dims, Depth);

		if(TrjData != 0) delete[] TrjData;
		TrjData = 0;
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFocusPotent(int ElemKey, double* StPoi, long lenStPoi, double* FiPoi, long lenFiPoi, int Np)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hg.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		TVector3d StPoiVect;
		if(!ValidateVector3d(StPoi, lenStPoi, &StPoiVect)) return;
		TVector3d FiPoiVect;
		if(!ValidateVector3d(FiPoi, lenFiPoi, &FiPoiVect)) return;

		if(Np<2) { Send.ErrorMessage("Radia::Error042"); return;}

		radTPrtclTrj PrtclTrj(SourcePtr, CompCriterium);
		double FocPot = PrtclTrj.FocusingPotential(StPoiVect, FiPoiVect, Np);

		if(SendingIsRequired) Send.Double(FocPot);
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

//void radTApplication::ComputeFocusKickPer(int ElemKey, double* P1, double* Nlong, double per, int nper, double* N1, double r1, int np1, double r2, int np2, const char* StrComment, int nharm, int ns, double d1, double d2, const char* strKickUnit, double inEnergyGeV)
void radTApplication::ComputeFocusKickPer(int ElemKey, double* P1, double* Nlong, double per, double nper, double* N1, double r1, int np1, double r2, int np2, const char* StrComment, int nharm, int ns, double d1, double d2, const char* strKickUnit, double inEnergyGeV, const char* strOutFormat)
{
	double *pKickData1=0, *pKickData2=0, *pBtE2Int=0, *pCoordDir1=0, *pCoordDir2=0, *pAuxFlatArr=0;
	char *pStrReport=0;

	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hg.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		TVector3d P1Vect, NlongVect, N1Vect;
		if(!ValidateVector3d(P1, 3, &P1Vect)) return;
		if(!ValidateVector3d(Nlong, 3, &NlongVect)) return;
		if((NlongVect.x == 0) && (NlongVect.y == 0) && (NlongVect.z == 0)) { Send.ErrorMessage("Radia::Error082"); return;}
		if(!ValidateVector3d(N1, 3, &N1Vect)) return;
		if((N1Vect.x == 0) && (N1Vect.y == 0) && (N1Vect.z == 0)) { Send.ErrorMessage("Radia::Error082"); return;}

		double TolPar = 1.e-08;
		if(PracticallyParallel(NlongVect, N1Vect, TolPar)) { Send.ErrorMessage("Radia::Error087"); return;}

		if(per <= 0) { Send.ErrorMessage("Radia::Error076"); return;}
		if(nper <= 0) { Send.ErrorMessage("Radia::Error077"); return;}
		if((r1 <= 0) || (r2 <= 0)) { Send.ErrorMessage("Radia::Error078"); return;}
		if((np1 <= 0) || (np2 <= 0)) { Send.ErrorMessage("Radia::Error079"); return;}
		if(ns < 8) { Send.ErrorMessage("Radia::Error080"); return;}
		if(nharm <= 0) { Send.ErrorMessage("Radia::Error081"); return;}

		//Validate strKickUnit:
		char cKickUnit = 1;		
		if(!strcmp(strKickUnit, "rad")) cKickUnit = 3;
		else if((!strcmp(strKickUnit, "microrad")) || (!strcmp(strKickUnit, "micro-rad"))) cKickUnit = 2;
		else if(strcmp(strKickUnit, "T2m2")) { Send.ErrorMessage("Radia::Error501"); return;}

		//Validate strOutFormat:
		char cOutFormat = 1; //fixed width
		if((!strcmp(strOutFormat, "tab")) || (!strcmp(strOutFormat, "TAB"))) cOutFormat = 2;
		else if(strcmp(strOutFormat, "fix") && strcmp(strOutFormat, "FIX")) { Send.ErrorMessage("Radia::Error502"); return;}

		long NpKick = np1*np2;
		pKickData1 = new double[NpKick]; if(pKickData1==0) { Send.ErrorMessage("Radia::Error900"); return;}
		pKickData2 = new double[NpKick]; if(pKickData2==0) { Send.ErrorMessage("Radia::Error900"); return;}
		pBtE2Int = new double[NpKick]; if(pBtE2Int==0) { Send.ErrorMessage("Radia::Error900"); return;}
		pCoordDir1 = new double[np1]; if(pCoordDir1==0) { Send.ErrorMessage("Radia::Error900"); return;}
		pCoordDir2 = new double[np2]; if(pCoordDir2==0) { Send.ErrorMessage("Radia::Error900"); return;}

		radTPrtclTrj PrtclTrj(SourcePtr, CompCriterium, inEnergyGeV);
		//PrtclTrj.ComputeSecondOrderKickPer(P1Vect, NlongVect, per, nper, N1Vect, r1, np1, r2, np2, nharm, ns, d1, d2, pCoordDir1, pCoordDir2, pKickData1, pKickData2, pBtE2Int);
		PrtclTrj.ComputeSecondOrderKickPer(P1Vect, NlongVect, per, nper, N1Vect, r1, np1, r2, np2, nharm, ns, d1, d2, pCoordDir1, pCoordDir2, pKickData1, pKickData2, pBtE2Int, cKickUnit);

		char EmptyStr[] = "\0";
		//radTPrtclTrj::ComposeStrReportSecondOrderKickPer(StrComment, per*nper, np1, np2, pCoordDir1, pCoordDir2, pKickData1, pKickData2, pBtE2Int, pStrReport, cKickUnit);
		radTPrtclTrj::ComposeStrReportSecondOrderKickPer(StrComment, per*nper, np1, np2, pCoordDir1, pCoordDir2, pKickData1, pKickData2, pBtE2Int, pStrReport, cKickUnit, cOutFormat);

		if(SendingIsRequired) 
		{
#if defined(__MATHEMATICA__)
			Send.InitOutList(6);

			int Depth = 2;
			int Dims[] = { np1, np2 };
			Send.ArbNestedArrays(pKickData1, Dims, Depth);
			Send.ArbNestedArrays(pKickData2, Dims, Depth);
			Send.ArbNestedArrays(pBtE2Int, Dims, Depth); //OC100509
			Send.DoubleList(pCoordDir1, np1);
			Send.DoubleList(pCoordDir2, np2);

			char *pStrToSend = pStrReport;
			if(pStrToSend == 0) pStrToSend = EmptyStr;
			Send.String(pStrToSend);
#else
			//long LenArr = 2*np1*np2 + np1 + np2 + 1;
			long LenArr = 3*np1*np2 + np1 + np2 + 1;
			pAuxFlatArr = new double[LenArr]; 
			double *tAuxFlatArr = pAuxFlatArr;
			double *tKickData1 = pKickData1;
			double *tKickData2 = pKickData2;
			double *tBtE2Int = pBtE2Int;
			double *tCoordDir1 = pCoordDir1;
			double *tCoordDir2 = pCoordDir2;

			for(long i=0; i<NpKick; i++) *(tAuxFlatArr++) = *(tKickData1++);
			for(long j=0; j<NpKick; j++) *(tAuxFlatArr++) = *(tKickData2++);
			for(long ii=0; ii<NpKick; ii++) *(tAuxFlatArr++) = *(tBtE2Int++);

			for(long k=0; k<np1; k++) *(tAuxFlatArr++) = *(tCoordDir1++);
			for(long m=0; m<np2; m++) *(tAuxFlatArr++) = *(tCoordDir2++);

			*tAuxFlatArr += strlen(pStrReport);

			Send.DoubleList(pAuxFlatArr, LenArr);
			
			if(pAuxFlatArr != 0) { delete[] pAuxFlatArr; pAuxFlatArr = 0;}
#endif
		}

		if(pKickData1 != 0) { delete[] pKickData1; pKickData1 = 0;}
		if(pKickData2 != 0) { delete[] pKickData2; pKickData2 = 0;}
		if(pBtE2Int != 0) { delete[] pBtE2Int; pBtE2Int = 0;}
		if(pCoordDir1 != 0) { delete[] pCoordDir1; pCoordDir1 = 0;}
		if(pCoordDir2 != 0) { delete[] pCoordDir2; pCoordDir2 = 0;}
		if(pStrReport != 0) { delete[] pStrReport; pStrReport = 0;}
		if(pAuxFlatArr != 0) { delete[] pAuxFlatArr; pAuxFlatArr = 0;}
	}
	catch(...) 
	{ 
		if(pKickData1 != 0) delete[] pKickData1;
		if(pKickData2 != 0) delete[] pKickData2;
		if(pBtE2Int != 0) delete[] pBtE2Int;
		if(pCoordDir1 != 0) delete[] pCoordDir1;
		if(pCoordDir2 != 0) delete[] pCoordDir2;
		if(pStrReport != 0) delete[] pStrReport;
		if(pAuxFlatArr != 0) delete[] pAuxFlatArr;
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComposeFocusKickPerFormStrRep(double* pKickData1, double* pKickData2, double* pBtE2Int, double* pCoordDir1, double* pCoordDir2, int np1, int np2, double per, int nper, const char* StrComment)
{//Auxiliary function, added for DLL
	char EmptyStr[] = "\0";
	char *pStrReport=0;
	radTPrtclTrj::ComposeStrReportSecondOrderKickPer(StrComment, per*nper, np1, np2, pCoordDir1, pCoordDir2, pKickData1, pKickData2, pBtE2Int, pStrReport);

	char *pStrToSend = pStrReport;
	if(pStrToSend == 0) pStrToSend = EmptyStr;
	Send.String(pStrToSend);

	if(pStrReport != 0) delete[] pStrReport;
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFocusKick(int ElemKey, double* P1, double* Nlong, double* ArrLongDist, int lenArrLongDist, int ns, double* Ntr1, double r1, int np1, double r2, int np2, const char* StrComment, double d1, double d2)
{
	double *pKickData1=0, *pKickData2=0, *pIBe2=0, *pCoordDir1=0, *pCoordDir2=0;
	char *pStrReport=0;

	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hg.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		TVector3d P1Vect, NlongVect, N1Vect;
		if(!ValidateVector3d(P1, 3, &P1Vect)) return;
		if(!ValidateVector3d(Nlong, 3, &NlongVect)) return;
		if((NlongVect.x == 0) && (NlongVect.y == 0) && (NlongVect.z == 0)) { Send.ErrorMessage("Radia::Error082"); return;}
		if(!ValidateVector3d(Ntr1, 3, &N1Vect)) return;
		if((N1Vect.x == 0) && (N1Vect.y == 0) && (N1Vect.z == 0)) { Send.ErrorMessage("Radia::Error082"); return;}

		double TolPar = 1.e-08;
		if(PracticallyParallel(NlongVect, N1Vect, TolPar)) { Send.ErrorMessage("Radia::Error087"); return;}

		if((ArrLongDist == 0) || (lenArrLongDist <= 0)) { Send.ErrorMessage("Radia::Error083"); return;}
		if(lenArrLongDist >= 50) { Send.ErrorMessage("Radia::Error085"); return;}
		double PrevLongDist = *ArrLongDist;
		double *tArrLongDist = ArrLongDist + 1;
		for(int k=1; k<lenArrLongDist; k++)
		{
			double NewLongDist = *(tArrLongDist++);
			if(NewLongDist <= PrevLongDist) { Send.ErrorMessage("Radia::Error086"); return;}
			PrevLongDist = NewLongDist;
		}

		if(ns <= 2) { Send.ErrorMessage("Radia::Error042"); return;}
		if((r1 <= 0) || (r2 <= 0)) { Send.ErrorMessage("Radia::Error078"); return;}
		if((np1 <= 0) || (np2 <= 0)) { Send.ErrorMessage("Radia::Error079"); return;}
		if((np1 >= 100) || (np2 >= 100)) { Send.ErrorMessage("Radia::Error084"); return;}

		int TotAmOfPairsOfMatr = lenArrLongDist;
		if(TotAmOfPairsOfMatr > 1) TotAmOfPairsOfMatr++;
		long NpKick = np1*np2*TotAmOfPairsOfMatr;
		pKickData1 = new double[NpKick]; if(pKickData1==0) { Send.ErrorMessage("Radia::Error900"); return;}
		pKickData2 = new double[NpKick]; if(pKickData2==0) { Send.ErrorMessage("Radia::Error900"); return;}
		pIBe2 = new double[NpKick]; if(pIBe2==0) { Send.ErrorMessage("Radia::Error900"); return;}
		pCoordDir1 = new double[np1]; if(pCoordDir1==0) { Send.ErrorMessage("Radia::Error900"); return;}
		pCoordDir2 = new double[np2]; if(pCoordDir2==0) { Send.ErrorMessage("Radia::Error900"); return;}

		radTPrtclTrj PrtclTrj(SourcePtr, CompCriterium);
		PrtclTrj.ComputeSecondOrderKick(P1Vect, NlongVect, ArrLongDist, lenArrLongDist, ns, N1Vect, r1, np1, r2, np2, d1, d2, pCoordDir1, pCoordDir2, pKickData1, pKickData2, pIBe2);

		char EmptyStr[] = "\0";
		PrtclTrj.ComposeStrReportSecondOrderKick(StrComment, ArrLongDist, lenArrLongDist, np1, np2, pCoordDir1, pCoordDir2, pKickData1, pKickData2, pIBe2, pStrReport);

		if(SendingIsRequired) 
		{
			Send.InitOutList(4);

			Send.InitOutList(TotAmOfPairsOfMatr);
			int Depth = 2;
			int Dims[] = { np1, np2 };
			double *tKickData1 = pKickData1;
			double *tKickData2 = pKickData2;
			double *tIBe2 = pIBe2;
			long dataPer = np1*np2;
			for(int i=0; i<TotAmOfPairsOfMatr; i++)
			{
				//Send.InitOutList(2);
				Send.InitOutList(3);
				Send.ArbNestedArrays(tKickData1, Dims, Depth);
				Send.ArbNestedArrays(tKickData2, Dims, Depth);
				Send.ArbNestedArrays(tIBe2, Dims, Depth);
				tKickData1 += dataPer;
				tKickData2 += dataPer;
				tIBe2 += dataPer;
			}

			Send.DoubleList(pCoordDir1, np1);
			Send.DoubleList(pCoordDir2, np2);

			char *pStrToSend = pStrReport;
			if(pStrToSend == 0) pStrToSend = EmptyStr;
			Send.String(pStrToSend);
		}

		if(pKickData1 != 0) { delete[] pKickData1; pKickData1 = 0;}
		if(pKickData2 != 0) { delete[] pKickData2; pKickData2 = 0;}
		if(pIBe2 != 0) { delete[] pIBe2; pIBe2 = 0;}
		if(pCoordDir1 != 0) { delete[] pCoordDir1; pCoordDir1 = 0;}
		if(pCoordDir2 != 0) { delete[] pCoordDir2; pCoordDir2 = 0;}
		if(pStrReport != 0) { delete[] pStrReport; pStrReport = 0;}
	}
	catch(...) 
	{ 
		if(pKickData1 != 0) { delete[] pKickData1; pKickData1 = 0;}
		if(pKickData2 != 0) { delete[] pKickData2; pKickData2 = 0;}
		if(pIBe2 != 0) { delete[] pIBe2; pIBe2 = 0;}
		if(pCoordDir1 != 0) { delete[] pCoordDir1; pCoordDir1 = 0;}
		if(pCoordDir2 != 0) { delete[] pCoordDir2; pCoordDir2 = 0;}
		if(pStrReport != 0) { delete[] pStrReport; pStrReport = 0;}

		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeShimSignature(int ElemKey, char* FldID, double* Disp, double* StPoi, double* FiPoi, int Np, double* Vi)
{
	radTField **arr_pField = 0, *arr_resField = 0;
	int TwoNp = Np << 1;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hg.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radTFieldKey FieldKey;
		if(FldID == 0) { Send.ErrorMessage("Radia::Error096"); return;}
		if((*FldID == 'i') || (*FldID == 'I'))
		{
			//if(!ValidateFieldIntChar(FldID, "inf", &FieldKey)) 
			if(!ValidateFieldIntChar(FldID, (char*)"inf", &FieldKey)) //OC01052013 
			{
				Send.ErrorMessage("Radia::Error096"); return;
			}
		}
		else if(!ValidateFieldChar(FldID, &FieldKey, false))
		{
			Send.ErrorMessage("Radia::Error096"); return;
		}

		TVector3d vDisp;
		if(!ValidateVector3d(Disp, 3, &vDisp)) return;
		TVector3d vStPoi;
		if(!ValidateVector3d(StPoi, 3, &vStPoi)) return;
		TVector3d vFiPoi;
		if(!ValidateVector3d(FiPoi, 3, &vFiPoi)) return;
		TVector3d vVi;
		if(!ValidateVector3d(Vi, 3, &vVi)) return;

		if(Np <= 0) { Send.ErrorMessage("Radia::Error079"); return;}

		TVector3d ZeroVect(0.,0.,0.);
		arr_pField = new radTField*[TwoNp];
		arr_resField = new radTField[Np];

        radTField **tField = arr_pField;
		for(int s=0; s<TwoNp; s++) *(tField++) = 0;

		TVector3d vTransl = ((Np > 1)? (1./double(Np - 1)) : 1)*(vFiPoi - vStPoi);

		short prevSendingIsRequired = SendingIsRequired;
		SendingIsRequired = false;

		int indTranslat = SetTranslation(Disp, 3);
		int indSrcDisp = DuplicateElement_g3d(ElemKey, 0, 0, 0);
		indSrcDisp = ApplySymmetry(indSrcDisp, indTranslat, 1);
		radThg hgSrcDisp;
		if(!ValidateElemKey(indSrcDisp, hgSrcDisp)) return;
		radTg3d* SourceDispPtr = Cast.g3dCast(hgSrcDisp.rep);
		if(SourceDispPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radTg3d *arrSourceDispPtr[] = {SourcePtr, SourceDispPtr};

        tField = arr_pField;
		for(int k=0; k<2; k++)
		{
            radTg3d *curSourceDispPtr = arrSourceDispPtr[k];
            TVector3d vP1 = vStPoi;
			for(int i=0; i<Np; i++)
			{
				if(FieldKey.Ib_	|| FieldKey.Ih_)
				{
					TVector3d vP2 = vP1 + vVi;
                    *tField = new radTField(FieldKey, CompCriterium, vP1, vP2, ZeroVect, ZeroVect);
				}
				else
				{
					*tField = new radTField(FieldKey, CompCriterium, vP1, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
				}

				curSourceDispPtr->B_genComp(*tField);

				if(k == 1) arr_resField[i] = **tField - *(arr_pField[i]);

				tField++;
                vP1 += vTransl;
			}
			//vStPoi += vDisp;
			curSourceDispPtr++;
		}

		DeleteElement(indSrcDisp);
		DeleteElement(indTranslat);

		SendingIsRequired = prevSendingIsRequired;

		if(SendingIsRequired)
		{
			if(FieldKey.Ib_	|| FieldKey.Ih_)
			{
				Send.OutFieldIntCompRes(FldID, arr_resField, NULL, Np);
				//Send.OutFieldIntCompRes(FldID, &Field);
			}
			else
			{
                Send.OutFieldCompRes(FldID, arr_resField, NULL, Np);
			}
		}

		if(arr_pField != 0) 
		{
			for(int i=0; i<TwoNp; i++) if(arr_pField[i] != 0) delete arr_pField[i];
			delete[] arr_pField;
		}
		if(arr_resField != 0) delete[] arr_resField;
	}
	catch(...) 
	{ 
		if(arr_pField != 0) 
		{
			for(int i=0; i<TwoNp; i++) if(arr_pField[i] != 0) delete arr_pField[i];
			delete[] arr_pField;
		}
		if(arr_resField != 0) delete[] arr_resField;
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetAndShowPhysUnits()
{
	try
	{
		// Implement mechanism for setting different units in future
		char LengthUnitID[] = "mm";
		char AngleUnitID[] = "radian";
		char CurrentUnitID[] = "A";
		char CurrentDensityUnitID[] = "A/mm^2";
		char MagnetizationUnitID[] = "Tesla";
		char FieldInductionUnitID[] = "Tesla";
		char FieldStrengthUnitID[] = "Tesla";
		char FieldForceUnitID[] = "Newton";

//#ifdef __GCC__
//		ostrstream OutUnitsInfoStream;
//#else
		ostringstream OutUnitsInfoStream; // Porting
//#endif

		OutUnitsInfoStream << "Physical units currently in use:" << endl;
		OutUnitsInfoStream << "Length:  " << LengthUnitID << endl;
		OutUnitsInfoStream << "Angle:  " << AngleUnitID << endl;
		OutUnitsInfoStream << "Electric current:  " << CurrentUnitID << endl;
		OutUnitsInfoStream << "Electric current density:  " << CurrentDensityUnitID << endl;
		OutUnitsInfoStream << "Magnetization:  " << MagnetizationUnitID << endl;
		OutUnitsInfoStream << "Magnetic field induction:  " << FieldInductionUnitID << endl;
		OutUnitsInfoStream << "Magnetic field strength:  " << FieldInductionUnitID << endl;
		OutUnitsInfoStream << "Magnetic field integral along a line:  " << FieldInductionUnitID << " " << LengthUnitID << endl;
		OutUnitsInfoStream << "Magnetic field vector-potential:  " << FieldInductionUnitID << " " << LengthUnitID << endl;
		OutUnitsInfoStream << "Force:  " << FieldForceUnitID << endl;

		OutUnitsInfoStream << ends;

//#ifdef __GCC__
//		if(SendingIsRequired) Send.String(OutUnitsInfoStream.str());
//#else
		if(SendingIsRequired) Send.String(OutUnitsInfoStream.str().c_str()); // Porting
//#endif

		return 1;
	}
	catch(...) 
	{ 
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTApplication::ReplaceInAllGroups(radThg& OldHandle, radThg& NewHandle)
{
	radTg3d* g3dOldPtr = (radTg3d*)(OldHandle.rep);

	try
	{
		for(radTmhg::iterator GenIter = GlobalMapOfHandlers.begin(); GenIter != GlobalMapOfHandlers.end(); ++GenIter)
		{
			radTg* gP = ((*GenIter).second).rep;
			radTg3d* g3dP = Cast.g3dCast(gP); 
			if(g3dP != 0) 
			{
				radTGroup* GroupP = Cast.GroupCast(g3dP); 
				if(GroupP != 0) 
					for(radTmhg::iterator GroupIter = GroupP->GroupMapOfHandlers.begin();
						GroupIter != GroupP->GroupMapOfHandlers.end(); ++GroupIter)
						if(((*GroupIter).second).rep == g3dOldPtr) (*GroupIter).second = NewHandle;
			}
		}
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ReplaceInGlobalMap(radThg& OldHandle, radThg& NewHandle)
{
	radTg3d* g3dOldPtr = (radTg3d*)(OldHandle.rep);

	try
	{
		for(radTmhg::iterator GenIter = GlobalMapOfHandlers.begin(); GenIter != GlobalMapOfHandlers.end(); ++GenIter)
		{
			if(((*GenIter).second).rep == g3dOldPtr) 
			{
				(*GenIter).second = NewHandle;
			}
		}
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SubdivideElement_g3d(int ElemKey, double* SubdivArray, long lenSubdivArray, char TypeExtraSpec, double* ExtraSpec, long lenExtraSpec, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	radTCylindricSubdivSpec* pCylindricSubdivSpec = 0;

	try
	{
		radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
		if(iter == GlobalMapOfHandlers.end()) { Send.ErrorMessage("Radia::Error002"); return 0;}
		radTg3d* g3dPtr = Cast.g3dCast((*iter).second.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		if(lenSubdivArray != 6) { Send.ErrorMessage("Radia::Error052"); return 0;}
		if((SubdivArray[0] <= 0) || (SubdivArray[1] <= 0) || (SubdivArray[2] <= 0)
			|| (SubdivArray[3] <= 0) || (SubdivArray[4] <= 0) || (SubdivArray[5] <= 0))
		{
			Send.ErrorMessage("Radia::Error053"); return 0;
		}

		char SubdivisionFrame = 0; // 0- Local; 1- Laboratory, each Group Member separately; 2- Laboratory, all Group as whole;
		char SubdivisionParamCode = 0; // 0- kx,ky,kz are subdiv. numbers; 1- kx,ky,kz are average sizes of pieces;
		char SubdivideCoils = 0; // 0- do not subdivide coils; 1- subdivide coils;

		radTOptionNames OptNam;
		const char** BufNameString = OptionNames;
		const char** BufValString = OptionValues;

		for(int i=0; i<OptionCount; i++)
		{
			if(!strcmp(*BufNameString, OptNam.Frame))
			{
				if(!strcmp(*BufValString, (OptNam.FrameValues)[0])) SubdivisionFrame = 0;
				else if(!strcmp(*BufValString, (OptNam.FrameValues)[1])) SubdivisionFrame = 1;
				else if(!strcmp(*BufValString, (OptNam.FrameValues)[2])) SubdivisionFrame = 2;
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else if(!strcmp(*BufNameString, OptNam.SubdParamCode))
			{
				if(!strcmp(*BufValString, (OptNam.SubdParamCodeValues)[0])) SubdivisionParamCode = 0;
				else if(!strcmp(*BufValString, (OptNam.SubdParamCodeValues)[1])) SubdivisionParamCode = 1;
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else if(!strcmp(*BufNameString, OptNam.SubdCoils))
			{
				if((!strcmp(*BufValString, (OptNam.SubdCoilsValues)[0])) || (!strcmp(*BufValString, (OptNam.SubdCoilsValues)[2]))) SubdivideCoils = 0;
				else if((!strcmp(*BufValString, (OptNam.SubdCoilsValues)[1])) || (!strcmp(*BufValString, (OptNam.SubdCoilsValues)[3]))) SubdivideCoils = 1;
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
			BufNameString++; BufValString++;
		}

		if(SubdivisionParamCode == 0)
		{// Important: setting q=1 if k=1
			const double SubdZeroTol = 1.E-12;
			for(int kk=0; kk<3; kk++)
			{
				int TwoKk = kk*2;
				if(SubdivArray[TwoKk] < 1.) SubdivArray[TwoKk] = 1.;
				if(fabs(SubdivArray[TwoKk]-1.) < SubdZeroTol) SubdivArray[TwoKk+1] = 1.;
			}
		}

		double NewSubdivArray[15];
		for(int ii=0; ii<lenSubdivArray; ii++) NewSubdivArray[ii] = SubdivArray[ii];

		radTSubdivOptions SubdivOptions;
		SubdivOptions.SubdivisionFrame = SubdivisionFrame;
		SubdivOptions.SubdivisionParamCode = SubdivisionParamCode;
		SubdivOptions.SubdivideCoils = SubdivideCoils;
		SubdivOptions.PutNewStuffIntoGenCont = 1;
		SubdivOptions.ReplaceOldStuff = 1;

		char CylindricSubdivision = (TypeExtraSpec == 1)? 1 : 0;
		if(CylindricSubdivision)
		{
			if((ExtraSpec[3]==0.) && (ExtraSpec[4]==0.) && (ExtraSpec[5]==0.)) { Send.ErrorMessage("Radia::Error066"); return 0;}
			if((lenExtraSpec > 6) && (ExtraSpec[9]<=0.)) { Send.ErrorMessage("Radia::Error067"); return 0;}

			pCylindricSubdivSpec = new radTCylindricSubdivSpec(ExtraSpec, lenExtraSpec);
			if(pCylindricSubdivSpec == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

			SubdivOptions.MethForRadialSegmAtEllCylSubd = 0; // Make accessible from outside if necessary
		}

		char SubdivisionByParPlanes = (TypeExtraSpec == 2)? 1 : 0;
		int AmOfDir;
		if(SubdivisionByParPlanes)
		{
			AmOfDir = int((lenExtraSpec + 1.E-06)/3.);
			if(AmOfDir < 1) { Send.ErrorMessage("Radia::Error069"); return 0;}

			for(int i=0; i<AmOfDir; i++)
			{
				int im5 = i*5;
				int im2 = i*2;
				int im3 = i*3;
				NewSubdivArray[im5] = ExtraSpec[im3];
				NewSubdivArray[im5+1] = ExtraSpec[im3+1];
				NewSubdivArray[im5+2] = ExtraSpec[im3+2];

				NewSubdivArray[im5+3] = SubdivArray[im2];
				NewSubdivArray[im5+4] = SubdivArray[im2+1];

				if((NewSubdivArray[im5]==0.) && (NewSubdivArray[im5+1]==0.) && (NewSubdivArray[im5+2]==0.)) { Send.ErrorMessage("Radia::Error069"); return 0;}
				if((NewSubdivArray[im5+3]<=0.) || (NewSubdivArray[im5+4]<=0.)) { Send.ErrorMessage("Radia::Error053"); return 0;}
			}
			//double NewSubdivArray[] = {n1x,n1y,n1z,k1x,q1x, n2x,n2y,n2z,k2x,q2x, n3x,n3y,n3z,k3x,q3x};
		}

		radThg& hgNew = (*iter).second;
		radThg hgOld = hgNew;

		if(CylindricSubdivision)
		{
			if(!g3dPtr->SubdivideItselfByEllipticCylinder(NewSubdivArray, pCylindricSubdivSpec, hgNew, this, &SubdivOptions)) return 0;
		}
		else if(SubdivisionByParPlanes)
		{
			if(!g3dPtr->SubdivideItselfByParPlanes(NewSubdivArray, AmOfDir, hgNew, this, &SubdivOptions)) return 0;
		}
		else if(!g3dPtr->SubdivideItself(NewSubdivArray, hgNew, this, &SubdivOptions)) return 0;

		if(g3dPtr->IsGroupMember) ReplaceInAllGroups(hgOld, hgNew);
		if(SendingIsRequired) Send.Int(ElemKey);

		if(pCylindricSubdivSpec != 0) { delete pCylindricSubdivSpec; pCylindricSubdivSpec = 0;}
		return ElemKey;
	}
	catch(...)
	{
        if(pCylindricSubdivSpec != 0) delete pCylindricSubdivSpec;
        Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SubdivideElement_g3dByParPlanes(int ElemKey, double* SubdivArray, int AmOfSubdivDirections, const char* LocOrLabFrame)
{
	try
	{
		radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
		if(iter == GlobalMapOfHandlers.end()) { Send.ErrorMessage("Radia::Error002"); return 0;}
		radTg3d* g3dPtr = Cast.g3dCast((*iter).second.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		char SubdivisionFrame = 0; // 0- Local; 1- Laboratory;
		char SubdivisionParamCode = 0; // 0- kx,ky,kz are subdiv. numbers; 1- kx,ky,kz are average sizes of pieces;
		char SubdivideCoils = 0; // 0- do not subdivide coils; 1- subdivide coils;

		if((!strcmp(LocOrLabFrame, "loc")) || (!strcmp(LocOrLabFrame, "Loc")) || (!strcmp(LocOrLabFrame, "LOC"))) SubdivisionFrame = 0;
		else if((!strcmp(LocOrLabFrame, "lab")) || (!strcmp(LocOrLabFrame, "Lab")) || (!strcmp(LocOrLabFrame, "LAB"))) SubdivisionFrame = 1;
		else { Send.ErrorMessage("Radia::Error054"); return 0;}

		if(SubdivisionParamCode == 0)
		{// Important: setting q=1 if k=1
			const double SubdZeroTol = 1.E-12;
			for(int kk=0; kk<AmOfSubdivDirections; kk++)
			{
				int BufInd = 5*kk+4;
				if(fabs(SubdivArray[BufInd] - 1.) < SubdZeroTol) SubdivArray[BufInd] = 1.;
			}
		}

		radTSubdivOptions SubdivOptions;
		SubdivOptions.SubdivisionFrame = SubdivisionFrame;
		SubdivOptions.SubdivisionParamCode = SubdivisionParamCode;
		SubdivOptions.SubdivideCoils = SubdivideCoils;
		SubdivOptions.PutNewStuffIntoGenCont = 1;

		radThg& hgNew = (*iter).second;
		radThg hgOld = hgNew;

		int SubdOK = g3dPtr->SubdivideItselfByParPlanes(SubdivArray, AmOfSubdivDirections, hgNew, this, &SubdivOptions);
		if(!SubdOK) return 0;
		if(g3dPtr->IsGroupMember) ReplaceInAllGroups(hgOld, hgNew);

		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::CutElement_g3d(int ElemKey, double* PointOnPlane, long lenPointOnPlane, double* PlaneNormal, long lenPlaneNormal, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	try
	{
		radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
		if(iter == GlobalMapOfHandlers.end()) { Send.ErrorMessage("Radia::Error002"); return 0;}
		radTg3d* g3dPtr = Cast.g3dCast((*iter).second.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		TVector3d PointOnPlaneVect, PlaneNormalVect;
		if(!ValidateVector3d(PointOnPlane, lenPointOnPlane, &PointOnPlaneVect)) return 0;
		if(!ValidateVector3d(PlaneNormal, lenPlaneNormal, &PlaneNormalVect)) return 0;

		TVector3d PlaneSpecification[] = {PointOnPlaneVect, PlaneNormalVect};

		radTOptionNames OptNam;
		const char* OptNamesToFind[] = {OptNam.Frame, OptNam.SubdCoils};
		char OptValsFoundParsed[] = {0, 0};
		char &SubdivisionFrame = OptValsFoundParsed[0]; // 0- Local; 1- LabTot; 2- Lab;
		char &CutCoils = OptValsFoundParsed[1]; // 0- No; 1- Yes;

		if(!OptNam.findParseOptionValues(OptionNames, OptionValues, OptionCount, OptNamesToFind, 2, OptValsFoundParsed, 0, 0))
		{
			Send.ErrorMessage("Radia::Error062"); return 0;
		}

		//char SubdivisionFrame = 0; // 0- Local; 1- Laboratory;
		//const char** BufNameString = OptionNames;
		//const char** BufValString = OptionValues;
		//for(int i=0; i<OptionCount; i++)
		//{
		//	if(!strcmp(*BufNameString, OptNam.Frame))
		//	{
		//		if(!strcmp(*BufValString, (OptNam.FrameValues)[0])) SubdivisionFrame = 0;
		//		else if(!strcmp(*BufValString, (OptNam.FrameValues)[1])) SubdivisionFrame = 1;
		//		else if(!strcmp(*BufValString, (OptNam.FrameValues)[2])) SubdivisionFrame = 2;
		//		else { Send.ErrorMessage("Radia::Error062"); return 0;}
		//	}
		//	else { Send.ErrorMessage("Radia::Error062"); return 0;}
		//	BufNameString++; BufValString++;
		//}

		radTSubdivOptions SubdivOptions;
		SubdivOptions.SubdivisionFrame = SubdivisionFrame;
		SubdivOptions.SubdivisionParamCode = 0; // Is not used by Cut
		SubdivOptions.SubdivideCoils = CutCoils; //0;
		SubdivOptions.PutNewStuffIntoGenCont = 1;
		SubdivOptions.ReplaceOldStuff = 0;
		SubdivOptions.SeparatePiecesAtCutting = 1;
		SubdivOptions.MapInternalFacesAfterCut = 0;

		radThg hgOld = (*iter).second;
		radThg hgNewTot = hgOld;
		radTPair_int_hg NewLowerPair_int_hg, NewUpperPair_int_hg;

		if(!g3dPtr->CutItself(PlaneSpecification, hgNewTot, NewLowerPair_int_hg, NewUpperPair_int_hg, this, &SubdivOptions)) return 0;

		int LowerElemKey = 0, UpperElemKey = 0;
		if(NewLowerPair_int_hg.Handler_g.rep != 0)
		{
			LowerElemKey = NewLowerPair_int_hg.m;
			if(LowerElemKey == 0) LowerElemKey = AddElementToContainer(NewLowerPair_int_hg.Handler_g);
		}
		if(NewUpperPair_int_hg.Handler_g.rep != 0) 
		{
			UpperElemKey = NewUpperPair_int_hg.m;
			if(UpperElemKey == 0) UpperElemKey = AddElementToContainer(NewUpperPair_int_hg.Handler_g);
		}

		if(SendingIsRequired) 
		{
			//char LowerElemKeyIsNotZero = (LowerElemKey != 0);
			//char UpperElemKeyIsNotZero = (UpperElemKey != 0);
			//int NumberOfElem = (LowerElemKeyIsNotZero && UpperElemKeyIsNotZero)? 2 : 1;
			int NumberOfElem = 0; //OC290908

			int OutArray[2];
			int *pOutArray = OutArray;
			if(LowerElemKey != 0) { *(pOutArray++) = LowerElemKey; NumberOfElem++;}
			if(UpperElemKey != 0) { *pOutArray = UpperElemKey; NumberOfElem++;}

			if(NumberOfElem <= 0) { *pOutArray = ElemKey; NumberOfElem++;} //OC290908

			Send.IntList(OutArray, NumberOfElem);
		}
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::FieldCompMethForSubdividedRecMag(int ElemKey, int InFldCmpMeth, int SubLevel)
{
	radThg hg;
	if(!ValidateElemKey(ElemKey, hg)) return 0;
	radTg3d* g3dP = Cast.g3dCast(hg.rep); if(g3dP==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
	radTGroup* GroupP = Cast.GroupCast(g3dP); if(GroupP==0) { Send.ErrorMessage("Radia::Error004"); return 0;}
	radTSubdividedRecMag* SubdividedRecMagP = Cast.SubdividedRecMagCast(GroupP); 
	if(SubdividedRecMagP==0) { Send.ErrorMessage("Radia::Error037"); return 0;}

	int Probe = SubdividedRecMagP->SetupFldCmpData((short)InFldCmpMeth, SubLevel);
	if(Probe==0) { Send.ErrorMessage("Radia::Error040"); return 0;}
	else if(Probe==-38) { Send.ErrorMessage("Radia::Error038"); return 0;}

	if(SendingIsRequired) Send.Int(ElemKey);
	return ElemKey;
}

//-------------------------------------------------------------------------

int radTApplication::SetLocMgnInSbdRecMag(int ElemKey, TVector3d* ArrayOfVectIndx, TVector3d* ArrayOfMagn, int Len)
{
	radThg hg;
	if(!ValidateElemKey(ElemKey, hg)) return 0;
	radTg3d* g3dP = Cast.g3dCast(hg.rep); if(g3dP==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
	radTGroup* GroupP = Cast.GroupCast(g3dP); if(GroupP==0) { Send.ErrorMessage("Radia::Error004"); return 0;}
	radTSubdividedRecMag* SubdividedRecMagP = Cast.SubdividedRecMagCast(GroupP); 
	if(SubdividedRecMagP==0) { Send.ErrorMessage("Radia::Error037"); return 0;}
	
	for(int i=0; i<Len; i++)
	{
		int ix = int(ArrayOfVectIndx[i].x) - 1;
		if((ix >= int(SubdividedRecMagP->kx)) || (ix < 0)) { Send.ErrorMessage("Radia::Error039"); return 0;}
		int iy = int(ArrayOfVectIndx[i].y) - 1;
		if((iy >= int(SubdividedRecMagP->ky)) || (iy < 0)) { Send.ErrorMessage("Radia::Error039"); return 0;}
		int iz = int(ArrayOfVectIndx[i].z) - 1;
		if((iz >= int(SubdividedRecMagP->kz)) || (iz < 0)) { Send.ErrorMessage("Radia::Error039"); return 0;}

		int SubElemNo = (ix*int(SubdividedRecMagP->ky) + iy)*int(SubdividedRecMagP->kz) + iz;

		SubdividedRecMagP->SetLocalMagn(SubElemNo, ArrayOfMagn[i]);
	}

	if(SendingIsRequired) Send.Int(ElemKey);
	return ElemKey;
}

//-------------------------------------------------------------------------

int radTApplication::DuplicateElement_g3d(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		radTrans* TransPtr = 0; //OC270307: added duplication of transformations
		if(g3dPtr==0) 
		{
			TransPtr = Cast.TransCast(hg.rep); 
			if(TransPtr == 0) 
			{
				Send.ErrorMessage("Radia::Error003"); return 0;
			}
		}

		if(TransPtr != 0) //OC270307: added duplication of transformations
		{
			radTrans *pNewTrans = new radTrans(*TransPtr);
			radThg hg(pNewTrans);
			int NewElemKey = AddElementToContainer(hg);
			if(SendingIsRequired) Send.Int(NewElemKey);
			return NewElemKey;
		}

		radTOptionNames OptNam;
		const char** BufNameString = OptionNames;
		const char** BufValString = OptionValues;

		char ReleaseSym = 0;
		for(int i=0; i<OptionCount; i++)
		{
			if(!strcmp(*BufNameString, OptNam.FreeSym))
			{
				if(!strcmp(*BufValString, (OptNam.FreeSymValues)[0])) ReleaseSym = 0;
				else if(!strcmp(*BufValString, (OptNam.FreeSymValues)[1])) ReleaseSym = 1;
				else if(!strcmp(*BufValString, (OptNam.FreeSymValues)[2])) ReleaseSym = 0;
				else if(!strcmp(*BufValString, (OptNam.FreeSymValues)[3])) ReleaseSym = 1;
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
			BufNameString++; BufValString++;
		}

		char PutNewStuffIntoGenCont = 1;
		if(ReleaseSym)
		{
			if(!g3dPtr->CreateFromSym(hg, this, PutNewStuffIntoGenCont)) return 0;
		}
		else
		{
			if(!g3dPtr->DuplicateItself(hg, this, PutNewStuffIntoGenCont)) return 0;
		}

		int NewElemKey = AddElementToContainer(hg);
		CopyDrawAttr(ElemKey, NewElemKey);

		if(SendingIsRequired) Send.Int(NewElemKey);
		return NewElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::CreateFromObj_g3dWithSym(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		char PutNewStuffIntoGenCont = 1;
		if(!g3dPtr->CreateFromSym(hg, this, PutNewStuffIntoGenCont)) return 0;

		if(hg.rep != g3dPtr)
		{
			int NewElemKey = AddElementToContainer(hg);
			CopyDrawAttr(ElemKey, NewElemKey);

			if(SendingIsRequired) Send.Int(NewElemKey);
			return NewElemKey;
		}
		else
		{
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

int radTApplication::DeleteAllElements(int DeletionMethNo)
{
	try
	{
		if(DeletionMethNo == 1) GlobalMapOfHandlers.erase(GlobalMapOfHandlers.begin(), GlobalMapOfHandlers.end());
		else if(DeletionMethNo == 2)
		{
			radTmhg::iterator IterStartDel = GlobalMapOfHandlers.end(); --IterStartDel;
			for(radTmhg::iterator iter = IterStartDel; iter != GlobalMapOfHandlers.begin(); --iter)
				GlobalMapOfHandlers.erase(iter);
			GlobalMapOfHandlers.erase(GlobalMapOfHandlers.begin());
		}
		GlobalUniqueMapKey = 1;

		MapOfDrawAttr.erase(MapOfDrawAttr.begin(), MapOfDrawAttr.end());

		if(SendingIsRequired) Send.Int(0);
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ComputeGeometricalVolume(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		if(SendingIsRequired) Send.Double(g3dPtr->VolumeWithSym());
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ComputeGeometricalLimits(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		double LimArr[6];
		g3dPtr->Limits(0, LimArr);

		if(SendingIsRequired) Send.DoubleList(LimArr, 6);
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTApplication::ReturnInput(double Input, int NumTimes)
{
	double* OutArray = NULL;

	try
	{
		OutArray = new double[NumTimes];

		for(int i=0; i<NumTimes; ++i) 
		{
			OutArray[i] = Input;
		}
		Send.DoubleList(OutArray, NumTimes);
		if(OutArray != NULL) { delete[] OutArray; OutArray = NULL;}
	}
	catch(...)
	{
		if(OutArray != NULL) { delete[] OutArray; OutArray = NULL;}
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetMemAllocMethForIntrctMatr(char* TotOrParts)
{
	short InMemAllocForIntrctMatrTotAtOnce = 0;
	if((!strcmp(TotOrParts, "tot")) || (!strcmp(TotOrParts, "Tot")) || (!strcmp(TotOrParts, "TOT"))) InMemAllocForIntrctMatrTotAtOnce = 1;
	else if((!strcmp(TotOrParts, "parts")) || (!strcmp(TotOrParts, "Parts")) || (!strcmp(TotOrParts, "PARTS"))) InMemAllocForIntrctMatrTotAtOnce = 0;
	else { Send.ErrorMessage("Radia::Error046"); return 0;}

	MemAllocForIntrctMatrTotAtOnce = InMemAllocForIntrctMatrTotAtOnce;

	if(SendingIsRequired) Send.Int(1);
	return 1;
}

//-------------------------------------------------------------------------
