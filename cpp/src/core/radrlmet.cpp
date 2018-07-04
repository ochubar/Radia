/*-------------------------------------------------------------------------
*
* File name:      radrlmet.cpp
*
* Project:        RADIA
*
* Description:    Relaxation methods
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radrlmet.h"
#include "radyield.h"

#include <time.h>

//-------------------------------------------------------------------------

extern radTYield radYield;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTIterativeRelaxMeth::ComputeRelaxStatusParam(const TVector3d* NewMagnArray, const TVector3d* OldMagnArray, const TVector3d* NewFieldArray)
{
	double BufMisfitM, BufMaxModM, BufMaxModH, TestBufMaxModM, TestBufMaxModH;
	BufMisfitM=0.;
	BufMaxModM=BufMaxModH=TestBufMaxModM=TestBufMaxModH=1.E-17;
	TVector3d Mnew_mi_MoldVect;

	radTRelaxStatusParam& RelStatParR = IntrctPtr->RelaxStatusParam;

	for(int i=0; i<IntrctPtr->AmOfMainElem; i++)
	{
		if(RelStatParR.MisfitM >= 0.)
		{
			Mnew_mi_MoldVect = NewMagnArray[i] - OldMagnArray[i];
			BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y 
						+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;
		}
		if(RelStatParR.MaxModM >= 0.)
		{
			TestBufMaxModM = sqrt(NewMagnArray[i].x*NewMagnArray[i].x 
								+ NewMagnArray[i].y*NewMagnArray[i].y 
								+ NewMagnArray[i].z*NewMagnArray[i].z);
			if(TestBufMaxModM > BufMaxModM) BufMaxModM = TestBufMaxModM;
		}
		if(RelStatParR.MaxModH >= 0.)
		{
			TestBufMaxModH = sqrt(NewFieldArray[i].x*NewFieldArray[i].x 
								+ NewFieldArray[i].y*NewFieldArray[i].y 
								+ NewFieldArray[i].z*NewFieldArray[i].z);
			if(TestBufMaxModH > BufMaxModH) BufMaxModH = TestBufMaxModH;
		}
	}
	if(RelStatParR.MisfitM >= 0.) RelStatParR.MisfitM = sqrt(BufMisfitM/IntrctPtr->AmOfMainElem);
	if(RelStatParR.MaxModM >= 0.) RelStatParR.MaxModM = BufMaxModM;
	if(RelStatParR.MaxModH >= 0.) RelStatParR.MaxModH = BufMaxModH;
}

//-------------------------------------------------------------------------

void radTIterativeRelaxMeth::MakeN_iter(int IterNum)
{
	for(int i=0; i<(IterNum-1); i++)
	{
		DefineNewMagnetizations(); 

		if(radYield.Check()==0) return; // To allow multitasking on Mac: consider better places for this
	}

	//radTSend Send;
	TVector3d* OldMagnArray = NULL;
	//try
	//{
		OldMagnArray = new TVector3d[IntrctPtr->AmOfMainElem];
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	Send.ErrorMessage(radExceptionPtr->what());	return;
	//}
	//catch (...)
	//{
	//	Send.ErrorMessage("Radia::Error999");	return;
	//}

	for(int k=0; k<IntrctPtr->AmOfMainElem; k++) OldMagnArray[k] = (IntrctPtr->g3dRelaxPtrVect[k])->Magn;
	DefineNewMagnetizations();
	for(int q=0; q<IntrctPtr->AmOfMainElem; q++) 
		IntrctPtr->NewMagnArray[q] = (IntrctPtr->g3dRelaxPtrVect[q])->Magn;

	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, OldMagnArray, IntrctPtr->NewFieldArray);
	delete[] OldMagnArray;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTSimpleRelaxation::DefineNewMagnetizations()
{
	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;
	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d H_atElemStrNo(0.,0.,0.);
		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
			H_atElemStrNo += (IntrctPtr->InteractMatrix[StrNo][ColNo])*(IntrctPtr->NewMagnArray[ColNo]);

		IntrctPtr->NewFieldArray[StrNo] = H_atElemStrNo + IntrctPtr->ExternFieldArray[StrNo];
	}

	double One_mi_RelaxParam = 1.- RelaxParam;
	radTg3dRelax* g3dRelaxPtr = NULL;
	for(int StNo=0; StNo<LocAmOfMainElem; StNo++)
	{
		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StNo];

		IntrctPtr->NewMagnArray[StNo] = One_mi_RelaxParam*g3dRelaxPtr->Magn;
		g3dRelaxPtr->Magn = ((radTMaterial*)(g3dRelaxPtr->MaterHandle.rep))->M(IntrctPtr->NewFieldArray[StNo]);
		IntrctPtr->NewMagnArray[StNo] += RelaxParam*g3dRelaxPtr->Magn;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRelaxationMethNo_2::DefineNewMagnetizations()
{
	TMatrix3d InstantKsiTensor;
	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;
	int AmOfMainElem_mi_One = LocAmOfMainElem - 1;

	TVector3d* OldField = IntrctPtr->NewMagnArray;
	TVector3d* NewField = IntrctPtr->NewFieldArray;
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d H_atElemStrNo(0.,0.,0.);
		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
			H_atElemStrNo += (IntrcMat[StrNo][ColNo])*((IntrctPtr->g3dRelaxPtrVect[ColNo])->Magn);

		OldField[StrNo] = NewField[StrNo];
		NewField[StrNo] = H_atElemStrNo + IntrctPtr->ExternFieldArray[StrNo];
	}

	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.), MagnFromMaterRel, InstantMr; // The later is not actually used here
	TMatrix3d E(E_Str0, E_Str1, E_Str2), mi_Eta, E_pl_Eta, InvE_pl_Eta;
	radTg3dRelax* g3dRelaxPtr = NULL;
	double One_mi_RelaxParam = 1.- RelaxParam;
	for(int StNo=0; StNo<LocAmOfMainElem; StNo++)
	{
		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StNo];

		((radTMaterial*)(g3dRelaxPtr->MaterHandle.rep))->
			DefineInstantKsiTensor(OldField[StNo], InstantKsiTensor, InstantMr);
		mi_Eta = InstantKsiTensor*IntrcMat[StNo][StNo];

		E_pl_Eta = E - mi_Eta;
		Matrix3d_inv(E_pl_Eta, InvE_pl_Eta);

		MagnFromMaterRel = ((radTMaterial*)(g3dRelaxPtr->MaterHandle.rep))->M(NewField[StNo]);
		g3dRelaxPtr->Magn = RelaxParam*(InvE_pl_Eta*(MagnFromMaterRel - (mi_Eta*g3dRelaxPtr->Magn)))
							+ One_mi_RelaxParam*g3dRelaxPtr->Magn;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRelaxationMethNo_3::DefineNewMagnetizations()
{
	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
	TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr;
	TMatrix3d MultByInstKsi;
	TVector3d MultByInstMr, Mnew_mi_MoldVect;
	double BufMisfitM=0.;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = NULL;
	radTMaterial* MaterPtr = NULL;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;
	int AmOfMainElem_mi_One = LocAmOfMainElem - 1; // What for?
	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
			if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += IntrcMat[StrNo][ColNo] * MagnAr[ColNo];
		QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		MaterPtr->MultMatrByInstKsiAndMr(NewFieldAr[StrNo], IntrcMat[StrNo][StrNo], MultByInstKsi, MultByInstMr);

		BufMatr = E - MultByInstKsi;
		Matrix3d_inv(BufMatr, InvBufMatr);
		NewFieldAr[StrNo] = InvBufMatr * (MultByInstMr + QuasiExtFieldAtElemStrNo);

		MagnAr[StrNo] = MaterPtr->M(NewFieldAr[StrNo]);

		Mnew_mi_MoldVect = MagnAr[StrNo] - g3dRelaxPtr->Magn;
		BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y 
					+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;

		g3dRelaxPtr->Magn = MagnAr[StrNo];
	}
	InstMisfitM = sqrt(BufMisfitM/LocAmOfMainElem);
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_3::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM(); // Consider removing
	}

	int IterCount = 0;
	while(InstMisfitM > PrecOnMagnetiz)
	{
		if(++IterCount > MaxIterNumber) break;
		DefineNewMagnetizations();

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
	}

	IntrctPtr->RelaxStatusParam.MisfitM = -1.;
	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, NULL, IntrctPtr->NewFieldArray);
	IntrctPtr->RelaxStatusParam.MisfitM = InstMisfitM;

	return IterCount-1;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTRelaxationMethNo_a5::radTRelaxationMethNo_a5(radTInteraction* InInteractionPtr) : radTIterativeRelaxMeth(InInteractionPtr) 
{ 
	IntrctPtr = InInteractionPtr; InstMisfitM = 1.E+23;

	radTRelaxSubInterval* TmpSubIntervArray;
	//try
	//{
		TmpSubIntervArray = new radTRelaxSubInterval[IntrctPtr->AmOfRelaxSubInterv];

		int RelaxTogetherCount = 0;
		int MaxSize = 0;

		for(int i=0; i<IntrctPtr->AmOfRelaxSubInterv; i++)
		{
			radTRelaxSubInterval& LocSubInterv = IntrctPtr->RelaxSubIntervArray[i];
			if(LocSubInterv.SubIntervalID == RelaxTogether) 
			{
				TmpSubIntervArray[RelaxTogetherCount] = LocSubInterv;
				RelaxTogetherCount++;

				int CurSize = LocSubInterv.FinNo - LocSubInterv.StartNo + 1;
				if(CurSize>MaxSize) MaxSize = CurSize;
			}
		}
		AmOfRelaxTogether = RelaxTogetherCount;

		SizeOfAuxs = 3*MaxSize;

		MathMethPtr = new radTMathLinAlgEq(SizeOfAuxs);

		if(AmOfRelaxTogether != 0)
		{
			AuxMatr1 = new double*[SizeOfAuxs];
			AuxMatr2 = new double*[SizeOfAuxs];
			
			AuxArray = new double[SizeOfAuxs];

			for(int m=0; m<SizeOfAuxs; m++)
			{
				AuxMatr1[m] = new double[SizeOfAuxs];
				AuxMatr2[m] = new double[SizeOfAuxs];
			}
		}
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	Send.ErrorMessage(radExceptionPtr->what());	return;
	//}
	//catch (...)
	//{
	//	Send.ErrorMessage("Radia::Error999");	return;
	//}

	delete[] TmpSubIntervArray; 
}

//-------------------------------------------------------------------------

radTRelaxationMethNo_a5::~radTRelaxationMethNo_a5()
{
	for(int j=0; j<SizeOfAuxs; j++)
	{
		delete[] (AuxMatr1[j]);
		delete[] (AuxMatr2[j]);
	}
	delete[] AuxMatr1;
	delete[] AuxMatr2;
	delete[] AuxArray;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_a5::DefineNewMagnetizations()
{

	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.), ZeroVect(0.,0.,0.);
	TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr, ZeroMatr(ZeroVect, ZeroVect, ZeroVect);
	TMatrix3d MultByInstKsi;
	TVector3d MultByInstMr, Mnew_mi_MoldVect;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = NULL;
	radTMaterial* MaterPtr = NULL;

	double BufMisfitM=0.;

	int StrNo = 0;
	int RelaxTogetherCount = -1;

	for(int IntrvNo=0; IntrvNo<IntrctPtr->AmOfRelaxSubInterv; IntrvNo++)
	{
		radTRelaxSubInterval& CurrentSubInterv = IntrctPtr->RelaxSubIntervArray[IntrvNo];

		if(CurrentSubInterv.SubIntervalID == RelaxTogether)
		{
			RelaxTogetherCount++;
			for(StrNo = CurrentSubInterv.StartNo; StrNo <= CurrentSubInterv.FinNo; StrNo++)
			{
				TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
				int ColNo=0;
				for(ColNo = 0; ColNo < CurrentSubInterv.StartNo; ColNo++)
					QuasiExtFieldAtElemStrNo += IntrcMat[StrNo][ColNo] * MagnAr[ColNo];
				for(ColNo = CurrentSubInterv.FinNo+1; ColNo < LocAmOfMainElem; ColNo++)
					QuasiExtFieldAtElemStrNo += IntrcMat[StrNo][ColNo] * MagnAr[ColNo];
				QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

				int AuxMatrStrNo = 3*(StrNo - CurrentSubInterv.StartNo);
				int AuxMatrStrNo_p1=AuxMatrStrNo+1, AuxMatrStrNo_p2=AuxMatrStrNo+2;

				double* AuxMatr1StrNoPtr = AuxMatr1[AuxMatrStrNo];
				double* AuxMatr1StrNo_p1Ptr = AuxMatr1[AuxMatrStrNo_p1];
				double* AuxMatr1StrNo_p2Ptr = AuxMatr1[AuxMatrStrNo_p2];

				g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
				MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

				TVector3d ContribFromMr(0.,0.,0.);
				for(ColNo = CurrentSubInterv.StartNo; ColNo <= CurrentSubInterv.FinNo; ColNo++)
				{
					MaterPtr->MultMatrByInstKsiAndMr(NewFieldAr[ColNo], IntrcMat[StrNo][ColNo], MultByInstKsi, MultByInstMr);

					if(ColNo == StrNo) BufMatr = E - MultByInstKsi;
					else BufMatr = ZeroMatr - MultByInstKsi;
					
					ContribFromMr += MultByInstMr;

					int AuxMatrColNo = 3*(ColNo - CurrentSubInterv.StartNo);
					int AuxMatrColNo_p1=AuxMatrColNo+1, AuxMatrColNo_p2=AuxMatrColNo+2;

					TVector3d& BufMatrStr0 = BufMatr.Str0;
					AuxMatr1StrNoPtr[AuxMatrColNo] = BufMatrStr0.x; AuxMatr1StrNoPtr[AuxMatrColNo_p1] = BufMatrStr0.y; AuxMatr1StrNoPtr[AuxMatrColNo_p2] = BufMatrStr0.z;
					TVector3d& BufMatrStr1 = BufMatr.Str1;
					AuxMatr1StrNo_p1Ptr[AuxMatrColNo] = BufMatrStr1.x; AuxMatr1StrNo_p1Ptr[AuxMatrColNo_p1] = BufMatrStr1.y; AuxMatr1StrNo_p1Ptr[AuxMatrColNo_p2] = BufMatrStr1.z;
					TVector3d& BufMatrStr2 = BufMatr.Str2;
					AuxMatr1StrNo_p2Ptr[AuxMatrColNo] = BufMatrStr2.x; AuxMatr1StrNo_p2Ptr[AuxMatrColNo_p1] = BufMatrStr2.y; AuxMatr1StrNo_p2Ptr[AuxMatrColNo_p2] = BufMatrStr2.z;
				}

				QuasiExtFieldAtElemStrNo += ContribFromMr;
				AuxArray[AuxMatrStrNo] = QuasiExtFieldAtElemStrNo.x;
				AuxArray[AuxMatrStrNo_p1] = QuasiExtFieldAtElemStrNo.y;
				AuxArray[AuxMatrStrNo_p2] = QuasiExtFieldAtElemStrNo.z;
			}
			int SizeMatr = 3*(CurrentSubInterv.FinNo - CurrentSubInterv.StartNo + 1);
			MathMethPtr->InverseMatrix(AuxMatr1, SizeMatr, AuxMatr2);

			int TriplCount = -1;
			int NewFieldArIndx = CurrentSubInterv.StartNo;

			for(int i=0; i<SizeMatr; i++)
			{
				TriplCount++;

				double Sum=0.;
				double* AuxMatr2_i = AuxMatr2[i];
				for(int j=0; j<SizeMatr; j++) Sum += AuxMatr2_i[j] * AuxArray[j];

				if(TriplCount==0) NewFieldAr[NewFieldArIndx].x = Sum;
				else if(TriplCount==1) NewFieldAr[NewFieldArIndx].y = Sum;
				else
				{
					NewFieldAr[NewFieldArIndx].z = Sum;

					NewFieldArIndx++;
					TriplCount = -1;
				}
			}
			for(StrNo = CurrentSubInterv.StartNo; StrNo <= CurrentSubInterv.FinNo; StrNo++)
			{
				g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
				MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

				MagnAr[StrNo] = MaterPtr->M(NewFieldAr[StrNo]);
				Mnew_mi_MoldVect = MagnAr[StrNo] - g3dRelaxPtr->Magn;
				BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y 
							+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;
				g3dRelaxPtr->Magn = MagnAr[StrNo];
			}
		}
		
		if(CurrentSubInterv.SubIntervalID == RelaxApart)
		{
			for(StrNo = CurrentSubInterv.StartNo; StrNo <= CurrentSubInterv.FinNo; StrNo++)
			{
				TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
				for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
					if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += IntrcMat[StrNo][ColNo] * MagnAr[ColNo];
				QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

				g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
				MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
				MaterPtr->MultMatrByInstKsiAndMr(NewFieldAr[StrNo], IntrcMat[StrNo][StrNo], MultByInstKsi, MultByInstMr);

				BufMatr = E - MultByInstKsi;
				Matrix3d_inv(BufMatr, InvBufMatr);
				NewFieldAr[StrNo] = InvBufMatr * (MultByInstMr + QuasiExtFieldAtElemStrNo);

				MagnAr[StrNo] = MaterPtr->M(NewFieldAr[StrNo]);

				Mnew_mi_MoldVect = MagnAr[StrNo] - g3dRelaxPtr->Magn;
				BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y 
							+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;
				g3dRelaxPtr->Magn = MagnAr[StrNo];
			}
		}
	}
	InstMisfitM = sqrt(BufMisfitM/LocAmOfMainElem);
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_a5::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{ // Absolutely the same as for MethNo_3

	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();  // Consider removing
	}

	int IterCount = 0;
	while(InstMisfitM > PrecOnMagnetiz)
	{
		if(++IterCount > MaxIterNumber) break;
		DefineNewMagnetizations();

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
	}

	IntrctPtr->RelaxStatusParam.MisfitM = -1.;
	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, NULL, IntrctPtr->NewFieldArray);
	IntrctPtr->RelaxStatusParam.MisfitM = InstMisfitM;

	return IterCount-1;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// Good operational version; failed to converge in case of fine mesh of Soleil booster quads
// restored 150505, because another version appears to fail to converge even in the Radia dipole magnet example
void radTRelaxationMethNo_4::DefineNewMagnetizations()
{
	//TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
	//TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr;
	//TMatrix3d MultByInstKsi;
	//TVector3d MultByInstMr, Mnew_mi_MoldVect;
	TVector3d Mnew_mi_MoldVect;
	double BufMisfitM=0.;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504

	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = NULL;
	radTMaterial* MaterPtr = NULL;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	double NormFact = 1./double(LocAmOfMainElem);
	double BestPrecMagnE2 =  DesiredPrecOnMagnetizE2*NormFact;
	double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.25*NormFact*InstMisfitMe2;
	//double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.1*NormFact*InstMisfitMe2;

	if(LocPrecMagnE2 < BestPrecMagnE2) LocPrecMagnE2 = BestPrecMagnE2;
	radTRelaxAuxData *tRelaxAuxData = mpRelaxAuxData; //OC06112003
	const int MaxConseqBadPasses = 1; //OC06112003

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
		TMatrix3df* MatrArrayPtr = IntrcMat[StrNo]; //OC250504
		//TMatrix3d* MatrArrayPtr = IntrcMat[StrNo]; //OC250504

		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
		{
			if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += MatrArrayPtr[ColNo] * MagnAr[ColNo];
		}
		QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d& InstantH = NewFieldAr[StrNo];
		TVector3d& InstantM = MagnAr[StrNo];

		TVector3d PrevH = InstantH;

		MaterPtr->FindNewH(InstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2);
		//MaterPtr->FindNewH(InstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2, g3dRelaxPtr);
		//MaterPtr->FindNewH(InstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2, g3dRelaxPtr, gpRelaxAuxData + StrNo); //OC140103

		//InstantH = (tRelaxAuxData->mRelaxPar)*InstantH + (1. - tRelaxAuxData->mRelaxPar)*PrevH; //OC131103
		//OC: commented out 150304

		//InstantM = MaterPtr->M(InstantH);
		
		TVector3d PureNewM = MaterPtr->M(InstantH);
		InstantM = PureNewM;

		Mnew_mi_MoldVect = PureNewM - g3dRelaxPtr->Magn;
		double NewDifMe2 = Mnew_mi_MoldVect.AmpE2(); //OC06112003
        BufMisfitM += NewDifMe2;

		//tRelaxAuxData->Update(NewDifMe2, MaxConseqBadPasses, mRelaxParModFact, mRelaxParMin); //OC06112003
		//OC: commented out 150304

		//InstantM = mRelaxPar*PureNewM + (1. - mRelaxPar)*(g3dRelaxPtr->Magn); //OC041103
		//InstantM = (tRelaxAuxData->mRelaxPar)*PureNewM + (1. - tRelaxAuxData->mRelaxPar)*(g3dRelaxPtr->Magn); //OC041103
		//InstantM = PureNewM;

		//tRelaxAuxData++; //OC041103
		//OC: commented out 150304

		//Mnew_mi_MoldVect = InstantM - g3dRelaxPtr->Magn;
		//BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y + Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;

		//double CurRelaxPar = (gpRelaxAuxData + StrNo)->RelaxPar; //OC140103
		//BufMisfitM += (Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y + Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z)/(CurRelaxPar*CurRelaxPar); //OC140103

		g3dRelaxPtr->Magn = InstantM; 
	}
	double NewInstMisfitMe2 = BufMisfitM/LocAmOfMainElem;

	//if(NewInstMisfitMe2 > mMisfitE2RatToStartModifRelaxPar*InstMisfitMe2) //OC041103
	//{
	//	if(mRelaxPar > mRelaxParMin) mRelaxPar *= mRelaxParModFact;
	//	int Aha = 1;
	//}
	InstMisfitMe2 = NewInstMisfitMe2;
}

//-------------------------------------------------------------------------
// Strange version; failed to converge in case of Radia dipole magnet example
// commented out 150505
/**
void radTRelaxationMethNo_4::DefineNewMagnetizations()
{
	TVector3d Mnew_mi_MoldVect;
	double BufMisfitM=0.;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;

	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = NULL;
	radTMaterial* MaterPtr = NULL;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	double NormFact = 1./double(LocAmOfMainElem);
	double BestPrecMagnE2 =  DesiredPrecOnMagnetizE2*NormFact;
	double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.25*NormFact*InstMisfitMe2;

	if(LocPrecMagnE2 < BestPrecMagnE2) LocPrecMagnE2 = BestPrecMagnE2;
	//radTRelaxAuxData *tRelaxAuxData = mpRelaxAuxData; //OC06112003
	//const int MaxConseqBadPasses = 1; //OC06112003
	const int MaxInitialBadPasses = 10; //OC06112003

	//if((!mKeepPrevOldValues) && ((mNumConvergPasses > 0) || ((mNumConvergPasses == 0) && (mNumDivergPasses == 0)) || ((mIterCount < MaxInitialBadPasses)))) IntrctPtr->StoreAuxOldArrays(); //OC300504
	if(!mKeepPrevOldValues) IntrctPtr->StoreAuxOldArrays(); //OC300504

	double *tElemVolumes = mElemVolumeArray; //OC010604
	const double ConstForM = -1./(400.*(3.14159265358979));
	double SumEnergy = 0;

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TMatrix3df* MatrArrayPtr = IntrcMat[StrNo]; //OC250504
		TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);

		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
		{
			if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += MatrArrayPtr[ColNo] * MagnAr[ColNo];
		}
		QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d& InstantH = NewFieldAr[StrNo];
		TVector3d& InstantM = MagnAr[StrNo];
		TVector3d PrevH = InstantH;


		if(StrNo == 42)
		{
			int aha = 1;

		}

		MaterPtr->FindNewH(InstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2);

		TVector3d PureNewM = MaterPtr->M(InstantH);
		Mnew_mi_MoldVect = PureNewM - g3dRelaxPtr->Magn;
		double NewDifMe2 = Mnew_mi_MoldVect.AmpE2(); //OC06112003
        BufMisfitM += NewDifMe2;

        //InstantH = mRelaxPar*InstantH + (1. - mRelaxPar)*PrevH; //OC260504
		//InstantM = MaterPtr->M(InstantH);

		//InstantM = mRelaxPar*PureNewM + (1. - mRelaxPar)*(g3dRelaxPtr->Magn); //OC041103
		//InstantM = 0.5*(MaterPtr->M(InstantH) + (mRelaxPar*PureNewM + (1. - mRelaxPar)*(g3dRelaxPtr->Magn))); //OC280504

		//Mnew_mi_MoldVect = InstantM - g3dRelaxPtr->Magn;
		//BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y + Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;

		InstantM = PureNewM;
		g3dRelaxPtr->Magn = InstantM; 

		SumEnergy += (InstantM*(InstantH + InstantM))*(*(tElemVolumes++)); //OC010604
	}
	SumEnergy *= ConstForM; //OC010604
	
	//InstMisfitMe2 = BufMisfitM*NormFact;
    double NewInstMisfitMe2 = BufMisfitM*NormFact;

	//if((NewInstMisfitMe2 > InstMisfitMe2) && (mIterCount > MaxInitialBadPasses)) //OC053004
	//{
	//       //if(mNumConvergPasses > 0) IntrctPtr->StoreAuxOldArrays();
	//       if(mRelaxPar > mRelaxParMin) 
	//	{
	//		mRelaxPar *= mRelaxParModFact; //OC300504
	//           mKeepPrevOldValues = true;
	//	}
	//	else mKeepPrevOldValues = false;

	//       mNumConvergPasses = 0;
	//       mNumDivergPasses++;

	////	//if(mNumDivergPasses > mNumDivergPassesMax)
	////	//{
	////       //    if(mRelaxPar > mRelaxParMin) mRelaxPar *= mRelaxParModFact;
	////	//           int Aha = 1;
	////	//		////test
	////	//		//char ErrorMesTitle[] = "Radia Debug";
	////	//		//char ErrorStr[100];
	////	//		//int j = sprintf(ErrorStr, "mRelaxPar: %g ", mRelaxPar);
	////	//		//j += sprintf(ErrorStr + j, "          NewInstMisfitMe2: %g   ", NewInstMisfitMe2);
	////	//		//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	////	//		//int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
	////	//		////end test
	////	//}
	//}
	//else
	//{
	//       mNumDivergPasses = 0;
	//       mNumConvergPasses++;
	//	mKeepPrevOldValues = false;
	//	mSysEnergyMin = SumEnergy;

	//	//mRelaxPar = 1.; //?????

	////	//      if(mNumConvergPasses > mNumConvergPassesMax)
	////	//{
	////	//          if(mRelaxPar < 1.) 
	////	//	{
	////	//		mRelaxPar /= mRelaxParModFact;
	////	//		if(mRelaxPar > 1.) mRelaxPar = 1.;
	////	//		mNumConvergPasses = 0;
	////	//	}
	////	//          //int Aha = 1;
	////	//}
	////	//      InstMisfitMe2 = NewInstMisfitMe2;
	//}

	if((NewInstMisfitMe2 > InstMisfitMe2) && (mIterCount > MaxInitialBadPasses))
	{
        mNumConvergPasses = 0;
        mNumDivergPasses++;
		if(mNumDivergPasses >= MaxInitialBadPasses) mBadConverg = true;
	}
	else
	{
        mNumDivergPasses = 0;
        mNumConvergPasses++;
	}

	if(mBadConverg)
	{
		if(SumEnergy < mSysEnergyMin)
		{
			StoreOptimValuesFromOldArrays();
			mRelaxPar = 1.; //?????
			mKeepPrevOldValues = false;
			mSysEnergyMin = SumEnergy;
		}
		else
		{
			if(SumEnergy < mSysEnergy)
			{
				mKeepPrevOldValues = false;
			}
			else
			{
				//if(mRelaxPar > mRelaxParMin)
				//{
				mRelaxPar *= mRelaxParModFact;
				mKeepPrevOldValues = true;
				RestoreOptimValuesToOldArrays();
				//}
				//else
				//{
				//	mKeepPrevOldValues = false;
				//}
			}
		}
		if(mRelaxPar < 1.) CorrectMagnAndFieldArraysWithRelaxPar();
	}

	mSysEnergy = SumEnergy;
    InstMisfitMe2 = NewInstMisfitMe2;
}
**/
//-------------------------------------------------------------------------

int radTRelaxationMethNo_4::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	DesiredPrecOnMagnetizE2 = PrecOnMagnetiz * PrecOnMagnetiz;

	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	//SetupAuxArrays(); //OC140103
	//SetupElemVolumeArray(); //OC150505 //OC010604
	//SetupOptimValuesArrays(); //OC150505 //OC020604

	mMethNo = 1;
	mIterCount = 0;
	mSysEnergy = 0;
	mSysEnergyMin = 1e+23;
	mKeepPrevOldValues = false;
	mBadConverg = false;

	double MinInstMisfitMe2 = 1.e+30;
	while(InstMisfitMe2 > DesiredPrecOnMagnetizE2)
	{
		if(++mIterCount > MaxIterNumber) break;
		DefineNewMagnetizations();

		//if(MinInstMisfitMe2 > InstMisfitMe2) 
		//{
		//	MinInstMisfitMe2 = InstMisfitMe2;
		//	IntrctPtr->StoreAuxOldArrays();
		//}

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this

		//test
		//mRelaxPar = 1./pow((double)IterCount + 1., 0.35);
		//mRelaxPar *= pow(((double)mIterCount)/((double)mIterCount + 1.), 0.4);
		//end test
	}

	//if(mIterCount > MaxIterNumber)
	//{
    //	IntrctPtr->RestoreAuxOldArrays();
	//	InstMisfitMe2 = MinInstMisfitMe2;
	//}

	IntrctPtr->RelaxStatusParam.MisfitM = -1.;
	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, NULL, IntrctPtr->NewFieldArray);
	
	IntrctPtr->RelaxStatusParam.MisfitM = sqrt(InstMisfitMe2);

	//DeleteAuxArrays(); //OC140103
	//DeleteElemVolumeArray(); //OC150505 //OC010604
	//DeleteOptimValuesArrays(); //OC150505 //OC020604

	return mIterCount-1;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_4::CorrectMagnAndFieldArraysWithRelaxPar() //OC300504
{
	if((mRelaxPar < 0.) || (mRelaxPar > 1.)) return;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	TVector3d *tMagnAr = IntrctPtr->NewMagnArray;
	TVector3d *tFieldAr = IntrctPtr->NewFieldArray;
	TVector3d *tOldMagnAr = IntrctPtr->AuxOldMagnArray;
	TVector3d *tOldFieldAr = IntrctPtr->AuxOldFieldArray;

	radTg3dRelax* g3dRelaxPtr = NULL;
	radTMaterial* MaterPtr = NULL;

    for(int k=0; k<LocAmOfMainElem; k++)
	{
		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[k];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

        *tFieldAr = mRelaxPar*(*tFieldAr) + (1. - mRelaxPar)*(*tOldFieldAr);
		*tMagnAr = MaterPtr->M(*tFieldAr);
        //*tMagnAr = mRelaxPar*(*tMagnAr) + (1. - mRelaxPar)*(*tOldMagnAr);

		g3dRelaxPtr->Magn = *tMagnAr;
		tMagnAr++; tFieldAr++; tOldMagnAr++; tOldFieldAr++;
	}
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_4::DefineNewMagnetizationsTest()
{
	double MaxNumDivergPassesToSwitchMeth = 5;

	TVector3d Mnew_mi_MoldVect;
	double BufMisfitM=0.;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	//TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = NULL;
	radTMaterial* MaterPtr = NULL;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	double NormFact = 1./double(LocAmOfMainElem);
	//double BestPrecMagnE2 =  DesiredPrecOnMagnetizE2*NormFact;
	//double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.25*NormFact*InstMisfitMe2;
	//if(LocPrecMagnE2 < BestPrecMagnE2) LocPrecMagnE2 = BestPrecMagnE2;

	if((mMethNo == 1) && (mNumDivergPasses >= MaxNumDivergPassesToSwitchMeth)) mMethNo = 2;

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d& InstantH = NewFieldAr[StrNo];
		TVector3d PrevH = InstantH;

        InstantH = FindNewFieldTwoSteps(StrNo, mMethNo);

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d& InstantM = MagnAr[StrNo];

		InstantM = MaterPtr->M(InstantH);
		Mnew_mi_MoldVect = InstantM - g3dRelaxPtr->Magn;
		double NewDifMe2 = Mnew_mi_MoldVect.AmpE2();
        BufMisfitM += NewDifMe2;

		g3dRelaxPtr->Magn = InstantM; 
	}
	double NewInstMisfitMe2 = BufMisfitM*NormFact;

	if(NewInstMisfitMe2 > InstMisfitMe2)
	{
		mNumConvergPasses = 0;
		mNumDivergPasses++;
		//if(mNumDivergPasses > mNumDivergPassesMax)
		//{
        //          if(mRelaxPar > mRelaxParMin) mRelaxPar *= mRelaxParModFact;
        //          int Aha = 1;
		//	////test
		//	//char ErrorMesTitle[] = "Radia Debug";
		//	//char ErrorStr[100];
		//	//int j = sprintf(ErrorStr, "mRelaxPar: %g ", mRelaxPar);
		//	//j += sprintf(ErrorStr + j, "          NewInstMisfitMe2: %g   ", NewInstMisfitMe2);
		//	//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
		//	//int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
		//	////end test
		//}
	}
	else
	{
		mNumDivergPasses = 0;
		mNumConvergPasses++;
		//if(mNumConvergPasses > mNumConvergPassesMax)
		//{
        //          if(mRelaxPar < 1.) 
		//	{
		//		mRelaxPar /= mRelaxParModFact;
		//		if(mRelaxPar > 1.) mRelaxPar = 1.;
		//	}
        //          int Aha = 1;
		//}
	}

	InstMisfitMe2 = BufMisfitM*NormFact;
}

//-------------------------------------------------------------------------

TVector3d radTRelaxationMethNo_4::FindNewFieldTwoSteps(int i, int MethNo)
{
	const double DesiredPrecOnMagn = 1.e-08;
	//const double DesiredPrecOnMagnE2 = DesiredPrecOnMagn*DesiredPrecOnMagn;
	int MaxLoopPasses = 10;
	double AbsMisfitM = 5;

    TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
    TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr, LocBufMatr, InvLocBufMatr, ExtBufMatr, InvExtBufMatr;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
    TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

    radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
	radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

	int N = IntrctPtr->AmOfMainElem;
	TVector3d& Hi = NewFieldAr[i];
	TVector3d& Mi = MagnAr[i];
    TMatrix3d Qii = IntrcMat[i][i];

	TMatrix3d SumMatr1;
	TVector3d SumVect1 = ExternFieldAr[i];
	for(int j=0; j<N; j++)
	{
		if(j == i) continue;

		TMatrix3d Qij = IntrcMat[i][j];

		if(MethNo == 1)
		{
            SumVect1 += Qij*MagnAr[j];
		}
		else if(MethNo == 2)
		{
			TMatrix3d Qjj = IntrcMat[j][j];
			TMatrix3d Qji = IntrcMat[j][i];

			TVector3d &Hj = NewFieldAr[j], Mr_j;
			TMatrix3d Ksi_j;
			radTg3dRelax* Loc_g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[j];
			radTMaterial* Loc_MaterPtr = (radTMaterial*)(Loc_g3dRelaxPtr->MaterHandle.rep);
			Loc_MaterPtr->DefineInstantKsiTensor(Hj, Ksi_j, Mr_j);

			//LocBufMatr = E - Qjj*Ksi_j;
			//Matrix3d_inv(LocBufMatr, InvLocBufMatr);
			//TMatrix3d Ksi_j_InvLocBufMatr = Ksi_j*InvLocBufMatr;
			//TMatrix3d Qij_Ksi_j_InvLocBufMatr = Qij*Ksi_j_InvLocBufMatr;
			//SumMatr1 += Qij_Ksi_j_InvLocBufMatr*Qji;

			SumMatr1 += (Qij*Ksi_j)*Qji;

			TVector3d LocSumVect(0,0,0);
			for(int k=0; k<N; k++)
			{
				if((k == j) || (k == i)) continue;
				LocSumVect += (IntrcMat[j][k])*MagnAr[k];
			}
			TVector3d Qjj_Mr_j_p_Qji_Mr_i_p_He_j_p_LocSumVect = (Qjj*Mr_j) + ExternFieldAr[j] + LocSumVect;

			//TVector3d Qjj_Mr_j_p_Qji_Mr_i_p_He_j_p_LocSumVect = (Qjj*Mr_j) + (Qji*Mr_i) + ExternFieldAr[j] + LocSumVect;
			//SumVect1 += Qij*((Ksi_j_InvLocBufMatr*Qjj_Mr_j_p_Qji_Mr_i_p_LocSumVect) + Mr_j);

			SumVect1 += Qij*((Ksi_j*Qjj_Mr_j_p_Qji_Mr_i_p_He_j_p_LocSumVect) + Mr_j);
		}
	}

	for(int p=0; p<MaxLoopPasses; p++)
	{
		TMatrix3d Ksi_i;
		TVector3d Mr_i;
		MaterPtr->DefineInstantKsiTensor(Hi, Ksi_i, Mr_i);
        TVector3d PrevHi = Hi;

		TVector3d AuxVect1 = SumVect1 + (Qii*Mr_i);

		if(MethNo == 1)
		{
            ExtBufMatr = E - (Qii*Ksi_i);
		}
		else if(MethNo == 2)
		{
            AuxVect1 += (SumMatr1*Mr_i);
            ExtBufMatr = E - ((Qii + SumMatr1)*Ksi_i);
		}

        Matrix3d_inv(ExtBufMatr, InvExtBufMatr);
		Hi = InvExtBufMatr*AuxVect1;

		Mi = MaterPtr->M(Hi);
		TVector3d MiLin = (Ksi_i*Hi) + Mr_i;
		TVector3d Mi_mi_MiLin = Mi - MiLin;
		double AbsNewMisfitM = sqrt(Mi_mi_MiLin.AmpE2());

		if(AbsNewMisfitM <= DesiredPrecOnMagn) break;

		double Alpha = AbsMisfitM/(AbsMisfitM + AbsNewMisfitM);
		//AbsInstantH = Alpha*NewAbsInstantH + (1 - Alpha)*AbsInstantH;

		AbsMisfitM = AbsNewMisfitM;

		Hi = Alpha*Hi + (1 - Alpha)*PrevHi;
	}

	return Hi;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_4::LpTau(int i, double& q)
{/**
	double a = i;
	int m = 1 + int(log(a)/0.693147);

	double s = 0.;
	for(int k = 1; k <= m; k++)
	{
		int ns = 0;
		for(int l = k; l <= m; l++)
		{
			ns += int(2*D(a/pow(2,l)))*int(2*D(1./pow(2,l+1-k)));
		}
		s += D(0.5*ns)/pow(2,k-1);
	}
	q = s;
**/
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTRelaxationMethNo_8::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsMotNeeded)
{
	if(IntrctPtr == 0) return 0;
	mDesiredPrecOnMagnetizE2 = PrecOnMagnetiz * PrecOnMagnetiz;

	if(!MagnResetIsMotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	double MinInstMisfitMe2 = 1.e+30;
	int ItCnt=0;
	for(ItCnt=0; ItCnt<MaxIterNumber; ItCnt++)
	{
		if(ItCnt > MaxIterNumber) break;
		DefineNewMagnetizations();

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
	}

/**
	//SetupAuxArrays(); //OC140103
	SetupElemVolumeArray(); //OC010604
	SetupOptimValuesArrays(); //OC020604

	mMethNo = 1;
	mSysEnergy = 0;
	mSysEnergyMin = 1e+23;
	mKeepPrevOldValues = false;
	mBadConverg = false;

	while(InstMisfitMe2 > DesiredPrecOnMagnetizE2)
	{
		if(++mIterCount > MaxIterNumber) break;
		DefineNewMagnetizations();

		//if(MinInstMisfitMe2 > InstMisfitMe2) 
		//{
		//	MinInstMisfitMe2 = InstMisfitMe2;
		//	IntrctPtr->StoreAuxOldArrays();
		//}

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this

		//test
		//mRelaxPar = 1./pow((double)IterCount + 1., 0.35);
		//mRelaxPar *= pow(((double)mIterCount)/((double)mIterCount + 1.), 0.4);
		//end test
	}

	//if(mIterCount > MaxIterNumber)
	//{
    //	IntrctPtr->RestoreAuxOldArrays();
	//	InstMisfitMe2 = MinInstMisfitMe2;
	//}

	IntrctPtr->RelaxStatusParam.MisfitM = -1.;
	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, NULL, IntrctPtr->NewFieldArray);
	
	IntrctPtr->RelaxStatusParam.MisfitM = sqrt(InstMisfitMe2);

	//DeleteAuxArrays(); //OC140103
	DeleteElemVolumeArray(); //OC010604
	DeleteOptimValuesArrays(); //OC020604

**/

	return ItCnt;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_8::DefineNewMagnetizations()
{
	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.), MatrElemByInstMr, Mnew_mi_MoldVect;
	TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr, MatrElemByInstKsi;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = NULL;
	radTMaterial* MaterPtr = NULL;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;
	double NormFact = 1./double(LocAmOfMainElem);
    double BufMisfitM=0.;

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
		TMatrix3df* MatrArrayPtr = IntrcMat[StrNo];

		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
		{
			if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += MatrArrayPtr[ColNo] * MagnAr[ColNo];
		}
		QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d& InstantH = NewFieldAr[StrNo];
		TVector3d& InstantM = MagnAr[StrNo];

		MaterPtr->MultMatrByInstKsiAndMr(InstantH, MatrArrayPtr[StrNo], MatrElemByInstKsi, MatrElemByInstMr);

        BufMatr = E - MatrElemByInstKsi;
        Matrix3d_inv(BufMatr, InvBufMatr);
        InstantH = InvBufMatr*(QuasiExtFieldAtElemStrNo + MatrElemByInstMr);

        InstantM = MaterPtr->M(InstantH);

	//BufMatr = E - Matr*InstantKsiTensor;
	//Matrix3d_inv(BufMatr, InvBufMatr);
	//H = InvBufMatr*(H_Ext + Matr*RemMagn);

		Mnew_mi_MoldVect = InstantM - g3dRelaxPtr->Magn;
		double NewDifMe2 = Mnew_mi_MoldVect.AmpE2();
        BufMisfitM += NewDifMe2;

		g3dRelaxPtr->Magn = InstantM; 
	}
	mInstMisfitMe2 = BufMisfitM*NormFact;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRelaxationMethNo_6::SetupInteractionMatrices(const radThg& hg, const radTCompCriterium& CompCrit)
{
	mAmOfParts = 0;

	radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
	if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); throw 0;}
	radTGroup* GroupPtr = Cast.GroupCast(g3dPtr); 
	if(GroupPtr==0) { Send.ErrorMessage("Radia::Error091"); throw 0;}
	mhGroup = hg;

	radThg hEmpty;
	mAmOfParts = (int)(GroupPtr->GroupMapOfHandlers.size());
	IntrctPtr = new radTInteraction[mAmOfParts];
	radTInteraction *tIntrct = IntrctPtr;

	int LocMapCount = 0;
	for(radTmhg::const_iterator iter = GroupPtr->GroupMapOfHandlers.begin(); iter != GroupPtr->GroupMapOfHandlers.end(); ++iter)
	{
		tIntrct->Setup((*iter).second, hEmpty, CompCrit, 0, 1, 1);
		if(tIntrct->SomethingIsWrong)
		{
			delete[] IntrctPtr; IntrctPtr = NULL;
            Send.ErrorMessage("Radia::Error116"); throw 0;
		}

		radThg hgGroup(GroupPtr->CreateGroupIncludingAllMembersExceptIt(iter));
		mMapOfPartHandlers[LocMapCount++] = hgGroup;
		tIntrct++;
	}
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_6::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, double* RelaxStatusParamArray)
{
	int TotAmOfRelaxObj = 0;
	radTInteraction *tIntrct = IntrctPtr;
	for(radTmhg::const_iterator iter = mMapOfPartHandlers.begin(); iter != mMapOfPartHandlers.end(); ++iter)
	{
        tIntrct->ResetM();
        tIntrct->ResetAuxParam();
		tIntrct->AddMoreExternField((*iter).second);
		TotAmOfRelaxObj += tIntrct->OutAmOfRelaxObjs();
		tIntrct++;
	}

    radTGroup* GroupPtr = Cast.GroupCast(Cast.g3dCast(mhGroup.rep)); 
	int ActualIterNum = 0, ActualOuterIterNum = 0;
	double DesiredPrecOnMagnetizE2 = PrecOnMagnetiz*PrecOnMagnetiz;
	double TotMisfitMe2 = 1.E+23;

	for(int i=0; i<MaxIterNumber; i++)
	{
		tIntrct = IntrctPtr;
		int PartCount = 0;
		double BufTotMisfitMe2 = 0;

		for(radTmhg::const_iterator it = GroupPtr->GroupMapOfHandlers.begin(); it != GroupPtr->GroupMapOfHandlers.end(); ++it)
		{
			if(tIntrct->OutAmOfRelaxObjs() > 0)
			{
				//tIntrct->StoreAuxOldMagnArray();
				tIntrct->StoreAuxOldArrays();
				radTRelaxationMethNo_4 RelaxMethNo_4(tIntrct);
				ActualIterNum = RelaxMethNo_4.AutoRelax(PrecOnMagnetiz, MaxIterNumber, 1);
				if(ActualIterNum >= MaxIterNumber) { Send.WarningMessage("Radia::Warning015");}

				tIntrct->SubstractOldMagn();
				UpdateExternFiledInAllIntrctExceptOne(tIntrct, (*it).second);
				tIntrct->AddOldMagn();
				BufTotMisfitMe2 += tIntrct->CalcQuadNewOldMagnDif();
			}
			tIntrct++;
		}
		TotMisfitMe2 = BufTotMisfitMe2/TotAmOfRelaxObj;
		if(TotMisfitMe2 <= DesiredPrecOnMagnetizE2)
		{
			ActualOuterIterNum = i; break;
		}
	}
	if(ActualOuterIterNum == 0) ActualOuterIterNum = MaxIterNumber;

	if(RelaxStatusParamArray != 0)
	{
		RelaxStatusParamArray[0] = sqrt(TotMisfitMe2);

		double MaxModM = 0, MaxModH = 0;
        tIntrct = IntrctPtr;
        for(int j=0; j<mAmOfParts; j++)
		{
			double LocMaxModM = 0, LocMaxModH = 0;
            (tIntrct++)->FindMaxModMandH(LocMaxModM, LocMaxModH);
			if(MaxModM < LocMaxModM) MaxModM = LocMaxModM;
			if(MaxModH < LocMaxModH) MaxModH = LocMaxModH;
		}
		RelaxStatusParamArray[1] = MaxModM;
        RelaxStatusParamArray[2] = MaxModH;
	}
	return ActualOuterIterNum;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_6::UpdateExternFiledInAllIntrctExceptOne(radTInteraction* pIntrctToSkip, const radThg& hgSrc)
{
	if((mAmOfParts <= 0) || (IntrctPtr == 0)) return;

	radTInteraction *tIntrct = IntrctPtr;
    for(int j=0; j<mAmOfParts; j++)
	{
		if(tIntrct == pIntrctToSkip) 
		{
			tIntrct++; continue;
		}
        (tIntrct++)->AddMoreExternField(hgSrc);
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::SetupMainInteractionData(const radThg& hg, const radTCompCriterium& CompCrit)
{
	radTg3d* g3dPtr = radTCast::g3dCast(hg.rep); 
	if(g3dPtr==0) { radTSend::ErrorMessage("Radia::Error003"); throw 0;}

	mArrAuxQuasiExtField = NULL;

	DeleteInterMatrData();
	radThg hEmpty;
	IntrctPtr = new radTInteraction(hg, hEmpty, CompCrit, 1, 1, 1);
	if(IntrctPtr->SomethingIsWrong)
	{
        delete[] IntrctPtr; IntrctPtr = NULL;
	}
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_7::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, double* RelaxStatusParamArray)
{
	if(IntrctPtr == NULL) return 0;
	int AmOfRelaxElem = IntrctPtr->OutAmOfRelaxObjs();
	if(AmOfRelaxElem <= 0) return 0;

	int AmOfSubMatr = 0;
	int *TotArrSubMatrNos = new int[AmOfRelaxElem];
	int *SubMatrLengths=NULL;
	int AmOfInitIter = FillInSubMatrixArrays(PrecOnMagnetiz, TotArrSubMatrNos, SubMatrLengths, AmOfSubMatr);
	if((TotArrSubMatrNos == NULL) || (SubMatrLengths == NULL) || (AmOfSubMatr <= 0)) return 0;

	mArrAuxQuasiExtField = new TVector3d[AmOfRelaxElem];

	CalcQuasiExtFieldForAll(TotArrSubMatrNos, SubMatrLengths, AmOfSubMatr);

	double DesiredPrecOnMagnetizE2 = PrecOnMagnetiz*PrecOnMagnetiz, InstMisfitMe2 = 1.e+30;
	int IterCount=0;
	try
	{
		double InvNumRelaxElem = 1./((double)AmOfRelaxElem);
		for(IterCount=0; IterCount<MaxIterNumber; IterCount++)
		{
			//IntrctPtr->StoreAuxOldMagnArray();
			IntrctPtr->StoreAuxOldArrays();
            int *tSubMatrLengths = SubMatrLengths;
            int *tTotArrSubMatrNos = TotArrSubMatrNos;
			int OffsetCurrentSubMatr = 0;

			for(int i=0; i<AmOfSubMatr; i++)
			{
				int SizeCurrentSubMatr = *(tSubMatrLengths++);
				int CurIterNum = RelaxCurrentSubMatrix(tTotArrSubMatrNos, SizeCurrentSubMatr, DesiredPrecOnMagnetizE2, MaxIterNumber);
				if(CurIterNum >= MaxIterNumber) radTSend::WarningMessage("Radia::Warning015");
				UpdateQuasiExtFieldFromCurrentSubMatrix(TotArrSubMatrNos, OffsetCurrentSubMatr, SizeCurrentSubMatr);

				OffsetCurrentSubMatr += SizeCurrentSubMatr;
                tTotArrSubMatrNos += SizeCurrentSubMatr;
			}
			InstMisfitMe2 = (IntrctPtr->CalcQuadNewOldMagnDif())*InvNumRelaxElem;
			if(InstMisfitMe2 <= DesiredPrecOnMagnetizE2) break;
			if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
		}

		if(RelaxStatusParamArray != 0)
		{
			RelaxStatusParamArray[0] = sqrt(InstMisfitMe2);
			IntrctPtr->FindMaxModMandH(RelaxStatusParamArray[1], RelaxStatusParamArray[2]);
		}
	}
	catch(int ErrNo)
	{
		DeleteArrays_AutoRelax(TotArrSubMatrNos, SubMatrLengths);
		throw ErrNo;
	}

	DeleteArrays_AutoRelax(TotArrSubMatrNos, SubMatrLengths);
	return IterCount;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::CalcQuasiExtFieldForAll(int* TotArrSubMatrNos, int* SubMatrLengths, int AmOfSubMatr)
{
	if((TotArrSubMatrNos == NULL) || (SubMatrLengths == NULL) || (AmOfSubMatr == 0) || (mArrAuxQuasiExtField == NULL)) return;

	int *tTotArrSubMatrNos = TotArrSubMatrNos;
	int *tSubMatrLengths = SubMatrLengths;

	int StartOffset = 0;
	for(int i=0; i<AmOfSubMatr; i++)
	{
		int CurSubMatrSize = *(tSubMatrLengths++);
		for(int j=0; j<CurSubMatrSize; j++)
		{
            int CurInd = *(tTotArrSubMatrNos++);
			CalcQuasiExtFieldForOneElem(CurInd, TotArrSubMatrNos, StartOffset, CurSubMatrSize);
		}
		StartOffset += CurSubMatrSize;
	}
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::CalcQuasiExtFieldForOneElem(int CurInd, int* TotArrSubMatrNos, int StartOffsetSkip, int SkipLength)
{
    int AmOfRelaxElem = IntrctPtr->OutAmOfRelaxObjs();
	if(StartOffsetSkip > AmOfRelaxElem) StartOffsetSkip = AmOfRelaxElem;

	TVector3d &QuasiExtField = *(mArrAuxQuasiExtField + CurInd);
	QuasiExtField.Zero();
	TMatrix3df *MatrArrayPtr = (IntrctPtr->InteractMatrix)[CurInd]; //OC250504
	//TMatrix3d *MatrArrayPtr = (IntrctPtr->InteractMatrix)[CurInd]; //OC250504
	TVector3d *MagnAr = IntrctPtr->NewMagnArray;

    int *tTotArrSubMatrNos = TotArrSubMatrNos;
    for(int i=0; i<StartOffsetSkip; i++)
	{
		int CurElemInd = *(tTotArrSubMatrNos++);
        QuasiExtField += MatrArrayPtr[CurElemInd] * MagnAr[CurElemInd];
	}
	int OffsetStart = StartOffsetSkip + SkipLength;
	if(OffsetStart > AmOfRelaxElem) OffsetStart = AmOfRelaxElem;
    tTotArrSubMatrNos = TotArrSubMatrNos + OffsetStart;
    for(int j=OffsetStart; j<AmOfRelaxElem; j++)
	{
		int CurElemInd = *(tTotArrSubMatrNos++);
        QuasiExtField += MatrArrayPtr[CurElemInd] * MagnAr[CurElemInd];
	}
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::AddQuasiExtFieldFromOneElem(int SrcElemInd, int* TotArrSubMatrNos, int StartOffsetSkip, int SkipLength)
{
    int AmOfRelaxElem = IntrctPtr->OutAmOfRelaxObjs();
	if(StartOffsetSkip > AmOfRelaxElem) StartOffsetSkip = AmOfRelaxElem;

	TMatrix3df **InteractMatr = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d **InteractMatr = IntrctPtr->InteractMatrix; //OC250504
	TVector3d SrcMagn = *(IntrctPtr->NewMagnArray + SrcElemInd);
    int *tTotArrSubMatrNos = TotArrSubMatrNos;

    for(int i=0; i<StartOffsetSkip; i++)
	{
		int CurElemInd = *(tTotArrSubMatrNos++);
        mArrAuxQuasiExtField[CurElemInd] += (*(InteractMatr + CurElemInd))[SrcElemInd] * SrcMagn;
	}
	int OffsetStart = StartOffsetSkip + SkipLength;
	if(OffsetStart > AmOfRelaxElem) OffsetStart = AmOfRelaxElem;
    tTotArrSubMatrNos = TotArrSubMatrNos + OffsetStart;
    for(int j=OffsetStart; j<AmOfRelaxElem; j++)
	{
		int CurElemInd = *(tTotArrSubMatrNos++);
        mArrAuxQuasiExtField[CurElemInd] += (*(InteractMatr + CurElemInd))[SrcElemInd] * SrcMagn;
	}
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_7::RelaxCurrentSubMatrix(int* pTotArrSubMatrNos, int SubMatrSize, double PrecOnMagnetizE2, int MaxIterNumForSubMatr)
{
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

	radTg3dRelax* g3dRelaxPtr = NULL;
	radTMaterial* MaterPtr = NULL;
	TVector3d Mnew_mi_MoldVect, QuasiExtFieldAtElemStrNo;

	double NormFact = 1./double(SubMatrSize);
	double BestPrecMagnE2 =  PrecOnMagnetizE2*NormFact;
	double InstMisfitMe2 = 1.E+23;
    int IterCount;
	for(IterCount = 0; IterCount < MaxIterNumForSubMatr; IterCount++)
	{
        double BufMisfitM=0.;
        double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.25*NormFact*InstMisfitMe2;
        if(LocPrecMagnE2 < BestPrecMagnE2) LocPrecMagnE2 = BestPrecMagnE2;

		int *tTotArrSubMatrNos = pTotArrSubMatrNos;
		for(int RelStrNo=0; RelStrNo<SubMatrSize; RelStrNo++)
		{
			int StrNo = *(tTotArrSubMatrNos++);
            QuasiExtFieldAtElemStrNo.Zero();
            TMatrix3df* MatrArrayPtr = IntrcMat[StrNo]; //OC250504
            //TMatrix3d* MatrArrayPtr = IntrcMat[StrNo]; //OC250504

			int *tColTotArrSubMatrNos = pTotArrSubMatrNos;
            for(int RelColNo=0; RelColNo<SubMatrSize; RelColNo++)
			{
				int ColNo = *(tColTotArrSubMatrNos++);
				if(RelColNo != RelStrNo) QuasiExtFieldAtElemStrNo += (MatrArrayPtr[ColNo] * MagnAr[ColNo]);
			}
            QuasiExtFieldAtElemStrNo += (ExternFieldAr[StrNo] + mArrAuxQuasiExtField[StrNo]);

			g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
			MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

            TVector3d *pInstantH = NewFieldAr + StrNo;
			MaterPtr->FindNewH(*pInstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2);

            TVector3d *pInstantM = MagnAr + StrNo;
            *pInstantM = MaterPtr->M(*pInstantH);

            Mnew_mi_MoldVect = *pInstantM - g3dRelaxPtr->Magn;
            BufMisfitM += Mnew_mi_MoldVect.AmpE2();
            g3dRelaxPtr->Magn = *pInstantM; 
		}
		InstMisfitMe2 = BufMisfitM/SubMatrSize;
		if(InstMisfitMe2 <= PrecOnMagnetizE2) break;
		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
	}
	return IterCount;
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_7::FillInSubMatrixArrays(double PrecOnMagnetiz, int*& TotArrSubMatrNos, int*& SubMatrLengths, int& AmOfSubMatr)
{
	if(IntrctPtr == NULL) return 0;
	int AmOfRelaxElem = IntrctPtr->OutAmOfRelaxObjs();
	if(AmOfRelaxElem <= 0) return 0;

	int AmOfInitIter = 10; //to tune
	double ApproxPortionElemInSubMatr = 0.1; //0.03 //to tune

	radTRelaxationMethNo_4 RelaxMethNo_4(IntrctPtr);
    int ActualIterNum = RelaxMethNo_4.AutoRelax(PrecOnMagnetiz, AmOfInitIter);
	if(ActualIterNum < AmOfInitIter) return ActualIterNum;

	radTlAuxIndNorm *ArrAuxIndNorm = new radTlAuxIndNorm[AmOfRelaxElem];
	radTlAuxIndNorm *tAuxIndNorm = ArrAuxIndNorm;

    //radTg3dRelax *g3dRelaxPtr = NULL;
    TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
    //TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	TMatrix3df *MatrArrayPtr = NULL; //OC250504
	//TMatrix3d *MatrArrayPtr = NULL; //OC250504

    TVector3d ContribH;
	for(int StrNo=0; StrNo<AmOfRelaxElem; StrNo++)
	{
		MatrArrayPtr = IntrcMat[StrNo];
		TVector3d *tMagnAr = IntrctPtr->NewMagnArray;

		for(int ColNo=0; ColNo<AmOfRelaxElem; ColNo++)
		{
            //g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[ColNo];
			if(ColNo != StrNo)
			{
                ContribH = (*MatrArrayPtr) * (*tMagnAr);
				radTAuxIndNorm AuxIndNorm(ColNo, ContribH.AmpE2());
				tAuxIndNorm->push_back(AuxIndNorm);
			}
			tMagnAr++;
            MatrArrayPtr++;
		}
		tAuxIndNorm->sort(radTAuxIndNorm::greater);

		//DEBUG
		//for(radTlAuxIndNorm::iterator it = tAuxIndNorm->begin(); it != tAuxIndNorm->end(); ++it)
		//{
		//	radTAuxIndNorm test = *it;
		//	int aha = 1;
		//}
		//END DBUG

		tAuxIndNorm++;
	}

	int ApproxAmOfElemInSubMatr = (int)(AmOfRelaxElem*ApproxPortionElemInSubMatr);
	if(ApproxAmOfElemInSubMatr < 1) ApproxAmOfElemInSubMatr = 1;

	AmOfSubMatr = (int)(AmOfRelaxElem/ApproxAmOfElemInSubMatr);
	if(AmOfSubMatr <= 0) AmOfSubMatr = 1;

	radTvInt* ArrVectSubMatrNos = new radTvInt[AmOfSubMatr];

	radTlAuxIndNorm::iterator it1, it2;
	tAuxIndNorm = ArrAuxIndNorm;
	int SubMatrixCount = 0;

	int SubMatrCount = 0;
	radTvInt *tArrVectSubMatrNos = ArrVectSubMatrNos;

	for(int i=0; i<AmOfRelaxElem; i++)
	{
		AddElemToCurrentSubMatrixIfNecessary(i, ArrVectSubMatrNos, SubMatrCount + 1);

        int LocElemCount = 0;
		for(it1 = tAuxIndNorm->begin(); it1 != tAuxIndNorm->end(); ++it1)
		{
			if(++LocElemCount > ApproxAmOfElemInSubMatr) break;
			int IndToCheck = (*it1).mInd;

			radTlAuxIndNorm *pAuxIndNormToCheck = ArrAuxIndNorm + IndToCheck;
            int LocElemCount2 = 0;
			for(it2 = pAuxIndNormToCheck->begin(); it2 != pAuxIndNormToCheck->end(); ++it2)
			{
				if(++LocElemCount2 > ApproxAmOfElemInSubMatr) break;
				if((*it2).mInd == i)
				{
					AddElemToCurrentSubMatrixIfNecessary(IndToCheck, ArrVectSubMatrNos, SubMatrCount + 1);
                    break;
				}
			}
		}
		tAuxIndNorm++;

		int AmOfElemsInCurSubMatr = (int)(tArrVectSubMatrNos->size());
		if(AmOfElemsInCurSubMatr >= ApproxAmOfElemInSubMatr)
		{
			if(SubMatrCount < AmOfSubMatr)
			{
				tArrVectSubMatrNos++;
				SubMatrCount++;
			}
		}
	}
    AmOfSubMatr = SubMatrCount;
	if(!tArrVectSubMatrNos->empty()) AmOfSubMatr++;

	tAuxIndNorm = ArrAuxIndNorm;
	for(int j=0; j<AmOfRelaxElem; j++)
	{
		if(CheckIfElemIsPresentInAnySubMatr(j, ArrVectSubMatrNos, AmOfSubMatr)) continue;
		AddElemToAppropriateSubMatrix(j, tAuxIndNorm, ApproxAmOfElemInSubMatr, ArrVectSubMatrNos, AmOfSubMatr);
		tAuxIndNorm++;
	}

	SubMatrLengths = new int[AmOfSubMatr];
	CopyVectSubMatrDataToArrays(ArrVectSubMatrNos, AmOfSubMatr, TotArrSubMatrNos, SubMatrLengths);

	if(ArrAuxIndNorm != NULL)
	{
		for(int k=0; k<AmOfRelaxElem; k++) ArrAuxIndNorm[k].erase(ArrAuxIndNorm[k].begin(), ArrAuxIndNorm[k].end());
		delete[] ArrAuxIndNorm; ArrAuxIndNorm = NULL;
	}
	if((ArrVectSubMatrNos != NULL) && (AmOfSubMatr > 0))
	{
        radTvInt *tArrVectSubMatrNos = ArrVectSubMatrNos;
		for(int i=0; i<AmOfSubMatr; i++) 
		{
			tArrVectSubMatrNos->erase(tArrVectSubMatrNos->begin(), tArrVectSubMatrNos->end());
			tArrVectSubMatrNos++;
		}
		delete[] ArrVectSubMatrNos;
		ArrVectSubMatrNos = NULL;
	}
	return ActualIterNum;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::FindSubMatricesToWhichElemCanBeAdded(radTlAuxIndNorm* pAuxIndNorm, int ApproxAmOfElemInSubMatr, radTvInt* ArrVectSubMatrNos, int AmOfSubMatr, radTvInt& VectIndPossibleSubMatr)
{
	if((ArrVectSubMatrNos == 0) || (AmOfSubMatr == 0) || (ApproxAmOfElemInSubMatr <= 0)) return;
	if(pAuxIndNorm == 0) return;

	radTlAuxIndNorm::iterator it;
	int ElemCount = 0;
	for(it = pAuxIndNorm->begin(); it != pAuxIndNorm->end(); ++it)
	{
		if(++ElemCount > ApproxAmOfElemInSubMatr) break;
		int IndToCheck = (*it).mInd;

		radTvInt *tArrVectSubMatrNos = ArrVectSubMatrNos;
		for(int i=0; i<AmOfSubMatr; i++)
		{
			if(CheckIfElemIsPresentInThisSubMatr(IndToCheck, tArrVectSubMatrNos))
			{
				VectIndPossibleSubMatr.push_back(i);
				break;
			}
			tArrVectSubMatrNos++;
		}
	}
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_7::FindSubMatrWithSmallestNumOfElem(radTvInt& VectIndPossibleSubMatr, radTvInt* ArrVectSubMatrNos, int AmOfSubMatr)
{
	if((ArrVectSubMatrNos == 0) || (AmOfSubMatr == 0)) return -1;

    int AmOfSubMatrToWhichElemCanBeAdded = (int)(VectIndPossibleSubMatr.size());
	if(AmOfSubMatrToWhichElemCanBeAdded == 0) return FindSubMatrWithSmallestNumOfElem(ArrVectSubMatrNos, AmOfSubMatr);

	int IndMinSize = -1, MinSize = IntrctPtr->OutAmOfRelaxObjs();
	radTvInt::iterator it;
	for(it = VectIndPossibleSubMatr.begin(); it != VectIndPossibleSubMatr.end(); ++it)
	{
		int IndSubMatr = (*it);
		if(IndSubMatr >= AmOfSubMatr) continue;
		radTvInt *tSubMatr = ArrVectSubMatrNos + IndSubMatr;
		int CurSize = (int)(tSubMatr->size());
        if(MinSize > CurSize) { MinSize = CurSize; IndMinSize = IndSubMatr;}
	}
	return IndMinSize;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
