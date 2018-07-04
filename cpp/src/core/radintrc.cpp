/*-------------------------------------------------------------------------
*
* File name:      radintrc.cpp
*
* Project:        RADIA
*
* Description:    Magnetic interaction between "relaxable" field source objects
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radintrc.h"
#include "radsbdrc.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTInteraction::radTInteraction(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char ExtraExternFieldArrayIsNeeded, char KeepTransData)
{
	if(!Setup(In_hg, In_hgMoreExtSrc, InCompCriterium, InMemAllocTotAtOnce, ExtraExternFieldArrayIsNeeded, KeepTransData)) 
	{
		SomethingIsWrong = 1;
		Send.ErrorMessage("Radia::Error118");
		throw 0;
	}
}

//-------------------------------------------------------------------------

radTInteraction::radTInteraction()
{
	AmOfMainElem = 0;
	AmOfExtElem = 0;
	InteractMatrix = NULL;
	ExternFieldArray = NULL;
	AuxOldMagnArray = NULL;
	AuxOldFieldArray = NULL;

	NewMagnArray = NULL;
	NewFieldArray = NULL;
	IdentTransPtr = NULL;

	RelaxSubIntervArray = NULL; // New
	mKeepTransData = 0;
}

//-------------------------------------------------------------------------

int radTInteraction::Setup(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char AuxOldMagnArrayIsNeeded, char KeepTransData)
{
	SomethingIsWrong = 0;

	AmOfMainElem = 0;
	AmOfExtElem = 0;
	InteractMatrix = NULL;
	ExternFieldArray = NULL;
	AuxOldMagnArray = NULL;
	AuxOldFieldArray = NULL;

	NewMagnArray = NULL;
	NewFieldArray = NULL;
	IdentTransPtr = NULL;

	RelaxSubIntervArray = NULL; // New
	AmOfRelaxSubInterv = 0; // New

	SourceHandle = In_hg;
	CompCriterium = InCompCriterium;
	FillInMainTransOnly = 0;
	RelaxationStarted = 0;

	MoreExtSourceHandle = In_hgMoreExtSrc;

	MemAllocTotAtOnce = InMemAllocTotAtOnce;

	IdentTransPtr = new radIdentTrans();

	radTlphgPtr NewListOfTransPtr;
	CountMainRelaxElems((radTg3d*)(SourceHandle.rep), &NewListOfTransPtr);

	if(!NotEmpty()) return 0;
	AllocateMemory(AuxOldMagnArrayIsNeeded);
	if(SomethingIsWrong)
	{
		EmptyVectOfPtrToListsOfTrans(); return 0;
	}

	FillInRelaxSubIntervArray(); // New

	FillInMainTransPtrArray();
	SetupInteractMatrix();
	SetupExternFieldArray();

	AddExternFieldFromMoreExtSource();
	//ZeroAuxOldMagnArray();
	ZeroAuxOldArrays();

	mKeepTransData = KeepTransData;
	if(!KeepTransData) //OC021103
	{
        DestroyMainTransPtrArray();
        EmptyVectOfPtrToListsOfTrans();
	}

	//ResetM();
	InitAuxArrays();

	return 1;
}

//-------------------------------------------------------------------------

radTInteraction::~radTInteraction()
{
	if(MemAllocTotAtOnce)
	{
		if(InteractMatrix != NULL)
		{
			if(InteractMatrix[0] != NULL) delete[] (InteractMatrix[0]);
			delete[] InteractMatrix;
		}
	}
	else
	{
		if(InteractMatrix != NULL)
		{
			for(int i=0; i<AmOfMainElem; i++)
			{
				TMatrix3df* Matrix3dPtr = InteractMatrix[i]; //OC250504
				//TMatrix3d* Matrix3dPtr = InteractMatrix[i]; //OC250504
				if(Matrix3dPtr != NULL) delete[] Matrix3dPtr;
			}
			delete[] InteractMatrix;
		}
	}

	g3dExternPtrVect.erase(g3dExternPtrVect.begin(), g3dExternPtrVect.end()); //OC240408, to enable current scaling/update

	if(ExternFieldArray != NULL) delete[] ExternFieldArray;
	if(AuxOldMagnArray != NULL) delete[] AuxOldMagnArray;
	if(AuxOldFieldArray != NULL) delete[] AuxOldFieldArray;

	if(NewMagnArray != NULL) delete[] NewMagnArray;
	if(NewFieldArray != NULL) delete[] NewFieldArray;

	if(RelaxSubIntervArray != NULL) delete[] RelaxSubIntervArray;

	if(mKeepTransData) //OC021103
	{
        DestroyMainTransPtrArray();
        EmptyVectOfPtrToListsOfTrans();
	}
	if(IdentTransPtr != NULL) delete IdentTransPtr; //required by EmptyVectOfPtrToListsOfTrans();
}

//-------------------------------------------------------------------------

void radTInteraction::CountMainRelaxElems(radTg3d* g3dPtr, radTlphgPtr* CurrListOfTransPtrPtr)
{
	radTGroup* GroupPtr = Cast.GroupCast(g3dPtr);
	if(GroupPtr == 0)
	{
		radTg3dRelax* g3dRelaxPtr = Cast.g3dRelaxCast(g3dPtr);
		if((g3dRelaxPtr != 0) && (g3dRelaxPtr->MaterHandle.rep != 0))
		{
			g3dRelaxPtrVect.push_back(g3dRelaxPtr);
			AmOfMainElem++;

			radTlphgPtr* TotalListOfElemTransPtrPtr = new radTlphgPtr(*CurrListOfTransPtrPtr);
			PushFrontNativeElemTransList(g3dRelaxPtr, TotalListOfElemTransPtrPtr);
			IntVectOfPtrToListsOfTransPtr.push_back(TotalListOfElemTransPtrPtr);
		}
		else 
		{
			g3dExternPtrVect.push_back(g3dPtr);
			AmOfExtElem++;

			radTlphgPtr* TotalListOfElemTransPtrPtr	= new radTlphgPtr(*CurrListOfTransPtrPtr);
			PushFrontNativeElemTransList(g3dPtr, TotalListOfElemTransPtrPtr);
			ExtVectOfPtrToListsOfTransPtr.push_back(TotalListOfElemTransPtrPtr);
		}
	}
	else
	{
		//--New
		radTSubdividedRecMag* SubdividedRecMagPtr = Cast.SubdividedRecMagCast(GroupPtr);
		if(SubdividedRecMagPtr != 0)
		{
			radTg3dRelax* g3dRelaxFromSbdRecMagPtr = (radTg3dRelax*)SubdividedRecMagPtr;

			radTRecMag* SubElRecMagPtr = Cast.RecMagCast((radTg3dRelax*)((*(SubdividedRecMagPtr->GroupMapOfHandlers.begin())).second.rep));

			if((g3dRelaxFromSbdRecMagPtr->MaterHandle.rep != 0) && (SubElRecMagPtr != 0))
			{
				int SubIntervStart = AmOfMainElem;
				if(SubdividedRecMagPtr->FldCmpMeth==1)
				{
					for(int ix=0; ix<int(SubdividedRecMagPtr->kx); ix++)
						for(int iy=0; iy<int(SubdividedRecMagPtr->ky); iy++)
							for(int iz=0; iz<int(SubdividedRecMagPtr->kz); iz++)
							{
								g3dRelaxPtrVect.push_back(g3dRelaxFromSbdRecMagPtr);
								AmOfMainElem++;

								radTlphgPtr* TotalListOfElemTransPtrPtr = new radTlphgPtr(*CurrListOfTransPtrPtr);
								PushFrontNativeElemTransList(g3dRelaxFromSbdRecMagPtr, TotalListOfElemTransPtrPtr);
								IntVectOfPtrToListsOfTransPtr.push_back(TotalListOfElemTransPtrPtr);
							}
				}
				int SubIntervFin = SubIntervStart + (int)(SubdividedRecMagPtr->GroupMapOfHandlers.size()) - 1;

				if(RelaxSubIntervConstrVect.empty())
				{
					radTRelaxSubInterval RlxSbIntrv(SubIntervStart, SubIntervFin, RelaxTogether);
					RelaxSubIntervConstrVect.push_back(RlxSbIntrv);
				}
				else
				{
					radTRelaxSubInterval& LastEnteredSubIntrv = RelaxSubIntervConstrVect.back();
					if((SubIntervStart != LastEnteredSubIntrv.StartNo) && (SubIntervFin != LastEnteredSubIntrv.FinNo))
					{
						radTRelaxSubInterval RlxSbIntrv(SubIntervStart, SubIntervFin, RelaxTogether);
						RelaxSubIntervConstrVect.push_back(RlxSbIntrv);
					}
				}
			}
		}
		if((SubdividedRecMagPtr == 0) || ((SubdividedRecMagPtr != 0) && (SubdividedRecMagPtr->FldCmpMeth != 1)))
		{
		//--EndNew
			radTlphgPtr* LocListOfTransPtrPtr = CurrListOfTransPtrPtr;
			
			short GroupListOfTransIsNotEmpty = 1;
			if(GroupPtr->g3dListOfTransform.empty()) GroupListOfTransIsNotEmpty = 0;

			if(GroupListOfTransIsNotEmpty) 
			{
				LocListOfTransPtrPtr = new radTlphgPtr(*CurrListOfTransPtrPtr);
				PushFrontNativeElemTransList(GroupPtr, LocListOfTransPtrPtr);
			}

			for(radTmhg::iterator iter = GroupPtr->GroupMapOfHandlers.begin();
				iter != GroupPtr->GroupMapOfHandlers.end(); ++iter) 
				CountMainRelaxElems((radTg3d*)((*iter).second.rep), LocListOfTransPtrPtr);

			if(GroupListOfTransIsNotEmpty) delete LocListOfTransPtrPtr;
		//--New
		}
		//--EndNew
	}
}

//-------------------------------------------------------------------------

void radTInteraction::FillInRelaxSubIntervArray() // New
{
	if(RelaxSubIntervConstrVect.size() == 0) return;

	int CurrentStartNo = 0;
	int PlainCount = -1;

#ifdef __GCC__
	vector<radTRelaxSubInterval>::iterator Iter;
#else
	vector<radTRelaxSubInterval, allocator<radTRelaxSubInterval> >::iterator Iter;
#endif

	for(Iter = RelaxSubIntervConstrVect.begin(); Iter != RelaxSubIntervConstrVect.end(); ++Iter)
	{
		int LocStartNo = (*Iter).StartNo;
		if(LocStartNo != CurrentStartNo)
		{
			RelaxSubIntervArray[++PlainCount] = radTRelaxSubInterval(CurrentStartNo, LocStartNo-1, RelaxApart);
		}
		RelaxSubIntervArray[++PlainCount] = *Iter;
		CurrentStartNo = (*Iter).FinNo + 1;
	}
	if(CurrentStartNo != AmOfMainElem)
		RelaxSubIntervArray[++PlainCount] = radTRelaxSubInterval(CurrentStartNo, AmOfMainElem-1, RelaxApart);
	
	AmOfRelaxSubInterv = ++PlainCount;

	RelaxSubIntervConstrVect.erase(RelaxSubIntervConstrVect.begin(), RelaxSubIntervConstrVect.end());
}

//-------------------------------------------------------------------------

void radTInteraction::AllocateMemory(char AuxOldMagnArrayIsNeeded)
{
	//try
	//{
		ExternFieldArray = new TVector3d[AmOfMainElem];
		if(AuxOldMagnArrayIsNeeded) 
		{
			AuxOldMagnArray = new TVector3d[AmOfMainElem];
			AuxOldFieldArray = new TVector3d[AmOfMainElem];
		}

		NewMagnArray = new TVector3d[AmOfMainElem];
		NewFieldArray = new TVector3d[AmOfMainElem];
		InteractMatrix = new TMatrix3df*[AmOfMainElem]; //OC250504
		//InteractMatrix = new TMatrix3d*[AmOfMainElem]; //OC250504

		for(int k=0; k<AmOfMainElem; k++) InteractMatrix[k] = NULL;
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	Send.ErrorMessage(radExceptionPtr->what());	return;
	//}
	//catch (...)
	//{
	//	Send.ErrorMessage("Radia::Error999"); return;
	//}

	if(MemAllocTotAtOnce)
	{
		TMatrix3df* GenMatrPtr = 0; //OC250504
		//TMatrix3d* GenMatrPtr = 0; //OC250504
		//try
		//{
			GenMatrPtr = new TMatrix3df[AmOfMainElem*AmOfMainElem]; //OC250504
			//GenMatrPtr = new TMatrix3d[AmOfMainElem*AmOfMainElem]; //OC250504
		//}
		//catch (radTException* radExceptionPtr)
		//{
		//	InteractMatrix[0] = NULL;
		//	SomethingIsWrong = 1;
		//	Send.ErrorMessage(radExceptionPtr->what());	return;
		//}
		//catch (...)
		//{
		//	Send.ErrorMessage("Radia::Error999"); return;
		//}

		if(GenMatrPtr != 0) // Check for allocation failure
			for(int i=0; i<AmOfMainElem; i++) InteractMatrix[i] = &(GenMatrPtr[i*AmOfMainElem]);
		else
		{
			InteractMatrix[0] = NULL;
			SomethingIsWrong = 1;
			Send.ErrorMessage("Radia::Error900"); return;
		}
	}
	else
	{
		for(int i=0; i<AmOfMainElem; i++)
		{
			InteractMatrix[i] = new TMatrix3df[AmOfMainElem]; //OC250504
			//InteractMatrix[i] = new TMatrix3d[AmOfMainElem]; //OC250504
			if(InteractMatrix[i] == 0) // Check for allocation failure
			{
				for(int k=0; k<i; k++) delete[] (InteractMatrix[i]);
				delete[] InteractMatrix;

				SomethingIsWrong = 1;
				Send.ErrorMessage("Radia::Error900"); return;
			}
		}
	}

	int MaxSubIntervArraySize = 2 * ((int)(RelaxSubIntervConstrVect.size())) + 1; // New
	//try
	//{
		if(MaxSubIntervArraySize > 1) RelaxSubIntervArray = new radTRelaxSubInterval[MaxSubIntervArraySize]; // New
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	Send.ErrorMessage(radExceptionPtr->what());	return;
	//}
	//catch (...)
	//{
	//	Send.ErrorMessage("Radia::Error999"); return;
	//}
}

//-------------------------------------------------------------------------

void radTInteraction::NestedFor_Trans(radTrans* BaseTransPtr, const radTlphgPtr::const_iterator& Iter, int ElemLocInd, char I_or_E)
{
	radTrans* TransPtr = (radTrans*)(((**Iter).Handler_g).rep);
	radTrans* LocTotTransPtr = BaseTransPtr;
	radTrans LocTotTrans;

	radTlphgPtr::const_iterator LocalNextIter = Iter;
	LocalNextIter++;
	int Mult = (**Iter).m;

	if(Mult == 1)
	{
		TrProduct(LocTotTransPtr, TransPtr, LocTotTrans);
		AddTransOrNestedFor(&LocTotTrans, LocalNextIter, ElemLocInd, I_or_E);
	}
	else
	{
		AddTransOrNestedFor(LocTotTransPtr, LocalNextIter, ElemLocInd, I_or_E);
		if(FillInMainTransOnly) return;
		for(int km = 1; km < Mult; km++)
		{
			TrProduct(LocTotTransPtr, TransPtr, LocTotTrans);
			LocTotTransPtr = &LocTotTrans;
			AddTransOrNestedFor(LocTotTransPtr, LocalNextIter, ElemLocInd, I_or_E);
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::FillInMainTransPtrArray()
{
	MainTransPtrArray = new radTrans*[AmOfMainElem];
	FillInMainTransOnly = 1;

	for(int i=0; i<AmOfMainElem; i++)
	{
		FillInTransPtrVectForElem(i, 'I');
		if(Cast.IdentTransCast(TransPtrVect[0]) == 0) 
		{
			MainTransPtrArray[i] = new radTrans(*(TransPtrVect[0]));
		}
		else MainTransPtrArray[i] = IdentTransPtr;
		EmptyTransPtrVect();
	}
	FillInMainTransOnly = 0;
}

//-------------------------------------------------------------------------

int radTInteraction::CountRelaxElemsWithSym()
{
	int AmOfElemWithSym = 0;

	for(int i=0; i<AmOfMainElem; i++)
	{
		radTlphgPtr& Loc_lphgPtr = *(IntVectOfPtrToListsOfTransPtr[i]);
		int LocTotMult = 1;

		for(radTlphgPtr::iterator TrIter = Loc_lphgPtr.begin();	
			TrIter != Loc_lphgPtr.end(); ++TrIter)
		{
			LocTotMult *= (**TrIter).m;
		}
		AmOfElemWithSym += LocTotMult;
	}
	return AmOfElemWithSym;
}

//-------------------------------------------------------------------------

void radTInteraction::SetupInteractMatrix()
{
	radTFieldKey FieldKeyInteract; FieldKeyInteract.B_=FieldKeyInteract.H_=FieldKeyInteract.PreRelax_=1;
	TVector3d ZeroVect(0.,0.,0.);

	//--New
	int AmOfElemWithSym = CountRelaxElemsWithSym();
	//--EndNew

	for(int ColNo=0; ColNo<AmOfMainElem; ColNo++)
	{
		FillInTransPtrVectForElem(ColNo, 'I');
		radTg3dRelax* g3dRelaxPtrColNo = g3dRelaxPtrVect[ColNo];

		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
		{
			TVector3d InitObsPoiVect = MainTransPtrArray[StrNo]->TrPoint((g3dRelaxPtrVect[StrNo])->ReturnCentrPoint());

			TMatrix3d SubMatrix(ZeroVect, ZeroVect, ZeroVect), BufSubMatrix;
			for(unsigned i=0; i<TransPtrVect.size(); i++)
			{
				TVector3d ObsPoiVect = TransPtrVect[i]->TrPoint_inv(InitObsPoiVect);

				radTField Field(FieldKeyInteract, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
				Field.AmOfIntrctElemWithSym = AmOfElemWithSym; // New, may be changed later

				g3dRelaxPtrColNo->B_comp(&Field);

				BufSubMatrix.Str0 = Field.B; 
				BufSubMatrix.Str1 = Field.H; 
				BufSubMatrix.Str2 = Field.A; 

				TransPtrVect[i]->TrMatrix(BufSubMatrix);
				SubMatrix += BufSubMatrix;
			}
			MainTransPtrArray[StrNo]->TrMatrix_inv(SubMatrix);
			InteractMatrix[StrNo][ColNo] = SubMatrix;
		}
		EmptyTransPtrVect();
	}

	//--New
	for(int ClNo=0; ClNo<AmOfMainElem; ClNo++)
	{
		radTg3dRelax* g3dRelaxPtrClNo = g3dRelaxPtrVect[ClNo];
		g3dRelaxPtrVect[ClNo] = g3dRelaxPtrClNo->FormalIntrctMemberPtr();
	}
	//--EndNew
}

//-------------------------------------------------------------------------

void radTInteraction::SetupExternFieldArray()
{
	radTFieldKey FieldKeyExtern; FieldKeyExtern.H_=1;
	TVector3d ZeroVect(0.,0.,0.), InitObsPoiVect(0.,0.,0.), ObsPoiVect(0.,0.,0.);

	for(int k=0; k<AmOfMainElem; k++) ExternFieldArray[k] = ZeroVect;

	for(int ExtElNo=0; ExtElNo<AmOfExtElem; ExtElNo++)
	{
		FillInTransPtrVectForElem(ExtElNo, 'E');
		radTg3d* ExtElPtr = g3dExternPtrVect[ExtElNo];

		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
		{
			InitObsPoiVect = MainTransPtrArray[StrNo]->TrPoint((g3dRelaxPtrVect[StrNo])->CentrPoint);
			TVector3d BufVect(0.,0.,0.);
			for(unsigned i=0; i<TransPtrVect.size(); i++)
			{
				TVector3d ObsPoiVect = TransPtrVect[i]->TrPoint_inv(InitObsPoiVect);
				radTField Field(FieldKeyExtern, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.); // Improve
				ExtElPtr->B_comp(&Field);
				BufVect += TransPtrVect[i]->TrVectField(Field.H);
			}
			ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrVectField_inv(BufVect);
		}
		EmptyTransPtrVect();
	}
	//g3dExternPtrVect.erase(g3dExternPtrVect.begin(), g3dExternPtrVect.end()); //OC240408, to enable current scaling/update
}

//-------------------------------------------------------------------------

void radTInteraction::AddExternFieldFromMoreExtSource()
{
	if(MoreExtSourceHandle.rep != 0)
	{
		radTFieldKey FieldKeyExtern; FieldKeyExtern.H_=1;
		TVector3d ZeroVect(0.,0.,0.), InitObsPoiVect(0.,0.,0.);

		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
		{
			radTrans* ATransPtr = MainTransPtrArray[StrNo];

			InitObsPoiVect = MainTransPtrArray[StrNo]->TrPoint((g3dRelaxPtrVect[StrNo])->CentrPoint);
			radTField Field(FieldKeyExtern, CompCriterium, InitObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.); // Improve

			((radTg3d*)(MoreExtSourceHandle.rep))->B_genComp(&Field);

			//TVector3d BufVect = ExternFieldArray[StrNo];

			ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrVectField_inv(Field.H);
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::AddMoreExternField(const radThg& hExtraExtSrc)
{
	if(hExtraExtSrc.rep == 0) return;

	radTg3d* pExtraExtSrc = (radTg3d*)(hExtraExtSrc.rep);

	radTFieldKey FieldKeyExtern; FieldKeyExtern.H_=1;
	TVector3d ZeroVect(0.,0.,0.), InitObsPoiVect(0.,0.,0.);

	for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
	{
		radTrans* aTransPtr = MainTransPtrArray[StrNo];
        InitObsPoiVect = MainTransPtrArray[StrNo]->TrPoint((g3dRelaxPtrVect[StrNo])->CentrPoint);

        radTField Field(FieldKeyExtern, CompCriterium, InitObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.); // Improve
        pExtraExtSrc->B_genComp(&Field);

        ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrVectField_inv(Field.H);
	}
}

//-------------------------------------------------------------------------

void radTInteraction::ZeroAuxOldArrays()
{
	if(AmOfMainElem <= 0) return;

	if(AuxOldMagnArray != NULL)
	{
		TVector3d *tAuxOldMagn = AuxOldMagnArray;
		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
		{
			tAuxOldMagn->x = 0;
			tAuxOldMagn->y = 0;
			(tAuxOldMagn++)->z = 0;
		}
	}
	if(AuxOldFieldArray != NULL)
	{
		TVector3d *tAuxOldField = AuxOldFieldArray;
		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
		{
			tAuxOldField->x = 0;
			tAuxOldField->y = 0;
			(tAuxOldField++)->z = 0;
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::SubstractOldMagn()
{
	if((AuxOldMagnArray == NULL) || (AmOfMainElem <= 0)) return;

	TVector3d *tAuxOldMagn = AuxOldMagnArray;
	for(int StNo=0; StNo<AmOfMainElem; StNo++)
	{
		TVector3d &M = (g3dRelaxPtrVect[StNo])->Magn;
		M -= *(tAuxOldMagn++); 
    }
}

//-------------------------------------------------------------------------

void radTInteraction::AddOldMagn()
{
	if((AuxOldMagnArray == NULL) || (AmOfMainElem <= 0)) return;

	TVector3d *tAuxOldMagn = AuxOldMagnArray;
	for(int StNo=0; StNo<AmOfMainElem; StNo++)
	{
		TVector3d &M = (g3dRelaxPtrVect[StNo])->Magn;
		M += *(tAuxOldMagn++); 
    }
}

//-------------------------------------------------------------------------

double radTInteraction::CalcQuadNewOldMagnDif()
{
	if((AuxOldMagnArray == NULL) || (AmOfMainElem <= 0)) return 0;

	double SumE2 = 0;
	TVector3d *tAuxOldMagn = AuxOldMagnArray;
	for(int StNo=0; StNo<AmOfMainElem; StNo++)
	{
		TVector3d CurDifM = (g3dRelaxPtrVect[StNo])->Magn - *(tAuxOldMagn++); 
		SumE2 += CurDifM.AmpE2(); //CurDifM*CurDifM;
    }
	return SumE2;
}

//-------------------------------------------------------------------------

void radTInteraction::FindMaxModMandH(double& MaxModM, double& MaxModH)
{
	double BufMaxModMe2, BufMaxModHe2, TestBufMaxModMe2, TestBufMaxModHe2;
	BufMaxModMe2 = BufMaxModHe2 = TestBufMaxModMe2 = TestBufMaxModHe2 = 1.E-17;

	for(int i=0; i<AmOfMainElem; i++)
	{
		TVector3d &NewMagn = NewMagnArray[i];
        TestBufMaxModMe2 = NewMagn.x*NewMagn.x + NewMagn.y*NewMagn.y + NewMagn.z*NewMagn.z;
        if(BufMaxModMe2 < TestBufMaxModMe2) BufMaxModMe2 = TestBufMaxModMe2;

		TVector3d &NewField = NewFieldArray[i];
		TestBufMaxModHe2 = NewField.x*NewField.x + NewField.y*NewField.y + NewField.z*NewField.z;
        if(BufMaxModHe2 < TestBufMaxModHe2) BufMaxModHe2 = TestBufMaxModHe2;
	}
	MaxModM = sqrt(BufMaxModMe2);
	MaxModH = sqrt(BufMaxModHe2);
}

//-------------------------------------------------------------------------

void radTInteraction::DumpBinVectOfPtrToListsOfTransPtr(CAuxBinStrVect& oStr, radVectPtr_lphgPtr& VectOfPtrToListsOfTransPtr, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers)
{
	int sizeVectOfPtrToListsOfTransPtr = (int)VectOfPtrToListsOfTransPtr.size();
	oStr << sizeVectOfPtrToListsOfTransPtr;
	for(int i=0; i<sizeVectOfPtrToListsOfTransPtr; i++)
	{
		radTlphgPtr* curListOfElemTransPtrPtr = VectOfPtrToListsOfTransPtr[i];
		int size_curListOfElemTransPtr = 0;
		if(curListOfElemTransPtrPtr != 0) size_curListOfElemTransPtr = (int)curListOfElemTransPtrPtr->size();
		
		oStr << size_curListOfElemTransPtr;
		if(size_curListOfElemTransPtr > 0)
		{
			for(radTlphgPtr::iterator TrIter = curListOfElemTransPtrPtr->begin();	TrIter != curListOfElemTransPtrPtr->end(); ++TrIter)
			{
				radTPair_int_hg *p_m_hg = *TrIter;
				//int mult = 0;

				if(p_m_hg != 0) 
				{
					int mult = p_m_hg->m;
					radThg &hg = p_m_hg->Handler_g;

					int existKey = 0;
					for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
					{
						if(mit->second == hg) { existKey = mit->first; break;}
					}
					oStr << mult;
					oStr << existKey;
				}
				else oStr << (int)0;
			}
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//radThg SourceHandle;
	int existKeySource = 0;
	if(SourceHandle.rep != 0)
	{
		//oStr << (char)1;
		//int existKey = 0;
		//const radThg &cur_hg = iter->second;
		for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
		{
			if(mit->second == SourceHandle) { existKeySource = mit->first; break;}
		}
		if(existKeySource == 0)
		{
			existKeySource = gUniqueMapKey; 
			gMapOfHandlers[gUniqueMapKey++] = SourceHandle;
		}
		int indExist = CAuxParse::FindElemInd(existKeySource, vElemKeysOut);
		if(indExist < 0) SourceHandle.rep->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, existKeySource);
	}
	//else oStr << (char)0;

	//radThg MoreExtSourceHandle;
	int existKeyMoreExtSource = 0;
	if(MoreExtSourceHandle.rep != 0)
	{
		//oStr << (char)1;
		//int existKey = 0;
		//const radThg &cur_hg = iter->second;
		for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
		{
			if(mit->second == MoreExtSourceHandle) { existKeyMoreExtSource = mit->first; break;}
		}
		if(existKeyMoreExtSource == 0)
		{
			existKeyMoreExtSource = gUniqueMapKey; 
			gMapOfHandlers[gUniqueMapKey++] = MoreExtSourceHandle;
		}
		int indExist = CAuxParse::FindElemInd(existKeyMoreExtSource, vElemKeysOut);
		if(indExist < 0) MoreExtSourceHandle.rep->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, existKeyMoreExtSource);
	}
	//else oStr << (char)0;

	//radTVectPtrg3dRelax g3dRelaxPtrVect;
	vector<int> vInd_g3dRelax;
	int size_g3dRelaxPtrVect = (int)g3dRelaxPtrVect.size();
	//oStr << size_g3dRelaxPtrVect;
	for(int i=0; i<size_g3dRelaxPtrVect; i++)
	{
		radTg3dRelax *p_g3dRelax = g3dRelaxPtrVect[i];
		if(p_g3dRelax != 0)
		{
			radTg *p_g = (radTg*)p_g3dRelax;
			//try to find element in the global map by pointer
			int oldKey = 0;
			for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
			{
				if(mit->second.rep == p_g) { oldKey = mit->first; break;}
			}
			if(oldKey == 0)
			{
				oldKey = gUniqueMapKey;
				radThg hg(p_g3dRelax);
				gMapOfHandlers[gUniqueMapKey++] = hg;
			}
			int indExist = CAuxParse::FindElemInd(oldKey, vElemKeysOut);
			if(indExist < 0) p_g3dRelax->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, oldKey);

			vInd_g3dRelax.push_back(oldKey);
		}
	}

	//radTVectPtr_g3d g3dExternPtrVect;
	vector<int> vInd_g3dExternPtrVect;
	int size_g3dExternPtrVect = (int)g3dExternPtrVect.size();
	for(int i=0; i<size_g3dExternPtrVect; i++)
	{
		radTg3d *p_g3d = g3dExternPtrVect[i];
		if(p_g3d != 0)
		{
			radTg *p_g = (radTg*)p_g3d;

			//try to find this element in the global map by pointer
			int oldKey = 0;
			for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
			{
				if(mit->second.rep == p_g) { oldKey = mit->first; break;}
			}
			if(oldKey == 0)
			{
				oldKey = gUniqueMapKey;
				radThg hg(p_g3d);
				gMapOfHandlers[gUniqueMapKey++] = hg;
			}
			int indExist = CAuxParse::FindElemInd(oldKey, vElemKeysOut);
			if(indExist < 0) p_g3d->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, oldKey);

			vInd_g3dExternPtrVect.push_back(oldKey);
		}
	}

	//radTVectPtrTrans TransPtrVect; //not required?
	vector<int> vIndTransPtrVect;
	int size_TransPtrVect = (int)TransPtrVect.size();
	for(int i=0; i<size_TransPtrVect; i++)
	{
		radTrans *pTrans = TransPtrVect[i];
		if(pTrans != 0)
		{
			if(Cast.IdentTransCast(pTrans))
			{
				vIndTransPtrVect.push_back(-1); //indicator of IdentTrans
			}
			else
			{
				radTrans *pTransCopy = new radTrans(*pTrans);

				radThg hg(pTransCopy);
				int oldKey = gUniqueMapKey;
				gMapOfHandlers[gUniqueMapKey++] = hg;
				
				pTransCopy->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, oldKey);
				vIndTransPtrVect.push_back(oldKey);
			}
		}
		else vIndTransPtrVect.push_back(0);
	}

	//radTrans** MainTransPtrArray; //required
	vector<int> vIndMainTrans;
	if(mKeepTransData && (MainTransPtrArray != 0))
	{
		for(int i=0; i<AmOfMainElem; i++)
		{
			radTrans *pTrans = MainTransPtrArray[i];
			if(pTrans != 0)
			{
				if(Cast.IdentTransCast(pTrans))
				{
					vIndTransPtrVect.push_back(-1); //indicator of IdentTrans
				}
				else
				{
					radTrans *pTransCopy = new radTrans(*pTrans);

					radThg hg(pTransCopy);
					int oldKey = gUniqueMapKey;
					gMapOfHandlers[gUniqueMapKey++] = hg;
				
					pTransCopy->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, oldKey);
					vIndMainTrans.push_back(oldKey);
				}
			}
			else vIndMainTrans.push_back(0);
		}
	}

	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)0;
	oStr << (char)0;
	oStr << (char)0;
	oStr << (char)0;

	//int AmOfMainElem;
	oStr << AmOfMainElem;

	//int AmOfExtElem;
	oStr << AmOfExtElem;

	//radThg SourceHandle;
	oStr << existKeySource;

	//radThg MoreExtSourceHandle;
	oStr << existKeyMoreExtSource;

	//radTVectPtrg3dRelax g3dRelaxPtrVect;
	int size_vInd_g3dRelax = (int)vInd_g3dRelax.size();
	oStr << size_vInd_g3dRelax;
	for(int i=0; i<size_vInd_g3dRelax; i++) oStr << vInd_g3dRelax[i];

	//radTVectPtr_g3d g3dExternPtrVect;
	int size_vInd_g3dExternPtrVect = (int)vInd_g3dExternPtrVect.size();
	oStr << size_vInd_g3dExternPtrVect;
	for(int i=0; i<size_vInd_g3dExternPtrVect; i++) oStr << vInd_g3dExternPtrVect[i];

	//radTVectPtrTrans TransPtrVect; //not required?
	int size_vIndTransPtrVect = (int)vIndTransPtrVect.size();
	oStr << size_vIndTransPtrVect;
	for(int i=0; i<size_vIndTransPtrVect; i++) oStr << vIndTransPtrVect[i];

	//radTCompCriterium CompCriterium;
	//short BasedOnPrecLevel; // Actually this is used nowhere at the moment
	oStr << CompCriterium.BasedOnPrecLevel;
	//double AbsPrecB;
	oStr << CompCriterium.AbsPrecB;
	//double AbsPrecA;
	oStr << CompCriterium.AbsPrecA;
	//double AbsPrecB_int;
	oStr << CompCriterium.AbsPrecB_int;
	//double AbsPrecForce;
	oStr << CompCriterium.AbsPrecForce;
	//double AbsPrecTorque;
	oStr << CompCriterium.AbsPrecTorque;
	//double AbsPrecEnergy;
	oStr << CompCriterium.AbsPrecTorque;
	//double AbsPrecTrjCoord;
	oStr << CompCriterium.AbsPrecTrjCoord;
	//double AbsPrecTrjAngle;
	oStr << CompCriterium.AbsPrecTrjAngle;
	//double MltplThresh[4]; // Threshold ratios for 4 diff. orders of multipole approx. at field computation
	oStr << CompCriterium.MltplThresh[0] << CompCriterium.MltplThresh[1] << CompCriterium.MltplThresh[2] << CompCriterium.MltplThresh[3];
	//double WorstRelPrec;
	oStr << CompCriterium.WorstRelPrec;
	//char BasedOnWorstRelPrec; // Used at energy - force computation
	oStr << CompCriterium.BasedOnWorstRelPrec;

	//radTRelaxStatusParam RelaxStatusParam;
	//double MisfitM, MaxModM, MaxModH;
	oStr << RelaxStatusParam.MisfitM;
	oStr << RelaxStatusParam.MaxModM;
	oStr << RelaxStatusParam.MaxModH;

	//short RelaxationStarted;
	oStr << RelaxationStarted;

	//TMatrix3df** InteractMatrix; //OC250504
	////TMatrix3d** InteractMatrix; //OC250504
	if(InteractMatrix != NULL)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++)
		{
			TMatrix3df *pLineInteractMatrix = InteractMatrix[i];
			if(pLineInteractMatrix != NULL)
			{
				oStr << (char)1;
				for(int j=0; j<AmOfMainElem; j++)
				{
					oStr << pLineInteractMatrix[j];
				}
			}
			else oStr << (char)0;
		}
	}
	else oStr << (char)0;

	//TVector3d* ExternFieldArray;
	if(ExternFieldArray != NULL)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << ExternFieldArray[i];
	}
	else oStr << (char)0;

	//TVector3d* NewMagnArray;
	if(NewMagnArray != NULL)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << NewMagnArray[i];
	}
	else oStr << (char)0;

	//TVector3d* NewFieldArray;
	if(NewFieldArray != NULL)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << NewFieldArray[i];
	}
	else oStr << (char)0;

	//TVector3d* AuxOldMagnArray;
	if(AuxOldMagnArray != NULL)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << AuxOldMagnArray[i];
	}
	else oStr << (char)0;

	//TVector3d* AuxOldFieldArray;
	if(AuxOldFieldArray != NULL)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << AuxOldFieldArray[i];
	}
	else oStr << (char)0;

	//radTVectRelaxSubInterval RelaxSubIntervConstrVect; // New
	int sizeRelaxSubIntervConstrVect = (int)RelaxSubIntervConstrVect.size();
	oStr << sizeRelaxSubIntervConstrVect;	
	if(sizeRelaxSubIntervConstrVect > 0)
	{
		for(int i=0; i<sizeRelaxSubIntervConstrVect; i++)
		{
			radTRelaxSubInterval &relaxSubInterval = RelaxSubIntervConstrVect[i];
			oStr << relaxSubInterval.StartNo;
			oStr << relaxSubInterval.FinNo;
			oStr << (int)(relaxSubInterval.SubIntervalID);
		}

		//radTRelaxSubInterval* RelaxSubIntervArray; // New 
		if(RelaxSubIntervArray != NULL)
		{
			int MaxSubIntervArraySize = 2*sizeRelaxSubIntervConstrVect + 1;
			oStr << (int)MaxSubIntervArraySize;
			radTRelaxSubInterval *t_RelaxSubIntervArray = RelaxSubIntervArray;
			for(int i=0; i<MaxSubIntervArraySize; i++)
			{
				oStr << (t_RelaxSubIntervArray->StartNo);
				oStr << (t_RelaxSubIntervArray->FinNo);
				oStr << (int)(t_RelaxSubIntervArray->SubIntervalID);
				t_RelaxSubIntervArray++;
			}
		}
		else oStr << (int)0;
	}

	//radVectPtr_lphgPtr IntVectOfPtrToListsOfTransPtr; //required
	DumpBinVectOfPtrToListsOfTransPtr(oStr, IntVectOfPtrToListsOfTransPtr, gMapOfHandlers);

	//radVectPtr_lphgPtr ExtVectOfPtrToListsOfTransPtr; //required
	DumpBinVectOfPtrToListsOfTransPtr(oStr, ExtVectOfPtrToListsOfTransPtr, gMapOfHandlers);

	//radIdentTrans* IdentTransPtr; //required, but doesn't need to be saved
	//radTCast Cast; //no members?
	//radTSend Send; //no members?

	//short FillInMainTransOnly;
	oStr << FillInMainTransOnly;

	//char mKeepTransData;
	oStr << mKeepTransData;

	//radTrans** MainTransPtrArray; //required
	int size_vIndMainTrans = (int)vIndMainTrans.size();
	oStr << size_vIndMainTrans;
	for(int k=0; k<size_vIndMainTrans; k++) oStr << vIndMainTrans[k];
	
	//int AmOfRelaxSubInterv;
	oStr << AmOfRelaxSubInterv;

	//short SomethingIsWrong;
	oStr << SomethingIsWrong;

	//short MemAllocTotAtOnce;
	oStr << MemAllocTotAtOnce;
}

//-------------------------------------------------------------------------

//void radTInteraction::DumpBinParseSourceHandle(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers, bool do_g3dCast, bool do_g3dRelaxCast, radThg& out_hg)
int radTInteraction::DumpBinParseSourceHandle(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers, bool do_g3dCast, bool do_g3dRelaxCast, radThg& out_hg)
{//move to g3d?
	int oldKey = 0;
	inStr >> oldKey;
	if(oldKey > 0)
	{
		map<int, int>::const_iterator itKey = mKeysOldNew.find(oldKey);
		if(itKey != mKeysOldNew.end())
		{
			int newKey = itKey->second;
			if(newKey > 0)
			{
				radTmhg::const_iterator iter = gMapOfHandlers.find(newKey);
				if(iter != gMapOfHandlers.end())
				{
					radThg hg = (*iter).second;
					if(hg.rep != 0)
					{
						if(do_g3dCast || do_g3dRelaxCast)
						{
							radTg3d *g3dPtr = radTCast::g3dCast(hg.rep);
							if(g3dPtr != 0)
							{
								if(do_g3dRelaxCast)
								{
									if(radTCast::g3dRelaxCast(g3dPtr) != 0) out_hg = hg;
								}
								else out_hg = hg;
							}
						}
						else out_hg = hg;
					}
				}
			}
		}
	}
	return oldKey;
}

//-------------------------------------------------------------------------

void radTInteraction::DumpBinParseVectOfPtrToListsOfTransPtr(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers, radVectPtr_lphgPtr& VectOfPtrToListsOfTransPtr)
{
	int sizeVectOfPtrToListsOfTransPtr = 0;
	inStr >> sizeVectOfPtrToListsOfTransPtr;

	for(int i=0; i<sizeVectOfPtrToListsOfTransPtr; i++)
	{
		int size_curListOfElemTransPtr = 0;
		inStr >> size_curListOfElemTransPtr;

		if(size_curListOfElemTransPtr > 0)
		{
			radTlphgPtr *pCurListOfElemTransPtr = new radTlphgPtr();
			for(int j=0; j<size_curListOfElemTransPtr; j++)
			{
				int mult = 0;
				inStr >> mult;
				if(mult > 0)
				{
					radThg hg;
					DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, false, false, hg);
					pCurListOfElemTransPtr->push_back(new radTPair_int_hg(mult, hg));
				}
			}
			VectOfPtrToListsOfTransPtr.push_back(pCurListOfElemTransPtr);
		}
	}
}

//-------------------------------------------------------------------------

radTInteraction::radTInteraction(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{
	//radIdentTrans* IdentTransPtr; //required
	IdentTransPtr = new radIdentTrans();

	//int AmOfMainElem;
	inStr >> AmOfMainElem;

	//int AmOfExtElem;
	inStr >> AmOfExtElem;

	//radThg SourceHandle;
	DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, true, false, SourceHandle);

	//radThg MoreExtSourceHandle;
	DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, true, false, MoreExtSourceHandle);

	//radTVectPtrg3dRelax g3dRelaxPtrVect;
	int size_g3dRelaxPtrVect = 0;
	inStr >> size_g3dRelaxPtrVect;
	if(g3dRelaxPtrVect.size() > 0) g3dRelaxPtrVect.erase(g3dRelaxPtrVect.begin(), g3dRelaxPtrVect.end()); //?
	for(int i=0; i<size_g3dRelaxPtrVect; i++)
	{
		radThg hg;
		DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, true, true, hg);
		if(hg.rep != 0) g3dRelaxPtrVect.push_back((radTg3dRelax*)((radTg3d*)hg.rep));
	}

	//radTVectPtr_g3d g3dExternPtrVect;
	int size_g3dExternPtrVect = 0;
	inStr >> size_g3dExternPtrVect;
	if(g3dExternPtrVect.size() > 0) g3dExternPtrVect.erase(g3dExternPtrVect.begin(), g3dExternPtrVect.end()); //?
	for(int i=0; i<size_g3dExternPtrVect; i++)
	{
		radThg hg;
		DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, true, false, hg);
		if(hg.rep != 0) g3dExternPtrVect.push_back((radTg3d*)hg.rep);
	}

	//radTVectPtrTrans TransPtrVect; //not required?
	int sizeTransPtrVect = 0;
	inStr >> sizeTransPtrVect;
	if(TransPtrVect.size() > 0) TransPtrVect.erase(TransPtrVect.begin(), TransPtrVect.end()); //?
	for(int i=0; i<sizeTransPtrVect; i++)
	{
		radThg hg;
		int oldKey = DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, false, false, hg);
		if(oldKey < 0) TransPtrVect.push_back(IdentTransPtr);
		else if(hg.rep != 0) TransPtrVect.push_back(new radTrans(*((radTrans*)hg.rep))); //will be deleted at distraction
	}

	//radTCompCriterium CompCriterium;
	//short BasedOnPrecLevel; // Actually this is used nowhere at the moment
	inStr >> CompCriterium.BasedOnPrecLevel;
	//double AbsPrecB;
	inStr >> CompCriterium.AbsPrecB;
	//double AbsPrecA;
	inStr >> CompCriterium.AbsPrecA;
	//double AbsPrecB_int;
	inStr >> CompCriterium.AbsPrecB_int;
	//double AbsPrecForce;
	inStr >> CompCriterium.AbsPrecForce;
	//double AbsPrecTorque;
	inStr >> CompCriterium.AbsPrecTorque;
	//double AbsPrecEnergy;
	inStr >> CompCriterium.AbsPrecTorque;
	//double AbsPrecTrjCoord;
	inStr >> CompCriterium.AbsPrecTrjCoord;
	//double AbsPrecTrjAngle;
	inStr >> CompCriterium.AbsPrecTrjAngle;
	//double MltplThresh[4]; // Threshold ratios for 4 diff. orders of multipole approx. at field computation
	inStr >> CompCriterium.MltplThresh[0];
	inStr >> CompCriterium.MltplThresh[1];
	inStr >> CompCriterium.MltplThresh[2];
	inStr >> CompCriterium.MltplThresh[3];
	//double WorstRelPrec;
	inStr >> CompCriterium.WorstRelPrec;
	//char BasedOnWorstRelPrec; // Used at energy - force computation
	inStr >> CompCriterium.BasedOnWorstRelPrec;

	//radTRelaxStatusParam RelaxStatusParam;
	//double MisfitM, MaxModM, MaxModH;
	inStr >> RelaxStatusParam.MisfitM;
	inStr >> RelaxStatusParam.MaxModM;
	inStr >> RelaxStatusParam.MaxModH;

	//short RelaxationStarted;
	inStr >> RelaxationStarted;

	//TMatrix3df** InteractMatrix;
	char matrixExists = 0;
	inStr >> matrixExists;
	if(matrixExists && (AmOfMainElem > 0))
	{
		InteractMatrix = new TMatrix3df*[AmOfMainElem];
		TMatrix3df **pLineInteractMatrix = InteractMatrix;

		for(int i=0; i<AmOfMainElem; i++)
		{
			char matrixRowExists = 0;
			*pLineInteractMatrix = NULL;

			inStr >> matrixRowExists;
			if(matrixRowExists)
			{
				*pLineInteractMatrix = new TMatrix3df[AmOfMainElem];
				TMatrix3df *tLine = *(pLineInteractMatrix++);
				for(int j=0; j<AmOfMainElem; j++)
				{
					inStr >> *(tLine++);
				}
			}
		}
	}

	//TVector3d* ExternFieldArray;
	char externFieldArrayExists = 0;
	ExternFieldArray = 0;
	inStr >> externFieldArrayExists;
	if(externFieldArrayExists && (AmOfMainElem > 0))
	{
		ExternFieldArray = new TVector3d[AmOfMainElem];
		for(int i=0; i<AmOfMainElem; i++) inStr >> ExternFieldArray[i];
	}

	//TVector3d* NewMagnArray;
	char newMagnArrayExists = 0;
	NewMagnArray = 0;
	inStr >> newMagnArrayExists;
	if(newMagnArrayExists && (AmOfMainElem > 0))
	{
		NewMagnArray = new TVector3d[AmOfMainElem];
		for(int i=0; i<AmOfMainElem; i++) inStr >> NewMagnArray[i];
	}

	//TVector3d* NewFieldArray;
	char newFieldArrayExists = 0;
	NewFieldArray = 0;
	inStr >> newFieldArrayExists;
	if(newFieldArrayExists && (AmOfMainElem > 0))
	{
		NewFieldArray = new TVector3d[AmOfMainElem];
		for(int i=0; i<AmOfMainElem; i++) inStr >> NewFieldArray[i];
	}

	//TVector3d* AuxOldMagnArray;
	char auxOldMagnArrayExists = 0;
	AuxOldMagnArray = 0;
	inStr >> auxOldMagnArrayExists;
	if(auxOldMagnArrayExists && (AmOfMainElem > 0))
	{
		AuxOldMagnArray = new TVector3d[AmOfMainElem];
		for(int i=0; i<AmOfMainElem; i++) inStr >> AuxOldMagnArray[i];
	}

	//TVector3d* AuxOldFieldArray;
	char auxOldFieldArrayExists = 0;
	AuxOldFieldArray = 0;
	inStr >> auxOldFieldArrayExists;
	if(auxOldFieldArrayExists && (AmOfMainElem > 0))
	{
		AuxOldFieldArray = new TVector3d[AmOfMainElem];
		for(int i=0; i<AmOfMainElem; i++) inStr >> AuxOldFieldArray[i];
	}

	//radTVectRelaxSubInterval RelaxSubIntervConstrVect; // New
	int sizeRelaxSubIntervConstrVect = 0;
	RelaxSubIntervArray = 0;
	inStr >> sizeRelaxSubIntervConstrVect;
	if(sizeRelaxSubIntervConstrVect > 0)
	{
		for(int i=0; i<sizeRelaxSubIntervConstrVect; i++)
		{
			radTRelaxSubInterval relaxSubInterval;
			inStr >> relaxSubInterval.StartNo;
			inStr >> relaxSubInterval.FinNo;
			int subIntervalID = 0;
			inStr >> subIntervalID;
			relaxSubInterval.SubIntervalID = (TRelaxSubIntervalID)subIntervalID;

			RelaxSubIntervConstrVect.push_back(relaxSubInterval);
		}

		//radTRelaxSubInterval* RelaxSubIntervArray; // New 
		int MaxSubIntervArraySize = 0;
		inStr >> MaxSubIntervArraySize;
		if(MaxSubIntervArraySize > 0)
		{
			RelaxSubIntervArray = new radTRelaxSubInterval[MaxSubIntervArraySize];
			radTRelaxSubInterval *t_RelaxSubIntervArray = RelaxSubIntervArray;
			for(int i=0; i<MaxSubIntervArraySize; i++)
			{
				inStr >> (t_RelaxSubIntervArray->StartNo);
				inStr >> (t_RelaxSubIntervArray->FinNo);
				int subIntervalID = 0;
				inStr >> subIntervalID;
				t_RelaxSubIntervArray->SubIntervalID = (TRelaxSubIntervalID)subIntervalID;
				t_RelaxSubIntervArray++;
			}
		}
	}

	//radVectPtr_lphgPtr IntVectOfPtrToListsOfTransPtr; //required
	DumpBinParseVectOfPtrToListsOfTransPtr(inStr, mKeysOldNew, gMapOfHandlers, IntVectOfPtrToListsOfTransPtr);

	//radVectPtr_lphgPtr ExtVectOfPtrToListsOfTransPtr; //required
	DumpBinParseVectOfPtrToListsOfTransPtr(inStr, mKeysOldNew, gMapOfHandlers, ExtVectOfPtrToListsOfTransPtr);

	//radTCast Cast; //no members?
	//radTSend Send; //no members?

	//short FillInMainTransOnly;
	inStr >> FillInMainTransOnly;

	//char mKeepTransData;
	inStr >> mKeepTransData;

	//radTrans** MainTransPtrArray; //required
	MainTransPtrArray= 0;
	int size_vIndMainTrans = 0;
	inStr >> size_vIndMainTrans;
	if(size_vIndMainTrans > 0)
	{
		MainTransPtrArray = new radTrans*[AmOfMainElem];

		for(int i=0; i<AmOfMainElem; i++)
		{
			radThg hg;
			int oldKey = DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, false, false, hg);
			if(oldKey < 0) MainTransPtrArray[i] = IdentTransPtr;
			else if(hg.rep != 0) MainTransPtrArray[i] = new radTrans(*((radTrans*)hg.rep)); //will be deleted at distraction
		}
	}

	//int AmOfRelaxSubInterv;
	inStr >> AmOfRelaxSubInterv;

	//short SomethingIsWrong;
	inStr >> SomethingIsWrong;

	//short MemAllocTotAtOnce;
	inStr >> MemAllocTotAtOnce;
}

//-------------------------------------------------------------------------
