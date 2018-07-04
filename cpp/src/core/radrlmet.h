/*-------------------------------------------------------------------------
*
* File name:      radrlmet.h
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

#ifndef __RADRLMET_H
#define __RADRLMET_H

#include "radintrc.h"
#include "radmamet.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTIterativeRelaxMeth {
protected:
	radTInteraction* IntrctPtr;

public:
	radTIterativeRelaxMeth(radTInteraction* InIntrctPtr) { IntrctPtr = InIntrctPtr;}
	radTIterativeRelaxMeth() { IntrctPtr = 0;}

	virtual void DefineNewMagnetizations() {}
	
	void MakeN_iter(int);
	void ComputeRelaxStatusParam(const TVector3d*, const TVector3d*, const TVector3d*);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTSimpleRelaxation : public radTIterativeRelaxMeth {
	double RelaxParam;
	
public:
	radTSimpleRelaxation(radTInteraction* InInteractionPtr, double InRelaxParam) 
	: radTIterativeRelaxMeth(InInteractionPtr) { RelaxParam = InRelaxParam;}

	void DefineNewMagnetizations();
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTRelaxationMethNo_2 : public radTIterativeRelaxMeth {
	double RelaxParam;

public:
	radTRelaxationMethNo_2(radTInteraction* InInteractionPtr, double InRelaxParam)
	: radTIterativeRelaxMeth(InInteractionPtr) { RelaxParam = InRelaxParam;}

	void DefineNewMagnetizations();
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTRelaxationMethNo_3 : public radTIterativeRelaxMeth {
	double InstMisfitM;

public:
	radTRelaxationMethNo_3(radTInteraction* InInteractionPtr) 
	: radTIterativeRelaxMeth(InInteractionPtr) 
	{ 
		IntrctPtr = InInteractionPtr; InstMisfitM = 1.E+23;
	}

	void DefineNewMagnetizations();
	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTRelaxAuxData { //OC140103
	float mAbsDelt;
	float mRelaxPar;
	int mBadPassCounts;

	radTRelaxAuxData()
	{
		Init();
	}
	void Init()
	{
		mAbsDelt = (float)(1.E+23);
		mRelaxPar = 1;
		mBadPassCounts = 0;
	}
	void Update(double NewDelt, int MaxConseqBadCounts, double RelaxParMod, double RelaxParMin)
	{
		float PrevAbsDelt = mAbsDelt;
		mAbsDelt = (float)NewDelt;
		if(NewDelt < PrevAbsDelt) 
		{
			mBadPassCounts = 0; return;
		}
		mBadPassCounts++;
		if((mBadPassCounts >= MaxConseqBadCounts) && (mRelaxPar > RelaxParMin))
		//if(mBadPassCounts > MaxConseqBadCounts)
		{
            mRelaxPar *= (float)RelaxParMod;
            //mRelaxPar = 0.5;
			mBadPassCounts = 0;
		}
	}
};

//-------------------------------------------------------------------------

class radTRelaxationMethNo_4 : public radTIterativeRelaxMeth {
	
	double InstMisfitMe2, DesiredPrecOnMagnetizE2;
	radTRelaxAuxData *mpRelaxAuxData; //OC140103

	double mRelaxPar, mRelaxParMin, mRelaxParModFact, mMisfitE2RatToStartModifRelaxPar, mSysEnergy, mSysEnergyMin;
	int mNumConvergPasses, mNumConvergPassesMax, mNumDivergPasses, mNumDivergPassesMax, mMethNo, mIterCount;
	bool mKeepPrevOldValues, mBadConverg;

	double *mElemVolumeArray; //OC010604
	TVector3d *mOptMagnArray, *mOptFieldArray;

public:
	radTRelaxationMethNo_4(radTInteraction* InInteractionPtr) 
	: radTIterativeRelaxMeth(InInteractionPtr) 
	{
		double DesiredPrecOnMagnetiz = 1.E-03;
		DesiredPrecOnMagnetizE2 = DesiredPrecOnMagnetiz * DesiredPrecOnMagnetiz;
		IntrctPtr = InInteractionPtr; InstMisfitMe2 = 1.E+23;
		
		mRelaxPar = 1; //OC041103
		mRelaxParMin = 1.e-04;
		mRelaxParModFact = 0.985; //0.96;
		mMisfitE2RatToStartModifRelaxPar = 1.; //0.99;

        mNumConvergPasses = 0;
		mNumConvergPassesMax = 30; //to steer
        mNumDivergPasses = 0;
		mNumDivergPassesMax = 5; //to steer
        mMethNo = 1;
		mIterCount = 0;
		mSysEnergy = mSysEnergyMin = 0;
		mKeepPrevOldValues = false;
		mBadConverg = false;

		mpRelaxAuxData = 0; //OC140103
		mElemVolumeArray = 0; //OC010604

		mOptMagnArray = mOptFieldArray = 0;
	}

	~radTRelaxationMethNo_4() 
	{
		DeleteAuxArrays(); //OC140103
	}

	void DefineNewMagnetizations();
	void DefineNewMagnetizationsTest();
	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

	void CorrectMagnAndFieldArraysWithRelaxPar();

    TVector3d FindNewFieldTwoSteps(int StrNo, int MethNo);

	void LpTau(int i, double& q);
	double D(double x) { return x - int(x);}


	void SetupAuxArrays() //OC140103
	{
		DeleteAuxArrays();

		long AmOfMainEl = IntrctPtr->AmOfMainElem;
		mpRelaxAuxData = new radTRelaxAuxData[AmOfMainEl];

		radTRelaxAuxData *tRelaxAuxData = mpRelaxAuxData;
		for(long i=0; i<AmOfMainEl; i++) 
		{
			(tRelaxAuxData++)->Init();
		}
	}
	void DeleteAuxArrays() //OC140103
	{
		if(mpRelaxAuxData != 0) { delete[] mpRelaxAuxData; mpRelaxAuxData = 0;}
	}

	void SetupElemVolumeArray()
	{
		DeleteElemVolumeArray();

		if(IntrctPtr != NULL)
		{
            long AmOfMainEl = IntrctPtr->AmOfMainElem;
			mElemVolumeArray = new double[AmOfMainEl];

			double *tElemVolumeArray = mElemVolumeArray;
			radTg3dRelax* g3dRelaxPtr = NULL;
			for(long i=0; i<AmOfMainEl; i++) 
			{
				g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
				*(tElemVolumeArray++) = g3dRelaxPtr->Volume(); //or VolumeWithSymmetries() ?
			}
		}
	}
	void DeleteElemVolumeArray() //OC140103
	{
		if(mElemVolumeArray != 0) { delete[] mElemVolumeArray; mElemVolumeArray = 0;}
	}

	void SetupOptimValuesArrays()
	{
		DeleteOptimValuesArrays();

		if(IntrctPtr != NULL)
		{
            long AmOfMainEl = IntrctPtr->AmOfMainElem;
			mOptMagnArray = new TVector3d[AmOfMainEl];
			mOptFieldArray = new TVector3d[AmOfMainEl];
		}
	}
	void StoreOptimValuesFromOldArrays()
	{
        TVector3d *tOptMagnArray = mOptMagnArray, *tOptFieldArray = mOptFieldArray;
        TVector3d *tOldMagnAr = IntrctPtr->AuxOldMagnArray, *tOldFieldAr = IntrctPtr->AuxOldFieldArray;
		//TVector3d *tMagnAr = IntrctPtr->NewMagnArray, *tFieldAr = IntrctPtr->NewFieldArray;
        long AmOfMainEl = IntrctPtr->AmOfMainElem;
		for(long i=0; i<AmOfMainEl; i++) 
		{
			//*(tOptMagnArray++) = *(tMagnAr++);
			//*(tOptFieldArray++) = *(tFieldAr++);
			*(tOptMagnArray++) = *(tOldMagnAr++);
			*(tOptFieldArray++) = *(tOldFieldAr++);
		}
	}
	void RestoreOptimValuesToOldArrays()
	{
        TVector3d *tOptMagnArray = mOptMagnArray, *tOptFieldArray = mOptFieldArray;
        TVector3d *tOldMagnAr = IntrctPtr->AuxOldMagnArray, *tOldFieldAr = IntrctPtr->AuxOldFieldArray;
        long AmOfMainEl = IntrctPtr->AmOfMainElem;
		for(long i=0; i<AmOfMainEl; i++) 
		{
			*(tOldMagnAr++) = *(tOptMagnArray++);
			*(tOldFieldAr++) = *(tOptFieldArray++);
		}
	}
	void DeleteOptimValuesArrays()
	{
		if(mOptMagnArray != 0) { delete[] mOptMagnArray; mOptMagnArray = 0;}
		if(mOptFieldArray != 0) { delete[] mOptFieldArray; mOptFieldArray = 0;}
	}
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTRelaxationMethNo_8 : public radTIterativeRelaxMeth {
	
	double mInstMisfitMe2, mDesiredPrecOnMagnetizE2;
	//radTRelaxAuxData *mpRelaxAuxData; //OC140103

	//double mRelaxPar, mRelaxParMin, mRelaxParModFact, mMisfitE2RatToStartModifRelaxPar, mSysEnergy, mSysEnergyMin;
	//int mNumConvergPasses, mNumConvergPassesMax, mNumDivergPasses, mNumDivergPassesMax, mMethNo, mIterCount;
	//bool mKeepPrevOldValues, mBadConverg;

	//double *mElemVolumeArray; //OC010604
	//TVector3d *mOptMagnArray, *mOptFieldArray;

public:

	radTRelaxationMethNo_8(radTInteraction* InInteractionPtr) : radTIterativeRelaxMeth(InInteractionPtr) 
	{
		double DesiredPrecOnMagnetiz = 1.E-03;
		mDesiredPrecOnMagnetizE2 = DesiredPrecOnMagnetiz*DesiredPrecOnMagnetiz;
		IntrctPtr = InInteractionPtr; 
		mInstMisfitMe2 = 1.E+23;

		//mRelaxPar = 1; //OC041103
		//mRelaxParMin = 1.e-04;
		//mRelaxParModFact = 0.985; //0.96;
		//mMisfitE2RatToStartModifRelaxPar = 1.; //0.99;

        //mNumConvergPasses = 0;
		//mNumConvergPassesMax = 30; //to steer
        //mNumDivergPasses = 0;
		//mNumDivergPassesMax = 5; //to steer
        //mMethNo = 1;
		//mIterCount = 0;
		//mSysEnergy = mSysEnergyMin = 0;
		//mKeepPrevOldValues = false;
		//mBadConverg = false;

		//mpRelaxAuxData = 0; //OC140103
		//mElemVolumeArray = 0; //OC010604

		//mOptMagnArray = mOptFieldArray = 0;
	}

	~radTRelaxationMethNo_8() 
	{
		//DeleteAuxArrays(); //OC140103
	}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);
	void DefineNewMagnetizations();


/**
	void DefineNewMagnetizationsTest();

	void CorrectMagnAndFieldArraysWithRelaxPar();

    TVector3d FindNewFieldTwoSteps(int StrNo, int MethNo);

	void LpTau(int i, double& q);
	double D(double x) { return x - int(x);}


	void SetupAuxArrays() //OC140103
	{
		DeleteAuxArrays();

		long AmOfMainEl = IntrctPtr->AmOfMainElem;
		mpRelaxAuxData = new radTRelaxAuxData[AmOfMainEl];

		radTRelaxAuxData *tRelaxAuxData = mpRelaxAuxData;
		for(long i=0; i<AmOfMainEl; i++) 
		{
			(tRelaxAuxData++)->Init();
		}
	}
	void DeleteAuxArrays() //OC140103
	{
		if(mpRelaxAuxData != 0) { delete[] mpRelaxAuxData; mpRelaxAuxData = 0;}
	}

	void SetupElemVolumeArray()
	{
		DeleteElemVolumeArray();

		if(IntrctPtr != NULL)
		{
            long AmOfMainEl = IntrctPtr->AmOfMainElem;
			mElemVolumeArray = new double[AmOfMainEl];

			double *tElemVolumeArray = mElemVolumeArray;
			radTg3dRelax* g3dRelaxPtr = NULL;
			for(long i=0; i<AmOfMainEl; i++) 
			{
				g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
				*(tElemVolumeArray++) = g3dRelaxPtr->Volume(); //or VolumeWithSymmetries() ?
			}
		}
	}
	void DeleteElemVolumeArray() //OC140103
	{
		if(mElemVolumeArray != 0) { delete[] mElemVolumeArray; mElemVolumeArray = 0;}
	}

	void SetupOptimValuesArrays()
	{
		DeleteOptimValuesArrays();

		if(IntrctPtr != NULL)
		{
            long AmOfMainEl = IntrctPtr->AmOfMainElem;
			mOptMagnArray = new TVector3d[AmOfMainEl];
			mOptFieldArray = new TVector3d[AmOfMainEl];
		}
	}
	void StoreOptimValuesFromOldArrays()
	{
        TVector3d *tOptMagnArray = mOptMagnArray, *tOptFieldArray = mOptFieldArray;
        TVector3d *tOldMagnAr = IntrctPtr->AuxOldMagnArray, *tOldFieldAr = IntrctPtr->AuxOldFieldArray;
		//TVector3d *tMagnAr = IntrctPtr->NewMagnArray, *tFieldAr = IntrctPtr->NewFieldArray;
        long AmOfMainEl = IntrctPtr->AmOfMainElem;
		for(long i=0; i<AmOfMainEl; i++) 
		{
			//*(tOptMagnArray++) = *(tMagnAr++);
			//*(tOptFieldArray++) = *(tFieldAr++);
			*(tOptMagnArray++) = *(tOldMagnAr++);
			*(tOptFieldArray++) = *(tOldFieldAr++);
		}
	}
	void RestoreOptimValuesToOldArrays()
	{
        TVector3d *tOptMagnArray = mOptMagnArray, *tOptFieldArray = mOptFieldArray;
        TVector3d *tOldMagnAr = IntrctPtr->AuxOldMagnArray, *tOldFieldAr = IntrctPtr->AuxOldFieldArray;
        long AmOfMainEl = IntrctPtr->AmOfMainElem;
		for(long i=0; i<AmOfMainEl; i++) 
		{
			*(tOldMagnAr++) = *(tOptMagnArray++);
			*(tOldFieldAr++) = *(tOptFieldArray++);
		}
	}
	void DeleteOptimValuesArrays()
	{
		if(mOptMagnArray != 0) { delete[] mOptMagnArray; mOptMagnArray = 0;}
		if(mOptFieldArray != 0) { delete[] mOptFieldArray; mOptFieldArray = 0;}
	}
**/
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTRelaxationMethNo_a5 : public radTIterativeRelaxMeth {
	double InstMisfitM;
	double** AuxMatr1;
	double** AuxMatr2;

	double* AuxArray;

	int AmOfRelaxTogether;
	int SizeOfAuxs;

	radTMathLinAlgEq* MathMethPtr;

public:
	radTRelaxationMethNo_a5(radTInteraction*); 
	~radTRelaxationMethNo_a5(); 

	void DefineNewMagnetizations();
	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTRelaxationMethNo_6 : public radTIterativeRelaxMeth {

	int mAmOfParts;
	radTmhg mMapOfPartHandlers;
	radThg mhGroup;

    radTCast Cast;
	radTSend Send;

public:
	radTRelaxationMethNo_6(const radThg& hObj, const radTCompCriterium& CompCrit) 
	{
		SetupInteractionMatrices(hObj, CompCrit);
	}; 
	~radTRelaxationMethNo_6() 
	{
		if(IntrctPtr != NULL) delete[] IntrctPtr;
		mMapOfPartHandlers.erase(mMapOfPartHandlers.begin(), mMapOfPartHandlers.end());
	}; 

	void SetupInteractionMatrices(const radThg& hObj, const radTCompCriterium& CompCrit);
	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, double* RelaxStatusParamArray);
	void UpdateExternFiledInAllIntrctExceptOne(radTInteraction* pIntrctToSkip, const radThg& hg);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTAuxIndNorm {
	int mInd;
	float mNorm;

	radTAuxIndNorm(int Ind, double Norm)
	{
        mInd = Ind; mNorm = (float)Norm;
	}

	static bool greater(const radTAuxIndNorm&, const radTAuxIndNorm&);

	//friend int operator <(const radTAuxIndNorm&, const radTAuxIndNorm&);
	//friend bool operator >(const radTAuxIndNorm&, const radTAuxIndNorm&);
    //friend int operator ==(const radTAuxIndNorm&, const radTAuxIndNorm&);
    //friend bool operator !=(const radTAuxIndNorm&, const radTAuxIndNorm&);
};

inline bool radTAuxIndNorm::greater(const radTAuxIndNorm& P1, const radTAuxIndNorm& P2)
{
	return (P1.mNorm > P2.mNorm);
}
//inline bool operator <(const radTAuxIndNorm& P1, const radTAuxIndNorm& P2)
//{
//	return (P1.mNorm < P2.mNorm);
//}
//inline bool operator >(const radTAuxIndNorm& P1, const radTAuxIndNorm& P2)
//{
//	return (P1.mNorm > P2.mNorm);
//}
//inline bool operator ==(const radTAuxIndNorm& P1, const radTAuxIndNorm& P2)
//{
//	return ((P1.mNorm == P2.mNorm) && (P1.mInd == P2.mInd));
//}
//inline bool operator !=(const radTAuxIndNorm& P1, const radTAuxIndNorm& P2)
//{
//	return ((P1.mNorm != P2.mNorm) || (P1.mInd != P2.mInd));
//}

//-------------------------------------------------------------------------

#ifdef __GCC__
typedef list <radTAuxIndNorm> radTlAuxIndNorm;
typedef vector <int> radTvInt;
//typedef vector <radTvInt> radTvvInt;
#else
typedef list <radTAuxIndNorm, allocator<radTAuxIndNorm> > radTlAuxIndNorm;
typedef vector <int, allocator<int> > radTvInt;
//typedef vector <radTvInt, allocator<radTvInt> > radTvvInt;
#endif

//-------------------------------------------------------------------------

class radTRelaxationMethNo_7 : public radTIterativeRelaxMeth {

	TVector3d* mArrAuxQuasiExtField;

public:
	radTRelaxationMethNo_7(const radThg& hObj, const radTCompCriterium& CompCrit) : radTIterativeRelaxMeth()
	{
		SetupMainInteractionData(hObj, CompCrit);
	}; 
	~radTRelaxationMethNo_7() 
	{
        DeleteInterMatrData();
	}; 

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, double* RelaxStatusParamArray);
	void SetupMainInteractionData(const radThg& hg, const radTCompCriterium& CompCrit);
	int FillInSubMatrixArrays(double PrecOnMagnetiz, int*& TotArrSubMatrNos, int*& SubMatrLengths, int& AmOfSubMatr);
    void FindSubMatricesToWhichElemCanBeAdded(radTlAuxIndNorm* pAuxIndNorm, int ApproxAmOfElemInSubMatr, radTvInt* ArrVectSubMatrNos, int AmOfSubMatr, radTvInt& VectIndPossibleSubMatr);
    int FindSubMatrWithSmallestNumOfElem(radTvInt& VectIndPossibleSubMatr, radTvInt* ArrVectSubMatrNos, int AmOfSubMatr);
	int RelaxCurrentSubMatrix(int* pTotArrSubMatrNos, int SubMatrSize, double PrecOnMagnetizE2, int MaxIterNumForSubMatr);
	void CalcQuasiExtFieldForAll(int* TotArrSubMatrNos, int* SubMatrLengths, int AmOfSubMatr);
	void CalcQuasiExtFieldForOneElem(int CurInd, int* TotArrSubMatrNos, int StartOffsetSkip, int SkipLength);
    void AddQuasiExtFieldFromOneElem(int SrcElemInd, int* TotArrSubMatrNos, int StartOffsetSkip, int SkipLength);

	void AddElemToAppropriateSubMatrix(int IndElemToAdd, radTlAuxIndNorm* pAuxIndNorm, int ApproxAmOfElemInSubMatr, radTvInt* ArrVectSubMatrNos, int AmOfSubMatr)
	{
		if((ArrVectSubMatrNos == 0) || (AmOfSubMatr == 0)) return;

		radTvInt VectIndPossibleSubMatr;
		FindSubMatricesToWhichElemCanBeAdded(pAuxIndNorm, ApproxAmOfElemInSubMatr, ArrVectSubMatrNos, AmOfSubMatr, VectIndPossibleSubMatr);

		int IndMatrToAddElem = FindSubMatrWithSmallestNumOfElem(VectIndPossibleSubMatr, ArrVectSubMatrNos, AmOfSubMatr);
		if((IndMatrToAddElem < 0) || (IndMatrToAddElem >= AmOfSubMatr)) return;

		(ArrVectSubMatrNos + IndMatrToAddElem)->push_back(IndElemToAdd);
	}
	bool CheckIfElemIsPresentInThisSubMatr(int IndElem, radTvInt* pCurSubMatr)
    {
        if(pCurSubMatr == 0) return false;
        radTvInt::iterator it;
		for(it = pCurSubMatr->begin(); it != pCurSubMatr->end(); ++it) 
		{
			if(*it == IndElem) return true;
		}
		return false;
	}
    void AddElemToCurrentSubMatrixIfNecessary(int IndElemToInsert, radTvInt* ArrSubMatr, int AmOfSubMatr)
	{
		if((ArrSubMatr == 0) || (AmOfSubMatr <= 0)) return;
		if(!CheckIfElemIsPresentInAnySubMatr(IndElemToInsert, ArrSubMatr, AmOfSubMatr)) 
		{
			(ArrSubMatr + (AmOfSubMatr - 1))->push_back(IndElemToInsert);
		}
	}
	bool CheckIfElemIsPresentInAnySubMatr(int ElemInd, radTvInt* ArrVectSubMatrNos, int AmOfSubMatr)
	{
		if((ArrVectSubMatrNos == 0) || (AmOfSubMatr <= 0)) return false;
		radTvInt *tArrVectSubMatrNos = ArrVectSubMatrNos;
		for(int i=0; i<AmOfSubMatr; i++)
		{
			if(CheckIfElemIsPresentInThisSubMatr(ElemInd, tArrVectSubMatrNos++)) return true;
		}
		return false;
	}
	int FindSubMatrWithSmallestNumOfElem(radTvInt* ArrVectSubMatrNos, int AmOfSubMatr)
	{
		if((ArrVectSubMatrNos == 0) || (AmOfSubMatr <= 0)) return -1;
		int IndMinSize = -1, MinSize = IntrctPtr->OutAmOfRelaxObjs();
		radTvInt *tArrVectSubMatrNos = ArrVectSubMatrNos;
		for(int i=0; i<AmOfSubMatr; i++)
		{
			int CurSize = (int)(tArrVectSubMatrNos->size());
			if(MinSize > CurSize) { MinSize = CurSize; IndMinSize = i;}
		}
		return IndMinSize;
	}
	void CopyVectSubMatrDataToArrays(radTvInt* ArrVectSubMatrNos, int AmOfSubMatr, int* TotArrSubMatrNos, int* SubMatrSizes)
	{
		if((ArrVectSubMatrNos == 0) || (AmOfSubMatr == 0) || (TotArrSubMatrNos == 0) || (SubMatrSizes == 0)) return;

		radTvInt::iterator it;
		radTvInt *tArrVectSubMatrNos = ArrVectSubMatrNos;
		int *tTotArrSubMatrNos = TotArrSubMatrNos;
		int *tSubMatrSizes = SubMatrSizes;

		for(int i=0; i<AmOfSubMatr; i++)
		{
			*(tSubMatrSizes++) = (int)(tArrVectSubMatrNos->size());
			for(it = tArrVectSubMatrNos->begin(); it != tArrVectSubMatrNos->end(); ++it) *(tTotArrSubMatrNos++) = *it;
			tArrVectSubMatrNos++;
		}
	}
	void SubstractOldMagnFromCurrentSubMatrixArray(int* pCurrentSubMatrNos, int SizeCurrentSubMatr)
	{//this does not modify IntrctPtr->g3dRelaxPtrVect !
		if((pCurrentSubMatrNos == 0) || (SizeCurrentSubMatr <= 0)) return;
		TVector3d* MagnAr = IntrctPtr->NewMagnArray;
		TVector3d* OldMagnAr = IntrctPtr->AuxOldMagnArray;
		int *tCurrentSubMatrNos = pCurrentSubMatrNos;
		for(int i=0; i<SizeCurrentSubMatr; i++)
		{
			int CurElemInd = *(tCurrentSubMatrNos++);
			MagnAr[CurElemInd] -= OldMagnAr[CurElemInd];
		}
	}
	void AddOldMagnFromCurrentSubMatrixArray(int* pCurrentSubMatrNos, int SizeCurrentSubMatr)
	{//this does not modify IntrctPtr->g3dRelaxPtrVect !
		if((pCurrentSubMatrNos == 0) || (SizeCurrentSubMatr <= 0)) return;
		TVector3d* MagnAr = IntrctPtr->NewMagnArray;
		TVector3d* OldMagnAr = IntrctPtr->AuxOldMagnArray;
		int *tCurrentSubMatrNos = pCurrentSubMatrNos;
		for(int i=0; i<SizeCurrentSubMatr; i++)
		{
			int CurElemInd = *(tCurrentSubMatrNos++);
			MagnAr[CurElemInd] += OldMagnAr[CurElemInd];
		}
	}
	void UpdateQuasiExtFieldFromCurrentSubMatrix(int* TotArrSubMatrNos, int OffsetCurrentSubMatr, int SizeCurrentSubMatr)
	{
		int *pCurrentSubMatrNos = TotArrSubMatrNos + OffsetCurrentSubMatr;
		SubstractOldMagnFromCurrentSubMatrixArray(pCurrentSubMatrNos, SizeCurrentSubMatr);

		int *tTotArrSubMatrNos = TotArrSubMatrNos + OffsetCurrentSubMatr;
		for(int i=0; i<SizeCurrentSubMatr; i++)
		{
			int CurElemInd = *(tTotArrSubMatrNos++);
			AddQuasiExtFieldFromOneElem(CurElemInd, TotArrSubMatrNos, OffsetCurrentSubMatr, SizeCurrentSubMatr);
		}
		AddOldMagnFromCurrentSubMatrixArray(pCurrentSubMatrNos, SizeCurrentSubMatr);
	}

	void DeleteArrays_AutoRelax(int* TotArrSubMatrNos, int* SubMatrLengths)
	{
		if(TotArrSubMatrNos != NULL)
		{
			delete[] TotArrSubMatrNos; TotArrSubMatrNos = NULL;
		}
		if(SubMatrLengths != NULL)
		{
			delete[] SubMatrLengths; SubMatrLengths = NULL;
		}
		if(mArrAuxQuasiExtField != NULL)
		{
			delete[] mArrAuxQuasiExtField; mArrAuxQuasiExtField = NULL;
		}
	}

	void DeleteInterMatrData()
	{
		if(IntrctPtr != NULL) { delete IntrctPtr; IntrctPtr = NULL;}
	}

};

//-------------------------------------------------------------------------

#endif
