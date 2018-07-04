
#ifndef __IGINTERF_H
#define __IGINTERF_H

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
//#include "XOPStructureAlignmentReset.h" //OC080613 (moved to XOPSupport 6)

#ifndef __IGINTERFSTR_H
#include "iginterfstr.h"
#endif
#ifndef __IGERCODE_H
#include "igercode.h"
#endif

//#ifndef __STLSTART_H
//#include "stlstart.h"
//#endif
#include <vector>
using namespace std;

//*************************************************************************

class srTIgorInterf {

	static char gProgressIndicatorWinName[256];
	static long gProgressIndicatorScale;
	static long gRadIndices[MAX_DIMENSIONS+1];

public:

	static void Message(const char* Text);
	static void WarningMessage(const char* WarnText);
	static void ErrorMessage(const char* ErrText);

	static int GetStringFromTextWave1D(waveHndl wavH, int IndNo, string& sValue);
	static int GetDoubleFromTextWave1D(waveHndl wavH, int IndNo, double& Value);
	static int UpdateTextWave1D(waveHndl wavH, vector<string>* pVectStr);
	static int UpdateTextWaveMD(waveHndl wavH, char** arStrFlat, int* arNp, int numDim);
	static int UpdateTextWaveMD_NestedFor(waveHndl wavH, char** arStrFlat, int* arNp, int numDim);

	static int GetArrDoubleFromNumWave1D(waveHndl wavH, long MaxNp, double*& pData, long& Np);
	static int GetArrDoubleFromCmplxWave1D(waveHndl wavH, long MaxNp, double*& pData, long& Np);
    static int GetArrIntFromNumWave1D(waveHndl wavH, long MaxNp, int*& pData, long& Np);
	static int GetArrDoubleScaledAndNpFromNumWave1D(waveHndl wavH, long MaxNp, double*& pData, long& Np, double& ArgStart, double& ArgStep);
	static int GetArrDoubleAndNpFromNumWave1D(waveHndl wavH, int DeclarNp, int NumCoordPerPoint, double*& pData, int& Np);
 	static int Get4ElemArrDoubleFromNumWave1D(waveHndl wavH, double* pData);
	static int Get3ElemArrDoubleFromNumWave1D(waveHndl wavH, double* pData);
    static int Get2ElemArrDoubleFromNumWave1D(waveHndl wavH, double* pData);
    static int GetDataPtrFromWaveDouble(waveHndl wavH, double*& pData, int& hState);
    static int GetDataPtrFromWaveDoubleOrFloat(waveHndl wavH, double*& pdData, float*& pfData, int& hState);
	static int GetDataPtrFromWaveDoubleOrFloat1D(waveHndl wavH, double*& pdData, float*& pfData, long& LenData, int& hState);
    static int ExtractCStrFromHandle(Handle hStr, int MaxLenStr, char* Str);
	static int GetArrDoubleFromNumWaveN3(waveHndl wavH, double*& pData, int& Np);
	static int GetArrDoubleFromNumWave2D(waveHndl wavH, int& NumRows, int& NumCols, double*& pData);
	//static int GetArrIntFromNumWave2D(waveHndl wavH, int& NumRows, int& NumCols, int*& pData);
	static int GetArrDoubleFromNumWave3D(waveHndl wavH, int& NumRows, int& NumCols, int& NumChanks, double*& pData);
	static int GetArrDoubleScaledFromNumWave2D(waveHndl wavH, int& NumRows, int& NumCols, double*& pData, double* ArgStart, double* ArgStep);
	static int GetArrDoubleScaledFromNumWave3D(waveHndl wavH, int& NumRows, int& NumCols, int& NumChanks, double*& pData, double* ArgStart, double* ArgStep);
	static int GetArrStringFromTextWave1D(waveHndl wavH, int& NumRows, char**& pStrData);
	static int GetArrStringFromTextWave2D(waveHndl wavH, int& NumRows, int& NumCols, char**& pStrData);
	//static int GetArrDoubleScaledFromNumWave4D(waveHndl wavH, double*& pData, int* arNumPt, double* arArgStart, double* arArgStep, char** arStrUnits);
	static int GetArrFloatScaledFromNumWave4D(waveHndl wavH, float*& pData, int* arNumPt, double* arArgStart, double* arArgStep, char* arStrUnits[], bool* pDataIsNotZero =0);
	static int GetTotNumPointsInWave4D(waveHndl wavH, long& totNumPt);

	static int GetSizeNumWave1D(waveHndl wavH, long& size);
	static int GetDimSizesWaveMD(waveHndl wavH, long* arDimSizes, long& numDims);

	static int GetWaveData(waveHndl wavH, TDataWaveMD* pWaveData, int& hState);

	static int ReDimNumWave1D(waveHndl wavH, int NumRows);
	static int ReDimNumWave2D(waveHndl wavH, int NumRows, int NumCols);
	static int ReDimNumWave(waveHndl wavH, int* dimSizes, int numDims, int newDataType);
	static int ReDimNumWave(waveHndl wavH, long* dimSizes, long numDims, int newDataType);
	static int ReDimWave3D(waveHndl wavH, int NumRows, int NumCols, int numChunks);

	static int SetDataInNumWave(waveHndl wavH, double* pTempCalcData, long Nb, char* DataUnitStr);
	static int SetDataInNumWave(waveHndl wavH, float* pTempCalcData, long Nb, char* DataUnitStr);
	static int SetDataInNumWave(waveHndl wavH, int* pData, long Np, char* DataUnitStr);

	static int SetScaleNumWave(waveHndl wavH, int nDim, double* ArgStarts, double* ArgSteps, char** ArgValUnits);
	static int SetDataInTextWave2D(waveHndl wavH, char** arrStr, int* arrLen, int numArr);
	static int SetDataInTextWave3D(waveHndl wavH, char** arrStr, int numRows, int numCols, int numChunks);
	static int SetZeroTextWave2D(waveHndl wavH);
	static int SetZeroTextWave3D(waveHndl wavH);

    static int ReleaseWave(waveHndl wavH, int hState);

	static int ShowCompProgress(double CurNormVal);
    static int InitCompProgressIndicator();
    static int UpdateCompProgressIndicator(double CurNormVal);
    static void DestroyCompProgressIndicator();

    static int SetNumWavePointValue1D(waveHndl wavH, int ind, DOUBLE val);
    static int UpdateScalingAndFinishWorkingWithWave(TDataWaveMD* pWaveData, waveHndl wavH, int hState);

	template<class T> static int GetArrFromNumWave2D(waveHndl wavH, int& NumRows, int& NumCols, T*& pData)
	{
		if(wavH == NIL) return NOWAV;
		int waveType = WaveType(wavH);
		if((waveType != NT_FP64) && (waveType != NT_FP32) && (waveType != NT_I8) && (waveType != NT_I16) && (waveType != NT_I32)) return NUMERIC_WAVE_REQUIRED;

		//long numDimensions;
		int numDimensions; //OC080613 (moving to XOPSupport 6)
		long dimensionSizes[MAX_DIMENSIONS+1];
		int result;

		if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;

		long nRows = dimensionSizes[0];
		long nCols = dimensionSizes[1];
		if(nCols <= 0) nCols = 1;

		if(NumRows > 0)
		{
	      if(nRows < NumRows) return TOO_SHORT_WAVE;
		}
		else NumRows = nRows;

		if(NumCols > 0)
		{
	      if(nCols < NumCols) return TOO_SHORT_WAVE;
		}
		else NumCols = nCols;

		long Np = NumRows*NumCols;
		if(Np <= 0) return 0; //allow empty waves
	
		if(pData == 0) pData = new T[Np];
	    T *tData = pData;

		long dataOffset;
		if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
		//int hState = MoveLockHandle(wavH); //OC080613 (moving to XOPSupport 6)
		char* dataStartPtr = (char*)(*wavH) + dataOffset;

		if((waveType == NT_I8) || (waveType == NT_I16) || (waveType == NT_I32)) //to check / correct / make separate type conversion
		{
			int* fp = (int*)dataStartPtr;
			for(long i=0; i<Np; i++) *(tData++) = (T)(*(fp++));
		}
		else if(waveType == NT_FP32)
		{
			float* fp = (float*)dataStartPtr;
			for(long i=0; i<Np; i++) *(tData++) = (T)(*(fp++));
		}
		else if(waveType == NT_FP64)
		{
			DOUBLE* dp = (DOUBLE*)dataStartPtr;
			for(long i=0; i<Np; i++) *(tData++) = (T)(*(dp++));
		}
		//HSetState((Handle)wavH, hState); //OC080613 (moving to XOPSupport 6)
		return 0;
	}

	//static int GetElecBeamThin(waveHndl wavH, double& I, double* Mom1, double& s);
	//static int SetElecBeamThin(waveHndl wavH, double I, double* Mom1, double s);
	//static int GetElecBeamThick(waveHndl wavH, double& I, double* pMom1, double* pMom2, int& nMom2, double& s0, int& TypeDistrTransverse, int& TypeDistrLongitudinal, double& ShortNoiseFactor);
    //static int SetElecBeamThick(waveHndl wavH, double I, double* pMom1, double* pMom2, double s0, int TypeDistrTransverse, int TypeDistrLongitudinal, double ShortNoiseFactor);
	//static int GetElecBeamTwiss(waveHndl wavH, double* pHorTwiss, double* pVertTwiss);
	//static int SetElecBeamEmitAndTwiss(waveHndl wavH, double HorEmit, double VertEmit, double* pHorTwiss, double* pVertTwiss);
	//static int GetMagFieldTransvUnif(waveHndl wavH, double FieldAbsZeroTol, double& sStart, double& sStep, long& np, double*& pBH, bool& BHIsZero, double*& pBV, bool& BVIsZero);
	//static int GetMagFieldPeriodic(waveHndl wavH, double& PerLength, double& TotLength, int& AmOfHarm, srTMagFldHarm*& MagFldHarmArr, int& TypeOfUnd, double& TaperPar_TU, double& PhaseSh_OK, int& FldErrTypeSASE, double& FldErrRMS, double& NatFocNxSASE, double& NatFocNySASE, int& TaperTypeSASE, double& TaperStartSASE, double& TaperRelFldChgSASE);
	//static int GetMagFieldPeriodic(waveHndl wavH, double& PerLength, double& TotLength, int& AmOfHarm, int*& ArrHarmNo, char*& ArrXorZ, double*& ArrK, double*& ArrPhase, int& TypeOfUnd, double& TaperPar_TU, double& PhaseSh_OK, int& FldErrTypeSASE, double& FldErrRMS, double& NatFocNxSASE, double& NatFocNySASE, int& TaperTypeSASE, double& TaperStartSASE, double& TaperRelFldChgSASE);
	//static int SetMagFieldPeriodic(waveHndl wavH, double PerLength, double TotLength, int AmOfHarm, int* ArrHarmNo, char* ArrXorZ, double* ArrK, double* ArrPhase, int TypeOfUnd, double TaperPar_TU, double PhaseSh_OK, double Fund_keV_per_GeV2);
	//static int GetMagFieldConstant(waveHndl wavH, double& Bx, double& Bz);
	//static int IdentifyMagFieldType(waveHndl wavH, int& MagFieldType);
	//static int IdentifyMagFieldTypeFromName(char* MagFieldName);
	//static int GetTrjDataPointers(waveHndl wOutHorAng, waveHndl wOutHorCoor, waveHndl wOutVerAng, waveHndl wOutVerCoor, double*& pOutBtxData, double*& pOutXData, double*& pOutBtzData, double*& pOutZData, int& hStateOutBtxData, int& hStateOutXData, int& hStateOutBtzData, int& hStateOutZData);
	//static int GetPrecParamWfrSamplingForPropag(waveHndl wavH, bool& AllowAutoChoiceOfNxNzForPropag, double& NxNzOversamplingParam);
	//static int GetWfrSampling(waveHndl wavH, double& s, double& zSt, double& zFi, int& nz, double& xSt, double& xFi, int& nx, double& eSt, double& eFi, int& ne, char* PhotEnUnits);
	//static int GetPrecParamElectricFieldComp(waveHndl wavH, int& IntegMeth, double& RelPrecOrStep, double& sStart, double& sEnd);
	//static int GetPrecParamStokesPerComp(waveHndl wavH, int& InitHarm, int& FinHarm, double& Kns, double& Knphi, char& IntensityOrFlux);
	//static int GetPrecParamPowDensComp(waveHndl wavH, double& PrecFact, int& Method);
	//static int GetPrecParamMagArb2Per(waveHndl wavH, double& RelPrec, int& MaxAmOfHarm);
	//static int GetPrecParamStokesArbComp(waveHndl wavH, srTParPrecStokesArb* pPrecStokesArb);
	//static int GetSRWRadInData(waveHndl wavH, srTSRWRadInData* pSRWRadInData);
	//static int GetSRWStokesInData(waveHndl wavH, srTSRWStokesInData* pStokesAccessData);
	//static int GetSRWPowDensInData(waveHndl wavH, srTSRWPowDensInData* pSRWPowDensInData);
	//static int GetIntensExtractParam(waveHndl wExtractParam, int& PolarizCompon, int& Int_or_Phase, int& PlotType, int& TransvPres, double& ePh, double& x, double& z);
 //   static int GetIntensExtractData(waveHndl wExtractedData, int& hStateExtractedData, char*& pExtractedData);
	//static int GetVectorOfStrings(waveHndl wavH, vector<string>*);
 //   static int GetVectorOfStrings(const char* StructName, vector<string>*);
	//static int GetGaussianBeam(waveHndl& wavH, double& Phot, int& Polar, double& SigmaX, double& SigmaZ, int& mx, int& mz, double* Mom1Arr, double& s0, double& SigmaT, int& TypeDistrT);

 //   static int GetAndSetElecBeamThick(int* piElecBeam, waveHndl& wavH);
	//static int GetAndSetElecBeamThin(int* piElecBeam, waveHndl& wavH);
	//static int GetAndSetGaussianBeam(int* piGsnBeam, waveHndl& wavH);
 //   static int GetAndSetMagFieldGen(int* piMagFld, waveHndl& wavH);
	//static int GetAndSetWfrSampling(int* piWfrSmp, waveHndl& wavH);

	//static int WfrModify(int ActionNo, srTSRWRadInData* pNewRadStructAccessData, char PolarizComp);
	//static int WfrCreateNew(srTSRWRadInData* pNewRadStructAccessData);
	//static int WfrModifyNeNxNz(srTSRWRadInData* pNewRadStructAccessData, char PolarizComp);
	//static int WfrDelete(srTSRWRadInData* pNewRadStructAccessData);
	//static int WfrRename(srTSRWRadInData* pNewRadStructAccessData);
	//static int WfrGetNames(srTSRWRadInData* pNewRadStructAccessData);

	//static int SetupExtractedWaveData(srTSRWRadInData* pSRWRadInData, int Int_or_Phase, int PlotType, int TransvPres, char* pExtrData, waveHndl wavH, int hStateExtractedData, srTIgorWaveAccessData* pExtrWaveData);
	//static int FinishWorkingWithSRWRadStruct(srTSRWRadInData* pSRWRadInData);
	//static int FinishWorkingWithSRWStokesStruct(srTSRWStokesInData* pStokesAccessData);
	//static int FinishWorkingWithSRWPowDensStruct(srTSRWPowDensInData* pPowDensAccessData);
	//static int FinishWorkingWithWave(srTIgorWaveAccessData* pExtrWaveData);
	//static int FinishWorkingWithTrjDataPointers(waveHndl wOutHorAng, waveHndl wOutHorCoor, waveHndl wOutVerAng, waveHndl wOutVerCoor, int hStateHorAng, int hStateHorCoor, int hStateVerAng, int hStateVerCoor);

private:

	//static int UpdateNumberPositionInSRWRad(srTSRWRadInData* pSRWRadStructAccessData, long* RadIndices, char* CharBuf, int AmOfBytes);
	//static int UpdateSrwWfrAuxData(srTSRWRadInData* pSRWRadInData);
	//static int SetupSrwWfrAuxData(srTSRWRadInData* pSRWRadInData);
	//static int UpdateTextPositionInSRWRad(srTSRWRadInData* pRad, int Index, char* Text);

};

//*************************************************************************

struct srTIgorAccessNumFunc2Waves {
	char FuncName[200];
	char ResVarName[200];
	int NumVars;
	waveHndl wArg1; //first argument
	waveHndl wArg2; //second argument
	char ArgWaveName1[100];
	char ArgWaveName2[100];
	char CommandStr[2048];

	srTIgorAccessNumFunc2Waves() { wArg1 = NIL; wArg2 = NIL; *CommandStr = '\0'; NumVars = 0;};
	~srTIgorAccessNumFunc2Waves() {};

	int Initialize(char* InFuncName, int NumVar);
	double CalcFunc(double* pArgs, double* pArgsMin, double FMin);
	//template<class T> double CalcFunc(T* pArgs, T* pArgsMin, double FMin);
};

//*************************************************************************

#endif