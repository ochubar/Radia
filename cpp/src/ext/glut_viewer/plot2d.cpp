
#ifndef __PLOT2D_H
#include "plot2d.h"
#endif
//#ifndef __GMMETH_H
//#include "gmmeth.h"
//#endif

//#include <math.h>

//---------------------------------------------------------------------------------------------------

CWinContPlot2D::CWinContPlot2D(double** FuncValues, double** ArgValues, double* ArgStart, double* ArgStep, long* Size, double** CurveOptions, int CurveNumber, char** Units, char** Labels, double* GraphOptions, char* WinTitle)
{
    if((FuncValues == 0) || (Size == 0) || (CurveNumber <= 0)) throw 1;

	bool* ArrArgStartStepDefined = new bool[CurveNumber];
	bool* ArrArgTableDefined = new bool[CurveNumber];
    if(!CheckHowArgIsDefined(ArgValues, ArgStart, ArgStep, CurveNumber, ArrArgStartStepDefined, ArrArgTableDefined))
	{
        delete[] ArrArgStartStepDefined;
        delete[] ArrArgTableDefined;
		throw 2;
	}

	

/**
	double*	m_LnVertCoord;
	int	m_LnNv;
	int* m_LnVertInd;
	int* m_LnLen;
	float* m_LnColors;
	int	m_LnN;
	bool m_LnWereDefined; // = false;
**/

    if(WinTitle != 0) 
	{
		int LenWinTitle = (int)strlen(WinTitle);
		if(LenWinTitle > 500) LenWinTitle = 500;
		strncpy(m_WinTitle, WinTitle, LenWinTitle);
		m_WinTitle[LenWinTitle] = '\0';
	}
	else strcpy(m_WinTitle, CSimpleGraph::sm_DefWinTitle);

	if(ArrArgStartStepDefined != 0) delete[] ArrArgStartStepDefined;
	if(ArrArgTableDefined != 0) delete[] ArrArgTableDefined;
}

//---------------------------------------------------------------------------------------------------

bool CWinContPlot2D::CheckHowArgIsDefined(double** ArgValues, double* ArgStart, double* ArgStep, int CurveNumber, bool* ArrArgStartStepDefined, bool* ArrArgTableDefined)
{
	if((ArrArgStartStepDefined == 0) || (ArrArgTableDefined == 0) || (CurveNumber <= 0)) return false;

	for(int i=0; i<CurveNumber; i++)
	{
        ArrArgStartStepDefined[i] = false;
		ArrArgTableDefined[i] = false;

		if(ArgValues != 0)
		{
			if(ArgValues[i] != 0) ArrArgTableDefined[i] = true;
		}
		if(!ArrArgTableDefined[i])
		{
			if((ArgStart != 0) && (ArgStep != 0))
			{
				ArrArgStartStepDefined[i] = true;
			}
		}
		if((!ArrArgStartStepDefined[i]) && (!ArrArgTableDefined[i])) return false;
	}
	return true;
}

//---------------------------------------------------------------------------------------------------

void CWinContPlot2D::FindCurveAbsMinAndMax(double** FuncValues, double** ArgValues, double* ArgStart, double* ArgStep, long* Size, bool* ArrArgTableDefined, bool* ArrArgStartStepDefined, int AmOfCurves, double& MinArg, double& MaxArg, double& MinVal, double& MaxVal)
{
    MinArg = MaxArg = MinVal = MaxVal = 0;
	if((FuncValues == 0) || (Size == 0) || (AmOfCurves <= 0) || (ArrArgStartStepDefined == 0) || (ArrArgTableDefined == 0)) return;

	MinArg = MinVal = 1.E+53;
	MaxArg = MaxVal = -1.E+53;

	for(int i=0; i<AmOfCurves; i++)
	{
		double* CurArgValues = 0;
		double CurArgStart = 0, CurArgStep = 0;
		if(ArrArgTableDefined[i]) CurArgValues = ArgValues[i];
		else if(ArrArgStartStepDefined[i]) { CurArgStart = ArgStart[i]; CurArgStep = ArgStep[i];}

		double LocMinArg = 1.E+53, LocMaxArg = -1.E+53, LocMinVal = 1.E+53, LocMaxVal = -1.E+53;
		FindOneCurveMinMax(FuncValues[i], CurArgValues, CurArgStart, CurArgStep, Size[i], LocMinArg, LocMaxArg, LocMinVal, LocMaxVal);

		//dddddddddddddddddddddddddddddddddddddd
		//if(mCurveAbsMinVal > tCurvesExt->mAbsMinVal) mCurveAbsMinVal = tCurvesExt->mAbsMinVal;
		//if(mCurveAbsMaxVal < tCurvesExt->mAbsMaxVal) mCurveAbsMaxVal = tCurvesExt->mAbsMaxVal;
		//if(mCurveAbsMinArg > tCurvesExt->mAbsMinArg) mCurveAbsMinArg = tCurvesExt->mAbsMinArg;
		//if(mCurveAbsMaxArg < tCurvesExt->mAbsMaxArg) mCurveAbsMaxArg = tCurvesExt->mAbsMaxArg;
		//tCurvesExt++;
	}
/*

	//Shifting limits to include 0, if they are close to it
	double RelZeroTol = 0.01;
	if(mCurveAbsMaxArg != mCurveAbsMinArg)
	{
		double ArgRange = mCurveAbsMaxArg - mCurveAbsMinArg;
		if(mCurveAbsMinArg > 0) 
		{
			if(mCurveAbsMinArg/ArgRange < RelZeroTol) mCurveAbsMinArg = 0;
		}
		if(mCurveAbsMaxArg < 0)
		{
			if(-mCurveAbsMaxArg/ArgRange < RelZeroTol) mCurveAbsMaxArg = 0;
		}
	}
	if(mCurveAbsMaxVal != mCurveAbsMinVal)
	{
		double ValRange = mCurveAbsMaxVal - mCurveAbsMinVal;
		if(mCurveAbsMinVal > 0) 
		{
			if(mCurveAbsMinVal/ValRange < RelZeroTol) mCurveAbsMinVal = 0;
		}
		if(mCurveAbsMaxVal < 0)
		{
			if(-mCurveAbsMaxVal/ValRange < RelZeroTol) mCurveAbsMaxVal = 0;
		}
	}
*/
}

//---------------------------------------------------------------------------------------------------

void CWinContPlot2D::FindOneCurveMinMax(double* FuncValues, double* ArgValues, double ArgStart, double ArgStep, long Size, double& MinArg, double& MaxArg, double& MinVal, double& MaxVal)
{
    MinArg = MaxArg = MinVal = MaxVal = 0;
	if((FuncValues == 0) || (Size == 0)) return;

	bool ArgTableDefined = (ArgValues != 0);
	if(!ArgTableDefined) 
	{
		MinArg = ArgStart; MaxArg = MinArg + (Size - 1)*ArgStep;
		if(MaxArg < MinArg)
		{
			double Buf = MaxArg;
            MaxArg = MinArg; MinArg = Buf;
		}
	}

	MinArg = 1.E+53; MaxArg = -1.E+53; MinVal = 1.E+53; MaxVal = -1.E+53;
	for(long i=0; i<Size; i++)
	{
		double CurFuncVal = FuncValues[i];
		if(MinVal > CurFuncVal) MinVal = CurFuncVal;
		if(MaxVal < CurFuncVal) MaxVal = CurFuncVal;
		if(ArgTableDefined)
		{
            double CurArgVal = FuncValues[i];
			if(MinArg > CurArgVal) MinArg = CurArgVal;
			if(MaxArg < CurArgVal) MaxArg = CurArgVal;
		}
	}
}

//---------------------------------------------------------------------------------------------------

