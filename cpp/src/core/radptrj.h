/*-------------------------------------------------------------------------
*
* File name:      radptrj.h
*
* Project:        RADIA
*
* Description:    Charged particle trajectory / dynamics in 3D magnetic field
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADPTRJ_H
#define __RADPTRJ_H

#include "radg3d.h"
#include <math.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTMagHarm {

	int HarmNo;
	char XorZ; //'x' or 'z'
	double B;
	double Phase;

	radTMagHarm(int In_HarmNo, char In_XorZ, double In_B, double In_Phase)
	{
		HarmNo = In_HarmNo; XorZ = In_XorZ; B = In_B; Phase = In_Phase;
	}
	radTMagHarm() 
	{
		HarmNo = 0; XorZ = 0; B = 0; Phase = 0;
	}
};


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTPrtclTrj {
	double Energy; // In GeV
	double ChargeToMomentum; // In SI units

	radTg3d* FldSrcPtr;
	radTField Field;
	TVector3d ZeroVect;

	double *dym_rk4, *dyt_rk4, *yt_rk4;
	double *dysav_rks5, *ysav_rks5, *ytemp_rks5;
	double *y_ap, *dydx_ap;
	int kmax_ap, count_ap;
	double *xp_ap, **yp_ap, dxsav_ap; 

	double* Y;
	double* dYdx;

	radTSend Send;
	int AmOfEq;

	short OnPrc;
	double *PrecArray, EpsTol;
	int MaxAutoStp;

public:
	radTPrtclTrj(double, radTg3d*, const radTCompCriterium&, short, double* =NULL, double =1., int =5000);
	radTPrtclTrj(radTg3d* InFldSourcePtr, const radTCompCriterium& InCompCrit, double inEnergyGeV =0) // For Focusing Potential
	{
		FldSrcPtr = InFldSourcePtr;
		
		Field.CompCriterium = InCompCrit; //? OC050413
		Energy = inEnergyGeV;

		dym_rk4 = dyt_rk4 = yt_rk4 = Y = dYdx = NULL;
		OnPrc = 0;
	}
	~radTPrtclTrj();

	void RungeKutta4(double*, double*, double, double, double*);
	int RungeKuttaStep5(double*, double*, double*, double, double*, double*);
	int AutoPropagate(double*, double, double, double, double, int*, int*);
	void Tabulation(double*, double, double, int, double*);

	void TrjEqs(double s, double* F, double* dFds)
	{
		double xd=F[1], zd=F[3];
		double xdxd=xd*xd, zdzd=zd*zd, xdzd=xd*zd;
		double Buf = ChargeToMomentum * sqrt(1.+xdxd+zdzd);

		Field.P = TVector3d(F[0], s, F[2]);
		Field.B = ZeroVect;
		FldSrcPtr->B_genComp(&Field); TVector3d& rB = Field.B;

		dFds[0] = xd;
		dFds[1] = -Buf*(zd*rB.y - (1.+xdxd)*rB.z + xdzd*rB.x);
		dFds[2] = zd;
		dFds[3] = Buf*(xd*rB.y - (1.+zdzd)*rB.x + xdzd*rB.z);
	}
	void DefOrtComponent(const TVector3d& vUn, const TVector3d& A, TVector3d& OrtCompA)
	{
		OrtCompA = A - ((vUn*A)*vUn);
	}

	double FocusingPotential(const TVector3d&, const TVector3d&, int);

	//void ComputeSecondOrderKickPer(const TVector3d& P1Vect, const TVector3d& NlongVect, double per, int nper, const TVector3d& N1Vect, double r1, int np1, double r2, int np2, int nharm, int ns, double d1, double d2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pBtE2Int);
	void ComputeSecondOrderKickPer(const TVector3d& P1Vect, const TVector3d& NlongVect, double per, double nper, const TVector3d& N1Vect, double r1, int np1, double r2, int np2, int nharm, int ns, double d1, double d2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pBtE2Int, char cKickUnit); //050413
    //void ComputeFocusPotentPerOnMesh2D(const TVector3d& P1Vect, const TVector3d& NlongVect, double per, int nper, const TVector3d& N1Vect, int np1, double step1, double d1, int np2, double step2, double d2, int nharm, int ns, double* pFocPotData, double* pBtE2Int);
    void ComputeFocusPotentPerOnMesh2D(const TVector3d& P1Vect, const TVector3d& NlongVect, double per, double nper, const TVector3d& N1Vect, int np1, double step1, double d1, int np2, double step2, double d2, int nharm, int ns, double* pFocPotData, double* pBtE2Int);
    //void ComputeFocusPotentPerInOneOrFourPoints(const TVector3d& P0, const TVector3d& dsNlongU, const TVector3d& d1N1U, const TVector3d& d2N2U, const TVector3d& N1U, const TVector3d& N2U, char UseSpecStepForDeriv1, char UseSpecStepForDeriv2, double per, int nper, int nharm, int ns, double* pAuxBuffer, double* pFocPotData);
    void ComputeFocusPotentPerInOneOrFourPoints(const TVector3d& P0, const TVector3d& dsNlongU, const TVector3d& d1N1U, const TVector3d& d2N2U, const TVector3d& N1U, const TVector3d& N2U, char UseSpecStepForDeriv1, char UseSpecStepForDeriv2, double per, double nper, int nharm, int ns, double* pAuxBuffer, double* pFocPotData);
    //double FocusPotentPerInOnePoint(const TVector3d& P0, const TVector3d& dsNlongU, const TVector3d& N1U, const TVector3d& N2U, double per, int nper, int nharm, int ns, double* pAuxBuffer);
    double FocusPotentPerInOnePoint(const TVector3d& P0, const TVector3d& dsNlongU, const TVector3d& N1U, const TVector3d& N2U, double per, double nper, int nharm, int ns, double* pAuxBuffer);
	//void ComputeFocusPotentDerivOnMesh2D(double* pFocPotData, int np1, double d1, char UseSpecStepForDeriv1, int np2, double d2, char UseSpecStepForDeriv2, double* pKickData1, double* pKickData2);
	void ComputeFocusPotentDerivOnMesh2D(double* pFocPotData, int np1, double d1, char UseSpecStepForDeriv1, int np2, double d2, char UseSpecStepForDeriv2, double* pKickData1, double* pKickData2, char cKickUnit); //050413
	
	void ComputeSecondOrderKick(const TVector3d& P1Vect, const TVector3d& NlongVect, double* ArrLongDist, int lenArrLongDist, int ns, const TVector3d& N1Vect, double r1, int np1, double r2, int np2, double d1, double d2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pBtE2Int);
	void ComputeFocusPotentOnMesh3D(const TVector3d& P1Vect, const TVector3d& NlongVect, double* ArrLongDist, int AmOfLongGridPts, int ns, const TVector3d& N1Vect, int np1, double step1, double d1, int np2, double step2, double d2, char UseSpecStepForDeriv1, char UseSpecStepForDeriv2, double* pFocPotData, double* pIntBtrE2Data);
	void SetupTransverseMeshForKickCalc(double r1, int np1, double r2, int np2, double& d1, double& d2, double& step1, double& step2, char& UseSpecStepForDeriv1, char& UseSpecStepForDeriv2, double* pCoordDir1, double* pCoordDir2, long& AmOfTransvPtsFocPot);
	void SetupFourOrOneTansvInitPt(const TVector3d& P0, const TVector3d& d1N1U, const TVector3d& d2N2U, char UseSpecStepForDeriv1, char UseSpecStepForDeriv2, TVector3d* InitTransvP)
	{
		if(UseSpecStepForDeriv1)
		{
			InitTransvP[0] = P0 - d1N1U;
			InitTransvP[1] = P0 + d1N1U;
			if(UseSpecStepForDeriv2)
			{
				InitTransvP[2] = P0 - d2N2U;
				InitTransvP[3] = P0 + d2N2U;
			}
		}
		else
		{
			if(UseSpecStepForDeriv2)
			{
				InitTransvP[0] = P0 - d2N2U;
				InitTransvP[1] = P0 + d2N2U;
			}
			else
			{
				InitTransvP[0] = P0;
			}
		}
	}
    void SetupLongNumPointsArray(double* ArrLongDist, int AmOfLongGridPts, int ns, int* ArrLongNumPts);
	void ComputeFocusPotentArrayForOneTransvPoint(const TVector3d& TrInitP, const TVector3d& NlongU, double* ArrLongDist, int AmOfLongGridPts, int* ArrLongNumPts, const TVector3d& N1U, const TVector3d& N2U, double* tFocPotData, double* pIntBtrE2);
	void ComputeFocusPotentDerivOnMesh3D(double* pFocPotData, int AmOfLongGridPts, int np1, double d1, char UseSpecStepForDeriv1, int np2, double d2, char UseSpecStepForDeriv2, double* pKickData1, double* pKickData2);
	void SumUpPartialSecondOrderKicks(long AmOfLongGridPts, long AmOfTransvPts, double* pKickData1, double* pKickData2, double* pIBe2, double* pAuxTotKickData1, double* pAuxTotKickData2, double* pAuxTotIBe2);
	void ArrangeAllSecondOrderKickData(long AmOfLongGridPts, long AmOfTransvPts, double* pTotKickData1, double* pTotKickData2, double* pTotIBe2, double* pKickData1, double* pKickData2, double* pIBe2);

    void AnalyzeForHarmonics(double* pB, int AmOfPts, double Per, double RelPrec, char XorZ, int& AmOfHarm, radTMagHarm*& MagHarmArr);
	void AnalyzeForHarmonics_DeleteAuxArrays(float*& AuxDataContIn, float*& AuxDataContOut, double*& CkArr, double*& PhikArr, int*& HarmNoArr)
	{
		if(AuxDataContIn != 0) { delete[] AuxDataContIn; AuxDataContIn = 0;}
		if(AuxDataContOut != 0) { delete[] AuxDataContOut; AuxDataContOut = 0;}
		if(CkArr != 0) { delete[] CkArr; CkArr = 0;}
		if(PhikArr != 0) { delete[] PhikArr; PhikArr = 0;}
		if(HarmNoArr != 0) { delete[] HarmNoArr; HarmNoArr = 0;}
	}


	static void ComposeStrReportSecondOrderKickPer(const char* StrComment, double len, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pBtE2Int, char*& StrReport, char cKickUnit=1, char cOutFormat=1);
//#ifdef __GCC__
//	static void ComposeStrReportSecondOrderKickPer_AddMainData(ostrstream& OutStream, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData);
//#else
	//static void ComposeStrReportSecondOrderKickPer_AddMainData(ostringstream& OutStream, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData);
	static void ComposeStrReportSecondOrderKickPer_AddMainData(ostringstream& OutStream, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData, char cOutFormat=1);
//#endif
	static void ComposeStrReportSecondOrderKick(const char* StrComment, double* ArrLongDist, int lenArrLongDist, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pIBe2, char*& pStrReport);

	static double ApplyZeroTol(double v, double tol) { return (::fabs(v) <= tol)? 0. : v;}
};

//-------------------------------------------------------------------------

#endif
