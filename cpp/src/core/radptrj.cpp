/*-------------------------------------------------------------------------
*
* File name:      radptrj.cpp
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

#include "radptrj.h"
#include "gmfft.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTPrtclTrj::radTPrtclTrj(double InEnergy, radTg3d* InFldSrcPtr, const radTCompCriterium& InCompCriterium, 
						   short InOnPrc, double* InPrecArray, double InEpsTol, int InMaxAutoStp) 
{
	Energy = InEnergy;
	ChargeToMomentum = -2.99792458E+05/(1.E+09 * Energy); // Electron has negative charge
	FldSrcPtr = InFldSrcPtr;

	ZeroVect.x = ZeroVect.y = ZeroVect.z = 0.;
	radTFieldKey LocFieldKey; LocFieldKey.B_= 1;
	Field.FieldKey = LocFieldKey;
	Field.CompCriterium = InCompCriterium;

	OnPrc = InOnPrc; AmOfEq = 4;
	
	dym_rk4 = new double[AmOfEq];
	dyt_rk4 = new double[AmOfEq];
	yt_rk4 = new double[AmOfEq];

	Y = new double[AmOfEq];
	dYdx = new double[AmOfEq];

	if(OnPrc)
	{
		PrecArray = new double[AmOfEq];
		for(int i=0; i<AmOfEq; i++) 
		{
			PrecArray[i] = InPrecArray[i];
			if(PrecArray[i]==0.) PrecArray[i] = 1.E+23;
		}
		EpsTol = InEpsTol;
		MaxAutoStp = InMaxAutoStp;

		dysav_rks5 = new double[AmOfEq];
		ysav_rks5 = new double[AmOfEq];
		ytemp_rks5 = new double[AmOfEq];

		y_ap = new double[AmOfEq];
		dydx_ap = new double[AmOfEq];

		kmax_ap = count_ap = 0;                    // To allow storage of intermediate results in AutoPropagate,
		xp_ap = NULL; yp_ap = NULL; dxsav_ap = 0.; // set this to desired values and allocate memory
	}
}

//-------------------------------------------------------------------------

radTPrtclTrj::~radTPrtclTrj()
{
	if(dym_rk4!=NULL) delete[] dym_rk4;
	if(dyt_rk4!=NULL) delete[] dyt_rk4;
	if(yt_rk4!=NULL) delete[] yt_rk4;

	if(Y!=NULL) delete[] Y; 
	if(dYdx!=NULL) delete[] dYdx;

	if(OnPrc)
	{
		delete[] PrecArray;

		delete[] dysav_rks5;
		delete[] ysav_rks5;
		delete[] ytemp_rks5;

		delete[] y_ap;
		delete[] dydx_ap;
	}
}

//-------------------------------------------------------------------------

void radTPrtclTrj::RungeKutta4(double* y, double* dydx, double x, double h, double* yout)
{
	double xh, hh=0.5*h, h6=h/6.;
	xh=x+hh;

	int i;
	for(i=0; i<AmOfEq; i++) yt_rk4[i] = y[i] + hh*dydx[i];
	TrjEqs(xh, yt_rk4, dyt_rk4);
	for(i=0; i<AmOfEq; i++) yt_rk4[i] = y[i] + hh*dyt_rk4[i];
	TrjEqs(xh, yt_rk4, dym_rk4);
	for(i=0; i<AmOfEq; i++)
	{
		yt_rk4[i] = y[i] + h*dym_rk4[i];
		dym_rk4[i] += dyt_rk4[i];
	}
	TrjEqs(x+h, yt_rk4, dyt_rk4);
	for(i=0; i<AmOfEq; i++) yout[i] = y[i] + h6*(dydx[i]+dyt_rk4[i]+2.*dym_rk4[i]);
}

//-------------------------------------------------------------------------

void radTPrtclTrj::Tabulation(double* F0, double Xmin, double Xmax, int Np, double* Solution)
{
	double Step_x = (Xmax-Xmin)/double(Np-1);
	double x = Xmin;

	double MinStepAllowed = 1.E-10 * Step_x;
	int Nok=0, Nbad=0;

	int k;
	for(k=0; k<AmOfEq; k++) Y[k] = F0[k];

	int AmOfEq_p_1 = AmOfEq+1;

	int Np_mi_1 = Np-1;
	for(int i=0; i<Np; i++)
	{
		if(!OnPrc) TrjEqs(x, Y, dYdx);

		int i_AmOfEq_p_1 = i*AmOfEq_p_1, i_AmOfEq_p_1_p_1 = i_AmOfEq_p_1 + 1;
		Solution[i_AmOfEq_p_1] = x;
		for(k=0; k<AmOfEq; k++) Solution[i_AmOfEq_p_1_p_1 + k] = Y[k];

		if(i != Np_mi_1)
		{
			double x1 = x+Step_x;
			if(OnPrc) AutoPropagate(Y, x, x1, Step_x, MinStepAllowed, &Nok, &Nbad);
			else RungeKutta4(Y, dYdx, x, Step_x, Y);
			x = x1;
		}
	}
}

//-------------------------------------------------------------------------

int radTPrtclTrj::RungeKuttaStep5(double* y, double* dydx, double* x, double htry, double* hdid, double* hnext)
{
	int i;
	double xsav= *x, hh, h= htry, temp, errmax;

	const double PGROW = -0.2;
	const double PSHRNK = -0.25;
	const double FCOR = 1./15.;
	const double SAFETY = 0.9;
	const double ERRCON = 6.0E-04;

	for(i=0; i<AmOfEq; i++)
	{
		ysav_rks5[i] = y[i];
		dysav_rks5[i] = dydx[i];
	}

	for(;;)
	{
		hh = 0.5*h;
		RungeKutta4(ysav_rks5, dysav_rks5, xsav, hh, ytemp_rks5);
		*x = xsav+hh;
		TrjEqs(*x, ytemp_rks5, dydx);
		RungeKutta4(ytemp_rks5, dydx, *x, hh, y);
		*x = xsav+h;
		if(*x == xsav) { Send.ErrorMessage("Radia::Error200"); return 0;}
		RungeKutta4(ysav_rks5, dysav_rks5, xsav, h, ytemp_rks5);

		errmax = 0.;
		for(i=0; i<AmOfEq; i++)
		{
			ytemp_rks5[i] = y[i]-ytemp_rks5[i];
			temp = fabs(ytemp_rks5[i]/PrecArray[i]);
			if(errmax < temp) errmax = temp;
		}
		errmax /= EpsTol;
		if(errmax <= 1.)
		{
			*hdid = h;
			*hnext = (errmax > ERRCON ? SAFETY*h*exp(PGROW*log(errmax)) : 4.*h);
			break;
		}
		h *= SAFETY*exp(PSHRNK*log(errmax));
	}
	for(i=0; i<AmOfEq; i++) y[i] += ytemp_rks5[i]*FCOR;
	return 1;
}

//-------------------------------------------------------------------------

int radTPrtclTrj::AutoPropagate(double* ystart, double x1, double x2, double h1, double hmin, int* nok, int* nbad)
{
	int nstp, i;
	double xsav, x=x1, hnext, hdid, h;

	h = (x2>x1) ? fabs(h1) : -fabs(h1);
	*nok = (*nbad) = count_ap = 0;
	for(i=0; i<AmOfEq; i++) y_ap[i] = ystart[i];

	if(kmax_ap>0) xsav = x-2.*dxsav_ap;
	for(nstp=1; nstp<=MaxAutoStp; nstp++)
	{
		TrjEqs(x, y_ap, dydx_ap);
		if(kmax_ap>0)
			if(fabs(x-xsav) > fabs(dxsav_ap))
				if(count_ap < kmax_ap-1)
				{
					xp_ap[++count_ap] = x;
					for(i=0; i<AmOfEq; i++) yp_ap[i][count_ap] = y_ap[i];
					xsav = x;
				}
		if((x+h-x2)*(x+h-x1) > 0.) h = x2-x;
		if(!RungeKuttaStep5(y_ap, dydx_ap, &x, h, &hdid, &hnext)) return 0;
		if(hdid==h) ++(*nok);
		else ++(*nbad);
		if((x-x2)*(x2-x1) >= 0.)
		{
			for(i=0; i<AmOfEq; i++) ystart[i] = y_ap[i];
			if(kmax_ap>0)
			{
				xp_ap[++count_ap] = x;
				for(i=0; i<AmOfEq; i++) yp_ap[i][count_ap] = y_ap[i];
			}
			return 1;
		}
		if(fabs(hnext) <= hmin) { Send.ErrorMessage("Radia::Error200"); return 0;}
		h = hnext;
	}
	Send.ErrorMessage("Radia::Error201"); return 0;
}

//-------------------------------------------------------------------------

double radTPrtclTrj::FocusingPotential(const TVector3d& InStPoi, const TVector3d& InFiPoi, int Np)
{
	TVector3d Vpar = InFiPoi - InStPoi;
	double AbsVpar = sqrt(Vpar.x*Vpar.x + Vpar.y*Vpar.y + Vpar.z*Vpar.z);
	TVector3d VparUn = (1./AbsVpar)*Vpar;
	double tSt = 1./(Np-1);
	double Weight = 0.5*AbsVpar*tSt;
	TVector3d dVpar = tSt*Vpar;

	TVector3d ZeroVect(0.,0.,0.), OrtCompB, PrOrtCompB(0.,0.,0.);

	radTFieldKey FieldKey;
	FieldKey.B_ = 1;
	radTField Field(FieldKey, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect);
	TVector3d& CurP = Field.P;
	TVector3d& B = Field.B;

	CurP = InStPoi; B = ZeroVect;
	FldSrcPtr->B_genComp(&Field);
	DefOrtComponent(VparUn, B, OrtCompB);
	TVector3d SumOrtCompB = OrtCompB;

	double GenSum=0., GenIntgrnd=0., PrGenIntgrnd=0.;

	for(int i=1; i<Np; i++)
	{
		CurP += dVpar; B = ZeroVect;
		FldSrcPtr->B_genComp(&Field);
		DefOrtComponent(VparUn, B, OrtCompB);

		SumOrtCompB += (OrtCompB + PrOrtCompB);
		PrOrtCompB = OrtCompB;

		GenIntgrnd = (SumOrtCompB.x*SumOrtCompB.x + SumOrtCompB.y*SumOrtCompB.y + SumOrtCompB.z*SumOrtCompB.z);

		GenSum += (GenIntgrnd + PrGenIntgrnd);
		PrGenIntgrnd = GenIntgrnd;
	}
	return GenSum*Weight*Weight*Weight;
}

//-------------------------------------------------------------------------

//void radTPrtclTrj::ComputeSecondOrderKickPer(const TVector3d& P1Vect, const TVector3d& NlongVect, double per, int nper, const TVector3d& N1Vect, double r1, int np1, double r2, int np2, int nharm, int ns, double d1, double d2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pBtE2Int)
void radTPrtclTrj::ComputeSecondOrderKickPer(const TVector3d& P1Vect, const TVector3d& NlongVect, double per, double nper, const TVector3d& N1Vect, double r1, int np1, double r2, int np2, int nharm, int ns, double d1, double d2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pBtE2Int, char cKickUnit) //OC050413
{
	double step1 = 0., step2 = 0.;
	char UseSpecStepForDeriv1 = 0, UseSpecStepForDeriv2 = 0;
	long AmOfPtsFocPot = 0;
	SetupTransverseMeshForKickCalc(r1, np1, r2, np2, d1, d2, step1, step2, UseSpecStepForDeriv1, UseSpecStepForDeriv2, pCoordDir1, pCoordDir2, AmOfPtsFocPot);

	double *pFocPotData = new double[AmOfPtsFocPot];
	ComputeFocusPotentPerOnMesh2D(P1Vect, NlongVect, per, nper, N1Vect, np1, step1, d1, np2, step2, d2, nharm, ns, pFocPotData, pBtE2Int);
	//ComputeFocusPotentDerivOnMesh2D(pFocPotData, np1, d1, UseSpecStepForDeriv1, np2, d2, UseSpecStepForDeriv2, pKickData1, pKickData2);
	ComputeFocusPotentDerivOnMesh2D(pFocPotData, np1, d1, UseSpecStepForDeriv1, np2, d2, UseSpecStepForDeriv2, pKickData1, pKickData2, cKickUnit); //OC050413

	if(pFocPotData != 0) delete[] pFocPotData;
}

//-------------------------------------------------------------------------

//void radTPrtclTrj::ComputeFocusPotentPerOnMesh2D(const TVector3d& P1Vect, const TVector3d& NlongVect, double per, int nper, const TVector3d& N1Vect, int np1, double step1, double d1, int np2, double step2, double d2, int nharm, int ns, double* pFocPotData, double* pBtE2Int)
void radTPrtclTrj::ComputeFocusPotentPerOnMesh2D(const TVector3d& P1Vect, const TVector3d& NlongVect, double per, double nper, const TVector3d& N1Vect, int np1, double step1, double d1, int np2, double step2, double d2, int nharm, int ns, double* pFocPotData, double* pBtE2Int)
{
	const double AbsTolStep = 1.0e-10;

	char UseSpecStepForDeriv1 = (::fabs(step1 - d1) > AbsTolStep)? 1 : 0;
	char UseSpecStepForDeriv2 = (::fabs(step2 - d2) > AbsTolStep)? 1 : 0;

	double InvLenNlong = 1./sqrt(NlongVect.x*NlongVect.x + NlongVect.y*NlongVect.y + NlongVect.z*NlongVect.z);
	TVector3d NlongU = InvLenNlong*NlongVect;
	double InvLenN1 = 1./sqrt(N1Vect.x*N1Vect.x + N1Vect.y*N1Vect.y + N1Vect.z*N1Vect.z);
	TVector3d N1U = InvLenN1*N1Vect;
	TVector3d N2U = N1U^NlongU;

	//double StepLong = per;
	//if(ns > 1) StepLong = per/(ns - 1);
    double StepLong = per/ns;

    TVector3d dsNlongU = StepLong*NlongU;

	TVector3d d1N1U, d2N2U;

	int NumFocPotValPerPoint = 1;
	if(UseSpecStepForDeriv1) 
	{ 
		NumFocPotValPerPoint <<= 1;
        d1N1U = d1*N1U;
	}
	if(UseSpecStepForDeriv2) 
	{ 
		NumFocPotValPerPoint <<= 1;
        d2N2U = d2*N2U;
	}

	int twoNs = ns << 1;
	if(twoNs < 4) twoNs = 4;

	double *pAuxBuffer = new double[twoNs];
	double *tFocPotData = pFocPotData;

	double Arg2Min = -0.5*step2*(np2 - 1), Arg1Min = -0.5*step1*(np1 - 1);
	if(np2 == 1) Arg2Min = -0.5*step2;
	if(np1 == 1) Arg1Min = -0.5*step1;

	double *tBtE2Int = pBtE2Int;
	double arg2 = Arg2Min;
	for(int j=0; j<np2; j++)
	{
        TVector3d CurOffsetN2 = arg2*N2U;

        double arg1 = Arg1Min;
		for(int i=0; i<np1; i++)
		{
            TVector3d CurOffsetN1 = arg1*N1U;
            TVector3d CurP = P1Vect + (CurOffsetN2 + CurOffsetN1);

			ComputeFocusPotentPerInOneOrFourPoints(CurP, dsNlongU, d1N1U, d2N2U, N1U, N2U, UseSpecStepForDeriv1, UseSpecStepForDeriv2, per, nper, nharm, ns, pAuxBuffer, tFocPotData);
			
			tFocPotData += NumFocPotValPerPoint;
			arg1 += step1;
			*(tBtE2Int++) = *pAuxBuffer; //OC100509
		}
		arg2 += step2;
	}
	if(pAuxBuffer != 0) delete[] pAuxBuffer;
}

//-------------------------------------------------------------------------

//void radTPrtclTrj::ComputeFocusPotentPerInOneOrFourPoints(const TVector3d& P0, const TVector3d& dsNlongU, const TVector3d& d1N1U, const TVector3d& d2N2U, const TVector3d& N1U, const TVector3d& N2U, char UseSpecStepForDeriv1, char UseSpecStepForDeriv2, double per, int nper, int nharm, int ns, double* pAuxBuffer, double* pFocPotData)
void radTPrtclTrj::ComputeFocusPotentPerInOneOrFourPoints(const TVector3d& P0, const TVector3d& dsNlongU, const TVector3d& d1N1U, const TVector3d& d2N2U, const TVector3d& N1U, const TVector3d& N2U, char UseSpecStepForDeriv1, char UseSpecStepForDeriv2, double per, double nper, int nharm, int ns, double* pAuxBuffer, double* pFocPotData)
{//at return: pFocPotData contains Focusing Potential data; *pAuxBuffer contains integral of squared transverse magnetic field (always one point only)
	TVector3d InitTransvP[4];
	int Np = 1;

	if(UseSpecStepForDeriv1)
	{
        InitTransvP[0] = P0 - d1N1U;
        InitTransvP[1] = P0 + d1N1U;
		Np = 2;

		if(UseSpecStepForDeriv2)
		{
            InitTransvP[2] = P0 - d2N2U;
            InitTransvP[3] = P0 + d2N2U;
			Np = 4;
		}
	}
	else
	{
		if(UseSpecStepForDeriv2)
		{
            InitTransvP[0] = P0 - d2N2U;
            InitTransvP[1] = P0 + d2N2U;
			Np = 2;
		}
		else
		{
			InitTransvP[0] = P0;
			Np = 1;
		}
	}

	double *tFocPotData = pFocPotData, sumValsIntBe2 = 0;
	for(int i=0; i<Np; i++)
	{
		*(tFocPotData++) = FocusPotentPerInOnePoint(InitTransvP[i], dsNlongU, N1U, N2U, per, nper, nharm, ns, pAuxBuffer);
		sumValsIntBe2 += *pAuxBuffer;
	}
	//return only one point of integral over B^2
	*pAuxBuffer = sumValsIntBe2/Np;
}

//-------------------------------------------------------------------------

//double radTPrtclTrj::FocusPotentPerInOnePoint(const TVector3d& P0, const TVector3d& dsNlongU, const TVector3d& N1U, const TVector3d& N2U, double per, int nper, int nharm, int ns, double* pAuxBuffer)
double radTPrtclTrj::FocusPotentPerInOnePoint(const TVector3d& P0, const TVector3d& dsNlongU, const TVector3d& N1U, const TVector3d& N2U, double per, double nper, int nharm, int ns, double* pAuxBuffer)
{//assumes distances in [mm] at input; returns focusing potential in [T^2*mm^3]
 //also, returns IntBtr^2 in pAuxBuffer[0] in [T^2*m]
	double *ArrB1 = pAuxBuffer;
	double *ArrB2 = pAuxBuffer + ns;

    TVector3d ZeroVect(0.,0.,0.);
	radTFieldKey FieldKey;
	FieldKey.B_ = 1;

	radTField Field(FieldKey, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect);
	TVector3d& CurP = Field.P;
	TVector3d& CurB = Field.B;
	//CurP = P0 - dsNlongU;
	CurP = P0;
	for(int i=0; i<ns; i++)
	{
        CurB = ZeroVect;
        FldSrcPtr->B_genComp(&Field);
		ArrB1[i] = N1U*CurB;
		ArrB2[i] = N2U*CurB;
		CurP += dsNlongU;
	}

	const double RelPrec = 0.001;

    radTMagHarm *MagHarmArr1 = 0, *MagHarmArr2 = 0;
	int NumHarm1 = nharm, NumHarm2 = nharm;

	//double StepLong = per;
	//if(ns > 3) StepLong = per/(ns - 3);
	double StepLong = per/ns;

	//FindFieldHarmonics(ArrB1, ns, -StepLong, StepLong, per, 0., RelPrec, 'x', NumHarm1, MagHarmArr1);
	//FindFieldHarmonics(ArrB2, ns, -StepLong, StepLong, per, 0., RelPrec, 'z', NumHarm2, MagHarmArr2);
	AnalyzeForHarmonics(ArrB1, ns, per, RelPrec, 'x', NumHarm1, MagHarmArr1);
	AnalyzeForHarmonics(ArrB2, ns, per, RelPrec, 'z', NumHarm2, MagHarmArr2);

	int AmOfHarm = NumHarm1;
	if(AmOfHarm < NumHarm2) AmOfHarm = NumHarm2;

	double SumBe2Norm = 0., SumBe2 = 0.;
	double B1, B2, B1e2, B2e2, k1, k2;
	for(int i=0; i<AmOfHarm; i++)
	{
		double CurTerm = 0;
		if(i < NumHarm1) 
		{
			radTMagHarm *pMagHarm1 = MagHarmArr1 + i;
			B1 = pMagHarm1->B;
			B1e2 = B1*B1;
			k1 = pMagHarm1->HarmNo;
			CurTerm += B1e2/(k1*k1);
			SumBe2 += B1e2;
		}
		if(i < NumHarm2) 
		{
			radTMagHarm *pMagHarm2 = MagHarmArr2 + i;
			B2 = pMagHarm2->B;
			B2e2 = B2*B2;
			k2 = pMagHarm2->HarmNo;
			CurTerm += B2e2/(k2*k2);
			SumBe2 += B2e2;
		}
		SumBe2Norm += CurTerm;
	}

	if(MagHarmArr1 != 0) delete[] MagHarmArr1;
	if(MagHarmArr2 != 0) delete[] MagHarmArr2;

	*pAuxBuffer = 0.0005*nper*per*SumBe2; //integral over transverse B^2 in [T^2*m], added OC090509

	const double Pi	= 3.141592653589793238;
	const double Inv_8PiE2 = 1./(8.*Pi*Pi);

	return (Inv_8PiE2*nper*per*per*per)*SumBe2Norm;
}

//-------------------------------------------------------------------------

void radTPrtclTrj::AnalyzeForHarmonics(double* pB, int AmOfPts, double Per, double RelPrec, char XorZ, int& AmOfHarm, radTMagHarm*& MagHarmArr)
{
	if((pB == 0) || (AmOfPts <= 0) || (Per <= 0) || (RelPrec <= 0)) return;

	float *AuxDataContIn=0, *AuxDataContOut=0;
	double *CkArr=0, *PhikArr=0;
	int *HarmNoArr=0;

	AuxDataContIn = new float[AmOfPts << 1];
	if(AuxDataContIn == 0) throw 0;

	AuxDataContOut = new float[AmOfPts << 1];
	if(AuxDataContOut == 0) throw 0;
	double *tB = pB;
	float *tIn = AuxDataContIn;
	double MaxAbsB = 0;
	for(int i=0; i<AmOfPts; i++) 
	{
		double CurB = *(tB++);
		*(tIn++) = (float)CurB; *(tIn++) = 0.;
		
		double CurAbsB = ::fabs(CurB);
		if(MaxAbsB < CurAbsB) MaxAbsB = CurAbsB;
	}
	if(MaxAbsB <= 0)
	{
	    AnalyzeForHarmonics_DeleteAuxArrays(AuxDataContIn, AuxDataContOut, CkArr, PhikArr, HarmNoArr);
		return;
	}

	double Step = Per/double(AmOfPts);
	double Start = -0.5*Per;

	//srTFFT1DInfo FFT1DInfo;
	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.pInData = AuxDataContIn;
	FFT1DInfo.pOutData = AuxDataContOut;
	FFT1DInfo.Dir = -1;
	FFT1DInfo.xStep = Step;
	FFT1DInfo.xStart = Start;
	FFT1DInfo.Nx = AmOfPts;
	FFT1DInfo.HowMany = 1;
	FFT1DInfo.UseGivenStartTrValue = 0;

	//srTFFT1D FFT1D;
	CGenMathFFT1D FFT1D;
	FFT1D.Make1DFFT(FFT1DInfo);

	int HalfAmOfPts = AmOfPts >> 1;
	int MaxAmOfHarm = AmOfHarm;
	if(MaxAmOfHarm >= HalfAmOfPts) MaxAmOfHarm = HalfAmOfPts - 1;

	CkArr = new double[MaxAmOfHarm];
	if(CkArr == 0) throw MEMORY_ALLOCATION_FAILURE;

	PhikArr = new double[MaxAmOfHarm];
	if(PhikArr == 0) throw MEMORY_ALLOCATION_FAILURE;

	HarmNoArr = new int[MaxAmOfHarm];
	if(HarmNoArr == 0) throw MEMORY_ALLOCATION_FAILURE;

	double CoefMult = 2./Per;
	double AbsThreshold = RelPrec*MaxAbsB/CoefMult;

	float *tOutFFT = AuxDataContOut + AmOfPts + 2;
	double *tCkArr = CkArr, *tPhikArr = PhikArr;
	int *tHarmNoArr = HarmNoArr;
	int HarmCount = 0;
	for(int j=0; j<MaxAmOfHarm; j++)
	{
		double Ak = *(tOutFFT++);
		double Bk = *(tOutFFT++);
		if((::fabs(Ak) < AbsThreshold) && (::fabs(Bk) < AbsThreshold)) continue;

		double Ck = sqrt(Ak*Ak + Bk*Bk);
		if(Ck < AbsThreshold) continue;

		*(tCkArr++) = CoefMult*Ck;
		*(tPhikArr++) = TMathFunctions::Argument(Ak, -Bk);
		*(tHarmNoArr++) = (j + 1);
		HarmCount++;
	}
	if(HarmCount <= 0) 
	{
		AnalyzeForHarmonics_DeleteAuxArrays(AuxDataContIn, AuxDataContOut, CkArr, PhikArr, HarmNoArr);
		return;
	}

	MagHarmArr = new radTMagHarm[HarmCount];
	radTMagHarm *tMagHarmArr = MagHarmArr;
	tCkArr = CkArr;
	tPhikArr = PhikArr;
    tHarmNoArr = HarmNoArr;
	for(int k=0; k<HarmCount; k++)
	{
        tMagHarmArr->HarmNo = *(tHarmNoArr++);
		tMagHarmArr->XorZ = XorZ;
		tMagHarmArr->B = (*(tCkArr++));
		tMagHarmArr->Phase = *(tPhikArr++);
		tMagHarmArr++;
	}
	AmOfHarm = HarmCount;
    AnalyzeForHarmonics_DeleteAuxArrays(AuxDataContIn, AuxDataContOut, CkArr, PhikArr, HarmNoArr);
}

//-------------------------------------------------------------------------

//void radTPrtclTrj::ComputeFocusPotentDerivOnMesh2D(double* pFocPotData, int np1, double d1, char UseSpecStepForDeriv1, int np2, double d2, char UseSpecStepForDeriv2, double* pKickData1, double* pKickData2)
void radTPrtclTrj::ComputeFocusPotentDerivOnMesh2D(double* pFocPotData, int np1, double d1, char UseSpecStepForDeriv1, int np2, double d2, char UseSpecStepForDeriv2, double* pKickData1, double* pKickData2, char cKickUnit)
{//computes KickData in [T^2*m^2], assuming FocPotData in [T^2*mm^3]
	double *tFocPotData = pFocPotData;
    double *tKickData1 = pKickData1, *tKickData2 = pKickData2;

	//double InvStep1 = -0.5*0.000001/d1;
	//double InvStep2 = -0.5*0.000001/d2;
	
	double unitConst = 0.5e-06; //to get [T^2*m^2] //OC050413
	if((cKickUnit == 2) || (cKickUnit == 3)) //to get [micro-rad] (2) or [rad] (3)
	{
		//*=(0.299792458/ener)^2
		double mult = 0.299792458;
		if(Energy > 0) mult /= Energy;
		unitConst *= (mult*mult);

		if(cKickUnit == 2) unitConst *= 1e+06; //to get [micro-rad]
	}

	double InvStep1 = -unitConst/d1;
	double InvStep2 = -unitConst/d2;

	double HalfInvStep1 = 0.5*InvStep1;
	double QuartInvStep1 = 0.25*InvStep1;

	double HalfInvStep2 = 0.5*InvStep2;
	double QuartInvStep2 = 0.25*InvStep2;

	double a1, a2;
	int Off1m, Off1p, Off2m, Off2p;
	int np1_mi_1 = np1 - 1, np2_mi_1 = np2 - 1;

	if((!UseSpecStepForDeriv1) && (!UseSpecStepForDeriv2))
	{
		for(int j=0; j<np2; j++)
		{
			a2 = HalfInvStep2;
			Off2m = -np1; Off2p = np1;
			if(j == 0) { Off2m = 0; a2 = InvStep2;}
			else if(j == np2_mi_1) { Off2p = 0; a2 = InvStep2;}

			for(int i=0; i<np1; i++)
			{
                a1 = HalfInvStep1;
                Off1m = -1; Off1p = 1;
				if(i == 0) { Off1m = 0; a1 = InvStep1;}
				else if(i == np1_mi_1) { Off1p = 0; a1 = InvStep1;}

				*(tKickData1++) = ((*(tFocPotData + Off1p)) - (*(tFocPotData + Off1m)))*a1;
				*(tKickData2++) = ((*(tFocPotData + Off2p)) - (*(tFocPotData + Off2m)))*a2;
				tFocPotData++;
			}
		}
	}
	else if(UseSpecStepForDeriv1 && (!UseSpecStepForDeriv2))
	{
        a1 = HalfInvStep1;
		for(int j=0; j<np2; j++)
		{
			a2 = QuartInvStep2;
			Off2m = -(np1 << 1); Off2p = (np1 << 1);
			if(j == 0) { Off2m = 0; a2 = HalfInvStep2;}
			else if(j == np2_mi_1) { Off2p = 0; a2 = HalfInvStep2;}

			for(int i=0; i<np1; i++)
			{
                *(tKickData1++) = a1*((*(tFocPotData + 1)) - (*tFocPotData));
				*(tKickData2++) = a2*((*(tFocPotData + Off2p)) - (*(tFocPotData + Off2m)) + (*(tFocPotData + (1 + Off2p))) - (*(tFocPotData + (1 + Off2m))));
				tFocPotData += 2;
			}
		}
	}
	else if((!UseSpecStepForDeriv1) && UseSpecStepForDeriv2)
	{
		a2 = HalfInvStep2;
		for(int j=0; j<np2; j++)
		{
			for(int i=0; i<np1; i++)
			{
                a1 = QuartInvStep1;
                Off1m = -2; Off1p = 2;
                if(i == 0) { Off1m = 0; a1 = HalfInvStep1;}
                else if(i == np1_mi_1) { Off1p = 0; a1 = HalfInvStep1;}

                *(tKickData1++) = a1*((*(tFocPotData + Off1p)) - (*(tFocPotData + Off1m)) + (*(tFocPotData + (1 + Off1p))) - (*(tFocPotData + (1 + Off1m))));
                *(tKickData2++) = a2*((*(tFocPotData + 1)) - *(tFocPotData));
                tFocPotData += 2;
			}
		}
	}
	else if(UseSpecStepForDeriv1 && UseSpecStepForDeriv2)
	{
		a1 = HalfInvStep1;
		a2 = HalfInvStep2;
		for(int j=0; j<np2; j++)
		{
			for(int i=0; i<np1; i++)
			{
                *(tKickData1++) = a1*((*(tFocPotData + 1)) - (*tFocPotData));
                *(tKickData2++) = a2*((*(tFocPotData + 3)) - *(tFocPotData + 2));
                tFocPotData += 4;
			}
		}
	}
}

//-------------------------------------------------------------------------

void radTPrtclTrj::SetupTransverseMeshForKickCalc(double r1, int np1, double r2, int np2, double& d1, double& d2, double& step1, double& step2, char& UseSpecStepForDeriv1, char& UseSpecStepForDeriv2, double* pCoordDir1, double* pCoordDir2, long& AmOfTransvPtsFocPot)
{
	const double DefaultStepForDeriv = 0.1; //[mm] to tune

	step1 = r1;
	step2 = r2;
	if(np1 > 1) step1 = r1/(np1 - 1);
	if(np2 > 1) step2 = r2/(np2 - 1);

	double arg = -0.5*r1;
	double *tCoordDir = pCoordDir1;
	for(int i=0; i<np1; i++) { *(tCoordDir++) = arg; arg += step1;}

	arg = -0.5*r2;
	tCoordDir = pCoordDir2;
	for(int j=0; j<np2; j++) { *(tCoordDir++) = arg; arg += step2;}

	UseSpecStepForDeriv1 = 1;
	UseSpecStepForDeriv2 = 1;

	AmOfTransvPtsFocPot = (np1*np2) << 2;
	if(d1 <= 0.)
	{
		if(np1 <= 1) d1 = DefaultStepForDeriv;
		else
		{
			d1 = step1;
			UseSpecStepForDeriv1 = 0;
			AmOfTransvPtsFocPot >>= 1;
		}
	}
	if(d2 <= 0.)
	{
		if(np2 <= 1) d2 = DefaultStepForDeriv;
		else
		{
			d2 = step2;
			UseSpecStepForDeriv2 = 0;
			AmOfTransvPtsFocPot >>= 1;
		}
	}
}

//-------------------------------------------------------------------------

void radTPrtclTrj::ComputeSecondOrderKick(const TVector3d& P1Vect, const TVector3d& NlongVect, double* ArrLongDist, int AmOfLongGridPts, int ns, const TVector3d& N1Vect, double r1, int np1, double r2, int np2, double d1, double d2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pIBe2)
{
	double TrGridStep1 = 0., TrGridStep2 = 0.;
	char UseSpecStepForDeriv1 = 0, UseSpecStepForDeriv2 = 0;
	long AmOfTrPtsFocPot = 0;
	SetupTransverseMeshForKickCalc(r1, np1, r2, np2, d1, d2, TrGridStep1, TrGridStep2, UseSpecStepForDeriv1, UseSpecStepForDeriv2, pCoordDir1, pCoordDir2, AmOfTrPtsFocPot);

	double *pFocPotData = new double[AmOfTrPtsFocPot*AmOfLongGridPts];
	ComputeFocusPotentOnMesh3D(P1Vect, NlongVect, ArrLongDist, AmOfLongGridPts, ns, N1Vect, np1, TrGridStep1, d1, np2, TrGridStep2, d2, UseSpecStepForDeriv1, UseSpecStepForDeriv2, pFocPotData, pIBe2);
	ComputeFocusPotentDerivOnMesh3D(pFocPotData, AmOfLongGridPts, np1, d1, UseSpecStepForDeriv1, np2, d2, UseSpecStepForDeriv2, pKickData1, pKickData2);

	if(AmOfLongGridPts > 1)
	{
        long AmOfTransvPts = np1*np2;
        double *pAuxTotKickData1 = new double[AmOfTransvPts];
        double *pAuxTotKickData2 = new double[AmOfTransvPts];
        double *pAuxTotIBe2 = new double[AmOfTransvPts];

        SumUpPartialSecondOrderKicks(AmOfLongGridPts, AmOfTransvPts, pKickData1, pKickData2, pIBe2, pAuxTotKickData1, pAuxTotKickData2, pAuxTotIBe2);
		ArrangeAllSecondOrderKickData(AmOfLongGridPts, AmOfTransvPts, pAuxTotKickData1, pAuxTotKickData2, pAuxTotIBe2, pKickData1, pKickData2, pIBe2);

		if(pAuxTotKickData1 != 0) delete[] pAuxTotKickData1;
		if(pAuxTotKickData2 != 0) delete[] pAuxTotKickData2;
		if(pAuxTotIBe2 != 0) delete[] pAuxTotIBe2;
	}
	if(pFocPotData != 0) delete[] pFocPotData;
}

//-------------------------------------------------------------------------

void radTPrtclTrj::ComputeFocusPotentOnMesh3D(const TVector3d& P1Vect, const TVector3d& NlongVect, double* ArrLongDist, int AmOfLongGridPts, int ns, const TVector3d& N1Vect, int np1, double TrGridStep1, double d1, int np2, double TrGridStep2, double d2, char UseSpecStepForDeriv1, char UseSpecStepForDeriv2, double* pFocPotData, double* pIntBtrE2Data)
{// produces pFocPotData aligned as:
 // external loop - z; middle loop - x; internal loop - longitud.
	double InvLenNlong = 1./sqrt(NlongVect.x*NlongVect.x + NlongVect.y*NlongVect.y + NlongVect.z*NlongVect.z);
	TVector3d NlongU = InvLenNlong*NlongVect;
	double InvLenN1 = 1./sqrt(N1Vect.x*N1Vect.x + N1Vect.y*N1Vect.y + N1Vect.z*N1Vect.z);
	TVector3d N1U = InvLenN1*N1Vect;
	TVector3d N2U = N1U^NlongU;

	TVector3d d1N1U, d2N2U;

	int NumFocPotArrPerTrPoint = 1;
	if(UseSpecStepForDeriv1) 
	{ 
		NumFocPotArrPerTrPoint <<= 1;
        d1N1U = d1*N1U;
	}
	if(UseSpecStepForDeriv2) 
	{ 
		NumFocPotArrPerTrPoint <<= 1;
        d2N2U = d2*N2U;
	}
	double invNumFocPotArrPerTrPoint = 1./NumFocPotArrPerTrPoint;

	int *ArrLongNumPts = new int[AmOfLongGridPts];
	double *auxArIntBtrE2 = new double[AmOfLongGridPts];
	SetupLongNumPointsArray(ArrLongDist, AmOfLongGridPts, ns, ArrLongNumPts);

	double Arg2Min = -0.5*TrGridStep2*(np2 - 1), Arg1Min = -0.5*TrGridStep1*(np1 - 1);
	if(np2 == 1) Arg2Min = -0.5*TrGridStep2;
	if(np1 == 1) Arg1Min = -0.5*TrGridStep1;

	double arg2 = Arg2Min;
	TVector3d ArrTrInitP[4];
	double *tFocPotData = pFocPotData;
	double *tIntBtrE2Data = pIntBtrE2Data;

	for(int j=0; j<np2; j++)
	{
       TVector3d CurOffsetN2 = arg2*N2U;
        double arg1 = Arg1Min;

		for(int i=0; i<np1; i++)
		{
			TVector3d CurOffsetN1 = arg1*N1U;
			TVector3d CurP = P1Vect + (CurOffsetN2 + CurOffsetN1);

			SetupFourOrOneTansvInitPt(CurP, d1N1U, d2N2U, UseSpecStepForDeriv1, UseSpecStepForDeriv2, ArrTrInitP);

			for(int kk=0; kk<AmOfLongGridPts; kk++) tIntBtrE2Data[kk] = 0.;

			for(int k=0; k<NumFocPotArrPerTrPoint; k++)
			{
				ComputeFocusPotentArrayForOneTransvPoint(ArrTrInitP[k], NlongU, ArrLongDist, AmOfLongGridPts, ArrLongNumPts, N1U, N2U, tFocPotData, auxArIntBtrE2);
				tFocPotData += AmOfLongGridPts;

				for(int kk=0; kk<AmOfLongGridPts; kk++) tIntBtrE2Data[kk] += auxArIntBtrE2[kk];
			}
			if(NumFocPotArrPerTrPoint > 1)
			{
				for(int kk=0; kk<AmOfLongGridPts; kk++) tIntBtrE2Data[kk] *= invNumFocPotArrPerTrPoint;
			}
			tIntBtrE2Data += AmOfLongGridPts;

			arg1 += TrGridStep1;
		}
		arg2 += TrGridStep2;
	}
	delete[] ArrLongNumPts;
	delete[] auxArIntBtrE2;
}

//-------------------------------------------------------------------------

void radTPrtclTrj::SetupLongNumPointsArray(double* ArrLongDist, int AmOfLongGridPts, int ns, int* ArrLongNumPts)
{
	if(AmOfLongGridPts == 1) { *ArrLongNumPts = ns; return;}

	double MaxDist = ArrLongDist[AmOfLongGridPts - 1];
	double PrevDist = 0.;
	for(int i=0; i<AmOfLongGridPts; i++)
	{
		double CurDist = ArrLongDist[i];
		double RelNumPts = (CurDist - PrevDist)/MaxDist;
        PrevDist = CurDist;
		int CurNumPts = (int)(ns*RelNumPts);
		if(CurNumPts <= 2) CurNumPts = 3;
		ArrLongNumPts[i] = CurNumPts;
	}
}

//-------------------------------------------------------------------------

void radTPrtclTrj::ComputeFocusPotentArrayForOneTransvPoint(const TVector3d& InitP, const TVector3d& NlongU, double* ArrLongDist, int AmOfLongGridPts, int* ArrLongNumPts, const TVector3d& N1U, const TVector3d& N2U, double* pFocPotData, double* pIntBtrE2)
{//assumes distances in [mm] at input; returns focusing potential in [T^2*mm^3], longitudinally-integrated squared transverse magnetic field in [T^2*m]
	
	//TVector3d Nlong = N1U^N2U;
	//Nlong.Normalize();
	
	TVector3d ZeroVect(0.,0.,0.);
	radTFieldKey FieldKey;
	FieldKey.B_ = 1;
    radTField Field(FieldKey, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect);
    TVector3d& CurP = Field.P;
    TVector3d& CurB = Field.B;
    CurP = InitP;
    CurB = ZeroVect;
    
	FldSrcPtr->B_genComp(&Field);
	double CurB_N1U = CurB*N1U, CurB_N2U = CurB*N2U;
	TVector3d PrevBtr = (CurB_N1U*N1U) + (CurB_N2U*N2U);
	double prevBtrE2 = CurB_N1U*CurB_N1U + CurB_N2U*CurB_N2U;

	TVector3d IntegBtr(0,0,0);
	double PrevIntegBtrE2 = 0.;
	double PrevDist = 0.;

	for(int i=0; i<AmOfLongGridPts; i++)
	{
		int CurLongNumPts = ArrLongNumPts[i];
		double tSt = 1./(CurLongNumPts - 1);
		double CurDist = ArrLongDist[i];
		
		double CurStep = tSt*(CurDist - PrevDist);
		PrevDist = CurDist;
		double CurWeight = 0.5*CurStep;
        TVector3d dNlongU = CurStep*NlongU;

		double ExternSum = 0, intBtrE2 = 0.;
		for(int j=1; j<CurLongNumPts; j++)
		{
            CurP += dNlongU; CurB = ZeroVect;
            FldSrcPtr->B_genComp(&Field);

			CurB_N1U = CurB*N1U;
			CurB_N2U = CurB*N2U;
			TVector3d CurBtr = (CurB_N1U*N1U) + ((CurB*N2U)*N2U);
			
			double curBtrE2 = CurB_N1U*CurB_N1U + CurB_N2U*CurB_N2U;
			intBtrE2 += CurWeight*(prevBtrE2 + curBtrE2);
			prevBtrE2 = curBtrE2;

			IntegBtr += CurWeight*(PrevBtr + CurBtr);
			PrevBtr = CurBtr;

			double CurIntegBtrE2 = IntegBtr.x*IntegBtr.x + IntegBtr.y*IntegBtr.y + IntegBtr.z*IntegBtr.z;
			ExternSum += (PrevIntegBtrE2 + CurIntegBtrE2);
			PrevIntegBtrE2 = CurIntegBtrE2;
		}
		pFocPotData[i] = CurWeight*ExternSum;
		pIntBtrE2[i] = 0.001*intBtrE2; //to have it in [T^2*m]

		//IntegBtr = ZeroVect; //To leave or to remove???
	}
}

//-------------------------------------------------------------------------

void radTPrtclTrj::ComputeFocusPotentDerivOnMesh3D(double* pFocPotData, int AmOfLongGridPts, int np1, double d1, char UseSpecStepForDeriv1, int np2, double d2, char UseSpecStepForDeriv2, double* pKickData1, double* pKickData2)
{//computes KickData in [T^2*m^2], assuming FocPotData in [T^2*mm^3]
 //produces pKickData1, pKickData2 aligned as:
 //external loop - z; middle loop - x; internal loop - longitudinal

	double *pFocPotData_StartLong = pFocPotData;
    double *tKickData1 = pKickData1, *tKickData2 = pKickData2;

	double InvStep1 = -0.5*0.000001/d1;
	double HalfInvStep1 = 0.5*InvStep1;
	double QuartInvStep1 = 0.25*InvStep1;

	double InvStep2 = -0.5*0.000001/d2;
	double HalfInvStep2 = 0.5*InvStep2;
	double QuartInvStep2 = 0.25*InvStep2;

	double a1, a2;
	int Off1m, Off1p, Off2m, Off2p;
	int np1_mi_1 = np1 - 1, np2_mi_1 = np2 - 1;

	if((!UseSpecStepForDeriv1) && (!UseSpecStepForDeriv2))
	{
		long PerFocPotData = AmOfLongGridPts;
		for(int j=0; j<np2; j++)
		{
			a2 = HalfInvStep2;
			Off2m = -np1; Off2p = np1;
			if(j == 0) { Off2m = 0; a2 = InvStep2;}
			else if(j == np2_mi_1) { Off2p = 0; a2 = InvStep2;}

			Off2m *= AmOfLongGridPts; Off2p *= AmOfLongGridPts;

			for(int i=0; i<np1; i++)
			{
                a1 = HalfInvStep1;
                Off1m = -1; Off1p = 1;
				if(i == 0) { Off1m = 0; a1 = InvStep1;}
				else if(i == np1_mi_1) { Off1p = 0; a1 = InvStep1;}

				Off1m *= AmOfLongGridPts; Off1p *= AmOfLongGridPts;

				double *tFocPotData = pFocPotData_StartLong;
				for(int k=0; k<AmOfLongGridPts; k++)
                {
                    *(tKickData1++) = ((*(tFocPotData + Off1p)) - (*(tFocPotData + Off1m)))*a1;
                    *(tKickData2++) = ((*(tFocPotData + Off2p)) - (*(tFocPotData + Off2m)))*a2;
                    tFocPotData++;
				}
				pFocPotData_StartLong += PerFocPotData;
			}
		}
	}
	else if(UseSpecStepForDeriv1 && (!UseSpecStepForDeriv2))
	{
		long PerFocPotData = AmOfLongGridPts << 1;
        a1 = HalfInvStep1;
		for(int j=0; j<np2; j++)
		{
			a2 = QuartInvStep2;
			Off2m = -(np1 << 1); Off2p = (np1 << 1);
			if(j == 0) { Off2m = 0; a2 = HalfInvStep2;}
			else if(j == np2_mi_1) { Off2p = 0; a2 = HalfInvStep2;}

			Off2m *= AmOfLongGridPts; Off2p *= AmOfLongGridPts;

			for(int i=0; i<np1; i++)
			{
				double *tFocPotData = pFocPotData_StartLong;
                for(int k=0; k<AmOfLongGridPts; k++)
                {
                    *(tKickData1++) = a1*((*(tFocPotData + AmOfLongGridPts)) - (*tFocPotData));
                    *(tKickData2++) = a2*((*(tFocPotData + Off2p)) - (*(tFocPotData + Off2m)) + (*(tFocPotData + (AmOfLongGridPts + Off2p))) - (*(tFocPotData + (AmOfLongGridPts + Off2m))));
                    tFocPotData++;
				}
				pFocPotData_StartLong += PerFocPotData;
			}
		}
	}
	else if((!UseSpecStepForDeriv1) && UseSpecStepForDeriv2)
	{
		long PerFocPotData = AmOfLongGridPts << 1;
		a2 = HalfInvStep2;
		for(int j=0; j<np2; j++)
		{
			for(int i=0; i<np1; i++)
			{
                a1 = QuartInvStep1;
                Off1m = -2; Off1p = 2;
                if(i == 0) { Off1m = 0; a1 = HalfInvStep1;}
                else if(i == np1_mi_1) { Off1p = 0; a1 = HalfInvStep1;}

                Off1m *= AmOfLongGridPts; Off1p *= AmOfLongGridPts;

				double *tFocPotData = pFocPotData_StartLong;
                for(int k=0; k<AmOfLongGridPts; k++)
                {
                    *(tKickData1++) = a1*((*(tFocPotData + Off1p)) - (*(tFocPotData + Off1m)) + (*(tFocPotData + (AmOfLongGridPts + Off1p))) - (*(tFocPotData + (AmOfLongGridPts + Off1m))));
                    *(tKickData2++) = a2*((*(tFocPotData + AmOfLongGridPts)) - *(tFocPotData));
                    tFocPotData++;
				}
				pFocPotData_StartLong += PerFocPotData;
			}
		}
	}
	else if(UseSpecStepForDeriv1 && UseSpecStepForDeriv2)
	{
		long PerFocPotData = AmOfLongGridPts << 2;
		long TwoAmOfLongGridPts = AmOfLongGridPts << 1;
		long ThreeAmOfLongGridPts = AmOfLongGridPts*3;

		a1 = HalfInvStep1;
		a2 = HalfInvStep2;

		for(int j=0; j<np2; j++)
		{
			for(int i=0; i<np1; i++)
			{
				double *tFocPotData = pFocPotData_StartLong;
                for(int k=0; k<AmOfLongGridPts; k++)
                {
                    *(tKickData1++) = a1*((*(tFocPotData + AmOfLongGridPts)) - (*tFocPotData));
                    *(tKickData2++) = a2*((*(tFocPotData + ThreeAmOfLongGridPts)) - *(tFocPotData + TwoAmOfLongGridPts));
                    tFocPotData++;
				}
				pFocPotData_StartLong += PerFocPotData;
			}
		}
	}
}

//-------------------------------------------------------------------------

void radTPrtclTrj::SumUpPartialSecondOrderKicks(long AmOfLongGridPts, long AmOfTransvPts, double* pKickData1, double* pKickData2, double* pIBe2Data, double* pTotKickData1, double* pTotKickData2, double* pTotIBe2Data)
{
	double *tKickData1 = pKickData1;
	double *tKickData2 = pKickData2;
	double *tIBe2Data = pIBe2Data;

	double *tTotKickData1 = pTotKickData1;
	double *tTotKickData2 = pTotKickData2;
	double *tTotIBe2Data = pTotIBe2Data;

	for(long i=0; i<AmOfTransvPts; i++)
	{
		double SumKick1 = 0., SumKick2 = 0., SumIBe2 = 0.;
        for(long j=0; j<AmOfLongGridPts; j++)
		{
			SumKick1 += *(tKickData1++);
			SumKick2 += *(tKickData2++);
			SumIBe2 += *(tIBe2Data++);
		}
		*(tTotKickData1++) = SumKick1;
		*(tTotKickData2++) = SumKick2;
		*(tTotIBe2Data++) = SumIBe2;
	}
}

//-------------------------------------------------------------------------

void radTPrtclTrj::ArrangeAllSecondOrderKickData(long AmOfLongGridPts, long AmOfTransvPts, double* pTotKickData1, double* pTotKickData2, double* pTotIBe2, double* pKickData1, double* pKickData2, double* pIBe2)
{
	long TotAmOfData = (AmOfLongGridPts + 1)*AmOfTransvPts;
	double *pAuxContKick1 = new double[TotAmOfData];
	double *pAuxContKick2 = new double[TotAmOfData];
	double *pAuxContIBe2 = new double[TotAmOfData];

	double *tAuxContKick1 = pAuxContKick1;
	double *tAuxContKick2 = pAuxContKick2;
	double *tAuxContIBe2 = pAuxContIBe2;

	double *tTotKickData1 = pTotKickData1;
	double *tTotKickData2 = pTotKickData2;
	double *tTotIBe2 = pTotIBe2;

	for(long i=0; i<AmOfTransvPts; i++) //Total kicks first
	{
		*(tAuxContKick1++) = *(tTotKickData1++);
		*(tAuxContKick2++) = *(tTotKickData2++);
		*(tAuxContIBe2++) = *(tTotIBe2++);
	}

    for(long j=0; j<AmOfLongGridPts; j++) //then partial kicks
	{
		long OffsetKickData = j;
        for(long i=0; i<AmOfTransvPts; i++)
        {
            *(tAuxContKick1++) = *(pKickData1 + OffsetKickData);
            *(tAuxContKick2++) = *(pKickData2 + OffsetKickData);
            *(tAuxContIBe2++) = *(pIBe2 + OffsetKickData);
			OffsetKickData += AmOfLongGridPts;
		}
	}

	tAuxContKick1 = pAuxContKick1;
	tAuxContKick2 = pAuxContKick2;
	tAuxContIBe2 = pAuxContIBe2;

	double *tKickData1 = pKickData1;
	double *tKickData2 = pKickData2;
	double *tIBe2 = pIBe2;

    for(long k=0; k<=AmOfLongGridPts; k++)
	{
        for(long i=0; i<AmOfTransvPts; i++)
        {
			*(tKickData1++) = *(tAuxContKick1++);
			*(tKickData2++) = *(tAuxContKick2++);
			*(tIBe2++) = *(tAuxContIBe2++);
		}
	}

	if(pAuxContKick1 != 0) delete[] pAuxContKick1;
	if(pAuxContKick2 != 0) delete[] pAuxContKick2;
	if(pAuxContIBe2 != 0) delete[] pAuxContIBe2;
}

//-------------------------------------------------------------------------

void radTPrtclTrj::ComposeStrReportSecondOrderKickPer(const char* StrComment, double UndLen, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pBtE2Int, char*& pStrReport, char cKickUnit, char cOutFormat)
{// coordinates expected in mm at input

	const char sKickUnit_T2m2[] = "T2m2\0"; 
	const char sKickUnit_microrad[] = "micro-rad\0"; 
	const char sKickUnit_rad[] = "rad\0"; 
	const char *psKickUnit = sKickUnit_T2m2;
	if(cKickUnit == 2) psKickUnit = sKickUnit_microrad;
	else if(cKickUnit == 3) psKickUnit = sKickUnit_rad;

//#ifdef __GCC__
//	ostrstream OutStream;
//#else
	ostringstream OutStream;
//#endif

	OutStream << "# Author : Radia User" << endl;
	OutStream << "# ";
	if(StrComment != 0) OutStream << StrComment;
	OutStream << endl;
	OutStream << "# Undulator Length [m]" << endl;
	OutStream << UndLen*0.001 << endl;
	OutStream << "# Number of Horizontal Points" << endl;
	OutStream << np1 << endl;
	OutStream << "# Number of Vertical Points" << endl;
	OutStream << np2 << endl;

	//OutStream << "# Horizontal 2nd Order Kick [T2m2]" << endl;
	OutStream << "# Horizontal Kick [" << psKickUnit << "]" << endl; //OC050413
	OutStream << "START" << endl;
	//ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pKickData1);
	ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pKickData1, cOutFormat);

	//OutStream << "# Vertical 2nd Order Kick [T2m2]" << endl;
	OutStream << "# Vertical Kick [" << psKickUnit << "]" << endl;
	OutStream << "START" << endl;
	//ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pKickData2);
	ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pKickData2, cOutFormat);

	OutStream << "# Longitudinally Integrated Squared Transverse Magnetic Field [T2m]" << endl;
	OutStream << "START" << endl;
	//ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pBtE2Int);
	ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pBtE2Int, cOutFormat);

	OutStream << ends;

	long LenStrReport = (long)strlen(OutStream.str().c_str()) + 2;
	pStrReport = new char[LenStrReport];
	strcpy(pStrReport, OutStream.str().c_str());
	*(pStrReport + (LenStrReport - 1)) = '\0';
}

//-------------------------------------------------------------------------

//#ifdef __GCC__
//void radTPrtclTrj::ComposeStrReportSecondOrderKickPer_AddMainData(ostrstream& OutStream, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData)
//#else
//void radTPrtclTrj::ComposeStrReportSecondOrderKickPer_AddMainData(ostringstream& OutStream, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData)
void radTPrtclTrj::ComposeStrReportSecondOrderKickPer_AddMainData(ostringstream& OutStream, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData, char cOutFormat)
//#endif
{// coordinates expected in mm at input
	double AbsZeroTol = 1.e-09;
	char BufStr[16];
	char EmptyField[] = "              "; //14*" "
	char FormatStr[] = "% 14.5e";

	if(cOutFormat == 2)
	{//tab-delimited format
		EmptyField[1] = '\0';
	}

	OutStream << EmptyField;
	if(cOutFormat == 2) OutStream << '\t';

	int i, j;
	double dVal;
	int np1_mi_1 = np1 - 1;

	for(i=0; i<np1; i++)
	{
		dVal = 0.001*ApplyZeroTol(pCoordDir1[i], AbsZeroTol);
		if(cOutFormat == 1)
		{
			//sprintf(BufStr, FormatStr, 0.001*ApplyZeroTol(pCoordDir1[i], AbsZeroTol));
			sprintf(BufStr, FormatStr, dVal);
			OutStream << BufStr;
		}
		else if(cOutFormat == 2)
		{
			OutStream << dVal;
			if(i < np1_mi_1) OutStream << '\t';
		}
	}
	OutStream << endl;

	int np2_mi_1 = np2 - 1;
	for(j=0; j<np2; j++)
	{
		int IndRow = np2_mi_1 - j;
		dVal = 0.001*ApplyZeroTol(pCoordDir2[IndRow], AbsZeroTol);

		if(cOutFormat == 1)
		{
			//sprintf(BufStr, FormatStr, 0.001*ApplyZeroTol(pCoordDir2[IndRow], AbsZeroTol));
			sprintf(BufStr, FormatStr, dVal);
			OutStream << BufStr;
		}
		else if(cOutFormat == 2)
		{
			OutStream << dVal << '\t';
		}

        double *tKickData = pKickData + IndRow*np1;
        for(int k=0; k<np1; k++)
        {
			dVal = ApplyZeroTol(tKickData[k], AbsZeroTol);

			if(cOutFormat == 1)
			{
				//sprintf(BufStr, FormatStr, ApplyZeroTol(tKickData[k], AbsZeroTol));
				sprintf(BufStr, FormatStr, dVal);
				OutStream << BufStr;
			}
			else if(cOutFormat == 2)
			{
				OutStream << dVal;
				if(k < np1_mi_1) OutStream << '\t';
			}
        }
        OutStream << endl;
	}
}

//-------------------------------------------------------------------------

void radTPrtclTrj::ComposeStrReportSecondOrderKick(const char* StrComment, double* ArrLongDist, int lenArrLongDist, int np1, int np2, double* pCoordDir1, double* pCoordDir2, double* pKickData1, double* pKickData2, double* pIBe2Data, char*& pStrReport)
{// coordinates expected in mm at input
//#ifdef __GCC__
//	ostrstream OutStream;
//#else
	ostringstream OutStream;
//#endif

	char FormatStr[] = "% 14.5e";
	double AbsZeroTol = 1.e-09;

	OutStream << "# Author : Radia User" << endl;
	OutStream << "# ";
	if(StrComment != 0) OutStream << StrComment;
	OutStream << endl;
	OutStream << "# Total Length of Longitudinal Interval [m]" << endl;
	OutStream << ArrLongDist[lenArrLongDist - 1]*0.001 << endl;
	OutStream << "# Number of Horizontal Points" << endl;
	OutStream << np1 << endl;
	OutStream << "# Number of Vertical Points" << endl;
	OutStream << np2 << endl;

	OutStream << "# Total Horizontal 2nd Order Kick [T2m2]" << endl;
	OutStream << "START" << endl;
	ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pKickData1);
	
	OutStream << "# Total Vertical 2nd Order Kick [T2m2]" << endl;
	OutStream << "START" << endl;
	ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pKickData2);

	OutStream << "# Total Longitudinally Integrated Squared Transverse Magnetic Field [T2m]" << endl;
	OutStream << "START" << endl;
	ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, pIBe2Data);

	if(lenArrLongDist > 1)
	{
        OutStream << "# " << endl;
        OutStream << "# Partial 2nd Order Kicks and Longitudinally Integrated Squared Transverse Magnetic Field for Longitudinal Sub-Intervals [m]" << endl;

		long AmOfTransvPts = np1*np2;

		double *tKickData1 = pKickData1 + AmOfTransvPts;
		double *tKickData2 = pKickData2 + AmOfTransvPts;
		double *tIBe2Data = pIBe2Data + AmOfTransvPts;

		char BufStr[16];
		double PrevPos = 0;
		for(int i=0; i<lenArrLongDist; i++)
		{
	        OutStream << "# " << endl;

			double CurPos = ArrLongDist[i];
			sprintf(BufStr, FormatStr, 0.001*ApplyZeroTol(PrevPos, AbsZeroTol));
            OutStream << BufStr;
			sprintf(BufStr, FormatStr, 0.001*ApplyZeroTol(CurPos, AbsZeroTol));
            OutStream << BufStr;
            OutStream << endl;
			PrevPos = CurPos;

            OutStream << "# Horizontal 2nd Order Kick [T2m2]" << endl;
            OutStream << "START" << endl;
            ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, tKickData1);

			OutStream << "# Vertical 2nd Order Kick [T2m2]" << endl;
            OutStream << "START" << endl;
            ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, tKickData2);

			OutStream << "# Longitudinally Integrated Squared Transverse Magnetic Field [T2m]" << endl;
			OutStream << "START" << endl;
			ComposeStrReportSecondOrderKickPer_AddMainData(OutStream, np1, np2, pCoordDir1, pCoordDir2, tIBe2Data);

			tKickData1 += AmOfTransvPts;
			tKickData2 += AmOfTransvPts;
			tIBe2Data += AmOfTransvPts;
		}
	}
	OutStream << ends;

	long LenStrReport = (long)strlen(OutStream.str().c_str()) + 1;
	pStrReport = new char[LenStrReport];
	strcpy(pStrReport, OutStream.str().c_str());
}

//-------------------------------------------------------------------------
