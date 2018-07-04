/*-------------------------------------------------------------------------
*
* File name:      radmamet.h
*
* Project:        RADIA
*
* Description:    Some numerical methods
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADMAMET_H
#define __RADMAMET_H

#include "radsend.h"
#include "gmvect.h"
//#include "srercode.h"

#include <time.h>
#include <stdlib.h>

#ifndef NEGATIVE_NUM_INTEG_OF_MODIF_BESSEL_FUNC //in line with SRW
#define NEGATIVE_NUM_INTEG_OF_MODIF_BESSEL_FUNC 67 + 23000
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTMathLinAlgEq {
	int SizeOfMatr;

	double* vvLU_Dcmp;
	double* colInverseMatrix;
	int* indxInverseMatrix;

public:
	radTMathLinAlgEq(int);
	radTMathLinAlgEq() {}
	~radTMathLinAlgEq();

	void LU_Dcmp(double**, int*, double*);
	void LU_BkSb(double**, int*, double*);
	void InverseMatrix(double**, int, double**);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
/**
template<class T> class radTMathOrdDifEq {
	
	double *dym_rk4, *dyt_rk4, *yt_rk4;
	double *dysav_rks5, *ysav_rks5, *ytemp_rks5;
	double *y_ap, *dydx_ap;
	int kmax_ap, count_ap;
	double *xp_ap, **yp_ap, dxsav_ap; 

	radTSend Send;
	T* PtrT;
	void (T::*FunDerivs)(double,double*,double*);
	int AmOfEq;

	short OnPrc;
	double *PrecArray, EpsTol;
	int MaxAutoStp;

public:
	radTMathOrdDifEq(int, T*, void (T::*)(double,double*,double*), short, double* =NULL, double =1., int =5000);
	~radTMathOrdDifEq();

	void RungeKutta4(double*, double*, double, double, double*);
	int RungeKuttaStep5(double*, double*, double*, double, double*, double*);
	int AutoPropagate(double*, double, double, double, double, int*, int*);
	void Tabulation(double*, double, double, int, double*);
};

//-------------------------------------------------------------------------

template <class T> radTMathOrdDifEq<T>::
	radTMathOrdDifEq(int InAmOfEq, T* InPtrT, void (T::*InFunDerivs)(double,double*,double*), short InOnPrc, 
					 double* InPrecArray, double InEpsTol, int InMaxAutoStp)
{
	OnPrc = InOnPrc; AmOfEq = InAmOfEq;
	PtrT = InPtrT; FunDerivs = InFunDerivs;
	
	dym_rk4 = new double[AmOfEq];
	dyt_rk4 = new double[AmOfEq];
	yt_rk4 = new double[AmOfEq];

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

template <class T> radTMathOrdDifEq<T>::~radTMathOrdDifEq()
{
	delete[] dym_rk4;
	delete[] dyt_rk4;
	delete[] yt_rk4;

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

template <class T> void radTMathOrdDifEq<T>::
	RungeKutta4(double* y, double* dydx, double x, double h, double* yout)
{
	double xh, hh=0.5*h, h6=h/6.;
	xh=x+hh;

	int i;
	for(i=0; i<AmOfEq; i++) yt_rk4[i] = y[i] + hh*dydx[i];
	(PtrT->*FunDerivs)(xh, yt_rk4, dyt_rk4);
	for(i=0; i<AmOfEq; i++) yt_rk4[i] = y[i] + hh*dyt_rk4[i];
	(PtrT->*FunDerivs)(xh, yt_rk4, dym_rk4);
	for(i=0; i<AmOfEq; i++)
	{
		yt_rk4[i] = y[i] + h*dym_rk4[i];
		dym_rk4[i] += dyt_rk4[i];
	}
	(PtrT->*FunDerivs)(x+h, yt_rk4, dyt_rk4);
	for(i=0; i<AmOfEq; i++) yout[i] = y[i] + h6*(dydx[i]+dyt_rk4[i]+2.*dym_rk4[i]);
}

//-------------------------------------------------------------------------

template <class T> void radTMathOrdDifEq<T>::
	Tabulation(double* F0, double Xmin, double Xmax, int Np, double* Solution)
{
	double Step_x = (Xmax-Xmin)/double(Np-1);
	double x = Xmin;
	double* Y = new double[AmOfEq];
	double* dYdx;
	if(!OnPrc) dYdx = new double[AmOfEq];

	double MinStepAllowed = 1.E-10 * Step_x;
	int Nok=0, Nbad=0;

	int k;
	for(k=0; k<AmOfEq; k++) Y[k] = F0[k];

	int AmOfEq_p_1 = AmOfEq+1;

	int Np_mi_1 = Np-1;
	for(int i=0; i<Np; i++)
	{
		if(!OnPrc) (PtrT->*FunDerivs)(x, Y, dYdx);

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
	delete[] Y; 
	if(!OnPrc) delete[] dYdx;
}

//-------------------------------------------------------------------------

template <class T> int radTMathOrdDifEq<T>::
	RungeKuttaStep5(double* y, double* dydx, double* x, double htry, double* hdid, double* hnext)
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
		(PtrT->*FunDerivs)(*x, ytemp_rks5, dydx);
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

template <class T> int radTMathOrdDifEq<T>::
	AutoPropagate(double* ystart, double x1, double x2, double h1, double hmin, int* nok, int* nbad)
{
	int nstp, i;
	double xsav, x=x1, hnext, hdid, h;

	h = (x2>x1) ? fabs(h1) : -fabs(h1);
	*nok = (*nbad) = count_ap = 0;
	for(i=0; i<AmOfEq; i++) y_ap[i] = ystart[i];

	if(kmax_ap>0) xsav = x-2.*dxsav_ap;
	for(nstp=1; nstp<=MaxAutoStp; nstp++)
	{
		(PtrT->*FunDerivs)(x, y_ap, dydx_ap);
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

**/
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

inline double Max(double x1, double x2) { return (x1<x2)? x2 : x1;}

//-------------------------------------------------------------------------

inline double Abs(double x) { return (x<0)? -x : x;}

//-------------------------------------------------------------------------

inline double Norm(double x) { return (x<0)? -x : x;}

//-------------------------------------------------------------------------

inline double Norm(const TVector3d& Vect) { return Max( Max( Abs(Vect.x), Abs(Vect.y)), Abs(Vect.z));}

//-------------------------------------------------------------------------

inline void MakeItLarge(double& x) { x = 1.E+23;}

//-------------------------------------------------------------------------

inline void MakeItLarge(TVector3d& Vect) { Vect.x = Vect.y = Vect.z = 1.E+23;}

//-------------------------------------------------------------------------

inline void MakeItZero(double& x) { x = 0.;}

//-------------------------------------------------------------------------

inline void MakeItZero(TVector3d& Vect) { Vect.x = Vect.y = Vect.z = 0.;}

//-------------------------------------------------------------------------

template <class T, class DoubleOrTVector3d> void FormalOneFoldInteg(
	T* PtrT, void (T::*FunPtr)(double,DoubleOrTVector3d*), 
	int LenVal, double* AbsPrecAndLimitsArray, short* ElemCompNotFinished, DoubleOrTVector3d** IntegVal)
{// This uses Newton method (n=3)
 // This function requies memory allocation from outside:
 //
 // double* AbsPrecAndLimitsArray = new double[LenVal+2];
 //	DoubleOrTVector3d* IntegVal[6];
 //	for(int j=0; j<6; j++) IntegVal[j] = new DoubleOrTVector3d[LenVal];
 // short* ElemCompNotFinished = new short[LenVal];

	const double IntegWeight[] = {3./8., 9./8., 9./8., 3./4.};
	DoubleOrTVector3d LargeEntity, ZeroEntity; 
	MakeItLarge(LargeEntity); MakeItZero(ZeroEntity);

	DoubleOrTVector3d* LocSumArray = IntegVal[1];
	DoubleOrTVector3d* GenSumArray = IntegVal[2];
	DoubleOrTVector3d* InstIntegVal = IntegVal[3];
	DoubleOrTVector3d* PrInstIntegVal = IntegVal[4];
	DoubleOrTVector3d* ValF = IntegVal[5];
	for(int k=0; k<LenVal; k++)
	{
		PrInstIntegVal[k] = LargeEntity; 
		LocSumArray[k] = GenSumArray[k] = (IntegVal[0])[k] = ZeroEntity;
		ElemCompNotFinished[k] = 1;
	}

	double ArgStrt = AbsPrecAndLimitsArray[LenVal];
	double ArgFin = AbsPrecAndLimitsArray[LenVal+1];

	double StepArg, Arg;
	double IntervLength = ArgFin - ArgStrt;

	short IndForWeight, IndForPass;
	short NotFirstPass = 0;
	short MorePassesNeeded = 1;

	int AmOfPoi = 4;
	int AmOfPoi_mi_1;

	while(MorePassesNeeded)
	{
		AmOfPoi_mi_1 = AmOfPoi - 1; StepArg = IntervLength/AmOfPoi_mi_1;
		Arg = ArgStrt;
		IndForWeight = IndForPass = 0;

		for(int i=0; i<AmOfPoi; i++)
		{
			if(IndForPass==2) IndForPass = 0;
			if(IndForWeight==4) IndForWeight = 1;
			if(!(NotFirstPass && (IndForPass==0)))
			{
				if(i==AmOfPoi_mi_1) IndForWeight = 0;

				(PtrT->*FunPtr)(Arg, ValF);

				for(int kk=0; kk<LenVal; kk++)
					if(ElemCompNotFinished[kk]) LocSumArray[kk] += IntegWeight[IndForWeight]*ValF[kk];
			}
			IndForPass++; IndForWeight++; Arg += StepArg;
		}
		short MorePassesNeededProbe = 0;
		for(int kkk=0; kkk<LenVal; kkk++)
			if(ElemCompNotFinished[kkk])
			{
				GenSumArray[kkk] += LocSumArray[kkk];
				InstIntegVal[kkk] = StepArg*GenSumArray[kkk];
				if(Norm(InstIntegVal[kkk] - PrInstIntegVal[kkk]) < AbsPrecAndLimitsArray[kkk]) ElemCompNotFinished[kkk] = 0;
				else
				{
					MorePassesNeededProbe = 1;
					PrInstIntegVal[kkk] = InstIntegVal[kkk]; LocSumArray[kkk] = ZeroEntity;
				}
			}
		MorePassesNeeded = MorePassesNeededProbe;

		AmOfPoi = AmOfPoi_mi_1*2 + 1;
		NotFirstPass = 1;
	}
	for(int kkkk=0; kkkk<LenVal; kkkk++) (IntegVal[0])[kkkk] = InstIntegVal[kkkk];
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTRandGenLPTau {

	long NR[3][20];
	long IV[3][20];
	long NG, IQ;
	long NA, NB, NC, ND, NT;

	long iQ, mQ;
	long MaskArr[30], *MaskArrTrav;

public:
	radTRandGenLPTau()
	{
		long LocNR1[] = {1,3,5,15,17,51,85,255,257,771,1285,3855,4369,13107,21845,65535,65537,196611,327685,983055};
		long LocNR2[] = {1,1,7,11,13,61,67,79,465,721,823,4091,4125,4141,28723,45311,53505,250113,276231,326411};
		long *IV0Trav = IV[0], *IV1Trav = IV[1], *IV2Trav = IV[2], 
			 *LocNR1Trav = LocNR1, *LocNR2Trav = LocNR2,
			 *NR0Trav = NR[0], *NR1Trav = NR[1], *NR2Trav = NR[2];
		for(int i=0; i<20; i++)
		{
			*(NR0Trav++) = 1; 
			*(NR1Trav++) = *LocNR1Trav;
			*(NR2Trav++) = *LocNR2Trav;

			long OneBuf = 1 << 21;
			*(IV0Trav++) = OneBuf; 
			*(IV1Trav++) = (*(LocNR1Trav++)) << 20;
			*(IV2Trav++) = (*(LocNR2Trav++)) << 19;
		}
		NG = IQ = 0;
		NA = 436207616;
		NB = 872415232;
		NC = 4194303;
		ND = 536870912;
		NT = 64;

		MaskArrTrav = MaskArr; 
		long BufNumb = 1;
		for(int k=0; k<29; k++) { BufNumb <<= 1; *(MaskArrTrav++) = BufNumb;}
		MaskArrTrav = MaskArr;
		iQ = 0; mQ = 1;

		srand(1);
	}

	void LPTauSlow(int n, double* Q)
	{
		if((++iQ) >= *MaskArrTrav) { mQ++; MaskArrTrav++;}
		for(int j=0; j<n; j++)
		{
			long* NRj = NR[j];
			double s = 0.;
			for(int k=0; k<mQ; k++)
			{
				long ns = 0;
				long *NRjL_Trav = NRj, *LocMaskArrTrav1 = &(MaskArr[k]), *LocMaskArrTrav2 = MaskArr;
				for(int L=k; L<mQ; L++)
					ns += int(2*D(double(iQ)/(*(LocMaskArrTrav1++))) + 1E-08)*int(2*D(double(*(NRjL_Trav++))/(*(LocMaskArrTrav2++))) + 1E-08);
				s += D(0.5*ns)/(1 << k);
			}
			Q[j] = s;
		}
	}
	double D(double x) { return x - long(x);}

	void LPTauQuick(long i, int n, double q)
	{
	}

	void SimpleRand(int n, double* Q)
	{
		double D_RAND_MAX = double(RAND_MAX);
		for(int i=0; i<n; i++) *(Q++) = double(rand())/D_RAND_MAX;
	}
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class TMathFunctions {

public:

	/** Computes the mofified Bessel Function Kn(x)
	@see Method from V.O. Kostrum, NIM 172 (1980) p371-374
    nint=0 : computes k= Kn(x)
    nint=1 : computes k= ¤ Kn(x) dx from x to infinity
    nint=p : computes k= ¤¤..¤ Kn(x) dx from x to infinity p times
    prec : precision parameter (1e-5 is advised)
    h arbitrary parameter related to the cpu time (0.5 is advised) */
	static int Kmu(int ni, double mu, double x, double& f)
	{
		const double h = 0.5;
		const double prec = 1.E-05;

		if(ni < 0) return NEGATIVE_NUM_INTEG_OF_MODIF_BESSEL_FUNC;
		if((x > 100.) && (ni == 0)) return 0;
		double k = 0., c3 = 1E+10;
		long rr = 0;

		while(c3 > prec)
		{
			double rrh = (++rr)*h;
			double c1 = cosh(rrh);
			double c2 = cosh(mu*rrh);

			double c1eni = 1.;
			if(ni > 0) for(int i=0; i<ni; i++) c1eni *= c1;
			
			c3 = exp(-x*c1)*c2/c1eni;
			k += c3;
		}
		f = h*(0.5*exp(-x) + k);
		return 0;
	}

	/** Calculates principal value of argument of a complex number (-Pi < Phi <= Pi)  
	@param [out] x real part
	@param [out] y imaginary part
 	@return	calculated argument value
 	@see		... */
	static double Argument(double x, double y)
	{
		const double Pi = 3.1415926535897932;
		if(x == 0)
		{
			if(y < 0) return -0.5*Pi;
			else if(y == 0) return 0;
			else return 0.5*Pi;
		}
		if(y == 0)
		{
			if(x >= 0) return 0.;
			else return Pi;
		}
		if(y < 0)
		{
			if(x < 0) return -Pi + atan(y/x);
			else return atan(y/x); // x > 0
		}
		else // y > 0
		{
			if(x < 0) return Pi + atan(y/x);
			else return atan(y/x); // x > 0
		}
	}

	/** Searches for minimum and maximum value of an array */
	static void FindMinMax(double* pData, int Np, double& Min, int& IndMin, double& Max, int& IndMax)
	{
		IndMin = IndMax = -1;
		Min = 1.E+30; Max = -1.E+30;
		if((pData == 0) || (Np <= 0)) return;
		double* tData = pData;
		for(int i=0; i<Np; i++)
		{
			double CurVal = *(tData++);
			if(Min > CurVal) { Min = CurVal; IndMin = i;}
			if(Max < CurVal) { Max = CurVal; IndMax = i;}
		}
	}
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class TMathInterpol1D {

	double* AllCf;
	double** PlnCf;
	double Orig_sStart, Orig_sStep, InvStep;
	int OrigNp;

public:

	TMathInterpol1D(double* OrigF, int InOrigNp, double InOrig_sStart, double InOrig_sStep)
	{
		if((OrigF == 0) || (InOrigNp == 0)) throw 0;

        AllCf = 0;
		PlnCf = 0;
        OrigNp = InOrigNp;
        Orig_sStart = InOrig_sStart;
		Orig_sStep = InOrig_sStep;
        InvStep = 1./Orig_sStep;

		SetupPolinomCfs(OrigF);
	}
	TMathInterpol1D()
	{
		AllCf = 0;
		PlnCf = 0;
		OrigNp = 0;
	}
	~TMathInterpol1D()
	{
		DeallocateMemoryForCfs();
	}

	static void CubicPlnCfs(double f1, double f2, double fpr1, double fpr2, double sStep, double* aa)
	{
		double f1mf2_d_s1ms2 = (f2 - f1)/sStep;
		aa[0] = f1;
		aa[1] = fpr1;
		aa[2] = (3.*f1mf2_d_s1ms2 - 2.*fpr1 - fpr2)/sStep;
		aa[3] = (-2.*f1mf2_d_s1ms2 + fpr1 + fpr2)/(sStep*sStep);
	}
	void Interpolate(double sSt, double sStp, int Np, double* pInterpData)
	{
		double s = sSt;
		for(long i=0; i<Np; i++)
		{
			int Indx = int((s - Orig_sStart)/Orig_sStep); 
			if(Indx >= OrigNp - 1) Indx = OrigNp - 2;

			double sb = Orig_sStart + Indx*Orig_sStep;
			double smsb = s - sb;
			double *B_CfP = PlnCf[Indx];
			*(pInterpData++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			s += sStp;
		}
	}
	static double Derivative(double* f,	double h, int PoIndx, int AmOfPo=5)
	{
		if(AmOfPo==5)
		{
			if(PoIndx==2) return 0.08333333333333*(f[0]	- 8.*f[1] +	8.*f[3]	- f[4])/h;
			else if(PoIndx==1) return 0.08333333333333*(-3.*f[0] - 10.*f[1]	+ 18.*f[2] - 6.*f[3] + f[4])/h;
			else if(PoIndx==3) return 0.08333333333333*(-f[0] +	6.*f[1]	- 18.*f[2] + 10.*f[3] +	3.*f[4])/h;
			else if(PoIndx==0) return 0.5*(-3.*f[0]	+ 4.*f[1] -	f[2])/h;
			else if(PoIndx==4) return 0.5*(f[2]	- 4.*f[3] +	3.*f[4])/h;
			else return	1.E+23;
		}
		else if(AmOfPo==4)
		{
			if(PoIndx==1) return 0.5*(-f[0]	+ f[2])/h;
			else if(PoIndx==2) return 0.5*(-f[1] + f[3])/h;
			else if(PoIndx==0) return 0.5*(-3.*f[0]	+ 4.*f[1] -	f[2])/h;
			else if(PoIndx==3) return 0.5*(f[1]	- 4.*f[2] +	3.*f[3])/h;
			else return	1.E+23;
		}
		else if(AmOfPo==3)
		{
			if(PoIndx==1) return 0.5*(-f[0]	+ f[2])/h;
			else if(PoIndx==0) return 0.5*(-3.*f[0]	+ 4.*f[1] -	f[2])/h;
			else if(PoIndx==2) return 0.5*(f[0]	- 4.*f[1] +	3.*f[2])/h;
			else return	1.E+23;
		}
		else if(AmOfPo==2) return (-f[0] + f[1])/h;
		else return	1.E+23;
	}

private:

	void CompDerivForOrigData(double* OrigF, double* DerF)
	{
		if((OrigF == 0) || (DerF == 0) || (OrigNp <= 0)) throw 0;

		double SubArray[5];

		for(int k=0; k<5; k++) 
		{
			SubArray[k] = OrigF[k];
		}
		DerF[0] = Derivative(SubArray, Orig_sStep, 0);
		DerF[1] = Derivative(SubArray, Orig_sStep, 1);
		DerF[2] = Derivative(SubArray, Orig_sStep, 2);

		int LenFieldData_m_2 = OrigNp - 2;
		for(int i=3; i<LenFieldData_m_2; i++)
		{
			int i_m_2 = i - 2;
			for(int k=0; k<5; k++) 
			{
				SubArray[k] = OrigF[i_m_2 + k];
			}
			DerF[i] = Derivative(SubArray, Orig_sStep, 2);
		}

		DerF[LenFieldData_m_2] = Derivative(SubArray, Orig_sStep, 3);
		DerF[OrigNp - 1] = Derivative(SubArray, Orig_sStep, 4);
	}

	void SetupPolinomCfs(double* OrigF)
	{
		if((OrigF == 0) || (OrigNp <= 0)) throw 0;
		AllocateMemoryForCfs();

		double* DerF = new double[OrigNp];
		CompDerivForOrigData(OrigF, DerF);
		CalcPlnCfs(OrigF, DerF);
		if(DerF != 0) delete[] DerF;
	}
	void CalcPlnCfs(double* OrigF, double* DerF)
	{
		if((OrigF == 0) || (DerF == 0) || (OrigNp <= 0)) throw 0;
        int LenFieldData_m_1 = OrigNp - 1;
        double f1 = OrigF[0], f2;
        double fpr1 = DerF[0], fpr2;
        for(int is=1; is<OrigNp; is++)
        {
            f2 = OrigF[is];
            fpr2 = DerF[is];
            CubicPlnCfs(f1, f2, fpr1, fpr2, Orig_sStep, PlnCf[is-1]);
            f1 = f2; fpr1 = fpr2;
		}
	}
    void AllocateMemoryForCfs()
    {
		if(OrigNp <= 0) throw 0;
        DeallocateMemoryForCfs();
        int LenFieldData_m_1 = OrigNp - 1;

        PlnCf = new double*[LenFieldData_m_1];
        if(PlnCf == 0) throw 0;
        AllCf = new double[LenFieldData_m_1*4];
        if(AllCf == 0) { delete[] PlnCf; throw 0;}
        double* tAllCf = AllCf;
        for(int i=0; i<LenFieldData_m_1; i++) { PlnCf[i] = tAllCf; tAllCf += 4;}
    }
    void DeallocateMemoryForCfs()
    {
        if(AllCf != 0) { delete[] AllCf; AllCf = 0;}
        if(PlnCf != 0) { delete[] PlnCf; PlnCf = 0;}
    }
};

//-------------------------------------------------------------------------

#endif
