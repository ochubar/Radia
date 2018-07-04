/*-------------------------------------------------------------------------
*
* File name:      radmater.cpp
*
* Project:        RADIA
*
* Description:    Material relations and auxiliary functions for relaxation
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radmater.h"
#include "radmtra1.h"
#include "radg3d.h"
#include "auxparse.h"
#include "radrlmet.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
/**
void radTMaterial::SteerNewH(TVector3d& PrevH, TVector3d& InstantH, void* pvAuxRelax) //OC140103
{// may modify InstantH
	const float MinRelaxPar = (float)0.1; //(float)0.0001; //(float)0.24; //to tune
	const float RelaxParReduceCoef = (float)0.9; //to tune
	const double RelDifE2 = 0.99; //1E-05; //to tune
	const int AmOfBadPassesToAllow = 1; //1; //to tune

	if(pvAuxRelax == 0) return;

	//double AbsInstantH_E2 = InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z;
	double dx = InstantH.x - PrevH.x, dy = InstantH.y - PrevH.y, dz = InstantH.z - PrevH.z;
	double AbsDeltH_E2 = dx*dx + dy*dy + dz*dz;

	radTRelaxAuxData *pCurRelaxAuxData = (radTRelaxAuxData*)pvAuxRelax;

	float &PrevAbsDeltH_E2 = pCurRelaxAuxData->AbsDeltH;
	float &RelaxPar = pCurRelaxAuxData->RelaxPar;
	int &BadPassCount = pCurRelaxAuxData->BadPassCounts;

	//if(AbsDeltH_E2 >= PrevAbsDeltH_E2 - AbsInstantH_E2*RelMinDifE2)
	if(AbsDeltH_E2 >= RelDifE2*PrevAbsDeltH_E2)
	{
		if(BadPassCount >= AmOfBadPassesToAllow)
		{
			RelaxPar *= RelaxParReduceCoef;
			if(RelaxPar < MinRelaxPar) RelaxPar = MinRelaxPar;
			BadPassCount = 0;
		}
		else 
		{
			BadPassCount++;
		}
	}
	else
	{
		BadPassCount = 0;
	}
	PrevAbsDeltH_E2 = (float)AbsDeltH_E2;
	InstantH = RelaxPar*InstantH + (1. - RelaxPar)*PrevH;
}
**/
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//void radTNonlinearIsotropMaterial::FindNewH(TVector3d& InstantH, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnE2, radTg3dRelax* pMag, void* p) //OC140103
void radTNonlinearIsotropMaterial::FindNewH(TVector3d& InstantH, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnE2) //OC140103
{
	const double AbsHZeroTol = 1.E-10;
	const double AbsMZeroTol = 1.E-10;
	const int MaxIterToFindH = 15; //50; //15;

	TMatrix3d E; E.Str0.x = 1.; E.Str1.y = 1.; E.Str2.z = 1.;
	TMatrix3d BufMatr, InvBufMatr;

	double MisfitM = 3;
	double AbsInstantH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
	double f, InstKsi=0.;

	for(int i=0; i<MaxIterToFindH; i++)
	{
		f=0.;
		if(gLenArrayHM == 0)
		{
			if(AbsInstantH <= AbsHZeroTol)
			{
				for(int j=0; j<lenMs_ks; j++) InstKsi += ks[j];
				f = 0;
			}
			else
			{
				for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*AbsInstantH/Ms[i]);
				InstKsi = f/AbsInstantH;
			}
		}
		else
		{
			if(AbsInstantH <= AbsHZeroTol) 
			{
				InstKsi = *gdMdH; f = 0;
			}
			else
			{
				f = AbsMvsAbsH_Interpol(AbsInstantH, gArrayHM, gdMdH, gLenArrayHM);
				InstKsi = f/AbsInstantH;
			}
		}

		BufMatr = E - InstKsi*Matr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		InstantH = InvBufMatr*H_Ext;

		double NewAbsInstantH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
		double NewMisfitM = f - NewAbsInstantH*InstKsi;

		double AbsMisfitM = ::fabs(MisfitM), AbsNewMisfitM = ::fabs(NewMisfitM);

		double ProbNew = AbsMisfitM;
		if(NewMisfitM*MisfitM > 0) ProbNew += 0.5*AbsNewMisfitM;
		double ProbOld = AbsNewMisfitM;
		double Alpha = ProbNew/(ProbNew + ProbOld);
		AbsInstantH = Alpha*NewAbsInstantH + (1 - Alpha)*AbsInstantH;

		MisfitM = f - AbsInstantH*InstKsi;

/**
		TVector3d CurM = M(InstantH);
		f = CurM.Abs();

		TMatrix3d InstantKsiTensor;
		TVector3d InstMr, PrevInstantH = InstantH;
		DefineInstantKsiTensor(InstantH, InstantKsiTensor, InstMr);
        BufMatr = E - Matr*InstantKsiTensor;
        Matrix3d_inv(BufMatr, InvBufMatr);
		InstantH = InvBufMatr*(H_Ext + Matr*InstMr);
		TVector3d NewM = InstantKsiTensor*InstantH + InstMr;
        double NewMisfitM = f - NewM.Abs();

		double NewAbsInstantH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
		double AbsMisfitM = ::fabs(MisfitM), AbsNewMisfitM = ::fabs(NewMisfitM);

		if(AbsNewMisfitM < AbsMZeroTol) break;

		double ProbNew = AbsMisfitM;
		if(NewMisfitM*MisfitM > 0) ProbNew += 0.5*AbsNewMisfitM;
		double ProbOld = AbsNewMisfitM;
		double Alpha = ProbNew/(ProbNew + ProbOld);
		InstantH = Alpha*InstantH + (1. - Alpha)*PrevInstantH;
		AbsInstantH = InstantH.Abs();

		NewM = InstantKsiTensor*InstantH + InstMr;
		MisfitM = f - NewM.Abs();
**/

		if(MisfitM*MisfitM <= DesiredPrecOnMagnE2) break;
	}
}

//-------------------------------------------------------------------------

double radTNonlinearIsotropMaterial::FuncNewAbsH(double AbsH, const TMatrix3d& Matr, const TVector3d& H_Ext)
{
	const double AbsHZeroTol = 1.E-10;
	double f=0., InstKsi=0.;

	if(gLenArrayHM == 0)
	{
		if(AbsH <= AbsHZeroTol)
		{
			for(int j=0; j<lenMs_ks; j++) InstKsi += ks[j];
			f = 0;
		}
		else
		{
			for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*AbsH/Ms[i]);
			InstKsi = f/AbsH;
		}
	}
	else
	{
		if(AbsH <= AbsHZeroTol) 
		{
			InstKsi = *gdMdH; f = 0;
		}
		else
		{
			f = AbsMvsAbsH_Interpol(AbsH, gArrayHM, gdMdH, gLenArrayHM);
			InstKsi = f/AbsH;
		}
	}
	TMatrix3d E; E.Str0.x = 1.; E.Str1.y = 1.; E.Str2.z = 1.;
	TMatrix3d BufMatr, InvBufMatr;
    BufMatr = E - InstKsi*Matr;
	Matrix3d_inv(BufMatr, InvBufMatr);
	TVector3d AuxVectInstH = InvBufMatr*H_Ext;
	return AuxVectInstH.Abs();
}

//-------------------------------------------------------------------------

void radTNonlinearIsotropMaterial::DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstKsi, TVector3d& InstMr)
{
	double H0 = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);

	double f=0, dfdH=0;
	DefineScalarM_dMdH(H0, f, dfdH);

	double AbsZeroTolH = 1.e-10;

	if(::fabs(H0) < AbsZeroTolH)
	{
		InstKsi.Str0.x = dfdH; InstKsi.Str0.y = 0; InstKsi.Str0.z = 0; 
		InstKsi.Str1.x = 0; InstKsi.Str1.y = dfdH; InstKsi.Str1.z = 0; 
		InstKsi.Str2.x = 0; InstKsi.Str2.y = 0; InstKsi.Str2.z = dfdH;
		InstMr = RemMagn;
		return;
	}

	double H0e3 = H0*H0*H0;
	double InvH0e3 = 1./H0e3;
	double Hx0 = InstantH.x, Hy0 = InstantH.y, Hz0 = InstantH.z;
	double Hx0e2 = Hx0*Hx0, Hy0e2 = Hy0*Hy0, Hz0e2 = Hz0*Hz0;

	double H0_dfdH_mi_f = H0*dfdH - f;
	double BufOffDiag = InvH0e3*H0_dfdH_mi_f;
	double H0_dfdH = H0*dfdH;
	double f_mi_H0_dfdH_d_H0 = -H0_dfdH_mi_f/H0;

	InstKsi.Str0.x = InvH0e3*((Hy0e2 + Hz0e2)*f + Hx0e2*H0_dfdH);
	InstKsi.Str0.y = BufOffDiag*Hx0*Hy0; //InvH0e3*Hx0*Hy0*(H0*dfdH - f);
	InstKsi.Str0.z = BufOffDiag*Hx0*Hz0; //InvH0e3*Hx0*Hz0*(H0*dfdH - f);
	InstKsi.Str1.x = InstKsi.Str0.y;
	InstKsi.Str1.y = InvH0e3*((Hx0e2 + Hz0e2)*f + Hy0e2*H0_dfdH);
	InstKsi.Str1.z = BufOffDiag*Hy0*Hz0; //InvH0e3*Hy0*Hz0*(H0*dfdH - f);
	InstKsi.Str2.x = InstKsi.Str0.z;
	InstKsi.Str2.y = InstKsi.Str1.z;
	InstKsi.Str2.z = InvH0e3*((Hx0e2 + Hy0e2)*f + Hz0e2*H0_dfdH);

	InstMr.x = Hx0*f_mi_H0_dfdH_d_H0 + RemMagn.x;
	InstMr.y = Hy0*f_mi_H0_dfdH_d_H0 + RemMagn.y;
	InstMr.z = Hz0*f_mi_H0_dfdH_d_H0 + RemMagn.z;

/**
	double Der, f, dInstKsi;
	Der = f = dInstKsi = 0.;

	if(gLenArrayHM == 0)
	{
		if(H0 == 0.)
		{
			for(int j=0; j<lenMs_ks; j++) Der += ks[j];
			dInstKsi = Der;
		}
		else
		{
			for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*H0/Ms[i]);
			dInstKsi = f/H0;
		}
	}
	else
	{
		if(H0 == 0.) dInstKsi = *gdMdH;
		else dInstKsi = AbsMvsAbsH_Interpol(H0, gArrayHM, gdMdH, gLenArrayHM)/H0;
	}

	TVector3d Ksi_Str0(dInstKsi,0.,0.), Ksi_Str1(0.,dInstKsi,0.), Ksi_Str2(0.,0.,dInstKsi);
	InstKsi.Str0 = Ksi_Str0; InstKsi.Str1 = Ksi_Str1; InstKsi.Str2 = Ksi_Str2;
	InstMr = RemMagn;
**/
}

//-------------------------------------------------------------------------

void radTNonlinearIsotropMaterial::DefineScalarM_dMdH(double AbsH, double& f, double& dfdH)
{
	double AbsZeroTolH = 1.e-10;
	if(AbsH < 0.) AbsH = 0.;

	f = dfdH = 0.;
	if(gLenArrayHM == 0)
	{
        for(int i=0; i<lenMs_ks; i++) 
		{
			double ms_i = Ms[i], ks_i = ks[i];
			if(ms_i == 0.) continue;

			double arg = ks_i*AbsH/ms_i;
			f += ms_i*tanh(arg);

			double cosh_arg = cosh(arg);
			dfdH += ks_i/(cosh_arg*cosh_arg);
		}
	}
	else
	{
		if(AbsH < AbsZeroTolH) { f = 0; dfdH = *gdMdH; return;}
        AbsMvsAbsH_FuncAndDer_Interpol(AbsH, gArrayHM, gdMdH, gLenArrayHM, f, dfdH);
	}
}

//-------------------------------------------------------------------------

void radTNonlinearIsotropMaterial::DefineScalarM(double AbsInstantH, double& f, double& InstKsi)
{
	const double AbsHZeroTol = 1.E-09;
	if(gLenArrayHM == 0)
	{
		if(AbsInstantH <= AbsHZeroTol)
		{
			for(int j=0; j<lenMs_ks; j++) InstKsi += ks[j];
			f = 0;
		}
		else
		{
			for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*AbsInstantH/Ms[i]);
			InstKsi = f/AbsInstantH;
		}
	}
	else
	{
		if(AbsInstantH <= AbsHZeroTol) 
		{
			InstKsi = *gdMdH; f = 0;
		}
		else
		{
			//f = AbsMvsAbsH_Interpol(AbsInstantH);
			f = AbsMvsAbsH_Interpol(AbsInstantH, gArrayHM, gdMdH, gLenArrayHM);
			InstKsi = f/AbsInstantH;
		}
	}
}

//-------------------------------------------------------------------------

double radTNonlinearIsotropMaterial::Derivative5(TVector2d* f, int PoIndx)
{
	double x0 = f[0].x, x1 = f[1].x, x2 = f[2].x, x3 = f[3].x, x4 = f[4].x;
	x1 -= x0; x2 -= x0; x3 -= x0; x4 -= x0; 
	double f0 = f[0].y, f1 = f[1].y, f2 = f[2].y, f3 = f[3].y, f4 = f[4].y;

	char IsIncreasing = ((x1 > 0.) && (x2 > x1) && (x3 > x2) && (x4 > x3)) && 
						((f0 <= f1) && (f1 <= f2) && (f2 <= f3) && (f3 <= f4));

	double x1mix2 = (x1-x2), x1mix3 = (x1-x3), x1mix4 = (x1-x4), x2mix3 = (x2-x3), x2mix4 = (x2-x4), x3mix4 = (x3-x4);
	double x1e2 = x1*x1, x2e2 = x2*x2, x3e2 = x3*x3, x4e2 = x4*x4;
	if(PoIndx==0)
	{
		double x1mix3e2 = x1mix3*x1mix3;
		double Der = (f4*x1e2*x1mix2*x2e2*x1mix3*x2mix3*x3e2 + 
					(-f3*x1e2*x1mix2*x2e2*x1mix4*x2mix4 + 
					x3e2*(f2*x1e2*x1mix3*x1mix4 - f1*x2e2*x2mix3*x2mix4)*
					x3mix4)*x4e2 - f0*x1mix2*x1mix3*x2mix3*x1mix4*
					x2mix4*x3mix4*(x1*x2*x3 + x2*x3*x4 + x1*(x2 + x3)*x4))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx==1)
	{
		double x1e3 = x1e2*x1;
		double x1mix2e2 = x1mix2*x1mix2, x1mix3e2 = x1mix3*x1mix3, x1mix4e2 = x1mix4*x1mix4;
		double x1mix3e4 = x1mix3e2*x1mix3e2;
		double Der = (-f4*x1e2*x1mix2e2*x2*x1mix3e2*x2mix3*x3 + 
					f0*x1mix2e2*x1mix3e2*x2mix3*x1mix4e2*x2mix4*x3mix4 + 
					x4*(f3*x1e2*x1mix2e2*x2*x1mix4e2*x2mix4 - 
					x3*x3mix4*(f2*x1e2*x1mix3e2*x1mix4e2 + f1*x2*x2mix3*x2mix4*
					(-4*x1e3 + x2*x3*x4 + 3*x1e2*(x2 + x3 + x4) - 
					2*x1*(x3*x4 + x2*(x3 + x4))))))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx==2)
	{
		double x1mix2e2 = x1mix2*x1mix2, x2mix3e2 = x2mix3*x2mix3, x1mix3e2 = x1mix3*x1mix3, x2mix4e2 = x2mix4*x2mix4;
		double Der = (f4*x1*x1mix2e2*x2e2*x1mix3*x2mix3e2*x3 - 
					f0*x1mix2e2*x1mix3*x2mix3e2*x1mix4*x2mix4e2*x3mix4 + 
					x4*(-f3*x1*x1mix2e2*x2e2*x1mix4*x2mix4e2 + 
					x3*x3mix4*(f1*x2e2*x2mix3e2*x2mix4e2 + f2*x1*x1mix3*x1mix4*
					(x1*(3*x2e2 + x3*x4 - 2*x2*(x3 + x4)) + 
					x2*(-4*x2e2 - 2*x3*x4 + 3*x2*(x3 + x4))))))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx==3)
	{
		double x1mix3e2 = x1mix3*x1mix3, x2mix3e2 = x2mix3*x2mix3, x3mix4e2 = x3mix4*x3mix4;
		double x1mix3e3 = x1mix3e2*x1mix3;
		double Der = (-f4*x1*x1mix2*x2*x1mix3e2*x2mix3e2*x3e2 + 
					f0*x1mix2*x1mix3e2*x2mix3e2*x1mix4*x2mix4*x3mix4e2 + 
					x4*(f3*x1*x1mix2*x2*x1mix4*x2mix4*
					(x3*(2*x1*x2 - 3*(x1 + x2)*x3 + 4*x3e2) - 
					(x1*x2 - 2*(x1 + x2)*x3 + 3*x3e2)*x4) + 
					x3e2*x3mix4e2*(f2*x1*x1mix3e2*x1mix4 + 
					f1*x2*x2mix3e2*(-x2 + x4))))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx==4)
	{
		double x4e3 = x4e2*x4;
		double x2mix4e2 = x2mix4*x2mix4, x3mix4e2 = x3mix4*x3mix4, x1mix3e2 = x1mix3*x1mix3, x1mix4e2 = x1mix4*x1mix4;
		double x1mix3e3 = x1mix3e2*x1mix3;
		double Der = (-f0*x1mix2*x1mix3*x2mix3*x1mix4e2*x2mix4e2*
					x3mix4e2 + (f3*x1*x1mix2*x2*x1mix4e2*x2mix4e2 + 
					x3*(-f2*x1*x1mix3*x1mix4e2 + f1*x2*x2mix3*x2mix4e2)*
					x3mix4e2)*x4e2 + f4*x1*x1mix2*x2*x1mix3*x2mix3*x3*
					(x1*x2*x3 - 2*(x2*x3 + x1*(x2 + x3))*x4 + 3*(x1 + x2 + x3)*x4e2 - 4*x4e3))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else return 0.;
}

//-------------------------------------------------------------------------

double radTNonlinearIsotropMaterial::Derivative3(TVector2d* f, int PoIndx)
{
	double x0 = f[0].x, x1 = f[1].x, x2 = f[2].x;
	x1 -= x0; x2 -= x0;
	double f0 = f[0].y, f1 = f[1].y, f2 = f[2].y;

	char IsIncreasing = ((x1 > 0.) && (x2 > x1)) && ((f0 <= f1) && (f1 <= f2));

	double x1e2 = x1*x1, x2e2 = x2*x2;
	double x1mix2 = (x1-x2);
	if(PoIndx == 0)
	{
		double Der = (f0*x1e2 - f2*x1e2 - f0*x2e2 + f1*x2e2)/(-x1*x2*x1mix2);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx == 1)
	{
		double Der = (-f2*x1e2 + f0*x1mix2*x1mix2 + f1*(2*x1 - x2)*x2)/(x1*x1mix2*x2);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx == 2)
	{
		double Der = (f2*x1*(x1 - 2*x2) - f0*x1mix2*x1mix2 + f1*x2e2)/(x1*x1mix2*x1mix2*x2);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else return 0.;
}

//-------------------------------------------------------------------------

//void radTNonlinearIsotropMaterial::Compute_dMdH(double& MaxKsi)
void radTNonlinearIsotropMaterial::Compute_dMdH(TVector2d* ArrayHM, double* dMdH, int LenArrayHM, double& MaxKsi)
{
	if((ArrayHM == 0) || (dMdH == 0) || (LenArrayHM <= 0)) return;
	double *tdMdH = dMdH;
	TVector2d *tHM = ArrayHM;

	MaxKsi = 0;

	int LenArrayHM_mi_4 = LenArrayHM - 4;
	int LenArrayHM_mi_3 = LenArrayHM - 3;
	int LenArrayHM_mi_2 = LenArrayHM - 2;
	int LenArrayHM_mi_1 = LenArrayHM - 1;

	*tdMdH = Derivative3(tHM, 0);
	if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	tdMdH++;

	*tdMdH = Derivative3(tHM, 1);
	if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	tdMdH++; tHM++;

	*tdMdH = Derivative3(tHM, 1);
	if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	tdMdH++;

	for(int i=3; i<LenArrayHM_mi_4; i++)
	{
		*tdMdH = Derivative5(tHM, 2);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
		tdMdH++; tHM++;
	}

	if(LenArrayHM >= 7) //OC130604
	{
		*tdMdH = Derivative5(tHM, 2);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
		tdMdH++; tHM += 2;

		*tdMdH = Derivative3(tHM, 1);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
		tdMdH++; tHM++;

		*tdMdH = Derivative3(tHM, 1);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
		tdMdH++; tHM++;

		*tdMdH = ((tHM+1)->y - tHM->y)/((tHM+1)->x - tHM->x);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	}
	else if(LenArrayHM >= 6)
	{
		*tdMdH = Derivative5(tHM, 2);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
		tdMdH++; tHM += 2;

		*tdMdH = Derivative3(tHM, 1);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
		tdMdH++; tHM++;

		*tdMdH = ((tHM+1)->y - tHM->y)/((tHM+1)->x - tHM->x);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	}
	else if(LenArrayHM >= 5)
	{
		*tdMdH = Derivative3(tHM, 1);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
		tdMdH++; tHM++;

		*tdMdH = ((tHM+1)->y - tHM->y)/((tHM+1)->x - tHM->x);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	}
	else if(LenArrayHM >= 4)
	{
		*tdMdH = ((tHM+1)->y - tHM->y)/((tHM+1)->x - tHM->x);
		if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	}

	//*tdMdH = Derivative3(tHM, 1);
	//if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	//tdMdH++; tHM++;

	//*tdMdH = Derivative3(tHM, 1);
	//if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;
	//tdMdH++; tHM++;

	//*tdMdH = ((tHM+1)->y - tHM->y)/((tHM+1)->x - tHM->x);
	//if(MaxKsi < *tdMdH) MaxKsi = *tdMdH;

	//CheckAndCorrect_dMdH();
	CheckAndCorrect_dMdH(ArrayHM, dMdH, LenArrayHM);
}

//-------------------------------------------------------------------------

//void radTNonlinearIsotropMaterial::CheckAndCorrect_dMdH()
void radTNonlinearIsotropMaterial::CheckAndCorrect_dMdH(TVector2d* ArrayHM, double* dMdH, int LenArrayHM)
{
	if((ArrayHM == 0) || (dMdH == 0) || (LenArrayHM <= 0)) return;

	int LenArrayHM_mi_1 = LenArrayHM - 1;
	double *tdMdH = dMdH + LenArrayHM_mi_1;
	TVector2d *tHM = ArrayHM + LenArrayHM_mi_1;

	for(int i=0; i<LenArrayHM_mi_1; i++)
	{
		tdMdH--; tHM--;
		double x2 = (tHM+1)->x - tHM->x;
		double f1 = tHM->y, f2 = (tHM+1)->y;
		double fd1 = *tdMdH, fd2 = *(tdMdH+1);

		double x2e2 = x2*x2, f1mif2 = f1 - f2;
		double x2e3 = x2e2*x2;
		
		double a0 = f1;
		double a1 = fd1;
		double a2 = -((3*f1mif2 + (2*fd1 + fd2)*x2)/(x2e2));
		double a3 = -((-2*f1mif2 - (fd1 + fd2)*x2)/x2e3);

		double D = a2*a2 - 3.*a1*a3;
		double xc1 = 1.E+23, xc2 = 1.E+23;
		if(a3 != 0.)
		{
			if(D >= 0.)
			{
				double R = sqrt(D), Buf = 1./(3.*a3);
				xc1 = (-a2 - R)*Buf;
				xc2 = (-a2 + R)*Buf;
			}
		}

		char IsIncreasing = ((f1 < f2) && (fd1 > 0.) && (fd2 > 0.));
		char IsDecreasing = ((f1 > f2) && (fd1 < 0.) && (fd2 < 0.));
		char IsMonotone = (IsIncreasing || IsDecreasing);

		char xc1IsInside = 0, xc2IsInside = 0;
		if(a3 != 0.)
		{
			xc1IsInside = ((0. < xc1) && (xc1 < x2));
			xc2IsInside = ((0. < xc2) && (xc2 < x2));
		}

		char CorrectionNeeded = (IsMonotone && (D > 0.) && (xc1IsInside || xc2IsInside));
		if(CorrectionNeeded)
		{
			double LocD = -3*fd2*x2e3*(4*f1mif2 + fd2*x2);
			if(LocD >= 0.)
			{
				double LocR = sqrt(LocD), LocBuf1 = 1./(2*x2e2), LocBuf2 = -x2*(6*f1mif2 + fd2*x2);
				double fd1a = (LocBuf2 - LocR)*LocBuf1, fd1b = (LocBuf2 + LocR)*LocBuf1;
				double ra = fabs(fd1 - fd1a), rb = fabs(fd1 - fd1b);
				*tdMdH = (ra < rb)? fd1a : fd1b;
			}
		}
	}
}

//-------------------------------------------------------------------------

//double radTNonlinearIsotropMaterial::AbsMvsAbsH_Interpol(double AbsH)
double radTNonlinearIsotropMaterial::AbsMvsAbsH_Interpol(double AbsH, TVector2d* ArrayHM, double* dMdH, int LenArrayHM)
{
	TVector2d *tHM = ArrayHM;
	int Indx = 0;
	for(int i=0; i<LenArrayHM; i++)
	{
		if(tHM->x > AbsH) break;
		tHM++; Indx++;
	}
	tHM--; Indx--; 

	if(Indx < 0) return ArrayHM->y;
	if(Indx >= (LenArrayHM - 1)) 
	{
		TVector2d *pHM = ArrayHM + (LenArrayHM - 1);

		double Arg = AbsH - pHM->x;
		return pHM->y + Arg*dMdH[LenArrayHM - 1];
	}

	double Arg = AbsH - tHM->x;
	double Step = (tHM + 1)->x - tHM->x;
	double f1 = tHM->y, f2 = (tHM + 1)->y;
	double fpr1 = dMdH[Indx], fpr2 = dMdH[Indx + 1];
	double a[4];
	CubPln(Step, f1, f2, fpr1, fpr2, a);
	return a[0] + Arg*(a[1] + Arg*(a[2] + Arg*a[3]));
}

//-------------------------------------------------------------------------

//double radTNonlinearIsotropMaterial::AbsHvsAbsM_Interpol(double AbsM)
double radTNonlinearIsotropMaterial::AbsHvsAbsM_Interpol(double AbsM, TVector2d* ArrayHM, double* dMdH, int LenArrayHM)
{
	TVector2d *tHM = ArrayHM;
	int Indx = 0;
	for(int i=0; i<LenArrayHM; i++)
	{
		if(tHM->y > AbsM) break;
		tHM++; Indx++;
	}
	tHM--; Indx--; 

	if(Indx < 0) return ArrayHM->x;
	if(Indx >= (LenArrayHM - 1)) return ArrayHM[LenArrayHM - 1].x;

	double Arg = AbsM - tHM->y;
	double Step = (tHM + 1)->y - tHM->y;
	double f1 = tHM->x, f2 = (tHM + 1)->x;
	double fpr1 = 1./dMdH[Indx], fpr2 = 1./dMdH[Indx + 1];

	double a[4];
	CubPln(Step, f1, f2, fpr1, fpr2, a);
	return a[0] + Arg*(a[1] + Arg*(a[2] + Arg*a[3]));
}

//-------------------------------------------------------------------------

//void radTNonlinearIsotropMaterial::AbsMvsAbsH_FuncAndDer_Interpol(double AbsH, double& f, double& fDer)
void radTNonlinearIsotropMaterial::AbsMvsAbsH_FuncAndDer_Interpol(double AbsH, TVector2d* ArrayHM, double* dMdH, int LenArrayHM, double& f, double& fDer)
{
	TVector2d *tHM = ArrayHM;
	int Indx = 0;
	for(int i=0; i<LenArrayHM; i++)
	{
		if(tHM->x > AbsH) break;
		tHM++; Indx++;
	}
	tHM--; Indx--; 

	if(Indx < 0) 
	{
		f = ArrayHM->y;
		fDer = *dMdH;
		return;
	}
	if(Indx >= (LenArrayHM - 1)) 
	{
		f = ArrayHM[LenArrayHM - 1].y;
		fDer = 0.;
		return;
	}

	double Arg = AbsH - tHM->x;
	double Step = (tHM + 1)->x - tHM->x;
	double f1 = tHM->y, f2 = (tHM + 1)->y;
	double fpr1 = dMdH[Indx], fpr2 = dMdH[Indx + 1];
	double a[4];
	CubPln(Step, f1, f2, fpr1, fpr2, a);
	f = a[0] + Arg*(a[1] + Arg*(a[2] + Arg*a[3]));
	fDer = a[1] + Arg*(2*a[2] + Arg*3*a[3]);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTHMatDBVect radTMaterial::SetupMaterDB()
{
	radTHMatDBVect LocDB;
	radTHMatDB h001(new radTLinearAnisotropMaterialDB("NdFeB", 0.06, 0.17, 1.2)); LocDB.push_back(h001);
	radTHMatDB h002(new radTLinearAnisotropMaterialDB("SmCo5", 0.005,0.04, 0.85)); LocDB.push_back(h002);
	radTHMatDB h003(new radTLinearAnisotropMaterialDB("Sm2Co17", 0.005, 0.04, 1.05)); LocDB.push_back(h003);
	radTHMatDB h004(new radTLinearAnisotropMaterialDB("Ferrite", 0.07, 0.2, 0.35)); LocDB.push_back(h004);

	radTHMatDB h101(new radTNonlinearIsotropMaterialDB("Xc06", 1.362,0.2605,0.4917, 2118.,63.06,17.138)); LocDB.push_back(h101);
	radTHMatDB h102(new radTNonlinearIsotropMaterialDB("XcO6", 1.362,0.2605,0.4917, 2118.,63.06,17.138)); LocDB.push_back(h102);
	radTHMatDB h103(new radTNonlinearIsotropMaterialDB("Steel37", 1.1488,0.4268,0.4759, 1596.3,133.11,18.713)); LocDB.push_back(h103);
	radTHMatDB h104(new radTNonlinearIsotropMaterialDB("Steel42", 1.441,0.2912,0.3316, 968.66,24.65,8.3)); LocDB.push_back(h104);
	radTHMatDB h105(new radTNonlinearIsotropMaterialDB("AFK502", 1.788,0.437,0.115, 10485.,241.5,7.43)); LocDB.push_back(h105);
	radTHMatDB h106(new radTNonlinearIsotropMaterialDB("AFK1", 1.704,0.493,0.152, 2001.,38.56,1.24)); LocDB.push_back(h106);
	
	return LocDB;
}

radTHMatDBVect radTMaterial::MaterDB = SetupMaterDB();

//-------------------------------------------------------------------------

int radTMaterial::GetMaterIndexDB(const char* InName)
{
	if(InName == 0) return -1;

	char cLocBuffer[100], cInLocBuffer[100];
	cLocBuffer[0] = '\0';
	strcpy(cInLocBuffer, InName);
	CAuxParse::toUpperCase(cInLocBuffer);

	for(int i=0; i<(int)(MaterDB.size()); i++)
	{
		strcpy(cLocBuffer, MaterDB[i].rep->Name);
		CAuxParse::toUpperCase(cLocBuffer);
		if(strcmp(cInLocBuffer, cLocBuffer) == 0) return i;
	}
	return -1;
}

//-------------------------------------------------------------------------

radTMaterial* radTMaterial::SetupStandardMater(const char* InName, double InMr)
{
	int IndMater = GetMaterIndexDB(InName);
	if(IndMater < 0) return 0;

	radTMaterialDB *pMatDB = MaterDB[IndMater].rep;
	if(pMatDB == 0) return 0;

	return pMatDB->SetupMater(InMr);
}

//-------------------------------------------------------------------------

radTMaterial* radTLinearAnisotropMaterialDB::SetupMater(double InMr)
{
	if(InMr <= 0.) InMr = Mr;
	return new radTLinearAnisotropMaterial(KsiPar, KsiPerp, InMr);
}

//-------------------------------------------------------------------------

radTMaterial* radTNonlinearIsotropMaterialDB::SetupMater(double InMr)
{
	return new radTNonlinearIsotropMaterial(ms, ks, 3);
}

//-------------------------------------------------------------------------
