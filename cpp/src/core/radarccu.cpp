/*-------------------------------------------------------------------------
*
* File name:      radarccu.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 rectangular cross-section arc with azimuthal current
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radappl.h"
#include "radarccu.h"
#include "radg3dgr.h"
#include "radsbdac.h"

#include <math.h>
#include <sstream>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTArcCur::B_compWithTrapeth(radTField* FieldPtr)
{
	const double Pi = 3.141592653589793238;
	const double ZeroTol = 1.E-10;
	const double RelZeroTolToDecompose = 1.E-07;

	TVector3d P_mi_CenPo = FieldPtr->P - CircleCentrPoint;

	const double SmallPositive = 1.E-10;
	double r = sqrt(P_mi_CenPo.x*P_mi_CenPo.x + P_mi_CenPo.y*P_mi_CenPo.y + SmallPositive);
	double phi = ((P_mi_CenPo.y < 0)? (2*Pi-acos(P_mi_CenPo.x/r)) : (acos(P_mi_CenPo.x/r)));
	double z1 = -Height/2. - P_mi_CenPo.z;
	double z2 = z1 + Height;

	double z1z1 = z1*z1;
	double z2z2 = z2*z2;
	double rr = r*r;
	double r1r1 = R_min*R_min;
	double r2r2 = R_max*R_max;
	double rr_pl_r1r1_pl_z1z1 = rr + r1r1 + z1z1;
	double rr_pl_r1r1_pl_z2z2 = rr + r1r1 + z2z2;
	double rr_pl_r2r2_pl_z1z1 = rr + r2r2 + z1z1;
	double rr_pl_r2r2_pl_z2z2 = rr + r2r2 + z2z2;
	double two_r = 2.*r;
	double two_rr1 = two_r*R_min;
	double two_rr2 = two_r*R_max;

	double CosDelPhi, AbsSinDelPhi, D11, D12, D21, D22, rCos, r1_mi_rCos, r2_mi_rCos,
		   ln_r1mrCospD11_di_r2mrCospD21, ln_r2mrCospD22_di_r1mrCospD12, rAbsSin, 
		   Absr1mrCos_di_rAbsSin, Absr2mrCos_di_rAbsSin, PiMult1, PiMult2, Buf,
		   z1pD11_di_z2pD12, z2pD22_di_z1pD21, ln_z1pD11diz2pD12, ln_z2pD22diz1pD21,
		   phiPrime, Cos_phiPrime, Sin_phiPrime, rr_mi_two_rrCosCos, FA, FBxy, FBz;

	int AmOfPoi = NumberOfSectors + 1;
	double Step_phiPrime = (Phi_max - Phi_min + SmallPositive)/(AmOfPoi-1);
	double phiPrime_mi_phi = Phi_min - phi + SmallPositive;

	double S_forAx, S_forAy, S_forBx, S_forBy, S_forBz;
	S_forAx = S_forAy = S_forBx = S_forBy = S_forBz = 0.;

	int AmOfPoi_mi_1 = AmOfPoi - 1;
	for(int i=0; i<AmOfPoi; i++)
	{
		CosDelPhi = cos(phiPrime_mi_phi);
		AbsSinDelPhi = Abs(sin(phiPrime_mi_phi));

		double two_rr1_CosDelPhi = two_rr1*CosDelPhi;
		double two_rr2_CosDelPhi = two_rr2*CosDelPhi;

		double BufVal11 = fabs(rr_pl_r1r1_pl_z1z1 - two_rr1_CosDelPhi); if(BufVal11<ZeroTol) BufVal11=ZeroTol;
		D11 = sqrt(BufVal11);
		D12 = sqrt(fabs(rr_pl_r1r1_pl_z2z2 - two_rr1*CosDelPhi));
		double BufVal21 = fabs(rr_pl_r2r2_pl_z1z1 - two_rr2_CosDelPhi);	if(BufVal21<ZeroTol) BufVal21=ZeroTol;
		D21 = sqrt(BufVal21);
		D22 = sqrt(fabs(rr_pl_r2r2_pl_z2z2 - two_rr2*CosDelPhi));
		rCos = r*CosDelPhi;
		r1_mi_rCos = R_min - rCos;
		r2_mi_rCos = R_max - rCos;

		double D11_RelZeroTolToDecomp = D11*RelZeroTolToDecompose; //OC190504
		double D12_RelZeroTolToDecomp = D12*RelZeroTolToDecompose;
		double D21_RelZeroTolToDecomp = D21*RelZeroTolToDecompose;
		double D22_RelZeroTolToDecomp = D22*RelZeroTolToDecompose;

		double re2SinE2 = rr*AbsSinDelPhi*AbsSinDelPhi; //OC180504
		double z1z1pre2SinE2 = z1z1 + re2SinE2;

		double r1_mi_rCosE2 = r1_mi_rCos*r1_mi_rCos; //OC180504
		double r1_mi_rCos_p_D11 = r1_mi_rCos + D11;
		if((Abs(r1_mi_rCos_p_D11) < D11_RelZeroTolToDecomp) && (z1z1pre2SinE2 < r1_mi_rCosE2*RelZeroTolToDecompose))
		{
			r1_mi_rCos_p_D11 = 0.5*z1z1pre2SinE2/Abs(r1_mi_rCos);
		}

		double r1_mi_rCos_p_D12 = r1_mi_rCos + D12; //OC180504
		if((Abs(r1_mi_rCos_p_D12) < D12_RelZeroTolToDecomp) && (z1z1pre2SinE2 < r1_mi_rCosE2*RelZeroTolToDecompose))
		{
			r1_mi_rCos_p_D12 = 0.5*z1z1pre2SinE2/Abs(r1_mi_rCos);
		}

		double r2_mi_rCosE2 = r2_mi_rCos*r2_mi_rCos; //OC180504
		double r2_mi_rCos_p_D21 = r2_mi_rCos + D21;
		if((Abs(r2_mi_rCos_p_D21) < D21_RelZeroTolToDecomp) && (z1z1pre2SinE2 < r2_mi_rCosE2*RelZeroTolToDecompose))
		{
			r2_mi_rCos_p_D21 = 0.5*z1z1pre2SinE2/Abs(r2_mi_rCos);
		}
		double r2_mi_rCos_p_D22 = r2_mi_rCos + D22; //OC180504
		if((Abs(r2_mi_rCos_p_D22) < D22_RelZeroTolToDecomp) && (z1z1pre2SinE2 < r2_mi_rCosE2*RelZeroTolToDecompose))
		{
			r2_mi_rCos_p_D22 = 0.5*z1z1pre2SinE2/Abs(r2_mi_rCos);
		}

		ln_r1mrCospD11_di_r2mrCospD21 = log((r1_mi_rCos_p_D11)/(r2_mi_rCos_p_D21));
		ln_r2mrCospD22_di_r1mrCospD12 = log((r2_mi_rCos_p_D22)/(r1_mi_rCos_p_D12));

		rAbsSin = r*AbsSinDelPhi;
		Absr1mrCos_di_rAbsSin = Abs(r1_mi_rCos)/rAbsSin;
		Absr2mrCos_di_rAbsSin = Abs(r2_mi_rCos)/rAbsSin;
	
		PiMult1 = PiMult2 = 0.;
		Buf = rAbsSin*(Sign(r1_mi_rCos)*(atan(TransAtans(z2*Absr1mrCos_di_rAbsSin/D12, -z1*Absr1mrCos_di_rAbsSin/D11, PiMult1)) + Pi*PiMult1)
					   +Sign(r2_mi_rCos)*(atan(TransAtans(z1*Absr2mrCos_di_rAbsSin/D21, -z2*Absr2mrCos_di_rAbsSin/D22, PiMult2)) + Pi*PiMult2))
			  +z1*ln_r1mrCospD11_di_r2mrCospD21 + z2*ln_r2mrCospD22_di_r1mrCospD12;

		double z1_p_D11 = z1 + D11; //OC190504
		double rr_p_r1r1_mi_2rr1Cos = rr + r1r1 - two_rr1_CosDelPhi;

		double z1z1_RelZeroTolToDecompose = z1z1*RelZeroTolToDecompose;
		double z2z2_RelZeroTolToDecompose = z2z2*RelZeroTolToDecompose;

		if((Abs(z1_p_D11) < D11_RelZeroTolToDecomp) && (rr_p_r1r1_mi_2rr1Cos < z1z1_RelZeroTolToDecompose))
		{
			z1_p_D11 = 0.5*rr_p_r1r1_mi_2rr1Cos/Abs(z1);
		}

		double z1_p_D21 = z1 + D21; //OC190504
		double rr_p_r2r2_mi_2rr2Cos = rr + r2r2 - two_rr2_CosDelPhi;
		if((Abs(z1_p_D21) < D21_RelZeroTolToDecomp) && (rr_p_r2r2_mi_2rr2Cos < z1z1_RelZeroTolToDecompose))
		{
            z1_p_D21 = 0.5*rr_p_r2r2_mi_2rr2Cos/Abs(z1);
		}

		double z2_p_D12 = z2 + D12; //OC190504
		if((Abs(z2_p_D12) < D12_RelZeroTolToDecomp) && (rr_p_r1r1_mi_2rr1Cos < z2z2_RelZeroTolToDecompose))
		{
			z2_p_D12 = 0.5*rr_p_r1r1_mi_2rr1Cos/Abs(z2);
		}

		double z2_p_D22 = z2 + D22; //OC190504
		if((Abs(z2_p_D22) < D22_RelZeroTolToDecomp) && (rr_p_r2r2_mi_2rr2Cos < z2z2_RelZeroTolToDecompose))
		{
			z2_p_D22 = 0.5*rr_p_r2r2_mi_2rr2Cos/Abs(z2);
		}

		z1pD11_di_z2pD12 = z1_p_D11/z2_p_D12;
		z2pD22_di_z1pD21 = z2_p_D22/z1_p_D21;

		phiPrime = phiPrime_mi_phi+phi;
		Cos_phiPrime = cos(phiPrime);
		Sin_phiPrime = sin(phiPrime);

		if((i==0) || (i==AmOfPoi_mi_1))
		{
			if(FieldPtr->FieldKey.A_)
			{
				ln_z1pD11diz2pD12 = log(z1pD11_di_z2pD12);
				ln_z2pD22diz1pD21 = log(z2pD22_di_z1pD21);
				rr_mi_two_rrCosCos = rr*(1.-2.*CosDelPhi*CosDelPhi);
				FA = rCos*Buf + 0.5*(z1*(D11-D21)+z2*(D22-D12)+(rr_mi_two_rrCosCos+r1r1)*ln_z1pD11diz2pD12+(rr_mi_two_rrCosCos+r2r2)*ln_z2pD22diz1pD21);
				S_forAx += 0.5*(-Sin_phiPrime*FA);
				S_forAy += 0.5*Cos_phiPrime*FA;
			}
			if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
			{
				FBxy = D11-D21+D22-D12 + rCos*(ln_r1mrCospD11_di_r2mrCospD21+ln_r2mrCospD22_di_r1mrCospD12);
				if(FieldPtr->FieldKey.A_) FBz = Buf - rCos*(ln_z1pD11diz2pD12+ln_z2pD22diz1pD21);
				else FBz = Buf - rCos*log(z1pD11_di_z2pD12*z2pD22_di_z1pD21);
				S_forBx += 0.5*Cos_phiPrime*FBxy;
				S_forBy += 0.5*Sin_phiPrime*FBxy;
				S_forBz += 0.5*FBz;
			}
		}
		else
		{
			if(FieldPtr->FieldKey.A_)
			{
				ln_z1pD11diz2pD12 = log(z1pD11_di_z2pD12);
				ln_z2pD22diz1pD21 = log(z2pD22_di_z1pD21);
				rr_mi_two_rrCosCos = rr*(1.-2.*CosDelPhi*CosDelPhi);
				FA = rCos*Buf + 0.5*(z1*(D11-D21)+z2*(D22-D12)+(rr_mi_two_rrCosCos+r1r1)*ln_z1pD11diz2pD12+(rr_mi_two_rrCosCos+r2r2)*ln_z2pD22diz1pD21);
				S_forAx += -Sin_phiPrime*FA;
				S_forAy += Cos_phiPrime*FA;
			}
			if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
			{
				FBxy = D11-D21+D22-D12 + rCos*(ln_r1mrCospD11_di_r2mrCospD21+ln_r2mrCospD22_di_r1mrCospD12);
				if(FieldPtr->FieldKey.A_) FBz = Buf - rCos*(ln_z1pD11diz2pD12+ln_z2pD22diz1pD21);
				else FBz = Buf - rCos*log(z1pD11_di_z2pD12*z2pD22_di_z1pD21);
				S_forBx += Cos_phiPrime*FBxy;
				S_forBy += Sin_phiPrime*FBxy;
				S_forBz += FBz;
			}
		}
		phiPrime_mi_phi += Step_phiPrime;
	}

	const double ConstForJ = 0.0001;
	double Fact = ConstForJ*J_azim;
	if(FieldPtr->FieldKey.A_)
	{
		double IntForAx = Step_phiPrime*S_forAx; 
		double IntForAy = Step_phiPrime*S_forAy;

		TVector3d BufA(IntForAx, IntForAy, 0.);
		FieldPtr->A += Fact*BufA;
	}
	if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
	{
		double IntForBx = Step_phiPrime*S_forBx; 
		double IntForBy = Step_phiPrime*S_forBy; 
		double IntForBz = Step_phiPrime*S_forBz;

		TVector3d BufB(IntForBx, IntForBy, IntForBz);
		BufB = Fact*BufB;
		FieldPtr->B += BufB;
		FieldPtr->H += BufB;
	}
}

//-------------------------------------------------------------------------

void radTArcCur::B_compWithNewtonCotes4(radTField* FieldPtr)
{
	const double Pi = 3.141592653589793238;
	const double ConstForJ = 0.0001;
	double Fact = ConstForJ*J_azim;

	TVector3d P_mi_CenPo = FieldPtr->P - CircleCentrPoint;

	const double SmallPositive = 1.E-10;
	double r = sqrt(P_mi_CenPo.x*P_mi_CenPo.x + P_mi_CenPo.y*P_mi_CenPo.y + SmallPositive);
	double phi = ((P_mi_CenPo.y < 0)? (2*Pi-acos(P_mi_CenPo.x/r)) : (acos(P_mi_CenPo.x/r)));
	double z1 = -Height/2. - P_mi_CenPo.z;
	double z2 = z1 + Height;

	double z1z1 = z1*z1;
	double z2z2 = z2*z2;
	double rr = r*r;
	double r1r1 = R_min*R_min;
	double r2r2 = R_max*R_max;
	double rr_pl_r1r1_pl_z1z1 = rr + r1r1 + z1z1;
	double rr_pl_r1r1_pl_z2z2 = rr + r1r1 + z2z2;
	double rr_pl_r2r2_pl_z1z1 = rr + r2r2 + z1z1;
	double rr_pl_r2r2_pl_z2z2 = rr + r2r2 + z2z2;
	double two_r = 2.*r;
	double two_rr1 = two_r*R_min;
	double two_rr2 = two_r*R_max;

	double phiPrime_mi_phi_min = Phi_min - phi + SmallPositive;

	const double IntegWeight[] = {14./45., 64./45., 24./45., 64./45., 28./45.};
	double IntForAx, IntForAy, IntForBx, IntForBy, IntForBz,
		   PrIntForAx, PrIntForAy, PrIntForBx, PrIntForBy, PrIntForBz,
		   GenS_forAx, GenS_forAy, GenS_forBx, GenS_forBy, GenS_forBz,
		   S_forAx, S_forAy, S_forBx, S_forBy, S_forBz;
    IntForAx = IntForAy = IntForBx = IntForBy = IntForBz = 1.E+23;
	GenS_forAx = GenS_forAy = GenS_forBx = GenS_forBy = GenS_forBz = 0.;

	double Step_phiPrime, phiPrime_mi_phi;
	short IndForWeight, IndForPass;

	double CosDelPhi, AbsSinDelPhi, D11, D12, D21, D22, rCos, r1_mi_rCos, r2_mi_rCos,
		   ln_r1mrCospD11_di_r2mrCospD21, ln_r2mrCospD22_di_r1mrCospD12, rAbsSin, 
		   Absr1mrCos_di_rAbsSin, Absr2mrCos_di_rAbsSin, PiMult1, PiMult2, Buf,
		   z1pD11_di_z2pD12, z2pD22_di_z1pD21, ln_z1pD11diz2pD12, ln_z2pD22diz1pD21,
		   phiPrime, Cos_phiPrime, Sin_phiPrime, rr_mi_two_rrCosCos, FA, FBxy, FBz;

	int AmOfPoi = 5;
	short NotFirstPass = 0;
	double PrecParamB = 1.E+23; 
	double PrecParamA = 1.E+23; 

	while((PrecParamB > FieldPtr->CompCriterium.AbsPrecB) || 
		  (PrecParamA > FieldPtr->CompCriterium.AbsPrecA))
	{
		Step_phiPrime = (Phi_max - Phi_min + SmallPositive)/(AmOfPoi-1);
		phiPrime_mi_phi = phiPrime_mi_phi_min;

		PrIntForAx = IntForAx; PrIntForAy = IntForAy;
		PrIntForBx = IntForBx; PrIntForBy = IntForBy; PrIntForBz = IntForBz;

		IndForWeight = IndForPass = 0;
		S_forAx = S_forAy = S_forBx = S_forBy = S_forBz = 0.;

		int AmOfPoi_mi_1 = AmOfPoi - 1;
		for(int i=0; i<AmOfPoi; i++)
		{
			if(IndForPass==3) IndForPass = 0;
			if(IndForWeight==5) IndForWeight = 1;
			if(NotFirstPass && (IndForPass==0)) goto BottomOfThisLoop;
			if(i==AmOfPoi_mi_1) IndForWeight = 0;

			CosDelPhi = cos(phiPrime_mi_phi);
			AbsSinDelPhi = Abs(sin(phiPrime_mi_phi));
			D11 = sqrt(rr_pl_r1r1_pl_z1z1 - two_rr1*CosDelPhi);
			D12 = sqrt(rr_pl_r1r1_pl_z2z2 - two_rr1*CosDelPhi);
			D21 = sqrt(rr_pl_r2r2_pl_z1z1 - two_rr2*CosDelPhi);
			D22 = sqrt(rr_pl_r2r2_pl_z2z2 - two_rr2*CosDelPhi);
			rCos = r*CosDelPhi;
			r1_mi_rCos = R_min - rCos;
			r2_mi_rCos = R_max - rCos;
			ln_r1mrCospD11_di_r2mrCospD21 = log((r1_mi_rCos + D11)/(r2_mi_rCos + D21));
			ln_r2mrCospD22_di_r1mrCospD12 = log((r2_mi_rCos + D22)/(r1_mi_rCos + D12));
			rAbsSin = r*AbsSinDelPhi;
			Absr1mrCos_di_rAbsSin = Abs(r1_mi_rCos)/rAbsSin;
			Absr2mrCos_di_rAbsSin = Abs(r2_mi_rCos)/rAbsSin;

			PiMult1 = PiMult2 = 0.;
			Buf = rAbsSin*(Sign(r1_mi_rCos)*(atan(TransAtans(z2*Absr1mrCos_di_rAbsSin/D12, -z1*Absr1mrCos_di_rAbsSin/D11, PiMult1)) + Pi*PiMult1)
						   +Sign(r2_mi_rCos)*(atan(TransAtans(z1*Absr2mrCos_di_rAbsSin/D21, -z2*Absr2mrCos_di_rAbsSin/D22, PiMult2)) + Pi*PiMult2))
				  +z1*ln_r1mrCospD11_di_r2mrCospD21 + z2*ln_r2mrCospD22_di_r1mrCospD12;

			z1pD11_di_z2pD12 = (z1 + D11)/(z2 + D12);
			z2pD22_di_z1pD21 = (z2 + D22)/(z1 + D21);
			phiPrime = phiPrime_mi_phi+phi;
			Cos_phiPrime = cos(phiPrime);
			Sin_phiPrime = sin(phiPrime);
			if(FieldPtr->FieldKey.A_)
			{
				ln_z1pD11diz2pD12 = log(z1pD11_di_z2pD12);
				ln_z2pD22diz1pD21 = log(z2pD22_di_z1pD21);
				rr_mi_two_rrCosCos = rr*(1.-2.*CosDelPhi*CosDelPhi);
				FA = rCos*Buf + 0.5*(z1*(D11-D21)+z2*(D22-D12)+(rr_mi_two_rrCosCos+r1r1)*ln_z1pD11diz2pD12+(rr_mi_two_rrCosCos+r2r2)*ln_z2pD22diz1pD21);
				S_forAx += IntegWeight[IndForWeight]*(-Sin_phiPrime*FA);
				S_forAy += IntegWeight[IndForWeight]*Cos_phiPrime*FA;
			}
			if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
			{
				FBxy = D11-D21+D22-D12 + rCos*(ln_r1mrCospD11_di_r2mrCospD21+ln_r2mrCospD22_di_r1mrCospD12);
				if(FieldPtr->FieldKey.A_) FBz = Buf - rCos*(ln_z1pD11diz2pD12+ln_z2pD22diz1pD21);
				else FBz = Buf - rCos*log(z1pD11_di_z2pD12*z2pD22_di_z1pD21);
				S_forBx += IntegWeight[IndForWeight]*Cos_phiPrime*FBxy;
				S_forBy += IntegWeight[IndForWeight]*Sin_phiPrime*FBxy;
				S_forBz += IntegWeight[IndForWeight]*FBz;
			}
BottomOfThisLoop:
			IndForPass++; IndForWeight++;
			phiPrime_mi_phi += Step_phiPrime;
		}

		PrecParamA = PrecParamB = 0.;

		if(FieldPtr->FieldKey.A_)
		{
			GenS_forAx += S_forAx; GenS_forAy += S_forAy;
			IntForAx = Step_phiPrime*GenS_forAx; IntForAy = Step_phiPrime*GenS_forAy;

			PrecParamA = Fact * Max(Abs(IntForAx-PrIntForAx), Abs(IntForAy-PrIntForAy));
		}
		if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
		{
			GenS_forBx += S_forBx; GenS_forBy += S_forBy; GenS_forBz += S_forBz;
			IntForBx = Step_phiPrime*GenS_forBx; 
			IntForBy = Step_phiPrime*GenS_forBy; 
			IntForBz = Step_phiPrime*GenS_forBz;

			PrecParamB = Fact * Max( Max( Abs(IntForBx-PrIntForBx), Abs(IntForBy-PrIntForBy)), Abs(IntForBz-PrIntForBz));
		}
		AmOfPoi = (AmOfPoi-1)*3+1;
		NotFirstPass = 1;
	}
	
	if(FieldPtr->FieldKey.A_)
	{
		TVector3d BufA(IntForAx, IntForAy, 0.);
		FieldPtr->A += Fact*BufA;
	}
	if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
	{
		TVector3d BufB(IntForBx, IntForBy, IntForBz);
		BufB = Fact*BufB;
		FieldPtr->B += BufB;
		FieldPtr->H += BufB;
	}
}

//-------------------------------------------------------------------------

void radTArcCur::B_intUtilSpecCaseZeroVxVy(double r0, double r1, double r2, double ph0, double ph1, double ph2, double h, double& Iz)
{
	const double Pi = 3.141592653589793238;
	const double SmallPositive = 1.E-10;

	double ph0_mi_ph1 = ph0-ph1;
	double ph0_mi_ph2 = ph0-ph2;
	double a1 = tan(0.5*ph0_mi_ph1);
	double a2 = tan(0.5*ph0_mi_ph2);
	double a1a1 = a1*a1;
	double a2a2 = a2*a2;
	double a1a1_pl_1 = a1a1+1.;
	double a2a2_pl_1 = a2a2+1.;
	double a1a1_mi_1 = a1a1-1.;
	double a2a2_mi_1 = a2a2-1.;
	double r1_mi_r2 = r1-r2;

CorrectionStart:
	double r0_mi_r1 = r0-r1; if(r0_mi_r1==0.) { r0+=SmallPositive; goto CorrectionStart;}
	double r0_mi_r2 = r0-r2; if(r0_mi_r2==0.) { r0+=SmallPositive; goto CorrectionStart;}
	double r1_pl_r2 = r1+r2;
	double r0_pl_r1 = r0+r1;
	double r0_pl_r2 = r0+r2;
	double r0plr1_di_r0mir1 = r0_pl_r1/r0_mi_r1;
	double r0plr2_di_r0mir2 = r0_pl_r2/r0_mi_r2;
	double Two_a2r0 = 2.*a2*r0;
	double Two_a1r0 = 2.*a1*r0;
	double a2a2mi1_mu_r0 = a2a2_mi_1*r0;
	double a1a1mi1_mu_r0 = a1a1_mi_1*r0;
	double a2a2mi1_r0_pl_a2a2pl1_r1 = a2a2mi1_mu_r0+a2a2_pl_1*r1;
	double a2a2mi1_r0_pl_a2a2pl1_r2 = a2a2mi1_mu_r0+a2a2_pl_1*r2; 
	double a1a1mi1_r0_pl_a1a1pl1_r1 = a1a1mi1_mu_r0+a1a1_pl_1*r1;
	double a1a1mi1_r0_pl_a1a1pl1_r2 = a1a1mi1_mu_r0+a1a1_pl_1*r2;
	if(a2a2mi1_r0_pl_a2a2pl1_r1==0. || 
	   a2a2mi1_r0_pl_a2a2pl1_r2==0. ||
	   a1a1mi1_r0_pl_a1a1pl1_r1==0. ||
	   a1a1mi1_r0_pl_a1a1pl1_r2==0.) { r0+=SmallPositive; goto CorrectionStart;}

	double a2a2mi1_r0_di_a2a2pl1 = a2a2mi1_mu_r0/a2a2_pl_1;
	double a1a1mi1_r0_di_a1a1pl1 = a1a1mi1_mu_r0/a1a1_pl_1;

	double PiMult1, PiMult2, PiMult3, PiMult4, Sz;
	PiMult1 = PiMult2 = PiMult3 = PiMult4 = 0.;

	double SzL1 = 0.5*(ph2-ph1)*(r2-r1);
	double SzL2 = Pi*(Step(ph0_mi_ph1-Pi)*Step(Pi-ph0_mi_ph2) + Step(Pi+ph0_mi_ph1)*Step(-Pi-ph0_mi_ph2))*(r1_mi_r2*(Step(r0_mi_r2)*Step(r0_mi_r1)-Step(-r0_mi_r1)*Step(-r0_mi_r2)) + (r1_pl_r2-2.*r0)*Step(r0_mi_r1)*Step(-r0_mi_r2));
	double SzL3 = Pi*r0*Step(r0_mi_r1)*Step(-r0_mi_r2)*(Sign(a2)-Sign(a1));
	double SzL4 = r1*(atan(TransAtans(a1*r0plr1_di_r0mir1, -a2*r0plr1_di_r0mir1, PiMult1)) + Pi*PiMult1);
	double SzL5 = r2*(atan(TransAtans(a2*r0plr2_di_r0mir2, -a1*r0plr2_di_r0mir2, PiMult2)) + Pi*PiMult2);
	double SumAtan1 = atan(TransAtans(Two_a2r0/a2a2mi1_r0_pl_a2a2pl1_r1, -Two_a2r0/a2a2mi1_r0_pl_a2a2pl1_r2, PiMult3)) + Pi*PiMult3;
	double PiTerm1 = Pi*Sign(a2)*Step(-a2a2mi1_r0_di_a2a2pl1-r1)*Step(r2+a2a2mi1_r0_di_a2a2pl1);
	double SzL6 = a2a2mi1_r0_di_a2a2pl1*(SumAtan1 + PiTerm1);
	double SumAtan2 = atan(TransAtans(Two_a1r0/a1a1mi1_r0_pl_a1a1pl1_r1, -Two_a1r0/a1a1mi1_r0_pl_a1a1pl1_r2, PiMult4)) + Pi*PiMult4;
	double PiTerm2 = Pi*Sign(a1)*Step(-a1a1mi1_r0_di_a1a1pl1-r1)*Step(r2+a1a1mi1_r0_di_a1a1pl1);
	double SzL7 = -a1a1mi1_r0_di_a1a1pl1*(SumAtan2 + PiTerm2);
	double SzL8 = (a2*r0/a2a2_pl_1)*log((r0_mi_r1*r0_mi_r1 + a2a2*r0_pl_r1*r0_pl_r1)/(r0_mi_r2*r0_mi_r2 + a2a2*r0_pl_r2*r0_pl_r2));
	double SzL9 = (a1*r0/a1a1_pl_1)*log((r0_mi_r2*r0_mi_r2 + a1a1*r0_pl_r2*r0_pl_r2)/(r0_mi_r1*r0_mi_r1 + a1a1*r0_pl_r1*r0_pl_r1));
	Sz = SzL1 + SzL2 + SzL3 + SzL4 + SzL5 + SzL6 + SzL7 + SzL8 + SzL9;
	Iz = 2.*h*Sz;
}

//-------------------------------------------------------------------------

void radTArcCur::B_intCompWithTrapeth(radTField* FieldPtr)
{
	const double Pi = 3.141592653589793238;
	const double ConForJ = 1.E-04;
	const double SmallPositive = 1.E-12;
	const double VxVyCaseZeroToler = 1.E-05; // SmallPositive < VxVyCaseZeroToler !!!
	double Fact = ConForJ*J_azim;

	TVector3d VectV = FieldPtr->NextP - FieldPtr->P;
	double ModV = sqrt(VectV.x*VectV.x + VectV.y*VectV.y + VectV.z*VectV.z);
	VectV.x = VectV.x/ModV; VectV.y = VectV.y/ModV; VectV.z = VectV.z/ModV;
	if(VectV.x==0. || VectV.x==-1.) VectV.x += SmallPositive;
	else if(VectV.x==1.) VectV.x -= SmallPositive;
	if(VectV.y==0. || VectV.y==-1.) VectV.y += SmallPositive;
	else if(VectV.y==1.) VectV.y -= SmallPositive;

	TVector3d IntForB(0.,0.,0.);

	TVector3d P_mi_CPoi = FieldPtr->P - CircleCentrPoint;
	double r0 = sqrt(P_mi_CPoi.x*P_mi_CPoi.x + P_mi_CPoi.y*P_mi_CPoi.y); if(r0==0.) r0 = SmallPositive;
	double PmiCPoix_di_r0 = P_mi_CPoi.x/r0;
	double ph0 = (P_mi_CPoi.y<0)? 2.*Pi-acos(PmiCPoix_di_r0) : acos(PmiCPoix_di_r0);
	double z1 = -0.5*Height - P_mi_CPoi.z;
	double z2 = z1 + Height;
	double r1 = R_min;
	double r2 = R_max;

// Check for Special Case
	if(Abs(VectV.x)<VxVyCaseZeroToler && Abs(VectV.y)<VxVyCaseZeroToler) 
	{
		B_intUtilSpecCaseZeroVxVy(r0, r1, r2, ph0, Phi_min, Phi_max, Height, IntForB.z);
		IntForB.x = IntForB.y = 0.;
		goto FinalDefinitionOfFieldIntegrals;
	}

// General Case: numerical integration over Phi. This uses Newton method (n=3).
	{
		double Sinph0 = sin(ph0);
		double Cosph0 = cos(ph0);
		double Vph0x = Sinph0*VectV.y + Cosph0*VectV.x;
		double Vph0y = Cosph0*VectV.y - Sinph0*VectV.x;
		double VzVz = VectV.z*VectV.z;
		double VxVx = VectV.x*VectV.x;
		double VyVy = VectV.y*VectV.y;
		double VxVxpVyVy = VxVx+VyVy;

			if(VxVxpVyVy==0.) VxVxpVyVy = SmallPositive;

		double Vzz1 = VectV.z*z1;
		double Vzz2 = VectV.z*z2;
		double VzVzr0 = VzVz*r0;
		double r0r0 = r0*r0;
		double r1r1 = r1*r1;
		double r2r2 = r2*r2;
		double r2mir1 = r2-r1;
		double z2miz1 = z2-z1;
		double Vph0xr0 = Vph0x*r0;
		double Vph0xr0_mi_Vzz1 = Vph0xr0-Vzz1;
		double Vph0xr0_mi_Vzz2 = Vph0xr0-Vzz2;
		double C1 = z1*z1+r0r0 - Vph0xr0_mi_Vzz1*Vph0xr0_mi_Vzz1;
		double C2 = z2*z2+r0r0 - Vph0xr0_mi_Vzz2*Vph0xr0_mi_Vzz2;
		double T = -Vph0y*r0;
		double R1 = z1 + VectV.z*Vph0xr0_mi_Vzz1;
		double R2 = z2 + VectV.z*Vph0xr0_mi_Vzz2;

		double ph1 = Phi_min + SmallPositive;
		double ph2 = Phi_max - SmallPositive;
		double PhMax_mi_PhMin = ph2 - ph1;

		double Sinph, Cosph, Cosphmiph0, Vphx, Vphy, VphyT, R1Vphy, R2Vphy, Vphyr1pT, Vphyr2pT, 
			   A, PreB, B1, B2, B1B1, B2B2, P, PT, K;
		double ArgLnZ1, ArgLnZ2, T_di_Vphy, Pr1, Pr2, PR1_pl_VphyT, PR2_pl_VphyT, R1Vphy_mi_PT, R2Vphy_mi_PT,
			   Kr1, Kr2, Kr1_pl_PR1_pl_VphyT, Kr2_pl_PR1_pl_VphyT, Kr1_pl_PR2_pl_VphyT, Kr2_pl_PR2_pl_VphyT, 
			   PR1pVphyT_di_K, PR2pVphyT_di_K, One_di_AA, One_di_KK, Radical1, Radical2, TwoAr1, TwoAr2, 
			   AC1, AC2, Ar1r1_pl_B1r1_pl_C1, Ar2r2_pl_B1r2_pl_C1, Ar1r1_pl_B2r1_pl_C2, Ar2r2_pl_B2r2_pl_C2,
			   Ar1r1, Ar2r2;
		double PiMult1, PiMult2, PiMult3, PiMult4;
		PiMult1 = PiMult2 = PiMult3 = PiMult4 = 0.;
		double Wz, WzL1, WzL2, WzL3, WzL4, WzL5, WzL6, Uz, UzL1, UzL2, UzL3, UzL4;

		int AmOfPoi = NumberOfSectors + 1;
		int AmOfPoi_mi_1 = NumberOfSectors;
		double Step_ph = PhMax_mi_PhMin/AmOfPoi_mi_1;
		double ph = ph1;

		TVector3d ZeroVect(0.,0.,0.), S_forIntB(0.,0.,0.), Func(0.,0.,0.);

		for(int i=0; i<AmOfPoi; i++)
		{
			Sinph = sin(ph); Cosph = cos(ph);
			Cosphmiph0 = cos(ph-ph0);
			Vphx = Sinph*VectV.y + Cosph*VectV.x;
			Vphy = Cosph*VectV.y - Sinph*VectV.x;
			R1Vphy = R1*Vphy;
			R2Vphy = R2*Vphy;
			Vphyr1pT = Vphy*r1+T;
			Vphyr2pT = Vphy*r2+T;
			A = 1.- Vphx*Vphx; 
			
//TO FIX: loss of accuracy if A is 0 or very small
//make special case? 

				if(A==0.) A = SmallPositive;

			One_di_AA = 1./(A*A);
			TwoAr1 = 2.*A*r1;
			TwoAr2 = 2.*A*r2;
			AC1 = A*C1;
			AC2 = A*C2;
			PreB = T*Vphy - Cosphmiph0*VzVzr0;
			B1 = 2.*(PreB - Vzz1*Vphx);
			B1B1 = B1*B1;
			B2 = 2.*(PreB - Vzz2*Vphx);
			B2B2 = B2*B2;
			P = -VectV.z*Vphx;
			PT = P*T;
			R1Vphy_mi_PT = R1Vphy-PT;
			R2Vphy_mi_PT = R2Vphy-PT;
			Pr1 = P*r1;
			Pr2 = P*r2;
			VphyT = Vphy*T;
			PR1_pl_VphyT = P*R1+VphyT;
			PR2_pl_VphyT = P*R2+VphyT;
			K = VxVxpVyVy*A;
			One_di_KK = 1./K/K;
			PR1pVphyT_di_K = PR1_pl_VphyT/K;
			PR2pVphyT_di_K = PR2_pl_VphyT/K;
			Kr1 = K*r1;
			Kr2 = K*r2;
			Kr1_pl_PR1_pl_VphyT = Kr1+PR1_pl_VphyT;
			Kr1_pl_PR2_pl_VphyT = Kr1+PR2_pl_VphyT;
			Kr2_pl_PR1_pl_VphyT = Kr2+PR1_pl_VphyT;
			Kr2_pl_PR2_pl_VphyT = Kr2+PR2_pl_VphyT;
			Ar1r1 = A*r1r1;
			Ar2r2 = A*r2r2;
			Ar1r1_pl_B1r1_pl_C1 = Ar1r1+B1*r1+C1;
			Ar2r2_pl_B1r2_pl_C1 = Ar2r2+B1*r2+C1;
			Ar1r1_pl_B2r1_pl_C2 = Ar1r1+B2*r1+C2;
			Ar2r2_pl_B2r2_pl_C2 = Ar2r2+B2*r2+C2;
			ArgLnZ1 = Ar1r1_pl_B1r1_pl_C1/Ar2r2_pl_B1r2_pl_C1;
			ArgLnZ2 = Ar1r1_pl_B2r1_pl_C2/Ar2r2_pl_B2r2_pl_C2;
			T_di_Vphy = T/Vphy;

			//Radical1 = sqrt(4.*A*C1-B1B1);
			//Radical2 = sqrt(4.*A*C2-B2B2);
			Radical1 = sqrt(::fabs(4.*A*C1-B1B1)); //OC fix 06/03: make a more accurate treatment
			Radical2 = sqrt(::fabs(4.*A*C2-B2B2));

			WzL1 = -(Vphy/A)*r2mir1*z2miz1;
			WzL2 = Pi*T_di_Vphy*T_di_Vphy*Step(-T_di_Vphy-r1)*Step(r2+T_di_Vphy)*(Sign(PT-R1Vphy)-Sign(PT-R2Vphy));
			WzL3 = r1r1*(atan(TransAtans((Pr1+R2)/Vphyr1pT, -(Pr1+R1)/Vphyr1pT, PiMult1)) + Pi*PiMult1);
			WzL4 = r2r2*(atan(TransAtans((Pr2+R1)/Vphyr2pT, -(Pr2+R2)/Vphyr2pT, PiMult2)) + Pi*PiMult2);
			WzL5 = One_di_KK*((PR1_pl_VphyT*PR1_pl_VphyT - R1Vphy_mi_PT*R1Vphy_mi_PT)*((atan(TransAtans(R1Vphy_mi_PT/Kr1_pl_PR1_pl_VphyT, -R1Vphy_mi_PT/Kr2_pl_PR1_pl_VphyT, PiMult3)) + Pi*PiMult3) + Pi*Sign(R1Vphy_mi_PT)*Step(-PR1pVphyT_di_K-r1)*Step(r2+PR1pVphyT_di_K)) 
				             +(PR2_pl_VphyT*PR2_pl_VphyT - R2Vphy_mi_PT*R2Vphy_mi_PT)*((atan(TransAtans(-R2Vphy_mi_PT/Kr1_pl_PR2_pl_VphyT, R2Vphy_mi_PT/Kr2_pl_PR2_pl_VphyT, PiMult4)) + Pi*PiMult4) - Pi*Sign(R2Vphy_mi_PT)*Step(-PR2pVphyT_di_K-r1)*Step(r2+PR2pVphyT_di_K)));
			WzL6 = One_di_KK*(PR1_pl_VphyT*R1Vphy_mi_PT*log(ArgLnZ1) - PR2_pl_VphyT*R2Vphy_mi_PT*log(ArgLnZ2));
			Wz = WzL1+WzL2+WzL3+WzL4+WzL5+WzL6;
			UzL1 = (B1-B2)*r2mir1/A;
			UzL2 = One_di_AA*(B1*Radical1*(atan(TransAtans((TwoAr1+B1)/Radical1, -(TwoAr2+B1)/Radical1, PiMult1)) + Pi*PiMult1)
							 +B2*Radical2*(atan(TransAtans(-(TwoAr1+B2)/Radical2, (TwoAr2+B2)/Radical2, PiMult2)) + Pi*PiMult2));
			UzL3 = 0.5*One_di_AA*(-(2.*AC1-B1B1)*log(ArgLnZ1) + (2.*AC2-B2B2)*log(ArgLnZ2));
			UzL4 = r1r1*log(Ar1r1_pl_B2r1_pl_C2/Ar1r1_pl_B1r1_pl_C1) + r2r2*log(Ar2r2_pl_B1r2_pl_C1/Ar2r2_pl_B2r2_pl_C2);
			Uz = UzL1+UzL2+UzL3+UzL4;

			Func.x = Cosph*Uz;
			Func.y = Sinph*Uz;
			Func.z = VectV.z*Vphx*Uz - 2.*Vphy*Wz;

			if((i==0) || (i==AmOfPoi_mi_1)) S_forIntB += 0.5*Func;
			else S_forIntB += Func;
			ph += Step_ph;
		}
		IntForB = (0.5*Step_ph)*S_forIntB;
	}

FinalDefinitionOfFieldIntegrals:
	TVector3d BufIb = Fact*IntForB;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIb;
	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIb;
}

//-------------------------------------------------------------------------

void radTArcCur::B_intCompWithNewton3(radTField* FieldPtr)
{
	const double Pi = 3.141592653589793238;
	const double ConForJ = 1.E-04;
	const double SmallPositive = 1.E-10;
	const double VxVyCaseZeroToler = 1.E-05; // SmallPositive < VxVyCaseZeroToler !!!
	double Fact = ConForJ*J_azim;

	TVector3d VectV = FieldPtr->NextP - FieldPtr->P;
	double ModV = sqrt(VectV.x*VectV.x + VectV.y*VectV.y + VectV.z*VectV.z);
	VectV.x = VectV.x/ModV; VectV.y = VectV.y/ModV; VectV.z = VectV.z/ModV;
	if(VectV.x==0. || VectV.x==-1.) VectV.x += SmallPositive;
	else if(VectV.x==1.) VectV.x -= SmallPositive;
	if(VectV.y==0. || VectV.y==-1.) VectV.y += SmallPositive;
	else if(VectV.y==1.) VectV.y -= SmallPositive;

	TVector3d IntForB(1.E+23, 1.E+23, 1.E+23);

	TVector3d P_mi_CPoi = FieldPtr->P - CircleCentrPoint;
	double r0 = sqrt(P_mi_CPoi.x*P_mi_CPoi.x + P_mi_CPoi.y*P_mi_CPoi.y); if(r0==0.) r0 = SmallPositive;
	double PmiCPoix_di_r0 = P_mi_CPoi.x/r0;
	double ph0 = (P_mi_CPoi.y<0)? 2.*Pi-acos(PmiCPoix_di_r0) : acos(PmiCPoix_di_r0);
	double z1 = -0.5*Height - P_mi_CPoi.z;
	double z2 = z1 + Height;
	double r1 = R_min;
	double r2 = R_max;

// Check for Special Case
	if(Abs(VectV.x)<VxVyCaseZeroToler && Abs(VectV.y)<VxVyCaseZeroToler) 
	{
		B_intUtilSpecCaseZeroVxVy(r0, r1, r2, ph0, Phi_min, Phi_max, Height, IntForB.z);
		IntForB.x = IntForB.y = 0.;
		goto FinalDefinitionOfFieldIntegrals;
	}

// General Case: numerical integration over Phi. This uses Newton method (n=3).
	{
		double Sinph0 = sin(ph0);
		double Cosph0 = cos(ph0);
		double Vph0x = Sinph0*VectV.y + Cosph0*VectV.x;
		double Vph0y = Cosph0*VectV.y - Sinph0*VectV.x;
		double VzVz = VectV.z*VectV.z;
		double VxVx = VectV.x*VectV.x;
		double VyVy = VectV.y*VectV.y;
		double VxVxpVyVy = VxVx+VyVy;
		double Vzz1 = VectV.z*z1;
		double Vzz2 = VectV.z*z2;
		double VzVzr0 = VzVz*r0;
		double r0r0 = r0*r0;
		double r1r1 = r1*r1;
		double r2r2 = r2*r2;
		double r2mir1 = r2-r1;
		double z2miz1 = z2-z1;
		double Vph0xr0 = Vph0x*r0;
		double Vph0xr0_mi_Vzz1 = Vph0xr0-Vzz1;
		double Vph0xr0_mi_Vzz2 = Vph0xr0-Vzz2;
		double C1 = z1*z1+r0r0 - Vph0xr0_mi_Vzz1*Vph0xr0_mi_Vzz1;
		double C2 = z2*z2+r0r0 - Vph0xr0_mi_Vzz2*Vph0xr0_mi_Vzz2;
		double T = -Vph0y*r0;
		double R1 = z1 + VectV.z*Vph0xr0_mi_Vzz1;
		double R2 = z2 + VectV.z*Vph0xr0_mi_Vzz2;

		double ph1 = Phi_min + SmallPositive;
		double ph2 = Phi_max - SmallPositive;
		double PhMax_mi_PhMin = ph2 - ph1;

		double Sinph, Cosph, Cosphmiph0, Vphx, Vphy, VphyT, R1Vphy, R2Vphy, Vphyr1pT, Vphyr2pT, 
			   A, PreB, B1, B2, B1B1, B2B2, P, PT, K;
		double ArgLnZ1, ArgLnZ2, T_di_Vphy, Pr1, Pr2, PR1_pl_VphyT, PR2_pl_VphyT, R1Vphy_mi_PT, R2Vphy_mi_PT,
			   Kr1, Kr2, Kr1_pl_PR1_pl_VphyT, Kr2_pl_PR1_pl_VphyT, Kr1_pl_PR2_pl_VphyT, Kr2_pl_PR2_pl_VphyT, 
			   PR1pVphyT_di_K, PR2pVphyT_di_K, One_di_AA, One_di_KK, Radical1, Radical2, TwoAr1, TwoAr2, 
			   AC1, AC2, Ar1r1_pl_B1r1_pl_C1, Ar2r2_pl_B1r2_pl_C1, Ar1r1_pl_B2r1_pl_C2, Ar2r2_pl_B2r2_pl_C2,
			   Ar1r1, Ar2r2;
		double PiMult1, PiMult2, PiMult3, PiMult4;
		PiMult1 = PiMult2 = PiMult3 = PiMult4 = 0.;
		double Wz, WzL1, WzL2, WzL3, WzL4, WzL5, WzL6, Uz, UzL1, UzL2, UzL3, UzL4;

		const double IntegWeight[] = {3./8., 9./8., 9./8., 3./4.};
		TVector3d ZeroVect(0.,0.,0.), S_forB, GenS_forB(0.,0.,0.), PrIntForB, Func(0.,0.,0.);

		double Step_ph, ph;
		short IndForWeight, IndForPass;
		short NotFirstPass = 0;

		int AmOfPoi = 4;
		int AmOfPoi_mi_1;
		double PrecParamInt = 1.E+23; 

		while(PrecParamInt > FieldPtr->CompCriterium.AbsPrecB_int)
		{
			AmOfPoi_mi_1 = AmOfPoi - 1;
			Step_ph = PhMax_mi_PhMin/AmOfPoi_mi_1;
			ph = ph1;

			PrIntForB = IntForB;
			S_forB = ZeroVect;
			IndForWeight = IndForPass = 0;

			for(int i=0; i<AmOfPoi; i++)
			{
				if(IndForPass==2) IndForPass = 0;
				if(IndForWeight==4) IndForWeight = 1;
				if(NotFirstPass && (IndForPass==0)) goto BottomOfThisLoop;
				if(i==AmOfPoi_mi_1) IndForWeight = 0;

				Sinph = sin(ph); Cosph = cos(ph);
				Cosphmiph0 = cos(ph-ph0);
				Vphx = Sinph*VectV.y + Cosph*VectV.x;
				Vphy = Cosph*VectV.y - Sinph*VectV.x;
				R1Vphy = R1*Vphy;
				R2Vphy = R2*Vphy;
				Vphyr1pT = Vphy*r1+T;
				Vphyr2pT = Vphy*r2+T;
				A = 1.- Vphx*Vphx;
				One_di_AA = 1./(A*A);
				TwoAr1 = 2.*A*r1;
				TwoAr2 = 2.*A*r2;
				AC1 = A*C1;
				AC2 = A*C2;
				PreB = T*Vphy - Cosphmiph0*VzVzr0;
				B1 = 2.*(PreB - Vzz1*Vphx);
				B1B1 = B1*B1;
				B2 = 2.*(PreB - Vzz2*Vphx);
				B2B2 = B2*B2;
				P = -VectV.z*Vphx;
				PT = P*T;
				R1Vphy_mi_PT = R1Vphy-PT;
				R2Vphy_mi_PT = R2Vphy-PT;
				Pr1 = P*r1;
				Pr2 = P*r2;
				VphyT = Vphy*T;
				PR1_pl_VphyT = P*R1+VphyT;
				PR2_pl_VphyT = P*R2+VphyT;
				K = VxVxpVyVy*A;
				One_di_KK = 1./K/K;
				PR1pVphyT_di_K = PR1_pl_VphyT/K;
				PR2pVphyT_di_K = PR2_pl_VphyT/K;
				Kr1 = K*r1;
				Kr2 = K*r2;
				Kr1_pl_PR1_pl_VphyT = Kr1+PR1_pl_VphyT;
				Kr1_pl_PR2_pl_VphyT = Kr1+PR2_pl_VphyT;
				Kr2_pl_PR1_pl_VphyT = Kr2+PR1_pl_VphyT;
				Kr2_pl_PR2_pl_VphyT = Kr2+PR2_pl_VphyT;
				Ar1r1 = A*r1r1;
				Ar2r2 = A*r2r2;
				Ar1r1_pl_B1r1_pl_C1 = Ar1r1+B1*r1+C1;
				Ar2r2_pl_B1r2_pl_C1 = Ar2r2+B1*r2+C1;
				Ar1r1_pl_B2r1_pl_C2 = Ar1r1+B2*r1+C2;
				Ar2r2_pl_B2r2_pl_C2 = Ar2r2+B2*r2+C2;
				ArgLnZ1 = Ar1r1_pl_B1r1_pl_C1/Ar2r2_pl_B1r2_pl_C1;
				ArgLnZ2 = Ar1r1_pl_B2r1_pl_C2/Ar2r2_pl_B2r2_pl_C2;
				T_di_Vphy = T/Vphy;
				Radical1 = sqrt(4.*A*C1-B1B1);
				Radical2 = sqrt(4.*A*C2-B2B2);

				WzL1 = -(Vphy/A)*r2mir1*z2miz1;
				WzL2 = Pi*T_di_Vphy*T_di_Vphy*Step(-T_di_Vphy-r1)*Step(r2+T_di_Vphy)*(Sign(PT-R1Vphy)-Sign(PT-R2Vphy));
				WzL3 = r1r1*(atan(TransAtans((Pr1+R2)/Vphyr1pT, -(Pr1+R1)/Vphyr1pT, PiMult1)) + Pi*PiMult1);
				WzL4 = r2r2*(atan(TransAtans((Pr2+R1)/Vphyr2pT, -(Pr2+R2)/Vphyr2pT, PiMult2)) + Pi*PiMult2);
				WzL5 = One_di_KK*((PR1_pl_VphyT*PR1_pl_VphyT - R1Vphy_mi_PT*R1Vphy_mi_PT)*((atan(TransAtans(R1Vphy_mi_PT/Kr1_pl_PR1_pl_VphyT, -R1Vphy_mi_PT/Kr2_pl_PR1_pl_VphyT, PiMult3)) + Pi*PiMult3) + Pi*Sign(R1Vphy_mi_PT)*Step(-PR1pVphyT_di_K-r1)*Step(r2+PR1pVphyT_di_K)) 
					             +(PR2_pl_VphyT*PR2_pl_VphyT - R2Vphy_mi_PT*R2Vphy_mi_PT)*((atan(TransAtans(-R2Vphy_mi_PT/Kr1_pl_PR2_pl_VphyT, R2Vphy_mi_PT/Kr2_pl_PR2_pl_VphyT, PiMult4)) + Pi*PiMult4) - Pi*Sign(R2Vphy_mi_PT)*Step(-PR2pVphyT_di_K-r1)*Step(r2+PR2pVphyT_di_K)));
				WzL6 = One_di_KK*(PR1_pl_VphyT*R1Vphy_mi_PT*log(ArgLnZ1) - PR2_pl_VphyT*R2Vphy_mi_PT*log(ArgLnZ2));
				Wz = WzL1+WzL2+WzL3+WzL4+WzL5+WzL6;
				UzL1 = (B1-B2)*r2mir1/A;
				UzL2 = One_di_AA*(B1*Radical1*(atan(TransAtans((TwoAr1+B1)/Radical1, -(TwoAr2+B1)/Radical1, PiMult1)) + Pi*PiMult1)
								 +B2*Radical2*(atan(TransAtans(-(TwoAr1+B2)/Radical2, (TwoAr2+B2)/Radical2, PiMult2)) + Pi*PiMult2));
				UzL3 = 0.5*One_di_AA*(-(2.*AC1-B1B1)*log(ArgLnZ1) + (2.*AC2-B2B2)*log(ArgLnZ2));
				UzL4 = r1r1*log(Ar1r1_pl_B2r1_pl_C2/Ar1r1_pl_B1r1_pl_C1) + r2r2*log(Ar2r2_pl_B1r2_pl_C1/Ar2r2_pl_B2r2_pl_C2);
				Uz = UzL1+UzL2+UzL3+UzL4;

				Func.x = Cosph*Uz;
				Func.y = Sinph*Uz;
				Func.z = VectV.z*Vphx*Uz - 2.*Vphy*Wz;
				
				S_forB += IntegWeight[IndForWeight] * Func;

BottomOfThisLoop:
				IndForPass++; IndForWeight++;
				ph += Step_ph;
			}
			S_forB = 0.5 * S_forB;
			S_forB.z = S_forB.z/VxVxpVyVy;

			GenS_forB += S_forB; 
			IntForB = Step_ph * GenS_forB;
			PrecParamInt = Fact * Max( Max( Abs(IntForB.x-PrIntForB.x), Abs(IntForB.y-PrIntForB.y)), Abs(IntForB.z-PrIntForB.z));

			AmOfPoi = AmOfPoi_mi_1*2+1;
			NotFirstPass = 1;
		}
	}

FinalDefinitionOfFieldIntegrals:
	TVector3d BufIb = Fact * IntForB;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIb;
	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIb;
}

//-------------------------------------------------------------------------

void radTArcCur::Dump(std::ostream& o, int ShortSign)
{
	radTg3d::Dump(o);
	DumpPureObjInfo(o, ShortSign);
	if(ShortSign==1) return;

	DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTArcCur::DumpPureObjInfo(std::ostream& o, int ShortSign)
{
	o << "Current carrying: ";
	o << "ArcCur";

	if(ShortSign==1) return;

	o << endl;
	o << "   {x,y,z}= {" << CircleCentrPoint.x << ','
						 << CircleCentrPoint.y << ','
						 << CircleCentrPoint.z << "}" << endl
	  << "   {rmin,rmax}= {" << R_min << ',' << R_max << "}" << endl
	  << "   {phimin,phimax}= {" << Phi_min << ',' << Phi_max << "}" << endl
	  << "   h= " << Height << endl
	  << "   nseg= " << NumberOfSectors << endl
	  << "   j= " << J_azim << endl
	  << "   Field computation mode: " << (BasedOnPrecLevel? "\"auto\"" : "\"man\"");
}

//-------------------------------------------------------------------------

void radTArcCur::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//Dumping objects that may be used by this object
	vector<pair<int, int> > vTrfKeys;
	DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys);

	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)Type_g3d();
	oStr << (char)0;
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTg3d
	DumpBin_g3d(oStr, vTrfKeys);

	//Members of radTArcCur
	//TVector3d CircleCentrPoint;
	oStr << CircleCentrPoint;

	//double R_min, R_max;
	oStr << R_min << R_max;

	//double Phi_min, Phi_max;
	oStr << Phi_min << Phi_max;

	//double Height;
	oStr << Height;

	//double J_azim;
	oStr << J_azim;

	//int NumberOfSectors;
	oStr << NumberOfSectors;

	//short BasedOnPrecLevel;
	oStr << BasedOnPrecLevel;

	//short InternalFacesAfterCut;
	oStr << InternalFacesAfterCut;

	//char J_IsNotZero;
	oStr << J_IsNotZero;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTArcCur::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTArcCurGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

int radTArcCur::SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	if(!SubdivideCoils) return 1;
	radTSend Send;
	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty())) 
	{
		Send.ErrorMessage("Radia::Error108"); return 0;
	}

	const double ZeroTol = 1.E-10;

	double kr = SubdivArray[0], kPhi = SubdivArray[2], kz = SubdivArray[4];
	double qr = SubdivArray[1], qPhi = SubdivArray[3], qz = SubdivArray[5];

	double DelPhi = Phi_max - Phi_min;
	double Del_r = R_max - R_min;

	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		double DelPhiL = DelPhi*0.666666666667*(R_max*R_max + R_max*R_min + R_min*R_min)/(R_max + R_min);
		kPhi = (kPhi < DelPhiL)? Round(DelPhiL/kPhi) : 1.;

		kr = (kr < Del_r)? Round(Del_r/kr) : 1.;
		kz = (kz < Height)? Round(Height/kz) : 1.;
	}

	if((fabs(kPhi-1.)<ZeroTol) && (fabs(kr-1.)<ZeroTol) && (fabs(kz-1.)<ZeroTol)) return 1;

	radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedArcCur(this);
	radThg NewHandle(GroupInPlaceOfThisPtr);

	const double AbsZeroTol = 5.E-12;
	double q0Phi = (fabs(kPhi-1.)>AbsZeroTol)? pow(qPhi, 1./(kPhi-1.)) : qPhi;
	double q0r = (fabs(kr-1.)>AbsZeroTol)? pow(qr, 1./(kr-1.)) : qr;
	double q0z = (fabs(kz-1.)>AbsZeroTol)? pow(qz, 1./(kz-1.)) : qz;
	double BufPhi = qPhi*q0Phi - 1., BufR = qr*q0r - 1., BufZ = qz*q0z - 1.;

	double a1Phi = (fabs(BufPhi) > AbsZeroTol)? DelPhi*(q0Phi - 1.)/BufPhi : DelPhi/kPhi;
	double a1r = (fabs(BufR) > AbsZeroTol)? Del_r*(q0r - 1.)/BufR : Del_r/kr;
	double a1z = (fabs(BufZ) > AbsZeroTol)? Height*(q0z - 1.)/BufZ : Height/kz;

	TVector3d InitNewDims(a1Phi, a1r, a1z);
	TVector3d NewDims = InitNewDims;

	short NewFacesState[6], ParentFacesState[6];
	ListFacesInternalAfterCut(ParentFacesState);

	int kPhiInt = int(kPhi), krInt = int(kr), kzInt = int(kz);
	int kPhi_mi_1 = kPhiInt-1, kr_mi_1 = krInt-1, kz_mi_1 = kzInt-1;

	TVector3d InitNewCircleCenPoi = TVector3d(CircleCentrPoint.x, CircleCentrPoint.y, CircleCentrPoint.z - 0.5*(Height - a1z));
	TVector3d NewCircleCenPoi = InitNewCircleCenPoi;

	double NewAngles[2], NewRadii[2], NewHeight = a1z;

	double &NewStartAngle = *NewAngles, &NewFinAngle = *(NewAngles+1), &NewStartRad = *NewRadii, &NewFinRad = *(NewRadii+1);
	NewStartAngle = Phi_min;
	NewFinAngle = Phi_min + a1Phi;
	double SmallDelPhi = a1Phi;

	NewStartRad = R_min;
	NewFinRad = R_min + a1r;
	double SmallDel_r = a1r;

	int NewStuffCounter = 0;
	for(int iPhi=0; iPhi<kPhiInt; iPhi++)
	{
		NewFacesState[0] = NewFacesState[1] = 1;
		if(iPhi==0) NewFacesState[0] = ParentFacesState[0];
		if(iPhi==kPhi_mi_1) NewFacesState[1] = ParentFacesState[1];

		float FloatNumOfSect = (float)((SmallDelPhi/DelPhi)*NumberOfSectors);
		int IntNumOfSect = int(FloatNumOfSect);
		int NewNumberOfSectors = ((FloatNumOfSect - IntNumOfSect) < 0.5)? IntNumOfSect : IntNumOfSect + 1;

		if(NewNumberOfSectors < 1) NewNumberOfSectors = 1;

		for(int ir=0; ir<krInt; ir++)
		{
			NewFacesState[2] = NewFacesState[3] = 1;
			if(ir==0) NewFacesState[2] = ParentFacesState[2];
			if(ir==kr_mi_1) NewFacesState[3] = ParentFacesState[3];

			for(int iz=0; iz<kzInt; iz++)
			{
				NewFacesState[4] = NewFacesState[5] = 1;
				if(iz==0) NewFacesState[4] = ParentFacesState[4];
				if(iz==kz_mi_1) NewFacesState[5] = ParentFacesState[5];

				radTArcCur* ArcCurPtr = new radTArcCur(NewCircleCenPoi, NewRadii, NewAngles, NewHeight, J_azim, NewNumberOfSectors, BasedOnPrecLevel);
				if(ArcCurPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

				radThg hg(ArcCurPtr);
				if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(hg), hg);
				else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, hg);

				ArcCurPtr->SetFacesInternalAfterCut(NewFacesState);

				NewCircleCenPoi.z += 0.5*NewHeight;
				NewHeight *= q0z;
				NewCircleCenPoi.z += 0.5*NewHeight;
			}
			NewStartRad = NewFinRad;
			SmallDel_r *= q0r;
			NewFinRad += SmallDel_r;

			NewCircleCenPoi.z = InitNewCircleCenPoi.z;
			NewHeight = a1z;
		}
		NewStartAngle = NewFinAngle;
		SmallDelPhi *= q0Phi;
		NewFinAngle += SmallDelPhi;

		NewStartRad = R_min;
		NewFinRad = R_min + a1r;
		SmallDel_r = a1r;
	}
	In_hg = NewHandle;
	return 1;
}

//-------------------------------------------------------------------------

void radTArcCur::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
	if(pFieldKey->M_) return;
	else
	{
		TVector3d CP = CentrPoint;
		radTrans* pTrans = (g3dListOfTransform.empty())? 0 : (radTrans*)((*(g3dListOfTransform.begin())).Handler_g.rep);

		radTrans TotTrans;
		if(pTrans != 0)
		{
			if(pBaseTrans != 0) 
			{
				TrProduct(pBaseTrans, pTrans, TotTrans);
				pTrans = &TotTrans;
			}
		}
		else
		{
			if(pBaseTrans != 0) pTrans = pBaseTrans;
		}

		if(pTrans != 0) CP = pTrans->TrPoint(CP);
		radTPairOfVect3d Pair(CP);
		if(pFieldKey->J_)
		{
			double Phic = 0.5*(Phi_min + Phi_max);
			TVector3d J(-J_azim*sin(Phic), J_azim*cos(Phic), 0.);
			if(pTrans != 0) J = pTrans->TrVectField(J);
			Pair.V2 = J;
		}
		else
		{
			radTCompCriterium CompCriterium;
			TVector3d ZeroVect(0.,0.,0.);
			radTField Field(*pFieldKey, CompCriterium, CP, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
			g3dSrcPtr->B_genComp(&Field);
			Pair.V2 = (pFieldKey->M_)? Field.M : ((pFieldKey->B_)? Field.B : ((pFieldKey->H_)? Field.H : ((pFieldKey->A_)? Field.A : ZeroVect)));
		}
		pVectPairOfVect3d->push_back(Pair);
	}
}

//-------------------------------------------------------------------------

void radTArcCur::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique) 
{//adds vertices of sectors 
	
	double HalfHeight = 0.5*Height;
	double dPhi = 0.;
	if(NumberOfSectors > 0) dPhi = (Phi_max - Phi_min)/NumberOfSectors;
	double CurPhi = Phi_min;
	
	for(int i=0; i<=NumberOfSectors; i++)
	{
		double CosPhi = cos(CurPhi), SinPhi = sin(CurPhi);
		double RelInnerX = R_min*CosPhi, RelOuterX = R_max*CosPhi;
		double RelInnerY = R_min*SinPhi, RelOuterY = R_max*SinPhi;

		double x = CircleCentrPoint.x + RelInnerX;
		double y = CircleCentrPoint.y + RelInnerY;
		double z = CircleCentrPoint.z - HalfHeight;
		TVector3d LowerInner(x, y, z);

		TVector3d UpperInner = LowerInner;
		UpperInner.z += Height;

		x = CircleCentrPoint.x + RelOuterX;
		y = CircleCentrPoint.y + RelOuterY;
		TVector3d LowerOuter(x, y, z);

		TVector3d UpperOuter = LowerOuter;
		UpperOuter.z += Height;

		OutVect.push_back(LowerInner);
		OutVect.push_back(UpperInner);
		OutVect.push_back(LowerOuter);
		OutVect.push_back(UpperOuter);

		CurPhi += dPhi;
	}
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTg3dGraphPresent* radTBackgroundFieldSource::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTBackgroundFldSrcGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------
