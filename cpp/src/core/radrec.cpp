/*-------------------------------------------------------------------------
*
* File name:      radrec.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 rectangular parallelepiped with constant magnetization 
*                 or currect density
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

//-------------------------------------------------------------------------
// Implementation of class radTRecMag - a class of objects of rectangular
// parallelipipedic shape capable to generate magnetic field.
// RecMag is derived from radTg3d.
//-------------------------------------------------------------------------

#include "radrec.h"
#include "radg3dgr.h"
#include "radgroup.h"
#include "radappl.h"
#include "radg3da1.h"

#ifndef _INC_MATH
#include <math.h>
#endif

//#ifdef __GCC__
//#include <strstream.h>
//#else
#include <sstream>
//#endif

//-------------------------------------------------------------------------

extern radTYield radYield;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRecMag::B_comp(radTField* FieldPtr)
{
	const double ZeroToler = 1.E-20;

	short J_is_Zero = 1;
	if(J.x!=0. || J.y!=0. || J.z!=0.) J_is_Zero = 0;
	TVector3d P_min_CenPo = FieldPtr->P - CentrPoint;

	double AlreadyComputedStuff[10];
	AlreadyComputedStuff[0] = P_min_CenPo.x; 
	AlreadyComputedStuff[1] = P_min_CenPo.y; 
	AlreadyComputedStuff[2] = P_min_CenPo.z;
	AlreadyComputedStuff[3] = P_min_CenPo.x*P_min_CenPo.x; 
	AlreadyComputedStuff[4] = P_min_CenPo.y*P_min_CenPo.y; 
	AlreadyComputedStuff[5] = P_min_CenPo.z*P_min_CenPo.z;
	AlreadyComputedStuff[6] = Dimensions.x*Dimensions.x; 
	AlreadyComputedStuff[7] = Dimensions.y*Dimensions.y; 
	AlreadyComputedStuff[8] = Dimensions.z*Dimensions.z;
	AlreadyComputedStuff[9] = (AlreadyComputedStuff[6] + AlreadyComputedStuff[7] + AlreadyComputedStuff[8])/
							  (AlreadyComputedStuff[3] + AlreadyComputedStuff[4] + AlreadyComputedStuff[5]);
	double SquaredMltpolCritRatio = AlreadyComputedStuff[9];
	double* MltThr = FieldPtr->CompCriterium.MltplThresh;
	if(J_is_Zero && !FieldPtr->FieldKey.A_ && !FieldPtr->FieldKey.Phi_ && 
	   (SquaredMltpolCritRatio<MltThr[0] || SquaredMltpolCritRatio<MltThr[1] || SquaredMltpolCritRatio<MltThr[2] || SquaredMltpolCritRatio<MltThr[3]))
	{
		if(SquaredMltpolCritRatio<MltThr[0]) return;
		else B_compMultipole(FieldPtr, AlreadyComputedStuff); return;
	}

	if(radYield.Check()==0) return; // To allow multitasking on Mac: consider better places for this

	TVector3d T0(0,0,0), T1(0,0,0), T(0,0,0), S(1.,1.,1.);
	const double dConst1 = 0.5;
	TVector3d HalfDim = dConst1*Dimensions;

	struct { double x[2],y[2],z[2];} BfSt;
	for(int ii=0; ii<=1; ii++)
	{
		int Eps=ii*2-1;
		BfSt.x[ii] = -P_min_CenPo.x+Eps*HalfDim.x;
		BfSt.y[ii] = -P_min_CenPo.y+Eps*HalfDim.y;
		BfSt.z[ii] = -P_min_CenPo.z+Eps*HalfDim.z;

// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(BfSt.x[ii]==0.) BfSt.x[ii] = radCR.AbsRandMagnitude(HalfDim.x);
		if(BfSt.y[ii]==0.) BfSt.y[ii] = radCR.AbsRandMagnitude(HalfDim.y);
		if(BfSt.z[ii]==0.) BfSt.z[ii] = radCR.AbsRandMagnitude(HalfDim.z);
	}

	double x0 = BfSt.x[0], x1 = BfSt.x[1];
	double y0 = BfSt.y[0], y1 = BfSt.y[1];
	double z0 = BfSt.z[0], z1 = BfSt.z[1];

	if(FieldPtr->FieldKey.M_)
	{
		if((x0*x1<0) && (y0*y1<0) && (z0*z1<0)) FieldPtr->M += Magn;
		//if(!FieldPtr->FieldKey.B_ && !FieldPtr->FieldKey.H_ && !FieldPtr->FieldKey.A_ && !FieldPtr->FieldKey.Phi_) return;
	}
	if(FieldPtr->FieldKey.J_) //OC061008
	{
		if(!J_is_Zero) //to drop this test?
		{
			if((x0*x1<0) && (y0*y1<0) && (z0*z1<0)) FieldPtr->J += J;
		}
	}
	if((!FieldPtr->FieldKey.B_) && (!FieldPtr->FieldKey.H_) && (!FieldPtr->FieldKey.A_) && (!FieldPtr->FieldKey.Phi_) && (!FieldPtr->FieldKey.PreRelax_)) return;

	double x0e2 = x0*x0, x1e2 = x1*x1;
	double y0e2 = y0*y0, y1e2 = y1*y1;
	double z0e2 = z0*z0, z1e2 = z1*z1;

	double D000 = sqrt(x0e2+y0e2+z0e2);	
	double D100 = sqrt(x1e2+y0e2+z0e2);	
	double D010 = sqrt(x0e2+y1e2+z0e2);	
	double D110 = sqrt(x1e2+y1e2+z0e2);	
	double D001 = sqrt(x0e2+y0e2+z1e2);	
	double D101 = sqrt(x1e2+y0e2+z1e2);	
	double D011 = sqrt(x0e2+y1e2+z1e2);	
	double D111 = sqrt(x1e2+y1e2+z1e2);	

	const double Pi = 3.141592653589793238;
	double PiMult1, PiMult2, PiMult3;
	PiMult1 = PiMult2 = PiMult3 = 0.;

	T0.x = atan(TransAtans(TransAtans(y0*z0/(x0*D000), -y0*z1/(x0*D001), PiMult1), 
						   TransAtans(-y1*z0/(x0*D010), y1*z1/(x0*D011), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T1.x = atan(TransAtans(TransAtans(-y0*z0/(x1*D100), y0*z1/(x1*D101), PiMult1), 
						   TransAtans(y1*z0/(x1*D110), -y1*z1/(x1*D111), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T0.y = atan(TransAtans(TransAtans(x0*z0/(y0*D000), -x0*z1/(y0*D001), PiMult1), 
						   TransAtans(-x1*z0/(y0*D100), x1*z1/(y0*D101), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T1.y = atan(TransAtans(TransAtans(-x0*z0/(y1*D010), x0*z1/(y1*D011), PiMult1), 
						   TransAtans(x1*z0/(y1*D110), -x1*z1/(y1*D111), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T0.z = atan(TransAtans(TransAtans(x0*y0/(z0*D000), -x1*y0/(z0*D100), PiMult1), 
						   TransAtans(-x0*y1/(z0*D010), x1*y1/(z0*D110), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T1.z = atan(TransAtans(TransAtans(-x0*y0/(z1*D001), x1*y0/(z1*D101), PiMult1), 
						   TransAtans(x0*y1/(z1*D011), -x1*y1/(z1*D111), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);

	double AbsRandD000 = 10.*radCR.AbsRandMagnitude(D000);
	double AbsRandD010 = 10.*radCR.AbsRandMagnitude(D010);
	double AbsRandD001 = 10.*radCR.AbsRandMagnitude(D001);
	double AbsRandD011 = 10.*radCR.AbsRandMagnitude(D011);
	double AbsRandD100 = 10.*radCR.AbsRandMagnitude(D100);
	double AbsRandD110 = 10.*radCR.AbsRandMagnitude(D110);
	double AbsRandD101 = 10.*radCR.AbsRandMagnitude(D101);
	double AbsRandD111 = 10.*radCR.AbsRandMagnitude(D111);

	double z0plD100 = z0+D100; if(z0plD100 < AbsRandD100) z0plD100 = 0.5*(x1e2 + y0e2)/Abs(z0);
	double z1plD101 = z1+D101; if(z1plD101 < AbsRandD101) z1plD101 = 0.5*(x1e2 + y0e2)/Abs(z1);
	double z1plD001 = z1+D001; if(z1plD001 < AbsRandD001) z1plD001 = 0.5*(x0e2 + y0e2)/Abs(z1);
	double z0plD000 = z0+D000; if(z0plD000 < AbsRandD000) z0plD000 = 0.5*(x0e2 + y0e2)/Abs(z0);
	double z0plD010 = z0+D010; if(z0plD010 < AbsRandD010) z0plD010 = 0.5*(x0e2 + y1e2)/Abs(z0);
	double z1plD011 = z1+D011; if(z1plD011 < AbsRandD011) z1plD011 = 0.5*(x0e2 + y1e2)/Abs(z1);
	double z1plD111 = z1+D111; if(z1plD111 < AbsRandD111) z1plD111 = 0.5*(x1e2 + y1e2)/Abs(z1);
	double z0plD110 = z0+D110; if(z0plD110 < AbsRandD110) z0plD110 = 0.5*(x1e2 + y1e2)/Abs(z0);

	double y0plD100 = y0+D100; if(y0plD100 < AbsRandD100) y0plD100 = 0.5*(x1e2 + z0e2)/Abs(y0);
	double y1plD110 = y1+D110; if(y1plD110 < AbsRandD110) y1plD110 = 0.5*(x1e2 + z0e2)/Abs(y1);
	double y1plD010 = y1+D010; if(y1plD010 < AbsRandD010) y1plD010 = 0.5*(x0e2 + z0e2)/Abs(y1);
	double y0plD000 = y0+D000; if(y0plD000 < AbsRandD000) y0plD000 = 0.5*(x0e2 + z0e2)/Abs(y0);
	double y0plD001 = y0+D001; if(y0plD001 < AbsRandD001) y0plD001 = 0.5*(x0e2 + z1e2)/Abs(y0);
	double y1plD011 = y1+D011; if(y1plD011 < AbsRandD011) y1plD011 = 0.5*(x0e2 + z1e2)/Abs(y1);
	double y1plD111 = y1+D111; if(y1plD111 < AbsRandD111) y1plD111 = 0.5*(x1e2 + z1e2)/Abs(y1);
	double y0plD101 = y0+D101; if(y0plD101 < AbsRandD101) y0plD101 = 0.5*(x1e2 + z1e2)/Abs(y0);

	double x0plD010 = x0+D010; if(x0plD010 < AbsRandD010) x0plD010 = 0.5*(y1e2 + z0e2)/Abs(x0);
	double x1plD110 = x1+D110; if(x1plD110 < AbsRandD110) x1plD110 = 0.5*(y1e2 + z0e2)/Abs(x1);
	double x1plD100 = x1+D100; if(x1plD100 < AbsRandD100) x1plD100 = 0.5*(y0e2 + z0e2)/Abs(x1);
	double x0plD000 = x0+D000; if(x0plD000 < AbsRandD000) x0plD000 = 0.5*(y0e2 + z0e2)/Abs(x0);
	double x0plD001 = x0+D001; if(x0plD001 < AbsRandD001) x0plD001 = 0.5*(y0e2 + z1e2)/Abs(x0);
	double x1plD101 = x1+D101; if(x1plD101 < AbsRandD101) x1plD101 = 0.5*(y0e2 + z1e2)/Abs(x1);
	double x1plD111 = x1+D111; if(x1plD111 < AbsRandD111) x1plD111 = 0.5*(y1e2 + z1e2)/Abs(x1);
	double x0plD011 = x0+D011; if(x0plD011 < AbsRandD011) x0plD011 = 0.5*(y1e2 + z1e2)/Abs(x0);

	double z0plD100_di_z1plD101 = z0plD100/z1plD101;
	double z1plD001_di_z0plD000 = z1plD001/z0plD000;
	double z0plD010_di_z1plD011 = z0plD010/z1plD011;
	double z1plD111_di_z0plD110 = z1plD111/z0plD110;
	double y0plD100_di_y1plD110 = y0plD100/y1plD110;
	double y1plD010_di_y0plD000 = y1plD010/y0plD000;
	double y0plD001_di_y1plD011 = y0plD001/y1plD011;
	double y1plD111_di_y0plD101 = y1plD111/y0plD101;
	double x0plD010_di_x1plD110 = x0plD010/x1plD110;
	double x1plD100_di_x0plD000 = x1plD100/x0plD000;
	double x0plD001_di_x1plD101 = x0plD001/x1plD101;
	double x1plD111_di_x0plD011 = x1plD111/x0plD011;

	const double dConst2 = 1./4./Pi;
	const double ConstForJ = 0.0001;

	double ln_z0plD100_di_z1plD101, ln_z1plD001_di_z0plD000, ln_z0plD010_di_z1plD011, ln_z1plD111_di_z0plD110,
		   ln_y0plD100_di_y1plD110, ln_y1plD010_di_y0plD000, ln_y0plD001_di_y1plD011, ln_y1plD111_di_y0plD101,
		   ln_x0plD010_di_x1plD110, ln_x1plD100_di_x0plD000, ln_x0plD001_di_x1plD101, ln_x1plD111_di_x0plD011;

	if(FieldPtr->FieldKey.A_ || FieldPtr->FieldKey.Phi_ || !J_is_Zero)
	{
		ln_z0plD100_di_z1plD101 = log(z0plD100_di_z1plD101);
		ln_z1plD001_di_z0plD000 = log(z1plD001_di_z0plD000);
		ln_z0plD010_di_z1plD011 = log(z0plD010_di_z1plD011);
		ln_z1plD111_di_z0plD110 = log(z1plD111_di_z0plD110);
		ln_y0plD100_di_y1plD110 = log(y0plD100_di_y1plD110);
		ln_y1plD010_di_y0plD000 = log(y1plD010_di_y0plD000);
		ln_y0plD001_di_y1plD011 = log(y0plD001_di_y1plD011);
		ln_y1plD111_di_y0plD101 = log(y1plD111_di_y0plD101);
		ln_x0plD010_di_x1plD110 = log(x0plD010_di_x1plD110);
		ln_x1plD100_di_x0plD000 = log(x1plD100_di_x0plD000);
		ln_x0plD001_di_x1plD101 = log(x0plD001_di_x1plD101);
		ln_x1plD111_di_x0plD011 = log(x1plD111_di_x0plD011);

		TVector3d BufVect(x0*T0.x + x1*T1.x
						 +y0*(ln_z0plD100_di_z1plD101+ln_z1plD001_di_z0plD000)
						 +y1*(ln_z0plD010_di_z1plD011+ln_z1plD111_di_z0plD110)
						 +z0*(ln_y0plD100_di_y1plD110+ln_y1plD010_di_y0plD000)
						 +z1*(ln_y0plD001_di_y1plD011+ln_y1plD111_di_y0plD101),
						  y0*T0.y + y1*T1.y
						 +x0*(ln_z0plD010_di_z1plD011+ln_z1plD001_di_z0plD000)
						 +x1*(ln_z0plD100_di_z1plD101+ln_z1plD111_di_z0plD110)
						 +z0*(ln_x0plD010_di_x1plD110+ln_x1plD100_di_x0plD000)
						 +z1*(ln_x0plD001_di_x1plD101+ln_x1plD111_di_x0plD011),
						  z0*T0.z + z1*T1.z
						 +y0*(ln_x0plD001_di_x1plD101+ln_x1plD100_di_x0plD000)
						 +y1*(ln_x0plD010_di_x1plD110+ln_x1plD111_di_x0plD011)
						 +x0*(ln_y0plD001_di_y1plD011+ln_y1plD010_di_y0plD000)
						 +x1*(ln_y0plD100_di_y1plD110+ln_y1plD111_di_y0plD101));
		if(J_is_Zero)
		{
			if(FieldPtr->FieldKey.A_)
			{
				TVector3d BufForA(Magn.y*BufVect.z-Magn.z*BufVect.y, 
								  Magn.z*BufVect.x-Magn.x*BufVect.z, 
								  Magn.x*BufVect.y-Magn.y*BufVect.x);
				FieldPtr->A += dConst2*BufForA;
			}
			if(FieldPtr->FieldKey.Phi_)	FieldPtr->Phi += dConst2*(Magn*BufVect);
		}
		else
		{
			if(FieldPtr->FieldKey.A_)
			{
				FieldPtr->A +=(ConstForJ
							  *((x0*x0*T0.x+x1*x1*T1.x
							   +y0*y0*T0.y+y1*y1*T1.y
							   +z0*z0*T0.z+z1*z1*T1.z)/2.
							   +x0*y0*ln_z1plD001_di_z0plD000
							   +x1*y0*ln_z0plD100_di_z1plD101
							   +x0*y1*ln_z0plD010_di_z1plD011
							   +x1*y1*ln_z1plD111_di_z0plD110
							   +x0*z0*ln_y1plD010_di_y0plD000
							   +x1*z0*ln_y0plD100_di_y1plD110
							   +x0*z1*ln_y0plD001_di_y1plD011
							   +x1*z1*ln_y1plD111_di_y0plD101
							   +y0*z0*ln_x1plD100_di_x0plD000
							   +y1*z0*ln_x0plD010_di_x1plD110
							   +y0*z1*ln_x0plD001_di_x1plD101
							   +y1*z1*ln_x1plD111_di_x0plD011))*J;
			}
			if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
			{
				TVector3d BufForB(J.y*BufVect.z-J.z*BufVect.y, 
								  J.z*BufVect.x-J.x*BufVect.z, 
								  J.x*BufVect.y-J.y*BufVect.x);
				FieldPtr->B += ConstForJ*BufForB;
				FieldPtr->H += ConstForJ*BufForB;
			}
		}
	}
	
	if(((FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_) && J_is_Zero) || FieldPtr->FieldKey.PreRelax_)
	{
		T = T0 + T1;
		if(!(FieldPtr->FieldKey.A_ || FieldPtr->FieldKey.Phi_))
		{
			S.x = -log(x0plD010_di_x1plD110*x1plD100_di_x0plD000*x0plD001_di_x1plD101*x1plD111_di_x0plD011);
			S.y = -log(y0plD100_di_y1plD110*y1plD010_di_y0plD000*y0plD001_di_y1plD011*y1plD111_di_y0plD101);
			S.z = -log(z0plD100_di_z1plD101*z1plD001_di_z0plD000*z0plD010_di_z1plD011*z1plD111_di_z0plD110);
		}
		else
		{
			S.x = -(ln_x0plD010_di_x1plD110+ln_x1plD100_di_x0plD000+ln_x0plD001_di_x1plD101+ln_x1plD111_di_x0plD011);
			S.y = -(ln_y0plD100_di_y1plD110+ln_y1plD010_di_y0plD000+ln_y0plD001_di_y1plD011+ln_y1plD111_di_y0plD101);
			S.z = -(ln_z0plD100_di_z1plD101+ln_z1plD001_di_z0plD000+ln_z0plD010_di_z1plD011+ln_z1plD111_di_z0plD110);
		}
		T=dConst2*T; S=dConst2*S;
		
		if(FieldPtr->FieldKey.PreRelax_)
		{
			TVector3d& RefB = FieldPtr->B;
			TVector3d& RefH = FieldPtr->H;
			TVector3d& RefA = FieldPtr->A;
			RefB.x = T.x; RefB.y = -S.z; RefB.z = -S.y;
			RefH.x = -S.z; RefH.y = T.y; RefH.z = -S.x;
			RefA.x = -S.y; RefA.y = -S.x; RefA.z = T.z;
			return;
		}

		TVector3d Str0(T.x,-S.z,-S.y), Str1(-S.z,T.y,-S.x), Str2(-S.y,-S.x,T.z);
		TMatrix3d MatrF(Str0, Str1, Str2);
		TVector3d BufH=MatrF*Magn;
		if(FieldPtr->FieldKey.B_)
		{
			TVector3d BufB=BufH;
			if((x0*x1<0) && (y0*y1<0) && (z0*z1<0)) BufB+=Magn;
			FieldPtr->B += BufB;
		}
		FieldPtr->H += BufH;
	}
}

//-------------------------------------------------------------------------

void radTRecMag::B_compMultipole(radTField* FieldPtr, double* AlreadyComputedStuff)
{
	double* MltThr = FieldPtr->CompCriterium.MltplThresh;
	double SquaredMltpolCritRatio = AlreadyComputedStuff[9];
	TVector3d T(0.,0.,0.), S(0.,0.,0.);
	double t1, t2, t3, t4, t5, t6, t7, t9, t10, t12, 
		   t13, t29, t32, t35, t36, t37, t39, t40, t42, t43, 
		   t44, t45, t46, t47, t50, t51, t53, t54, t58, t60, 
		   t64, t66, t68, t69, t70, t71, t72, t74, t75, t86, 
		   t88, t90, t91, t92, t93, t95, t110, t111, t112, t114, 
		   t116, t117, t118, t119, t120, t121, t122, t123, t124, t125, 
		   t128, t129, t131, t142, t144, t145, t146, t147, t148, t152, 
		   t153, t166, t168, t169, t170, t171, t175, t190, t191, t195, 
		   t196, t208, t212, t216, t230, t231, t233,
		   t9t13, t10t13, t10t13mumi3, t54t40, t66t72, t88t93, t120mu180, t120mu90, t122mu6, t119mu101, 
		   t121mu101, t117mu101, t124mu101, t118mu101, t123mu101, t117mu116, t118mu116, t119mu116, 
		   t123mu116, t121mu116, t124mu116,
		   t116mu8, t122mu8, t125mu8, t42mu8, t46mu8, t47mu8, t125mu6, t116mu6, t43mu23, t44mu23, t45mu23, 
		   t44mu7, t45mu7, t43mu7, t46mu2, t42mu2, t47mu2,
		   t1pt2, t1pt3, t2pt3, t42pt46, t44pt45, t42pt47, t43pt45, t46pt47, t43pt44,
		   t121pt123, t118pt124, t119pt124, t117pt123, t117pt121, t118pt119,
		   t121pt123mu11, t121pt123mu4, t118pt124mu11, t118pt124mu4, t119pt124mu11, t119pt124mu4, t117pt123mu11, t117pt123mu4,
		   t117pt121mu11, t117pt121mu4, t118pt119mu11, t118pt119mu4;
	const double c21d128 = 21.0/128.0;
	const double c3d128 = 3.0/128.0;
	const double c35d64 = 35.0/64.0;
	const double c5d64 = 5.0/64.0;
	const double c105d64 = 105.0/64.0;
	double c21d128t131t114, c21d128t142t148, c21d128t166t171, c35d64t191, c5d64t208, 
		   c35d64t212, c35d64t231, c3d128t110, c5d64t110;
	double x, y, z, Lx, Ly, Lz;

	x = AlreadyComputedStuff[0]; y = AlreadyComputedStuff[1]; z = AlreadyComputedStuff[2];
	Lx = Dimensions.x; Ly = Dimensions.y; Lz = Dimensions.z;

	t1 = AlreadyComputedStuff[3]; t2 = AlreadyComputedStuff[4]; t3 = AlreadyComputedStuff[5];
	t1pt2 = t1+t2; t1pt3 = t1+t3; t2pt3 = t2+t3;
	t4 = t1pt2+t3; t5 = sqrt(t4); t6 = t4*t4; t7 = t6*t4; t9 = t5/t7; t10 = t9*y; t12 = Lx*Ly; t13 = t12*Lz;
	t9t13 = t9*t13; t10t13 = t10*t13; t10t13mumi3 = -3.*t10t13;
	t53 = x*z;
	S.x = t10t13mumi3*z;
	T.x = (2.0*t1-t2pt3)*t9t13;
	S.y = -3.0*t9t13*t53;
	T.y = (2.0*t2-t1pt3)*t9t13;
	S.z = t10t13mumi3*x;
	T.z = (2.0*t3-t1pt2)*t9t13;
	if(SquaredMltpolCritRatio<MltThr[1]) goto FinalFieldDefinition;

	t29 = y*z; t32 = t6*t6;
	t35 = t5/t32/t4;
	t36 = AlreadyComputedStuff[6];
	t37 = t36*Lx; t39 = Ly*Lz;
	t40 = t35*t37*t39; t42 = t1*t1; t43 = t1*t2; t44 = t1*t3; t45 = t2*t3; t46 = t2*t2; t47 = t3*t3;
	t50 = t37*Ly; t51 = t50*Lz; 
	t54 = 3.0*t2pt3-4.0*t1;
	t42pt46 = t42+t46; t44pt45 = t44+t45; t42pt47 = t42+t47; t43pt45 = t43+t45; 
	t46pt47 = t46+t47; t43pt44 = t43+t44;
	t42mu8 = 8.0*t42; t46mu8 = 8.0*t46; t47mu8 = 8.0*t47; 
	t58 = (4.0*t42pt46-27.0*t43+3.0*t44pt45-t47)*t35; t60 = y*x;
	t64 = (4.0*t42pt47-27.0*t44+3.0*t43pt45-t46)*t35;
	t66 = 3.0*t1-4.0*t2+3.0*t3; t68 = t35*Lx; t69 = AlreadyComputedStuff[7];
	t70 = t69*Ly; t71 = t70*Lz; t72 = t68*t71; t74 = Lx*t70; t75 = t74*Lz;
	t86 = (4.0*t46pt47-27.0*t45+3.0*t43pt44-t42)*t35;
	t88 = -4.0*t3+3.0*t1pt2;
	t90 = AlreadyComputedStuff[8]; t91 = t90*Lz; t92 = Ly*t91; t93 = t68*t92; t95 = t12*t91;
	t54t40 = t54*t40; t66t72 = t66*t72; t88t93 = t88*t93;
	S.x += (t66t72+t88t93-(6.0*t1-t2pt3)*t40)*0.625*t29;
	S.y += (t54t40+t88t93-(6.0*t2-t1pt3)*t72)*0.625*t53;
	S.z += (t54t40+t66t72-(6.0*t3-t1pt2)*t93)*0.625*t60;
	T.x += ((t42mu8-24.0*t43pt44+6.0*t45+3.0*t46pt47)*t35*t51-t58*t75-t64*t95)*0.125;
	T.y += ((t46mu8-24.0*t43pt45+6.0*t44+3.0*t42pt47)*t35*t75-t58*t51-t86*t95)*0.125;
	T.z += ((t47mu8-24.0*t44pt45+6.0*t43+3.0*t42pt46)*t35*t95-t64*t51-t86*t75)*0.125;
	if(SquaredMltpolCritRatio<MltThr[2]) goto FinalFieldDefinition;

	t110 = t5/t32/t7; t111 = t36*t36; t112 = t111*Lx; t114 = t110*t112*t39; t116 = t42*t1; t117 = t42*t2;
	t118 = t42*t3; t119 = t1*t46; t120 = t43*t3; t121 = t1*t47; t122 = t46*t2; t123 = t46*t3; t124 = t2*t47;
	t125 = t47*t3; t128 = t112*Ly; t129 = t128*Lz;
	t131 = 5.0*t46pt47-20.0*t43pt44+10.0*t45+t42mu8;
	t142 = 5.0*t42pt47-20.0*t43pt45+10.0*t44+t46mu8;
	t144 = t110*Lx; t145 = t69*t69; t146 = t145*Ly; t147 = t146*Lz; t148 = t144*t147; t152 = Lx*t146; t153 = t152*Lz;
	t166 = t47mu8-20.0*t44pt45+5.0*t42pt46+10.0*t43;
	t168 = t90*t90; t169 = t168*Lz; t170 = Ly*t169; t171 = t144*t170; t175 = t12*t169; t190 = t110*t37;
	t191 = t190*t71; t195 = t37*t70; t196 = t195*Lz;
	t121pt123 = t121+t123; t118pt124 = t118+t124; t119pt124 = t124+t119; t117pt123 = t117+t123; t117pt121 = t117+t121; t118pt119 = t118+t119;
	t120mu180 = 180.0*t120;
	t208 = (2.0*(t116+t125+t122)-15.0*(t118pt124+t117pt123+t119+t121)+t120mu180)*t110;
	t212 = t190*t92; t216 = t50*t91; t230 = t70*t91; t231 = t144*t230; t233 = t74*t91;
	c3d128t110 = c3d128*t110; c5d64t110 = c5d64*t110;
	c21d128t131t114 = c21d128*t131*t114; c21d128t142t148 = c21d128*t142*t148; c21d128t166t171 = c21d128*t166*t171;
	c35d64t191 = c35d64*t191; c5d64t208 = c5d64*t208; c35d64t212 = c35d64*t212; c35d64t231 = c35d64*t231;
	t120mu90 = 90.0*t120; t122mu6 = 6.0*t122; t119mu101 = 101.0*t119; t121mu101 = 101.0*t121;
	t117mu101 = 101.0*t117; t124mu101 = 101.0*t124; t123mu101 = 101.0*t123; t118mu101 = 101.0*t118; t117mu116 = 116.0*t117; 
	t118mu116 = 116.0*t118; t119mu116 = 116.0*t119; t123mu116 = 116.0*t123; t121mu116 = 116.0*t121; t124mu116 = 116.0*t124; 
	t116mu8 = 8.0*t116; t122mu8 = 8.0*t122; t125mu8 = 8.0*t125; 
	t125mu6 = 6.0*t125; t116mu6 = 6.0*t116; t43mu23 = 23.0*t43; t44mu23 = 23.0*t44; t45mu23 = 23.0*t45;
	t44mu7 = 7.0*t44; t45mu7 = 7.0*t45; t43mu7 = 7.0*t43; t46mu2 = 2.0*t46; t42mu2 = 2.0*t42; t47mu2 = 2.0*t47;
	t121pt123mu11 = 11.0*t121pt123; t121pt123mu4 = 4.0*t121pt123; t118pt124mu11 = 11.0*t118pt124; t118pt124mu4 = 4.0*t118pt124;
	t119pt124mu11 = 11.0*t119pt124; t119pt124mu4 = 4.0*t119pt124; t117pt123mu11 = 11.0*t117pt123; t117pt123mu4 = 4.0*t117pt123;
	t117pt121mu11 = 11.0*t117pt121; t117pt121mu4 = 4.0*t117pt121; t118pt119mu11 = 11.0*t118pt119; t118pt119mu4 = 4.0*t118pt119;
	S.x += t29*(c35d64t191*(-t43mu23+t46mu2+t45+t42mu8+t44mu7-t47)
				+c35d64t212*(t42mu8+t43mu7-t44mu23-t46+t45+t47mu2)
				-c21d128*(-16.0*(t43pt44-t42)+t46+2.0*t45+t47)*t114
				-c105d64*(t42-t43pt44-2.0*t46pt47+t45mu7)*t231-c21d128t142t148-c21d128t166t171);
	S.y += t53*(c35d64t191*(-t43mu23+t42mu2+t44+t46mu8+t45mu7-t47)
				+c35d64t231*(t46mu8+t43mu7-t45mu23-t42+t44+t47mu2)
				-c21d128*(-16.0*(t43pt45-t46)+t42+2.0*t44+t47)*t148
				-c105d64*(t46-t43pt45-2.0*t42pt47+t44mu7)*t212-c21d128t131t114-c21d128t166t171);
	S.z += t60*(c35d64t212*(t47mu8+t45mu7-t44mu23-t46+t43+t42mu2)
				+c35d64t231*(-t45mu23+t46mu2+t43+t47mu8+t44mu7-t42)
				-c21d128*(-16.0*(t44pt45-t47)+t46+2.0*t43+t42)*t171
				-c105d64*(t47-t44pt45-2.0*t42pt46+t43mu7)*t191-c21d128t131t114-c21d128t142t148);
	T.x += c3d128t110*(t129*(16.0*t116-120.0*(t117+t118)+90.0*(t119+t121)+t120mu180-5.0*(t122+t125)-15.0*(t123+t124))
					   +t153*(t116mu6+t118pt124mu11-t117mu101+t121pt123mu4+t119mu116-t120mu90-t122mu8-t125)
					   +t175*(t116mu6+t117pt123mu11-t118mu101+t119pt124mu4-t120mu90+t121mu116-t122-t125mu8))
		   -c5d64t110*(t196*(t116mu8-t117mu116-t118pt124mu4+t120mu90-t121pt123mu11+t119mu101+t125-t122mu6)
					   +t216*(t116mu8-t118mu116-t117pt123mu4-t119pt124mu11+t121mu101+t120mu90-t125mu6+t122))+c5d64t208*t233;
	T.y += c3d128t110*(t129*(t122mu6+t121pt123mu11-t119mu101+t118pt124mu4+t117mu116-t120mu90-t116mu8-t125)
					   +t153*(16.0*t122-120.0*(t119+t123)+90.0*(t117+t124)+t120mu180-5.0*(t116+t125)-15.0*(t118+t121))
					   +t175*(t122mu6+t118pt119mu11-t123mu101+t117pt121mu4-t120mu90+t124mu116-t116-t125mu8))
		   -c5d64t110*(t196*(t122mu8-t119mu116-t121pt123mu4+t120mu90-t118pt124mu11+t117mu101+t125-t116mu6)
					   +t233*(t122mu8-t123mu116-t118pt119mu4-t117pt121mu11+t124mu101+t120mu90-t125mu6+t116))+c5d64t208*t216;
	T.z += c3d128t110*(t129*(t125mu6+t119pt124mu11-t121mu101+t117pt123mu4-t120mu90+t118mu116-t122-t116mu8)
					   +t153*(t125mu6+t117pt121mu11-t124mu101+t118pt119mu4+t123mu116-t120mu90-t122mu8-t116)
					   +t175*(16.0*t125-120.0*(t124+t121)+90.0*(t123+t118)+t120mu180-5.0*(t122+t116)-15.0*(t119+t117)))
		   -c5d64t110*(t216*(t125mu8-t121mu116-t119pt124mu4-t117pt123mu11+t118mu101+t120mu90-t116mu6+t122)
					   +t233*(t125mu8-t124mu116-t117pt121mu4+t120mu90-t118pt119mu11+t123mu101+t116-t122mu6))+c5d64t208*t196;

FinalFieldDefinition:
	const double Pi = 3.141592653589793238;
	const double dConst2 = 1./4./Pi;
	T=dConst2*T; S=dConst2*S;
	if(FieldPtr->FieldKey.PreRelax_)
	{
		FieldPtr->B = T; FieldPtr->H = S; return;
	}
	TVector3d Str0(T.x,-S.z,-S.y), Str1(-S.z,T.y,-S.x), Str2(-S.y,-S.x,T.z);
	TMatrix3d MatrF(Str0, Str1, Str2);
	TVector3d BufH=MatrF*Magn;
	if(FieldPtr->FieldKey.B_)
	{
		TVector3d BufB=BufH;
		if((Abs(x)<0.5*Lx) && (Abs(y)<0.5*Ly) && (Abs(z)<0.5*Lz)) BufB+=Magn;
		FieldPtr->B += BufB;
	}
	FieldPtr->H += BufH;
}

//-------------------------------------------------------------------------

void radTRecMag::B_intComp(radTField* FieldPtr)
{
	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

// An analytical algorithm for infinite Field Integrals:
	TVector3d CenPo_mi_StPo = CentrPoint - FieldPtr->P;
	TVector3d HalfDim = 0.5 * Dimensions;
	TVector3d P1 = CenPo_mi_StPo - HalfDim;
	TVector3d P2 = CenPo_mi_StPo + HalfDim;
	TVector3d V = FieldPtr->NextP - FieldPtr->P;

	double ModV = sqrt(V.x*V.x + V.y*V.y + V.z*V.z);
	V.x /= ModV;  V.y /= ModV;  V.z /= ModV;

	const double Pi = 3.141592653589793238;
	const double ZeroToler = 1.E-06; // Relative tolerance to switch to special cases
	const double SmallestRelTolerV = 1.E-12; // Relative tolerance to repair trapping V.i to zero at general case

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(CentrPoint.z);

	TMatrix3d F;
	TVector3d G;

	short J_is_Zero = 1;
	if(J.x!=0. || J.y!=0. || J.z!=0.) J_is_Zero = 0;

// Tests for special cases:
	double AbsVx = Abs(V.x);
	double AbsVy = Abs(V.y);
	double AbsVz = Abs(V.z);
	if(AbsVx<ZeroToler && AbsVy<ZeroToler)
	{
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(P1.x==0.) P1.x = AbsRandX;
		if(P1.y==0.) P1.y = AbsRandY;
		if(P2.x==0.) P2.x = AbsRandX;
		if(P2.y==0.) P2.y = AbsRandY;

		B_intUtilSpecCaseZeroVxVy(P1, P2, J_is_Zero, F, G);
		goto FinalDefinitionOfFieldIntegrals;
	}
	if(AbsVx<ZeroToler && AbsVz<ZeroToler) 
	{
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(P1.x==0.) P1.x = AbsRandX;
		if(P1.z==0.) P1.z = AbsRandZ;
		if(P2.x==0.) P2.x = AbsRandX;
		if(P2.z==0.) P2.z = AbsRandZ;

		TVector3d LocP1(P1.x, P1.z, P1.y), LocP2(P2.x, P2.z, P2.y), LocG;
		TMatrix3d LocF;
		B_intUtilSpecCaseZeroVxVy(LocP1, LocP2, J_is_Zero, LocF, LocG);
		TVector3d& F_str0 = F.Str0;
		TVector3d& F_str1 = F.Str1;
		TVector3d& F_str2 = F.Str2;
		TVector3d& LocF_str0 = LocF.Str0;
		TVector3d& LocF_str1 = LocF.Str1;
		TVector3d& LocF_str2 = LocF.Str2;
		F_str0.x = LocF_str0.y; F_str0.y = LocF_str0.x; F_str0.z = LocF_str0.z;
		F_str1.x = LocF_str2.y; F_str1.y = LocF_str2.x; F_str1.z = LocF_str2.z;
		F_str2.x = LocF_str1.y; F_str2.y = LocF_str1.x; F_str2.z = LocF_str1.z;
		G.x = LocG.x; G.y = LocG.z; G.z = LocG.y;
		goto FinalDefinitionOfFieldIntegrals;
	}
	if(AbsVy<ZeroToler && AbsVz<ZeroToler) 
	{
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(P1.y==0.) P1.y = AbsRandY;
		if(P1.z==0.) P1.z = AbsRandZ;
		if(P2.y==0.) P2.y = AbsRandY;
		if(P2.z==0.) P2.z = AbsRandZ;

		TVector3d LocP1(P1.z, P1.y, P1.x), LocP2(P2.z, P2.y, P2.x), LocG;
		TMatrix3d LocF;
		B_intUtilSpecCaseZeroVxVy(LocP1, LocP2, J_is_Zero, LocF, LocG);
		TVector3d& F_str0 = F.Str0;
		TVector3d& F_str1 = F.Str1;
		TVector3d& F_str2 = F.Str2;
		TVector3d& LocF_str0 = LocF.Str0;
		TVector3d& LocF_str1 = LocF.Str1;
		TVector3d& LocF_str2 = LocF.Str2;
		F_str0.x = LocF_str2.z; F_str0.y = LocF_str2.y; F_str0.z = LocF_str2.x;
		F_str1.x = LocF_str1.z; F_str1.y = LocF_str1.y; F_str1.z = LocF_str1.x;
		F_str2.x = LocF_str0.z; F_str2.y = LocF_str0.y; F_str2.z = LocF_str0.x;
		G.x = LocG.z; G.y = LocG.y; G.z = LocG.x;
		goto FinalDefinitionOfFieldIntegrals;
	}

	{
// Prevent trapping each of V.x,V.y,V.z to Zero separately???
		if(AbsVx<SmallestRelTolerV) V.x = SmallestRelTolerV;
		if(AbsVy<SmallestRelTolerV) V.y = SmallestRelTolerV;
		if(AbsVz<SmallestRelTolerV) V.z = SmallestRelTolerV;
//GeneralCaseStart:
		double vxvx = V.x*V.x;
		double vyvy = V.y*V.y;
		double vzvz = V.z*V.z;
		double vxvx_p_vyvy = vxvx + vyvy;
		double vxvx_p_vzvz = vxvx + vzvz;
		double vyvy_p_vzvz = vyvy + vzvz;
		double vzvz_m_vyvy = vzvz - vyvy;
		double vzvz_m_vxvx = vzvz - vxvx;
		double vxvx_m_vyvy = vxvx - vyvy;
		double vyvy_p_vzvz_p_2vxvx = vyvy_p_vzvz + 2.*vxvx;
		double vxvx_p_vzvz_p_2vyvy = vxvx_p_vzvz + 2.*vyvy;
		double vxvx_p_vyvy_p_2vzvz = vxvx_p_vyvy + 2.*vzvz;
		double One_d_vxvxpvzvz = 1./vxvx_p_vzvz;
		double One_d_vxvxpvyvy = 1./vxvx_p_vyvy;
		double One_d_vyvypvzvz = 1./vyvy_p_vzvz;
		double One_d_vxvxpvzvz_d_vyvypvzvz = One_d_vxvxpvzvz*One_d_vyvypvzvz;
		double One_d_vxvxpvyvy_d_vxvxpvzvz = One_d_vxvxpvzvz*One_d_vxvxpvyvy;
		double One_d_vxvxpvyvy_d_vyvypvzvz = One_d_vxvxpvyvy*One_d_vyvypvzvz;
		double vxx1 = V.x*P1.x;
		double vxy1 = V.x*P1.y;
		double vxz1 = V.x*P1.z;
		double vyx1 = V.y*P1.x;
		double vyy1 = V.y*P1.y;
		double vyz1 = V.y*P1.z;
		double vzx1 = V.z*P1.x;
		double vzy1 = V.z*P1.y;
		double vzz1 = V.z*P1.z;
		double vxx2 = V.x*P2.x;
		double vxy2 = V.x*P2.y;
		double vxz2 = V.x*P2.z;
		double vyx2 = V.y*P2.x;
		double vyy2 = V.y*P2.y;
		double vyz2 = V.y*P2.z;
		double vzx2 = V.z*P2.x;
		double vzy2 = V.z*P2.y;
		double vzz2 = V.z*P2.z;   		  // Remove this ???
		double vzy1_m_vyz1 = vzy1 - vyz1; if(vzy1_m_vyz1==0.) vzy1_m_vyz1 = SmallestRelTolerV*AbsRandY;
		double vzx1_m_vxz1 = vzx1 - vxz1; if(vzx1_m_vxz1==0.) vzx1_m_vxz1 = SmallestRelTolerV*AbsRandX;
		double vxy1_m_vyx1 = vxy1 - vyx1; if(vxy1_m_vyx1==0.) vxy1_m_vyx1 = SmallestRelTolerV*AbsRandY;
		double vzy2_m_vyz1 = vzy2 - vyz1; if(vzy2_m_vyz1==0.) vzy2_m_vyz1 = SmallestRelTolerV*AbsRandY;
		double vzx1_m_vxz2 = vzx1 - vxz2; if(vzx1_m_vxz2==0.) vzx1_m_vxz2 = SmallestRelTolerV*AbsRandX;
		double vxy1_m_vyx2 = vxy1 - vyx2; if(vxy1_m_vyx2==0.) vxy1_m_vyx2 = SmallestRelTolerV*AbsRandY;
		double vzy1_m_vyz2 = vzy1 - vyz2; if(vzy1_m_vyz2==0.) vzy1_m_vyz2 = SmallestRelTolerV*AbsRandY;
		double vzx2_m_vxz1 = vzx2 - vxz1; if(vzx2_m_vxz1==0.) vzx2_m_vxz1 = SmallestRelTolerV*AbsRandX;
		double vxy2_m_vyx1 = vxy2 - vyx1; if(vxy2_m_vyx1==0.) vxy2_m_vyx1 = SmallestRelTolerV*AbsRandY;
		double vzy2_m_vyz2 = vzy2 - vyz2; if(vzy2_m_vyz2==0.) vzy2_m_vyz2 = SmallestRelTolerV*AbsRandY;
		double vzx2_m_vxz2 = vzx2 - vxz2; if(vzx2_m_vxz2==0.) vzx2_m_vxz2 = SmallestRelTolerV*AbsRandX;
		double vxy2_m_vyx2 = vxy2 - vyx2; if(vxy2_m_vyx2==0.) vxy2_m_vyx2 = SmallestRelTolerV*AbsRandY;
		double ArgAtanX111 = (V.z*vzx1_m_vxz1 - V.y*vxy1_m_vyx1)/vzy1_m_vyz1;
		double ArgAtanX211 = (V.z*vzx2_m_vxz1 - V.y*vxy1_m_vyx2)/vzy1_m_vyz1;
		double ArgAtanX121 = (V.z*vzx1_m_vxz1 - V.y*vxy2_m_vyx1)/vzy2_m_vyz1;
		double ArgAtanX221 = (V.z*vzx2_m_vxz1 - V.y*vxy2_m_vyx2)/vzy2_m_vyz1;
		double ArgAtanX112 = (V.z*vzx1_m_vxz2 - V.y*vxy1_m_vyx1)/vzy1_m_vyz2;
		double ArgAtanX212 = (V.z*vzx2_m_vxz2 - V.y*vxy1_m_vyx2)/vzy1_m_vyz2;
		double ArgAtanX122 = (V.z*vzx1_m_vxz2 - V.y*vxy2_m_vyx1)/vzy2_m_vyz2;
		double ArgAtanX222 = (V.z*vzx2_m_vxz2 - V.y*vxy2_m_vyx2)/vzy2_m_vyz2;
		double ArgAtanY111 = (V.z*vzy1_m_vyz1 + V.x*vxy1_m_vyx1)/vzx1_m_vxz1;
		double ArgAtanY211 = (V.z*vzy1_m_vyz1 + V.x*vxy1_m_vyx2)/vzx2_m_vxz1;
		double ArgAtanY121 = (V.z*vzy2_m_vyz1 + V.x*vxy2_m_vyx1)/vzx1_m_vxz1;
		double ArgAtanY221 = (V.z*vzy2_m_vyz1 + V.x*vxy2_m_vyx2)/vzx2_m_vxz1;
		double ArgAtanY112 = (V.z*vzy1_m_vyz2 + V.x*vxy1_m_vyx1)/vzx1_m_vxz2;
		double ArgAtanY212 = (V.z*vzy1_m_vyz2 + V.x*vxy1_m_vyx2)/vzx2_m_vxz2;
		double ArgAtanY122 = (V.z*vzy2_m_vyz2 + V.x*vxy2_m_vyx1)/vzx1_m_vxz2;
		double ArgAtanY222 = (V.z*vzy2_m_vyz2 + V.x*vxy2_m_vyx2)/vzx2_m_vxz2;
		double ArgAtanZ111 = (V.x*vzx1_m_vxz1 + V.y*vzy1_m_vyz1)/vxy1_m_vyx1;
		double ArgAtanZ211 = (V.x*vzx2_m_vxz1 + V.y*vzy1_m_vyz1)/vxy1_m_vyx2;
		double ArgAtanZ121 = (V.x*vzx1_m_vxz1 + V.y*vzy2_m_vyz1)/vxy2_m_vyx1;
		double ArgAtanZ221 = (V.x*vzx2_m_vxz1 + V.y*vzy2_m_vyz1)/vxy2_m_vyx2;
		double ArgAtanZ112 = (V.x*vzx1_m_vxz2 + V.y*vzy1_m_vyz2)/vxy1_m_vyx1;
		double ArgAtanZ212 = (V.x*vzx2_m_vxz2 + V.y*vzy1_m_vyz2)/vxy1_m_vyx2;
		double ArgAtanZ122 = (V.x*vzx1_m_vxz2 + V.y*vzy2_m_vyz2)/vxy2_m_vyx1;
		double ArgAtanZ222 = (V.x*vzx2_m_vxz2 + V.y*vzy2_m_vyz2)/vxy2_m_vyx2;
		double PiMult = 0.;
		double SumAtanXy1z1 = atan(TransAtans(ArgAtanX111, -ArgAtanX211, PiMult)) + Pi*PiMult;
		double SumAtanXy2z1 = atan(TransAtans(ArgAtanX121, -ArgAtanX221, PiMult)) + Pi*PiMult;
		double SumAtanXy1z2 = atan(TransAtans(ArgAtanX112, -ArgAtanX212, PiMult)) + Pi*PiMult;
		double SumAtanXy2z2 = atan(TransAtans(ArgAtanX122, -ArgAtanX222, PiMult)) + Pi*PiMult;
		double SumAtanYx1z1 = atan(TransAtans(ArgAtanY111, -ArgAtanY121, PiMult)) + Pi*PiMult;
		double SumAtanYx2z1 = atan(TransAtans(ArgAtanY211, -ArgAtanY221, PiMult)) + Pi*PiMult;
		double SumAtanYx1z2 = atan(TransAtans(ArgAtanY112, -ArgAtanY122, PiMult)) + Pi*PiMult;
		double SumAtanYx2z2 = atan(TransAtans(ArgAtanY212, -ArgAtanY222, PiMult)) + Pi*PiMult;
		double SumAtanZx1y1 = atan(TransAtans(ArgAtanZ111, -ArgAtanZ112, PiMult)) + Pi*PiMult;
		double SumAtanZx2y1 = atan(TransAtans(ArgAtanZ211, -ArgAtanZ212, PiMult)) + Pi*PiMult;
		double SumAtanZx1y2 = atan(TransAtans(ArgAtanZ121, -ArgAtanZ122, PiMult)) + Pi*PiMult;
		double SumAtanZx2y2 = atan(TransAtans(ArgAtanZ221, -ArgAtanZ222, PiMult)) + Pi*PiMult;
		double vxy1mvyx1_mu_vxy1mvyx1 = vxy1_m_vyx1*vxy1_m_vyx1;
		double vxy1mvyx2_mu_vxy1mvyx2 = vxy1_m_vyx2*vxy1_m_vyx2;
		double vxy2mvyx1_mu_vxy2mvyx1 = vxy2_m_vyx1*vxy2_m_vyx1;
		double vxy2mvyx2_mu_vxy2mvyx2 = vxy2_m_vyx2*vxy2_m_vyx2;
		double vzy1mvyz1_mu_vzy1mvyz1 = vzy1_m_vyz1*vzy1_m_vyz1;
		double vzy2mvyz1_mu_vzy2mvyz1 = vzy2_m_vyz1*vzy2_m_vyz1;
		double vzy1mvyz2_mu_vzy1mvyz2 = vzy1_m_vyz2*vzy1_m_vyz2;
		double vzy2mvyz2_mu_vzy2mvyz2 = vzy2_m_vyz2*vzy2_m_vyz2;
		double vzx1mvxz1_mu_vzx1mvxz1 = vzx1_m_vxz1*vzx1_m_vxz1;
		double vzx2mvxz1_mu_vzx2mvxz1 = vzx2_m_vxz1*vzx2_m_vxz1;
		double vzx1mvxz2_mu_vzx1mvxz2 = vzx1_m_vxz2*vzx1_m_vxz2;
		double vzx2mvxz2_mu_vzx2mvxz2 = vzx2_m_vxz2*vzx2_m_vxz2;
		double Log111 = log(vxy1mvyx1_mu_vxy1mvyx1 + vzy1mvyz1_mu_vzy1mvyz1 + vzx1mvxz1_mu_vzx1mvxz1);
		double Log211 = log(vxy1mvyx2_mu_vxy1mvyx2 + vzy1mvyz1_mu_vzy1mvyz1 + vzx2mvxz1_mu_vzx2mvxz1);
		double Log121 = log(vxy2mvyx1_mu_vxy2mvyx1 + vzy2mvyz1_mu_vzy2mvyz1 + vzx1mvxz1_mu_vzx1mvxz1);
		double Log221 = log(vxy2mvyx2_mu_vxy2mvyx2 + vzy2mvyz1_mu_vzy2mvyz1 + vzx2mvxz1_mu_vzx2mvxz1);
		double Log112 = log(vxy1mvyx1_mu_vxy1mvyx1 + vzy1mvyz2_mu_vzy1mvyz2 + vzx1mvxz2_mu_vzx1mvxz2);
		double Log212 = log(vxy1mvyx2_mu_vxy1mvyx2 + vzy1mvyz2_mu_vzy1mvyz2 + vzx2mvxz2_mu_vzx2mvxz2);
		double Log122 = log(vxy2mvyx1_mu_vxy2mvyx1 + vzy2mvyz2_mu_vzy2mvyz2 + vzx1mvxz2_mu_vzx1mvxz2);
		double Log222 = log(vxy2mvyx2_mu_vxy2mvyx2 + vzy2mvyz2_mu_vzy2mvyz2 + vzx2mvxz2_mu_vzx2mvxz2);
	
		if(J_is_Zero)
		{
			double CommonLogTermZ = 0.5*(P1.z*(Log111-Log211+Log221-Log121) + P2.z*(Log212-Log112+Log122-Log222));
			double CommonLogTermX = 0.5*(P1.x*(Log111-Log121+Log122-Log112) + P2.x*(Log221-Log211+Log212-Log222));
			double CommonLogTermY = 0.5*(P1.y*(Log111-Log211+Log212-Log112) + P2.y*(Log221-Log121+Log122-Log222));
/* Fxyx*/	F.Str0.x = One_d_vxvxpvzvz*(vzx1_m_vxz1*SumAtanYx1z1 - vzx2_m_vxz1*SumAtanYx2z1 - vzx1_m_vxz2*SumAtanYx1z2 + vzx2_m_vxz2*SumAtanYx2z2)
					 + 0.5*V.y*One_d_vxvxpvzvz*((vzz1+vxx1)*(Log121-Log111) + (vzz1+vxx2)*(Log211-Log221) + (vzz2+vxx1)*(Log112-Log122) + (vzz2+vxx2)*(Log222-Log212))
					 + CommonLogTermY;
/* Fxzx*/	F.Str0.y = One_d_vxvxpvyvy*(-vxy1_m_vyx1*SumAtanZx1y1 + vxy1_m_vyx2*SumAtanZx2y1 + vxy2_m_vyx1*SumAtanZx1y2 - vxy2_m_vyx2*SumAtanZx2y2)
					 + 0.5*V.z*One_d_vxvxpvyvy*((vyy1+vxx1)*(Log112-Log111) + (vyy1+vxx2)*(Log211-Log212) + (vyy2+vxx1)*(Log121-Log122) + (vyy2+vxx2)*(Log222-Log221))
					 + CommonLogTermZ;
/* -Fyz0*/	F.Str0.z = One_d_vxvxpvyvy_d_vxvxpvzvz*((vzvz_m_vyvy*vxx1+vxvx_p_vzvz*vyy1)*SumAtanZx1y1 - (vzvz_m_vyvy*vxx1+vxvx_p_vzvz*vyy2)*SumAtanZx1y2 - (vzvz_m_vyvy*vxx2+vxvx_p_vzvz*vyy1)*SumAtanZx2y1 + (vzvz_m_vyvy*vxx2+vxvx_p_vzvz*vyy2)*SumAtanZx2y2)
					 + V.z*One_d_vxvxpvzvz*(P1.z*(SumAtanYx1z1-SumAtanYx2z1) + P2.z*(SumAtanYx2z2-SumAtanYx1z2))
					 - 0.5*V.z*One_d_vxvxpvyvy_d_vxvxpvzvz*((vyvy_p_vzvz_p_2vxvx*vyx1-vxvx_p_vzvz*vxy1)*(Log112-Log111) + (vyvy_p_vzvz_p_2vxvx*vyx1-vxvx_p_vzvz*vxy2)*(Log121-Log122) + (vyvy_p_vzvz_p_2vxvx*vyx2-vxvx_p_vzvz*vxy1)*(Log211-Log212) + (vyvy_p_vzvz_p_2vxvx*vyx2-vxvx_p_vzvz*vxy2)*(Log222-Log221))
					 - V.x*V.y*One_d_vxvxpvzvz*CommonLogTermZ
					 - (Pi*vzvz*One_d_vxvxpvzvz/V.x)*(P1.x*Step(vzx1/V.x-P1.z)*Step(P2.z-vzx1/V.x)*(Sign(vxy1_m_vyx1)-Sign(vxy2_m_vyx1)) + P2.x*Step(vzx2/V.x-P1.z)*Step(P2.z-vzx2/V.x)*(Sign(vxy2_m_vyx2)-Sign(vxy1_m_vyx2)));
/* Fxyy*/	F.Str1.x = One_d_vyvypvzvz*(vzy1_m_vyz1*SumAtanXy1z1 - vzy2_m_vyz1*SumAtanXy2z1 - vzy1_m_vyz2*SumAtanXy1z2 + vzy2_m_vyz2*SumAtanXy2z2)
					 + 0.5*V.x*One_d_vyvypvzvz*((vzz1+vyy1)*(Log211-Log111) + (vzz1+vyy2)*(Log121-Log221) + (vzz2+vyy1)*(Log112-Log212) + (vzz2+vyy2)*(Log222-Log122))
					 + CommonLogTermX;
/* -Fxz0*/	F.Str1.y =-One_d_vxvxpvyvy_d_vyvypvzvz*((vzvz_m_vxvx*vyy1+vyvy_p_vzvz*vxx1)*SumAtanZx1y1 - (vzvz_m_vxvx*vyy1+vyvy_p_vzvz*vxx2)*SumAtanZx2y1 - (vzvz_m_vxvx*vyy2+vyvy_p_vzvz*vxx1)*SumAtanZx1y2 + (vzvz_m_vxvx*vyy2+vyvy_p_vzvz*vxx2)*SumAtanZx2y2)
					 + V.z*One_d_vyvypvzvz*(P1.z*(SumAtanXy1z1-SumAtanXy2z1) + P2.z*(SumAtanXy2z2-SumAtanXy1z2))
					 - 0.5*V.z*One_d_vxvxpvyvy_d_vyvypvzvz*((vxvx_p_vzvz_p_2vyvy*vxy1-vyvy_p_vzvz*vyx1)*(Log112-Log111) + (vxvx_p_vzvz_p_2vyvy*vxy1-vyvy_p_vzvz*vyx2)*(Log211-Log212) + (vxvx_p_vzvz_p_2vyvy*vxy2-vyvy_p_vzvz*vyx1)*(Log121-Log122) + (vxvx_p_vzvz_p_2vyvy*vxy2-vyvy_p_vzvz*vyx2)*(Log222-Log221))
					 - V.x*V.y*One_d_vyvypvzvz*CommonLogTermZ
					 - (Pi*vzvz*One_d_vyvypvzvz/V.y)*(P1.y*Step(vzy1/V.y-P1.z)*Step(P2.z-vzy1/V.y)*(Sign(vxy1_m_vyx2)-Sign(vxy1_m_vyx1)) + P2.y*Step(vzy2/V.y-P1.z)*Step(P2.z-vzy2/V.y)*(Sign(vxy2_m_vyx1)-Sign(vxy2_m_vyx2)));
/* Fyzy*/	F.Str1.z = F.Str0.y;
/* -Fxy0*/	F.Str2.x =-One_d_vxvxpvzvz_d_vyvypvzvz*((vxvx_m_vyvy*vzz1+vxvx_p_vzvz*vyy1)*SumAtanXy1z1 - (vxvx_m_vyvy*vzz1+vxvx_p_vzvz*vyy2)*SumAtanXy2z1 - (vxvx_m_vyvy*vzz2+vxvx_p_vzvz*vyy1)*SumAtanXy1z2 + (vxvx_m_vyvy*vzz2+vxvx_p_vzvz*vyy2)*SumAtanXy2z2)
					 - V.x*One_d_vxvxpvzvz*(P1.x*(SumAtanYx1z1-SumAtanYx1z2) + P2.x*(SumAtanYx2z2-SumAtanYx2z1))
					 - 0.5*V.x*One_d_vxvxpvzvz_d_vyvypvzvz*((vxvx_p_vyvy_p_2vzvz*vyz1-vxvx_p_vzvz*vzy1)*(Log211-Log111) + (vxvx_p_vyvy_p_2vzvz*vyz1-vxvx_p_vzvz*vzy2)*(Log121-Log221) + (vxvx_p_vyvy_p_2vzvz*vyz2-vxvx_p_vzvz*vzy1)*(Log112-Log212) + (vxvx_p_vyvy_p_2vzvz*vyz2-vxvx_p_vzvz*vzy2)*(Log222-Log122))
					 - V.y*V.z*One_d_vxvxpvzvz*CommonLogTermX
					 - (Pi*vxvx*One_d_vxvxpvzvz/V.z)*(P1.z*Step(vxz1/V.z-P1.x)*Step(P2.x-vxz1/V.z)*(Sign(vzy1_m_vyz1)-Sign(vzy2_m_vyz1)) + P2.z*Step(vxz2/V.z-P1.x)*Step(P2.x-vxz2/V.z)*(Sign(vzy2_m_vyz2)-Sign(vzy1_m_vyz2)));
/* Fxzz*/	F.Str2.y = F.Str1.x;
/* Fyzz*/	F.Str2.z = F.Str0.x;
		}
		else
		{
			double vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx = V.y*One_d_vxvxpvzvz*vyvy_p_vzvz_p_2vxvx;
			double vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy = V.x*One_d_vyvypvzvz*vxvx_p_vzvz_p_2vyvy;
			double vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz = V.y*One_d_vxvxpvzvz*vxvx_p_vyvy_p_2vzvz;
			double x1x1 = P1.x*P1.x;
			double x2x2 = P2.x*P2.x;
			double y1y1 = P1.y*P1.y;
			double y2y2 = P2.y*P2.y;
			double z1z1 = P1.z*P1.z;
			double z2z2 = P2.z*P2.z;
			double TwoGx, TwoGy, TwoGz;

			TwoGx =-One_d_vxvxpvyvy_d_vxvxpvzvz*((vxvx_p_vzvz*(vxy1_m_vyx1-vyx1)*P1.y-vzvz_m_vyvy*vxx1*P1.x)*SumAtanZx1y1 - (vxvx_p_vzvz*(vxy2_m_vyx1-vyx1)*P2.y-vzvz_m_vyvy*vxx1*P1.x)*SumAtanZx1y2 - (vxvx_p_vzvz*(vxy1_m_vyx2-vyx2)*P1.y-vzvz_m_vyvy*vxx2*P2.x)*SumAtanZx2y1 + (vxvx_p_vzvz*(vxy2_m_vyx2-vyx2)*P2.y-vzvz_m_vyvy*vxx2*P2.x)*SumAtanZx2y2)
				  + One_d_vxvxpvzvz*(P1.z*(vzx1_m_vxz1+vzx1)*SumAtanYx1z1 - P2.z*(vzx1_m_vxz2+vzx1)*SumAtanYx1z2 - P1.z*(vzx2_m_vxz1+vzx2)*SumAtanYx2z1 + P2.z*(vzx2_m_vxz2+vzx2)*SumAtanYx2z2)
				  + 0.5*V.y*One_d_vxvxpvzvz*((2.*vxx1+vzz1)*P1.z*(Log121-Log111) + (2.*vxx1+vzz2)*P2.z*(Log112-Log122) + (2.*vxx2+vzz1)*P1.z*(Log211-Log221) + (2.*vxx2+vzz2)*P2.z*(Log222-Log212))
				  + 0.5*V.z*One_d_vxvxpvyvy*((vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx*x1x1-(2.*vxx1+vyy1)*P1.y)*(Log111-Log112) + (-vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx*x1x1+(2.*vxx1+vyy2)*P2.y)*(Log121-Log122) + (vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx*x2x2-(2.*vxx2+vyy1)*P1.y)*(Log212-Log211) + (-vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx*x2x2+(2.*vxx2+vyy2)*P2.y)*(Log222-Log221))
				  + P1.y*P1.z*(Log111-Log211) + P1.y*P2.z*(Log212-Log112) + P2.y*P1.z*(Log221-Log121) + P2.y*P2.z*(Log122-Log222)
				  + (Pi*V.z*V.z*One_d_vxvxpvzvz/V.x)*(x1x1*Step(vzx1/V.x-P1.z)*Step(P2.z-vzx1/V.x)*(Sign(vxy2_m_vyx1)-Sign(vxy1_m_vyx1)) + x2x2*Step(vzx2/V.x-P1.z)*Step(P2.z-vzx2/V.x)*(Sign(vxy1_m_vyx2)-Sign(vxy2_m_vyx2)));
			TwoGy = One_d_vxvxpvyvy_d_vyvypvzvz*(-(vyvy_p_vzvz*(vxy1_m_vyx1+vxy1)*P1.x+vzvz_m_vxvx*vyy1*P1.y)*SumAtanZx1y1 + (vyvy_p_vzvz*(vxy1_m_vyx2+vxy1)*P2.x+vzvz_m_vxvx*vyy1*P1.y)*SumAtanZx2y1 + (vyvy_p_vzvz*(vxy2_m_vyx1+vxy2)*P1.x+vzvz_m_vxvx*vyy2*P2.y)*SumAtanZx1y2 - (vyvy_p_vzvz*(vxy2_m_vyx2+vxy2)*P2.x+vzvz_m_vxvx*vyy2*P2.y)*SumAtanZx2y2)
				  + One_d_vyvypvzvz*(P1.z*(vzy1_m_vyz1+vzy1)*SumAtanXy1z1 - P2.z*(vzy1_m_vyz2+vzy1)*SumAtanXy1z2 - P1.z*(vzy2_m_vyz1+vzy2)*SumAtanXy2z1 + P2.z*(vzy2_m_vyz2+vzy2)*SumAtanXy2z2)
				  + 0.5*V.x*One_d_vyvypvzvz*((2.*vyy1+vzz1)*P1.z*(Log211-Log111) + (2.*vyy1+vzz2)*P2.z*(Log112-Log212) + (2.*vyy2+vzz1)*P1.z*(Log121-Log221) + (2.*vyy2+vzz2)*P2.z*(Log222-Log122))
				  + 0.5*V.z*One_d_vxvxpvyvy*((vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy*y1y1-(2.*vyy1+vxx1)*P1.x)*(Log111-Log112) + (-vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy*y1y1+(2.*vyy1+vxx2)*P2.x)*(Log211-Log212) + (vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy*y2y2-(2.*vyy2+vxx1)*P1.x)*(Log122-Log121) + (-vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy*y2y2+(2.*vyy2+vxx2)*P2.x)*(Log222-Log221))
				  + P1.x*P1.z*(Log111-Log121) + P1.x*P2.z*(Log122-Log112) + P2.x*P1.z*(Log221-Log211) + P2.x*P2.z*(Log212-Log222)
				  + (Pi*V.z*V.z*One_d_vyvypvzvz/V.y)*(y1y1*Step(vzy1/V.y-P1.z)*Step(P2.z-vzy1/V.y)*(Sign(vxy1_m_vyx1)-Sign(vxy1_m_vyx2)) + y2y2*Step(vzy2/V.y-P1.z)*Step(P2.z-vzy2/V.y)*(Sign(vxy2_m_vyx2)-Sign(vxy2_m_vyx1)));
			TwoGz = One_d_vxvxpvzvz_d_vyvypvzvz*((vxvx_p_vzvz*(vzy1_m_vyz1-vyz1)*P1.y-vxvx_m_vyvy*vzz1*P1.z)*SumAtanXy1z1 - (vxvx_p_vzvz*(vzy2_m_vyz1-vyz1)*P2.y-vxvx_m_vyvy*vzz1*P1.z)*SumAtanXy2z1 - (vxvx_p_vzvz*(vzy1_m_vyz2-vyz2)*P1.y-vxvx_m_vyvy*vzz2*P2.z)*SumAtanXy1z2 + (vxvx_p_vzvz*(vzy2_m_vyz2-vyz2)*P2.y-vxvx_m_vyvy*vzz2*P2.z)*SumAtanXy2z2)
				  + One_d_vxvxpvzvz*(P1.x*(vzx1_m_vxz1-vxz1)*SumAtanYx1z1 - P2.x*(vzx2_m_vxz1-vxz1)*SumAtanYx2z1 - P1.x*(vzx1_m_vxz2-vxz2)*SumAtanYx1z2 + P2.x*(vzx2_m_vxz2-vxz2)*SumAtanYx2z2)
				  + 0.5*V.y*One_d_vxvxpvzvz*((2.*vzz1+vxx1)*P1.x*(Log121-Log111) + (2.*vzz1+vxx2)*P2.x*(Log211-Log221) + (2.*vzz2+vxx1)*P1.x*(Log112-Log122) + (2.*vzz2+vxx2)*P2.x*(Log222-Log212))
				  + 0.5*V.x*One_d_vyvypvzvz*((vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz*z1z1-(2.*vzz1+vyy1)*P1.y)*(Log111-Log211) + (-vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz*z1z1+(2.*vzz1+vyy2)*P2.y)*(Log121-Log221) + (vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz*z2z2-(2.*vzz2+vyy1)*P1.y)*(Log212-Log112) + (-vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz*z2z2+(2.*vzz2+vyy2)*P2.y)*(Log222-Log122))
				  + P1.x*P1.y*(Log111-Log112) + P2.x*P1.y*(Log212-Log211) + P1.x*P2.y*(Log122-Log121) + P2.x*P2.y*(Log221-Log222)
				  + (Pi*V.x*V.x*One_d_vxvxpvzvz/V.z)*(z1z1*Step(vxz1/V.z-P1.x)*Step(P2.x-vxz1/V.z)*(Sign(vzy2_m_vyz1)-Sign(vzy1_m_vyz1)) + z2z2*Step(vxz2/V.z-P1.x)*Step(P2.x-vxz2/V.z)*(-Sign(vzy2_m_vyz2)+Sign(vzy1_m_vyz2)));
			G.x = 0.5*TwoGx; G.y = 0.5*TwoGy; G.z = 0.5*TwoGz;
		}
	}

FinalDefinitionOfFieldIntegrals:
	if(J_is_Zero)
	{
		const double ConstForM = 1./2./Pi;
		TVector3d ConByM = ConstForM * Magn;
		if(FieldPtr->FieldKey.Ib_)
		{
			TVector3d BufIb(-(F.Str2.x+F.Str1.y)*ConByM.x + F.Str1.z*ConByM.y + F.Str2.z*ConByM.z,
							 F.Str0.y*ConByM.x - (F.Str2.x+F.Str0.z)*ConByM.y + F.Str2.y*ConByM.z,
							 F.Str0.x*ConByM.x + F.Str1.x*ConByM.y - (F.Str1.y+F.Str0.z)*ConByM.z);
			FieldPtr->Ib += BufIb;
		}
		if(FieldPtr->FieldKey.Ih_)
		{
			TVector3d BufIh(F.Str0.z*ConByM.x + F.Str0.y*ConByM.y + F.Str0.x*ConByM.z, 
							F.Str1.z*ConByM.x + F.Str1.y*ConByM.y + F.Str1.x*ConByM.z,
							F.Str2.z*ConByM.x + F.Str2.y*ConByM.y + F.Str2.x*ConByM.z);
			FieldPtr->Ih += BufIh;
		}
	}
	else
	{
		const double ConstForJ = 0.0002;
		TVector3d BufI((J.y*G.z-J.z*G.y)*ConstForJ, (J.z*G.x-J.x*G.z)*ConstForJ, (J.x*G.y-J.y*G.x)*ConstForJ);
		if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufI;
		if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufI;
	}
}

//-------------------------------------------------------------------------

void radTRecMag::B_intUtilSpecCaseZeroVxVy(const TVector3d& P1, const TVector3d& P2, short J_is_Zero, TMatrix3d& F, TVector3d& G)
{
	double z2_m_z1 = P2.z - P1.z;

	const double Pi = 3.141592653589793238;
	double PiMult1, PiMult2, PiMult3;
	PiMult1 = PiMult2 = PiMult3 = 0.;

	double x1x1 = P1.x*P1.x;
	double y1y1 = P1.y*P1.y;
	double x2x2 = P2.x*P2.x;
	double y2y2 = P2.y*P2.y;

	if(J_is_Zero)
	{
		F.Str0.x = F.Str1.x = F.Str2.x = F.Str2.y = F.Str2.z = 0.;
		double SumAtan1 = atan(TransAtans(TransAtans(-P2.x/P1.y, P1.x/P1.y, PiMult1), TransAtans(P2.x/P2.y, -P1.x/P2.y, PiMult2), PiMult3)) 
			            + Pi*(PiMult1+PiMult2+PiMult3);
		F.Str1.y = -z2_m_z1*SumAtan1;
		double SumAtan2 = atan(TransAtans(TransAtans(-P2.y/P1.x, P1.y/P1.x, PiMult1), TransAtans(P2.y/P2.x, -P1.y/P2.x, PiMult2), PiMult3)) 
			            + Pi*(PiMult1+PiMult2+PiMult3);
		F.Str0.z = -z2_m_z1*SumAtan2;
		F.Str0.y = F.Str1.z = 0.5*z2_m_z1*log(((x1x1+y2y2)*(x2x2+y1y1))/((x1x1+y1y1)*(x2x2+y2y2)));
	}
	else
	{
		double SumAtan_x1 = atan(TransAtans(-P1.y/P1.x, P2.y/P1.x, PiMult1)) + Pi*PiMult1;
		double SumAtan_x2 = atan(TransAtans(-P2.y/P2.x, P1.y/P2.x, PiMult1)) + Pi*PiMult1;
		G.x = z2_m_z1*(P1.x*SumAtan_x1 + P2.x*SumAtan_x2 + 0.5*(P1.y*log((x2x2+y1y1)/(x1x1+y1y1)) + P2.y*log((x1x1+y2y2)/(x2x2+y2y2))));
		double SumAtan_y1 = atan(TransAtans(-P1.x/P1.y, P2.x/P1.y, PiMult1)) + Pi*PiMult1;
		double SumAtan_y2 = atan(TransAtans(-P2.x/P2.y, P1.x/P2.y, PiMult1)) + Pi*PiMult1;
		G.y = z2_m_z1*(P1.y*SumAtan_y1 + P2.y*SumAtan_y2 + 0.5*(P1.x*log((x1x1+y2y2)/(x1x1+y1y1)) + P2.x*log((x2x2+y1y1)/(x2x2+y2y2))));
		G.z = 0.;
	}
}

//-------------------------------------------------------------------------

void radTRecMag::FunForOuterIntAtSurfInt(double Arg, TVector3d* VectArray)
{
	const double PrecEnhFact = 1.; // Don't make it >1 : it's dangerous for convergence of outer itegral !
	double* OuterIntPrecArray = SurfIntDataPtr->Field.ShapeIntDataPtr->AbsPrecArray;

	double SmallPositive = 1.E-10;

	if(SurfIntDataPtr->SurfBoundInd==1 || SurfIntDataPtr->SurfBoundInd==2)
	{
		for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++)
			(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[i] = PrecEnhFact*OuterIntPrecArray[i]/Dimensions.y;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen] = CentrPoint.x - 0.5*Dimensions.x + SmallPositive;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen + 1] = CentrPoint.x + 0.5*Dimensions.x;
		SurfIntDataPtr->PointOnSurface.y = Arg;
	}
	else if(SurfIntDataPtr->SurfBoundInd==3 || SurfIntDataPtr->SurfBoundInd==4)
	{
		for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++)
			(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[i] = PrecEnhFact*OuterIntPrecArray[i]/Dimensions.z;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen] = CentrPoint.x - 0.5*Dimensions.x + SmallPositive;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen + 1] = CentrPoint.x + 0.5*Dimensions.x;
		SurfIntDataPtr->PointOnSurface.z = Arg;
	}
	else if(SurfIntDataPtr->SurfBoundInd==5 || SurfIntDataPtr->SurfBoundInd==6)
	{
		for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++)
			(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[i] = PrecEnhFact*OuterIntPrecArray[i]/Dimensions.z;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen] = CentrPoint.y - 0.5*Dimensions.y + SmallPositive;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen + 1] = CentrPoint.y + 0.5*Dimensions.y;
		SurfIntDataPtr->PointOnSurface.z = Arg;
	}
	FormalOneFoldInteg(this, &radTRecMag::FunForInnerIntAtSurfInt, SurfIntDataPtr->IntegrandLen, 
					   SurfIntDataPtr->InnerAbsPrecAndLimitsArray, 
					   SurfIntDataPtr->InnerElemCompNotFinished, SurfIntDataPtr->InnerIntegVal);

	for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++) VectArray[i] = ((SurfIntDataPtr->InnerIntegVal)[0])[i];
}

//-------------------------------------------------------------------------

void radTRecMag::IntOverSurf(radTField* FieldPtr)
{
	int LenVal = FieldPtr->ShapeIntDataPtr->IntegrandLength;
	int LenValp2 = LenVal+2;
	TVector3d ZeroVect(0.,0.,0.);

	TVector3d* InnerIntegVal[6]; 
	TVector3d* OuterIntegVal[6];

	short* InnerElemCompNotFinished;
	short* OuterElemCompNotFinished;
	double* InnerAbsPrecAndLimitsArray;
	double* OuterAbsPrecAndLimitsArray;
	TVector3d* LocalVectArray;

	int j;
	for(j=0; j<6; j++)
	{
		InnerIntegVal[j] = new TVector3d[LenVal];
		OuterIntegVal[j] = new TVector3d[LenVal];
	}
	InnerElemCompNotFinished = new short[LenVal];
	OuterElemCompNotFinished = new short[LenVal];
	InnerAbsPrecAndLimitsArray = new double[LenValp2];
	OuterAbsPrecAndLimitsArray = new double[LenValp2];
	LocalVectArray = new TVector3d[LenVal];
	SurfIntDataPtr = new radTParallelepSurfIntData();

	SurfIntDataPtr->IntegrandLen = LenVal;
	SurfIntDataPtr->IntegrandFunPtr = FieldPtr->ShapeIntDataPtr->IntegrandFunPtr;
	SurfIntDataPtr->InnerAbsPrecAndLimitsArray = InnerAbsPrecAndLimitsArray;
	SurfIntDataPtr->InnerElemCompNotFinished = InnerElemCompNotFinished;
	SurfIntDataPtr->InnerIntegVal = InnerIntegVal;
	
	SurfIntDataPtr->Field = *FieldPtr;
	radTStructForShapeInt LocShapeIntData = *(FieldPtr->ShapeIntDataPtr);

	TVector3d* InputFieldPtrVectArrayPtr = FieldPtr->ShapeIntDataPtr->VectArray;

	LocShapeIntData.VectArray = LocalVectArray;
	SurfIntDataPtr->Field.ShapeIntDataPtr = &LocShapeIntData;

	int i;
	for(i=0; i<LenVal; i++)
	{
		OuterAbsPrecAndLimitsArray[i] = (FieldPtr->ShapeIntDataPtr->AbsPrecArray)[i];
	}
	TVector3d HalfDim = 0.5*Dimensions;

	TVector3d* OutVectArray = InputFieldPtrVectArrayPtr;

	double SmallPositive = 1.E-10;

//For lower and upper bounds
	OuterAbsPrecAndLimitsArray[LenVal] = CentrPoint.y - HalfDim.y + SmallPositive;
	OuterAbsPrecAndLimitsArray[LenVal+1] = CentrPoint.y + HalfDim.y;
//Integration over lower bound
	SurfIntDataPtr->SurfBoundInd = 1;
	SurfIntDataPtr->PointOnSurface.z = CentrPoint.z - HalfDim.z + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,0.,-1.);
	FormalOneFoldInteg(this, &radTRecMag::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];
//Integration over upper bound
	SurfIntDataPtr->SurfBoundInd = 2;
	SurfIntDataPtr->PointOnSurface.z = CentrPoint.z + HalfDim.z + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,0.,1.);
	FormalOneFoldInteg(this, &radTRecMag::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];

//For left, right, back and front bounds
	OuterAbsPrecAndLimitsArray[LenVal] = CentrPoint.z - HalfDim.z + SmallPositive;
	OuterAbsPrecAndLimitsArray[LenVal+1] = CentrPoint.z + HalfDim.z;
//Integration over left bound
	SurfIntDataPtr->SurfBoundInd = 3;
	SurfIntDataPtr->PointOnSurface.y = CentrPoint.y - HalfDim.y + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,-1.,0.);
	FormalOneFoldInteg(this, &radTRecMag::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];
//Integration over right bound
	SurfIntDataPtr->SurfBoundInd = 4;
	SurfIntDataPtr->PointOnSurface.y = CentrPoint.y + HalfDim.y + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,1.,0.);
	FormalOneFoldInteg(this, &radTRecMag::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];
//Integration over back bound
	SurfIntDataPtr->SurfBoundInd = 5;
	SurfIntDataPtr->PointOnSurface.x = CentrPoint.x - HalfDim.x + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(-1.,0.,0.);
	FormalOneFoldInteg(this, &radTRecMag::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];
//Integration over right bound
	SurfIntDataPtr->SurfBoundInd = 6;
	SurfIntDataPtr->PointOnSurface.x = CentrPoint.x + HalfDim.x + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(1.,0.,0.);
	FormalOneFoldInteg(this, &radTRecMag::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];

// Delete everything created above:
	delete[] InnerElemCompNotFinished;
	delete[] OuterElemCompNotFinished;
	delete[] InnerAbsPrecAndLimitsArray;
	delete[] OuterAbsPrecAndLimitsArray;
	for(j=0; j<6; j++)
	{
		delete[] InnerIntegVal[j];
		delete[] OuterIntegVal[j];
	}
	delete[] LocalVectArray;
	delete SurfIntDataPtr;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTRecMag::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTRecMagGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

int radTRecMag::SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions) 
{
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	const double ZeroTol = 1.E-10;
	if((!SubdivideCoils) && (J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero)) return 1;

	double LocSubdivArray[15];
	for(int jj=0; jj<15; jj++) LocSubdivArray[jj] = SubdivArray[jj];

	radTSend Send;
	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty())) 
	{
		if(J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero) { Send.ErrorMessage("Radia::Error108"); return 0;}
		char ConversionToPolyhedronIsNeeded = 0;
		CheckAxesExchangeForSubdInLabFrame(LocSubdivArray, ConversionToPolyhedronIsNeeded);
		if(ConversionToPolyhedronIsNeeded)
		{
			radThg& NewHandle = In_hg;
			radThg OldHandle = In_hg;
			if(!((radTg3d*)(OldHandle.rep))->ConvertToPolyhedron(NewHandle, radPtr, pSubdivOptions->PutNewStuffIntoGenCont)) // PutNewStuffIntoGenCont is only necessary for "to replace in all groups"
			{
				Send.ErrorMessage("Radia::Error108"); return 0;
			}
			radThg OldNewHandle = NewHandle;
			int SubdOK = ((radTPolyhedron*)(OldNewHandle.rep))->SubdivideItself(LocSubdivArray, NewHandle, radPtr, pSubdivOptions);
			if(!SubdOK) return 0;
			return 1;
		}
	}
	double &kx = LocSubdivArray[0], &ky = LocSubdivArray[2], &kz = LocSubdivArray[4];
	double &qx = LocSubdivArray[1], &qy = LocSubdivArray[3], &qz = LocSubdivArray[5];
	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		kx = (kx < Dimensions.x)? Round(Dimensions.x/kx) : 1.;
		ky = (ky < Dimensions.y)? Round(Dimensions.y/ky) : 1.;
		kz = (kz < Dimensions.z)? Round(Dimensions.z/kz) : 1.;
	}
	if((fabs(kx-1.)<ZeroTol) && (fabs(ky-1.)<ZeroTol) && (fabs(kz-1.)<ZeroTol)) return 1;

	radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedRecMag(this, LocSubdivArray);
	radThg NewHandle(GroupInPlaceOfThisPtr);

	const double AbsZeroTol = 5.E-12;
	double q0x = (fabs(kx-1.)>AbsZeroTol)? pow(qx, 1./(kx-1.)) : qx;
	double q0y = (fabs(ky-1.)>AbsZeroTol)? pow(qy, 1./(ky-1.)) : qy;
	double q0z = (fabs(kz-1.)>AbsZeroTol)? pow(qz, 1./(kz-1.)) : qz;
	double BufX = qx*q0x - 1., BufY = qy*q0y - 1., BufZ = qz*q0z - 1.;
	double a1x = (fabs(BufX) > AbsZeroTol)? Dimensions.x*(q0x - 1.)/BufX : Dimensions.x/kx;
	double a1y = (fabs(BufY) > AbsZeroTol)? Dimensions.y*(q0y - 1.)/BufY : Dimensions.y/ky;
	double a1z = (fabs(BufZ) > AbsZeroTol)? Dimensions.z*(q0z - 1.)/BufZ : Dimensions.z/kz;

	TVector3d InitNewDims(a1x, a1y, a1z);
	TVector3d NewDims = InitNewDims;
	TVector3d InitNewCenPoi = CentrPoint - 0.5*(Dimensions - NewDims);
	TVector3d NewCenPoi = InitNewCenPoi;

	short NewFacesState[6], ParentFacesState[6];
	ListFacesInternalAfterCut(ParentFacesState);

	int kxInt = int(kx), kyInt = int(ky), kzInt = int(kz);
	int kx_mi_1 = kxInt-1, ky_mi_1 = kyInt-1, kz_mi_1 = kzInt-1;

	int NewStuffCounter = 0;
	for(int ix=0; ix<kxInt; ix++)
	{
		NewFacesState[0] = NewFacesState[1] = 1;
		if(ix==0) NewFacesState[0] = ParentFacesState[0];
		if(ix==kx_mi_1) NewFacesState[1] = ParentFacesState[1];

		for(int iy=0; iy<kyInt; iy++)
		{
			NewFacesState[2] = NewFacesState[3] = 1;
			if(iy==0) NewFacesState[2] = ParentFacesState[2];
			if(iy==ky_mi_1) NewFacesState[3] = ParentFacesState[3];

			for(int iz=0; iz<kzInt; iz++)
			{
				NewFacesState[4] = NewFacesState[5] = 1;
				if(iz==0) NewFacesState[4] = ParentFacesState[4];
				if(iz==kz_mi_1) NewFacesState[5] = ParentFacesState[5];

				radTRecMag* RecMagPtr = new radTRecMag(NewCenPoi, NewDims, Magn, J, MaterHandle, J_IsNotZero);
				if(RecMagPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
				radThg hg(RecMagPtr);
				if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(hg), hg);
				else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, hg);

				RecMagPtr->SetFacesInternalAfterCut(NewFacesState);

				NewCenPoi.z += 0.5*NewDims.z;
				NewDims.z *= q0z;
				NewCenPoi.z += 0.5*NewDims.z;
			}
			NewCenPoi.y += 0.5*NewDims.y;
			NewDims.y *= q0y;
			NewCenPoi.y += 0.5*NewDims.y;

			NewCenPoi.z = InitNewCenPoi.z;
			NewDims.z = InitNewDims.z;
		}
		NewCenPoi.x += 0.5*NewDims.x;
		NewDims.x *= q0x;
		NewCenPoi.x += 0.5*NewDims.x;

		NewCenPoi.y = InitNewCenPoi.y;
		NewDims.y = InitNewDims.y;
	}
	In_hg = NewHandle;
	return 1;
}

//-------------------------------------------------------------------------

void radTRecMag::Dump(std::ostream& o, int ShortSign) // Porting
{
	radTg3dRelax::Dump(o);
	DumpPureObjInfo(o, ShortSign);
	if(ShortSign==1) return;

	if(!J_IsNotZero) DumpMaterApplied(o);
	DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTRecMag::DumpPureObjInfo(std::ostream& o, int ShortSign)
{
	if(J_IsNotZero) 
	{ 
		o << "Current carrying: ";
		o << "RecCur";
	}
	else 
	{ 
		o << "Relaxable: ";
		o << "RecMag";
	}
	if(ShortSign==1) return;
	o << endl;
	o << "   {x,y,z}= {" << CentrPoint.x << ','
						 << CentrPoint.y << ','
						 << CentrPoint.z << "}" << endl
	  << "   {wx,wy,wz}= {" << Dimensions.x << ','
							<< Dimensions.y << ','
							<< Dimensions.z << "}" << endl;
	if(J_IsNotZero) 
	{ 
	  o << "   {jx,jy,jz}= {" << J.x << ','
						      << J.y << ','
							  << J.z << "}";
	}
	else
	{
	  o << "   {mx,my,mz}= {" << Magn.x << ','
							  << Magn.y << ','
							  << Magn.z << "}";
	}
}

//-------------------------------------------------------------------------

void radTRecMag::DumpBin_RecMag(CAuxBinStrVect& oStr)
{
	//radTParallelepSurfIntData* SurfIntDataPtr; /is not dumped
	//TVector3d Dimensions;
	oStr << Dimensions;
	
	//TVector3d J;
	oStr << J;

	//short J_IsNotZero;
	oStr << J_IsNotZero;

	//short InternalFacesAfterCut;
	oStr << InternalFacesAfterCut;
}

//-------------------------------------------------------------------------

void radTRecMag::DumpBinParse_RecMag(CAuxBinStrVect& inStr)
{
	//radTParallelepSurfIntData* SurfIntDataPtr; /is not dumped
	//TVector3d Dimensions;
	inStr >> Dimensions;

	//TVector3d J;
	inStr >> J;

	//short J_IsNotZero;
	inStr >> J_IsNotZero;

	//short InternalFacesAfterCut;
	inStr >> InternalFacesAfterCut;
}

//-------------------------------------------------------------------------

//void radTRecMag::DumpBin(CAuxBinStrVect& oStr, radTmhg& mEl, radThg& hg)
void radTRecMag::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//Dumping objects that may be used by this object
	vector<pair<int, int> > vTrfKeys;
	//DumpBin_g3d_TreatTrfs(oStr, mEl, vTrfKeys);
	DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys);
	int matKey=0;
	DumpBin_g3dRelax_TreatMat(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, matKey);

	//int newKey = (int)mEl.size() + 1;
	//mEl[newKey] = hg;
	 
	//Start dumping this object
	//oStr << newKey;

	//elemCount++;
	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)Type_g3d();
	oStr << (char)Type_g3dRelax();
	oStr << (char)Type_RecMag();
	oStr << (char)0;

	//Members of radTg3d
	DumpBin_g3d(oStr, vTrfKeys);

	//Members of radTg3dRelax
	DumpBin_g3dRelax(oStr, matKey);

	//Members of radTRecMag
	DumpBin_RecMag(oStr);
}

//-------------------------------------------------------------------------

int radTRecMag::ConvertToPolyhedron(radThg& In_hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{
	if(J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero) return 1;

	TVector3d hDims = 0.5*Dimensions;
	double xMin = CentrPoint.x - hDims.x, xMax = CentrPoint.x + hDims.x;
	double yMin = CentrPoint.y - hDims.y, yMax = CentrPoint.y + hDims.y;
	double zMin = CentrPoint.z - hDims.z, zMax = CentrPoint.z + hDims.z;

	TVector3d ArrayOfPoints[8];
	ArrayOfPoints[0] = TVector3d(xMin, yMin, zMin);
	ArrayOfPoints[1] = TVector3d(xMax, yMin, zMin);
	ArrayOfPoints[2] = TVector3d(xMax, yMax, zMin);
	ArrayOfPoints[3] = TVector3d(xMin, yMax, zMin);
	ArrayOfPoints[4] = TVector3d(xMin, yMin, zMax);
	ArrayOfPoints[5] = TVector3d(xMax, yMin, zMax);
	ArrayOfPoints[6] = TVector3d(xMax, yMax, zMax);
	ArrayOfPoints[7] = TVector3d(xMin, yMax, zMax);
	int* ArrayOfFaces[6];

	int Face0[] = { 1,5,8,4 }; ArrayOfFaces[0] = Face0;
	int Face1[] = { 2,3,7,6 }; ArrayOfFaces[1] = Face1;
	int Face2[] = { 1,2,6,5 }; ArrayOfFaces[2] = Face2;
	int Face3[] = { 3,4,8,7 }; ArrayOfFaces[3] = Face3;
	int Face4[] = { 4,3,2,1 }; ArrayOfFaces[4] = Face4;
	int Face5[] = { 5,6,7,8 }; ArrayOfFaces[5] = Face5;

	int ArrayOfLengths[] = { 4,4,4,4,4,4 };

	if(ConsiderOnlyWithTrans)
	{
		radTrans ResTransf;
		short SomethingFound = 0;
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
		if(SomethingFound) 
		{
			for(int i=0; i<8; i++)
			{
				ArrayOfPoints[i] = ResTransf.TrPoint(ArrayOfPoints[i]);
			}
		}
	}

	radTSend Send;
	radTPolyhedron* PolyhedronPtr = new radTPolyhedron(ArrayOfPoints, 8, ArrayOfFaces, ArrayOfLengths, 6, Magn);
	if(PolyhedronPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
	PolyhedronPtr->MaterHandle = MaterHandle;
	PolyhedronPtr->IsGroupMember = IsGroupMember;

	PolyhedronPtr->g3dListOfTransform = g3dListOfTransform;
	if(ConsiderOnlyWithTrans) PolyhedronPtr->EraseInnerTransform();
	PolyhedronPtr->ConsiderOnlyWithTrans = 0;

	PolyhedronPtr->HandleAuxCompData = HandleAuxCompData;
	PolyhedronPtr->MessageChar = MessageChar;

	short InternalFacesState[6];
	ListFacesInternalAfterCut(InternalFacesState);
	for(int k=0; k<6; k++)
		PolyhedronPtr->VectHandlePgnAndTrans[k].FaceIsInternalAfterCut = InternalFacesState[k]? true : false;

	radThg NewHandle(PolyhedronPtr);
	In_hg = NewHandle;
	return 1;
}

//-------------------------------------------------------------------------

int radTRecMag::SubdivideItselfByOneSetOfParPlanes(
	TVector3d& InPlanesNormal, TVector3d* InPointsOnCuttingPlanes, int AmOfPieces_mi_1, 
	radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions, radTvhg* pVectOfHgChanged)
{
	if((!pSubdivOptions->SubdivideCoils) && (J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero)) { return 1;}
	radTSend Send;
	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty()) && (J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero)) 
	{
		Send.ErrorMessage("Radia::Error108"); return 0;
	}

	TVector3d PlanesNormal, *PointsOnCuttingPlanes;
	if(!TransferSubdivisionStructToLocalFrame(InPlanesNormal, InPointsOnCuttingPlanes, AmOfPieces_mi_1, pSubdivOptions, PlanesNormal, PointsOnCuttingPlanes)) return 0;

	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double RelTol = RelAbsTol[0];

	short DirNo = -1; //1- "x", 2- "y", 3- "z"
	double AbsNy = Abs(PlanesNormal.y), AbsNz = Abs(PlanesNormal.z);
	if((AbsNy<RelTol) && (AbsNz<RelTol)) DirNo = 1;
	else
	{
		double AbsNx = Abs(PlanesNormal.x);
		if((AbsNx<RelTol) && (AbsNz<RelTol)) DirNo = 2;
		else
		{
			if((AbsNx<RelTol) && (AbsNy<RelTol)) DirNo = 3;
		}
	}
	if(DirNo > 0) 
	{
		radThg hgOld = In_hg;
		int SubdOK = SubdivideItselfByPlanesParToFace(DirNo, PointsOnCuttingPlanes, AmOfPieces_mi_1, In_hg, radPtr, pSubdivOptions->SubdivideCoils, pSubdivOptions->PutNewStuffIntoGenCont);
		if(SubdOK) 
		{
			if(hgOld.rep != In_hg.rep) if(MessageChar==0) pVectOfHgChanged->push_back(In_hg);
		}
		return SubdOK;
	}

	radThg hgBuf = In_hg;
	if(!ConvertToPolyhedron(hgBuf, radPtr, pSubdivOptions->PutNewStuffIntoGenCont)) return 0;
	radThg hgBuf1 = hgBuf;
	if(!((radTg3d*)(hgBuf.rep))->SubdivideItselfByOneSetOfParPlanes(InPlanesNormal, InPointsOnCuttingPlanes, AmOfPieces_mi_1, hgBuf1, radPtr, pSubdivOptions, pVectOfHgChanged)) return 0;
	In_hg = hgBuf1;

	if(PointsOnCuttingPlanes != InPointsOnCuttingPlanes) delete[] PointsOnCuttingPlanes;
	return 1;
}

//-------------------------------------------------------------------------

int radTRecMag::SubdivideItselfByPlanesParToFace(
	short DirNo, TVector3d* PointsOnCuttingPlanes, int AmOfPieces_mi_1, 
	radThg& In_hg, radTApplication* radPtr, char SubdivideCoils, char PutNewStuffIntoGenCont)
{
	radTSend Send;
	double* CutCoords = new double[AmOfPieces_mi_1];
	if(CutCoords==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double AbsTol = RelAbsTol[1];

	double StEdge, FiEdge;
	if(DirNo==1) 
	{
		double HalfDim = 0.5*Dimensions.x;
		StEdge = CentrPoint.x - HalfDim;
		FiEdge = CentrPoint.x + HalfDim;
	}
	else if(DirNo==2)
	{
		double HalfDim = 0.5*Dimensions.y;
		StEdge = CentrPoint.y - HalfDim;
		FiEdge = CentrPoint.y + HalfDim;
	}
	else
	{
		double HalfDim = 0.5*Dimensions.z;
		StEdge = CentrPoint.z - HalfDim;
		FiEdge = CentrPoint.z + HalfDim;
	}

	int AmOfActualCuts = 0;
	for(int k=0; k<AmOfPieces_mi_1; k++)
	{
		TVector3d& CurrentPo = PointsOnCuttingPlanes[k];
		double CurrentCutCoord = (DirNo==1)? CurrentPo.x : ((DirNo==2)? CurrentPo.y : CurrentPo.z);
		if(((CurrentCutCoord - AbsTol)>StEdge) && ((CurrentCutCoord + AbsTol)<FiEdge))
		{
			CutCoords[AmOfActualCuts] = CurrentCutCoord; 
			AmOfActualCuts++;
		}
	}
	if(AmOfActualCuts==0) return 1;

	if(AmOfActualCuts > 1) if(CutCoords[1] < CutCoords[0])
	{
		for(int i=0; i<0.5*AmOfActualCuts; i++)
		{
			double Buf = CutCoords[i];
			double& Back = CutCoords[AmOfActualCuts - i - 1];
			CutCoords[i] = Back; Back = Buf;
		}
	}

	radThg Buf_hg = In_hg;
	char LocPutNewStuffIntoGenCont = 0;
	if(!DuplicateItself(Buf_hg, radPtr, LocPutNewStuffIntoGenCont)) return 0;
	radThg Buf_hgOld = Buf_hg;

	((radTg3d*)(Buf_hgOld.rep))->ConsiderOnlyWithTrans = 0;
	if(!((radTg3d*)(Buf_hgOld.rep))->ConvertToPolyhedron(Buf_hg, radPtr, LocPutNewStuffIntoGenCont)) return 0;
	radTPolyhedron* pPolyhdrCopy = (radTPolyhedron*)(Buf_hg.rep);
	pPolyhdrCopy->ConsiderOnlyWithTrans = ConsiderOnlyWithTrans;

	radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedPolyhedron(pPolyhdrCopy);
	radThg NewHandle(GroupInPlaceOfThisPtr);
	
	int LocElemCount = 0;
	int AmOfActualCuts_mi_1 = AmOfActualCuts-1;
	double CurCutCoord = StEdge, NextCutCoord = CutCoords[0];

	short NewFacesState[6], ParentFacesState[6];
	ListFacesInternalAfterCut(ParentFacesState);

	int StInd = 2*(DirNo-1);
	int FiInd = StInd+1;

	for(int j=0; j<AmOfActualCuts+1; j++)
	{
		short *pNew = NewFacesState, *pParent = ParentFacesState;
		for(int ii=0; ii<6; ii++) *(pNew++) = *(pParent++);

		double CurDim = NextCutCoord - CurCutCoord;
		double CurCen = CurCutCoord + 0.5*CurDim;

		TVector3d NewCenPoi = CentrPoint;
		TVector3d NewDims = Dimensions;
		if(DirNo==1) { NewCenPoi.x = CurCen; NewDims.x = CurDim;}
		else if(DirNo==2) { NewCenPoi.y = CurCen; NewDims.y = CurDim;}
		else { NewCenPoi.z = CurCen; NewDims.z = CurDim;}

		radTRecMag* NewRecMagPtr = new radTRecMag(NewCenPoi, NewDims, Magn, J, MaterHandle, J_IsNotZero);
		if(NewRecMagPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(NewRecMagPtr);
		GroupInPlaceOfThisPtr->AddElement(++LocElemCount, hg);

		if(j==0) NewFacesState[FiInd] = 1;
		else if(j==AmOfActualCuts) NewFacesState[StInd] = 1;
		else NewFacesState[StInd] = NewFacesState[FiInd] = 1;

		NewRecMagPtr->SetFacesInternalAfterCut(NewFacesState);
		NewRecMagPtr->SetMessageChar(1);

		CurCutCoord = NextCutCoord;
		NextCutCoord = (j>=AmOfActualCuts_mi_1)? FiEdge : CutCoords[j+1];
	}

	In_hg = NewHandle;

	if(CutCoords != 0) delete[] CutCoords;
	return 1;
}

//-------------------------------------------------------------------------

int radTRecMag::CutItself(TVector3d* InCuttingPlane, radThg& In_hg, radTPair_int_hg& LowerNewPair_int_hg, radTPair_int_hg& UpperNewPair_int_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char AddNewElemsToGenCont = pSubdivOptions->PutNewStuffIntoGenCont;
	char SeparatePiecesAtCutting = pSubdivOptions->SeparatePiecesAtCutting;
	char MapInternalFacesAfterCut = pSubdivOptions->MapInternalFacesAfterCut;

	double RelAbsTol[2];
	DefineRelAndAbsTol(RelAbsTol);
	double RelTol = RelAbsTol[0];
	double AbsTol = RelAbsTol[1];
	radTSend Send;

	if(J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero)
	{
		if(!SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error108"); return 0;}
	}

	TVector3d CuttingPlane[] = { *InCuttingPlane, InCuttingPlane[1] };
	TVector3d& PointOnCutPlane = *CuttingPlane;
	TVector3d& CutPlaneNormal = CuttingPlane[1];
	double SqLen = CutPlaneNormal.x*CutPlaneNormal.x + CutPlaneNormal.y*CutPlaneNormal.y + CutPlaneNormal.z*CutPlaneNormal.z;
	if(fabs(SqLen - 1.) > RelTol) CutPlaneNormal = (1./sqrt(SqLen))*CutPlaneNormal;

	radTrans ResTransf;
	short SomethingFound = 0;
	if(!SubdInLocFrame)
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	else if(ConsiderOnlyWithTrans)
	{
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
	}
	if(SomethingFound) 
	{
		CutPlaneNormal = ResTransf.TrBiPoint_inv(CutPlaneNormal);
		PointOnCutPlane = ResTransf.TrPoint_inv(PointOnCutPlane);
	}

	short NormDirForward;
	short DirNo = -1; //1- "x", 2- "y", 3- "z"
	double AbsNy = Abs(CutPlaneNormal.y), AbsNz = Abs(CutPlaneNormal.z);
	if((AbsNy < RelTol) && (AbsNz < RelTol))
	{
		DirNo = 1; NormDirForward = (Abs(CutPlaneNormal.x - 1.) < RelTol);
	}
	else
	{
		double AbsNx = Abs(CutPlaneNormal.x);
		if((AbsNx < RelTol) && (AbsNz < RelTol))
		{
			DirNo = 2; NormDirForward = (Abs(CutPlaneNormal.y - 1.) < RelTol);
		}
		else
		{
			if((AbsNx < RelTol) && (AbsNy < RelTol))
			{
				DirNo = 3; NormDirForward = (Abs(CutPlaneNormal.z - 1.) < RelTol);
			}
		}
	}
	if(DirNo < 0)
	{
		char VertexPtsPositionsChar;
		CheckVertexPtsPositionsWithRespectToPlane(CuttingPlane, VertexPtsPositionsChar);

		if(VertexPtsPositionsChar == 'I')
		{
			radThg hgBuf = In_hg;
			if(!ConvertToPolyhedron(hgBuf, radPtr, AddNewElemsToGenCont)) return 0;
			radThg hgBuf1 = hgBuf;
			if(!((radTg3d*)(hgBuf.rep))->CutItself(InCuttingPlane, hgBuf1, LowerNewPair_int_hg, UpperNewPair_int_hg, radPtr, pSubdivOptions)) return 0;
			In_hg = hgBuf1;
		}
		else if(VertexPtsPositionsChar == 'L')
		{
			if(SeparatePiecesAtCutting)
			{
				LowerNewPair_int_hg.m = radPtr->RetrieveElemKey(this);
				LowerNewPair_int_hg.Handler_g = In_hg;
			}
		}
		else if(VertexPtsPositionsChar == 'H')
		{
			if(SeparatePiecesAtCutting)
			{
				UpperNewPair_int_hg.m = radPtr->RetrieveElemKey(this);
				UpperNewPair_int_hg.Handler_g = In_hg;
			}
		}
		return 1;
	}

	double StEdge, FiEdge;
	if(DirNo==1) 
	{
		double HalfDim = 0.5*Dimensions.x;
		StEdge = CentrPoint.x - HalfDim;
		FiEdge = CentrPoint.x + HalfDim;
	}
	else if(DirNo==2)
	{
		double HalfDim = 0.5*Dimensions.y;
		StEdge = CentrPoint.y - HalfDim;
		FiEdge = CentrPoint.y + HalfDim;
	}
	else
	{
		double HalfDim = 0.5*Dimensions.z;
		StEdge = CentrPoint.z - HalfDim;
		FiEdge = CentrPoint.z + HalfDim;
	}

	TVector3d& CurrentPo = CuttingPlane[0];
	double CurrentCutCoord = (DirNo==1)? CurrentPo.x : ((DirNo==2)? CurrentPo.y : CurrentPo.z);

	char ThisItselfIsLowerElem = ((CurrentCutCoord+AbsTol) > FiEdge);
	char ThisItselfIsUpperElem = ((CurrentCutCoord-AbsTol) < StEdge);

	int StInd = 2*(DirNo-1);
	int FiInd = StInd+1;

	if(ThisItselfIsLowerElem)
	{
		if(SeparatePiecesAtCutting)
		{
			if(!MapInternalFacesAfterCut)
			{
				char CoinsidesWithFace = (CurrentCutCoord-AbsTol) < FiEdge;
				if(CoinsidesWithFace)
				{
					MapFaceAsExternal(FiInd + 1);
				}
			}
			LowerNewPair_int_hg.m = radPtr->RetrieveElemKey(this);
			LowerNewPair_int_hg.Handler_g = In_hg;
		}
		return 1;
	}
	if(ThisItselfIsUpperElem)
	{
		if(SeparatePiecesAtCutting)
		{
			if(!MapInternalFacesAfterCut)
			{
				char CoinsidesWithFace = (CurrentCutCoord+AbsTol) > StEdge;
				if(CoinsidesWithFace)
				{
					MapFaceAsExternal(StInd + 1);
				}
			}
			UpperNewPair_int_hg.m = radPtr->RetrieveElemKey(this);
			UpperNewPair_int_hg.Handler_g = In_hg;
		}
		return 1;
	}

	radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedRecMag(this);
	radThg NewHandle(GroupInPlaceOfThisPtr);

	short NewFacesState[6], ParentFacesState[6];
	ListFacesInternalAfterCut(ParentFacesState);
	short *pNew = NewFacesState, *pParent = ParentFacesState;
	for(int ii=0; ii<6; ii++) *(pNew++) = *(pParent++);

	double CurDim = CurrentCutCoord - StEdge;
	double CurCen = StEdge + 0.5*CurDim;

	TVector3d NewCenPoi = CentrPoint;
	TVector3d NewDims = Dimensions;
	if(DirNo==1) { NewCenPoi.x = CurCen; NewDims.x = CurDim;}
	else if(DirNo==2) { NewCenPoi.y = CurCen; NewDims.y = CurDim;}
	else { NewCenPoi.z = CurCen; NewDims.z = CurDim;}

	radTRecMag* NewRecMagPtr = new radTRecMag(NewCenPoi, NewDims, Magn, J, MaterHandle, J_IsNotZero);
	if(NewRecMagPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

	if(MapInternalFacesAfterCut) NewFacesState[FiInd] = 1;
	else NewFacesState[FiInd] = 0;

	NewRecMagPtr->SetFacesInternalAfterCut(NewFacesState);
	NewRecMagPtr->J_IsNotZero = J_IsNotZero; // To remove?
	if(ConsiderOnlyWithTrans)
	{
		NewRecMagPtr->ConsiderOnlyWithTrans = 1;
		NewRecMagPtr->g3dListOfTransform = g3dListOfTransform;
	}
	radThg hg1(NewRecMagPtr);

	NewFacesState[FiInd] = ParentFacesState[FiInd];

	CurDim = FiEdge - CurrentCutCoord;
	CurCen = CurrentCutCoord + 0.5*CurDim;

	if(DirNo==1) { NewCenPoi.x = CurCen; NewDims.x = CurDim;}
	else if(DirNo==2) { NewCenPoi.y = CurCen; NewDims.y = CurDim;}
	else { NewCenPoi.z = CurCen; NewDims.z = CurDim;}

	NewRecMagPtr = new radTRecMag(NewCenPoi, NewDims, Magn, J, MaterHandle, J_IsNotZero);
	if(NewRecMagPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

	if(MapInternalFacesAfterCut) NewFacesState[StInd] = 1; 
	else NewFacesState[StInd] = 0; 

	NewFacesState[FiInd] = ParentFacesState[FiInd];
	NewRecMagPtr->SetFacesInternalAfterCut(NewFacesState);
	NewRecMagPtr->J_IsNotZero = J_IsNotZero; // To remove?
	if(ConsiderOnlyWithTrans)
	{
		NewRecMagPtr->ConsiderOnlyWithTrans = 1;
		NewRecMagPtr->g3dListOfTransform = g3dListOfTransform;
	}
	radThg hg2(NewRecMagPtr);

	int LowerElemKey, UpperElemKey;
	radThg hgLowerNew, hgUpperNew;
	if(NormDirForward)
	{
		LowerElemKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hg1) : 1;
		hgLowerNew = hg1;
		UpperElemKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hg2) : 2;
		hgUpperNew = hg2;
	}
	else
	{
		LowerElemKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hg2) : 1;
		hgLowerNew = hg2;
		UpperElemKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hg1) : 2;
		hgUpperNew = hg1;
	}
		
	GroupInPlaceOfThisPtr->AddElement(LowerElemKey, hgLowerNew);
	GroupInPlaceOfThisPtr->AddElement(UpperElemKey, hgUpperNew);

	In_hg = NewHandle;

	if(SeparatePiecesAtCutting)
	{
		if(!g3dListOfTransform.empty())
		{
			radThg hgNewDpl = hgLowerNew;
			if(!(hgLowerNew.rep)->DuplicateItself(hgNewDpl, radPtr, AddNewElemsToGenCont)) return 0;
			((radTg3d*)(hgNewDpl.rep))->g3dListOfTransform = g3dListOfTransform;
			LowerElemKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hgNewDpl) : 1;
			hgLowerNew = hgNewDpl;

			hgNewDpl = hgUpperNew;
			if(!(hgUpperNew.rep)->DuplicateItself(hgNewDpl, radPtr, AddNewElemsToGenCont)) return 0;
			((radTg3d*)(hgNewDpl.rep))->g3dListOfTransform = g3dListOfTransform;
			UpperElemKey = AddNewElemsToGenCont? radPtr->AddElementToContainer(hgNewDpl) : 2;
			hgUpperNew = hgNewDpl;
		}
		LowerNewPair_int_hg.m = LowerElemKey;
		LowerNewPair_int_hg.Handler_g = hgLowerNew;

		UpperNewPair_int_hg.m = UpperElemKey;
		UpperNewPair_int_hg.Handler_g = hgUpperNew;
	
		if(AddNewElemsToGenCont)
		{
			int NewGroupElemKey = radPtr->RetrieveElemKey(GroupInPlaceOfThisPtr);
			radPtr->CopyDrawAttr(NewGroupElemKey, LowerElemKey);
			radPtr->CopyDrawAttr(NewGroupElemKey, UpperElemKey);
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTRecMag::FindLowestAndUppestVertices(TVector3d& PlanesNormal, radTSubdivOptions* pSubdivOptions, 
	TVector3d& LowestVertexPoint, TVector3d& UppestVertexPoint, radTrans& Trans, char& TransWasSet, char& Ignore)
{
	if(!pSubdivOptions->SubdivideCoils) if(J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero) { Ignore = 1; return 1;}
	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;

	Ignore = 0;
	TransWasSet = 0;
	TVector3d& ActualPlanesNormal = PlanesNormal;

	radTrans ResTransf;
	short SomethingFound = 0;
	if(!SubdInLocFrame)
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	else if(ConsiderOnlyWithTrans)
	{
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
	}
	if(SomethingFound) 
	{
		ActualPlanesNormal = ResTransf.TrBiPoint_inv(ActualPlanesNormal);
		Trans = ResTransf;
		TransWasSet = 1;
	}

	TVector3d LowestPo, UppestPo;

	TVector3d Vertices[8];
	TVector3d HalfDim = 0.5*Dimensions;
	TVector3d DimX(Dimensions.x, 0., 0.), DimY(0., Dimensions.y, 0.), DimZ(0., 0., Dimensions.z);
	*Vertices = CentrPoint - HalfDim;
	Vertices[1] = *Vertices + DimX;
	Vertices[2] = Vertices[1] + DimY;
	Vertices[3] = *Vertices + DimY;
	Vertices[4] = *Vertices + DimZ;
	Vertices[5] = Vertices[1] + DimZ;
	Vertices[6] = Vertices[2] + DimZ;
	Vertices[7] = Vertices[3] + DimZ;

	char Starting = 1;
	for(int i=0; i<8; i++)
	{
		TVector3d& TestPo = Vertices[i];
		if(Starting)
		{
			LowestPo = TestPo; UppestPo = TestPo; Starting = 0;
		}
		else
		{
			TVector3d TestLoV = TestPo - LowestPo;
			TVector3d TestUpV = TestPo - UppestPo;
			if(TestLoV*ActualPlanesNormal < 0.) LowestPo = TestPo;
			if(TestUpV*ActualPlanesNormal > 0.) UppestPo = TestPo;
		}
	}
	LowestVertexPoint = LowestPo;
	UppestVertexPoint = UppestPo;
	return 1;
}

//-------------------------------------------------------------------------

void radTRecMag::CheckVertexPtsPositionsWithRespectToPlane(TVector3d* Plane, char& VertexPtsPositionsChar)
{
	TVector3d Vertices[8];
	TVector3d HalfDim = 0.5*Dimensions;
	*Vertices = CentrPoint - HalfDim;
	Vertices[1] = *Vertices; Vertices[1].x += Dimensions.x;
	Vertices[2] = Vertices[1]; Vertices[2].y += Dimensions.y;
	Vertices[3] = *Vertices; Vertices[3].y += Dimensions.y;
	Vertices[4] = *Vertices; Vertices[4].z += Dimensions.z;
	Vertices[5] = Vertices[1]; Vertices[5].z += Dimensions.z;
	Vertices[6] = Vertices[2]; Vertices[6].z += Dimensions.z;
	Vertices[7] = Vertices[3]; Vertices[7].z += Dimensions.z;

	TVector3d& PtOnPlane = *Plane;
	TVector3d& PlaneNormal = Plane[1];
	TVector3d V;

	char AllPtsAreHigher = 1, AllPtsAreLower = 1;
	for(int i=0; i<8; i++)
	{
		V = Vertices[i] - PtOnPlane;
		double ScalProd = V*PlaneNormal;
		if(ScalProd >= 0.) AllPtsAreLower = 0;
		else AllPtsAreHigher = 0;
	}

	if(AllPtsAreHigher) VertexPtsPositionsChar = 'H';
	else if(AllPtsAreLower) VertexPtsPositionsChar = 'L';
	else VertexPtsPositionsChar = 'I';
}

//-------------------------------------------------------------------------

void radTRecMag::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
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

	if(pFieldKey->M_)
	{
		if(J_IsNotZero) return;
		Pair.V2 = (pTrans == 0)? Magn : pTrans->TrVectField(Magn);
	}
	else if(pFieldKey->J_)
	{
		if(!J_IsNotZero) return;
		Pair.V2 = (pTrans == 0)? J : pTrans->TrVectField(J);
	}
	else
	{
		radTCompCriterium CompCriterium;
		TVector3d ZeroVect(0.,0.,0.);
		radTField Field(*pFieldKey, CompCriterium, CP, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		g3dSrcPtr->B_genComp(&Field);
		Pair.V2 = (pFieldKey->B_)? Field.B : ((pFieldKey->H_)? Field.H : ((pFieldKey->A_)? Field.A : ZeroVect));
	}
	pVectPairOfVect3d->push_back(Pair);
}

//-------------------------------------------------------------------------

void radTRecMag::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique)
{
	TVector3d HalfDim = 0.5*Dimensions;
	double xMin = CentrPoint.x - HalfDim.x, xMax = CentrPoint.x + HalfDim.x;
	double yMin = CentrPoint.y - HalfDim.y, yMax = CentrPoint.y + HalfDim.y;
	double zMin = CentrPoint.z - HalfDim.z, zMax = CentrPoint.z + HalfDim.z;

	TVector3d P1(xMin, yMin, zMin); OutVect.push_back(P1);
	TVector3d P2(xMax, yMin, zMin); OutVect.push_back(P2);
	TVector3d P3(xMin, yMax, zMin); OutVect.push_back(P3);
	TVector3d P4(xMax, yMax, zMin); OutVect.push_back(P4);

	TVector3d P5(xMin, yMin, zMax); OutVect.push_back(P5);
	TVector3d P6(xMax, yMin, zMax); OutVect.push_back(P6);
	TVector3d P7(xMin, yMax, zMax); OutVect.push_back(P7);
	TVector3d P8(xMax, yMax, zMax); OutVect.push_back(P8);
}

//-------------------------------------------------------------------------
