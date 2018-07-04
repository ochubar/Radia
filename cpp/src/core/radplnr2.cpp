/*-------------------------------------------------------------------------
*
* File name:      radplnr2.cpp
*
* Project:        RADIA
*
* Description:    Auxiliary 2D objects
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radplnr.h"
#include "radappl.h"

//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;
extern radTYield radYield;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTPolygon::B_comp(radTField* FieldPtr)
{
// Orientation: The Polygon normal parallel to vertical ort !!!
	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;

	const double Max_k = 1.E+09; //1.E+08; // To skip segments in general field computation loop.

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	TVector2d& FirstPo2d = *BaseIter;

	double AbsRandX = radCR.AbsRandMagnitude(FirstPo2d.x - CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(FirstPo2d.y - CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(CoordZ);
	if(AbsRandZ < AbsRandX) AbsRandZ = AbsRandX;
	if(AbsRandZ < AbsRandY) AbsRandZ = AbsRandY;
	double RelRandMagn = radCR.AbsRandMagnitude(1.);

	if(radYield.Check()==0) return; // To allow multitasking on Mac: consider better places for this

	TVector3d& ObsPo = FieldPtr->P;

	double z = CoordZ - ObsPo.z;  // If this is zero, we get error !
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(z==0.) 
	{
		z = AbsRandZ;
        ObsPo.z -= AbsRandZ; //OC040504 test
	}
	double ze2 = z*z;
	double absz = Abs(z);

	TVector2d First2d(FirstPo2d.x - ObsPo.x, FirstPo2d.y - ObsPo.y);
	TVector2d Vect2dToAdd(-ObsPo.x, -ObsPo.y);

	double x1 = First2d.x;  // If this is zero, we get error ?
	double y1 = First2d.y;  // If this is zero, we get error ?
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(x1==0.) 
	{
		x1 = AbsRandX;

        First2d.x += AbsRandX; //OC040504 test
		Vect2dToAdd.x += AbsRandX; //OC040504 test
		ObsPo.x -= AbsRandX; //OC040504 test
	}
	if(y1==0.) 
	{
		y1 = AbsRandY;

        First2d.y += AbsRandY; //OC040504 test
		Vect2dToAdd.y += AbsRandY; //OC040504 test
		ObsPo.y -= AbsRandY; //OC040504 test
	}
	double x1e2 = x1*x1, y1e2 = y1*y1;
	double x2, y2, x2e2, y2e2;

	short A_CompNeeded = FieldPtr->FieldKey.A_;
	short B_orH_CompNeeded = FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_ || FieldPtr->FieldKey.PreRelax_;

	int AmOfEdgePoints_mi_1 = AmOfEdgePoints - 1;

	double Sx=0., Sy=0., Sz=0., AS=0.;
	double ArgSumAtans1=0., PiMultSumAtans1=0.;
	double ArgSumLogs2=1.;

	for(int i=0; i<AmOfEdgePoints; i++)
	{
		++BaseIter;
		if(i!=AmOfEdgePoints_mi_1)
		{
			TVector2d& BufPo2d = *BaseIter;
			x2 = BufPo2d.x + Vect2dToAdd.x;  // If this is zero, we get error !
			y2 = BufPo2d.y + Vect2dToAdd.y;  // If this is zero, we get error !
		}
		else
		{
			x2 = First2d.x;  // If this is zero, we get error !
			y2 = First2d.y;  // If this is zero, we get error !
		}

		// Artificial shift of an observation point a bit right of the block's border
		// if the point is exactly on the boarder (to avoid "divide by zero" error):
		// Removing this may be dangerous for the Checking-If-Inside
		if(x2==0.) 
		{
			x2 = AbsRandX;
            
			First2d.x += AbsRandX; //OC040504 test
            Vect2dToAdd.x += AbsRandX; //OC040504 test
            ObsPo.x -= AbsRandX; //OC040504 test
		}
		if(y2==0.) 
		{
			y2 = AbsRandY;

            First2d.y += AbsRandY; //OC040504 test
            Vect2dToAdd.y += AbsRandY; //OC040504 test
            ObsPo.y -= AbsRandY; //OC040504 test
		}

		x2e2 = x2*x2; y2e2 = y2*y2;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = y2my1/x2mx1; // x2mx1 != 0 here
			double b = y1 - k*x1;  // If b is zero, we get error !
			if(b==0.) 
			{
				b = AbsRandY; //OC040504 test

				y1 += AbsRandY;
				y1e2 = y1*y1;
				y2my1 -= AbsRandY;
				abs_y2my1 = Abs(y2my1);
				k = y2my1/x2mx1;

				First2d.y += AbsRandY; //OC040504 test
				Vect2dToAdd.y += AbsRandY; //OC040504 test
				ObsPo.y -= AbsRandY; //OC040504 test
			}

			double bk = b*k, ke2 = k*k, be2 = b*b, twob = 2.*b;
			double ke2p1 = ke2+1.;
			if(ke2p1==0.) ke2p1 = RelRandMagn;
			double sqrtke2p1 = sqrt(ke2p1);

			double kx1 = k*x1, kx2 = k*x2;
			double bpkx1 = y1, bpkx2 = y2;
			double bpkx1e2 = bpkx1*bpkx1, bpkx2e2 = bpkx2*bpkx2;
			double kx1mb = -b+kx1, kx2mb = -b+kx2; 
			//double R1 = radCR.DoublePlus(sqrt(x1e2 + bpkx1e2 + ze2)), R2 = radCR.DoublePlus(sqrt(x2e2 + bpkx2e2 + ze2));
			double R1 = sqrt(x1e2 + bpkx1e2 + ze2), R2 = sqrt(x2e2 + bpkx2e2 + ze2); //OC040504

			double x1e2pze2 = x1e2+ze2, x2e2pze2 = x2e2+ze2;
			double bkpx1pke2x1 = bk+ke2p1*x1, bkpx2pke2x2 = bk+ke2p1*x2;  // If this is zero, we get error !
			double kze2 = k*ze2;
			double ke2ze2 = k*kze2;
			double ke2ze2mbe2 = ke2ze2-be2, ke2ze2pbe2 = ke2ze2+be2;
			double bx1 = b*x1, bx2 = b*x2;
			double R1pbpkx1 = bpkx1+R1, R2pbpkx2 = bpkx2+R2;

			const double MaxRelTolToSwitch = 1.E-07; //1.E-05;
            //double AbsRandR1 = 10.*radCR.AbsRandMagnitude(R1);
			//double AbsRandR2 = 10.*radCR.AbsRandMagnitude(R2);

			double AbsRandR2 = 100.*radCR.AbsRandMagnitude(R2);
			double AbsRandR1 = 100.*radCR.AbsRandMagnitude(R1);
			double MaxAbsRandR1 = MaxRelTolToSwitch*R1;
			double MaxAbsRandR2 = MaxRelTolToSwitch*R2;
			if(AbsRandR1 > MaxAbsRandR1) AbsRandR1 = MaxAbsRandR1;
			if(AbsRandR2 > MaxAbsRandR2) AbsRandR2 = MaxAbsRandR2;
			//double AbsRandR1 = MaxRelTolToSwitch*R1; //OC040504
			//double AbsRandR2 = MaxRelTolToSwitch*R2; //OC040504

			//if(R1pbpkx1 < AbsRandR1) R1pbpkx1 = 0.5*(x1e2 + ze2)/Abs(bpkx1);
			//if(R2pbpkx2 < AbsRandR2) R2pbpkx2 = 0.5*(x2e2 + ze2)/Abs(bpkx2);
			//if((R1pbpkx1 < AbsRandR1) && (bpkx1 != 0.)) R1pbpkx1 = 0.5*(x1e2 + ze2)/Abs(bpkx1); //OC040504
			//if((R2pbpkx2 < AbsRandR2) && (bpkx2 != 0.)) R2pbpkx2 = 0.5*(x2e2 + ze2)/Abs(bpkx2); //OC040504
			if((R1pbpkx1 < AbsRandR1) && (R1 > 100.*AbsRandR1) && ((x1e2 + ze2) < bpkx1e2*MaxRelTolToSwitch)) R1pbpkx1 = 0.5*(x1e2 + ze2)/Abs(bpkx1); //OC170504
			if((R2pbpkx2 < AbsRandR2) && (R2 > 100.*AbsRandR2) && ((x2e2 + ze2) < bpkx2e2*MaxRelTolToSwitch)) R2pbpkx2 = 0.5*(x2e2 + ze2)/Abs(bpkx2); //OC170504


			double FlpRep1ForSumAtans1 = 0.;

			double four_be2ke2 = 4.*be2*ke2; 
			double four_be2be2ke2 = be2*four_be2ke2;
			double be2mke2ze2 = be2-ke2ze2, be2pke2ze2 = be2+ke2ze2;
			double be2mke2ze2e2 = be2mke2ze2*be2mke2ze2, be2pke2ze2e2 = be2pke2ze2*be2pke2ze2;
			double DFlipRepSumAtans1 = (be2+ke2p1*ze2)*(four_be2ke2*(be2+ke2ze2)-be2mke2ze2e2);
			double BufDen = four_be2be2ke2-ke2p1*be2mke2ze2e2;

			if((DFlipRepSumAtans1 >= 0.) && (BufDen != 0.))
			{
				double Buf1Num = bk*be2pke2ze2e2;
				double Buf2Num = be2mke2ze2*sqrt(DFlipRepSumAtans1);
				double xFlp1 = (Buf1Num - Buf2Num)/BufDen;
				double xFlp2 = (Buf1Num + Buf2Num)/BufDen;

				double xFlp = xFlp1;
				if((x1<x2)? ((xFlp>x1) && (xFlp<x2)) : ((xFlp<x1) && (xFlp>x2)))
				{
					double xFlpe2 = xFlp*xFlp, kxFlp = k*xFlp;
					double kxFlppb = kxFlp+b, kxFlpmb = kxFlp-b;
					double kxFlppbe2 = kxFlppb*kxFlppb;
					double SqRoot = sqrt(xFlpe2+kxFlppbe2+ze2);

					if(Sign((xFlpe2+ze2)*(-be2mke2ze2) + (-be2+ke2*xFlpe2)*be2pke2ze2) == Sign(-kxFlpmb)) // RootIsReal?
					{
						double DenomDerivSign = Sign(-2.*xFlp*be2mke2ze2 + kxFlpmb*be2pke2ze2*(k+(bk+ke2p1*xFlp)/SqRoot) + k*be2pke2ze2*(kxFlppb + SqRoot));
						//double DenomDerivSign = Sign(-2.*xFlp*be2mke2ze2*SqRoot + kxFlpmb*be2pke2ze2*(k*SqRoot + (bk+ke2p1*xFlp)) + k*be2pke2ze2*(kxFlppb + SqRoot)*SqRoot); //OC040504
						
						double NumSign = Sign((2.*bk*ze2*(xFlpe2+ze2) + (b*xFlp+kze2)*be2pke2ze2*(kxFlppb + SqRoot))/z);
						//double NumSign = Sign((2.*bk*ze2*(xFlpe2+ze2) + (b*xFlp+kze2)*be2pke2ze2*(kxFlppb + SqRoot))*z); //OC040504
						double Buf = -DenomDerivSign*NumSign*Sign(x2mx1);

						FlpRep1ForSumAtans1 += Buf;
					}
				}
				xFlp = xFlp2;
				if((x1<x2)? ((xFlp>x1) && (xFlp<x2)) : ((xFlp<x1) && (xFlp>x2)))
				{
					double xFlpe2 = xFlp*xFlp, kxFlp = k*xFlp;
					double kxFlppb = kxFlp+b, kxFlpmb = kxFlp-b;
					double kxFlppbe2 = kxFlppb*kxFlppb;
					double SqRoot = sqrt(xFlpe2+kxFlppbe2+ze2);
					if(Sign((xFlpe2+ze2)*(-be2mke2ze2) + (-be2+ke2*xFlpe2)*be2pke2ze2) == Sign(-kxFlpmb))
					{
						double DenomDerivSign = Sign(-2.*xFlp*be2mke2ze2 + kxFlpmb*be2pke2ze2*(k+(bk+ke2p1*xFlp)/SqRoot) + k*be2pke2ze2*(kxFlppb + SqRoot));
						//double DenomDerivSign = Sign(-2.*xFlp*be2mke2ze2*SqRoot + kxFlpmb*be2pke2ze2*(k*SqRoot + (bk+ke2p1*xFlp)) + k*be2pke2ze2*(kxFlppb + SqRoot)*SqRoot); //OC040504

						double NumSign = Sign((2.*bk*ze2*(xFlpe2+ze2) + (b*xFlp+kze2)*be2pke2ze2*(kxFlppb + SqRoot))/z);
						//double NumSign = Sign((2.*bk*ze2*(xFlpe2+ze2) + (b*xFlp+kze2)*be2pke2ze2*(kxFlppb + SqRoot))*z);
						double Buf = -DenomDerivSign*NumSign*Sign(x2mx1);

						FlpRep1ForSumAtans1 += Buf;
					}
				}
			}

			double Arg1ForSumAtans1 = -(ke2ze2pbe2*(bx1 + kze2)*R1pbpkx1 + kze2*twob*x1e2pze2);
			double Arg2ForSumAtans1 = (ke2ze2pbe2*kx1mb*R1pbpkx1 + ke2ze2mbe2*x1e2pze2)*z;
			double Arg3ForSumAtans1 = ke2ze2pbe2*(bx2 + kze2)*R2pbpkx2 + kze2*twob*x2e2pze2;
			double Arg4ForSumAtans1 = (ke2ze2pbe2*kx2mb*R2pbpkx2 + ke2ze2mbe2*x2e2pze2)*z;

			double PiMult1=0., PiMult2=0.;
			double CurArgSumAtans1 = TransAtans(Arg1ForSumAtans1/Arg2ForSumAtans1, Arg3ForSumAtans1/Arg4ForSumAtans1, PiMult1);

			//double DivVal = Arg2ForSumAtans1*Arg4ForSumAtans1 - Arg1ForSumAtans1*Arg3ForSumAtans1; //OC040504
            //if(DivVal == 0.) DivVal = 1.e-50; //OC040504
			//double CurArgSumAtans1 = (Arg1ForSumAtans1*Arg4ForSumAtans1 + Arg3ForSumAtans1*Arg2ForSumAtans1)/DivVal; //OC040504
			//if(Arg2ForSumAtans1*Arg4ForSumAtans1 != 0.) //OC040504
			//{
			//	double x = Arg1ForSumAtans1/Arg2ForSumAtans1, y = Arg3ForSumAtans1/Arg4ForSumAtans1;
            //  double Buf = 1.-x*y;
            //             PiMult1 = (((Buf > 0)? 0.:1.) * ((x < 0)? -1.:1.));
			//}
			//else //OC040504
			//{
			//	if(Arg1ForSumAtans1*Arg3ForSumAtans1 < 0)
			//	{
			//		PiMult1 = 0;
			//	}
			//	else
			//	{
			//      if(Arg2ForSumAtans1 < 0.)
			//		{
			//			if(Arg1ForSumAtans1 > 0.)
			//			{
			//                         PiMult1 = -1;
			//			}
			//			else PiMult1 = 1;
			//		}
			//		else
			//		{
			//			if(Arg1ForSumAtans1 > 0.)
			//			{
			//                         PiMult1 = 1;
			//			}
			//			else PiMult1 = -1;
			//		}
			//	}
			//}

			ArgSumAtans1 = TransAtans(ArgSumAtans1, CurArgSumAtans1, PiMult2);
			PiMultSumAtans1 += PiMult1+PiMult2 + FlpRep1ForSumAtans1;

			double bkpx1pke2x1dsqrtke2p1pR1 = bkpx1pke2x1/sqrtke2p1 + R1; // sqrtke2p1 > 0 always
			double bkpx2pke2x2dsqrtke2p1pR2 = bkpx2pke2x2/sqrtke2p1 + R2; // sqrtke2p1 > 0 always
			//if(bkpx1pke2x1dsqrtke2p1pR1 < AbsRandR1) bkpx1pke2x1dsqrtke2p1pR1 = 0.5*(be2 + ze2)/(Abs(x1)*sqrtke2p1);
			//if(bkpx2pke2x2dsqrtke2p1pR2 < AbsRandR2) bkpx2pke2x2dsqrtke2p1pR2 = 0.5*(be2 + ze2)/(Abs(x2)*sqrtke2p1);
			//if((bkpx1pke2x1dsqrtke2p1pR1 < AbsRandR1) && (x1 != 0.)) bkpx1pke2x1dsqrtke2p1pR1 = 0.5*(be2 + ze2)/(Abs(x1)*sqrtke2p1);
			//if((bkpx2pke2x2dsqrtke2p1pR2 < AbsRandR2) && (x2 != 0.)) bkpx2pke2x2dsqrtke2p1pR2 = 0.5*(be2 + ze2)/(Abs(x2)*sqrtke2p1);

			double be2pze2 = be2 + ze2;
			if((bkpx1pke2x1dsqrtke2p1pR1 < AbsRandR1) && (R1 > 100.*AbsRandR1) && ((be2pze2 + 2*bk*x1) < x1e2*ke2p1*MaxRelTolToSwitch)) bkpx1pke2x1dsqrtke2p1pR1 = 0.5*be2pze2/(Abs(x1)*sqrtke2p1);
			if((bkpx2pke2x2dsqrtke2p1pR2 < AbsRandR2) && (R2 > 100.*AbsRandR2) && ((be2pze2 + 2*bk*x2) < x2e2*ke2p1*MaxRelTolToSwitch)) bkpx2pke2x2dsqrtke2p1pR2 = 0.5*be2pze2/(Abs(x2)*sqrtke2p1);

			if(bkpx1pke2x1dsqrtke2p1pR1 == 0.) bkpx1pke2x1dsqrtke2p1pR1 = 1.e-50; //OC040504
			if(bkpx2pke2x2dsqrtke2p1pR2 == 0.) bkpx2pke2x2dsqrtke2p1pR2 = 1.e-50; //OC040504

			double SumLogs1 = log(bkpx2pke2x2dsqrtke2p1pR2/bkpx1pke2x1dsqrtke2p1pR1);
			double SumLogs1dsqrtke2p1 = SumLogs1/sqrtke2p1; // sqrtke2p1 > 0 always

			if(B_orH_CompNeeded)
			{
				if(R1pbpkx1 == 0.) R1pbpkx1 = 1.e-50; //OC040504
				ArgSumLogs2 *= (R2pbpkx2/R1pbpkx1);
				
				Sx += -k*SumLogs1dsqrtke2p1;
				Sy += SumLogs1dsqrtke2p1;
			}
			if(A_CompNeeded)
			{
				if(R2pbpkx2 <= 0.) R2pbpkx2 = 1.e-50; //OC040504
				if(R1pbpkx1 <= 0.) R1pbpkx1 = 1.e-50; //OC040504
				AS += -b*SumLogs1dsqrtke2p1 - (x2*log(R2pbpkx2) - x1*log(R1pbpkx1));
			}
		}
		x1 = x2; y1 = y2;
		x1e2 = x2e2; y1e2 = y2e2;
	}

	FieldPtr->PointIsInsideFrame = (z > 0.);

	Sz = atan(ArgSumAtans1) + PiMultSumAtans1*PI;
	if(B_orH_CompNeeded)
	{
		if(ArgSumLogs2 <= 0.) ArgSumLogs2 = 1.e-50; //OC040504
		Sx += log(ArgSumLogs2);

		if(FieldPtr->FieldKey.PreRelax_)
		{
			TVector3d St0(0., 0., -ConstForH*Sx); 
			TVector3d St1(0., 0., -ConstForH*Sy); 
			TVector3d St2(0., 0., -ConstForH*Sz);
			FieldPtr->B += St0;
			FieldPtr->H += St1;
			FieldPtr->A += St2;
			return;
		}
		TVector3d BufH(Sx, Sy, Sz);
		FieldPtr->H += (-ConstForH*Magn.z)*BufH; // Check if "-" is necessary
	}
	if(A_CompNeeded)
	{
		AS += -z*Sz;
		TVector3d BufA(-Magn.y, Magn.x, 0.);
		FieldPtr->A += (-AS*ConstForH)*BufA; // Check if "-" is necessary
	}
}

//-------------------------------------------------------------------------

void radTPolygon::B_intComp(radTField* FieldPtr)
{
// Orientation: The Polygon normal parallel to vertical ort !!!
// Field Integrals are correct only if line does not cross the prism body.
// This was noticed to have the "divide by zero" problem if (Vx*Vx+Vy*Vy==0)||(Vx*Vx+Vz*Vz==0)||(Vy*Vy+Vz*Vz==0)

	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;
	const double Max_k = 1.E+08; // To skip segments in general field int. computation loop.
	const double ZeroToler = 1.E-06; // This is to switch to Special Cases
	const double SmallestRelTolerV = 1.E-12; // Relative tolerance to repair trapping V.i to zero in general case

	TVector3d V = FieldPtr->NextP - FieldPtr->P;
	double InvAbsV = 1./sqrt(V.x*V.x + V.y*V.y + V.z*V.z);
	V = InvAbsV*V;

	short InitVxIsZero = (Abs(V.x) <= ZeroToler);
	short InitVyIsZero = (Abs(V.y) <= ZeroToler);
	short InitVzIsZero = (Abs(V.z) <= ZeroToler);

// This is the attempt to avoid "divide by zero" problem
	TSpecCaseID SpecCaseID = NoSpecCase;
	if(InitVxIsZero && InitVyIsZero) SpecCaseID = ZeroVxVy;
	else if(InitVxIsZero && InitVzIsZero) SpecCaseID = ZeroVxVz;
	else if(InitVyIsZero && InitVzIsZero) SpecCaseID = ZeroVyVz;
	if(SpecCaseID==ZeroVxVy || SpecCaseID==ZeroVxVz || SpecCaseID==ZeroVyVz) { B_intCompSpecCases(FieldPtr, SpecCaseID); return;}

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(CoordZ);

	double Vx = V.x; if(Vx==0.) Vx = SmallestRelTolerV;
	double Vy = V.y; if(Vy==0.) Vy = SmallestRelTolerV;
	double Vz = V.z; if(Vz==0.) Vz = SmallestRelTolerV;

	double Vxe2 = Vx*Vx, Vye2 = Vy*Vy, Vze2 = Vz*Vz;
	double Vxe2pVye2 = Vxe2+Vye2, Vxe2pVze2 = Vxe2+Vze2, Vye2pVze2 = Vye2+Vze2;
	double Vye2pVze2Vy = Vye2pVze2*Vy, Vx1pVye2 = Vx*(1.+Vye2);
	double VxVy = Vx*Vy, VyVz = Vy*Vz;

	double p2dVxe2pVze2 = 2./Vxe2pVze2;
	TVector3d& StPo = FieldPtr->P;

	double z = CoordZ - StPo.z;
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(z==0.) z = AbsRandZ;
	double Vxz = Vx*z, Vyz = Vy*z, VyVzz = VyVz*z, Vzz = Vz*z;

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	int AmOfEdgePoints_mi_1 = AmOfEdgePoints - 1;

	TVector2d& FirstPo2d = *BaseIter;
	TVector2d First2d(FirstPo2d.x - StPo.x, FirstPo2d.y - StPo.y);
	TVector2d Vect2dToAdd(-StPo.x, -StPo.y);

	double x1 = First2d.x;
	double y1 = First2d.y;
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(x1==0.) x1 = AbsRandY;
	if(y1==0.) y1 = AbsRandZ;
	double x2, y2;

	double Sx=0., Sy=0., Sz=0.;

	for(int i=0; i<AmOfEdgePoints; i++)
	{
		++BaseIter;
		if(i!=AmOfEdgePoints_mi_1)
		{
			TVector2d& BufPo2d = *BaseIter;
			x2 = BufPo2d.x + Vect2dToAdd.x;
			y2 = BufPo2d.y + Vect2dToAdd.y;
		}
		else
		{
			x2 = First2d.x;
			y2 = First2d.y;
		}
		// Artificial shift of an observation point a bit right of the block's border
		// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(x2==0.) x2 = AbsRandX;
		if(y2==0.) y2 = AbsRandY;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = y2my1/x2mx1, b = y1 - k*x1;
			double ke2 = k*k;
			double ke2p1 = ke2 + 1.;
			double kVx = k*Vx, kVy = k*Vy, kVz = k*Vz;
			double kVxmVy = kVx - Vy, kVypVx = kVy + Vx;
			double kVypVxVz = kVypVx*Vz;
			double kVxmVye2p1pke2Vze2 = kVxmVy*kVxmVy + ke2p1*Vze2;
			double bVx = b*Vx, bVz = b*Vz;
			double kVxe2pVze2mVxVy = k*Vxe2pVze2 - VxVy;
			double bkVxe2pVze2mVxVy = b*kVxe2pVze2mVxVy;
			double kVxmVyx1 = kVxmVy*x1, kVxmVyx2 = kVxmVy*x2;
			double kVxmVye2p1pke2Vze2x1 = kVxmVye2p1pke2Vze2*x1, kVxmVye2p1pke2Vze2x2 = kVxmVye2p1pke2Vze2*x2;
			double kVypVxVzz = kVypVxVz*z, kVxmVyz = kVxmVy*z;
			double bVxpkVxmVyx1 = bVx + kVxmVyx1, bVxpkVxmVyx2 = bVx + kVxmVyx2;
			double Vzx1 = Vz*x1, Vzx2 = Vz*x2;
			double Vzx1mVxz = Vzx1-Vxz, Vzx2mVxz = Vzx2-Vxz;
			double kVzx1 = kVz*x1, kVzx2 = kVz*x2;
			double bVzpkVzx1mVyz = bVz + kVzx1 - Vyz, bVzpkVzx2mVyz = bVz + kVzx2 - Vyz;
			double bVzpkVxmVyz = bVz + kVxmVyz; 
			double bVxpkVxmVyx1e2 = bVxpkVxmVyx1*bVxpkVxmVyx1, bVxpkVxmVyx2e2 = bVxpkVxmVyx2*bVxpkVxmVyx2;
			double bVxe2pVze2 = b*Vxe2pVze2;
			double kVxe2pVze2mVxVyx1 = kVxe2pVze2mVxVy*x1, kVxe2pVze2mVxVyx2 = kVxe2pVze2mVxVy*x2;

			double PiMult1=0.;
			double SumAtans1 = atan(TransAtans(-(bkVxe2pVze2mVxVy + kVxmVye2p1pke2Vze2x1 - kVypVxVzz)/bVzpkVxmVyz, (bkVxe2pVze2mVxVy + kVxmVye2p1pke2Vze2x2 - kVypVxVzz)/bVzpkVxmVyz, PiMult1));
			SumAtans1 += PiMult1*PI;

			double AtanX1 = atan((bVxe2pVze2 + kVxe2pVze2mVxVyx1 - VyVzz)/(Vxz - Vzx1));
			double AtanX2 = atan((bVxe2pVze2 + kVxe2pVze2mVxVyx2 - VyVzz)/(Vxz - Vzx2));

			double LogX1 = log(bVxpkVxmVyx1e2 + bVzpkVzx1mVyz*bVzpkVzx1mVyz + Vzx1mVxz*Vzx1mVxz);
			double LogX2 = log(bVxpkVxmVyx2e2 + bVzpkVzx2mVyz*bVzpkVzx2mVyz + Vzx2mVxz*Vzx2mVxz);

			double kVypVxVzVz = kVypVxVz*Vz;
			double kVxVy = k*VxVy;
			double BufLogMult1 = (kVxVy-Vye2pVze2)*bVxe2pVze2 + Vzz*(Vye2pVze2Vy-k*Vx1pVye2);

			double BufLogMult2 = (kVypVxVzz - b*kVxe2pVze2mVxVy)/kVxmVye2p1pke2Vze2;
			double BufLogMult3 = kVypVxVz*bVxe2pVze2 + (Vx*kVxmVy - Vy*kVypVxVzVz)*z;

			double p2dVxe2pVze2dkVxmVye2p1pke2Vze2 = p2dVxe2pVze2/kVxmVye2p1pke2Vze2;

			Sx += p2dVxe2pVze2dkVxmVye2p1pke2Vze2*((kVz*bVxe2pVze2 + (VxVy*kVxmVy - kVypVxVzVz)*z)*SumAtans1 + Vz*(kVxmVye2p1pke2Vze2x2*AtanX2 - kVxmVye2p1pke2Vze2x1*AtanX1) 
				+ 0.5*((VxVy*kVxmVye2p1pke2Vze2x2 + BufLogMult1)*LogX2 - (VxVy*kVxmVye2p1pke2Vze2x1 + BufLogMult1)*LogX1));

			Sy += -2.*(bVzpkVxmVyz/kVxmVye2p1pke2Vze2)*SumAtans1 + (BufLogMult2-x2)*LogX2 - (BufLogMult2-x1)*LogX1;

			Sz += p2dVxe2pVze2dkVxmVye2p1pke2Vze2*(-(b*kVxmVy*Vxe2pVze2 + (Vye2-Vxe2-2.*kVxVy)*Vzz)*SumAtans1 - Vx*(kVxmVye2p1pke2Vze2x2*AtanX2 - kVxmVye2p1pke2Vze2x1*AtanX1)
				+ 0.5*((VyVz*kVxmVye2p1pke2Vze2x2 + BufLogMult3)*LogX2 - (VyVz*kVxmVye2p1pke2Vze2x1 + BufLogMult3)*LogX1));
		}
		x1 = x2; y1 = y2;
	}
	TVector3d BufIh(Sx, Sy, Sz);
	double MultIh = -ConstForH*Magn.z; // Check if "-" is necessary
	BufIh = MultIh*BufIh;

	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIh;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIh;
}

//-------------------------------------------------------------------------

void radTPolygon::B_intCompSpecCases(radTField* FieldPtr, const TSpecCaseID& SpecCaseID)
{
	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;

	const double Max_k = 1.E+08;

	TVector3d& StPo = FieldPtr->P;

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(CoordZ);

	double z = CoordZ - StPo.z;
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(z==0.) z = AbsRandZ;
	double ze2 = z*z;

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	int AmOfEdgePoints_mi_1 = AmOfEdgePoints - 1;

	TVector2d& FirstPo2d = *BaseIter;
	TVector2d First2d(FirstPo2d.x - StPo.x, FirstPo2d.y - StPo.y);
	TVector2d Vect2dToAdd(-StPo.x, -StPo.y);

	double x1 = First2d.x;
	double y1 = First2d.y;
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(x1==0.) x1 = AbsRandX;
	if(y1==0.) y1 = AbsRandY;
	double x2, y2;
	double x1e2 = x1*x1, x2e2, y1e2 = y1*y1, y2e2;

	double Sx=0., Sy=0., Sz=0.;
	double PiMult1=0.;

	for(int i=0; i<AmOfEdgePoints; i++)
	{
		++BaseIter;
		if(i!=AmOfEdgePoints_mi_1)
		{
			TVector2d& BufPo2d = *BaseIter;
			x2 = BufPo2d.x + Vect2dToAdd.x;
			y2 = BufPo2d.y + Vect2dToAdd.y;
		}
		else
		{
			x2 = First2d.x;
			y2 = First2d.y;
		}
		// Artificial shift of an observation point a bit right of the block's border
		// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(x2==0.) x2 = AbsRandX;
		if(y2==0.) y2 = AbsRandY;
		x2e2 = x2*x2; y2e2 = y2*y2;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			if(SpecCaseID==ZeroVxVy)
			{
				double k = y2my1/x2mx1, b = y1 - k*x1;
				double ke2 = k*k, bk = b*k;
				double ke2p1 = ke2 + 1.;

				double AtanX1 = atan(k+b/x1);
				double AtanX2 = atan(k+b/x2);
				double SumAtans1 = atan(TransAtans((bk+ke2p1*x2)/b, -(bk+ke2p1*x1)/b, PiMult1));

				SumAtans1 += PiMult1*PI;
				double Log1 = log(x1e2 + y1e2);
				double Log2 = log(x2e2 + y2e2);
				double SumLogs1 = Log2 - Log1;
				double bdke2p1 = b/ke2p1;
				double bkdke2p1 = bdke2p1*k;

				Sx += -2.*((x2*AtanX2 - x1*AtanX1) - bkdke2p1*SumAtans1) - bdke2p1*SumLogs1;
				Sy += -2.*bdke2p1*SumAtans1 - ((bkdke2p1+x2)*Log2 - (bkdke2p1+x1)*Log1);
			}
			else if(SpecCaseID==ZeroVxVz)
			{
				double k = y2my1/x2mx1, b = y1 - k*x1;
				double kz = k*z;
				double SumAtans1 = atan(TransAtans(x2/z, -x1/z, PiMult1));
				SumAtans1 += PiMult1*PI;
				double SumLogs1 = log((x2e2 + ze2)/(x1e2 + ze2));

				Sx += 2.*(-y2my1 + kz*SumAtans1) - b*SumLogs1;
				Sz += -2.*b*SumAtans1 - kz*SumLogs1;
			}
		}
		if(abs_y2my1*Max_k > abs_x2mx1)
		{
			if(SpecCaseID==ZeroVyVz)
			{
				double k1 = x2mx1/y2my1, b1 = x1 - k1*y1;
				double k1z = k1*z;
				double SumAtans1 = atan(TransAtans(y2/z, -y1/z, PiMult1));
				SumAtans1 += PiMult1*PI;
				double SumLogs1 = log((y2e2 + ze2)/(y1e2 + ze2));

				Sy -= 2.*(-x2mx1 + k1z*SumAtans1) - b1*SumLogs1;
				Sz -= -2.*b1*SumAtans1 - k1z*SumLogs1;
			}
		}
		x1 = x2; y1 = y2;
		x1e2 = x2e2; y1e2 = y2e2;
	}
	TVector3d BufIh(Sx, Sy, Sz);
	double MultIh = -ConstForH*Magn.z; // Check if "-" is necessary
	BufIh = MultIh*BufIh;

	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIh;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIh;
}

//-------------------------------------------------------------------------

//void radTPolygon::B_comp_frJ(radTField* pField)
//{
//	radTVect2dVect::iterator baseIter = EdgePointsVector.begin();
//
//	for(int i=0; i<AmOfEdgePoints; i++)
//	{
//		++baseIter;
//	}
//}

//-------------------------------------------------------------------------
