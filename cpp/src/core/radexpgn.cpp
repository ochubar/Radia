/*-------------------------------------------------------------------------
*
* File name:      radexpgn.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source: extruded polygon (prism)
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
#include "radexpgn.h"
#include "radg3dgr.h"
#include "radsbdep.h"
#include "radg3da1.h"

//-------------------------------------------------------------------------

extern radTYield radYield;
extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTExtrPolygon::Dump(std::ostream& o, int ShortSign) // Porting
{
	radTg3dRelax::Dump(o);
	DumpPureObjInfo(o, ShortSign);
	if(ShortSign==1) return;

	DumpMaterApplied(o);
	DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTExtrPolygon::DumpPureObjInfo(std::ostream& o, int ShortSign)
{
	o << "Relaxable: ";
	o << "ThckPgn";

	if(ShortSign==1) return;

	o << endl;
	o << "   {x,y,z}= {" << CentrPoint.x << ',' << CentrPoint.y << ',' << CentrPoint.z << "}" << endl;
	o << "   lx= " << Thickness << endl;

	o << "   Pgn= {";
	radTPolygon* BasePgnPtr = (radTPolygon*)(BasePolygonHandle.rep);
	int AmOfEdgePoints_m_1 = BasePgnPtr->AmOfEdgePoints - 1;
	for(int i=0; i<=AmOfEdgePoints_m_1; i++)
	{
		o << "{" << BasePgnPtr->EdgePointsVector[i].x << ',' << BasePgnPtr->EdgePointsVector[i].y << "}";
		if(i!=AmOfEdgePoints_m_1) { o << ",";}
	}
	o << "}" << endl;

	o << "   {mx,my,mz}= {" << Magn.x << ',' << Magn.y << ',' << Magn.z << "}";
}

//-------------------------------------------------------------------------

void radTExtrPolygon::DumpBin_ExtrPolygon(CAuxBinStrVect& oStr)
{
	//TVector3d FirstPoint;
	oStr << FirstPoint;

	//TAxisOrient AxOrnt;
	char cAxOrnt = AxOrnt;
	oStr << cAxOrnt;

	//radThg BasePolygonHandle;
	radTPolygon* pBasePgn = (radTPolygon*)(BasePolygonHandle.rep);
	char cBasePgnDefined = (pBasePgn != 0);
	oStr << cBasePgnDefined;
	if(cBasePgnDefined) pBasePgn->DumpBin_Polygon(oStr);

	//double Thickness;
	oStr << Thickness;
}

//-------------------------------------------------------------------------

void radTExtrPolygon::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//Dumping objects that may be used by this object
	vector<pair<int, int> > vTrfKeys;
	DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys);

	int matKey=0;
	DumpBin_g3dRelax_TreatMat(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, matKey);

	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)Type_g3d();
	oStr << (char)Type_g3dRelax();
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTg3d
	DumpBin_g3d(oStr, vTrfKeys);

	//Members of radTg3dRelax
	DumpBin_g3dRelax(oStr, matKey);

	//Members of radTExtrPolygon
	DumpBin_ExtrPolygon(oStr);
}

//-------------------------------------------------------------------------

void radTExtrPolygon::DumpBinParse_ExtrPolygon(CAuxBinStrVect& inStr)
{
	//TVector3d FirstPoint;
	inStr >> FirstPoint;

	//TAxisOrient AxOrnt;
	char cAxOrnt = 0;
	inStr >> cAxOrnt;
	if(cAxOrnt == 0) AxOrnt = ParallelToX;
	else if(cAxOrnt == 1) AxOrnt = ParallelToY;
	else if(cAxOrnt == 2) AxOrnt = ParallelToZ;

	//radThg BasePolygonHandle;
	char cBasePgnDefined = 0;
	inStr >> cBasePgnDefined;
	if(cBasePgnDefined)
	{
		radThg hg(new radTPolygon(inStr));
		BasePolygonHandle = hg;
	}

	//double Thickness;
	inStr >> Thickness;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTExtrPolygon::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTExtrPolygonGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

void radTExtrPolygon::B_comp(radTField* FieldPtr) 
{
	// Orientation: The prism exis parallel to X !!!
	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;

	const double Max_k = 1.E+08; // To skip segments in general field computation loop.

	double AbsRandX = radCR.AbsRandMagnitude(FirstPoint.x - CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(FirstPoint.y - CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(FirstPoint.z - CentrPoint.z);
	double RelRandMagn = radCR.AbsRandMagnitude(1.);

	if(radYield.Check()==0) return; // To allow multitasking on Mac: consider better places for this

	radTPolygon* BasePolygonPtr = (radTPolygon*)(BasePolygonHandle.rep);
	int AmOfEdPoInBase = BasePolygonPtr->AmOfEdgePoints;

	TVector3d& ObsPo = FieldPtr->P;

	double z1 = FirstPoint.x - ObsPo.x;  // If this is zero, we get error !
	double z2 = z1 + Thickness;  // If this is zero, we get error !

// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(z1==0.) z1 = AbsRandX;
	if(z2==0.) z2 = AbsRandX;

	double z1e2 = z1*z1, z2e2 = z2*z2;
	double absz1 = Abs(z1), absz2 = Abs(z2);

#ifdef __GCC__
	vector<TVector2d>::iterator BaseIter = (BasePolygonPtr->EdgePointsVector).begin();
#else
	vector<TVector2d, allocator<TVector2d> >::iterator BaseIter = (BasePolygonPtr->EdgePointsVector).begin();
#endif
	
	TVector2d First2d(FirstPoint.y - ObsPo.y, FirstPoint.z - ObsPo.z);
	TVector2d Vect2dToAdd(First2d.x - (*BaseIter).x, First2d.y - (*BaseIter).y);

	double x1 = First2d.x;  // If this is zero, we get error ?
	double y1 = First2d.y;  // If this is zero, we get error ?

// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(x1==0.) x1 = AbsRandY;
	if(y1==0.) y1 = AbsRandZ;

	double x1e2 = x1*x1, y1e2 = y1*y1;
	double x2, y2, x2e2, y2e2;

	short A_CompNeeded = FieldPtr->FieldKey.A_;
	short B_orH_CompNeeded = FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_;
	short FieldCompNeeded = FieldPtr->FieldKey.PreRelax_ || B_orH_CompNeeded || A_CompNeeded;

	short InsideZ = (z1*z2<0);
	int X_In_Count=0, X_Out_Count=0;
	int Y_In_Count=0, Y_Out_Count=0;

	short OutsideGateX=1, OutsideGateY=1;

	short MagnCompNeeded = (FieldPtr->FieldKey.M_ || FieldPtr->FieldKey.B_) && InsideZ;
	// This is not final yet!

	int AmOfEdPoInBase_mi_1 = AmOfEdPoInBase - 1;

	double S11=0., S12=0., S13=0., S22=0., S23=0., S33=0.;
	double ArgSumAtans1=0., PiMultSumAtans1=0.;
	double ArgSumLogs2=1., ArgSumLogs4=1.;

	double SumLogs1, SumLogs3;
	double SumLogs4A, SumLogs5A, SumLogs6A, SumLogs7A;

	double AS12, AS13, AS23;
	double ArgSumAtans3A, ArgSumAtans4A, ArgSumAtans5A, ArgSumAtans6A, ArgSumAtans8A, ArgSumAtans9A,
		   PiMultSumAtans3A, PiMultSumAtans4A, PiMultSumAtans5A, PiMultSumAtans6A, PiMultSumAtans8A, PiMultSumAtans9A,
		   PhCorrSumAtans3A, PhCorrSumAtans4A;
	if(A_CompNeeded)
	{
		AS12 = AS13 = AS23 = 0.;
		ArgSumAtans3A = ArgSumAtans4A = ArgSumAtans5A = ArgSumAtans6A = ArgSumAtans8A = ArgSumAtans9A = 0.;
		PiMultSumAtans3A = PiMultSumAtans4A = PiMultSumAtans5A = PiMultSumAtans6A = PiMultSumAtans8A = PiMultSumAtans9A = 0.;
		PhCorrSumAtans3A = PhCorrSumAtans4A = 0.;
	}

	double four_be2ke2, four_be2be2ke2, be2mke2z2e2, be2pke2z2e2, be2mke2z2e2e2, be2pke2z2e2e2, 
		   DFlipRepSumAtans1, BufDen,
		   be2mke2z1e2, be2pke2z1e2, be2mke2z1e2e2, be2pke2z1e2e2,
		   Buf1Num, Buf2Num, xFlp1, xFlp2, xFlp,
		   xFlpe2, kxFlp, kxFlppb, kxFlpmb, kxFlppbe2, SqRoot;

	for(int i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x + Vect2dToAdd.x;  // If this is zero, we get error !
			y2 = (*BaseIter).y + Vect2dToAdd.y;  // If this is zero, we get error !
		}
		else
		{
			x2 = First2d.x;  // If this is zero, we get error !
			y2 = First2d.y;  // If this is zero, we get error !
		}

		// Artificial shift of an observation point a bit right of the block's border
		// if the point is exactly on the boarder (to avoid "divide by zero" error):
		// Removing this may be dangerous for the Checking-If-Inside
		if(x2==0.) x2 = AbsRandY;
		if(y2==0.) y2 = AbsRandZ;

		x2e2 = x2*x2; y2e2 = y2*y2;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = y2my1/x2mx1, b = y1 - k*x1;  // If b is zero, we get error !
			if(b==0.) b = AbsRandZ;

			if(MagnCompNeeded)
			{
				if(x1*x2 <= 0.)
				{
					short LocInside = ((x2>x1)? (b<=0) : (b>=0));
					X_In_Count += LocInside? 1 : 0;
					X_Out_Count += (!LocInside)? 1 : 0;

					OutsideGateX = 0;
				}
				if(y1*y2 <= 0.)	OutsideGateY = 0;
			}

			if(FieldCompNeeded)
			{
				double bk = b*k, ke2 = k*k, be2 = b*b, twob = 2.*b;
				double ke2p1 = ke2+1.;
				if(ke2p1==0.) ke2p1 = RelRandMagn;
				double sqrtke2p1 = sqrt(ke2p1);

				double kx1 = k*x1, kx2 = k*x2;
				double bpkx1 = b+kx1, bpkx2 = b+kx2;
				if(bpkx1==0.) bpkx1 = AbsRandZ;
				if(bpkx2==0.) bpkx2 = AbsRandZ;

				double bpkx1e2 = bpkx1*bpkx1, bpkx2e2 = bpkx2*bpkx2;
				double kx1mb = -b+kx1, kx2mb = -b+kx2; 
				//double R11 = sqrt(x1e2 + bpkx1e2 + z1e2), R12 = sqrt(x1e2 + bpkx1e2 + z2e2), R21 = sqrt(x2e2 + bpkx2e2 + z1e2), R22 = sqrt(x2e2 + bpkx2e2 + z2e2);
				double R11 = radCR.DoublePlus(sqrt(x1e2 + bpkx1e2 + z1e2)), R12 = radCR.DoublePlus(sqrt(x1e2 + bpkx1e2 + z2e2)), R21 = radCR.DoublePlus(sqrt(x2e2 + bpkx2e2 + z1e2)), R22 = radCR.DoublePlus(sqrt(x2e2 + bpkx2e2 + z2e2));

				double x1e2pz1e2 = x1e2+z1e2, x1e2pz2e2 = x1e2+z2e2, x2e2pz1e2 = x2e2+z1e2, x2e2pz2e2 = x2e2+z2e2;
				double bkpx1pke2x1 = bk+ke2p1*x1, bkpx2pke2x2 = bk+ke2p1*x2;  // If this is zero, we get error !
				if(bkpx1pke2x1==0.) bkpx1pke2x1 = AbsRandY;
				if(bkpx2pke2x2==0.) bkpx2pke2x2 = AbsRandY;

				double kz1e2 = k*z1e2, kz2e2 = k*z2e2;
				double ke2z1e2 = k*kz1e2, ke2z2e2 = k*kz2e2;
				double ke2z1e2mbe2 = ke2z1e2-be2, ke2z2e2mbe2 = ke2z2e2-be2, ke2z1e2pbe2 = ke2z1e2+be2, ke2z2e2pbe2 = ke2z2e2+be2;
				double bdz1 = b/z1, bdz2 = b/z2, bx1 = b*x1, bx2 = b*x2;
				double R11pbpkx1 = bpkx1+R11, R12pbpkx1 = bpkx1+R12, R21pbpkx2 = bpkx2+R21, R22pbpkx2 = bpkx2+R22;

				double AbsRandR11 = 100.*radCR.AbsRandMagnitude(R11);
				double AbsRandR12 = 100.*radCR.AbsRandMagnitude(R12);
				double AbsRandR21 = 100.*radCR.AbsRandMagnitude(R21);
				double AbsRandR22 = 100.*radCR.AbsRandMagnitude(R22);

				if(R11pbpkx1 < AbsRandR11) R11pbpkx1 = 0.5*(x1e2 + z1e2)/Abs(bpkx1);
				if(R12pbpkx1 < AbsRandR12) R12pbpkx1 = 0.5*(x1e2 + z2e2)/Abs(bpkx1);
				if(R21pbpkx2 < AbsRandR21) R21pbpkx2 = 0.5*(x2e2 + z1e2)/Abs(bpkx2);
				if(R22pbpkx2 < AbsRandR22) R22pbpkx2 = 0.5*(x2e2 + z2e2)/Abs(bpkx2);

				double bdz1R11 = bdz1*R11, bdz2R12 = bdz2*R12, bdz1R21 = bdz1*R21, bdz2R22 = bdz2*R22;
				double R11pz1 = R11+z1, R12pz2 = R12+z2, R21pz1 = R21+z1, R22pz2 = R22+z2;

				if(R11pz1 < AbsRandR11) R11pz1 = 0.5*(x1e2 + bpkx1e2)/Abs(z1);
				if(R12pz2 < AbsRandR12) R12pz2 = 0.5*(x1e2 + bpkx1e2)/Abs(z2);
				if(R21pz1 < AbsRandR21) R21pz1 = 0.5*(x2e2 + bpkx2e2)/Abs(z1);
				if(R22pz2 < AbsRandR22) R22pz2 = 0.5*(x2e2 + bpkx2e2)/Abs(z2);

				double bkpx1pke2x1dsqrtke2p1 = bkpx1pke2x1/sqrtke2p1;
				double R11_p_bkpx1pke2x1dsqrtke2p1 = R11 + bkpx1pke2x1dsqrtke2p1;
				double R12_p_bkpx1pke2x1dsqrtke2p1 = R12 + bkpx1pke2x1dsqrtke2p1;
				double bkpx2pke2x2dsqrtke2p1 = bkpx2pke2x2/sqrtke2p1;
				double R21_p_bkpx2pke2x2dsqrtke2p1 = R21 + bkpx2pke2x2dsqrtke2p1;
				double R22_p_bkpx2pke2x2dsqrtke2p1 = R22 + bkpx2pke2x2dsqrtke2p1;

				if(R11_p_bkpx1pke2x1dsqrtke2p1 < AbsRandR11) R11_p_bkpx1pke2x1dsqrtke2p1 = 0.5*(be2 + z1e2)/(Abs(x1)*sqrtke2p1);
				if(R12_p_bkpx1pke2x1dsqrtke2p1 < AbsRandR12) R12_p_bkpx1pke2x1dsqrtke2p1 = 0.5*(be2 + z2e2)/(Abs(x1)*sqrtke2p1);
				if(R21_p_bkpx2pke2x2dsqrtke2p1 < AbsRandR21) R21_p_bkpx2pke2x2dsqrtke2p1 = 0.5*(be2 + z1e2)/(Abs(x2)*sqrtke2p1);
				if(R22_p_bkpx2pke2x2dsqrtke2p1 < AbsRandR22) R22_p_bkpx2pke2x2dsqrtke2p1 = 0.5*(be2 + z2e2)/(Abs(x2)*sqrtke2p1);

				double FlpRep1ForSumAtans1 = 0., FlpRep1z1ForSumAtansA1 = 0., FlpRep1z2ForSumAtansA1 = 0.;

				four_be2ke2 = 4.*be2*ke2; four_be2be2ke2 = be2*four_be2ke2;

				be2mke2z2e2 = be2-ke2z2e2; be2pke2z2e2 = be2+ke2z2e2;
				be2mke2z2e2e2 = be2mke2z2e2*be2mke2z2e2; be2pke2z2e2e2 = be2pke2z2e2*be2pke2z2e2;
				DFlipRepSumAtans1 = (be2+ke2p1*z2e2)*(four_be2ke2*(be2+ke2z2e2)-be2mke2z2e2e2);
				BufDen = four_be2be2ke2-ke2p1*be2mke2z2e2e2;

				if((DFlipRepSumAtans1 >= 0.) && BufDen!=0.)
				{
					Buf1Num = bk*be2pke2z2e2e2;
					Buf2Num = be2mke2z2e2*sqrt(DFlipRepSumAtans1);
					xFlp1 = (Buf1Num - Buf2Num)/BufDen;
					xFlp2 = (Buf1Num + Buf2Num)/BufDen;

					xFlp = xFlp1;
					if((x1<x2)? ((xFlp>x1) && (xFlp<x2)) : ((xFlp<x1) && (xFlp>x2)))
					{
						xFlpe2 = xFlp*xFlp; kxFlp = k*xFlp;
						kxFlppb = kxFlp+b; kxFlpmb = kxFlp-b;
						kxFlppbe2 = kxFlppb*kxFlppb;
						SqRoot = sqrt(xFlpe2+kxFlppbe2+z2e2);
						if(Sign((xFlpe2+z2e2)*(-be2mke2z2e2) + (-be2+ke2*xFlpe2)*be2pke2z2e2) == Sign(-kxFlpmb)) // RootIsReal?
						{
							double DenomDerivSign = Sign(-2.*xFlp*be2mke2z2e2 + kxFlpmb*be2pke2z2e2*(k+(bk+ke2p1*xFlp)/SqRoot) + k*be2pke2z2e2*(kxFlppb + SqRoot));
							double NumSign = Sign((2.*bk*z2e2*(xFlpe2+z2e2) + (b*xFlp+kz2e2)*be2pke2z2e2*(kxFlppb + SqRoot))/z2);
							double Buf = DenomDerivSign*NumSign*Sign(x2mx1);

							FlpRep1ForSumAtans1 += Buf;
							FlpRep1z2ForSumAtansA1 -= Buf;
						}
					}

					xFlp = xFlp2;
					if((x1<x2)? ((xFlp>x1) && (xFlp<x2)) : ((xFlp<x1) && (xFlp>x2)))
					{
						xFlpe2 = xFlp*xFlp; kxFlp = k*xFlp;
						kxFlppb = kxFlp+b; kxFlpmb = kxFlp-b;
						kxFlppbe2 = kxFlppb*kxFlppb;
						SqRoot = sqrt(xFlpe2+kxFlppbe2+z2e2);
						if(Sign((xFlpe2+z2e2)*(-be2mke2z2e2) + (-be2+ke2*xFlpe2)*be2pke2z2e2) == Sign(-kxFlpmb))
						{
							double DenomDerivSign = Sign(-2.*xFlp*be2mke2z2e2 + kxFlpmb*be2pke2z2e2*(k+(bk+ke2p1*xFlp)/SqRoot) + k*be2pke2z2e2*(kxFlppb + SqRoot));
							double NumSign = Sign((2.*bk*z2e2*(xFlpe2+z2e2) + (b*xFlp+kz2e2)*be2pke2z2e2*(kxFlppb + SqRoot))/z2);
							double Buf = DenomDerivSign*NumSign*Sign(x2mx1);

							FlpRep1ForSumAtans1 += Buf;
							FlpRep1z2ForSumAtansA1 -= Buf;
						}
					}
				}

				be2mke2z1e2 = be2-ke2z1e2; be2pke2z1e2 = be2+ke2z1e2;
				be2mke2z1e2e2 = be2mke2z1e2*be2mke2z1e2; be2pke2z1e2e2 = be2pke2z1e2*be2pke2z1e2;
				DFlipRepSumAtans1 = (be2+ke2p1*z1e2)*(four_be2ke2*(be2+ke2z1e2)-be2mke2z1e2e2);
				BufDen = four_be2be2ke2-ke2p1*be2mke2z1e2e2;

				if((DFlipRepSumAtans1 >= 0.) && BufDen!=0.)
				{
					Buf1Num = bk*be2pke2z1e2e2;
					Buf2Num = be2mke2z1e2*sqrt(DFlipRepSumAtans1);
					xFlp1 = (Buf1Num - Buf2Num)/BufDen;
					xFlp2 = (Buf1Num + Buf2Num)/BufDen;

					xFlp = xFlp1;
					if((x1<x2)? ((xFlp>x1) && (xFlp<x2)) : ((xFlp<x1) && (xFlp>x2)))
					{
						xFlpe2 = xFlp*xFlp; kxFlp = k*xFlp;
						kxFlppb = kxFlp+b; kxFlpmb = kxFlp-b;
						kxFlppbe2 = kxFlppb*kxFlppb;
						SqRoot = sqrt(xFlpe2+kxFlppbe2+z1e2);
						if(Sign((xFlpe2+z1e2)*(-be2mke2z1e2) + (-be2+ke2*xFlpe2)*be2pke2z1e2) == Sign(-kxFlpmb))
						{
							double DenomDerivSign = Sign(-2.*xFlp*be2mke2z1e2 + kxFlpmb*be2pke2z1e2*(k+(bk+ke2p1*xFlp)/SqRoot) + k*be2pke2z1e2*(kxFlppb + SqRoot));
							double NumSign = Sign((2.*bk*z1e2*(xFlpe2+z1e2) + (b*xFlp+kz1e2)*be2pke2z1e2*(kxFlppb + SqRoot))/z1);
							double Buf = DenomDerivSign*NumSign*Sign(x2mx1);

							FlpRep1ForSumAtans1 -= Buf;
							FlpRep1z1ForSumAtansA1 += Buf;
						}
					}

					xFlp = xFlp2;
					if((x1<x2)? ((xFlp>x1) && (xFlp<x2)) : ((xFlp<x1) && (xFlp>x2)))
					{
						xFlpe2 = xFlp*xFlp; kxFlp = k*xFlp;
						kxFlppb = kxFlp+b; kxFlpmb = kxFlp-b;
						kxFlppbe2 = kxFlppb*kxFlppb;
						SqRoot = sqrt(xFlpe2+kxFlppbe2+z1e2);
						if(Sign((xFlpe2+z1e2)*(-be2mke2z1e2) + (-be2+ke2*xFlpe2)*be2pke2z1e2) == Sign(-kxFlpmb))
						{
							double DenomDerivSign = Sign(-2.*xFlp*be2mke2z1e2 + kxFlpmb*be2pke2z1e2*(k+(bk+ke2p1*xFlp)/SqRoot) + k*be2pke2z1e2*(kxFlppb + SqRoot));
							double NumSign = Sign((2.*bk*z1e2*(xFlpe2+z1e2) + (b*xFlp+kz1e2)*be2pke2z1e2*(kxFlppb + SqRoot))/z1);
							double Buf = DenomDerivSign*NumSign*Sign(x2mx1);

							FlpRep1ForSumAtans1 -= Buf;
							FlpRep1z1ForSumAtansA1 += Buf;
						}
					}
				}

				double mbkdke2p1 = -bk/ke2p1;
				double FlpRep2ForSumAtans1 = (1./ke2p1)*Sign(b)*(Sign(z2)-Sign(z1))*(Step(mbkdke2p1-x1)*Step(x2-mbkdke2p1)*Step(x2mx1) - Step(x1-mbkdke2p1)*Step(mbkdke2p1-x2)*Step(-x2mx1));
				double FlpRep3ForSumAtans1 = (Step(-x1)*Step(x2)*Step(x2mx1) - Step(x1)*Step(-x2)*Step(-x2mx1))*Step(-z1)*Step(z2)*((-2.- Sign(b)*(Sign(z2)-Sign(z1)))*(Step(-y1)*Step(y2)*Step(y2my1) + Step(y1)*Step(-y2)*Step(-y2my1)) - 4.*Step((y1<y2)? y1 : y2));

				double PiMult1=0., PiMult2=0., PiMult3=0., PiMult4=0.;

				if(A_CompNeeded)
				{
					SumLogs4A = log(R12pz2/R11pz1);
					SumLogs5A = log(R21pz1/R22pz2);

					SumLogs6A = log(R21_p_bkpx2pke2x2dsqrtke2p1/R11_p_bkpx1pke2x1dsqrtke2p1);
					SumLogs7A = log(R12_p_bkpx1pke2x1dsqrtke2p1/R22_p_bkpx2pke2x2dsqrtke2p1);

					SumLogs1 = SumLogs4A + SumLogs5A;
					SumLogs3 = -SumLogs6A - SumLogs7A;
				}
				else if(B_orH_CompNeeded)
				{
					SumLogs1 = log((R21pz1*R12pz2)/(R11pz1*R22pz2));
					SumLogs3 = log(R11_p_bkpx1pke2x1dsqrtke2p1*R22_p_bkpx2pke2x2dsqrtke2p1/(R12_p_bkpx1pke2x1dsqrtke2p1*R21_p_bkpx2pke2x2dsqrtke2p1));
				}

				double Arg1ForSumAtans1 = -(ke2z1e2pbe2*(bx1 + kz1e2)*R11pbpkx1 + kz1e2*twob*x1e2pz1e2);
				double Arg2ForSumAtans1 = (ke2z1e2pbe2*kx1mb*R11pbpkx1 + ke2z1e2mbe2*x1e2pz1e2)*z1;
				double Arg3ForSumAtans1 = ke2z1e2pbe2*(bx2 + kz1e2)*R21pbpkx2 + kz1e2*twob*x2e2pz1e2;
				double Arg4ForSumAtans1 = (ke2z1e2pbe2*kx2mb*R21pbpkx2 + ke2z1e2mbe2*x2e2pz1e2)*z1;
				double Arg5ForSumAtans1 = ke2z2e2pbe2*(bx1 + kz2e2)*R12pbpkx1 + kz2e2*twob*x1e2pz2e2;
				double Arg6ForSumAtans1 = (ke2z2e2pbe2*kx1mb*R12pbpkx1 + ke2z2e2mbe2*x1e2pz2e2)*z2;
				double Arg7ForSumAtans1 = -(ke2z2e2pbe2*(bx2 + kz2e2)*R22pbpkx2 + kz2e2*twob*x2e2pz2e2);
				double Arg8ForSumAtans1 = (ke2z2e2pbe2*kx2mb*R22pbpkx2 + ke2z2e2mbe2*x2e2pz2e2)*z2;

				if(B_orH_CompNeeded)
				{
					double BufTransAtans1 = TransAtans(Arg1ForSumAtans1/Arg2ForSumAtans1, Arg3ForSumAtans1/Arg4ForSumAtans1, PiMult1);
					double BufTransAtans2 = TransAtans(Arg5ForSumAtans1/Arg6ForSumAtans1, Arg7ForSumAtans1/Arg8ForSumAtans1, PiMult2);
					double CurArgSumAtans1 = TransAtans(BufTransAtans1, BufTransAtans2, PiMult3);

					ArgSumAtans1 = TransAtans(ArgSumAtans1, CurArgSumAtans1, PiMult4);
					PiMultSumAtans1 += PiMult1+PiMult2+PiMult3+PiMult4 + FlpRep1ForSumAtans1;

					double SumAtans2 = atan(TransAtans(TransAtans(-bdz1R11/bkpx1pke2x1, bdz2R12/bkpx1pke2x1, PiMult1), TransAtans(bdz1R21/bkpx2pke2x2, -bdz2R22/bkpx2pke2x2, PiMult2), PiMult3));
					SumAtans2 += (PiMult1+PiMult2+PiMult3)*PI;
					double PhCorrSumAtans2 = PhCorrForArcTanTwo(bkpx1pke2x1, -bdz1R11) + PhCorrForArcTanTwo(bkpx1pke2x1, bdz2R12) + PhCorrForArcTanTwo(bkpx2pke2x2, bdz1R21) + PhCorrForArcTanTwo(bkpx2pke2x2, -bdz2R22);
					double SumAtans3 = atan(TransAtans(TransAtans(bkpx1pke2x1/bdz1R11, -bkpx1pke2x1/bdz2R12, PiMult1), TransAtans(-bkpx2pke2x2/bdz1R21, bkpx2pke2x2/bdz2R22, PiMult2), PiMult3));
					SumAtans3 += (PiMult1+PiMult2+PiMult3)*PI;

					ArgSumLogs2 *= R11pbpkx1*R22pbpkx2/(R12pbpkx1*R21pbpkx2);
					ArgSumLogs4 *= x1e2pz2e2*x2e2pz1e2/(x1e2pz1e2*x2e2pz2e2);

					S11 += SumAtans2/ke2p1 + k*SumLogs1/ke2p1 + PI*(FlpRep2ForSumAtans1 + FlpRep3ForSumAtans1);
					S12 += (k*SumAtans3 - SumLogs1)/ke2p1;
					S13 -= k*SumLogs3/sqrtke2p1;
					S22 += (-(SumAtans2 + PhCorrSumAtans2) - k*SumLogs1)/ke2p1;
					S23 += SumLogs3/sqrtke2p1;
				}
				if(A_CompNeeded)
				{
					double PiMult1=0., PiMult2=0.;
					double CurArgSumAtans3A = TransAtans(-Arg1ForSumAtans1/Arg2ForSumAtans1, -Arg3ForSumAtans1/Arg4ForSumAtans1, PiMult1);
					ArgSumAtans3A = (ArgSumAtans3A==0.)? CurArgSumAtans3A : TransAtans(ArgSumAtans3A, CurArgSumAtans3A, PiMult2);
					PiMultSumAtans3A += PiMult1 + PiMult2 + FlpRep1z1ForSumAtansA1;

					PiMult1=0., PiMult2=0.;
					double CurArgSumAtans4A = TransAtans(-Arg5ForSumAtans1/Arg6ForSumAtans1, -Arg7ForSumAtans1/Arg8ForSumAtans1, PiMult1);
					ArgSumAtans4A = (ArgSumAtans4A==0.)? CurArgSumAtans4A : TransAtans(ArgSumAtans4A, CurArgSumAtans4A, PiMult2);
					PiMultSumAtans4A += PiMult1 + PiMult2 + FlpRep1z2ForSumAtansA1;

					double SumLogs1A = B_orH_CompNeeded? SumLogs3 : log(R11_p_bkpx1pke2x1dsqrtke2p1*R22_p_bkpx2pke2x2dsqrtke2p1/(R12_p_bkpx1pke2x1dsqrtke2p1*R21_p_bkpx2pke2x2dsqrtke2p1));

					double SumLogs2A = log(R11pbpkx1/R12pbpkx1);
					double SumLogs3A = log(R22pbpkx2/R21pbpkx2);

					AS12 += b*SumLogs1A/sqrtke2p1 + x1*SumLogs2A + x2*SumLogs3A;

					double FirstArgForSumAtansA = bkpx1pke2x1/bdz1;
					double SecondArgForSumAtansA = -R11;
					double ThirdArgForSumAtansA = bkpx2pke2x2/bdz1;
					double FourthArgForSumAtansA = R21;
					double CurArgSumAtans7A1 = TransAtans(SecondArgForSumAtansA/FirstArgForSumAtansA, FourthArgForSumAtansA/ThirdArgForSumAtansA, PiMult1);
					double CurPhCorrSumAtans7A1 = PhCorrForArcTanTwo(FirstArgForSumAtansA, SecondArgForSumAtansA) + PhCorrForArcTanTwo(ThirdArgForSumAtansA, FourthArgForSumAtansA);
					FirstArgForSumAtansA = bkpx1pke2x1/bdz2;
					SecondArgForSumAtansA = R12;
					ThirdArgForSumAtansA = bkpx2pke2x2/bdz2;
					FourthArgForSumAtansA = -R22;
					double CurArgSumAtans7A2 = TransAtans(SecondArgForSumAtansA/FirstArgForSumAtansA, FourthArgForSumAtansA/ThirdArgForSumAtansA, PiMult2);
					double CurPhCorrSumAtans7A2 = PhCorrForArcTanTwo(FirstArgForSumAtansA, SecondArgForSumAtansA) + PhCorrForArcTanTwo(ThirdArgForSumAtansA, FourthArgForSumAtansA);
					double SumAtans7A = atan(TransAtans(CurArgSumAtans7A1, CurArgSumAtans7A2, PiMult3));
					SumAtans7A += (PiMult1 + PiMult2 + PiMult3)*PI + CurPhCorrSumAtans7A1 + CurPhCorrSumAtans7A2;

					AS13 += (b*SumAtans7A + bkpx1pke2x1*SumLogs4A + bkpx2pke2x2*SumLogs5A)/ke2p1
						  + (z1*SumLogs6A + z2*SumLogs7A)/sqrtke2p1;
				}
			}
		}
		else
		{
			if(MagnCompNeeded) 
				if(y1*y2 <= 0.)
				{
					short LocInside = ((y2>y1)? (x1>=0) : (x1<=0));
					Y_In_Count += LocInside? 1 : 0;
					Y_Out_Count += (!LocInside)? 1 : 0;

					OutsideGateY = 0;
				}
		}

		if(abs_y2my1*Max_k > abs_x2mx1)
		{
			if(A_CompNeeded)
			{
				double k1 = (x2-x1)/(y2-y1), b1 = x1 - k1*y1;  // If b1 is zero, we get error !
				if(b1==0.) b1 = AbsRandY;

				double b1k1 = b1*k1, k1e2p1 = k1*k1+1;
				if(k1e2p1==0.) k1e2p1 = RelRandMagn;

				double sqrtk1e2p1 = sqrt(k1e2p1);
				double b1k1py1pk1e2y1 = b1k1+k1e2p1*y1, b1k1py2pk1e2y2 = b1k1+k1e2p1*y2;  // If this is zero, we get error !
				double b1pk1y1 = b1+k1*y1, b1pk1y2 = b1+k1*y2;
				double b1pk1y1e2 = b1pk1y1*b1pk1y1, b1pk1y2e2 = b1pk1y2*b1pk1y2;
				double RY11 = radCR.DoublePlus(sqrt(y1e2 + b1pk1y1e2 + z1e2)), RY12 = radCR.DoublePlus(sqrt(y1e2 + b1pk1y1e2 + z2e2)), RY21 = radCR.DoublePlus(sqrt(y2e2 + b1pk1y2e2 + z1e2)), RY22 = radCR.DoublePlus(sqrt(y2e2 + b1pk1y2e2 + z2e2));
				double RY11pz1 = RY11+z1, RY12pz2 = RY12+z2, RY21pz1 = RY21+z1, RY22pz2 = RY22+z2;

				double AbsRandRY11 = 100.*radCR.AbsRandMagnitude(RY11);
				double AbsRandRY12 = 100.*radCR.AbsRandMagnitude(RY12);
				double AbsRandRY21 = 100.*radCR.AbsRandMagnitude(RY21);
				double AbsRandRY22 = 100.*radCR.AbsRandMagnitude(RY22);
				
				if(RY11pz1 < AbsRandRY11) RY11pz1 = 0.5*(y1e2 + b1pk1y1e2)/Abs(z1);
				if(RY12pz2 < AbsRandRY12) RY12pz2 = 0.5*(y1e2 + b1pk1y1e2)/Abs(z2);
				if(RY21pz1 < AbsRandRY21) RY21pz1 = 0.5*(y2e2 + b1pk1y2e2)/Abs(z1);
				if(RY22pz2 < AbsRandRY22) RY22pz2 = 0.5*(y2e2 + b1pk1y2e2)/Abs(z2);

				double b1k1py1pk1e2y1dsqrtk1e2p1 = b1k1py1pk1e2y1/sqrtk1e2p1;
				double RY11_p_b1k1py1pk1e2y1dsqrtk1e2p1 = RY11 + b1k1py1pk1e2y1dsqrtk1e2p1;
				double RY12_p_b1k1py1pk1e2y1dsqrtk1e2p1 = RY12 + b1k1py1pk1e2y1dsqrtk1e2p1;
				double b1k1py2pk1e2y2dsqrtk1e2p1 = b1k1py2pk1e2y2/sqrtk1e2p1;
				double RY21_p_b1k1py2pk1e2y2dsqrtk1e2p1 = RY21 + b1k1py2pk1e2y2dsqrtk1e2p1;
				double RY22_p_b1k1py2pk1e2y2dsqrtk1e2p1 = RY22 + b1k1py2pk1e2y2dsqrtk1e2p1;

				if(RY11_p_b1k1py1pk1e2y1dsqrtk1e2p1 < AbsRandRY11) RY11_p_b1k1py1pk1e2y1dsqrtk1e2p1 = 0.5*(b1*b1 + z1e2)/(Abs(y1)*sqrtk1e2p1);
				if(RY12_p_b1k1py1pk1e2y1dsqrtk1e2p1 < AbsRandRY12) RY12_p_b1k1py1pk1e2y1dsqrtk1e2p1 = 0.5*(b1*b1 + z2e2)/(Abs(y1)*sqrtk1e2p1);
				if(RY21_p_b1k1py2pk1e2y2dsqrtk1e2p1 < AbsRandRY21) RY21_p_b1k1py2pk1e2y2dsqrtk1e2p1 = 0.5*(b1*b1 + z1e2)/(Abs(y2)*sqrtk1e2p1);
				if(RY22_p_b1k1py2pk1e2y2dsqrtk1e2p1 < AbsRandRY22) RY22_p_b1k1py2pk1e2y2dsqrtk1e2p1 = 0.5*(b1*b1 + z2e2)/(Abs(y2)*sqrtk1e2p1);

				double b1dz1 = b1/z1, b1dz2 = b1/z2;
				double PiMult1=0., PiMult2=0., PiMult3=0.;

				double FirstArgForSumAtansA = b1k1py1pk1e2y1/b1dz1;
				double SecondArgForSumAtansA = RY11;
				double ThirdArgForSumAtansA = b1k1py2pk1e2y2/b1dz1;
				double FourthArgForSumAtansA = -RY21;
				double CurArgSumAtans10A1 = TransAtans(SecondArgForSumAtansA/FirstArgForSumAtansA, FourthArgForSumAtansA/ThirdArgForSumAtansA, PiMult1);
				double CurPhCorrSumAtans10A1 = PhCorrForArcTanTwo(FirstArgForSumAtansA, SecondArgForSumAtansA) + PhCorrForArcTanTwo(ThirdArgForSumAtansA, FourthArgForSumAtansA);
				FirstArgForSumAtansA = b1k1py1pk1e2y1/b1dz2;
				SecondArgForSumAtansA = -RY12;
				ThirdArgForSumAtansA = b1k1py2pk1e2y2/b1dz2;
				FourthArgForSumAtansA = RY22;
				double CurArgSumAtans10A2 = TransAtans(SecondArgForSumAtansA/FirstArgForSumAtansA, FourthArgForSumAtansA/ThirdArgForSumAtansA, PiMult2);
				double CurPhCorrSumAtans10A2 = PhCorrForArcTanTwo(FirstArgForSumAtansA, SecondArgForSumAtansA) + PhCorrForArcTanTwo(ThirdArgForSumAtansA, FourthArgForSumAtansA);
				double SumAtans10A = atan(TransAtans(CurArgSumAtans10A1, CurArgSumAtans10A2, PiMult3));
				SumAtans10A += (PiMult1 + PiMult2 + PiMult3)*PI + CurPhCorrSumAtans10A1 + CurPhCorrSumAtans10A2;

				double SumLogs8A = log((RY11pz1)/(RY12pz2));
				double SumLogs9A = log((RY22pz2)/(RY21pz1));
				double SumLogs10A = log(RY11_p_b1k1py1pk1e2y1dsqrtk1e2p1/RY21_p_b1k1py2pk1e2y2dsqrtk1e2p1);
				double SumLogs11A = log(RY22_p_b1k1py2pk1e2y2dsqrtk1e2p1/RY12_p_b1k1py1pk1e2y1dsqrtk1e2p1);

				AS23 += (b1*SumAtans10A + b1k1py1pk1e2y1*SumLogs8A + b1k1py2pk1e2y2*SumLogs9A)/k1e2p1
					  + (z1*SumLogs10A + z2*SumLogs11A)/sqrtk1e2p1;
			}
		}
		x1 = x2; y1 = y2;
		x1e2 = x2e2; y1e2 = y2e2;
	}
	if(FieldCompNeeded)
	{
		if(B_orH_CompNeeded)
		{
			double SumAtans1 = atan(ArgSumAtans1) + PiMultSumAtans1*PI;
			S11 += SumAtans1;
			S13 += log(ArgSumLogs2) + 0.5*log(ArgSumLogs4);
			S33 -= SumAtans1;
		}
		if(A_CompNeeded)
		{
			AS12 += z1*(atan(ArgSumAtans3A) + PiMultSumAtans3A*PI) 
				  + z2*(atan(ArgSumAtans4A) + PiMultSumAtans4A*PI);
		}
	}

	int X_In = X_In_Count - X_Out_Count;
	int Y_In = X_In_Count - X_Out_Count;
	short InsideBlock = InsideZ? (CheckIfPosEven(X_In) && CheckIfPosEven(Y_In) && (!OutsideGateX) && (!OutsideGateY)) : 0;

	MagnCompNeeded = MagnCompNeeded && InsideBlock; // This is final

	if(MagnCompNeeded) FieldPtr->M += Magn;
	if(FieldCompNeeded)
	{
		double& mx = Magn.y;
		double& my = Magn.z;
		double& mz = Magn.x;

		if(FieldPtr->FieldKey.PreRelax_)
		{
			TVector3d Str0(S33, S13, S23);
			TVector3d Str1(S13, S11, S12);
			TVector3d Str2(S23, S12, S22);

			FieldPtr->B = (-ConstForH)*Str0; 
			FieldPtr->H = (-ConstForH)*Str1; 
			FieldPtr->A = (-ConstForH)*Str2; 
			return;
		}
		if(B_orH_CompNeeded)
		{
			TVector3d BufH(-ConstForH*(S13*mx+S23*my+S33*mz),
						   -ConstForH*(S11*mx+S12*my+S13*mz),
						   -ConstForH*(S12*mx+S22*my+S23*mz));

			if(FieldPtr->FieldKey.H_) FieldPtr->H += BufH;
			if(FieldPtr->FieldKey.B_)
			{
				FieldPtr->B += BufH;
				if(InsideBlock) FieldPtr->B += Magn;
			}
		}
		if(FieldPtr->FieldKey.A_)
		{
			TVector3d BufA;

			double& ax = BufA.y;
			double& ay = BufA.z;
			double& az = BufA.x;
			
			ax = ConstForH*(-AS12*my-AS13*mz);
			ay = ConstForH*(AS12*mx+AS23*mz);
			az = ConstForH*(AS13*mx-AS23*my);

			FieldPtr->A += BufA;
		}
	}
}

//-------------------------------------------------------------------------

void radTExtrPolygon::B_intComp(radTField* FieldPtr) 
{
// Orientation: The prism exis parallel to X !!!
// Field Integrals are correct only if line does not cross the prism body.
// This was noticed the "divide by zero" problem to take place if (Vx*Vx+Vy*Vy==0)||(Vx*Vx+Vz*Vz==0)||(Vy*Vy+Vz*Vz==0)

	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;

	const double Max_k = 1.E+08;

	const double ZeroToler = 1.E-06; // This is to switch to Special Cases
	const double SmallestRelTolerV = 1.E-12; // Relative tolerance to repair trapping V.i to zero at general case

	TVector3d V = FieldPtr->NextP - FieldPtr->P;

	double InvAbsV = 1./sqrt(V.x*V.x + V.y*V.y + V.z*V.z);
	V = InvAbsV*V;

	short InitVxIsZero = (Abs(V.x) <= ZeroToler);
	short InitVyIsZero = (Abs(V.y) <= ZeroToler);
	short InitVzIsZero = (Abs(V.z) <= ZeroToler);

// This is the attempt to avoid "divide by zero" problem
	TSpecCaseID SpecCaseID = NoSpecCase;
	if(InitVxIsZero && InitVyIsZero) SpecCaseID = ZeroVxVz;
	else if(InitVxIsZero && InitVzIsZero) SpecCaseID = ZeroVyVz;
	else if(InitVyIsZero && InitVzIsZero) SpecCaseID = ZeroVxVy;
	if(SpecCaseID==ZeroVxVy || SpecCaseID==ZeroVxVz || SpecCaseID==ZeroVyVz) { B_intCompSpecCases(FieldPtr, SpecCaseID); return;}

	double AbsRandX = radCR.AbsRandMagnitude(FirstPoint.x - CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(FirstPoint.y - CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(FirstPoint.z - CentrPoint.z);

	double Vx = V.y; if(Vx==0.) Vx = SmallestRelTolerV;
	double Vy = V.z; if(Vy==0.) Vy = SmallestRelTolerV;
	double Vz = V.x; if(Vz==0.) Vz = SmallestRelTolerV;

	double Vxe2 = Vx*Vx, Vye2 = Vy*Vy, Vze2 = Vz*Vz;
	double Vxe2pVye2 = Vxe2+Vye2, Vxe2pVze2 = Vxe2+Vze2, Vye2pVze2 = Vye2+Vze2;
	double VyVze2mVxe2 = Vy*(Vze2-Vxe2), VxVze2mVye2 = Vx*(Vze2-Vye2), Vxe2pVze2Vy = Vxe2pVze2*Vy;
	double VxVy = Vx*Vy, VxVz = Vx*Vz, VyVz = Vy*Vz;
	double VzdVxe2pVye2 = Vz/Vxe2pVye2, p2dVxe2pVye2 = 2./Vxe2pVye2;
	double p2Vxe2Vy = 2*Vxe2*Vy, m2Vz = -2.*Vz;

	TVector3d& StPo = FieldPtr->P;

	double z1 = FirstPoint.x - StPo.x;
	double z2 = z1 + Thickness;

// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(z1==0.) z1 = AbsRandX;
	if(z2==0.) z2 = AbsRandX;

	double Vxz1 = Vx*z1, Vxz2 = Vx*z2, Vyz1 = Vy*z1, Vyz2 = Vy*z2;
	double m2Vzz1 = m2Vz*z1, m2Vzz2 = m2Vz*z2;

	radTPolygon* BasePolygonPtr = (radTPolygon*)(BasePolygonHandle.rep);

#ifdef __GCC__
	vector<TVector2d>::iterator BaseIter = (BasePolygonPtr->EdgePointsVector).begin();
#else
	vector<TVector2d, allocator<TVector2d> >::iterator BaseIter = (BasePolygonPtr->EdgePointsVector).begin();
#endif

	int AmOfEdPoInBase = BasePolygonPtr->AmOfEdgePoints;
	int AmOfEdPoInBase_mi_1 = AmOfEdPoInBase - 1;
	
	TVector2d First2d(FirstPoint.y - StPo.y, FirstPoint.z - StPo.z);
	TVector2d Vect2dToAdd(First2d.x - (*BaseIter).x, First2d.y - (*BaseIter).y);

	double x1 = First2d.x;
	double y1 = First2d.y;

// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(x1==0.) x1 = AbsRandY;
	if(y1==0.) y1 = AbsRandZ;

	double x2, y2;

	double S11=0., S12=0., S13=0., S22=0., S23=0., S33=0.;
	double PiMult1=0.;

	for(int i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x + Vect2dToAdd.x;
			y2 = (*BaseIter).y + Vect2dToAdd.y;
		}
		else
		{
			x2 = First2d.x;
			y2 = First2d.y;
		}

		// Artificial shift of an observation point a bit right of the block's border
		// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(x2==0.) x2 = AbsRandY;
		if(y2==0.) y2 = AbsRandZ;

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
			double bVx = b*Vx, bVy = b*Vy, bVz = b*Vz;
			double bkVxVxe2pVze2pVyVze2mVxe2 = b*(kVx*Vxe2pVze2 + VyVze2mVxe2);
			double bVyVz = b*VyVz;
			double bkVxe2pVze2mVxVy = b*(k*Vxe2pVze2 - VxVy);
			double kVxmVye2p1pke2Vze2Vx = kVxmVye2p1pke2Vze2*Vx;
			double bVxmkVxe2pVze2mVxVyVy = b*(Vx - (k*Vxe2pVze2 - VxVy)*Vy);
			double kVxmVye2p1pke2Vze2Vy = kVxmVye2p1pke2Vze2*Vy, kVxe2pVze2mVxVy = k*Vxe2pVze2 - VxVy;
			double kVypVxVzx1 = kVypVxVz*x1, kVypVxVzx2 = kVypVxVz*x2;
			double kVxmVyx1 = kVxmVy*x1, kVxmVyx2 = kVxmVy*x2;
			double kVxmVye2p1pke2Vze2x1 = kVxmVye2p1pke2Vze2*x1, kVxmVye2p1pke2Vze2x2 = kVxmVye2p1pke2Vze2*x2;
			double kVypVxVzz1 = kVypVxVz*z1, kVypVxVzz2 = kVypVxVz*z2, kVxmVyz1 = kVxmVy*z1, kVxmVyz2 = kVxmVy*z2;
			double bVxpkVxmVyx1 = bVx + kVxmVyx1, bVxpkVxmVyx2 = bVx + kVxmVyx2;
			double Vzx1 = Vz*x1, Vzx2 = Vz*x2, Vxz1 = Vx*z1, Vxz2 = Vx*z2;
			double Vzx1mVxz1 = Vzx1-Vxz1, Vzx1mVxz2 = Vzx1-Vxz2, Vzx2mVxz1 = Vzx2-Vxz1, Vzx2mVxz2 = Vzx2-Vxz2;
			double kVzx1 = kVz*x1, kVzx2 = kVz*x2, Vyz1 = Vy*z1, Vyz2 = Vy*z2;
			double bVzpkVzx1mVyz1 = bVz + kVzx1 - Vyz1, bVzpkVzx1mVyz2 = bVz + kVzx1 - Vyz2, bVzpkVzx2mVyz1 = bVz + kVzx2 - Vyz1, bVzpkVzx2mVyz2 = bVz + kVzx2 - Vyz2;
			double bVzpkVxmVyz1 = bVz + kVxmVyz1, bVzpkVxmVyz2 = bVz + kVxmVyz2; 
			double kVxe2pVze2mVxVyz1 = kVxe2pVze2mVxVy*z1, kVxe2pVze2mVxVyz2 = kVxe2pVze2mVxVy*z2;
			double bVxpkVxmVyx1e2 = bVxpkVxmVyx1*bVxpkVxmVyx1, bVxpkVxmVyx2e2 = bVxpkVxmVyx2*bVxpkVxmVyx2;
			double BufForLogsX1 = VzdVxe2pVye2*(bVxmkVxe2pVze2mVxVyVy - kVxmVye2p1pke2Vze2Vy*x1);
			double BufForLogsX2 = VzdVxe2pVye2*(bVxmkVxe2pVze2mVxVyVy - kVxmVye2p1pke2Vze2Vy*x2);
			double bkVxmVyVxe2pVze2 = b*kVxmVy*Vxe2pVze2, bVxe2pVze2 = b*Vxe2pVze2;
			double ke2p1Vye2mkVypVxe2Vz = (ke2p1*Vye2 - kVypVx*kVypVx)*Vz;
			double kVxe2pVze2mVxVyx1 = kVxe2pVze2mVxVy*x1, kVxe2pVze2mVxVyx2 = kVxe2pVze2mVxVy*x2;
			double VyVzz1 = VyVz*z1, VyVzz2 = VyVz*z2;
			double kVxmVyVxmkVypVxVyVze2 = kVxmVy*Vx - kVypVx*Vy*Vze2;
			double BufForLogsZ1 = (bVxe2pVze2*kVypVxVz + kVxmVyVxmkVypVxVyVze2*z1)/kVxmVye2p1pke2Vze2;
			double BufForLogsZ2 = (bVxe2pVze2*kVypVxVz + kVxmVyVxmkVypVxVyVze2*z2)/kVxmVye2p1pke2Vze2;
			double VyVzx1 = VyVz*x1, VyVzx2 = VyVz*x2;

			double SumAtans1 = atan(TransAtans((-bVyVz - kVypVxVzx1 + Vxe2pVye2*z1)/bVxpkVxmVyx1, (bVyVz + kVypVxVzx1 - Vxe2pVye2*z2)/bVxpkVxmVyx1, PiMult1));
			SumAtans1 += PiMult1*PI;
			double SumAtans2 = atan(TransAtans((-bVyVz - kVypVxVzx2 + Vxe2pVye2*z2)/bVxpkVxmVyx2, (bVyVz + kVypVxVzx2 - Vxe2pVye2*z1)/bVxpkVxmVyx2, PiMult1));
			SumAtans2 += PiMult1*PI;
			double SumAtans3 = atan(TransAtans((bkVxe2pVze2mVxVy + kVxmVye2p1pke2Vze2x1 - kVypVxVzz1)/bVzpkVxmVyz1, -(bkVxe2pVze2mVxVy + kVxmVye2p1pke2Vze2x2 - kVypVxVzz1)/bVzpkVxmVyz1, PiMult1));
			SumAtans3 += PiMult1*PI;
			double SumAtans4 = atan(TransAtans((bkVxe2pVze2mVxVy + kVxmVye2p1pke2Vze2x2 - kVypVxVzz2)/bVzpkVxmVyz2, -(bkVxe2pVze2mVxVy + kVxmVye2p1pke2Vze2x1 - kVypVxVzz2)/bVzpkVxmVyz2, PiMult1));
			SumAtans4 += PiMult1*PI;
			double SumAtans5 = atan(TransAtans((bVxe2pVze2 + kVxe2pVze2mVxVyx1 - VyVzz1)/(Vxz1 - Vzx1), -(bVxe2pVze2 + kVxe2pVze2mVxVyx1 - VyVzz2)/(Vxz2 - Vzx1), PiMult1));
			SumAtans5 += PiMult1*PI;
			double SumAtans6 = atan(TransAtans((bVxe2pVze2 + kVxe2pVze2mVxVyx2 - VyVzz2)/(Vxz2 - Vzx2), -(bVxe2pVze2 + kVxe2pVze2mVxVyx2 - VyVzz1)/(Vxz1 - Vzx2), PiMult1));
			SumAtans6 += PiMult1*PI;

			double Log1 = log(bVxpkVxmVyx1e2 + bVzpkVzx1mVyz1*bVzpkVzx1mVyz1 + Vzx1mVxz1*Vzx1mVxz1);
			double Log2 = log(bVxpkVxmVyx1e2 + bVzpkVzx1mVyz2*bVzpkVzx1mVyz2 + Vzx1mVxz2*Vzx1mVxz2);
			double Log3 = log(bVxpkVxmVyx2e2 + bVzpkVzx2mVyz1*bVzpkVzx2mVyz1 + Vzx2mVxz1*Vzx2mVxz1);
			double Log4 = log(bVxpkVxmVyx2e2 + bVzpkVzx2mVyz2*bVzpkVzx2mVyz2 + Vzx2mVxz2*Vzx2mVxz2);

			S22 += (-p2dVxe2pVye2*((bkVxVxe2pVze2pVyVze2mVxe2 + kVxmVye2p1pke2Vze2Vx*x1)*SumAtans1 + (bkVxVxe2pVze2pVyVze2mVxe2 + kVxmVye2p1pke2Vze2Vx*x2)*SumAtans2) + m2Vzz1*SumAtans3 + m2Vzz2*SumAtans4
					- (BufForLogsX1 + kVxe2pVze2mVxVyz1)*Log1 + (BufForLogsX1 + kVxe2pVze2mVxVyz2)*Log2 + (BufForLogsX2 + kVxe2pVze2mVxVyz1)*Log3 - (BufForLogsX2 + kVxe2pVze2mVxVyz2)*Log4)/kVxmVye2p1pke2Vze2;

			S23 += (2.*(bVz*(SumAtans1 + SumAtans2) - kVxmVyz1*SumAtans3 - kVxmVyz2*SumAtans4)
					+ (kVypVxVzz1-bkVxe2pVze2mVxVy-kVxmVye2p1pke2Vze2x1)*Log1 - (kVypVxVzz2-bkVxe2pVze2mVxVy-kVxmVye2p1pke2Vze2x1)*Log2 - (kVypVxVzz1-bkVxe2pVze2mVxVy-kVxmVye2p1pke2Vze2x2)*Log3 + (kVypVxVzz2-bkVxe2pVze2mVxVy-kVxmVye2p1pke2Vze2x2)*Log4)/kVxmVye2p1pke2Vze2;

			S33 += (-2.*(((bkVxmVyVxe2pVze2 + ke2p1Vye2mkVypVxe2Vz*z1)*SumAtans3 + (bkVxmVyVxe2pVze2 + ke2p1Vye2mkVypVxe2Vz*z2)*SumAtans4)/kVxmVye2p1pke2Vze2 + Vx*(x1*SumAtans5 + x2*SumAtans6))
					+ (BufForLogsZ1+VyVzx1)*Log1 - (BufForLogsZ2+VyVzx1)*Log2 - (BufForLogsZ1+VyVzx2)*Log3 + (BufForLogsZ2+VyVzx2)*Log4)/Vxe2pVze2;
		}

		if(abs_y2my1*Max_k > abs_x2mx1)
		{
			double k1 = x2mx1/y2my1, b1 = x1 - k1*y1;

			double k1e2 = k1*k1; 
			double k1Vx = k1*Vx, k1Vy = k1*Vy, k1Vz = k1*Vz;
			double Vxmk1Vy = Vx - k1Vy, k1VxpVy = k1Vx + Vy;
			double k1VxpVyVz = k1VxpVy*Vz;
			double Vxmk1Vye2p1pk1e2Vze2 = Vxmk1Vy*Vxmk1Vy + (1.+k1e2)*Vze2;
			double b1Vx = b1*Vx, b1Vy = b1*Vy, b1Vz = b1*Vz, b1VxVz = b1*VxVz;
			double b1k1VyVye2pVze2pVxVze2mVye2 = b1*(k1Vy*Vye2pVze2 + VxVze2mVye2);
			double b1k1VxVye2pVze2mVxe2pVze2Vy = b1*(k1Vx*Vye2pVze2-Vxe2pVze2Vy);
			double b1k1Vye2pVze2mVxVy = b1*(k1*Vye2pVze2-VxVy);
			double Vxmk1Vye2p1pk1e2Vze2Vy = Vxmk1Vye2p1pk1e2Vze2*Vy;
			double mb12Vxe2Vymk1VxmVyVye2pVze2 = -b1*(p2Vxe2Vy - (k1Vx - Vy)*Vye2pVze2);
			double VxVxmk1Vye2p1pk1e2Vze2 = Vx*Vxmk1Vye2p1pk1e2Vze2;
			double Vxmk1VyVymk1Vze2 = Vxmk1Vy*Vy - k1*Vze2;
			double VxVxmk1VypVze2 = Vx*Vxmk1Vy + Vze2;
			double b1VymVxmk1Vyy1 = b1Vy - Vxmk1Vy*y1;
			double b1VymVxmk1Vyy2 = b1Vy - Vxmk1Vy*y2;
			double k1Vzy1 = k1Vz*y1, k1Vzy2 = k1Vz*y2;
			double b1Vzpk1Vzy1mVxz1 = b1Vz + k1Vzy1 - Vxz1;
			double b1Vzpk1Vzy1mVxz2 = b1Vz + k1Vzy1 - Vxz2;
			double b1Vzpk1Vzy2mVxz1 = b1Vz + k1Vzy2 - Vxz1;
			double b1Vzpk1Vzy2mVxz2 = b1Vz + k1Vzy2 - Vxz2;
			double Vzy1 = Vz*y1, Vzy2 = Vz*y2;
			double Vzy1mVyz1 = Vzy1 - Vyz1, Vzy1mVyz2 = Vzy1 - Vyz2, Vzy2mVyz1 = Vzy2 - Vyz1, Vzy2mVyz2 = Vzy2 - Vyz2;
			double k1VxpVyVzy1 = k1VxpVyVz*y1, k1VxpVyVzy2 = k1VxpVyVz*y2, Vxe2pVye2z1 = Vxe2pVye2*z1, Vxe2pVye2z2 = Vxe2pVye2*z2;
			double Vxmk1Vye2p1pk1e2Vze2y1 = Vxmk1Vye2p1pk1e2Vze2*y1, Vxmk1Vye2p1pk1e2Vze2y2 = Vxmk1Vye2p1pk1e2Vze2*y2;
			double k1VxpVyVzz1 = k1VxpVyVz*z1, k1VxpVyVzz2 = k1VxpVyVz*z2;
			double b1VzmVxmk1Vyz1 = b1Vz - Vxmk1Vy*z1, b1VzmVxmk1Vyz2 = b1Vz - Vxmk1Vy*z2;
			double Vxmk1Vye2p1pk1e2Vze2Vyy1 = Vxmk1Vye2p1pk1e2Vze2Vy*y1, Vxmk1Vye2p1pk1e2Vze2Vyy2 = Vxmk1Vye2p1pk1e2Vze2Vy*y2;
			double VxVxmk1Vye2p1pk1e2Vze2y1 = VxVxmk1Vye2p1pk1e2Vze2*y1, VxVxmk1Vye2p1pk1e2Vze2y2 = VxVxmk1Vye2p1pk1e2Vze2*y2;
			double Vxmk1VyVymk1Vze2z1 = Vxmk1VyVymk1Vze2*z1, Vxmk1VyVymk1Vze2z2 = Vxmk1VyVymk1Vze2*z2;
			double b1VymVxmk1Vyy1e2 = b1VymVxmk1Vyy1*b1VymVxmk1Vyy1, b1VymVxmk1Vyy2e2 = b1VymVxmk1Vyy2*b1VymVxmk1Vyy2;
			double BufForLogsY1 = (mb12Vxe2Vymk1VxmVyVye2pVze2 + VxVxmk1Vye2p1pk1e2Vze2y1)*VzdVxe2pVye2;
			double BufForLogsY2 = (mb12Vxe2Vymk1VxmVyVye2pVze2 + VxVxmk1Vye2p1pk1e2Vze2y2)*VzdVxe2pVye2;
			double p2k1Vz = 2.*k1Vz;
			double b1VxVxe2pVze2pk1VyVye2pVze2 = b1*(Vx*Vxe2pVze2 + k1Vy*Vye2pVze2);
			double VxVxmk1VypVze2z1 = VxVxmk1VypVze2*z1, VxVxmk1VypVze2z2 = VxVxmk1VypVze2*z2;
			double Buf1ForLogsY1 = (b1VxVxe2pVze2pk1VyVye2pVze2 + Vxmk1Vye2p1pk1e2Vze2Vyy1)*VzdVxe2pVye2;
			double Buf1ForLogsY2 = (b1VxVxe2pVze2pk1VyVye2pVze2 + Vxmk1Vye2p1pk1e2Vze2Vyy2)*VzdVxe2pVye2;

			double SumAtans1 = atan(TransAtans((b1VxVz + k1VxpVyVzy1 - Vxe2pVye2z1)/b1VymVxmk1Vyy1, -(b1VxVz + k1VxpVyVzy1 - Vxe2pVye2z2)/b1VymVxmk1Vyy1, PiMult1));
			SumAtans1 += PiMult1*PI;
			double SumAtans2 = atan(TransAtans((b1VxVz + k1VxpVyVzy2 - Vxe2pVye2z2)/b1VymVxmk1Vyy2, -(b1VxVz + k1VxpVyVzy2 - Vxe2pVye2z1)/b1VymVxmk1Vyy2, PiMult1));
			SumAtans2 += PiMult1*PI;
			double SumAtans3 = atan(TransAtans((b1k1Vye2pVze2mVxVy - k1VxpVyVzz1 + Vxmk1Vye2p1pk1e2Vze2y1)/b1VzmVxmk1Vyz1, -(b1k1Vye2pVze2mVxVy - k1VxpVyVzz1 + Vxmk1Vye2p1pk1e2Vze2y2)/b1VzmVxmk1Vyz1, PiMult1));
			SumAtans3 += PiMult1*PI;
			double SumAtans4 = atan(TransAtans((b1k1Vye2pVze2mVxVy - k1VxpVyVzz2 + Vxmk1Vye2p1pk1e2Vze2y2)/b1VzmVxmk1Vyz2, -(b1k1Vye2pVze2mVxVy - k1VxpVyVzz2 + Vxmk1Vye2p1pk1e2Vze2y1)/b1VzmVxmk1Vyz2, PiMult1));
			SumAtans4 += PiMult1*PI;

			double Log1 = log(b1VymVxmk1Vyy1e2 + b1Vzpk1Vzy1mVxz1*b1Vzpk1Vzy1mVxz1 + Vzy1mVyz1*Vzy1mVyz1);
			double Log2 = log(b1VymVxmk1Vyy1e2 + b1Vzpk1Vzy1mVxz2*b1Vzpk1Vzy1mVxz2 + Vzy1mVyz2*Vzy1mVyz2);
			double Log3 = log(b1VymVxmk1Vyy2e2 + b1Vzpk1Vzy2mVxz1*b1Vzpk1Vzy2mVxz1 + Vzy2mVyz1*Vzy2mVyz1);
			double Log4 = log(b1VymVxmk1Vyy2e2 + b1Vzpk1Vzy2mVxz2*b1Vzpk1Vzy2mVxz2 + Vzy2mVyz2*Vzy2mVyz2);

			S11 += (p2dVxe2pVye2*((b1k1VyVye2pVze2pVxVze2mVye2 + Vxmk1Vye2p1pk1e2Vze2Vyy1)*SumAtans1 + (b1k1VyVye2pVze2pVxVze2mVye2 + Vxmk1Vye2p1pk1e2Vze2Vyy2)*SumAtans2) + m2Vzz1*SumAtans3 + m2Vzz2*SumAtans4
					+ (Vxmk1VyVymk1Vze2z1 + BufForLogsY1)*Log1 - (Vxmk1VyVymk1Vze2z2 + BufForLogsY1)*Log2 - (Vxmk1VyVymk1Vze2z1 + BufForLogsY2)*Log3 + (Vxmk1VyVymk1Vze2z2 + BufForLogsY2)*Log4)/Vxmk1Vye2p1pk1e2Vze2;

			S12 += (-p2dVxe2pVye2*((b1k1VxVye2pVze2mVxe2pVze2Vy + VxVxmk1Vye2p1pk1e2Vze2y1)*SumAtans1 + (b1k1VxVye2pVze2mVxe2pVze2Vy + VxVxmk1Vye2p1pk1e2Vze2y2)*SumAtans2) + p2k1Vz*z1*SumAtans3 + p2k1Vz*z2*SumAtans4
					+ (Buf1ForLogsY1 - VxVxmk1VypVze2z1)*Log1 - (Buf1ForLogsY1 - VxVxmk1VypVze2z2)*Log2 - (Buf1ForLogsY2 - VxVxmk1VypVze2z1)*Log3 + (Buf1ForLogsY2 - VxVxmk1VypVze2z2)*Log4)/Vxmk1Vye2p1pk1e2Vze2;

			S13 += (2.*(-b1Vz*(SumAtans1 + SumAtans2) + Vxmk1Vy*(z1*SumAtans3 + z2*SumAtans4))
					- (b1k1Vye2pVze2mVxVy - k1VxpVyVzz1 + Vxmk1Vye2p1pk1e2Vze2y1)*Log1 + (b1k1Vye2pVze2mVxVy - k1VxpVyVzz2 + Vxmk1Vye2p1pk1e2Vze2y1)*Log2 + (b1k1Vye2pVze2mVxVy - k1VxpVyVzz1 + Vxmk1Vye2p1pk1e2Vze2y2)*Log3 - (b1k1Vye2pVze2mVxVy - k1VxpVyVzz2 + Vxmk1Vye2p1pk1e2Vze2y2)*Log4)/Vxmk1Vye2p1pk1e2Vze2;
		}
		x1 = x2; y1 = y2;
	}
//FinalDefinitionOfFieldIntegrals:

	double mx = Magn.y;
	double my = Magn.z;
	double mz = Magn.x;

	TVector3d BufIh;
	double& Ihx = BufIh.y;
	double& Ihy = BufIh.z;
	double& Ihz = BufIh.x;
	Ihx = ConstForH*(S11*mx+S12*my+S13*mz);
	Ihy = ConstForH*(S12*mx-S22*my-S23*mz);
	Ihz = ConstForH*(S13*mx-S23*my-S33*mz);

	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIh;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIh;
}

//-------------------------------------------------------------------------

void radTExtrPolygon::B_intCompSpecCases(radTField* FieldPtr, const TSpecCaseID& SpecCaseID)
{
	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;

	const double Max_k = 1.E+08;

	TVector3d& StPo = FieldPtr->P;

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(CentrPoint.z);

	double z1 = FirstPoint.x - StPo.x;
	double z2 = z1 + Thickness;

// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(z1==0.) z1 = AbsRandX;
	if(z2==0.) z2 = AbsRandX;

	double z1e2 = z1*z1, z2e2 = z2*z2;
	double z2mz1 = z2 - z1;

	radTPolygon* BasePolygonPtr = (radTPolygon*)(BasePolygonHandle.rep);

#ifdef __GCC__
	vector<TVector2d>::iterator BaseIter = (BasePolygonPtr->EdgePointsVector).begin();
#else
	vector<TVector2d, allocator<TVector2d> >::iterator BaseIter = (BasePolygonPtr->EdgePointsVector).begin();
#endif

	int AmOfEdPoInBase = BasePolygonPtr->AmOfEdgePoints;
	int AmOfEdPoInBase_mi_1 = AmOfEdPoInBase - 1;
	
	TVector2d First2d(FirstPoint.y - StPo.y, FirstPoint.z - StPo.z);
	TVector2d Vect2dToAdd(First2d.x - (*BaseIter).x, First2d.y - (*BaseIter).y);

	double x1 = First2d.x;
	double y1 = First2d.y;

// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(x1==0.) x1 = AbsRandY;
	if(y1==0.) y1 = AbsRandZ;

	double x2, y2;
	double x1e2 = x1*x1, x2e2, y1e2 = y1*y1, y2e2;

	double S11=0., S12=0., S13=0., S22=0., S23=0., S33=0.;
	double PiMult1=0.;

	for(int i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x + Vect2dToAdd.x;
			y2 = (*BaseIter).y + Vect2dToAdd.y;
		}
		else
		{
			x2 = First2d.x;
			y2 = First2d.y;
		}

		// Artificial shift of an observation point a bit right of the block's border
		// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(x2==0.) x2 = AbsRandY;
		if(y2==0.) y2 = AbsRandZ;

		x2e2 = x2*x2; y2e2 = y2*y2;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			if(SpecCaseID==ZeroVxVz)
			{
				double k = y2my1/x2mx1, b = y1 - k*x1;
				
				double SumAtans1 = atan(TransAtans(x1/z1, -x2/z1, PiMult1));
				SumAtans1 += PiMult1*PI;
				double SumAtans2 = atan(TransAtans(x2/z2, -x1/z2, PiMult1));
				SumAtans1 += PiMult1*PI;

				double SumLogs1 = log((x1e2 + z1e2)/(x2e2 + z1e2));
				double SumLogs2 = log((x2e2 + z2e2)/(x1e2 + z2e2));

				double Buf1 = 2.*b*(SumAtans1 + SumAtans2) + k*(z1*SumLogs1 + z2*SumLogs2);

				S11 -= Buf1;
				S13 -= 2.*k*(z1*SumAtans1 + z2*SumAtans2) - b*(SumLogs1 + SumLogs2);
				S33 -= Buf1;
			}
		}
		if(abs_y2my1*Max_k > abs_x2mx1)
		{
			if(SpecCaseID==ZeroVxVy)
			{
				double k1 = x2mx1/y2my1, b1 = x1 - k1*y1;

				double k1e2 = k1*k1;
				double k1e2p1 = k1e2+1;
				double b1k1 = b1*k1;

				double b1db1k1pk1e2p1y1 = b1/(k1*x1+y1);   /* b1/(b1k1 + k1e2p1*y1); */
				double b1db1k1pk1e2p1y2 = b1/(k1*x2+y2);   /* b1/(b1k1 + k1e2p1*y2); */
				double b1pk1y1 = x1;   /* b1+k1*y1; */
				double b1pk1y2 = x2;   /* b1+k1*y2; */

				double SumAtans1 = atan(TransAtans(b1db1k1pk1e2p1y1, -b1db1k1pk1e2p1y2, PiMult1));
				double yFlp = -b1k1/k1e2p1;
				double FlpRep = Sign(b1)*(Step(yFlp-y1)*Step(y2-yFlp)*Step(y2my1) - Step(yFlp-y2)*Step(y1-yFlp)*Step(-y2my1));
				SumAtans1 += (PiMult1 + FlpRep)*PI;
				double SumLogs1 = log((y1e2 + b1pk1y1*b1pk1y1)/(y2e2 + b1pk1y2*b1pk1y2));
				double Buf1 = -z2mz1*(2.*SumAtans1 - k1*SumLogs1)/k1e2p1;

				S11 += Buf1;
				S12 += z2mz1*(2.*k1*SumAtans1 + SumLogs1)/k1e2p1;
				S22 += Buf1;
			}
			else if(SpecCaseID==ZeroVyVz)
			{
				double k1 = x2mx1/y2my1, b1 = x1 - k1*y1;

				double SumAtans1 = atan(TransAtans(y1/z1, -y2/z1, PiMult1));
				SumAtans1 += PiMult1*PI;
				double SumAtans2 = atan(TransAtans(y2/z2, -y1/z2, PiMult1));
				SumAtans1 += PiMult1*PI;

				double SumLogs1 = log((y1e2 + z1e2)/(y2e2 + z1e2));
				double SumLogs2 = log((y2e2 + z2e2)/(y1e2 + z2e2));

				double Buf1 = 2.*b1*(SumAtans1 + SumAtans2) + k1*(z1*SumLogs1 + z2*SumLogs2);

				S22 -= Buf1;
				S23 -= 2.*k1*(z1*SumAtans1 + z2*SumAtans2) - b1*(SumLogs1 + SumLogs2);
				S33 += Buf1;
			}
		}
		x1 = x2; y1 = y2;
		x1e2 = x2e2; y1e2 = y2e2;
	}
//FinalDefinitionOfFieldIntegrals:

	double mx = Magn.y;
	double my = Magn.z;
	double mz = Magn.x;

	TVector3d BufIh;
	double& Ihx = BufIh.y;
	double& Ihy = BufIh.z;
	double& Ihz = BufIh.x;
	Ihx = ConstForH*(S11*mx+S12*my+S13*mz);
	Ihy = ConstForH*(S12*mx-S22*my-S23*mz);
	Ihz = ConstForH*(S13*mx-S23*my-S33*mz);

	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIh;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIh;
}

//-------------------------------------------------------------------------

int radTExtrPolygon::SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	double LocSubdivArray[15];
	for(int jj=0; jj<15; jj++) LocSubdivArray[jj] = SubdivArray[jj];

	const double ZeroTol = 1.E-10;
	if((pSubdivOptions->SubdivisionParamCode == 0) && (fabs(LocSubdivArray[0]-1.)<ZeroTol) && (fabs(LocSubdivArray[2]-1.)<ZeroTol) && (fabs(LocSubdivArray[4]-1.)<ZeroTol)) return 1;

	radTSend Send;
	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty())) 
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

	double &kx = LocSubdivArray[0];
	double &qx = LocSubdivArray[1];
	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		kx = (kx < Thickness)? Round(Thickness/kx) : 1.;
	}

	radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedExtrPolygon(this, LocSubdivArray);
	radThg NewHandle(GroupInPlaceOfThisPtr);

	TVector2d OldBaseFirstEdgePoi = ((radTPolygon*)(BasePolygonHandle.rep))->EdgePointsVector[0];

	const double AbsZeroTol = 5.E-12;
	double q0x = (fabs(kx-1.)>AbsZeroTol)? pow(qx, 1./(kx-1.)) : qx;
	double BufX = qx*q0x - 1.;
	double a1x = (fabs(BufX) > AbsZeroTol)? Thickness*(q0x - 1.)/BufX : Thickness/kx;

	double SubdArrayForBase[] = { LocSubdivArray[2], LocSubdivArray[3], LocSubdivArray[4], LocSubdivArray[5]};

	radThg GrInPlaceOfBasePgnHandle = BasePolygonHandle;
// Maybe duplicate the base first ?
	int SubdOK = ((radTPolygon*)(BasePolygonHandle.rep))->SubdivideItself(SubdArrayForBase, GrInPlaceOfBasePgnHandle, radPtr, pSubdivOptions);
	if(!SubdOK) return 0;

	short BaseReallySubdivided = (GrInPlaceOfBasePgnHandle.rep != BasePolygonHandle.rep)? 1 : 0;

	int NewStuffCounter = 0;
	if(BaseReallySubdivided)
	{
		radTGroup* GrInPlaceOfBasePgnPtr = (radTGroup*)(GrInPlaceOfBasePgnHandle.rep);
		for(radTmhg::const_iterator iter = GrInPlaceOfBasePgnPtr->GroupMapOfHandlers.begin();
			iter != GrInPlaceOfBasePgnPtr->GroupMapOfHandlers.end(); ++iter)
		{
			radThg NewBasePgnHandle = (*iter).second;
			radTPolygon* NewBasePgnPtr = (radTPolygon*)(NewBasePgnHandle.rep);

			TVector2d AddForFirstPoi = NewBasePgnPtr->EdgePointsVector[0] - OldBaseFirstEdgePoi;
			TVector3d NewFirstPoint(FirstPoint.x, FirstPoint.y + AddForFirstPoi.x, FirstPoint.z + AddForFirstPoi.y);

			double NewThickness = a1x;

			for(int ix=0; ix<int(kx); ix++)
			{
				if(!NewBasePgnPtr->DuplicateItself(NewBasePgnHandle, radPtr, PutNewStuffIntoGenCont)) return 0;

				radThg hg(new radTExtrPolygon(NewFirstPoint, ParallelToX, NewThickness, NewBasePgnHandle, Magn, MaterHandle));
				if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(hg), hg);
				else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, hg);

				NewFirstPoint.x += NewThickness;
				NewThickness *= q0x;
			}
			if(PutNewStuffIntoGenCont)
			{
				short OldSendingIsRequired = radPtr->SendingIsRequired;
				radPtr->SendingIsRequired = 0;
				radPtr->DeleteElement(radPtr->RetrieveElemKey(NewBasePgnPtr));
				radPtr->SendingIsRequired = OldSendingIsRequired;
			}
		}
	}
	else
	{
		radThg NewBasePgnHandle = BasePolygonHandle;
		radTPolygon* NewBasePgnPtr = (radTPolygon*)(NewBasePgnHandle.rep);
		TVector2d AddForFirstPoi = NewBasePgnPtr->EdgePointsVector[0] - OldBaseFirstEdgePoi;

		TVector3d NewFirstPoint(FirstPoint.x, FirstPoint.y + AddForFirstPoi.x, FirstPoint.z + AddForFirstPoi.y);
		double NewThickness = a1x;

		for(int ix=0; ix<int(kx); ix++)
		{
			if(!NewBasePgnPtr->DuplicateItself(NewBasePgnHandle, radPtr, PutNewStuffIntoGenCont)) return 0;

			radThg hg(new radTExtrPolygon(NewFirstPoint, ParallelToX, NewThickness, NewBasePgnHandle, Magn, MaterHandle));
			if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(hg), hg);
			else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, hg);

			NewFirstPoint.x += NewThickness;
			NewThickness *= q0x;
		}
	}
	((radTSubdividedExtrPolygon*)GroupInPlaceOfThisPtr)->AmOfSubElem = (int)(GroupInPlaceOfThisPtr->GroupMapOfHandlers.size());

	In_hg = NewHandle;
	return 1;
}

//-------------------------------------------------------------------------

int radTExtrPolygon::ConvertToPolyhedron(radThg& In_hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{// This does not handle internal faces after cut,
 // so be careful applying this to Extr. Polygons obtained after subdivision !
	radTSend Send;

	radTPolygon* BasePgnPtr = (radTPolygon*)(BasePolygonHandle.rep);
	int AmOfEdgePoInBase = BasePgnPtr->AmOfEdgePoints;
	radTVect2dVect& BasePgnEdgePointsVector = BasePgnPtr->EdgePointsVector;

	int AmOfVertexPoints = 2*AmOfEdgePoInBase;

	TVector3d* ArrayOfPoints = new TVector3d[AmOfVertexPoints];
	if(ArrayOfPoints == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
	TVector3d* ArrayOfPointsTravers1 = ArrayOfPoints;
	TVector3d* ArrayOfPointsTravers2 = &(ArrayOfPoints[AmOfEdgePoInBase]);

	int AmOfFaces = AmOfEdgePoInBase+2;
	int** ArrayOfFaces = new int*[AmOfFaces];
	if(ArrayOfFaces == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

	int* FacesIndexes = new int[AmOfEdgePoInBase*7+2];

	int* BaseFace1 = &(FacesIndexes[0]);
	int* BaseFace2 = &(FacesIndexes[AmOfEdgePoInBase]);

	int** ArrayOfFacesTravers = ArrayOfFaces;
	*(ArrayOfFacesTravers++) = BaseFace1;
	*(ArrayOfFacesTravers++) = BaseFace2;

	int* ArrayOfLengths = &(FacesIndexes[AmOfEdgePoInBase*6]);
	int* ArrayOfLengthsTravers = ArrayOfLengths;
	*(ArrayOfLengthsTravers++) = AmOfEdgePoInBase;
	*(ArrayOfLengthsTravers++) = AmOfEdgePoInBase;

	TVector2d OldBaseFirstEdgePoi = ((radTPolygon*)(BasePolygonHandle.rep))->EdgePointsVector[0]; //OC160406

	int AmOfEdgePoInBase_mi_1 = AmOfEdgePoInBase - 1;
	for(int k=0; k<AmOfEdgePoInBase; k++)
	{
		int kp1 = k+1, kp2 = k+2;
		*(BaseFace1++) = kp1;
		*(BaseFace2++) = AmOfEdgePoInBase+kp1;
		int* CurrentMantleFace = &(FacesIndexes[AmOfVertexPoints + k*4]);
		*(ArrayOfFacesTravers++) = CurrentMantleFace;
		*(CurrentMantleFace++) = kp1;
		if(k==AmOfEdgePoInBase_mi_1)
		{
			*(CurrentMantleFace++) = 1;
			*(CurrentMantleFace++) = AmOfEdgePoInBase+1;
		}
		else
		{
			*(CurrentMantleFace++) = kp2;
			*(CurrentMantleFace++) = AmOfEdgePoInBase+kp2;
		}
		*CurrentMantleFace = AmOfEdgePoInBase+kp1;
		*(ArrayOfLengthsTravers++) = 4;

		TVector2d& CurBasePo = BasePgnEdgePointsVector[k];
		ArrayOfPointsTravers1->x = FirstPoint.x; ArrayOfPointsTravers2->x = FirstPoint.x + Thickness; 

		//ArrayOfPointsTravers1->y = ArrayOfPointsTravers2->y = CurBasePo.x;
		//ArrayOfPointsTravers1->z = ArrayOfPointsTravers2->z = CurBasePo.y;
		ArrayOfPointsTravers1->y = ArrayOfPointsTravers2->y = FirstPoint.y + (CurBasePo.x - OldBaseFirstEdgePoi.x); //OC160406
		ArrayOfPointsTravers1->z = ArrayOfPointsTravers2->z = FirstPoint.z + (CurBasePo.y - OldBaseFirstEdgePoi.y);

		ArrayOfPointsTravers1++; ArrayOfPointsTravers2++; 
	}

	if(ConsiderOnlyWithTrans)
	{
		radTrans ResTransf;
		short SomethingFound = 0;
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
		if(SomethingFound) 
		{
			for(int i=0; i<AmOfVertexPoints; i++)
			{
				ArrayOfPoints[i] = ResTransf.TrPoint(ArrayOfPoints[i]);
			}
		}
	}

	radTPolyhedron* PolyhedronPtr = new radTPolyhedron(ArrayOfPoints, AmOfVertexPoints, ArrayOfFaces, ArrayOfLengths, AmOfFaces, Magn);
	if(PolyhedronPtr == 0) return 0;
	PolyhedronPtr->MaterHandle = MaterHandle;
	PolyhedronPtr->IsGroupMember = IsGroupMember;

	PolyhedronPtr->g3dListOfTransform = g3dListOfTransform;
	if(ConsiderOnlyWithTrans) PolyhedronPtr->EraseInnerTransform();
	PolyhedronPtr->ConsiderOnlyWithTrans = 0;

	PolyhedronPtr->CentrPoint = CentrPoint; // For full compatibility

	PolyhedronPtr->HandleAuxCompData = HandleAuxCompData;
	PolyhedronPtr->MessageChar = MessageChar;

	radThg NewHandle(PolyhedronPtr);
	In_hg = NewHandle;

	if(ArrayOfPoints != 0) delete[] ArrayOfPoints;
	if(ArrayOfFaces != 0) delete[] ArrayOfFaces;
	if(FacesIndexes != 0) delete[] FacesIndexes;
	return 1;
}

//-------------------------------------------------------------------------

int radTExtrPolygon::FindLowestAndUppestVertices(TVector3d& PlanesNormal, radTSubdivOptions* pSubdivOptions, 
	TVector3d& LowestVertexPoint, TVector3d& UppestVertexPoint, radTrans& Trans, char& TransWasSet, char& Ignore)
{
	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;

	Ignore = 0;
	TransWasSet = 0;
	TVector3d& ActualPlanesNormal = PlanesNormal;
	if(!SubdInLocFrame)
	{
		radTrans ResTransf;
		short SomethingFound = 0;
		FindResTransfWithMultOne(ResTransf, SomethingFound);
		if(SomethingFound) 
		{
			ActualPlanesNormal = ResTransf.TrBiPoint_inv(ActualPlanesNormal);
			Trans = ResTransf;
			TransWasSet = 1;
		}
	}
	else if(ConsiderOnlyWithTrans)
	{
		radTrans ResTransf;
		short SomethingFound = 0;
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
		if(SomethingFound) 
		{
			ActualPlanesNormal = ResTransf.TrBiPoint_inv(ActualPlanesNormal);
			Trans = ResTransf;
			TransWasSet = 1;
		}
	}

	TVector3d LowestPo = (1.E+23)*PlanesNormal, UppestPo = (-1.E+23)*PlanesNormal;

	radTPolygon* BasePgnPtr = (radTPolygon*)(BasePolygonHandle.rep);
	int AmOfEdgePoInBase = BasePgnPtr->AmOfEdgePoints;
	radTVect2dVect& BasePgnEdgePointsVector = BasePgnPtr->EdgePointsVector;

	TVector3d TestLoV, TestUpV;
	int AmOfEdgePoInBase_mi_1 = AmOfEdgePoInBase - 1;
	for(int k=0; k<AmOfEdgePoInBase; k++)
	{
		TVector2d& CurBasePo = BasePgnEdgePointsVector[k];
		TVector3d TestPo1(FirstPoint.x, CurBasePo.x, CurBasePo.y);
		TVector3d TestPo2(FirstPoint.x + Thickness, CurBasePo.x, CurBasePo.y);

		TestLoV = TestPo1 - LowestPo;
		TestUpV = TestPo1 - UppestPo;
		if(TestLoV*ActualPlanesNormal < 0.) LowestPo = TestPo1;
		if(TestUpV*ActualPlanesNormal > 0.) UppestPo = TestPo1;

		TestLoV = TestPo2 - LowestPo;
		TestUpV = TestPo2 - UppestPo;
		if(TestLoV*ActualPlanesNormal < 0.) LowestPo = TestPo2;
		if(TestUpV*ActualPlanesNormal > 0.) UppestPo = TestPo2;
	}
	LowestVertexPoint = LowestPo;
	UppestVertexPoint = UppestPo;
	return 1;
}

//-------------------------------------------------------------------------

void radTExtrPolygon::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique)
{
	radTPolygon* BasePgnPtr = (radTPolygon*)(BasePolygonHandle.rep);
	int AmOfEdgePoInBase = BasePgnPtr->AmOfEdgePoints;

	radTVect2dVect& BasePgnEdgePointsVector = BasePgnPtr->EdgePointsVector;

	int AmOfEdgePoInBase_mi_1 = AmOfEdgePoInBase - 1;
	for(int k=0; k<AmOfEdgePoInBase; k++)
	{
		TVector2d& CurBasePo = BasePgnEdgePointsVector[k];
		TVector3d P1(FirstPoint.x, CurBasePo.x, CurBasePo.y);
		TVector3d P2(FirstPoint.x + Thickness, CurBasePo.x, CurBasePo.y);

		OutVect.push_back(P1);
		OutVect.push_back(P2);
	}
}

//-------------------------------------------------------------------------

