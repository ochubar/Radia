/*-------------------------------------------------------------------------
*
* File name:      radflm.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source: filament conductor
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
#include "radflm.h"
#include "radg3dgr.h"

#include <math.h>
#include <sstream>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTFlmLinCur::radTFlmLinCur(const TVector3d& InStartPoint, const TVector3d& InEndPoint, double InI)
{
	I = InI; StartPoint = InStartPoint; EndPoint = InEndPoint;

	TVector3d LinVect = EndPoint - StartPoint;
	double SqLength = LinVect.x*LinVect.x + LinVect.y*LinVect.y + LinVect.z*LinVect.z;
	Length = sqrt(SqLength);

	if(LinVect.y==0. && LinVect.z==0.)
	{
		TVector3d St0(1.,0.,0.);
		TVector3d St1(0.,1.,0.);
		TVector3d St2(0.,0.,1.);
		TMatrix3d M(St0, St1, St2);
		TVector3d ZeroVect(0.,0.,0.);
		NativeRotation = radTrans(M, M, ZeroVect, 1., 1.); // Identity
		return;
	}

	TVector3d LinVectProto(Length, 0., 0.);
	TVector3d RotAx(LinVectProto.y*LinVect.z-LinVectProto.z*LinVect.y,
					LinVectProto.z*LinVect.x-LinVectProto.x*LinVect.z,
					LinVectProto.x*LinVect.y-LinVectProto.y*LinVect.x);
	double cosAngle = (LinVectProto * LinVect)/SqLength;
	double Angle = acos(cosAngle);
	SetNativeRotation(RotAx, Angle);

	TVector3d TestVect = NativeRotation.TrBiPoint(LinVectProto);
	const double SmallPositive = 1.E-10;
	const double TwoPi = 2.*3.141592653589793238;
	if((Abs(TestVect.x-LinVect.x)/Length > SmallPositive) ||
	   (Abs(TestVect.y-LinVect.y)/Length > SmallPositive) ||
	   (Abs(TestVect.z-LinVect.z)/Length > SmallPositive)) SetNativeRotation(RotAx, TwoPi-Angle);
}

//-------------------------------------------------------------------------

void radTFlmLinCur::SetNativeRotation(const TVector3d& InAxVect, double Angle)
{
	double NormFact = 1./sqrt(InAxVect.x*InAxVect.x+InAxVect.y*InAxVect.y+InAxVect.z*InAxVect.z);
	TVector3d AxVect = NormFact*InAxVect;
	double VxVx, VyVy, VzVz;
	VxVx=AxVect.x*AxVect.x; VyVy=AxVect.y*AxVect.y; VzVz=AxVect.z*AxVect.z;

	double cosAng, sinAng, One_m_cos;
	cosAng = cos(Angle); sinAng = sin(Angle); One_m_cos = 1. - cosAng;
	double One_m_cosVxVy, One_m_cosVxVz, One_m_cosVyVz, sinVx, sinVy, sinVz;
	One_m_cosVxVy = One_m_cos*AxVect.x*AxVect.y;
	One_m_cosVxVz = One_m_cos*AxVect.x*AxVect.z;
	One_m_cosVyVz = One_m_cos*AxVect.y*AxVect.z;
	sinVx = sinAng*AxVect.x; sinVy = sinAng*AxVect.y; sinVz = sinAng*AxVect.z;

	TVector3d St0(VxVx+cosAng*(VyVy+VzVz), One_m_cosVxVy-sinVz, One_m_cosVxVz+sinVy);
	TVector3d St1(One_m_cosVxVy+sinVz, VyVy+cosAng*(VxVx+VzVz), One_m_cosVyVz-sinVx);
	TVector3d St2(One_m_cosVxVz-sinVy, One_m_cosVyVz+sinVx, VzVz+cosAng*(VxVx+VyVy));
	TMatrix3d M(St0, St1, St2);
	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
	TMatrix3d M0(St00, St01, St02);

	NativeRotation = radTrans(M, M0*StartPoint, 1., 1.);
}

//-------------------------------------------------------------------------

void radTFlmLinCur::B_comp(radTField* FieldPtr)
{
	TVector3d BufP = NativeRotation.TrPoint_inv(FieldPtr->P);
	TVector3d V0 = StartPoint - BufP;
	double x1 = V0.x + Length;

	const double SmallPositive = 1.E-23;
	double y0y0_p_z0z0 = V0.y*V0.y + V0.z*V0.z;
	if(y0y0_p_z0z0 == 0.) y0y0_p_z0z0 = SmallPositive;

	double SqRt0 = sqrt(V0.x*V0.x + y0y0_p_z0z0);
	double SqRt1 = sqrt(x1*x1 + y0y0_p_z0z0);

	double ConI = 1.E-04 * I;
	double ComMult;

	if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
	{
		ComMult = ConI*(V0.x/SqRt0 - x1/SqRt1)/y0y0_p_z0z0;
		TVector3d BufB(0., -ComMult*V0.z, ComMult*V0.y);
		BufB = NativeRotation.TrVectField(BufB);
		if(FieldPtr->FieldKey.B_) FieldPtr->B += BufB;
		if(FieldPtr->FieldKey.H_) FieldPtr->H += BufB;
	}
	if(FieldPtr->FieldKey.A_)
	{
		double V0x_p_SqRt0 = V0.x+SqRt0;
		if(V0x_p_SqRt0==0.) V0x_p_SqRt0 = SmallPositive;
		ComMult = ConI*log((x1+SqRt1)/V0x_p_SqRt0);
		TVector3d BufA(ComMult, 0., 0.);
		FieldPtr->A += NativeRotation.TrVectPoten(BufA);
	}
}

//-------------------------------------------------------------------------

void radTFlmLinCur::B_intComp(radTField* FieldPtr)
{
	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

// An analytical algorithm for infinite Field Integrals:
	TVector3d BufP = NativeRotation.TrPoint_inv(FieldPtr->P);
	TVector3d BufNextP = NativeRotation.TrPoint_inv(FieldPtr->NextP);

	TVector3d VV = StartPoint - BufP;
	double x2 = VV.x + Length;

	TVector3d v = BufNextP - BufP;
	double Mod_v = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	v.x /= Mod_v; v.y /= Mod_v; v.z /= Mod_v;

	double ConI = 2.E-04 * I;
	TVector3d BufIntB(0.,0.,0.);
	double ComMult;

	const double SpecCaseZeroToler = 1.E-12;
	double vyvy_p_vzvz = v.y*v.y + v.z*v.z;
	if(vyvy_p_vzvz < SpecCaseZeroToler)
	{
		ComMult = ConI*(x2-VV.x)/(VV.z*VV.z + VV.y*VV.y);
		BufIntB.y = ComMult*VV.z;
		BufIntB.z = -ComMult*VV.y;
		goto FinalDefinitionOfFieldIntegrals;
	}
	{
		const double Pi = 3.141592653589793238;
		const double SmallestZeroToler = 1.E-14;
		double vzY_mi_vyZ = v.z*VV.y - v.y*VV.z; if(vzY_mi_vyZ==0.) vzY_mi_vyZ = SmallestZeroToler;
		double vxZ = v.x*VV.z;
		double vzX1_mi_vxZ = v.z*VV.x - vxZ;
		double vzX2_mi_vxZ = v.z*x2 - vxZ;
		double vxY = v.x*VV.y;
		double vxY_mi_vyX1 = vxY - v.y*VV.x;
		double vxY_mi_vyX2 = vxY - v.y*x2;
		double vxvyY_p_vxvzZ = v.x*(v.y*VV.y+v.z*VV.z);
		double vzYmivyZ_mu_vzYmivyZ = vzY_mi_vyZ*vzY_mi_vyZ;

		double PiMult = 0.;
		double F = (atan(TransAtans((vyvy_p_vzvz*VV.x-vxvyY_p_vxvzZ)/vzY_mi_vyZ, -(vyvy_p_vzvz*x2-vxvyY_p_vxvzZ)/vzY_mi_vyZ, PiMult)) + Pi*PiMult)/vyvy_p_vzvz;
		double G = 0.5*v.x*log((vzYmivyZ_mu_vzYmivyZ + vzX1_mi_vxZ*vzX1_mi_vxZ + vxY_mi_vyX1*vxY_mi_vyX1)/(vzYmivyZ_mu_vzYmivyZ + vzX2_mi_vxZ*vzX2_mi_vxZ + vxY_mi_vyX2*vxY_mi_vyX2))/vyvy_p_vzvz;

		BufIntB.y = ConI*(v.y*F+v.z*G);
		BufIntB.z = ConI*(v.z*F-v.y*G);
	}
FinalDefinitionOfFieldIntegrals:
	BufIntB = NativeRotation.TrVectField(BufIntB);
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIntB;
	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIntB;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTFlmLinCur::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTFlmLinCurGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

void radTFlmLinCur::Dump(std::ostream& o, int ShortSign) // Porting
{
	radTg3d::Dump(o);
	o << "FlmLinCur";
	if(ShortSign==1) return;
	o << endl;

	o << "   {x1,y1,z1}= {" << StartPoint.x << ','
							<< StartPoint.y << ','
							<< StartPoint.z << "}" << endl
	  << "   {x2,y2,z2}= {" << EndPoint.x << ','
							<< EndPoint.y << ','
							<< EndPoint.z << "}" << endl
	  << "   i= " << I;

	DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTFlmLinCur::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
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

	//radTrans NativeRotation;
	NativeRotation.DumpBin_Trans(oStr);

	//double Length;
	oStr << Length;

	//double I;
	oStr << I;

	//TVector3d StartPoint, EndPoint;
	oStr << StartPoint << EndPoint;
}

//-------------------------------------------------------------------------

int radTFlmLinCur::SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	if(!SubdivideCoils) return 1;
	radTSend Send;
	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty())) 
	{
		Send.ErrorMessage("Radia::Error108"); return 0;
	}

	double k = SubdivArray[0], q = SubdivArray[1];

	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		k = (k < Length)? Round(Length/k) : 1.;
	}

	const double ZeroTol = 1.E-10;
	if(fabs(k-1.)<ZeroTol) return 1;

	radTGroup* GroupInPlaceOfThisPtr = new radTGroup();
	IsGroupMember = GroupInPlaceOfThisPtr->IsGroupMember;
	g3dListOfTransform = GroupInPlaceOfThisPtr->g3dListOfTransform;

	radThg NewHandle(GroupInPlaceOfThisPtr);

	const double AbsZeroTol = 5.E-12;
	double q0 = (fabs(k-1.)>AbsZeroTol)? pow(q, 1./(k-1.)) : q;
	double Buf = q*q0 - 1.;
	double a1 = (fabs(Buf) > AbsZeroTol)? Length*(q0 - 1.)/Buf : Length/k;

	int kInt = int(k), k_mi_1 = kInt-1;

	TVector3d UnitVect = (1./Length)*(EndPoint - StartPoint);
	TVector3d NewStartPoint = StartPoint, NewEndPoint = StartPoint + a1*UnitVect;
	double NewLength = a1;

	int NewStuffCounter = 0;
	for(int i=0; i<kInt; i++)
	{
		radTFlmLinCur* FlmLinCurPtr = new radTFlmLinCur(NewStartPoint, NewEndPoint, I);
		if(FlmLinCurPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(FlmLinCurPtr);
		if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(hg), hg);
		else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, hg);

		NewLength *= q0;
		NewStartPoint = NewEndPoint;
		NewEndPoint = NewStartPoint + NewLength*UnitVect;
	}
	In_hg = NewHandle;
	return 1;
}

//-------------------------------------------------------------------------
