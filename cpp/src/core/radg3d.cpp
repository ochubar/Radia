/*-------------------------------------------------------------------------
*
* File name:      radg3d.cpp
*
* Project:        RADIA
*
* Description:    Base class for 3D objects - magnetic field sources
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
#include "radg3d.h"
#include "radg3da1.h"
//#include "radtrans.h"
#include "gmtrans.h"
#include "radcast.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

extern radTApplication rad;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTg3d::NestedFor_B(radTField* FieldPtr, const radTlphg::iterator& Iter)
{
	radTrans* TransPtr = (radTrans*)(((*Iter).Handler_g).rep);
	radTlphg::iterator LocalNextIter = Iter;
	LocalNextIter++;

	TVector3d ZeroVect(0.,0.,0.);

	short FldIntNeeded = FieldPtr->FieldKey.Ib_ || FieldPtr->FieldKey.Ih_; // Plus this string

	if((*Iter).m == 1)
	{
		radTField BufField = TransPtr->TrField_inv(*FieldPtr);
		BufField.P = TransPtr->TrPoint_inv(FieldPtr->P);
		if(FldIntNeeded) BufField.NextP = TransPtr->TrPoint_inv(FieldPtr->NextP); // Plus this string

		B_comp_Or_NestedFor(&BufField, LocalNextIter);

		BufField.P = FieldPtr->P;
		if(FldIntNeeded) BufField.NextP = FieldPtr->NextP; // Plus this string

		*FieldPtr = TransPtr->TrField(BufField);
	}
	else
	{
		radTField BufField(FieldPtr->FieldKey, FieldPtr->CompCriterium, FieldPtr->P, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		if(FldIntNeeded) // Plus this
		{
			BufField.Ib = ZeroVect; BufField.Ih = ZeroVect; BufField.NextP = FieldPtr->NextP;
		}

		B_comp_Or_NestedFor(&BufField, LocalNextIter);
		radTField BufField1 = BufField;

		BufField = radTField(FieldPtr->FieldKey, FieldPtr->CompCriterium, TransPtr->TrPoint_inv(BufField.P), ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		if(FldIntNeeded) // Plus this
		{
			BufField.Ib = ZeroVect; BufField.Ih = ZeroVect; 
			BufField.NextP = TransPtr->TrPoint_inv(BufField1.NextP);
		}

		int Mult = (*Iter).m;
		for(int km = 1; km < Mult-1; km++)
		{
			B_comp_Or_NestedFor(&BufField, LocalNextIter);
			BufField = TransPtr->TrField_inv(BufField);
			BufField.P = TransPtr->TrPoint_inv(BufField.P);
			if(FldIntNeeded) BufField.NextP = TransPtr->TrPoint_inv(BufField.NextP); // Plus this string
		}
		B_comp_Or_NestedFor(&BufField, LocalNextIter);
		radTField BufField2 = BufField;

		for(int km1 = 1; km1 < Mult; km1++)	BufField2 = TransPtr->TrField(BufField2);

		BufField1 += BufField2;
		*FieldPtr += BufField1;
	}
}

//-------------------------------------------------------------------------

void radTg3d::NestedFor_Energy(radTField* FieldPtr, const radTlphg::iterator& Iter)
{
	radTrans* TransPtr = (radTrans*)(((*Iter).Handler_g).rep);
	radTlphg::iterator LocalNextIter = Iter;
	LocalNextIter++;

	radTg3d* SourcePtr = (radTg3d*)(FieldPtr->HandleEnergyForceTorqueCompData.rep->hSource.rep);
	radTrans* InvTransPtr = new radTrans(*TransPtr); // This is to let handler work correctly
	InvTransPtr->Invert();
	radThg hInvTrans(InvTransPtr);

	double& LocEnergy = FieldPtr->Energy;
	TVector3d ZeroVect(0.,0.,0.);

	radTFieldKey& FieldPtrFieldKey = FieldPtr->FieldKey;

	if((*Iter).m == 1)
	{
		SourcePtr->AddTransform(1, hInvTrans);
		Energy_Or_NestedFor(FieldPtr, LocalNextIter);
		SourcePtr->EraseOuterTransform();
	}
	else
	{
		TVector3d LocForceMult1(0.,0.,0.), LocTorqueMult1(0.,0.,0.);
		double LocEnergyMult1 =0.;

		radTField BufField(FieldPtrFieldKey, FieldPtr->CompCriterium, FieldPtr->HandleEnergyForceTorqueCompData);
		BufField.Energy = FieldPtr->Energy;

		Energy_Or_NestedFor(&BufField, LocalNextIter);
		LocEnergyMult1 = BufField.Energy;

		BufField = radTField(FieldPtrFieldKey, FieldPtr->CompCriterium, FieldPtr->HandleEnergyForceTorqueCompData);
		BufField.Energy = FieldPtr->Energy;

		double &BufFieldEnergy = BufField.Energy;
		SourcePtr->AddTransform(1, hInvTrans);

		int Mult = (*Iter).m;
		for(int km = 1; km < Mult-1; km++)
		{
			Energy_Or_NestedFor(&BufField, LocalNextIter);
			SourcePtr->AddTransform(1, hInvTrans);
		}
		Energy_Or_NestedFor(&BufField, LocalNextIter);

		for(int km1 = 1; km1 < Mult; km1++) SourcePtr->EraseOuterTransform();

		LocEnergy += LocEnergyMult1 + BufFieldEnergy;
	}
}

//-------------------------------------------------------------------------

void radTg3d::NestedFor_IntOverShape(radTField* FieldPtr, const radTlphg::iterator& Iter)
{
	radTrans* TransPtr = (radTrans*)(((*Iter).Handler_g).rep);
	radTlphg::iterator LocalNextIter = Iter;
	LocalNextIter++;

	radTg3d* SourcePtr = (radTg3d*)((FieldPtr->ShapeIntDataPtr->HandleOfSource).rep);
	radTrans* InvTransPtr = new radTrans(*TransPtr); // This is to let handler work correctly
	InvTransPtr->Invert();
	radThg hInvTrans((radTg*)InvTransPtr);

	int LocLenVal = FieldPtr->ShapeIntDataPtr->IntegrandLength;
	TVector3d* LocVectArray = FieldPtr->ShapeIntDataPtr->VectArray;
	char* LocVectTypeArray = FieldPtr->ShapeIntDataPtr->VectTypeArray;
	TVector3d ZeroVect(0.,0.,0.);

	if((*Iter).m == 1)
	{
		for(int i=0; i<LocLenVal; i++)
		{
			if(LocVectTypeArray[i]=='r') LocVectArray[i] = TransPtr->TrVectField_inv(LocVectArray[i]);
			if(LocVectTypeArray[i]=='a') LocVectArray[i] = TransPtr->TrVectPoten_inv(LocVectArray[i]);
		}
		SourcePtr->AddTransform(1, hInvTrans);

		IntOverShape_Or_NestedFor(FieldPtr, LocalNextIter);

		for(int ii=0; ii<LocLenVal; ii++)
		{
			if(LocVectTypeArray[ii]=='r') LocVectArray[ii] = TransPtr->TrVectField(LocVectArray[ii]);
			if(LocVectTypeArray[ii]=='a') LocVectArray[ii] = TransPtr->TrVectPoten(LocVectArray[ii]);
		}
		SourcePtr->EraseOuterTransform();
	}
	else
	{
		radTStructForShapeInt LocStructForShapeInt1 = *(FieldPtr->ShapeIntDataPtr);
		radTStructForShapeInt LocStructForShapeInt = *(FieldPtr->ShapeIntDataPtr);
		TVector3d* LocVectArrayMult1 = new TVector3d[LocLenVal];
		TVector3d* LocVectArrayMult = new TVector3d[LocLenVal];
		LocStructForShapeInt1.VectArray = LocVectArrayMult1;
		LocStructForShapeInt.VectArray = LocVectArrayMult;
		for(int i=0; i<LocLenVal; i++)
		{
			LocVectArrayMult1[i] = LocVectArrayMult[i] = ZeroVect;
		}
		radTField BufField(FieldPtr->FieldKey, FieldPtr->CompCriterium, &LocStructForShapeInt1);

		IntOverShape_Or_NestedFor(&BufField, LocalNextIter);

		BufField = radTField(FieldPtr->FieldKey, FieldPtr->CompCriterium, &LocStructForShapeInt);
		SourcePtr->AddTransform(1, hInvTrans);

		int Mult = (*Iter).m;
		for(int km = 1; km < Mult-1; km++)
		{
			IntOverShape_Or_NestedFor(&BufField, LocalNextIter);

			for(int ii=0; ii<LocLenVal; ii++)
			{
				if(LocVectTypeArray[ii]=='r') LocVectArrayMult[ii] = TransPtr->TrVectField_inv(LocVectArrayMult[ii]);
				if(LocVectTypeArray[ii]=='a') LocVectArrayMult[ii] = TransPtr->TrVectPoten_inv(LocVectArrayMult[ii]);
			}
			SourcePtr->AddTransform(1, hInvTrans);
		}

		IntOverShape_Or_NestedFor(&BufField, LocalNextIter);
		for(int km1 = 1; km1 < Mult; km1++)
		{
			for(int iii=0; iii<LocLenVal; iii++)
			{
				if(LocVectTypeArray[iii]=='r') LocVectArrayMult[iii] = TransPtr->TrVectField(LocVectArrayMult[iii]);
				if(LocVectTypeArray[iii]=='a') LocVectArrayMult[iii] = TransPtr->TrVectPoten(LocVectArrayMult[iii]);
			}
			SourcePtr->EraseOuterTransform();
		}
		for(int iiii=0; iiii<LocLenVal; iiii++) LocVectArray[iiii] += LocVectArrayMult1[iiii] + LocVectArrayMult[iiii];
		
		delete[] LocVectArrayMult1;
		delete[] LocVectArrayMult;
	}
}

//-------------------------------------------------------------------------

void radTg3d::B_intCompFinNum(radTField* FieldPtr)
{// This uses Newton method (n=3)
	const double IntegWeight[] = {3./8., 9./8., 9./8., 3./4.};

	TVector3d VectV = FieldPtr->NextP - FieldPtr->P;
	double Fact = sqrt(VectV.x*VectV.x + VectV.y*VectV.y + VectV.z*VectV.z);

	radTFieldKey LocFieldKey;
	short LocIb_ = FieldPtr->FieldKey.Ib_;
	short LocIh_ = FieldPtr->FieldKey.Ih_;
	LocFieldKey.B_ = LocIb_;
	LocFieldKey.H_ = LocIh_;

	radTCompCriterium LocCompCriterium;
	LocCompCriterium = FieldPtr->CompCriterium;

	TVector3d ZeroVect(0.,0.,0.), S_forB, S_forH, GenS_forB(0.,0.,0.), GenS_forH(0.,0.,0.),
			  IntForB(1.E+23, 1.E+23, 1.E+23), IntForH(1.E+23, 1.E+23, 1.E+23),
			  PrIntForB, PrIntForH;

	radTField LocField(LocFieldKey, LocCompCriterium, FieldPtr->P, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);

	double t_min = 0.;
	double Step_t, t;
	short IndForWeight, IndForPass;
	short NotFirstPass = 0;

	int AmOfPoi = 4;
	int AmOfPoi_mi_1;
	double PrecParamB_int, PrecParamH_int, PrecParamInt;
	PrecParamB_int = PrecParamH_int = 0.;
	PrecParamInt = 1.E+23; 

	while(PrecParamInt > FieldPtr->CompCriterium.AbsPrecB_int)
	{
		AmOfPoi_mi_1 = AmOfPoi - 1;
		Step_t = 1./AmOfPoi_mi_1;
		t = t_min;

		PrIntForB = IntForB; PrIntForH = IntForH;
		S_forB = ZeroVect; S_forH = ZeroVect;

		IndForWeight = IndForPass = 0;

		for(int i=0; i<AmOfPoi; i++)
		{
			if(IndForPass==2) IndForPass = 0;
			if(IndForWeight==4) IndForWeight = 1;
			if(NotFirstPass && (IndForPass==0)) goto BottomOfThisLoop;
			if(i==AmOfPoi_mi_1) IndForWeight = 0;

			LocField.P = FieldPtr->P + (t * VectV);
			if(LocIb_) LocField.B = ZeroVect;
			if(LocIh_) LocField.H = ZeroVect;

			B_comp(&LocField);
						
			if(LocIb_) S_forB += IntegWeight[IndForWeight] * LocField.B;
			if(LocIh_) S_forH += IntegWeight[IndForWeight] * LocField.H;

BottomOfThisLoop:
			IndForPass++; IndForWeight++;
			t += Step_t;
		}

		if(LocIb_)
		{
			GenS_forB += S_forB; 
			IntForB = Step_t * GenS_forB;
			PrecParamB_int = Fact * Max( Max( Abs(IntForB.x-PrIntForB.x), Abs(IntForB.y-PrIntForB.y)), Abs(IntForB.z-PrIntForB.z));
		}
		if(LocIh_)
		{
			GenS_forH += S_forH; 
			IntForH = Step_t * GenS_forH;
			PrecParamH_int = Fact * Max( Max( Abs(IntForH.x-PrIntForH.x), Abs(IntForH.y-PrIntForH.y)), Abs(IntForH.z-PrIntForH.z));
		}

		PrecParamInt = Max(PrecParamB_int, PrecParamH_int);
		AmOfPoi = AmOfPoi_mi_1 * 2 + 1;
		NotFirstPass = 1;
	}

	if(LocIb_) FieldPtr->Ib += Fact * IntForB;
	if(LocIh_) FieldPtr->Ih += Fact * IntForH;
}

//-------------------------------------------------------------------------

void radTg3d::NormStressTensor(radTField* FieldPtr)
{
	TVector3d ZeroVect(0.,0.,0.);
	FieldPtr->FieldKey.Force_= 0; 
	short PrevB_= FieldPtr->FieldKey.B_; FieldPtr->FieldKey.B_= 1;
	FieldPtr->B = ZeroVect;

	((radTg3d*)(FieldPtr->ShapeIntDataPtr->HandleOfSource.rep))->B_genComp(FieldPtr);

	const double ConForStrTensInSI = 1.E-06/(4*3.14159265358979*1.E-07); // (10^(-3))^2/mu0
	TVector3d LocB = FieldPtr->B;

	//Out normal projection of the Maxwell Stress Tensor
	*(FieldPtr->ShapeIntDataPtr->VectArray) =
		ConForStrTensInSI*((LocB*FieldPtr->ShapeIntDataPtr->Normal)*LocB
		-(0.5*(LocB*LocB))*FieldPtr->ShapeIntDataPtr->Normal);

	FieldPtr->FieldKey.Force_= 1; FieldPtr->FieldKey.B_= PrevB_;
}

//-------------------------------------------------------------------------

void radTg3d::DumpTransApplied(std::ostream& o) // Porting
{
	o << endl;
	o << "   Transformations applied: ";

	if(g3dListOfTransform.empty()) { o << "None";}
	else
	{
		for(radTlphg::reverse_iterator iter = g3dListOfTransform.rbegin();
			iter != g3dListOfTransform.rend(); ++iter)
		{
			o << endl;

			long TrElemKey = rad.RetrieveElemKey((*iter).Handler_g.rep);
			if(TrElemKey > 0)
			{
				o << "      Index " << TrElemKey << ": ";
			}
			else
			{
				o << "      Index " << "n/a" << ": ";
			}
			((*iter).Handler_g.rep)->Dump(o, 1);
			o << ";  Multiplicity: " << (*iter).m;
		}
	}
}

//-------------------------------------------------------------------------

void radTg3d::ActualEnergyForceTorqueComp(radTField* FieldPtr)
{
	TVector3d ZeroVect(0.,0.,0.);
	const double DeltaL = 1.E-01; // mm
	//const double DeltaL = 1.; // mm
	const double DeltaTeta = 1.E-03; // rad

	radTField LocField(FieldPtr->FieldKey, FieldPtr->CompCriterium, FieldPtr->HandleEnergyForceTorqueCompData);
	double& Enr = LocField.Energy;
	double Ex1, Ex2, Ey1, Ey2, Ez1, Ez2;

	if(FieldPtr->FieldKey.Energy_)
	{
		Enr = 0.;
		(g3dListOfTransform.empty())? SimpleEnergyComp(&LocField) : NestedFor_Energy(&LocField, g3dListOfTransform.begin());
		FieldPtr->Energy += LocField.Energy;
	}
	if(FieldPtr->FieldKey.ForceEnr_)
	{
		LocField.CompCriterium.BasedOnWorstRelPrec = 1;

		TVector3d DeltaX(-DeltaL/2., 0., 0.), DeltaY(0., -DeltaL/2., 0.),  DeltaZ(0., 0., -DeltaL/2.);
		radTrans *SmallTranslXPtr = new radTrans(), *SmallTranslYPtr = new radTrans(), *SmallTranslZPtr = new radTrans();
		radThg hTranslX(SmallTranslXPtr), hTranslY(SmallTranslYPtr), hTranslZ(SmallTranslZPtr);
		SmallTranslXPtr->SetupTranslation(DeltaX); SmallTranslYPtr->SetupTranslation(DeltaY); SmallTranslZPtr->SetupTranslation(DeltaZ);

		AddTransform(1, hTranslX); 
		//AddTransform_OtherSide(1, hTranslX);
		Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ex1 = Enr;
		SmallTranslXPtr->Invert(); Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ex2 = Enr;
		EraseOuterTransform(); 
		//EraseInnerTransform(); 

		AddTransform(1, hTranslY);
		//AddTransform_OtherSide(1, hTranslY);
		Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ey1 = Enr;
		SmallTranslYPtr->Invert(); Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ey2 = Enr;
		EraseOuterTransform(); 
		//EraseInnerTransform(); 

		AddTransform(1, hTranslZ); 
		//AddTransform_OtherSide(1, hTranslZ);
		Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ez1 = Enr;
		SmallTranslZPtr->Invert(); Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ez2 = Enr;
		EraseOuterTransform();
		//EraseInnerTransform(); 

		double Const = 1000./DeltaL;
		TVector3d LocForce(Const*(Ex1 - Ex2), Const*(Ey1 - Ey2), Const*(Ez1 - Ez2));
		if(FieldPtr->FieldKey.ForceEnr_) FieldPtr->Force += LocForce;
	}
	if(FieldPtr->FieldKey.Torque_)
	{
		LocField.CompCriterium.BasedOnWorstRelPrec = 1;

		radTrans *SmallRotXaxPtr = new radTrans(), *SmallRotYaxPtr = new radTrans(), *SmallRotZaxPtr = new radTrans();
		radThg hRotX(SmallRotXaxPtr), hRotY(SmallRotYaxPtr), hRotZ(SmallRotZaxPtr);

		TVector3d OrtX(1.,0.,0.), OrtY(0.,1.,0.), OrtZ(0.,0.,1.);
		SmallRotXaxPtr->SetupRotation(FieldPtr->P, OrtX, -DeltaTeta/2.);
		SmallRotYaxPtr->SetupRotation(FieldPtr->P, OrtY, -DeltaTeta/2.);
		SmallRotZaxPtr->SetupRotation(FieldPtr->P, OrtZ, -DeltaTeta/2.);

		AddTransform(1, hRotX); 
		//AddTransform_OtherSide(1, hRotX);
		Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ex1 = Enr;
		SmallRotXaxPtr->Invert(); Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ex2 = Enr;
		EraseOuterTransform(); 
		//EraseInnerTransform(); 

		AddTransform(1, hRotY); 
		//AddTransform_OtherSide(1, hRotY);
		Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ey1 = Enr;
		SmallRotYaxPtr->Invert(); Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ey2 = Enr;
		EraseOuterTransform(); 
		//EraseInnerTransform(); 

		AddTransform(1, hRotZ); 
		//AddTransform_OtherSide(1, hRotZ);
		Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ez1 = Enr;
		SmallRotZaxPtr->Invert(); Enr = 0.; NestedFor_Energy(&LocField, g3dListOfTransform.begin()); Ez2 = Enr;
		EraseOuterTransform();
		//EraseInnerTransform(); 

		double Const = 1000./DeltaTeta;
		TVector3d LocTorque(Const*(Ex1 - Ex2), Const*(Ey1 - Ey2), Const*(Ez1 - Ez2));
		if(FieldPtr->FieldKey.Torque_) FieldPtr->Torque += LocTorque;
	}
}

//-------------------------------------------------------------------------

void radTg3d::EnergyForceTorqueComp(radTField* FieldPtr)
{
	const double RelTol = 1.E-09;
	radTStructForEnergyForceTorqueComp* EnFrcTrqCompDataPtr = FieldPtr->HandleEnergyForceTorqueCompData.rep;

	if(EnFrcTrqCompDataPtr->AutoDestSubdivision)
	{
		EnergyForceTorqueCompAutoDestSubd(FieldPtr);
	}
	else
	{
		if(!DuplicateItself(EnFrcTrqCompDataPtr->hDest, EnFrcTrqCompDataPtr->radPtr, 0))
		{ 
			EnFrcTrqCompDataPtr->SomethingIsWrong = 1; return;
		}
		radTg3d* DuplDestPtr = (radTg3d*)(EnFrcTrqCompDataPtr->hDest.rep);

		if((fabs(EnFrcTrqCompDataPtr->DestSubdivArray[0]-1.)<RelTol) && (fabs(EnFrcTrqCompDataPtr->DestSubdivArray[2]-1.)<RelTol) && (fabs(EnFrcTrqCompDataPtr->DestSubdivArray[4]-1.)<RelTol))
		{
			DuplDestPtr->ActualEnergyForceTorqueComp(FieldPtr); return;
		}
		else
		{
			radTSubdivOptions SubdivOptions;
			SubdivOptions.SubdivisionFrame = 0;
			SubdivOptions.SubdivisionParamCode = 0;
			SubdivOptions.SubdivideCoils = 1;
			SubdivOptions.PutNewStuffIntoGenCont = 0;

			if(!DuplDestPtr->SubdivideItself(EnFrcTrqCompDataPtr->DestSubdivArray, EnFrcTrqCompDataPtr->hDest, EnFrcTrqCompDataPtr->radPtr, &SubdivOptions))
			{ 
				EnFrcTrqCompDataPtr->SomethingIsWrong = 1; return;
			}
			((radTg3d*)(EnFrcTrqCompDataPtr->hDest.rep))->ActualEnergyForceTorqueComp(FieldPtr); 
		}
	}
}

//-------------------------------------------------------------------------

char radTg3d::CheckIfMoreEnrFrcTrqCompNeededAndUpdate(radTField* WorkFieldPtr, radTField* TotFieldPtr)
{
	radTFieldKey &FieldKey = WorkFieldPtr->FieldKey;
	char EnergyCompNotNeeded = 1, ForceCompNotNeeded = 1, TorqueCompNotNeeded = 1;
	
	radTg3d* LocDestPtr = (radTg3d*)(WorkFieldPtr->HandleEnergyForceTorqueCompData.rep->hDest.rep);
	radTAuxCompDataG3D* AuxCompDataPtr = LocDestPtr->HandleAuxCompData.rep;

	if(FieldKey.Energy_)
	{
		double &Old = AuxCompDataPtr->Energy, &New = WorkFieldPtr->Energy;
		double Diff = New - Old;
		EnergyCompNotNeeded = (Abs(Diff) <= WorkFieldPtr->CompCriterium.AbsPrecEnergy)? 1 : 0;
		if(!EnergyCompNotNeeded)
		{
			Old = New; // Is it actually necessary?
		}
		TotFieldPtr->Energy += Diff;
	}
	if(FieldKey.ForceEnr_)
	{
		TVector3d &Old = AuxCompDataPtr->Force, &New = WorkFieldPtr->Force;
		TVector3d Diff = New - Old;
		ForceCompNotNeeded = (NormAbs(Diff) <= WorkFieldPtr->CompCriterium.AbsPrecForce)? 1 : 0;
		if(!ForceCompNotNeeded)
		{
			Old = New; // Is it actually necessary?
		}
		TotFieldPtr->Force += Diff;
	}
	if(FieldKey.Torque_)
	{
		TVector3d &Old = AuxCompDataPtr->Torque, &New = WorkFieldPtr->Torque;
		TVector3d Diff = New - Old;
		TorqueCompNotNeeded = (NormAbs(Diff) <= WorkFieldPtr->CompCriterium.AbsPrecTorque)? 1 : 0;
		if(!TorqueCompNotNeeded)
		{
			Old = New; // Is it actually necessary?
		}
		TotFieldPtr->Torque += Diff;
	}
	return !(EnergyCompNotNeeded && ForceCompNotNeeded && TorqueCompNotNeeded);
}

//-------------------------------------------------------------------------

void radTg3d::ActualEnergyForceTorqueCompWithAdd(radTField* FieldPtr)
{// This is for everything except for Groups and their childs.
	TVector3d ZeroVect(0.,0.,0.);
	radTField LocField = *FieldPtr;
	LocField.Force = LocField.Torque = ZeroVect; LocField.Energy = 0.;
	
	ActualEnergyForceTorqueComp(&LocField);

	if(HandleAuxCompData.rep == 0) CreateAuxCompData();
	HandleAuxCompData.rep->StoreDataFromField(&LocField);

	radTFieldKey &FieldKey = FieldPtr->FieldKey;
	if(FieldKey.Energy_) FieldPtr->Energy += LocField.Energy;
	if(FieldKey.ForceEnr_) FieldPtr->Force += LocField.Force;
	if(FieldKey.Torque_) FieldPtr->Torque += LocField.Torque;
}

//-------------------------------------------------------------------------

int radTg3d::ProceedNextStepEnergyForceTorqueComp(double* SubdArr, radThg& HandleOfThis, radTField* LocFieldPtr, radTField* FieldPtr, char& SubdNeed, char XorYorZ)
{
	radThg hOld = HandleOfThis;

	radTSubdivOptions SubdivOptions;
	SubdivOptions.SubdivisionFrame = 0;
	SubdivOptions.SubdivisionParamCode = 0;
	SubdivOptions.SubdivideCoils = 1;
	SubdivOptions.PutNewStuffIntoGenCont = 0;

	if(!((radTg3d*)(hOld.rep))->SubdivideItself(SubdArr, HandleOfThis, FieldPtr->HandleEnergyForceTorqueCompData.rep->radPtr, &SubdivOptions))
	{
		FieldPtr->HandleEnergyForceTorqueCompData.rep->SomethingIsWrong = 1; return 0;
	}
	TVector3d &LocForce = LocFieldPtr->Force, &LocTorque = LocFieldPtr->Torque;
	LocFieldPtr->Energy = LocForce.x = LocForce.y = LocForce.z = LocTorque.x = LocTorque.y = LocTorque.z = 0.;

	LocFieldPtr->HandleEnergyForceTorqueCompData.rep->hDest = HandleOfThis;

	radTg3d* NewGroupPtr = (radTg3d*)(HandleOfThis.rep);
	NewGroupPtr->ActualEnergyForceTorqueCompWithAdd(LocFieldPtr);
	NewGroupPtr->HandleAuxCompData = ((radTg3d*)(hOld.rep))->HandleAuxCompData; // Move to constructors of subdivided items

	NewGroupPtr->SetupFurtherSubdInd(NewGroupPtr->HandleAuxCompData.rep->SubdNeedInd);

	SubdNeed = CheckIfMoreEnrFrcTrqCompNeededAndUpdate(LocFieldPtr, FieldPtr);
	if(!SubdNeed) HandleOfThis = hOld;

	radTg3d* GroupInPlaceOfThis_Or_This = (radTg3d*)(HandleOfThis.rep);
	GroupInPlaceOfThis_Or_This->MarkFurtherSubdNeed1D(SubdNeed, XorYorZ);
	return 1;
}

//-------------------------------------------------------------------------

int radTg3d::NextStepEnergyForceTorqueComp(double* TotSubdArr, radThg& HandleOfThis, radTField* FieldPtr, char& MoreSubdNeeded)
{
	char SubdNeedX, SubdNeedY, SubdNeedZ;
	HandleAuxCompData.rep->ShowSubdNeed(SubdNeedX, SubdNeedY, SubdNeedZ);
	if(!(SubdNeedX || SubdNeedY || SubdNeedZ)) return 1;

	double SubdArrX[] = {TotSubdArr[0], TotSubdArr[1], 1.,1.,1.,1.};
	double SubdArrY[] = {1.,1., TotSubdArr[2], TotSubdArr[3], 1.,1.};
	double SubdArrZ[] = {1.,1.,1.,1., TotSubdArr[4], TotSubdArr[5]};

	radTStructForEnergyForceTorqueComp* EnFrcTrqCompDataPtr = FieldPtr->HandleEnergyForceTorqueCompData.rep;
	char SubdivideCoils =1, PutNewStuffIntoGenCont =0;

	radTField LocField = *FieldPtr;
	LocField.HandleEnergyForceTorqueCompData = radTHandleStructForEnergyForceTorqueComp(new radTStructForEnergyForceTorqueComp(*(FieldPtr->HandleEnergyForceTorqueCompData.rep)));

	if(SubdNeedX)
		if(!ProceedNextStepEnergyForceTorqueComp(SubdArrX, HandleOfThis, &LocField, FieldPtr, SubdNeedX, 'x')) return 0;
	if(SubdNeedY)
	{
		radThg hgOld = HandleOfThis;
		radTg3d* DestPtr = (radTg3d*)(hgOld.rep);
		if(!DestPtr->ProceedNextStepEnergyForceTorqueComp(SubdArrY, HandleOfThis, &LocField, FieldPtr, SubdNeedY, 'y')) return 0;
	}
	if(SubdNeedZ)
	{
		radThg hgOld = HandleOfThis;
		radTg3d* DestPtr = (radTg3d*)(hgOld.rep);
		if(!DestPtr->ProceedNextStepEnergyForceTorqueComp(SubdArrZ, HandleOfThis, &LocField, FieldPtr, SubdNeedZ, 'z')) return 0;
	}
	MoreSubdNeeded = SubdNeedX || SubdNeedY || SubdNeedZ;
	return 1;
}

//-------------------------------------------------------------------------

void radTg3d::EnergyForceTorqueCompAutoDestSubd(radTField* FieldPtr)
{
	radTStructForEnergyForceTorqueComp* EnFrcTrqCompDataPtr = FieldPtr->HandleEnergyForceTorqueCompData.rep;

	if(!DuplicateItself(EnFrcTrqCompDataPtr->hDest, EnFrcTrqCompDataPtr->radPtr, 0))
	{ 
		EnFrcTrqCompDataPtr->SomethingIsWrong = 1; return;
	}

	TVector3d ZeroVect(0.,0.,0.);
	double TotSubdArray[] = {2.,1.,2.,1.,2.,1}; 

	radTField LocField = *FieldPtr;
	LocField.Force = LocField.Torque = ZeroVect; LocField.Energy = 0.;

	radTg3d* DestPtr = (radTg3d*)(EnFrcTrqCompDataPtr->hDest.rep);
	DestPtr->ActualEnergyForceTorqueCompWithAdd(&LocField);
	DestPtr->MarkFurtherSubdNeed(1, 1, 1);

	char MoreSubdNeeded = 1;
	while(MoreSubdNeeded)
	{
		radThg hDestOld = LocField.HandleEnergyForceTorqueCompData.rep->hDest; // To prevent from automatic deletion through handle
		if(!((radTg3d*)(hDestOld.rep))->NextStepEnergyForceTorqueComp(TotSubdArray, EnFrcTrqCompDataPtr->hDest, &LocField, MoreSubdNeeded))
		{
			EnFrcTrqCompDataPtr->SomethingIsWrong = 1; return;
		}
	}
	radTFieldKey &FieldKey = FieldPtr->FieldKey;
	if(FieldKey.Energy_) FieldPtr->Energy += LocField.Energy;
	if(FieldKey.ForceEnr_) FieldPtr->Force += LocField.Force;
	if(FieldKey.Torque_) FieldPtr->Torque += LocField.Torque;
}

//-------------------------------------------------------------------------

void radTg3d::CheckAxesExchangeForSubdInLabFrame(double* SubdivArray, char& ConversionToPolyhedronsIsNeeded)
{
	ConversionToPolyhedronsIsNeeded = 0;

	radTrans ResTransf;
	short SomethingFound = 0;
	FindResTransfWithMultOne(ResTransf, SomethingFound);
	if(SomethingFound) 
	{
		const double ZeroTol = 1.E-13;
		TVector3d ex(1.,0.,0.), ey(0.,1.,0.), ez(0.,0.,1.);
		TVector3d ex1 = ResTransf.TrBiPoint_inv(ex);
		TVector3d ey1 = ResTransf.TrBiPoint_inv(ey);
		TVector3d ez1 = ResTransf.TrBiPoint_inv(ez);

		if(PracticallyEqualOrAnti(ex, ex1, ZeroTol))
		{
			if(!(PracticallyEqualOrAnti(ey, ey1, ZeroTol) || PracticallyEqualOrAnti(ey, ez1, ZeroTol)))
				ConversionToPolyhedronsIsNeeded = 1;
		}
		else if(PracticallyEqualOrAnti(ex, ey1, ZeroTol))
		{
			if(!(PracticallyEqualOrAnti(ey, ex1, ZeroTol) || PracticallyEqualOrAnti(ey, ez1, ZeroTol)))
				ConversionToPolyhedronsIsNeeded = 1;
		}
		else if(PracticallyEqualOrAnti(ex, ez1, ZeroTol))
		{
			if(!(PracticallyEqualOrAnti(ey, ex1, ZeroTol) || PracticallyEqualOrAnti(ey, ey1, ZeroTol)))
				ConversionToPolyhedronsIsNeeded = 1;
		}
		else ConversionToPolyhedronsIsNeeded = 1;

		if(!ConversionToPolyhedronsIsNeeded)
		{
			double &kx = SubdivArray[0], &ky = SubdivArray[2],  &kz = SubdivArray[4];
			double &qx = SubdivArray[1], &qy = SubdivArray[3],  &qz = SubdivArray[5];

			TVector3d kx_ex1 = kx*ex1, ky_ey1 = ky*ey1, kz_ez1 = kz*ez1;
			TVector3d k1Tot = kx_ex1 + ky_ey1 + kz_ez1;
			TVector3d qx_ex1 = qx*ex1, qy_ey1 = qy*ey1, qz_ez1 = qz*ez1;
			TVector3d q1Tot = qx_ex1 + qy_ey1 + qz_ez1;

			kx = ex*k1Tot; qx = ex*q1Tot;
			if(kx < 0.) { kx = -kx; qx = -1./qx;}
			ky = ey*k1Tot; qy = ey*q1Tot;
			if(ky < 0.) { ky = -ky; qy = -1./qy;}
			kz = ez*k1Tot; qz = ez*q1Tot;
			if(kz < 0.) { kz = -kz; qz = -1./qz;}
		}
	}
}

//-------------------------------------------------------------------------

int radTg3d::TransferSubdivisionStructToLocalFrame(TVector3d& InNormal, TVector3d* InPoints, int AmOfPoints, radTSubdivOptions* pSubdivOptions, TVector3d& OutNormal, TVector3d*& OutPoints)
{// AmOfPoints = AmOfPieces_mi_1
	OutNormal = InNormal;
	OutPoints = InPoints; // This is important
		
	radTrans ResTransf;
	short SomethingFound = 0;

	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty()))
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	else if(ConsiderOnlyWithTrans && (!g3dListOfTransform.empty()))
	{
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
	}

	radTSend Send;
	if(SomethingFound) 
	{
		OutNormal = ResTransf.TrBiPoint_inv(InNormal);

		OutPoints = new TVector3d[AmOfPoints];
		if(OutPoints == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		TVector3d *tOut = OutPoints, *tIn = InPoints;
		for(int k=0; k<AmOfPoints; k++) *(tOut++) = ResTransf.TrPoint_inv(*(tIn++));
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTg3d::FindEllipticCoordOfPoint(radTCylindricSubdivSpec* pCylSpec, TVector3d& Point, double& a, double& Phi)
{
	TVector3d& CylAxVect = pCylSpec->CylAxVect;
	TVector3d& PointOnCylAx = pCylSpec->PointOnCylAx;
	TVector3d& PointOnEllAx = pCylSpec->PointOnEllAx;
	double Ratio = 1./pCylSpec->EllAxRatio;

	TVector3d PointOnEllAx_mi_PointOnCylAx = PointOnEllAx - PointOnCylAx;
	TVector3d Ex = PointOnEllAx_mi_PointOnCylAx - (PointOnEllAx_mi_PointOnCylAx*CylAxVect)*CylAxVect;
	double InvLenEx = 1./sqrt(Ex.x*Ex.x + Ex.y*Ex.y + Ex.z*Ex.z);
	Ex = InvLenEx*Ex;
	TVector3d Ey = CylAxVect^Ex;

	TVector3d Point_mi_PointOnCylAx = Point - PointOnCylAx;
	TVector3d PointOrtComp = Point_mi_PointOnCylAx - (Point_mi_PointOnCylAx*CylAxVect)*CylAxVect;

	double xLoc = PointOrtComp*Ex, yLoc = PointOrtComp*Ey;

	const double HalfPI = 1.5707963267949;
	const double TolForZero = 1.E-12;
	if(fabs(xLoc) < TolForZero)
	{
		if(yLoc > TolForZero) { Phi = HalfPI; a = Ratio*yLoc; return;}
		else if(yLoc < -TolForZero) { Phi = 3.*HalfPI; a = -Ratio*yLoc; return;}
		else { a = 0.; Phi = 0.; return;}
	}
	else if(xLoc > TolForZero)
	{
		if(yLoc > TolForZero) Phi = atan(Ratio*yLoc/xLoc);
		else if(yLoc < -TolForZero) Phi = 4.*HalfPI + atan(Ratio*yLoc/xLoc);
		else { a = xLoc; Phi = 0.; return;}
	}
	else if(xLoc < -TolForZero)
	{
		Phi = 2.*HalfPI + atan(Ratio*yLoc/xLoc);
	}

	a = xLoc/cos(Phi);
}

//-------------------------------------------------------------------------

double radTg3d::EstimateLengthAlongEllipse(double a, double Ratio, double PhMin, double PhMax)
{
	const double PointsPerRad = 4./1.5708;
	double DelPhi = PhMax - PhMin;
	int Np = int(DelPhi*PointsPerRad) + 1;

	double Buf = Np*0.5;
	if((Buf - int(Buf)) < 0.4) Np++;
	if(Np < 3) Np = 3;

	double Inv_k = Ratio;
	double Inv_ke2 = Inv_k*Inv_k;

	double SinPh = sin(PhMin), CosPh = cos(PhMin);
	double f0 = sqrt(SinPh*SinPh + Inv_ke2*CosPh*CosPh);
	SinPh = sin(PhMax), CosPh = cos(PhMax);
	double fn = sqrt(SinPh*SinPh + Inv_ke2*CosPh*CosPh);

	double StPhi = DelPhi/(Np - 1);
	double Phi = PhMin + StPhi;
	double Sum1=0., Sum2=0.;
	int nPass = (Np - 3) >> 1;
	for(int k=0; k<nPass; k++)
	{
		SinPh = sin(Phi), CosPh = cos(Phi);
		Sum1 += sqrt(SinPh*SinPh + Inv_ke2*CosPh*CosPh);
		Phi += StPhi;
		SinPh = sin(Phi), CosPh = cos(Phi);
		Sum2 += sqrt(SinPh*SinPh + Inv_ke2*CosPh*CosPh);
		Phi += StPhi;
	}
	SinPh = sin(Phi), CosPh = cos(Phi);
	Sum1 += sqrt(SinPh*SinPh + Inv_ke2*CosPh*CosPh);
	return 0.33333*StPhi*a*(f0 + fn + 4.*Sum1 + 2.*Sum2);
}

//-------------------------------------------------------------------------

int radTg3d::CreateFromSym(radThg& In_hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{
	int MaxMult = 0;
	for(radTlphg::const_iterator iter = g3dListOfTransform.begin(); iter != g3dListOfTransform.end(); ++iter)
		if((*iter).m > MaxMult) MaxMult = (*iter).m;

	if(MaxMult <= 1) 
	//if((MaxMult <= 1) && (g3dListOfTransform.size() <= 1)) //OC061007_BNL
	{
		return DuplicateItself(In_hg, radPtr, PutNewStuffIntoGenCont);
	}
	else
	{
		radTSend Send;
		radTGroup* pGroup = new radTGroup();
		if(pGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgGroup(pGroup);

		radTrans BaseTrans;
		BaseTrans.SetupIdent();
		if(!NestedFor_CreateFromSym(pGroup, radPtr, PutNewStuffIntoGenCont, &BaseTrans, g3dListOfTransform.begin())) return 0;

		if(!pGroup->GroupMapOfHandlers.empty()) 
		{
			pGroup->IsGroupMember = IsGroupMember;
			pGroup->ConsiderOnlyWithTrans = ConsiderOnlyWithTrans;
			pGroup->MessageChar = MessageChar;
			In_hg = hgGroup;
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTg3d::NestedFor_CreateFromSym(radTGroup* pGroup, radTApplication* radPtr, char PutNewStuffIntoGenCont, radTrans* BaseTransPtr, const radTlphg::iterator& Iter)
{
	radTrans* TransPtr = (radTrans*)(((*Iter).Handler_g).rep);
	radTlphg::iterator LocalNextIter = Iter;
	LocalNextIter++;

	radTrans LocTotTrans = *BaseTransPtr;

	if((*Iter).m == 1)
	{
		LocTotTrans = Product(LocTotTrans, *TransPtr);
		if(!CreateAndAddToGroupOrNestedFor(pGroup, radPtr, PutNewStuffIntoGenCont, &LocTotTrans, LocalNextIter)) return 0;
	}
	else
	{
		if(!CreateAndAddToGroupOrNestedFor(pGroup, radPtr, PutNewStuffIntoGenCont, &LocTotTrans, LocalNextIter)) return 0;
		int Mult = (*Iter).m;
		for(int km = 1; km < Mult; km++)
		{
			LocTotTrans = Product(LocTotTrans, *TransPtr);
			if(!CreateAndAddToGroupOrNestedFor(pGroup, radPtr, PutNewStuffIntoGenCont, &LocTotTrans, LocalNextIter)) return 0;
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTg3d::CreateAndAddToGroupOrNestedFor(radTGroup* pGroup, radTApplication* radPtr, char PutNewStuffIntoGenCont, radTrans* BaseTransPtr, const radTlphg::iterator& Iter)
{
	if(Iter == g3dListOfTransform.end())
	{
		radThg hgNew;
		if(!DuplicateItself(hgNew, radPtr, PutNewStuffIntoGenCont)) return 0;

		radTg3d* g3dDplPtr = (radTg3d*)(hgNew.rep);
		g3dDplPtr->EraseAllTransformations();

		double RelTol = 1.E-12;		
		if(!BaseTransPtr->IsIdent(RelTol))
		{
			radTSend Send;
			radTrans* pNewTrans = new radTrans(*BaseTransPtr);
			if(pNewTrans == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			radThg hgTrans(pNewTrans);
			int LocMult = 1;
			g3dDplPtr->AddTransform(LocMult, hgTrans);
		}

		if(PutNewStuffIntoGenCont)
		{
			pGroup->AddElement(radPtr->AddElementToContainer(hgNew), hgNew);
		}
		else
		{
			int LocKey = (int)(pGroup->GroupMapOfHandlers.size()) + 1;
			pGroup->AddElement(LocKey, hgNew);
		}
		return 1;
	}
	else return NestedFor_CreateFromSym(pGroup, radPtr, PutNewStuffIntoGenCont, BaseTransPtr, Iter);
}

//-------------------------------------------------------------------------

double radTg3d::VolumeWithSym()
{
	int TotMult = TotalMultiplicity();
	return TotMult*Volume();
}

//-------------------------------------------------------------------------

void radTg3d::Limits(radTrans* pExtTr, double* LimArr)
{
	LimArr[0] = 1E+23; LimArr[1] = -1E+23;
	LimArr[2] = 1E+23; LimArr[3] = -1E+23;
	LimArr[4] = 1E+23; LimArr[5] = -1E+23;

	if(g3dListOfTransform.empty()) 
	{
		LimitsAtTransform(pExtTr, LimArr);
		return;
	}
	else
	{
		radTvhg FlatTransforms;
		FlattenSpaceTransforms(FlatTransforms);
		int SizeFlatTransforms = (int)(FlatTransforms.size());

		double AuxLimArr[6];
		for(int k=0; k<SizeFlatTransforms; k++)
		{
			radTrans* pTr = (radTrans*)(FlatTransforms[k].rep);
			if(pExtTr != 0)
			{
				*pTr = Product(*pExtTr, *pTr);
			}

			LimitsAtTransform(pTr, AuxLimArr);
			
			if(LimArr[0] > AuxLimArr[0]) LimArr[0] = AuxLimArr[0];
			if(LimArr[1] < AuxLimArr[1]) LimArr[1] = AuxLimArr[1];
			if(LimArr[2] > AuxLimArr[2]) LimArr[2] = AuxLimArr[2];
			if(LimArr[3] < AuxLimArr[3]) LimArr[3] = AuxLimArr[3];
			if(LimArr[4] > AuxLimArr[4]) LimArr[4] = AuxLimArr[4];
			if(LimArr[5] < AuxLimArr[5]) LimArr[5] = AuxLimArr[5];
		}
	}
}

//-------------------------------------------------------------------------

void radTg3d::LimitsAtTransform(radTrans* pTr, double* LimArr)
{
	radTVectorOfVector3d LocVertices;
	VerticesInLocFrame(LocVertices, false);

	int AmOfVertices = (int)(LocVertices.size());
	if(AmOfVertices <= 0) return;

	LimArr[0] = 1E+23; LimArr[1] = -1E+23;
	LimArr[2] = 1E+23; LimArr[3] = -1E+23;
	LimArr[4] = 1E+23; LimArr[5] = -1E+23;

	for(int j=0; j<AmOfVertices; j++) 
	{
		TVector3d P = LocVertices[j];
		if(pTr != 0) P = pTr->TrPoint(P);

		if(LimArr[0] > P.x) LimArr[0] = P.x;
		if(LimArr[1] < P.x) LimArr[1] = P.x;
		if(LimArr[2] > P.y) LimArr[2] = P.y;
		if(LimArr[3] < P.y) LimArr[3] = P.y;
		if(LimArr[4] > P.z) LimArr[4] = P.z;
		if(LimArr[5] < P.z) LimArr[5] = P.z;
	}
}

//-------------------------------------------------------------------------
/**
void radTg3d::Limits(double* LimArr)
{
	for(int n=0; n<6; n++) LimArr[n] = 0;

	radTVectorOfVector3d LocVertices;
	VerticesInLocFrame(LocVertices);

	int AmOfVertices = LocVertices.size();
	if(AmOfVertices <= 0) return;

	TVector3d* VertArr = new TVector3d[AmOfVertices];
	for(int i=0; i<AmOfVertices; i++) VertArr[i] = LocVertices[i];

	if(g3dListOfTransform.empty()) 
	{
		DeterminePointsLimits(VertArr, AmOfVertices, LimArr);
		return;
	}
	else
	{
		radTvhg FlatTransforms;
		FlattenSpaceTransforms(FlatTransforms);
		int SizeFlatTransforms = FlatTransforms.size();

		double AuxLimArr[6];
		for(int k=0; k<SizeFlatTransforms; k++)
		{
			radTrans* pCurFlatTrans = (radTrans*)(FlatTransforms[k].rep);
			for(int j=0; j<AmOfVertices; j++) VertArr[j] = pCurFlatTrans->TrPoint(LocVertices[j]);

			DeterminePointsLimits(VertArr, AmOfVertices, AuxLimArr);

			if(LimArr[0] > AuxLimArr[0]) LimArr[0] = AuxLimArr[0];
			if(LimArr[1] < AuxLimArr[1]) LimArr[1] = AuxLimArr[1];
			if(LimArr[2] > AuxLimArr[2]) LimArr[2] = AuxLimArr[2];
			if(LimArr[3] < AuxLimArr[3]) LimArr[3] = AuxLimArr[3];
			if(LimArr[4] > AuxLimArr[4]) LimArr[4] = AuxLimArr[4];
			if(LimArr[5] < AuxLimArr[5]) LimArr[5] = AuxLimArr[5];
		}

	}
	delete[] VertArr;
}
**/
//-------------------------------------------------------------------------

void radTg3d::DeterminePointsLimits(TVector3d* PointsArr, int AmOfPoints, double* LimArr)
{
	if((PointsArr == 0) || (LimArr == 0) || (AmOfPoints <= 0)) return;

	double &xMin = LimArr[0], &xMax = LimArr[1];
	double &yMin = LimArr[2], &yMax = LimArr[3];
	double &zMin = LimArr[4], &zMax = LimArr[5];

	xMin = 1e+23; xMax = -1e+23;
	yMin = 1e+23; yMax = -1e+23;
	zMin = 1e+23; zMax = -1e+23;

	TVector3d* pPoint = PointsArr;
	for(int i=0; i<AmOfPoints; i++)
	{
		if(xMin > pPoint->x) xMin = pPoint->x;
		if(xMax < pPoint->x) xMax = pPoint->x;

		if(yMin > pPoint->y) yMin = pPoint->y;
		if(yMax < pPoint->y) yMax = pPoint->y;

		if(zMin > pPoint->z) zMin = pPoint->z;
		if(zMax < pPoint->z) zMax = pPoint->z;

		pPoint++;
	}
}

//-------------------------------------------------------------------------

int radTg3d::TotalMultiplicity()
{
	int TotMult = 1;
	if(g3dListOfTransform.empty()) return TotMult;

	for(radTlphg::reverse_iterator iter = g3dListOfTransform.rbegin(); iter != g3dListOfTransform.rend(); ++iter)
	{
		radTrans* pTrans = (radTrans*)(((*iter).Handler_g).rep);
		TotMult *= (*iter).m;
	}

	return TotMult;
}

//-------------------------------------------------------------------------

void radTg3d::FlattenSpaceTransforms(radTvhg& FlatTransforms)
{
	if(g3dListOfTransform.empty()) return;

	//radThg ihg(new radIdentTrans());
	radTrans *pOrigIdentTrf = new radTrans(); //OC061007
	pOrigIdentTrf->SetupIdent();
	radThg ihg(pOrigIdentTrf);
	FlatTransforms.push_back(ihg);

	for(radTlphg::reverse_iterator iter = g3dListOfTransform.rbegin(); iter != g3dListOfTransform.rend(); ++iter)
	{
		radTrans* pTrans = (radTrans*)(((*iter).Handler_g).rep);
		int mult = (*iter).m;

		int CurFlatTranSize = (int)(FlatTransforms.size());
		
		if(mult == 1)
		{
			for(int k=0; k<CurFlatTranSize; k++)
			{
				radTrans* pCurFlatTrans = (radTrans*)(FlatTransforms[k].rep);
				*pCurFlatTrans = Product(*pTrans, *pCurFlatTrans); //multiply flat trans from left
			}
			continue;
		}
		
		radTvhg AuxDuplVect;
		for(int k=0; k<CurFlatTranSize; k++)
		{
			radTrans* pCurFlatTrans = (radTrans*)(FlatTransforms[k].rep);
			radThg LocHg(new radTrans(*pCurFlatTrans));
			AuxDuplVect.push_back(LocHg);
		}

		for(int j=1; j<mult; j++)
		{
			for(int k=0; k<CurFlatTranSize; k++)
			{
				radTrans* pCurFlatTrans = (radTrans*)(AuxDuplVect[k].rep);
				*pCurFlatTrans = Product(*pTrans, *pCurFlatTrans); //multiply from left

				radThg LocHg(new radTrans(*pCurFlatTrans));
				FlatTransforms.push_back(LocHg);
			}
		}

		AuxDuplVect.erase(AuxDuplVect.begin(), AuxDuplVect.end());
	}
}

//-------------------------------------------------------------------------

void radTg3d::DumpBinParse_g3d(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{
	//radTlphg g3dListOfTransform;
	int nTrfs = 0;
	inStr >> nTrfs;
	int oldKey = 0, m = 0; 
	for(int iTr=0; iTr<nTrfs; iTr++)
	{
		inStr >> oldKey;
		inStr >> m;
		if(m <= 0) throw 0;

		map<int, int>::const_iterator itKey = mKeysOldNew.find(oldKey);
		int newKey = 0;
		if(itKey == mKeysOldNew.end()) throw 0;

		newKey = itKey->second;
		if(newKey > 0)
		{
			radTmhg::const_iterator iter = gMapOfHandlers.find(newKey);
			if(iter == gMapOfHandlers.end()) throw 0;
			radThg hg = (*iter).second;
			if(radTCast::TransCast(hg.rep)==0) throw 0;
			AddTransform(m, hg);
		}
	}

	//int IsGroupMember;
	inStr >> IsGroupMember;

	//char ConsiderOnlyWithTrans;
	inStr >> ConsiderOnlyWithTrans;

	//TVector3d CentrPoint;
	inStr >> CentrPoint;

	//char MessageChar; 
	inStr >> MessageChar;

	//radTHandleAuxCompDataG3D HandleAuxCompData;
	char swAuxCompData = 0;
	inStr >> swAuxCompData;
	if(swAuxCompData > 0)
	{
		CreateAuxCompData();
		radTAuxCompDataG3D *pAux = HandleAuxCompData.rep;
		if(pAux != 0)
		{
			//radTAuxCompDataG3D
			//TVector3d Force;
			inStr >> pAux->Force;

			//TVector3d Torque;
			inStr >> pAux->Torque;

			//double Energy;
			inStr >> pAux->Energy;

			//char SubdNeedInd;
			inStr >> pAux->SubdNeedInd;
		}
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTg3dRelax::DumpMaterApplied(std::ostream& o) // Porting
{
	o << endl;
	o << "   Material applied: ";
	if(MaterHandle.rep != 0)
	{
		o << endl;
		o << "      Index " << rad.RetrieveElemKey(MaterHandle.rep) << ": ";
		((radTMaterial*)(MaterHandle.rep))->Dump(o, 1);
	}
	else { o << "None";}
}

//-------------------------------------------------------------------------

void radTg3dRelax::DumpBinParse_g3dRelax(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{
	//radThg MaterHandle;
	int oldMatKey = 0;
	inStr >> oldMatKey;
	if(oldMatKey > 0)
	{
		map<int, int>::const_iterator itKey = mKeysOldNew.find(oldMatKey);
		int newMatKey = 0;
		if(itKey == mKeysOldNew.end()) throw 0;
		newMatKey = itKey->second;
		if(newMatKey > 0)
		{
			radTmhg::const_iterator iter = gMapOfHandlers.find(newMatKey);
			if(iter == gMapOfHandlers.end()) throw 0;
			radThg hg = (*iter).second;
			if(radTCast::MaterCast(hg.rep)==0) throw 0;

			MaterHandle = hg;
		}
	}

	//TMatrix3d* pM_LinCoef;
	char swM_LinCoef = 0;
	inStr >> swM_LinCoef;
	if(swM_LinCoef == 0) pM_LinCoef = 0;
	else
	{
		pM_LinCoef = new TMatrix3d();
		inStr >> *pM_LinCoef;
	}

	//TVector3d Magn;
	inStr >> Magn;

	//float AuxFloat1, AuxFloat2, AuxFloat3;
	inStr >> AuxFloat1;
	inStr >> AuxFloat2;
	inStr >> AuxFloat3;
}

//-------------------------------------------------------------------------

void radTg3dRelax::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
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

	if(pFieldKey->J_) return;
	else if(pFieldKey->M_) Pair.V2 = (pTrans == 0)? Magn : pTrans->TrVectField(Magn);
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
//-------------------------------------------------------------------------

radTStructForShapeInt::radTStructForShapeInt(const radTStructForShapeInt& InStr) 
{
	HandleOfSource = InStr.HandleOfSource;
	HandleOfShape = InStr.HandleOfShape;
	IntegrandLength = InStr.IntegrandLength;
	Normal = InStr.Normal;
	VectArray = InStr.VectArray;
	VectTypeArray = InStr.VectTypeArray;
	IntegrandFunPtr = InStr.IntegrandFunPtr;
	IntOverLine_ = InStr.IntOverLine_;
	IntOverSurf_ = InStr.IntOverSurf_;
	IntOverVol_ = InStr.IntOverVol_;
	AbsPrecArray = InStr.AbsPrecArray;
}

//-------------------------------------------------------------------------
