/*-------------------------------------------------------------------------
*
* File name:      radmater.h
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

#ifndef __RADMATER_H
#define __RADMATER_H

#include "radsend.h"
#include "gmvect.h"
#include "radg.h"

#include <math.h>
#include <sstream>

//-------------------------------------------------------------------------

class radTMaterial;

//-------------------------------------------------------------------------

class radTMaterialDB {
public:
	char Name[50];

	virtual radTMaterial* SetupMater(double InMr=0) { return 0;}
};

typedef radTHandle<radTMaterialDB> radTHMatDB;
typedef vector <radTHMatDB> radTHMatDBVect;

//-------------------------------------------------------------------------

class radTLinearAnisotropMaterialDB : public radTMaterialDB {
public:
	double KsiPar, KsiPerp, Mr;

	radTLinearAnisotropMaterialDB(const char* InName, double InKsiPar, double InKsiPerp, double InMr)
	{
		SetParam(InName, InKsiPar, InKsiPerp, InMr);
	}
	radTLinearAnisotropMaterialDB() {}

	void SetParam(const char* InName, double InKsiPar, double InKsiPerp, double InMr)
	{
		strcpy(Name, InName);
		KsiPar = InKsiPar; KsiPerp = InKsiPerp; Mr = InMr;
	}

	radTMaterial* SetupMater(double InMr=0); //virtual
};

//-------------------------------------------------------------------------

class radTNonlinearIsotropMaterialDB : public radTMaterialDB {
public:
	double ks[3], ms[3];

	radTNonlinearIsotropMaterialDB(const char* InName, double ms0, double ms1, double ms2, double ks0, double ks1, double ks2)
	{
		SetParam(InName, ms0, ms1, ms2, ks0, ks1, ks2);
	}
	radTNonlinearIsotropMaterialDB() {}

	void SetParam(const char* InName, double ms0, double ms1, double ms2, double ks0, double ks1, double ks2)
	{
		strcpy(Name, InName);
		ks[0] = ks0; ms[0] = ms0;
		ks[1] = ks1; ms[1] = ms1;
		ks[2] = ks2; ms[2] = ms2;
	}

	radTMaterial* SetupMater(double InMr=0); //virtual
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTg3dRelax;

//-------------------------------------------------------------------------

class radTMaterial : public radTg {
protected:
	static radTHMatDBVect MaterDB;

public:
	TVector3d RemMagn; // Don't make it private nor protected
	char EasyAxisDefined;

	radTMaterial(const TVector3d& InRemMagn, char InEasyAxisDefined) 
	{ 
		RemMagn = InRemMagn; EasyAxisDefined = InEasyAxisDefined;
	}
	radTMaterial() { EasyAxisDefined = 0;}

	int Type_g() { return 3;}
	virtual int Type_Material() { return 0;}

	static radTHMatDBVect SetupMaterDB();
	static int GetMaterIndexDB(const char*);
	static radTMaterial* SetupStandardMater(const char*, double Mr=0);

	virtual TVector3d M(const TVector3d& H) { return 0.*H;}
	virtual void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&) {}
	virtual void MultMatrByInstKsiAndMr(const TVector3d&, const TMatrix3d&, TMatrix3d&, TVector3d&) {}

	virtual void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double) {}
	//virtual void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double, radTg3dRelax*) {}
	//virtual void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double, radTg3dRelax*, void* p=0) {} //OC140103

	virtual int FinishSetup(TVector3d&) { return 1;}

	int FinishDuplication(radTMaterial* MatPtr, radThg& hg)
	{
		radTSend Send;
		if(MatPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(MatPtr); hg = hgLoc; return 1;
	}

	void Dump(std::ostream& o, int ShortSign =0) // Porting
	{
		radTg::Dump(o);
		o << "Magnetic material: ";
	}

	void DumpBin_Material(CAuxBinStrVect& oStr)
	{
		//static radTHMatDBVect MaterDB; //no need to dump static members
		//TVector3d RemMagn; // Don't make it private nor protected
		oStr << RemMagn;
		
		//char EasyAxisDefined;
		oStr << EasyAxisDefined;
	}

	void DumpBinParse_Material(CAuxBinStrVect& inStr)
	{
		//TVector3d RemMagn; // Don't make it private nor protected
		inStr >> RemMagn;

		//char EasyAxisDefined;
		inStr >> EasyAxisDefined;
	}

	//static void SteerNewH(TVector3d& PrevH, TVector3d& InstantH, void* pvAuxRelax); //OC140103
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTLinearAnisotropMaterial : public radTMaterial {
	double KsiPar, KsiPerp;
	TMatrix3d KsiTensor;

public:
	radTLinearAnisotropMaterial(const double* InKsiArray, const TVector3d& InRemMagn, char InEasyAxisDefined)
		: radTMaterial(InRemMagn, InEasyAxisDefined)
	{
		KsiPar = InKsiArray[0]; KsiPerp = InKsiArray[1];
		if(InEasyAxisDefined) SetupKsiTensor();
	}

	radTLinearAnisotropMaterial(double InKsiPar, double InKsiPerp, double InMr)
	{
		KsiPar = InKsiPar; KsiPerp = InKsiPerp;
		RemMagn.x = InMr; RemMagn.y = RemMagn.z = 0.;
	}

/*
	radTLinearAnisotropMaterial(int IndDB, const TVector3d& InRemMagn, char InEasyAxisDefined)
		: radTMaterial(InRemMagn, InEasyAxisDefined)
	{
		radTLinearAnisotropMaterialDB *pMat = (radTLinearAnisotropMaterialDB*)(MaterDB[IndDB].rep);

		KsiPar = pMat->KsiPar; KsiPerp = pMat->KsiPerp;
		if(InEasyAxisDefined) SetupKsiTensor();
		else
		{
			RemMagn.x = pMat->Mr;
			RemMagn.y = RemMagn.z = 0.;
		}
	}
*/

	radTLinearAnisotropMaterial(CAuxBinStrVect& inStr) //, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);

		//double KsiPar, KsiPerp;
		inStr >> KsiPar;
		inStr >> KsiPerp;

		//TMatrix3d KsiTensor;
		inStr >> KsiTensor;
	}

	radTLinearAnisotropMaterial() {}

	int Type_Material() { return 1;}
	
	inline void SetupKsiTensor();
	TVector3d M(const TVector3d& H) { return KsiTensor*H + RemMagn;}  // Should we add RemMagn here?
	void DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
	{
		InstantKsiTensor = KsiTensor; InstantMr = RemMagn;
	}
	void MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		MultByKsi = Matr * KsiTensor; MultByMr = Matr * RemMagn;
	}
	//void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2, radTg3dRelax* pMag, void* p=0) //OC140103
	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2)
	{
		TVector3d ESt1(1.,0.,0.), ESt2(0.,1.,0.), ESt3(0.,0.,1.);
		TMatrix3d E(ESt1, ESt2, ESt3);
		TMatrix3d BufMatr = E - Matr*KsiTensor;
		TMatrix3d InvBufMatr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		H = InvBufMatr*(H_Ext + Matr*RemMagn);
	}

	int FinishSetup(TVector3d& Magn)
	{
		if(!EasyAxisDefined)
		{
			double AbsLocMagn = sqrt(Magn.x*Magn.x + Magn.y*Magn.y + Magn.z*Magn.z);
		
			radTSend Send;
			const double AbsTol = 1.E-10;
			if(AbsLocMagn < AbsTol) { Send.ErrorMessage("Radia::Error107"); return 0;}

			RemMagn = (RemMagn.x/AbsLocMagn)*Magn;
			SetupKsiTensor();
			EasyAxisDefined = 1;
		}
		return 1;
	}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTLinearAnisotropMaterial(*this), hg);
	}

	inline void Dump(std::ostream& o, int ShortSign =0); // Porting

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
	//void DumpBin(CAuxBinStrVect& oStr, radTmhg& mEl, radThg& hg)
	{
		//int newKey = (int)mEl.size() + 1;
		//mEl[newKey] = hg;
		//Start dumping this object
		//oStr << newKey;

		vElemKeysOut.push_back(elemKey);
		oStr << elemKey;

		//Next 5 bytes define/encode element type:
		oStr << (char)Type_g();
		oStr << (char)Type_Material();
		oStr << (char)0;
		oStr << (char)0;
		oStr << (char)0;

		//Members of radTMaterial
		DumpBin_Material(oStr);

		//double KsiPar, KsiPerp;
		oStr << KsiPar << KsiPerp;

		//TMatrix3d KsiTensor;
		oStr << KsiTensor;
	}

	int SizeOfThis() { return sizeof(radTLinearAnisotropMaterial);}
};

//-------------------------------------------------------------------------

inline void radTLinearAnisotropMaterial::SetupKsiTensor()
{
	double AbsRemMagn = sqrt(RemMagn.x*RemMagn.x + RemMagn.y*RemMagn.y + RemMagn.z*RemMagn.z);
	TVector3d L = (1./AbsRemMagn)*RemMagn;
	double DeltaKsi = KsiPar-KsiPerp;
	double LxLx, LyLy, LzLz;
	LxLx=L.x*L.x; LyLy=L.y*L.y; LzLz=L.z*L.z;
	TVector3d Str0(KsiPar*LxLx+KsiPerp*(LyLy+LzLz), DeltaKsi*L.x*L.y, DeltaKsi*L.x*L.z);
	TVector3d Str1(Str0.y, KsiPar*LyLy+KsiPerp*(LxLx+LzLz), DeltaKsi*L.y*L.z);
	TVector3d Str2(Str0.z, Str1.z, KsiPar*LzLz+KsiPerp*(LxLx+LyLy));
	KsiTensor.Str0 = Str0; KsiTensor.Str1 = Str1; KsiTensor.Str2 = Str2;
}

//-------------------------------------------------------------------------

inline void radTLinearAnisotropMaterial::Dump(std::ostream& o, int ShortSign) //Porting
{
	radTMaterial::Dump(o);
	o << "Linear anisotropic";

	if(ShortSign==1) return;
	o << endl;
	o << "   {ksipar,ksiper}= {" << KsiPar << ',' << KsiPerp << "}" << endl;

	if(EasyAxisDefined)
		o << "   {mrx,mry,mrz}= {" << RemMagn.x << ',' << RemMagn.y << ',' << RemMagn.z << "}";
	else
		o << "   mr= " << RemMagn.x;

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTLinearIsotropMaterial : public radTMaterial {
	double Ksi;

public:
	radTLinearIsotropMaterial(double InKsi) { Ksi = InKsi;}
	radTLinearIsotropMaterial(const double* InKsiArray, const TVector3d& InRemMagn, char InEasyAxisDefined) 
		: radTMaterial(InRemMagn, InEasyAxisDefined) { Ksi = InKsiArray[0];}
	
	radTLinearIsotropMaterial(CAuxBinStrVect& inStr) //, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);
		//double Ksi;
		inStr << Ksi;
	}

	radTLinearIsotropMaterial() {}

	int Type_Material() { return 2;}

	TVector3d M(const TVector3d& H) { return Ksi * H + RemMagn;} // Should not we add RemMagn here?
	void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
	void MultMatrByInstKsiAndMr(const TVector3d&, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		MultByKsi = Ksi * Matr; MultByMr = Matr * RemMagn;
	}
	//void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2, radTg3dRelax* pMag, void* p=0) //OC140103
	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2) //OC140103
	{
		TVector3d ESt1(1.,0.,0.), ESt2(0.,1.,0.), ESt3(0.,0.,1.);
		TMatrix3d E(ESt1, ESt2, ESt3);
		TMatrix3d BufMatr = E - Ksi*Matr;
		TMatrix3d InvBufMatr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		H = InvBufMatr*(H_Ext + Matr*RemMagn);
	}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTLinearIsotropMaterial(*this), hg);
	}

	void Dump(std::ostream& o, int ShortSign =0) // Porting
	//inline void radTLinearIsotropMaterial::Dump(std::ostream& o, int ShortSign) // Porting
	{
		radTMaterial::Dump(o);
		o << "Linear isotropic";

		if(ShortSign==1) return;
		o << endl;
		o << "   ksi= " << Ksi;

		o << endl;
		o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
	{
		vElemKeysOut.push_back(elemKey);
		oStr << elemKey;

		//Next 5 bytes define/encode element type:
		oStr << (char)Type_g();
		oStr << (char)Type_Material();
		oStr << (char)0;
		oStr << (char)0;
		oStr << (char)0;

		//Members of radTMaterial
		DumpBin_Material(oStr);

		//double Ksi;
		oStr << Ksi;
	}

	int SizeOfThis() { return sizeof(radTLinearIsotropMaterial);}
};

//-------------------------------------------------------------------------

inline void radTLinearIsotropMaterial::DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
{
	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
	TMatrix3d E(E_Str0, E_Str1, E_Str2);
	InstantKsiTensor = Ksi*E; InstantMr = RemMagn;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTNonlinearIsotropMaterial : public radTMaterial {
	double Ms[3], ks[3];
	int lenMs_ks;

	TVector2d* gArrayHM;
	double* gdMdH;
	int gLenArrayHM;

	double gMaxKsi;

public:
	radTNonlinearIsotropMaterial(const double* InMsArray, const double* In_ksArray, int In_lenMs_ks)	
	{ 
		gMaxKsi = 0;
		lenMs_ks = In_lenMs_ks;
		for(int i=0; i<lenMs_ks; i++) { Ms[i]=InMsArray[i]; ks[i]=In_ksArray[i]; gMaxKsi+=ks[i];}

		gArrayHM = 0; gLenArrayHM = 0; gdMdH = 0;
	}
	radTNonlinearIsotropMaterial(TVector2d* InArrayHM, int InLenArrayHM)
	{
		gArrayHM = 0; gdMdH = 0; gLenArrayHM = 0;
		double ZeroTol = 1e-10;
		char PrependZero = 0;
		if((InArrayHM->x > ZeroTol) && (InArrayHM->y > ZeroTol))
		{
			InLenArrayHM++; PrependZero = 1;
		}

		gLenArrayHM = InLenArrayHM;
		AllocateArrays(InLenArrayHM);
		CopyArrayHM(gArrayHM, InArrayHM, InLenArrayHM, PrependZero);
		//Compute_dMdH(gMaxKsi);
		Compute_dMdH(gArrayHM, gdMdH, gLenArrayHM, gMaxKsi);
	}
	radTNonlinearIsotropMaterial(CAuxBinStrVect& inStr) //, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);

		//double Ms[3];
		inStr >> Ms[0]; inStr >> Ms[1]; inStr >> Ms[2];

		//double ks[3];
		inStr >> ks[0]; inStr >> ks[1]; inStr >> ks[2];

		//int lenMs_ks;
		inStr >> lenMs_ks;

		//int gLenArrayHM;
		inStr >> gLenArrayHM;

		//TVector2d* gArrayHM;
		gArrayHM = 0;
		char cTest=0;
		inStr >> cTest;
		if(cTest > 0)
		{
			gArrayHM = new TVector2d[gLenArrayHM];
			if(gArrayHM == 0) throw 0;
			TVector2d *t_gArrayHM = gArrayHM;
			for(int i=0; i<gLenArrayHM; i++) inStr >> (*(t_gArrayHM++));
		}
		//double* gdMdH;
		gdMdH = 0;
		inStr >> cTest;
		if(cTest > 0)
		{
			gdMdH = new double[gLenArrayHM];
			if(gdMdH == 0) throw 0;
			double *t_gdMdH = gdMdH;
			for(int i=0; i<gLenArrayHM; i++) inStr >> (*(t_gdMdH++));
		}

		//double gMaxKsi;
		inStr >> gMaxKsi;
	}

	radTNonlinearIsotropMaterial() { gArrayHM = 0; gLenArrayHM = 0; gdMdH = 0; gMaxKsi = 0;}
	~radTNonlinearIsotropMaterial() { DeallocateArrays();}

	int Type_Material() { return 3;}

	TVector3d M(const TVector3d& H);
	void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
	void MultMatrByInstKsiAndMr(const TVector3d&, const TMatrix3d&, TMatrix3d&, TVector3d&);
	//void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double, radTg3dRelax*, void*); //OC140103
	void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double); //OC140103

	inline void Dump(std::ostream& o, int ShortSign =0);
	inline void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey);

	//void Compute_dMdH(double& MaxKsi);
	static void Compute_dMdH(TVector2d* ArrayHM, double* dMdH, int LenArrayHM, double& MaxKsi);

	//void CheckAndCorrect_dMdH();
    static void CheckAndCorrect_dMdH(TVector2d* ArrayHM, double* dMdH, int LenArrayHM);

	static double Derivative5(TVector2d* f, int PoIndx);
	static double Derivative3(TVector2d* f, int PoIndx);

	//double AbsMvsAbsH_Interpol(double AbsH);
	static double AbsMvsAbsH_Interpol(double AbsH, TVector2d* ArrayHM, double* dMdH, int LenArrayHM);

	//double AbsHvsAbsM_Interpol(double AbsM);
	static double AbsHvsAbsM_Interpol(double AbsM, TVector2d* ArrayHM, double* dMdH, int LenArrayHM);

	//void AbsMvsAbsH_FuncAndDer_Interpol(double AbsH, double& f, double& fDer);
	static void AbsMvsAbsH_FuncAndDer_Interpol(double AbsH, TVector2d* ArrayHM, double* dMdH, int LenArrayHM, double& f, double& fDer);
	void DefineScalarM(double AbsInstantH, double& f, double& InstKsi);
	void DefineScalarM_dMdH(double AbsInstantH, double& f, double& dfdH);

	double FuncNewAbsH(double AbsH, const TMatrix3d& Matr, const TVector3d& H_Ext);
	double FuncToZero(double AbsH, const TMatrix3d& Matr, const TVector3d& H_Ext) 
	{
		return FuncNewAbsH(AbsH, Matr, H_Ext) - AbsH;
	}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{// Add more if new members!
		radTSend Send;
		radTNonlinearIsotropMaterial* pNewMater = 0;
		if((gArrayHM != 0) && (gdMdH != 0) && (gLenArrayHM != 0))
		{
			pNewMater = new radTNonlinearIsotropMaterial();
			if(pNewMater == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

			if(!pNewMater->AllocateArrays(gLenArrayHM)) { Send.ErrorMessage("Radia::Error900"); return 0;}
			CopyArrayHM(pNewMater->gArrayHM, gArrayHM, gLenArrayHM, 0);
			CopyArray_dMdH(pNewMater->gdMdH, gdMdH, gLenArrayHM);
		}
		else pNewMater = new radTNonlinearIsotropMaterial(*this);
		if(pNewMater == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		return FinishDuplication(pNewMater, hg);
	}

	int SizeOfThis() { return sizeof(radTNonlinearIsotropMaterial);}

	void FindNewH_FromKsi(TVector3d& InstantH, const TMatrix3d& Matr, const TVector3d& H_Ext, double gInstKsi) 
	{
		TMatrix3d E; E.Str0.x = 1.; E.Str1.y = 1.; E.Str2.z = 1.;
		TMatrix3d BufMatr = gInstKsi*Matr;
		BufMatr = E - BufMatr;
		TMatrix3d InvBufMatr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		InstantH = InvBufMatr*H_Ext;
	}

	int AllocateArrays(int InLenArrayHM)
	{
		DeallocateArrays();
		gLenArrayHM = InLenArrayHM;
		
		gArrayHM = new TVector2d[gLenArrayHM];
		if(gArrayHM == 0) return 0;

		gdMdH = new double[gLenArrayHM];
		if(gdMdH == 0) return 0;
		return 1;
	}
	void DeallocateArrays()
	{
		if(gArrayHM != 0) delete[] gArrayHM; gArrayHM = 0;
		if(gdMdH != 0) delete[] gdMdH; gdMdH = 0;
		gLenArrayHM = 0;
	}

	static void CopyArrayHM(TVector2d* Dst, TVector2d* Src, int InLenArrayHM, char PrependZero)
	{
		if((Dst == 0) || (Src == 0) || (InLenArrayHM <= 0)) return;

		TVector2d *tArrayHM = Dst, *tInArrayHM = Src;
		if(PrependZero)
		{
			tArrayHM->x = 0.; (tArrayHM++)->y = 0.;
			InLenArrayHM--;
		}
		for(int i=0; i<InLenArrayHM; i++) *(tArrayHM++) = *(tInArrayHM++);
	}
	void CopyArray_dMdH(double* Dst, double* Src, int InLenArrayHM)
	{
		if((Dst == 0) || (Src == 0) || (InLenArrayHM <= 0)) return;

		double *tDst = Dst, *tSrc = Src;
		for(int i=0; i<InLenArrayHM; i++) *(tDst++) = *(tSrc++);
	}
	
	static void CubPln(double Step, double f1, double f2, double fpr1, double fpr2, double* aa)
	{
		double InvStep = 1./Step;
		double f1mf2_d_s1ms2 = (f2 - f1)*InvStep;
		*(aa++) = f1;
		*(aa++) = fpr1;
		*(aa++) = (3.*f1mf2_d_s1ms2 - 2.*fpr1 - fpr2)*InvStep;
		*aa = (-2.*f1mf2_d_s1ms2 + fpr1 + fpr2)*InvStep*InvStep;
	}
};

//-------------------------------------------------------------------------

inline TVector3d radTNonlinearIsotropMaterial::M(const TVector3d& H) 
{
	double AbsH = sqrt(H.x*H.x + H.y*H.y + H.z*H.z);
	double AbsM = 0.;
	if(gLenArrayHM == 0)
	{
		for(int i=0; i<lenMs_ks; i++) 
			if(Ms[i]!=0.) AbsM += Ms[i]*tanh(ks[i]*AbsH/Ms[i]);
	}
	else AbsM = AbsMvsAbsH_Interpol(AbsH, gArrayHM, gdMdH, gLenArrayHM);

	if(AbsH!=0) return (AbsM/AbsH)*H + RemMagn;  // Should not we add RemMagn here?
	else return RemMagn;
}

//-------------------------------------------------------------------------

inline void radTNonlinearIsotropMaterial::MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
{
/**
	double AbsInstantH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
	double Der, f, InstKsi;
	Der = f = InstKsi = 0.;

	if(gLenArrayHM == 0)
	{
		if(AbsInstantH==0.)
		{
			for(int j=0; j<lenMs_ks; j++) Der += ks[j];
			InstKsi = Der;
		}
		else
		{
			for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*AbsInstantH/Ms[i]);
			InstKsi = f/AbsInstantH;
		}
	}
	else
	{
		if(AbsInstantH == 0.) InstKsi = *gdMdH;
		//else InstKsi = AbsMvsAbsH_Interpol(AbsInstantH)/AbsInstantH;
		else InstKsi = AbsMvsAbsH_Interpol(AbsInstantH, gArrayHM, gdMdH, gLenArrayHM)/AbsInstantH;
	}
	MultByKsi = InstKsi*Matr; MultByMr = Matr*RemMagn;
**/

	TMatrix3d InstantKsiTensor;
	TVector3d InstantMr;
	DefineInstantKsiTensor(InstantH, InstantKsiTensor, InstantMr);

	MultByKsi = Matr*InstantKsiTensor; 
	MultByMr = Matr*InstantMr;
}

//-------------------------------------------------------------------------

inline void radTNonlinearIsotropMaterial::Dump(std::ostream& o, int ShortSign) // Porting
{
	radTMaterial::Dump(o);
	o << "Nonlinear isotropic";

	if(ShortSign==1) return;

	o << endl;
	if((gArrayHM == 0) || (gLenArrayHM == 0))
	{
		o << "   {ms1,ms2,ms3}= {" << Ms[0] << ',' << Ms[1] << ',' << Ms[2] << "}" << endl;
		o << "   {ks1,ks2,ks3}= {" << ks[0] << ',' << ks[1] << ',' << ks[2] << "}";
	}
	else
	{
        o << "   M(H) defined by table of values";
	}

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

inline void radTNonlinearIsotropMaterial::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)Type_Material();
	oStr << (char)0;
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTMaterial
	DumpBin_Material(oStr);

	//double Ms[3];
	oStr << Ms[0] << Ms[1] << Ms[2];

	//double ks[3];
	oStr << ks[0] << ks[1] << ks[2];

	//int lenMs_ks;
	oStr << lenMs_ks;

	//int gLenArrayHM;
	oStr << gLenArrayHM;

	//TVector2d* gArrayHM;
	if((gLenArrayHM > 0) && (gArrayHM != 0))
	{
		oStr << (char)1;
		TVector2d *t_gArrayHM = gArrayHM;
		for(int i=0; i<gLenArrayHM; i++) oStr << (*(t_gArrayHM++));
	}
	else oStr << (char)0;

	//double* gdMdH;
	if((gLenArrayHM > 0) && (gdMdH != 0))
	{
		oStr << (char)1;
		double *t_gdMdH = gdMdH;
		for(int i=0; i<gLenArrayHM; i++) oStr << (*(t_gdMdH++));
	}
	else oStr << (char)0;

	//double gMaxKsi;
	oStr << gMaxKsi;
}

//-------------------------------------------------------------------------

#endif
