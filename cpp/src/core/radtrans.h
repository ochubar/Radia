
//-------------------------------------------------------------------------
// Definition of radTrans - a class of objects keeping information on
// symmetries.
//-------------------------------------------------------------------------

#ifndef __RADTRANS_H
#define __RADTRANS_H

//#ifndef __GMVECT_H
//#include "gmvect.h"
//#endif
#ifndef __GMTRANS_H
#include "gmtrans.h"
#endif
#ifndef __RADG3D_H
#include "radg3d.h"
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTrans : public gmTrans, public radTg {
protected:
	//TMatrix3d M, M_inv;
	//TVector3d V;
	//double detM, s;

public:
	//int ID_No; //OC060113 (it exists in gmTrans)

	radTrans(const TMatrix3d& InM, const TMatrix3d& InM_inv, const TVector3d& InV, double In_detM, double In_s, int InID_No =-1) 
		: gmTrans(InM, InM_inv, InV, In_detM, In_s, InID_No)
	{
		//M = InM; M_inv = InM_inv; V = InV; s = In_s; detM = In_detM; ID_No = InID_No;
	}
	radTrans(const TMatrix3d& InM, const TVector3d& InV, double In_detM, double In_s, int InID_No =-1)
		: gmTrans(InM, InV, In_detM, In_s, InID_No =-1)
	{
		//M = InM; V = InV; s = In_s; detM = In_detM; M_inv = Matrix3d_inv(InM); ID_No = InID_No;
	}
	radTrans(CAuxBinStrVect& inStr)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Trans(inStr);
	}
	radTrans() : gmTrans() {}

	int Type_g() { return 2;}
	virtual int Type_Trans() { return 0;}

	void Dump(std::ostream& o, int ShortSign) // Porting
	{
		radTg::Dump(o);
		o << "Transformation: ";
		if(ID_No == 1) o << "Translation";
		else if(ID_No == 2) o << "Rotation";
		else if(ID_No == 3) o << "Plane symmetry";
		else if(ID_No == 4) o << "Field inversion";
		else if(ID_No == 10) o << "Composite";

		if(ShortSign) return;

		o << endl;
		o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}

	void DumpBin_Trans(CAuxBinStrVect& oStr)
	{
		//TMatrix3d M, M_inv;
		oStr << M << M_inv;

		//TVector3d V;
		oStr << V;

		//double detM, s;
		oStr << detM << s;

		//int ID_No;
		oStr << ID_No;
	}

	void DumpBinParse_Trans(CAuxBinStrVect& inStr)
	{
		//TMatrix3d M, M_inv;
		inStr >> M;
		inStr >> M_inv;

		//TVector3d V;
		inStr >> V;

		//double detM, s;
		inStr >> detM;
		inStr >> s;

		//int ID_No;
		inStr >> ID_No;
	}

	//void DumpBin(CAuxBinStrVect& oStr, radTmhg& mEl, radThg& hg)
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
	{
		//int newKey = (int)mEl.size() + 1;
		//mEl[newKey] = hg;

		//Start dumping this object
		//oStr << newKey;
		//elemCount++;
		vElemKeysOut.push_back(elemKey);
		oStr << elemKey;

		//Next 5 bytes define/encode element type:
		oStr << (char)Type_g();
		oStr << (char)Type_Trans();
		oStr << (char)0;
		oStr << (char)0;
		oStr << (char)0;

		DumpBin_Trans(oStr);
	}

	//virtual TVector3d TrPoint(const TVector3d& P) { return M*P+V;}
	//virtual TVector3d TrBiPoint(const TVector3d& P) { return M*P;}
	//virtual TVector3d TrVectField(const TVector3d& B) { return s*(M*B);}
	//virtual TVector3d TrVectPoten(const TVector3d& A) { return /* s*detM*(M_inv*A); */ s*detM*(M*A);}
	//virtual TVector3d TrPoint_inv(const TVector3d& P) { return M_inv*(P-V);}
	//virtual TVector3d TrBiPoint_inv(const TVector3d& P) { return M_inv*P;}
	//virtual TVector3d TrVectField_inv(const TVector3d& B) { return s*(M_inv*B);}
	//virtual TVector3d TrVectPoten_inv(const TVector3d& A) { return /* (s/detM)*(M*A); */ (s/detM)*(M_inv*A);}

	//virtual TVector3d TrAxialVect(const TVector3d& A) { return detM*(M*A);}
	//virtual TVector3d TrAxialVect_inv(const TVector3d& A) { return (1./detM)*(M_inv*A);}

	//virtual void TrMatrix(TMatrix3d& Matrix) { Matrix = s*M*Matrix;}
	//virtual void TrMatrix_inv(TMatrix3d& Matrix) { Matrix = s*M_inv*Matrix;}

	//virtual void TrMatrixLeft(TMatrix3d& Matrix) { Matrix = s*Matrix*M;}
	//virtual void TrMatrixLeft_inv(TMatrix3d& Matrix) { Matrix = s*Matrix*M_inv;}

	//int ShowParity() { return int(detM);}

	radTField TrField(const radTField& InField)
	{
		radTField OutField(InField);
		
		if(!InField.FieldKey.Q_) //OC191005
		{
			if(InField.FieldKey.B_) OutField.B = TrVectField(InField.B);
			if(InField.FieldKey.H_) OutField.H = TrVectField(InField.H);
			if(InField.FieldKey.A_) OutField.A = TrVectPoten(InField.A);
		}
		else
		{
			TMatrix3d Q(InField.B, InField.H, InField.A);
			TrMatrix(Q);
			OutField.B = Q.Str0; OutField.H = Q.Str1; OutField.A = Q.Str2;
		}
		
		if(InField.FieldKey.M_) OutField.M = TrVectField(InField.M);
		if(InField.FieldKey.Ib_) OutField.Ib = TrVectField(InField.Ib);
		if(InField.FieldKey.Ih_) OutField.Ih = TrVectField(InField.Ih);
		return OutField;
	}
	radTField TrField_inv(const radTField& InField)
	{
		radTField OutField(InField);
		
		if(!InField.FieldKey.Q_) //OC191005
		{
			if(InField.FieldKey.B_) OutField.B = TrVectField_inv(InField.B);
			if(InField.FieldKey.H_) OutField.H = TrVectField_inv(InField.H);
			if(InField.FieldKey.A_) OutField.A = TrVectPoten_inv(InField.A);
		}
		else
		{
			TMatrix3d Q(InField.B, InField.H, InField.A);
			TrMatrix_inv(Q);
			OutField.B = Q.Str0; OutField.H = Q.Str1; OutField.A = Q.Str2;
		}
		
		if(InField.FieldKey.M_) OutField.M = TrVectField_inv(InField.M);    // Really _inv here ?
		if(InField.FieldKey.Ib_) OutField.Ib = TrVectField_inv(InField.Ib);
		if(InField.FieldKey.Ih_) OutField.Ih = TrVectField_inv(InField.Ih);
		return OutField;
	}
	//void Invert()
	//{
	//	TMatrix3d OldM = M;
	//	M = M_inv; M_inv = OldM;
	//	V = (-1.)*(M*V);
	//	detM = 1./detM;
	//}

	int SizeOfThis() { return sizeof(radTrans);}

	friend radTrans Product(const radTrans& Tr2, const radTrans& Tr1);
	friend void TrProduct(radTrans* Tr2Ptr, radTrans* Tr1Ptr, radTrans& ResTrPtr);

	//void SetupRotation(const TVector3d&, const TVector3d&, double);
	//void SetupTranslation(const TVector3d& InV)
	//{
	//	TVector3d St0(1.,0.,0.), St1(0.,1.,0.), St2(0.,0.,1.);
	//	M = TMatrix3d(St0, St1, St2); M_inv = M; V = InV;
	//	detM = s = 1.;
	//	ID_No = 1;
	//}
	//void SetupIdent()
	//{
	//	TVector3d St0(1.,0.,0.), St1(0.,1.,0.), St2(0.,0.,1.), Zero(0.,0.,0.);
	//	M = TMatrix3d(St0, St1, St2); M_inv = M; 
	//	V = Zero;
	//	detM = s = 1.;
	//	ID_No = 10;
	//}
	//char IsIdent(double RelTol)
	//{
	//	if((fabs(M.Str0.x - 1.) > RelTol) || (fabs(M.Str0.y) > RelTol) || (fabs(M.Str0.z) > RelTol)) return 0;
	//	if((fabs(M.Str1.y - 1.) > RelTol) || (fabs(M.Str1.x) > RelTol) || (fabs(M.Str1.z) > RelTol)) return 0;
	//	if((fabs(M.Str2.z - 1.) > RelTol) || (fabs(M.Str2.x) > RelTol) || (fabs(M.Str2.y) > RelTol)) return 0;
	//	if((fabs(V.x) > RelTol) || (fabs(V.y) > RelTol) || (fabs(V.z) > RelTol)) return 0;
	//	if(fabs(detM - 1.) > RelTol) return 0;
	//	return 1;
	//}
};

//-------------------------------------------------------------------------

inline radTrans Product(const radTrans& Tr2, const radTrans& Tr1)
{
	return radTrans(Tr2.M*Tr1.M, Tr1.M_inv*Tr2.M_inv, Tr2.M*Tr1.V+Tr2.V, Tr2.detM*Tr1.detM, Tr2.s*Tr1.s, 10); // ID_No = 10; // Composite
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radIdentTrans : public radTrans {
public:
	radIdentTrans() 
	{
		M.Str0.x = 1; M.Str0.y = 0; M.Str0.z = 0;
		M.Str1.x = 0; M.Str1.y = 1; M.Str1.z = 0;
		M.Str2.x = 0; M.Str2.y = 0; M.Str2.z = 1;
		M_inv = M; 
		V.x = V.y = V.z = 0; 
		s = 1; 
		detM = 1; 
	}

	TVector3d TrPoint(const TVector3d& P) { return P;}
	TVector3d TrBiPoint(const TVector3d& P) { return P;}
	TVector3d TrVectField(const TVector3d& B) { return B;}
	TVector3d TrVectPoten(const TVector3d& A) { return A;}
	TVector3d TrPoint_inv(const TVector3d& P) { return P;}
	TVector3d TrBiPoint_inv(const TVector3d& P) { return P;}
	TVector3d TrVectField_inv(const TVector3d& B) { return B;}
	TVector3d TrVectPoten_inv(const TVector3d& A) { return A;}
		
	TVector3d TrAxialVect(const TVector3d& A) { return A;}
	TVector3d TrAxialVect_inv(const TVector3d& A) { return A;}

	void TrMatrix(TMatrix3d& Matrix) {}
	void TrMatrix_inv(TMatrix3d& Matrix) {}

	int Type_Trans() { return 1;}
	int SizeOfThis() { return sizeof(radIdentTrans);}

	friend void TrProduct(radTrans* Tr2Ptr, radTrans* Tr1Ptr, radTrans& ResTrPtr);
};

//-------------------------------------------------------------------------

inline void TrProduct(radTrans* Tr2Ptr, radTrans* Tr1Ptr, radTrans& ResTr)
{
	radIdentTrans IdentTr;
	int IdentTrID = IdentTr.Type_Trans();
	int Tr1ID = Tr1Ptr->Type_Trans();
	int Tr2ID = Tr2Ptr->Type_Trans();

	if(Tr2ID==IdentTrID) 
	{
		ResTr = *Tr1Ptr;
	}
	else if(Tr1ID==IdentTrID)
	{
		ResTr = *Tr2Ptr;
	}
	else
	{
		ResTr.M = Tr2Ptr->M*Tr1Ptr->M;
		ResTr.M_inv = Tr1Ptr->M_inv*Tr2Ptr->M_inv;
		ResTr.V = Tr2Ptr->M*Tr1Ptr->V+Tr2Ptr->V;
		ResTr.detM = Tr2Ptr->detM*Tr1Ptr->detM;
		ResTr.s = Tr2Ptr->s*Tr1Ptr->s;
	}

	ResTr.ID_No = 10; // Composite
}

//-------------------------------------------------------------------------

#endif
