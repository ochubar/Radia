/*-------------------------------------------------------------------------
*
* File name:      radintrc.h
*
* Project:        RADIA
*
* Description:    Magnetic interaction between "relaxable" field source objects
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADINTRC_H
#define __RADINTRC_H

#include "radcast.h"
#include "gmvectf.h"
//#include "radtrans.h"
#include "gmtrans.h"
#include "radg3d.h"

#include <sstream>
#include <vector>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifdef __GCC__
typedef list <radTPair_int_hg*> radTlphgPtr;
#else
typedef list <radTPair_int_hg*, allocator<radTPair_int_hg*> > radTlphgPtr;
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTRelaxStatusParam {
	double MisfitM, MaxModM, MaxModH;
	radTRelaxStatusParam(double InMisfitM =0., double InMaxModM =0., double InMaxModH =0.) 
	{ 
		MisfitM=InMisfitM; MaxModM=InMaxModM; MaxModH=InMaxModH;
	}
};

//-------------------------------------------------------------------------

enum TRelaxSubIntervalID { RelaxTogether, RelaxApart };

//-------------------------------------------------------------------------

struct radTRelaxSubInterval {
	int StartNo, FinNo;
	TRelaxSubIntervalID SubIntervalID;

	radTRelaxSubInterval(int InStartNo, int InFinNo, TRelaxSubIntervalID InSubIntervalID)
	{
		StartNo = InStartNo; FinNo = InFinNo; SubIntervalID = InSubIntervalID;
	}
	radTRelaxSubInterval() {}

	inline friend int operator <(const radTRelaxSubInterval&, const radTRelaxSubInterval&);
	inline friend int operator ==(const radTRelaxSubInterval&, const radTRelaxSubInterval&);
};

//-------------------------------------------------------------------------

inline int operator <(const radTRelaxSubInterval&, const radTRelaxSubInterval&) { return 1;}

//-------------------------------------------------------------------------

inline int operator ==(const radTRelaxSubInterval& i1, const radTRelaxSubInterval& i2) 
{ 
	return (i1.StartNo == i2.StartNo) && (i1.FinNo == i2.FinNo) && (i1.SubIntervalID == i2.SubIntervalID);
}

//-------------------------------------------------------------------------

#ifdef __GCC__
typedef vector<radTg3dRelax*> radTVectPtrg3dRelax;
typedef vector<radTg3d*> radTVectPtr_g3d;
typedef vector<radTrans*> radTVectPtrTrans;
typedef vector<radTlphgPtr*> radVectPtr_lphgPtr;
typedef vector<radTRelaxSubInterval> radTVectRelaxSubInterval;
#else
typedef vector<radTg3dRelax*, allocator<radTg3dRelax*> > radTVectPtrg3dRelax;
typedef vector<radTg3d*, allocator<radTg3d*> > radTVectPtr_g3d;
typedef vector<radTrans*, allocator<radTrans*> > radTVectPtrTrans;
typedef vector<radTlphgPtr*, allocator<radTlphgPtr*> > radVectPtr_lphgPtr;
typedef vector<radTRelaxSubInterval, allocator<radTRelaxSubInterval> > radTVectRelaxSubInterval;
#endif

#ifdef __MWERKS__
/*
null_template
struct iterator_traits <radTg3dRelax**> {
     typedef ptrdiff_t difference_type;
     typedef radTg3dRelax* value_type;
     typedef radTg3dRelax** pointer;
     typedef radTg3dRelax*& reference;
     typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTg3d**> {
     typedef ptrdiff_t difference_type;
     typedef radTg3d* value_type;
     typedef radTg3d** pointer;
     typedef radTg3d*& reference;
     typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTrans**> {
     typedef ptrdiff_t difference_type;
     typedef radTrans* value_type;
     typedef radTrans** pointer;
     typedef radTrans*& reference;
     typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTlphgPtr**> {
     typedef ptrdiff_t difference_type;
     typedef radTlphgPtr* value_type;
     typedef radTlphgPtr** pointer;
     typedef radTlphgPtr*& reference;
     typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTRelaxSubInterval*> {
     typedef ptrdiff_t difference_type;
     typedef radTRelaxSubInterval value_type;
     typedef radTRelaxSubInterval* pointer;
     typedef radTRelaxSubInterval& reference;
     typedef random_access_iterator_tag iterator_category;
};
*/
#endif

//-------------------------------------------------------------------------

class radTInteraction : public radTg {

	int AmOfMainElem;
	int AmOfExtElem;
	radThg SourceHandle;
	radThg MoreExtSourceHandle;
	radTCompCriterium CompCriterium;
	radTRelaxStatusParam RelaxStatusParam;
	short RelaxationStarted;

	TMatrix3df** InteractMatrix; //OC250504
	//TMatrix3d** InteractMatrix; //OC250504

	TVector3d* ExternFieldArray;
	TVector3d* NewMagnArray;
	TVector3d* NewFieldArray;
	TVector3d* AuxOldMagnArray;
	TVector3d* AuxOldFieldArray;

	radTRelaxSubInterval* RelaxSubIntervArray; // New 
	radTVectPtrg3dRelax g3dRelaxPtrVect;
	radTVectPtr_g3d g3dExternPtrVect;
	radTVectPtrTrans TransPtrVect;
	radVectPtr_lphgPtr IntVectOfPtrToListsOfTransPtr;
	radVectPtr_lphgPtr ExtVectOfPtrToListsOfTransPtr;
	radTVectRelaxSubInterval RelaxSubIntervConstrVect; // New
	radTrans** MainTransPtrArray;

	radTCast Cast;
	radTSend Send;
	radIdentTrans* IdentTransPtr;
	short FillInMainTransOnly;
	char mKeepTransData;

public:

	int AmOfRelaxSubInterv;

	short SomethingIsWrong;
	short MemAllocTotAtOnce;

	radTInteraction(const radThg&, const radThg&, const radTCompCriterium&, short =0, char =0, char =0);
	radTInteraction(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);
	radTInteraction();
	~radTInteraction();

	int Setup(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char AuxOldMagnArrayIsNeeded, char KeepTransData);

	void CountMainRelaxElems(radTg3d*, radTlphgPtr*);
	void AllocateMemory(char ExtraExternFieldArrayIsNeeded);
	void SetupInteractMatrix();
	void SetupExternFieldArray();
	void AddExternFieldFromMoreExtSource();
	void AddMoreExternField(const radThg& hExtraExtSrc);
	//void ZeroAuxOldMagnArray();
	//void StoreAuxOldMagnArray();
	void SubstractOldMagn();
    void AddOldMagn();
	double CalcQuadNewOldMagnDif();
	int CountRelaxElemsWithSym();
	int OutAmOfRelaxObjs() { return AmOfMainElem;}
	void FindMaxModMandH(double& MaxModM, double& MaxModH);

	inline void PushFrontNativeElemTransList(radTg3d*, radTlphgPtr*);
	inline void EmptyVectOfPtrToListsOfTrans();

	inline void FillInTransPtrVectForElem(int, char);
	inline void EmptyTransPtrVect();

	void NestedFor_Trans(radTrans*, const radTlphgPtr::const_iterator&, int, char);
	inline void AddTransOrNestedFor(radTrans*, const radTlphgPtr::const_iterator&, int, char);

	void FillInMainTransPtrArray();
	inline void DestroyMainTransPtrArray();

	void FillInRelaxSubIntervArray(); //New

	int NotEmpty() { return (AmOfMainElem==0)? 0 : 1;}
	inline void Dump(std::ostream&, int =0); // Porting
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey);
	void DumpBinVectOfPtrToListsOfTransPtr(CAuxBinStrVect& oStr, radVectPtr_lphgPtr& VectOfPtrToListsOfTransPtr, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers);
	int DumpBinParseSourceHandle(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers, bool do_g3dCast, bool do_g3dRelaxCast, radThg& out_hg);
	void DumpBinParseVectOfPtrToListsOfTransPtr(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers, radVectPtr_lphgPtr& VectOfPtrToListsOfTransPtr);

	int Type_g() { return 4;}

	inline void ResetM();
	inline void ResetAuxParam();
	inline void InitAuxArrays();

	void ZeroAuxOldArrays(); //OC300504
	inline void StoreAuxOldArrays(); //OC300504
    inline void RestoreAuxOldArrays(); //OC300504

	inline void OutRelaxStatusParam(double*);
	inline void ShowInteractVector(char);
	inline void ShowInteractMatrix();

	inline int SizeOfThis();

	inline void UpdateExternalField();

	friend class radTIterativeRelaxMeth;
	friend class radTSimpleRelaxation;
	friend class radTRelaxationMethNo_2;
	friend class radTRelaxationMethNo_3;
	friend class radTRelaxationMethNo_4;
	friend class radTRelaxationMethNo_a5;
	friend class radTRelaxationMethNo_7;
	friend class radTRelaxationMethNo_8;
};

//-------------------------------------------------------------------------

inline void radTInteraction::PushFrontNativeElemTransList(radTg3d* g3dPtr, radTlphgPtr* ListOfPtrToTransPtr)
{
	for(radTlphg::iterator TrIter = g3dPtr->g3dListOfTransform.begin();	
		TrIter != g3dPtr->g3dListOfTransform.end(); ++TrIter)
		ListOfPtrToTransPtr->push_back(&(*TrIter)); // Improve dereferentiation?
}

//-------------------------------------------------------------------------

inline void radTInteraction::EmptyVectOfPtrToListsOfTrans()
{
	for(unsigned i=1; i<IntVectOfPtrToListsOfTransPtr.size(); i++) 
	{
		radTlphgPtr*& p_lphgPtr = IntVectOfPtrToListsOfTransPtr[i];
		if(p_lphgPtr != 0) delete p_lphgPtr;
		p_lphgPtr = 0;
	}
	IntVectOfPtrToListsOfTransPtr.erase(IntVectOfPtrToListsOfTransPtr.begin(), IntVectOfPtrToListsOfTransPtr.end());
	for(unsigned k=1; k<ExtVectOfPtrToListsOfTransPtr.size(); k++) 
	{
		radTlphgPtr*& p_lphgPtr = ExtVectOfPtrToListsOfTransPtr[k];
		if(p_lphgPtr != 0) delete p_lphgPtr;
		p_lphgPtr = 0;
	}
	ExtVectOfPtrToListsOfTransPtr.erase(ExtVectOfPtrToListsOfTransPtr.begin(), ExtVectOfPtrToListsOfTransPtr.end());
}

//-------------------------------------------------------------------------

inline void radTInteraction::FillInTransPtrVectForElem(int ElemLocInd, char I_or_E)
{
	radTlphgPtr* PtrToListOfPtrToTrans = NULL;
	if(I_or_E == 'I') PtrToListOfPtrToTrans = IntVectOfPtrToListsOfTransPtr[ElemLocInd];
	else PtrToListOfPtrToTrans = ExtVectOfPtrToListsOfTransPtr[ElemLocInd];

	if(PtrToListOfPtrToTrans->empty()) TransPtrVect.push_back(IdentTransPtr);
	else NestedFor_Trans(IdentTransPtr, PtrToListOfPtrToTrans->begin(), ElemLocInd, I_or_E);
}

//-------------------------------------------------------------------------

inline void radTInteraction::EmptyTransPtrVect()
{
	if(Cast.IdentTransCast(TransPtrVect[0])==0) delete TransPtrVect[0];
	for(unsigned i=1; i<TransPtrVect.size(); i++) delete TransPtrVect[i];
	TransPtrVect.erase(TransPtrVect.begin(), TransPtrVect.end());
}

//-------------------------------------------------------------------------

inline void radTInteraction::AddTransOrNestedFor(radTrans* BaseTransPtr, const radTlphgPtr::const_iterator& Iter, int ElemLocInd, char I_or_E)
{
	radTlphgPtr* PtrToListOfPtrToTrans = NULL;
	if(I_or_E == 'I') PtrToListOfPtrToTrans = IntVectOfPtrToListsOfTransPtr[ElemLocInd];
	else PtrToListOfPtrToTrans = ExtVectOfPtrToListsOfTransPtr[ElemLocInd];

	if(Iter == PtrToListOfPtrToTrans->end()) 
	{
		if(Cast.IdentTransCast(BaseTransPtr) == 0) TransPtrVect.push_back(new radTrans(*BaseTransPtr));
		else TransPtrVect.push_back(BaseTransPtr);
	}
	else NestedFor_Trans(BaseTransPtr, Iter, ElemLocInd, I_or_E);
}

//-------------------------------------------------------------------------

inline void radTInteraction::DestroyMainTransPtrArray()
{
	if(MainTransPtrArray == 0) return;

	for(int i=0; i<AmOfMainElem; i++)
	{
		radTrans* MainTransPtr = MainTransPtrArray[i];
		if(MainTransPtr != 0)
		{
			if(Cast.IdentTransCast(MainTransPtr)==0) 
			{
				delete (MainTransPtr);
				MainTransPtr = 0;
			}
		}
	}
	delete[] MainTransPtrArray;
	MainTransPtrArray = 0;
}

//-------------------------------------------------------------------------

inline void radTInteraction::ResetM()
{
	for(int i=0; i<AmOfMainElem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = g3dRelaxPtrVect[i];

		g3dRelaxPtr->Magn = ((radTMaterial*)(g3dRelaxPtrVect[i]->MaterHandle.rep))->RemMagn;
		NewMagnArray[i] = g3dRelaxPtr->Magn;
		NewFieldArray[i] = TVector3d(0.,0.,0.); // Or make it TVector3df
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::ResetAuxParam()
{
	for(int i=0; i<AmOfMainElem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = g3dRelaxPtrVect[i];

		g3dRelaxPtr->AuxFloat1 = 0;
		g3dRelaxPtr->AuxFloat2 = 0;
		g3dRelaxPtr->AuxFloat3 = 0;
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::InitAuxArrays()
{
	for(int i=0; i<AmOfMainElem; i++)
	{
		NewMagnArray[i] = g3dRelaxPtrVect[i]->Magn;
		NewFieldArray[i] = TVector3d(0.,0.,0.); // Or make it TVector3df
	}
}

//-------------------------------------------------------------------------

//inline void radTInteraction::StoreOldMagnData() //OC300504
//{
//	if(AuxOldMagnArray == NULL) return;
//
//    TVector3d *tAuxOldMagnArray = AuxOldMagnArray;
//	for(int i=0; i<AmOfMainElem; i++)
//	{
//		*(tAuxOldMagnArray++) = g3dRelaxPtrVect[i]->Magn;
//	}
//}

//-------------------------------------------------------------------------

inline void radTInteraction::StoreAuxOldArrays()
{
	if(AmOfMainElem <= 0) return;
	
	if((AuxOldMagnArray != NULL) && (AuxOldFieldArray != NULL))
	{
        TVector3d *tAuxOldMagn = AuxOldMagnArray;
		TVector3d *tAuxOldField = AuxOldFieldArray;
		TVector3d *tNewFieldArray = NewFieldArray;

        for(int StNo=0; StNo<AmOfMainElem; StNo++)
		{
			TVector3d &M = (g3dRelaxPtrVect[StNo])->Magn; 
			*(tAuxOldMagn++) = M;
			*(tAuxOldField++) = *(tNewFieldArray++);
		}
	}
	else
	{
		if(AuxOldMagnArray != NULL)
		{
			TVector3d *tAuxOldMagn = AuxOldMagnArray;
			for(int StNo=0; StNo<AmOfMainElem; StNo++)
			{
				TVector3d &M = (g3dRelaxPtrVect[StNo])->Magn; 
				*(tAuxOldMagn++) = M;
			}
		}
		if(AuxOldFieldArray != NULL)
		{
			TVector3d *tAuxOldField = AuxOldFieldArray;
			TVector3d *tNewFieldArray = NewFieldArray;
			for(int StNo=0; StNo<AmOfMainElem; StNo++)
			{
				*(tAuxOldField++) = *(tNewFieldArray++);
			}
		}
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::RestoreAuxOldArrays() //OC300504
{
	if((AuxOldMagnArray == NULL) && (AuxOldFieldArray == NULL)) return;

    TVector3d *tAuxOldMagnArray = AuxOldMagnArray;
	TVector3d *tNewMagnArray = NewMagnArray;
    TVector3d *tAuxOldFieldArray = AuxOldFieldArray;
    TVector3d *tNewFieldArray = NewFieldArray;

	if((AuxOldMagnArray != NULL) && (AuxOldFieldArray != NULL))
	{
		for(int i=0; i<AmOfMainElem; i++)
		{
			g3dRelaxPtrVect[i]->Magn = *tAuxOldMagnArray;
			*(tNewMagnArray++) = *(tAuxOldMagnArray++);
			*(tNewFieldArray++) = *(tAuxOldFieldArray++);
		}
	}
	else
	{
		if(AuxOldMagnArray != NULL)
		{
			for(int i=0; i<AmOfMainElem; i++)
			{
				g3dRelaxPtrVect[i]->Magn = *tAuxOldMagnArray;
				*(tNewMagnArray++) = *(tAuxOldMagnArray++);
			}
		}
		if(AuxOldFieldArray != NULL)
		{
			for(int i=0; i<AmOfMainElem; i++)
			{
                *(tNewFieldArray++) = *(tAuxOldFieldArray++);
			}
		}
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::ShowInteractVector(char Ch)
{
	TVector3d* Vect3dPtr = NULL;
	switch(Ch) 
	{
		case 'E':
			Vect3dPtr = ExternFieldArray; break;
		case 'T':
			Vect3dPtr = NewFieldArray; break;
		case 'M':
			Vect3dPtr = NewMagnArray; break;
		default :
			Vect3dPtr = NewFieldArray; break;
	}
	Send.ArrayOfVector3d(Vect3dPtr, AmOfMainElem);
}

//-------------------------------------------------------------------------

inline void radTInteraction::ShowInteractMatrix()
{
	Send.MatrixOfMatrix3d(InteractMatrix, AmOfMainElem, AmOfMainElem);
}

//-------------------------------------------------------------------------

inline void radTInteraction::Dump(std::ostream& o, int ShortSign) // Porting
{
	radTg::Dump(o);
	o << "Interaction: ";

	if(ShortSign) return;

	o << endl;
	o << "   Number of \"atomic\" relaxable objects: " << AmOfMainElem << endl;
	o << "   Total number of degrees of freedom to relax on: " << AmOfMainElem*3 << endl;
	o << "   Number of external field sources in the general container: " << AmOfExtElem;

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

inline int radTInteraction::SizeOfThis()
{
	long GenSize = sizeof(*this);
	GenSize += AmOfMainElem*AmOfMainElem*sizeof(TMatrix3d);
	GenSize += AmOfMainElem*sizeof(TMatrix3d*);
	GenSize += 3*AmOfMainElem*sizeof(TVector3d);
	GenSize += AmOfMainElem*sizeof(radTg3dRelax*);
	GenSize += AmOfRelaxSubInterv*sizeof(radTRelaxSubInterval);
	return GenSize;
}

//-------------------------------------------------------------------------

inline void radTInteraction::OutRelaxStatusParam(double* RelaxStatusParamArray)
{
	RelaxStatusParamArray[0] = RelaxStatusParam.MisfitM;
	RelaxStatusParamArray[1] = RelaxStatusParam.MaxModM;
	RelaxStatusParamArray[2] = RelaxStatusParam.MaxModH;
	// Add more members of radTRelaxStatusParam here, should they appear in Future
}

//-------------------------------------------------------------------------

inline void radTInteraction::UpdateExternalField()
{
	SetupExternFieldArray(); //zeros and then sets the ExternFieldArray from g3dExternPtrVect
	AddExternFieldFromMoreExtSource(); //adds field to ExternFieldArray from MoreExtSourceHandle
}

//-------------------------------------------------------------------------

#endif
