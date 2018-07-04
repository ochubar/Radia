/*-------------------------------------------------------------------------
*
* File name:      radapl2.cpp
*
* Project:        RADIA
*
* Description:    Wrapping RADIA application function calls
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radappl.h"
#include "radg3da1.h"
#include "radintrc.h"
#include "radmater.h"
#include "radrlmet.h"
#include "radptrj.h"
#include "radopnam.h"
#include "radmtra1.h"
#include "gmvbstr.h"

#include <math.h>
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTApplication::SetLinearMaterial(double* KsiArray, long lenKsiArray, double* RemMagnArray, long lenRemMagnArray)
{
	try
	{
		if(lenKsiArray != 2)
		{
			Send.ErrorMessage("Radia::Error022"); return 0;
		}
		TVector3d RemMagnVect;
		char EasyAxisDefined;
		if(lenRemMagnArray == 3)
		{
			EasyAxisDefined = 1;

			if(!ValidateVector3d(RemMagnArray, lenRemMagnArray, &RemMagnVect)) return 0;
			if((RemMagnVect.x==0) && (RemMagnVect.y==0) && (RemMagnVect.z==0) && (KsiArray[0]!=KsiArray[1]))
			{ Send.ErrorMessage("Radia::Error023"); return 0;}
		}
		else if(lenRemMagnArray == 1) 
		{
			EasyAxisDefined = 0;
			RemMagnVect.x = *RemMagnArray;
		}

		radTMaterial* MaterPtr = new radTLinearAnisotropMaterial(KsiArray, RemMagnVect, EasyAxisDefined);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg hg(MaterPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetMaterialStd(char* MatName, double Mr)
{
	try
	{
		if(MatName == 0)
		{
			Send.ErrorMessage("Radia::Error072"); return 0;
		}

		radTMaterial* MaterPtr = radTMaterial::SetupStandardMater(MatName, Mr);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error073"); return 0;}

		radThg hg(MaterPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetNonlinearIsotropMaterial(double* Ms, long lenMs, double* ks, long len_ks)
{
	try
	{
		if((lenMs != len_ks) || (lenMs > 3)) { Send.ErrorMessage("Radia::Error024"); return 0;}

		radTNonlinearIsotropMaterial* MaterPtr = new radTNonlinearIsotropMaterial(Ms, ks, (int)len_ks);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(MaterPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetNonlinearIsotropMaterial(TVector2d* ArrayHM, int LenArrayHM)
{
	try
	{
		if(!ValidateIsotropMaterDescrByPoints(ArrayHM, LenArrayHM)) return 0;

		radTNonlinearIsotropMaterial* MaterPtr = new radTNonlinearIsotropMaterial(ArrayHM, LenArrayHM);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(MaterPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ValidateIsotropMaterDescrByPoints(TVector2d* ArrayHM, int LenArrayHM)
{
	try
	{
		if((ArrayHM->x < 0.) || (ArrayHM->y < 0.)) { Send.ErrorMessage("Radia::Error071"); return 0;}

		TVector2d* tArrayHM = ArrayHM;
		double Hprev = -1, Mprev = -1;
		for(int i=0; i<LenArrayHM; i++)
		{
			if(tArrayHM->x < Hprev) { Send.ErrorMessage("Radia::Error071"); return 0;}
			if(tArrayHM->y < 0.95*Mprev) { Send.WarningMessage("Radia::Warning014");}

			Hprev = tArrayHM->x; Mprev = tArrayHM->y;
			tArrayHM++;
		}
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetNonlinearLaminatedMaterial(TVector2d* ArrayOfPoints2d, int lenArrayOfPoints2d, double PackFactor, double* dN)
{
	if(lenArrayOfPoints2d <= 0) { Send.ErrorMessage("Radia::Error000"); return 0;}
	if((PackFactor <= 0) || (PackFactor > 1)) { Send.ErrorMessage("Radia::Error074"); return 0;}

	radTMaterial *MaterPtr = 0;

	try
	{
		if(lenArrayOfPoints2d <= 3) 
		{
			double Ms[] = {0,0,0};
			double Ks[] = {0,0,0};
			int lenMs = lenArrayOfPoints2d;
			for(int i=0; i<lenMs; i++) 
			{
				Ks[i] = ArrayOfPoints2d[i].x;
				Ms[i] = ArrayOfPoints2d[i].y;
			}

			if((PackFactor <= 0.) || (PackFactor >= 1.)) MaterPtr = new radTNonlinearIsotropMaterial(Ms, Ks, lenMs);
			else MaterPtr = new radTNonlinearLaminatedMaterial(Ms, Ks, lenMs, PackFactor, dN);
		}
		else
		{
			if(!ValidateIsotropMaterDescrByPoints(ArrayOfPoints2d, lenArrayOfPoints2d)) return 0;

			if((PackFactor <= 0.) || (PackFactor >= 1.)) MaterPtr = new radTNonlinearIsotropMaterial(ArrayOfPoints2d, lenArrayOfPoints2d);
			else MaterPtr = new radTNonlinearLaminatedMaterial(ArrayOfPoints2d, lenArrayOfPoints2d, PackFactor, dN);
		}

		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error000"); return 0;}
		radThg hg(MaterPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Send.ErrorMessage("Radia::Error075"); 
		return 0;
	}
}

//-------------------------------------------------------------------------

//int radTApplication::SetNonlinearAnisotropMaterial(double** Ksi, double** Ms, double* Hc, char* DependenceIsNonlinear)
int radTApplication::SetNonlinearAnisotropMaterial(double** Ksi, double** Ms, double* Hc, int lenHc, char* DependenceIsNonlinear)
{
	try
	{
		radTMaterial* MaterPtr = 0;
		char MaterialIsIsotropic = 1;
		if(DependenceIsNonlinear[0] || DependenceIsNonlinear[1])
		{
			if(!(DependenceIsNonlinear[0] && DependenceIsNonlinear[1])) MaterialIsIsotropic = 0;
			else if((lenHc == 2) && (!(Hc[0]==0. && Hc[1]==0.))) MaterialIsIsotropic = 0;
			else if((lenHc == 4) && ((Hc[0]!=0.) || (Hc[1]!=0.) || (Hc[2]!=0.) || (Hc[3]!=0.))) MaterialIsIsotropic = 0;
			else
			{
				double *KsiPar = Ksi[0], *KsiPer = Ksi[1], *MsPar = Ms[0], *MsPer = Ms[1];
				for(int i=0; i<3; i++)
				{
					if(!((*(KsiPar++)==*(KsiPer++)) && (*(MsPar++)==*(MsPer++)))) { MaterialIsIsotropic = 0; break;}
				}
				if((*KsiPar != 0.) || (*KsiPer != 0.)) MaterialIsIsotropic = 0;
			}
			if(MaterialIsIsotropic) MaterPtr = new radTNonlinearIsotropMaterial(Ms[0], Ksi[0], 3);
			else 
			{
				double Hci[4];
				if(lenHc <= 2) 
				{
					for(int i=0; i<lenHc; i++) Hci[i] = Hc[0];
				}
				else
				{
					for(int i=0; i<lenHc; i++) Hci[i] = Hc[i];
				}
				//MaterPtr = new radTNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
				MaterPtr = new radTNonlinearAnisotropMaterial(Ksi, Ms, Hci, DependenceIsNonlinear);
			}
		}
		else { Send.ErrorMessage("Radia::Error060"); return 0;}

		if(DependenceIsNonlinear[1] && (Hc[1] != 0.) && (lenHc == 2)) { Send.ErrorMessage("Radia::Error061"); return 0;}

		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(MaterPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetNonlinearAnisotropMaterial0(double* pDataPar, int lenDataPar, double* pDataPer, int lenDataPer)
{
	try
	{
		radTMaterial* MaterPtr = 0;
		bool MaterialIsIsotropic = true;
		char DependenceIsNonlinear[] = {(lenDataPar > 1), (lenDataPer > 1)};

		double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hci[4];
		double *tKsiPar = KsiPar, *tKsiPer = KsiPer, *tMsPar = MsPar, *tMsPer = MsPer, *tHci = Hci;
		for(int j=0; j<3; j++)
		{
			*(tKsiPar++) = 0;
            *(tKsiPer++) = 0;
            *(tMsPar++) = 0;
			*(tMsPer++) = 0;
			*(tHci++) = 0;
		}
		*tKsiPar = 0; *tKsiPer = 0; *tHci = 0;

		double *Ksi[] = {KsiPar, KsiPer}, *Ms[] = {MsPar, MsPer};

		if(lenDataPar == 11) //{ksi1,ms1,hc1,ksi2,ms2,hc2,ksi3,ms3,hc3,ksi0,hc0}
		{
			double *tKsi = Ksi[0], *tMs = Ms[0], *tHci = Hci;
			double *tDataPar = pDataPar;
            for(int i=0; i<3; i++) 
			{
				*(tKsi++) = *(tDataPar++);
				*(tMs++) = *(tDataPar++);
				*(tHci++) = *(tDataPar++);
			}
            *tKsi = *(tDataPar++);
            *tHci = *tDataPar;
		}
		else if(lenDataPar == 8) //{ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0,hc}
		{
			double *tKsi = Ksi[0], *tMs = Ms[0];
			double *tDataPar = pDataPar;
            for(int i=0; i<3; i++) 
			{
				*(tKsi++) = *(tDataPar++);
				*(tMs++) = *(tDataPar++);
			}
            *tKsi = *(tDataPar++);
            Hci[0] = Hci[1] = Hci[2] = Hci[3] = *tDataPar;
		}
		else if(lenDataPar == 1) //{ksi0}
		{
			double *tKsi = Ksi[0]; //, *tMs = Ms[0];
			double *tDataPar = pDataPar;
			//for(int i=0; i<3; i++) 
			//{
			//	*(tKsi++) = 0;
			//	*(tMs++) = 0;
			//}
            *tKsi = *tDataPar;
            Hci[0] = Hci[1] = Hci[2] = Hci[3] = 0;
		}

		if(lenDataPer == 7) //{ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0}
		{
			double *tKsi = Ksi[1], *tMs = Ms[1];
			double *tDataPer = pDataPer;
            for(int i=0; i<3; i++) 
			{
				*(tKsi++) = *(tDataPer++);
				*(tMs++) = *(tDataPer++);
			}
            *tKsi = *tDataPer;
		}
		else if(lenDataPer == 1) //{ksi0}
		{
			double *tKsi = Ksi[1]; //, *tMs = Ms[1];
			double *tDataPer = pDataPer;
			//for(int i=0; i<3; i++) 
			//{
			//	*(tKsi++) = 0;
			//	*(tMs++) = 0;
			//}
            *tKsi = *tDataPer;
		}

		if(DependenceIsNonlinear[0] || DependenceIsNonlinear[1])
		{
			if(!(DependenceIsNonlinear[0] && DependenceIsNonlinear[1])) MaterialIsIsotropic = false;
			else if(!(Hci[0]==0. && Hci[1]==0. && Hci[2]==0. && Hci[3]==0.)) MaterialIsIsotropic = false;
			else
			{
				double *KsiPar = Ksi[0], *KsiPer = Ksi[1], *MsPar = Ms[0], *MsPer = Ms[1];
				for(int i=0; i<3; i++)
				{
					if(!((*(KsiPar++)==*(KsiPer++)) && (*(MsPar++)==*(MsPer++)))) { MaterialIsIsotropic = false; break;}
				}
				if((*KsiPar != 0.) || (*KsiPer != 0.)) MaterialIsIsotropic = false;
			}
			if(MaterialIsIsotropic) MaterPtr = new radTNonlinearIsotropMaterial(Ms[0], Ksi[0], 3);
			else 
			{
				MaterPtr = new radTNonlinearAnisotropMaterial(Ksi, Ms, Hci, DependenceIsNonlinear);
			}
		}
		else { Send.ErrorMessage("Radia::Error060"); return 0;}

		//if(DependenceIsNonlinear[1] && (Hc[1] != 0.)) { Send.ErrorMessage("Radia::Error061"); return 0;}

		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(MaterPtr);
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ApplyMaterial(int g3dElemKey, int MaterElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(g3dElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		radTg3dRelax* g3dRelaxPtr = Cast.g3dRelaxCast(g3dPtr); 
		if(g3dRelaxPtr==0) 
		{
			radTGroup* GroupPtr = Cast.GroupCast(g3dPtr);
			if(GroupPtr==0) { Send.ErrorMessage("Radia::Error015"); return 0;}
		}

		if(!ValidateElemKey(MaterElemKey, hg)) return 0;
		radTMaterial* MaterPtr = Cast.MaterCast(hg.rep);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error016"); return 0;}

		char PutNewStuffIntoGenCont = 1;
		if(!MaterPtr->EasyAxisDefined) 
		{
			if(!MaterPtr->DuplicateItself(hg, this, PutNewStuffIntoGenCont)) return 0;
			AddElementToContainer(hg); // Maybe not necessary
		}

		if(!g3dPtr->SetMaterial(hg, this)) return 0;

		if(SendingIsRequired) Send.Int(g3dElemKey);
		return g3dElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------
/**
void radTApplication::OutMagnetizCompRes(char* MagnChar, TVector3d& M_vect)
{
	char* BufChar = MagnChar;
	char* EqEmptyStr = "MxMyMz";

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			if((((*BufChar == 'M') || (*BufChar == 'm')) && 
			   (*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) ||
			   (*BufChar == 'X') || (*BufChar == 'x') ||
			   (*BufChar == 'Y') || (*BufChar == 'y') ||
			   (*BufChar == 'Z') || (*BufChar == 'z')) ItemCount++;
			BufChar++;
		}
		BufChar = MagnChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}

	if(ItemCount > 1) Send.InitOutList(ItemCount);

	while (*BufChar != '\0') 
	{
		if((*(BufChar)=='M') || (*(BufChar)=='m'))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) Send.Vector3d(&M_vect);
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) Send.Double(M_vect.x);
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) Send.Double(M_vect.y);
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) Send.Double(M_vect.z);

		BufChar++;
	}
}
**/
//-------------------------------------------------------------------------

void radTApplication::ComputeMvsH(int g3dRelaxOrMaterElemKey, char* MagnChar, double* H, long lenH)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(g3dRelaxOrMaterElemKey, hg)) return;

		radTMaterial* MaterPtr = NULL;

		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr!=NULL)
		{
			radTg3dRelax* g3dRelaxPtr = Cast.g3dRelaxCast(g3dPtr);
			if(g3dRelaxPtr!=NULL) 
			{
				MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
				if(MaterPtr==NULL) { Send.ErrorMessage("Radia::Error027"); return;}
			}
			else
			{
				radTGroup* GroupPtr = Cast.GroupCast(g3dPtr);
				if(GroupPtr!=NULL)
				{
					radTg3dRelax* g3dSubdRelaxPtr = NULL;

					radTSubdividedRecMag* SubdividedRecMagPtr = Cast.SubdividedRecMagCast(GroupPtr);
					if(SubdividedRecMagPtr!=NULL) g3dSubdRelaxPtr = (radTg3dRelax*)SubdividedRecMagPtr;
					else
					{
						radTSubdividedExtrPolygon* SubdividedExtrPolygonPtr = Cast.SubdExtrPolygonCastFromGroup(GroupPtr);
						if(SubdividedExtrPolygonPtr!=NULL) g3dSubdRelaxPtr = (radTg3dRelax*)SubdividedExtrPolygonPtr;
						else
						{
							radTSubdividedPolyhedron* SubdividedPolyhedronPtr = Cast.SubdPolyhedronCastFromGroup(GroupPtr);
							if(SubdividedPolyhedronPtr!=NULL) g3dSubdRelaxPtr = (radTg3dRelax*)SubdividedPolyhedronPtr;
						}
					}
					if(g3dSubdRelaxPtr!=NULL) 
					{
						MaterPtr = (radTMaterial*)(g3dSubdRelaxPtr->MaterHandle.rep);
						if(MaterPtr==NULL) { Send.ErrorMessage("Radia::Error027"); return;}
					}
				}
			}
		}
		if(MaterPtr==NULL)
		{
			MaterPtr = Cast.MaterCast(hg.rep);
			if(MaterPtr==NULL) { Send.ErrorMessage("Radia::Error025"); return;}
		}

		if(!ValidateMagnChar(MagnChar)) return;
		TVector3d H_vect;
		if(!ValidateVector3d(H, lenH, &H_vect)) return;

		TVector3d M_vect = MaterPtr->M(H_vect);
		if(SendingIsRequired) Send.OutMagnetizCompRes(MagnChar, M_vect);
	}
	catch(...)
	{
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//void radTApplication::DumpElem(int ElemKey)
void radTApplication::DumpElem(int* arKeys, int nElem, const char* strFormat, bool arKeysAllocInMathLink)
{
	try
	{
		if((strFormat == 0) || (strcmp(strFormat, "asc") == 0) || (strcmp(strFormat, "ascii") == 0) || (strcmp(strFormat, "ASC") == 0) || (strcmp(strFormat, "ASCII") == 0))
		{
			ostringstream OutDumpStream;
			int nElem_mi_1 = nElem - 1;
			for(int i=0; i<nElem; i++)
			{
				int elemKey = arKeys[i];
				radThg hg;
				if(!ValidateElemKey(elemKey, hg)) return;
				
				OutDumpStream << "Index " << elemKey << ": ";
				(hg.rep)->Dump(OutDumpStream);

				if(i < nElem_mi_1) OutDumpStream << endl;
			}
			
			OutDumpStream << ends;
			Send.String(OutDumpStream.str().c_str());
		}
		else if((strcmp(strFormat, "bin") == 0) || (strcmp(strFormat, "binary") == 0) || (strcmp(strFormat, "BIN") == 0) || (strcmp(strFormat, "BINARY") == 0))
		{
			//CAuxBinStr oStr;
			CAuxBinStrVect oStr;

			oStr << (char)arKeysAllocInMathLink; //indicates whether at parsing a list of elements should be expected

			oStr << nElem; //number of "directly listed" elements, to start with
			for(int j=0; j<nElem; j++) //first writing keys of directly-dumped elements 
			{
				oStr << arKeys[j];
			}
			oStr << nElem; //number of elements (it will be corrected later with actual number of elements)

			//radTmhg locElMap;
			//int elemCount=0;
			vector<int> vElemKeysOut;
			for(int i=0; i<nElem; i++)
			{
				int elemKey = arKeys[i];

				int indExist = CAuxParse::FindElemInd(elemKey, vElemKeysOut);
				if(indExist < 0)
				{//to avoid duplication of objects in output byte string
					radThg hg;
					if(!ValidateElemKey(elemKey, hg)) return;
					//(hg.rep)->DumpBin(oStr, locElMap, hg);
					(hg.rep)->DumpBin(oStr, vElemKeysOut, GlobalMapOfHandlers, GlobalUniqueMapKey, elemKey);
				}
			}

			//int elemCount = (int)locElMap.size();
			//oStr.setFromPos(0, elemCount);
			int elemCount = (int)vElemKeysOut.size();
			//oStr.setFromPos((long)((nElem + 1)*(sizeof(int))), elemCount);
			oStr.setFromPos((long)((nElem + 1)*(sizeof(int)) + 1), elemCount); //OC060713

			//Saving "Drawing Attributes" of objects
			long drwAttrOfst = oStr.getCurOfst();
			int nDrwAttrFound = 0;
			oStr << nDrwAttrFound;
			for(int j=0; j<elemCount; j++)
			{
				int elemKey = vElemKeysOut[j];
				radTMapOfDrawAttr::const_iterator itDrw = MapOfDrawAttr.find(elemKey);
				if(itDrw != MapOfDrawAttr.end())
				{
					const radTDrawAttr &drwAttr = itDrw->second;
					oStr << elemKey;
					//Members of radTDrawAttr
					//double Red, Green, Blue; 
					oStr << drwAttr.RGB_col.Red << drwAttr.RGB_col.Green << drwAttr.RGB_col.Blue; 
					//double LineThickness;
					oStr << drwAttr.LineThickness;
					nDrwAttrFound++;
				}
			}
			if(nDrwAttrFound > 0) oStr.setFromPos(drwAttrOfst, nDrwAttrFound);

			Send.ByteString(reinterpret_cast<const unsigned char*>(oStr.data()), (long)oStr.size());

			//DEBUG
			//CAuxBinStrVect inStr(reinterpret_cast<const unsigned char*>(oStr.data()), (long)oStr.size());
			//int i0, i1, i2, i3, i4, i5, i6;
			//inStr >> i0; inStr >> i1; inStr >> i2; inStr >> i3; inStr >> i4; inStr >> i5; inStr >> i6;
			//END DEBUG
		}
		else 
		{
			Send.ErrorMessage("Radia::Error000");
		}
	}
	catch(...)
	{
		Initialize();
	}
}

//-------------------------------------------------------------------------

int radTApplication::DumpElemParse(const unsigned char *bstr, int bstrLen)
{
	if((bstr == 0) || (bstrLen <= 0)) { Send.ErrorMessage("Radia::Error000"); return 0;}
	int *arDirElemOldKeys = 0;
	try
	{
		CAuxBinStrVect inStr(bstr, bstrLen);

		char listIsExpectedInOutput = 0;
		inStr >> listIsExpectedInOutput; //indicates whether at parsing a list of elements should be expected

		int nElemDir = 0;
		inStr >> nElemDir;

		arDirElemOldKeys = new int[nElemDir];
		int *t_arDirElemOldKeys = arDirElemOldKeys;
		for(int j=0; j<nElemDir; j++)
		{
			inStr >> (*(t_arDirElemOldKeys++));
		}

		int nElemTot = 0;
		inStr >> nElemTot;

		int oldKey;
		char cType1, cType2, cType3, cType4, cType5;
		map<int, int> mKeysOldNew;
		vector<int> vElemKeys;
		
		radTg3d g3d;
		radTArcCur arcCur;
		radTFlmLinCur flmCur;
		radTBackgroundFieldSource bkgFldSrc;
		radTGroup grp;
		radTg3dRelax g3dRelax;
		radTRecMag recMag;
		radTSubdividedRecMag sbdRecMag;
		radTExtrPolygon extPgn;
		radTSubdividedExtrPolygon sbdExtPgn;
		radTPolyhedron polyhdr;
		radTSubdividedPolyhedron sbdPolyhdr;
		radTrans tr;
		radTMaterial mat;
		radTLinearAnisotropMaterial matLinAniso;
		radTLinearIsotropMaterial matLinIso;
		radTNonlinearIsotropMaterial matNonLinIso;
		radTNonlinearAnisotropMaterial matNonLinAniso;
		radTNonlinearLaminatedMaterial matNonLinLam;
		radTInteraction intrc;

		for(int i=0; i<nElemTot; i++)
		{
			inStr >> oldKey;

			inStr >> cType1;
			inStr >> cType2;
			inStr >> cType3;
			inStr >> cType4;
			inStr >> cType5;

			radThg hg;

			if(cType1 == g3d.Type_g())
			{
				if(cType2 == g3dRelax.Type_g3d())
				{
					if(cType3 == recMag.Type_g3dRelax())
					{//Instantiate RecMag
						hg = radThg(new radTRecMag(inStr, mKeysOldNew, GlobalMapOfHandlers));
					}
					else if(cType3 == extPgn.Type_g3dRelax())
					{//Instantiate ExtrPolygon
						hg = radThg(new radTExtrPolygon(inStr, mKeysOldNew, GlobalMapOfHandlers));
					}
					else if(cType3 == polyhdr.Type_g3dRelax())
					{//Instantiate Polyhedron
						hg = radThg(new radTPolyhedron(inStr, mKeysOldNew, GlobalMapOfHandlers));
					}
				}
				else if(cType2 == grp.Type_g3d())
				{//Instantiate Group
					if(cType3 == grp.Type_Group())
					{//Instantiate Group
						hg = radThg(new radTGroup(inStr, mKeysOldNew, GlobalMapOfHandlers));
					}
					else if(cType3 == sbdRecMag.Type_Group())
					{
						hg = radThg((radTGroup*)(new radTSubdividedRecMag(inStr, mKeysOldNew, GlobalMapOfHandlers)));
					}
					else if(cType3 == sbdExtPgn.Type_Group())
					{//Instantiate Subdivided ExtrPolygon
						hg = radThg((radTGroup*)(new radTSubdividedExtrPolygon(inStr, mKeysOldNew, GlobalMapOfHandlers)));
					}
					else if(cType3 == sbdPolyhdr.Type_Group())
					{//Instantiate Subdivided Polyhedron
						hg = radThg((radTGroup*)(new radTSubdividedPolyhedron(inStr, mKeysOldNew, GlobalMapOfHandlers)));
					}
				}
				else if(cType2 == arcCur.Type_g3d())
				{//Instantiate ArcCur
					hg = radThg(new radTArcCur(inStr, mKeysOldNew, GlobalMapOfHandlers));
				}
				else if(cType2 == flmCur.Type_g3d())
				{//Instantiate FlmLinCur
					hg = radThg(new radTFlmLinCur(inStr, mKeysOldNew, GlobalMapOfHandlers));
				}
				else if(cType2 == bkgFldSrc.Type_g3d())
				{//Instantiate BackgroundFieldSource
					hg = radThg(new radTBackgroundFieldSource(inStr, mKeysOldNew, GlobalMapOfHandlers));
				}
			}
			else if(cType1 == tr.Type_g())
			{//Instantiate Transformation
				hg = radThg(new radTrans(inStr));
			}
			else if(cType1 == mat.Type_g())
			{
				if(cType2 == matLinIso.Type_Material())
				{//Instantiate Linear Isotropic Material
					hg = radThg(new radTLinearIsotropMaterial(inStr));
				}
				else if(cType2 == matLinAniso.Type_Material())
				{//Instantiate Linear Anisotropic Material
					hg = radThg(new radTLinearAnisotropMaterial(inStr));
				}
				else if(cType2 == matNonLinIso.Type_Material())
				{//Instantiate Non-Linear Isotropic Material
					hg = radThg(new radTNonlinearIsotropMaterial(inStr));
				}
				else if(cType2 == matNonLinAniso.Type_Material())
				{//Instantiate Non-Linear Anisotropic Material
					if(cType3 == matNonLinAniso.Type_NonlinearAnisotropMaterial())
					{
						hg = radThg(new radTNonlinearAnisotropMaterial(inStr));
					}
					else if(cType3 == matNonLinLam.Type_NonlinearAnisotropMaterial())
					{
						hg = radThg(new radTNonlinearLaminatedMaterial(inStr));
					}
				}
			}
			else if(cType1 == intrc.Type_g())
			{//Instantiate Interaction Matrix
				hg = radThg(new radTInteraction(inStr, mKeysOldNew, GlobalMapOfHandlers));
			}

			int elemKey = AddElementToContainer(hg);
			vElemKeys.push_back(elemKey);
			mKeysOldNew[oldKey] = elemKey;
		}

		//Drawing Attributes
		int nDrwAttrFound = 0;
		inStr >> nDrwAttrFound;
		//double red, green, blue, lineThick;
		for(int id=0; id<nDrwAttrFound; id++)
		{
			int oldElemKey=0;
			inStr >> oldElemKey;

			radTDrawAttr DrawAttr;
			inStr >> DrawAttr.RGB_col.Red;
			inStr >> DrawAttr.RGB_col.Green;
			inStr >> DrawAttr.RGB_col.Blue;
			inStr >> DrawAttr.LineThickness;

			int newElemKey=0;
			map<int, int>::const_iterator itOldNewKey = mKeysOldNew.find(oldElemKey);
			if(itOldNewKey != mKeysOldNew.end())
			{
				MapOfDrawAttr[itOldNewKey->second] = DrawAttr;
			}
		}

		int res = 0;
		//int trueNumElems = (int)(vElemKeys.size());
		//if(trueNumElems == 1)
		//if(nElemDir == 1)
		if((nElemDir == 1) && (!listIsExpectedInOutput))
		{
			//int elemKey = vElemKeys[0];
			map<int, int>::const_iterator itKeyOldNew = mKeysOldNew.find(*arDirElemOldKeys);
			int elemKey = 0;
			if(itKeyOldNew != mKeysOldNew.end())
			{
				elemKey = itKeyOldNew->second;
			}
			if(SendingIsRequired) Send.Int(elemKey);
			//return elemKey;
			res = elemKey;
		}
		//else if(trueNumElems > 1)
		//else if(nElemDir > 1)
		else if(listIsExpectedInOutput)
		{
			if(SendingIsRequired)
			{
				for(int j=0; j<nElemDir; j++)
				{
					int oldElemKey = arDirElemOldKeys[j];
					map<int, int>::const_iterator itKeyOldNew = mKeysOldNew.find(oldElemKey);
					int elemKey = 0;
					if(itKeyOldNew != mKeysOldNew.end())
					{
						elemKey = itKeyOldNew->second;
					}
					arDirElemOldKeys[j] = elemKey;
				}
				Send.IntList(arDirElemOldKeys, nElemDir);
			}
			//return 1;
			res = 1;
		}

		if(arDirElemOldKeys != 0) delete[] arDirElemOldKeys;
		return res;
	}
	catch(...)
	{
		Send.ErrorMessage("Radia::Error202");
		Initialize(); 
		if(arDirElemOldKeys != 0) delete[] arDirElemOldKeys;
		return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::RetrieveElemKey(const radTg* IngPtr)
{
	try
	{
		int ElemKey = 0;
		for(radTmhg::iterator GenIter = GlobalMapOfHandlers.begin();
			GenIter != GlobalMapOfHandlers.end(); ++GenIter)
			if(((*GenIter).second).rep == IngPtr) { ElemKey = (*GenIter).first; break;}

			return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::GenDump()
{
	try
	{
//#ifdef __GCC__
//		ostrstream OutDumpStream;
//#else
		ostringstream OutDumpStream; // Porting
//#endif

		OutDumpStream << "rad: Currently in memory:\n";
		int AmOfElem = (int)(GlobalMapOfHandlers.size());
		if(AmOfElem > 0)
		{
			for(radTmhg::const_iterator iter = GlobalMapOfHandlers.begin();
				iter != GlobalMapOfHandlers.end(); ++iter)
			{
				OutDumpStream << "  Elem. No.:" << (*iter).first << endl;
				(((*iter).second).rep)->Dump(OutDumpStream, 1);
			}
			OutDumpStream << ends;

//#ifdef __GCC__
//			Send.String(OutDumpStream.str());
//#else
			Send.String(OutDumpStream.str().c_str()); // Porting
//#endif
		}
		else Send.ErrorMessage("Radia::Error100");
	}
	catch(...)
	{
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTApplication::ApplyDrawAttrToElem_g3d(int ElemKey, double* RGB_col, long lenRGB_col, double InLineThickness)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
		TVector3d RGB_colVect;
		if(!ValidateVector3d(RGB_col, lenRGB_col, &RGB_colVect)) return 0;
		// May be not necessary?
		if((RGB_col[0]<0.) || (RGB_col[1]<0.) || (RGB_col[2]<0.) || 
		   (RGB_col[0]>1.) || (RGB_col[1]>1.) || (RGB_col[2]>1.))
		{
			Send.ErrorMessage("Radia::Error008"); return 0;
		}

		radTMapOfDrawAttr::iterator iter = MapOfDrawAttr.find(ElemKey);
		if(!(iter == MapOfDrawAttr.end())) MapOfDrawAttr.erase(iter);

		radRGB ColRGB(RGB_colVect.x, RGB_colVect.y, RGB_colVect.z);
		radTDrawAttr DrawAttr;
		DrawAttr.RGB_col = ColRGB;
		DrawAttr.LineThickness = (InLineThickness<0)? 0.001 : InLineThickness;

		MapOfDrawAttr[ElemKey] = DrawAttr;

		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::RemoveDrawAttrFromElem_g3d(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		radTMapOfDrawAttr::iterator iter = MapOfDrawAttr.find(ElemKey);
		if(iter == MapOfDrawAttr.end()) 
		{
			Send.ErrorMessage("Radia::Error013");
			return 0;
		}
		else
		{
			MapOfDrawAttr.erase(iter);
			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::GraphicsForElem_g3d(int ElemKey, int InShowSymmetryChilds, const char** arOptionNames, const char** arOptionValues, int numOptions)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		radTOptionNames OptNam;
		const char* OptNamesToFind[] = {OptNam.Debug};
		char OptValsFoundParsed[] = {0};
		char &doDebug = OptValsFoundParsed[0]; // 0- No; 1- Yes;
		if(!OptNam.findParseOptionValues(arOptionNames, arOptionValues, numOptions, OptNamesToFind, 1, OptValsFoundParsed, 0, 0))
		{
			Send.ErrorMessage("Radia::Error062"); return 0;
		}

		radGraphPresOptions InGraphPresOptions((char)InShowSymmetryChilds, doDebug);

		Send.GenInitDraw();

		radTg3dGraphPresent* g3dGraphPresentPtr = g3dPtr->CreateGraphPresent();

		g3dGraphPresentPtr->SetGraphPresOptions(InGraphPresOptions);
		g3dGraphPresentPtr->MapOfDrawAttrPtr = &MapOfDrawAttr;
		g3dGraphPresentPtr->RetrieveDrawAttr(ElemKey);

		g3dGraphPresentPtr->GenDraw();

		delete g3dGraphPresentPtr;
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::GraphicsForAll_g3d(int InShowSymmetryChilds)
{
	try
	{
		int TotalElem = (int)(GlobalMapOfHandlers.size());
		radTg3d** g3dPtrPtr = new radTg3d*[TotalElem];
		int* KeyPtr = new int[TotalElem];

		radGraphPresOptions InGraphPresOptions((char)InShowSymmetryChilds);

		int g3dPresElemCount = 0;
		for(radTmhg::const_iterator iter = GlobalMapOfHandlers.begin();
			iter != GlobalMapOfHandlers.end(); ++iter)
		{
			radTg* gPtr = ((*iter).second).rep;
			radTg3d g3d;
			if(gPtr->Type_g() == g3d.Type_g())
				if(!((radTg3d*)gPtr)->IsGroupMember) 
				{
					g3dPtrPtr[g3dPresElemCount] = (radTg3d*)gPtr;
					KeyPtr[g3dPresElemCount++] = (*iter).first;
				}
		}
		if(g3dPresElemCount != 0) 
		{
			Send.GenInitDraw();

			Send.InitOutList(g3dPresElemCount);
			for(int i = 0; i < g3dPresElemCount; i++)
			{
				radTg3dGraphPresent* g3dGraphPresentPtr = g3dPtrPtr[i]->CreateGraphPresent();

				g3dGraphPresentPtr->SetGraphPresOptions(InGraphPresOptions);
				g3dGraphPresentPtr->MapOfDrawAttrPtr = &MapOfDrawAttr;
				g3dGraphPresentPtr->RetrieveDrawAttr(KeyPtr[i]);
				g3dGraphPresentPtr->GenDraw();
				delete g3dGraphPresentPtr;
			}
		}
		else Send.ErrorMessage("Radia::Error101");
		delete[] g3dPtrPtr;
	}
	catch(...)
	{
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTApplication::PreRelax(int ElemKey, int SrcElemKey)
{
	radThg hg;
	if(!ValidateElemKey(ElemKey, hg)) return 0;
	radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
	if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

	radThg hgMoreExtSrc;
	if(SrcElemKey!=0)
	{
		if(!ValidateElemKey(SrcElemKey, hgMoreExtSrc)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hgMoreExtSrc.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
	}

	try
	{
		char AllocateExtraArray = 1; //OC300504
		char KeepTransData = 1; //OC240408 to enable update after scaling of currents

		radTInteraction* InteractionPtr = new radTInteraction(hg, hgMoreExtSrc, CompCriterium, MemAllocForIntrctMatrTotAtOnce, AllocateExtraArray, KeepTransData);

		if(InteractionPtr->SomethingIsWrong) 
		{ 
			delete InteractionPtr; 
			return 0;
		} // The message has already been sent
		else if(!(InteractionPtr->NotEmpty())) { delete InteractionPtr; Send.ErrorMessage("Radia::Error102"); return 0;}
		else
		{
			radThg InteractHandle(InteractionPtr);
			int InteractElemKey = AddElementToContainer(InteractHandle);
			if(SendingIsRequired) Send.Int(InteractElemKey);
			return InteractElemKey;
		}
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ShowInteractMatrix(int InteractElemKey)
{
	radThg hg;
	if(!ValidateElemKey(InteractElemKey, hg)) return;
	radTInteraction* InteractPtr = Cast.InteractCast(hg.rep); 
	if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return;}

	InteractPtr->ShowInteractMatrix();
}

//-------------------------------------------------------------------------

void radTApplication::ShowInteractVector(int InteractElemKey, char* FieldVectID)
{
	radThg hg;
	if(!ValidateElemKey(InteractElemKey, hg)) return;
	radTInteraction* InteractPtr = Cast.InteractCast(hg.rep); 
	if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return;}

	if(!strcmp(FieldVectID, "ext")) InteractPtr->ShowInteractVector('E');
	else if(!strcmp(FieldVectID, "tot")) InteractPtr->ShowInteractVector('T');
	else if(!strcmp(FieldVectID, "mag")) InteractPtr->ShowInteractVector('M');
	else { Send.ErrorMessage("Radia::Error020"); return;}
}

//-------------------------------------------------------------------------

int radTApplication::MakeManualRelax(int InteractElemKey, int MethNo, int IterNumber, double RelaxParam)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(InteractElemKey, hg)) return 0;
		radTInteraction* InteractPtr = Cast.InteractCast(hg.rep); 
		if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return 0;}

		if(MethNo<0 || MethNo>4) { Send.ErrorMessage("Radia::Error028"); return 0;}
		if(IterNumber<0) { Send.ErrorMessage("Radia::Error019"); return 0;}
		if((RelaxParam<0.) || (RelaxParam>1.)) { Send.ErrorMessage("Radia::Error018"); return 0;}

		switch(MethNo)
		{
		case 0:
			{ 
				InteractPtr->ResetM();
			}
			break;
		case 1:
			{
				radTSimpleRelaxation SimpleRelaxation(InteractPtr, RelaxParam);
				SimpleRelaxation.MakeN_iter(IterNumber);
			}
			break;
		case 2:
			{
				radTRelaxationMethNo_2 RelaxationMethNo_2(InteractPtr, RelaxParam);
				RelaxationMethNo_2.MakeN_iter(IterNumber);
			}
			break;
		case 3:
			{
				radTRelaxationMethNo_3 RelaxationMethNo_3(InteractPtr);
				RelaxationMethNo_3.MakeN_iter(IterNumber);
			}
			break;
		case 5:
			{
				if(InteractPtr->AmOfRelaxSubInterv == 0)
				{
					radTRelaxationMethNo_3 RelaxationMethNo_3(InteractPtr);
					RelaxationMethNo_3.MakeN_iter(IterNumber);
				}
				else
				{
					radTRelaxationMethNo_a5 RelaxationMethNo_a5(InteractPtr);
					RelaxationMethNo_a5.MakeN_iter(IterNumber);
				}
			}
			break;
		}

		int lenRelaxStatusParamArray = 3;
		double RelaxStatusParamArray[3];
		InteractPtr->OutRelaxStatusParam(RelaxStatusParamArray);

		if(SendingIsRequired) Send.DoubleList(RelaxStatusParamArray, lenRelaxStatusParamArray);
		return InteractElemKey;
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::MakeAutoRelax(int InteractElemKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, const char** arOptionNames, const char** arOptionValues, int numOptions)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(InteractElemKey, hg)) return 0;
		radTInteraction* InteractPtr = Cast.InteractCast(hg.rep); 
		if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return 0;}

		if(PrecOnMagnetiz <= 0.) { Send.ErrorMessage("Radia::Error030"); return 0;}
		if(MaxIterNumber <= 0) { Send.ErrorMessage("Radia::Error031"); return 0;}

		//if((MethNo < 3) || (MethNo > 5)) { Send.ErrorMessage("Radia::Error041"); return 0;}
		if(MethNo < 3) { Send.ErrorMessage("Radia::Error041"); return 0;}

		radTOptionNames OptNam;
		const char** BufNameString = arOptionNames;
		const char** BufValString = arOptionValues;
		char MagnResetIsNotNeeded = 0;
		for(int i=0; i<numOptions; i++)
		{
			if(!strcmp(*BufNameString, OptNam.ZeroM))
			{
				if(!strcmp(*BufValString, (OptNam.ZeroM_Values)[0])) MagnResetIsNotNeeded = 1; //no
				else if(!strcmp(*BufValString, (OptNam.ZeroM_Values)[1])) MagnResetIsNotNeeded = 0; //yes
				else if(!strcmp(*BufValString, (OptNam.ZeroM_Values)[2])) MagnResetIsNotNeeded = 1; //false
				else if(!strcmp(*BufValString, (OptNam.ZeroM_Values)[3])) MagnResetIsNotNeeded = 0; //true
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
			BufNameString++; BufValString++;
		}

		int ActualIterNum = 0;
		switch(MethNo)
		{
		case 3:
			{
				radTRelaxationMethNo_3 RelaxMethNo_3(InteractPtr);
				ActualIterNum = RelaxMethNo_3.AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
			}
			break;
		case 4:
			{
				radTRelaxationMethNo_4 RelaxMethNo_4(InteractPtr);
				ActualIterNum = RelaxMethNo_4.AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
			}
			break;
		case 5:
			{
				if(InteractPtr->AmOfRelaxSubInterv == 0)
				{
					radTRelaxationMethNo_3 RelaxMethNo_3(InteractPtr);
					ActualIterNum = RelaxMethNo_3.AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
				}
				else
				{
					radTRelaxationMethNo_a5 RelaxMethNo_a5(InteractPtr);
					ActualIterNum = RelaxMethNo_a5.AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
				}
			}
			break;
		case 8:
			{
				radTRelaxationMethNo_8 RelaxMethNo_8(InteractPtr);
				ActualIterNum = RelaxMethNo_8.AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
			}
			break;
		}

		int lenRelaxStatusParamArray = 3;
		double RelaxStatusParamArray[3];
		InteractPtr->OutRelaxStatusParam(RelaxStatusParamArray);

		if(ActualIterNum >= MaxIterNumber) { Send.WarningMessage("Radia::Warning015");}
		if(SendingIsRequired) Send.OutRelaxResultsInfo(RelaxStatusParamArray, lenRelaxStatusParamArray, ActualIterNum);

		return ActualIterNum;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::UpdateSourcesForRelax(int InteractElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(InteractElemKey, hg)) return 0;
		radTInteraction* InteractPtr = Cast.InteractCast(hg.rep); 
		if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return 0;}

		InteractPtr->UpdateExternalField();

		if(SendingIsRequired) Send.Int(InteractElemKey);
		return InteractElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SolveGen(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo)
{
	long ActualIterNum = 0;
	try
	{
		if((MethNo >= 6) && (MethNo <= 7))
		{
			radThg hg;
			if(!ValidateElemKey(ObjKey, hg)) return 0;
			radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
			if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
            radTGroup* GroupPtr = Cast.GroupCast(g3dPtr); 
			if(GroupPtr==0) { Send.ErrorMessage("Radia::Error091"); return 0;}

			double RelaxStatusParamArray[3];

			if(MethNo == 6)
			{
				radTRelaxationMethNo_6 SolveMethNo_6(hg, CompCriterium);
				ActualIterNum = SolveMethNo_6.AutoRelax(PrecOnMagnetiz, MaxIterNumber, RelaxStatusParamArray);
			}
			else if(MethNo == 7)
			{
				radTRelaxationMethNo_7 SolveMethNo_7(hg, CompCriterium);
				ActualIterNum = SolveMethNo_7.AutoRelax(PrecOnMagnetiz, MaxIterNumber, RelaxStatusParamArray);
			}

			if(ActualIterNum >= MaxIterNumber) 
			{ 
				Send.WarningMessage("Radia::Warning015");
			}
			if(SendingIsRequired) 
			{
				Send.OutRelaxResultsInfo(RelaxStatusParamArray, 3, ActualIterNum);
			}
		}
		else
		{
			short PrevSendingIsRequired = SendingIsRequired;
			SendingIsRequired = 0;

			int InteractElemKey = PreRelax(ObjKey, 0);
			if(InteractElemKey <= 0) return 0;

			SendingIsRequired = PrevSendingIsRequired;

			try
			{
				ActualIterNum = MakeAutoRelax(InteractElemKey, PrecOnMagnetiz, MaxIterNumber, MethNo);
			}
			catch(...)
			{
				SendingIsRequired = 0;
				DeleteElement(InteractElemKey);
				throw 0;
			}

			PrevSendingIsRequired = SendingIsRequired; SendingIsRequired = 0;
			DeleteElement(InteractElemKey);
			SendingIsRequired = PrevSendingIsRequired;
		}
	}
	catch(...) { Initialize(); return 0;}
	return ActualIterNum;
}

//-------------------------------------------------------------------------
/**
void radTApplication::OutFieldCompRes(char* FieldChar, radTField* FieldArray, double* ArgArray, int Np)
{
	char* BufChar = FieldChar;
	char* EqEmptyStr = "BHAM";

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b') || 
			   (*BufChar == 'H') || (*BufChar == 'h') ||
			   (*BufChar == 'A') || (*BufChar == 'a') ||
			   (*BufChar == 'M') || (*BufChar == 'm') ||
			   (*BufChar == 'P') || (*BufChar == 'p')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 4;
	}
	char* ActualInitCharPtr = BufChar;

	if(Np > 1) Send.InitOutList(Np);

	radTField* FieldPtr = FieldArray;
	for(int i=0; i<Np; i++)
	{
		if(ArgArray != NULL) // Argument Needed
		{
			Send.InitOutList(2);
			Send.Double(ArgArray[i]);
		}

		if(ItemCount > 1) Send.InitOutList(ItemCount);
		while (*BufChar != '\0') 
		{
			char* BufChar_p_1 = BufChar+1;
			if(*(BufChar)=='B' || *(BufChar)=='b')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->B.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->B.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->B.z);
				else Send.Vector3d(&(FieldPtr->B));
			}
			else if(*(BufChar)=='H' || *(BufChar)=='h')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->H.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->H.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->H.z);
				else Send.Vector3d(&(FieldPtr->H));
			}
			else if(*(BufChar)=='A' || *(BufChar)=='a')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->A.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->A.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->A.z);
				else Send.Vector3d(&(FieldPtr->A));
			}
			else if(*(BufChar)=='M' || *(BufChar)=='m')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->M.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->M.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->M.z);
				else Send.Vector3d(&(FieldPtr->M));
			}
			else if(*(BufChar)=='P' || *(BufChar)=='p')	Send.Double(FieldPtr->Phi);
			BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}
}
**/
//-------------------------------------------------------------------------

void radTApplication::OutFieldCompRes(char* FieldChar, radTField* FieldArray, long Np)
{
	char* BufChar = FieldChar;
	//char* EqEmptyStr = "BHAMJ";
	char EqEmptyStr[] = "BHAMJ"; //OC01052013

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b') || 
			   (*BufChar == 'H') || (*BufChar == 'h') ||
			   (*BufChar == 'A') || (*BufChar == 'a') ||
			   (*BufChar == 'M') || (*BufChar == 'm') ||
			   (*BufChar == 'J') || (*BufChar == 'j') ||
			   (*BufChar == 'P') || (*BufChar == 'p') ||
			   (*BufChar == 'Q') || (*BufChar == 'q')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 4;
	}
	char* ActualInitCharPtr = BufChar;

	double *TotOutArray = new double[Np*12];
	if(TotOutArray == 0) { Send.ErrorMessage("Radia::Error900"); return;}

	int InnerCount=0;
	radTField* FieldPtr = FieldArray;
	double *t = TotOutArray;
	for(int i=0; i<Np; i++)
	{
		InnerCount = 0;
		while (*BufChar != '\0') 
		{
			char* BufChar_p_1 = BufChar+1;
			if(*(BufChar)=='B' || *(BufChar)=='b')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->B.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->B.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->B.z; InnerCount++;}
				else { *(t++) = FieldPtr->B.x; *(t++) = FieldPtr->B.y; *(t++) = FieldPtr->B.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='H' || *(BufChar)=='h')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->H.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->H.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->H.z; InnerCount++;}
				else { *(t++) = FieldPtr->H.x; *(t++) = FieldPtr->H.y; *(t++) = FieldPtr->H.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='A' || *(BufChar)=='a')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->A.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->A.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->A.z; InnerCount++;}
				else { *(t++) = FieldPtr->A.x; *(t++) = FieldPtr->A.y; *(t++) = FieldPtr->A.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='M' || *(BufChar)=='m')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->M.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->M.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->M.z; InnerCount++;}
				else { *(t++) = FieldPtr->M.x; *(t++) = FieldPtr->M.y; *(t++) = FieldPtr->M.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='J' || *(BufChar)=='j')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->J.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->J.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->J.z; InnerCount++;}
				else { *(t++) = FieldPtr->J.x; *(t++) = FieldPtr->J.y; *(t++) = FieldPtr->J.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='Q' || *(BufChar)=='q') //OC161005
			{
				*(t++) = FieldPtr->B.x; *(t++) = FieldPtr->B.y; *(t++) = FieldPtr->B.z; InnerCount += 3;
				*(t++) = FieldPtr->H.x; *(t++) = FieldPtr->H.y; *(t++) = FieldPtr->H.z; InnerCount += 3;
				*(t++) = FieldPtr->A.x; *(t++) = FieldPtr->A.y; *(t++) = FieldPtr->A.z; InnerCount += 3;
			}
			
			BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}
	int Dims[] = { InnerCount, Np};
	Send.MultiDimArrayOfDouble(TotOutArray, Dims, 2);

	if(TotOutArray != 0) delete[] TotOutArray;
}

//-------------------------------------------------------------------------

void radTApplication::OutFieldCompRes(char* FieldChar, radTField* FieldArray, long Np, radTVectInputCell& VectInputCell)
{
	char* BufChar = FieldChar;
	//char* EqEmptyStr = "BHAMJ";
	char EqEmptyStr[] = "BHAMJ"; //OC01052013

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b') || 
			   (*BufChar == 'H') || (*BufChar == 'h') ||
			   (*BufChar == 'A') || (*BufChar == 'a') ||
			   (*BufChar == 'M') || (*BufChar == 'm') ||
			   (*BufChar == 'J') || (*BufChar == 'j') ||
			   (*BufChar == 'P') || (*BufChar == 'p')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 4;
	}
	char* ActualInitCharPtr = BufChar;

	radTField* tField = FieldArray;
	for(radTVectInputCell::iterator iterCell = VectInputCell.begin(); iterCell != VectInputCell.end(); ++iterCell)
	{
		radTInputCell& Cell = *iterCell;
		if(Cell.Type == 'L')
		{
			Send.InitOutList(Cell.AuxNum);
		}
		else if(Cell.Type == 'P')
		{
			ParseAndSendOneFieldValue(tField++, BufChar, ItemCount);
		}
	}
}

//-------------------------------------------------------------------------

void radTApplication::ParseAndSendOneFieldValue(radTField* FieldPtr, char* BufChar, int AmOfItem)
{
	if(AmOfItem > 1) Send.InitOutList(AmOfItem);
	while(*BufChar != '\0') 
	{
		char* BufChar_p_1 = BufChar+1;
		if(*(BufChar)=='B' || *(BufChar)=='b')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->B.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->B.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->B.z);
			else Send.Vector3d(&(FieldPtr->B));
		}
		else if(*(BufChar)=='H' || *(BufChar)=='h')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->H.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->H.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->H.z);
			else Send.Vector3d(&(FieldPtr->H));
		}
		else if(*(BufChar)=='A' || *(BufChar)=='a')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->A.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->A.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->A.z);
			else Send.Vector3d(&(FieldPtr->A));
		}
		else if(*(BufChar)=='M' || *(BufChar)=='m')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->M.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->M.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->M.z);
			else Send.Vector3d(&(FieldPtr->M));
		}
		else if(*(BufChar)=='J' || *(BufChar)=='j')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->J.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->J.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->J.z);
			else Send.Vector3d(&(FieldPtr->J));
		}
		else if(*(BufChar)=='P' || *(BufChar)=='p')	Send.Double(FieldPtr->Phi);
		BufChar++;
	}
}

//-------------------------------------------------------------------------
/**
void radTApplication::OutFieldIntCompRes(char* FieldIntChar, radTField* FieldPtr)
{
	char* BufChar = FieldIntChar;
	char* BufCharPrev = NULL;
	char* EqEmptyStr = "Ib";

	short I_used = 0;
	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if(((*BufChar == 'B') || (*BufChar == 'b') || 
			    (*BufChar == 'H') || (*BufChar == 'h')) ||
			   (((*BufChar == 'X') || (*BufChar == 'x') ||
			     (*BufChar == 'Y') || (*BufChar == 'y') ||
				 (*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufCharPrev != 'B') && (*BufCharPrev != 'b') &&
				(*BufCharPrev != 'H') && (*BufCharPrev != 'h'))) ItemCount++;

			if((*BufChar == 'I') || (*BufChar == 'i')) I_used = 1;
			BufCharPrev = BufChar;
			BufChar++;
		}
		BufChar = FieldIntChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 1;
	}
	if(I_used && (ItemCount == 0))
	{
		BufChar = EqEmptyStr;
		ItemCount = 1;
	}

	if(ItemCount > 1) Send.InitOutList(ItemCount);

	while (*BufChar != '\0') 
	{
		char* BufChar_pl_1 = BufChar+1;
		char* BufChar_mi_1 = BufChar-1;

		if((*BufChar =='I') || (*BufChar == 'i'))
		{
			if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Send.Double(FieldPtr->Ib.x);
			else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Send.Double(FieldPtr->Ib.y);
			else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Send.Double(FieldPtr->Ib.z);
			else if((*BufChar_pl_1 != 'B') && (*BufChar_pl_1 != 'b') &&
					(*BufChar_pl_1 != 'H') && (*BufChar_pl_1 != 'h') &&
					(*BufChar_pl_1 != 'X') && (*BufChar_pl_1 != 'x') &&
					(*BufChar_pl_1 != 'Y') && (*BufChar_pl_1 != 'y') &&
					(*BufChar_pl_1 != 'Z') && (*BufChar_pl_1 != 'z')) { Send.Vector3d(&(FieldPtr->Ib));	break;}
		}
		else if((*BufChar == 'B') || (*BufChar == 'b'))
		{
			if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Send.Double(FieldPtr->Ib.x);
			else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Send.Double(FieldPtr->Ib.y);
			else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Send.Double(FieldPtr->Ib.z);
			else Send.Vector3d(&(FieldPtr->Ib));
		}
		else if((*BufChar == 'H') || (*BufChar == 'h'))
		{
			if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Send.Double(FieldPtr->Ih.x);
			else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Send.Double(FieldPtr->Ih.y);
			else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Send.Double(FieldPtr->Ih.z);
			else Send.Vector3d(&(FieldPtr->Ih));
		}
		else if(((*BufChar == 'X') || (*BufChar == 'x')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
 				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Send.Double(FieldPtr->Ib.x);
		else if(((*BufChar == 'Y') || (*BufChar == 'y')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
 				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Send.Double(FieldPtr->Ib.y);
		else if(((*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
 				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Send.Double(FieldPtr->Ib.z);
		BufChar++;
	}
}
**/
//-------------------------------------------------------------------------

void radTApplication::OutFieldEnergyForceCompRes(char* FieldChar, radTField* FieldPtr)
{
	char* BufChar = FieldChar;
	//char* EqEmptyStr = "EFT";
	char EqEmptyStr[] = "EFT"; //OC01052013

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'E') || (*BufChar == 'e') || 
			   (*BufChar == 'F') || (*BufChar == 'f') ||
			   (*BufChar == 'T') || (*BufChar == 't')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}
	char* ActualInitCharPtr = BufChar;

	if(ItemCount > 1) Send.InitOutList(ItemCount);
	while (*BufChar != '\0') 
	{
		char* BufChar_p_1 = BufChar+1;

		if(*(BufChar)=='E' || *(BufChar)=='e') Send.Double(FieldPtr->Energy);
		else if(*(BufChar)=='F' || *(BufChar)=='f')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->Force.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->Force.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->Force.z);
			else Send.Vector3d(&(FieldPtr->Force));
		}
		else if(*(BufChar)=='T' || *(BufChar)=='t')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->Torque.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->Torque.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->Torque.z);
			else Send.Vector3d(&(FieldPtr->Torque));
		}
		BufChar++;
	}
}

//-------------------------------------------------------------------------

void radTApplication::OutCenFieldCompRes(radTVectPairOfVect3d* pVectPairOfVect3d)
{
	int AmOfPoints = (int)(pVectPairOfVect3d->size());
	radTSend Send;
	Send.InitOutList(AmOfPoints, 0);

	for(int i=0; i<AmOfPoints; i++)
	{
		Send.InitOutList(2, 0);
		radTPairOfVect3d& Pair = (*pVectPairOfVect3d)[i];
		Send.Vector3d(&(Pair.V1));
		Send.Vector3d(&(Pair.V2));
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetCompPrecisions(const char** ValNames, double* Values, int ValCount)
{
	try
	{
		radTOptionNames OptionNames;
		const char** BufString = ValNames;
		double* Ptr = Values;
		for(int i=0; i<ValCount; i++)
		{
			if(!strcmp(*BufString, OptionNames.B)) CompCriterium.AbsPrecB = *Ptr;
			else if(!strcmp(*BufString, OptionNames.A)) CompCriterium.AbsPrecA = *Ptr;
			else if(!strcmp(*BufString, OptionNames.BInt)) CompCriterium.AbsPrecB_int = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Force)) CompCriterium.AbsPrecForce = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Torque)) CompCriterium.AbsPrecTorque = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Energy)) CompCriterium.AbsPrecEnergy = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Coord)) CompCriterium.AbsPrecTrjCoord = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Angle)) CompCriterium.AbsPrecTrjAngle = *Ptr;
			else { Send.ErrorMessage("Radia::Error057"); return 0;}
			BufString++; Ptr++;
		}
		if(SendingIsRequired) Send.Int(1);
		return 1;
	}
	catch(...) { Initialize(); return 0;}
}

//-------------------------------------------------------------------------

int radTApplication::SetCompCriterium(double InAbsPrecB, double InAbsPrecA, double InAbsPrecB_int, double InAbsPrecFrc, double InAbsPrecTrjCoord, double InAbsPrecTrjAngle)
{
	short InBasedOnPrecFlag = 0;

	try
	{
		CompCriterium.BasedOnPrecLevel = InBasedOnPrecFlag; 
		CompCriterium.AbsPrecB = InAbsPrecB;
		CompCriterium.AbsPrecA = InAbsPrecA;
		CompCriterium.AbsPrecB_int = InAbsPrecB_int;
		CompCriterium.AbsPrecForce = InAbsPrecFrc;
		CompCriterium.AbsPrecTrjCoord = InAbsPrecTrjCoord;
		CompCriterium.AbsPrecTrjAngle = InAbsPrecTrjAngle;

		if(SendingIsRequired) Send.Int(1);
		return 1;
	}
	catch(...) { Initialize(); return 0;}
}

//-------------------------------------------------------------------------

int radTApplication::SetMltplThresh(double* InMltplThresh) // Maybe to be removed later
{
	for(int i=0; i<4; i++) CompCriterium.MltplThresh[i] = InMltplThresh[i]*InMltplThresh[i];
	if(SendingIsRequired) Send.Int(1);
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::SetTolForConvergence(double AbsRandMagnitude, double RelRandMagnitude, double ZeroRandMagnitude)
{
	CnRep.SwitchActOnDoubles(1, fabs(AbsRandMagnitude), fabs(RelRandMagnitude), fabs(ZeroRandMagnitude));
	if(SendingIsRequired) Send.Int(1);
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::RandomizationOnOrOff(char* OnOrOff)
{
	char SwitchOn;
	if((!strcmp(OnOrOff, "on")) || (!strcmp(OnOrOff, "On")) || (!strcmp(OnOrOff, "ON"))) SwitchOn = 1;
	else if((!strcmp(OnOrOff, "off")) || (!strcmp(OnOrOff, "Off")) || (!strcmp(OnOrOff, "OFF"))) SwitchOn = 0;
	else { Send.ErrorMessage("Radia::Error043"); return 0;}

	if(SwitchOn) CnRep = CnRepAux;
	else
	{
		CnRepAux = CnRep;
		CnRep.AbsRand = CnRep.RelRand = CnRep.ZeroRand = 0.;
	}
	if(SendingIsRequired) Send.Int(int(SwitchOn));
	return 1;
}

//-------------------------------------------------------------------------
