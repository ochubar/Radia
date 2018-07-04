
#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

#ifndef __IGINTERF_H
#include "iginterf.h"
#endif
#ifndef __SRINTERF_H
#include "radentry.h"
#endif
#ifndef __IGERCODE_H
#include "igercode.h"
#endif
#ifndef __AUXPARSE_H
#include "auxparse.h"
#endif

//*************************************************************************

// All structures are 2-byte-aligned.
#if GENERATINGPOWERPC
	#pragma options align=mac68k
#endif
#ifdef _WINDOWS_
	#pragma pack(2)
#endif

/* Global Variables (none) */
//void
//Go(IORecHandle ioRecHandle)
//{
//	HOST_IMPORT void main(IORecHandle);
//	
//	main(ioRecHandle);
//}
HOST_IMPORT void main(IORecHandle); //OC080613 (moving to XOPSupport 6)
void Go(IORecHandle ioRecHandle)
{
	//HOST_IMPORT void main(IORecHandle);
	main(ioRecHandle);
}

//*************************************************************************

static int CorrectErrorCode(int ErrCode)
{
	//ErrCode += (FIRST_XOP_ERR - FIRST_ERR_SRWDLL);
	if(ErrCode == 1) return SR_COMP_PROC_ABORTED; // 1 is returned when user presses Abort
	return ErrCode;
}

//*************************************************************************

static int CorrectErrorCodeRad(int ErrCode)
{
	const int FirstRadiaErrNo = 100;
	return FIRST_XOP_ERR + FirstRadiaErrNo + ErrCode;
}

//*************************************************************************

static void ProcErr(int ErrNo)
{
	if(ErrNo == 0) return;
	else if(ErrNo < 0) //warning
	{
		char* WarnStr = 0;
		int WarnLen = 0;
        RadWarGetSize(&WarnLen, ErrNo);
		if(WarnLen > 0) 
		{
			WarnStr = new char[WarnLen + 10];
            RadWarGetText(WarnStr, ErrNo);
	        srTIgorInterf::WarningMessage(WarnStr);
		}
		if(WarnStr != 0) delete[] WarnStr;
	}
	else throw CorrectErrorCode(ErrNo);
}

//*************************************************************************

static void ProcErrRad(int ErrNo)
{
	if(ErrNo == 0) return;
	if(ErrNo < 0) { ProcErr(ErrNo); return;}
	else throw CorrectErrorCodeRad(ErrNo);
}

//*************************************************************************

struct radTIgorRadObjRecMag {
    waveHndl wM;
    waveHndl wL;
    waveHndl wP;
	DOUBLE obj;
};
static int radObjRecMag(void* p_void)
{
	try
	{
		radTIgorRadObjRecMag *p = (radTIgorRadObjRecMag*)p_void;
		p->obj = 0;
		double P[3], L[3], M[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wL, L));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));
		int Obj=0;
		ProcErrRad(RadObjRecMag(&Obj, P, L, M));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjRecMagC {
    //waveHndl wM;
	DOUBLE mz;
	DOUBLE my;
	DOUBLE mx;
    //waveHndl wL;
	DOUBLE wz;
	DOUBLE wy;
	DOUBLE wx;
    //waveHndl wP;
	DOUBLE z;
	DOUBLE y;
	DOUBLE x;
	DOUBLE obj;
};
static int radObjRecMagC(void* p_void)
{
	try
	{
		radTIgorRadObjRecMagC *p = (radTIgorRadObjRecMagC*)p_void;
		p->obj = 0;
		//double P[3], L[3], M[3];
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wL, L));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));
		double P[] = {p->x, p->y, p->z};
		double L[] = {p->wx, p->wy, p->wz};
		double M[] = {p->mx, p->my, p->mz};

		int Obj=0;
		ProcErrRad(RadObjRecMag(&Obj, P, L, M));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjThckPgn {
    waveHndl wM;
    Handle a;
    waveHndl wVertCoord;
	DOUBLE lx;
	DOUBLE xc;
	DOUBLE obj;
};
static int radObjThckPgn(void* p_void)
{
	double* pCoord=0;
	try
	{
		radTIgorRadObjThckPgn *p = (radTIgorRadObjThckPgn*)p_void;
		p->obj = 0;

		int Nv = -1, DimP = 2;
		//ProcErr(srTIgorInterf::GetArrDoubleAndNpFromNumWave1D(p->wVertCoord, (int)(p->nv), 2, pCoord, Nv));
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wVertCoord, DimP, Nv, pCoord));

		double M[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));

		char aStr[3];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));

		int Obj=0;
		ProcErrRad(RadObjThckPgn(&Obj, (double)(p->xc), (double)(p->lx), pCoord, Nv, aStr[0], M));

		p->obj = Obj;
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
	}
	catch(int ErrNo) 
	{
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjPolyhdr {
    waveHndl wM;
	waveHndl wFacesLengths;
	waveHndl wFaces;
    waveHndl wVertCoord;
	DOUBLE obj;
};
static int radObjPolyhdr(void* p_void)
{
	double* pCoord=0;
	int *pFaces=0, *pFacesLengths=0;

	try
	{
		radTIgorRadObjPolyhdr *p = (radTIgorRadObjPolyhdr*)p_void;
		p->obj = 0;

		int Nv = -1, DimP = 3;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wVertCoord, DimP, Nv, pCoord));

		long N_Faces = 0, N_FacesLengths;
        ProcErr(srTIgorInterf::GetArrIntFromNumWave1D(p->wFaces, -1, pFaces, N_Faces));
        ProcErr(srTIgorInterf::GetArrIntFromNumWave1D(p->wFacesLengths, -1, pFacesLengths, N_FacesLengths));

		double M[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));

		int Obj=0;
		//ProcErrRad(RadObjPolyhdr(&Obj, pCoord, Nv, pFaces, pFacesLengths, N_FacesLengths, M));
		ProcErrRad(RadObjPolyhdr(&Obj, pCoord, Nv, pFaces, pFacesLengths, N_FacesLengths, M, 0, 0, 0)); //to update

		p->obj = Obj;
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pFaces != 0) { delete[] pFaces; pFaces = 0;}
		if(pFacesLengths != 0) { delete[] pFacesLengths; pFacesLengths = 0;}
	}
	catch(int ErrNo) 
	{
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pFaces != 0) { delete[] pFaces; pFaces = 0;}
		if(pFacesLengths != 0) { delete[] pFacesLengths; pFacesLengths = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjMltExtPgn {
    waveHndl wM;
	//DOUBLE ns;
    waveHndl wAttitudes;
	waveHndl wSlicesLen;
    waveHndl wVertCoord;
	DOUBLE obj;
};
static int radObjMltExtPgn(void* p_void)
{
	double *pCoord=0, *pAttitudes=0;
	int* pSlicesLen=0;
	try
	{
		radTIgorRadObjMltExtPgn *p = (radTIgorRadObjMltExtPgn*)p_void;
		p->obj = 0;

		//long nVertCoordAlloc=0, nSlicesLenAlloc=0, nAttitudesAlloc=0;
		//ProcErr(srTIgorInterf::GetArrDoubleFromNumWave1D(p->wVertCoord, -1, pCoord, nVertCoordAlloc));

		int DimP = 2, nVertices = -1;
		long nSlicesLenAlloc=0, nAttitudesAlloc=0;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wVertCoord, DimP, nVertices, pCoord));
		ProcErr(srTIgorInterf::GetArrIntFromNumWave1D(p->wSlicesLen, -1, pSlicesLen, nSlicesLenAlloc));
		ProcErr(srTIgorInterf::GetArrDoubleFromNumWave1D(p->wAttitudes, -1, pAttitudes, nAttitudesAlloc));

		//if(nAttitudesAlloc < (long)(p->ns)) throw TOO_SHORT_WAVE;
		if(nSlicesLenAlloc != nAttitudesAlloc) throw INCOMPATIBLE_WAVE_LENGTHS;

		double M[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));

		int Obj=0;
		//ProcErrRad(RadObjMltExtPgn(&Obj, pCoord, pSlicesLen, pAttitudes, (int)(p->ns), M));
		ProcErrRad(RadObjMltExtPgn(&Obj, pCoord, pSlicesLen, pAttitudes, nAttitudesAlloc, M));

		p->obj = Obj;
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pSlicesLen != 0) { delete[] pSlicesLen; pSlicesLen = 0;}
		if(pAttitudes != 0) { delete[] pAttitudes; pAttitudes = 0;}
	}
	catch(int ErrNo) 
	{
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pSlicesLen != 0) { delete[] pSlicesLen; pSlicesLen = 0;}
		if(pAttitudes != 0) { delete[] pAttitudes; pAttitudes = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjMltExtRtg {
    waveHndl wM;
	waveHndl wRtgSizes;
    waveHndl wCenPtsCoord;
	DOUBLE obj;
};
static int radObjMltExtRtg(void* p_void)
{
	double *pCenPtsCoord=0, *pRtgSizes=0;
	try
	{
		radTIgorRadObjMltExtRtg *p = (radTIgorRadObjMltExtRtg*)p_void;
		p->obj = 0;

		//long nCenPtsCoordAlloc=0, nRtgSizesAlloc=0;
		//ProcErr(srTIgorInterf::GetArrDoubleFromNumWave1D(p->wCenPtsCoord, -1, pCenPtsCoord, nCenPtsCoordAlloc));
		//ProcErr(srTIgorInterf::GetArrDoubleFromNumWave1D(p->wRtgSizes, -1, pRtgSizes, nRtgSizesAlloc));
		//if(nCenPtsCoordAlloc < (long)(3*(p->ns))) throw TOO_SHORT_WAVE;
		//if(nRtgSizesAlloc < (long)(2*(p->ns))) throw TOO_SHORT_WAVE;

		int DimP = 3, DimR = 2, nPts = -1, nRtg = -1;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wCenPtsCoord, DimP, nPts, pCenPtsCoord));
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wRtgSizes, DimR, nRtg, pRtgSizes));
		if(nPts != nRtg) throw INCOMPATIBLE_WAVE_LENGTHS;

		double M[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));

		int Obj=0;
		//ProcErrRad(RadObjMltExtRtg(&Obj, pCenPtsCoord, pRtgSizes, (int)(p->ns), M));
		ProcErrRad(RadObjMltExtRtg(&Obj, pCenPtsCoord, pRtgSizes, nPts, M));

		p->obj = Obj;
		if(pCenPtsCoord != 0) { delete[] pCenPtsCoord; pCenPtsCoord = 0;}
		if(pRtgSizes != 0) { delete[] pRtgSizes; pRtgSizes = 0;}
	}
	catch(int ErrNo) 
	{
		if(pCenPtsCoord != 0) { delete[] pCenPtsCoord; pCenPtsCoord = 0;}
		if(pRtgSizes != 0) { delete[] pRtgSizes; pRtgSizes = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjMltExtTri {
    Handle hOpt;
	waveHndl wM;
    Handle a;
    waveHndl wSbdPar;
    waveHndl wVertCoord;
	DOUBLE lx;
	DOUBLE xc;
	DOUBLE obj;
};
static int radObjMltExtTri(void* p_void)
{
	double *pCoord=0, *pSbd=0;
	try
	{
		radTIgorRadObjMltExtTri *p = (radTIgorRadObjMltExtTri*)p_void;
		p->obj = 0;

		int Nv = -1, DimP = 2;
		ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wVertCoord, DimP, Nv, pCoord));

		int InnerDimSizeSbd = -1;
		ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wSbdPar, InnerDimSizeSbd, Nv, pSbd));

		if(InnerDimSizeSbd != 2)
		{
			double *pSbdAux = new double[Nv << 1];
			double *tSbdAux = pSbdAux, *tSbd = pSbd;
			for(int i=0; i<Nv; i++)
			{
				*(tSbdAux++) = *(tSbd++);
				*(tSbdAux++) = 1.;
			}
			delete[] pSbd;
			pSbd = pSbdAux;
			pSbdAux = 0;
		}

		char aStr[3];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));

		double M[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));

        char sOpt[512];
		*sOpt = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt, 511, sOpt));

		int Obj=0;
		ProcErrRad(RadObjMltExtTri(&Obj, (double)(p->xc), (double)(p->lx), pCoord, pSbd, Nv, aStr[0], M, sOpt));
		p->obj = Obj;

		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pSbd != 0) { delete[] pSbd; pSbd = 0;}
	}
	catch(int ErrNo) 
	{
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pSbd != 0) { delete[] pSbd; pSbd = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

//struct radTIgorRadObjArcMag {
//    waveHndl wM;
//    Handle a;
//	DOUBLE nseg;
//	DOUBLE h;
//    waveHndl wPhi;
//    waveHndl wR;
//    waveHndl wP;
//	DOUBLE obj;
//};
//static int radObjArcMag(void* p_void)
//{
//	try
//	{
//		radTIgorRadObjArcMag *p = (radTIgorRadObjArcMag*)p_void;
//		p->obj = 0;
//		double P[3], M[3], R[2], Phi[2];
//		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
//		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));
//		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wR, R));
//		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wPhi, Phi));
//
//		char aStr[3];
//		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));
//
//		int Obj=0;
//		ProcErrRad(RadObjArcMag(&Obj, P, R, Phi, p->h, (int)(p->nseg), aStr[0], M));
//		p->obj = Obj;
//	}
//	catch(int ErrNo) 
//	{
//		return ErrNo;
//	}
//	return 0;
//}

//*************************************************************************

struct radTIgorRadObjArcPgnMag {
    waveHndl wM;
    Handle hSym;
	DOUBLE nseg;
    waveHndl wPhi;
    waveHndl wFlatVert;
    Handle ha;
    waveHndl wP;
	DOUBLE obj;
};
static int radObjArcPgnMag(void* p_void)
{
	double *pPtsCoord=0;
	try
	{
		radTIgorRadObjArcPgnMag *p = (radTIgorRadObjArcPgnMag*)p_void;
		p->obj = 0;
		double P[2], Phi[2], M[3];

		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wPhi, Phi));
        ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));

		char aStr[3], SymStr[10];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->ha, 2, aStr));
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hSym, 6, SymStr));

		int DimP = 2, nPts = -1;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wFlatVert, DimP, nPts, pPtsCoord));

		int Obj=0;
        ProcErrRad(RadObjArcPgnMag(&Obj, P, aStr[0], pPtsCoord, nPts, Phi, (int)(p->nseg), SymStr[0], M));

		if(pPtsCoord != 0) { delete[] pPtsCoord; pPtsCoord = 0;}
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		if(pPtsCoord != 0) { delete[] pPtsCoord; pPtsCoord = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjCylMag {
    waveHndl wM;
    Handle a;
	DOUBLE nseg;
	DOUBLE h;
	DOUBLE r;
    waveHndl wP;
	DOUBLE obj;
};
static int radObjCylMag(void* p_void)
{
	try
	{
		radTIgorRadObjCylMag *p = (radTIgorRadObjCylMag*)p_void;
		p->obj = 0;
		double P[3], M[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));

		char aStr[3];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));

		int Obj=0;
		ProcErrRad(RadObjCylMag(&Obj, P, p->r, p->h, (int)(p->nseg), aStr[0], M));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjCylMagC {
	//waveHndl wM;
	DOUBLE mz;
	DOUBLE my;
	DOUBLE mx;
	Handle a;
	DOUBLE nseg;
	DOUBLE h;
	DOUBLE r;
	//waveHndl wP;
	DOUBLE z;
	DOUBLE y;
	DOUBLE x;
	DOUBLE obj;
};
static int radObjCylMagC(void* p_void)
{
	try
	{
		radTIgorRadObjCylMagC *p = (radTIgorRadObjCylMagC*)p_void;
		p->obj = 0;
		double P[] = {p->x, p->y, p->z}, M[] = {p->mx, p->my, p->mz};
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));

		char aStr[3];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));

		int Obj=0;
		ProcErrRad(RadObjCylMag(&Obj, P, p->r, p->h, (int)(p->nseg), aStr[0], M));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjRecCur {
    waveHndl wJ;
    waveHndl wL;
    waveHndl wP;
	DOUBLE obj;
};
static int radObjRecCur(void* p_void)
{
	try
	{
		radTIgorRadObjRecCur *p = (radTIgorRadObjRecCur*)p_void;
		p->obj = 0;
		double P[3], L[3], J[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wL, L));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wJ, J));
		int Obj=0;
		ProcErrRad(RadObjRecCur(&Obj, P, L, J));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjRecCurC {
	//waveHndl wJ;
 	DOUBLE jz;
	DOUBLE jy;
	DOUBLE jx;
   //waveHndl wL;
	DOUBLE wz;
	DOUBLE wy;
	DOUBLE wx;
	//waveHndl wP;
	DOUBLE z;
	DOUBLE y;
	DOUBLE x;
	DOUBLE obj;
};
static int radObjRecCurC(void* p_void)
{
	try
	{
		radTIgorRadObjRecCurC *p = (radTIgorRadObjRecCurC*)p_void;
		p->obj = 0;
		double P[] = {p->x, p->y, p->z}, L[] = {p->wx, p->wy, p->wz}, J[] = {p->jx, p->jy, p->jz};
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wL, L));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wJ, J));
		int Obj=0;
		ProcErrRad(RadObjRecCur(&Obj, P, L, J));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjArcCur {
	DOUBLE j;
    Handle a;
    Handle hMode;
	DOUBLE nseg;
	DOUBLE h;
    waveHndl wPhi;
    waveHndl wR;
    waveHndl wP;
	DOUBLE obj;
};
static int radObjArcCur(void* p_void)
{
	try
	{
		radTIgorRadObjArcCur *p = (radTIgorRadObjArcCur*)p_void;
		p->obj = 0;
		double P[3], R[2], Phi[2];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wR, R));
		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wPhi, Phi));

        char StrMode[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hMode, 30, StrMode));

		char aStr[3];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));

		int Obj=0;
		ProcErrRad(RadObjArcCur(&Obj, P, R, Phi, p->h, (int)(p->nseg), StrMode[0], aStr[0], p->j));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjArcCurC {
	DOUBLE j;
    Handle a;
    Handle hMode;
	DOUBLE nseg;
	DOUBLE h;
	//waveHndl wPhi;
	DOUBLE phi2;
	DOUBLE phi1;
	//waveHndl wR;
	DOUBLE r2;
	DOUBLE r1;
	//waveHndl wP;
	DOUBLE z;
	DOUBLE y;
	DOUBLE x;
	DOUBLE obj;
};
static int radObjArcCurC(void* p_void)
{
	try
	{
		radTIgorRadObjArcCurC *p = (radTIgorRadObjArcCurC*)p_void;
		p->obj = 0;
		double P[] = {p->x, p->y, p->z}, R[] = {p->r1, p->r2}, Phi[] = {p->phi1, p->phi2};
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		//ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wR, R));
		//ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wPhi, Phi));

        char StrMode[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hMode, 30, StrMode));

		char aStr[3];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));

		int Obj=0;
		ProcErrRad(RadObjArcCur(&Obj, P, R, Phi, p->h, (int)(p->nseg), StrMode[0], aStr[0], p->j));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjRaceTrk {
	DOUBLE j;
    Handle a;
    Handle hMode;
	DOUBLE nseg;
	DOUBLE h;
    waveHndl wL;
    waveHndl wR;
    waveHndl wP;
	DOUBLE obj;
};
static int radObjRaceTrk(void* p_void)
{
	try
	{
		radTIgorRadObjRaceTrk *p = (radTIgorRadObjRaceTrk*)p_void;
		p->obj = 0;
		double P[3], R[2], L[2];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wR, R));
		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wL, L));

        char StrMode[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hMode, 30, StrMode));

		char aStr[3];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));

		int Obj=0;
		ProcErrRad(RadObjRaceTrk(&Obj, P, R, L, p->h, (int)(p->nseg), StrMode[0], aStr[0], p->j));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjRaceTrkC {
	DOUBLE j;
    Handle a;
    Handle hMode;
	DOUBLE nseg;
	DOUBLE h;
    //waveHndl wL;
	DOUBLE ly;
	DOUBLE lx;
	//waveHndl wR;
	DOUBLE r2;
	DOUBLE r1;
	//waveHndl wP;
 	DOUBLE z;
 	DOUBLE y;
	DOUBLE x;
	DOUBLE obj;
};
static int radObjRaceTrkC(void* p_void)
{
	try
	{
		radTIgorRadObjRaceTrkC *p = (radTIgorRadObjRaceTrkC*)p_void;
		p->obj = 0;
		double P[] = {p->x, p->y, p->z};
		double R[] = {p->r1, p->r2};
		double L[] = {p->lx, p->ly};
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		//ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wR, R));
		//ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wL, L));

        char StrMode[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hMode, 30, StrMode));

		char aStr[3];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->a, 2, aStr));

		int Obj=0;
		ProcErrRad(RadObjRaceTrk(&Obj, P, R, L, p->h, (int)(p->nseg), StrMode[0], aStr[0], p->j));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjFlmCur {
    DOUBLE i;
    //DOUBLE np;
    waveHndl wCoord;
	DOUBLE obj;
};
static int radObjFlmCur(void* p_void)
{
	double* pCoord=0;
	try
	{
        radTIgorRadObjFlmCur *p = (radTIgorRadObjFlmCur*)p_void;
		p->obj = 0;

        //int Np = 0;
		//ProcErr(srTIgorInterf::GetArrDoubleAndNpFromNumWave1D(p->wCoord, (int)(p->np), 3, pCoord, Np));
		int DimP = 3, Np = -1;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wCoord, DimP, Np, pCoord));

		int Obj=0;
        ProcErrRad(RadObjFlmCur(&Obj, pCoord, Np, (int)(p->i)));
		p->obj = Obj;
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
	}
	catch(int ErrNo) 
	{
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjBckg {
    waveHndl wB;
	DOUBLE obj;
};
static int radObjBckg(void* p_void)
{
	try
	{
		radTIgorRadObjBckg *p = (radTIgorRadObjBckg*)p_void;
		p->obj = 0;
		double B[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wB, B));
		int Obj=0;
		ProcErrRad(RadObjBckg(&Obj, B));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjCnt0 {
	DOUBLE obj;
};
static int radObjCnt0(void* p_void)
{
	int* pObjs=0;
	try
	{
        radTIgorRadObjCnt0 *p = (radTIgorRadObjCnt0*)p_void;
		p->obj = 0;

		int Obj=0;
        ProcErrRad(RadObjCnt(&Obj, 0, 0));

		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjCnt {
    //DOUBLE nobj;
    waveHndl wObjs;
	DOUBLE obj;
};
static int radObjCnt(void* p_void)
{
	int* pObjs=0;
	try
	{
        radTIgorRadObjCnt *p = (radTIgorRadObjCnt*)p_void;
		p->obj = 0;

        long nObjAlloc=0;
		ProcErr(srTIgorInterf::GetArrIntFromNumWave1D(p->wObjs, -1, pObjs, nObjAlloc));
		//if(nObjAlloc < (long)(p->nobj)) throw TOO_SHORT_WAVE;

		int Obj=0;
        //ProcErrRad(RadObjCnt(&Obj, pObjs, (int)(p->nobj)));
        ProcErrRad(RadObjCnt(&Obj, pObjs, (int)(nObjAlloc)));

		p->obj = Obj;
		if(pObjs != 0) { delete[] pObjs; pObjs = 0;}
	}
	catch(int ErrNo) 
	{
		if(pObjs != 0) { delete[] pObjs; pObjs = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjAddToCnt {
    //DOUBLE nobj;
    waveHndl wObjs;
	DOUBLE obj_in;
	DOUBLE obj;
};
static int radObjAddToCnt(void* p_void)
{
	int* pObjs=0;
	try
	{
        radTIgorRadObjAddToCnt *p = (radTIgorRadObjAddToCnt*)p_void;
		p->obj = 0;

        long nObjAlloc=0;
		ProcErr(srTIgorInterf::GetArrIntFromNumWave1D(p->wObjs, -1, pObjs, nObjAlloc));
		//if(nObjAlloc < (long)(p->nobj)) throw TOO_SHORT_WAVE;

        //ProcErrRad(RadObjAddToCnt((int)(p->obj_in), pObjs, (int)(p->nobj)));
        ProcErrRad(RadObjAddToCnt((int)(p->obj_in), pObjs, (int)(nObjAlloc)));

		p->obj = (int)(p->obj_in);
		if(pObjs != 0) { delete[] pObjs; pObjs = 0;}
	}
	catch(int ErrNo) 
	{
		if(pObjs != 0) { delete[] pObjs; pObjs = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjCntSize {
	DOUBLE obj_in;
	DOUBLE n;
};
static int radObjCntSize(void* p_void)
{
	try
	{
        radTIgorRadObjCntSize *p = (radTIgorRadObjCntSize*)p_void;
		p->n = 0;

		int CntLen=0;
        ProcErrRad(RadObjCntSize(&CntLen, (int)(p->obj_in)));
		p->n = CntLen;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjCntStuf {
	DOUBLE obj_in;
    waveHndl wObjs;
	DOUBLE obj;
};
static int radObjCntStuf(void* p_void)
{
	float *pf=0;
	double *pd=0;
	int hState, *piObjs=0;
    radTIgorRadObjCntStuf *p = (radTIgorRadObjCntStuf*)p_void;
	try
	{
		p->obj = 0;

		long nObjsAlloc=0;
		ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wObjs, pd, pf, nObjsAlloc, hState));
		if(nObjsAlloc <= 0) throw ZERO_NUMBER_OF_ELEM_IN_WAVE;

		int CntLen=0;
        ProcErrRad(RadObjCntSize(&CntLen, (int)(p->obj_in)));
		if(nObjsAlloc < CntLen) throw TOO_SHORT_WAVE;

		piObjs = new int[CntLen];
		int *tiObjs = piObjs;
		if(pd != 0)
		{
			double *td = pd;
			for(int i=0; i<CntLen; i++) *(tiObjs++) = (int)(*(td++));
			for(int k=CntLen; k<nObjsAlloc; k++) *(td++) = 0;
		}
		else if(pf != 0)
		{
			float *tf = pf;
			for(int i=0; i<CntLen; i++) *(tiObjs++) = (int)(*(tf++));
			for(int k=CntLen; k<nObjsAlloc; k++) *(tf++) = 0;
		}
		
        ProcErrRad(RadObjCntStuf(piObjs, (int)(p->obj_in)));
		p->obj = p->obj_in;

		if(piObjs != 0)
		{
			tiObjs = piObjs;
            if(pd != 0)
			{
                double *td = pd;
                for(int i=0; i<CntLen; i++) *(td++) = *(tiObjs++);
            }
            else if(pf != 0)
            {
                float *tf = pf;
                for(int i=0; i<CntLen; i++) *(tf++) = (float)(*(tiObjs++));
            }

			delete[] piObjs;
			piObjs = 0;
		}
		ProcErr(srTIgorInterf::ReleaseWave(p->wObjs, hState));
	}
	catch(int ErrNo) 
	{
		if(piObjs != 0) { delete[] piObjs; piObjs = 0;}
		srTIgorInterf::ReleaseWave(p->wObjs, hState);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjM {
	DOUBLE obj_in;
    waveHndl wM;
	DOUBLE obj;
};
static int radObjM(void* p_void)
{
	radTIgorRadObjM *p = (radTIgorRadObjM*)p_void;

	int hState;
	double *pd=0;
	float *pf=0;
	bool DoubleArrIsLocal = false;
	try
	{
		long LenM=0;
		ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wM, pd, pf, LenM, hState));
		if(LenM < 6) throw TOO_SHORT_WAVE;

		int CntLen=0;
        ProcErrRad(RadObjCntSize(&CntLen, (int)(p->obj_in)));
		int NecessaryLenM = 6*CntLen;
		if(LenM < NecessaryLenM) throw TOO_SHORT_WAVE;

        if((pd == 0) && (pf != 0) && (NecessaryLenM > 0))
		{
            pd = new double[NecessaryLenM];
			DoubleArrIsLocal = true;
		}

        ProcErrRad(RadObjM(pd, (int)(p->obj_in)));

		if(DoubleArrIsLocal && (pd != 0)) 
		{
			float *tf = pf;
			double *td = pd;
			for(long i=0; i<NecessaryLenM; i++) *(tf++) = (float)(*(td++));
			delete[] pd;
			pd = 0;
		}
		
		p->obj = p->obj_in;
		ProcErr(srTIgorInterf::ReleaseWave(p->wM, hState));
	}
	catch(int ErrNo) 
	{
		srTIgorInterf::ReleaseWave(p->wM, hState);
		if(DoubleArrIsLocal && (pd != 0)) delete[] pd;
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjDpl {
    Handle hOpt1;
    DOUBLE obj_in;
	DOUBLE obj;
};
static int radObjDpl(void* p_void)
{
	try
	{
        radTIgorRadObjDpl *p = (radTIgorRadObjDpl*)p_void;
        p->obj = 0;

        char StrOpt1[31];
		*StrOpt1 = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt1, 30, StrOpt1));

		int Obj=0;
        ProcErrRad(RadObjDpl(&Obj, (int)(p->obj_in), StrOpt1));
        p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjCutMag {
    Handle hOpt1;
    waveHndl wN;
    waveHndl wP;
	DOUBLE obj_in;
    waveHndl wObjs;
	DOUBLE nobj;
};
static int radObjCutMag(void* p_void)
{
	float *pf=0;
	double *pd=0;
	int hState, *piObjs=0;
    radTIgorRadObjCutMag *p = (radTIgorRadObjCutMag*)p_void;
	try
	{
		p->nobj = 0;

		long nObjsAlloc=0;
		ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wObjs, pd, pf, nObjsAlloc, hState));
		if(nObjsAlloc < 2) throw TOO_SHORT_WAVE;

		piObjs = new int[nObjsAlloc];
		int *tiObjs = piObjs;
		if(pd != 0)
		{
			double *td = pd;
			for(int i=0; i<nObjsAlloc; i++) *(tiObjs++) = (int)(*(td++));
		}
		else if(pf != 0)
		{
			float *tf = pf;
			for(int i=0; i<nObjsAlloc; i++) *(tiObjs++) = (int)(*(tf++));
		}

		double P[3], N[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wN, N));

        char StrOpt1[61];
		*StrOpt1 = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt1, 60, StrOpt1));

		int nObj=0;
        ProcErrRad(RadObjCutMag(piObjs, &nObj, (int)(p->obj_in), P, N, StrOpt1));
		p->nobj = nObj;

		if((piObjs != 0) && (nObj > 0))
		{
			tiObjs = piObjs;
            if(pd != 0)
			{
                double *td = pd;
                for(int i=0; i<nObj; i++) *(td++) = *(tiObjs++);
				if(nObjsAlloc > nObj)
				{
                    for(int k=nObj; k<nObjsAlloc; k++) *(td++) = 0;
				}
            }
            else if(pf != 0)
            {
                float *tf = pf;
                for(int i=0; i<nObj; i++) *(tf++) = (float)(*(tiObjs++));
				if(nObjsAlloc > nObj)
				{
                    for(int k=nObj; k<nObjsAlloc; k++) *(tf++) = 0;
				}
            }

			delete[] piObjs;
			piObjs = 0;
		}
		ProcErr(srTIgorInterf::ReleaseWave(p->wObjs, hState));
	}
	catch(int ErrNo) 
	{
		if(piObjs != 0) { delete[] piObjs; piObjs = 0;}
		srTIgorInterf::ReleaseWave(p->wObjs, hState);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjDivMagPln {
    Handle hOpt1;
    waveHndl wNorm;
    waveHndl wSbdPar;
	DOUBLE obj_in;
	DOUBLE obj;
};
static int radObjDivMagPln(void* p_void)
{
	double *pSbdPar=0, *pNorm=0;

	try
	{
        radTIgorRadObjDivMagPln *p = (radTIgorRadObjDivMagPln*)p_void;
		p->obj = 0;

		int Dim1 = -1, Dim2 = -1;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wSbdPar, Dim1, Dim2, pSbdPar));
		if(Dim1 < 0) throw ZERO_NUMBER_OF_ELEM_IN_WAVE;
		if(Dim2 < 0) Dim2 = 1;
        int nSbdPar = Dim1*Dim2;

        Dim1 = 3, Dim2 = 3;
		ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wNorm, Dim1, Dim2, pNorm));

        char StrOpt1[1024]; //, StrOpt2[61];
		*StrOpt1 = '\0'; //*StrOpt2 = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt1, 1023, StrOpt1));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt2, 60, StrOpt2));

		int Obj=0;
        //ProcErrRad(RadObjDivMagPln(&Obj, (int)(p->obj_in), pSbdPar, nSbdPar, pNorm, StrOpt1, StrOpt2));
        ProcErrRad(RadObjDivMagPln(&Obj, (int)(p->obj_in), pSbdPar, nSbdPar, pNorm, StrOpt1));
		p->obj = Obj;

		if(pSbdPar != 0) { delete[] pSbdPar; pSbdPar = 0;}
		if(pNorm != 0) { delete[] pNorm; pNorm = 0;}
	}
	catch(int ErrNo) 
	{
		if(pSbdPar != 0) { delete[] pSbdPar; pSbdPar = 0;}
		if(pNorm != 0) { delete[] pNorm; pNorm = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjDivMagCyl {
    Handle hOpt1;
	DOUBLE rat;
    waveHndl wCylPar;
    waveHndl wSbdPar;
	DOUBLE obj_in;
	DOUBLE obj;
};
static int radObjDivMagCyl(void* p_void)
{
	double *pSbdPar=0, *pCylPar=0;

	try
	{
        radTIgorRadObjDivMagCyl *p = (radTIgorRadObjDivMagCyl*)p_void;
		p->obj = 0;

		int Dim1 = -1, Dim2 = -1;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wSbdPar, Dim1, Dim2, pSbdPar));
		if(Dim1 < 0) throw ZERO_NUMBER_OF_ELEM_IN_WAVE;
		if(Dim2 < 0) Dim2 = 1;
        int nSbdPar = Dim1*Dim2;

        Dim1 = 3, Dim2 = 3;
		ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wCylPar, Dim1, Dim2, pCylPar));

        char StrOpt1[1061]; //, StrOpt2[61];
		*StrOpt1 = '\0'; //*StrOpt2 = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt1, 1060, StrOpt1));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt2, 60, StrOpt2));

		int Obj=0;
        //ProcErrRad(RadObjDivMagCyl(&Obj, (int)(p->obj_in), pSbdPar, nSbdPar, pCylPar, (double)(p->rat), StrOpt1, StrOpt2));
        ProcErrRad(RadObjDivMagCyl(&Obj, (int)(p->obj_in), pSbdPar, nSbdPar, pCylPar, (double)(p->rat), StrOpt1));

		p->obj = Obj;

		if(pSbdPar != 0) { delete[] pSbdPar; pSbdPar = 0;}
		if(pCylPar != 0) { delete[] pCylPar; pCylPar = 0;}
	}
	catch(int ErrNo) 
	{
		if(pSbdPar != 0) { delete[] pSbdPar; pSbdPar = 0;}
		if(pCylPar != 0) { delete[] pCylPar; pCylPar = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjGeoVol {
	DOUBLE obj_in;
	DOUBLE v;
};
static int radObjGeoVol(void* p_void)
{
	try
	{
        radTIgorRadObjGeoVol *p = (radTIgorRadObjGeoVol*)p_void;
		p->v = 0;

		double V=0;
        ProcErrRad(RadObjGeoVol(&V, (int)(p->obj_in)));
		p->v = V;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjGeoLim {
	DOUBLE obj_in;
    waveHndl wLim;
	DOUBLE obj;
};
static int radObjGeoLim(void* p_void)
{
	try
	{
        radTIgorRadObjGeoLim *p = (radTIgorRadObjGeoLim*)p_void;

		double Lim[] = {0,0,0,0,0,0};
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wLim, Lim));
        ProcErrRad(RadObjGeoLim(Lim, (int)(p->obj_in)));

		int NumCols = 2;
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wLim, 3, NumCols));
        ProcErr(srTIgorInterf::SetDataInNumWave(p->wLim, Lim, 6, 0));

		p->obj = p->obj_in;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjDrwOpenGL {
    Handle hOpt1;
    DOUBLE obj;
	DOUBLE res;
};
static int radObjDrwOpenGL(void* p_void)
{
	try
	{
        radTIgorRadObjDrwOpenGL *p = (radTIgorRadObjDrwOpenGL*)p_void;

        char StrOpt1[1024];
		*StrOpt1 = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt1, 1023, StrOpt1));

        ProcErrRad(RadObjDrwOpenGL((int)(p->obj), StrOpt1));
        p->res = p->obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjDrwQD3D {
    Handle hOpt;
    DOUBLE obj;
	DOUBLE res;
};
static int radObjDrwQD3D(void* p_void)
{
	try
	{
        radTIgorRadObjDrwQD3D *p = (radTIgorRadObjDrwQD3D*)p_void;

        char StrOpt[1024];
		*StrOpt = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt, 1023, StrOpt));

        ProcErrRad(RadObjDrwQD3D((int)(p->obj), StrOpt));

        p->res = p->obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjDrwAtr {
    DOUBLE thcn;
    waveHndl wRGB;
	DOUBLE obj_in;
	DOUBLE obj;
};
static int radObjDrwAtr(void* p_void)
{
	try
	{
		radTIgorRadObjDrwAtr *p = (radTIgorRadObjDrwAtr*)p_void;
		p->obj = 0;
		double RGB[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wRGB, RGB));

		int Obj=0;
		ProcErrRad(RadObjDrwAtr((int)(p->obj_in), RGB, (double)(p->thcn)));
		p->obj = p->obj_in;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadObjFullMag {
	waveHndl wRGB;
	DOUBLE mat;
	DOUBLE grp;
	waveHndl wK;
	waveHndl wM;
    waveHndl wL;
    waveHndl wP;
	DOUBLE obj;
};
static int radObjFullMag(void* p_void)
{
	double* pSbdPar = 0;
	try
	{
		radTIgorRadObjFullMag *p = (radTIgorRadObjFullMag*)p_void;
		p->obj = 0;

		double P[3], L[3], M[3], RGB[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wL, L));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wM, M));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wRGB, RGB));

		int Dim1 = -1, Dim2 = -1;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wK, Dim1, Dim2, pSbdPar));
		if(Dim1 < 0) throw ZERO_NUMBER_OF_ELEM_IN_WAVE;
		if(Dim2 < 0) Dim2 = 1;
        int nSbdPar = Dim1*Dim2;

		int Obj=0;
		ProcErrRad(RadObjFullMag(&Obj, P, L, M, pSbdPar, nSbdPar, (int)(p->grp), (int)(p->mat), RGB));
		p->obj = Obj;

		if(pSbdPar != 0) { delete[] pSbdPar; pSbdPar = 0;}
	}
	catch(int ErrNo) 
	{
		if(pSbdPar != 0) { delete[] pSbdPar; pSbdPar = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfPlSym {
    waveHndl wN;
    waveHndl wP;
	DOUBLE obj;
};
static int radTrfPlSym(void* p_void)
{
	try
	{
		radTIgorRadTrfPlSym *p = (radTIgorRadTrfPlSym*)p_void;
		p->obj = 0;
		double P[3], N[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wN, N));
		int Obj=0;
		ProcErrRad(RadTrfPlSym(&Obj, P, N));
		p->obj = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfRot {
	DOUBLE phi;
    waveHndl wV;
    waveHndl wP;
	DOUBLE trf;
};
static int radTrfRot(void* p_void)
{
	try
	{
		radTIgorRadTrfRot *p = (radTIgorRadTrfRot*)p_void;
		p->trf = 0;
		double P[3], V[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wV, V));

		int Trf=0;
		ProcErrRad(RadTrfRot(&Trf, P, V, (double)(p->phi)));
		p->trf = Trf;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfTrsl {
    waveHndl wV;
	DOUBLE trf;
};
static int radTrfTrsl(void* p_void)
{
	try
	{
		radTIgorRadTrfTrsl *p = (radTIgorRadTrfTrsl*)p_void;
		p->trf = 0;
		double V[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wV, V));

		int Trf=0;
		ProcErrRad(RadTrfTrsl(&Trf, V));
		p->trf = Trf;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfInv {
	DOUBLE trf;
};
static int radTrfInv(void* p_void)
{
	try
	{
		radTIgorRadTrfInv *p = (radTIgorRadTrfInv*)p_void;
		p->trf = 0;

		int Trf=0;
		ProcErrRad(RadTrfInv(&Trf));
		p->trf = Trf;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfCmbL {
	DOUBLE trf;
	DOUBLE origtrf;
	DOUBLE newtrf;
};
static int radTrfCmbL(void* p_void)
{
	try
	{
		radTIgorRadTrfCmbL *p = (radTIgorRadTrfCmbL*)p_void;
		p->newtrf = 0;

		int Trf=0;
		ProcErrRad(RadTrfCmbL(&Trf, (int)(p->origtrf), (int)(p->trf)));
		p->newtrf = Trf;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfCmbR {
	DOUBLE trf;
	DOUBLE origtrf;
	DOUBLE newtrf;
};
static int radTrfCmbR(void* p_void)
{
	try
	{
		radTIgorRadTrfCmbR *p = (radTIgorRadTrfCmbR*)p_void;
		p->newtrf = 0;

		int Trf=0;
		ProcErrRad(RadTrfCmbR(&Trf, (int)(p->origtrf), (int)(p->trf)));
		p->newtrf = Trf;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfMlt {
	DOUBLE mlt;
	DOUBLE trf;
	DOUBLE obj;
	DOUBLE objout;
};
static int radTrfMlt(void* p_void)
{
	try
	{
		radTIgorRadTrfMlt *p = (radTIgorRadTrfMlt*)p_void;
		p->objout = 0;

		int ObjOut=0;
		ProcErrRad(RadTrfMlt(&ObjOut, (int)(p->obj), (int)(p->trf), (int)(p->mlt)));
		p->objout = ObjOut;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfOrnt {
	DOUBLE trf;
	DOUBLE obj;
	DOUBLE objout;
};
static int radTrfOrnt(void* p_void)
{
	try
	{
		radTIgorRadTrfOrnt *p = (radTIgorRadTrfOrnt*)p_void;
		p->objout = 0;

		int ObjOut=0;
		ProcErrRad(RadTrfOrnt(&ObjOut, (int)(p->obj), (int)(p->trf)));
		p->objout = ObjOut;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfZerPara {
    waveHndl wN;
    waveHndl wP;
	DOUBLE obj;
	DOUBLE objout;
};
static int radTrfZerPara(void* p_void)
{
	try
	{
		radTIgorRadTrfZerPara *p = (radTIgorRadTrfZerPara*)p_void;
		p->objout = 0;
		double P[3], N[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wN, N));

		int Obj=0;
		ProcErrRad(RadTrfZerPara(&Obj, (int)(p->obj), P, N));
		p->objout = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadTrfZerPerp {
    waveHndl wN;
    waveHndl wP;
	DOUBLE obj;
	DOUBLE objout;
};
static int radTrfZerPerp(void* p_void)
{
	try
	{
		radTIgorRadTrfZerPerp *p = (radTIgorRadTrfZerPerp*)p_void;
		p->objout = 0;
		double P[3], N[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wN, N));

		int Obj=0;
		ProcErrRad(RadTrfZerPerp(&Obj, (int)(p->obj), P, N));
		p->objout = Obj;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadMatLin {
    //DOUBLE nMr;
    waveHndl wMr;
	waveHndl wKsi;
	DOUBLE mat;
};
static int radMatLin(void* p_void)
{
	try
	{
		radTIgorRadMatLin *p = (radTIgorRadMatLin*)p_void;
		p->mat = 0;
		double Ksi[2], Mr[3];
		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wKsi, Ksi));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wMr, Mr));

		long nMr = -1;
		double *pMr = Mr;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave1D(p->wMr, 3, pMr, nMr));

		int Mat=0;
		ProcErrRad(RadMatLin(&Mat, Ksi, Mr, (int)(nMr)));
		p->mat = Mat;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadMatSatIso {
    //DOUBLE np;
    Handle hDataType;
	waveHndl wMatData;
	DOUBLE mat;
};
static int radMatSatIso(void* p_void)
{
	double* MatData=0;

	try
	{
		radTIgorRadMatSatIso *p = (radTIgorRadMatSatIso*)p_void;
		p->mat = 0;

		//int NpDeclar = (int)(p->np);
        //int Np = 0;
		//ProcErr(srTIgorInterf::GetArrDoubleAndNpFromNumWave1D(p->wMatData, NpDeclar, 2, MatData, Np));
		//if(NpDeclar > Np) throw TOO_SHORT_WAVE;

		int Dim1 = 2, Dim2 = -1;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wMatData, Dim1, Dim2, MatData));
		if(Dim2 <= 0) throw ZERO_NUMBER_OF_ELEM_IN_WAVE;
		int Np = Dim2;

        char StrDataType[1025];
		*StrDataType = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hDataType, 1024, StrDataType));
		CAuxParse::toUpperCase(StrDataType);

		int Mat=0;
		if(strcmp(StrDataType, "TAB") == 0) ProcErrRad(RadMatSatIsoTab(&Mat, MatData, Np));
		else if(strcmp(StrDataType, "FRM") == 0)
		{
			double KsiMs1[] = {MatData[0], MatData[1]}, KsiMs2[] = {0,0}, KsiMs3[] = {0,0};
			if(Np > 1) { KsiMs2[0] = MatData[2]; KsiMs2[1] = MatData[3];}
			if(Np > 2) { KsiMs3[0] = MatData[4]; KsiMs3[1] = MatData[5];}
            ProcErrRad(RadMatSatIsoFrm(&Mat, KsiMs1, KsiMs2, KsiMs3));
		}
		p->mat = Mat;

		if(MatData != 0) { delete[] MatData; MatData = 0;}
	}
	catch(int ErrNo) 
	{
		if(MatData != 0) { delete[] MatData; MatData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadMatSatLam {
	waveHndl wN;
    DOUBLE stack;
    Handle hDataType;
	waveHndl wMatData;
	DOUBLE mat;
};
static int radMatSatLam(void* p_void)
{
	double* MatData=0;

	try
	{
		radTIgorRadMatSatLam *p = (radTIgorRadMatSatLam*)p_void;
		p->mat = 0;

		int Dim1 = 2, Dim2 = -1;
        ProcErr(srTIgorInterf::GetArrDoubleFromNumWave2D(p->wMatData, Dim1, Dim2, MatData));
		if(Dim2 <= 0) throw ZERO_NUMBER_OF_ELEM_IN_WAVE;
		int Np = Dim2;

        char StrDataType[1025];
		*StrDataType = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hDataType, 1024, StrDataType));
		CAuxParse::toUpperCase(StrDataType);

		double N[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wN, N));

		int Mat=0;
		if(strcmp(StrDataType, "TAB") == 0) ProcErrRad(RadMatSatLamTab(&Mat, MatData, Np, (double)(p->stack), N));
		else if(strcmp(StrDataType, "FRM") == 0)
		{
			double KsiMs1[] = {MatData[0], MatData[1]}, KsiMs2[] = {0,0}, KsiMs3[] = {0,0};
			if(Np > 1) { KsiMs2[0] = MatData[2]; KsiMs2[1] = MatData[3];}
			if(Np > 2) { KsiMs3[0] = MatData[4]; KsiMs3[1] = MatData[5];}
            ProcErrRad(RadMatSatLamFrm(&Mat, KsiMs1, KsiMs2, KsiMs3, (double)(p->stack), N));
		}
		p->mat = Mat;

		if(MatData != 0) { delete[] MatData; MatData = 0;}
	}
	catch(int ErrNo) 
	{
		if(MatData != 0) { delete[] MatData; MatData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadMatSatAniso {
    //DOUBLE nDataPer;
	waveHndl wDataPer;
    //DOUBLE nDataPar;
	waveHndl wDataPar;
	DOUBLE mat;
};
static int radMatSatAniso(void* p_void)
{
	double *DataPar=0, *DataPer=0;
    long Np = 0;

	try
	{
		radTIgorRadMatSatAniso *p = (radTIgorRadMatSatAniso*)p_void;
		p->mat = 0;

		//int nDataParDeclar = (int)(p->nDataPar);
		ProcErr(srTIgorInterf::GetArrDoubleFromNumWave1D(p->wDataPar, -1, DataPar, Np));
		//if(nDataParDeclar > Np) throw TOO_SHORT_WAVE;
		if(Np <= 0) throw ZERO_NUMBER_OF_ELEM_IN_WAVE;
		int nDataPar = (int)Np;

		//int nDataPerDeclar = (int)(p->nDataPer);
		ProcErr(srTIgorInterf::GetArrDoubleFromNumWave1D(p->wDataPer, -1, DataPer, Np));
		//if(nDataPerDeclar > Np) throw TOO_SHORT_WAVE;
		if(Np <= 0) throw ZERO_NUMBER_OF_ELEM_IN_WAVE;
		int nDataPer = (int)Np;

		int Mat=0;
		ProcErrRad(RadMatSatAniso(&Mat, DataPar, nDataPar, DataPer, nDataPer));
		p->mat = Mat;

		if(DataPar != 0) { delete[] DataPar; DataPar = 0;}
		if(DataPer != 0) { delete[] DataPer; DataPer = 0;}
	}
	catch(int ErrNo) 
	{
		if(DataPar != 0) { delete[] DataPar; DataPar = 0;}
		if(DataPer != 0) { delete[] DataPer; DataPer = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadMatStd {
    DOUBLE mr;
    Handle id;
    DOUBLE mat;
};
static int radMatStd(void* p_void)
{
	try
	{
        radTIgorRadMatStd *p = (radTIgorRadMatStd*)p_void;

        char MatID[2048];
		*MatID = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->id, 1024, MatID));

		int Mat=0;
        ProcErrRad(RadMatStd(&Mat, MatID, (double)(p->mr)));
		p->mat = Mat;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadMatApl {
	DOUBLE mat;
	DOUBLE obj;
	DOUBLE objout;
};
static int radMatApl(void* p_void)
{
	try
	{
		radTIgorRadMatApl *p = (radTIgorRadMatApl*)p_void;
		p->objout = 0;
		int ObjOut=0;
		ProcErrRad(RadMatApl(&ObjOut, (int)(p->obj), (int)(p->mat)));
		p->objout = ObjOut;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadMatMvsH {
    waveHndl wH;
	Handle id;
	DOUBLE obj;
    waveHndl wM;
	DOUBLE M0;
};
static int radMatMvsH(void* p_void)
{
	radTIgorRadMatMvsH *p = (radTIgorRadMatMvsH*)p_void;

	//int hState;
	//double *pdM=0;
	//float *pfM=0;
	//bool DoubleArrMIsLocal = false;
	try
	{
		double ComponH[3];
        ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wH, ComponH));

		//p->nM = 0;
		//long LenM=0;
		//ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wM, pdM, pfM, LenM, hState));
		//if((pdM == 0) && (pfM != 0) && (LenM > 0))
		//{
        //	pdM = new double[LenM];
		//	DoubleArrMIsLocal = true;
		//}

        char StrID[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->id, 30, StrID));

		p->M0 = 0;
		double TempCalcData[10];
        int Nm=0;
        ProcErrRad(RadMatMvsH(TempCalcData, &Nm, (int)(p->obj), StrID, ComponH));

		int NumCols = 0;
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wM, Nm, NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wM, TempCalcData, Nm, "T"));

		//if(DoubleArrMIsLocal && (pdM != 0)) 
		//{
		//	float *tfM = pfM;
		//	double *tdM = pdM;
		//	for(long i=0; i<Nm; i++) *(tfM++) = (float)(*(tdM++));
		//	delete[] pdM; pdM = 0;
		//}

		//p->nM = Nm;
		p->M0 = TempCalcData[0];

		//ProcErr(srTIgorInterf::ReleaseWave(p->wM, hState));
	}
	catch(int ErrNo) 
	{
		//srTIgorInterf::ReleaseWave(p->wM, hState);
		//if(DoubleArrMIsLocal && (pdM != 0)) delete[] pdM;
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadRlxPre {
	DOUBLE objsrc;
	DOUBLE obj;
	DOUBLE objout;
};
static int radRlxPre(void* p_void)
{
	try
	{
		radTIgorRadRlxPre *p = (radTIgorRadRlxPre*)p_void;
		p->objout = 0;
		int ObjOut=0;
		ProcErrRad(RadRlxPre(&ObjOut, (int)(p->obj), (int)(p->objsrc)));
		p->objout = ObjOut;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadRlxMan {
	DOUBLE rlxpar;
	DOUBLE iter;
	DOUBLE meth;
	DOUBLE intrc;
	waveHndl wD;
	DOUBLE nD;
};
static int radRlxMan(void* p_void)
{
	try
	{
		radTIgorRadRlxMan *p = (radTIgorRadRlxMan*)p_void;

		double ResArr[4];
		int nD;
		ProcErrRad(RadRlxMan(ResArr, &nD, (int)(p->intrc), (int)(p->meth), (int)(p->iter), (double)(p->rlxpar)));
        //p->nD = nD;
		p->nD = (int)(p->iter);

		for(int i=0; i<nD; i++) ProcErrRad(srTIgorInterf::SetNumWavePointValue1D(p->wD, i, ResArr[i]));
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadRlxAuto {
	DOUBLE meth;
	DOUBLE iter;
	DOUBLE prec;
	DOUBLE intrc;
	waveHndl wD;
	DOUBLE nD; //number of iterations made
};
static int radRlxAuto(void* p_void)
{
	try
	{
		radTIgorRadRlxAuto *p = (radTIgorRadRlxAuto*)p_void;

		double ResArr[4];
		int nD;
		ProcErrRad(RadRlxAuto(ResArr, &nD, (int)(p->intrc), (double)(p->prec), (int)(p->iter), (int)(p->meth), 0));

		//p->nD = nD;
		p->nD = (long)(ResArr[3]);

		for(int i=0; i<nD; i++) ProcErrRad(srTIgorInterf::SetNumWavePointValue1D(p->wD, i, ResArr[i]));
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadRlxAutoExt {
	Handle hOpt1; //relaxation options
	DOUBLE meth;
	DOUBLE iter;
	DOUBLE prec;
	DOUBLE intrc;
	waveHndl wD;
	DOUBLE nD; //output: number of iterations made
};
static int radRlxAutoExt(void* p_void)
{
	try
	{
		radTIgorRadRlxAutoExt *p = (radTIgorRadRlxAutoExt*)p_void;

		char sOpt1[1024];
		*sOpt1 = '\0';
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt1, 1023, sOpt1));

		double ResArr[4];
		int nD;
		ProcErrRad(RadRlxAuto(ResArr, &nD, (int)(p->intrc), (double)(p->prec), (int)(p->iter), (int)(p->meth), sOpt1));

		p->nD = (long)(ResArr[3]);

		for(int i=0; i<nD; i++) ProcErrRad(srTIgorInterf::SetNumWavePointValue1D(p->wD, i, ResArr[i]));
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadSolve {
	DOUBLE meth;
	DOUBLE iter;
	DOUBLE prec;
	DOUBLE obj;
	waveHndl wD;
	DOUBLE nD; //number of iterations made
};
static int radSolve(void* p_void)
{
	try
	{
		radTIgorRadSolve *p = (radTIgorRadSolve*)p_void;

		double ResArr[4];
		int nD;
		ProcErrRad(RadSolve(ResArr, &nD, (int)(p->obj), (double)(p->prec), (int)(p->iter), (int)(p->meth)));
        //p->nD = nD;
        p->nD = (long)(ResArr[3]);

        for(int i=0; i<nD; i++) ProcErrRad(srTIgorInterf::SetNumWavePointValue1D(p->wD, i, ResArr[i]));
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFld {
    waveHndl wCoord;
    Handle hID;
	DOUBLE obj;
    waveHndl wB;
	DOUBLE B0;
};
static int radFld(void* p_void)
{
	radTIgorRadFld *p = (radTIgorRadFld*)p_void;
    double *pTempCalcData = 0;
	double* pCoord=0;
	try
	{
        int Np = 0;
		ProcErr(srTIgorInterf::GetArrDoubleFromNumWaveN3(p->wCoord, pCoord, Np));
		if(Np < 1) throw THREE_ELEM_NUM_WAVE_REQUIRED;

        char StrID[256];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 256, StrID));

		const int MaxAmOfFieldValPerPoint = 30;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

		p->B0 = 0;
		int Nb=0; //total number of values produced
        ProcErrRad(RadFld(pTempCalcData, &Nb, (int)(p->obj), StrID, pCoord, Np));
		p->B0 = *pTempCalcData;

		//int NumCols = (int)(((double)Nb)/((double)Np) + 1E-06);
		//if(NumCols <= 1) NumCols = 0;
        //ProcErr(srTIgorInterf::ReDimNumWave2D(p->wB, Np, NumCols));

		int NumRows=0, NumCols=0;
		if(Nb <= Np)
		{
			NumRows = Np;
		}
		else if(Np == 1)
		{
			NumRows = Nb;
		}
		else
		{
            NumRows = (int)(((double)Nb)/((double)Np) + 1E-06);
			NumCols = Np;
		}
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wB, NumRows, NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wB, pTempCalcData, Nb, "T"));

		double ArgStarts[] = {0, 0}, ArgSteps[] = {1, 1};
		char* ArgValUnits[] = {"", "", "T"};
        ProcErr(srTIgorInterf::SetScaleNumWave(p->wB, 2, ArgStarts, ArgSteps, ArgValUnits));

		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
	}
	catch(int ErrNo) 
	{
		if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldC {
    //waveHndl wCoord;
	DOUBLE z;
	DOUBLE y;
	DOUBLE x;
    Handle hID;
	DOUBLE obj;
    waveHndl wB;
	DOUBLE B0;
};
static int radFldC(void* p_void)
{
	radTIgorRadFldC *p = (radTIgorRadFldC*)p_void;
    double *pTempCalcData = 0;
	//double* pCoord=0;
	try
	{
		//int Np = 0;
		//ProcErr(srTIgorInterf::GetArrDoubleFromNumWaveN3(p->wCoord, pCoord, Np));
		//if(Np < 1) throw THREE_ELEM_NUM_WAVE_REQUIRED;

		int Np = 1;
		double pCoord[] = {p->x, p->y, p->z};

        char StrID[256];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 256, StrID));

		const int MaxAmOfFieldValPerPoint = 30;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

		p->B0 = 0;
		int Nb=0; //total number of values produced
        ProcErrRad(RadFld(pTempCalcData, &Nb, (int)(p->obj), StrID, pCoord, Np));
		p->B0 = *pTempCalcData;

		int NumRows=0, NumCols=0;
		if(Nb <= Np)
		{
			NumRows = Np;
		}
		else if(Np == 1)
		{
			NumRows = Nb;
		}
		else
		{
            NumRows = (int)(((double)Nb)/((double)Np) + 1E-06);
			NumCols = Np;
		}
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wB, NumRows, NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wB, pTempCalcData, Nb, "T"));

		double ArgStarts[] = {0, 0}, ArgSteps[] = {1, 1};
		char* ArgValUnits[] = {"", "", "T"};
        ProcErr(srTIgorInterf::SetScaleNumWave(p->wB, 2, ArgStarts, ArgSteps, ArgValUnits));

		//if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
	}
	catch(int ErrNo) 
	{
		//if(pCoord != 0) { delete[] pCoord; pCoord = 0;}
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldLst {
	DOUBLE start;
	//Handle hArgOrNot;
	DOUBLE np;
    waveHndl wP2;
    waveHndl wP1;
    Handle hID;
	DOUBLE obj;
    waveHndl wB;
	DOUBLE B0;
};
static int radFldLst(void* p_void)
{
	radTIgorRadFldLst *p = (radTIgorRadFldLst*)p_void;
    double *pTempCalcData = 0;

	try
	{
		double P1[3], P2[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP1, P1));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP2, P2));

		p->B0 = 0;
		int Np = (int)(p->np);

        const int MaxAmOfFieldValPerPoint = 30;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        char StrID[31], StrArgOrNot[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 30, StrID));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hArgOrNot, 30, StrArgOrNot));
		strcpy(StrArgOrNot, "noarg");

        int Nb=0;
        ProcErrRad(RadFldLst(pTempCalcData, &Nb, (int)(p->obj), StrID, P1, P2, Np, StrArgOrNot, (double)(p->start)));
		p->B0 = *pTempCalcData;

		//int NumCols = (int)(((double)Nb)/((double)(p->np)) + 1E-06);
		//if(NumCols <= 1) NumCols = 0;
        //ProcErr(srTIgorInterf::ReDimNumWave2D(p->wB, (int)(p->np), NumCols));

		double dx = P2[0] - P1[0], dy = P2[1] - P1[1], dz = P2[2] - P1[2];
		double ArgRange = sqrt(dx*dx + dy*dy + dz*dz);
		double ArgStep = (Np <= 1)? 0 : ArgRange/(Np - 1);
		double ArgStarts[] = {0, 0}, ArgSteps[] = {1, 1};
		char* ArgValUnits[] = {"", "", "T"};

		//double ArgStarts[] = {(p->start)*0.001, 0}, ArgSteps[] = {ArgStep*0.001, 1};

		int NumRows=0, NumCols=0, IndDim = 0;
		if(Nb <= Np)
		{
			NumRows = Np;
			IndDim = 0;
		}
		else if(Np == 1)
		{
			NumRows = Nb;
			IndDim = 0;
		}
		else
		{
            NumRows = (int)(((double)Nb)/((double)Np) + 1E-06);
			NumCols = Np;
			IndDim = 1;
		}
        ArgStarts[IndDim] = (p->start)*0.001;
		ArgSteps[IndDim] = ArgStep*0.001;
		ArgValUnits[IndDim] = "m";

        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wB, NumRows, NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wB, pTempCalcData, Nb, "T"));
        ProcErr(srTIgorInterf::SetScaleNumWave(p->wB, 2, ArgStarts, ArgSteps, ArgValUnits));

		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
	}
	catch(int ErrNo) 
	{
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldLstC {
	DOUBLE start;
	DOUBLE np;
	DOUBLE z2;
	DOUBLE y2;
	DOUBLE x2;
	DOUBLE z1;
	DOUBLE y1;
	DOUBLE x1;
    Handle hID;
	DOUBLE obj;
    waveHndl wB;
	DOUBLE B0;
};
static int radFldLstC(void* p_void)
{
	radTIgorRadFldLstC *p = (radTIgorRadFldLstC*)p_void;
    double *pTempCalcData = 0;

	try
	{
		//double P1[3], P2[3];
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP1, P1));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP2, P2));

		double P1[] = {p->x1, p->y1, p->z1};
		double P2[] = {p->x2, p->y2, p->z2};

		p->B0 = 0;
		int Np = (int)(p->np);

        const int MaxAmOfFieldValPerPoint = 30;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        char StrID[31], StrArgOrNot[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 30, StrID));
		strcpy(StrArgOrNot, "noarg");

        int Nb=0;
        ProcErrRad(RadFldLst(pTempCalcData, &Nb, (int)(p->obj), StrID, P1, P2, Np, StrArgOrNot, (double)(p->start)));
		p->B0 = *pTempCalcData;

		double dx = P2[0] - P1[0], dy = P2[1] - P1[1], dz = P2[2] - P1[2];
		double ArgRange = sqrt(dx*dx + dy*dy + dz*dz);
		double ArgStep = (Np <= 1)? 0 : ArgRange/(Np - 1);
		double ArgStarts[] = {0, 0}, ArgSteps[] = {1, 1};
		char* ArgValUnits[] = {"", "", "T"};

		int NumRows=0, NumCols=0, IndDim = 0;
		if(Nb <= Np)
		{
			NumRows = Np;
			IndDim = 0;
		}
		else if(Np == 1)
		{
			NumRows = Nb;
			IndDim = 0;
		}
		else
		{
            NumRows = (int)(((double)Nb)/((double)Np) + 1E-06);
			NumCols = Np;
			IndDim = 1;
		}
        ArgStarts[IndDim] = (p->start)*0.001;
		ArgSteps[IndDim] = ArgStep*0.001;
		ArgValUnits[IndDim] = "m";

        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wB, NumRows, NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wB, pTempCalcData, Nb, "T"));
        ProcErr(srTIgorInterf::SetScaleNumWave(p->wB, 2, ArgStarts, ArgSteps, ArgValUnits));

		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
	}
	catch(int ErrNo) 
	{
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldInt {
    waveHndl wP2;
    waveHndl wP1;
    Handle hID;
    Handle hInfOrFin;
	DOUBLE obj;
    waveHndl wIb;
	DOUBLE Ib0; // first component of field integral calculated
};
static int radFldInt(void* p_void)
{
	radTIgorRadFldInt *p = (radTIgorRadFldInt*)p_void;
    double *pTempCalcData = 0;

	//int hState;
	//double *pdIb=0;
	//float *pfIb =0;
	//bool DoubleArrBIsLocal = false;

	try
	{
		double P1[3], P2[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP1, P1));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP2, P2));

		p->Ib0 = 0;
		//long LenB=0;
		//ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wIb, pdIb, pfIb, LenB, hState));
		//if((pdIb == 0) && (pfIb != 0) && (LenB > 0))
		//{
        //	pdIb = new double[LenB];
		//	DoubleArrBIsLocal = true;
		//}

		int Np = 1;
        const int MaxAmOfFieldValPerPoint = 30;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        char StrID[31], StrInfOrFin[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 30, StrID));
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hInfOrFin, 30, StrInfOrFin));

        int NIb=0;
        //ProcErrRad(RadFldInt(pdIb, &NIb, (int)(p->obj), StrInfOrFin, StrID, P1, P2));
		//p->Ib0 = *pdIb;

        ProcErrRad(RadFldInt(pTempCalcData, &NIb, (int)(p->obj), StrInfOrFin, StrID, P1, P2));
		p->Ib0 = *pTempCalcData;

		//int NumCols = (int)(((double)Nb)/((double)(p->np)) + 1E-06);
		//if(NumCols <= 1) NumCols = 0;
		int NumCols = 0;

        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wIb, NIb, NumCols));
        ProcErr(srTIgorInterf::SetDataInNumWave(p->wIb, pTempCalcData, NIb, "T*mm"));
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}

		//if(DoubleArrBIsLocal && (pdIb != 0)) 
		//{
		//	float *tfIb = pfIb;
		//	double *tdIb = pdIb;
		//	for(long i=0; i<NIb; i++) *(tfIb++) = (float)(*(tdIb++));
		//	delete[] pdIb;
		//	pdIb = 0;
		//}
		//ProcErr(srTIgorInterf::ReleaseWave(p->wIb, hState));
	}
	catch(int ErrNo) 
	{
		//srTIgorInterf::ReleaseWave(p->wIb, hState);
		//if(DoubleArrBIsLocal && (pdIb != 0)) delete[] pdIb;

        if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldIntC {
	DOUBLE z2;
	DOUBLE y2;
	DOUBLE x2;
	DOUBLE z1;
	DOUBLE y1;
	DOUBLE x1;
    Handle hID;
    Handle hInfOrFin;
	DOUBLE obj;
    waveHndl wIb;
	DOUBLE Ib0; // first component of field integral calculated
};
static int radFldIntC(void* p_void)
{
	radTIgorRadFldIntC *p = (radTIgorRadFldIntC*)p_void;
    double *pTempCalcData = 0;

	try
	{
		double P1[] = {p->x1, p->y1, p->z1};
		double P2[] = {p->x2, p->y2, p->z2};
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP1, P1));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP2, P2));

		p->Ib0 = 0;

		int Np = 1;
        const int MaxAmOfFieldValPerPoint = 30;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        char StrID[31], StrInfOrFin[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 30, StrID));
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hInfOrFin, 30, StrInfOrFin));

        int NIb=0;
        ProcErrRad(RadFldInt(pTempCalcData, &NIb, (int)(p->obj), StrInfOrFin, StrID, P1, P2));
		p->Ib0 = *pTempCalcData;

		int NumCols = 0;
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wIb, NIb, NumCols));
        ProcErr(srTIgorInterf::SetDataInNumWave(p->wIb, pTempCalcData, NIb, "T*mm"));
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
	}
	catch(int ErrNo) 
	{
        if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldPtcTrj {
	DOUBLE np;
    waveHndl wLongLim;
    waveHndl wInitCond;
	DOUBLE E;
	DOUBLE obj;
    waveHndl wF;
	DOUBLE nF; //total number of trajectory data components calculated
};
static int radFldPtcTrj(void* p_void)
{
	radTIgorRadFldPtcTrj *p = (radTIgorRadFldPtcTrj*)p_void;
    double *pTempCalcData = 0;

	//int hState;
	//double *pdF=0;
	//float *pfF =0;
	//bool DoubleArrBIsLocal = false;

	try
	{
		double InitCond[4], LongLim[2];
		ProcErr(srTIgorInterf::Get4ElemArrDoubleFromNumWave1D(p->wInitCond, InitCond));
		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wLongLim, LongLim));

		p->nF = 0;
		//long LenF=0;
		//ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wF, pdF, pfF, LenF, hState));
		//if((pdF == 0) && (pfF != 0) && (LenF > 0))
		//{
        //	pdF = new double[LenF];
		//	DoubleArrBIsLocal = true;
		//}

		int Np = (int)(p->np);
        const int MaxAmOfFieldValPerPoint = 10;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        int NF=0;
        //ProcErrRad(RadFldPtcTrj(pdF, &NF, (int)(p->obj), (double)(p->E), InitCond, LongLim, (int)(p->np)));
        ProcErrRad(RadFldPtcTrj(pTempCalcData, &NF, (int)(p->obj), (double)(p->E), InitCond, LongLim, (int)(p->np)));
		p->nF = NF;

		int NumCols = (int)(((double)NF)/((double)(p->np)) + 1E-06);
		if(NumCols <= 1) NumCols = 0;
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wF, (int)(p->np), NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wF, pTempCalcData, NF, 0));
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}

		//if(DoubleArrBIsLocal && (pdF != 0)) 
		//{
		//	float *tfF = pfF;
		//	double *tdF = pdF;
		//	for(long i=0; i<NF; i++) *(tfF++) = (float)(*(tdF++));
		//	delete[] pdF;
		//	pdF = 0;
		//}
		//ProcErr(srTIgorInterf::ReleaseWave(p->wF, hState));
	}
	catch(int ErrNo) 
	{
		//srTIgorInterf::ReleaseWave(p->wF, hState);
		//if(DoubleArrBIsLocal && (pdF != 0)) delete[] pdF;

		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldPtcTrjC {
	DOUBLE np;
    //waveHndl wLongLim;
	DOUBLE y2;
	DOUBLE y1;
	//waveHndl wInitCond;
	DOUBLE dzdy0;
	DOUBLE z0;
	DOUBLE dxdy0;
	DOUBLE x0;
	DOUBLE E;
	DOUBLE obj;
    waveHndl wF;
	DOUBLE nF; //total number of trajectory data components calculated
};
static int radFldPtcTrjC(void* p_void)
{
	radTIgorRadFldPtcTrjC *p = (radTIgorRadFldPtcTrjC*)p_void;
    double *pTempCalcData = 0;

	try
	{
		//double InitCond[4], LongLim[2];
		//ProcErr(srTIgorInterf::Get4ElemArrDoubleFromNumWave1D(p->wInitCond, InitCond));
		//ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wLongLim, LongLim));

		double InitCond[] = {p->x0, p->dxdy0, p->z0, p->dzdy0};
		double LongLim[] = {p->y1, p->y2};

		p->nF = 0;
		int Np = (int)(p->np);
        const int MaxAmOfFieldValPerPoint = 10;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        int NF=0;
        ProcErrRad(RadFldPtcTrj(pTempCalcData, &NF, (int)(p->obj), (double)(p->E), InitCond, LongLim, (int)(p->np)));
		p->nF = NF;

		int NumCols = (int)(((double)NF)/((double)(p->np)) + 1E-06);
		if(NumCols <= 1) NumCols = 0;
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wF, (int)(p->np), NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wF, pTempCalcData, NF, 0));
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
	}
	catch(int ErrNo) 
	{
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldFocPot {
	DOUBLE np;
    waveHndl wP2;
    waveHndl wP1;
	DOUBLE obj;
	DOUBLE d;
};
static int radFldFocPot(void* p_void)
{
	try
	{
        radTIgorRadFldFocPot *p = (radTIgorRadFldFocPot*)p_void;

		double P1[3], P2[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP1, P1));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP2, P2));

		double D;
		p->d = 0;
        ProcErrRad(RadFldFocPot(&D, (int)(p->obj), P1, P2, (int)(p->np)));
		p->d = D;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldFocPotC {
	DOUBLE np;
    //waveHndl wP2;
	DOUBLE z2;
	DOUBLE y2;
	DOUBLE x2;
    //waveHndl wP1;
	DOUBLE z1;
	DOUBLE y1;
	DOUBLE x1;
	DOUBLE obj;
	DOUBLE d;
};
static int radFldFocPotC(void* p_void)
{
	try
	{
        radTIgorRadFldFocPotC *p = (radTIgorRadFldFocPotC*)p_void;

		//double P1[3], P2[3];
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP1, P1));
		//ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP2, P2));
		double P1[] = {p->x1, p->y1, p->z1};
		double P2[] = {p->x2, p->y2, p->z2};

		double D;
		p->d = 0;
        ProcErrRad(RadFldFocPot(&D, (int)(p->obj), P1, P2, (int)(p->np)));
		p->d = D;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldFocKickPer {
	Handle hCom;
	DOUBLE nh;
	DOUBLE d2;
	DOUBLE np2;
	DOUBLE r2;
	DOUBLE d1;
	DOUBLE np1;
	DOUBLE r1;
    waveHndl wNtr;
	DOUBLE nps;
	DOUBLE nper;
	DOUBLE per;
    waveHndl wNs;
    waveHndl wP1;
	DOUBLE obj;
    waveHndl wArg2;
    waveHndl wArg1;
    waveHndl wIBe2;
    waveHndl wM2;
    waveHndl wM1;
	Handle hFormStr;
};
static int radFldFocKickPer(void* p_void)
{
    double *pTempM1=0, *pTempM2=0, *pTempIBe2=0, *pTempArg1=0, *pTempArg2=0;
	char *pTempStr=0;
    radTIgorRadFldFocKickPer *p = (radTIgorRadFldFocKickPer*)p_void;

	try
	{
		double P1[3], Ns[3], Ntr[3];
        ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP1, P1));
        ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wNs, Ns));
        ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wNtr, Ntr));

		int Np1 = (int)(p->np1);
		int Np2 = (int)(p->np2);
		long NpKick = Np1*Np2;

		pTempM1 = new double[NpKick];
		pTempM2 = new double[NpKick];
		pTempM2 = new double[NpKick];
		pTempIBe2 = new double[NpKick];
		pTempArg1 = new double[Np1];
		pTempArg2 = new double[Np2];

        char StrCom[2048];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hCom, 2045, StrCom));

		int LenStr = 0;
        ProcErrRad(RadFldFocKickPer(pTempM1, pTempM2, pTempIBe2, pTempArg1, pTempArg2, &LenStr, (int)(p->obj), P1, Ns, (double)(p->per), (int)(p->nper), (int)(p->nps), Ntr, (double)(p->r1), Np1, (double)(p->d1), (double)(p->r2), Np2, (double)(p->d2), (int)(p->nh), StrCom));

		int NumCols = Np2;
		if(NumCols <= 1) NumCols = 0;
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wM1, Np1, NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wM1, pTempM1, NpKick, 0));
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wM2, Np1, NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wM2, pTempM2, NpKick, 0));
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wArg1, Np1, 0));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wArg1, pTempArg1, Np1, 0));
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wArg2, Np2, 0));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wArg2, pTempArg2, Np2, 0));

		pTempStr = new char[LenStr + 5];
		ProcErrRad(RadFldFocKickPerFormStr(pTempStr, pTempM1, pTempM2, pTempIBe2, pTempArg1, pTempArg2, Np1, Np2, (double)(p->per), (int)(p->nper), StrCom));

        Handle hStr = NewHandle(LenStr);
        strncpy(*hStr, pTempStr, LenStr);
        p->hFormStr = hStr;

		if(pTempM1 != 0) { delete[] pTempM1; pTempM1 = 0;}
		if(pTempM2 != 0) { delete[] pTempM2; pTempM2 = 0;}
		if(pTempIBe2 != 0) { delete[] pTempIBe2; pTempIBe2 = 0;}
		if(pTempArg1 != 0) { delete[] pTempArg1; pTempArg1 = 0;}
		if(pTempArg2 != 0) { delete[] pTempArg2; pTempArg2 = 0;}
		if(pTempStr != 0) { delete[] pTempStr; pTempStr = 0;}
	}
	catch(int ErrNo) 
	{
		if(pTempM1 != 0) { delete[] pTempM1; pTempM1 = 0;}
		if(pTempM2 != 0) { delete[] pTempM2; pTempM2 = 0;}
		if(pTempIBe2 != 0) { delete[] pTempIBe2; pTempIBe2 = 0;}
		if(pTempArg1 != 0) { delete[] pTempArg1; pTempArg1 = 0;}
		if(pTempArg2 != 0) { delete[] pTempArg2; pTempArg2 = 0;}
		if(pTempStr != 0) { delete[] pTempStr; pTempStr = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldEnr {
    waveHndl wSbdPar;
	DOUBLE objsrc;
	DOUBLE objdst;
	DOUBLE d;
};
static int radFldEnr(void* p_void)
{
	try
	{
        radTIgorRadFldEnr *p = (radTIgorRadFldEnr*)p_void;

		double dSbdPar[] = {1,1,1};
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wSbdPar, dSbdPar));
		//int SbdPar[] = {dSbdPar[0], dSbdPar[1], dSbdPar[2]};
		int SbdPar[] = {(int)dSbdPar[0], (int)dSbdPar[1], (int)dSbdPar[2]};

		double D;
		p->d = 0;
        ProcErrRad(RadFldEnr(&D, (int)(p->objdst), (int)(p->objsrc), SbdPar));
		p->d = D;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldEnrFrc {
    waveHndl wSbdPar;
    Handle hID;
	DOUBLE objsrc;
	DOUBLE objdst;
    waveHndl wF;
	DOUBLE nF; // first component of the force calculated
};
static int radFldEnrFrc(void* p_void)
{
	radTIgorRadFldEnrFrc *p = (radTIgorRadFldEnrFrc*)p_void;
    double *pTempCalcData = 0;

	//int hState;
	//double *pdF=0;
	//float *pfF =0;
	//bool DoubleArrBIsLocal = false;

	try
	{
		double dSbdPar[] = {1,1,1};
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wSbdPar, dSbdPar));
		//int SbdPar[] = {dSbdPar[0], dSbdPar[1], dSbdPar[2]};
		int SbdPar[] = {(int)dSbdPar[0], (int)dSbdPar[1], (int)dSbdPar[2]};

        char StrID[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 30, StrID));

		p->nF = 0;
		//long LenF=0;
		//ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wF, pdF, pfF, LenF, hState));
		//if((pdF == 0) && (pfF != 0) && (LenF > 0))
		//{
        //	pdF = new double[LenF];
		//	DoubleArrBIsLocal = true;
		//}

		int Np = 1;
        const int MaxAmOfFieldValPerPoint = 10;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        int NF=0;
        //ProcErrRad(RadFldEnrFrc(pdF, &NF, (int)(p->objdst), (int)(p->objsrc), StrID, SbdPar));
		//p->nF = *pdF;

        ProcErrRad(RadFldEnrFrc(pTempCalcData, &NF, (int)(p->objdst), (int)(p->objsrc), StrID, SbdPar));
        p->nF = *pTempCalcData;

		//int NumCols = (int)(((double)NF)/((double)(p->np)) + 1E-06);
		//if(NumCols <= 1) NumCols = 0;

        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wF, NF, 0));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wF, pTempCalcData, NF, "N"));
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}

		//if(DoubleArrBIsLocal && (pdF != 0)) 
		//{
		//	float *tfF = pfF;
		//	double *tdF = pdF;
		//	for(long i=0; i<NF; i++) *(tfF++) = (float)(*(tdF++));
		//	delete[] pdF;
		//	pdF = 0;
		//}
		//ProcErr(srTIgorInterf::ReleaseWave(p->wF, hState));
	}
	catch(int ErrNo) 
	{
		//srTIgorInterf::ReleaseWave(p->wF, hState);
		//if(DoubleArrBIsLocal && (pdF != 0)) delete[] pdF;

		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldEnrTrq {
    waveHndl wSbdPar;
    waveHndl wP;
    Handle hID;
	DOUBLE objsrc;
	DOUBLE objdst;
    waveHndl wF;
	DOUBLE nF; // first component of the torque calculated
};
static int radFldEnrTrq(void* p_void)
{
	radTIgorRadFldEnrTrq *p = (radTIgorRadFldEnrTrq*)p_void;
    double *pTempCalcData = 0;

	//int hState;
	//double *pdF=0;
	//float *pfF =0;
	//bool DoubleArrBIsLocal = false;

	try
	{
		double dSbdPar[] = {1,1,1}, P[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wSbdPar, dSbdPar));
		//int SbdPar[] = {dSbdPar[0], dSbdPar[1], dSbdPar[2]};
		int SbdPar[] = {(int)dSbdPar[0], (int)dSbdPar[1], (int)dSbdPar[2]};

		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));

        char StrID[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 30, StrID));

		p->nF = 0;
		//long LenF=0;
		//ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wF, pdF, pfF, LenF, hState));
		//if((pdF == 0) && (pfF != 0) && (LenF > 0))
		//{
        //	pdF = new double[LenF];
		//	DoubleArrBIsLocal = true;
		//}

		int Np = 1;
        const int MaxAmOfFieldValPerPoint = 10;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        int NF=0;
        //ProcErrRad(RadFldEnrTrq(pdF, &NF, (int)(p->objdst), (int)(p->objsrc), StrID, P, SbdPar));
		//p->nF = NF;

        ProcErrRad(RadFldEnrTrq(pTempCalcData, &NF, (int)(p->objdst), (int)(p->objsrc), StrID, P, SbdPar));
        p->nF = *pTempCalcData;

        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wF, NF, 0));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wF, pTempCalcData, NF, "N*mm"));
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}

		//if(DoubleArrBIsLocal && (pdF != 0)) 
		//{
		//	float *tfF = pfF;
		//	double *tdF = pdF;
		//	for(long i=0; i<NF; i++) *(tfF++) = (float)(*(tdF++));
		//	delete[] pdF;
		//	pdF = 0;
		//}
        //ProcErr(srTIgorInterf::ReleaseWave(p->wF, hState));
	}
	catch(int ErrNo) 
	{
		//srTIgorInterf::ReleaseWave(p->wF, hState);
		//if(DoubleArrBIsLocal && (pdF != 0)) delete[] pdF;

		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldFrc {
	DOUBLE shape;
	DOUBLE obj;
    waveHndl wF;
	DOUBLE nF; // first component of the force calculated
};
static int radFldFrc(void* p_void)
{
	radTIgorRadFldFrc *p = (radTIgorRadFldFrc*)p_void;
    double *pTempCalcData = 0;

	//int hState;
	//double *pdF=0;
	//float *pfF =0;
	//bool DoubleArrBIsLocal = false;

	try
	{
		p->nF = 0;
		//long LenF=0;
		//ProcErr(srTIgorInterf::GetDataPtrFromWaveDoubleOrFloat1D(p->wF, pdF, pfF, LenF, hState));
		//if((pdF == 0) && (pfF != 0) && (LenF > 0))
		//{
        //	pdF = new double[LenF];
		//	DoubleArrBIsLocal = true;
		//}

        int NF=0;
        //ProcErrRad(RadFldFrc(pdF, &NF, (int)(p->obj), (int)(p->shape)));
		//p->nF = NF;

        ProcErrRad(RadFldFrc(pTempCalcData, &NF, (int)(p->obj), (int)(p->shape)));
		p->nF = *pTempCalcData;

        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wF, NF, 0));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wF, pTempCalcData, NF, "N"));
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}

		//if(DoubleArrBIsLocal && (pdF != 0)) 
		//{
		//	float *tfF = pfF;
		//	double *tdF = pdF;
		//	for(long i=0; i<NF; i++) *(tfF++) = (float)(*(tdF++));
		//	delete[] pdF;
		//	pdF = 0;
		//}
		//ProcErr(srTIgorInterf::ReleaseWave(p->wF, hState));
	}
	catch(int ErrNo) 
	{
		//srTIgorInterf::ReleaseWave(p->wF, hState);
		//if(DoubleArrBIsLocal && (pdF != 0)) delete[] pdF;

		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldShimSig {
    waveHndl wVi;
	DOUBLE np;
    waveHndl wP2;
    waveHndl wP1;
    waveHndl wV;
    Handle hID;
	DOUBLE obj;
    waveHndl wB;
	DOUBLE B0;
};
static int radFldShimSig(void* p_void)
{
	radTIgorRadFldShimSig *p = (radTIgorRadFldShimSig*)p_void;
    double *pTempCalcData = 0;

	try
	{
		double V[3], P1[3], P2[3], Vi[3];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wV, V));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP1, P1));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP2, P2));
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wVi, Vi));

		p->B0 = 0;
		int Np = (int)(p->np);

        const int MaxAmOfFieldValPerPoint = 30;
		pTempCalcData = new double[MaxAmOfFieldValPerPoint*Np];

        char StrID[31];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hID, 30, StrID));

        int Nb=0;
        ProcErrRad(RadFldShimSig(pTempCalcData, &Nb, (int)(p->obj), StrID, V, P1, P2, Np, Vi));
		p->B0 = *pTempCalcData;

		double dx = P2[0] - P1[0], dy = P2[1] - P1[1], dz = P2[2] - P1[2];
		double ArgRange = sqrt(dx*dx + dy*dy + dz*dz);
		double ArgStep = (Np <= 1)? 0 : ArgRange/(Np - 1);

		int NumRows=0, NumCols=0, IndDim = 0;
		if(Nb <= Np)
		{
			NumRows = Np;
			IndDim = 0;
		}
		else if(Np == 1)
		{
			NumRows = Nb;
			IndDim = 0;
		}
		else
		{
            NumRows = (int)(((double)Nb)/((double)Np) + 1E-06);
			NumCols = Np;
			IndDim = 1;
		}

		char strFldUnit[] = "T*mm\0";
		bool isFieldInt = ((*StrID == 'i') || (*StrID == 'I'));
		if(!isFieldInt) { strFldUnit[0] = 'T'; strFldUnit[1] = '\0';}

		double ArgStarts[] = {0, 0}, ArgSteps[] = {1, 1};
		char* ArgValUnits[] = {"", "", strFldUnit};
        ArgStarts[IndDim] = 0.; //(p->start)*0.001;
		ArgSteps[IndDim] = ArgStep*0.001;
		ArgValUnits[IndDim] = "m";
        ProcErr(srTIgorInterf::ReDimNumWave2D(p->wB, NumRows, NumCols));
		ProcErr(srTIgorInterf::SetDataInNumWave(p->wB, pTempCalcData, Nb, strFldUnit));
        ProcErr(srTIgorInterf::SetScaleNumWave(p->wB, 2, ArgStarts, ArgSteps, ArgValUnits));

		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
	}
	catch(int ErrNo) 
	{
		if(pTempCalcData != 0) { delete[] pTempCalcData; pTempCalcData = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldCmpCrt {
	DOUBLE prcTrjAng;
	DOUBLE prcTrjCrd;
	DOUBLE prcFrc;
	DOUBLE prcBInt;
	DOUBLE prcA;
	DOUBLE prcB;
	DOUBLE n;
};
static int radFldCmpCrt(void* p_void)
{
	try
	{
        radTIgorRadFldCmpCrt *p = (radTIgorRadFldCmpCrt*)p_void;

		int N;
		p->n = 0;
        ProcErrRad(RadFldCmpCrt(&N, (double)(p->prcB), (double)(p->prcA), (double)(p->prcBInt), (double)(p->prcFrc), (double)(p->prcTrjCrd), (double)(p->prcTrjAng)));
		p->n = N;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldCmpPrc {
	//Handle hOpt8;
	//Handle hOpt7;
	//Handle hOpt6;
	//Handle hOpt5;
	//Handle hOpt4;
	//Handle hOpt3;
	//Handle hOpt2;
	Handle hOpt1;
	DOUBLE n;
};
static int radFldCmpPrc(void* p_void)
{
	try
	{
        radTIgorRadFldCmpPrc *p = (radTIgorRadFldCmpPrc*)p_void;

        char StrOpt1[2512]; //, StrOpt2[512], StrOpt3[512], StrOpt4[512], StrOpt5[512], StrOpt6[512], StrOpt7[512], StrOpt8[512];
        *StrOpt1 = '\0'; //*StrOpt2 = '\0'; *StrOpt3 = '\0'; *StrOpt4 = '\0'; *StrOpt5 = '\0'; *StrOpt6 = '\0'; *StrOpt7 = '\0'; *StrOpt8 = '\0';

		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt1, 2511, StrOpt1));

		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt2, 512, StrOpt2));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt3, 512, StrOpt3));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt4, 512, StrOpt4));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt5, 512, StrOpt5));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt6, 512, StrOpt6));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt7, 512, StrOpt7));
		//ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOpt8, 512, StrOpt8));

		int N;
		p->n = 0;
        //ProcErrRad(RadFldCmpPrc(&N, StrOpt1, StrOpt2, StrOpt3, StrOpt4, StrOpt5, StrOpt6, StrOpt7, StrOpt8));
        ProcErrRad(RadFldCmpPrc(&N, StrOpt1));
		p->n = N;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldUnits {
	Handle result;
};
static int radFldUnits(void* p_void)
{
	try
	{
        radTIgorRadFldUnits *p = (radTIgorRadFldUnits*)p_void;
        char Str[2048];

        ProcErrRad(RadFldUnits(Str));
        long LenStr = (long)strlen(Str);
        Handle hStr = NewHandle(LenStr);
        strncpy(*hStr, Str, LenStr);
        p->result = hStr;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldUnitsSize {
	DOUBLE result;
};
static int radFldUnitsSize(void* p_void)
{
	try
	{
        radTIgorRadFldUnitsSize *p = (radTIgorRadFldUnitsSize*)p_void;

		p->result = 0;
		int Size=0;
		ProcErrRad(RadFldUnitsSize(&Size));
        p->result = Size;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldLenRndSw {
	Handle hOnOrOff;
	DOUBLE n;
};
static int radFldLenRndSw(void* p_void)
{
	try
	{
        radTIgorRadFldLenRndSw *p = (radTIgorRadFldLenRndSw*)p_void;

        char StrOnOrOff[512];
        *StrOnOrOff = '\0'; 
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hOnOrOff, 512, StrOnOrOff));

		p->n = 0;
		int N = 0;
		ProcErrRad(RadFldLenRndSw(&N, StrOnOrOff));
        p->n = N;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldLenTol {
	DOUBLE ZeroVal;
	DOUBLE RelVal;
	DOUBLE AbsVal;
	DOUBLE n;
};
static int radFldLenTol(void* p_void)
{
	try
	{
        radTIgorRadFldLenTol *p = (radTIgorRadFldLenTol*)p_void;

		p->n = 0;
		int N = 0;
		ProcErrRad(RadFldLenTol(&N, (double)(p->AbsVal), (double)(p->RelVal), (double)(p->ZeroVal)));
        p->n = N;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadFldFrcShpRtg {
	waveHndl wW;
	waveHndl wP;
	DOUBLE n;
};
static int radFldFrcShpRtg(void* p_void)
{
	try
	{
        radTIgorRadFldFrcShpRtg *p = (radTIgorRadFldFrcShpRtg*)p_void;

		double P[3], W[2];
		ProcErr(srTIgorInterf::Get3ElemArrDoubleFromNumWave1D(p->wP, P));
		ProcErr(srTIgorInterf::Get2ElemArrDoubleFromNumWave1D(p->wW, W));

		p->n = 0;
		int N = 0;
		ProcErrRad(RadFldFrcShpRtg(&N, P, W));
        p->n = N;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadUtiDel {
	DOUBLE obj;
	DOUBLE n;
};
static int radUtiDel(void* p_void)
{
	try
	{
        radTIgorRadUtiDel *p = (radTIgorRadUtiDel*)p_void;

		p->n = 0;
		int N = 0;
		ProcErrRad(RadUtiDel(&N, (int)(p->obj)));
        p->n = N;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadUtiDelAll {
	DOUBLE n;
};
static int radUtiDelAll(void* p_void)
{
	try
	{
        radTIgorRadUtiDelAll *p = (radTIgorRadUtiDelAll*)p_void;

		p->n = 0;
		int N = 0;
		ProcErrRad(RadUtiDelAll(&N));
        p->n = N;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadUtiDmp {
	Handle hAscOrBin;
	//DOUBLE obj;
	waveHndl wElems;
	Handle result;
};
static int radUtiDmp(void* p_void)
{
	int *arElems=0;
	try
	{
        radTIgorRadUtiDmp *p = (radTIgorRadUtiDmp*)p_void;

		char strAscOrBin[4];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hAscOrBin, 3, strAscOrBin));

        long nElems=0;
		ProcErr(srTIgorInterf::GetArrIntFromNumWave1D(p->wElems, -1, arElems, nElems));

		int Size=0;
		ProcErrRad(RadUtiDmpSize(&Size, arElems, nElems, strAscOrBin, false));

        //char Str[4096];
        ////ProcErrRad(RadUtiDmp(Str, (int)(p->obj)));
        //ProcErrRad(RadUtiDmp(Str, (int)(p->obj), strAscOrBin));

        //long LenStr = (long)strlen(Str);
        //Handle hStr = NewHandle(LenStr);
        //strncpy(*hStr, Str, LenStr);

        Handle hStr = NewHandle(Size);
		char *thStr = *hStr;
		ProcErrRad(RadUtiDmpRead(thStr, strAscOrBin));

        p->result = hStr;

		if(arElems != 0) { delete[] arElems; arElems = 0;}
	}
	catch(int ErrNo) 
	{
		if(arElems != 0) { delete[] arElems; arElems = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadUtiDmpSize {
	Handle hAscOrBin;
	//DOUBLE obj;
	waveHndl wElems;
	DOUBLE result;
};
static int radUtiDmpSize(void* p_void)
{//not necessary, since radUtiDmp gets the correct size
	int *arElems=0;
	try
	{
        radTIgorRadUtiDmpSize *p = (radTIgorRadUtiDmpSize*)p_void;

		char strAscOrBin[4];
		ProcErr(srTIgorInterf::ExtractCStrFromHandle(p->hAscOrBin, 3, strAscOrBin));

        long nElems=0;
		ProcErr(srTIgorInterf::GetArrIntFromNumWave1D(p->wElems, -1, arElems, nElems));

		p->result = 0;
		int Size=0;
		ProcErrRad(RadUtiDmpSize(&Size, arElems, nElems, strAscOrBin, false));

		//int Size=0;
		//ProcErrRad(RadUtiDmpSize(&Size, (int)(p->obj)));
        p->result = Size;

		if(arElems != 0) { delete[] arElems; arElems = 0;}
	}
	catch(int ErrNo) 
	{
		if(arElems != 0) { delete[] arElems; arElems = 0;}
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadUtiIntrptTim {
	DOUBLE t;
	DOUBLE d;
};
static int radUtiIntrptTim(void* p_void)
{
	try
	{
        radTIgorRadUtiIntrptTim *p = (radTIgorRadUtiIntrptTim*)p_void;

		p->d = 0;
		double D = 0;
		ProcErrRad(RadUtiIntrptTim(&D, (double)(p->t)));
        p->d = D;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct radTIgorRadUtiVer {
	DOUBLE d;
};
static int radUtiVer(void* p_void)
{
	try
	{
        radTIgorRadUtiVer *p = (radTIgorRadUtiVer*)p_void;

		p->d = 0;
		double D = 0;
		ProcErrRad(RadUtiVer(&D));
        p->d = D;
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

void radEnsureInit()
{
	char TestName[] = "RadVerProc";
	char StrDummy[50];
	int RadiaWasNotInit = FetchStrVar(TestName, StrDummy);
	if(RadiaWasNotInit)
	{
		char InitStr[] = "RadInit(1)";
		XOPCommand(InitStr);
	}
	SetXOPType(RESIDENT);
}

//*************************************************************************
/*	RegisterFunction()
	Igor calls this at startup time to find the address of the
	XFUNCs added by this XOP. See XOP manual regarding "Direct XFUNCs".
*/
static long RegisterFunction()
{
	int funcIndex;
	funcIndex = GetXOPItem(0);		// Which function is Igor asking about?
	switch (funcIndex) {
		case 0:
			return((long)radObjRecMag);
			break;
		case 1:
			return((long)radObjThckPgn);
			break;
		case 2:
			return((long)radObjPolyhdr);
			break;
		case 3:
			return((long)radObjMltExtPgn);
			break;
		case 4:
			return((long)radObjMltExtRtg);
			break;
		case 5:
			return((long)radObjMltExtTri);
			break;
		case 6:
			return((long)radObjArcPgnMag);
			break;
		case 7:
			return((long)radObjCylMag);
			break;
		case 8:
			return((long)radObjRecCur);
			break;
		case 9:
			return((long)radObjArcCur);
			break;
		case 10:
			return((long)radObjRaceTrk);
			break;
		case 11:
			return((long)radObjFlmCur);
			break;
		case 12:
			return((long)radObjBckg);
			break;
		case 13:
			return((long)radObjCnt);
			break;
		case 14:
			return((long)radObjAddToCnt);
			break;
		case 15:
			return((long)radObjCntSize);
			break;
		case 16:
			return((long)radObjCntStuf);
			break;
		case 17:
			return((long)radObjM);
			break;
		case 18:
			return((long)radObjDpl);
			break;
		case 19:
			return((long)radObjCutMag);
			break;
		case 20:
			return((long)radObjDivMagPln);
			break;
		case 21:
			return((long)radObjDivMagCyl);
			break;
		case 22:
			return((long)radObjGeoVol);
			break;
		case 23:
			return((long)radObjGeoLim);
			break;
		case 24:
			return((long)radObjDrwOpenGL);
			break;
		case 25:
			return((long)radObjDrwAtr);
			break;
		case 26:
			return((long)radTrfPlSym);
			break;
		case 27:
			return((long)radTrfRot);
			break;
		case 28:
			return((long)radTrfTrsl);
			break;
		case 29:
			return((long)radTrfInv);
			break;
		case 30:
			return((long)radTrfCmbL);
			break;
		case 31:
			return((long)radTrfCmbR);
			break;
		case 32:
			return((long)radTrfMlt);
			break;
		case 33:
			return((long)radTrfOrnt);
			break;
		case 34:
			return((long)radTrfZerPara);
			break;
		case 35:
			return((long)radTrfZerPerp);
			break;
		case 36:
			return((long)radMatLin);
			break;
		case 37:
			return((long)radMatSatIso);
			break;
		case 38:
			return((long)radMatSatLam);
			break;
		case 39:
			return((long)radMatSatAniso);
			break;
		case 40:
			return((long)radMatStd);
			break;
		case 41:
			return((long)radMatApl);
			break;
		case 42:
			return((long)radMatMvsH);
			break;
		case 43:
			return((long)radRlxPre);
			break;
		case 44:
			return((long)radRlxMan);
			break;
		case 45:
			return((long)radRlxAuto);
			break;
		case 46:
			return((long)radRlxAutoExt);
			break;
		case 47:
			return((long)radSolve);
			break;
		case 48:
			return((long)radFld);
			break;
		case 49:
			return((long)radFldLst);
			break;
		case 50:
			return((long)radFldInt);
			break;
		case 51:
			return((long)radFldPtcTrj);
			break;
		case 52:
			return((long)radFldFocPot);
			break;
		case 53:
			return((long)radFldFocKickPer);
			break;
		case 54:
			return((long)radFldEnr);
			break;
		case 55:
			return((long)radFldEnrFrc);
			break;
		case 56:
			return((long)radFldEnrTrq);
			break;
		case 57:
			return((long)radFldFrc);
			break;
		case 58:
			return((long)radFldCmpCrt);
			break;
		case 59:
			return((long)radFldCmpPrc);
			break;
		case 60:
			return((long)radFldUnits);
			break;
		//case 55:
		//	return((long)radFldUnitsSize);
		//	break;
		case 61:
			return((long)radFldLenRndSw);
			break;
		case 62:
			return((long)radFldLenTol);
			break;
		case 63:
			return((long)radFldFrcShpRtg);
			break;
		case 64:
			return((long)radUtiDel);
			break;
		case 65:
			return((long)radUtiDelAll);
			break;
		case 66:
			return((long)radUtiDmp);
			break;
		//case 62:
		//	return((long)radUtiDmpSize);
		//	break;
		case 67:
			return((long)radUtiIntrptTim);
			break;
		case 68:
			return((long)radUtiVer);
			break;
		case 69:
			return((long)radObjCnt0);
			break;
		case 70:
			return((long)radObjFullMag);
			break;			
		case 71:
			return((long)radObjDrwQD3D);
			break;			
		case 72:
			return((long)radFldShimSig);
			break;			
		case 73:
			return((long)radFldC);
			break;			
		case 74:
			return((long)radFldLstC);
			break;			
		case 75:
			return((long)radFldIntC);
			break;			
		case 76:
			return((long)radFldPtcTrjC);
			break;			
		case 77:
			return((long)radFldFocPotC);
			break;			
		case 78:
			return((long)radObjRecMagC);
			break;			
		case 79:
			return((long)radObjCylMagC);
			break;			
		case 80:
			return((long)radObjRecCurC);
			break;			
		case 81:
			return((long)radObjArcCurC);
			break;			
		case 82:
			return((long)radObjRaceTrkC);
			break;
	}
	return(NIL);
}

//*************************************************************************
/*	DoFunction()
	
	Igor calls this when the user invokes one if the XOP's XFUNCs
	if we returned NIL for the XFUNC from RegisterFunction. In this
	XOP, we always use the direct XFUNC method, so Igor will never call
	this function. See XOP manual regarding "Direct XFUNCs".
*/
static int DoFunction()
{
	int funcIndex;
	void *p;				// Pointer to structure containing function parameters and result.
	int err = 0;

	funcIndex = GetXOPItem(0);		// Which function is being invoked ?
	p = (void *)GetXOPItem(1);		// Get pointer to params/result.
	switch (funcIndex) {
		case 0:
			err = radObjRecMag(p);
			break;
		case 1:
			err = radObjThckPgn(p);
			break;
		case 2:
			err = radObjPolyhdr(p);
			break;
		case 3:
			err = radObjMltExtPgn(p);
			break;
		case 4:
			err = radObjMltExtRtg(p);
			break;
		case 5:
			err = radObjMltExtTri(p);
			break;
		case 6:
			err = radObjArcPgnMag(p);
			break;
		case 7:
			err = radObjCylMag(p);
			break;
		case 8:
			err = radObjRecCur(p);
			break;
		case 9:
			err = radObjArcCur(p);
			break;
		case 10:
			err = radObjRaceTrk(p);
			break;
		case 11:
			err = radObjFlmCur(p);
			break;
		case 12:
			err = radObjBckg(p);
			break;
		case 13:
			err = radObjCnt(p);
			break;
		case 14:
			err = radObjAddToCnt(p);
			break;
		case 15:
			err = radObjCntSize(p);
			break;
		case 16:
			err = radObjCntStuf(p);
			break;
		case 17:
			err = radObjM(p);
			break;
		case 18:
			err = radObjDpl(p);
			break;
		case 19:
			err = radObjCutMag(p);
			break;
		case 20:
			err = radObjDivMagPln(p);
			break;
		case 21:
			err = radObjDivMagCyl(p);
			break;
		case 22:
			err = radObjGeoVol(p);
			break;
		case 23:
			err = radObjGeoLim(p);
			break;
		case 24:
			err = radObjDrwOpenGL(p);
			break;
		case 25:
			err = radObjDrwAtr(p);
			break;
		case 26:
			err = radTrfPlSym(p);
			break;
		case 27:
			err = radTrfRot(p);
			break;
		case 28:
			err = radTrfTrsl(p);
			break;
		case 29:
			err = radTrfInv(p);
			break;
		case 30:
			err = radTrfCmbL(p);
			break;
		case 31:
			err = radTrfCmbR(p);
			break;
		case 32:
			err = radTrfMlt(p);
			break;
		case 33:
			err = radTrfOrnt(p);
			break;
		case 34:
			err = radTrfZerPara(p);
			break;
		case 35:
			err = radTrfZerPerp(p);
			break;
		case 36:
			err = radMatLin(p);
			break;
		case 37:
			err = radMatSatIso(p);
			break;
		case 38:
			err = radMatSatLam(p);
			break;
		case 39:
			err = radMatSatAniso(p);
			break;
		case 40:
			err = radMatStd(p);
			break;
		case 41:
			err = radMatApl(p);
			break;
		case 42:
			err = radMatMvsH(p);
			break;
		case 43:
			err = radRlxPre(p);
			break;
		case 44:
			err = radRlxMan(p);
			break;
		case 45:
			err = radRlxAuto(p);
			break;
		case 46:
			err = radRlxAutoExt(p);
			break;
		case 47:
			err = radSolve(p);
			break;
		case 48:
			err = radFld(p);
			break;
		//case 44:
		//	err = radFld1(p);
		//	break;
		case 49:
			err = radFldLst(p);
			break;
		case 50:
			err = radFldInt(p);
			break;
		case 51:
			err = radFldPtcTrj(p);
			break;
		case 52:
			err = radFldFocPot(p);
			break;
		case 53:
			err = radFldFocKickPer(p);
			break;
		case 54:
			err = radFldEnr(p);
			break;
		case 55:
			err = radFldEnrFrc(p);
			break;
		case 56:
			err = radFldEnrTrq(p);
			break;
		case 57:
			err = radFldFrc(p);
			break;
		case 58:
			err = radFldCmpCrt(p);
			break;
		case 59:
			err = radFldCmpPrc(p);
			break;
		case 60:
			err = radFldUnits(p);
			break;
		//case 55:
		//	err = radFldUnitsSize(p);
		//	break;
		case 61:
			err = radFldLenRndSw(p);
			break;
		case 62:
			err = radFldLenTol(p);
			break;
		case 63:
			err = radFldFrcShpRtg(p);
			break;
		case 64:
			err = radUtiDel(p);
			break;
		case 65:
			err = radUtiDelAll(p);
			break;
		case 66:
			err = radUtiDmp(p);
			break;
		//case 62:
		//	err = radUtiDmpSize(p);
		//	break;
		case 67:
			err = radUtiIntrptTim(p);
			break;
		case 68:
			err = radUtiVer(p);
			break;
		case 69:
			err = radObjCnt0(p);
			break;
		case 70:
			err = radObjFullMag(p);
			break;
		case 71:
			err = radObjDrwQD3D(p);
			break;
		case 72:
			err = radFldShimSig(p);
			break;
		case 73:
			err = radFldC(p);
			break;
		case 74:
			err = radFldLstC(p);
			break;
		case 75:
			err = radFldIntC(p);
			break;
		case 76:
			err = radFldPtcTrjC(p);
			break;
		case 77:
			err = radFldFocPotC(p);
			break;
		case 78:
			err = radObjRecMagC(p);
			break;
		case 79:
			err = radObjCylMagC(p);
			break;
		case 80:
			err = radObjRecCurC(p);
			break;
		case 81:
			err = radObjArcCurC(p);
			break;
		case 82:
			err = radObjRaceTrkC(p);
			break;
	}
	return(err);
}

//*************************************************************************
/*	XOPEntry()

	This is the entry point from the host application to the XOP for all messages after the
	INIT message.
*/
static void
XOPEntry(void)
{	
	long result = 0;

	switch (GetXOPMessage()) {
		case FUNCTION:						// Our external function being invoked ?
			result = DoFunction();
			break;

		case FUNCADDRS:
			result = RegisterFunction();
			break;

		case IDLE:
			radEnsureInit();
			break;
	}
	SetXOPResult(result);
}

//*************************************************************************
/*	main(ioRecHandle)

	This is the initial entry point at which the host application calls XOP.
	The message sent by the host must be INIT.
	main() does any necessary initialization and then sets the XOPEntry field of the
	ioRecHandle to the address to be called for future messages.
*/
HOST_IMPORT void main(IORecHandle ioRecHandle)
{
//#ifdef applec					/* for MPW C for 68K only */
//	void _DATAINIT(void);
//	_DATAINIT();				/* for MPW C only */
//	UnloadSeg(_DATAINIT);
//#endif
//
//#ifdef XOP_GLOBALS_ARE_A4_BASED
//	#ifdef __MWERKS__
//		SetCurrentA4();									/* Set up correct A4. This allows globals to work. */
//		SendXOPA4ToIgor(ioRecHandle, GetA4());			/* And communicate it to Igor. */
//	#endif
//#endif
//	
//	LoadXOPSegs();
//	XOPInit(ioRecHandle);							/* do standard XOP initialization */
//	SetXOPEntry(XOPEntry);							/* set entry point for future calls */
//	
//	if (igorVersion < 200)
//	{} //SetXOPResult(REQUIRES_IGOR_200);
//	else
//		SetXOPResult(0L);

	//#ifdef applec					// For MPW C for 68K only.
	//	void _DATAINIT(void);
	//	_DATAINIT();				// For MPW C only.
	//	UnloadSeg(_DATAINIT);
	//#endif
	
	#ifdef XOP_GLOBALS_ARE_A4_BASED
		#ifdef __MWERKS__
			SetCurrentA4();							// Set up correct A4. This allows globals to work.
			SendXOPA4ToIgor(ioRecHandle, GetA4());	// And communicate it to Igor.
		#endif
	#endif

	//LoadXOPSegs(); //OC210705
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
	
	//if(igorVersion < 200) SetXOPResult(OLD_IGOR);
	//else SetXOPResult(0L);
    SetXOPResult(0L);
	SetXOPType(RESIDENT | IDLES);

	RadUtiYeldFuncSet(&SpinProcess);
	double d = 0;
	RadUtiIntrptTim(&d, 0.5);


	//srSetCompProgressIndicFunc(&(srTIgorSend::ShowCompProgress));
	//srSetWfrExtModifFunc(&(srTIgorSend::WfrModify));
}

//*************************************************************************

// All structures are 2-byte-aligned.
#if GENERATINGPOWERPC
	#pragma options align=reset
#endif
#ifdef _WINDOWS_
	#pragma pack()
#endif
