
#include "radstlon.h"
#include "auxparse.h"

//#ifdef WIN32
//#define WIN32_LEAN_AND_MEAN
////#include <windows.h>
//#endif

#include "radentry.h"
#include "radiobuf.h"

//-------------------------------------------------------------------------

extern "C" {

void RecMag( double,double,double, double,double,double, double,double,double );
void ExtrudedPolygonDLL( double, double, double*, int, char, double* );
//void PolyhedronDLL( double*, int, int*, int*, int, double* );
void PolyhedronDLL( double*, int, int*, int*, int, double*, double*, double*, double* );
void MultGenExtrPolygonDLL( double*, int*, double*, int, double* );
void MultGenExtrRectangleDLL( double*, double*, int, double* );
//void MultGenExtrTriangleDLL( double, double, double*, double*, int, char, double*, const char*,const char*,const char* );
void MultGenExtrTriangleDLL( double, double, double*, double*, int, char, double*, const char*,const char*,const char*,const char* ); //OC30072018

void ArcMag( double,double,double, double,double, double,double, double, int, char*, double,double,double );
void ArcPolygon();
void ArcPolygonDLL( double,double, char, double*, int, double,double, int, char, double,double,double );
void CylMag( double,double,double, double, double, int, char*, double,double,double );
void RecCur( double,double,double, double,double,double, double,double,double );
void ArcCur( double,double,double, double,double, double,double, double, int, double, char*, char* );
void RaceTrack( double,double,double, double,double, double,double, double, int, double, char*, char* );
void FlmCurDLL( double*, int, double );
void ScaleCurInObj( int,double );
void BackgroundFieldSource( double,double,double );
void Rectngl( double,double,double, double,double );
void Group( int*, long );
void AddToGroup( int, int*, long );
void OutGroupSize( int );
void OutGroupSubObjectKeys( int );

void DuplicateElementG3DOpt( int, const char* );
void CutElementG3DOpt( int, double,double,double, double,double,double, const char* );
void SubdivideElementG3DOpt( int, double*, char, double*, int, const char*, const char*, const char* );
void GeometricalVolume( int );
void GeometricalLimits( int );
void NumberOfDegOfFreedom( int );

void MagnOfObj( int );
void ObjField( int, char* );
void SetObjMagn( int, double,double,double );

void Translation( double,double,double );
void Rotation( double,double,double, double,double,double, double );
void PlaneSym( double,double,double, double,double,double );
void FieldInversion();
void CombineTransformLeft( int, int );
void CombineTransformRight( int, int );
void TransformObject( int, int );
void ApplySymmetry( int, int, int );

void LinearMaterial( double,double, double,double,double );
void LinearMaterial2( double,double, double );
void MaterialStd( char*, double );

void NonlinearIsotropMaterial2( double,double, double,double, double,double );
void NonlinearIsotropMaterial3Opt( double**, long );
void NonlinearLaminatedMaterialFrm( double*,double*,double*, double, double* );
void NonlinearLaminatedMaterialTab( double*, int, double, double* );
void NonlinearAnisotropMaterialOpt0( double*, int, double*, int );
void NonlinearAnisotropMaterialOpt1( double**, double** );
void NonlinearAnisotropMaterialOpt2( double**, double );
void NonlinearAnisotropMaterialOpt3( double, double** );
void ApplyMaterial( int, int );
void MvsH( int, char*, double,double,double );

void PreRelax( int, int );
void ShowInteractMatrix(int);
void ShowInteractVector(int, char*);
void ManualRelax( int, int, int, double );
//void AutoRelax( int, double, int, int );
void AutoRelaxOpt( int, double, int, int, const char* );
void UpdateSourcesForRelax( int );
void SolveGen( int, double, int, int );

void FieldArbitraryPointsArray( long, const char*, double**, long );
void Field( int, char*, double,double,double, double,double,double, int, char*, double );
void FieldEnergy( int, int, int,int,int );
void FieldForce( int, int );
void FieldForceThroughEnergy( int, int, char*, int,int,int );
void FieldTorqueThroughEnergy( int, int, char*, double,double,double, int,int,int );
void FocusingPotential( int, double,double,double, double,double,double, int );
void FocusingKickPer( int, double,double,double, double,double,double, double,int, double,double,double, double,int,double,int, const char*, int,int,double,double );
void FocusingKickPerFormStrRep( double*,double*,double*,double*,double*, int,int, double, int, const char* );

void ParticleTrajectory( int, double, double,double,double,double, double,double, int );
void FieldInt( int, char*, char*, double,double,double, double,double,double );
void CompCriterium( double, double, double, double, double,double );
void CompPrecisionOpt( const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char* );
void PhysicalUnits();
void RandomizationOnOrOff( char* );
void TolForConvergence( double, double, double );
void ShimSignature( int, char*, double,double,double, double,double,double, double,double,double, int, double,double,double );

void QuickDraw3D_ViewerOpt( int, const char*, const char*, const char* );
void OpenGL_3D_ViewerOpt( int, const char*, const char*, const char* );

void ApplyDrawAttrToElem( int, double,double,double, double );

void DeleteElement( int );
void DeleteAllElements1();
void InterruptTime( double );
void RadiaVersion();
//void DumpElem( int );
void DumpElemOpt( int*, int, const char* );
void DumpElemParseOpt( const unsigned char*, int );

void GenDump();

}

//-------------------------------------------------------------------------

extern radTIOBuffer ioBuffer;
//radTIOBuffer ioBuffer; //OC, to place back!!!

//-------------------------------------------------------------------------

//#ifdef WIN32
//
//extern HINSTANCE hinstCurrentRadia;
//extern HINSTANCE hinstPreviousRadia;
//extern LPSTR lpszCmdLineRadia;
//extern int nCmdShowRadia;
//
//BOOL APIENTRY DllMain(HANDLE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
//{
//	hinstCurrentRadia = (HINSTANCE)hModule;
//    return TRUE;
//}
//
//#endif

//-------------------------------------------------------------------------

int (*pgRadYieldExternFunc)() = 0;

//-------------------------------------------------------------------------

int CALL RadUtiYeldFuncSet(int (*pExtFunc)())
{
	if(pExtFunc != 0) 
	{
		pgRadYieldExternFunc = pExtFunc;
	}
	return OK;
}

//-------------------------------------------------------------------------
// Copied from AlpDllEntry.cpp
const char* CALL RadErrGet(int er)
{
	return ioBuffer.GetError(er);
}

int CALL RadErrGetSize(int* siz,int er)
{
	*siz= ioBuffer.GetErrorSize(er);
	return OK;
}

int CALL RadErrGetText(char* t,int er)
{
	//strcpy(t,ioBuffer.GetError(er));
	//OC02102018 (to avoid a need to have a separate call-back function for warning in an interface):
	if(er > 0) strcpy(t,ioBuffer.GetError(er));
	else strcpy(t,ioBuffer.GetWarning(er));
	return OK;
}

//-------------------------------------------------------------------------
// Copied from AlpDllEntry.cpp
const char* CALL RadWarGet(int er)
{
	return ioBuffer.GetWarning(er);
}

int CALL RadWarGetSize(int* siz,int er)
{
	*siz= ioBuffer.GetWarningSize(er);
	return OK;
}

int CALL RadWarGetText(char* t,int er)
{
	strcpy(t,ioBuffer.GetWarning(er));
	return OK;
}

//-------------------------------------------------------------------------

int CALL RadObjRecMag(int* n, double* pP, double* pL, double* pM)
{
	RecMag(pP[0], pP[1], pP[2], pL[0], pL[1], pL[2], pM[0], pM[1], pM[2]);
	
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjThckPgn(int* n, double xc, double lx, double* pFlatVertices, int NumVertices, char a, double* pM)
{
	ExtrudedPolygonDLL(xc, lx, pFlatVertices, NumVertices, a, pM);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjPolyhdr(int* n, double* pFlatVertices, int NumVertices, int* pFlatFaces, int* pFacesLengths, int NumFaces, double* pM, double* pM_LinCoef, double* pJ, double* pJ_LinCoef)
{
	PolyhedronDLL(pFlatVertices, NumVertices, pFlatFaces, pFacesLengths, NumFaces, pM, pM_LinCoef, pJ, pJ_LinCoef);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjMltExtPgn(int* n, double* pFlatVertices, int* pLayerLengths, double* pAttitudes, int NumLayers, double* pM)
{// pFlatVertices - flat array of 2d points
	MultGenExtrPolygonDLL(pFlatVertices, pLayerLengths, pAttitudes, NumLayers, pM);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjMltExtRtg(int* n, double* pFlatCenPts, double* pFlatRtgSizes, int NumLayers, double* pM)
{
	MultGenExtrRectangleDLL(pFlatCenPts, pFlatRtgSizes, NumLayers, pM);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjMltExtTri(int* n, double xc, double lx, double* pFlatVert, double* pFlatSubd, int nv, char a, double* pM, char* sOpt)
{
	const char *sOpt1=0, *sOpt2=0, *sOpt3=0, *sOpt4=0;
	vector<string> AuxStrings;
	if(sOpt != 0)
	{
		//char *SepStrArr[] = {";", ","};
		//CAuxParse::StringSplit(sOpt, SepStrArr, 2, " ", AuxStrings);
		//OC30072018 
		int lenStrOpt = (int)strlen(sOpt);
		char *sOptLoc = new char[lenStrOpt + 1];
		CAuxParse::StringSymbolsRemove(sOpt, (char*)" ", sOptLoc);
		CAuxParse::StringSplitNested(sOptLoc,";,", AuxStrings);
		delete[] sOptLoc;

		int AmOfTokens = (int)AuxStrings.size();
		if(AmOfTokens > 0) 
		{
			sOpt1 = (AuxStrings[0]).c_str();
			if(AmOfTokens > 1) 
			{
				sOpt2 = (AuxStrings[1]).c_str();
				if(AmOfTokens > 2) 
				{
					sOpt3 = (AuxStrings[2]).c_str();
					if(AmOfTokens > 3) sOpt4 = (AuxStrings[3]).c_str();
				}
			}
		}
	}
	
	MultGenExtrTriangleDLL(xc, lx, pFlatVert, pFlatSubd, nv, a, pM, sOpt1, sOpt2, sOpt3, sOpt4);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

//int CALL RadObjArcMag(int* n, double* P, double* R, double* Phi, double h, int nseg, char a, double* M)
//{
//	ArcMag(P[0], P[1], P[2], R[0], R[1], Phi[0], Phi[1], h, nseg, &a, M[0], M[1], M[2]);
//
//	*n = ioBuffer.OutInt();
//	return ioBuffer.OutErrorStatus();
//}

//-------------------------------------------------------------------------

int CALL RadObjArcPgnMag(int* n, double* P, char a, double* pFlatVert, int nv, double* Phi, int nseg, char sym_no, double* M)
{
	if((P == 0) || (pFlatVert == 0) || (Phi == 0) || (M == 0)) { ioBuffer.StoreErrorMessage("Radia::Error000"); *n = ioBuffer.OutInt(); return ioBuffer.OutErrorStatus();}
	//consider puting this everywhere

    ArcPolygonDLL(P[0], P[1], a, pFlatVert, nv, Phi[0], Phi[1], nseg, sym_no, M[0], M[1], M[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCylMag(int* n, double* P, double r, double h, int nseg, char a, double* M)
{
	CylMag(P[0], P[1], P[2], r, h, nseg, &a, M[0], M[1], M[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjRecCur(int* n, double* pP, double* pL, double* pJ)
{
	RecCur(pP[0], pP[1], pP[2], pL[0], pL[1], pL[2], pJ[0], pJ[1], pJ[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

//int CALL RadObjArcCur(int* n, double* pP, double* pR, double* pPhi, double h, double j, int nseg, char* cManOrAuto)
int CALL RadObjArcCur(int* n, double* pP, double* pR, double* pPhi, double h, int nseg, char man_auto, char a, double j)
{
	char strManAuto[] = "man\0  ";
	man_auto = (char)toupper(man_auto);
	if(man_auto == 'A') strcpy(strManAuto, "auto\0");

	ArcCur(pP[0], pP[1], pP[2], pR[0], pR[1], pPhi[0], pPhi[1], h, nseg, j, strManAuto, &a);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

//int CALL RadObjRaceTrk(int* n, double* pP, double* pR, double* pL, double h, double j, int nseg, char* cManOrAuto)
int CALL RadObjRaceTrk(int* n, double* pP, double* pR, double* pL, double h, int nseg, char man_auto, char a, double j)
{
	char strManAuto[] = "man\0  ";
	man_auto = (char)toupper(man_auto);
	if(man_auto == 'A') strcpy(strManAuto, "auto\0");

	RaceTrack(pP[0], pP[1], pP[2], pR[0], pR[1], pL[0], pL[1], h, nseg, j, strManAuto, &a);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjFlmCur(int* n, double* pPts, int np, double i)
{
	FlmCurDLL(pPts, np, i);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjScaleCur(int n, double scaleCoef)
{
	ScaleCurInObj(n, scaleCoef);
	ioBuffer.OutInt(); // to clear buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjBckg(int* n, double* pB)
{
	BackgroundFieldSource(pB[0], pB[1], pB[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCnt(int* n, int* pKeys, int NumKeys)
{
	Group(pKeys, NumKeys);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjAddToCnt(int Cnt, int* pKeys, int NumKeys)
{
	AddToGroup(Cnt, pKeys, NumKeys);

	Cnt = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCntSize(int* n, int Cnt)
{
	OutGroupSize(Cnt);
	//OutGroupSize(Cnt, deep);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCntStuf(int* pCntIndexes, int Cnt)
{
	OutGroupSubObjectKeys(Cnt);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20], NumDims;
	ioBuffer.OutMultiDimArrayOfInt(pCntIndexes, Dims, NumDims);
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjDpl(int* n, int Obj, char* Opt1)
{
	DuplicateElementG3DOpt(Obj, Opt1);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjM(double* pM, int* arMesh, int Obj) //OC21092018
//int CALL RadObjM(int* arMesh, int Obj) //OC15092018
//int CALL RadObjM(double* pM, int Obj)
{
	arMesh[0] = 0; //OC15092018

	MagnOfObj(Obj);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	//ioBuffer.OutMultiDimArrayOfDouble(0, Dims, NumDims); //OC15092018
	ioBuffer.OutMultiDimArrayOfDouble(pM, Dims, NumDims); //OC27092018

	if(arMesh != 0)
	{
		arMesh[0] = NumDims; //OC15092018
		for(int i=0; i<NumDims; i++) arMesh[i+1] = Dims[i];
	}

	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjCenFld(double* pM, int* arMesh, int Obj, char type) //OC27092018
//int CALL RadObjCenFld(int* arMesh, int Obj, char type) //OC22092018
//int CALL RadObjCenFld(double* pM, int Obj, char type)
{
	ObjField(Obj, &type);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	//ioBuffer.OutMultiDimArrayOfDouble(0, Dims, NumDims); //OC22092018
	ioBuffer.OutMultiDimArrayOfDouble(pM, Dims, NumDims); //OC27092018

	if(arMesh != 0)
	{
		arMesh[0] = NumDims; //OC22092018
		for(int i=0; i<NumDims; i++) arMesh[i+1] = Dims[i];
	}

	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjSetM(int obj, double* pM)
{
	SetObjMagn(obj, pM[0],pM[1],pM[2]);

	ioBuffer.OutInt(); // to clear buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjCutMag(int* pIndexes, int* pAmOfIndexes, int Obj, double* pP, double* pN, char* Opt1)
{
	CutElementG3DOpt(Obj, pP[0], pP[1], pP[2], pN[0], pN[1], pN[2], Opt1);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20], NumDims;
	ioBuffer.OutMultiDimArrayOfInt(pIndexes, Dims, NumDims);
	*pAmOfIndexes = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjDivMagPln(int* n, int Obj, double* pSbdPar, int nSbdPar, double* pFlatNormals, char* sOpt) //OC22092018
//int CALL RadObjDivMagPln(int* n, int Obj, double* pSbdPar, int nSbdPar, double* pFlatNormals, char* Opt)
{
	const char *sOpt1=0, *sOpt2=0;
	//const char *Opt1=0, *Opt2=0;
	vector<string> AuxStrings;
	if((sOpt != 0) && (strlen(sOpt) > 0)) //OC22092018
	//if(Opt != 0)
	{
		//char *SepStrArr[] = {(char*)";", (char*)","}; //OC04082018 (to please GCC 4.9)
		//CAuxParse::StringSplit(Opt, SepStrArr, 2, (char*)" ", AuxStrings);
		//int AmOfTokens = (int)AuxStrings.size();
		//if(AmOfTokens > 0) Opt1 = (AuxStrings[0]).c_str();
		//if(AmOfTokens > 1) Opt2 = (AuxStrings[1]).c_str();
		//OC22092018
		int lenStrOpt = (int)strlen(sOpt);
		char *sOptLoc = new char[lenStrOpt + 1];
		CAuxParse::StringSymbolsRemove(sOpt, (char*)" ", sOptLoc);
		CAuxParse::StringSplitNested(sOptLoc,";,", AuxStrings);
		delete[] sOptLoc;
		int AmOfTokens = (int)AuxStrings.size();
		if(AmOfTokens > 0) 
		{
			sOpt1 = (AuxStrings[0]).c_str();
			if(AmOfTokens > 1) 
			{
				sOpt2 = (AuxStrings[1]).c_str();
			}
		}
	}

	double AuxSbdPar[6];
	if(nSbdPar == 3)
	{
		AuxSbdPar[0] = pSbdPar[0]; AuxSbdPar[1] = 1.;
		AuxSbdPar[2] = pSbdPar[1]; AuxSbdPar[3] = 1.;
		AuxSbdPar[4] = pSbdPar[2]; AuxSbdPar[5] = 1.;
	}
	else if(nSbdPar == 6)
	{
		for(int i=0; i<6; i++) AuxSbdPar[i] = pSbdPar[i];
	}
	char TypeExtraSpec = 2; //pln; = 1;//cyl
	int LenExtraSpec = 9;

	SubdivideElementG3DOpt(Obj, AuxSbdPar, TypeExtraSpec, pFlatNormals, LenExtraSpec, sOpt1, sOpt2, "\0");
	//SubdivideElementG3DOpt(Obj, AuxSbdPar, TypeExtraSpec, pFlatNormals, LenExtraSpec, Opt1, Opt2, "\0");
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjDivMagCyl(int* n, int Obj, double* pSbdPar, int nSbdPar, double* pFlatCylDefPts, double Rat, char* Opt)
{
	const char *Opt1=0, *Opt2=0;
	vector<string> AuxStrings;
	if(Opt != 0)
	{
		//char *SepStrArr[] = {";", ","};
		char *SepStrArr[] = {(char*)";", (char*)","};  //OC04082018 (to please GCC 4.9)
		CAuxParse::StringSplit(Opt, SepStrArr, 2, (char*)" ", AuxStrings);
		int AmOfTokens = (int)AuxStrings.size();
		if(AmOfTokens > 0) Opt1 = (AuxStrings[0]).c_str();
		if(AmOfTokens > 1) Opt2 = (AuxStrings[1]).c_str();
	}

	double AuxSbdPar[6];
	if(nSbdPar == 3)
	{
		AuxSbdPar[0] = pSbdPar[0]; AuxSbdPar[1] = 1.;
		AuxSbdPar[2] = pSbdPar[1]; AuxSbdPar[3] = 1.;
		AuxSbdPar[4] = pSbdPar[2]; AuxSbdPar[5] = 1.;
	}
	else if(nSbdPar == 6)
	{
		for(int i=0; i<6; i++) AuxSbdPar[i] = pSbdPar[i];
	}
	char TypeExtraSpec = 1; //pln; = 1;//cyl
	int LenExtraSpec = 10;
	double AuxExtraSpec[10];
	for(int i=0; i<(LenExtraSpec - 1); i++) AuxExtraSpec[i] = pFlatCylDefPts[i];
	AuxExtraSpec[LenExtraSpec - 1] = Rat;

	SubdivideElementG3DOpt(Obj, AuxSbdPar, TypeExtraSpec, AuxExtraSpec, LenExtraSpec, Opt1, Opt2, "\0");
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjGeoVol(double* Vol, int Obj)
{
	GeometricalVolume(Obj);

	*Vol = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjGeoLim(double* Lim, int Obj)
{
	GeometricalLimits(Obj);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(Lim, Dims, NumDims);
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjDegFre(int* Num, int Obj)
{//03102018
	NumberOfDegOfFreedom(Obj);

	*Num = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfCmbL(int* n, int OrigTrf, int trf)
{
	CombineTransformLeft(OrigTrf, trf);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfCmbR(int* n, int OrigTrf, int trf)
{
	CombineTransformRight(OrigTrf, trf);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfInv(int* n)
{
	FieldInversion();

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfMlt(int* n, int obj, int trf, int mlt)
{
	ApplySymmetry(obj, trf, mlt);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfOrnt(int* n, int obj, int trf)
{
	TransformObject(obj, trf);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfPlSym(int* n, double* pP, double* pN)
{
	PlaneSym(pP[0], pP[1], pP[2], pN[0], pN[1], pN[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfRot(int* n, double* pP, double* pV, double phi)
{
	Rotation(pP[0], pP[1], pP[2], pV[0], pV[1], pV[2], phi);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadTrfTrsl(int* n, double* pV)
{
	Translation(pV[0], pV[1], pV[2]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatApl(int* n, int Obj, int Mat) 
{
	ApplyMaterial(Obj, Mat);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatLin(int* n, double* pKsi, double* pMr, int LenMr)
{
	if(LenMr == 3) LinearMaterial(pKsi[0], pKsi[1], pMr[0] , pMr[1], pMr[2]);
	else if(LenMr == 1) LinearMaterial2(pKsi[0], pKsi[1], pMr[0]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatStd(int* n, char* Name, double m)
{
	MaterialStd(Name, m);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatMvsH(double* pM, int* pNm, int Obj, char* id, double* pH)
{
	MvsH(Obj, id, pH[0], pH[1], pH[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pM, Dims, NumDims);
	*pNm = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadMatSatIsoFrm(int* n, double* pKsiMs1, double* pKsiMs2, double* pKsiMs3)
{
	//NonlinearIsotropMaterial2(pKsiMs1[0],pKsiMs1[1], pKsiMs2[0],pKsiMs2[1], pKsiMs3[0],pKsiMs3[1]);
	//OC03102018
	double KsiMs1_0=0, KsiMs1_1=0, KsiMs2_0=0, KsiMs2_1=0, KsiMs3_0=0, KsiMs3_1=0;
	if(pKsiMs1 != 0)
	{
		KsiMs1_0 = *pKsiMs1; KsiMs1_1 = *(pKsiMs1+1);
	}
	if(pKsiMs2 != 0)
	{
		KsiMs2_0 = *pKsiMs2; KsiMs2_1 = *(pKsiMs2+1);
	}
	if(pKsiMs3 != 0)
	{
		KsiMs3_0 = *pKsiMs3; KsiMs3_1 = *(pKsiMs3+1);
	}
	NonlinearIsotropMaterial2(KsiMs1_0,KsiMs1_1, KsiMs2_0,KsiMs2_1, KsiMs3_0,KsiMs3_1);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSatIsoTab(int* n, double* pFlatMatDef, int AmOfPts)
{
	double **PointsArray = new double*[AmOfPts];
	if(PointsArray == 0) 
	{ 
		ioBuffer.StoreErrorMessage("Radia::Error900"); return ioBuffer.OutErrorStatus();
	}
	double **tPointsArray = PointsArray;
	double *tCoord = pFlatMatDef;
	for(int i=0; i<AmOfPts; i++)
	{
		//*(tPointsArray++) = tCoord;
		//tCoord += 2;

		*tPointsArray = new double[2];
		(*tPointsArray)[0] = *(tCoord++);
		(*tPointsArray)[1] = *(tCoord++);
		tPointsArray++;
	}

	NonlinearIsotropMaterial3Opt(PointsArray, (long)AmOfPts);

	*n = ioBuffer.OutInt();
	if(PointsArray != 0) 
	{
		for(int i=0; i<AmOfPts; i++) if(PointsArray[i] != 0) delete[] (PointsArray[i]);
		delete[] PointsArray;
	}
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSatLamFrm(int* n, double* pKsiMs1, double* pKsiMs2, double* pKsiMs3, double p, double* N)
{
	NonlinearLaminatedMaterialFrm(pKsiMs1, pKsiMs2, pKsiMs3, p, N);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSatLamTab(int* n, double* pFlatMatDef, int AmOfMatPts, double p, double* N)
{
	NonlinearLaminatedMaterialTab(pFlatMatDef, AmOfMatPts, p, N);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadMatSatAniso(int* n, double* pDataPar, int LenDataPar, double* pDataPer, int LenDataPer)
{
	double **ParArray = 0;
	double **PerArray = 0;

	//if(LenDataPar == 11) // {ksi1,ms1,hc1,ksi2,ms2,hc2,ksi3,ms3,hc3,ksi0,hc0}
	//{
	//}
	if(LenDataPar == 8) // {ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0,hc}
	{
		ParArray = new double*[4];
		double **tParArray = ParArray;
		double *tCoord = pDataPar;
		for(int i=0; i<4; i++)
		{
			*(tParArray++) = tCoord;
			tCoord += 2;
		}
	}

	//if(LenDataPer == 8)
	if(LenDataPer == 7)
	{
		PerArray = new double*[4];
		double **tPerArray = PerArray;
		double *tCoord = pDataPer;
		for(int i=0; i<4; i++)
		{
			*(tPerArray++) = tCoord;
			tCoord += 2;
		}
	}

	if(LenDataPar == 11) NonlinearAnisotropMaterialOpt0(pDataPar, LenDataPar, pDataPer, LenDataPer);
	//else if((LenDataPar == 8) && (LenDataPer == 8)) NonlinearAnisotropMaterialOpt1(ParArray, PerArray);
	else if((LenDataPar == 8) && (LenDataPer == 7)) NonlinearAnisotropMaterialOpt1(ParArray, PerArray);
	else if((LenDataPar == 8) && (LenDataPer == 1)) NonlinearAnisotropMaterialOpt2(ParArray, pDataPer[0]);
	//else if((LenDataPar == 1) && (LenDataPer == 8)) NonlinearAnisotropMaterialOpt3(pDataPar[0], PerArray);
	else if((LenDataPar == 1) && (LenDataPer == 7)) NonlinearAnisotropMaterialOpt3(pDataPar[0], PerArray);
	else
	{
        ioBuffer.StoreErrorMessage("Radia::Error000"); 
		return ioBuffer.OutErrorStatus();
	}

	*n = ioBuffer.OutInt();
	if(ParArray != 0) delete[] ParArray;
	if(PerArray != 0) delete[] PerArray;
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadRlxPre(int* n, int Obj, int SrcObj)
{
	PreRelax(Obj, SrcObj);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadRlxMan(double* dOut, int* nOut, int Intrc, int Meth, int IterNum, double RlxPar)
{
	ManualRelax(Intrc, Meth, IterNum, RlxPar);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(dOut, Dims, NumDims);
	*nOut = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadRlxAuto(double* dOut, int* nOut, int Intrc, double Prec, int MaxIter, int Meth, const char* Opt1)
{
	//AutoRelax(Intrc, Prec, MaxIter, Meth);
	AutoRelaxOpt(Intrc, Prec, MaxIter, Meth, Opt1); //OC240408

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(dOut, Dims, NumDims);
	*nOut = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

EXP int CALL RadRlxUpdSrc(int intrc)
{
	UpdateSourcesForRelax(intrc);

	ioBuffer.OutInt(); // to clean buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFld(double* pB, int* pNb, int Obj, char* ID, double* pCoord, int Np)
{
	double **PointsArray = new double*[Np];
	if(PointsArray == 0) 
	{ 
		ioBuffer.StoreErrorMessage("Radia::Error900"); return ioBuffer.OutErrorStatus();
	}
	double **tPointsArray = PointsArray;
	double *tCoord = pCoord;
	for(int i=0; i<Np; i++)
	{
		*(tPointsArray++) = tCoord;
		tCoord += 3;
	}

	FieldArbitraryPointsArray((long)Obj, ID, PointsArray, (long)Np);
	
	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0)
	{
        if(PointsArray != 0) delete[] PointsArray;
		return ErrStat;
	}

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pB, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNb = TotLen;

	if(PointsArray != 0) delete[] PointsArray;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldEnr(double* d, int objdst, int objsrc, int* SbdPar)
{
	int LocSbdArr[] = {0,0,0};
	if(SbdPar != 0) { LocSbdArr[0] = SbdPar[0]; LocSbdArr[1] = SbdPar[1]; LocSbdArr[2] = SbdPar[2];}

	FieldEnergy(objdst, objsrc, LocSbdArr[0], LocSbdArr[1], LocSbdArr[2]);

	*d = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldEnrFrc(double* pF, int* pNf, int objdst, int objsrc, char* id, int* SbdPar)
{
	int LocSbdArr[] = {0,0,0};
	if(SbdPar != 0) { LocSbdArr[0] = SbdPar[0]; LocSbdArr[1] = SbdPar[1]; LocSbdArr[2] = SbdPar[2];}

	FieldForceThroughEnergy(objdst, objsrc, id, LocSbdArr[0], LocSbdArr[1], LocSbdArr[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);
	*pNf = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldEnrTrq(double* pF, int* pNf, int objdst, int objsrc, char* id, double* pP, int* SbdPar)
{
	int LocSbdArr[] = {0,0,0};
	if(SbdPar != 0) { LocSbdArr[0] = SbdPar[0]; LocSbdArr[1] = SbdPar[1]; LocSbdArr[2] = SbdPar[2];}

	FieldTorqueThroughEnergy(objdst, objsrc, id, pP[0], pP[1], pP[2], LocSbdArr[0], LocSbdArr[1], LocSbdArr[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);
	*pNf = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldFrc(double* pF, int* pNf, int Obj, int Shape)
{
	FieldForce(Obj, Shape);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);
	*pNf = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldFocPot(double* d, int Obj, double* P1, double* P2, int np)
{
	FocusingPotential(Obj, P1[0], P1[1], P1[2], P2[0], P2[1], P2[2], np);

	*d = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldFocKickPer(double* pMatr1, double* pMatr2, double* pIntBtrE2, double* pArg1, double* pArg2, int* psize, int obj, double* P1, double* Ns, double per, int nper, int nps, double* Ntr, double r1, int np1, double d1, double r2, int np2, double d2, int nh, char* com)
{
    FocusingKickPer(obj, P1[0], P1[1], P1[2], Ns[0], Ns[1], Ns[2], per, nper, Ntr[0], Ntr[1], Ntr[2], r1, np1, r2, np2, com, nh, nps, d1, d2);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	long NpKick = np1*np2;
	//long LenArr = 2*np1*np2 + np1 + np2 + 1;
	long LenArr = 3*np1*np2 + np1 + np2 + 1;

	double* pAuxBuf = new double[LenArr]; 
	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pAuxBuf, Dims, NumDims);

	double *tAuxBuf = pAuxBuf;

	double *tMatr1 = pMatr1;
	double *tMatr2 = pMatr2;
	double *tIntBtrE2 = pIntBtrE2;

	double *tArg1 = pArg1;
	double *tArg2 = pArg2;

	for(long i=0; i<NpKick; i++) *(tMatr1++) = *(tAuxBuf++);
	for(long j=0; j<NpKick; j++) *(tMatr2++) = *(tAuxBuf++);
	for(long ii=0; ii<NpKick; ii++) *(tIntBtrE2++) = *(tAuxBuf++);

	for(long k=0; k<np1; k++) *(tArg1++) = *(tAuxBuf++);
	for(long m=0; m<np2; m++) *(tArg2++) = *(tAuxBuf++);
	*psize = (int)(*tAuxBuf + 1.E-10);

	if(pAuxBuf != 0) delete[] pAuxBuf;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldFocKickPerFormStr(char* OutStr, double* pMatr1, double* pMatr2, double* pIntBtrE2, double* pArg1, double* pArg2, int np1, int np2, double per, int nper, char* com)
{
	FocusingKickPerFormStrRep(pMatr1, pMatr2, pIntBtrE2, pArg1, pArg2, np1, np2, per, nper, com);
	
	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	ioBuffer.OutStringClean(OutStr);
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldInt(double* pF, int* pNf, int Obj, char* InfOrFin, char* id, double* P1, double* P2)
{
	FieldInt(Obj, InfOrFin, id, P1[0], P1[1], P1[2], P2[0], P2[1], P2[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNf = TotLen;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldLst(double* pF, int* pNf, int Obj, char* id, double* P1, double* P2, int np, char* ArgOrNoArg, double Strt)
{
	Field(Obj, id, P1[0], P1[1], P1[2], P2[0], P2[1], P2[2], np, ArgOrNoArg, Strt);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNf = TotLen;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldPtcTrj(double* pF, int* pNf, int Obj, double E, double* pIC, double* pIL, int np)
{
	ParticleTrajectory(Obj, E, pIC[0], pIC[1], pIC[2], pIC[3], pIL[0], pIL[1], np);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNf = TotLen;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldShimSig(double* pF, int* pNf, int obj, char* id, double* V, double* P1, double* P2, int np, double* inVi)
{
	double Vi[] = {0,0,0};
	if(inVi != 0) { Vi[0] = inVi[0]; Vi[1] = inVi[1]; Vi[2] = inVi[2];}

    ShimSignature(obj, id, V[0], V[1], V[2], P1[0], P1[1], P1[2], P2[0], P2[1], P2[2], np, Vi[0], Vi[1], Vi[2]);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(pF, Dims, NumDims);

	int TotLen = 1;
	for(int k=0; k<NumDims; k++) TotLen *= Dims[k];
	*pNf = TotLen;
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldCmpCrt(int* n, double PrcB, double PrcA, double PrcBInt, double PrcFrc, double PrcTrjCrd, double PrcTrjAng)
{
	CompCriterium(PrcB, PrcA, PrcBInt, PrcFrc, PrcTrjCrd, PrcTrjAng);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldCmpPrc(int* n, char* Opt)
{
	const char *Opt1=0, *Opt2=0, *Opt3=0, *Opt4=0, *Opt5=0, *Opt6=0, *Opt7=0, *Opt8=0;
	vector<string> AuxStrings;
	if(Opt != 0)
	{
		//char *SepStrArr[] = {";", ","};
		char *SepStrArr[] = {(char*)";", (char*)","}; //OC04082018 (to please GCC 4.9)
		CAuxParse::StringSplit(Opt, SepStrArr, 2, (char*)" ", AuxStrings);
		int AmOfTokens = (int)AuxStrings.size();
		if(AmOfTokens > 0) Opt1 = (AuxStrings[0]).c_str();
		if(AmOfTokens > 1) Opt2 = (AuxStrings[1]).c_str();
		if(AmOfTokens > 2) Opt3 = (AuxStrings[2]).c_str();
		if(AmOfTokens > 3) Opt4 = (AuxStrings[3]).c_str();
		if(AmOfTokens > 4) Opt5 = (AuxStrings[4]).c_str();
		if(AmOfTokens > 5) Opt6 = (AuxStrings[5]).c_str();
		if(AmOfTokens > 6) Opt7 = (AuxStrings[6]).c_str();
		if(AmOfTokens > 7) Opt8 = (AuxStrings[7]).c_str();
	}

	CompPrecisionOpt(Opt1, Opt2, Opt3, Opt4, Opt5, Opt6, Opt7, Opt8);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldUnits(char* OutStr)
{
	PhysicalUnits();

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	strcpy(OutStr, ioBuffer.OutStringPtr()); //OC27092018
	//strcpy(OutStr, ioBuffer.OutString());
	ioBuffer.EraseStringBuffer();
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldUnitsSize(int* OutSize)
{
	PhysicalUnits();

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	*OutSize = (int)strlen(ioBuffer.OutStringPtr()); //27092018
	//*OutSize = (int)strlen(ioBuffer.OutString());
	ioBuffer.EraseStringBuffer();
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadFldFrcShpRtg(int* n, double* pP, double* pW)
{
	Rectngl(pP[0], pP[1], pP[2], pW[0], pW[1]);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldLenRndSw(int* n, char* OnOrOff)
{
	RandomizationOnOrOff(OnOrOff);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadFldLenTol(int* n, double AbsVal, double RelVal, double ZeroVal)
{
	TolForConvergence(AbsVal, RelVal, ZeroVal);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjDrwQD3D(int Obj, char* Opt)
{
	const char *Opt1=0, *Opt2=0, *Opt3=0;
	vector<string> AuxStrings;
	if(Opt != 0)
	{
		//char *SepStrArr[] = {";", ","};
		char *SepStrArr[] = {(char*)";", (char*)","}; //OC04082018 (to please GCC 4.9)
		CAuxParse::StringSplit(Opt, SepStrArr, 2, (char*)" ", AuxStrings);
		int AmOfTokens = (int)AuxStrings.size();
		if(AmOfTokens > 0) Opt1 = (AuxStrings[0]).c_str();
		if(AmOfTokens > 1) Opt2 = (AuxStrings[1]).c_str();
		if(AmOfTokens > 2) Opt3 = (AuxStrings[2]).c_str();
	}

	QuickDraw3D_ViewerOpt(Obj, Opt1, Opt2, Opt3);
	ioBuffer.OutInt(); // to clear buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjDrwOpenGL(int Obj, char* Opt)
{
	const char *Opt1=0, *Opt2=0, *Opt3=0;
	vector<string> AuxStrings;
	if(Opt != 0)
	{
		char *SepStrArr[] = {(char*)";", (char*)","};
		CAuxParse::StringSplit(Opt, SepStrArr, 2, (char*)" ", AuxStrings);
		int AmOfTokens = (int)AuxStrings.size();
		if(AmOfTokens > 0) 
		{
			Opt1 = (AuxStrings[0]).c_str();
			if(AmOfTokens > 1) 
			{
				Opt2 = (AuxStrings[1]).c_str();
				if(AmOfTokens > 2) Opt3 = (AuxStrings[2]).c_str();
			}
		}
	}

	OpenGL_3D_ViewerOpt(Obj, Opt1, Opt2, Opt3);
	ioBuffer.OutInt(); // to clear buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadObjDrwAtr(int Obj, double* pRGB, double Thcn)
{
	ApplyDrawAttrToElem(Obj, *pRGB, *(pRGB + 1), *(pRGB + 2), Thcn);
	ioBuffer.OutInt(); // to clear buffer
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiDel(int* n, int Elem)
{
	DeleteElement(Elem);

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiDelAll(int* n)
{
	DeleteAllElements1();

	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiDmp(char* OutStr, int* pSize, int* arElem, int nElem, char* AscOrBin) //OC27092018
//int CALL RadUtiDmp(char* OutStr, int* arElem, int nElem, char* AscOrBin)
//int CALL RadUtiDmp(char* OutStr, int Elem)
{
	//DumpElem(Elem);
	DumpElemOpt(arElem, nElem, AscOrBin); //OC230713

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;


	if((strcmp(AscOrBin, "asc\0") == 0) || (strcmp(AscOrBin, "Asc\0") == 0) || (strcmp(AscOrBin, "ASC\0") == 0)) 
	{
		if(pSize != 0) *pSize = (int)strlen(ioBuffer.OutStringPtr()) + 1; //to include terminating '\0' (?)
		//if(pSize != 0) *pSize = (int)strlen(ioBuffer.OutStringPtr()); //OC27092018
		if(OutStr != 0) ioBuffer.OutStringClean(OutStr); 
		//strcpy(OutStr, ioBuffer.OutString());
	}
	else 
	{
		if(pSize != 0) *pSize = (int)ioBuffer.OutByteStringSize(); //27092018
		if(OutStr != 0) ioBuffer.OutByteStringClean(OutStr); //27092018
	}
	//{
	//	long sizeData = ioBuffer.OutByteStringSize();
	//	const char *tData = ioBuffer.OutByteStringPtr(); //27092018
	//	//const char *tData = ioBuffer.OutByteString();
	//	char *tOutStr = OutStr;
	//	for(long i=0; i<sizeData; i++) *(tOutStr++) = *(tData++);
	//}

	//ioBuffer.EraseStringBuffer(); //in any case
	return ErrStat;
}

//-------------------------------------------------------------------------

//int CALL RadUtiDmpRead(char* OutStr, char* AscOrBin)
//{
//	int ErrStat = ioBuffer.OutErrorStatus();
//	if(ErrStat > 0) return ErrStat;
//
//	if((strcmp(AscOrBin, "asc\0") == 0) || (strcmp(AscOrBin, "Asc\0") == 0) || (strcmp(AscOrBin, "ASC\0") == 0))
//		strcpy(OutStr, ioBuffer.OutStringPtr()); //27092018
//		//strcpy(OutStr, ioBuffer.OutString());
//	else
//	{
//		long sizeData = ioBuffer.OutByteStringSize();
//		const char *tData = ioBuffer.OutByteStringPtr(); //27092018
//		//const char *tData = ioBuffer.OutByteString();
//		char *tOutStr = OutStr;
//		for(long i=0; i<sizeData; i++) *(tOutStr++) = *(tData++);
//	}
//
//	ioBuffer.EraseStringBuffer(); //in any case
//	return ErrStat;
//}

//-------------------------------------------------------------------------

////int CALL RadUtiDmpSize(int* OutSize, int Elem)
//int CALL RadUtiDmpSize(int* OutSize, int* arElem, int nElem, char* AscOrBin, bool doEraseBuf)
//{
//	//DumpElem(Elem);
//	DumpElemOpt(arElem, nElem, AscOrBin); //OC230713
//
//	int ErrStat = ioBuffer.OutErrorStatus();
//	if(ErrStat > 0) return ErrStat;
//
//	if((strcmp(AscOrBin, "asc\0") == 0) || (strcmp(AscOrBin, "Asc\0") == 0) || (strcmp(AscOrBin, "ASC\0") == 0))
//		*OutSize = (int)strlen(ioBuffer.OutStringPtr()); //OC27092018
//		//*OutSize = (int)strlen(ioBuffer.OutString());
//	else *OutSize = (int)ioBuffer.OutByteStringSize();
//
//	if(doEraseBuf) ioBuffer.EraseStringBuffer();
//	//leaving buffer not erased may be interesting if immediately after this RadUtiDmp will be called (then the DumpElemOpt(..) call won't need to be repeated)
//	return ErrStat;
//}

//-------------------------------------------------------------------------

int CALL RadUtiDmpPrs(int* arElem, int* pnElem, unsigned char* sBytes, int nBytes)
{//OC01102018
	DumpElemParseOpt(sBytes, nBytes);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	bool resIsList = (bool)sBytes[0];
	if(pnElem != 0)
	{
		*pnElem = 0;
		if(resIsList) 
		{
			int arDims[20], nDims=0;
			if(arElem == 0) ioBuffer.OutMultiDimArrayOfIntDims(arDims, nDims);
			else ioBuffer.OutMultiDimArrayOfInt(arElem, arDims, nDims);
			*pnElem = arDims[0]; //output can only be 1D array in this case
		}
		else 
		{
			*pnElem = 1;
			if(arElem != 0)
			{
				*arElem = ioBuffer.OutInt();
			}
		}
	}
	else if(arElem != 0)
	{
		if(resIsList) 
		{
			int arDims[20], nDims=0;
			ioBuffer.OutMultiDimArrayOfInt(arElem, arDims, nDims);
		}
		else *arElem = ioBuffer.OutInt();
	}
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadUtiIntrptTim(double* d, double t)
{
	InterruptTime(t);

	*d = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------

int CALL RadUtiVer(double* d)
{
	RadiaVersion();

	*d = ioBuffer.OutDouble();
	return ioBuffer.OutErrorStatus();
}

//-------------------------------------------------------------------------
// "Secondary" functions
//-------------------------------------------------------------------------

int CALL RadTrfZerPara(int* n, int obj, double* P, double* N)
{
	*n = 0;
	int OutRes = 0, LocRes = 0;
	int sym, inv, trf;

	if((LocRes = RadTrfPlSym(&sym, P, N)) > 0) return LocRes;
	if(LocRes < 0) OutRes = LocRes;

	if((LocRes = RadTrfInv(&inv)) > 0) return LocRes;
	if((LocRes < 0) && (OutRes == 0)) OutRes = LocRes;

	if((LocRes = RadTrfCmbL(&trf, sym, inv)) > 0) return LocRes;
	if((LocRes < 0) && (OutRes == 0)) OutRes = LocRes;

	if((LocRes = RadTrfMlt(n, obj, trf, 2)) > 0) return LocRes;
	if((LocRes < 0) && (OutRes == 0)) OutRes = LocRes;

	return OutRes;
}

//-------------------------------------------------------------------------

int CALL RadTrfZerPerp(int* n, int obj, double* P, double* N)
{
	*n = 0;
	int OutRes = 0, LocRes = 0;
	int sym;

	if((LocRes = RadTrfPlSym(&sym, P, N)) > 0) return LocRes;
	if(LocRes < 0) OutRes = LocRes;

	if((LocRes = RadTrfMlt(n, obj, sym, 2)) > 0) return LocRes;
	if((LocRes < 0) && (OutRes == 0)) OutRes = LocRes;

	return OutRes;
}

//-------------------------------------------------------------------------

int CALL RadSolve(double* dOut, int* nOut, int obj, double prec, int iter, int meth)
{
	SolveGen(obj, prec, iter, meth);

	int ErrStat = ioBuffer.OutErrorStatus();
	if(ErrStat > 0) return ErrStat;

	int Dims[20];
	int NumDims;
	ioBuffer.OutMultiDimArrayOfDouble(dOut, Dims, NumDims);
	*nOut = Dims[0];
	return ErrStat;
}

//-------------------------------------------------------------------------

int CALL RadObjFullMag(int* n, double* pP, double* pL, double* pM, double* SbdPar, int nSbdPar, int grp, int mat, double* pRGB)
{
	int OutRes = 0, LocRes = 0;

	if(LocRes = RadObjRecMag(n, pP, pL, pM)) return LocRes;
	if(LocRes < 0) OutRes = LocRes;

	if(LocRes = RadObjDrwAtr(*n, pRGB, 0.001)) return LocRes;
	if(LocRes < 0) OutRes = LocRes;

	if(mat > 0)
	{
		if(LocRes = RadMatApl(n, *n, mat)) return LocRes;
		if(LocRes < 0) OutRes = LocRes;
	}

	double FlatNorm[] = {1,0,0,0,1,0,0,0,1};
	if(LocRes = RadObjDivMagPln(n, *n, SbdPar, nSbdPar, FlatNorm, 0)) return LocRes;
	if(LocRes < 0) OutRes = LocRes;

	if(grp > 0)
	{
		int KeysArr[] = {*n};
		if(LocRes = RadObjAddToCnt(grp, KeysArr, 1)) return LocRes;
		if(LocRes < 0) OutRes = LocRes;
	}

	return OutRes;
}

//-------------------------------------------------------------------------

int CALL RadUtiDataGet(char* pcData, const char typeData[3], long key) //OC04102018
//int CALL RadUtiDataGet(char* pcData, char typeData[3], long key) //OC27092018
//int CALL RadUtiDataGet(double* pData, long key) //OC15092018
{
	ioBuffer.OutDataClean(pcData, typeData, key);

	int ErrStat = ioBuffer.OutErrorStatus();
	return ErrStat;
}

//-------------------------------------------------------------------------
