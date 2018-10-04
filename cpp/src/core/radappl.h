/*-------------------------------------------------------------------------
*
* File name:      radapl.h
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

#ifndef __RADAPPL_H
#define __RADAPPL_H

#include "radsend.h"
#include "radg.h"
#include "radcast.h"
#include "radg3dgr.h"
#include "radyield.h"
#include "radcnvrg.h"
#include "radapl1.h"
#include "radauxst.h"

#include <sstream> // Porting

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct TVector2d;

//-------------------------------------------------------------------------

class radTApplication {

	radTmhg GlobalMapOfHandlers;
	int GlobalUniqueMapKey;
	radTMapOfDrawAttr MapOfDrawAttr;
	radTCast Cast;
	radTCompCriterium CompCriterium;

	char QD3D_ViewerWasInitialized;

public:

	short SendingIsRequired;

	radTSend Send;
	radTConvergRepair CnRep, CnRepAux;
	short TreatRecMagsAsExtrPolygons, TreatRecMagsAsPolyhedrons, RecognizeRecMagsInPolyhedrons, TreatExtrPgnsAsPolyhedrons;
	short MemAllocForIntrctMatrTotAtOnce;

	radTApplication()
	{
		Initialize();
	}
	~radTApplication() {}

	void Initialize()
	{
		GlobalUniqueMapKey = 1;
		CompCriterium.BasedOnPrecLevel = 0;
		SendingIsRequired = 1;
		TreatRecMagsAsExtrPolygons = TreatRecMagsAsPolyhedrons = TreatExtrPgnsAsPolyhedrons = 0;
		RecognizeRecMagsInPolyhedrons = 1; // If possible, of course
		MemAllocForIntrctMatrTotAtOnce = 0;

		QD3D_ViewerWasInitialized = 0;
	}

	int ValidateVector3d(double* ArrayToCheck, long LenArray, TVector3d* VectorPtr);
	int ValidateVector2d(double* ArrayToCheck, long LenArray, TVector2d* VectorPtr);
	int ValidateMatrix3d(double* arToCheck, long LenAr, TMatrix3d* MatrixPtr);

	int ValidateElemKey(long ElemKey, radThg& hg);
	int ValidateFieldChar(char* FieldChar, radTFieldKey* FieldKeyPtr, bool LocSendRequired = true);
	int ValidateFieldIntChar(char* FieldIntChar, char* FinOrInfChar, radTFieldKey* FieldKeyPtr, bool LocSendRequired = true);
	int ValidateFieldEnergyForceChar(char* ComponIDChar, radTFieldKey* FieldKeyPtr);
	int ValidateMagnChar(char* MagnChar);
	int ValidateForceChar(char* ForceChar);
	int ValidateTorqueChar(char* TorqueChar);
	int ValidateIsotropMaterDescrByPoints(TVector2d* ArrayHM, int LenArrayArrayHM);

	inline int AddElementToContainer(radThg& hg);

	int SetRecMag(double* CPoi, long lenCPoi, double* Dims, long lenDims, double* Magn, long lenMagn, double* J, long lenJ, short J_IsZero);
	int SetArcCur(double* CPoi, long lenCPoi, double* Radii, long lenRadii, double* Angles, long lenAngles, double InHeight, double InJ_azim, int NumberOfSegm, char* ManOrAuto, char* Orient);
	int SetArcMag(double* CPoi, long lenCPoi, double* Radii, long lenRadii, double* Angles, long lenAngles, double InHeight, int InNumberOfSegm, double* Magn, long lenMagn, char* Orient);
	int SetCylMag(double* CPoi, long lenCPoi, double r, double h, int NumberOfSegm, double* Magn, long lenMagn, char* Orient);
	//int OrientObjAlongMainAxis(int ObjInd, double* CPoi, char DefOrient, char Orient);
	int FindSpaceTransToOrientObjAlongMainAxis(double* CPoi, char DefOrient, char Orient);
	void TransformBackMagnOrCurDensArr(int IndTr, double* Magn, long lenMagn);
	void TransformBackPointArr(int IndTr, double* arP, long lenP);

	int SetExtrudedPolygon(double* FirstPoi, long lenFirstPoi, double, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* Magn, long lenMagn, const char* OrientStr);
	inline int CheckIfExtrudedPolygonIsRecMag(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d);
	int SetPlanarPolygon(double CoordZ, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* Magn, long lenMagn);
	
	int SetPolyhedron1(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfNumOfPoInFaces, int lenArrayOfFaces, double* Magn, double* arM_LinCoef=0, double* J=0, double* arJ_LinCoef=0, const char** OptionNames=0, const char** OptionValues=0, int OptionCount=0);
	//int SetPolyhedron1(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfNumOfPoInFaces, int lenArrayOfFaces, double* Magn, long lenMagn);
	
	int SetPolyhedron2(TVector3d** ArrayOfFaces, int* ArrayOfNumOfPoInFaces, long lenArrayOfFaces, double* Magn, long lenMagn);
	int SetArcPolygon(double* CenP, const char* OrientStr, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* Angles, int NumberOfSegm, const char* SymOrNoSymStr, double* Magn);

	int SetMultGenExtrPolygon(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons, double* Magn, long lenMagn);
	//int SetMultGenExtrPolygonCur(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, const char** arOptionNames=0, const char** arOptionValues=0, int numOptions=0);
	int SetMultGenExtrPolygonCur(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double* arSubdData, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInSteps, const char** arOptionNames=0, const char** arOptionValues=0, int numOptions=0);
	
	int SetMultGenExtrRectangle(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect, double* Magn, long lenMagn);
	int SetMultGenExtrTriangle(double* FirstPoi, long lenFirstPoi, double Lx, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* arSubdData, double* Magn, long lenMagn, const char* OrientStr, const char** OptionNames, const char** OptionValues, int OptionCount);
	//int TriangulatePolygon(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* arSubdData, const char** OptionNames, const char** OptionValues, int OptionCount, TVector2d*& arTriVertPt, int& numTriVertPt, int*& arTriVertInd, int& numTri);
	int TriangulatePolygon(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* arSubdData, char triSubdParamBorderCode, double triAngMin, double triAreaMax, const char* sTriExtOpt, TVector2d*& arTriVertPt, int& numTriVertPt, int*& arTriVertInd, int& numTri);
	
	//int SetUpPolyhedronsFromBaseFacePolygons(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, char frame, radThg& hgOut);
	int SetUpPolyhedronsFromBaseFacePolygons(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInSteps, char frame, radThg& hgOut);
	int SetUpPolyhedronsFromBaseFacePolygonsTri(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, TVector2d* arTriVertPt, int numTriVertPt, int* arTriVertInd, int numTri, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInExtrSteps, char frame, radThg& hgOut);

	int SetUpPolyhedronsFromLayerPolygons(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons, TVector3d& Magn, radThg& hg);
	int SetUpPolyhedronsFromLayerRectangles(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect, TVector3d& MagnVect, radThg& hg);
	int CheckLayerPolygonStructures(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons);
	int CheckLayerRectangleStructures(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect);
	int SetUpOnePolyhedronSegment(radTPtrsToPgnAndVect2d* pPtrsToPgnAndVect2d, double* z1z2, char StageChar, double* RelAbsTol, TVector3d& Magn, radTVectVect3d* pVertexPointsVect, radTVectIntPtrAndInt* pFacesVect, radThg& hgLoc);
	int CheckIfGroupIsNeeded(radThg& In_hg);
	int SetUpTetrahedronBasedOnTwoLinSegm(radTVect2dVect* pFirstVect2dVect, double z1, radTVect2dVect* pSecondVect2dVect, double z2, double* RelAbsTol, TVector3d& Magn, radThg& In_hg);
	int FindLowestPoint(radTVect2dVect* pVectP2d, TVector2d& V, double* RelAbsTol, int& LowestPointInd, char& AmOfPo);
	int FindTwoAdjacentFaces(int OneVertPoInd, int AnotherVertPoInd, radTVectIntPtrAndInt* pFacesVect, int& OneFaceInd, int& IndOfPoOnOneFace, int& AnotherFaceInd, int& IndOfPoOnAnotherFace);
	int NextCircularNumber(int CurrentNo, int Total) { return (CurrentNo == Total-1)? 0 : CurrentNo + 1;}
	int ShiftVertexPointNumbersInFaces(radTVectIntPtrAndInt* pFacesVect, int AmOfPoToDelete, char);

	int RecMagsAsExtrPolygons(char* OnOrOff);
	int RecMagsAsPolyhedrons(char* OnOrOff);
	int RecognizeRecMags(char* OnOrOff);
	int ExtPgnsAsPolyhedrons(char* OnOrOff);

	int SetGroup(int* ArrayOfKeys, long lenArrayOfKeys);
	int AddToGroup(int GroupKey, int* ArrayOfKeys, long lenArrayOfKeys);
	int OutGroupSize(int ElemKey);
	//int OutGroupSize(int ElemKey, char deep=0);
	int OutGroupSubObjectKeys(int ElemKey);

	int SetRaceTrack(double* CPoi, long lenCPoi, double* Radii, long lenRadii, double* StrPartDims, long lenStrPartDims, double InHeight, double InJ_azim, int NumberOfSegm, char* ManOrAuto, char* Orient);
	int SetFlmCur(double I, TVector3d* ArrayOfPoints, int lenArrayOfPoints);
	int SetRectangle(double* CPoi, long lenCPoi, double* Dims, long lenDims);
	int SetBackgroundFieldSource(double* B, long lenB);

	int ComputeNumberOfDegOfFreedom(int ElemKey);
	int ComputeGeometricalVolume(int ElemKey);
	int ComputeGeometricalLimits(int ElemKey);
	//void ComputeMagnInCenter(int ElemKey);
	void ComputeMagnOrJ_InCenter(int ElemKey, char MorJ);
	int ScaleCurrent(int ElemKey, double scaleCoef);
	int SetObjMagn(int ElemKey, double Mx, double My, double Mz);
	void OutCenFieldCompRes(radTVectPairOfVect3d*);

	int FieldCompMethForSubdividedRecMag(int ElemKey, int Switch, int SubLevel);
	int SetLocMgnInSbdRecMag(int ElemKey, TVector3d* ArrayOfVectIndx, TVector3d* ArrayOfMagn, int Len);

	int SetTranslation(double* Transl, long lenTransl);
	int SetRotation(double* PoiOnAx, long lenPoiOnAx, double* AxVect, long lenAxVect, double Angle);
	int SetPlaneSym(double* PoiOnPlane, long lenPoiOnPlane, double* PlaneNormal, long lenPlaneNormal, int s);
	int SetFieldInversion();

	int CombineTransformations(int ThisElemKey, int AnotherElemKey, char L_or_R);
	int ApplySymmetry(int g3dElemKey, int TransElemKey, int Multiplicity);

	int SetLinearMaterial(double* KsiArray, long lenKsiArray, double* RemMagnArray, long lenRemMagnArray);
	int SetMaterialStd(char* MatName, double Mr=0);

	int SetNonlinearIsotropMaterial(double* Ms, long lenMs, double* ks, long len_ks);
	int SetNonlinearIsotropMaterial(TVector2d* ArrayHM, int LenArrayArrayHM);
    int SetNonlinearLaminatedMaterial(TVector2d* ArrayOfPoints2d, int lenArrayOfPoints2d, double PackFactor, double* dN);

	int SetNonlinearAnisotropMaterial(double** Ksi, double** Ms, double* Hc, int lenHc, char* DependenceIsNonlinear);
	int SetNonlinearAnisotropMaterial0(double* pDataPar, int lenDataPar, double* pDataPer, int lenDataPer);

	int ApplyMaterial(int g3dRelaxElemKey, int MaterElemKey);
	void ComputeMvsH(int g3dRelaxOrMaterElemKey, char* MagnChar, double* H, long lenH);
	//void OutMagnetizCompRes(char* MagnChar, TVector3d& M_vect);

	int PreRelax(int ElemKey, int SrcElemKey);
	void ShowInteractMatrix(int InteractElemKey);
	void ShowInteractVector(int InteractElemKey, char* FieldVectID);
	int MakeManualRelax(int InteractElemKey, int MethNo, int IterNumber, double RelaxParam);
	int MakeAutoRelax(int InteractElemKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, const char** arOptionNames=0, const char** arOptionValues=0, int numOptions=0);
	int UpdateSourcesForRelax(int InteractElemKey);
	int SolveGen(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo);

	void ComputeField(int ElemKey, char* FieldChar, double* StObsPoi, long lenStObsPoi, double* FiObsPoi, long lenFiObsPoi, int Np, char* ShowArgFlag, double StrtArg);
	void ComputeField(int ElemKey, char* FieldChar, radTVectorOfVector3d& VectorOfVector3d, radTVectInputCell& VectInputCell);
	void ComputeField(int ElemKey, char* FieldChar, double** Points, long LenPoints);

	void ComputeFieldInt(int ElemKey, char* IntID, char* FieldIntChar, double* StPoi, long lenStPoi, double* FiPoi, long lenFiPoi);
	void ComputeFieldForce(int ElemKey, int ShapeElemKey);
	void ComputeFieldEnergy(int DestElemKey, int SourceElemKey, int* SubdArray, long lenSubdivArray);
	void ComputeFieldForceThroughEnergy(int DestElemKey, int SourceElemKey, char* ForceComponID, int* SubdArray, long lenSubdivArray);
	void ComputeFieldTorqueThroughEnergy(int DestElemKey, int SourceElemKey, char* TorqueComponID, int* SubdArray, long lenSubdivArray, double* TorqueCenPo, long lenTorqueCenPo);
	//void OutFieldForceOrTorqueThroughEnergyCompRes(char* ForceComponID, TVector3d& Force, char ID);
	inline char CheckForAutoDestSubdivision(double* SubdivArray);

	void ComputeParticleTrajectory(int ElemKey, double E, double x0, double dxdy0, double z0, double dzdy0, double y0, double y1, int Np);
	void ComputeFocusPotent(int ElemKey, double* StPoi, long lenStPoi, double* FiPoi, long lenFiPoi, int Np);
	//void ComputeFocusKickPer(int ElemKey, double* P1, double* Nlong, double per, int nper, double* N1, double r1, int np1, double r2, int np2, const char* Comment, int nharm, int ns, double d1, double d2, const char* strKickUnit, double inEnergyGeV=0);
	void ComputeFocusKickPer(int ElemKey, double* P1, double* Nlong, double per, double nper, double* N1, double r1, int np1, double r2, int np2, const char* Comment, int nharm, int ns, double d1, double d2, const char* strKickUnit, double inEnergyGeV=0, const char* strOutFormat=0);
	void ComposeFocusKickPerFormStrRep(double* pKickData1, double* pKickData2, double* pBtE2Int, double* pCoordDir1, double* pCoordDir2, int np1, int np2, double per, int nper, const char* Comment);
	void ComputeFocusKick(int ElemKey, double* P1, double* Nlong, double* ArrLongDist, int lenArrLongDist, int ns, double* Ntr1, double r1, int np1, double r2, int np2, const char* StrComment, double d1, double d2);
	void ComputeShimSignature(int ElemKey, char* FldID, double* V, double* StPoi, double* FiPoi, int Np, double* Vi);

	//void OutFieldCompRes(char* FieldChar, radTField* FieldPtr, double* ArgArray, int Np);
	void OutFieldCompRes(char* FieldChar, radTField* FieldArray, long Np, radTVectInputCell& VectInputCell);
	void OutFieldCompRes(char* FieldChar, radTField* FieldArray, long Np);
	void ParseAndSendOneFieldValue(radTField* tField, char* BufChar, int AmOfItem);

	//void OutFieldIntCompRes(char* FieldIntChar, radTField* FieldPtr);
	void OutFieldEnergyForceCompRes(char* ComponIDChar, radTField* FieldPtr);

	int SetCompPrecisions(const char** ValNames, double* Values, int ValCount);
	int SetCompCriterium(double InAbsPrecB, double InAbsPrecA, double InAbsPrecB_int, double InAbsPrecFrc, double InAbsPrecTrjCoord, double InAbsPrecTrjAngle);
	int SetMltplThresh(double* InMltplThresh); // Maybe to be removed later
	int SetTolForConvergence(double AbsRandMagnitude, double RelRandMagnitude, double ZeroRandMagnitude);
	int RandomizationOnOrOff(char* OnOrOff);

	int SetAndShowPhysUnits();

	//void DumpElem(int ElemKey);
	void DumpElem(int* arKeys, int nElem, const char* strFormat, bool arKeysAllocInMathLink=false);
	int DumpElemParse(const unsigned char *bstr, int bstrLen);
	void GenDump();
	int RetrieveElemKey(const radTg* gPtr);

	inline double ReturnVersionID();
	void ReturnInput(double Input, int NumTimes);
	int SetMemAllocMethForIntrctMatr(char* TotOrParts);

	int ApplyDrawAttrToElem_g3d(int ElemKey, double* RGB_col, long lenRGB_col, double LineThickness =-1.);
	int RemoveDrawAttrFromElem_g3d(int ElemKey);
	int GraphicsForElem_g3d(int ElemKey, int InShowSymmetryChilds, const char** arOptionNames =0, const char** arOptionValues =0, int numOptions =0);

	void GraphicsForAll_g3d(int InShowSymmetryChilds);

	int GoQuickDraw3D_Viewer(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount);
	int GoOpenGL_3D_Viewer(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount);
	int DecodeViewingOptions(const char** OptionNames, const char** OptionValues, int OptionCount, char* OptBits);
	void PrepareGeomPolygDataForViewing(radTVectGeomPolygon& GeomPolygons, double*& VertCoord, int& Nv, int*& VertInd, int*& PgLen, float*& PgColors, int& Npg);
	void DeallocateAuxPgnViewData(double** dArr, int** iArr1, int** iArr2, float** fArr);

	int SubdivideElement_g3d(int ElemKey, double* SubdivArray, long lenSubdivArray, char TypeExtraSpec, double* ExtraSpec, long lenExtraSpec, const char** OptionNames, const char** OptionValues, int OptionCount);
	int CutElement_g3d(int ElemKey, double* PointOnPlane, long lenPointOnPlane, double* PlaneNormal, long lenPlaneNormal, const char** OptionNames, const char** OptionValues, int OptionCount);
	int SubdivideElement_g3dByParPlanes(int ElemKey, double* SubdivArray, int AmOfSubdivDirections, const char* LabOrLocFrame);
	int DuplicateElement_g3d(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount);

	int CreateFromObj_g3dWithSym(int ElemKey);
	inline int DeleteElement(int ElemKey);
	int DeleteAllElements(int DeletionMethNo);

	void ReplaceInAllGroups(radThg& OldHandle, radThg& NewHandle);
	void ReplaceInGlobalMap(radThg& OldHandle, radThg& NewHandle);
	inline void CopyDrawAttr(int OldElemKey, int NewElemKey);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

inline double radTApplication::ReturnVersionID()
{
	//double VersionID = 4.29; // Modified May 12, 2009
	//double VersionID = 4.30; // Modified June 24, 2012
	//double VersionID = 4.31; // Modified July 07, 2013
	double VersionID = 4.32; // Modified March 04, 2017

	if(SendingIsRequired) Send.Double(VersionID);
	return VersionID;
}

//-------------------------------------------------------------------------

inline int radTApplication::AddElementToContainer(radThg& hg)
{
	//CheckMemoryAvailable();
	GlobalMapOfHandlers[GlobalUniqueMapKey++] = hg;
	return GlobalUniqueMapKey - 1;
}

//-------------------------------------------------------------------------

inline int radTApplication::DeleteElement(int ElemKey)
{
	radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
	if(iter == GlobalMapOfHandlers.end())
	{
		if(SendingIsRequired) Send.ErrorMessage("Radia::Error002"); 
		return 0;
	}
	GlobalMapOfHandlers.erase(iter);

	radTMapOfDrawAttr::iterator iterDrawAttr = MapOfDrawAttr.find(ElemKey);
	if(iterDrawAttr != MapOfDrawAttr.end()) MapOfDrawAttr.erase(iterDrawAttr);

	if(SendingIsRequired) Send.Int(0);
	return 1;
}

//-------------------------------------------------------------------------

inline void radTApplication::CopyDrawAttr(int OldElemKey, int NewElemKey)
{
	radTMapOfDrawAttr::const_iterator iter = MapOfDrawAttr.find(OldElemKey);
	if(iter != MapOfDrawAttr.end()) MapOfDrawAttr[NewElemKey] = (*iter).second;
}

//-------------------------------------------------------------------------

inline char radTApplication::CheckForAutoDestSubdivision(double* SubdivArray)
{
	const double RelTol = 1.E-09;
	if((fabs(SubdivArray[0])<RelTol) && (fabs(SubdivArray[2])<RelTol) && (fabs(SubdivArray[4])<RelTol)) return 1;
	else return 0;
}

//-------------------------------------------------------------------------

inline int radTApplication::CheckIfExtrudedPolygonIsRecMag(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d)
{
	if((ArrayOfPoints2d == 0) || (lenArrayOfPoints2d != 4)) return 0;

	TVector2d &P2d_0 = ArrayOfPoints2d[0], &P2d_1 = ArrayOfPoints2d[1], &P2d_2 = ArrayOfPoints2d[2], &P2d_3 = ArrayOfPoints2d[3];
	double x0 = P2d_0.x, y0 = P2d_0.y;
	double x1 = P2d_1.x, y1 = P2d_1.y;
	double x2 = P2d_2.x, y2 = P2d_2.y;
	double x3 = P2d_3.x, y3 = P2d_3.y;
	
	double x2_mi_x0 = x2 - x0, y2_mi_y0 = y2 - y0;
	double charactLen1 = sqrt(x2_mi_x0*x2_mi_x0 + y2_mi_y0*y2_mi_y0);
	double x3_mi_x1 = x3 - x1, y3_mi_y1 = y3 - y1;
	double charactLen2 = sqrt(x3_mi_x1*x3_mi_x1 + y3_mi_y1*y3_mi_y1);
	double charactLen = (charactLen1 > charactLen2)? charactLen1 : charactLen2;
	double absTol = 2*radCR.AbsRandMagnitude(charactLen); //*2 precaution
	if(fabs(charactLen2 - charactLen1) > absTol) return 0;

	double abs_x10 = fabs(x1 - x0), abs_y10 = fabs(y1 - y0);
	double abs_x21 = fabs(x2 - x1), abs_y21 = fabs(y2 - y1);
	double abs_x32 = fabs(x3 - x2), abs_y32 = fabs(y3 - y2);
	double abs_x03 = fabs(x0 - x3), abs_y03 = fabs(y0 - y3);

	if(((abs_x10 < absTol) && (abs_y10 > absTol) &&
	    (abs_x21 > absTol) && (abs_y21 < absTol) &&
	    (abs_x32 < absTol) && (abs_y32 > absTol) &&
	    (abs_x03 > absTol) && (abs_y03 < absTol)) ||
	   ((abs_x10 > absTol) && (abs_y10 < absTol) &&
	    (abs_x21 < absTol) && (abs_y21 > absTol) &&
		(abs_x32 > absTol) && (abs_y32 < absTol) &&
		(abs_x03 < absTol) && (abs_y03 > absTol))) return 1;

	return 0;
}

//-------------------------------------------------------------------------

#endif
