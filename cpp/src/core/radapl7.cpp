/*-------------------------------------------------------------------------
*
* File name:      radapl7.cpp
*
* Project:        RADIA
*
* Description:    RADIA function calls
*
* Author(s):      Oleg Chubar
*
* First release:  2007
* 
* Copyright (C):  2007 by Synchrotron SOLEIL & European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radappl.h"
#include "radopnam.h"
#include "auxparse.h"

#define REAL double
extern "C" {
	#include "triangle.h"
	void triangulate(char *, struct triangulateio *, struct triangulateio *, struct triangulateio *);
	int gErrorTRIANGLE;
}
//extern void triangulate(char *, struct triangulateio *, struct triangulateio *, struct triangulateio *);

//-------------------------------------------------------------------------

int radTApplication::SetMultGenExtrTriangle(double* FirstPoi, long lenFirstPoi, double Lx, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* arSubdData, double* Magn, long lenMagn, const char* OrientStr, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	//SetMultGenExtrTriangle(FirstPoi, 3, radCR.Double(Lx), ArrayOfPoints2d, lenArrayOfPoints2d, AuxSubdDataPtr, Magn, 3, OrientStr, OptionNames, OptionValues, OptionCount);
	//Data structures for TRIANGLE
	struct triangulateio trIn, trOut, trVorOut;
	trIn.pointlist = NULL;
	trIn.pointattributelist = NULL;
	trIn.pointmarkerlist = NULL;
	trIn.segmentlist = NULL;
	trIn.segmentmarkerlist = NULL;
	//trIn.regionlist = NULL;
	
	trOut.pointlist = (REAL*) NULL; //Not needed if -N switch used.
	trOut.pointattributelist = (REAL*) NULL; //Not needed if -N switch used or number of point attributes is zero:
	trOut.pointmarkerlist = (int*) NULL; //Not needed if -N or -B switch used.
	trOut.trianglelist = (int*) NULL; //Not needed if -E switch used.
	trOut.triangleattributelist = (REAL*) NULL; //Not needed if -E switch used or number of triangle attributes is zero:
	trOut.neighborlist = (int*) NULL; //Needed only if -n switch used.
	trOut.segmentlist = (int*) NULL; //Needed only if segments are output (-p or -c) and -P not used:
	trOut.segmentmarkerlist = (int*) NULL; //Needed only if segments are output (-p or -c) and -P and -B not used:
	trOut.edgelist = (int*) NULL; //Needed only if -e switch used.
	trOut.edgemarkerlist = (int*) NULL; //Needed if -e used and -B not used.

	trVorOut.pointlist = (REAL*) NULL; //Needed only if -v switch used.
	trVorOut.pointattributelist = (REAL*) NULL; //Needed only if -v switch used and number of attributes is not zero:
	trVorOut.edgelist = (int*) NULL; //Needed only if -v switch used.
	trVorOut.normlist = (REAL*) NULL; //Needed only if -v switch used.

	short PrevSendingIsRequired = SendingIsRequired; 
	try
	{
		double CPoi[] = {0, 0, 0};

		SendingIsRequired = 0; //OC301207
		int IndTr = FindSpaceTransToOrientObjAlongMainAxis(CPoi, 'X', *OrientStr);
		if((Magn != 0) && (lenMagn >= 3) && (IndTr != 0)) TransformBackMagnOrCurDensArr(IndTr, Magn, lenMagn);
		
		TransformBackPointArr(IndTr, FirstPoi, lenFirstPoi); //OC090106
		
		TVector3d FirstPoiVect, MagnVect; //OC090106
		if(!ValidateVector3d(FirstPoi, lenFirstPoi, &FirstPoiVect)) throw 0;
		if(!ValidateVector3d(Magn, lenMagn, &MagnVect)) throw 0;

		//ki->Numb|Size,TriAngMin->...,TriAreaMax->...,ExtOpt->\"...\"
		char SubdParamBorderCode=0; // 0- ki are subdiv. numbers; 1- ki are average sizes of pieces;
		double triAngMin=0, triAreaMax=0;
		const char *sTriExtOpt=0;

		radTOptionNames OptNames;
		const char** BufNameString = OptionNames;
		const char** BufValString = OptionValues;

		for(int i=0; i<OptionCount; i++)
		{
			if(!strcmp(*BufNameString, OptNames.SubdParamBorderCode))
			{
				if(!strcmp(*BufValString, (OptNames.SubdParamCodeValues)[0])) SubdParamBorderCode = 0;
				else if(!strcmp(*BufValString, (OptNames.SubdParamCodeValues)[1])) SubdParamBorderCode = 1;
				else { Send.ErrorMessage("Radia::Error062"); throw 0;}
			}
			else if(!strcmp(*BufNameString, OptNames.TriAngMin))
			{
				triAngMin = atof(*BufValString);
			}
			else if(!strcmp(*BufNameString, OptNames.TriAreaMax))
			{
				triAreaMax = atof(*BufValString);
			}
			else if(!strcmp(*BufNameString, OptNames.TriExtOpt))
			{
				sTriExtOpt = *BufValString;
			}
			else { Send.ErrorMessage("Radia::Error062"); throw 0;}

			BufNameString++; BufValString++;
		}

		if(triAngMin > 35.) { Send.ErrorMessage("Radia::Error098"); throw 0;} //keep the constraint coherent with the error message!

	//Setting up final 2D points array taking into account subdivision
		const double AbsZeroTol = 5.E-12;
		int lenArrayOfPoints2d_mi_1 = lenArrayOfPoints2d - 1;

		vector<TVector2d> vResBorderPoints;
		TVector2d *tP1 = ArrayOfPoints2d, *tP2 = ArrayOfPoints2d + 1;

		double *t_arSubdData = arSubdData;
		for(int j=0; j<lenArrayOfPoints2d; j++)
		{
			if(j == lenArrayOfPoints2d_mi_1) tP2 = ArrayOfPoints2d;

			TVector2d P1P2 = (*tP2) - (*tP1);
			double absP1P2 = sqrt(P1P2*P1P2);
			double invAbsP1P2 = 0.;
			if(absP1P2 != 0) invAbsP1P2 = 1./absP1P2;

			TVector2d vE = invAbsP1P2*P1P2;

			double dSubdPar = *(t_arSubdData++);
			double q = *(t_arSubdData++);
			if((SubdParamBorderCode == 1) && (dSubdPar != 0.)) //ki are average sizes of pieces
			{
				dSubdPar = absP1P2/dSubdPar;
			}
			int numSegmParts = radTg3d::Round(dSubdPar);
			if(numSegmParts <= 0) numSegmParts = 1;

			double q0 = (fabs(numSegmParts - 1.) > AbsZeroTol)? pow(q, 1./(numSegmParts - 1.)) : q;
			double Buf = q*q0 - 1.;
			double a = (fabs(Buf) > AbsZeroTol)? absP1P2*(q0 - 1.)/Buf : absP1P2/numSegmParts;

			TVector2d SubP = *tP1;
			vResBorderPoints.push_back(SubP);
			for(int k=1; k<numSegmParts; k++)
			{
				SubP += a*vE;
				vResBorderPoints.push_back(SubP);
				a *= q0;
			}
			tP1++; tP2++;
		}

	//Setting up data for TRIANGLE
		int numBorderPoints = (int)vResBorderPoints.size();
		int numBorderPoints_mi_1 = numBorderPoints - 1;

		trIn.numberofpoints = numBorderPoints;
		trIn.numberofpointattributes = 1;
		//trIn.pointlist = (REAL*) malloc(trIn.numberofpoints * 2 * sizeof(REAL));
		trIn.pointlist = new REAL[numBorderPoints << 1];
		//trIn.pointattributelist = (REAL*) malloc(trIn.numberofpoints*trIn.numberofpointattributes*sizeof(REAL));
		trIn.pointattributelist = new REAL[numBorderPoints*trIn.numberofpointattributes];
		//trIn.pointmarkerlist = (int*) malloc(trIn.numberofpoints * sizeof(int));
		trIn.pointmarkerlist = new int[numBorderPoints];

		trIn.numberofsegments = numBorderPoints;
		trIn.segmentlist = new int[numBorderPoints << 1];
		trIn.segmentmarkerlist = new int[numBorderPoints];

		//trIn.numberofregions = 1;
		//trIn.regionlist = new REAL[trIn.numberofregions << 2];
		//trIn.regionlist[0] = 2.5;
		//trIn.regionlist[1] = 2.5;
		//trIn.regionlist[2] = 7.0; //Regional attribute (for whole mesh).
		//trIn.regionlist[3] = 0.1; //Area constraint that will not be used.

		//to modify:
		//TVector2d *tArrayOfPoints2d = ArrayOfPoints2d;
		REAL *t_pointlist = trIn.pointlist, *t_pointattributelist = trIn.pointattributelist;
		int *t_pointmarkerlist = trIn.pointmarkerlist;
		int *t_segmentlist = trIn.segmentlist, *t_segmentmarkerlist = trIn.segmentmarkerlist;
		for(int j=0; j<numBorderPoints; j++)
		{
			TVector2d &curPoint = vResBorderPoints[j];
			*(t_pointlist++) = curPoint.x;
			*(t_pointlist++) = curPoint.y;
			*(t_pointattributelist++) = 0.;
			*(t_pointmarkerlist++) = 0;

			*(t_segmentmarkerlist++) = 0;
			*t_segmentlist = j + 1;
			*(t_segmentlist + 1) = (j == numBorderPoints_mi_1)? 1 : j + 2;
			t_segmentlist += 2;
		}

		//trIn.numberofsegments = 0;
		//trIn.segmentlist = NULL;
		//trIn.segmentmarkerlist = NULL;

		trIn.trianglelist = NULL;
		trIn.triangleattributelist = NULL;
		trIn.trianglearealist = NULL;
		trIn.neighborlist = NULL;
		trIn.numberoftriangles = 0;
		trIn.numberofcorners = 0;
		trIn.numberoftriangleattributes = 0;

		trIn.holelist = NULL;
		trIn.numberofholes = 0;

		trIn.regionlist = NULL;
		trIn.numberofregions = 0;

		ostringstream osTriSwitches;
		osTriSwitches << "Qpq";
		if(triAngMin > 0.)
		{
			osTriSwitches << triAngMin;
		}
		if(triAreaMax > 0.)
		{
			osTriSwitches << "a" << triAreaMax;
		}
		if(sTriExtOpt != 0)
		{
			if(strlen(sTriExtOpt) > 0)
			{
				char sBufTriExtOpt[1024];
				//CAuxParse::StringSymbolsRemove(sTriExtOpt, "-", sBufTriExtOpt);
				CAuxParse::StringSymbolsRemove(sTriExtOpt, (char*)"-", sBufTriExtOpt); //OC01052013
				osTriSwitches << sBufTriExtOpt;
			}
		}
		osTriSwitches << ends;

		string sTriSwitches = osTriSwitches.str();
		const char *strTriSwitches = sTriSwitches.c_str();
		int lenStrTriSwitches = (int)strlen(strTriSwitches);
		if(lenStrTriSwitches > 255) lenStrTriSwitches = 255;
		char trSwitches[256];
		strncpy(trSwitches, strTriSwitches, lenStrTriSwitches);
		trSwitches[lenStrTriSwitches] = '\0';

	//Triangulate:
		gErrorTRIANGLE = 0;
		triangulate(trSwitches, &trIn, &trOut, &trVorOut);

		int numTri = trOut.numberoftriangles;
		int numPts = trOut.numberofpoints;
		if((gErrorTRIANGLE != 0) || (trOut.trianglelist == NULL) || (trOut.pointlist == NULL) || (numTri <= 0) || (numPts <= 0))
		{
			Send.ErrorMessage("Radia::Error119"); throw 0;
		}

		radTGroup* pGroup = new radTGroup();
		if(pGroup == 0) { Send.ErrorMessage("Radia::Error900"); throw 0;}

		REAL *arVertexCoord = trOut.pointlist;
		int *t_trianglelist = trOut.trianglelist;
		TVector2d arThreePts[3];
		TVector2d &trP1 = arThreePts[0], &trP2 = arThreePts[1], &trP3 = arThreePts[2];
		TAxisOrient AxOrnt = ParallelToX;
		for(int k=0; k<numTri; k++)
		{
			int indCoordTriVertex1 = (*(t_trianglelist++) - 1) << 1;
			int indCoordTriVertex2 = (*(t_trianglelist++) - 1) << 1;
			int indCoordTriVertex3 = (*(t_trianglelist++) - 1) << 1;

			trP1.x = arVertexCoord[indCoordTriVertex1];
			trP1.y = arVertexCoord[indCoordTriVertex1 + 1];
			trP2.x = arVertexCoord[indCoordTriVertex2];
			trP2.y = arVertexCoord[indCoordTriVertex2 + 1];
			trP3.x = arVertexCoord[indCoordTriVertex3];
			trP3.y = arVertexCoord[indCoordTriVertex3 + 1];

			FirstPoiVect.y = trP1.x;
			FirstPoiVect.z = trP1.y;

			radTExtrPolygon* pExtrPgn = new radTExtrPolygon(FirstPoiVect, AxOrnt, Lx, arThreePts, 3, MagnVect);
			if(pExtrPgn == 0) { Send.ErrorMessage("Radia::Error900"); throw 0;}
			if(((radTPolygon*)(pExtrPgn->BasePolygonHandle.rep))->SomethingIsWrong)
			{
				delete pExtrPgn; throw 0;
			}

			radThg hPgn(pExtrPgn);
			if(TreatExtrPgnsAsPolyhedrons) if(!pExtrPgn->ConvertToPolyhedron(hPgn, this, 1)) throw 0;

			int pgnKey = AddElementToContainer(hPgn);
            if(IndTr != 0) pgnKey = ApplySymmetry(pgnKey, IndTr, 1);
			pGroup->AddElement(pgnKey, hPgn);
		}

		radThg hGroup(pGroup);
		int resGrpKey = AddElementToContainer(hGroup);

		SendingIsRequired = PrevSendingIsRequired;
		if(SendingIsRequired) Send.Int(resGrpKey);

		if(trIn.pointlist != NULL) delete[] trIn.pointlist;
		if(trIn.pointattributelist != NULL) delete[] trIn.pointattributelist;
		if(trIn.pointmarkerlist != NULL) delete[] trIn.pointmarkerlist;
		if(trIn.segmentlist != NULL) delete[] trIn.segmentlist;
		if(trIn.segmentmarkerlist != NULL) delete[] trIn.segmentmarkerlist;
		//if(trIn.regionlist != NULL) delete[] trIn.regionlist;

		if(trOut.pointlist != NULL) free(trOut.pointlist); //Not needed if -N switch used.
		if(trOut.pointattributelist != NULL) free(trOut.pointattributelist); //Not needed if -N switch used or number of point attributes is zero:
		if(trOut.pointmarkerlist != NULL) free(trOut.pointmarkerlist); //Not needed if -N or -B switch used.
		if(trOut.trianglelist != NULL) free(trOut.trianglelist); //Not needed if -E switch used.
		if(trOut.triangleattributelist != NULL) free(trOut.triangleattributelist); //Not needed if -E switch used or number of triangle attributes is zero:
		if(trOut.neighborlist != NULL) free(trOut.neighborlist); //Needed only if -n switch used.
		if(trOut.segmentlist != NULL) free(trOut.segmentlist); //Needed only if segments are output (-p or -c) and -P not used:
		if(trOut.segmentmarkerlist != NULL) free(trOut.segmentmarkerlist); //Needed only if segments are output (-p or -c) and -P and -B not used:
		if(trOut.edgelist != NULL) free(trOut.edgelist); //Needed only if -e switch used.
		if(trOut.edgemarkerlist != NULL) free(trOut.edgemarkerlist); //Needed if -e used and -B not used.

		if(trVorOut.pointlist != NULL) free(trVorOut.pointlist); //Needed only if -v switch used.
		if(trVorOut.pointattributelist != NULL) free(trVorOut.pointattributelist); //Needed only if -v switch used and number of attributes is not zero:
		if(trVorOut.edgelist != NULL) free(trVorOut.edgelist); //Needed only if -v switch used.
		if(trVorOut.normlist != NULL) free(trVorOut.normlist); //Needed only if -v switch used.

		return resGrpKey;
	}
	catch(...)
	{
		SendingIsRequired = PrevSendingIsRequired;

		if(trIn.pointlist != NULL) delete[] trIn.pointlist;
		if(trIn.pointattributelist != NULL) delete[] trIn.pointattributelist;
		if(trIn.pointmarkerlist != NULL) delete[] trIn.pointmarkerlist;
		if(trIn.segmentlist != NULL) delete[] trIn.segmentlist;
		if(trIn.segmentmarkerlist != NULL) delete[] trIn.segmentmarkerlist;
		//if(trIn.regionlist != NULL) delete[] trIn.regionlist;

		if(trOut.pointlist != NULL) free(trOut.pointlist); //Not needed if -N switch used.
		if(trOut.pointattributelist != NULL) free(trOut.pointattributelist); //Not needed if -N switch used or number of point attributes is zero:
		if(trOut.pointmarkerlist != NULL) free(trOut.pointmarkerlist); //Not needed if -N or -B switch used.
		if(trOut.trianglelist != NULL) free(trOut.trianglelist); //Not needed if -E switch used.
		if(trOut.triangleattributelist != NULL) free(trOut.triangleattributelist); //Not needed if -E switch used or number of triangle attributes is zero:
		if(trOut.neighborlist != NULL) free(trOut.neighborlist); //Needed only if -n switch used.
		if(trOut.segmentlist != NULL) free(trOut.segmentlist); //Needed only if segments are output (-p or -c) and -P not used:
		if(trOut.segmentmarkerlist != NULL) free(trOut.segmentmarkerlist); //Needed only if segments are output (-p or -c) and -P and -B not used:
		if(trOut.edgelist != NULL) free(trOut.edgelist); //Needed only if -e switch used.
		if(trOut.edgemarkerlist != NULL) free(trOut.edgemarkerlist); //Needed if -e used and -B not used.

		if(trVorOut.pointlist != NULL) free(trVorOut.pointlist); //Needed only if -v switch used.
		if(trVorOut.pointattributelist != NULL) free(trVorOut.pointattributelist); //Needed only if -v switch used and number of attributes is not zero:
		if(trVorOut.edgelist != NULL) free(trVorOut.edgelist); //Needed only if -v switch used.
		if(trVorOut.normlist != NULL) free(trVorOut.normlist); //Needed only if -v switch used.

		Initialize(); 
		return 0;
	}
}

//-------------------------------------------------------------------------

//int radTApplication::TriangulatePolygon(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* arSubdData, const char** OptionNames, const char** OptionValues, int OptionCount, TVector2d*& arTriVertPt, int& numTriVertPt, int*& arTriVertInd, int& numTri)
int radTApplication::TriangulatePolygon(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* arSubdData, char triSubdParamBorderCode, double triAngMin, double triAreaMax, const char* sTriExtOpt, TVector2d*& arTriVertPt, int& numTriVertPt, int*& arTriVertInd, int& numTri)
{
	//Data structures for TRIANGLE
	struct triangulateio trIn, trOut, trVorOut;
	trIn.pointlist = NULL;
	trIn.pointattributelist = NULL;
	trIn.pointmarkerlist = NULL;
	trIn.segmentlist = NULL;
	trIn.segmentmarkerlist = NULL;
	//trIn.regionlist = NULL;
	
	trOut.pointlist = (REAL*) NULL; //Not needed if -N switch used.
	trOut.pointattributelist = (REAL*) NULL; //Not needed if -N switch used or number of point attributes is zero:
	trOut.pointmarkerlist = (int*) NULL; //Not needed if -N or -B switch used.
	trOut.trianglelist = (int*) NULL; //Not needed if -E switch used.
	trOut.triangleattributelist = (REAL*) NULL; //Not needed if -E switch used or number of triangle attributes is zero:
	trOut.neighborlist = (int*) NULL; //Needed only if -n switch used.
	trOut.segmentlist = (int*) NULL; //Needed only if segments are output (-p or -c) and -P not used:
	trOut.segmentmarkerlist = (int*) NULL; //Needed only if segments are output (-p or -c) and -P and -B not used:
	trOut.edgelist = (int*) NULL; //Needed only if -e switch used.
	trOut.edgemarkerlist = (int*) NULL; //Needed if -e used and -B not used.

	trVorOut.pointlist = (REAL*) NULL; //Needed only if -v switch used.
	trVorOut.pointattributelist = (REAL*) NULL; //Needed only if -v switch used and number of attributes is not zero:
	trVorOut.edgelist = (int*) NULL; //Needed only if -v switch used.
	trVorOut.normlist = (REAL*) NULL; //Needed only if -v switch used.

	//short PrevSendingIsRequired = SendingIsRequired; 
	try
	{
		//ki->Numb|Size,TriAngMin->...,TriAreaMax->...,ExtOpt->\"...\"
		//char SubdParamBorderCode=0; // 0- ki are subdiv. numbers; 1- ki are average sizes of pieces;
		//double triAngMin=0, triAreaMax=0;
		//const char *sTriExtOpt=0;
		//radTOptionNames OptNames;
		//const char** BufNameString = OptionNames;
		//const char** BufValString = OptionValues;
		//for(int i=0; i<OptionCount; i++)
		//{
		//	if(!strcmp(*BufNameString, OptNames.SubdParamBorderCode))
		//	{
		//		if(!strcmp(*BufValString, (OptNames.SubdParamCodeValues)[0])) SubdParamBorderCode = 0;
		//		else if(!strcmp(*BufValString, (OptNames.SubdParamCodeValues)[1])) SubdParamBorderCode = 1;
		//		else { Send.ErrorMessage("Radia::Error062"); throw 0;}
		//	}
		//	else if(!strcmp(*BufNameString, OptNames.TriAngMin))
		//	{
		//		triAngMin = atof(*BufValString);
		//	}
		//	else if(!strcmp(*BufNameString, OptNames.TriAreaMax))
		//	{
		//		triAreaMax = atof(*BufValString);
		//	}
		//	else if(!strcmp(*BufNameString, OptNames.TriExtOpt))
		//	{
		//		sTriExtOpt = *BufValString;
		//	}
		//	else { Send.ErrorMessage("Radia::Error062"); throw 0;}
		//	BufNameString++; BufValString++;
		//}
		//if(triAngMin > 35.) { Send.ErrorMessage("Radia::Error098"); throw 0;} //keep the constraint coherent with the error message!

	//Setting up final 2D points array (on boundary) taking into account subdivision
		const double AbsZeroTol = 5.E-12;
		int lenArrayOfPoints2d_mi_1 = lenArrayOfPoints2d - 1;

		vector<TVector2d> vResBorderPoints;
		TVector2d *tP1 = ArrayOfPoints2d, *tP2 = ArrayOfPoints2d + 1;

		double *t_arSubdData = arSubdData;
		for(int j=0; j<lenArrayOfPoints2d; j++)
		{
			if(j == lenArrayOfPoints2d_mi_1) tP2 = ArrayOfPoints2d;

			TVector2d P1P2 = (*tP2) - (*tP1);
			double absP1P2 = sqrt(P1P2*P1P2);
			double invAbsP1P2 = 0.;
			if(absP1P2 != 0) invAbsP1P2 = 1./absP1P2;

			TVector2d vE = invAbsP1P2*P1P2;

			double dSubdPar = *(t_arSubdData++);
			double q = *(t_arSubdData++);
			if((triSubdParamBorderCode == 1) && (dSubdPar != 0.)) //ki are average sizes of pieces
			{
				dSubdPar = absP1P2/dSubdPar;
			}
			int numSegmParts = radTg3d::Round(dSubdPar);
			if(numSegmParts <= 0) numSegmParts = 1;

			double q0 = (fabs(numSegmParts - 1.) > AbsZeroTol)? pow(q, 1./(numSegmParts - 1.)) : q;
			double Buf = q*q0 - 1.;
			double a = (fabs(Buf) > AbsZeroTol)? absP1P2*(q0 - 1.)/Buf : absP1P2/numSegmParts;

			TVector2d SubP = *tP1;
			vResBorderPoints.push_back(SubP);
			for(int k=1; k<numSegmParts; k++)
			{
				SubP += a*vE;
				vResBorderPoints.push_back(SubP);
				a *= q0;
			}
			tP1++; tP2++;
		}

	//Setting up data for TRIANGLE
		int numBorderPoints = (int)vResBorderPoints.size();
		int numBorderPoints_mi_1 = numBorderPoints - 1;

		trIn.numberofpoints = numBorderPoints;
		trIn.numberofpointattributes = 1;
		//trIn.pointlist = (REAL*) malloc(trIn.numberofpoints * 2 * sizeof(REAL));
		trIn.pointlist = new REAL[numBorderPoints << 1];
		//trIn.pointattributelist = (REAL*) malloc(trIn.numberofpoints*trIn.numberofpointattributes*sizeof(REAL));
		trIn.pointattributelist = new REAL[numBorderPoints*trIn.numberofpointattributes];
		//trIn.pointmarkerlist = (int*) malloc(trIn.numberofpoints * sizeof(int));
		trIn.pointmarkerlist = new int[numBorderPoints];

		trIn.numberofsegments = numBorderPoints;
		trIn.segmentlist = new int[numBorderPoints << 1];
		trIn.segmentmarkerlist = new int[numBorderPoints];

		//trIn.numberofregions = 1;
		//trIn.regionlist = new REAL[trIn.numberofregions << 2];
		//trIn.regionlist[0] = 2.5;
		//trIn.regionlist[1] = 2.5;
		//trIn.regionlist[2] = 7.0; //Regional attribute (for whole mesh).
		//trIn.regionlist[3] = 0.1; //Area constraint that will not be used.

		//to modify:
		//TVector2d *tArrayOfPoints2d = ArrayOfPoints2d;
		REAL *t_pointlist = trIn.pointlist, *t_pointattributelist = trIn.pointattributelist;
		int *t_pointmarkerlist = trIn.pointmarkerlist;
		int *t_segmentlist = trIn.segmentlist, *t_segmentmarkerlist = trIn.segmentmarkerlist;
		for(int j=0; j<numBorderPoints; j++)
		{
			TVector2d &curPoint = vResBorderPoints[j];
			*(t_pointlist++) = curPoint.x;
			*(t_pointlist++) = curPoint.y;
			*(t_pointattributelist++) = 0.;
			*(t_pointmarkerlist++) = 0;

			*(t_segmentmarkerlist++) = 0;
			*t_segmentlist = j + 1;
			*(t_segmentlist + 1) = (j == numBorderPoints_mi_1)? 1 : j + 2;
			t_segmentlist += 2;
		}

		//trIn.numberofsegments = 0;
		//trIn.segmentlist = NULL;
		//trIn.segmentmarkerlist = NULL;

		trIn.trianglelist = NULL;
		trIn.triangleattributelist = NULL;
		trIn.trianglearealist = NULL;
		trIn.neighborlist = NULL;
		trIn.numberoftriangles = 0;
		trIn.numberofcorners = 0;
		trIn.numberoftriangleattributes = 0;

		trIn.holelist = NULL;
		trIn.numberofholes = 0;

		trIn.regionlist = NULL;
		trIn.numberofregions = 0;

		ostringstream osTriSwitches;
		osTriSwitches << "Qpq";
		if(triAngMin > 0.)
		{
			osTriSwitches << triAngMin;
		}
		if(triAreaMax > 0.)
		{
			osTriSwitches << "a" << triAreaMax;
		}
		if(sTriExtOpt != 0)
		{
			if(strlen(sTriExtOpt) > 0)
			{
				char sBufTriExtOpt[1024];
				//CAuxParse::StringSymbolsRemove(sTriExtOpt, "-", sBufTriExtOpt);
				CAuxParse::StringSymbolsRemove(sTriExtOpt, (char*)"-", sBufTriExtOpt); //OC01052013
				osTriSwitches << sBufTriExtOpt;
			}
		}
		osTriSwitches << ends;

		string sTriSwitches = osTriSwitches.str();
		const char *strTriSwitches = sTriSwitches.c_str();
		int lenStrTriSwitches = (int)strlen(strTriSwitches);
		if(lenStrTriSwitches > 255) lenStrTriSwitches = 255;
		char trSwitches[256];
		strncpy(trSwitches, strTriSwitches, lenStrTriSwitches);
		trSwitches[lenStrTriSwitches] = '\0';

	//Triangulate:
		gErrorTRIANGLE = 0;
		triangulate(trSwitches, &trIn, &trOut, &trVorOut);

		//int numTri = trOut.numberoftriangles;
		//int numPts = trOut.numberofpoints;
		numTri = trOut.numberoftriangles;
		numTriVertPt = trOut.numberofpoints;

		if((gErrorTRIANGLE != 0) || (trOut.trianglelist == NULL) || (trOut.pointlist == NULL) || (numTri <= 0) || (numTriVertPt <= 0))
		{
			Send.ErrorMessage("Radia::Error119"); throw 0;
		}

		arTriVertPt = new TVector2d[numTriVertPt];
		TVector2d *t_arTriVert = arTriVertPt;
		REAL *arVertexCoord = trOut.pointlist;
		REAL *t_arVertexCoord = arVertexCoord;
		for(int k=0; k<numTriVertPt; k++)
		{
			t_arTriVert->x = *(t_arVertexCoord++);
			(t_arTriVert++)->y = *(t_arVertexCoord++);
		}

		arTriVertInd = new int[numTri*3];
		int *t_arTriVertInd = arTriVertInd;
		int *t_trianglelist = trOut.trianglelist;
		for(int j=0; j<numTri*3; j++)
		{
			*(t_arTriVertInd++) = *(t_trianglelist++) - 1;

			//int indCoordTriVertex1 = (*(t_trianglelist++) - 1) << 1;
			//int indCoordTriVertex2 = (*(t_trianglelist++) - 1) << 1;
			//int indCoordTriVertex3 = (*(t_trianglelist++) - 1) << 1;
			//t_arTriVert->x = arVertexCoord[indCoordTriVertex1];
			//(t_arTriVert++)->y = arVertexCoord[indCoordTriVertex1 + 1];
			//t_arTriVert->x = arVertexCoord[indCoordTriVertex2];
			//(t_arTriVert++)->y = arVertexCoord[indCoordTriVertex2 + 1];
			//t_arTriVert->x = arVertexCoord[indCoordTriVertex3];
			//(t_arTriVert++)->y = arVertexCoord[indCoordTriVertex3 + 1];
		}

		//SendingIsRequired = PrevSendingIsRequired;
		//if(SendingIsRequired) Send.Int(resGrpKey);

		if(trIn.pointlist != NULL) delete[] trIn.pointlist;
		if(trIn.pointattributelist != NULL) delete[] trIn.pointattributelist;
		if(trIn.pointmarkerlist != NULL) delete[] trIn.pointmarkerlist;
		if(trIn.segmentlist != NULL) delete[] trIn.segmentlist;
		if(trIn.segmentmarkerlist != NULL) delete[] trIn.segmentmarkerlist;
		//if(trIn.regionlist != NULL) delete[] trIn.regionlist;

		if(trOut.pointlist != NULL) free(trOut.pointlist); //Not needed if -N switch used.
		if(trOut.pointattributelist != NULL) free(trOut.pointattributelist); //Not needed if -N switch used or number of point attributes is zero:
		if(trOut.pointmarkerlist != NULL) free(trOut.pointmarkerlist); //Not needed if -N or -B switch used.
		if(trOut.trianglelist != NULL) free(trOut.trianglelist); //Not needed if -E switch used.
		if(trOut.triangleattributelist != NULL) free(trOut.triangleattributelist); //Not needed if -E switch used or number of triangle attributes is zero:
		if(trOut.neighborlist != NULL) free(trOut.neighborlist); //Needed only if -n switch used.
		if(trOut.segmentlist != NULL) free(trOut.segmentlist); //Needed only if segments are output (-p or -c) and -P not used:
		if(trOut.segmentmarkerlist != NULL) free(trOut.segmentmarkerlist); //Needed only if segments are output (-p or -c) and -P and -B not used:
		if(trOut.edgelist != NULL) free(trOut.edgelist); //Needed only if -e switch used.
		if(trOut.edgemarkerlist != NULL) free(trOut.edgemarkerlist); //Needed if -e used and -B not used.

		if(trVorOut.pointlist != NULL) free(trVorOut.pointlist); //Needed only if -v switch used.
		if(trVorOut.pointattributelist != NULL) free(trVorOut.pointattributelist); //Needed only if -v switch used and number of attributes is not zero:
		if(trVorOut.edgelist != NULL) free(trVorOut.edgelist); //Needed only if -v switch used.
		if(trVorOut.normlist != NULL) free(trVorOut.normlist); //Needed only if -v switch used.

		return 1;
	}
	catch(...)
	{
		//SendingIsRequired = PrevSendingIsRequired;

		if(trIn.pointlist != NULL) delete[] trIn.pointlist;
		if(trIn.pointattributelist != NULL) delete[] trIn.pointattributelist;
		if(trIn.pointmarkerlist != NULL) delete[] trIn.pointmarkerlist;
		if(trIn.segmentlist != NULL) delete[] trIn.segmentlist;
		if(trIn.segmentmarkerlist != NULL) delete[] trIn.segmentmarkerlist;
		//if(trIn.regionlist != NULL) delete[] trIn.regionlist;

		if(trOut.pointlist != NULL) free(trOut.pointlist); //Not needed if -N switch used.
		if(trOut.pointattributelist != NULL) free(trOut.pointattributelist); //Not needed if -N switch used or number of point attributes is zero:
		if(trOut.pointmarkerlist != NULL) free(trOut.pointmarkerlist); //Not needed if -N or -B switch used.
		if(trOut.trianglelist != NULL) free(trOut.trianglelist); //Not needed if -E switch used.
		if(trOut.triangleattributelist != NULL) free(trOut.triangleattributelist); //Not needed if -E switch used or number of triangle attributes is zero:
		if(trOut.neighborlist != NULL) free(trOut.neighborlist); //Needed only if -n switch used.
		if(trOut.segmentlist != NULL) free(trOut.segmentlist); //Needed only if segments are output (-p or -c) and -P not used:
		if(trOut.segmentmarkerlist != NULL) free(trOut.segmentmarkerlist); //Needed only if segments are output (-p or -c) and -P and -B not used:
		if(trOut.edgelist != NULL) free(trOut.edgelist); //Needed only if -e switch used.
		if(trOut.edgemarkerlist != NULL) free(trOut.edgemarkerlist); //Needed if -e used and -B not used.

		if(trVorOut.pointlist != NULL) free(trVorOut.pointlist); //Needed only if -v switch used.
		if(trVorOut.pointattributelist != NULL) free(trVorOut.pointattributelist); //Needed only if -v switch used and number of attributes is not zero:
		if(trVorOut.edgelist != NULL) free(trVorOut.edgelist); //Needed only if -v switch used.
		if(trVorOut.normlist != NULL) free(trVorOut.normlist); //Needed only if -v switch used.

		Initialize(); 
		return 0;
	}
}

//-------------------------------------------------------------------------
