/*-------------------------------------------------------------------------
*
* File name:      radapl6.cpp
*
* Project:        RADIA
*
* Description:    Wrapping RADIA application function calls: GLUT-OpenGL
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifdef _WITH_GLUT

#include "radappl.h"

#ifndef __GLUTLINK_H
#ifdef ALPHA__DLL__
#undef ALPHA__DLL__
#define ALPHA_NONE
#define RESET_ALPHA__DLL__
#endif

#include "glutlink.h"

#ifdef RESET_ALPHA__DLL__
#define ALPHA__DLL__
#undef ALPHA_NONE
#undef RESET_ALPHA__DLL__
#endif
#endif

COpenGLViewer glViewer;

//-------------------------------------------------------------------------

int radTApplication::GoOpenGL_3D_Viewer(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		bool SendingWasAlreadyDone = false;

		char OptBits[4];
		char& DoShowLines = OptBits[0];
		char& DoShowFaces = OptBits[1];
		char& DoShowFrameAxes = OptBits[2];
		char& DoShowSymChilds = OptBits[3];
		if(!DecodeViewingOptions(OptionNames, OptionValues, OptionCount, OptBits)) return 0;

		//short ShowSymmetryChilds = 1;
		char ShowSymmetryChilds = DoShowSymChilds;

		radGraphPresOptions InGraphPresOptions(ShowSymmetryChilds);
		radTg3dGraphPresent* g3dGraphPresentPtr = g3dPtr->CreateGraphPresent();

		char DrawFacilityInd = 2; // OpenGL Draw facility index
		g3dGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;

		radTg3dGraphPresent::Send = Send;
		//radTSend& PreSend = radTg3dGraphPresent::Send;

		//PreSend.ShowLines = g3dGraphPresentPtr->ShowEdgeLines = DoShowLines;
		//PreSend.ShowFaces = g3dGraphPresentPtr->ShowFaces = DoShowFaces;
		//g3dGraphPresentPtr->SetGraphPresOptions(InGraphPresOptions, DoShowLines, DoShowFaces);
		g3dGraphPresentPtr->SetGraphPresOptionsExt(InGraphPresOptions, DoShowLines, DoShowFaces);
		g3dGraphPresentPtr->MapOfDrawAttrPtr = &MapOfDrawAttr;
		g3dGraphPresentPtr->RetrieveDrawAttr(ElemKey);

		//PreSend.GeomPolygons.reserve(100);
		//PreSend.GeomLines.reserve(100);
		//PreSend.InitLimits3D();

		g3dGraphPresentPtr->GenDraw();
		if(DoShowFrameAxes) g3dGraphPresentPtr->DrawFrameLines();

		//if(DoShowFrameAxes) 
		//{
		//	PreSend.ShowLines = 1;
		//	PreSend.FrameLines(DrawFacilityInd);
		//}

		int AmOfPolygons = (int)((radTg3dGraphPresent::Send).GeomPolygons.size());
		int AmOfLines = (int)((radTg3dGraphPresent::Send).GeomLines.size());

		if((AmOfPolygons > 0) || (AmOfLines > 0))
		{
			double *PgnVertCoord=0, *LineVertCoord=0;
			int PgNv=0, LineNv=0;
			int *PgVertInd=0, *LineVertInd=0;
			int *PgLen=0, *LineLen=0;
			float *PgColors=0, *LineColors=0;
			int Npg=0, Nln=0;

			if(AmOfPolygons > 0)
			{
				PrepareGeomPolygDataForViewing((radTg3dGraphPresent::Send).GeomPolygons, PgnVertCoord, PgNv, PgVertInd, PgLen, PgColors, Npg);
			}
			if(AmOfLines > 0)
			{
				PrepareGeomPolygDataForViewing((radTg3dGraphPresent::Send).GeomLines, LineVertCoord, LineNv, LineVertInd, LineLen, LineColors, Nln);
			}

			double* VertCoordArr[] = {PgnVertCoord, LineVertCoord};
			int NvArr[] = {PgNv, LineNv};
			int* VertIndArr[] = {PgVertInd, LineVertInd};
			int* LenArr[] = {PgLen, LineLen};
			float* ColorsArr[] = {PgColors, LineColors};
			int NpArr[] = {Npg, Nln};

			//double* VertCoordArr[] = {PgnVertCoord, 0};
			//int NvArr[] = {PgNv, 0};
			//int* VertIndArr[] = {PgVertInd, 0};
			//int* LenArr[] = {PgLen, 0};
			//float* ColorsArr[] = {PgColors, 0};
			//int NpArr[] = {Npg, 0};

			//test
			//char ErrorMesTitle[] = "SRW Debug";
			//char ErrorStr[100];
			//int j = sprintf(ErrorStr, "AmOfLines: %d", AmOfLines);
			//j += sprintf(ErrorStr + j, "          Nln: %d", Nln);
			//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
			//int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
			//end test

			char ViewerStartMode = 1; //multi-thread
				//TEST
				//char ViewerStartMode = 0; //single-thread 
				//ENDTEST

			if(glViewer.viewScene3D(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr, "Radia 3D Viewer", ViewerStartMode, 0))
			{
				Send.ErrorMessage("Radia::Error115");
				SendingWasAlreadyDone = true;
			}

			(radTg3dGraphPresent::Send).DeallocateGeomPolygonData();
			DeallocateAuxPgnViewData(VertCoordArr, VertIndArr, LenArr, ColorsArr);
		}

		delete g3dGraphPresentPtr;

		if(SendingIsRequired && (!SendingWasAlreadyDone)) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

#else

//-------------------------------------------------------------------------

#ifndef __RADAPPL_H
#include "radappl.h"
#endif

//-------------------------------------------------------------------------

int radTApplication::GoOpenGL_3D_Viewer(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	Send.ErrorMessage("Radia::Error114"); return 0;
	return ElemKey;
}

//-------------------------------------------------------------------------

#endif //_WITH_GLUT

