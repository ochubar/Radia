/*-------------------------------------------------------------------------
*
* File name:      radapl5.cpp
*
* Project:        RADIA
*
* Description:    Wrapping RADIA application function calls: QD3D
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifdef _WITH_QD3D

#include "radq3ld.h"

#ifdef _WINDOWS
extern "C" {
#include "resource.h" // Windows resource IDs

#ifndef __RADQ3VW_H
#include "radq3vw.h" // QD3D Viewer functionality
#endif

#include "QD3DShader.h"

int LoadQD3DViewerDLL();
int LoadQD3D_DLL();

BOOL QD3D_InitApplication(HINSTANCE);
BOOL QD3D_InitInstance(HINSTANCE, int);
void ErrorHandler(TQ3Error, TQ3Error, long);

// Function pointers in QD3DViewer DLL:
extern radTFQ3WinViewerDraw radQ3WinViewerDraw;    
extern radTFQ3WinViewerUseGroup radQ3WinViewerUseGroup;
// Function pointers in QD3D DLL:
extern radTFQ3Error_Register radQ3Error_Register;
extern radTFQ3Object_Dispose radQ3Object_Dispose;
extern radTFQ3Group_AddObject radQ3Group_AddObject;
extern radTFQ3PhongIllumination_New radQ3PhongIllumination_New;
extern radTFQ3Group_New radQ3Group_New;

//Test
	extern radTFQ3Polygon_New radQ3Polygon_New;
	extern radTFQ3AttributeSet_New radQ3AttributeSet_New;
	extern radTFQ3AttributeSet_Add radQ3AttributeSet_Add;
//End Test
}

extern HINSTANCE hinstCurrentRadia;
extern HINSTANCE hinstPreviousRadia;
extern LPSTR lpszCmdLineRadia;
extern int nCmdShowRadia;

extern HWND gHwnd; // hwnd of QD3D Viewer
#endif

//-------------------------------------------------------------------------

#ifdef _MAC_OS
extern "C" {
#include "QD3DViewer.h"
//#include "QD3DGroup.h"

extern Boolean SupportsQuickDraw3D(void);
extern Boolean SupportsQuickDraw3DViewer(void);
extern void radStartQD3D_ViewerOnMac();
extern void radRunQD3D_ViewerEventLoopOnMac();
}
#endif

//-------------------------------------------------------------------------

extern TQ3ViewerObject gViewer;

//-------------------------------------------------------------------------

#ifndef __RADAPPL_H
#include "radappl.h"
#endif
#ifndef __RADOPNAM_H
#include "radopnam.h"
#endif

//-------------------------------------------------------------------------

int radTApplication::GoQuickDraw3D_Viewer(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		/**
		radTOptionNames OptNam;
		const char** BufNameString = OptionNames;
		const char** BufValString = OptionValues;

		char DoShowLines = 1;
		char DoShowFaces = 1;
		char DoShowFrameAxes = 1;
		for(int i=0; i<OptionCount; i++)
		{
		if(!strcmp(*BufNameString, OptNam.ShowLines))
		{
		if(!strcmp(*BufValString, (OptNam.ShowLinesValues)[0])) DoShowLines = 0;
		else if(!strcmp(*BufValString, (OptNam.ShowLinesValues)[1])) DoShowLines = 1;
		else if(!strcmp(*BufValString, (OptNam.ShowLinesValues)[2])) DoShowLines = 0;
		else if(!strcmp(*BufValString, (OptNam.ShowLinesValues)[3])) DoShowLines = 1;
		else { Send.ErrorMessage("Radia::Error062"); return 0;}
		}
		if(!strcmp(*BufNameString, OptNam.ShowFaces))
		{
		if(!strcmp(*BufValString, (OptNam.ShowFacesValues)[0])) DoShowFaces = 0;
		else if(!strcmp(*BufValString, (OptNam.ShowFacesValues)[1])) DoShowFaces = 1;
		else if(!strcmp(*BufValString, (OptNam.ShowFacesValues)[2])) DoShowFaces = 0;
		else if(!strcmp(*BufValString, (OptNam.ShowFacesValues)[3])) DoShowFaces = 1;
		else { Send.ErrorMessage("Radia::Error062"); return 0;}
		}
		else if(!strcmp(*BufNameString, OptNam.ShowFrameAxes))
		{
		if(!strcmp(*BufValString, (OptNam.ShowFrameAxesValues)[0])) DoShowFrameAxes = 0;
		else if(!strcmp(*BufValString, (OptNam.ShowFrameAxesValues)[1])) DoShowFrameAxes = 1;
		else if(!strcmp(*BufValString, (OptNam.ShowFrameAxesValues)[2])) DoShowFrameAxes = 0;
		else if(!strcmp(*BufValString, (OptNam.ShowFrameAxesValues)[3])) DoShowFrameAxes = 1;
		else { Send.ErrorMessage("Radia::Error062"); return 0;}
		}
		else { Send.ErrorMessage("Radia::Error062"); return 0;}
		BufNameString++; BufValString++;
		}
		**/

		char OptBits[4];
		char& DoShowLines = OptBits[0];
		char& DoShowFaces = OptBits[1];
		char& DoShowFrameAxes = OptBits[2];
		char& DoShowSymChilds = OptBits[3];

		if(!DecodeViewingOptions(OptionNames, OptionValues, OptionCount, OptBits)) return 0;

		//short ShowSymmetryChilds = DoShowSymChilds;
		//radGraphPresOptions InGraphPresOptions(ShowSymmetryChilds);
		radGraphPresOptions InGraphPresOptions(DoShowSymChilds);
		radTg3dGraphPresent* g3dGraphPresentPtr = g3dPtr->CreateGraphPresent();

#ifdef _WINDOWS
		if(!QD3D_ViewerWasInitialized)
		{
			if(!LoadQD3D_DLL()) 
			{
				Send.ErrorMessage("Radia::Error113"); return 0;
			}
			if(!LoadQD3DViewerDLL()) { Send.ErrorMessage("Radia::Error113"); return 0;}
		}
#endif
#ifdef _MAC_OS
		if(!QD3D_ViewerWasInitialized)
		{
			if(!SupportsQuickDraw3D()) { Send.ErrorMessage("Radia::Error113"); return 0;}
			if(!SupportsQuickDraw3DViewer()) { Send.ErrorMessage("Radia::Error113"); return 0;}

			TQ3Status myStatus = Q3Initialize();
			if(myStatus == kQ3Failure) { Send.ErrorMessage("Radia::Error113"); return 0;}

			QD3D_ViewerWasInitialized = 1;
		}
#endif

		Send.QD3D_GroupToDraw = radQ3Group_New();
		if(Send.QD3D_GroupToDraw == 0) { Send.ErrorMessage("Radia::Error113"); return 0;}
		Send.InitLimits3D();
		if(!Send.InitSmallRotForQD3D()) return 0;

		radTg3dGraphPresent::Send = Send;

		char DrawFacilityInd = 1; // QD3D Draw facility index
		g3dGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;

		if(DoShowLines)
		{
			(radTg3dGraphPresent::Send).ShowLines = g3dGraphPresentPtr->ShowEdgeLines = 1;
		}
		else
		{
			(radTg3dGraphPresent::Send).ShowLines = g3dGraphPresentPtr->ShowEdgeLines = 0;
		}

		if(DoShowFaces)
		{
			(radTg3dGraphPresent::Send).ShowFaces = g3dGraphPresentPtr->ShowFaces = 1;
		}
		else
		{
			(radTg3dGraphPresent::Send).ShowFaces = g3dGraphPresentPtr->ShowFaces = 0;
		}

		g3dGraphPresentPtr->SetGraphPresOptions(InGraphPresOptions);
		g3dGraphPresentPtr->MapOfDrawAttrPtr = &MapOfDrawAttr;
		g3dGraphPresentPtr->RetrieveDrawAttr(ElemKey);

		g3dGraphPresentPtr->GenDraw();

		Send = radTg3dGraphPresent::Send;
		if(DoShowFrameAxes) 
		{
			//Send.ShowLines = 1;
			//Send.FrameLines(DrawFacilityInd);

			g3dGraphPresentPtr->DrawFrameLines();
		}
		Send.DelSmallRotForQD3D(); // No more sending polygons nor lines after this !

		if(Send.QD3D_GroupToDraw != 0)
		{
#ifdef _WINDOWS

			TQ3Status aStatus = kQ3Failure;
			LPSTR lpCmdLine = NULL;
			int QD3D_nCmdShow = SW_SHOW; /* | SW_SHOWNORMAL; */

			if(!QD3D_ViewerWasInitialized)
			{
				//if(!QD3D_InitApplication(hinstCurrentRadia)) return (FALSE);
				QD3D_InitApplication(hinstCurrentRadia);
				QD3D_ViewerWasInitialized = 1;
			}
			if(!QD3D_InitInstance(hinstCurrentRadia, QD3D_nCmdShow)) return (FALSE);
			radQ3Error_Register(ErrorHandler, 0);
			HACCEL hAccelTable = LoadAccelerators(hinstCurrentRadia, MAKEINTRESOURCE(IDR_GENERIC));

			radQ3WinViewerUseGroup(gViewer, Send.QD3D_GroupToDraw);
			radQ3WinViewerDraw(gViewer);

			MSG msg;
			while(GetMessage(&msg, NULL, 0, 0)) 
			{
				if(!TranslateAccelerator(gHwnd, hAccelTable, &msg)) 
				{
					TranslateMessage(&msg);
					DispatchMessage(&msg);
				}
			}
#endif
#ifdef _MAC_OS
			radStartQD3D_ViewerOnMac();
			if(gViewer != 0)
			{
				/* test
				Q3ViewerUseGroup(gViewer, Send.QD3D_GroupToDraw);
				Q3ViewerDraw(gViewer);
				*/
				radRunQD3D_ViewerEventLoopOnMac();
			}
#endif
		}
		TQ3Status DeletedOK = radQ3Object_Dispose(Send.QD3D_GroupToDraw);
		(radTg3dGraphPresent::Send).QD3D_GroupToDraw = Send.QD3D_GroupToDraw = 0;

		delete g3dGraphPresentPtr;

		if(SendingIsRequired) Send.Int(ElemKey);
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

int radTApplication::GoQuickDraw3D_Viewer(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	Send.ErrorMessage("Radia::Error114"); return 0;
	return ElemKey;
}

//-------------------------------------------------------------------------

#endif //_WITH_QD3D
