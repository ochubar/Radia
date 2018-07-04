/*
 * This file automatically produced by C:\SoftwareDevelopments\Radia_Dev\cpp\vc\..\..\ext_lib\mathlink\bin\mprep.exe from:
 *	C:\SoftwareDevelopments\Radia_Dev\cpp\vc\..\src\clients\mathlink\radia.tm
 * mprep Revision 16 Copyright (c) Wolfram Research, Inc. 1990-2009
 */

#define MPREP_REVISION 16


#include "mathlink.h"

int MLAbort = 0;
int MLDone  = 0;
long MLSpecialCharacter = '\0';
HANDLE MLInstance = (HANDLE)0;
HWND MLIconWindow = (HWND)0;

MLINK stdlink = 0;
MLEnvironment stdenv = 0;
#if MLINTERFACE >= 3
MLYieldFunctionObject stdyielder = (MLYieldFunctionObject)0;
MLMessageHandlerObject stdhandler = (MLMessageHandlerObject)0;
#else
MLYieldFunctionObject stdyielder = 0;
MLMessageHandlerObject stdhandler = 0;
#endif /* MLINTERFACE >= 3 */

#include <windows.h>

#if defined(__GNUC__)

#	ifdef TCHAR
#		undef TCHAR
#	endif
#	define TCHAR char

#	ifdef PTCHAR
#		undef PTCHAR
#	endif
#	define PTCHAR char *

#	ifdef __TEXT
#		undef __TEXT
#	endif
#	define __TEXT(arg) arg

#	ifdef _tcsrchr
#		undef _tcsrchr
#	endif
#	define _tcsrchr strrchr

#	ifdef _tcscat
#		undef _tcscat
#	endif
#	define _tcscat strcat

#	ifdef _tcsncpy
#		undef _tcsncpy
#	endif
#	define _tcsncpy _fstrncpy
#else
#	include <tchar.h>
#endif

#include <stdlib.h>
#include <string.h>
#if (WIN32_MATHLINK || WIN64_MATHLINK || __GNUC__) && !defined(_fstrncpy)
#       define _fstrncpy strncpy
#endif

#ifndef CALLBACK
#define CALLBACK FAR PASCAL
typedef LONG LRESULT;
typedef unsigned int UINT;
typedef WORD WPARAM;
typedef DWORD LPARAM;
#endif


LRESULT CALLBACK MLEXPORT
IconProcedure( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

LRESULT CALLBACK MLEXPORT
IconProcedure( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch( msg){
	case WM_CLOSE:
		MLDone = 1;
		MLAbort = 1;
		break;
	case WM_QUERYOPEN:
		return 0;
	}
	return DefWindowProc( hWnd, msg, wParam, lParam);
}


#ifdef _UNICODE
#define _APISTR(i) L ## #i
#else
#define _APISTR(i) #i
#endif

#define APISTR(i) _APISTR(i)

HWND MLInitializeIcon( HINSTANCE hInstance, int nCmdShow)
{
	TCHAR path_name[260];
	PTCHAR icon_name;

	WNDCLASS  wc;
	HMODULE hdll;

	MLInstance = hInstance;
	if( ! nCmdShow) return (HWND)0;

	hdll = GetModuleHandle( __TEXT("ml32i" APISTR(MLINTERFACE)));

	(void)GetModuleFileName( hInstance, path_name, 260);

	icon_name = _tcsrchr( path_name, '\\') + 1;
	*_tcsrchr( icon_name, '.') = '\0';


	wc.style = 0;
	wc.lpfnWndProc = IconProcedure;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = hInstance;

	if( hdll)
		wc.hIcon = LoadIcon( hdll, __TEXT("MLIcon"));
	else
		wc.hIcon = LoadIcon( NULL, IDI_APPLICATION);

	wc.hCursor = LoadCursor( NULL, IDC_ARROW);
	wc.hbrBackground = (HBRUSH)GetStockObject( WHITE_BRUSH);
	wc.lpszMenuName =  NULL;
	wc.lpszClassName = __TEXT("mprepIcon");
	(void)RegisterClass( &wc);

	MLIconWindow = CreateWindow( __TEXT("mprepIcon"), icon_name,
			WS_OVERLAPPEDWINDOW | WS_MINIMIZE, CW_USEDEFAULT,
			CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
			(HWND)0, (HMENU)0, hInstance, (void FAR*)0);

	if( MLIconWindow){
		ShowWindow( MLIconWindow, SW_MINIMIZE);
		UpdateWindow( MLIconWindow);
	}
	return MLIconWindow;
}


#if __BORLANDC__
#pragma argsused
#endif

#if MLINTERFACE >= 3
MLYDEFN( int, MLDefaultYielder, ( MLINK mlp, MLYieldParameters yp))
#else
MLYDEFN( devyield_result, MLDefaultYielder, ( MLINK mlp, MLYieldParameters yp))
#endif /* MLINTERFACE >= 3 */
{
	MSG msg;

#if !__BORLANDC__
	mlp = mlp; /* suppress unused warning */
	yp = yp; /* suppress unused warning */
#endif

	if( PeekMessage( &msg, (HWND)0, 0, 0, PM_REMOVE)){
		TranslateMessage( &msg);
		DispatchMessage( &msg);
	}
	return MLDone;
}


/********************************* end header *********************************/


# line 3 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//:Evaluate:      BeginPackage["Radia`"]


//-------------------------------------------------------------------------


# line 191 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 34 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//:Evaluate:      radObjDrwQD3D::usage = "radObjDrwQD3D[obj,EdgeLines->True|False,Faces->True|False,Axes->True|False] starts an application for viewing 3D geometry of the object obj in interactive mode. The feature is implemented using QuickDraw 3D graphics library from Apple Computer. The option EdgeLines->True|False (default EdgeLines->True) highlights the edge lines of objects; the option Faces->True|False (default Faces->True) shows faces of the objects; the option Axes->True|False (default Axes->True) shows the Cartesian frame axes."
# line 196 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 107 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//:Evaluate:      Begin["Radia`Private`"]


//-------------------------------------------------------------------------


void RecMag P(( double,double,double, double,double,double, double,double,double ));

# line 208 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 124 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ExtrudedPolygon P(( ));

# line 214 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 135 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ExtrudedPolygon2 P(( ));

# line 220 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 146 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//:Begin:
//:Function:      ExtrudedPolygon2
//:Pattern:       radObjThckPgnMag[xCoordin_, lxWidth_, ListOf2dPoints_, Orient_String:"x", Magnetiz_List:{0.,0.,0.}]
//:Arguments:     { N[xCoordin], N[lxWidth], N[ListOf2dPoints], Orient, N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void PlanarPolygon P(( ));

# line 235 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 166 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void Polyhedron1 P(( ));

# line 241 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 176 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//:Begin:
//:Function:      Polyhedron1
//:Pattern:       radObjPolyhdr[ListOfPoints_, ListOfListOfIndexes_, Magnetiz_List:{0,0,0}, OptPar1_:0]
//:Arguments:     { N[ListOfPoints], Round[ListOfListOfIndexes], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]], N[OptPar1_] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:

//:Begin:
//:Function:      Polyhedron1
//:Pattern:       radObjPolyhdrMag[ListOfPoints_, ListOfListOfIndexes_, Magnetiz_List:{0,0,0}]
//:Arguments:     { N[ListOfPoints], Round[ListOfListOfIndexes], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void RecMagsAsExtrPolygons P(( const char* ));

# line 264 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 204 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void RecMagsAsPolyhedrons P(( const char* ));

# line 270 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 215 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ExtPgnsAsPolyhedrons P(( const char* ));

# line 276 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 226 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void RecognizeRecMags P(( const char* ));

# line 282 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 237 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void MultGenExtrPolygon P(( ));

# line 288 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 248 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//:Begin:
//:Function:      MultGenExtrPolygon
//:Pattern:       radObjMltExtPgnMag[ListOfLayerPolygons_, Magnetiz_List:{0,0,0}]
//:Arguments:     { N[ListOfLayerPolygons], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void MultGenExtrPolygonCur P(( ));
//radObjMltExtPgnCur[z:0,a:"z",{{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...},I]

# line 304 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 269 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void MultGenExtrPolygonMag P(( ));
//radObjMltExtPgnMag[z:0,a:"z",{{x1,y1},{x2,y2},...},{{k1,q1},{k2,q2},...},{{R1,T1,H1},{R2,T2,H2},...},{{mx1,my1,mz1},{mx2,my2,mz2},...}]

# line 311 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 281 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void MultGenExtrRectangle P(( ));

# line 317 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 292 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//:Begin:
//:Function:      MultGenExtrRectangle
//:Pattern:       radObjMltExtRtgMag[ListOfLayerRectangles_, Magnetiz_List:{0,0,0}]
//:Arguments:     { N[ListOfLayerRectangles], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void MultGenExtrTriangle P(( ));

# line 332 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 312 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void ArcMag P(( double,double,double, double,double, double,double, double, int, const char*, double,double,double ));
//
//:Begin:
//:Function:      ArcMag
//:Pattern:       radObjArcMag[{xCoordin_,yCoordin_,zCoordin_}, {rminWidth_,rmaxWidth_}, {phiminWidth_,phimaxWidth_}, HeightWidth_, SectN_, Orient_String:"z", Mv_List:{0,0,0}]
//:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[rminWidth],N[rmaxWidth], N[phiminWidth],N[phimaxWidth], N[HeightWidth], Round[SectN], Orient, N[Mv[[1]]],N[Mv[[2]]],N[Mv[[3]]] }
//:ArgumentTypes: { Real,Real,Real, Real,Real, Real,Real, Real, Integer, String, Real,Real,Real }
//:ReturnType:    Manual
//:End:


void ArcPolygon P(( ));

# line 349 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 334 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void CylMag P(( double,double,double, double, double, int, double,double,double, const char* ));

//:Begin:
//:Function:      CylMag
//:Pattern:       radObjCylMag[{xCoordin_,yCoordin_,zCoordin_}, Rad_, HeightWidth_, SectN_, Mv_List:{0,0,0}, Orient_String:"z"]
//:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[Rad], N[HeightWidth], Round[SectN], N[Mv[[1]]],N[Mv[[2]]],N[Mv[[3]]], Orient }
//:ArgumentTypes: { Real,Real,Real, Real, Real, Integer, Real,Real,Real, String }
//:ReturnType:    Manual
//:End:


void CylMag P(( double,double,double, double, double, int, const char*, double,double,double ));

# line 366 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 356 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void RecCur P(( double,double,double, double,double,double, double,double,double ));

# line 372 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 367 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ArcCur P(( double,double,double, double,double, double,double, double, int, double, const char*, const char* ));

# line 378 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 378 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void RaceTrack P(( double,double,double, double,double, double,double, double, int, double, const char*, const char* ));

# line 384 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 389 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//:Begin:
//:Function:      RaceTrack
//:Pattern:       radObjRaceTrkCur[{xCoordin_,yCoordin_,zCoordin_}, {rminWidth_,rmaxWidth_}, {lxWidth_,lyWidth_}, HeightWidth_, SectN_, Jaz_, ManOrAuto_String:"man", Orient_String:"z"]
//:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[rminWidth],N[rmaxWidth], N[lxWidth],N[lyWidth], N[HeightWidth], Round[SectN], N[Jaz], ManOrAuto, Orient }
//:ArgumentTypes: { Real,Real,Real, Real,Real, Real,Real, Real, Integer, Real, String, String }
//:ReturnType:    Manual
//:End:


void FlmCur P(( ));

# line 399 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 409 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void Rectngl P(( double,double,double, double,double ));

# line 405 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 420 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void Group P(( int*, long ));

# line 411 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 431 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void AddToGroup P(( int, int*, long ));

# line 417 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 442 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void OutGroupSubObjectKeys P(( int ));

# line 423 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 453 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void BackgroundFieldSource P(( double,double,double ));

# line 429 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 464 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void SubdivideElementG3D P(( ));

# line 435 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 475 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void CutElementG3D P(( ));

# line 441 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 486 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void SubdivideElementG3DByParPlanes P(( ));

# line 447 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 497 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void DuplicateElementG3D P(( ));

# line 453 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 508 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void CreateFromG3DObjectWithSymmetries P(( int ));

# line 459 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 519 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void NumberOfDegOfFreedom P(( int ));

# line 465 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 530 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void MagnOfObj P(( int ));

# line 471 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 541 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ObjField P(( int, const char* ));

# line 477 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 552 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ScaleCurInObj P(( int,double ));

# line 483 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 563 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void GeometricalVolume P(( int ));

# line 489 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 582 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void GeometricalLimits P(( int ));

# line 495 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 593 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FldCmpMetForSubdRecMag P(( int, int, int ));

# line 501 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 604 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void SetLocMgnInSbdRecMag P(( ));

# line 507 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 615 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//-------------------------------------------------------------------------


void Translation P(( double,double,double ));

# line 516 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 629 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void Rotation P(( double,double,double, double,double,double, double ));

# line 522 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 640 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void PlaneSym P(( double,double,double, double,double,double ));

# line 528 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 651 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FieldInversion P(( ));

# line 534 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 662 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void TransformObject P(( int, int ));

# line 540 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 673 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ApplySymmetry P(( int, int, int ));

# line 546 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 684 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void CombineTransformLeft P(( int, int ));

# line 552 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 695 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void CombineTransformRight P(( int, int ));

# line 558 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 706 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//-------------------------------------------------------------------------


void LinearMaterial P(( double,double, double,double,double ));

# line 567 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 720 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void LinearMaterial2 P(( double,double, double ));

# line 573 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 731 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void NonlinearIsotropMaterial P(( double,double,double, double,double,double ));

# line 579 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 742 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void NonlinearIsotropMaterial2 P(( double,double, double,double, double,double ));

# line 585 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 753 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void NonlinearIsotropMaterial3 P(( ));

# line 591 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 764 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void NonlinearLaminatedMaterialML P(( ));

# line 597 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 775 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void NonlinearAnisotropMaterial P(( ));

# line 603 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 786 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ApplyMaterial P(( int, int ));

# line 609 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 797 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void MvsH P(( int, const char*, double,double,double ));

# line 615 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 808 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//-------------------------------------------------------------------------


void PreRelax P(( int, int ));

# line 624 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 822 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ShowInteractMatrix P(( int ));

# line 630 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 833 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ShowInteractVector P(( int, const char* ));

# line 636 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 844 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ManualRelax P(( int, int, int, double ));

# line 642 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 855 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void AutoRelax P(( int, double, int, int ));
void AutoRelax P(( ));

# line 649 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 867 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void UpdateSourcesForRelax P(( int ));

# line 655 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 878 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void SolveGen P(( int, double, int, int ));

# line 661 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 889 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void CompCriterium P(( double, double, double, double, double,double ));

# line 667 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 900 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void CompPrecision P(( ));

# line 673 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 911 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void MultipoleThresholds P(( double, double, double, double )); // Maybe to be removed later

# line 679 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 922 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void Field P(( int, const char*, double,double,double, double,double,double, int, const char*, double ));

//:Begin:
//:Function:      Field
//:Pattern:       radFld[ElemKey_, FieldChar_String, {xCoordin_,yCoordin_,zCoordin_}, x2Coordin_:10.^23,y2Coordin_:10.^23,z2Coordin_:10.^23, np_:1, ShowArgFlag_String:"noarg", strtarg_:0.]
//:Arguments:     { Round[ElemKey], FieldChar, N[xCoordin],N[yCoordin],N[zCoordin], N[x2Coordin],N[y2Coordin],N[z2Coordin], Round[np], ShowArgFlag, N[strtarg]}
//:ArgumentTypes: { Integer, String, Real,Real,Real, Real,Real,Real, Integer, String, Real }
//:ReturnType:    Manual
//:End:


# line 694 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 942 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ShimSignature P(( int, const char*, double,double,double, double,double,double, double,double,double, int, double,double,double ));

# line 700 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 953 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FieldArbitraryPointsStruct P(( ));

//:Begin:
//:Function:      FieldArbitraryPointsStruct
//:Pattern:       radFld[ElemKey_, FieldChar_String, PointsStructure_]
//:Arguments:     { Round[ElemKey], FieldChar, N[PointsStructure] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:

# line 714 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 972 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FieldInt P(( int, const char*, const char*, double,double,double, double,double,double ));

# line 720 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 983 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FieldForce P(( int, int ));

# line 726 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 994 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FieldEnergy P(( int, int, int,int,int ));

# line 732 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1005 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FieldForceThroughEnergy P(( int, int, const char*, int,int,int ));

# line 738 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1016 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FieldTorqueThroughEnergy P(( int, int, const char*, double,double,double, int,int,int ));

# line 744 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1027 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void ParticleTrajectory P(( int, double, double,double,double,double, double,double, int ));

# line 750 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1038 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void FocusingPotential P(( int, double,double,double, double,double,double, int ));

# line 756 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1049 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void FocusingKickPer P(( int, double,double,double, double,double,double, double,int, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double ));
void FocusingKickPer P(( int, double,double,double, double,double,double, double,double, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double, const char* ));

# line 763 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1061 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void FocusingKick P(( int, double,double,double, double,double,double, double*,long,int, double,double,double, double,int,double,int, const char*, double,double ));
void FocusingKickML P(( ));

# line 770 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1073 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void PhysicalUnits P(( ));

# line 776 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1084 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void OffsetForConvergence P(( char*, double,double,double, double,double,double ));

//:Begin:
//:Function:      OffsetForConvergence
//:Pattern:       radFldOfst[PointsOrDims_String, {axReg_,ayReg_,azReg_}, {axRnd_,ayRnd_,azRnd_}]
//:Arguments:     { PointsOrDims, N[axReg],N[ayReg],N[azReg], N[axRnd],N[ayRnd],N[azRnd] }
//:ArgumentTypes: { String, Real,Real,Real, Real,Real,Real }
//:ReturnType:    Manual
//:End:


void TolForConvergence P(( double, double, double ));

# line 793 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1106 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void RandomizationOnOrOff P(( const char* ));

# line 799 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1117 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//-------------------------------------------------------------------------


void ApplyDrawAttrToElem P(( int, double,double,double, double ));

# line 808 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1131 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void ApplyColorToElem P(( int, double,double,double ));

//:Begin:
//:Function:      ApplyColorToElem
//:Pattern:       radObjCol[ElemKey_, {r_,g_,b_}]
//:Arguments:     { Round[ElemKey], N[r],N[g],N[b] }
//:ArgumentTypes: { Integer, Real,Real,Real }
//:ReturnType:    Manual
//:End:


void RemoveDrawAttrFromElem P(( int ));

# line 825 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1153 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void GraphicsForElemWithSymChilds P(( int ));

//:Begin:
//:Function:      GraphicsForElemWithSymChilds
//:Pattern:       radObjDrw[ElemKey_]
//:Arguments:     { Round[ElemKey] }
//:ArgumentTypes: { Integer }
//:ReturnType:    Manual
//:End:


void GraphicsForElemWithSymChildsExt P(( ));

# line 842 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1175 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void GraphicsForAllWithSymChilds P(( ));

//:Begin:
//:Function:      GraphicsForAllWithSymChilds
//:Pattern:       radObjDrwAll
//:Arguments:     { }
//:ArgumentTypes: { }
//:ReturnType:    Manual
//:End:


void GraphicsForElemWithoutSymChilds P(( int ));

# line 859 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1197 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void GraphicsForAllWithoutSymChilds P(( ));

//:Begin:
//:Function:      GraphicsForAllWithoutSymChilds
//:Pattern:       radObjDrwAllWithoutTrfMlt
//:Arguments:     { }
//:ArgumentTypes: { }
//:ReturnType:    Manual
//:End:


//void QuickDraw3D_Viewer P(( ));

//:Begin:
//:Function:      QuickDraw3D_Viewer
//:Pattern:       radObjDrwQD3D[ElemKey_, OptPar1_:0, OptPar2_:0, OptPar3_:0]
//:Arguments:     { Round[ElemKey], OptPar1, OptPar2, OptPar3 }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void OpenGL_3D_Viewer P(( ));

# line 887 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1230 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//-------------------------------------------------------------------------


void DeleteElement P(( int ));

# line 896 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1244 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void DeleteAllElements1 P(( ));

# line 902 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1255 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void DeleteAllElements2 P(( ));

# line 908 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1266 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void DumpElem P(( int ));

//:Begin:
//:Function:      DumpElem
//:Pattern:       radUtiDmp[ElemKey_, OutFormat_String:"asc"]
//:Arguments:     { Round[ElemKey], OutFormat }
//:ArgumentTypes: { Integer, String }
//:ReturnType:    Manual
//:End:

void DumpElem P(( ));

# line 924 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1286 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void GenDump P(( ));

//:Begin:
//:Function:      GenDump
//:Pattern:       radUtiDmpAll
//:Arguments:     { }
//:ArgumentTypes: { }
//:ReturnType:    Manual
//:End:


void DumpElemParse P(( ));

# line 941 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1308 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void RadiaVersion P(());

# line 947 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1319 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//void OutCommandsInfo P(( ));

//:Begin:
//:Function:      OutCommandsInfo
//:Pattern:       radUtiInfo
//:Arguments:     { }
//:ArgumentTypes: { }
//:ReturnType:    Manual
//:End:


void ReturnInput P(( double, int ));

# line 964 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1341 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void MemAllocMethForIntrctMatr P(( const char* ));

# line 970 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1352 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//------------ P ELLEAUME -------------------------------------------------


void StartProf P(( int, int, int ));

# line 979 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1366 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void StopProf P(());

# line 985 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1377 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
void InterruptTime P(( double ));

# line 991 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


# line 1388 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia.tm"
//-------------------------------------------------------------------------


//:Evaluate:      End[]


//-------------------------------------------------------------------------


//:Evaluate:      EndPackage[]
# line 1005 "C:\\SoftwareDevelopments\\Radia_Dev\\cpp\\vc\\..\\src\\clients\\mathlink\\radia_tm_win.c"


void RecMag P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9));

#if MLPROTOTYPES
static int _tr0( MLINK mlp)
#else
static int _tr0(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLNewPacket(mlp) ) goto L9;

	RecMag(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9);

	res = 1;
L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr0 */


void ExtrudedPolygon P(( void));

#if MLPROTOTYPES
static int _tr1( MLINK mlp)
#else
static int _tr1(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	ExtrudedPolygon();

	res = 1;

	return res;
} /* _tr1 */


void ExtrudedPolygon2 P(( void));

#if MLPROTOTYPES
static int _tr2( MLINK mlp)
#else
static int _tr2(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	ExtrudedPolygon2();

	res = 1;

	return res;
} /* _tr2 */


void PlanarPolygon P(( void));

#if MLPROTOTYPES
static int _tr3( MLINK mlp)
#else
static int _tr3(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	PlanarPolygon();

	res = 1;

	return res;
} /* _tr3 */


void Polyhedron1 P(( void));

#if MLPROTOTYPES
static int _tr4( MLINK mlp)
#else
static int _tr4(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	Polyhedron1();

	res = 1;

	return res;
} /* _tr4 */


void RecMagsAsExtrPolygons P(( const char * _tp1));

#if MLPROTOTYPES
static int _tr5( MLINK mlp)
#else
static int _tr5(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	RecMagsAsExtrPolygons(_tp1);

	res = 1;
L1:	MLDisownString(mlp, _tp1);

L0:	return res;
} /* _tr5 */


void RecMagsAsPolyhedrons P(( const char * _tp1));

#if MLPROTOTYPES
static int _tr6( MLINK mlp)
#else
static int _tr6(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	RecMagsAsPolyhedrons(_tp1);

	res = 1;
L1:	MLDisownString(mlp, _tp1);

L0:	return res;
} /* _tr6 */


void ExtPgnsAsPolyhedrons P(( const char * _tp1));

#if MLPROTOTYPES
static int _tr7( MLINK mlp)
#else
static int _tr7(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	ExtPgnsAsPolyhedrons(_tp1);

	res = 1;
L1:	MLDisownString(mlp, _tp1);

L0:	return res;
} /* _tr7 */


void RecognizeRecMags P(( const char * _tp1));

#if MLPROTOTYPES
static int _tr8( MLINK mlp)
#else
static int _tr8(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	RecognizeRecMags(_tp1);

	res = 1;
L1:	MLDisownString(mlp, _tp1);

L0:	return res;
} /* _tr8 */


void MultGenExtrPolygon P(( void));

#if MLPROTOTYPES
static int _tr9( MLINK mlp)
#else
static int _tr9(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	MultGenExtrPolygon();

	res = 1;

	return res;
} /* _tr9 */


void MultGenExtrPolygonCur P(( void));

#if MLPROTOTYPES
static int _tr10( MLINK mlp)
#else
static int _tr10(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	MultGenExtrPolygonCur();

	res = 1;

	return res;
} /* _tr10 */


void MultGenExtrPolygonMag P(( void));

#if MLPROTOTYPES
static int _tr11( MLINK mlp)
#else
static int _tr11(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	MultGenExtrPolygonMag();

	res = 1;

	return res;
} /* _tr11 */


void MultGenExtrRectangle P(( void));

#if MLPROTOTYPES
static int _tr12( MLINK mlp)
#else
static int _tr12(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	MultGenExtrRectangle();

	res = 1;

	return res;
} /* _tr12 */


void MultGenExtrTriangle P(( void));

#if MLPROTOTYPES
static int _tr13( MLINK mlp)
#else
static int _tr13(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	MultGenExtrTriangle();

	res = 1;

	return res;
} /* _tr13 */


void ArcPolygon P(( void));

#if MLPROTOTYPES
static int _tr14( MLINK mlp)
#else
static int _tr14(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	ArcPolygon();

	res = 1;

	return res;
} /* _tr14 */


void CylMag P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, int _tp6, const char * _tp7, double _tp8, double _tp9, double _tp10));

#if MLPROTOTYPES
static int _tr15( MLINK mlp)
#else
static int _tr15(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	int _tp6;
	const char * _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLNewPacket(mlp) ) goto L10;

	CylMag(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10);

	res = 1;
L10: L9: L8: L7:	MLDisownString(mlp, _tp7);
L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr15 */


void RecCur P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9));

#if MLPROTOTYPES
static int _tr16( MLINK mlp)
#else
static int _tr16(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLNewPacket(mlp) ) goto L9;

	RecCur(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9);

	res = 1;
L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr16 */


void ArcCur P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, int _tp9, double _tp10, const char * _tp11, const char * _tp12));

#if MLPROTOTYPES
static int _tr17( MLINK mlp)
#else
static int _tr17(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	int _tp9;
	double _tp10;
	const char * _tp11;
	const char * _tp12;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetString( mlp, &_tp12) ) goto L11;
	if ( ! MLNewPacket(mlp) ) goto L12;

	ArcCur(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12);

	res = 1;
L12:	MLDisownString(mlp, _tp12);
L11:	MLDisownString(mlp, _tp11);
L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr17 */


void RaceTrack P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, int _tp9, double _tp10, const char * _tp11, const char * _tp12));

#if MLPROTOTYPES
static int _tr18( MLINK mlp)
#else
static int _tr18(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	int _tp9;
	double _tp10;
	const char * _tp11;
	const char * _tp12;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetString( mlp, &_tp12) ) goto L11;
	if ( ! MLNewPacket(mlp) ) goto L12;

	RaceTrack(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12);

	res = 1;
L12:	MLDisownString(mlp, _tp12);
L11:	MLDisownString(mlp, _tp11);
L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr18 */


void FlmCur P(( void));

#if MLPROTOTYPES
static int _tr19( MLINK mlp)
#else
static int _tr19(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	FlmCur();

	res = 1;

	return res;
} /* _tr19 */


void Rectngl P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr20( MLINK mlp)
#else
static int _tr20(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	Rectngl(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr20 */


void Group P(( int * _tp1, long _tpl1));

#if MLPROTOTYPES
static int _tr21( MLINK mlp)
#else
static int _tr21(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int * _tp1;
	long _tpl1;
	if ( ! MLGetIntegerList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	Group(_tp1, _tpl1);

	res = 1;
L1:	MLDisownIntegerList( mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr21 */


void AddToGroup P(( int _tp1, int * _tp2, long _tpl2));

#if MLPROTOTYPES
static int _tr22( MLINK mlp)
#else
static int _tr22(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int * _tp2;
	long _tpl2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetIntegerList( mlp, &_tp2, &_tpl2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	AddToGroup(_tp1, _tp2, _tpl2);

	res = 1;
L2:	MLDisownIntegerList( mlp, _tp2, _tpl2);
L1: 
L0:	return res;
} /* _tr22 */


void OutGroupSubObjectKeys P(( int _tp1));

#if MLPROTOTYPES
static int _tr23( MLINK mlp)
#else
static int _tr23(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	OutGroupSubObjectKeys(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr23 */


void BackgroundFieldSource P(( double _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr24( MLINK mlp)
#else
static int _tr24(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	BackgroundFieldSource(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr24 */


void SubdivideElementG3D P(( void));

#if MLPROTOTYPES
static int _tr25( MLINK mlp)
#else
static int _tr25(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	SubdivideElementG3D();

	res = 1;

	return res;
} /* _tr25 */


void CutElementG3D P(( void));

#if MLPROTOTYPES
static int _tr26( MLINK mlp)
#else
static int _tr26(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	CutElementG3D();

	res = 1;

	return res;
} /* _tr26 */


void SubdivideElementG3DByParPlanes P(( void));

#if MLPROTOTYPES
static int _tr27( MLINK mlp)
#else
static int _tr27(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	SubdivideElementG3DByParPlanes();

	res = 1;

	return res;
} /* _tr27 */


void DuplicateElementG3D P(( void));

#if MLPROTOTYPES
static int _tr28( MLINK mlp)
#else
static int _tr28(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	DuplicateElementG3D();

	res = 1;

	return res;
} /* _tr28 */


void CreateFromG3DObjectWithSymmetries P(( int _tp1));

#if MLPROTOTYPES
static int _tr29( MLINK mlp)
#else
static int _tr29(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	CreateFromG3DObjectWithSymmetries(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr29 */


void NumberOfDegOfFreedom P(( int _tp1));

#if MLPROTOTYPES
static int _tr30( MLINK mlp)
#else
static int _tr30(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	NumberOfDegOfFreedom(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr30 */


void MagnOfObj P(( int _tp1));

#if MLPROTOTYPES
static int _tr31( MLINK mlp)
#else
static int _tr31(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	MagnOfObj(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr31 */


void ObjField P(( int _tp1, const char * _tp2));

#if MLPROTOTYPES
static int _tr32( MLINK mlp)
#else
static int _tr32(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	const char * _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	ObjField(_tp1, _tp2);

	res = 1;
L2:	MLDisownString(mlp, _tp2);
L1: 
L0:	return res;
} /* _tr32 */


void ScaleCurInObj P(( int _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr33( MLINK mlp)
#else
static int _tr33(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	ScaleCurInObj(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr33 */


void GeometricalVolume P(( int _tp1));

#if MLPROTOTYPES
static int _tr34( MLINK mlp)
#else
static int _tr34(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	GeometricalVolume(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr34 */


void GeometricalVolume P(( int _tp1));

#if MLPROTOTYPES
static int _tr35( MLINK mlp)
#else
static int _tr35(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	GeometricalVolume(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr35 */


void GeometricalLimits P(( int _tp1));

#if MLPROTOTYPES
static int _tr36( MLINK mlp)
#else
static int _tr36(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	GeometricalLimits(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr36 */


void FldCmpMetForSubdRecMag P(( int _tp1, int _tp2, int _tp3));

#if MLPROTOTYPES
static int _tr37( MLINK mlp)
#else
static int _tr37(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	FldCmpMetForSubdRecMag(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr37 */


void SetLocMgnInSbdRecMag P(( void));

#if MLPROTOTYPES
static int _tr38( MLINK mlp)
#else
static int _tr38(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	SetLocMgnInSbdRecMag();

	res = 1;

	return res;
} /* _tr38 */


void Translation P(( double _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr39( MLINK mlp)
#else
static int _tr39(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	Translation(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr39 */


void Rotation P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7));

#if MLPROTOTYPES
static int _tr40( MLINK mlp)
#else
static int _tr40(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLNewPacket(mlp) ) goto L7;

	Rotation(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7);

	res = 1;
L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr40 */


void PlaneSym P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6));

#if MLPROTOTYPES
static int _tr41( MLINK mlp)
#else
static int _tr41(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLNewPacket(mlp) ) goto L6;

	PlaneSym(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6);

	res = 1;
L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr41 */


void FieldInversion P(( void));

#if MLPROTOTYPES
static int _tr42( MLINK mlp)
#else
static int _tr42(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	FieldInversion();

	res = 1;

L0:	return res;
} /* _tr42 */


void TransformObject P(( int _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr43( MLINK mlp)
#else
static int _tr43(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	TransformObject(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr43 */


void ApplySymmetry P(( int _tp1, int _tp2, int _tp3));

#if MLPROTOTYPES
static int _tr44( MLINK mlp)
#else
static int _tr44(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	ApplySymmetry(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr44 */


void CombineTransformLeft P(( int _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr45( MLINK mlp)
#else
static int _tr45(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	CombineTransformLeft(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr45 */


void CombineTransformRight P(( int _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr46( MLINK mlp)
#else
static int _tr46(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	CombineTransformRight(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr46 */


void LinearMaterial P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr47( MLINK mlp)
#else
static int _tr47(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	LinearMaterial(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr47 */


void LinearMaterial2 P(( double _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr48( MLINK mlp)
#else
static int _tr48(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	LinearMaterial2(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr48 */


void NonlinearIsotropMaterial P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6));

#if MLPROTOTYPES
static int _tr49( MLINK mlp)
#else
static int _tr49(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLNewPacket(mlp) ) goto L6;

	NonlinearIsotropMaterial(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6);

	res = 1;
L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr49 */


void NonlinearIsotropMaterial2 P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6));

#if MLPROTOTYPES
static int _tr50( MLINK mlp)
#else
static int _tr50(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLNewPacket(mlp) ) goto L6;

	NonlinearIsotropMaterial2(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6);

	res = 1;
L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr50 */


void NonlinearIsotropMaterial3 P(( void));

#if MLPROTOTYPES
static int _tr51( MLINK mlp)
#else
static int _tr51(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	NonlinearIsotropMaterial3();

	res = 1;

	return res;
} /* _tr51 */


void NonlinearLaminatedMaterialML P(( void));

#if MLPROTOTYPES
static int _tr52( MLINK mlp)
#else
static int _tr52(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	NonlinearLaminatedMaterialML();

	res = 1;

	return res;
} /* _tr52 */


void NonlinearAnisotropMaterial P(( void));

#if MLPROTOTYPES
static int _tr53( MLINK mlp)
#else
static int _tr53(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	NonlinearAnisotropMaterial();

	res = 1;

	return res;
} /* _tr53 */


void ApplyMaterial P(( int _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr54( MLINK mlp)
#else
static int _tr54(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	ApplyMaterial(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr54 */


void MvsH P(( int _tp1, const char * _tp2, double _tp3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr55( MLINK mlp)
#else
static int _tr55(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	const char * _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	MvsH(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
L5: L4: L3: L2:	MLDisownString(mlp, _tp2);
L1: 
L0:	return res;
} /* _tr55 */


void PreRelax P(( int _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr56( MLINK mlp)
#else
static int _tr56(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	PreRelax(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr56 */


void ShowInteractMatrix P(( int _tp1));

#if MLPROTOTYPES
static int _tr57( MLINK mlp)
#else
static int _tr57(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	ShowInteractMatrix(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr57 */


void ShowInteractVector P(( int _tp1, const char * _tp2));

#if MLPROTOTYPES
static int _tr58( MLINK mlp)
#else
static int _tr58(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	const char * _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	ShowInteractVector(_tp1, _tp2);

	res = 1;
L2:	MLDisownString(mlp, _tp2);
L1: 
L0:	return res;
} /* _tr58 */


void ManualRelax P(( int _tp1, int _tp2, int _tp3, double _tp4));

#if MLPROTOTYPES
static int _tr59( MLINK mlp)
#else
static int _tr59(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	double _tp4;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLNewPacket(mlp) ) goto L4;

	ManualRelax(_tp1, _tp2, _tp3, _tp4);

	res = 1;
L4: L3: L2: L1: 
L0:	return res;
} /* _tr59 */


void AutoRelax P(( void));

#if MLPROTOTYPES
static int _tr60( MLINK mlp)
#else
static int _tr60(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	AutoRelax();

	res = 1;

	return res;
} /* _tr60 */


void UpdateSourcesForRelax P(( int _tp1));

#if MLPROTOTYPES
static int _tr61( MLINK mlp)
#else
static int _tr61(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	UpdateSourcesForRelax(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr61 */


void SolveGen P(( int _tp1, double _tp2, int _tp3, int _tp4));

#if MLPROTOTYPES
static int _tr62( MLINK mlp)
#else
static int _tr62(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	int _tp3;
	int _tp4;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLNewPacket(mlp) ) goto L4;

	SolveGen(_tp1, _tp2, _tp3, _tp4);

	res = 1;
L4: L3: L2: L1: 
L0:	return res;
} /* _tr62 */


void CompCriterium P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6));

#if MLPROTOTYPES
static int _tr63( MLINK mlp)
#else
static int _tr63(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLNewPacket(mlp) ) goto L6;

	CompCriterium(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6);

	res = 1;
L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr63 */


void CompPrecision P(( void));

#if MLPROTOTYPES
static int _tr64( MLINK mlp)
#else
static int _tr64(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	CompPrecision();

	res = 1;

	return res;
} /* _tr64 */


void MultipoleThresholds P(( double _tp1, double _tp2, double _tp3, double _tp4));

#if MLPROTOTYPES
static int _tr65( MLINK mlp)
#else
static int _tr65(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLNewPacket(mlp) ) goto L4;

	MultipoleThresholds(_tp1, _tp2, _tp3, _tp4);

	res = 1;
L4: L3: L2: L1: 
L0:	return res;
} /* _tr65 */


void Field P(( int _tp1, const char * _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, int _tp9, const char * _tp10, double _tp11));

#if MLPROTOTYPES
static int _tr66( MLINK mlp)
#else
static int _tr66(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	const char * _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	int _tp9;
	const char * _tp10;
	double _tp11;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLNewPacket(mlp) ) goto L11;

	Field(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11);

	res = 1;
L11: L10:	MLDisownString(mlp, _tp10);
L9: L8: L7: L6: L5: L4: L3: L2:	MLDisownString(mlp, _tp2);
L1: 
L0:	return res;
} /* _tr66 */


void ShimSignature P(( int _tp1, const char * _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, int _tp12, double _tp13, double _tp14, double _tp15));

#if MLPROTOTYPES
static int _tr67( MLINK mlp)
#else
static int _tr67(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	const char * _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	int _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLNewPacket(mlp) ) goto L15;

	ShimSignature(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15);

	res = 1;
L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2:	MLDisownString(mlp, _tp2);
L1: 
L0:	return res;
} /* _tr67 */


void FieldArbitraryPointsStruct P(( int _tp1, const char * _tp2));

#if MLPROTOTYPES
static int _tr68( MLINK mlp)
#else
static int _tr68(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	const char * _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;

	FieldArbitraryPointsStruct(_tp1, _tp2);

	res = 1;
	MLDisownString(mlp, _tp2);
L1: 
L0:	return res;
} /* _tr68 */


void FieldInt P(( int _tp1, const char * _tp2, const char * _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9));

#if MLPROTOTYPES
static int _tr69( MLINK mlp)
#else
static int _tr69(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	const char * _tp2;
	const char * _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLNewPacket(mlp) ) goto L9;

	FieldInt(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9);

	res = 1;
L9: L8: L7: L6: L5: L4: L3:	MLDisownString(mlp, _tp3);
L2:	MLDisownString(mlp, _tp2);
L1: 
L0:	return res;
} /* _tr69 */


void FieldForce P(( int _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr70( MLINK mlp)
#else
static int _tr70(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	FieldForce(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr70 */


void FieldEnergy P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5));

#if MLPROTOTYPES
static int _tr71( MLINK mlp)
#else
static int _tr71(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	FieldEnergy(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr71 */


void FieldForceThroughEnergy P(( int _tp1, int _tp2, const char * _tp3, int _tp4, int _tp5, int _tp6));

#if MLPROTOTYPES
static int _tr72( MLINK mlp)
#else
static int _tr72(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	const char * _tp3;
	int _tp4;
	int _tp5;
	int _tp6;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLNewPacket(mlp) ) goto L6;

	FieldForceThroughEnergy(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6);

	res = 1;
L6: L5: L4: L3:	MLDisownString(mlp, _tp3);
L2: L1: 
L0:	return res;
} /* _tr72 */


void FieldTorqueThroughEnergy P(( int _tp1, int _tp2, const char * _tp3, double _tp4, double _tp5, double _tp6, int _tp7, int _tp8, int _tp9));

#if MLPROTOTYPES
static int _tr73( MLINK mlp)
#else
static int _tr73(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	const char * _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLNewPacket(mlp) ) goto L9;

	FieldTorqueThroughEnergy(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9);

	res = 1;
L9: L8: L7: L6: L5: L4: L3:	MLDisownString(mlp, _tp3);
L2: L1: 
L0:	return res;
} /* _tr73 */


void ParticleTrajectory P(( int _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, int _tp9));

#if MLPROTOTYPES
static int _tr74( MLINK mlp)
#else
static int _tr74(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	int _tp9;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLNewPacket(mlp) ) goto L9;

	ParticleTrajectory(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9);

	res = 1;
L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr74 */


void FocusingPotential P(( int _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, int _tp8));

#if MLPROTOTYPES
static int _tr75( MLINK mlp)
#else
static int _tr75(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	int _tp8;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLNewPacket(mlp) ) goto L8;

	FocusingPotential(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8);

	res = 1;
L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr75 */


void FocusingKickPer P(( int _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, int _tp14, double _tp15, int _tp16, const char * _tp17, int _tp18, int _tp19, double _tp20, double _tp21, const char * _tp22, double _tp23, const char * _tp24));

#if MLPROTOTYPES
static int _tr76( MLINK mlp)
#else
static int _tr76(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	double _tp15;
	int _tp16;
	const char * _tp17;
	int _tp18;
	int _tp19;
	double _tp20;
	double _tp21;
	const char * _tp22;
	double _tp23;
	const char * _tp24;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetString( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetString( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetString( mlp, &_tp24) ) goto L23;
	if ( ! MLNewPacket(mlp) ) goto L24;

	FocusingKickPer(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24);

	res = 1;
L24:	MLDisownString(mlp, _tp24);
L23: L22:	MLDisownString(mlp, _tp22);
L21: L20: L19: L18: L17:	MLDisownString(mlp, _tp17);
L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr76 */


void FocusingKickML P(( void));

#if MLPROTOTYPES
static int _tr77( MLINK mlp)
#else
static int _tr77(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	FocusingKickML();

	res = 1;

	return res;
} /* _tr77 */


void PhysicalUnits P(( void));

#if MLPROTOTYPES
static int _tr78( MLINK mlp)
#else
static int _tr78(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	PhysicalUnits();

	res = 1;

L0:	return res;
} /* _tr78 */


void TolForConvergence P(( double _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr79( MLINK mlp)
#else
static int _tr79(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	TolForConvergence(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr79 */


void RandomizationOnOrOff P(( const char * _tp1));

#if MLPROTOTYPES
static int _tr80( MLINK mlp)
#else
static int _tr80(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	RandomizationOnOrOff(_tp1);

	res = 1;
L1:	MLDisownString(mlp, _tp1);

L0:	return res;
} /* _tr80 */


void ApplyDrawAttrToElem P(( int _tp1, double _tp2, double _tp3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr81( MLINK mlp)
#else
static int _tr81(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	ApplyDrawAttrToElem(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr81 */


void RemoveDrawAttrFromElem P(( int _tp1));

#if MLPROTOTYPES
static int _tr82( MLINK mlp)
#else
static int _tr82(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	RemoveDrawAttrFromElem(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr82 */


void GraphicsForElemWithSymChildsExt P(( void));

#if MLPROTOTYPES
static int _tr83( MLINK mlp)
#else
static int _tr83(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	GraphicsForElemWithSymChildsExt();

	res = 1;

	return res;
} /* _tr83 */


void GraphicsForElemWithoutSymChilds P(( int _tp1));

#if MLPROTOTYPES
static int _tr84( MLINK mlp)
#else
static int _tr84(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	GraphicsForElemWithoutSymChilds(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr84 */


void OpenGL_3D_Viewer P(( void));

#if MLPROTOTYPES
static int _tr85( MLINK mlp)
#else
static int _tr85(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	OpenGL_3D_Viewer();

	res = 1;

	return res;
} /* _tr85 */


void DeleteElement P(( int _tp1));

#if MLPROTOTYPES
static int _tr86( MLINK mlp)
#else
static int _tr86(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	DeleteElement(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr86 */


void DeleteAllElements1 P(( void));

#if MLPROTOTYPES
static int _tr87( MLINK mlp)
#else
static int _tr87(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	DeleteAllElements1();

	res = 1;

L0:	return res;
} /* _tr87 */


void DeleteAllElements2 P(( void));

#if MLPROTOTYPES
static int _tr88( MLINK mlp)
#else
static int _tr88(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	DeleteAllElements2();

	res = 1;

L0:	return res;
} /* _tr88 */


void DumpElem P(( void));

#if MLPROTOTYPES
static int _tr89( MLINK mlp)
#else
static int _tr89(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	DumpElem();

	res = 1;

	return res;
} /* _tr89 */


void DumpElemParse P(( void));

#if MLPROTOTYPES
static int _tr90( MLINK mlp)
#else
static int _tr90(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	DumpElemParse();

	res = 1;

	return res;
} /* _tr90 */


void RadiaVersion P(( void));

#if MLPROTOTYPES
static int _tr91( MLINK mlp)
#else
static int _tr91(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	RadiaVersion();

	res = 1;

L0:	return res;
} /* _tr91 */


void ReturnInput P(( double _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr92( MLINK mlp)
#else
static int _tr92(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	int _tp2;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	ReturnInput(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr92 */


void MemAllocMethForIntrctMatr P(( const char * _tp1));

#if MLPROTOTYPES
static int _tr93( MLINK mlp)
#else
static int _tr93(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	MemAllocMethForIntrctMatr(_tp1);

	res = 1;
L1:	MLDisownString(mlp, _tp1);

L0:	return res;
} /* _tr93 */


void StartProf P(( int _tp1, int _tp2, int _tp3));

#if MLPROTOTYPES
static int _tr94( MLINK mlp)
#else
static int _tr94(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	StartProf(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr94 */


void StopProf P(( void));

#if MLPROTOTYPES
static int _tr95( MLINK mlp)
#else
static int _tr95(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	StopProf();

	res = 1;

L0:	return res;
} /* _tr95 */


void InterruptTime P(( double _tp1));

#if MLPROTOTYPES
static int _tr96( MLINK mlp)
#else
static int _tr96(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	InterruptTime(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr96 */


static struct func {
	int   f_nargs;
	int   manual;
	int   (*f_func)P((MLINK));
	const char  *f_name;
	} _tramps[97] = {
		{ 9, 0, _tr0, "RecMag" },
		{ 0, 2, _tr1, "ExtrudedPolygon" },
		{ 0, 2, _tr2, "ExtrudedPolygon2" },
		{ 0, 2, _tr3, "PlanarPolygon" },
		{ 0, 2, _tr4, "Polyhedron1" },
		{ 1, 0, _tr5, "RecMagsAsExtrPolygons" },
		{ 1, 0, _tr6, "RecMagsAsPolyhedrons" },
		{ 1, 0, _tr7, "ExtPgnsAsPolyhedrons" },
		{ 1, 0, _tr8, "RecognizeRecMags" },
		{ 0, 2, _tr9, "MultGenExtrPolygon" },
		{ 0, 2, _tr10, "MultGenExtrPolygonCur" },
		{ 0, 2, _tr11, "MultGenExtrPolygonMag" },
		{ 0, 2, _tr12, "MultGenExtrRectangle" },
		{ 0, 2, _tr13, "MultGenExtrTriangle" },
		{ 0, 2, _tr14, "ArcPolygon" },
		{10, 0, _tr15, "CylMag" },
		{ 9, 0, _tr16, "RecCur" },
		{12, 0, _tr17, "ArcCur" },
		{12, 0, _tr18, "RaceTrack" },
		{ 0, 2, _tr19, "FlmCur" },
		{ 5, 0, _tr20, "Rectngl" },
		{ 1, 0, _tr21, "Group" },
		{ 2, 0, _tr22, "AddToGroup" },
		{ 1, 0, _tr23, "OutGroupSubObjectKeys" },
		{ 3, 0, _tr24, "BackgroundFieldSource" },
		{ 0, 2, _tr25, "SubdivideElementG3D" },
		{ 0, 2, _tr26, "CutElementG3D" },
		{ 0, 2, _tr27, "SubdivideElementG3DByParPlanes" },
		{ 0, 2, _tr28, "DuplicateElementG3D" },
		{ 1, 0, _tr29, "CreateFromG3DObjectWithSymmetries" },
		{ 1, 0, _tr30, "NumberOfDegOfFreedom" },
		{ 1, 0, _tr31, "MagnOfObj" },
		{ 2, 0, _tr32, "ObjField" },
		{ 2, 0, _tr33, "ScaleCurInObj" },
		{ 1, 0, _tr34, "GeometricalVolume" },
		{ 1, 0, _tr35, "GeometricalVolume" },
		{ 1, 0, _tr36, "GeometricalLimits" },
		{ 3, 0, _tr37, "FldCmpMetForSubdRecMag" },
		{ 0, 2, _tr38, "SetLocMgnInSbdRecMag" },
		{ 3, 0, _tr39, "Translation" },
		{ 7, 0, _tr40, "Rotation" },
		{ 6, 0, _tr41, "PlaneSym" },
		{ 0, 0, _tr42, "FieldInversion" },
		{ 2, 0, _tr43, "TransformObject" },
		{ 3, 0, _tr44, "ApplySymmetry" },
		{ 2, 0, _tr45, "CombineTransformLeft" },
		{ 2, 0, _tr46, "CombineTransformRight" },
		{ 5, 0, _tr47, "LinearMaterial" },
		{ 3, 0, _tr48, "LinearMaterial2" },
		{ 6, 0, _tr49, "NonlinearIsotropMaterial" },
		{ 6, 0, _tr50, "NonlinearIsotropMaterial2" },
		{ 0, 2, _tr51, "NonlinearIsotropMaterial3" },
		{ 0, 2, _tr52, "NonlinearLaminatedMaterialML" },
		{ 0, 2, _tr53, "NonlinearAnisotropMaterial" },
		{ 2, 0, _tr54, "ApplyMaterial" },
		{ 5, 0, _tr55, "MvsH" },
		{ 2, 0, _tr56, "PreRelax" },
		{ 1, 0, _tr57, "ShowInteractMatrix" },
		{ 2, 0, _tr58, "ShowInteractVector" },
		{ 4, 0, _tr59, "ManualRelax" },
		{ 0, 2, _tr60, "AutoRelax" },
		{ 1, 0, _tr61, "UpdateSourcesForRelax" },
		{ 4, 0, _tr62, "SolveGen" },
		{ 6, 0, _tr63, "CompCriterium" },
		{ 0, 2, _tr64, "CompPrecision" },
		{ 4, 0, _tr65, "MultipoleThresholds" },
		{11, 0, _tr66, "Field" },
		{15, 0, _tr67, "ShimSignature" },
		{ 2, 2, _tr68, "FieldArbitraryPointsStruct" },
		{ 9, 0, _tr69, "FieldInt" },
		{ 2, 0, _tr70, "FieldForce" },
		{ 5, 0, _tr71, "FieldEnergy" },
		{ 6, 0, _tr72, "FieldForceThroughEnergy" },
		{ 9, 0, _tr73, "FieldTorqueThroughEnergy" },
		{ 9, 0, _tr74, "ParticleTrajectory" },
		{ 8, 0, _tr75, "FocusingPotential" },
		{24, 0, _tr76, "FocusingKickPer" },
		{ 0, 2, _tr77, "FocusingKickML" },
		{ 0, 0, _tr78, "PhysicalUnits" },
		{ 3, 0, _tr79, "TolForConvergence" },
		{ 1, 0, _tr80, "RandomizationOnOrOff" },
		{ 5, 0, _tr81, "ApplyDrawAttrToElem" },
		{ 1, 0, _tr82, "RemoveDrawAttrFromElem" },
		{ 0, 2, _tr83, "GraphicsForElemWithSymChildsExt" },
		{ 1, 0, _tr84, "GraphicsForElemWithoutSymChilds" },
		{ 0, 2, _tr85, "OpenGL_3D_Viewer" },
		{ 1, 0, _tr86, "DeleteElement" },
		{ 0, 0, _tr87, "DeleteAllElements1" },
		{ 0, 0, _tr88, "DeleteAllElements2" },
		{ 0, 2, _tr89, "DumpElem" },
		{ 0, 2, _tr90, "DumpElemParse" },
		{ 0, 0, _tr91, "RadiaVersion" },
		{ 2, 0, _tr92, "ReturnInput" },
		{ 1, 0, _tr93, "MemAllocMethForIntrctMatr" },
		{ 3, 0, _tr94, "StartProf" },
		{ 0, 0, _tr95, "StopProf" },
		{ 1, 0, _tr96, "InterruptTime" }
		};

static const char* evalstrs[] = {
	"radObjRecMag::usage = \"radObjRecMag[{x,y,z},{wx,wy,wz},{mx,my,mz",
	"}:{0,0,0}] creates a rectangular parallelepipedic block with cen",
	"ter point {x,y,z}, dimensions {wx,wy,wz}, and magnetization {mx,",
	"my,mz}.\"",
	(const char*)0,
	"radObjThckPgn::usage = \"radObjThckPgn[x,lx,{{y1,z1},{y2,z2},...}",
	",a:\\\"x\\\",{mx,my,mz}:{0,0,0}] creates an extruded polygon block; x ",
	"is the position of the block's center of gravity in the extrusio",
	"n direction, lx is the thickness, {{y1,z1},{y2,z2},...} is a lis",
	"t of points describing the polygon in 2D; the extrusion directio",
	"n is defined by the character a (which can be \\\"x\\\", \\\"y\\\" or \\\"z\\\"), ",
	"{mx,my,mz} is the block magnetization.\"",
	(const char*)0,
	"radObjPolyhdr::usage = \"radObjPolyhdr[{{x1,y1,z1},{x2,y2,z2},...",
	"},{{f1i1,f1i2,...},{f2i1,f2i2,...},...},{mx,my,mz}:{0,0,0},J->{j",
	"x,jy,jz}|{{jx,jy,jz},{{djxdy,djxdy,djxdz},{djydy,djydy,djydz},{d",
	"jzdy,djzdy,djzdz}}},Lin->Rel|Abs] creates a polyhedron. {{x1,y1,",
	"z1},{x2,y2,z2},...} is a list of the polyhedron vertex points, {",
	"{f1n1,f1n2,...},{f2n1,f2n2,...},...} is a list of lists of index",
	"es of vertex points defining the polyhedron faces, {mx,my,mz} is",
	" magnetization inside the polyhedron. The option J->... can be u",
	"sed to define constant {jx,jy,jz} or linearly-varying current de",
	"nsity vector inside the polyhedron; the linear dependence can be",
	" defined through 3x3 matrix of coefficients {{djxdy,djxdy,djxdz}",
	",{djydy,djydy,djydz},{djzdy,djzdy,djzdz}}; depending on the valu",
	"e of the option Lin->... this linear dependence is treated with ",
	"respect to the polyhedron center (Lin->Rel, default) or with res",
	"pect to the origin of the Cartesian frame (Lin->Abs).\"",
	(const char*)0,
	"radObjMltExtPgn::usage = \"radObjMltExtPgn[{{{{x11,y11},{x12,y12}",
	",...},z1},{{{x21,y21},{x22,y22},...},z2},...},{mx,my,mz}:{0,0,0}",
	"] attempts to create a convex polyhedron or a set of convex poly",
	"hedrons based on horizontal slices, which should be convex plana",
	"r polygons. The slice polygons are defined by the nested list {{",
	"{{x11,y11},{x12,y12},...},z1},{{{x21,y21},{x22,y22},...},z2},...",
	"}, with {{x11,y11},{x12,y12},...},... describing the polygons in",
	" 2D, and z1, z2,... giving their attitudes (vertical coordinates",
	"). {mx,my,mz} is the magnetization inside the polyhedron(s) crea",
	"ted.\"",
	(const char*)0,
	"radObjMltExtPgnCur::usage = \"radObjMltExtPgnCur[z:0,a:\\\"z\\\",{{{x1,",
	"y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...}},I,Frame->Loc|Lab] ",
	"attempts to create a set of current-carrying convex polyhedron o",
	"bjects by applying a generalized extrusion to the initial planar",
	" convex polygon. The initial polygon is defined for the \\\"attitud",
	"e\\\" z (0 by default) by the list of 2D points {{x1,y1},{x2,y2},..",
	".}, with the  a  character specifying orientation of this polygo",
	"n normal in 3D space: if a = \\\"z\\\" (default orientation), the poly",
	"gon is assumed to be parallel to XY plane of the laboratory fram",
	"e (\\\"y\\\" for ZX plane, \\\"x\\\" for YZ plane). The extrusion can consis",
	"t of a number of \\\"steps\\\", with each step creating one convex pol",
	"yhedron defined optionally by one (or combination of) rotation(s",
	"), and/or translation(s), and/or one homothety: {Rk,Tk,Hk}, k = ",
	"1,2,..., applied to the base polygon (i.e. either the initial ba",
	"se polygon, or the polygon obtained by previous extrusion step).",
	" In case if k-th extrusion step contains one rotation Rk about a",
	"n axis, it is defined as {{xRk,yRk,zRk},{vxRk,vyRk,vzRk},phRk}},",
	" where {xRk,yRk,zRk} and {vxRk,vyRk,vzRk} are respectively 3D co",
	"ordinates of a point and a vector difining the rotation axis, an",
	"d phRk the rotation angle in radians; in case if Rk is a combina",
	"tion of \\\"atomic\\\" rotations about different axes, it should be de",
	"fined as list: {Rk1,Rk2,...}. If k-th extrusion step includes tr",
	"anslation Tk, it must be defined by vector {vxTk,vyTk,vzTk}; opt",
	"ional homothety with respect to the base polygon center of gravi",
	"ty should be defined either by two different coefficients {pxHk,",
	"pyHk} with respect to two orthogonal axes of the base polygon lo",
	"cal frame, or by nested list {{pxHk,pyHk},phHk}, where phHk is r",
	"otation angle of the two homothety axes in radians. A real numbe",
	"r I defines average current in Amperes along the extrusion path.",
	" The Frame->Loc or Frame->Lab option specifies whether the trans",
	"formations at each step of the extrusion path are defined in the",
	" frame of the previous base polygon (Frame->Loc, default), or al",
	"l the transformations are defined in the laboratory frame (Frame",
	"->Lab).\"",
	(const char*)0,
	"radObjMltExtPgnMag::usage = \"radObjMltExtPgnMag[z:0,a:\\\"z\\\",{{{x1,",
	"y1},{x2,y2},...},{{k1,q1},{k2,q2},...}:{{1,1},{1,1},...},{{R1,T1",
	",H1},{R2,T2,H2},...}},{{mx1,my1,mz1},{mx2,my2,mz2},...}:{{0,0,0}",
	",{0,0,0},...},Frame->Loc|Lab,ki->Numb|Size,TriAngMin->...,TriAre",
	"aMax->...,TriExtOpt->\\\"...\\\"] attempts to create a set of uniforml",
	"y magnetized polyhedron objects by applying a generalized extrus",
	"ion to the initial planar convex polygon. The initial polygon is",
	" defined for the \\\"attitude\\\" z (0 by default) by the list of 2D p",
	"oints {{x1,y1},{x2,y2},...}, with the  a  character specifying o",
	"rientation of this polygon normal in 3D space: if a = \\\"z\\\" (defau",
	"lt orientation), the polygon is assumed to be parallel to XY pla",
	"ne of the laboratory frame (\\\"y\\\" for ZX plane, \\\"x\\\" for YZ plane).",
	" The extrusion can consist of a number of \\\"steps\\\", with each ste",
	"p creating one convex polyhedron defined optionally by one (or c",
	"ombination of) rotation(s), and/or translation(s), and/or one ho",
	"mothety: {Rk,Tk,Hk}, k = 1,2,..., applied to the base polygon (i",
	".e. either the initial base polygon, or the polygon obtained by ",
	"previous extrusion step). In case if k-th extrusion step contain",
	"s one rotation Rk about an axis, it is defined as {{xRk,yRk,zRk}",
	",{vxRk,vyRk,vzRk},phRk}}, where {xRk,yRk,zRk} and {vxRk,vyRk,vzR",
	"k} are respectively 3D coordinates of a point and a vector difin",
	"ing the rotation axis, and phRk the rotation angle in radians; i",
	"n case if Rk is a combination of \\\"atomic\\\" rotations about differ",
	"ent axes, it should be defined as a list: {Rk1,Rk2,...}. If k-th",
	" extrusion step includes translation Tk, it must be defined by v",
	"ector {vxTk,vyTk,vzTk}; optional homothety with respect to the b",
	"ase polygon center of gravity should be defined either by two di",
	"fferent coefficients {pxHk,pyHk} with respect to two orthogonal ",
	"axes of the base polygon local frame, or by nested list {{pxHk,p",
	"yHk},phHk}, where phHk is rotation angle of the two homothety ax",
	"es in radians. Optional list {{mx1,my1,mz1},{mx2,my2,mz2},...} d",
	"efines magnetization vectors in each of polyhedrons to be create",
	"d. The Frame->Loc or Frame->Lab option specifies whether the tra",
	"nsformations at each step of the extrusion path are defined in t",
	"he frame of the previous base polygon (Frame->Loc, default), or ",
	"all the transformations are defined in the laboratory frame (Fra",
	"me->Lab). Optionally, the object can be subdivided by (extruded)",
	" triangulation at its creation; this occurs if {{k1,q1},{k2,q2},",
	"...} subdivision (triangulation) parameters for each segment of ",
	"the base polygon border are defined; the meaning of k1, k2,... d",
	"epends on value of the option ki: if ki->Numb (default), then k1",
	", k2,... are subdivision numbers; if ki->Size, they are average ",
	"sizes of sub-segments to be produced; q1, q2,... are ratios of t",
	"he last-to-first sub-segment lengths; the TriAngMin option defin",
	"es minimal angle of triangles to be produced (in degrees, defaul",
	"t is 20); the TriAreaMax option defines maximal area of traingle",
	"s to be produced (in mm^2, not defined by default); the ExtOpt o",
	"ption allows to specify additional parameters for triangulation ",
	"function in a string.\"",
	(const char*)0,
	"radObjMltExtRtg::usage = \"radObjMltExtRtg[{{{x1,y1,z1},{wx1,wy1}",
	"},{{x2,y2,z2},{wx2,wy2}},...},{mx,my,mz}:{0,0,0}] attempts to cr",
	"eate a convex polyhedron or a set of convex polyhedrons based on",
	" horizontal slices of rectangular shape. The slice rectangles ar",
	"e defined by the nested list {{{x1,y1,z1},{wx1,wy1}},{{x2,y2,z2}",
	",{wx2,wy2}},...}, with {x1,y1,z1}, {x2,y2,z2},... being center p",
	"oints of the rectangles, and {wx1,wy1}, {wx2,wy2},... their dime",
	"nsions. {mx,my,mz} is the magnetization inside the polyhedron(s)",
	" created.\"",
	(const char*)0,
	"radObjMltExtTri::usage = \"radObjMltExtTri[x,lx,{{y1,z1},{y2,z2},",
	"...},{{k1,q1},{k2,q2},...},a:\\\"x\\\",{mx,my,mz}:{0,0,0},ki->Numb|Siz",
	"e,TriAngMin->...,TriAreaMax->...,TriExtOpt->\\\"...\\\"] creates trian",
	"gulated extruded polygon block, i.e. a straight prism with its b",
	"ases being (possibly non-convex) polygons subdivided by triangul",
	"ation; x is the position of the block's center of gravity in the",
	" extrusion direction, lx is the thickness, {{y1,z1},{y2,z2},...}",
	" is a list of points describing the polygon in 2D; {{k1,q1},{k2,",
	"q2},...} are subdivision (triangulation) parameters for each seg",
	"ment of the base polygon border; the meaning of k1, k2,... depen",
	"ds on value of the option ki: if ki->Numb (default), then k1, k2",
	",... are subdivision numbers; if ki->Size, they are average size",
	"s of sub-segments to be produced; q1, q2,... are ratios of the l",
	"ast-to-first sub-segment lengths; the extrusion direction is def",
	"ined by the character a (which can be \\\"x\\\", \\\"y\\\" or \\\"z\\\"); {mx,my,m",
	"z} is magnetization inside the block. The TriAngMin option defin",
	"es minimal angle of triangles to be produced (in degrees, defaul",
	"t is 20); the TriAreaMax option defines maximal area of traingle",
	"s to be produced (in mm^2, not defined by default); the ExtOpt o",
	"ption allows to specify additional parameters for triangulation ",
	"function in a string.\"",
	(const char*)0,
	"radObjArcPgnMag::usage = \"radObjArcPgnMag[{x,y},a,{{r1,z1},{r2,z",
	"2},...},{phimin,phimax},nseg,\\\"sym|nosym\\\":\\\"nosym\\\",{mx,my,mz}:{0,0",
	",0}] creates a finite-length arc of polygonal cross-section with",
	" the position and orientation of the rotation axis defined by pa",
	"ir of coordinates {x,y} and character a (which can be \\\"x\\\", \\\"y\\\" o",
	"r \\\"z\\\"), the cross-section 2D polygon {{r1,z1},{r2,z2},...}, init",
	"ial and final rotation angles {phimin,phimax}, number of sectors",
	" (segments) vs azimuth nseg, and magnetization vector {mx,my,mz}",
	". Depending on the value of the \\\"sym|nosym\\\" switch, the magnetiz",
	"ation vectors in nseg sector polyhedrons are either assumed to s",
	"atisfy rotational symmetry conditions (\\\"sym\\\"), or are assumed in",
	"dependent (i.e. will be allowed to vary independently at further",
	" relaxation).\"",
	(const char*)0,
	"radObjCylMag::usage = \"radObjCylMag[{x,y,z},r,h,nseg,a:\\\"z\\\",{mx,m",
	"y,mz}:{0,0,0}] creates a cylindrical magnet approximated by a ri",
	"ght polygon with center point {x,y,z}, base radius r, height h, ",
	"number of segments nseg, orientation of the rotation axis define",
	"d by character a (which can be \\\"x\\\", \\\"y\\\" or \\\"z\\\"), and magnetizati",
	"on vector {mx,my,mz}.\"",
	(const char*)0,
	"radObjRecCur::usage = \"radObjRecCur[{x,y,z},{wx,wy,wz},{jx,jy,jz",
	"}] creates a current carrying rectangular parallelepipedic block",
	" with center point {x,y,z}, dimensions {wx,wy,wz}, and current d",
	"ensity {jx,jy,jz}.\"",
	(const char*)0,
	"radObjArcCur::usage = \"radObjArcCur[{x,y,z},{rmin,rmax},{phimin,",
	"phimax},h,nseg,j,\\\"man|auto\\\":\\\"man\\\",a:\\\"z\\\"] creates a current-carry",
	"ing finite-length arc of rectangular cross-section, center point",
	" {x,y,z}, inner and outer radii {rmin,rmax}, initial and final a",
	"ngles {phimin,phimax}, height h, number of segments nseg, and az",
	"imuthal current density j. According to the value of the \\\"man|au",
	"to\\\" switch, the field from the arc is computed based on the numb",
	"er of segments nseg (\\\"man\\\"), or on the general absolute precisio",
	"n level specified by the function radFldCmpCrt (\\\"auto\\\"). The ori",
	"entation of the rotation axis is defined by the character a (whi",
	"ch can be either \\\"x\\\", \\\"y\\\" or \\\"z\\\").\"",
	(const char*)0,
	"radObjRaceTrk::usage = \"radObjRaceTrk[{x,y,z},{rmin,rmax},{lx,ly",
	"},h,nseg,j,\\\"man|auto\\\":\\\"man\\\",a:\\\"z\\\"] creates a current carrying ra",
	"cetrack coil consisting of four 90-degree bents connected by fou",
	"r straight parts of rectangular straight section, center  point ",
	"{x,y,z}, inner and outer bent radii {rmin,rmax}, straight sectio",
	"n lengths {lx,ly}, height h, number of segments in bents nseg, a",
	"nd azimuthal current density j. According to the value of the \\\"m",
	"an|auto\\\" switch, the field from the bents is computed based on t",
	"he number of segments nseg (\\\"man\\\"), or on the general absolute p",
	"recision level specified by the function radFldCmpCrt (\\\"auto\\\"). ",
	"The orientation of the racetrack axis is defined by the characte",
	"r a (which can be either \\\"x\\\", \\\"y\\\" or \\\"z\\\").\"",
	(const char*)0,
	"radObjFlmCur::usage = \"radObjFlmCur[{{x1,y1,z1},{x2,y2,z2},...},",
	"i] creates a filament polygonal line conductor defined by the se",
	"quence of points {{x1,y1,z1},{x2,y2,z2},...} with current i.\"",
	(const char*)0,
	"radObjCnt::usage = \"radObjCnt[{obj1,obj2,...}] creates a contain",
	"er object for the objects {obj1,obj2,...}.\"",
	(const char*)0,
	"radObjAddToCnt::usage = \"radObjAddToCnt[cnt,{obj1,obj2,...}] add",
	"s objects {obj1,obj2,...} to the container object cnt.\"",
	(const char*)0,
	"radObjCntStuf::usage = \"radObjCntStuf[obj] gives a list of gener",
	"al indexes of the objects present in container if obj is a conta",
	"iner; or returns {obj} if obj is not a container.\"",
	(const char*)0,
	"radObjBckg::usage = \"radObjBckg[{bx,by,bz}] creates a source of ",
	"uniform background magnetic field {bx,by,bz}.\"",
	(const char*)0,
	"radObjDivMag::usage = \"radObjDivMag[obj,{{k1,q1},{k2,q2},{k3,q3}",
	"},{\\\"pln\\\",{n1x,n1y,n1z},{n2x,n2y,n2z},{n3x,n3y,n3z}}|{\\\"cyl\\\",{{ax,",
	"ay,az},{vx,vy,vz}},{px,py,pz},rat},kxkykz->Numb|Size,Frame->Loc|",
	"Lab|LabTot] subdivides the object obj. The main subdivision para",
	"meters are defined by the list {{k1,q1},{k2,q2},{k3,q3}}. The me",
	"aning of k1, k2 and k3 depends on the value of the option kxkykz",
	": if kxkykz->Numb (default), then k1, k2 and k3 are subdivision ",
	"numbers; if kxkykz->Size, they are average sizes of the sub-obje",
	"cts to be produced. q1, q2 and q3 in any case are ratios of the ",
	"last-to-first sub-object sizes. The third variable is optional. ",
	"If it is not used, the subdivision is performed in directions X,",
	" Y and Z. If it is used in the form {\\\"pln\\\",...}, the function pe",
	"rforms subdivision by three sets of parallel planes normal to th",
	"e vectors {n1x,n1y,n1z}, {n2x,n2y,n2z} and {n3x,n3y,n3z}. The di",
	"stances between the parallel planes are defined by the parameter",
	"s {k1,q1},{k2,q2} and {k3,q3}. If the third variable is used in ",
	"the form {\\\"cyl\\\",...}, the function performs subdivision by a sys",
	"tem of coaxial elliptic cylinders. The cylinder axis is defined ",
	"by the point {ax,ay,az} and vector {vx,vy,vz}. One of two axes o",
	"f the cylinder base ellipses is exactly the perpendicular from t",
	"he point {px,py,pz} to the cylinder axis; rat is the ratio of th",
	"e ellipse axes lengths. In the case of the subdivision by ellipt",
	"ic cylinders, the parameters {k1,q1},{k2,q2} and {k3,q3} corresp",
	"ond to radial, azimuthal, and axial directions respectively. If ",
	"the option Frame is set to Frame->Loc (default), the subdivision",
	" is performed in local frame(s) of the object(s); if Frame->Lab ",
	"or Frame->LabTot, the subdivision is performed in the laboratory",
	" frame. The actions of Frame->Lab and Frame->LabTot differ for c",
	"ontainers only: Frame->Lab means that each of the objects in the",
	" container is subdivided separately; Frame->LabTot means that th",
	"e objects in the container are subdivided as one object, by the ",
	"same planes.\"",
	(const char*)0,
	"radObjCutMag::usage = \"radObjCutMag[obj,{x,y,z},{nx,ny,nz},Frame",
	"->Loc|Lab|LabTot] cuts the object obj by a plane normal to the v",
	"ector {nx,ny,nz} and passing through the point {x,y,z}. The Fram",
	"e->Loc, Frame->Lab or Frame->LabTot option specifies whether the",
	" cuting plane is defined in the local frame of the object obj or",
	" in the laboratory frame (default). The actions of Frame->Lab an",
	"d Frame->LabTot differ for containers only: Frame->Lab means tha",
	"t each of the objects in the container is cut separately; Frame-",
	">LabTot means that the objects in the container are cut as one o",
	"bject, by the same plane. The function returns a list of indexes",
	" of the objects produced by the cutting.\"",
	(const char*)0,
	"radObjGeoVol::usage = \"radObjGeoVol[obj] computes geometrical vo",
	"lume of the object obj.\"",
	(const char*)0,
	"radObjGeoLim::usage = \"radObjGeoLim[obj] computes geometrical li",
	"mits of the object obj in the laboratory frame. Returns {xmin, x",
	"max, ymin, ymax, zmin, zmax}.\"",
	(const char*)0,
	"radObjDpl::usage = \"radObjDpl[obj,FreeSym->False|True] duplicate",
	"s the object obj. The option FreeSym->False|True specifies wheth",
	"er the symmetries (transformations with multiplicity more than o",
	"ne) previously applied to the object obj should be simply copied",
	" at the duplication (FreeSym->False, default), or a container of",
	" new independent objects should be created in place of any symme",
	"try previously applied to the object obj. In both cases the fina",
	"l object created by the duplication has exactly the same geometr",
	"y as the initial object obj.\"",
	(const char*)0,
	"radObjDrwAtr::usage = \"radObjDrwAtr[obj,{r,g,b},thcn] applies RG",
	"B color {r,g,b} and line thickness thcn to object obj.\"",
	(const char*)0,
	"radObjDrw::usage = \"radObjDrw[obj] prepares a set of 3D graphica",
	"l primitives representing object obj in space.\"",
	(const char*)0,
	"radObjDrwOpenGL::usage = \"radObjDrwOpenGL[obj,EdgeLines->True|Fa",
	"lse,Faces->True|False,Axes->True|False] starts an application fo",
	"r viewing 3D geometry of the object obj. The viewer is based on ",
	"the GLUT / OpenGL graphics library. The option EdgeLines->True|F",
	"alse (default EdgeLines->True) highlights the edge lines of obje",
	"cts; the option Faces->True|False (default Faces->True) shows fa",
	"ces of the objects; the option Axes->True|False (default Axes->T",
	"rue) shows the Cartesian frame axes.\"",
	(const char*)0,
	"radObjDegFre::usage = \"radObjDegFre[obj] gives number of degrees",
	" of freedom for the relaxation of the object obj.\"",
	(const char*)0,
	"radObjM::usage = \"radObjM[obj] returns coordinates of geometrica",
	"l center point and magnetization vector of the object obj at tha",
	"t point. If obj is a container, a list of the container members'",
	" center points and their magnetizations is returned.\"",
	(const char*)0,
	"radObjCenFld::usage = \"radObjCenFld[obj,\\\"A|B|H|J|M\\\"] returns coo",
	"rdinates of geometrical center point of the object obj and a fie",
	"ld characteristic vector at that point. The type of field charac",
	"teristic is defined by the second parameter (character); it can ",
	"be one of the following: \\\"A\\\" for vector potential, \\\"B\\\" for magne",
	"tic field induction, \\\"H\\\" for magnetic field strength, \\\"J\\\" for cu",
	"rrent density, \\\"M\\\" for magnetization. If obj is a container, a l",
	"ist of the container members' center points and their field char",
	"acteristics is returned.\"",
	(const char*)0,
	"radObjScaleCur::usage = \"radObjScaleCur[obj,k] scales current (d",
	"ensity) in the obj by multiplying it by k (if obj is a current-c",
	"arying object). If obj is a container, the current (density) sca",
	"ling applies to all its members.\"",
	(const char*)0,
	"radTrfTrsl::usage = \"radTrfTrsl[{vx,vy,vz}] creates a translatio",
	"n with vector {vx,vy,vz}.\"",
	(const char*)0,
	"radTrfRot::usage = \"radTrfRot[{x,y,z},{vx,vy,vz},phi] creates a ",
	"rotation of angle phi around the axis defined by the point {x,y,",
	"z} and the vector {vx,vy,vz}.\"",
	(const char*)0,
	"radTrfPlSym::usage = \"radTrfPlSym[{x,y,z},{nx,ny,nz}] creates a ",
	"symmetry with respect to plane defined by point {x,y,z} and norm",
	"al vector {nx,ny,nz}.\"",
	(const char*)0,
	"radTrfInv::usage = \"radTrfInv[] creates a field inversion.\"",
	(const char*)0,
	"radTrfCmbL::usage = \"radTrfCmbL[OrigTrf,trf] multiplies original",
	" transformation OrigTrf by another transformation trf from left.",
	"\"",
	(const char*)0,
	"radTrfCmbR::usage = \"radTrfCmbR[OrigTrf,trf] multiplies original",
	" transformation OrigTrf by another transformation trf from right",
	".\"",
	(const char*)0,
	"radTrfOrnt::usage = \"radTrfOrnt[obj,trf] orients object obj by a",
	"pplying transformation trf to it once.\"",
	(const char*)0,
	"radTrfMlt::usage = \"radTrfMlt[obj,trf,mlt] creates mlt-1 objects",
	". Each object is derived from the previous one by applying the t",
	"ransformation trf. Following this, the object obj becomes equiva",
	"lent to mlt different objects.\"",
	(const char*)0,
	"radMatLin::usage = \"radMatLin[{ksipar,ksiper},mr] or radMatLin[{",
	"ksipar,ksiper},{mrx,mry,mrz}] creates a linear anisotropic magne",
	"tic material with susceptibilities parallel (perpendicular) to t",
	"he easy magnetization axis given by ksipar (ksiper). In the firs",
	"t form of the function, mr is the magnitude of the remanent magn",
	"etization vector; the direction of the easy magnetisation axis i",
	"s set up by the magnetization vector in the object to which the ",
	"material is applied (the magnetization vector is specified at th",
	"e object creation). In the second form, {mrx,mry,mrz} is the rem",
	"anent magnetization vector explicitly defining the direction of ",
	"the easy magnetization axis in any object to which the material ",
	"is later applied.\"",
	(const char*)0,
	"radMatSatIso::usage = \"radMatSatIso[{ksi1,ms1},{ksi2,ms2},{ksi3,",
	"ms3}] creates a nonlinear isotropic magnetic material with the m",
	"agnetization magnitude M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*",
	"H/ms2) + ms3*tanh(ksi3*H/ms3) where H is the magnitude of the ma",
	"gnetic field strength vector. radMatSatIso[{{H1,M1},{H2,M2},...}",
	"] creates a nonlinear isotropic magnetic material with the M ver",
	"sus H curve defined by the list of pairs {{H1,M1},{H2,M2},...} i",
	"n Tesla.\"",
	(const char*)0,
	"radMatSatLam::usage = \"radMatSatLam[data,p,{nx,ny,nz}:0] creates",
	" laminated nonlinear anisotropic magnetic material with stacking",
	" factor p and vector normal to the lamination planes given by {n",
	"x,ny,nz}. data is a list of pairs of numbers defining the materi",
	"al dependence M(H) as for isotropic case. If this list contains ",
	"only 3 elements, it is interpreted as {{ksi1,ms1},{ksi2,ms2},{ks",
	"i3,ms3}} and M(H) is then computed by the formula M = ms1*tanh(k",
	"si1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3); if it ",
	"contains more than 3 elements, it is interpreted as tabulated de",
	"pendence M(H): {{H1,M1},{H2,M2},...} in Tesla. If the vector nor",
	"mal {nx,ny,nz} is not supplied, the lamination planes are assume",
	"d to be perpendicular to the magnetization vector in the object ",
	"to which the material is applied (the magnetization vector is sp",
	"ecified at the object creation).\"",
	(const char*)0,
	"radMatSatAniso::usage = \"radMatSatAniso[datapar,dataper] where d",
	"atapar can be {{ksi1,ms1,hc1},{ksi2,ms2,hc2},{ksi3,ms3,hc3},{ksi",
	"0,hc0}} or ksi0, and dataper can be {{ksi1,ms1},{ksi2,ms2},{ksi3",
	",ms3},ksi0} or ksi0 - creates a nonlinear anisotropic magnetic m",
	"aterial. If the first argument is set to {{ksi1,ms1,hc1},{ksi2,m",
	"s2,hc2},{ksi3,ms3,hc3},{ksi0,hc0}}, the magnetization vector com",
	"ponent parallel to the easy axis is computed as ms1*tanh(ksi1*(h",
	"pa-hc1)/ms1) + ms2*tanh(ksi2*(hpa-hc2)/ms2) + ms3*tanh(ksi3*(hpa",
	"-hc3)/ms3) + ksi0*(hpa-hc0), where hpa is the field strength vec",
	"tor component parallel to the easy axis. If the second argument ",
	"is set to {{ksi1,ms1},{ksi2,ms2},{ksi3,ms3},ksi0}, the magnetiza",
	"tion vector component perpendicular to the easy axis is computed",
	" as ms1*tanh(ksi1*hpe/ms1) + ms2*tanh(ksi2*hpe/ms2) + ms3*tanh(k",
	"si3*hpe/ms3) + ksi0*hpe, where hpe is the field strength vector ",
	"component perpendicular to the easy axis. If the first or second",
	" argument is set to ksi0, the magnetization component parallel (",
	"perpendicular) to the easy axis is computed by ksi0*hp, where hp",
	" is the corresponding component of field strength vector. At lea",
	"st one of the magnetization vector components should non-linearl",
	"y depend on the field strength. The direction of the easy magnet",
	"isation axis is set up by the magnetization vector in the object",
	" to which the material is later applied.\"",
	(const char*)0,
	"radMatApl::usage = \"radMatApl[obj,mat] applies material mat to o",
	"bject obj.\"",
	(const char*)0,
	"radMatMvsH::usage = \"radMatMvsH[obj,\\\"mx|my|mz\\\"|\\\"\\\",{hx,hy,hz}] co",
	"mputes magnetization from magnetic field strength vector (hx,hy,",
	"hz) for the material of the object obj; the magnetization compon",
	"ents are specified by the second argument.\"",
	(const char*)0,
	"radRlxPre::usage = \"radRlxPre[obj,srcobj:0] builds an interactio",
	"n matrix for the object obj, treating the object srcobj as addit",
	"ional external field source.\"",
	(const char*)0,
	"radRlxMan::usage = \"radRlxMan[intrc,meth,iternum,rlxpar] execute",
	"s manual relaxation procedure for interaction matrix intrc using",
	" method number meth, by making iternum iterations with relaxatio",
	"n parameter value rlxpar.\"",
	(const char*)0,
	"radRlxAuto::usage = \"radRlxAuto[intrc,prec,maxiter,meth:4,ZeroM-",
	">True|False] executes automatic relaxation procedure with the in",
	"teraction matrix intrc using the method number meth. Relaxation ",
	"stops whenever the change in magnetization (averaged over all su",
	"b-elements) between two successive iterations is smaller than pr",
	"ec or the number of iterations is larger than maxiter. The optio",
	"n value ZeroM->True (default) starts the relaxation by setting t",
	"he magnetization values in all paricipating objects to zero; Zer",
	"oM->False starts the relaxation with the existing magnetization ",
	"values in the sub-volumes.\"",
	(const char*)0,
	"radRlxUpdSrc::usage = \"radRlxUpdSrc[intrc] updates external fiel",
	"d data for the relaxation (to take into account e.g. modificatio",
	"n of currents in coils, if any) without rebuilding the interacti",
	"on matrix.\"",
	(const char*)0,
	"radSolve::usage = \"radSolve[obj,prec,maxiter,meth:4] builds inte",
	"raction matrix for the object obj and executes relaxation proced",
	"ure using the method number meth. Relaxation stops whenever the ",
	"change in magnetization (averaged over all sub-elements) between",
	" two successive iterations is smaller than prec or the number of",
	" iterations is larger than maxiter.\"",
	(const char*)0,
	"radFldCmpCrt::usage = \"radFldCmpCrt[prcB,prcA,prcBint,prcFrc,prc",
	"TrjCrd,prcTrjAng] sets general absolute accuracy levels for comp",
	"utation of field induction (prcB), vector potential (prcA), indu",
	"ction integrals along straight line (prcBint), field force (prcF",
	"rc), relativistic particle trajectory coordinates (prcTrjCrd) an",
	"d angles (prcTrjAng).\"",
	(const char*)0,
	"radFldCmpPrc::usage = \"radFldCmpPrc[PrcB->prb,PrcA->pra,PrcBInt-",
	">prbint,PrcForce->prfrc,PrcTorque->prtrq,PrcEnergy->pre,PrcCoord",
	"->prcrd,PrcAngle->prang] sets general absolute accuracy levels f",
	"or computation of magnetic field induction, vector potential, in",
	"duction integral along straight line, field force, torque, energ",
	"y; relativistic charged particle trajectory coordinates and angl",
	"es. The function works in line with the Mathematica mechanism of",
	" Options. PrcB, PrcA, PrcBInt, PrcForce, PrcTorque, PrcEnergy, P",
	"rcCoord, PrcAngle are names of the options; prb, pra, prbint, pr",
	"frc, prtrq, pre, prcrd, prang are the corresponding values (real",
	" numbers specifying the accuracy levels).\"",
	(const char*)0,
	"radFld::usage = \"radFld[obj,\\\"bx|by|bz|hx|hy|hz|ax|ay|az|mx|my|mz",
	"\\\"|\\\"\\\",{x,y,z}|{{x1,y1,z1},{x2,y2,z2},...}] computes magnetic fiel",
	"d created by the object obj in point(s) {x,y,z} ({x1,y1,z1},{x2,",
	"y2,z2},...). The field component is specified by the second inpu",
	"t variable. The function accepts a list of 3D points of arbitrar",
	"y nestness: in this case it returns the corresponding list of ma",
	"gnetic field values.\"",
	(const char*)0,
	"radFldLst::usage = \"radFldLst[obj,\\\"bx|by|bz|hx|hy|hz|ax|ay|az|mx",
	"|my|mz\\\"|\\\"\\\",{x1,y1,z1},{x2,y2,z2},np,\\\"arg|noarg\\\":\\\"noarg\\\",strt:0.]",
	" computes magnetic field created by object obj in np equidistant",
	" points along a line segment from {x1,y1,z1} to {x2,y2,z2}; the ",
	"field component is specified by the second input variable; the \\\"",
	"arg|noarg\\\" string variable specifies whether to output a longitu",
	"dinal position for each point where the field is computed, and s",
	"trt gives the start-value for the longitudinal position.\"",
	(const char*)0,
	"radFldInt::usage = \"radFldInt[obj,\\\"inf|fin\\\",\\\"ibx|iby|ibz\\\"|\\\"\\\",{x1",
	",y1,z1},{x2,y2,z2}] computes magnetic field induction integral p",
	"roduced by the object obj along a straight line specified by poi",
	"nts {x1,y1,z1} and {x2,y2,z2}; depending on the second variable ",
	"value, the integral is infinite (\\\"inf\\\") or finite from {x1,y1,z1",
	"} to {x2,y2,z2} (\\\"fin\\\"); the field integral component is specifi",
	"ed by the third input variable. The unit is Tesla x millimeter.\"",
	(const char*)0,
	"radFldFrc::usage = \"radFldFrc[obj,shape] computes force of the f",
	"ield produced by the object obj into a shape defined by shape. s",
	"hape can be the result of radObjRecMag[..] (parallelepiped) or r",
	"adFldFrcShpRtg[..] (rectangular surface).\"",
	(const char*)0,
	"radFldFrcShpRtg::usage = \"radFldFrcShpRtg[{x,y,z},{wx,wy}] creat",
	"es a rectangle with central point {x,y,z} and dimensions {wx,wy}",
	"; to be used for field force computation.\"",
	(const char*)0,
	"radFldEnr::usage = \"radFldEnr[objdst,objsrc] or radFldEnr[objdst",
	",objsrc,{kx,ky,kz}] computes potential energy (in Joule) of the ",
	"object objdst in the field created by the object objsrc. The fir",
	"st form of the function performes the computation based on absol",
	"ute accuracy value for the energy (by default 10 Joule; can be m",
	"odified by the function radFldCmpPrc). The second form performs ",
	"the computation based on the destination object subdivision numb",
	"ers {kx,ky,kz}.\"",
	(const char*)0,
	"radFldEnrFrc::usage = \"radFldEnrFrc[objdst,objsrc,\\\"fx|fy|fz|\\\"\\\"] ",
	"or radFldEnrFrc[objdst,objsrc,\\\"fx|fy|fz|\\\"\\\",{kx,ky,kz}] computes ",
	"force (in Newton) acting on the object objdst in the field produ",
	"ced by the object objsrc. The first form of the function perform",
	"es the computation based on absolute accuracy value for the forc",
	"e (by default 10 Newton; can be modified by the function radFldC",
	"mpPrc). The second form performs the computation based on the de",
	"stination object subdivision numbers {kx,ky,kz}.\"",
	(const char*)0,
	"radFldEnrTrq::usage = \"radFldEnrTrq[objdst,objsrc,\\\"tx|ty|tz|\\\"\\\",{",
	"x,y,z}] or radFldEnrTrq[objdst,objsrc,\\\"tx|ty|tz|\\\"\\\",{x,y,z},{kx,k",
	"y,kz}] computes torque (in Newton*mm) with respect to point {x,y",
	",z}, acting on the object objdst in the field produced by the ob",
	"ject objsrc. The first form of the function performes the comput",
	"ation based on absolute accuracy value for the torque (by defaul",
	"t 10 Newton*mm; can be modified by the function radFldCmpPrc). T",
	"he second form performs the computation based on the destination",
	" object subdivision numbers {kx,ky,kz}.\"",
	(const char*)0,
	"radFldPtcTrj::usage = \"radFldPtcTrj[obj,E,{x0,dxdy0,z0,dzdy0},{y",
	"0,y1},np] computes transverse coordinates and its derivatives (a",
	"ngles) of a relativistic charged particle trajectory in 3D magne",
	"tic field produced by the object obj, using the 4th order Runge-",
	"Kutta integration. The particle energy is E [GeV], initial trans",
	"verse coordinates and derivatives are {x0,dxdy0,z0,dzdy0}; the l",
	"ongitudinal coordinate y is varied from y0 to y1 in np steps. Al",
	"l positions are in millimeters and angles in radians.\"",
	(const char*)0,
	"radFldFocPot::usage = \"radFldFocPot[obj,{x1,y1,z1},{x2,y2,z2},np",
	"] computes the potential for trajectory of relativistic charged ",
	"particle in magnetic field produced by the object obj. The integ",
	"ration is made from {x1,y1,z1} to {x2,y2,z2} with np equidistant",
	" points.\"",
	(const char*)0,
	"radFldFocKickPer::usage = \"radFldFocKickPer[obj,{x1,y1,z1},{nsx,",
	"nsy,nsz},per,nper,{n1x,n1y,n1z},r1,np1,r2,np2,com:\\\"\\\",{nh:1,nps:8",
	",d1:0,d2:0},\\\"T2m2|rad|microrad\\\":\\\"T2m2\\\",en:1,\\\"fix|tab\\\":\\\"fix\\\"] com",
	"putes matrices of 2nd order kicks of trajectory of relativistic ",
	"charged particle in periodic magnetic field produced by the obje",
	"ct obj. The longitudinal integration along one period starts at ",
	"point {x1,y1,z1} and is done along direction pointed by vector {",
	"nsx,nsy,nsz}; per is period length, nper is number of full perio",
	"ds; one direction of the transverse grid is pointed by vector {n",
	"1x,n1y,n1z}, the other transverse direction is given by vector p",
	"roduct of {n1x,n1y,n1z} and {nsx,nsy,nsz}; r1 and r2 are ranges ",
	"of the transverse grid, np1 and np2 are corresponding numbers of",
	" points; com is arbitrary string comment; nh is maximum number o",
	"f magnetic field harmonics to treat (default 1), nps is number o",
	"f longitudinal points (default 8), d1 and d2 are steps of transv",
	"erse differentiation (by default equal to the steps of the trans",
	"verse grid); the \\\"T2m2|rad|microrad\\\" string variable specifies t",
	"he units for the resulting 2nd order kick values (default \\\"T2m2\\\"",
	"); en is electron elergy in GeV (optional, required only if unit",
	"s are \\\"rad\\\" or \\\"microrad\\\"); the \\\"fix|tab\\\":\\\"fix\\\" string variable ",
	"specifies the format of the output data string, \\\"fix\\\" for fixed-",
	"width, \\\"tab\\\" for tab-delimited (i.e. element [[6]] of the output",
	" list, default \\\"fix\\\"). Returns list containing: [[1]]- matrix of",
	" kick values in the first transverse direction, [[2]]- matrix of",
	" kick values in the second transverse direction, [[3]]- matrix o",
	"f longitudinally-integrated squared transverse magnetic field ca",
	"lculated on same transverse mesh as kicks, [[4]],[[5]]- lists of",
	" positions defining the transverse grid, [[6]]- formatted string",
	" containing the computed results (for saving into a text file).\"",
	(const char*)0,
	"radFldFocKick::usage = \"radFldFocKick[obj,{x1,y1,z1},{nsx,nsy,ns",
	"z},{ds1,ds2,ds3,...},nps,{n1x,n1y,n1z},r1,np1,r2,np2,com:\\\"\\\",{d1:",
	"0,d2:0}] computes matrices of 2nd order kicks of trajectory of r",
	"elativistic charged particle in arbitrary magnetic field produce",
	"d by the object obj. PLEASE NOTE that this is a time-consuming c",
	"alculation; therefore try using the radFldFocKickPer function we",
	"re applicable. The longitudinal integration starts at point {x1,",
	"y1,z1} and is done along direction pointed by vector {nsx,nsy,ns",
	"z}; the transverse matrices of the kick values are computed for ",
	"the planes located at distances {ds1,ds2,ds3,...} from {x1,y1,z1",
	"}; nps is total number of points for longitudinal integration; o",
	"ne direction of the transverse grid is pointed by vector {n1x,n1",
	"y,n1z}, the other transverse direction is given by vector produc",
	"t of {n1x,n1y,n1z} and {nsx,nsy,nsz}; r1 and r2 are ranges of th",
	"e transverse grid, np1 and np2 are corresponding numbers of poin",
	"ts; com is arbitrary string comment; d1 and d2 are steps of tran",
	"sverse differentiation (by default equal to the steps of the tra",
	"nsverse grid). Returns list containing: [[1]]- list of triplets ",
	"of 2D matrices representing kick values in the first and second ",
	"transverse directions and the longitudinally-integrated squared ",
	"transverse magnetic field, including the total matrices computed",
	" for whole range of longitudinal position and partial matrices c",
	"orresponding to given longitudinal positions, [[2]],[[3]]- lists",
	" of positions defining the transverse grid, [[4]]- formatted str",
	"ing containing all the computed results (for saving into a text ",
	"file).\"",
	(const char*)0,
	"radFldShimSig::usage = \"radFldShimSig[obj,\\\"bx|by|bz|hx|hy|hz|ibx",
	"|iby|ibz\\\"|\\\"\\\",{dx,dy,dz},{x1,y1,z1},{x2,y2,z2},np,{vix,viy,viz}:{",
	"0,0,0}] computes virtual \\\"shim signature\\\", i.e. variation of mag",
	"netic field component defined by the second variable, introduced",
	" by displacement {dx,dy,dz} of magnetic field source object obj.",
	" The field variation is computed at np equidistant points along ",
	"a line segment from {x1,y1,z1} to {x2,y2,z2}; the vector {vix,vi",
	"y,viz} is taken into account if a field integral variation shoul",
	"d be computed: in this case, it defines orientation of the integ",
	"ration line.\"",
	(const char*)0,
	"radFldLenTol::usage = \"radFldLenTol[abs,rel,zero:0] sets absolut",
	"e and relative randomization magnitudes for all the length value",
	"s, including coordinates and dimensions of the objects producing",
	" magnetic field, and coordinates of points where the field is co",
	"mputed. Optimal values of the variables can be: rel=10^(-11), ab",
	"s=L*rel, zero=abs, where L is the distance scale value (in mm) f",
	"or the problem to be solved. Too small randomization magnitudes ",
	"can result in run-time code errors.\"",
	(const char*)0,
	"radFldLenRndSw::usage = \"radFldLenRndSw[\\\"on|off\\\"] switches on or",
	" off the randomization of all the length values. The randomizati",
	"on magnitude can be set by the function radFldLenTol.\"",
	(const char*)0,
	"radFldUnits::usage = \"radFldUnits[] shows the physical units cur",
	"rently in use.\"",
	(const char*)0,
	"radUtiDel::usage = \"radUtiDel[elem] deletes element elem.\"",
	(const char*)0,
	"radUtiDelAll::usage = \"radUtiDelAll[] deletes all previously cre",
	"ated elements.\"",
	(const char*)0,
	"radUtiDmp::usage = \"radUtiDmp[elem,\\\"asc|bin\\\":\\\"asc\\\"] outputs info",
	"rmation about elem, which can be either one element or list of e",
	"lements; second argument specifies whether the output should be ",
	"in ASCII (\\\"asc\\\", default) or in Binary (\\\"bin\\\") format.\"",
	(const char*)0,
	"radUtiDmpPrs::usage = \"radUtiDmpPrs[bstr] parses byte-string bst",
	"r produced by radUtiDmp[elem,\\\"bin\\\"] and attempts to instantiate ",
	"elem objects(s); returns either index of one instantiated object",
	" (if elem was an index of one object) or a list of indexes of in",
	"stantiated objects (if elem was a list of objects).\"",
	(const char*)0,
	"radUtiIntrptTim::usage = \"radUtiIntrptTim[t] sets interruption t",
	"ime quanta in seconds for platforms with no preemptive multitask",
	"ing.\"",
	(const char*)0,
	"radUtiVer::usage = \"radUtiVer[] returns the version number of th",
	"e Radia executable file.\"",
	(const char*)0,
	"rAdObjPgn::usage = \"rAdObjPgn[z,{{x1,y1},{x2,y2},...},{mx,my,mz}",
	":{0,0,0}] \"",
	(const char*)0,
	"rAdObjSetLocMgn::usage = \"rAdObjSetLocMgn[obj,{{{ix1,iy1,iz1},{m",
	"x1,my1,mz1}},{{ix2,iy2,iz2},{mx2,my2,mz2}},...}] \"",
	(const char*)0,
	"rAdObjDrwWithoutTrfMlt::usage = \"rAdObjDrwWithoutTrfMlt[obj] \"",
	(const char*)0,
	"rAdObjDelDrwAtr::usage = \"rAdObjDelDrwAtr[obj] \"",
	(const char*)0,
	"rAdObjGeoVol::usage = \"rAdObjGeoVol[obj] \"",
	(const char*)0,
	"rAdObjRecMagsAsExtPgns::usage = \"rAdObjRecMagsAsExtPgns[\\\"on|off\\\"",
	"] \"",
	(const char*)0,
	"rAdObjRecMagsAsPolyhdrs::usage = \"rAdObjRecMagsAsPolyhdrs[\\\"on|of",
	"f\\\"] \"",
	(const char*)0,
	"rAdObjExtPgnsAsPolyhdrs::usage = \"rAdObjExtPgnsAsPolyhdrs[\\\"on|of",
	"f\\\"] \"",
	(const char*)0,
	"rAdObjRcgnRecMags::usage = \"rAdObjRcgnRecMags[\\\"on|off\\\"] \"",
	(const char*)0,
	"rAdObjDivByPlns::usage = \"rAdObjDivByPlns[obj,{{{n1x,n1y,n1z},k1",
	",q1},{{n2x,n2y,n2z},k2,q2},...},\\\"lab|loc\\\":\\\"lab\\\"] \"",
	(const char*)0,
	"rAdObjDplFrSym::usage = \"rAdObjDplFrSym[obj] duplicates the geom",
	"etry of the object obj in such a way that a container of new ind",
	"ependent objects appears instead of any symmetry (i.e., transfor",
	"mation with multiplicity more than one) previously applied to th",
	"e object obj.\"",
	(const char*)0,
	"rAdRlxOutIntrcMatr::usage = \"rAdRlxOutIntrcMatr[intrc] \"",
	(const char*)0,
	"rAdRlxOutIntrcVect::usage = \"rAdRlxOutIntrcVect[intrc,\\\"ext|tot|m",
	"ag\\\"] \"",
	(const char*)0,
	"rAdFldMltpolThr::usage = \"rAdFldMltpolThr[{a0,a1,a2,a3}] \"",
	(const char*)0,
	"rAdFldCmpMetNxNyNz::usage = \"rAdFldCmpMetNxNyNz[obj,1|0] \"",
	(const char*)0,
	"rAdUtiDelAll::usage = \"rAdUtiDelAll[] \"",
	(const char*)0,
	"rAdUtiRlxMemAllocMet::usage = \"rAdUtiRlxMemAllocMet[\\\"tot|parts\\\"]",
	" \"",
	(const char*)0,
	"rAdUtiRetInp::usage = \"rAdUtiRetInp[inp,ntimes] \"",
	(const char*)0,
	"rAdUtiStrtProf::usage = \"rAdUtiStrtProf[Flag,NumFuncs,Depth] \"",
	(const char*)0,
	"rAdUtiStopProf::usage = \"rAdUtiStopProf \"",
	(const char*)0,
	(const char*)0
};
#define CARDOF_EVALSTRS 93

static int _definepattern P(( MLINK, char*, char*, int));

static int _doevalstr P(( MLINK, int));

int  _MLDoCallPacket P(( MLINK, struct func[], int));


#if MLPROTOTYPES
int MLInstall( MLINK mlp)
#else
int MLInstall(mlp) MLINK mlp;
#endif
{
	int _res;
	_res = MLConnect(mlp);
	if (_res) _res = _doevalstr( mlp, 0);
	if (_res) _res = _doevalstr( mlp, 1);
	if (_res) _res = _doevalstr( mlp, 2);
	if (_res) _res = _doevalstr( mlp, 3);
	if (_res) _res = _doevalstr( mlp, 4);
	if (_res) _res = _doevalstr( mlp, 5);
	if (_res) _res = _doevalstr( mlp, 6);
	if (_res) _res = _doevalstr( mlp, 7);
	if (_res) _res = _doevalstr( mlp, 8);
	if (_res) _res = _doevalstr( mlp, 9);
	if (_res) _res = _doevalstr( mlp, 10);
	if (_res) _res = _doevalstr( mlp, 11);
	if (_res) _res = _doevalstr( mlp, 12);
	if (_res) _res = _doevalstr( mlp, 13);
	if (_res) _res = _doevalstr( mlp, 14);
	if (_res) _res = _doevalstr( mlp, 15);
	if (_res) _res = _doevalstr( mlp, 16);
	if (_res) _res = _doevalstr( mlp, 17);
	if (_res) _res = _doevalstr( mlp, 18);
	if (_res) _res = _doevalstr( mlp, 19);
	if (_res) _res = _doevalstr( mlp, 20);
	if (_res) _res = _doevalstr( mlp, 21);
	if (_res) _res = _doevalstr( mlp, 22);
	if (_res) _res = _doevalstr( mlp, 23);
	if (_res) _res = _doevalstr( mlp, 24);
	if (_res) _res = _doevalstr( mlp, 25);
	if (_res) _res = _doevalstr( mlp, 26);
	if (_res) _res = _doevalstr( mlp, 27);
	if (_res) _res = _doevalstr( mlp, 28);
	if (_res) _res = _doevalstr( mlp, 29);
	if (_res) _res = _doevalstr( mlp, 30);
	if (_res) _res = _doevalstr( mlp, 31);
	if (_res) _res = _doevalstr( mlp, 32);
	if (_res) _res = _doevalstr( mlp, 33);
	if (_res) _res = _doevalstr( mlp, 34);
	if (_res) _res = _doevalstr( mlp, 35);
	if (_res) _res = _doevalstr( mlp, 36);
	if (_res) _res = _doevalstr( mlp, 37);
	if (_res) _res = _doevalstr( mlp, 38);
	if (_res) _res = _doevalstr( mlp, 39);
	if (_res) _res = _doevalstr( mlp, 40);
	if (_res) _res = _doevalstr( mlp, 41);
	if (_res) _res = _doevalstr( mlp, 42);
	if (_res) _res = _doevalstr( mlp, 43);
	if (_res) _res = _doevalstr( mlp, 44);
	if (_res) _res = _doevalstr( mlp, 45);
	if (_res) _res = _doevalstr( mlp, 46);
	if (_res) _res = _doevalstr( mlp, 47);
	if (_res) _res = _doevalstr( mlp, 48);
	if (_res) _res = _doevalstr( mlp, 49);
	if (_res) _res = _doevalstr( mlp, 50);
	if (_res) _res = _doevalstr( mlp, 51);
	if (_res) _res = _doevalstr( mlp, 52);
	if (_res) _res = _doevalstr( mlp, 53);
	if (_res) _res = _doevalstr( mlp, 54);
	if (_res) _res = _doevalstr( mlp, 55);
	if (_res) _res = _doevalstr( mlp, 56);
	if (_res) _res = _doevalstr( mlp, 57);
	if (_res) _res = _doevalstr( mlp, 58);
	if (_res) _res = _doevalstr( mlp, 59);
	if (_res) _res = _doevalstr( mlp, 60);
	if (_res) _res = _doevalstr( mlp, 61);
	if (_res) _res = _doevalstr( mlp, 62);
	if (_res) _res = _doevalstr( mlp, 63);
	if (_res) _res = _doevalstr( mlp, 64);
	if (_res) _res = _doevalstr( mlp, 65);
	if (_res) _res = _doevalstr( mlp, 66);
	if (_res) _res = _doevalstr( mlp, 67);
	if (_res) _res = _doevalstr( mlp, 68);
	if (_res) _res = _doevalstr( mlp, 69);
	if (_res) _res = _doevalstr( mlp, 70);
	if (_res) _res = _doevalstr( mlp, 71);
	if (_res) _res = _doevalstr( mlp, 72);
	if (_res) _res = _doevalstr( mlp, 73);
	if (_res) _res = _doevalstr( mlp, 74);
	if (_res) _res = _doevalstr( mlp, 75);
	if (_res) _res = _doevalstr( mlp, 76);
	if (_res) _res = _doevalstr( mlp, 77);
	if (_res) _res = _doevalstr( mlp, 78);
	if (_res) _res = _doevalstr( mlp, 79);
	if (_res) _res = _doevalstr( mlp, 80);
	if (_res) _res = _doevalstr( mlp, 81);
	if (_res) _res = _doevalstr( mlp, 82);
	if (_res) _res = _doevalstr( mlp, 83);
	if (_res) _res = _doevalstr( mlp, 84);
	if (_res) _res = _doevalstr( mlp, 85);
	if (_res) _res = _doevalstr( mlp, 86);
	if (_res) _res = _doevalstr( mlp, 87);
	if (_res) _res = _doevalstr( mlp, 88);
	if (_res) _res = _doevalstr( mlp, 89);
	if (_res) _res = _doevalstr( mlp, 90);
	if (_res) _res = _doevalstr( mlp, 91);
	if (_res) _res = _doevalstr( mlp, 92);
	if (_res) _res = _definepattern(mlp, (char *)"radObjRecMag[{xCoordin_,yCoordin_,zCoordin_}, {wxWidth_,wyWidth_,wzWidth_}, Magnetiz_List:{0.,0.,0.}]", (char *)"{ N[xCoordin],N[yCoordin],N[zCoordin], N[wxWidth],N[wyWidth],N[wzWidth], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }", 0);
	if (_res) _res = _definepattern(mlp, (char *)"radObjThckPgn[xCoordin_, lxWidth_, ListOf2dPoints_, Magnetiz_List:{0.,0.,0.}, Orient_String:\"x\"]", (char *)"{ N[xCoordin], N[lxWidth], N[ListOf2dPoints], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]], Orient }", 1);
	if (_res) _res = _definepattern(mlp, (char *)"radObjThckPgn[xCoordin_, lxWidth_, ListOf2dPoints_, Orient_String:\"x\", Magnetiz_List:{0.,0.,0.}]", (char *)"{ N[xCoordin], N[lxWidth], N[ListOf2dPoints], Orient, N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }", 2);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjPgn[zCoordin_, ListOf2dPoints_, Magnetiz_List:{0.,0.,0.}]", (char *)"{ N[zCoordin], N[ListOf2dPoints], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }", 3);
	if (_res) _res = _definepattern(mlp, (char *)"radObjPolyhdr[ListOfPoints_, ListOfListOfIndexes_, OptPar1_:0, OptPar2_:0, OptPar3_:0]", (char *)"{ N[ListOfPoints], Round[ListOfListOfIndexes], N[OptPar1], N[OptPar2], N[OptPar3] }", 4);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjRecMagsAsExtPgns[OnOrOff_String]", (char *)"{ OnOrOff }", 5);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjRecMagsAsPolyhdrs[OnOrOff_String]", (char *)"{ OnOrOff }", 6);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjExtPgnsAsPolyhdrs[OnOrOff_String]", (char *)"{ OnOrOff }", 7);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjRcgnRecMags[OnOrOff_String]", (char *)"{ OnOrOff }", 8);
	if (_res) _res = _definepattern(mlp, (char *)"radObjMltExtPgn[ListOfLayerPolygons_, Magnetiz_List:{0,0,0}]", (char *)"{ N[ListOfLayerPolygons], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }", 9);
	if (_res) _res = _definepattern(mlp, (char *)"radObjMltExtPgnCur[zCoordin_, Orient_, ListOf2dPoints_:0, ExtrusionPath_:0, AvgCurrent_:0, OptPar1_:0]", (char *)"{ N[zCoordin], N[Orient], N[ListOf2dPoints], N[ExtrusionPath], N[AvgCurrent], N[OptPar1] }", 10);
	if (_res) _res = _definepattern(mlp, (char *)"radObjMltExtPgnMag[zCoordin_, Orient_, ListOf2dPoints_:0, ListOfSubdivParam_:0, ExtrusionPath_:0, Magnetiz_:0, OptPar1_:0, OptPar2_:0, OptPar3_:0, OptPar4_:0, OptPar5_:0]", (char *)"{ N[zCoordin], N[Orient], N[ListOf2dPoints], N[ListOfSubdivParam], N[ExtrusionPath], N[Magnetiz], OptPar1, OptPar2, OptPar3, OptPar4, OptPar5  }", 11);
	if (_res) _res = _definepattern(mlp, (char *)"radObjMltExtRtg[ListOfLayerRectangles_, Magnetiz_List:{0,0,0}]", (char *)"{ N[ListOfLayerRectangles], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }", 12);
	if (_res) _res = _definepattern(mlp, (char *)"radObjMltExtTri[xCoordin_, lxWidth_, ListOf2dPoints_, ListOfSubdivParam_, orient_String:\"x\", Magnetiz_List:{0,0,0}, OptPar1_:0, OptPar2_:0, OptPar3_:0, OptPar4_:0]", (char *)"{ N[xCoordin], N[lxWidth], N[ListOf2dPoints], N[ListOfSubdivParam], orient, N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]], OptPar1, OptPar2, OptPar3, OptPar4 }", 13);
	if (_res) _res = _definepattern(mlp, (char *)"radObjArcPgnMag[{xCoordin_,yCoordin_}, Orient_String, ListOf2dPoints_, {phiminWidth_,phimaxWidth_}, SectN_, SymOrNoSym_String:\"nosym\", Magn_List:{0.,0.,0.}]", (char *)"{ N[xCoordin],N[yCoordin], Orient, N[ListOf2dPoints], N[phiminWidth],N[phimaxWidth], Round[SectN], SymOrNoSym, N[Magn[[1]]],N[Magn[[2]]],N[Magn[[3]]] }", 14);
	if (_res) _res = _definepattern(mlp, (char *)"radObjCylMag[{xCoordin_,yCoordin_,zCoordin_}, Rad_, HeightWidth_, SectN_, Orient_String:\"z\", Mv_List:{0,0,0}]", (char *)"{ N[xCoordin],N[yCoordin],N[zCoordin], N[Rad], N[HeightWidth], Round[SectN], Orient, N[Mv[[1]]],N[Mv[[2]]],N[Mv[[3]]] }", 15);
	if (_res) _res = _definepattern(mlp, (char *)"radObjRecCur[{xCoordin_,yCoordin_,zCoordin_}, {wxWidth_,wyWidth_,wzWidth_}, {jxCurrDens_,jyCurrDens_,jzCurrDens_}]", (char *)"{ N[xCoordin],N[yCoordin],N[zCoordin], N[wxWidth],N[wyWidth],N[wzWidth], N[jxCurrDens],N[jyCurrDens],N[jzCurrDens] }", 16);
	if (_res) _res = _definepattern(mlp, (char *)"radObjArcCur[{xCoordin_,yCoordin_,zCoordin_}, {rminWidth_,rmaxWidth_}, {phiminWidth_,phimaxWidth_}, HeightWidth_, SectN_, Jaz_, ManOrAuto_String:\"man\", Orient_String:\"z\"]", (char *)"{ N[xCoordin],N[yCoordin],N[zCoordin], N[rminWidth],N[rmaxWidth], N[phiminWidth],N[phimaxWidth], N[HeightWidth], Round[SectN], N[Jaz], ManOrAuto, Orient }", 17);
	if (_res) _res = _definepattern(mlp, (char *)"radObjRaceTrk[{xCoordin_,yCoordin_,zCoordin_}, {rminWidth_,rmaxWidth_}, {lxWidth_,lyWidth_}, HeightWidth_, SectN_, Jaz_, ManOrAuto_String:\"man\", Orient_String:\"z\"]", (char *)"{ N[xCoordin],N[yCoordin],N[zCoordin], N[rminWidth],N[rmaxWidth], N[lxWidth],N[lyWidth], N[HeightWidth], Round[SectN], N[Jaz], ManOrAuto, Orient }", 18);
	if (_res) _res = _definepattern(mlp, (char *)"radObjFlmCur[ListOfPoints_, TotalCurrent_]", (char *)"{ N[ListOfPoints], N[TotalCurrent] }", 19);
	if (_res) _res = _definepattern(mlp, (char *)"radFldFrcShpRtg[{xCoordin_,yCoordin_,zCoordin_},{wxWidth_,wyWidth_}]", (char *)"{ N[xCoordin],N[yCoordin],N[zCoordin], N[wxWidth],N[wyWidth] }", 20);
	if (_res) _res = _definepattern(mlp, (char *)"radObjCnt[ArrayOfKeys_List]", (char *)"{ ArrayOfKeys }", 21);
	if (_res) _res = _definepattern(mlp, (char *)"radObjAddToCnt[GroupKey_, ArrayOfKeys_List]", (char *)"{ Round[GroupKey], ArrayOfKeys }", 22);
	if (_res) _res = _definepattern(mlp, (char *)"radObjCntStuf[GroupKey_]", (char *)"{ Round[GroupKey] }", 23);
	if (_res) _res = _definepattern(mlp, (char *)"radObjBckg[{bxField_,byField_,bzField_}]", (char *)"{ N[bxField],N[byField],N[bzField] }", 24);
	if (_res) _res = _definepattern(mlp, (char *)"radObjDivMag[ElemKey_, SubdivisionStructure_, OptAdditionalSpecifications_:0, OptPar1_:0, OptPar2_:0, OptPar3_:0]", (char *)"{ Round[ElemKey], N[SubdivisionStructure], N[OptAdditionalSpecifications], OptPar1, OptPar2, OptPar3 }", 25);
	if (_res) _res = _definepattern(mlp, (char *)"radObjCutMag[ElemKey_, PointOnThePlane_, NormalVector_, OptPar1_:0, OptPar2_:0]", (char *)"{ Round[ElemKey], N[PointOnThePlane], N[NormalVector], OptPar1, OptPar2 }", 26);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjDivByPlns[ElemKey_, SubdivisionStructure_, FrameFlag_String:\"lab\"]", (char *)"{ Round[ElemKey], N[SubdivisionStructure], FrameFlag }", 27);
	if (_res) _res = _definepattern(mlp, (char *)"radObjDpl[ElemKey_, OptPar1_:0]", (char *)"{ Round[ElemKey], OptPar1 }", 28);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjDplFrSym[ElemKey_]", (char *)"{ Round[ElemKey] }", 29);
	if (_res) _res = _definepattern(mlp, (char *)"radObjDegFre[ElemKey_]", (char *)"{ Round[ElemKey] }", 30);
	if (_res) _res = _definepattern(mlp, (char *)"radObjM[ElemKey_]", (char *)"{ Round[ElemKey] }", 31);
	if (_res) _res = _definepattern(mlp, (char *)"radObjCenFld[ElemKey_, FldCh_]", (char *)"{ Round[ElemKey], FldCh }", 32);
	if (_res) _res = _definepattern(mlp, (char *)"radObjScaleCur[ElemKey_, scaleCoef_]", (char *)"{ Round[ElemKey], scaleCoef }", 33);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjGeoVol[ElemKey_]", (char *)"{ Round[ElemKey] }", 34);
	if (_res) _res = _definepattern(mlp, (char *)"radObjGeoVol[ElemKey_]", (char *)"{ Round[ElemKey] }", 35);
	if (_res) _res = _definepattern(mlp, (char *)"radObjGeoLim[ElemKey_]", (char *)"{ Round[ElemKey] }", 36);
	if (_res) _res = _definepattern(mlp, (char *)"rAdFldCmpMetNxNyNz[ElemKey_, MethNo_, SubLevel_:0]", (char *)"{ Round[ElemKey], Round[MethNo], Round[SubLevel] }", 37);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjSetLocMgn[ElemKey_, ListOfMagnetiz_]", (char *)"{ Round[ElemKey], ListOfMagnetiz }", 38);
	if (_res) _res = _definepattern(mlp, (char *)"radTrfTrsl[{vxVector_,vyVector_,vzVector_}]", (char *)"{ N[vxVector],N[vyVector],N[vzVector] }", 39);
	if (_res) _res = _definepattern(mlp, (char *)"radTrfRot[{xCoordin_,yCoordin_,zCoordin_}, {vxVector_,vyVector_,vzVector_}, Angle_]", (char *)"{ N[xCoordin],N[yCoordin],N[zCoordin], N[vxVector],N[vyVector],N[vzVector], N[Angle] }", 40);
	if (_res) _res = _definepattern(mlp, (char *)"radTrfPlSym[{xCoordin_,yCoordin_,zCoordin_}, {nxVector_,nyVector_,nzVector_}]", (char *)"{ N[xCoordin],N[yCoordin],N[zCoordin], N[nxVector],N[nyVector],N[nzVector] }", 41);
	if (_res) _res = _definepattern(mlp, (char *)"radTrfInv[]", (char *)"{ }", 42);
	if (_res) _res = _definepattern(mlp, (char *)"radTrfOrnt[g3dElemKey_, TransElemKey_]", (char *)"{ Round[g3dElemKey], Round[TransElemKey] }", 43);
	if (_res) _res = _definepattern(mlp, (char *)"radTrfMlt[g3dElemKey_, TransElemKey_, Multiplicity_]", (char *)"{ Round[g3dElemKey], Round[TransElemKey], Round[Multiplicity] }", 44);
	if (_res) _res = _definepattern(mlp, (char *)"radTrfCmbL[ThisElemKey_, AnotherElemKey_]", (char *)"{ Round[ThisElemKey], Round[AnotherElemKey] }", 45);
	if (_res) _res = _definepattern(mlp, (char *)"radTrfCmbR[ThisElemKey_, AnotherElemKey_]", (char *)"{ Round[ThisElemKey], Round[AnotherElemKey] }", 46);
	if (_res) _res = _definepattern(mlp, (char *)"radMatLin[{ksiparMagneticStuff_,ksiperMagneticStuff_}, {mrxMagneticStuff_,mryMagneticStuff_,mrzMagneticStuff_}]", (char *)"{ N[ksiparMagneticStuff],N[ksiperMagneticStuff], N[mrxMagneticStuff],N[mryMagneticStuff],N[mrzMagneticStuff] }", 47);
	if (_res) _res = _definepattern(mlp, (char *)"radMatLin[{ksiparMagneticStuff_,ksiperMagneticStuff_}, RemanentMagnetizMagnitude_]", (char *)"{ N[ksiparMagneticStuff],N[ksiperMagneticStuff], N[RemanentMagnetizMagnitude] }", 48);
	if (_res) _res = _definepattern(mlp, (char *)"radMatSatIso[{ms1MagneticStuff_,ms2MagneticStuff_,ms3MagneticStuff_},{ks1MagneticStuff_,ks2MagneticStuff_,ks3MagneticStuff_}]", (char *)"{ N[ms1MagneticStuff],N[ms2MagneticStuff],N[ms3MagneticStuff], N[ks1MagneticStuff],N[ks2MagneticStuff],N[ks3MagneticStuff] }", 49);
	if (_res) _res = _definepattern(mlp, (char *)"radMatSatIso[{ks1MagneticStuff_, ms1MagneticStuff_}, {ks2MagneticStuff_, ms2MagneticStuff_}, {ks3MagneticStuff_, ms3MagneticStuff_}]", (char *)"{ N[ks1MagneticStuff],N[ms1MagneticStuff], N[ks2MagneticStuff],N[ms2MagneticStuff], N[ks3MagneticStuff],N[ms3MagneticStuff] }", 50);
	if (_res) _res = _definepattern(mlp, (char *)"radMatSatIso[ListOfPairsHM_]", (char *)"{ N[ListOfPairsHM] }", 51);
	if (_res) _res = _definepattern(mlp, (char *)"radMatSatLam[MagData_, LamCoef_, LamNormal_:0]", (char *)"{ N[MagData], N[LamCoef], N[LamNormal] }", 52);
	if (_res) _res = _definepattern(mlp, (char *)"radMatSatAniso[MagDataPara_, MagDataPerp_]", (char *)"{ N[MagDataPara], N[MagDataPerp] }", 53);
	if (_res) _res = _definepattern(mlp, (char *)"radMatApl[g3dRelaxElemKey_, MaterElemKey_]", (char *)"{ Round[g3dRelaxElemKey], Round[MaterElemKey] }", 54);
	if (_res) _res = _definepattern(mlp, (char *)"radMatMvsH[g3dRelaxOrMaterElemKey_, MagnChar_String, {hxMagneticStuff_,hyMagneticStuff_,hzMagneticStuff_}]", (char *)"{ Round[g3dRelaxOrMaterElemKey], MagnChar, N[hxMagneticStuff],N[hyMagneticStuff],N[hzMagneticStuff] }", 55);
	if (_res) _res = _definepattern(mlp, (char *)"radRlxPre[ElemKey_, SrcElemKey_:0]", (char *)"{ Round[ElemKey], Round[SrcElemKey] }", 56);
	if (_res) _res = _definepattern(mlp, (char *)"rAdRlxOutIntrcMatr[InteractElemKey_]", (char *)"{ Round[InteractElemKey] }", 57);
	if (_res) _res = _definepattern(mlp, (char *)"rAdRlxOutIntrcVect[InteractElemKey_, InteractFieldID_String]", (char *)"{ Round[InteractElemKey], InteractFieldID }", 58);
	if (_res) _res = _definepattern(mlp, (char *)"radRlxMan[InteractElemKey_, MethNo_, IterNum_, RelaxPar_]", (char *)"{ Round[InteractElemKey], Round[MethNo], Round[IterNum], N[RelaxPar] }", 59);
	if (_res) _res = _definepattern(mlp, (char *)"radRlxAuto[InteractElemKey_, Prec_, MaxIterNum_, MethNo_:-1, OptPar1_:-1]", (char *)"{ Round[InteractElemKey], N[Prec], Round[MaxIterNum], MethNo, OptPar1 }", 60);
	if (_res) _res = _definepattern(mlp, (char *)"radRlxUpdSrc[InteractElemKey_]", (char *)"{ Round[InteractElemKey] }", 61);
	if (_res) _res = _definepattern(mlp, (char *)"radSolve[ObjKey_, Prec_, MaxIterNum_, MethNo_:4]", (char *)"{ Round[ObjKey], N[Prec], Round[MaxIterNum], Round[MethNo] }", 62);
	if (_res) _res = _definepattern(mlp, (char *)"radFldCmpCrt[AbsPrecB_, AbsPrecA_, AbsPrecBInt_, AbsPrecFrc_:1., AbsPrecTrjCrd_:-1, AbsPrecTrjAng_:-1]", (char *)"{ N[AbsPrecB], N[AbsPrecA], N[AbsPrecBInt], N[AbsPrecFrc], N[AbsPrecTrjCrd],N[AbsPrecTrjAng] }", 63);
	if (_res) _res = _definepattern(mlp, (char *)"radFldCmpPrc[Prc1_, Prc2_:0, Prc3_:0, Prc4_:0, Prc5_:0, Prc6_:0, Prc7_:0, Prc8_:0, Prc9_:0]", (char *)"{ Prc1, Prc2, Prc3, Prc4, Prc5, Prc6, Prc7, Prc8, Prc9 }", 64);
	if (_res) _res = _definepattern(mlp, (char *)"rAdFldMltpolThr[{a0Threshold_,a1Threshold_,a2Threshold_,a3Threshold_}]", (char *)"{ N[a0Threshold],N[a1Threshold],N[a2Threshold],N[a3Threshold] }", 65);
	if (_res) _res = _definepattern(mlp, (char *)"radFldLst[ElemKey_, FieldChar_String, {xCoordin_,yCoordin_,zCoordin_}, {x2Coordin_,y2Coordin_,z2Coordin_}, np_:1, ShowArgFlag_String:\"noarg\", strtarg_:0.]", (char *)"{ Round[ElemKey], FieldChar, N[xCoordin],N[yCoordin],N[zCoordin], N[x2Coordin],N[y2Coordin],N[z2Coordin], Round[np], ShowArgFlag, N[strtarg]}", 66);
	if (_res) _res = _definepattern(mlp, (char *)"radFldShimSig[ElemKey_, FieldChar_String, {dxCoordin_,dyCoordin_,dzCoordin_}, {xCoordin_,yCoordin_,zCoordin_}, {x2Coordin_,y2Coordin_,z2Coordin_}, np_, {vxCoordin_:0.,vyCoordin_:0.,vzCoordin_:0.}]", (char *)"{ Round[ElemKey], FieldChar, N[dxCoordin],N[dyCoordin],N[dzCoordin], N[xCoordin],N[yCoordin],N[zCoordin], N[x2Coordin],N[y2Coordin],N[z2Coordin], Round[np], N[vxCoordin],N[vyCoordin],N[vzCoordin]}", 67);
	if (_res) _res = _definepattern(mlp, (char *)"radFld[ElemKey_, FieldChar_String, PointsStructure_]", (char *)"{ Round[ElemKey], FieldChar, N[PointsStructure] }", 68);
	if (_res) _res = _definepattern(mlp, (char *)"radFldInt[ElemKey_, CondChar_String, FieldIntChar_String, {x1Coordin_,y1Coordin_,z1Coordin_}, {x2Coordin_,y2Coordin_,z2Coordin_}]", (char *)"{ Round[ElemKey], CondChar, FieldIntChar, N[x1Coordin],N[y1Coordin],N[z1Coordin], N[x2Coordin],N[y2Coordin],N[z2Coordin] }", 69);
	if (_res) _res = _definepattern(mlp, (char *)"radFldFrc[ObjElemKey_, ShapeElemKey_]", (char *)"{ Round[ObjElemKey], Round[ShapeElemKey] }", 70);
	if (_res) _res = _definepattern(mlp, (char *)"radFldEnr[DestObjElemKey_, SourceObjElemKey_, SubdivParam_List:{0,0,0}]", (char *)"{ Round[DestObjElemKey], Round[SourceObjElemKey], Round[SubdivParam[[1]]],Round[SubdivParam[[2]]],Round[SubdivParam[[3]]] }", 71);
	if (_res) _res = _definepattern(mlp, (char *)"radFldEnrFrc[DestObjElemKey_, SourceObjElemKey_, ComponIDChar_String, SubdivParam_List:{0,0,0}]", (char *)"{ Round[DestObjElemKey], Round[SourceObjElemKey], ComponIDChar, Round[SubdivParam[[1]]],Round[SubdivParam[[2]]],Round[SubdivParam[[3]]] }", 72);
	if (_res) _res = _definepattern(mlp, (char *)"radFldEnrTrq[DestObjElemKey_, SourceObjElemKey_, ComponIDChar_String, TorqueCentrPo_List:{0.,0.,0.}, SubdivParam_List:{0,0,0}]", (char *)"{ Round[DestObjElemKey], Round[SourceObjElemKey], ComponIDChar, N[TorqueCentrPo[[1]]],N[TorqueCentrPo[[2]]],N[TorqueCentrPo[[3]]], Round[SubdivParam[[1]]],Round[SubdivParam[[2]]],Round[SubdivParam[[3]]] }", 73);
	if (_res) _res = _definepattern(mlp, (char *)"radFldPtcTrj[ElemKey_,Energy_,{x0Coordin_,dxdy0_,z0Coordin_,dzdy0_},{y0Coordin_,y1Coordin_},np_]", (char *)"{ Round[ElemKey], N[Energy], N[x0Coordin],N[dxdy0],N[z0Coordin],N[dzdy0], N[y0Coordin],N[y1Coordin], Round[np] }", 74);
	if (_res) _res = _definepattern(mlp, (char *)"radFldFocPot[ElemKey_,{xxStart_,yyStart_,zzStart_},{xxFin_,yyFin_,zzFin_},NumPo_]", (char *)"{ Round[ElemKey], N[xxStart],N[yyStart],N[zzStart], N[xxFin],N[yyFin],N[zzFin], Round[NumPo] }", 75);
	if (_res) _res = _definepattern(mlp, (char *)"radFldFocKickPer[ElemKey_,{x0_,y0_,z0_},{nsx_,nsy_,nsz_},per_,nper_,{n1x_,n1y_,n1z_},r1_,np1_,r2_,np2_,comment_String:\"\",{nh_:1,ns_:8,d1_:0,d2_:0},kickUnit_String:\"T2m2\",en_:1,format_String:\"fix\"]", (char *)"{ Round[ElemKey], N[x0],N[y0],N[z0], N[nsx],N[nsy],N[nsz], N[per],N[nper], N[n1x],N[n1y],N[n1z], N[r1],Round[np1],N[r2],Round[np2], comment, Round[nh],Round[ns],N[d1],N[d2], kickUnit, N[en], format }", 76);
	if (_res) _res = _definepattern(mlp, (char *)"radFldFocKick[ElemKey_,P1_List,Nlong_List,ArrayOfLongPos_List,nps_,Ntr1_List,r1_,np1_,r2_,np2_,comment_String:\"\",{d1_:0,d2_:0}]", (char *)"{ Round[ElemKey], N[P1], N[Nlong], N[ArrayOfLongPos],Round[nps], N[Ntr1], N[r1],Round[np1],N[r2],Round[np2], comment, N[d1],N[d2] }", 77);
	if (_res) _res = _definepattern(mlp, (char *)"radFldUnits[]", (char *)"{ }", 78);
	if (_res) _res = _definepattern(mlp, (char *)"radFldLenTol[AbsRandMagnitude_,RelRandMagnitude_,ZeroRandMagnitude_:0]", (char *)"{ N[AbsRandMagnitude], N[RelRandMagnitude], N[ZeroRandMagnitude] }", 79);
	if (_res) _res = _definepattern(mlp, (char *)"radFldLenRndSw[OnOrOff_String]", (char *)"{ OnOrOff }", 80);
	if (_res) _res = _definepattern(mlp, (char *)"radObjDrwAtr[ElemKey_, {rColorAttrib_,gColorAttrib_,bColorAttrib_}, LineThickness_: -1.]", (char *)"{ Round[ElemKey], N[rColorAttrib],N[gColorAttrib],N[bColorAttrib], N[LineThickness] }", 81);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjDelDrwAtr[ElemKey_]", (char *)"{ Round[ElemKey] }", 82);
	if (_res) _res = _definepattern(mlp, (char *)"radObjDrw[ElemKey_, OptPar1_:0]", (char *)"{ Round[ElemKey], OptPar1 }", 83);
	if (_res) _res = _definepattern(mlp, (char *)"rAdObjDrwWithoutTrfMlt[ElemKey_]", (char *)"{ Round[ElemKey] }", 84);
	if (_res) _res = _definepattern(mlp, (char *)"radObjDrwOpenGL[ElemKey_, OptPar1_:0, OptPar2_:0, OptPar3_:0]", (char *)"{ Round[ElemKey], OptPar1, OptPar2, OptPar3 }", 85);
	if (_res) _res = _definepattern(mlp, (char *)"radUtiDel[ElemKey_]", (char *)"{ Round[ElemKey] }", 86);
	if (_res) _res = _definepattern(mlp, (char *)"radUtiDelAll[]", (char *)"{ }", 87);
	if (_res) _res = _definepattern(mlp, (char *)"rAdUtiDelAll[]", (char *)"{ }", 88);
	if (_res) _res = _definepattern(mlp, (char *)"radUtiDmp[ElemKey_, OutFormat_String:\"asc\"]", (char *)"{ Round[ElemKey], OutFormat }", 89);
	if (_res) _res = _definepattern(mlp, (char *)"radUtiDmpPrs[bstr_]", (char *)"{ bstr }", 90);
	if (_res) _res = _definepattern(mlp, (char *)"radUtiVer[]", (char *)"{ }", 91);
	if (_res) _res = _definepattern(mlp, (char *)"rAdUtiRetInp[Input_, NumTimes_]", (char *)"{ N[Input], Round[NumTimes] }", 92);
	if (_res) _res = _definepattern(mlp, (char *)"rAdUtiRlxMemAllocMet[TotOrParts_String]", (char *)"{ TotOrParts }", 93);
	if (_res) _res = _definepattern(mlp, (char *)"rAdUtiStrtProf[Flag_, NumFunc_, Depth_]", (char *)"{ Round[Flag], Round[NumFunc], Round[Depth] }", 94);
	if (_res) _res = _definepattern(mlp, (char *)"rAdUtiStopProf", (char *)"{ }", 95);
	if (_res) _res = _definepattern(mlp, (char *)"radUtiIntrptTim[tTimeQuanta_:1.0]", (char *)"{ N[tTimeQuanta] }", 96);
	if (_res) _res = MLPutSymbol( mlp, "End");
	if (_res) _res = MLFlush( mlp);
	return _res;
} /* MLInstall */


#if MLPROTOTYPES
int MLDoCallPacket( MLINK mlp)
#else
int MLDoCallPacket( mlp) MLINK mlp;
#endif
{
	return _MLDoCallPacket( mlp, _tramps, 97);
} /* MLDoCallPacket */

/******************************* begin trailer ********************************/

#ifndef EVALSTRS_AS_BYTESTRINGS
#	define EVALSTRS_AS_BYTESTRINGS 1
#endif

#if CARDOF_EVALSTRS
static int  _doevalstr( MLINK mlp, int n)
{
	long bytesleft, charsleft, bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
	long charsnow;
#endif
	char **s, **p;
	char *t;

	s = (char **)evalstrs;
	while( n-- > 0){
		if( *s == 0) break;
		while( *s++ != 0){}
	}
	if( *s == 0) return 0;
	bytesleft = 0;
	charsleft = 0;
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = (long)(t - *p);
		bytesleft += bytesnow;
		charsleft += bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
		t = *p;
		charsleft -= MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
#endif
		++p;
	}


	MLPutNext( mlp, MLTKSTR);
#if EVALSTRS_AS_BYTESTRINGS
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = (long)(t - *p);
		bytesleft -= bytesnow;
		MLPut8BitCharacters( mlp, bytesleft, (unsigned char*)*p, bytesnow);
		++p;
	}
#else
	MLPut7BitCount( mlp, (long_st)charsleft, (long_st)bytesleft);

	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		t = *p;
		charsnow = bytesnow - MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
		charsleft -= charsnow;
		MLPut7BitCharacters(  mlp, charsleft, *p, bytesnow, charsnow);
		++p;
	}
#endif
	return MLError( mlp) == MLEOK;
}
#endif /* CARDOF_EVALSTRS */


static int  _definepattern( MLINK mlp, char *patt, char *args, int func_n)
{
	MLPutFunction( mlp, "DefineExternal", (long)3);
	  MLPutString( mlp, patt);
	  MLPutString( mlp, args);
	  MLPutInteger( mlp, func_n);
	return !MLError(mlp);
} /* _definepattern */


int _MLDoCallPacket( MLINK mlp, struct func functable[], int nfuncs)
{
	long len;
	int n, res = 0;
	struct func* funcp;

	if( ! MLGetInteger( mlp, &n) ||  n < 0 ||  n >= nfuncs) goto L0;
	funcp = &functable[n];

	if( funcp->f_nargs >= 0
	&& ( ! MLCheckFunction(mlp, "List", &len)
	     || ( !funcp->manual && (len != funcp->f_nargs))
	     || (  funcp->manual && (len <  funcp->f_nargs))
	   )
	) goto L0;

	stdlink = mlp;
	res = (*funcp->f_func)( mlp);

L0:	if( res == 0)
		res = MLClearError( mlp) && MLPutSymbol( mlp, "$Failed");
	return res && MLEndPacket( mlp) && MLNewPacket( mlp);
} /* _MLDoCallPacket */


mlapi_packet MLAnswer( MLINK mlp)
{
	mlapi_packet pkt = 0;

	while( !MLDone && !MLError(mlp)
	&& (pkt = MLNextPacket(mlp), pkt) && pkt == CALLPKT){
		MLAbort = 0;
		if( !MLDoCallPacket(mlp)) pkt = 0;
	}
	MLAbort = 0;
	return pkt;
}



/*
	Module[ { me = $ParentLink},
		$ParentLink = contents of RESUMEPKT;
		Message[ MessageName[$ParentLink, "notfe"], me];
		me]
*/

static int refuse_to_be_a_frontend( MLINK mlp)
{
	int pkt;

	MLPutFunction( mlp, "EvaluatePacket", 1);
	  MLPutFunction( mlp, "Module", 2);
	    MLPutFunction( mlp, "List", 1);
		  MLPutFunction( mlp, "Set", 2);
		    MLPutSymbol( mlp, "me");
	        MLPutSymbol( mlp, "$ParentLink");
	  MLPutFunction( mlp, "CompoundExpression", 3);
	    MLPutFunction( mlp, "Set", 2);
	      MLPutSymbol( mlp, "$ParentLink");
	      MLTransferExpression( mlp, mlp);
	    MLPutFunction( mlp, "Message", 2);
	      MLPutFunction( mlp, "MessageName", 2);
	        MLPutSymbol( mlp, "$ParentLink");
	        MLPutString( mlp, "notfe");
	      MLPutSymbol( mlp, "me");
	    MLPutSymbol( mlp, "me");
	MLEndPacket( mlp);

	while( (pkt = MLNextPacket( mlp), pkt) && pkt != SUSPENDPKT)
		MLNewPacket( mlp);
	MLNewPacket( mlp);
	return MLError( mlp) == MLEOK;
}


#if MLINTERFACE >= 3
int MLEvaluate( MLINK mlp, char *s)
#else
int MLEvaluate( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
{
	if( MLAbort) return 0;
	return MLPutFunction( mlp, "EvaluatePacket", 1L)
		&& MLPutFunction( mlp, "ToExpression", 1L)
		&& MLPutString( mlp, s)
		&& MLEndPacket( mlp);
}


#if MLINTERFACE >= 3
int MLEvaluateString( MLINK mlp, char *s)
#else
int MLEvaluateString( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
{
	int pkt;
	if( MLAbort) return 0;
	if( MLEvaluate( mlp, s)){
		while( (pkt = MLAnswer( mlp), pkt) && pkt != RETURNPKT)
			MLNewPacket( mlp);
		MLNewPacket( mlp);
	}
	return MLError( mlp) == MLEOK;
} /* MLEvaluateString */


#if __BORLANDC__
#pragma argsused
#endif

#if MLINTERFACE >= 3
MLMDEFN( void, MLDefaultHandler, ( MLINK mlp, int message, int n))
#else
MLMDEFN( void, MLDefaultHandler, ( MLINK mlp, unsigned long message, unsigned long n))
#endif /* MLINTERFACE >= 3 */
{
#if !__BORLANDC__
	mlp = (MLINK)0; /* suppress unused warning */
	n = 0;          /* suppress unused warning */
#endif

	switch (message){
	case MLTerminateMessage:
		MLDone = 1;
	case MLInterruptMessage:
	case MLAbortMessage:
		MLAbort = 1;
	default:
		return;
	}
}



#if MLINTERFACE >= 3
static int _MLMain( char **argv, char **argv_end, char *commandline)
#else
static int _MLMain( charpp_ct argv, charpp_ct argv_end, charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
{
	MLINK mlp;
#if MLINTERFACE >= 3
	int err;
#else
	long err;
#endif /* MLINTERFACE >= 3 */

	if( !stdenv)
		stdenv = MLInitialize( (MLParametersPointer)0);
	if( stdenv == (MLEnvironment)0) goto R0;

	if( !stdyielder)
#if MLINTERFACE >= 3
		stdyielder = (MLYieldFunctionObject)MLDefaultYielder;
#else
		stdyielder = MLCreateYieldFunction( stdenv,
			NewMLYielderProc( MLDefaultYielder), 0);
#endif /* MLINTERFACE >= 3 */


#if MLINTERFACE >= 3
	if( !stdhandler)
		stdhandler = (MLMessageHandlerObject)MLDefaultHandler;
#else
	if( !stdhandler)
		stdhandler = MLCreateMessageHandler( stdenv,
			NewMLHandlerProc( MLDefaultHandler), 0);
#endif /* MLINTERFACE >= 3 */


	mlp = commandline
		? MLOpenString( stdenv, commandline, &err)
#if MLINTERFACE >= 3
		: MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
#else
		: MLOpenArgv( stdenv, argv, argv_end, &err);
#endif
	if( mlp == (MLINK)0){
		MLAlert( stdenv, MLErrorString( stdenv, err));
		goto R1;
	}

	if( MLIconWindow){
#define TEXTBUFLEN 64
		TCHAR textbuf[TEXTBUFLEN];
		PTCHAR tmlname;
		const char *mlname;
		size_t namelen, i;
		int len;
		len = GetWindowText(MLIconWindow, textbuf, 62 );
		mlname = MLName(mlp);
		namelen = strlen(mlname);
		tmlname = (PTCHAR)malloc((namelen + 1)*sizeof(TCHAR));
		if(tmlname == NULL) goto R2;

		for(i = 0; i < namelen; i++){
			tmlname[i] = mlname[i];
		}
		tmlname[namelen] = '\0';
		
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
		_tcscat_s( textbuf + len, TEXTBUFLEN - len, __TEXT("("));
		_tcsncpy_s(textbuf + len + 1, TEXTBUFLEN - len - 1, tmlname, TEXTBUFLEN - len - 3);
		textbuf[TEXTBUFLEN - 2] = '\0';
		_tcscat_s(textbuf, TEXTBUFLEN, __TEXT(")"));
#else
		_tcscat( textbuf + len, __TEXT("("));
		_tcsncpy( textbuf + len + 1, tmlname, TEXTBUFLEN - len - 3);
		textbuf[TEXTBUFLEN - 2] = '\0';
		_tcscat( textbuf, __TEXT(")"));
#endif
		textbuf[len + namelen + 2] = '\0';
		free(tmlname);
		SetWindowText( MLIconWindow, textbuf);
	}

	if( MLInstance){
		if( stdyielder) MLSetYieldFunction( mlp, stdyielder);
		if( stdhandler) MLSetMessageHandler( mlp, stdhandler);
	}

	if( MLInstall( mlp))
		while( MLAnswer( mlp) == RESUMEPKT){
			if( ! refuse_to_be_a_frontend( mlp)) break;
		}

R2:	MLClose( mlp);
R1:	MLDeinitialize( stdenv);
	stdenv = (MLEnvironment)0;
R0:	return !MLDone;
} /* _MLMain */


#if MLINTERFACE >= 3
int MLMainString( char *commandline)
#else
int MLMainString( charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
{
#if MLINTERFACE >= 3
	return _MLMain( (char **)0, (char **)0, commandline);
#else
	return _MLMain( (charpp_ct)0, (charpp_ct)0, commandline);
#endif /* MLINTERFACE >= 3 */
}

int MLMainArgv( char** argv, char** argv_end) /* note not FAR pointers */
{   
	static char FAR * far_argv[128];
	int count = 0;
	
	while(argv < argv_end)
		far_argv[count++] = *argv++;
		 
#if MLINTERFACE >= 3
	return _MLMain( far_argv, far_argv + count, (char *)0);
#else
	return _MLMain( far_argv, far_argv + count, (charp_ct)0);
#endif /* MLINTERFACE >= 3 */

}

#if MLINTERFACE >= 3
int MLMain( int argc, char **argv)
#else
int MLMain( int argc, charpp_ct argv)
#endif /* MLINTERFACE >= 3 */
{
#if MLINTERFACE >= 3
 	return _MLMain( argv, argv + argc, (char *)0);
#else
 	return _MLMain( argv, argv + argc, (charp_ct)0);
#endif /* MLINTERFACE >= 3 */
}
 
