
#ifndef __GLUTLINK_H
#define __GLUTLINK_H

#ifndef __GLUTVIEW_H
#include "glutview.h"
#endif

//-------------------------------------------------------------------------
// This is dedicated for run-time linking of the GlutViewer DLL
//-------------------------------------------------------------------------

#ifdef _WINDOWS

#include <windows.h>		// required for all Windows applications

//typedef int (CALL* CPViewSceneOpenGL)(double**, int*, int**, int**, float**, int*);
//typedef int (CALL* CPViewPolygonsOpenGL)(double*, int, int*, int*, float*, int);
typedef int (CALL* CPViewScene3D)(double**, int*, int**, int**, float**, int*, char*, char, char*);
typedef int (CALL* CPViewPolygons3D)(double*, int, int*, int*, float*, int, char*, char, char*);

struct COpenGLViewer {
	HMODULE hInstDLL;
	//CPViewSceneOpenGL pViewSceneOpenGL;
	//CPViewPolygonsOpenGL pViewPolygonsOpenGL;
	CPViewScene3D pViewScene3D;
	CPViewPolygons3D pViewPolygons3D;
	char mLibType; // 's' or 'd' (static or dynamic)

	COpenGLViewer(char s_or_d)
	{
		mLibType = 's';
		if(s_or_d == 'd') //"dynamic"
		{
            ZeroFuncPointers();
            AssignFuncPointers();
			mLibType = 'd';
		}
	}
	COpenGLViewer() 
	{
		mLibType = 's';
	}

	~COpenGLViewer()
	{
		if(hInstDLL != NULL) FreeLibrary(hInstDLL); 
	}

	void ZeroFuncPointers()
	{
		hInstDLL = NULL;
		pViewScene3D = NULL;
		pViewPolygons3D = NULL;
	}

	void AssignFuncPointers()
	{
		if(hInstDLL == NULL) hInstDLL = LoadLibrary("GlutViewer");
		if(hInstDLL == NULL) return;

		if(pViewScene3D == NULL) pViewScene3D = (CPViewScene3D)SetFuncPointer("ViewScene3D", "_ViewScene3D@24");
		if(pViewPolygons3D == NULL) pViewPolygons3D = (CPViewPolygons3D)SetFuncPointer("ViewPolygons3D", "_ViewPolygons3D@24");
	}

	FARPROC SetFuncPointer(char* FuncNameA, char* FuncNameB)
	{
		if((hInstDLL == NULL) || ((FuncNameA == NULL) && (FuncNameB == NULL))) return NULL;

		FARPROC pOut = NULL;
		if(FuncNameA != NULL) pOut = GetProcAddress(hInstDLL, FuncNameA);
		if(pOut != NULL) return pOut;
		if(FuncNameB != NULL) pOut = GetProcAddress(hInstDLL, FuncNameB);
		if(pOut != NULL) return pOut;
		FreeLibrary(hInstDLL); 
		return NULL;
	}

	//int viewSceneOpenGL(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr)
	int viewScene3D(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr, char* WinTitle, char ViewerStartMode, char* ErrWarnText)
	{
		if(mLibType == 'd')
		{
            if(pViewScene3D == NULL) AssignFuncPointers();
            if(pViewScene3D == NULL) return -1;
            return pViewScene3D(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr, WinTitle, ViewerStartMode, ErrWarnText);
		}
		//return ViewSceneOpenGL(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr);
		return ViewScene3D(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr, WinTitle, ViewerStartMode, ErrWarnText);
	}
	//int viewPolygonsOpenGL(double* VertCoord, int Nv, int* VertInd, int* LenInd, float* Colors, int Npg)
	int viewPolygons3D(double* VertCoord, int Nv, int* VertInd, int* LenInd, float* Colors, int Npg, char* WinTitle, char ViewerStartMode, char* ErrWarnText)
	{
		if(mLibType == 'd')
		{
            if(pViewPolygons3D == NULL) AssignFuncPointers();
            if(pViewPolygons3D == NULL) return -1;
            return pViewPolygons3D(VertCoord, Nv, VertInd, LenInd, Colors, Npg, WinTitle, ViewerStartMode, ErrWarnText);
		}
        //return ViewPolygonsOpenGL(VertCoord, Nv, VertInd, LenInd, Colors, Npg);
        return ViewPolygons3D(VertCoord, Nv, VertInd, LenInd, Colors, Npg, WinTitle, ViewerStartMode, ErrWarnText);
	}
};

#else

struct COpenGLViewer {
	//int viewSceneOpenGL(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr)
	int viewScene3D(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr, char* ErrWarnText)
	{
		//return viewSceneOpenGL(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr);
		return ViewScene3D(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr, ErrWarnText);
	}
};

#endif

//-------------------------------------------------------------------------

#endif


