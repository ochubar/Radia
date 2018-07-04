
#include "glutview.h"

#ifndef __VIEWER3D_H
#include "viewer3d.h"
#endif
#ifndef __VIEWERROR_H
#include "viewerr.h"
#endif
#ifndef __PLOT2D_H
#include "plot2d.h"
#endif

//-------------------------------------------------------------------------

int CALL ViewGetErrorText(int ErrNo, char* ErrText)
{
	if(ErrText == 0) return 0;
	CErrWarn::CopyErrorText(ErrNo, ErrText);
	return 0;
}

//-------------------------------------------------------------------------

int CALL ViewGetWarningText(int WarnNo, char* WarnText)
{
	if(WarnText == 0) return 0;
	CErrWarn::CopyWarningText(WarnNo, WarnText);
	return 0;
}

//-------------------------------------------------------------------------

int CALL ViewGetAllWarningText(char* WarnText)
{
	CErrWarn::CopyAllWarningText(WarnText);
	return 0;
}

//-------------------------------------------------------------------------

int CALL ViewPolygons3D(double* VertCoord, int Nv, int* PgVertInd, int* PgLen, float* PgColors, int Npg, char* WinTitle, char StartMode, char* ErrWarnText)
{
	try
	{
		if((VertCoord == 0) || (Nv == 0) || (PgVertInd == 0) || (PgLen == 0)) throw 1;

		double* VertCoordArr[] = {VertCoord, 0};
		int NvArr[] = {Nv, 0};
		int* VertIndArr[] = {PgVertInd, 0};
		int* LenArr[] = {PgLen, 0};
		float* ColorsArr[] = {PgColors, 0};
		int NpArr[] = {Npg, 0};

		CHWinCont hWinCont(new CWinContViewer3D(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr, WinTitle));
		CSimpleGraph::StartNewWinApp(hWinCont, StartMode);
	}
	catch(int ErrNo) 
	{
		if(ErrWarnText != 0) CErrWarn::CopyErrorText(ErrNo, ErrWarnText);
		return ErrNo;
	}
	if(ErrWarnText != 0) CErrWarn::CopyAllWarningText(ErrWarnText);
	return 0;
}

//-------------------------------------------------------------------------

int CALL ViewScene3D(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr, char* WinTitle, char StartMode, char* ErrWarnText)
{
	try
	{
		if((VertCoordArr == 0) || (NvArr == 0) || (VertIndArr == 0) || (LenArr == 0) || (NpArr == 0)) throw 1;

		CHWinCont hWinCont(new CWinContViewer3D(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr, WinTitle));
		CSimpleGraph::StartNewWinApp(hWinCont, StartMode);
	}
	catch(int ErrNo) 
	{
		if(ErrWarnText != 0) CErrWarn::CopyErrorText(ErrNo, ErrWarnText);
		return ErrNo;
	}
	if(ErrWarnText != 0) CErrWarn::CopyAllWarningText(ErrWarnText);
	return 0;
}

//-------------------------------------------------------------------------

int CALL ViewPlot2D(double** FuncValues, double** ArgValues, double* ArgStart, double* ArgStep, long* Size, double** CurveOptions, int CurveNumber, char** Units, char** Labels, double* GraphOptions, char* WinTitle, char StartMode, char* ErrWarnText)
{
	try
	{
		if((FuncValues == 0) || (Size == 0) || (CurveNumber <= 0)) throw 1;

		CHWinCont hWinCont(new CWinContPlot2D(FuncValues, ArgValues, ArgStart, ArgStep, Size, CurveOptions, CurveNumber, Units, Labels, GraphOptions, WinTitle));
		CSimpleGraph::StartNewWinApp(hWinCont, StartMode);
	}
	catch(int ErrNo) 
	{
		if(ErrWarnText != 0) CErrWarn::CopyErrorText(ErrNo, ErrWarnText);
		return ErrNo;
	}
	if(ErrWarnText != 0) CErrWarn::CopyAllWarningText(ErrWarnText);
	return 0;
}

//-------------------------------------------------------------------------