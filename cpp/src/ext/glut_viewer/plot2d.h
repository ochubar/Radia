
#ifndef __PLOT2D_H
#define __PLOT2D_H

#ifndef __SIMPLEGRAPH_H
#include "simplegraph.h"
#endif
//#ifndef __GMTRANS_H
//#include "gmtrans.h"
//#endif
//#include <GL/glut.h>

//-------------------------------------------------------------------------
// The base class of the OpenGL-GLUT 3D Viewer and Grapher
//-------------------------------------------------------------------------

class CWinContPlot2D : public CWinCont {

	//gmTrans m_RotInit, m_TrfAux1, m_TrfAux2, m_RotCur, m_TrfCur; //Initial rotation - permutation of axes
	//static TVector3d m_Vz, m_ZeroV;

	//float m_Aspect;
	//float m_DepthNear, m_DepthFar;
	//float m_Depth, m_DepthMotionQuanta;
	//int m_DownX, m_DownY;
    //bool m_LeftButton, m_MiddleButton;
    //int m_DisplayMode;
    //int m_MotionType;

	//bool m_DrawFaceNorms, m_Antialias, m_EnvMap;
	//int m_MenuDistance, m_MenuMotion, m_MenuMain;

public:

	CWinContPlot2D(double** FuncValues, double** ArgValues, double* ArgStart, double* ArgStep, long* Size, double** CurveOptions, int CurveNumber, char** Units, char** Labels, double* GraphOptions, char* WinTitle);

	bool CheckHowArgIsDefined(double** ArgValues, double* ArgStart, double* ArgStep, int CurveNumber, bool* ArrArgStartStepDefined, bool* ArrArgTableDefined);
    void FindCurveAbsMinAndMax(double** FuncValues, double** ArgValues, double* ArgStart, double* ArgStep, long* Size, bool* ArrArgTableDefined, bool* ArrArgStartStepDefined, int AmOfCurves, double& MinArg, double& MaxArg, double& MinVal, double& MaxVal);
    void FindOneCurveMinMax(double* FuncValues, double* ArgValues, double ArgStart, double ArgStep, long Size, double& MinArg, double& MaxArg, double& MinVal, double& MaxVal);

/**
	float m_LightPosition[4];

	static const int sm_SMALL, sm_MEDIUM, sm_LARGE, sm_XLARGE;
	static const int sm_WIREFRAME, sm_HIDDENLINE, sm_FLATSHADED, sm_SMOOTHSHADED, sm_TEXTURED;
	static const int sm_FULLSCREEN, sm_FACENORMALS, sm_ANTIALIAS, sm_ENVMAP;
	static const int sm_ROTATION3D, sm_ROTATION2D, sm_TRANSLATION_DEPTH;

	CWinContViewer3D() 
	{
		InitViewingParams();
		SetupRotInit();
		SetupTrfCur(m_RotInit);
	}
	~CWinContViewer3D() {}

	static void DetermineAverageCoordAndLimits(double* VertCoord, int Nv, double* LnVertCoord, int LnNv, double* CenVert, double* ObjLimits, double& MaxSize);
    void CalculatePolygonNormals();

	void ProcDisplayContent(); // virtual
	void ProcReshapeContent(int Width, int Height); // virtual
	void ProcMouseContent(int button, int state, int x, int y); // virtual
	void ProcMotionContent(int x, int y); // virtual
	void ProcKeyboardContent(unsigned char ch, int x, int y); // virtual

	static void DrawArbPolygObjects(GLenum Mode, double* VertCoord, int Nv, int* VertInd, int* Lengths, float* Colors, int Npg, float* NormCoord, bool& NormalsWereDefined);
	//static void DrawArbPolygObjects(GLenum Mode, double* VertCoord, int Nv, int* VertInd, int* Lengths, float* Colors, int Npg, float* NormCoord, bool NormalsWereDefined);

	void SetViewingParamsGL(); // virtual
	void SetupMenuGLUT(); // virtual

	void DoSetSizeContent(int SizeCase); // virtual
	void DoSetMotionContent(int MotionCase); // virtual
	//void DoSetMainContent(int MainCase); // virtual

	//void SetViewSize(int SizeCase);
    void SetDisplayMode(int InMode);
    void SetOtherParams(int InPar);

	void InitViewingParams()
	{
		m_Aspect = 5.0/4.0;
		m_DepthNear = 15.0; m_DepthFar = 100.0;
		m_Depth = 0; m_DepthMotionQuanta = 0;
		//m_Phi = 270.0, m_Theta = 90;
		m_DisplayMode = sm_WIREFRAME;
		m_MotionType = sm_ROTATION3D;

		m_LightPosition[0] = 0; m_LightPosition[1] = 0; m_LightPosition[2] = 1; m_LightPosition[3] = 1; 

		m_DrawFaceNorms = m_Antialias = m_EnvMap = false;

		m_DownX = m_DownY = 0;
		m_LeftButton = m_MiddleButton = false;
	}
	void SetupRotInit()
	{
		const double PI = 3.1415926;
		const double Angle = 2*PI/3.;
		TVector3d ZeroPoiOnAxVect(0, 0, 0), AxVect(-1, -1, -1);
		m_RotInit.SetupRotation(ZeroPoiOnAxVect, AxVect, Angle);
		m_RotCur = m_RotInit;
	}
	void SetupTrfCur(gmTrans& RotCur)
	{
        TVector3d TrV1(-m_CenVert[0], -m_CenVert[1], -m_CenVert[2]);
        m_TrfAux1.SetupTranslation(TrV1);
		TrProduct(&RotCur, &m_TrfAux1, m_TrfAux2);
        TVector3d TrV2(0, 0, -m_Depth);
        m_TrfAux1.SetupTranslation(TrV2);
		TrProduct(&m_TrfAux1, &m_TrfAux2, m_TrfCur);
	}
	void SetupViewMatrArrGL(GLfloat *ArrMatrGL)
	{
		if(ArrMatrGL == 0) return;
		TMatrix3d TrfCurM = m_TrfCur.OutMatrix();
		TVector3d TrfCurV = m_TrfCur.OutVector();
		TVector3d &MStr0 = TrfCurM.Str0, &MStr1 = TrfCurM.Str1, &MStr2 = TrfCurM.Str2;

        ArrMatrGL[0] = MStr0.x; ArrMatrGL[4] = MStr0.y; ArrMatrGL[8] = MStr0.z; ArrMatrGL[12] = TrfCurV.x;
        ArrMatrGL[1] = MStr1.x; ArrMatrGL[5] = MStr1.y; ArrMatrGL[9] = MStr1.z; ArrMatrGL[13] = TrfCurV.y;
        ArrMatrGL[2] = MStr2.x; ArrMatrGL[6] = MStr2.y; ArrMatrGL[10] = MStr2.z; ArrMatrGL[14] = TrfCurV.z;
        ArrMatrGL[3] = 0; ArrMatrGL[7] = 0; ArrMatrGL[11] = 0; ArrMatrGL[15] = 1;
	}
**/
};

//-------------------------------------------------------------------------

#endif
