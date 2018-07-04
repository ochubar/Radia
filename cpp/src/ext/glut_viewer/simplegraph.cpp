
#ifndef __SIMPLEGRAPH_H
#include "simplegraph.h"
#endif

//#ifdef WIN32
//#include <process.h>
//#include <windows.h>
//#endif

//#ifdef GLUT_BUILDING_LIB
//#undef GLUT_BUILDING_LIB //to allow compilimg and linking with glut.h
//#endif
//#ifdef _DLL
//#undef _DLL
//#endif

//#ifndef __glut_h__
//#include <GL/glut.h>
//#endif

#ifndef GLAPI
#define GLAPI
#endif
#ifndef GLAPIENTRY
#define GLAPIENTRY __stdcall
#endif

//#include "GL/osmesa.h" //OC090105

#ifndef PNG_H
#include "png.h"
#endif
//#ifndef __IJL_H__
//#include "ijl.h"
//#endif

#include <string.h>
#include <math.h>
#include <stdio.h>

//-------------------------------------------------------------------------

const int CWinCont::sm_EXIT = 50;
const int CWinCont::sm_SAVE = 51;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//enum {WIREFRAME, HIDDENLINE, FLATSHADED, SMOOTHSHADED, TEXTURED};
//enum {FULLSCREEN, FACENORMALS, ANTIALIAS, ENVMAP};
//enum {SMALL, MEDIUM, LARGE, XLARGE};
//
//void MyTest_setSize(int value)
//{
//		CHWinCont hCurCont;
//		CSimpleGraph::GetCurHWinCont(hCurCont);
//		CWinCont *pc = hCurCont.rep;
//
//
//	double xDim = pc->m_ObjLimits[1] - pc->m_ObjLimits[0];
//	double yDim = pc->m_ObjLimits[3] - pc->m_ObjLimits[2];
//	double zDim = pc->m_ObjLimits[5] - pc->m_ObjLimits[4];
//
//	double TransSize = yDim;
//	if(TransSize < zDim) TransSize = zDim;
//
//	double LocGrid;
//    switch(value) 
//    {
//        case SMALL : LocGrid = 2*TransSize/4; break;
//        case MEDIUM: LocGrid = 2*TransSize/2; break;
//        case LARGE : LocGrid = 2*TransSize/1.5; break;
//        case XLARGE : LocGrid = 2*TransSize; break;
//    }
//
//	LocGrid += 0.5*xDim;
//
//    pc->m_DepthNear= LocGrid/20.0;
//    pc->m_DepthFar= LocGrid*30.0;
//
//    pc->m_Depth = (5.0/4.0)*LocGrid;
//	pc->m_DepthMotionQuanta = 0.05*pc->m_Depth;
//
//	pc->m_Phi=270.0;
//	pc->m_Theta=90; //45.0;
//
//    glutPostRedisplay();
//}
//
//void MyTest_setDisplay(int value)
//{
//			CHWinCont hCurCont;
//			CSimpleGraph::GetCurHWinCont(hCurCont);
//			CWinCont *pc = hCurCont.rep;
//
//    pc->m_DisplayMode = value;
//    switch(value) 
//    {
//        case WIREFRAME   : 
//            glShadeModel(GL_FLAT); 
//            glDisable(GL_LIGHTING);
//            break;
//        case HIDDENLINE: 
//            glShadeModel(GL_FLAT); 
//            glDisable(GL_LIGHTING);
//            break;
//        case FLATSHADED  : 
//            glShadeModel(GL_FLAT); 
//            glEnable(GL_LIGHTING);
//            break;
//        case SMOOTHSHADED: 
//            glShadeModel(GL_SMOOTH); 
//            glEnable(GL_LIGHTING);
//            break;
//        case TEXTURED: 
//            glShadeModel(GL_SMOOTH); 
//            glEnable(GL_LIGHTING);
//            break;
//    }
//    glutPostRedisplay();
//}
//
//void MyTest_setOther(int value)
//{
//			CHWinCont hCurCont;
//			CSimpleGraph::GetCurHWinCont(hCurCont);
//			CWinCont *pc = hCurCont.rep;
//
//    switch (value)
//    {
//        case FULLSCREEN: 
//            glutFullScreen();
//            break;
//        case FACENORMALS: 
//            pc->m_DrawFaceNorms = !pc->m_DrawFaceNorms;
//            break;
//        case ANTIALIAS: 
//            pc->m_Antialias = !pc->m_Antialias;
//            if (pc->m_Antialias)
//            {
//                glEnable(GL_BLEND);
//                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//                glEnable(GL_LINE_SMOOTH);
//                glLineWidth(1.5);
//            }
//            else
//            {
//                glDisable(GL_BLEND);
//                glDisable(GL_LINE_SMOOTH);
//                glLineWidth(1.0);
//            }
//            break;
//        case ENVMAP: 
//            pc->m_EnvMap = !pc->m_EnvMap;
//            if (pc->m_EnvMap)
//            {
//                //glBindTexture(GL_TEXTURE_2D, texId2);
//                glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
//                glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
//                glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
//                glEnable(GL_TEXTURE_GEN_S);
//                glEnable(GL_TEXTURE_GEN_T);
//            }
//            else
//            {
//                //glBindTexture(GL_TEXTURE_2D, texId1);
//                glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
//                glDisable(GL_TEXTURE_GEN_S);
//                glDisable(GL_TEXTURE_GEN_T);
//            }
//            break;
//    }
//    glutPostRedisplay();
//}

//-------------------------------------------------------------------------

//bool MyTest_CheckIfVectorsAreCollinear(double xV1, double yV1, double zV1, double xV2, double yV2, double zV2)
//{
//	const double RelPrec = 1.e-5;
//
//	double V1scV2 = xV1*xV2 + yV1*yV2 + zV1*zV2;
//	double V1scV2e2 = V1scV2*V1scV2;
//
//	double V1e2 = xV1*xV1 + yV1*yV1 + zV1*zV1;
//	double V2e2 = xV2*xV2 + yV2*yV2 + zV2*zV2;
//	double V1e2V2e2 = V1e2*V2e2;
//
//	if(fabs(fabs(V1scV2e2) - fabs(V1e2V2e2)) < RelPrec*V1e2V2e2) return true;
//	return false;
//}
//
//void MyTest_DrawArbPolygObjects(GLenum Mode, double* VertCoord, int Nv, int* VertInd, int* Lengths, float* Colors, int Npg, float* NormCoord, bool& NormalsWereDefined)
//{
//	bool ColorsShouldBeTreated = ((Colors != 0) && (Npg > 0));
//	bool NormalsShouldBeTreated = ((NormCoord != 0) && (Npg > 0));	
//
//	float xP1, yP1, zP1, xP2, yP2, zP2, xN, yN, zN;
//	int *tVertInd = VertInd;
//
//	for(int i=0; i<Npg; i++)
//	{
//		int AmOfPtInPg = Lengths[i];
//		int IndPgColorAndNormCoord = i*3;
//
//		if(ColorsShouldBeTreated)
//		{
//			float* tColor = Colors + IndPgColorAndNormCoord;
//			float r = *(tColor++);
//			float g = *(tColor++);
//			float b = *tColor;
//			if((r >= 0.) && (g >= 0.) && (b >= 0.))
//			{
//				glColor3f(r, g, b);
//			}
//		}
//
//		glBegin(Mode);
//
//		bool FirstPointDefined = false;
//		bool SecondPointDefined = false;
//		bool ThirdPointDefined = false;
//
//		for(int j=0; j<AmOfPtInPg; j++)
//		{
//			int VertInd = *(tVertInd++);
//			double *tVertCoord = VertCoord + ((VertInd - 1)*3);
//
//			float xVert = *(tVertCoord++);
//			float yVert = *(tVertCoord++);
//			float zVert = *tVertCoord;
//
//			if(!FirstPointDefined)
//			{
//				xP1 = xVert; yP1 = yVert; zP1 = zVert; 
//				FirstPointDefined = true;
//			}
//			else if(!SecondPointDefined)
//			{
//				xP2 = xVert; yP2 = yVert; zP2 = zVert; 
//				SecondPointDefined = true;
//
//				if(AmOfPtInPg == 2)
//				{
//					glVertex3f(xP1, yP1, zP1);
//					glVertex3f(xP2, yP2, zP2);
//				}
//			}
//			else if(!ThirdPointDefined)
//			{
//				bool CurNormalIsDefined = false;
//
//				if(NormalsShouldBeTreated && (!NormalsWereDefined))
//				{
//					double xV1 = xP2 - xP1, yV1 = yP2 - yP1, zV1 = zP2 - zP1;
//					double xV2 = xVert - xP2, yV2 = yVert - yP2, zV2 = zVert - zP2;
//
//					if(!MyTest_CheckIfVectorsAreCollinear(xV1, yV1, zV1, xV2, yV2, zV2))
//					{
//						xN = yV1*zV2 - zV1*yV2;
//						yN = zV1*xV2 - xV1*zV2;
//						zN = xV1*yV2 - yV1*xV2;
//
//						float Ne2 = xN*xN + yN*yN + zN*zN;
//						if(Ne2 == 0.) Ne2 = 1e-7;
//						float InvSqrtN = 1./sqrt(Ne2);
//						xN *= InvSqrtN; yN *= InvSqrtN; zN *= InvSqrtN;
//
//						float *tNormCoord = NormCoord + IndPgColorAndNormCoord;
//						*(tNormCoord++) = xN;
//						*(tNormCoord++) = yN;
//						*tNormCoord = zN;
//
//						CurNormalIsDefined = true;
//						ThirdPointDefined = true;
//					}
//				}
//				else if(NormalsShouldBeTreated && NormalsWereDefined)
//				{
//					float *tNormCoord = NormCoord + IndPgColorAndNormCoord;
//					xN = *(tNormCoord++);
//					yN = *(tNormCoord++);
//					zN = *tNormCoord;
//
//					CurNormalIsDefined = true;
//					ThirdPointDefined = true;
//				}
//				else 
//				{
//					ThirdPointDefined = true;
//				}
//
//				if(ThirdPointDefined)
//				{
//					if(CurNormalIsDefined) glNormal3f(xN, yN, zN);
//					glVertex3f(xP1, yP1, zP1);
//					glVertex3f(xP2, yP2, zP2);
//					glVertex3f(xVert, yVert, zVert);
//				}
//			}
//			else
//			{
//				glVertex3f(xVert, yVert, zVert);
//			}
//		}
//		glEnd();
//	}
//	if(NormalsShouldBeTreated) NormalsWereDefined = true;
//}

//-------------------------------------------------------------------------

//void MyTest_reshape(int width, int height)
//{
//			CHWinCont hCurCont;
//			CSimpleGraph::GetCurHWinCont(hCurCont);
//			CWinCont *pc = hCurCont.rep;
//
//    pc->m_xSize = width; 
//    pc->m_ySize = height;
//    pc->m_Aspect = (float)(pc->m_xSize)/(float)(pc->m_ySize);
//    glViewport(0, 0, pc->m_xSize, pc->m_ySize);
//
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    glutPostRedisplay();
//}
//
//void MyTest_display(void) 
//{
//			CHWinCont hCurCont;
//			CSimpleGraph::GetCurHWinCont(hCurCont);
//			CWinCont *pc = hCurCont.rep;
//
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(64.0, pc->m_Aspect, pc->m_DepthNear, pc->m_DepthFar);
//    glMatrixMode(GL_MODELVIEW);
//    glLoadIdentity(); 
//
//	//improve these movements
//    glTranslatef(0.0,0.0,-(pc->m_Depth));
//    glRotatef(-(pc->m_Theta), 1.0, 0.0, 0.0);
//    glRotatef(pc->m_Phi, 0.0, 0.0, 1.0);
//    glTranslatef(-(pc->m_CenVert[0]), -(pc->m_CenVert[1]), -(pc->m_CenVert[2]));
//
//	//DrawObjects();
//	if(pc->m_PgWereDefined)
//		MyTest_DrawArbPolygObjects(GL_POLYGON, pc->m_PgVertCoord, pc->m_PgNv, pc->m_PgVertInd, pc->m_PgLen, pc->m_PgColors, pc->m_PgN, pc->m_PgNormCoord, pc->m_PgNormalsWereDefined);
//
//    glutSwapBuffers();
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//}
//
//void MyTest_keyboard(unsigned char ch, int x, int y)
//{
//			CHWinCont hCurCont;
//			CSimpleGraph::GetCurHWinCont(hCurCont);
//			CWinCont *pc = hCurCont.rep;
//
//    switch (ch) 
//    {
//        case '-': 
//			if(pc->m_Depth < pc->m_DepthFar) pc->m_Depth += pc->m_DepthMotionQuanta; 
//			break;
//        case '+': 
//			if(pc->m_Depth > pc->m_DepthNear) pc->m_Depth -= pc->m_DepthMotionQuanta; 
//			break;
//        case 27: exit(0); break;
//        //case 27: { glutDestroyWindow(gActiveWindowID); gActiveWindowID = 0; glutHackStopMainLoop();} break;
//    }
//    glutPostRedisplay();
//}
//
//void MyTest_mouse(int button, int state, int x, int y)
//{
//			CHWinCont hCurCont;
//			CSimpleGraph::GetCurHWinCont(hCurCont);
//			CWinCont *pc = hCurCont.rep;
//
//    pc->m_DownX = x;
//    pc->m_DownY = y;
//    pc->m_LeftButton = ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN));
//    pc->m_MiddleButton = ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN));
//}
//
//void MyTest_motion(int x, int y)
//{//to improve
//			CHWinCont hCurCont;
//			CSimpleGraph::GetCurHWinCont(hCurCont);
//			CWinCont *pc = hCurCont.rep;
//
//    if(pc->m_LeftButton)
//    {
//        pc->m_Theta += (float)(0.25*(pc->m_DownY - y));
//
//		if(pc->m_Theta >= 180.) pc->m_Theta = pc->m_Theta - 360.;
//		else if(pc->m_Theta < -180.) pc->m_Theta = pc->m_Theta + 360.;
//
//		double PhiMult = 0.25;
//		if(pc->m_Theta <= 0.) PhiMult = -PhiMult;
//
//        pc->m_Phi += (float)(PhiMult*(x - pc->m_DownX));
//
//		//int Aha = 1;
//    }
//    if(pc->m_MiddleButton)
//    {
//        pc->m_Depth += (float)(pc->m_DownY - y) / 10.0;
//    }
//    pc->m_DownX = x;
//    pc->m_DownY = y;
//    glutPostRedisplay();
//}

//-------------------------------------------------------------------------

//GLfloat MyTest_lightPosition[] = { 0.0, 0.0, 1.0, 1.0}; 
//
//int MyTest_StartViewer()
//{
//    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
//    glutInitWindowSize(500, 500);
//    
//	int CurWindowID = glutCreateWindow("3D Viewer");
//		if(CurWindowID != 0) CSimpleGraph::sm_WinContMap[CurWindowID] = CSimpleGraph::sm_hFirstWinCont;
//  
//    glEnable(GL_DEPTH_TEST);
//    glDepthFunc(GL_LEQUAL);
//    glClearColor(1.0, 1.0, 1.0, 0.0);
//    glPolygonOffset(1.0, 1.0);
//    glEnable(GL_CULL_FACE);
//    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
//    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
//    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
//    glEnable(GL_COLOR_MATERIAL);
//    glColorMaterial(GL_FRONT, GL_DIFFUSE);
//    glLightfv (GL_LIGHT0, GL_POSITION, MyTest_lightPosition);
//    glEnable(GL_LIGHT0);
//    
//    MyTest_setSize(MEDIUM);
//    MyTest_setDisplay(SMOOTHSHADED);
//    MyTest_setOther(ENVMAP);
//
//	glutReshapeFunc(MyTest_reshape);
//    glutDisplayFunc(MyTest_display);
//	//glutVisibilityFunc(visibility);
//
//	glutKeyboardFunc(MyTest_keyboard);
//    glutMouseFunc(MyTest_mouse);
//	glutMotionFunc(MyTest_motion);
//
//    glutMainLoop();
//
//    return 0;             /* ANSI C requires main to return int. */
//}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void CWinCont::SetupPgAndLnData(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr, double* CenVert, double* ObjLim, double MaxSize)
{
	m_PgVertCoord = 0;
	m_PgNv = 0;
	m_PgVertInd = 0;
	m_PgLen = 0;
	m_PgN = 0;
	m_PgColors = 0;
	m_PgNormCoord = 0;
	m_PgWereDefined = false;

	m_LnVertCoord = 0;
	m_LnNv = 0;
	m_LnVertInd = 0;
	m_LnLen = 0;
	m_LnN = 0;
	m_LnColors = 0;
	m_LnWereDefined = false;

	if((VertCoordArr[0] != 0) && (NvArr[0] > 0) && (VertIndArr[0] != 0) && (LenArr[0] != 0) && (NpArr[0] != 0))
	{
        m_PgNv = NvArr[0]; //Number of Vertex points
		m_PgN = NpArr[0]; //Number of Polygons

		if(m_PgNv > 0)
		{
			long TotAmOfDouble = 3*m_PgNv;
            m_PgVertCoord = new double[TotAmOfDouble];

			double *tPgVertCoord = m_PgVertCoord, *tInVertCoord = VertCoordArr[0];
			for(long i=0; i<TotAmOfDouble; i++) *(tPgVertCoord++) = *(tInVertCoord++);
		}
		if(m_PgN > 0)
		{
            m_PgLen = new int[m_PgN];

			long TotNumPointsInPolygons = 0;
			int *tPgLen = m_PgLen, *tInLen = LenArr[0];
            for(long j=0; j<m_PgN; j++) 
            {
				TotNumPointsInPolygons += *tInLen;
                *(tPgLen++) = *(tInLen++);
            }

            m_PgVertInd = new int[TotNumPointsInPolygons];

            int *tPgVertInd = m_PgVertInd, *tInVertInd = VertIndArr[0];
            for(long k=0; k<TotNumPointsInPolygons; k++) *(tPgVertInd++) = *(tInVertInd++);

			if(ColorsArr[0] != 0)
			{
				long TotAmOfColorNum = 3*m_PgN;
				m_PgColors = new float[TotAmOfColorNum];

				float *tPgColors = m_PgColors, *tInColors = ColorsArr[0];
				for(long m=0; m<TotAmOfColorNum; m++) *(tPgColors++) = *(tInColors++);
			}
		}

        //m_PgNormCoord = new float[3*m_PgN];
        m_PgNormCoord = new double[3*m_PgN]; //OC050109
        m_PgNormalsWereDefined = false;

		m_PgWereDefined = true;
	}

	if((VertCoordArr[1] != 0) && (NvArr[1] > 0) && (VertIndArr[1] != 0) && (LenArr[1] != 0) && (NpArr[1] != 0))
	{
		m_LnNv = NvArr[1]; //Number of Vertex points
		m_LnN = NpArr[1]; //Number of Lines

		if(m_LnNv > 0)
		{
			long TotAmOfDouble = 3*m_LnNv;
            m_LnVertCoord = new double[TotAmOfDouble];

			double *tLnVertCoord = m_LnVertCoord, *tInVertCoord = VertCoordArr[1];
			for(long i=0; i<TotAmOfDouble; i++) *(tLnVertCoord++) = *(tInVertCoord++);
		}
		if(m_LnN > 0)
		{
            m_LnLen = new int[m_LnN];

			long TotNumPointsInLines = 0;
			int *tLnLen = m_LnLen, *tInLen = LenArr[1];
            for(long j=0; j<m_LnN; j++) 
            {
				TotNumPointsInLines += *tInLen;
                *(tLnLen++) = *(tInLen++);
            }
			
			m_LnVertInd = new int[TotNumPointsInLines];

			int *tLnVertInd = m_LnVertInd, *tInVertInd = VertIndArr[1];
            for(long k=0; k<TotNumPointsInLines; k++) *(tLnVertInd++) = *(tInVertInd++);

			if(ColorsArr[1] != 0)
			{
                long TotAmOfColorNum = 3*m_LnN;
                m_LnColors = new float[TotAmOfColorNum];

                float *tLnColors = m_LnColors, *tInColors = ColorsArr[1];
                for(long m=0; m<TotAmOfColorNum; m++) *(tLnColors++) = *(tInColors++);
			}

            m_LnWereDefined = true;
		}
	}

	if(CenVert != 0)
	{
		for(int i=0; i<3; i++) m_CenVert[i] = CenVert[i];
	}
	if(ObjLim != 0)
	{
		for(int i=0; i<6; i++) m_ObjLimits[i] = ObjLim[i];
	}
    m_MaxSize = MaxSize;
}

//-------------------------------------------------------------------------

CWinCont::~CWinCont()
{
	if(m_PgVertCoord != 0) delete[] m_PgVertCoord;
	if(m_PgLen != 0) delete[] m_PgLen;
	if(m_PgVertInd != 0) delete[] m_PgVertInd;
	if(m_PgColors != 0) delete[] m_PgColors;
	if(m_PgNormCoord != 0) delete[] m_PgNormCoord;

	if(m_LnVertCoord != 0) delete[] m_LnVertCoord;
    if(m_LnLen != 0) delete[] m_LnLen;
    if(m_LnVertInd != 0) delete[] m_LnVertInd;
    if(m_LnColors != 0) delete[] m_LnColors;
}

//-------------------------------------------------------------------------

void CWinCont::ProcDisplayContentForView()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	ProcDisplayContent();

    glutSwapBuffers();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

//-------------------------------------------------------------------------

void CWinCont::ProcReshapeContentForView(int Width, int Height)
{
    ProcReshapeContent(Width, Height);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glutPostRedisplay();
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

bool CSimpleGraph::sm_WasStarted = false;
bool CSimpleGraph::sm_StartingInProcess = false;

//bool CSimpleGraph::sm_WinShouldBeCreated = false;
//vector <int> CSimpleGraph::sm_WinIDsVect;
CHWinContMap CSimpleGraph::sm_WinContMap;
CHWinContList CSimpleGraph::sm_WinToCreateContList;
bool CSimpleGraph::sm_SomethingToCreate = false;
//bool CSimpleGraph::sm_NoMoreWindows = true;
CHWinCont CSimpleGraph::sm_hFirstWinCont;

int CSimpleGraph::sm_DefWinPosX = 10;
int CSimpleGraph::sm_DefWinPosY = 10;
int CSimpleGraph::sm_DefWinSizeX = 500;
int CSimpleGraph::sm_DefWinSizeY = 500;
char CSimpleGraph::sm_DefWinTitle[] = "Simple Graph\0";

//GLRGBQUAD* CSimpleGraph::sm_pDataToSaveToFile = NULL;
char CSimpleGraph::sm_GraphFileName[1024];

//#ifdef WIN32
//struct M__pixelformat__ M_pix[] =
//{
//    /* Double Buffer, alpha */
//    {	{	sizeof(PIXELFORMATDESCRIPTOR),	1,
//        PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_GENERIC_FORMAT|PFD_DOUBLEBUFFER|PFD_SWAP_COPY,
//        PFD_TYPE_RGBA,
//        24,	8,	0,	8,	8,	8,	16,	8,	24,
//        0,	0,	0,	0,	0,	16,	8,	0,	0,	0,	0,	0,	0 },
//        GL_TRUE
//    },
//    /* Single Buffer, alpha */
//    {	{	sizeof(PIXELFORMATDESCRIPTOR),	1,
//        PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_GENERIC_FORMAT,
//        PFD_TYPE_RGBA,
//        24,	8,	0,	8,	8,	8,	16,	8,	24,
//        0,	0,	0,	0,	0,	16,	8,	0,	0,	0,	0,	0,	0 },
//        GL_FALSE
//    },
//    /* Double Buffer, no alpha */
//    {	{	sizeof(PIXELFORMATDESCRIPTOR),	1,
//        PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_GENERIC_FORMAT|PFD_DOUBLEBUFFER|PFD_SWAP_COPY,
//        PFD_TYPE_RGBA,
//        24,	8,	0,	8,	8,	8,	16,	0,	0,
//        0,	0,	0,	0,	0,	16,	8,	0,	0,	0,	0,	0,	0 },
//        GL_TRUE
//    },
//    /* Single Buffer, no alpha */
//    {	{	sizeof(PIXELFORMATDESCRIPTOR),	1,
//        PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_GENERIC_FORMAT,
//        PFD_TYPE_RGBA,
//        24,	8,	0,	8,	8,	8,	16,	0,	0,
//        0,	0,	0,	0,	0,	16,	8,	0,	0,	0,	0,	0,	0 },
//        GL_FALSE
//    },
//};
//#endif

//-------------------------------------------------------------------------

void CSimpleGraph::Init()
{
    *sm_GraphFileName = '\0';
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    //glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH);
	//make basic initialisations here
}

//-------------------------------------------------------------------------

void CSimpleGraph::WrapForMT_glutMainLoop(void* InData)
{//wrap for Multi-Threading entry (on Win32)

    StartNewWin(sm_hFirstWinCont); //move args to input
    
	//glutIdleFunc(ProcIdle);

	int CallInterv_ms = 400; //to tune
    glutTimerFunc(CallInterv_ms, ProcTimer, CallInterv_ms);
    sm_WasStarted = true;

    glutMainLoop();
}

//-------------------------------------------------------------------------

void CSimpleGraph::StartNewWinApp(CHWinCont& hWinCont, char StartMode)
{
    //sm_hFirstWinCont = hWinCont;
	//MyTest_StartViewer();

	if(StartMode == 1) //multi-thread
	{
		if(sm_StartingInProcess)
		{
			while(sm_StartingInProcess)
			{
#ifdef WIN32
				Sleep(100); 
#else
		//to program MT for Linux and MacOS X
#endif
				continue;
			}
		}
	}

	if(!sm_WasStarted) // or sm_WasStarted
	{
		sm_StartingInProcess = true;

        Init();
		sm_hFirstWinCont = hWinCont;

		if(StartMode == 1)
		{
#ifdef WIN32
			//uintptr_t res = _beginthread(CSimpleGraph::WrapForMT_glutMainLoop, 0, NULL);
			_beginthread(CSimpleGraph::WrapForMT_glutMainLoop, 0, NULL);
#else
			//to program MT for Linux and MacOS X
#endif
		}
		else
		{
            WrapForMT_glutMainLoop(NULL);
		}
		//CSimpleGraph::WrapForMT_glutMainLoop(NULL);
	}
	else
	{
		CHWinCont hLocWinCont = hWinCont;
		sm_WinToCreateContList.push_back(hLocWinCont);
		sm_SomethingToCreate = true;

		//OC041107
		int CallInterv_ms = 400; //to tune
		glutTimerFunc(CallInterv_ms, ProcTimer, CallInterv_ms);
	}
}

//-------------------------------------------------------------------------

//void CSimpleGraph::StartNewWin(CHWinCont& hWinCont, int xPos, int yPos, int xSize, int ySize, char* strTitle)
void CSimpleGraph::StartNewWin(CHWinCont& hWinCont, int xPos, int yPos, int xSize, int ySize)
{
	char StrTitleBuf[1024];
	*StrTitleBuf = '\0';

	//if((strTitle == 0) || (*strTitle == '\0')) strcpy(StrTitleBuf, sm_DefWinTitle);
	//else strcpy(StrTitleBuf, strTitle);
	if(hWinCont.rep != 0)
	{
		if(hWinCont.rep->m_WinTitle[0] != '\0') strcpy(StrTitleBuf, hWinCont.rep->m_WinTitle);
		else strcpy(StrTitleBuf, sm_DefWinTitle);
	}

	glutInitWindowPosition(xPos, yPos);
    glutInitWindowSize(xSize, ySize);

	int CurWindowID = glutCreateWindow(StrTitleBuf);
	if(CurWindowID != 0) sm_WinContMap[CurWindowID] = hWinCont; //add to container

	if(hWinCont.rep != 0)
	{
        hWinCont.rep->SetViewingParamsGL();
        hWinCont.rep->SetupMenuGLUT();
	}

	glutReshapeFunc(ProcReshape);
	glutDisplayFunc(ProcDisplay);
	//glutVisibilityFunc(ProcVisibility);

    glutKeyboardFunc(ProcKeyboard);
	glutMouseFunc(ProcMouse);
	glutMotionFunc(ProcMotion);
	//glutPassiveMotionFunc(ProcPassiveMotion);
	//glutWindowStatusFunc(ProcWinStatus);
}

//void CSimpleGraph::StartNewWin(CHWinCont& hWinCont, int xPos, int yPos, int xSize, int ySize)
//{
//    StartNewWin(hWinCont, xPos, yPos, xSize, ySize, sm_DefWinTitle);
//}

void CSimpleGraph::StartNewWin(CHWinCont& hWinCont)
{
	const int MinPosX = -2000;
	const int MinPosY = -2000;
	int PosX = sm_DefWinPosX, PosY = sm_DefWinPosY;
	int SizeX = sm_DefWinSizeX, SizeY = sm_DefWinSizeY;
	if(hWinCont.rep != 0)
	{
		if(hWinCont.rep->m_xPos >= 0) PosX = hWinCont.rep->m_xPos;
		if(hWinCont.rep->m_yPos >= 0) PosY = hWinCont.rep->m_yPos;
		if(hWinCont.rep->m_xSize > 0) SizeX = hWinCont.rep->m_xSize;
		if(hWinCont.rep->m_ySize > 0) SizeY = hWinCont.rep->m_ySize;
	}
    StartNewWin(hWinCont, PosX, PosY, SizeX, SizeY);
}

//-------------------------------------------------------------------------

void CSimpleGraph::StartNewWindowIfNecessary()
{//this function starts new windows after creation of the first window
	if(!sm_WinToCreateContList.empty())
	{
		int AmOfWinToCreate = (int)(sm_WinToCreateContList.size());
		if(AmOfWinToCreate > 0)
		{
			for(int i=0; i<AmOfWinToCreate; i++)
			{
				CHWinCont hCont = *(sm_WinToCreateContList.begin());
				StartNewWin(hCont);
				sm_WinToCreateContList.pop_front();
			}
			if(!sm_WinToCreateContList.empty()) sm_WinToCreateContList.clear();
			sm_SomethingToCreate = false;
		}
	}
}

//-------------------------------------------------------------------------

void CSimpleGraph::ProcReshape(int Width, int Height)
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	hCurCont.rep->ProcReshapeContentForView(Width, Height);
	hCurCont.rep->ProcDisplayContentForView(); // to ensure continious update at resizing
	if(sm_StartingInProcess) sm_StartingInProcess = false;
}

//-------------------------------------------------------------------------

void CSimpleGraph::ProcDisplay()
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	hCurCont.rep->ProcDisplayContentForView();

	if(*sm_GraphFileName != '\0') FinishSavingWinContentToGraphFile(hCurCont);
	if(sm_StartingInProcess) sm_StartingInProcess = false;
}

//-------------------------------------------------------------------------

void CSimpleGraph::ProcKeyboard(unsigned char ch, int x, int y)
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	hCurCont.rep->ProcKeyboardContent(ch, x, y);
	hCurCont.rep->ProcDisplayContentForView(); // to ensure continious update
	//ProcDisplay(); // to ensure continious update
}

//-------------------------------------------------------------------------

void CSimpleGraph::ProcMouse(int button, int state, int x, int y)
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	hCurCont.rep->ProcMouseContent(button, state, x, y);
}

//-------------------------------------------------------------------------

void CSimpleGraph::ProcMotion(int x, int y)
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	hCurCont.rep->ProcMotionContent(x, y);
	hCurCont.rep->ProcDisplayContentForView(); // to ensure continious update
	//glutPostRedisplay();
}

//-------------------------------------------------------------------------

//void CSimpleGraph::ProcPassiveMotion(int x, int y)
//{
//	CHWinCont hCurCont;
//	if(!GetCurHWinCont(hCurCont)) return;
//	if(hCurCont.rep == 0) return;
//
//	//hCurCont.rep->ProcMotionContent(x, y);
//	//glutPostRedisplay();
//}

//-------------------------------------------------------------------------

//void CSimpleGraph::ProcWinStatus(int StatusID)
//{
//	CHWinCont hCurCont;
//	if(!GetCurHWinCont(hCurCont)) return;
//	if(hCurCont.rep == 0) return;
//
//    //glutPostRedisplay();
//	ProcDisplay();
//}

//-------------------------------------------------------------------------

void CSimpleGraph::DoSetSize(int SizeCase)
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	hCurCont.rep->DoSetSizeContent(SizeCase);
	hCurCont.rep->ProcDisplayContentForView(); // to ensure continious update
}

//-------------------------------------------------------------------------

void CSimpleGraph::DoSetMotion(int MotionCase)
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	hCurCont.rep->DoSetMotionContent(MotionCase);
	hCurCont.rep->ProcDisplayContentForView(); // to ensure continious update
}

//-------------------------------------------------------------------------

void CSimpleGraph::DoSetMain(int MainCase)
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	//hCurCont.rep->DoSetMainContent(MainCase);
	hCurCont.rep->ProcDisplayContentForView(); // to ensure continious update

	if(MainCase == CWinCont::sm_EXIT)
	{
		int CurWinID = glutGetWindow();
		if(CurWinID != 0) DestroyWindow(CurWinID);
		//exit(0);
	}
	else if(MainCase == CWinCont::sm_SAVE)
	{
		int CurWinID = glutGetWindow();
		if(CurWinID != 0) SaveWinContentToGraphFile(CurWinID);
	}
}

//-------------------------------------------------------------------------

void CSimpleGraph::DestroyWindow(int WinID)
{
    if(WinID == 0) return;
    glutDestroyWindow(WinID);

	if(sm_WinContMap.empty()) return;
	CHWinContMap::iterator iter = sm_WinContMap.find(WinID);
	if(iter == sm_WinContMap.end()) return;
	sm_WinContMap.erase(iter);
}

//-------------------------------------------------------------------------

bool CSimpleGraph::ExtractDispContentData(int WinID, CHWinCont& hCont)
{
	if(sm_WinContMap.empty()) return false;
	CHWinContMap::const_iterator iter = sm_WinContMap.find(WinID);
	if(iter == sm_WinContMap.end())
	{
		//any other action ?
		return false;
	}
	hCont = (*iter).second;
	return true;
}

//-------------------------------------------------------------------------

bool CSimpleGraph::GetCurHWinCont(CHWinCont& hCurCont)
{
	int CurWinID = glutGetWindow();
	if(!ExtractDispContentData(CurWinID, hCurCont)) return false;
	return true;
}

//-------------------------------------------------------------------------

void CSimpleGraph::SaveWinContentToGraphFile(int WinID) 
{
	int DefaultFileToScreenResolRatio = 2;
	if(WinID == 0) return;

	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	//hCurCont.rep->ProcDisplayContentForView(); // to ensure continious update at resizing

	hCurCont.rep->m_xPos = glutGet(GLUT_WINDOW_X);
	hCurCont.rep->m_yPos = glutGet(GLUT_WINDOW_Y);

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	int CurWidth = viewport[2], CurHeight = viewport[3], CurColorDepth = 32;

	int SaveWidth = DefaultFileToScreenResolRatio*CurWidth;
	int SaveHeight = DefaultFileToScreenResolRatio*CurHeight;

	char sFileName[1024];
	sFileName[0] = '\0';
	GetFileName(sFileName);
	if(sFileName[0] == '\0') return;

	strcpy(sm_GraphFileName, sFileName); //since sm_GraphFileName != '\0' FinishSavingWinContentToGraphFile will be called 

/**
//The code below worked last time with OSMesa off-screen rendering to improve resolution

	hCurCont.rep->ProcDisplayContentForView(); // to ensure continious update at resizing
    glFinish(); //wait until concurrent OpenGL processing ends

	GLRGBQUAD *pData = new GLRGBQUAD[SaveWidth*SaveHeight];
    if(!pData) throw "Memory Allocation Failure";
	bool SwapRB = false;

	//CRGBTri *pData3 = 0;

//TO IMPLEMENT: Off-screen rendering to improve resolution of saved files!!!
//#ifdef WIN32
//
//    PIXELFORMATDESCRIPTOR pfdAuxInf;
//	int res;
//
//	//HDC Old_hDC = wglGetCurrentDC();
//	//HGLRC Old_hGLRC = wglGetCurrentContext();
//
//		HDC Old_hDC = CreateDC("DISPLAY", NULL, NULL, NULL); 
//		//HDC hdcCompatible = CreateCompatibleDC(Old_hDC); 
//
//	int SaveInd = SaveDC(Old_hDC);
//
//	HDC hMemDC = CreateCompatibleDC(Old_hDC);
////	//HDC hDC = Old_hDC;
//
//	HBITMAP Aux_hBMP = CreateCompatibleBitmap(Old_hDC, SaveWidth, SaveHeight);
//	HGDIOBJ res00 = SelectObject(hMemDC, Aux_hBMP);
//
//	int nExistingPixFormat = GetPixelFormat(Old_hDC);
//	if(nExistingPixFormat != 0)
//	{
//		DescribePixelFormat(Old_hDC, nExistingPixFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfdAuxInf);
//		//if(!(pfdAuxInf.dwFlags & PFD_NEED_PALETTE)) return;
//
//		if(pfdAuxInf.dwFlags & PFD_DRAW_TO_WINDOW) pfdAuxInf.dwFlags ^= PFD_DRAW_TO_WINDOW;
//		pfdAuxInf.dwFlags |= PFD_DRAW_TO_BITMAP;
//	}
//	else
//	{
//		memset(&pfdAuxInf, 0, sizeof(PIXELFORMATDESCRIPTOR));
//		pfdAuxInf.nSize = sizeof(PIXELFORMATDESCRIPTOR);
//		pfdAuxInf.nVersion = 1;
//		pfdAuxInf.dwFlags = PFD_DRAW_TO_BITMAP | PFD_SUPPORT_OPENGL | PFD_SUPPORT_GDI; // replaces PFD_DRAW_TO_WINDOW
//		//pfdAuxInf.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_SUPPORT_GDI | PFD_DOUBLEBUFFER; // replaces PFD_DRAW_TO_WINDOW
//		pfdAuxInf.iPixelType = PFD_TYPE_RGBA; 
//		pfdAuxInf.cColorBits = 24; //32; //24; //32; //8;
//		pfdAuxInf.cDepthBits = 32; //8; //16;
//		pfdAuxInf.iLayerType = PFD_MAIN_PLANE; 
//	}
//
////	int Old_iPixFormat = GetPixelFormat(Old_hDC);
////	res = DescribePixelFormat(Old_hDC, Old_iPixFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfdAuxInf);
////
////	//int iPixFormat = GetPixelFormat(hdcCompatible);
////	//res = DescribePixelFormat(hdcCompatible, iPixFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfdAuxInf);
//
//	BITMAPINFO bmi;
//	int SizeBITMAPINFOHEADER = 40; //sizeof(BITMAPINFOHEADER); 
//	long SizeImage = SaveWidth*SaveHeight*4;
//	memset(&bmi, 0, sizeof(BITMAPINFO));
//    bmi.bmiHeader.biSize = SizeBITMAPINFOHEADER;
//    bmi.bmiHeader.biWidth = SaveWidth;
//    bmi.bmiHeader.biHeight = SaveHeight;
//	bmi.bmiHeader.biBitCount = 24; //32; 
//	bmi.bmiHeader.biPlanes = 1;
//	bmi.bmiHeader.biCompression = BI_RGB;
//	bmi.bmiHeader.biSizeImage = SizeImage;
//
//	//// Rearrange RGB component storage from BGR to RGB.
//	//bmi.bmiColors[0].rgbBlue  = 0xff ; // Store red in blue's normal position.
//	//bmi.bmiColors[1].rgbGreen = 0xff ; // Green stays the same.
//	//bmi.bmiColors[2].rgbRed   = 0xff ; // Store blue in red's normal position.
//
//	//PIXELFORMATDESCRIPTOR pfd;
// //   memset(&pfd,0,sizeof(PIXELFORMATDESCRIPTOR));
// //   pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR);
// //   pfd.nVersion = 1;
// //   //pfd.dwFlags = PFD_DRAW_TO_BITMAP | PFD_SUPPORT_OPENGL | PFD_SUPPORT_GDI; // replaces PFD_DRAW_TO_WINDOW
//	//pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_SUPPORT_GDI | PFD_DOUBLEBUFFER; // replaces PFD_DRAW_TO_WINDOW
// //   pfd.iPixelType = PFD_TYPE_RGBA; 
// //   pfd.cColorBits = 24; //32; //24; //32; //8;
// //   pfd.cDepthBits = 32; //8; //16;
// //   pfd.iLayerType = PFD_MAIN_PLANE; 
//
////PIXELFORMATDESCRIPTOR pfd = { 
////    sizeof(PIXELFORMATDESCRIPTOR),  //  size of this pfd 
////    1,                     // version number 
////    PFD_DRAW_TO_WINDOW |   // support window 
////    PFD_SUPPORT_OPENGL |   // support OpenGL 
////    PFD_DOUBLEBUFFER,      // double buffered 
////    PFD_TYPE_RGBA,         // RGBA type 
////    24,                    // 24-bit color depth 
////    0, 0, 0, 0, 0, 0,      // color bits ignored 
////    0,                     // no alpha buffer 
////    0,                     // shift bit ignored 
////    0,                     // no accumulation buffer 
////    0, 0, 0, 0,            // accum bits ignored 
////    32,                    // 32-bit z-buffer     
////    0,                     // no stencil buffer 
////    0,                     // no auxiliary buffer 
////    PFD_MAIN_PLANE,        // main layer 
////    0,                     // reserved 
////    0, 0, 0                // layer masks ignored 
////};
//
////    HDC  hdc;
////    int  iPixelFormat; 
//// 
////iPixelFormat = ChoosePixelFormat(hdc, &pfd);
//
////    //bool res0 = SetPixelFormat(hDC, Old_iPixFormat, &pfd);
////
////	SwapRB = true;
//	void *BitData=0;
//	HBITMAP hBMP = CreateDIBSection(hMemDC, &bmi, DIB_RGB_COLORS, &BitData, NULL, NULL);
//	HGDIOBJ Old_Obj = SelectObject(hMemDC, hBMP);
//	pData = (GLRGBQUAD*)BitData;
//
////	pData3 = (CRGBTri*)BitData;
////
////
////	iPixFormat = GetPixelFormat(hDC);
////	res = DescribePixelFormat(hDC, iPixFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfdAuxInf);
////
////    //hbrush = GetCurrentObject(hDC, OBJ_BITMAP); 
//// 
////    // Retrieve a LOGBRUSH structure that contains the 
////    // current brush attributes. 
//// 
////    //GetObject(hbrush, sizeof(LOGBRUSH), &lb); 
////
////	int hdcPixFormat = GetPixelFormat(hDC);
////	if(hdcPixFormat != 0)
////	{
////		res = DescribePixelFormat(hDC, hdcPixFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfdAuxInf);
////		res += 1;
////	}
////
////	int hdcOldPixFormat = GetPixelFormat(Old_hDC);
////	res = DescribePixelFormat(Old_hDC, hdcOldPixFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfdAuxInf);
//
//    //int nPixelFormat = ChoosePixelFormat(hMemDC, &pfd);
//    int nPixelFormat = ChoosePixelFormat(hMemDC, &pfdAuxInf);
//	//int nPixelFormat = M_wglChoosePixelFormat(hMemDC, &pfdAuxInf);
//
//    if(nPixelFormat == 0) 
//	{
//		int ErrNo = GetLastError();
//		return; //throw error here
//	}
//
//    //bool bResult = SetPixelFormat(hMemDC, nPixelFormat, &pfd); // Set the pixel format.
//    bool bResult = SetPixelFormat(hMemDC, nPixelFormat, &pfdAuxInf); // Set the pixel format.
////    //bool bResult = SetPixelFormat(hDC, 54, &pfd); // Set the pixel format.
////
////	if(hdcPixFormat != 0)
////	{
////        bool bResult = SetPixelFormat(hDC, hdcPixFormat, &pfd); // Set the pixel format.
////	}
////
////
////    //HGLRC hGLRC = wglCreateContext(hDC);
//    HGLRC hGLRC = wglCreateContext(hMemDC);
//
//	if(hGLRC == 0)
//	{
//			int ErrNo = GetLastError();	
//			return; //throw error here
//	}
//
//    wglMakeCurrent(hMemDC, hGLRC);
////#endif

	OSMesaContext ctx = OSMesaCreateContextExt(OSMESA_RGBA, CurColorDepth, 0, 0, NULL); //OC090105
	if(!OSMesaMakeCurrent(ctx, pData, GL_UNSIGNED_BYTE, SaveWidth, SaveHeight)) return; //OC090105

	//glutReshapeWindow(SaveWidth, SaveHeight);
	//hCurCont.rep->ProcReshapeContent(SaveWidth, SaveHeight);

	RenderImageToSave(WinID, SaveWidth, SaveHeight);

	//sm_pDataToSaveToFile = pData;
	//strcpy(sm_GraphFileName, sFileName);
	//hCurCont.rep->m_xSizeOrig = CurWidth;
	//hCurCont.rep->m_ySizeOrig = CurHeight;

	//glFlush();
    //glFinish();

	//GLboolean Param;
	//glGetBooleanv(GL_RGBA_MODE, &Param);

	//glReadBuffer(GL_FRONT_LEFT);
	//glPixelStorei(GL_PACK_ALIGNMENT, 4);
	//glPixelStorei(GL_PACK_ROW_LENGTH, 0);
	//glPixelStorei(GL_PACK_SKIP_ROWS, 0);
	//glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
    //glReadPixels(0, 0, SaveWidth, SaveHeight, GL_RGBA, GL_UNSIGNED_BYTE, pData);
	
	//AuxMakeWhite(pData, SaveWidth);

    SaveImageToGraphFile(pData, SaveWidth, SaveHeight, SwapRB, sFileName);
	delete[] pData;

    OSMesaDestroyContext(ctx); //OC090105

//#ifdef WIN32
//
//    RestoreDC(Old_hDC, SaveInd);
//
//
//	if(!wglMakeCurrent(Old_hDC, Old_hGLRC)) {};
//
//    if(!wglDeleteContext(hGLRC)) {}
//	//if(!DeleteDC(hDC)) {}
//	if(!DeleteObject(hBMP)) {}
//#endif

//	glutSetWindow(WinID);
	//RenderImageToSave(WinID, CurWidth, CurHeight);
//    hCurCont.rep->ProcDisplayContentForView();

	hCurCont.rep->m_xSize = CurWidth;
	hCurCont.rep->m_ySize = CurHeight;
	DestroyWindow(WinID);
	StartNewWinApp(hCurCont);
**/
}

//-------------------------------------------------------------------------

void CSimpleGraph::FinishSavingWinContentToGraphFile(CHWinCont& hCurCont)
{
	if(hCurCont.rep == 0) return;
	long SaveWidth = hCurCont.rep->m_xSize;
	long SaveHeight = hCurCont.rep->m_ySize;
	if((*sm_GraphFileName == '\0') || (SaveWidth == 0) || (SaveHeight == 0)) return;

	glFlush();
    glFinish();

	GLboolean Param;
	glGetBooleanv(GL_RGBA_MODE, &Param);

	glReadBuffer(GL_FRONT_LEFT);
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
	glPixelStorei(GL_PACK_ROW_LENGTH, 0);
	glPixelStorei(GL_PACK_SKIP_ROWS, 0);
	glPixelStorei(GL_PACK_SKIP_PIXELS, 0);

	GLRGBQUAD *pData = new GLRGBQUAD[SaveWidth*SaveHeight];
    if(!pData) throw "Memory Allocation Failure";

    glReadPixels(0, 0, SaveWidth, SaveHeight, GL_RGBA, GL_UNSIGNED_BYTE, pData);
	if((hCurCont.rep->m_xSizeOrig > 0) && (hCurCont.rep->m_ySizeOrig > 0)) 
		glutReshapeWindow(hCurCont.rep->m_xSizeOrig, hCurCont.rep->m_ySizeOrig);

    AuxMakeWhite(pData, SaveWidth); //first line - to correct for glut problem

	bool SwapRB = false;
    SaveImageToGraphFile(pData, SaveWidth, SaveHeight, SwapRB, sm_GraphFileName);
	delete[] pData;
	*sm_GraphFileName = '\0';
}

//-------------------------------------------------------------------------

bool CSimpleGraph::SaveImageToGraphFile(GLRGBQUAD *pData, int Width, int Height, bool SwapRB, const char *sFileName)
{
	if((pData == 0) || (Width == 0) || (Height == 0)) return false;

	const char *ext = &sFileName[strlen(sFileName) - 3]; //get extension: to improve - search for dot !
	if(SwapRB) BGRA2RGBA(pData, Width, Height);

	if(strcmp(ext, "bmp") == 0) return SaveImageToBMP(pData, Width, Height, sFileName);
	//if(strcmp(ext, "jpg") == 0) return SaveImageToJPG(pData, Width, Height, sFileName);
	if(strcmp(ext, "png") == 0) return SaveImageToPNG(pData, Width, Height, sFileName);
	//if(strcmp(ext, "tga") == 0) return SaveImageToTGA(sFileName);
	return false;
}

//-------------------------------------------------------------------------

bool CSimpleGraph::SaveImageToBMP(GLRGBQUAD *pData, int Width, int Height, const char *filename)
{
	const int BM = 19778;

	BITMAPFILEHEADER_HACK fHeader;
	BITMAPINFOHEADER_HACK iHeader;
	FILE *bmpFile;
    
    long length = Width*Height*4;

    if(!pData) return false;
    if(!(bmpFile = fopen(filename, "wb"))) return false;

	int SizeBITMAPFILEHEADER = 14; //sizeof(BITMAPFILEHEADER_HACK);
	int SizeBITMAPINFOHEADER = 40;

    fHeader.bfType = BM;
	fHeader.bfSize = SizeBITMAPFILEHEADER + SizeBITMAPINFOHEADER + length;
    fHeader.bfReserved1 = 0; 
    fHeader.bfReserved2 = 0; 

    // This sets the distance from the start of the file to the start
    // of the bitmaps color data
    fHeader.bfOffBits = SizeBITMAPFILEHEADER + SizeBITMAPINFOHEADER;
	
	//long bytesWritten = fwrite(&fHeader, 1, SizeBITMAPFILEHEADER, bmpFile);
	//if(bytesWritten != SizeBITMAPFILEHEADER) { fclose(bmpFile); return false;}
	//long bytesWritten = fwrite(&fHeader, SizeBITMAPFILEHEADER, 1, bmpFile);
	//if(bytesWritten != SizeBITMAPFILEHEADER) { fclose(bmpFile); return false;}

	long bytesWritten = (long)fwrite(&fHeader.bfType, 1, 2, bmpFile); if(bytesWritten != 2) { fclose(bmpFile); return false;}
	bytesWritten = (long)fwrite(&fHeader.bfSize, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
	bytesWritten = (long)fwrite(&fHeader.bfReserved1, 1, 2, bmpFile); if(bytesWritten != 2) { fclose(bmpFile); return false;}
	bytesWritten = (long)fwrite(&fHeader.bfReserved2, 1, 2, bmpFile); if(bytesWritten != 2) { fclose(bmpFile); return false;}
	bytesWritten = (long)fwrite(&fHeader.bfOffBits, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}

	//memset(&iHeader, 0, sizeof(iHeader));
	memset(&iHeader, 0, SizeBITMAPINFOHEADER);

    // Set up the information header
    iHeader.biSize = SizeBITMAPINFOHEADER; //sizeof(BITMAPINFOHEADER_HACK);
    iHeader.biWidth = Width;         // Current width
    iHeader.biHeight = Height;       // Current height
    iHeader.biPlanes = 1;              // Number of planes, must be set to 1
    iHeader.biBitCount = 32;           // Current color depth
    iHeader.biCompression = BI_RGB;    // No compression
    iHeader.biSizeImage = length;      // Number of bytes in bitmap

	//bytesWritten = fwrite(&iHeader, 1, SizeBITMAPINFOHEADER, bmpFile);
	//if(bytesWritten != SizeBITMAPINFOHEADER) { fclose(bmpFile); return false;}
    //bytesWritten = fwrite(&iHeader, SizeBITMAPINFOHEADER, 1, bmpFile);
	//if(bytesWritten != SizeBITMAPINFOHEADER) { fclose(bmpFile); return false;}

    bytesWritten = (long)fwrite(&iHeader.biSize, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biWidth, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biHeight, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biPlanes, 1, 2, bmpFile); if(bytesWritten != 2) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biBitCount, 1, 2, bmpFile); if(bytesWritten != 2) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biCompression, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biSizeImage, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biXPelsPerMeter, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biYPelsPerMeter, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biClrUsed, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}
    bytesWritten = (long)fwrite(&iHeader.biClrImportant, 1, 4, bmpFile); if(bytesWritten != 4) { fclose(bmpFile); return false;}

    AuxSwapRB(pData, Width, Height);

	bytesWritten = (long)fwrite(pData, 1, length, bmpFile);
	//bytesWritten = fwrite(pData, 4, Width*Height, bmpFile);
	if(bytesWritten != length) { fclose(bmpFile); return false;}

    fclose(bmpFile);
    return true;
}

//-------------------------------------------------------------------------

void CSimpleGraph::AuxSwapRB(GLRGBQUAD *pixel, int Width, int Height)
{
    for(long i = 0; i < (Width*Height); i++)
    {
        unsigned char tmp = pixel->blue;
        pixel->blue = pixel->red;
        pixel->red = tmp;
		pixel++;
    }
}

//-------------------------------------------------------------------------

bool CSimpleGraph::SaveImageToPNG(GLRGBQUAD *pData, int Width, int Height, const char *filename)
{
	png_structp png = 0;
	png_infop pngInfo = 0;
	FILE *pngFile;
	png_byte **rowPtrs;
	int i;

	if(!pData) return false;
    RotateBitmapVertically(pData, Width, Height);
	//CorrectAlphaInBitmapRGBA(pData, Width, Height);

	long NumPix = Width*Height;
	CRGBTri *LocArrRGB = new CRGBTri[NumPix];
	CopyRGBA2RBG(pData, LocArrRGB, Width, Height);

	if(!(pngFile = fopen(filename, "wb"))) return false;

	png = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	if(!png) { fclose(pngFile); return false;}

	pngInfo = png_create_info_struct(png);
	if(!pngInfo)
	{
		png_destroy_read_struct(&png, 0, 0);
		fclose(pngFile);
	}

	png_init_io(png, pngFile);

    png_set_compression_level(png, Z_BEST_COMPRESSION);

	//png_set_IHDR(png, pngInfo, Width, Height, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_set_IHDR(png, pngInfo, Width, Height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

	png_write_info(png, pngInfo);

	rowPtrs = new png_byte * [Height];
	if(!rowPtrs)
	{
		fclose(pngFile);
		png_destroy_read_struct(&png, &pngInfo, 0);
		return false;
	}

	//for(i = 0; i < Height; i++) 
	//{
	//	rowPtrs[i] = (unsigned char *)pData + (i * Width * 4);
	//}
	for(i = 0; i < Height; i++) 
	{
		rowPtrs[i] = (unsigned char *)LocArrRGB + (i * Width * 3);
	}

	png_write_image(png, rowPtrs);
	png_write_end(png, pngInfo);
	png_write_flush(png);

	delete[] rowPtrs;

	png_destroy_read_struct(&png, &pngInfo, 0);
	fclose(pngFile);

	if(LocArrRGB != 0) delete[] LocArrRGB;
	return true;
}

//-------------------------------------------------------------------------

bool CSimpleGraph::SaveImageToPNG(CRGBTri *pData, int Width, int Height, const char *filename)
{
	png_structp png = 0;
	png_infop pngInfo = 0;
	FILE *pngFile;
	png_byte **rowPtrs;
	int i;

	if(!pData) return false;
    //RotateBitmapVertically(pData, Width, Height);

	long NumPix = Width*Height;
	//CRGBTri *LocArrRGB = new CRGBTri[NumPix];
	//CopyRGBA2RBG(pData, LocArrRGB, Width, Height);

	if(!(pngFile = fopen(filename, "wb"))) return false;

	png = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	if(!png) { fclose(pngFile); return false;}

	pngInfo = png_create_info_struct(png);
	if(!pngInfo)
	{
		png_destroy_read_struct(&png, 0, 0);
		fclose(pngFile);
	}

	png_init_io(png, pngFile);

    png_set_compression_level(png, Z_BEST_COMPRESSION);

	//png_set_IHDR(png, pngInfo, Width, Height, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_set_IHDR(png, pngInfo, Width, Height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

	png_write_info(png, pngInfo);

	rowPtrs = new png_byte * [Height];
	if(!rowPtrs)
	{
		fclose(pngFile);
		png_destroy_read_struct(&png, &pngInfo, 0);
		return false;
	}

	//for(i = 0; i < Height; i++) 
	//{
	//	rowPtrs[i] = (unsigned char *)pData + (i * Width * 4);
	//}
	for(i = 0; i < Height; i++) 
	{
		rowPtrs[i] = (unsigned char *)pData + (i * Width * 3);
	}

	png_write_image(png, rowPtrs);
	png_write_end(png, pngInfo);
	png_write_flush(png);

	delete[] rowPtrs;

	png_destroy_read_struct(&png, &pngInfo, 0);
	fclose(pngFile);

	//if(LocArrRGB != 0) delete[] LocArrRGB;
	return true;
}

//-------------------------------------------------------------------------

void CSimpleGraph::RotateBitmapVertically(GLRGBQUAD *pData, int Width, int Height)
{
	if((pData == 0) || (Width <= 0) || (Height <= 0)) return;

	//Rotating the bitmap
	for(int k=0; k<(Height >> 1); k++)
	{
		GLRGBQUAD *pCurRow = pData + k*Width;
		GLRGBQUAD *pCurRowCompl = pData + (Height - k - 1)*Width;
		GLRGBQUAD *tCurRow = pCurRow, *tCurRowCompl = pCurRowCompl;

		for(int j=0; j<Width; j++) 
		{
			GLRGBQUAD AuxPix = *tCurRow;
			*tCurRow = *tCurRowCompl;
			*tCurRowCompl = AuxPix;
			tCurRow++; tCurRowCompl++;
		}
	}
}

//-------------------------------------------------------------------------

void CSimpleGraph::CorrectAlphaInBitmapRGBA(GLRGBQUAD *pData, int Width, int Height)
{
	if((pData == 0) || (Width <= 0) || (Height <= 0)) return;

	GLRGBQUAD *tPix = pData;
	for(int k=0; k<Height; k++)
	{
		for(int j=0; j<Width; j++) 
		{
			if((tPix->alpha != 0) && (tPix->alpha != 255))
			{
				if((tPix->red == 255) && (tPix->green == 255) && (tPix->blue == 255)) tPix->alpha = 0;
				else tPix->alpha = 255;
			}
			tPix++;
		}
	}
}

//-------------------------------------------------------------------------

void CSimpleGraph::CopyRGBA2RBG(GLRGBQUAD *pRGBA, CRGBTri *pRGB, int Width, int Height)
{
	if((pRGBA == 0) || (pRGB == 0) || (Width == 0) || (Height == 0)) return;

	GLRGBQUAD *tRGBA = pRGBA;
    CRGBTri *tRGB = pRGB;
	for(int k=0; k<Height; k++)
	{
		for(int j=0; j<Width; j++) 
		{
			tRGB->red = tRGBA->red;
			tRGB->green = tRGBA->green;
			tRGB->blue = tRGBA->blue;
			tRGB++; tRGBA++;
		}
	}
}

//-------------------------------------------------------------------------

void CSimpleGraph::BGRA2RGBA(GLRGBQUAD *pData, int Width, int Height)
{
	if((pData == 0) || (Width == 0) || (Height == 0)) return;

	GLRGBQUAD *tData = pData;
    unsigned char buf = 0;	
	for(int k=0; k<Height; k++)
	{
		for(int j=0; j<Width; j++) 
		{
			buf = tData->red;
			tData->red = tData->blue;
            tData->blue = buf;
			tData++;
		}
	}
}

//-------------------------------------------------------------------------

void CSimpleGraph::AuxMakeWhite(GLRGBQUAD *pData, long Size)
{
	if((pData == 0) || (Size == 0)) return;
	GLRGBQUAD *tData = pData;
	for(long j=0; j<Size; j++) 
	{
		tData->red = 255;
		tData->green = 255;
		tData->blue = 255;
		tData++;
	}
}

//-------------------------------------------------------------------------
/**
bool CSimpleGraph::SaveImageToJPG(GLRGBQUAD *pData, int Width, int Height, const char *filename)
{
	int jpgQuality = 100; //bw 1 and 100
    int jErr;
    JPEG_CORE_PROPERTIES jpgProps;
	PFNIJLINIT ijlInit;
    PFNIJLWRITE ijlWrite;

	jErr = ijlInit(&jpgProps);
	if(jErr	!= IJL_OK) return false;

	jpgProps.DIBWidth =	Width;
	jpgProps.DIBHeight = -Height;
	jpgProps.DIBBytes =	(unsigned char *)pData;
	jpgProps.DIBPadBytes = 0;
	jpgProps.DIBChannels = 4;
	jpgProps.DIBColor =	IJL_RGB;

	jpgProps.JPGFile = filename;
	jpgProps.JPGWidth =	Width;
	jpgProps.JPGHeight = Height;
	jpgProps.JPGChannels = 3;
	jpgProps.JPGColor =	IJL_YCBCR;
	jpgProps.JPGSubsampling	= IJL_411;
    jpgProps.jquality =	jpgQuality;	   
  
	jErr = ijlWrite(&jpgProps, IJL_JFILE_WRITEWHOLEIMAGE);

  if (jErr != IJL_OK)
  {
    ijlFree(&jpgProps);
    return false;
  }

  ijlFree(&jpgProps);

    return true;
}
**/
//-------------------------------------------------------------------------

void CSimpleGraph::GetFileName(char* sFileName)
{
#ifdef WIN32

	char szFile[1024], szFileTitle[1024];
	szFile[0] = '\0', szFileTitle[0] = '\0';

	OPENFILENAME ofn;
    ofn.lStructSize = sizeof(OPENFILENAME); 
    ofn.hwndOwner = NULL; //GetActiveWindow(); //NULL; //hWnd; 
    ofn.lpstrCustomFilter = NULL;
    ofn.nFilterIndex = 0;
    ofn.lpstrFile = szFile; 
    ofn.nMaxFile = sizeof(szFile)/sizeof(*szFile); 
    ofn.lpstrFileTitle = szFileTitle; 
    ofn.lpstrInitialDir = (LPSTR)NULL; 
    ofn.lpstrTitle = NULL; //"Enter File Name"; //szTitle; 
    ofn.Flags = OFN_EXPLORER | OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST; //OFN_SHOWHELP | OFN_OVERWRITEPROMPT | OFN_EXPLORER; 
    ofn.lpstrDefExt = "png";
    ofn.lpstrFilter = "PNG Files (*.png)\0*.png\0BMP Files (*.bmp)\0*.bmp\0\0";
    //ofn.lpstrFilter = "PNG Files (*.png)\0*.png\0BMP Files (*.bmp)\0*.bmp\0All Files (*.*)\0*.*\0";
    ofn.nMaxFileTitle = sizeof(szFileTitle); 

	if(!GetSaveFileName(&ofn)) return;
	strcpy(sFileName, szFile);

/**
   fn.lpstrFile         = DefaultFile;
   fn.nMaxFile          = sizeof(DefaultFile);
   fn.lpstrFileTitle    = NULL;
   fn.lpstrInitialDir   = NULL;
   fn.lpstrTitle        = NULL;
   fn.Flags             = OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT |
                          OFN_PATHMUSTEXIST;
   fn.lpstrDefExt       = DefaultExt;
   // construct file name filter string and default file name string
   lstrcpy(DefaultFile,FileName);
   if (lstrcmp(DefaultExt,"bmp") == 0)
   {
      fn.lpstrFilter = "BMP files (*.bmp)\0*.BMP\0\0";
      lstrcpy(strrchr(DefaultFile,'.'),".bmp");
   }
   else if (lstrcmp(DefaultExt,"pcx") == 0)
   {
      fn.lpstrFilter = "PCX files (*.pcx)\0*.PCX\0\0";
      lstrcpy(strrchr(DefaultFile,'.'),".pcx");
   }
   else
   {
      fn.lpstrFilter = NULL;
   }
   // activate the Open File dialog box
   if (GetSaveFileName(&fn))
   {
      lstrcpy(FileName,DefaultFile);
      return OK;
   }
   else
      return ERR;
**/

#endif
}

//-------------------------------------------------------------------------

void CSimpleGraph::RenderImageToSave(int WinID, int SaveWidth, int SaveHeight)
{
	CHWinCont hCurCont;
	if(!GetCurHWinCont(hCurCont)) return;
	if(hCurCont.rep == 0) return;

	hCurCont.rep->SetViewingParamsGL(); //because it enforces black image when saving without off-screen rendering
	hCurCont.rep->ProcReshapeContent(SaveWidth, SaveHeight);
	hCurCont.rep->ProcDisplayContent(); //to ensure continious update at resizing
	//hCurCont.rep->ProcReshapeContent(SaveWidth, SaveHeight);

	glFinish(); // This is very important!!! Make sure buffered commands are finished!!!
}

//-------------------------------------------------------------------------
