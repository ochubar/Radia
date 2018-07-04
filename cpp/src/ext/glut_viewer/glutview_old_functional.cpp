
#include "glutview.h"

#include <math.h>
#include <stdlib.h>

//#ifdef GLUT_BUILDING_LIB
//#undef GLUT_BUILDING_LIB //to allow compilimg and linking with glut.h
//#endif
//#ifdef _DLL
//#undef _DLL
//#endif

#include <GL/glut.h>

#ifndef __VIEWER3D_H
#include "viewer3d.h"
#endif

/* Grid */
enum {WIREFRAME, HIDDENLINE, FLATSHADED, SMOOTHSHADED, TEXTURED};
enum {FULLSCREEN, FACENORMALS, ANTIALIAS, ENVMAP};
//enum {WEAK, NORMAL, STRONG};
enum {SMALL, MEDIUM, LARGE, XLARGE};
//enum {CURRENT, FLAT, SPIKE, DIAGONALWALL, SIDEWALL, HOLE, 
//      MIDDLEBLOCK, DIAGONALBLOCK, CORNERBLOCK, HILL, HILLFOUR};
int displayMode = WIREFRAME;
//int resetMode = DIAGONALBLOCK;

//#define MAXGRID 5 //63

//int grid = 17;
//float dt = 0.004;

bool waving = false, editing = false, 
     drawFaceNorms = false, antialias = false,
     envMap = false;
//#define SQRTOFTWOINV 1.0 / 1.414213562

/* Some <math.h> files do not define M_PI... */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Viewing */
float gPhi=270.0, gTheta=90; //45.0;

float gDepth = 0;
float gDepthMotionQuanta = 0;
float gDepthNear=15.0, gDepthFar=100.0;

float aspect = 5.0/4.0;
long xsize, ysize;
int downX, downY;
bool leftButton = false, middleButton = false;

GLfloat lightPosition[] = { 0.0, 0.0, 1.0, 1.0}; 
//int displayMenu, otherMenu, sizeMenu, mainMenu;

//-------------------------------------------------------------------------

double* gPgVertCoord;
int gPgNv;
int* gPgVertInd;
int* gPgLen;
float* gPgColors;
int gPgN;
float* gPgNormCoord = 0;
bool gPgNormalsWereDefined = false;
bool gPgWereDefined = false;

double* gLnVertCoord;
int gLnNv;
int* gLnVertInd;
int* gLnLen;
float* gLnColors;
int gLnN;
bool gLnWereDefined = false;

double gCenVert[3];
double gObjLimits[6];
double gMaxSize;
int gActiveWindowID = 0;
//int gLastWinWidth = 500;
//int gLastWinHeight = 500;

//-------------------------------------------------------------------------

//void wave(void)
//{
//    if(waving)
//    {
//        //getforce();
//        //getvelocity();
//        //getposition();
//        glutPostRedisplay();
//    }
//}
//void go(void)
//{
//    waving = true;
//    editing = false;
//    glutIdleFunc(wave);
//}

//void stop(void)
//{
//    waving = false;
//    glutIdleFunc(NULL);
//}

//void edit(void)
//{
//    stop();
//    editing = true;
//}

void setSize(int value)
{
	double xDim = gObjLimits[1] - gObjLimits[0];
	double yDim = gObjLimits[3] - gObjLimits[2];
	double zDim = gObjLimits[5] - gObjLimits[4];

	double TransSize = yDim;
	if(TransSize < zDim) TransSize = zDim;

	double LocGrid;
    switch(value) 
    {
        case SMALL : LocGrid = 2*TransSize/4; break;
        case MEDIUM: LocGrid = 2*TransSize/2; break;
        case LARGE : LocGrid = 2*TransSize/1.5; break;
        case XLARGE : LocGrid = 2*TransSize; break;
    }

	LocGrid += 0.5*xDim;

    gDepthNear= LocGrid/20.0;
    gDepthFar= LocGrid*30.0;

    gDepth = (5.0/4.0)*LocGrid;
	gDepthMotionQuanta = 0.05*gDepth;

	gPhi=270.0;
	gTheta=90; //45.0;

    glutPostRedisplay();
}

void setDisplay(int value)
{
    displayMode = value;
    switch(value) 
    {
        case WIREFRAME   : 
            glShadeModel(GL_FLAT); 
            glDisable(GL_LIGHTING);
            break;
        case HIDDENLINE: 
            glShadeModel(GL_FLAT); 
            glDisable(GL_LIGHTING);
            break;
        case FLATSHADED  : 
            glShadeModel(GL_FLAT); 
            glEnable(GL_LIGHTING);
            break;
        case SMOOTHSHADED: 
            glShadeModel(GL_SMOOTH); 
            glEnable(GL_LIGHTING);
            break;
        case TEXTURED: 
            glShadeModel(GL_SMOOTH); 
            glEnable(GL_LIGHTING);
            break;
    }
    glutPostRedisplay();
}

void setOther(int value)
{
    switch (value)
    {
        case FULLSCREEN: 
            glutFullScreen();
            break;
        case FACENORMALS: 
            drawFaceNorms = !drawFaceNorms;
            break;
        case ANTIALIAS: 
            antialias = !antialias;
            if (antialias)
            {
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                glEnable(GL_LINE_SMOOTH);
                glLineWidth(1.5);
            }
            else
            {
                glDisable(GL_BLEND);
                glDisable(GL_LINE_SMOOTH);
                glLineWidth(1.0);
            }
            break;
        case ENVMAP: 
            envMap = !envMap;
            if (envMap)
            {
                //glBindTexture(GL_TEXTURE_2D, texId2);
                glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
                glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
                glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
                glEnable(GL_TEXTURE_GEN_S);
                glEnable(GL_TEXTURE_GEN_T);
            }
            else
            {
                //glBindTexture(GL_TEXTURE_2D, texId1);
                glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
                glDisable(GL_TEXTURE_GEN_S);
                glDisable(GL_TEXTURE_GEN_T);
            }
            break;
    }
    glutPostRedisplay();
}

void reshape(int width, int height)
{
    xsize = width; 
    ysize = height;
    aspect = (float)xsize/(float)ysize;
    glViewport(0, 0, xsize, ysize);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glutPostRedisplay();
}

//-------------------------------------------------------------------------

bool CheckIfVectorsAreCollinear(double xV1, double yV1, double zV1, double xV2, double yV2, double zV2)
{
	const double RelPrec = 1.e-5;

	double V1scV2 = xV1*xV2 + yV1*yV2 + zV1*zV2;
	double V1scV2e2 = V1scV2*V1scV2;

	double V1e2 = xV1*xV1 + yV1*yV1 + zV1*zV1;
	double V2e2 = xV2*xV2 + yV2*yV2 + zV2*zV2;
	double V1e2V2e2 = V1e2*V2e2;

	if(fabs(fabs(V1scV2e2) - fabs(V1e2V2e2)) < RelPrec*V1e2V2e2) return true;
	return false;
}

//-------------------------------------------------------------------------

void DrawArbPolygObjects(GLenum Mode, double* VertCoord, int Nv, int* VertInd, int* Lengths, float* Colors, int Npg, float* NormCoord, bool& NormalsWereDefined)
{
	bool ColorsShouldBeTreated = ((Colors != 0) && (Npg > 0));
	bool NormalsShouldBeTreated = ((NormCoord != 0) && (Npg > 0));	

	float xP1, yP1, zP1, xP2, yP2, zP2, xN, yN, zN;
	int *tVertInd = VertInd;

	for(int i=0; i<Npg; i++)
	{
		int AmOfPtInPg = Lengths[i];
		int IndPgColorAndNormCoord = i*3;

		if(ColorsShouldBeTreated)
		{
			float* tColor = Colors + IndPgColorAndNormCoord;
			float r = *(tColor++);
			float g = *(tColor++);
			float b = *tColor;
			if((r >= 0.) && (g >= 0.) && (b >= 0.))
			{
				glColor3f(r, g, b);
			}
		}

		glBegin(Mode);

		bool FirstPointDefined = false;
		bool SecondPointDefined = false;
		bool ThirdPointDefined = false;

		for(int j=0; j<AmOfPtInPg; j++)
		{
			int VertInd = *(tVertInd++);
			double *tVertCoord = VertCoord + ((VertInd - 1)*3);

			float xVert = *(tVertCoord++);
			float yVert = *(tVertCoord++);
			float zVert = *tVertCoord;

			if(!FirstPointDefined)
			{
				xP1 = xVert; yP1 = yVert; zP1 = zVert; 
				FirstPointDefined = true;
			}
			else if(!SecondPointDefined)
			{
				xP2 = xVert; yP2 = yVert; zP2 = zVert; 
				SecondPointDefined = true;

				if(AmOfPtInPg == 2)
				{
					glVertex3f(xP1, yP1, zP1);
					glVertex3f(xP2, yP2, zP2);
				}
			}
			else if(!ThirdPointDefined)
			{
				bool CurNormalIsDefined = false;

				if(NormalsShouldBeTreated && (!NormalsWereDefined))
				{
					double xV1 = xP2 - xP1, yV1 = yP2 - yP1, zV1 = zP2 - zP1;
					double xV2 = xVert - xP2, yV2 = yVert - yP2, zV2 = zVert - zP2;

					if(!CheckIfVectorsAreCollinear(xV1, yV1, zV1, xV2, yV2, zV2))
					{
						xN = yV1*zV2 - zV1*yV2;
						yN = zV1*xV2 - xV1*zV2;
						zN = xV1*yV2 - yV1*xV2;

						float Ne2 = xN*xN + yN*yN + zN*zN;
						if(Ne2 == 0.) Ne2 = 1e-7;
						float InvSqrtN = 1./sqrt(Ne2);
						xN *= InvSqrtN; yN *= InvSqrtN; zN *= InvSqrtN;

						float *tNormCoord = NormCoord + IndPgColorAndNormCoord;
						*(tNormCoord++) = xN;
						*(tNormCoord++) = yN;
						*tNormCoord = zN;

						CurNormalIsDefined = true;
						ThirdPointDefined = true;
					}
				}
				else if(NormalsShouldBeTreated && NormalsWereDefined)
				{
					float *tNormCoord = NormCoord + IndPgColorAndNormCoord;
					xN = *(tNormCoord++);
					yN = *(tNormCoord++);
					zN = *tNormCoord;

					CurNormalIsDefined = true;
					ThirdPointDefined = true;
				}
				else 
				{
					ThirdPointDefined = true;
				}

				if(ThirdPointDefined)
				{
					if(CurNormalIsDefined) glNormal3f(xN, yN, zN);
					glVertex3f(xP1, yP1, zP1);
					glVertex3f(xP2, yP2, zP2);
					glVertex3f(xVert, yVert, zVert);
				}
			}
			else
			{
				glVertex3f(xVert, yVert, zVert);
			}
		}
		glEnd();
	}
	if(NormalsShouldBeTreated) NormalsWereDefined = true;
}

//-------------------------------------------------------------------------

void DrawObjects()
{
	if(gPgWereDefined)
		DrawArbPolygObjects(GL_POLYGON, gPgVertCoord, gPgNv, gPgVertInd, gPgLen, gPgColors, gPgN, gPgNormCoord, gPgNormalsWereDefined);
	//if(gLnWereDefined)
	//	DrawArbPolygObjects(GL_LINE_STRIP, gLnVertCoord, gLnNv, gLnVertInd, gLnLen, gLnColors, gLnN, 0, gPgNormalsWereDefined);
}

//-------------------------------------------------------------------------

void display(void) 
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(64.0, aspect, gDepthNear, gDepthFar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity(); 

	//improve these movements
    glTranslatef(0.0,0.0,-gDepth);
    glRotatef(-gTheta, 1.0, 0.0, 0.0);
    glRotatef(gPhi, 0.0, 0.0, 1.0);
    glTranslatef(-gCenVert[0], -gCenVert[1], -gCenVert[2]);

	DrawObjects();

    glutSwapBuffers();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

//void visibility(int state)
//{
//    if((state == GLUT_VISIBLE) && waving) go();
//    else stop();
//}

void keyboard(unsigned char ch, int x, int y)
{
    switch (ch) 
    {
        case '-': 
			if(gDepth < gDepthFar) gDepth += gDepthMotionQuanta; 
			break;
        case '+': 
			if(gDepth > gDepthNear) gDepth -= gDepthMotionQuanta; 
			break;
        case 27: exit(0); break;
        //case 27: { glutDestroyWindow(gActiveWindowID); gActiveWindowID = 0; glutHackStopMainLoop();} break;
    }
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
    downX = x;
    downY = y;
    leftButton = ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN));
    middleButton = ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN));
}

void motion(int x, int y)
{//to improve
    if(leftButton)
    {
        gTheta += (float)(0.25*(downY - y));

		if(gTheta >= 180.) gTheta = gTheta - 360.;
		else if(gTheta < -180.) gTheta = gTheta + 360.;

		double PhiMult = 0.25;
		if(gTheta <= 0.) PhiMult = -PhiMult;

        gPhi += (float)(PhiMult*(x - downX));

		//int Aha = 1;
    }
    if(middleButton)
    {
        gDepth += (float)(downY - y) / 10.0;
    }
    downX = x;
    downY = y;
    glutPostRedisplay();
}

//void setMain(int value)
//{
//    switch(value) 
//    {
//        case 1: edit();    break;
//        case 2: go();      break; /* set idle func to something */
//        case 3: stop();    break; /* set idle func to null */
//        //case 4: reverse(); break;
//        case 5: exit(0);   break;
//        //case 5: 
//			//{ glutDestroyWindow(gActiveWindowID); gActiveWindowID = 0; glutHackStopMainLoop();} 
//			//break;
//        case 6: //Full Screen
//			{
//				gLastWinWidth = glutGet(GLUT_WINDOW_WIDTH);
//				gLastWinHeight = glutGet(GLUT_WINDOW_HEIGHT);
//
//				glutRemoveMenuItem(2); 
//				glutRemoveMenuItem(2); 
//
//				glutAddMenuEntry("Window View", 7);
//				glutAddMenuEntry("Exit", 5);
//
//				glutFullScreen(); 
//			}
//			break;
//		case 7: //Window View
//			{ 
//				glutRemoveMenuItem(2); 
//				glutRemoveMenuItem(2); 
//
//				glutAddMenuEntry("Full Screen", 6);
//				glutAddMenuEntry("Exit", 5);
//
//				glutReshapeWindow(gLastWinWidth, gLastWinHeight);
//				glutPositionWindow(10, 30); 
//			}
//			break;
//    }
//}

//-------------------------------------------------------------------------

int StartViewer()
{
	//try 
	//{ //GLUT_DISABLE_ATEXIT_HACK,
	//	if(gActiveWindowID > 0) 
	//	{
	//		glutDestroyWindow(gActiveWindowID);
	//		gActiveWindowID = 0;
	//	}
	//}
	//catch(...) { gActiveWindowID = 0;}
 
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    gActiveWindowID = glutCreateWindow("3D Viewer");
    
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glPolygonOffset(1.0, 1.0);
    glEnable(GL_CULL_FACE);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glLightfv (GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_LIGHT0);
    
    setSize(MEDIUM);
    setDisplay(SMOOTHSHADED);
    setOther(ENVMAP);

	glutReshapeFunc(reshape);
    glutDisplayFunc(display);
//    glutVisibilityFunc(visibility);

	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);

    glutMainLoop();

    return 0;             /* ANSI C requires main to return int. */
}

//-------------------------------------------------------------------------

void DetermineAverageCoordAndLimits(double* VertCoord, int Nv, double* LnVertCoord, int LnNv, double* CenVert, double* ObjLimits, double& MaxSize)
{
	double &xMin = ObjLimits[0], &xMax = ObjLimits[1];
	double &yMin = ObjLimits[2], &yMax = ObjLimits[3];
	double &zMin = ObjLimits[4], &zMax = ObjLimits[5];
	xMin = yMin = zMin = 1e+23;
	xMax = yMax = zMax = -1e+23;
	MaxSize = 0;

	if(((VertCoord == 0) || (Nv == 0)) && ((LnVertCoord == 0) || (LnNv == 0))) return;

	double xSum = 0, ySum = 0, zSum = 0;
	int TotNp = 0;

	if((VertCoord != 0) && (Nv > 0))
	{
		double *tVertCoord = VertCoord;
		for(int i=0; i<Nv; i++)
		{
			double x = *(tVertCoord++);
			double y = *(tVertCoord++);
			double z = *(tVertCoord++);
			
			if(xMin > x) xMin = x;
			if(xMax < x) xMax = x;
			if(yMin > y) yMin = y;
			if(yMax < y) yMax = y;
			if(zMin > z) zMin = z;
			if(zMax < z) zMax = z;
			
			xSum += x;
			ySum += y;
			zSum += z;
		}
		TotNp += Nv;
	}
	if((LnVertCoord != 0) && (LnNv > 0))
	{
		double *tVertCoord = LnVertCoord;
		for(int i=0; i<LnNv; i++)
		{
			double x = *(tVertCoord++);
			double y = *(tVertCoord++);
			double z = *(tVertCoord++);
			
			if(xMin > x) xMin = x;
			if(xMax < x) xMax = x;
			if(yMin > y) yMin = y;
			if(yMax < y) yMax = y;
			if(zMin > z) zMin = z;
			if(zMax < z) zMax = z;
			
			xSum += x;
			ySum += y;
			zSum += z;
		}
		TotNp += LnNv;
	}

	//double InvNv = 1./TotNp;
	//CenVert[0] = xSum*InvNv;
	//CenVert[1] = ySum*InvNv;
	//CenVert[2] = zSum*InvNv;

	CenVert[0] = 0.5*(xMax + xMin);
	CenVert[1] = 0.5*(yMax + yMin);
	CenVert[2] = 0.5*(zMax + zMin);

	double xDim = xMax - xMin;
	double yDim = yMax - yMin;
	double zDim = zMax - zMin;

	if(MaxSize < xDim) MaxSize = xDim;
	if(MaxSize < yDim) MaxSize = yDim;
	if(MaxSize < zDim) MaxSize = zDim;
}

//-------------------------------------------------------------------------

int CALL ViewPolygonsOpenGL(double* VertCoord, int Nv, int* PgVertInd, int* PgLen, float* PgColors, int Npg)
{
	gPgVertCoord = VertCoord;
	gPgNv = Nv;
	gPgVertInd = PgVertInd;
	gPgLen = PgLen;
	gPgColors = PgColors;
	gPgN = Npg;
	if(gPgN > 0) gPgNormCoord = new float[3*gPgN];
	gPgNormalsWereDefined = false;
	if(PgColors != 0) gPgColors = PgColors;
	gPgWereDefined = true;

	gLnVertCoord = 0;
	gLnNv = 0;
	gLnVertInd = 0;
	gLnLen = 0;
	gLnN = 0;
	gLnWereDefined = false;

	DetermineAverageCoordAndLimits(gPgVertCoord, gPgNv, 0, 0, gCenVert, gObjLimits, gMaxSize);

	int res = StartViewer();

	if(gPgWereDefined && (gPgNormCoord != 0)) { delete[] gPgNormCoord; gPgNormCoord = 0;}
	return res;
}

//-------------------------------------------------------------------------

int CALL ViewSceneOpenGL(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr)
{
	if((VertCoordArr == 0) || (NvArr == 0) || (VertIndArr == 0) || (LenArr == 0) || (NpArr == 0)) return 0;

    CHWinCont hWinCont(new CWinContViewer3D(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr));
    CSimpleGraph::StartNewWinApp(hWinCont);
	return 0;



	//gPgVertCoord = 0;
	//gPgNv = 0;
	//gPgVertInd = 0;
	//gPgLen = 0;
	//gPgN = 0;
	//gPgNormCoord = 0;
	//gPgWereDefined = false;

	//gLnVertCoord = 0;
	//gLnNv = 0;
	//gLnVertInd = 0;
	//gLnLen = 0;
	//gLnN = 0;
	//gLnWereDefined = false;

	//if((VertCoordArr[0] != 0) && (NvArr[0] > 0) && (VertIndArr[0] != 0) && (LenArr[0] != 0) && (NpArr[0] != 0))
	//{// closed polygons
	//	gPgVertCoord = VertCoordArr[0];
	//	gPgNv = NvArr[0];
	//	gPgVertInd = VertIndArr[0];
	//	gPgLen = LenArr[0];
	//	gPgN = NpArr[0];
	//	gPgNormCoord = new float[3*gPgN];
	//	gPgNormalsWereDefined = false;
	//	if((ColorsArr != 0) && (ColorsArr[0] != 0)) gPgColors = ColorsArr[0];
	//	gPgWereDefined = true;
	//}

	//if((VertCoordArr[1] != 0) && (NvArr[1] > 0) && (VertIndArr[1] != 0) && (LenArr[1] != 0) && (NpArr[1] != 0))
	//{// lines
	//	gLnVertCoord = VertCoordArr[1];
	//	gLnNv = NvArr[1];
	//	gLnVertInd = VertIndArr[1];
	//	gLnLen = LenArr[1];
	//	gLnN = NpArr[1];
	//	if((ColorsArr != 0) && (ColorsArr[1] != 0)) gLnColors = ColorsArr[1];
	//	gLnWereDefined = true;

	//	//char ErrorMesTitle[] = "SRW Debug";
 //       //char ErrorStr[100];
 //       //sprintf(ErrorStr, "%d", gLnNv);
 //       //UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
 //       //int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
	//}

	//DetermineAverageCoordAndLimits(gPgVertCoord, gPgNv, gLnVertCoord, gLnNv, gCenVert, gObjLimits, gMaxSize);

	//int res = StartViewer();


			////int res = 0;
			////CHWinCont hWinCont(new CWinContViewer3D(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr));
			////CSimpleGraph::StartNewWinApp(hWinCont);


	//if(gPgWereDefined) { delete[] gPgNormCoord; gPgNormCoord = 0;}
	//return res;

}

//-------------------------------------------------------------------------

