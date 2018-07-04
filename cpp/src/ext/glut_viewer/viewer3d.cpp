
#ifndef __VIEWER3D_H
#include "viewer3d.h"
#endif
#ifndef __GMMETH_H
#include "gmmeth.h"
#endif

#include <math.h>

//-------------------------------------------------------------------------

const int CWinContViewer3D::sm_SMALL = 0;
const int CWinContViewer3D::sm_MEDIUM = 1;
const int CWinContViewer3D::sm_LARGE = 2;
const int CWinContViewer3D::sm_XLARGE = 3;

const int CWinContViewer3D::sm_ROTATION3D = 0;
const int CWinContViewer3D::sm_ROTATION2D = 1;
const int CWinContViewer3D::sm_TRANSLATION_DEPTH = 2;

const int CWinContViewer3D::sm_WIREFRAME = 0;
const int CWinContViewer3D::sm_HIDDENLINE = 1;
const int CWinContViewer3D::sm_FLATSHADED = 2;
const int CWinContViewer3D::sm_SMOOTHSHADED = 3;
const int CWinContViewer3D::sm_TEXTURED = 4;

const int CWinContViewer3D::sm_FULLSCREEN = 0;
const int CWinContViewer3D::sm_FACENORMALS = 1;
const int CWinContViewer3D::sm_ANTIALIAS = 2;
const int CWinContViewer3D::sm_ENVMAP = 3;

TVector3d CWinContViewer3D::m_Vz(0, 0, 1);
TVector3d CWinContViewer3D::m_ZeroV(0, 0, 0);

//-------------------------------------------------------------------------

CWinContViewer3D::CWinContViewer3D(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr, char* WinTitle)
{
	InitViewingParams();
	if((VertCoordArr == 0) || (NvArr == 0) || (VertIndArr == 0) || (LenArr == 0) || (NpArr == 0)) return;

	double *LocPgVertCoord = 0, *LocLnVertCoord = 0;
	int LocPgNv = 0, LocLnNv = 0;
	if((VertCoordArr[0] != 0) && (NvArr[0] > 0))
	{
        LocPgVertCoord = VertCoordArr[0]; LocPgNv = NvArr[0];
	}
	if((VertCoordArr[1] != 0) && (NvArr[1] > 0))
	{
        LocLnVertCoord = VertCoordArr[1]; LocLnNv = NvArr[1];
	}
	double LocCenVert[3], LocObjLimits[6], LocMaxSize;
	DetermineAverageCoordAndLimits(LocPgVertCoord, LocPgNv, LocLnVertCoord, LocLnNv, LocCenVert, LocObjLimits, LocMaxSize);

	SetupPgAndLnData(VertCoordArr, NvArr, VertIndArr, LenArr, ColorsArr, NpArr, LocCenVert, LocObjLimits, LocMaxSize);
    //CalculatePolygonNormals();

	DoSetSizeContent(sm_MEDIUM); // to setup m_Depth
	SetupRotInit();
	SetupTrfCur(m_RotInit);

    if(WinTitle != 0) 
	{
		int LenWinTitle = (int)strlen(WinTitle);
		if(LenWinTitle > 500) LenWinTitle = 500;
		strncpy(m_WinTitle, WinTitle, LenWinTitle);
		m_WinTitle[LenWinTitle] = '\0';
	}
	else strcpy(m_WinTitle, CSimpleGraph::sm_DefWinTitle);
}

//-------------------------------------------------------------------------

void CWinContViewer3D::CalculatePolygonNormals()
{
    m_PgNormalsWereDefined = false;
	if((m_PgVertCoord == 0) || (m_PgVertInd == 0) || (m_PgLen == 0) || (m_PgNormCoord == 0) || (m_PgNv <= 0) || (m_PgN <= 0)) return;

	m_PgNormalsWereDefined = false;

	//const double RelPrec = 1.e-05;
	const double RelPrec = 1.e-09; //OC050109

	double xP1, yP1, zP1, xP2, yP2, zP2, xN, yN, zN;
	int *tVertInd = m_PgVertInd;
	//float *tNormCoord = m_PgNormCoord;
	double *tNormCoord = m_PgNormCoord; //OC050109

	for(long i=0; i<m_PgN; i++)
	{
        int AmOfPtInPg = m_PgLen[i];

		bool FirstPointDefined = false;
		bool SecondPointDefined = false;
		bool ThirdPointDefined = false;

		for(long j=0; j<AmOfPtInPg; j++)
		{
			int VertInd = *(tVertInd++);
			//double *tVertCoord = pViewCont->m_PgVertCoord + ((VertInd - 1)*3);
			double *tVertCoord = m_PgVertCoord + (VertInd*3); //to check

			double xVert = *(tVertCoord++);
			double yVert = *(tVertCoord++);
			double zVert = *tVertCoord;

			if(!FirstPointDefined)
			{
				xP1 = xVert; yP1 = yVert; zP1 = zVert; 
				FirstPointDefined = true;
			}
			else if(!SecondPointDefined)
			{
				xP2 = xVert; yP2 = yVert; zP2 = zVert; 
				SecondPointDefined = true;
			}
			else if(!ThirdPointDefined)
			{
                double xV1 = xP2 - xP1, yV1 = yP2 - yP1, zV1 = zP2 - zP1;
				double xV2 = xVert - xP2, yV2 = yVert - yP2, zV2 = zVert - zP2;
				
				if(!CGenMathMeth::VectCheckIfCollinear(xV1, yV1, zV1, xV2, yV2, zV2, RelPrec))
				{
                    xN = yV1*zV2 - zV1*yV2;
                    yN = zV1*xV2 - xV1*zV2;
                    zN = xV1*yV2 - yV1*xV2;

                    double Ne2 = xN*xN + yN*yN + zN*zN;
                    if(Ne2 == 0.) Ne2 = RelPrec*RelPrec;

                    double InvSqrtN = 1./sqrt(Ne2);
                    xN *= InvSqrtN; yN *= InvSqrtN; zN *= InvSqrtN;
                    //*(tNormCoord++) = (float)xN;
                    //*(tNormCoord++) = (float)yN;
                    //*(tNormCoord++) = (float)zN;
                    *(tNormCoord++) = xN; //OC050109
                    *(tNormCoord++) = yN;
                    *(tNormCoord++) = zN;

                    ThirdPointDefined = true;
					break;
				}
			}
		}
		if(!ThirdPointDefined)
		{
            *(tNormCoord++) = 0;
            *(tNormCoord++) = 0;
			*(tNormCoord++) = 0;
		}
	}
	m_PgNormalsWereDefined = true;
}

//-------------------------------------------------------------------------

void CWinContViewer3D::DetermineAverageCoordAndLimits(double* VertCoord, int Nv, double* LnVertCoord, int LnNv, double* CenVert, double* ObjLimits, double& MaxSize)
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

void CWinContViewer3D::ProcReshapeContent(int Width, int Height) // virtual
{
	m_xSize = Width;
    m_ySize = Height;
    //m_Aspect = (float)m_xSize/(float)m_ySize; //OCtest
    m_Aspect = (float)((double)m_xSize/(double)m_ySize); //OCtest
	glViewport(0, 0, m_xSize, m_ySize); //OCtest

	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glutPostRedisplay();
}

//-------------------------------------------------------------------------

void CWinContViewer3D::ProcMouseContent(int button, int state, int x, int y) // virtual
{
	m_DownX = x;
	m_DownY = y;
	m_LeftButton = ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN));
	m_MiddleButton = ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN));
	//glutPostRedisplay();
}

//-------------------------------------------------------------------------

void CWinContViewer3D::ProcMotionContent(int x, int y) // virtual
{//to improve
	const double MultTransvFact = 1.3;
	const double MultLongitFact = 0.5;
	if(m_LeftButton)
	{
		//m_Theta += (float)(0.25*(m_DownY - y));
		//if(m_Theta >= 180.) m_Theta = m_Theta - 360.;
		//else if(m_Theta < -180.) m_Theta = m_Theta + 360.;
		//double PhiMult = 0.25;
		//if(m_Theta <= 0.) PhiMult = -PhiMult;
		//m_Phi += (float)(PhiMult*(x - m_DownX));

		if((x != m_DownX) || (y != m_DownY))
		{
			if(m_MotionType == sm_ROTATION3D)
			{
				//double ApproxTransvFactY = MultTransvFact*m_Depth/((float)m_ySize);
				double ApproxTransvFactY = MultTransvFact*m_Depth/((double)m_ySize);

				//double ApproxTransvFactX = ApproxTransvFactY*m_Aspect;
				double ApproxTransvFactX = ApproxTransvFactY;
				double RelObjCoordX = (x - 0.5*m_xSize)*ApproxTransvFactX;
				double RelObjCoordY = (0.5*m_ySize - y)*ApproxTransvFactY;
				double PrevRelObjCoordX = (m_DownX - 0.5*m_xSize)*ApproxTransvFactX;
				double PrevRelObjCoordY = (0.5*m_ySize - m_DownY)*ApproxTransvFactY;

				TVector3d ProjVect(RelObjCoordX, RelObjCoordY, m_Depth); 
				TVector3d DispVect(RelObjCoordX - PrevRelObjCoordX, RelObjCoordY - PrevRelObjCoordY, 0); 

				//TVector3d DispVect(x - m_DownX, m_DownY - y, 0); 
				//TVector3d NewAxVect = m_Vz^DispVect;
				//double ExtraAngle = AngleMultQuantaRad*DispVect.Abs();

				TVector3d NewAxVect = ProjVect^DispVect;
				double ExtraAngle = DispVect.Abs()/(::fabs(m_Depth));

				m_TrfAux1.SetupRotation(m_ZeroV, NewAxVect, ExtraAngle);
				m_TrfAux2 = m_RotCur;

				TrProduct(&m_TrfAux1, &m_TrfAux2, m_RotCur);
				SetupTrfCur(m_RotCur);
			}
			else if(m_MotionType == sm_ROTATION2D)
			{
				//double ApproxTransvFactY = MultTransvFact*m_Depth/((float)m_ySize);
				double ApproxTransvFactY = MultTransvFact*m_Depth/((double)m_ySize);

				double ApproxTransvFactX = ApproxTransvFactY;
				double RelObjCoordX = (x - 0.5*m_xSize)*ApproxTransvFactX;
				double RelObjCoordY = (0.5*m_ySize - y)*ApproxTransvFactY;
				double PrevRelObjCoordX = (m_DownX - 0.5*m_xSize)*ApproxTransvFactX;
				double PrevRelObjCoordY = (0.5*m_ySize - m_DownY)*ApproxTransvFactY;

				TVector3d ProjVect(RelObjCoordX, RelObjCoordY, 0); 
				double LenProjVect = ProjVect.Abs();
				if(LenProjVect > 0)
				{
                    TVector3d DispVect(RelObjCoordX - PrevRelObjCoordX, RelObjCoordY - PrevRelObjCoordY, 0); 
                    double ExtraAngle = DispVect.Abs()/LenProjVect;

					TVector3d TestNewAxVect = ProjVect^DispVect;
					double TestScalProd = TestNewAxVect*m_Vz;
					if(TestScalProd < 0) ExtraAngle = -ExtraAngle;

                    m_TrfAux1.SetupRotation(m_ZeroV, m_Vz, ExtraAngle);
                    m_TrfAux2 = m_RotCur;
                    TrProduct(&m_TrfAux1, &m_TrfAux2, m_RotCur);
                    SetupTrfCur(m_RotCur);
				}
			}
			else if(m_MotionType == sm_TRANSLATION_DEPTH)
			{
                //m_Depth += (float)(m_DownY - y)/10.0;
                //m_Depth += MultLongitFact*((float)(m_DownY - y));
                m_Depth += MultLongitFact*((double)(m_DownY - y));

                SetupTrfCur(m_RotCur);
			}
		}
	}
	if(m_MiddleButton)
	{
		//m_Depth += (float)(m_DownY - y)/10.0;
		//m_Depth += MultLongitFact*((float)(m_DownY - y));
		m_Depth += MultLongitFact*((double)(m_DownY - y));
        SetupTrfCur(m_RotCur);
	}

	m_DownX = x;
	m_DownY = y;
}

//-------------------------------------------------------------------------

void CWinContViewer3D::ProcKeyboardContent(unsigned char ch, int x, int y) // virtual
{
    switch(ch) 
    {
        case '-': 
			if(m_Depth < m_DepthFar) m_Depth += m_DepthMotionQuanta; 
			break;
        case '+': 
			if(m_Depth > m_DepthNear) m_Depth -= m_DepthMotionQuanta; 
			break;
        //case 27: exit(0); break;
        case 27: 
			{ 
				//glutDestroyWindow(gActiveWindowID); gActiveWindowID = 0; glutHackStopMainLoop();
				int CurWinID = glutGetWindow();
				if(CurWinID != 0) CSimpleGraph::DestroyWindow(CurWinID);
			} 
			break;
    }
}

//-------------------------------------------------------------------------

void CWinContViewer3D::ProcDisplayContent() // virtual
{//update display of Viewer3D

//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(64.0, m_Aspect, m_DepthNear, m_DepthFar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity(); 

	//improve these movements
    //glTranslatef(0.0,0.0,-m_Depth);
    //glRotatef(-m_Theta, 1.0, 0.0, 0.0);
    //glRotatef(m_Phi, 0.0, 0.0, 1.0);
    //glTranslatef(-m_CenVert[0], -m_CenVert[1], -m_CenVert[2]);

	GLfloat ArrMatrGL[16];
	SetupViewMatrArrGL(ArrMatrGL);
    glLoadMatrixf(ArrMatrGL);

	bool DummyNormalsWereDefined = false;
	if(m_PgWereDefined)
		DrawArbPolygObjects(GL_POLYGON, m_PgVertCoord, m_PgNv, m_PgVertInd, m_PgLen, m_PgColors, m_PgN, m_PgNormCoord, m_PgNormalsWereDefined);
	if(m_LnWereDefined)
        DrawArbPolygObjects(GL_LINE_STRIP, m_LnVertCoord, m_LnNv, m_LnVertInd, m_LnLen, m_LnColors, m_LnN, 0, DummyNormalsWereDefined);

//    glutSwapBuffers();
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

//-------------------------------------------------------------------------

//void CWinContViewer3D::DrawArbPolygObjects(GLenum Mode, double* VertCoord, int Nv, int* VertInd, int* Lengths, float* Colors, int Npg, float* NormCoord, bool& NormalsWereDefined)
void CWinContViewer3D::DrawArbPolygObjects(GLenum Mode, double* VertCoord, int Nv, int* VertInd, int* Lengths, float* Colors, int Npg, double* NormCoord, bool& NormalsWereDefined) //OC050109
{
	//const double RelPrec = 1.e-05;
	const double RelPrec = 1.e-09; //OC050109

	bool ColorsShouldBeTreated = ((Colors != 0) && (Npg > 0));
	bool NormalsShouldBeTreated = ((NormCoord != 0) && (Npg > 0));	

	double xP1, yP1, zP1, xP2, yP2, zP2, xN, yN, zN;
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
			else glColor3f(1, 1, 1); //OC101008 re-setting default color
		}

		glBegin(Mode);

		bool FirstPointDefined = false;
		bool SecondPointDefined = false;
		bool ThirdPointDefined = false;

		for(int j=0; j<AmOfPtInPg; j++)
		{
			int VertInd = *(tVertInd++);
			//double *tVertCoord = VertCoord + ((VertInd - 1)*3);
			double *tVertCoord = VertCoord + (VertInd*3);

			double xVert = *(tVertCoord++);
			double yVert = *(tVertCoord++);
			double zVert = *tVertCoord;

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
					//glVertex3f((float)xP1, (float)yP1, (float)zP1);
					//glVertex3f((float)xP2, (float)yP2, (float)zP2);
					glVertex3d(xP1, yP1, zP1); //OC050109
					glVertex3d(xP2, yP2, zP2);
				}
			}
			else if(!ThirdPointDefined)
			{
				bool CurNormalIsDefined = false;

				if(NormalsShouldBeTreated && (!NormalsWereDefined))
				{
					double xV1 = xP2 - xP1, yV1 = yP2 - yP1, zV1 = zP2 - zP1;
					double xV2 = xVert - xP2, yV2 = yVert - yP2, zV2 = zVert - zP2;

					if(!CGenMathMeth::VectCheckIfCollinear(xV1, yV1, zV1, xV2, yV2, zV2, RelPrec))
					{
						xN = yV1*zV2 - zV1*yV2;
						yN = zV1*xV2 - xV1*zV2;
						zN = xV1*yV2 - yV1*xV2;

						double Ne2 = xN*xN + yN*yN + zN*zN;
						if(Ne2 == 0.) Ne2 = 1e-7;
						double InvSqrtN = 1./sqrt(Ne2);
						xN *= InvSqrtN; yN *= InvSqrtN; zN *= InvSqrtN;

						//float *tNormCoord = NormCoord + IndPgColorAndNormCoord;
						double *tNormCoord = NormCoord + IndPgColorAndNormCoord; //OC050109
						*(tNormCoord++) = xN; //OC050109
						*(tNormCoord++) = yN;
						*tNormCoord = zN;

						CurNormalIsDefined = true;
						ThirdPointDefined = true;
					}
				}
				else if(NormalsShouldBeTreated && NormalsWereDefined)
				{
					//float *tNormCoord = NormCoord + IndPgColorAndNormCoord;
					double *tNormCoord = NormCoord + IndPgColorAndNormCoord; //OC050109
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
					//if(CurNormalIsDefined) glNormal3f((float)xN, (float)yN, (float)zN);
					//glVertex3f((float)xP1, (float)yP1, (float)zP1);
					//glVertex3f((float)xP2, (float)yP2, (float)zP2);
					//glVertex3f((float)xVert, (float)yVert, (float)zVert);
					if(CurNormalIsDefined) glNormal3d(xN, yN, zN); //OC050109
					glVertex3d(xP1, yP1, zP1);
					glVertex3d(xP2, yP2, zP2);
					glVertex3d(xVert, yVert, zVert);
				}
			}
			else
			{
				//glVertex3f((float)xVert, (float)yVert, (float)zVert);
				glVertex3d(xVert, yVert, zVert); //OC050109
			}
		}
		glEnd();
	}
	if(NormalsShouldBeTreated) NormalsWereDefined = true;
}

//-------------------------------------------------------------------------

//void CWinContViewer3D::DrawArbPolygObjects(GLenum Mode, double* VertCoord, int Nv, int* VertInd, int* Lengths, float* Colors, int Npg, float* NormCoord, bool NormalsWereDefined)
//{
//	const double RelPrec = 1.e-05;
//
//	bool ColorsShouldBeTreated = ((Colors != 0) && (Npg > 0));
//	bool NormalsShouldBeTreated = ((NormCoord != 0) && (Npg > 0));	
//
//	//float xP1, yP1, zP1, xP2, yP2, zP2, xN, yN, zN;
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
//		//bool FirstPointDefined = false;
//		//bool SecondPointDefined = false;
//		//bool ThirdPointDefined = false;
//
//		if(NormalsWereDefined) 
//		{
//			float* tNormCoord = NormCoord + IndPgColorAndNormCoord;
//			float xN = *(tNormCoord++);
//			float yN = *(tNormCoord++);
//			float zN = *tNormCoord;
//			if((xN != 0.) || (yN != 0.) || (zN != 0.))
//			{
//                glNormal3f(xN, yN, zN);
//			}
//		}
//
//		for(int j=0; j<AmOfPtInPg; j++)
//		{
//			int VertInd = *(tVertInd++);
//			//double *tVertCoord = VertCoord + ((VertInd - 1)*3);
//			double *tVertCoord = VertCoord + (VertInd*3);
//
//			float xVert = *(tVertCoord++);
//			float yVert = *(tVertCoord++);
//			float zVert = *tVertCoord;
//
//			glVertex3f(xVert, yVert, zVert);
//		}
//		glEnd();
//	}
//	//if(NormalsShouldBeTreated) NormalsWereDefined = true;
//}

//-------------------------------------------------------------------------

void CWinContViewer3D::SetViewingParamsGL() // virtual
{//to make different for Viewer and Graph

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

	GLfloat LocLightPosition[] = { 0.0, 0.0, 1.0, 1.0};
	for(int i=0; i<4; i++) LocLightPosition[i] = m_LightPosition[i];
	glLightfv(GL_LIGHT0, GL_POSITION, LocLightPosition);
	glEnable(GL_LIGHT0);

	CSimpleGraph::DoSetSize(CWinContViewer3D::sm_MEDIUM);
	SetDisplayMode(CWinContViewer3D::sm_SMOOTHSHADED);
	SetOtherParams(CWinContViewer3D::sm_ENVMAP);
}

//-------------------------------------------------------------------------

void CWinContViewer3D::SetupMenuGLUT()
{
	m_MenuMotion = glutCreateMenu(CSimpleGraph::DoSetMotion);
    glutAddMenuEntry("3D Rotation", sm_ROTATION3D);
	glutAddMenuEntry("2D Rotation", sm_ROTATION2D);
	glutAddMenuEntry("In-Depth Translation", sm_TRANSLATION_DEPTH);

	m_MenuDistance = glutCreateMenu(CSimpleGraph::DoSetSize);
    glutAddMenuEntry("Small", sm_SMALL);
    glutAddMenuEntry("Normal", sm_MEDIUM);
    glutAddMenuEntry("Large", sm_LARGE);
    glutAddMenuEntry("Very Large", sm_XLARGE);

	//add more sub-menus here

	m_MenuMain = glutCreateMenu(CSimpleGraph::DoSetMain);
    glutAddSubMenu("Motion", m_MenuMotion);
    glutAddSubMenu("Distance", m_MenuDistance);

    glutAddMenuEntry("Save As...", sm_SAVE);
    glutAddMenuEntry("Exit", sm_EXIT);

	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

//-------------------------------------------------------------------------

void CWinContViewer3D::DoSetSizeContent(int SizeCase) // virtual
{
	double xDim = m_ObjLimits[1] - m_ObjLimits[0];
	double yDim = m_ObjLimits[3] - m_ObjLimits[2];
	double zDim = m_ObjLimits[5] - m_ObjLimits[4];

	double TransSize = yDim;
	if(TransSize < zDim) TransSize = zDim;

	double LocGrid;
	if(SizeCase == sm_SMALL) LocGrid = 2*TransSize/5;
	else if(SizeCase == sm_MEDIUM) LocGrid = 2*TransSize/3;
	else if(SizeCase == sm_LARGE) LocGrid = 2*TransSize/1.5;
	else if(SizeCase == sm_XLARGE) LocGrid = 4*TransSize;
	else LocGrid = TransSize;
    //switch(SizeCase) 
    //{
    //    case sm_SMALL : LocGrid = 2*TransSize/4; break;
    //    case sm_MEDIUM: LocGrid = 2*TransSize/2; break;
    //    case sm_LARGE : LocGrid = 2*TransSize/1.5; break;
    //    case sm_XLARGE : LocGrid = 2*TransSize; break;
    //}

	LocGrid += 0.5*xDim;

    m_DepthNear= LocGrid/20.0;
    m_DepthFar= LocGrid*30.0;

	//if(m_Depth <= 0) m_Depth = (5.0/4.0)*LocGrid;
	////m_DepthMotionQuanta = 0.05*m_Depth;
	//if(m_DepthMotionQuanta <= 0) m_DepthMotionQuanta = 0.2*m_Depth;

	m_Depth = (5.0/4.0)*LocGrid;
	m_DepthMotionQuanta = 0.2*m_Depth;

	SetupTrfCur(m_RotCur);
}

//-------------------------------------------------------------------------

void CWinContViewer3D::DoSetMotionContent(int MotionCase) // virtual
{
	m_MotionType = MotionCase;

	//if(MotionCase == sm_ROTATION3D)
	//{
	//}
	//else if(MotionCase == sm_ROTATION2D)
	//{
	//}
	//else if(MotionCase == sm_TRANSLATION_DEPTH)
	//{
	//}
}

//-------------------------------------------------------------------------

void CWinContViewer3D::SetDisplayMode(int InMode)
{
    m_DisplayMode = InMode;
	if(InMode == sm_WIREFRAME)
	{
        glShadeModel(GL_FLAT); 
        glDisable(GL_LIGHTING);
	}
	else if(InMode == sm_HIDDENLINE)
	{
        glShadeModel(GL_FLAT); 
        glDisable(GL_LIGHTING);
	}
	else if(InMode == sm_FLATSHADED)
	{
        glShadeModel(GL_FLAT); 
        glEnable(GL_LIGHTING);
	}
	else if(InMode == sm_SMOOTHSHADED)
	{
        glShadeModel(GL_SMOOTH); 
        glEnable(GL_LIGHTING);
	}
	else if(InMode == sm_TEXTURED)
	{
        glShadeModel(GL_SMOOTH); 
        glEnable(GL_LIGHTING);
	}
    //switch(value) 
    //{
    //    case WIREFRAME   : 
    //        glShadeModel(GL_FLAT); 
    //        glDisable(GL_LIGHTING);
    //        break;
    //    case HIDDENLINE: 
    //        glShadeModel(GL_FLAT); 
    //        glDisable(GL_LIGHTING);
    //        break;
    //    case FLATSHADED  : 
    //        glShadeModel(GL_FLAT); 
    //        glEnable(GL_LIGHTING);
    //        break;
    //    case SMOOTHSHADED: 
    //        glShadeModel(GL_SMOOTH); 
    //        glEnable(GL_LIGHTING);
    //        break;
    //    case TEXTURED: 
    //        glShadeModel(GL_SMOOTH); 
    //        glEnable(GL_LIGHTING);
    //        break;
    //}
    //glutPostRedisplay();
}

//-------------------------------------------------------------------------

void CWinContViewer3D::SetOtherParams(int InPar)
{
	if(InPar == sm_FULLSCREEN)
	{
        glutFullScreen();
	}
	else if(InPar == sm_FACENORMALS)
	{
        m_DrawFaceNorms = !m_DrawFaceNorms;
	}
	else if(InPar == sm_ANTIALIAS)
	{
		m_Antialias = !m_Antialias;
		if (m_Antialias)
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
	}
	else if(InPar == sm_ENVMAP)
	{
		m_EnvMap = !m_EnvMap;
		if(m_EnvMap)
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
	}

    //switch(InPar)
    //{
    //    case FULLSCREEN: 
    //        glutFullScreen();
    //        break;
    //    case FACENORMALS: 
    //        drawFaceNorms = !drawFaceNorms;
    //        break;
    //    case ANTIALIAS: 
    //        antialias = !antialias;
    //        if (antialias)
    //        {
    //            glEnable(GL_BLEND);
    //            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //            glEnable(GL_LINE_SMOOTH);
    //            glLineWidth(1.5);
    //        }
    //        else
    //        {
    //            glDisable(GL_BLEND);
    //            glDisable(GL_LINE_SMOOTH);
    //            glLineWidth(1.0);
    //        }
    //        break;
    //    case ENVMAP: 
    //        envMap = !envMap;
    //        if (envMap)
    //        {
    //            //glBindTexture(GL_TEXTURE_2D, texId2);
    //            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    //            glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    //            glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    //            glEnable(GL_TEXTURE_GEN_S);
    //            glEnable(GL_TEXTURE_GEN_T);
    //        }
    //        else
    //        {
    //            //glBindTexture(GL_TEXTURE_2D, texId1);
    //            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    //            glDisable(GL_TEXTURE_GEN_S);
    //            glDisable(GL_TEXTURE_GEN_T);
    //        }
    //        break;
    //}
	//glutPostRedisplay();
}

//-------------------------------------------------------------------------
