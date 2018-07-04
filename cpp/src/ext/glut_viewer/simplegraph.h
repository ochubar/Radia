
#ifndef __SIMPLEGRAPH_H
#define __SIMPLEGRAPH_H

#ifndef __OBJHNDL_H
#include "objhndl.h"
#endif

#include <vector>
#include <list>
#include <map>
using namespace std;

#ifdef WIN32
#include <process.h>
#include <Windows.h>
#endif

#ifndef __glut_h__
#include <GL/glut.h>
#endif

//-------------------------------------------------------------------------
// The structure keeping information about content of one window
//-------------------------------------------------------------------------

class CWinCont {
protected:
//public:

	double*	m_PgVertCoord;
	int	m_PgNv;
	int* m_PgVertInd;
	int* m_PgLen;
	float* m_PgColors;
	int	m_PgN;
	//float* m_PgNormCoord; // = 0;
	double* m_PgNormCoord; // = 0; //OC050109
	bool m_PgWereDefined; // = false;

	double*	m_LnVertCoord;
	int	m_LnNv;
	int* m_LnVertInd;
	int* m_LnLen;
	float* m_LnColors;
	int	m_LnN;
	bool m_LnWereDefined; // = false;

	double m_CenVert[3];
	double m_ObjLimits[6];
	double m_MaxSize;

public:

	long m_xSize, m_ySize;
	long m_xSizeOrig, m_ySizeOrig; //used for saving only
	long m_xPos, m_yPos;

	bool m_PgNormalsWereDefined; // = false;

	char m_WinTitle[512];

	static const int sm_EXIT;
	static const int sm_SAVE;

	CWinCont()
	{
		m_xSize = m_ySize = 0;
		m_xSizeOrig = m_ySizeOrig = 0;
		m_xPos = m_yPos = 0;

		m_PgNormCoord = 0;
		m_PgNormalsWereDefined = false;
		m_PgWereDefined = false;
		m_LnWereDefined = false;
		m_WinTitle[0] = '\0';
	}
	virtual ~CWinCont();

	void SetupPgAndLnData(double** VertCoordArr, int* NvArr, int** VertIndArr, int** LenArr, float** ColorsArr, int* NpArr, double* CenVert, double* ObjLim, double MaxSize);
	virtual void SetViewingParamsGL() {};
	virtual void SetupMenuGLUT() {};

	void ProcDisplayContentForView();
	void ProcReshapeContentForView(int Width, int Height);

	virtual void ProcDisplayContent() {};
    virtual void ProcReshapeContent(int Width, int Height) {};
    virtual void ProcMouseContent(int button, int state, int x, int y) {};
    virtual void ProcMotionContent(int x, int y) {};
	virtual void ProcKeyboardContent(unsigned char ch, int x, int y) {};

	virtual void DoSetSizeContent(int SizeCase) {};
    virtual void DoSetMotionContent(int MotionCase) {};
    //virtual void DoSetMainContent(int MainCase) {};
};

//-------------------------------------------------------------------------
// Container of contents of all windows
//-------------------------------------------------------------------------

typedef CHandle<CWinCont> CHWinCont;

typedef map <int, CHWinCont, less<int> > CHWinContMap;
#ifdef __GCC__
typedef list <CHWinCont> CHWinContList;
typedef vector <CHWinCont > CHWinContVect;
#else
typedef list <CHWinCont, allocator<CHWinCont> > CHWinContList;
typedef vector <CHWinCont, allocator<CHWinCont> > CHWinContVect;
#endif

//-------------------------------------------------------------------------
// Auxiliary structures
//-------------------------------------------------------------------------

typedef struct
{
	unsigned char red;
	unsigned char green;
	unsigned char blue;
	unsigned char alpha;
} GLRGBQUAD;

typedef struct
{
	unsigned char red;
	unsigned char green;
	unsigned char blue;
} CRGBTri;

typedef struct tagBITMAPFILEHEADER_HACK
{//derived from wingdi.h
        //WORD    bfType;
        //DWORD   bfSize;
        //WORD    bfReserved1;
        //WORD    bfReserved2;
        //DWORD   bfOffBits;
	unsigned short bfType; //2 bytes
	unsigned long bfSize; //4 bytes
	unsigned short bfReserved1;
	unsigned short bfReserved2;
	unsigned long bfOffBits;
} BITMAPFILEHEADER_HACK; //, FAR *LPBITMAPFILEHEADER, *PBITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER_HACK
{//derived from wingdi.h
        //DWORD      biSize;
        //LONG       biWidth;
        //LONG       biHeight;
        //WORD       biPlanes;
        //WORD       biBitCount;
        //DWORD      biCompression;
        //DWORD      biSizeImage;
        //LONG       biXPelsPerMeter;
        //LONG       biYPelsPerMeter;
        //DWORD      biClrUsed;
        //DWORD      biClrImportant;
	unsigned long biSize; //4 bytes
	long biWidth; //4 bytes
	long biHeight; //4 bytes
	unsigned short biPlanes; //2 bytes
	unsigned short biBitCount; //2 bytes
	unsigned long biCompression; //4 bytes
	unsigned long biSizeImage; //4 bytes
	long biXPelsPerMeter; //4 bytes
	long biYPelsPerMeter; //4 bytes
	unsigned long biClrUsed; //4 bytes
	unsigned long biClrImportant; //4 bytes
} BITMAPINFOHEADER_HACK; //, FAR *LPBITMAPINFOHEADER, *PBITMAPINFOHEADER;

//#ifdef WIN32
//
//#ifndef GLboolean
//typedef unsigned char GLboolean;
//#endif
//#ifndef GL_TRUE
//#define GL_TRUE 1
//#endif
//#ifndef GL_FALSE
//#define GL_FALSE 0
//#endif

//#ifndef M__pixelformat__
//struct M__pixelformat__
//{
//    PIXELFORMATDESCRIPTOR	pfd;
//    GLboolean doubleBuffered;
//};
//#endif
//#endif
////int qt_pix = sizeof(pix) / sizeof(pix[0]);


//-------------------------------------------------------------------------
// The base class of the OpenGL-GLUT 3D Viewer and Grapher
//-------------------------------------------------------------------------

class CSimpleGraph {

	//static GLRGBQUAD* sm_pDataToSaveToFile;
	static char sm_GraphFileName[1024];

public:

	static bool sm_WasStarted;
	static bool sm_StartingInProcess;
    //static bool sm_WinShouldBeCreated;
	//static vector <int> sm_WinIDsVect;
	static CHWinContMap sm_WinContMap;
	static CHWinContList sm_WinToCreateContList;
	static bool sm_SomethingToCreate; //, sm_NoMoreWindows;
	static CHWinCont sm_hFirstWinCont;

	static int sm_DefWinPosX, sm_DefWinPosY;
    static int sm_DefWinSizeX, sm_DefWinSizeY;
    static char sm_DefWinTitle[];

//public:

	CSimpleGraph() {}
	~CSimpleGraph() {}

	static void Init();
	static void StartNewWinApp(CHWinCont& hWinCont, char StartMode);
	//static void StartNewWin(CHWinCont& hWinCont, int xPos, int yPos, int xSize, int ySize, char* strTitle);
	static void StartNewWin(CHWinCont& hWinCont, int xPos, int yPos, int xSize, int ySize);
	static void StartNewWin(CHWinCont& hWinCont);
	static void StartNewWindowIfNecessary();

	static void DestroyWindow(int WinID);

	static bool ExtractDispContentData(int WinID, CHWinCont& hCont);

	static void WrapForMT_glutMainLoop(void*); //called at creation of new thread
	//static void RunNewGraphCreator(void*); //called at creation of new thread

	static bool GetCurHWinCont(CHWinCont& hCurCont);

	//static void ProcIdle()
	//{
	//	if(sm_SomethingToCreate)
	//	{
	//		StartNewWindowIfNecessary();
	//	}
	//#ifdef WIN32
	//	//Sleep(5); // To release CPU resources; 100 ms will reduce viewing performance
    //  Sleep(1); // To release CPU resources; 100 ms will reduce viewing performance
	//#endif
	//}
	static void ProcTimer(int Value)
	{
		if(sm_SomethingToCreate)
		{
			StartNewWindowIfNecessary();
		}
		glutTimerFunc(Value, ProcTimer, Value);
	}

	static void ProcDisplay();
	static void ProcReshape(int width, int height);
	static void ProcKeyboard(unsigned char ch, int x, int y);
    static void ProcMouse(int button, int state, int x, int y);
	static void ProcMotion(int x, int y);
	//static void ProcPassiveMotion(int x, int y);
    //static void ProcWinStatus(int StatusID);

	static void DoSetSize(int SizeCase);
	static void DoSetMotion(int MotionCase);
	static void DoSetMain(int MainCase);
	
	static void	SaveWinContentToGraphFile(int WinID);
	static void FinishSavingWinContentToGraphFile(CHWinCont& hCurCont);

	static void	GetFileName(char* sFileName);
	static bool	SaveImageToGraphFile(GLRGBQUAD *pData, int Width, int Height, bool SwapRB, const char* sFileName);
    static bool SaveImageToPNG(GLRGBQUAD *pData, int Width, int Height, const char *sFileName);
	static bool SaveImageToPNG(CRGBTri *pData, int Width, int Height, const char *filename);

    static bool SaveImageToBMP(GLRGBQUAD *pData, int Width, int Height, const char *sFileName);
    //static bool SaveImageToJPG(GLRGBQUAD *pData, int Width, int Height, const char *sFileName);
	static void AuxSwapRB(GLRGBQUAD *pixel, int Width, int Height);
	static void RotateBitmapVertically(GLRGBQUAD *pData, int Width, int Height);
    static void CorrectAlphaInBitmapRGBA(GLRGBQUAD *pData, int Width, int Height);
    static void CopyRGBA2RBG(GLRGBQUAD *pRGBA, CRGBTri *pRGB, int Width, int Height);
	static void BGRA2RGBA(GLRGBQUAD *pData, int Width, int Height);
	static void AuxMakeWhite(GLRGBQUAD *pData, long Size);

	static void RenderImageToSave(int WinID, int SaveWidth, int SaveHeight);

//#ifdef WIN32
//	static int M_wglChoosePixelFormat(HDC hdc, CONST PIXELFORMATDESCRIPTOR *ppfd);
//#endif
};

//-------------------------------------------------------------------------

#endif
