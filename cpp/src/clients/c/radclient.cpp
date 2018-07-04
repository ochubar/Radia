/////////////////////////////////////////////////////////////////////////////
// Name:        radclient.cpp
// Purpose:     Example C++ file illustrating use of Radia library
// Author:      O.C.
// Version:     4.115
// Modified:    30.05.04
/////////////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include <iostream>
#include <time.h>
using namespace std;

#include "radentry.h"

//----------------------------------------------------------------------
//Auxiliary function dedicated to process errors reported by DLL
//----------------------------------------------------------------------

void ErrProc(int er)
{
	char ErrorBuf[2048];

	if(er == 0) return;
	else if(er < 0)
	{
		RadWarGetText(ErrorBuf, er);
		cout << endl << "WARNING: " << ErrorBuf << endl;
		//return/send warning message in place of the above
	}
	else if(er > 0)
	{
		RadErrGetText(ErrorBuf, er);
		cout << endl << "ERROR: " << ErrorBuf << endl << endl;
		cout << "Press \"Enter\" to exit" << endl;
		getchar();
		exit(0);
		//return/send error message (e.g. display an error dialog in GUI) in place of the above
	}
}

//----------------------------------------------------------------------
//The main function iilustrates calls to radia.dll
//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
	cout << "///////////////////////////////////////////////////////////////" << endl; 
	cout << "Creating and solving a 3D geometry having rotational symmetry:" << endl; 
	cout << "///////////////////////////////////////////////////////////////" << endl; 
	cout << endl; 

	int IronPart1 = 0, IronPart2 = 0, IronPart3 = 0, IronGeom = 0; 
	int NumSegm = 15; //25; //Number of azimuthal segments
	double Pc[] = {0, 0}, Angles[] = {0, 2*3.14159}, M0[] = {0, 0, 0};
	
	cout << "Iron part 1 (exit 3D viewer to continue): "; 

	//First part of the iron geometry (NOTE: the cross-section should be a convex polygon !)
    //Flat 2D profile array {r1, z1, r2, z2, ...}, coordinates are in [mm]
	double FlatProfilePart1[] = {23, -50, 40, -50, 60, -36, 40, -24, 23, -36}; 

	int res = RadObjArcPgnMag(&IronPart1, Pc, 'z', FlatProfilePart1, 5, Angles, NumSegm, 's', M0);

	ErrProc(res);
	RadObjDrwOpenGL(IronPart1, 0); //Viewing the geometry in 3D (OpenGL) viewer

	cout << "done" << endl; 
	cout << "Iron part 2 (exit 3D viewer to continue): "; 

	//Second part of the iron geometry (NOTE: the cross-section should be a convex polygon !)
    double FlatProfilePart2[] = {40, -24, 60, -36, 60, 50, 40, 25}; //Flat 2D profile array {r1, z1, r2, z2, ...}
	ErrProc(RadObjArcPgnMag(&IronPart2, Pc, 'z', FlatProfilePart2, 4, Angles, NumSegm, 's', M0));
	RadObjDrwOpenGL(IronPart2, 0); //Viewing the geometry in 3D (OpenGL) viewer

	cout << "done" << endl; 
	cout << "Iron part 3 (exit 3D viewer to continue): "; 

	//Third part of the iron geometry (NOTE: the cross-section should be a convex polygon !)
    double FlatProfilePart3[] = {40, 25, 60, 50, 17, 50, 13, 46, 13, 33, 30, 25}; //Flat 2D profile array {r1, z1, r2, z2, ...}
	ErrProc(RadObjArcPgnMag(&IronPart3, Pc, 'z', FlatProfilePart3, 6, Angles, NumSegm, 's', M0));
	RadObjDrwOpenGL(IronPart3, 0); //Viewing the geometry in 3D (OpenGL) viewer

	cout << "done" << endl; 
	cout << "Iron geometry before segmentation: "; 

	//All iron geometry (container)
	int IronGeomElems[] = {IronPart1, IronPart2, IronPart3};
	ErrProc(RadObjCnt(&IronGeom, IronGeomElems, 3)); //Creating a container

	RadObjDpl(&IronGeom, IronGeom, "FreeSym->True"); //Releasing degrees of freedom for relaxation

	RadObjDrwOpenGL(IronGeom, 0); //Viewing the geometry in 3D (OpenGL) viewer

	cout << "done" << endl; 

	//Creating a standard (iron) material
	int IronMater = 0;
	//ErrProc(RadMatStd(&IronMater, "Xc06", 0));

	double MatData[] = {0.0000123, 0.0149877, 0.0001, 0.7999, 0.00015, 1.42, 0.0025, 2.098, 0.05, 2.29};
    ErrProc(RadMatSatIsoTab(&IronMater, MatData, 5));

	//Applying the magnetic material to the geometry
    ErrProc(RadMatApl(&IronGeom, IronGeom, IronMater));

	cout << "Iron geometry after complete segmentation (exit 3D viewer to continue): "; 

	//Segmeting the iron geometry over R and Z
	//int nr = 15, nz = 15; //Number of pieces in radial and axial directions
	int nr = 5, nz = 1; //Number of pieces in radial and axial directions

	double rmax = 60; //Maximal radial coordinate
	double SbdPar[] = {nr, 1, nz};
	double FlatCylPar[] = {Pc[0], Pc[1], 0, 0, 0, 1, rmax, 0, 0}; //Definition of the segmenting cylinders (see Radia DLL help for details)
    ErrProc(RadObjDivMagCyl(&IronGeom, IronGeom, SbdPar, 3, FlatCylPar, 1, "Frame->LabTot"));
	RadObjDrwOpenGL(IronGeom, 0); //Viewing the geometry in 3D (OpenGL) viewer

	cout << "done" << endl; 
	cout << "Coil (exit 3D viewer to continue): "; 

	//Creating a circular coil with rectangular cross-section
	int Coil = 0;
	double PcCoil[] = {Pc[0], Pc[1], 0}; //Center point
	double Radii[] = {30, 40};
	double height = 50;
	double j = 100; //Current density in [A/mm^2]
	double ColorRGB[] = {1., 0, 0};
    ErrProc(RadObjArcCur(&Coil, PcCoil, Radii, Angles, height, NumSegm, 'm', 'z', j));
	ErrProc(RadObjDrwAtr(Coil, ColorRGB, 0.001)); //Applying color to the coil
	RadObjDrwOpenGL(Coil, 0); //Viewing the geometry in 3D (OpenGL) viewer

	cout << "done" << endl; 
	cout << "Iron and coil together (exit 3D viewer to continue): "; 

	//Complete geometry: iron and coil together in a container
	int TotGeom = 0;
	int TotGeomElems[] = {IronGeom, Coil};
	ErrProc(RadObjCnt(&TotGeom, TotGeomElems, 2)); //Creating a container
	RadObjDrwOpenGL(TotGeom, 0); //Viewing the geometry in 3D (OpenGL) viewer

	cout << "done" << endl; 
	cout << "Started solving (relaxing) the geometry (please wait...): "; 

	clock_t t0 = clock();

	//Relaxing the geometry
	double Inf[10];
	int LenInf;
	double prec = 0.0001; //Accuracy level for magnetization in [T]
	int MaxIter = 2000; //Maximum number of iterations to be made
	ErrProc(RadSolve(Inf, &LenInf, TotGeom, prec, MaxIter, 0));

	int dt_s = (int)(0.001*(clock() - t0));

	cout << "Solved in " << dt_s << " s" << endl; 
	cout << "Stability of magnetization: " << Inf[0] << " T" << endl; 
	cout << "Number of iterations made: " << Inf[3] << endl; 
	cout << "Calculating axial magnetic field (please wait...): "; 

	t0 = clock();

	//Calculating magnetic field along straight line
	double Z_Bz[200];
	int NumB, NumP = 53;
	double zmin = -130, zmax = 130; //Min. and max axial positions
	//double StartP[] = {0, 0, zmin}, FinP[] = {0.1, 0.1, zmax}; //Initial and final observation points
	double StartP[] = {0.1, 0.1, zmin}, FinP[] = {0.1, 0.1, zmax}; //Initial and final observation points
	ErrProc(RadFldLst(Z_Bz, &NumB, TotGeom, "bz", StartP, FinP, NumP, "arg", zmin));

	dt_s = (int)(0.001*(clock() - t0));
	cout << "done" << endl << "Field in " << NumP << " points calculated in " << dt_s << " s" << endl << endl; 

	for(int i=0; i<NumP; i++)
	{
		cout << "z = " << Z_Bz[2*i] << " mm:   Bz = " << Z_Bz[2*i + 1] << " T" << endl; 
	}
	cout << endl; 

	cout << "///////////////////////////////////////////////////////////////" << endl; 
	cout << "Examples of other 3D geometries:" << endl; 
	cout << "///////////////////////////////////////////////////////////////" << endl; 
	cout << endl; 
	cout << "Press \"Enter\" to continue" << endl;
	getchar();

	cout << "- Rectangular Parallelepiped: (exit 3D viewer to continue) "; 

	double P[] = {0., 0., -3.}, L[] = {3., 3., 3.}, M[] = {0., 0., 1.};
	int mag = 0;
	ErrProc(RadObjRecMag(&mag, P, L, M)); //Creating a uniformly magnetised Rectangular Parallelepiped
	double MagColorRGB[] = {0, 0, 0.8};
	ErrProc(RadObjDrwAtr(mag, MagColorRGB, 0.001)); //Applying color to the coil
	RadObjDrwOpenGL(mag, 0); //Viewing the geometry in 3D (OpenGL) viewer

	cout << "done" << endl; 
	cout << "- Polyhedrons: (exit 3D viewer to continue) "; 

	double FV[] = {0.,0.,1.,1.,-1.,1.,-1.,-1.,1.,-1.,0.,0.,0.,1.,1.,1.,1.,0.};
	int SL[] = {1,4,1,3};
	double AT[] = {0., 1., 2., 3.};
	int ppoly;
	ErrProc(RadObjMltExtPgn(&ppoly,FV,SL,AT,4,M)); //Creating uniformly magnetised Polyhedrons
	ErrProc(RadObjDrwAtr(ppoly, MagColorRGB, 0.001)); //Applying color to the coil
	RadObjDrwOpenGL(ppoly, 0);

	cout << "done" << endl; 
	cout << "Press \"Enter\" to exit" << endl;

	getchar();
	return 0;
}

//----------------------------------------------------------------------
