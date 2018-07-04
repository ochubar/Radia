
#ifndef __GMINTERP_H
#define __GMINTERP_H

#ifndef _GM_WITHOUT_BASE
#ifndef __GMOBJ_H
#include "gmobj.h"
#endif
#endif

//-------------------------------------------------------------------------

class CGenMathInterp
#ifndef _GM_WITHOUT_BASE
	: public CGenMathObj
#endif
{
	double *mSplineY2Arr, *mSplineArgTabArr, *mSplineValTabArr;
	double mArgStep, mArgStart;
	int mSplineTabNp;
	int mMethNo;

public:

	CGenMathInterp(int MethNo, double *x, double *y, int np);
	CGenMathInterp(int MethNo, double xStart, double xStep, double *y, int np);

	CGenMathInterp()
	{
		mMethNo = 0;
        mSplineY2Arr = 0; mSplineArgTabArr = 0; mSplineValTabArr = 0;
		mArgStep = 0; mArgStart = 0;
		mSplineTabNp = 0;
	}
	~CGenMathInterp()
	{
		if(mSplineY2Arr != 0) { delete[] mSplineY2Arr; mSplineY2Arr = 0;}
		if(mSplineArgTabArr != 0) { delete[] mSplineArgTabArr; mSplineArgTabArr = 0;}
		if(mSplineValTabArr != 0) { delete[] mSplineValTabArr; mSplineValTabArr = 0;}
	}

	//static void CubSplinePrep(double *x, double *y, int n, double yp1, double ypn, double *y2);
	static void InterpCubicSplinePrep(double *x, double *y, int n, double *y2);
	static void InterpCubicSplinePrepU(double xStart, double xStep, double *y, int n, double *y2);
	static double InterpCubicSpline(double *xa, double *ya, double *y2a, int n, double x);
	static double Deriv1(double* f, double h, int PoIndx, int AmOfPo);

	void InitCubicSpline(double *x, double *y, int np);
	void InitCubicSplineU(double xStart, double xStep, double *y, int np);
	double Interp1D(double x)
	{
		if((mMethNo == 1) && (mSplineY2Arr != 0) && (mSplineArgTabArr != 0) && (mSplineValTabArr != 0) && (mSplineTabNp > 0))
		{
			return InterpCubicSpline(mSplineArgTabArr, mSplineValTabArr, mSplineY2Arr, mSplineTabNp, x);
		}
		return 0;
	}
	
	double InterpRelCubicSpline(double b, int i0, double curArgStep)
	{
		if((mSplineY2Arr == 0) || (mSplineValTabArr == 0)) return 0;

		//b = (x - xa[klo])/h; 
		double a = 1. - b; //(xa[khi] - x)/h; //Cubic spline polynomial is now evaluated.
		double ya_lo = mSplineValTabArr[i0], ya_hi = mSplineValTabArr[i0 + 1];
		double y2a_lo = mSplineY2Arr[i0], y2a_hi = mSplineY2Arr[i0 + 1];
		return a*ya_lo + b*ya_hi + ((a*a*a - a)*y2a_lo + (b*b*b - b)*y2a_hi)*(curArgStep*curArgStep)/6.0;
	}
	double InterpRelCubicSplineU(double b, int i0)
	{
		if((mSplineY2Arr == 0) || (mSplineValTabArr == 0) || (mArgStep == 0)) return 0;

		//b = (x - xa[klo])/h; 
		double a = 1. - b; //(xa[khi] - x)/h; //Cubic spline polynomial is now evaluated.
		double ya_lo = mSplineValTabArr[i0], ya_hi = mSplineValTabArr[i0 + 1];
		double y2a_lo = mSplineY2Arr[i0], y2a_hi = mSplineY2Arr[i0 + 1];
		return a*ya_lo + b*ya_hi + ((a*a*a - a)*y2a_lo + (b*b*b - b)*y2a_hi)*(mArgStep*mArgStep)/6.0;
	}

	static double Interp3dBilin(double* inP, double* inArrArgBounds, double* inArrFunc)
	{
	// The approach taken from Numerical Recipes (p. 104), extended to 3d
	// assumes:
	// inP[] = {x,y,z}; //point were the function should be computed
	// inArrArgBounds[] = {x0,x1,y0,y1,z0,z1}; //argument values at the corners of the cube
	// inArrFunc[] = {f(x0,y0,z0),f(x1,y0,z0),f(x0,y1,z0),f(x0,y0,z1),f(x1,y1,z0),f(x1,y0,z1),f(x0,y1,z1),f(x1,y1,z1)} //function values at the corners of the cube
		double xt = 0, yt = 0, zt = 0;
		double &x0 = inArrArgBounds[0], &x1 = inArrArgBounds[1], &y0 = inArrArgBounds[2], &y1 = inArrArgBounds[3], &z0 = inArrArgBounds[4], &z1 = inArrArgBounds[5];
		if(x1 != x0) xt = (inP[0] - x0)/(x1 - x0);
		if(y1 != y0) yt = (inP[1] - y0)/(y1 - y0);
		if(z1 != z0) zt = (inP[2] - z0)/(z1 - z0);
		double one_mi_xt = 1 - xt, one_mi_yt = 1 - yt, one_mi_zt = 1 - zt;
		return inArrFunc[0]*one_mi_xt*one_mi_yt*one_mi_zt
			+ inArrFunc[1]*xt*one_mi_yt*one_mi_zt
			+ inArrFunc[2]*one_mi_xt*yt*one_mi_zt
			+ inArrFunc[3]*one_mi_xt*one_mi_yt*zt
			+ inArrFunc[4]*xt*yt*one_mi_zt
			+ inArrFunc[5]*xt*one_mi_yt*zt
			+ inArrFunc[6]*one_mi_xt*yt*zt
			+ inArrFunc[7]*xt*yt*zt;
	}

	static double Interp3dBilinRel(double xt, double yt, double zt, double* inArFunc)
	{// The approach taken from Numerical Recipes (p. 104), extended to 3d
	 // assumes:
	 // inArFunc[] = {f(x0,y0,z0),f(x1,y0,z0),f(x0,y1,z0),f(x0,y0,z1),f(x1,y1,z0),f(x1,y0,z1),f(x0,y1,z1),f(x1,y1,z1)} //function values at the corners of the cube
		double one_mi_xt = 1.- xt, one_mi_yt = 1.- yt, one_mi_zt = 1.- zt;
		return inArFunc[0]*one_mi_xt*one_mi_yt*one_mi_zt
			+ inArFunc[1]*xt*one_mi_yt*one_mi_zt
			+ inArFunc[2]*one_mi_xt*yt*one_mi_zt
			+ inArFunc[3]*one_mi_xt*one_mi_yt*zt
			+ inArFunc[4]*xt*yt*one_mi_zt
			+ inArFunc[5]*xt*one_mi_yt*zt
			+ inArFunc[6]*one_mi_xt*yt*zt
			+ inArFunc[7]*xt*yt*zt;
	}

	static double Interp3dQuadRel(double xt, double yt, double zt, double* arF)
	{//assumes:
	 //double arF[] = {f(x0,y0,zm1),f(x0,ym1,z0),f(xm1,y0,z0),f(x0,y0,z0),f(x1,y0,z0),f(x0,y1,z0),f(x1,y1,z0),f(x0,y0,z1),f(x1,y0,z1),f(x0,y1,z1)};
		double ax2 = 0.5*(-2.*arF[3] + arF[4] + arF[2]); //1/2 (-2 f000 + f100 + fm100)
		double ay2 = 0.5*(-2.*arF[3] + arF[5] + arF[1]); //1/2 (-2 f000 + f010 + f0m10)
		double az2 = 0.5*(-2.*arF[3] + arF[7] + arF[0]); //1/2 (-2 f000 + f001 + f00m1)
		double axy = arF[3] - arF[5] - arF[4] + arF[6]; //f000 - f010 - f100 + f110
		double axz = arF[3] - arF[7] - arF[4] + arF[8]; //f000 - f001 - f100 + f101
		double ayz = arF[3] - arF[7] - arF[5] + arF[9]; //f000 - f001 - f010 + f011
		double ax = 0.5*(arF[4] - arF[2]); //(f100 - fm100)/2
		double ay = 0.5*(arF[5] - arF[1]); //(f010 - f0m10)/2
		double az = 0.5*(arF[7] - arF[0]); //(f001 - f00m1)/2
		return arF[3] + xt*(ax + xt*ax2 + yt*axy + zt*axz) + yt*(ay + yt*ay2 + zt*ayz) + zt*(az + zt*az2);
	}

	static double Interp3dCubicRel(double xt, double yt, double zt, double* arF)
	{//assumes:
	 //double arF[] = {f00m1,f0m10,fm100,f000,f100,f200,f010,f110,f210,f020,f120,f001,f101,f201,f011,f111,f021,f002,f102,f012};
		const double i6 = 1./6.;
		double ax3 = i6*(3.*(arF[3] - arF[4]) + arF[5] - arF[2]); //1/6 (3 f000 - 3 f100 + f200 - fm100)
		double ay3 = i6*(3.*(arF[3] - arF[6]) + arF[9] - arF[1]); //1/6 (3 f000 - 3 f010 + f020 - f0m10)
		double az3 = i6*(3.*(arF[3] - arF[11]) + arF[17] - arF[0]); //1/6 (3 f000 - 3 f001 + f002 - f00m1)
		double ax2y = 0.5*(-arF[3] + arF[6] + 2.*(arF[4] - arF[7]) - arF[5] + arF[8]); //1/2 (-f000 + f010 + 2 f100 - 2 f110 - f200 + f210)
		double ax2z = 0.5*(-arF[3] + arF[11] + 2.*(arF[4] - arF[12]) - arF[5] + arF[13]); //1/2 (-f000 + f001 + 2 f100 - 2 f101 - f200 + f201)
		double axy2 = 0.5*(-arF[3] + 2.*(arF[6] - arF[7]) - arF[9] + arF[4] + arF[10]); //1/2 (-f000 + 2 f010 - f020 + f100 - 2 f110 + f120)
		double axz2 = 0.5*(-arF[3] + 2.*arF[11] - arF[17] + arF[4] - 2.*arF[12] + arF[18]); //1/2 (-f000 + 2 f001 - f002 + f100 - 2 f101 + f102)
		double ay2z = 0.5*(-arF[3] + arF[11] + 2.*arF[6] - 2.*arF[14] - arF[9] + arF[16]); //1/2 (-f000 + f001 + 2 f010 - 2 f011 - f020 + f021)
		double ayz2 = 0.5*(-arF[3] + 2.*arF[11] - arF[17] + arF[6] - 2.*arF[14] + arF[19]); //1/2 (-f000 + 2 f001 - f002 + f010 - 2 f011 + f012)
		double axyz = -arF[3] + arF[11] + arF[6] - arF[14] + arF[4] - arF[12] - arF[7] + arF[15]; //-f000 + f001 + f010 - f011 + f100 - f101 - f110 + f111
		double ax2 = 0.5*(-2.*arF[3] + arF[4] + arF[2]); //1/2 (-2 f000 + f100 + fm100)
		double ay2 = 0.5*(-2.*arF[3] + arF[6] + arF[1]); //1/2 (-2 f000 + f010 + f0m10)
		double az2 = 0.5*(-2.*arF[3] + arF[11] + arF[0]); //1/2 (-2 f000 + f001 + f00m1)
		double axy = 0.5*(4.*arF[3] - 5.*arF[6] + arF[9] - 5.*arF[4] + 6.*arF[7] - arF[10] + arF[5] - arF[8]); //1/2 (4 f000 - 5 f010 + f020 - 5 f100 + 6 f110 - f120 + f200 - f210)
		double axz = 0.5*(4.*arF[3] - 5.*arF[11] + arF[17] - 5.*arF[4] + 6.*arF[12] - arF[18] + arF[5] - arF[13]); //1/2 (4 f000 - 5 f001 + f002 - 5 f100 + 6 f101 - f102 + f200 - f201)
		double ayz = 0.5*(4.*arF[3] - 5.*arF[11] + arF[17] - 5.*arF[6] + 6.*arF[14] - arF[19] + arF[9] - arF[16]); //1/2 (4 f000 - 5 f001 + f002 - 5 f010 + 6 f011 - f012 + f020 - f021)
		double ax = i6*(-3.*arF[3] + 6.*arF[4] - arF[5] - 2.*arF[2]); //1/6 (-3 f000 + 6 f100 - f200 - 2 fm100)
		double ay = i6*(-3.*arF[3] + 6.*arF[6] - arF[9] - 2.*arF[1]); //1/6 (-3 f000 + 6 f010 - f020 - 2 f0m10)
		double az = i6*(-3.*arF[3] + 6.*arF[11] - arF[17] - 2.*arF[0]); //1/6 (-3 f000 + 6 f001 - f002 - 2 f00m1)
		return arF[3] + xt*(ax + xt*(ax2 + xt*ax3 + yt*ax2y + zt*ax2z) + yt*(axy + yt*axy2 + zt*axyz) + zt*(axz + zt*axz2)) +
						yt*(ay + yt*(ay2 + yt*ay3 + xt*axy2 + zt*ay2z) + zt*(ayz + zt*ayz2)) + zt*(az + zt*(az2 + zt*az3));
	}

	static double Interp3dBiCubic32pRel(double xt, double yt, double zt, double* arF)
	{
		double *p = arF;
		double f00m1,f10m1,f01m1,f11m1,f0m10,f1m10,fm100,f000,f100,f200,fm110,f010,f110,f210,f020,f120,f0m11,f1m11,fm101,f001,f101,f201,fm111,f011,f111,f211,f021,f121,f002,f102,f012,f112;
		f00m1=*(p++); f10m1=*(p++); f01m1=*(p++); f11m1=*(p++);
		f0m10=*(p++); f1m10=*(p++); fm100=*(p++); f000=*(p++); f100=*(p++); f200=*(p++); fm110=*(p++); f010=*(p++); f110=*(p++); f210=*(p++); f020=*(p++); f120=*(p++); 
		f0m11=*(p++); f1m11=*(p++); fm101=*(p++); f001=*(p++); f101=*(p++); f201=*(p++); fm111=*(p++); f011=*(p++); f111=*(p++); f211=*(p++); f021=*(p++); f121=*(p++); 
		f002=*(p++); f102=*(p++); f012=*(p++); f112=*(p++);

		const double i6 = 1./6.;
		double b5 = 3.*(f000 - f001 - f010 + f011 - f100 + f101 + f110 - f111);
		double a311 = i6*(b5 + f200 - f201 - f210 + f211 - fm100 + fm101 + fm110 - fm111);
		double a131 = i6*(b5 + f020 - f021 - f0m10 + f0m11 - f120 + f121 + f1m10 - f1m11);
		double a113 = i6*(b5 + f002 - f00m1 - f102 - f012 + f01m1 + f10m1 + f112 - f11m1);
		double b4xy = 3.*(f010 - f000 + f100 - f110);
		double a310 = i6*(b4xy - f200 + f210 + fm100 - fm110);
		double a130 = i6*(b4xy - f020 + f0m10 + f120 - f1m10);
		double b4xz = 3.*(f001 - f000 + f100 - f101);
		double a301 = i6*(b4xz - f200 + f201 + fm100 - fm101);
		double a103 = i6*(b4xz - f002 + f00m1 + f102 - f10m1);
		double b4yz = 3.*(f001 - f000 + f010 - f011);
		double a031 = i6*(b4yz - f020 + f021 + f0m10 - f0m11);
		double a013 = i6*(b4yz - f002 + f00m1 + f012 - f01m1);
		double a211 = 0.5*(2.*(f001 - f000 + f010 - f011) + f100 - f101 - f110 + f111 + fm100 - fm101 - fm110 + fm111);
		double a121 = 0.5*(2.*(f001 - f000 + f100 - f101) + f010 - f011 + f0m10 - f0m11 - f110 + f111 - f1m10 + f1m11);
		double a112 = 0.5*(2.*(f010 - f000 + f100 - f110) + f001 + f00m1 - f011 - f01m1 - f101 - f10m1 + f111 + f11m1);
		double a300 = i6*(3.*(f000 - f100) + f200 - fm100);
		double a030 = i6*(3.*(f000 - f010) + f020 - f0m10);
		double a003 = i6*(3.*(f000 - f001) + f002 - f00m1);
		double a210 = 0.5*(2.*(f000 - f010) - f100 + f110 - fm100 + fm110); 
		double a201 = 0.5*(2.*(f000 - f001) - f100 + f101 - fm100 + fm101);
		double a120 = 0.5*(2.*(f000 - f100) - f010 - f0m10 + f110 + f1m10);
		double a021 = 0.5*(2.*(f000 - f001) - f010 + f011 - f0m10 + f0m11);
		double a102 = 0.5*(2.*(f000 - f100) - f001 - f00m1 + f101 + f10m1);
		double a012 = 0.5*(2.*(f000 - f010) - f001 - f00m1 + f011 + f01m1);
		double a111 = f111 + i6*(3.*(f000 - f011 - f101 - f110) + 2.*(f01m1 - f00m1 - f0m10 + f0m11 + f10m1 - f11m1 + f1m10 - f1m11 - fm100 + fm101 + fm110 - fm111) - f002 + f012 - f020 + f021 + f102 - f112 + f120 - f121 - f200 + f201 + f210 - f211);
		double a200 = 0.5*(f100 + fm100) - f000;
		double a020 = 0.5*(f010 + f0m10) - f000;
		double a002 = 0.5*(f001 + f00m1) - f000;
		double a110 = f110 + i6*(-3.*(f010 + f100) + 2.*(f0m10 - f1m10 + fm100 - fm110) + f020 - f120 + f200 - f210);
		double a101 = f101 + i6*(-3.*(f001 + f100) + 2.*(f00m1 - f10m1 + fm100 - fm101) + f002 - f102 + f200 - f201);
		double a011 = f011 + i6*(-3.*(f001 + f010) + 2.*(f00m1 - f01m1 + f0m10 - f0m11) + f002 - f012 + f020 - f021);
		double a100 = f100 + i6*(-3.*f000 - 2.*fm100 - f200);
		double a010 = f010 + i6*(-3.*f000 - 2.*f0m10 - f020);
		double a001 = f001 + i6*(-3.*f000 - 2.*f00m1 - f002);

		return f000 + xt*(a100 + xt*(a200 + xt*(a300 + yt*(a310 + a311*zt) + a301*zt) + yt*(a210 + a211*zt) + a201*zt) 
							   + yt*(a110 + yt*(a120 + yt*(a130 + a131*zt) + a121*zt) + zt*(a111 + zt*(a112 + a113*zt))) 
							   + zt*(a101 + zt*(a102 + a103*zt))) 
					+ yt*(a010 + yt*(a020 + yt*(a030 + a031*zt) + a021*zt) + zt*(a011 + zt*(a012 + a013*zt)))
					+ zt*(a001 + zt*(a002 + a003*zt));
	}

	static double Interp2dBiCubic12pRel(double xt, double yt, double* arF)
	{
		double *p = arF;
		double f0m1,f1m1,fm10,f00,f10,f20,fm11,f01,f11,f21,f02,f12;
		f0m1=*(p++); f1m1=*(p++); fm10=*(p++); f00=*(p++); f10=*(p++); f20=*(p++); fm11=*(p++); f01=*(p++); f11=*(p++); f21=*(p++); f02=*(p++); f12=*(p++); 

		const double i6 = 1./6.;
		double b4xy = 3.*(f01 - f00 + f10 - f11);
		double a31 = i6*(b4xy - f20 + f21 + fm10 - fm11);
		double a13 = i6*(b4xy - f02 + f0m1 + f12 - f1m1);
		double a30 = i6*(3.*(f00 - f10) + f20 - fm10);
		double a03 = i6*(3.*(f00 - f01) + f02 - f0m1);
		double a21 = 0.5*(2.*(f00 - f01) - f10 + f11 - fm10 + fm11); 
		double a12 = 0.5*(2.*(f00 - f10) - f01 - f0m1 + f11 + f1m1);
		double a20 = 0.5*(f10 + fm10) - f00;
		double a02 = 0.5*(f01 + f0m1) - f00;
		double a11 = f11 + i6*(-3.*(f01 + f10) + 2.*(f0m1 - f1m1 + fm10 - fm11) + f02 - f12 + f20 - f21);
		double a10 = f10 + i6*(-3.*f00 - 2.*fm10 - f20);
		double a01 = f01 + i6*(-3.*f00 - 2.*f0m1 - f02);

		return f00 + xt*(a10 + xt*(a20 + xt*(a30 + yt*a31) + yt*a21) + yt*(a11 + yt*(a12 + yt*a13))) + yt*(a01 + yt*(a02 + yt*a03));
	}

	static double InterpCubHalfStep(double* f, int i)
	{//interpolation by cubic poligon for a point in the middle between equidistant points for which the function f is defined 
		if(i < 0) return 0.0625*(5.*f[0] + 15.*f[1] - 5.*f[2] + f[3]);
		else if(i == 0) return 0.0625*(-f[0] + 9.*f[1] + 9.*f[2] - f[3]);
		else return 0.0625*(f[0] - 5.*f[1] + 15.*f[2] + 5.*f[3]);
	}
};

//-------------------------------------------------------------------------

#endif
