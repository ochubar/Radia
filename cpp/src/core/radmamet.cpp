/*-------------------------------------------------------------------------
*
* File name:      radmamet.cpp
*
* Project:        RADIA
*
* Description:    Some numerical methods
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radmamet.h"

#include <math.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTMathLinAlgEq::radTMathLinAlgEq(int InSize)
{
	SizeOfMatr = InSize;

	//radTSend Send;
	//try	
	//{
		indxInverseMatrix = new int[SizeOfMatr];
		colInverseMatrix = new double[SizeOfMatr];
		vvLU_Dcmp = new double[SizeOfMatr];
	//}
	//catch (radTException* radExceptionPtr) 
	//{
	//	Send.ErrorMessage(radExceptionPtr->what()); return;
	//}
	//catch (...)
	//{
	//	Send.ErrorMessage("Radia::Error999"); return;
	//}

	// Apply another allocation check procedure here !!!
}

//-------------------------------------------------------------------------

radTMathLinAlgEq::~radTMathLinAlgEq()
{
	delete[] indxInverseMatrix;
	delete[] colInverseMatrix;
	delete[] vvLU_Dcmp;
}

//-------------------------------------------------------------------------

void radTMathLinAlgEq::InverseMatrix(double** DirMatr, int Len, double** InvMatr)
{
	double d;
	double* col = colInverseMatrix;
	int* indx = indxInverseMatrix;
	SizeOfMatr = Len;

	LU_Dcmp(DirMatr, indx, &d);
	int i, j;
	for(j=0; j<Len; j++)
	{
		for(i=0; i<Len; i++) col[i]=0.;
		col[j]=1.;
		LU_BkSb(DirMatr, indx, col);
		for(i=0; i<Len; i++) InvMatr[i][j]=col[i];
	}
}

//-------------------------------------------------------------------------

void radTMathLinAlgEq::LU_Dcmp(double** a, int* indx, double* d)
{// LU decompisition of a square matrix. The source, with minor changes, taken from book
 // "Numerical Recipes in C" by W.H.Press, B.P.Flannery, S.A.Teukolsky, W.T.Vetterling, Cambridge University Press
 // The largest element in the matrix should not be Zero!!!
	const double TINY = 1.E-20;
	int i, imax, j, k;
	double big, dum, sum, temp;
	int n = SizeOfMatr;

	double* vv = vvLU_Dcmp;

	*d = 1.;
	for(i=0; i<n; i++)
	{
		big=0.;
		for(j=0; j<n; j++) if((temp=fabs(a[i][j])) > big) big=temp;
		vv[i]=1./big;
	}
	for(j=0; j<n; j++)
	{
		for(i=0; i<j; i++)
		{
			sum=a[i][j];
			for(k=0; k<i; k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.;
		for(i=j; i<n; i++)
		{
			sum=a[i][j];
			for(k=0; k<j; k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if((dum=vv[i]*fabs(sum)) >= big)
			{
				big=dum; imax=i;
			}
		}
		if(j!=imax)
		{
			for(k=0; k<n; k++)
			{
				dum=a[imax][k]; a[imax][k]=a[j][k]; a[j][k]=dum;
			}
			*d= -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if(a[j][j]==0.) a[j][j]=TINY;
		if(j!=(n-1))
		{
			dum=1./(a[j][j]);
			for(i=j+1; i<n; i++) a[i][j] *= dum;
		}
	}
}

//-------------------------------------------------------------------------

void radTMathLinAlgEq::LU_BkSb(double** a, int* indx, double* b)
{// Should be used, together with LU_Dcmp, for solving systems of linear algebr. equations
	int i, ii=-1, ip, j;
	double sum;
	int n = SizeOfMatr;

	for(i=0; i<n; i++)
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if(ii!=-1) for(j=ii; j<i; j++) sum -= a[i][j]*b[j];
		else if(sum) ii=i;
		b[i]=sum;
	}
	for(i=n-1; i>=0; i--)
	{
		sum=b[i];
		for(j=i+1; j<n; j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

//-------------------------------------------------------------------------
