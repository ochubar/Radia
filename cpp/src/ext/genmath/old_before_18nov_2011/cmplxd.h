
//-------------------------------------------------------------------------
//	Definition of simplest complex<double> structure
//-------------------------------------------------------------------------

#ifndef __CMPLXD_H
#define __CMPLXD_H

//#include <math.h>

//-------------------------------------------------------------------------

struct TComplexD {
	double x, y;

	TComplexD(double xx =0, double yy =0) { x=xx; y=yy;}
	TComplexD(double* dArray) { x=dArray[0]; y=dArray[1];}

	double AbsE2() { return (x*x + y*y);}
	double Abs() { return sqrt(x*x + y*y);}
	void Conjugate() { y = -y;}

	TComplexD& operator +=(const TComplexD& c)
	{
		x += c.x; y += c.y; return *this;
	}
	TComplexD& operator *=(const TComplexD& c)
	{
		double ax = x*c.x - y*c.y;
		double ay = x*c.y + y*c.x;
		x = ax; y = ay; 
		return *this;
	}

	inline friend TComplexD operator +(const TComplexD&, const TComplexD&);
	inline friend TComplexD operator +(const double, const TComplexD&);
	inline friend TComplexD operator +(const TComplexD&, const double);
	inline friend TComplexD operator -(const TComplexD&, const TComplexD&);
	inline friend TComplexD operator -(const double, const TComplexD&);
	inline friend TComplexD operator -(const TComplexD&, const double);
	inline friend TComplexD operator *(const TComplexD&, const TComplexD&);
	inline friend TComplexD operator *(const double, const TComplexD&);
	inline friend TComplexD operator *(const TComplexD&, const double);
	inline friend TComplexD operator /(const TComplexD&, const TComplexD&);
	inline friend TComplexD operator /(const double, const TComplexD&);
	inline friend TComplexD operator /(const TComplexD&, const double);

/**

	inline friend TVector3d operator ^(const TVector3d&, const TVector3d&); // Vector product
	inline friend TVector3d operator *(const TMatrix3d&, const TVector3d&);

	inline friend double Abs(const TVector3d&);
	inline friend double NormAbs(const TVector3d&);
	inline friend int operator <(const TVector3d&, const TVector3d&);
	inline friend int operator ==(const TVector3d&, const TVector3d&);
	inline friend bool operator !=(const TVector3d&, const TVector3d&);
	inline friend bool operator >(const TVector3d&, const TVector3d&);

	inline friend int PracticallyEqual(const TVector3d&, const TVector3d&, double);
**/
};

//-------------------------------------------------------------------------

inline TComplexD operator +(const TComplexD& P1, const TComplexD& P2)
{
	TComplexD Res(P1.x + P2.x, P1.y + P2.y);
	return Res;
}
inline TComplexD operator +(const double a, const TComplexD& P)
{
	TComplexD Res(a + P.x, P.y);
	return Res;
}
inline TComplexD operator +(const TComplexD& P, const double a)
{
	TComplexD Res(a + P.x, P.y);
	return Res;
}
inline TComplexD operator -(const TComplexD& P1, const TComplexD& P2)
{
	TComplexD Res(P1.x - P2.x, P1.y - P2.y);
	return Res;
}
inline TComplexD operator -(const double a, const TComplexD& P)
{
	TComplexD Res(a - P.x, -P.y);
	return Res;
}
inline TComplexD operator -(const TComplexD& P, const double a)
{
	TComplexD Res(P.x - a, P.y);
	return Res;
}
inline TComplexD operator *(const TComplexD& P1, const TComplexD& P2)
{
	TComplexD Res(P1.x*P2.x - P1.y*P2.y, P1.x*P2.y + P1.y*P2.x);
	return Res;
}
inline TComplexD operator *(const double a, const TComplexD& P)
{
	TComplexD Res(a*P.x, a*P.y);
	return Res;
}
inline TComplexD operator *(const TComplexD& P, const double a)
{
	TComplexD Res(a*P.x, a*P.y);
	return Res;
}
inline TComplexD operator /(const TComplexD& P1, const TComplexD& P2)
{
	double InvAbsP2E2 = 1./(P2.x*P2.x + P2.y*P2.y);
	TComplexD Res((P1.x*P2.x + P1.y*P2.y)*InvAbsP2E2, (P1.y*P2.x - P1.x*P2.y)*InvAbsP2E2);
	return Res;
}
inline TComplexD operator /(const double a, const TComplexD& P)
{
	double InvAbsP2E2 = 1./(P.x*P.x + P.y*P.y);
	TComplexD Res(a*P.x*InvAbsP2E2, -a*P.y*InvAbsP2E2);
	return Res;
}
inline TComplexD operator /(const TComplexD& P, const double a)
{
	double Inv_a = 1./a;
	TComplexD Res(Inv_a*P.x, -Inv_a*P.y);
	return Res;
}


/**
inline TVector3d operator +(const TVector3d& P1, const TVector3d& P2)
{
	// The following can cause problems with Code Warrior
	return TVector3d(P1.x+P2.x, P1.y+P2.y, P1.z+P2.z);
}
inline TVector3d Plus(const TVector3d& P1, const TVector3d& P2)
{
	// The following can cause problems with Code Warrior
	return TVector3d(P1.x+P2.x, P1.y+P2.y, P1.z+P2.z);
}

//-------------------------------------------------------------------------

inline TVector3d operator -(const TVector3d& P1, const TVector3d& P2)
{
	// The following can cause problems with Code Warrior
	return TVector3d(P1.x-P2.x, P1.y-P2.y, P1.z-P2.z);
}
inline TVector3d Minus(const TVector3d& P1, const TVector3d& P2)
{
	// The following can cause problems with Code Warrior
	return TVector3d(P1.x-P2.x, P1.y-P2.y, P1.z-P2.z);
}

//-------------------------------------------------------------------------

inline TVector3d operator *(const double D, const TVector3d& P)
{
	// The following can cause problems with Code Warrior
	return TVector3d(D*P.x, D*P.y, D*P.z);
}
inline TVector3d Mult(const double D, const TVector3d& P)
{
	// The following can cause problems with Code Warrior
	return TVector3d(D*P.x, D*P.y, D*P.z);
}

//-------------------------------------------------------------------------

inline double operator *(const TVector3d& P1, const TVector3d& P2)
{
	return P1.x*P2.x+P1.y*P2.y+P1.z*P2.z;
}
inline double Mult(const TVector3d& P1, const TVector3d& P2)
{
	return P1.x*P2.x+P1.y*P2.y+P1.z*P2.z;
}

//-------------------------------------------------------------------------

inline TVector3d operator ^(const TVector3d& v1, const TVector3d& v2)
{
	// The following can cause problems with Code Warrior
	return TVector3d(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

//-------------------------------------------------------------------------

inline int operator <(const TVector3d& P1, const TVector3d& P2)
{
	return ((P1.x*P1.x + P1.y*P1.y + P1.z*P1.z) < (P2.x*P2.x + P2.y*P2.y + P2.z*P2.z));
}

//-------------------------------------------------------------------------

inline int operator ==(const TVector3d& P1, const TVector3d& P2)
{
	return ((P1.x == P2.x) && (P1.y == P2.y) && (P1.z == P2.z));
}

//-------------------------------------------------------------------------

inline bool operator !=(const TVector3d& P1, const TVector3d& P2)
{
	return ((P1.x != P2.x) || (P1.y != P2.y) || (P1.z != P2.z));
}

//-------------------------------------------------------------------------

inline bool operator >(const TVector3d& P1, const TVector3d& P2)
{
	return ((P1.x*P1.x + P1.y*P1.y + P1.z*P1.z) > (P2.x*P2.x + P2.y*P2.y + P2.z*P2.z));
}

//-------------------------------------------------------------------------

inline double Abs(const TVector3d& v)
{
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

//-------------------------------------------------------------------------

inline double NormAbs(const TVector3d& v)
{
	double vx = v.x, vy = v.y, vz = v.z;
	double Abs_vx = (vx>=0.)? vx : -vx;
	double Abs_vy = (vy>=0.)? vy : -vy;
	double Abs_vz = (vz>=0.)? vz : -vz;
	return (Abs_vx > Abs_vy)? ((Abs_vx > Abs_vz)? Abs_vx : Abs_vz) : ((Abs_vy > Abs_vz)? Abs_vy : Abs_vz);
}

//-------------------------------------------------------------------------

inline int PracticallyEqual(const TVector3d& v1, const TVector3d& v2, double Tol)
{
	return (::fabs(v2.x-v1.x) < Tol) && (::fabs(v2.y-v1.y) < Tol) && (::fabs(v2.z-v1.z) < Tol);
}
**/
//-------------------------------------------------------------------------

#endif
