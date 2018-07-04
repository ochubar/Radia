
#include "radtrans.h"

//-------------------------------------------------------------------------

//void radTrans::SetupRotation(const TVector3d& PoiOnAxVect, const TVector3d& InAxVect, double Angle)
//{
//	double NormFact = 1./sqrt(InAxVect.x*InAxVect.x+InAxVect.y*InAxVect.y+InAxVect.z*InAxVect.z);
//	TVector3d AxVect = NormFact*InAxVect;
//	double VxVx, VyVy, VzVz;
//	VxVx=AxVect.x*AxVect.x; VyVy=AxVect.y*AxVect.y; VzVz=AxVect.z*AxVect.z;
//
//	double cosAng, sinAng, One_m_cos;
//	cosAng = cos(Angle); sinAng = sin(Angle); One_m_cos = 1. - cosAng;
//	double One_m_cosVxVy, One_m_cosVxVz, One_m_cosVyVz, sinVx, sinVy, sinVz;
//	One_m_cosVxVy = One_m_cos*AxVect.x*AxVect.y;
//	One_m_cosVxVz = One_m_cos*AxVect.x*AxVect.z;
//	One_m_cosVyVz = One_m_cos*AxVect.y*AxVect.z;
//	sinVx = sinAng*AxVect.x; sinVy = sinAng*AxVect.y; sinVz = sinAng*AxVect.z;
//
//	TVector3d St0(VxVx+cosAng*(VyVy+VzVz), One_m_cosVxVy-sinVz, One_m_cosVxVz+sinVy);
//	TVector3d St1(One_m_cosVxVy+sinVz, VyVy+cosAng*(VxVx+VzVz), One_m_cosVyVz-sinVx);
//	TVector3d St2(One_m_cosVxVz-sinVy, One_m_cosVyVz+sinVx, VzVz+cosAng*(VxVx+VyVy));
//	M = TMatrix3d(St0, St1, St2);
//	M_inv = Matrix3d_inv(M);
//
//	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
//	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
//	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
//	TMatrix3d M0(St00, St01, St02);
//	V = M0*PoiOnAxVect;
//	detM = s = 1.;
//	ID_No = 2;
//}

//-------------------------------------------------------------------------
