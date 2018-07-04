/************************************************************************//**
 * File: gmvbstr.h
 * Description: Binary/Byte Stream utilities for Vector/Matrix classes (header)
 * Project: 
 * First release: 2013
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __GMVBSTR_H
#define __GMVBSTR_H

#include "auxparse.h"
#include "gmvectf.h"

//-------------------------------------------------------------------------

class CAuxBinStrVect : public CAuxBinStr {

public:

	CAuxBinStrVect(const unsigned char *bStr, int len_bStr) : CAuxBinStr(bStr, len_bStr) {}
	CAuxBinStrVect() : CAuxBinStr() {}

	inline friend CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TVector3d& v);
	inline friend CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TVector3df& v);
	inline friend CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TMatrix3d& m);
	inline friend CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TMatrix3df& m);
	inline friend CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TVector2d& v);
	inline friend CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TMatrix2d& m);

	inline friend CAuxBinStrVect& operator>>(CAuxBinStrVect& str, TVector3d& v);
	inline friend CAuxBinStrVect& operator>>(CAuxBinStrVect& str, TVector3df& v);
	inline friend CAuxBinStrVect& operator>>(CAuxBinStrVect& str, TMatrix3d& m);
	inline friend CAuxBinStrVect& operator>>(CAuxBinStrVect& str, TMatrix3df& m);
};

inline CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TVector3d& v)
{
	str << v.x << v.y << v.z;
	return str;
};
inline CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TVector3df& v)
{
	str << v.x << v.y << v.z;
	return str;
};
inline CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TMatrix3d& m)
{
	str << m.Str0 << m.Str1 << m.Str2;
	return str;
};
inline CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TMatrix3df& m)
{
	str << m.Str0 << m.Str1 << m.Str2;
	return str;
};

inline CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TVector2d& v)
{
	str << v.x << v.y;
	return str;
};
inline CAuxBinStrVect& operator<<(CAuxBinStrVect& str, TMatrix2d& m)
{
	str << m.Str0 << m.Str1;
	return str;
};

inline CAuxBinStrVect& operator>>(CAuxBinStrVect& str, TVector3d& v)
{
	str >> v.x; str >> v.y; str >> v.z;
	return str;
};
inline CAuxBinStrVect& operator>>(CAuxBinStrVect& str, TVector3df& v)
{
	str >> v.x; str >> v.y; str >> v.z;
	return str;
};
inline CAuxBinStrVect& operator>>(CAuxBinStrVect& str, TMatrix3d& m)
{
	str >> m.Str0; str >> m.Str1; str >> m.Str2;
	return str;
};
inline CAuxBinStrVect& operator>>(CAuxBinStrVect& str, TMatrix3df& m)
{
	str >> m.Str0; str >> m.Str1; str >> m.Str2;
	return str;
};

//-------------------------------------------------------------------------

#endif
