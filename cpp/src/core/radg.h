/*-------------------------------------------------------------------------
*
* File name:      radg.h
*
* Project:        RADIA
*
* Description:    Base class for RADIA objects
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADG_H
#define __RADG_H

#include "radhandl.h"
#include "radauxst.h"
#include "gmvbstr.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTg;
class radTApplication;

//-------------------------------------------------------------------------

class radTg {
public:
	radTg() {}
	virtual ~radTg() {}

	virtual int Type_g() { return 0;}
	virtual void Dump(std::ostream& o, int ShortSign =0)
	{
		//o << "     Address: " << this << endl;
	}
	//virtual void DumpBin(CAuxBinStrVect& oStr, map<int, radTHandle<radTg>, less<int> >& mEl, radTHandle<radTg>& hg) {}
	//virtual void DumpBin(CAuxBinStrVect& oStr, map<int, radTHandle<radTg>, less<int> >& mEl, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey) {}
	virtual void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey) {}

	virtual int DuplicateItself(radTHandle<radTg>&, radTApplication*, char PutNewStuffIntoGenCont =1) { return 1;} // PutNewStuffIntoGenCont only for groups
	virtual int SizeOfThis() { return sizeof(*this);}
};

//-------------------------------------------------------------------------

typedef radTHandle<radTg> radThg;

#ifdef __GCC__
typedef map <int, radThg, less<int> > radTmhg;
typedef vector <radThg > radTvhg;
#else
typedef map <int, radThg, less<int> > radTmhg;
typedef vector <radThg, allocator<radThg> > radTvhg;
#endif

//-------------------------------------------------------------------------

struct radTPair_int_hg {
	int m;
	radThg Handler_g;

	radTPair_int_hg(int In_m, const radThg& hg) { m = In_m; Handler_g = hg; }
	radTPair_int_hg() {}

	int operator<(const radTPair_int_hg& r)
	{
		if(m<r.m) return 1;
		else return 0;
	}
	int operator==(const radTPair_int_hg& r)
	{
		if((m==r.m)&&(Handler_g==r.Handler_g)) return 1;
		else return 0;
	}

	inline friend bool operator!=(const radTPair_int_hg& r1, const radTPair_int_hg& r2);
	inline friend bool operator<(const radTPair_int_hg& r1, const radTPair_int_hg& r2);
	inline friend bool operator>(const radTPair_int_hg& r1, const radTPair_int_hg& r2);
};

//-------------------------------------------------------------------------

inline bool operator!=(const radTPair_int_hg& r1, const radTPair_int_hg& r2)
{
	return ((r1.m != r2.m) || (r1.Handler_g.rep != r2.Handler_g.rep));
}

//-------------------------------------------------------------------------

inline bool operator<(const radTPair_int_hg& r1, const radTPair_int_hg& r2)
{
	return (r1.m < r2.m);
}

//-------------------------------------------------------------------------

inline bool operator>(const radTPair_int_hg& r1, const radTPair_int_hg& r2)
{
	return (r1.m > r2.m);
}

//-------------------------------------------------------------------------

#endif

