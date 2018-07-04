/*-------------------------------------------------------------------------
*
* File name:      radexcep.h
*
* Project:        RADIA
*
* Description:    Auxiliary class: exceptions
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADEXCEP_H
#define __RADEXCEP_H

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTException {
public:	
	radTException() throw() {}
	radTException(const radTException&) throw() {}
	radTException& operator=(const radTException&) /*throw()*/ { return *this; }
//	virtual ~radTException() throw();
	virtual const char* what() const /*throw()*/{ return "";}
};
// later when all compilers support it radTBad_alloc will be derived from exception

//-------------------------------------------------------------------------

class radTBad_alloc : public radTException {
public:
	radTBad_alloc() throw() {}
	radTBad_alloc(const radTBad_alloc&) throw() {}
	radTBad_alloc& operator=(const radTBad_alloc&) /*throw()*/ { return *this; }
//	virtual ~radTBad_alloc() throw();
	virtual const char* what() const /*throw()*/{ return "Radia::Error900";}
};

//-------------------------------------------------------------------------

//static rad_bad_alloc nomem;

#ifdef _VC40
inline int radMemAllocHandler(size_t)
{
	radTBad_alloc Bad_alloc_exception;
	throw(&Bad_alloc_exception);
	return -1;
}
#else
inline void radMemAllocHandler()
{
	radTBad_alloc Bad_alloc_exception;
	throw(&Bad_alloc_exception);
}
#endif

//-------------------------------------------------------------------------

#endif
