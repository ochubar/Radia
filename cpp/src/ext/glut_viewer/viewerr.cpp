
#ifndef __VIEWERROR_H
#include "viewerr.h"
#endif

//-------------------------------------------------------------------------

int CErrWarn::m_BaseErrNo = 0;
int CErrWarn::m_BaseWarNo = 0;
vector<int> CErrWarn::m_WarnNosVect;

//-------------------------------------------------------------------------

string CErrWarn::m_error[] = {
	"Wrong error number.\0",
	"Incorrect function arguments.\0", //1
	"Function argument(s) are not defined correctly.\0", //2
};

//-------------------------------------------------------------------------

string CErrWarn::m_warning[] = {
	"Wrong warning number.\0",
};

//-------------------------------------------------------------------------
