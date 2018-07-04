
#ifndef __RADSTLON_H
#define __RADSTLON_H

//-------------------------------------------------------------------------

#if defined _MSC_VER 
#pragma warning(disable : 4786) // to suppress annoying warning from STL
#include <map>
//#include <strstream>
#include <sstream>
#include <iostream>
using namespace std;

#elif defined  __HP_aCC
#include <map>
#include <strstream.h>
#include <iostream.h>

#elif defined __MWERKS__
#include <map.h>
#include <strstream.h>
#include <iostream.h>

#else
#include <map>
//#include <strstream>
#include <sstream>
#include <iostream>
using namespace std;
#endif

//-------------------------------------------------------------------------

#endif
