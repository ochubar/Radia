
//+++++++++++++++++++++++++++++++++++++++
//
// Radia Globals
//
//++++++++++++++ Remarks++++++++++++++++++
//
// All global variables start with Rad
// All global functions and procs written in Igor start with Rad
// All functions written in C start with rad
//
//+++++++++++++++++++++++++++++++++++++++

#pragma rtGlobals=1		// Use modern global access method.

//+++++++++++++++++++++++++++++++++++++++
//
//Create and Intialize all global variables
//
//+++++++++++++++++++++++++++++++++++++++
proc RadInit(v)
Variable v=1
Silent 1				|	Initializing Radia  ...

String/G RadVerProc="4.23"
String VerExt=num2str(radUtiVer())

PathInfo Igor
newpath/O/Q RadUtiPath, S_path+"Radia"

DefaultFont "Times"

if(v==1)
	String info=IgorInfo(0)
	Print "Radia Procedure Version : ", RadVerProc
	Print "Radia DLL Version : ", VerExt
	Print ""
endif

if(v==1)
print " "
Print "Radia Examples are accessible from the \"Radia\" menu.\r  Reference to Radia functions is accessible via \"Radia->Help\" menu and via Igor Pro Help Browser (menu \"Help->Igor Help Browser\")."
print " "
endif

radFldLenTol(1.e-9, 1.e-9, 1.e-9)

end  // proc RadInit()
