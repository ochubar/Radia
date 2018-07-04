/*-------------------------------------------------------------------------
*
* File name:      radiobuf.cpp
*
* Project:        RADIA
*
* Description:    Input/output buffer; errors and warnings
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radstlon.h"
#include "radiobuf.h"

//-------------------------------------------------------------------------

int radTIOBuffer::AmOfErrors = 137; //modify this when adding new error !!!
string radTIOBuffer::err_ar[] = {

	"Radia::ErrorXXX::::Wrong Error Number.\0",
	"Radia::Error000::::Incorrect input: Function arguments.\0",
	"Radia::Error001::::Incorrect input: Block dimensions should be positive.\0",
	"Radia::Error002::::Incorrect input: Element with that key does not exist.\0",
	"Radia::Error003::::Incorrect input: Element with that key is not a 3D Object.\0",
	"Radia::Error004::::Incorrect input: Element with that key is not a Container.\0",
	"Radia::Error005::::Incorrect input: Container can not be added into itself.\0",
	"Radia::Error006::::Incorrect input: Element with that key is not a Transformation.\0",
	"Radia::Error007::::Incorrect input: Element with that key is not a Rectangular parallelepiped block.\0",
	"Radia::Error008::::Incorrect input: RGB color values should be between 0 and 1.\0",
	"Radia::Error009::::Incorrect input: Field identification string should contain either combination of characters B,b,H,h,A,a,M,m,x,y,z or nothing.\0",
	"Radia::Error010::::Incorrect input: Arc radii {rmin,rmax} should be two positive values, rmin < rmax.\0",
	"Radia::Error011::::Incorrect input: Arc angles {phimin,phimax} should be between 0 and 2*Pi, phimin < phimax.\0",
	"Radia::Error012::::Incorrect input: Height should be positive.\0",
	"Radia::Error013::::Incorrect input: Object with that key has no drawing attributes applied.\0",
	"Radia::Error014::::Incorrect input: Computation criterium string should be either  auto  or  man .\0",
	"Radia::Error015::::Incorrect input: Element with that key is neither a relaxable 3D object nor a Container.\0",
	"Radia::Error016::::Incorrect input: Element with that key is not a Material.\0",
	"Radia::Error017::::Incorrect input: Element with that key is not an Interaction.\0",
	"Radia::Error018::::Incorrect input: Relaxation Parameter value should be between 0 and 1.\0",
	"Radia::Error019::::Incorrect input: Number of iterations should be positive.\0",
	"Radia::Error020::::Incorrect input: Interaction vector identification string should be either  ext ,  tot  or  mag .\0",
	"Radia::Error021::::Incorrect input: Subdivision list should consist of 3 positive integer values: {kx,ky,kz}.\0",
	"Radia::Error022::::Incorrect input: Magnetic susceptibility list should consist of 2 real values: {ksipar,ksiper}.\0",
	"Radia::Error023::::Incorrect input: Non-zero Remanent Magnetization should be specified to determine the Easy Magnetization axis.\0",
	"Radia::Error024::::Incorrect input: Lists of the Material defining constants should be of equal lengths, not exceeding 3.\0",
	"Radia::Error025::::Incorrect input: Element with that key is neither Material nor 3D object with Material applied.\0",
	"Radia::Error026::::Incorrect input: Magnetization identification string should contain either combination of characters M,m,X,x,Y,y,Z,z or nothing.\0",
	"Radia::Error027::::Incorrect input: Material was not applied to that 3D object.\0",
	"Radia::Error028::::Incorrect input: Method number should be 0 or 1 or 2 or 3.\0",
	"Radia::Error029::::Incorrect input: Racetrack straight part dimensions should be set up by two positive values {lx,ly}.\0",
	"Radia::Error030::::Incorrect input: Absolute precision value should be positive.\0",
	"Radia::Error031::::Incorrect input: Maximum number of iterations should be positive.\0",
	"Radia::Error032::::Incorrect input: Field Integral component identification string should contain either combination of characters I,i,B,b,H,h,x,y,z or nothing.\0",
	"Radia::Error033::::Incorrect input: Integral identification string should be either  inf  or  fin .\0",
	"Radia::Error034::::Incorrect input: Argument flag string should be either  noarg  or  arg .\0",
	"Radia::Error035::::Incorrect input: List of the Filament Conductor edge points can not contain less than two points.\0",
	"Radia::Error036::::Incorrect input: At Force computation, the Shape element should be Rectangle or Rectangular Parallelepiped (RecMag).\0",
	"Radia::Error037::::Incorrect input: Element with that key is not a Subdivided Rectangular Parallelepiped.\0",
	"Radia::Error038::::Incorrect input: (At present,) This Field Computation Method is implemented only for subdivision not higher than {3,3,3}.\0",
	"Radia::Error039::::Incorrect input: Sub-element index out of subdivision limits.\0",
	"Radia::Error040::::Incorrect input: No subdivided rectangular parallelepipeds (RecMags) present at the specified subdivision level.\0",
	"Radia::Error041::::Incorrect input: Automatic Relaxation is only available with Methods 3 and 4.\0",
	"Radia::Error042::::Incorrect input: Number of points for Focusing Potential computation should be larger than 2.\0",
	"Radia::Error043::::Incorrect input: The function input string should be  on  or  off .\0",
	"Radia::Error044::::Incorrect input: Randomization type identification string should be  rel  or  abs .\0",
	"Radia::Error045::::Incorrect input: Field computation mode string should be either  man  or  auto .\0",
	"Radia::Error046::::Incorrect input: Memory allocation method identification string should be  tot  or  parts .\0",
	"Radia::Error047::::Incorrect input: The polyhedron's vertex points specified as defining one face do not belong to one plane.\0",
	"Radia::Error048::::Incorrect input: Polyhedron's face can not contain less than 3 vertex points.\0",
	"Radia::Error049::::Incorrect input: Can not find at least three vertex points not belonging to one line in the face of the polyhedron.\0",
	"Radia::Error050::::Incorrect input: Improper definition of the polyhedron: No face found with all the rest vertex points located from one side of it.\0",
	"Radia::Error051::::Incorrect input: Improper cutting plane definition. The plane should be defined by a point and a normal, in the form {{x,y,z},{nx,ny,nz}}.\0",
	"Radia::Error052::::Incorrect input: Improper definition of the subdivision parameters. See the radObjDivMag function template.\0",
	"Radia::Error053::::Incorrect input: Subdivision parameters should be positive.\0",
	"Radia::Error054::::Incorrect input: Subdivision frame identification string should be  loc  (for the subdivision in local frames of the 3D objects) or  lab  (for the subdivision in the laboratory frame).\0",
	"Radia::Error055::::Incorrect input: Coils subdivision flag string should be  SubdCoils  or  NotSubdCoils .\0",
	"Radia::Error056::::Incorrect input: Energy or Force (Torque) component identification string should contain either combination of characters E,e,F,f,T,t,x,y,z or nothing.\0",
	"Radia::Error057::::Incorrect input: Improper definition of absolute precision value(s). See the usage message for the function.\0",
	"Radia::Error058::::Incorrect input: Force component identification string should contain either combination of characters F,f,X,x,Y,y,Z,z or nothing.\0",
	"Radia::Error059::::Incorrect input: Torque component identification string should contain either combination of characters T,t,X,x,Y,y,Z,z or nothing.\0",
	"Radia::Error060::::Incorrect input: Improper Nonlinear Anisotropic material definition: No nonlinear dependence of M vs H detected. Use radMatLin function to define a linear material.\0",
	"Radia::Error061::::Incorrect input: Improper Nonlinear Anisotropic material definition: Coercivity parameter is not 0 for the component perpendicular to the easy magnetization axis.\0",
	"Radia::Error062::::Incorrect input: Improper definition of optional parameters of the function. See usage message or Reference Guide record for the function.\0",
	"Radia::Error063::::Multiple extruded polygon can not be generated from this input: incorrect definition of horizontal slice polygons.\0",
	"Radia::Error064::::Incorrect input: Rectangle dimensions should be positive.\0",
	"Radia::Error065::::Incorrect input: Improper definition of additional specifications for the subdivision. See usage message or Reference Guide record for the function.\0",
	"Radia::Error066::::Incorrect input: The vector defining cylindrical subdivision axis should not be zero.\0",
	"Radia::Error067::::Incorrect input: The subdivision ellipse exes ratio should be defined by a positive real number.\0",
	"Radia::Error068::::Incorrect input: Too small azimuthal subdivision number or too large averaged azimuthal size of pieces to be produced by cylindrical subdivision.\0",
	"Radia::Error069::::Incorrect input: Improper definition of directions for subdivision by parallel planes.\0",
	"Radia::Error070::::Incorrect input: Field component identification character should be one of the following: A,a,B,b,H,h,J,j,M,m.\0",
	"Radia::Error071::::Incorrect input: For proper definition of Nonlinear Isotropic material, the absolute magnetization should be tabulated vs absolute magnetic field strength (with no coercitivity), as {{h1,m1},{h2,m2},...}, where h1 < h2 < ...\0",
	"Radia::Error072::::Incorrect input: Material name was not specified.\0",
	"Radia::Error073::::Incorrect input: Material with this name was not found in the internal database.\0",
	"Radia::Error074::::Incorrect input: Laminated material stacking factor should be a positive number less than or equal to 1.\0",
	"Radia::Error075::::Incorrect input: Incomplete or incorrect definition of magnetic material.\0",
	"Radia::Error076::::Incorrect input: Period should be positive.\0",
	"Radia::Error077::::Incorrect input: Number or periods should be positive.\0",
	"Radia::Error078::::Incorrect input: Range of position should be positive.\0",
	"Radia::Error079::::Incorrect input: Number of points should be positive.\0",
	"Radia::Error080::::Incorrect input: Minimum number of points for harminic analysis should be 8.\0",
	"Radia::Error081::::Incorrect input: Number or magnetic field harmonics should be positive.\0",
	"Radia::Error082::::Incorrect input: Non-zero vector is expected.\0",
	"Radia::Error083::::Incorrect input: Non-zero array is expected.\0",
	"Radia::Error084::::Incorrect input: Number of transverse points should be less than 100.\0",
	"Radia::Error085::::Incorrect input: Number of longitudinal positions should be less than 50.\0",
	"Radia::Error086::::Incorrect input: Array of longitudinal positions should represent an increasing sequence.\0",
	"Radia::Error087::::Incorrect input: The two vectors can not be parallel.\0",
	"Radia::Error088::::Incorrect input: Nonlinear magnetic material definition requires more than 3 pairs of values (H,M).\0",
	"Radia::Error089::::Incorrect input: Radius should be positive.\0",
	"Radia::Error090::::Incorrect input: Number of segments should be larger than 2.\0",
	"Radia::Error091::::Incorrect input: This relaxation method requires a Container object.\0",
	"Radia::Error092::::Incorrect input: The character defining the object orientation should be either X or Y or Z.\0",
	"Radia::Error093::::Incorrect input: The string defining the symmetry state of the object should be either \"sym\" or \"nosym\".\0",
	"Radia::Error094::::Incorrect input: The angular size of one segment should be larger than Pi.\0",
	"Radia::Error095::::Incorrect input: Incorrect definition of radial coordinates.\0",
	"Radia::Error096::::Incorrect input: Field identification string should contain either combination of characters B,b,H,h,A,a,M,m,x,y,z or nothing; Field Integral component identification string should contain either combination of characters I,i,B,b,H,h,x,y,z or nothing.\0",
	"Radia::Error097::::Incorrect input: Subdivision parameter list length should be equal to the number of segments in the polygon border.\0",
	"Radia::Error098::::Incorrect input: The constraint on minimal angle in triangles should not exceed 35 degrees.\0",
	"Radia::Error099::::Incorrect input: Homothety coefficients must have same sign.\0",
	"Radia::Error100::::No elements are present in memory.\0",
	"Radia::Error101::::No 3D objects are present in memory.\0",
	"Radia::Error102::::No Relaxable 3D objects with Material applied are present in the entity.\0",
	"Radia::Error103::::Failed to subdivide the Polygon. Try to change the subdivision parameters or the Polygon shape.\0",
	"Radia::Error104::::Self-intersecting polygon encountered. Self-intersecting polygons are not supported by Radia.\0",
	"Radia::Error105::::Unclosed polyhedron encountered. Unclosed polyhedrons can not be used as basic 3D objects in Radia.\0",
	"Radia::Error106::::Non-convex polyhedron encountered. Non-convex polyhedrons can not be used as basic 3D objects in Radia. Try to represent the volume by a group of several convex polyhedrons.\0",
	"Radia::Error107::::The anisotropic material can not be applied to the object since neither the material nor the object have the easy magnetization axis defined.\0",
	"Radia::Error108::::Subdivision of this object can not be performed in the Laboratory frame. The object can only be subdivided in its Local frame.\0",
	"Radia::Error109::::The requested procedure is not implemented for coils.\0",
	"Radia::Error110::::Failed to generate convex polyhedron(s) from the given input.\0",
	"Radia::Error111::::Multiple extruded polygon can not be generated from this input: Non-convex horizontal slice polygon encountered.\0",
	"Radia::Error112::::Multiple extruded polygon can not be generated from this input: More than three points belonging to one line encountered in the horizontal slice polygon.\0",
	"Radia::Error113::::Can not start QuickDraw 3D library. Please check whether the QuickDraw 3D 1.5 or later from Apple Computer is properly installed on your system.\0",
	"Radia::Error114::::This function is not implemented on that platform.\0",
	"Radia::Error115::::Can not start Glut OpenGL shared library.\0",
	"Radia::Error116::::Failed to create Interaction Matrices for solving by parts.\0",
	"Radia::Error117::::Incorrectly defined polyhedron encountered. Note that non-convex polyhedrons are not supported in this version of the code.\0",
	"Radia::Error118::::Failed to create Interaction Matrix.\0",
	"Radia::Error119::::Triangulation failed.\0",
	"Radia::Error120::::In this Radia version, Magnetization and Current Density can not be defined simultaneously for the same 3D object.\0",
	"Radia::Error121::::In this Radia version, linearly varying Magnetization in a 3D object is not supported.\0",
	"Radia::Error122::::Failed to find finite Current Density for given polyhedron shape current-carrying oblect.\0",
	"Radia::Error123::::Multiple extruded polygon can not be generated from this input: incorrect definition of transformations at extrusion step(s).\0",
	"Radia::Error124::::Multiple extruded polygon can not be generated from this input: an extrusion step can not consist of a single homothety without any other transformations.\0",
	"Radia::Error125::::Failed to generate 3D object from the given input.\0",
	"Radia::Error200::::Step size is too small in automatic Runge-Kutta integration routine.\0",
	"Radia::Error201::::Maximum number of steps exceeded in automatic Runge-Kutta integration routine.\0",
	"Radia::Error202::::Failed to instantiate object(s).\0",
	"Radia::Error500::::Incorrect input: Byte string is expected.\0", //keep on adding new "Incorrect inputs" after this
	"Radia::Error501::::Incorrect input: Wrong / unsupported magnetic kick units.\0", //keep on adding new "Incorrect inputs" after this
	"Radia::Error502::::Incorrect input: Wrong / unsupported kick-map string/file format specification.\0", //keep on adding new "Incorrect inputs" after this
	"Radia::Error900::::Memory allocation failure.\0",
	"Radia::Error990::::Graphical presentation of this element is not available.\0",
	"Radia::Error998::::Execution aborted by User.\0",
	"Radia::Error999::::Unidentified error.\0",
	//insert new error here
};

//-------------------------------------------------------------------------

int radTIOBuffer::AmOfWarnings = 10; //modify this when adding new warning !!!
string radTIOBuffer::warn_ar[] = {

	"Radia::WarningXXX::::Wrong Warning Number.\0",
	"Radia::Warning010::::Polygon's center outside the polygon or improper polygon definition.\0",
	"Radia::Warning011::::Non-convex polygon or improper polygon definition.\0",
	"Radia::Warning012::::Non-convex polyhedron or improper polyhedron definition.\0",
	"Radia::Warning013::::Non-closed polyhedron or improper polyhedron definition.\0",
	"Radia::Warning014::::The dependence of magnetization vs field strength is not monotone in the material definition. This may affect the accuracy of relaxation.\0",
	"Radia::Warning015::::The specified accuracy was not reached. This may indicate diverging or incomplete relaxation.\0",
	"Radia::Warning016::::Polygon edge points seem to coinside. Certain algorithms may not work for such polygon definition.\0",
	"Radia::Warning017::::A symbol was supplied to field computation function, in place where a number was expected; the symbol was replaced by 0.\0",
	"Radia::Warning018::::Current density scaling was not performed because no current-carrying 3D objects were found.\0",
	//insert new warning here
};

//-------------------------------------------------------------------------

