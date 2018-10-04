/////////////////////////////////////////////////////////////////////////////
// Name:        radentry.h
// Purpose:     Radia DLL header file
// Authors:     O.Chubar, P.Elleaume
// Version:     DLL 4.115
// Modified:    30.05.04
/////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------
// Platform- and compiler-dependent macro definitions
//-------------------------------------------------------------------------
#if !(defined(ALPHA_NONE) || defined(ALPHA__LIB__))
/*---------------For CodeWarrior PowerMac---------------*/
#if defined __POWERPC__
#if defined ALPHA__DLL__ || defined MATLAB_MEX_FILE
#define EXP __declspec(export)
#endif
/*---------------For CodeWarrior PC and Visual C++---------------*/
#elif defined __INTEL__ || defined WIN32
#if defined ALPHA__DLL__ || defined MATLAB_MEX_FILE
#define EXP __declspec(dllexport)
#else
#define EXP __declspec(dllimport)
#endif
//#define CALL __stdcall
#define CALL __cdecl //use one of these calling convntions at preparing DLL
/*---------------For HP-UX, gcc---------------*/
#else
#endif
#endif /*ALPHA_NONE*/

#ifndef EXP
#define EXP
#endif

#ifndef CALL
#define CALL
#endif

//-------------------------------------------------------------------------

#ifdef __cplusplus  
extern "C" {
#endif

/* OK = 0 */
#define OK 0

/** Returns the error message associated to error number. 
This function cannot be called from Visual Basic For Applications. 
@param er [in] error number 
@return the chain of characters representing the error string
@author P. Elleaume 
@version 1.0 
@see GetErrorSize and GetErrorText 
*/ 
EXP const char* CALL RadErrGet(int er); 

/** Returns the length of the error message string not counting "\0". 
@param siz [out] length of the error message  
@param er [in] error number 
@return integer error code (0 : No Error, >0 : Error Number, <0 : Warning Number) 
@author P. Elleaume 
@version 1.0 
*/ 
EXP int CALL RadErrGetSize(int* siz, int er);

/** Returns the text of error message associated to error number. 
@param t [out] error message string
@param er [in] error number 
@return integer error code (0 : No Error, >0 : Error Number, <0 : Warning Number)
@author P. Elleaume 
@version 1.0 
*/ 
EXP int CALL RadErrGetText(char* t, int er);

/** Returns the warning message associated to warning number. 
This function cannot be called from Visual Basic For Applications. 
@param er [in] warning number 
@return the chain of characters representing the warning string
@author P. Elleaume 
@version 1.0
@see GetWarningSize and GetWarningText 
*/ 
EXP const char * CALL RadWarGet(int er); 

/** Returns the length of the warning message not counting "\0". 
@param siz [out] length of the warning message  
@param er [in] warning number 
@return	integer error code (0 : No Error, >0 : Error Number, <0 : Warning Number)
@author P. Elleaume 
@version 1.0 
@see dllGetWarningText
*/ 
EXP int CALL RadWarGetSize(int* siz, int er);

/** Returns the text of warning message associated to error number. 
@param t [out] warning message string
@param er [in] error number 
@return integer error code (0 : No Error, >0 : Error Number, <0 : Warning Number)
@author P. Elleaume 
@version 1.0 
@see GetWarningSize
*/ 
EXP int CALL RadWarGetText(char* t, int er);

/** Creates a uniformly magnetized rectangular parallelepiped.
The parallelepiped block is defined through its center point P[3], dimensions L[3], and magnetization M[3]."
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the block center of gravity
@param L [in] array of 3 dimensions of the block
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjRecMag(int* n, double* P, double* L, double* M);

/** Creates a uniformly magnetized extruded polygon.
The extrusion axis is directed along X axis of the laboratory frame.
@param n [out] reference number of the object created
@param xc [in] the horizontal coordinate of the block center of gravity
@param lx [in] the thickness (extrusion size)
@param FlatVert [in] flat array of y and z coordinates (y1, z1, y2, z2,...) of vertex points describing the polygon in 2D
@param nv [in] number of vertex points of the 2D polygon (the length of the FlatVert array is 2*nv)
@param a [in] character identifying extrusion direction (can be 'x', 'y' or 'z')
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjThckPgn(int* n, double xc, double lx, double* FlatVert, int nv, char a, double* M);

/** Creates a uniformly magnetized polyhedron (closed volume limited by planes).
@param n [out] reference number of the object created
@param FlatVert [in] flat array of x, y and z coordinates (x1, y1, z1, x2, y2, z2,...) of the polyhedron vertex points
@param nv [in] number of vertex points of the polyhedron (the length of the FlatVert array is 3*nv)
@param FlatFaces [in] flat array of indexes of vertex points defining the polyhedron faces (f1i1, f1i2,..., f2i1, f2i2,...)
@param FacesLen [in] array of integer numbers equal to the numbers of vertex points in each face of the polyhedron; the order of the faces is the same as in the FlatFaces array (the length of the FlatFaces array is equal to the sum of elements of the FacesLen array) 
@param nf [in] number of faces of the polyhedron (or the length of the FacesLen array)
@param M [in] array of 3 cartesian components of magnetization vector inside the block
@param M_LinCoef [in] array of 9 coefficients of linearly-varying magnetization vector inside the block
@param J [in] array of 3 cartesian components of current density vector inside the block
@param J_LinCoef [in] array of 9 coefficients of linearly-varying current density vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjPolyhdr(int* n, double* FlatVert, int nv, int* FlatFaces, int* FacesLen, int nf, double* M, double* M_LinCoef, double* J, double* J_LinCoef);

/**Creates a uniformly magnetized finite-length arc of polygonal cross-section.
@param n [out] reference number of the object created
@param P [in] array of 2 cartesian coordinates defining the position of the rotation axis in the plane perpendicular to this axis
@param a [in] character defining the orientation of the rotation axis in 3D space; it can be either to 'x', 'y' or 'z'
@param FlatVert [in] flat array of radial and axial coordinates (r1, z1, r2, z2,...) of vertex points of the cross-section polygon
@param nv [in] number of vertex points of the cross-section polygon (the length of the FlatVert array is 2*nv)
@param Phi [in] array of 2 numbers - initial and final azimuth angles
@param nseg [in] number of segments in the arc
@param sym_no [in] character which can be either 's' or 'n'; depending on the value of this switch, the magnetization vectors in nseg sector polyhedrons either will be linked by rotational symmetry ('s'), or will behave as independent magnetization vectors at any subsequent relaxation
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjArcPgnMag(int* n, double* P, char a, double* FlatVert, int nv, double* Phi, int nseg, char sym_no, double* M);

/** ?? TO IMPLEMENT ??
Creates a uniformly magnetized volume obtained by rotation of a planar convex polygon over 2 Pi around pre-defined axis.
@param n [out] reference number of the object created
@param P [in] array of 2 cartesian coordinates defining the position of the rotation axis in the plane perpendicular to this axis
@param FlatVert [in] flat array of radial and axial coordinates (r1, z1, r2, z2,...) of vertex points of the cross-section polygon
@param nv [in] number of vertex points of the cross-section polygon (the length of the FlatVert array is 2*nv)
@param nseg [in] number of segments in the arc
@param sym_no [in] character which can be either 's' or 'n'; depending on the value of this switch, the magnetization vectors in nseg sector polyhedrons either will be linked by rotational symmetry ('s'), or will behave as independent magnetization vectors at any subsequent relaxation
@param a [in] character defining the orientation of the arc rotation axis; it can be either to 'x', 'y' or 'z'
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
//EXP int CALL RadObjCircPgnMag(int* n, double* P, double* FlatVert, int nv, int nseg, char sym_no, char a, double* M);

/** Attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on slices.
The slices are assumed to be convex planar polygons parallel to the XY plane.
@param n [out] reference number of the object created
@param FlatVert [in] flat array of x and y coordinates (x11, y11, x12, y12,..., x21, y21, x22, y22,...) of vertex points of the slices (planar polygons)
@param SlicesLen [in] array of integer numbers equal to the numbers of vertex points in each slice; the order of the slices is the same as in the FlatVert array (the length of the FlatVert array is twice the sum of elements of the SlicesLen array) 
@param Attitudes [in] array of vertical coordinates (or attitudes) of the slices, in the ascending order (which is the same as for the SlicesLen array)
@param ns [in] number of slices (or the length of the Attitudes and SlicesLen arrays)
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjMltExtPgn(int* n, double* FlatVert, int* SlicesLen, double* Attitudes, int ns, double* M);

/** Attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on rectangular slices.
The rectangular slices are assumed to be parallel to the XY plane.
@param n [out] reference number of the object created
@param FlatCenPts [in] flat array of x, y and z coordinates (x1, y1, z1, x2, y2, z2,...) of the slices center points
@param FlatRtgSizes [in] flat array of sizes of the slice rectangles along x and y (wx1, wy1, wx2, wy2,...)
@param ns [in] number of slices (the length of the FlatCenPts array is 3*ns, the length of the FlatRtgSizes array is 2*ns)
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjMltExtRtg(int* n, double* FlatCenPts, double* FlatRtgSizes, int ns, double* M);

/** Creates triangulated extruded polygon block, i.e. an extruded polygon with its bases subdivided by triangulation.
@param n [out] reference number of the object created
@param xc [in] the horizontal coordinate of the block center of gravity
@param lx [in] the thickness (extrusion size)
@param FlatVert [in] flat array of y and z coordinates (y1, z1, y2, z2,...) of vertex points describing the polygon in 2D
@param FlatSubd [in] flat array of subdivision parameters for base polygon segments, two numbers for each segment: the first defining number of sub-segments, the second - "gradient" of the segmentation
@param nv [in] number of vertex points of the 2D polygon (the length of the FlatVert and FlatSubd arrays is 2*nv)
@param a [in] character identifying extrusion direction (can be 'x', 'y' or 'z')
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@param opt [in] pointer to options string, which can be e.g. "ki->...,TriAngMin->...,TriAreaMax->...,ExtOpt->..." or each of these tokens separately, or 0; default values are "ki->Numb" (rather than "ki->Size"), "TriAngMin->20"
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjMltExtTri(int* n, double xc, double lx, double* FlatVert, double* FlatSubd, int nv, char a, double* M, char* opt);

/** Creates a finite-length arc magnet of rectangular cross-section.
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the arc center point
@param R [in] array of 2 numbers - inner and outer radii
@param Phi [in] array of 2 numbers - initial and final azimuth angles
@param h [in] height
@param nseg [in] number of segments
@param a [in] character defining the orientation of the rotation axis of the arc; it can be either to 'x', 'y' or 'z'
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
//EXP int CALL RadObjArcMag(int* n, double* P, double* R, double* Phi, double h, int nseg, char a, double* M);

/** Creates a cylindrical magnet.
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the cylinder center point
@param r [in] cylinder radius
@param h [in] cylinder height
@param nseg [in] number of segments
@param a [in] character defining the orientation of the cylinder axis; it can be either to 'x', 'y' or 'z'
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjCylMag(int* n, double* P, double r, double h, int nseg, char a, double* M);

/** Creates a current carrying rectangular parallelepiped block.
The parallelepiped block is defined through its center point P[3], dimensions L[3], and current density vector J[3]."
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the block center of gravity
@param L [in] array of 3 dimensions of the block
@param J [in] array of 3 cartesian coordinates of the current density vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjRecCur(int* n, double* P, double* L, double* J);

/** Creates a current carrying finite-length arc of rectangular cross-section.
The arc rotation axis is directed along Z.
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the arc center point
@param R [in] array of 2 numbers - inner and outer radii
@param Phi [in] array of 2 numbers -  initial and final azimuth angles 
@param h [in] height
@param nseg [in] number of segments
@param man_auto [in] character which can be either 'm' or 'a'; the magnetic field from the arc is then computed based on the number of segments nseg ('m', i.e. "manual"), or on the general absolute precision level specified by the functions RadFldCmpCrt or RadFldCmpPrc ('a', i.e. "automatic")
@param a [in] character defining the orientation of the rotation axis of the arc; it can be either 'x', 'y' or 'z'
@param j [in] azimuthal current density
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjArcCur(int* n, double* P, double* R, double* Phi, double h, int nseg, char man_auto, char a, double j);

/** Creates a current carrying racetrack coil.
The coil consists of four 90-degree bents connected by four straight parts parallel to the XY plane.
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the racetrack center point
@param R [in] array of 2 numbers - inner and outer radii
@param L [in] array of 2 numbers - straight section lengths
@param h [in] height
@param nseg [in] number of segments
@param man_auto [in] character which can be either 'm' or 'a'; the magnetic field from the arc is then computed based on the number of segments nseg ('m', i.e. "manual"), or on the general absolute precision level specified by the functions RadFldCmpCrt or RadFldCmpPrc ('a', i.e. "automatic")
@param a [in] character defining the orientation of the rotation axis of the arc; it can be either 'x', 'y' or 'z'
@param j [in] azimuthal current density
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjRaceTrk(int* n, double* P, double* R, double* L, double h, int nseg, char man_auto, char a, double j);

/** Creates a filament polygonal line conductor with current.
The line conductor is defined by sequence of points in 3D space.
@param n [out] reference number of the object created
@param FlatPts [in] flat array of x, y and z coordinates of the points (x1, y1, z1, x2, y2, z2,...)
@param np [in] number of points (the length of the array FlatPts is 3*np)
@param i [in] current
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjFlmCur(int* n, double* FlatPts, int np, double i);

/** Scales current (density) in a 3D object by multiplying it by a constant.
@param n [out] reference number of the object with current (density) to be scaled
@param scaleCoef [in] scaling coefficient
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjScaleCur(int n, double scaleCoef);

/** Creates a source of uniform background magnetic field.
@param n [out] reference number of the object created
@param B [in] array of 3 cartesian coordinates of the magnetic field vector
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjBckg(int* n, double* B);

/** Creates a container of magnetic field sources.
@param n [out] reference number of the object created
@param Objs [in] array of reference numbers of the objects (n1, n2, n3,...) to put to the container
@param nobj [in] number of objects (the length of the array Elems)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjCnt(int* n, int* Objs, int nobj);

/** Adds objects to the container object cnt.
@param cnt [in] reference number of the container object
@param Objs [in] array of reference numbers of the objects (n1, n2, n3,...) to put to the container
@param nobj [in] number of objects to put to the container (the length of the array Elems)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjAddToCnt(int cnt, int* Objs, int nobj);

/** Calculates the number of objects in the container cnt.
@param n [out] calculated number of objects
@param cnt [in] reference number of the container object
@param deep [in] switch specifying whether all atomic elements of eventual member containers have to be counted (1) or not (0, default)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjCntSize(int* n, int cnt);

/** Fills-in an array with the reference numbers of objects present in the container cnt. The array should be allocated in the calling application. The necessary size of the array can be determined using the function RadObjCntSize.
@param Objs [out] array of reference numbers of the objects present in the container
@param cnt [in] reference number of the container object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjCntStuf(int* Objs, int cnt);

/** Duplicates the object obj. 
@param n [out] reference number of the object created
@param obj [in] reference number of the object to duplicate
@param opt [in] pointer to an option string, which can be "FreeSym->False", "FreeSym->True" or 0. This specifies whether the symmetries (transformations with multiplicity more than one) previously applied to the object obj should be simply copied at the duplication ("FreeSym->False" or 0), or a container of new independent objects should be created in place of any symmetry previously applied to the object obj. In both cases the final object created by the duplication has exactly the same geometry as the initial object obj.
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDpl(int* n, int obj, char* opt);

/** Provides coordinates of geometrical center point(s) and magnetization(s) of the object obj.
@param M [out] flat array of resulting magnetization array. If this pointer is 0 at input, the actual output data array can be obtained by calling RadUtiDataGet function. This aray will contain coordinates of geometrical center point(s) and magnetic field components. If obj is a container, the array will include the container members' center points and their magnetic field components. 
@param arMesh [out] flat array defining the structure of resulting magnetization array, i.e. number of dimensions (arMesh[0]) and number of values per dimension. The actual output data array can be obtained by calling RadUtiDataGet function. This aray will contain coordinates of geometrical center point(s) and magnetic field components. If obj is a container, the array will include the container members' center points and their magnetic field components. 
@param obj [in] reference number of a magnetic field source object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjM(double* M, int* arMesh, int obj); //OC27092018
//EXP int CALL RadObjM(double* M, int obj);

/** Provides coordinates of geometrical center point and magnetic field at that point.
@param B [out] flat array of resulting center point coordinates and field values. If this pointer is 0 at input, the actual output data array can be obtained by calling RadUtiDataGet function. This aray will contain coordinates of geometrical center point(s) and magnetic field components. If obj is a container, the array will include the container members' center points and their magnetic field components. 
@param arMesh [out] flat array defining the structure of resulting field array, i.e. number of dimensions (arMesh[0]) and number of values per dimension. The actual output data array can be obtained by calling RadUtiDataGet function. This aray will contain coordinates of geometrical center point(s) and magnetic field components. If obj is a container, the array will include the container members' center points and their magnetic field components. 
@param obj [in] reference number of a magnetic field source object
@param type [in] character identifying type of a magnetic field characteristic to return (can be 'A' or 'B' or 'H' or 'J' or 'M')
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjCenFld(double* B, int* arMesh, int obj, char type); //OC27092018
//EXP int CALL RadObjCenFld(double* B, int obj, char type);

/** Sets magnetization of the object obj.
@param obj [in] reference number of a magnetic field source object
@param M [in] flat array of 3 components of the magnetization vector
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjSetM(int obj, double* M);

/** Cuts the object obj by a plane passing through a given point normally to a given vector.
@param Objs [out] array of reference numbers of the objects produced by the cutting
@param nobj [out] amount of the objects produced by the cutting
@param obj [in] reference number of the object to cut
@param P [in] array of 3 cartesian coordinates of a point the cutting plane passes through
@param N [in] array of 3 cartesian coordinates of the vector normal to the cutting plane
@param opt [in] pointer to an option string, which can be "Frame->Lab", "Frame->Loc" or 0. This specifies whether the cuting plane is defined in the laboratory frame ("Frame->Lab" or 0) or in the local frame of the object obj. 
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjCutMag(int* Objs, int* nobj, int obj, double* P, double* N, char* opt);

/** Subdivides (segments) the object obj by 3 sets of parallel planes. 
@param n [out] reference number of the object created (as a rule, this is a container object)
@param obj [in] reference number of the object to subdivide
@param SbdPar [in] array of 3 (k1,k2,k3) or 6 (k1,q1,k2,q2,k3,q3) subdivision parameters. The meaning of k1, k2 and k3 depends on the value of the option kxkykz: if kxkykz->Numb (default), then k1, k2 and k3 are subdivision numbers; if kxkykz->Size, they are average sizes of the sub-objects to be produced; q1, q2 and q3 are ratios of the last-to-first sub-object sizes.
@param nSbdPar [in] number of subdivision parameters (length of the SbdPar array)
@param FlatNorm [in] array of 9 numbers specifying cartesian coordinates of 3 vectors normal to the subdivision planes
@param opt [in] pointer to options string, which can be "kxkykz->Numb" (default) or "kxkykz->Size" for the segmentation parameters to be interpreted as numbers of peices or their average dimensions; "Frame->Lab", "Frame->LabTot" or "Frame->Loc" for the subdivision to be performed in the laboratory frame or in local frame of the 3D object. The action of "Frame->Lab" and "Frame->LabTot" differs only for containers: "Frame->Lab" means that each of the objects in the container is subdivided separately; "Frame->LabTot" means that all objects in the container are subdivided as one object, by the same planes. opt can contain composition of these sub-strings separated by ";".
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDivMagPln(int* n, int obj, double* SbdPar, int nSbdPar, double* FlatNorm, char* opt);

/** Subdivides (segments) the object obj by a set of coaxial elliptic cylinders. 
@param n [out] reference number of the object created (as a rule, this is a container object)
@param obj [in] reference number of the object to subdivide
@param SbdPar [in] array of 3 (k1,k2,k3) or 6 (k1,q1,k2,q2,k3,q3) subdivision parameters. The meaning of k1, k2 and k3 depends on the value of the option kxkykz: if kxkykz->Numb (default), then k1, k2 and k3 are subdivision numbers; if kxkykz->Size, they are average sizes of the sub-objects to be produced; q1, q2 and q3 are ratios of the last-to-first sub-object sizes. The parameters (k1,q1),(k2,q2) and (k3,q3) correspond to radial, azimuthal, and axial directions respectively.
@param nSbdPar [in] number of subdivision parameters (length of the SbdPar array)
@param FlatCylPar [in] array of 9 numbers (ax,ay,az,vx,vy,vz,px,py,pz) specifying positions of subdividing coaxial elliptic cylinders in space. The cylinders axis is defined by the point (ax,ay,az) and vector (vx,vy,vz). One of two axes of the cylinder base ellipses is exactly the perpendicular from the point (px,py,pz) to the cylinder axis.
@param rat [in] the ratio of the ellipse axes lengths in the bases of subdividing coaxial elliptic cylinders
@param opt [in] pointer to options string, which can be "kxkykz->Numb" (default) or "kxkykz->Size" for the segmentation parameters to be interpreted as numbers of peices or their average dimensions; "Frame->Lab", "Frame->LabTot" or "Frame->Loc" for the subdivision to be performed in the laboratory frame or in local frame of the 3D object. The action of "Frame->Lab" and "Frame->LabTot" differs only for containers: "Frame->Lab" means that each of the objects in the container is subdivided separately; "Frame->LabTot" means that all objects in the container are subdivided as one object, by the same planes. opt can contain composition of these sub-strings separated by ";".
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDivMagCyl(int* n, int obj, double* SbdPar, int nSbdPar, double* FlatCylPar, double rat, char* opt);

/** Computes geometrical volume of a 3D object.
@param v [out] volume (in mm^3)
@param obj [in] reference number of a 3D object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjGeoVol(double* v, int obj);

/** Computes coordinates of object extrimities in laboratory frame.
@param L [out] array of 6 numbers representing cartesian coordinates of object extrimities (xmin, xmax, ymin, ymax, zmin, zmax)
@param obj [in] reference number of a 3D object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjGeoLim(double* L, int obj);

/** Gives number of degrees of freedom for the relaxation of an object.
@param num [out] number of degrees of freedom
@param obj [in] reference number of a 3D object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDegFre(int* num, int obj);

/** Starts an application for viewing of 3D geometry of the object obj. The viewer is based on the QuickDraw 3D graphics library. 
@param obj [in] reference number of the object to be viewed
@param opt [in] pointer to options string, which can be "Axes->True" (default) or "Axes->False" for showing or not the axes of the Cartesian laboratory frame; "Faces->True" (default) or "Faces->False" for showing or not visible faces of 3D objects; "EdgeLines->True" (default) or "EdgeLines->False" for highlighting or not the edge lines of 3D objects. opt can contain composition of these option sub-strings separated by ";".
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDrwQD3D(int obj, char* opt);

/** Starts an application for viewing of 3D geometry of the object obj. The viewer is based on the GLUT / OpenGL graphics library. 
@param obj [in] reference number of the object to be viewed
@param opt [in] pointer to options string, which can be "Axes->True" (default) or "Axes->False" for showing or not the axes of the Cartesian laboratory frame; "Faces->True" (default) or "Faces->False" for showing or not visible faces of 3D objects; "EdgeLines->True" (default) or "EdgeLines->False" for highlighting or not the edge lines of 3D objects. opt can contain composition of these option sub-strings separated by ";".
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDrwOpenGL(int obj, char* opt);

/** Applies drawing attributes - RGB color (r,g,b) and line thickness thcn - to object obj.
@param obj [in] reference number of the object to which drawing attributes should be applied
@param RGB [in] array of 3 numbers from 0 to 1 specifying intensities of red, green and blue colors
@param thcn [in] line thickness parameter
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDrwAtr(int obj, double* RGB, double thcn);

/** Creates a parallelepiped block with center point {P[0],P[1],P[2]}, dimensions {L[0],L[1],L[2]} and color {RGB[0],RGB[1],RGB[2]}. 
The block is magnetized according to {M[0],M[1],M[2]} then subdivided according to {K[0],K[1],K[2]} and added into the container grp. grp should be defined in advance by calling RadObjCnt().
@param n [out] reference number of the object created
@param P [in] three cartesian coordinates of the block center point
@param L [in] block dimensions
@param M [in] three components of the magnetization vector
@param K [in] array of subdivision parameters
@param nK [in] length of array of subdivision parameters
@param grp [in] reference number of the container object
@param mat [in] reference number of the magnetic material
@param RGB [in] array of 3 numbers from 0 to 1 specifying intensities of red, green and blue colors
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjFullMag(int* n, double* P, double* L, double* M, double* K, int nK, int grp, int mat, double* RGB);

/** Creates a symmetry with respect to plane defined by a point and a normal vector.
@param trf [out] reference number of the symmetry object created
@param P [in] array of 3 numbers representing cartesian coordinates of a point in the plane
@param N [in] array of 3 numbers representing cartesian coordinates of a vector normal to the plane
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfPlSym(int* trf, double* P, double* N);

/** Creates a rotation.
@param trf [out] reference number of the symmetry object created
@param P [in] array of 3 numbers representing cartesian coordinates of a point belonging to the rotation axis
@param V [in] array of 3 numbers representing cartesian coordinates of a vector parallel to the rotation axis
@param phi [in] rotation angle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfRot(int* trf, double* P, double* V, double phi);

/** Creates a translation.
@param trf [out] reference number of the symmetry object created
@param V [in] array of 3 numbers representing cartesian coordinates of the translation vector
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfTrsl(int* trf, double* V);

/** Creates a field inversion.
@param trf [out] reference number of the symmetry object created
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfInv(int* trf);

/** Multiplies original space transformation origtrf by another transformation trf from left.
@param fintrf [out] reference number of the final space transformation
@param origtrf [in] reference number of the original space transformation to be multiplied
@param trf [in] reference number of another space transformation
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfCmbL(int* fintrf, int origtrf, int trf);

/** Multiplies original space transformation origtrf by another transformation trf from right.
@param fintrf [out] reference number of the final space transformation
@param origtrf [in] reference number of the original space transformation to be multiplied
@param trf [in] reference number of another space transformation
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfCmbR(int* fintrf, int origtrf, int trf);

/** Creates mlt-1 symmetry objects of the object obj. Each symmetry object is derived from the previous one by applying the transformation trf to it. Following this, the object obj becomes equivalent to mlt different objects.
@param objout [out] reference number of the final object with symmetries applied
@param obj [in] reference number of the original object to which symmetries should be applied
@param trf [in] reference number of a space transformation
@param mlt [in] multiplicity of the space transformation
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfMlt(int* objout, int obj, int trf, int mlt);

/** Orients object obj by applying transformation trf to it once.
@param objout [out] reference number of the final object with space transformation applied
@param obj [in] reference number of the original object
@param trf [in] reference number of a space transformation
*/
EXP int CALL RadTrfOrnt(int* objout, int obj, int trf);

/** Creates an object mirror with respect to a plane. The object mirror possesses the same geometry as obj, but its magnetization and/or current densities are modified in such a way that the magnetic field produced by the obj and its mirror in the plane of mirroring is perpendicular to this plane.
@param objout [out] reference number of the final object
@param obj [in] an integer number referencing the original object
@param P [in] array of 3 cartesian coordinates of a pointg in the mirror plane
@param N [in] array of 3 cartesian coordinates of a vector normal to the mirror plane
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author P.E., O.C.
*/ 
EXP int CALL RadTrfZerPara(int* objout, int obj, double* P, double* N);

/** Creates an object mirror with respect to a plane. The object mirror possesses the same geometry as obj, but its magnetization and/or current densities are modified in such a way that the magnetic field produced by the obj and its mirror in the plane of mirroring is parallel to this plane.
@param objout [out] reference number of the object
@param obj [in] an integer number referencing the original object
@param P [in] array of 3 cartesian coordinates of a pointg in the mirror plane
@param N [in] array of 3 cartesian coordinates of a vector normal to the mirror plane
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author P.E., O.C.
*/ 
EXP int CALL RadTrfZerPerp(int* objout, int obj, double* P, double* N);

/** Applies material mat to object obj.
@param objout [out] reference number of the final object with material applied
@param obj [in] reference number of the original object
@param mat [in] reference number of the material to be applied
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatApl(int* objout, int obj, int mat);

/** Creates a linear anisotropic magnetic material.
@param mat [out] reference number of the material created
@param Ksi [in] array of 2 magnetic susceptibility values for the directions parallel and perpendicular to the easy magnetization axis
@param Mr [in] array of 3 cartesian coordinates of the remanent magnetization vector
@param nMr [in] number of components of the remanent magnetization vector. If nMr = 1, Mr[0] specifies absolute value of the remanent magnetization; the direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is applied (the magnetization vector is specified at the object creation).
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatLin(int* mat, double* Ksi, double* Mr, int nMr);

/** Creates a pre-defined magnetic material.
The material is identified by its name/formula (e.g. \"NdFeB\"). 
@param mat [out] reference number of the material created
@param id [in] null-terminated string identifying the material
@param m [in] amplitude of the remanent magnetization (for permanent-magnet type materials)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadMatStd(int* mat, char* id, double m);

/** Computes magnetization from magnetic field strength vector.
@param M [out] array of magnetization components calculated
@param nM [out] number of magnetization components calculated (length of the array M)
@param obj [in] reference number of a material or of an object with material applied
@param id [in] string specifying the magnetization components to be calculated (e.g. \"mz\")
@param H [in] magnetic field strength vector in Tesla (mu0*H)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatMvsH(double* M, int* nM, int obj, char* id, double* H);

/** Creates a nonlinear isotropic magnetic material with the magnetization magnitude equal M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3), where H is the magnitude of the magnetic field strength vector (in Tesla).
@param mat [out] reference number of the material created
@param KsiMs1 [in] array of 2 numbers specifying the parameters ms1 and ksi1 (see the formula above)
@param KsiMs2 [in] array of 2 numbers specifying the parameters ms2 and ksi2 (see the formula above)
@param KsiMs3 [in] array of 2 numbers specifying the parameters ms3 and ksi3 (see the formula above)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatIsoFrm(int* mat, double* KsiMs1, double* KsiMs2, double* KsiMs3);

/** Creates a nonlinear isotropic magnetic material with the M versus H curve defined by the list of pairs corresponding values of H and M (H1,M1,H2,M2,...).
@param mat [out] reference number of the material created
@param MatData [in] flat array of material data points (H1,M1,H2,M2,H3,M3,...)
@param np [in] number of material data points
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatIsoTab(int* mat, double* MatData, int np);

/** Creates laminated nonlinear anisotropic magnetic material with packing factor p and the lamination planes perpendicular to the vector N. The magnetization magnitude vs magnetic field strength for the corresponding isotropic material is defined by the formula M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3), where H is the magnitude of the magnetic field strength vector (in Tesla); ksi1, ms1, ksi2, ms2, ksi3, ms3 constants are given by parameters KsiMs1, KsiMs2, KsiMs3.
@param mat [out] reference number of the material created
@param KsiMs1 [in] array of 2 numbers specifying the parameters ms1 and ksi1 (see the formula above)
@param KsiMs2 [in] array of 2 numbers specifying the parameters ms2 and ksi2 (see the formula above)
@param KsiMs3 [in] array of 2 numbers specifying the parameters ms3 and ksi3 (see the formula above)
@param p [in] lamination stacking factor
@param N [in] array of 3 numbers specifying cartesian coordinates of a vector normal to the lamination planes; if the pointer N is 0, the lamination planes are assumed to be perpendicular to the magnetization vector in the object to which the material is applied (the magnetization vector should be specified at the object creation)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatLamFrm(int* mat, double* KsiMs1, double* KsiMs2, double* KsiMs3, double p, double* N);

/** Creates laminated nonlinear anisotropic magnetic material with packing factor p and the lamination planes perpendicular to the vector N. The magnetization magnitude vs magnetic field strength for the corresponding isotropic material is defined by pairs of values H, M in Tesla.
@param mat [out] reference number of the material created
@param MatData [in] flat array of material data points (H1,M1,H2,M2,H3,M3,...)
@param np [in] number of material data points
@param p [in] lamination stacking factor
@param N [in] array of 3 numbers specifying cartesian coordinates of a vector normal to the lamination planes; if the pointer N is 0, the lamination planes are assumed to be perpendicular to the magnetization vector in the object to which the material is applied (the magnetization vector should be specified at the object creation)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatLamTab(int* n, double* MatData, int np, double p, double* N);

/** Creates a nonlinear anisotropic magnetic material. The magnetization vector component parallel to the easy axis is computed either by the formula: ms1*tanh(ksi1*(hpa-hc1)/ms1) + ms2*tanh(ksi2*(hpa-hc2)/ms2) + ms3*tanh(ksi3*(hpa-hc3)/ms3) + ksi0*(hpa-hc0), where hpa is the field strength vector component parallel to the easy axis, or by ksi0*hpa. The magnetization vector component perpendicular to the easy axis is computed either by the formula: ms1*tanh(ksi1*hpe/ms1) + ms2*tanh(ksi2*hpe/ms2) + ms3*tanh(ksi3*hpe/ms3) + ksi0*hpe, where hpe is the field strength vector component perpendicular to the easy axis, or by ksi0*hpe. At least one of the magnetization components should non-linearly depend on the field strength. The direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is later applied.
@param mat [out] reference number of the material created
@param DataPar [in] flat array of constants defining the magnetic material behavior in the direction parallel to the easy magnetization axis. It can be {ksi1,ms1,hc1,ksi2,ms2,hc2,ksi3,ms3,hc3,ksi0,hc0} or {ksi0}. 
@param nDataPar [in] length of the array DataPar. Can be equal to 11 or 1, for the DataPar to be interpreted as {ksi1,ms1,hc1,ksi2,ms2,hc2,ksi3,ms3,hc3,ksi0,hc0} or {ksi0}.
@param DataPer [in] flat array of constants defining the magnetic material behavior in the direction perpendicular to the easy magnetization axis. It can be {ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0} or {ksi0}. 
@param nDataPer [in] length of the array DataPer. Can be equal to 7 or 1, for the DataPer to be interpreted as {ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0} or {ksi0}.
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatAniso(int* mat, double* DataPar, int nDataPar, double* DataPer, int nDataPer);

/** Builds interaction matrix for the object obj.
@param n [out] reference number of the interaction matrix created
@param obj [in] reference number of the object for which the interaction matrix should be created
@param srcobj [in] reference number of the object creating additional constant magnetic field 
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadRlxPre(int* n, int obj, int srcobj);

/** Executes manual relaxation procedure for the interaction matrix intrc.
@param D [out] an array of four numbers specifying: [0] average absolute change in magnetization after previous iteration over all the objects participating in the relaxation, [1] maximum absolute value of magnetization over all the objects participating in the relaxation, [2] maximum absolute value of magnetic field strength over central points of all the objects participating in the relaxation, and [3] actual number of iterations done. The values [0]-[2] are those of last iteration.
@param n [out] length of array D
@param intrc [in] an integer number referencing the interaction matrix
@param meth [in] an integer number specifying the method of relaxation to be used
@param iter [in] an integer number specifying number of iterations to be made
@param rlxpar [in] a floating point number between 0 and 1 specifying the relaxation parameter
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadRlxMan(double* D, int* n, int intrc, int meth, int iter, double rlxpar);

/** Executes automatic relaxation procedure for the interaction matrix intrc.
The relaxation stops whenever the change of magnetization (averaged over all sub-elements) between two successive iterations is smaller than prec or the number of iterations is larger than iter.
@param D [out] an array of four numbers specifying: [0] average absolute change in magnetization after previous iteration over all the objects participating in the relaxation, [1] maximum absolute value of magnetization over all the objects participating in the relaxation, [2] maximum absolute value of magnetic field strength over central points of all the objects participating in the relaxation, and [3] actual number of iterations done. The values [0]-[2] are those of last iteration.
@param n [out] length of array D
@param intrc [in] an integer number referencing the interaction matrix
@param prec [in] a real number specifying an absolute precision value for magnetization (in Tesla), to be reached by the end of the relaxation
@param iter [in] maximum number of iterations permitted to reach the specified precision
@param meth [in] an integer number specifying the method of relaxation to be used (values 0, 3 - 5 can be used; 0 means default method)
@param opt [in] pointer to an option string, which can be "ResetM->True" (default) or "ResetM->False"
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadRlxAuto(double* D, int* n, int intrc, double prec, int iter, int meth, const char* opt);

/** Updates external field data for the relaxation (to take into account e.g. modification of currents in coils, if any) without rebuilding the interaction matrix.
@param intrc [in] an integer number referencing the interaction object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadRlxUpdSrc(int intrc);

/** Builds an interaction matrix and performs a relaxation procedure. 
The relaxation stops whenever the change of magnetization (averaged over all sub-elements) between two successive iterations is smaller than prec or the number of iterations is larger than iter. The interaction matrix is deleted. 
@param D [out] an array of four numbers specifying: [0] average absolute change in magnetization after previous iteration over all the objects participating in the relaxation, [1] maximum absolute value of magnetization over all the objects participating in the relaxation, [2] maximum absolute value of magnetic field strength over central points of all the objects participating in the relaxation, and [3] actual number of iterations done. The values [0]-[2] are those of last iteration.
@param n [out] length of array D
@param obj [in] an integer number specifying the object to solve for magnetization
@param prec [in] a real number specifying an absolute precision value for magnetization (in Tesla), to be reached by the end of the relaxation
@param iter [in] maximum number of iterations permitted to reach the specified precision
@param meth [in] an integer number specifying the method of relaxation to be used (values 0, 3 - 5 can be used; 0 means default method)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author P.E., O.C.
*/ 
EXP int CALL RadSolve(double* D, int* n, int obj, double prec, int iter, int meth);

/** Computes magnetic field created by the object obj at one or many points.
@param B [out] flat array of all computed values of the magnetic field components (should be allocated by calling function)
@param nB [out] total number of calculated magnetic field component values
@param obj [in] reference number of the magnetic field source object
@param id [in] string identifying magnetic field components to be computed
@param Coords [in] flat array of coordinates of all points where the field should be computed (x1,y1,z1,x2,y2,z2,...)
@param np [in] number of points where the magnetic field should be calculated (the length of the array Coords is equal to 3*np)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFld(double* B, int* nB, int obj, char* id, double* Coords, int np);

/** Sets general absolute accuracy levels for computation of magnetic field induction (prcB), vector potential (prcA), induction integrals along straight line (prcBint), field force (prcFrc), relativistic particle trajectory coordinates (prcTrjCrd) and angles (prcTrjAng).
@param n [out] dummy
@param prcB [in] absolute accuracy level for magnetic field induction [T]
@param prcA [in] absolute accuracy level for vector potential
@param prcBInt [in] absolute accuracy level for magnetic field induction integrals along straight line [T*mm]
@param prcFrc [in] absolute accuracy level for the force [N]
@param prcTrjCrd [in] absolute accuracy level for particle trajectory coordinate
@param prcTrjAng [in] absolute accuracy level for particle trajectory angle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldCmpCrt(int* n, double prcB, double prcA, double prcBInt, double prcFrc, double prcTrjCrd, double prcTrjAng);

/** Sets general absolute accuracy levels for computation of magnetic field induction (PrcB), vector potential (PrcA), induction integral along straight line (PrcBInt), field force (PrcForce), torque (PrcTorque), energy (PrcEnergy); relativistic charged particle trajectory coordinates (PrcCoord) and angles (PrcAngle). The function works according to the mechanism of string options. The name(s) of the option(s) should be: PrcB, PrcA, PrcBInt, PrcForce, PrcTorque, PrcEnergy, PrcCoord, PrcAngle.
@param n [out] dummy
@param opt [in] pointer to an option string, which can be "PrcB->..." or "PrcA->..." or "PrcBInt->..." or "PrcForce->..." or "PrcTorque->..." or "PrcEnergy->..." or "PrcCoord->..." or "PrcAngle->...", where "..." should be replaced by the appropriate absolute accuracy level.
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldCmpPrc(int* n, char* opt);

/** Shows the physical units currently in use.
@param OutStr [out] string of physical units currently in use
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldUnits(char* OutStr);

/** Specifies the length of string about physical units currently in use.
@param size [out] length of string about physical units currently in use
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldUnitsSize(int* size);

/** Switches on or off the randomization of all the length values. The randomization magnitude can be set by the function radFldLenTol.
@param n [out] dummy
@param OnOrOff [in] string containing either "on" or "off"
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldLenRndSw(int* n, char* OnOrOff);

/** Sets absolute and relative randomization magnitudes for all the length values, including coordinates and dimensions of the objects producing magnetic field, and coordinates of points where the field is computed. Optimal values of the variables can be: RelVal=10^(-11), AbsVal=L*RelVal, ZeroVal=AbsVal, where L is the distance scale value (in mm) for the problem to be solved. Too small randomization magnitudes can result in run-time code errors.
@param n [out] dummy
@param AbsVal [in] absolute position/length randomization magnitude [mm]
@param RelVal [in] relative position/length randomization magnitude
@param ZeroVal [in] absolute zero tolerance [mm]
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldLenTol(int* n, double AbsVal, double RelVal, double ZeroVal);

/** Computes potential energy (in Joule) of the object objdst in the field created by the object objsrc. If SbdPar = 0, the function performes the computation based on absolute accuracy value for the energy (by default 10 Joule; can be modified by the function radFldCmpPrc). Otherwise, the computation is performed based on the destination object subdivision numbers (kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]).
@param d [out] computed energy value
@param objdst [in] reference number of the object, the energy of which should be computed
@param objsrc [in] reference number of the source object creating the magnetic field
@param SbdPar [in] array of 3 integer numbers specifying subdivision parameters for the energy computation (kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]). If SbdPar = 0, the function performes the computation based on absolute accuracy value for the energy (by default 10 Joule, can be modified by the function radFldCmpPrc). 
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldEnr(double* d, int objdst, int objsrc, int* SbdPar);

/** Computes force (in Newton) acting on the object objdst in the field produced by the object objsrc. If SbdPar = 0, the function performes the computation based on absolute accuracy value for the force (by default 10 Newton; can be modified by the function radFldCmpPrc). Otherwise, the computation is performed based on the destination object subdivision numbers {kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]}.
@param f [out] computed force component(s)
@param nf [out] number of force components computed
@param objdst [in] reference number of the object, on which the force is acting
@param objsrc [in] reference number of the source object creating the magnetic field
@param id [in] string identifying the force components to be computed (fx|fy|fz)
@param SbdPar [in] array of 3 integer numbers specifying subdivision parameters for the energy/force computation (kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]). If SbdPar = 0, the function performes the computation based on absolute accuracy value for the force (by default 10 Newton, can be modified by the function radFldCmpPrc). 
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldEnrFrc(double* f, int* nf, int objdst, int objsrc, char* id, int* SbdPar);

/** Computes torque (in Newton*mm) with respect to point P, acting on the object objdst in the field produced by the object objsrc. If SbdPar = 0, the function performes the computation based on absolute accuracy value for the torque (by default 10 Newton*mm; can be modified by the function radFldCmpPrc). Otherwise, the computation is performed based on the destination object subdivision numbers {kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]}.
@param f [out] computed torque component(s)
@param nf [out] number of torque components computed
@param objdst [in] reference number of the object on which the force is acting
@param objsrc [in] reference number of the source object creating the magnetic field
@param id [in] string identifying the force components to be computed (tx|ty|tz)
@param P [in] array of 3 real numbers specifying cartesian coordinates of the point with respect to which the torque should be computed  
@param SbdPar [in] array of 3 integer numbers specifying subdivision parameters for the energy/force computation (kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]). If SbdPar = 0, the function performes the computation based on absolute accuracy value for the torque (by default 10 Newton*mm, can be modified by the function radFldCmpPrc). 
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldEnrTrq(double* f, int* nf, int objdst, int objsrc, char* id, double* P, int* SbdPar);

/** Computes force of the field produced by the object obj into a shape defined by shape. shape can be the result of RadObjRecMag (parallelepiped) or RadFldFrcShpRtg (rectangular surface). This function uses the algorithm based on Maxwell tensor, which may not always provide high efficiency. We suggest to use the function RadFldEnrFrc instead of this function.
@param f [out] computed force component(s)
@param nf [out] number of force components computed
@param obj [in] reference number of the magnetic field source object
@param shape [in] reference number of the shape object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldFrc(double* f, int* nf, int obj, int shape);

/** Computes "focusing potential" for trajectory of relativistic charged particle in magnetic field produced by the object obj. The integration is made from P1 to P2 with np equidistant points.
@param d [out] computed focusing potential value
@param obj [in] reference number of the magnetic field source object
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of an edge point of the integration segment
@param P2 [in] array of 3 real numbers specifying cartesian coordinates of an edge point of the integration segment
@param np [in] number of points for the integration
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldFocPot(double* d, int obj, double* P1, double* P2, int np);

/** Computes magnetic field integral produced by the object obj along a straight line specified by points P1 and P2; depending on the InfOrFin variable value, the integral is infinite ("inf") or finite ("fin"), from P1 to P2; the field integral component is specified by the id input variable. The unit is T*mm.
@param f [out] computed field integral component(s)
@param nf [out] number of field integral components computed
@param obj [in] reference number of the object creating the magnetic field
@param InfOrFin [in] string specifying the type of field integral: finite ("fin") or infinite ("inf")
@param id [in] string identifying the field integral components to be computed (ibx|iby|ibz)
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of a point on the integration line
@param P2 [in] array of 3 real numbers specifying cartesian coordinates of another point on the integration line
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldInt(double* f, int* nf, int obj, char* InfOrFin, char* id, double* P1, double* P2);

/** Computes magnetic field created by object obj in np equidistant points along a line segment from P1 to P2; the field component is specified by the id input variable.
@param B [out] computed field value(s)
@param nB [out] number of field values computed
@param obj [in] reference number of the object creating the magnetic field
@param id [in] string identifying magnetic field components to be computed
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of an edge point of the line segment
@param P2 [in] array of 3 real numbers specifying cartesian coordinates of another edge point of the line segment
@param np [in] number of points where the magnetic field should be calculated
@param ArgOrNoArg [in] string specifying whether or not to output a longitudinal position for each point where the field is computed ("arg|noarg")
@param start [in] start value for the longitudinal position
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldLst(double* B, int* nB, int obj, char* id, double* P1, double* P2, int np, char* ArgOrNoArg, double start);

/** Computes transverse coordinates and their derivatives (angles) of a relativistic charged trajectory in the magnetic field produced by object obj, using a Runge-Kutta integration. The charge of the particle is that of electron. All positions are in millimeters and angles in radians.
@param f [out] computed tragectory parameters
@param nf [out] number of tragectory parameters computed
@param obj [in] reference number of the magnetic field source object
@param E [in] particle energy [GeV]
@param InitCond [in] array of 4 real numbers specifying initial transverse coordinates and angles of the trajectory (x0,dxdy0,z0,dzdy0, y is longitudinal coordinate)
@param LongLim [in] array of 2 real numbers specifying limits of the longitudinal position (y1,y2)
@param np [in] number of points where the trajectory should be calculated
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldPtcTrj(double* f, int* nf, int obj, double E, double* InitCond, double* LongLim, int np);

/** Computes matrices of 2nd order kicks for trajectory of relativistic charged particle in periodic magnetic field produced by the object obj. The computed kick matrices can be used in charged particle tracking codes.  The longitudinal integration along one period starts at point P1 and is done along direction pointed by vector Ns; one direction of the transverse grid is pointed by vector Ntr, the other transverse direction is given by vector product of Ntr and Ns.
@param M1 [out] flat array of real numbers representing the matrix of kick values in the first transverse direction given by the Ntr vector
@param M2 [out] flat array of real numbers representing the matrix of kick values in the second transverse direction given by vector product of Ntr and Ns
@param IntBtrE2 [out] flat array of real numbers representing the matrix of longitudinally-integrated squared transverse magnetic field
@param Arg1 [out] array of position values in the first transverse direction where kick matrix was computed
@param Arg2 [out] array of position values in the second transverse direction where kick matrix was computed
@param size [out] length of formatted string containing the computed results; such string can be prepared by the function RadFldFocKickPerFormStr
@param obj [in] reference number of the object creating the magnetic field
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of the integration start point
@param Ns [in] array of 3 real numbers specifying cartesian coordinates of a vector defining the longitudinal direction for the periodic magnetic field
@param per [in] period length
@param nper [in] number of full periods
@param nps [in] number of longitudinal points
@param Ntr [in] array of 3 real numbers specifying cartesian coordinates of a vector defining first transverse direction
@param r1 [in] range of the transverse grid along direction given by Ntr vector
@param np1 [in] number of points in the transverse grid along direction given by Ntr vector
@param d1 [in] step to calculate derivative along direction given by Ntr vector; d1=0 means that the step for derivative calculation is equal to the corresponding step of the transverse grid where kick matrix is computed
@param r2 [in] range of the transverse grid along direction given by vector product of Ntr and Ns
@param np2 [in] number of points in the transverse grid along direction given by vector product of Ntr and Ns
@param d2 [in] step to calculate derivative along direction given by vector product of Ntr and Ns; d2=0 means that the step for derivative calculation is equal to the corresponding step of the transverse grid where kick matrix is computed
@param nh [in] maximum number of magnetic field harmonics to take into account
@param com [in] arbitrary string comment
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@authors: P.E., O.C.
*/
EXP int CALL RadFldFocKickPer(double* M1, double* M2, double* IntBtrE2, double* Arg1, double* Arg2, int* size, int obj, double* P1, double* Ns, double per, int nper, int nps, double* Ntr, double r1, int np1, double d1, double r2, int np2, double d2, int nh, char* com);

/** Prepares a formatted string containing second-order kick matrices for trajectory of relativistic charged particle in periodic magnetic field.
@param str [out] C-string containing second-order kick matrices values in the format compatible with particle tracking computer codes (e.g. BETA, TRACY)
@param M1 [in] flat array of real numbers representing the matrix of kick values in one transverse direction
@param M2 [in] flat array of real numbers representing the matrix of kick values in another transverse direction
@param IntBtrE2 [in] flat array of real numbers representing the matrix of longitudinally-integrated squared transverse magnetic field
@param Arg1 [in] array of position values in the first transverse direction where kick matrix was computed
@param Arg2 [in] array of position values in the second transverse direction where kick matrix was computed
@param np1 [in] number of points in the transverse grid along the first transverse direction
@param np2 [in] number of points in the transverse grid along the second transverse direction
@param per [in] period length
@param nper [in] number of full periods
@param com [in] arbitrary string comment that will be included into the formatted string str
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@authors: O.C.
*/
EXP int CALL RadFldFocKickPerFormStr(char* str, double* M1, double* M2, double* IntBtrE2, double* Arg1, double* Arg2, int np1, int np2, double per, int nper, char* com);

/** Computes a virtual "shim signature", i.e. variation of a given magnetic field component introduced by given displacement of magnetic field source object.
@param f [out] computed array of field component values
@param nf [out] number of values in the field component array computed
@param obj [in] reference number of the magnetic field source object
@param id [in] string identifying magnetic field component to be computed (can be e.g. "bx", "bz", "ix", "iz",...)
@param V [in] array of 3 real numbers specifying cartesian coordinates of the field source object displacement
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of an edge point of the line segment along which the "shim signature" should be computed
@param P2 [in] array of 3 real numbers specifying cartesian coordinates of another edge point of the line segment
@param np [in] number of points where the magnetic field should be computed
@param Vi [in] array of 3 real numbers specifying cartesian coordinates of a vector defining the integration line (is taken into account only if id string specifies field integral component)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@authors: O.C.
*/
EXP int CALL RadFldShimSig(double* f, int* nf, int obj, char* id, double* V, double* P1, double* P2, int np, double* Vi);

/** Creates a rectangle with central point P and dimensions W (to be used for force computation via Maxwell tensor).
@param n [out] reference number of the object created
@param P [in] array of 3 real numbers specifying cartesian coordinates of the center point
@param W [in] array of 2 real numbers specifying dimensions of the rectangle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldFrcShpRtg(int* n, double* P, double* W);

/** Deletes object obj.
@param n [out] dummy
@param obj [in] reference number of the object to be deleted
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadUtiDel(int* n, int obj);

/** Deletes all previously created objects.
@param n [out] dummy
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadUtiDelAll(int* n);

/** Outputs information about object obj.
@param OutStr [out] string containing information about obj
@param arObj [in] array of object reference numbers to show information for
@param nObj [in] length of array of object reference numbers
@param AscOrBin [in] string specifying format of the output information string, can be \"asc\" (for ASCII) or \"bin\" (for Binary)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
//EXP int CALL RadUtiDmp(char* OutStr, int obj);
//EXP int CALL RadUtiDmp(char* OutStr, int* arObj, int nObj, char* AscOrBin);
EXP int CALL RadUtiDmp(char* OutStr, int* pSize, int* arObj, int nObj, char* AscOrBin); //OC01102018

/** Outputs information about object obj after reading it from internal buffer.
@param OutStr [out] string containing information about obj
@param AscOrBin [in] string specifying format of the output information string, can be \"asc\" (for ASCII) or \"bin\" (for Binary)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
//EXP int CALL RadUtiDmpRead(char* OutStr, char* AscOrBin);

/** Gives necessary length of the string to include dump information about object obj.
@param size [out] size of string containing information about obj
@param arObj [in] array of object reference numbers to show information for
@param nObj [in] length of array of object reference numbers
@param AscOrBin [in] string specifying format of the output information string, can be \"asc\" (for ASCII) or \"bin\" (for Binary)
@param doEraseBuf [in] switch specifying whether the output buffer has to be erased or not (leaving the buffer unerased allows not to repeat the dump operation at a subsequent call of RadUtiDmp)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
//EXP int CALL RadUtiDmpSize(int* size, int* arElem, int nElem, char* AscOrBin, bool doEraseBuf);
////EXP int CALL RadUtiDmpSize(int* size, int Elem);


EXP int CALL RadUtiDmpPrs(int* arElem, int* nElem, unsigned char* sBytes, int nBytes); //OC01102018

/** Sets interruption time quanta in seconds for platforms with no preemptive multitasking.
@param t [in] interruption time quanta [s]
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author P.E., O.C.
*/
EXP int CALL RadUtiIntrptTim(double* d, double t);

/** Returns data resulting from previous calculations in cases when the data size was not known 'a priori', e.g. after executing functions RadObjM, RadObjCenFld, RadUtiDmp,...
@param size [out] pointer to the resulting data (to be allocated in calling function)
@param typeData [in] string identifying type of the data: \"mad\" for multi-dim. array of double, \"mai\" for multi-dim. array of integer, \"bin\" for byte array, \"asc\" for ASCII string, \"d\" for double, \"i\" for integer
@param key [in] additional identifier of the data to be extracted, e.g. to ensure thread safety (not implemented yet)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadUtiDataGet(char* pcData, const char typeData[3], long key=0); //OC04102018
//EXP int CALL RadUtiDataGet(char* pcData, char typeData[3], long key=0); //OC27092018
//EXP int CALL RadUtiDataGet(double* pData, long key); //OC15092018

/** Identifies the version number of the Radia DLL.
@param d [out] version number
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author P.E., O.C.
*/ 
EXP int CALL RadUtiVer(double* d);

EXP int CALL RadUtiYeldFuncSet(int (*pExtFunc)());

#ifdef __cplusplus  
}
#endif
