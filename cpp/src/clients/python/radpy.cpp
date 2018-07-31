/************************************************************************//**
 * File: radpy.cpp
 * Description: Python binding
 * Project: Radia
 * First release: June 2018
 *
 * @authors O.Chubar (BNL), J.Edelen (RadiaSoft)
 * @version 0.02
 ***************************************************************************/

#include "radentry.h"
#include "pyparse.h"
#include "auxparse.h"

/************************************************************************//**
 * Error messages related to Python interface functions
 ***************************************************************************/
static const char strEr_BadFuncArg[] = "Incorrect function arguments";

///************************************************************************//**
// * Global objects to be used across different function calls
// ***************************************************************************/
static char g_strErTot[2000]; 
CPyParse g_pyParse(RadErrGetText);

/************************************************************************//**
 * Auxiliary function for combining error string out of pieces
 ***************************************************************************/
static char* CombErStr(const char* s1, const char* s2)
{
	return strcat(strcpy(g_strErTot, s1), s2);
}

/************************************************************************//**
 * Auxiliary function for parsing Orientation definition string
 ***************************************************************************/
static char ParseOrnt(PyObject* oOrnt, char aDef='x', const char* sFuncName=0)
{
	if(oOrnt == 0) return aDef;

	char sOrnt[256]; *sOrnt = '\0';
	CPyParse::CopyPyStringToC(oOrnt, sOrnt, 256);
	char a = *sOrnt;
	if((a != 'x') && (a != 'X') && (a != 'y') && (a != 'Y') && (a != 'z') && (a != 'Z')) 
	{
		const char sErCom[] = "orientation definition should be \'x\', \'y\' or \'z\'";
		char sAux[1024];
		strcpy(sAux, ": ");
		if(sFuncName == 0) 
		{
			strcat(sAux, sErCom);
		}
		else 
		{
			strcat(sAux, sFuncName);
			strcat(sAux, ", ");
			strcat(sAux, sErCom);
		}
		throw CombErStr(strEr_BadFuncArg, sAux);
	}
	return a;
}

/************************************************************************//**
 * Auxiliary function for parsing Magnetization vector (arM is expected to be allocated in calling function)
 ***************************************************************************/
static void ParseM(double arM[3], PyObject* oM, const char* sFuncName=0)
{
	if(oM != 0)
	{
		int lenM = 3;
		bool lenIsSmall = false;
		CPyParse::CopyPyListElemsToNumArray(oM, 'd', arM, lenM, lenIsSmall);
		if((lenM != 3) || lenIsSmall) 
		{
			const char sErCom[] = "incorrect definition of magnetization vector";
			char sAux[1024];
			strcpy(sAux, ": ");
			if(sFuncName == 0) 
			{
				strcat(sAux, sErCom);
			}
			else 
			{
				strcat(sAux, sFuncName);
				strcat(sAux, ", ");
				strcat(sAux, sErCom);
			}
			throw CombErStr(strEr_BadFuncArg, sAux);
		}
	}
	else
	{
		arM[0] = 0.; arM[1] = 0.; arM[2] = 0.;
	}
}

/************************************************************************//**
 * Magnetic Field Source: Rectangular Parallelepiped with Constant Magnetizatiom over volume
 ***************************************************************************/
static PyObject* radia_ObjRecMag(PyObject *self, PyObject *args)
{//The parallelepiped block is defined through its center point P[3], dimensions L[3], and magnetization M[3].

	PyObject *oP=0, *oL=0, *oM=0, *oResInd=0;

	try
	{
		if(!PyArg_ParseTuple(args, "OO|O:ObjRecMag", &oP, &oL, &oM)) throw CombErStr(strEr_BadFuncArg, ": ObjRecMag");
		if((oP == 0) || (oL == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjRecMag");

		double arP[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjRecMag, incorrect definition of center point"));

		double arL[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oL, 'd', arL, 3, CombErStr(strEr_BadFuncArg, ": ObjRecMag, incorrect definition of dimensions"));

		double arM[] = {0.,0.,0.};
		ParseM(arM, oM, "ObjRecMag");

		int ind = 0;
		g_pyParse.ProcRes(RadObjRecMag(&ind, arP, arL, arM));

		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
* Creates a uniformly magnetized extruded polygon.
***************************************************************************/
static PyObject* radia_ObjThckPgn(PyObject *self, PyObject *args)
{
	PyObject *oPgn=0, *oOrnt=0, *oM=0, *oResInd=0;
	double *arCrd=0;

	try
	{
		double xc=0, lx=0;

		if(!PyArg_ParseTuple(args, "ddO|OO:ObjThckPgn", &xc, &lx, &oPgn, &oOrnt, &oM)) throw CombErStr(strEr_BadFuncArg, ": ObjThckPgn");
		if((lx <= 0) || (oPgn == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjThckPgn");

		int nCrd=0;
		if(!(CPyParse::CopyPyNestedListElemsToNumAr(oPgn, 'd', arCrd, nCrd))) throw CombErStr(strEr_BadFuncArg, ": ObjThckPgn, incorrect polygon definition");
		int nv = nCrd >> 1;
		if(nCrd != (nv << 1)) throw CombErStr(strEr_BadFuncArg, ": ObjThckPgn, incorrect polygon definition (total number of coordinates should be even)");

		char a = ParseOrnt(oOrnt, 'x', "ObjThckPgn");

		double arM[] = {0.,0.,0.};
		ParseM(arM, oM, "ObjThckPgn");

		int ind = 0;
		g_pyParse.ProcRes(RadObjThckPgn(&ind, xc, lx, arCrd, nv, a, arM));

 		oResInd = Py_BuildValue("i", ind);
 		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arCrd != 0) delete[] arCrd;
	return oResInd;
}

/************************************************************************//**
* Creates a uniformly magnetized polyhedron (closed volume limited by planes).
***************************************************************************/
static PyObject* radia_ObjPolyhdr(PyObject *self, PyObject *args)
{
	PyObject *oVerts=0, *oFaces=0, *oM=0, *odMdr=0, *oJ=0, *odJdr=0, *oRelAbs=0, *oResInd=0;
	double *arCrd=0;
	int *arFaceInds=0, *arFaceLens=0;
	double *arM=0, *ar_dMdr=0, *arJ=0, *ar_dJdr=0;

	try
	{
		if(!PyArg_ParseTuple(args, "OO|OOOOO:ObjPolyhdr", &oVerts, &oFaces, &oM, &odMdr, &oJ, &odJdr, &oRelAbs)) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr");
		if((oVerts == 0) || (oFaces == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr");
		if((!PyList_Check(oVerts)) || (!PyList_Check(oFaces))) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr");

		int nVerts = (int)PyList_Size(oVerts);
		if(nVerts < 4) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, number of vertices should be >= 4");

		int nFaces=0;
		CPyParse::FindLengthsOfElemListsOrArrays(oFaces, arFaceLens, nFaces);
		if(nFaces < 4) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, number of faces should be >= 4");

		int nCrd=0;
		if(!(CPyParse::CopyPyNestedListElemsToNumAr(oVerts, 'd', arCrd, nCrd))) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of vertices");
		if(nCrd != nVerts*3) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of vertices (each vertex point should be defined by 3 cartesian coordinates)");

		int nFaceInds=0;
		if(!(CPyParse::CopyPyNestedListElemsToNumAr(oFaces, 'i', arFaceInds, nFaceInds))) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of faces");

		int len3 = 3, len9 = 9;
		bool lenIsSmall = false;
		if(oM != 0)
		{
			CPyParse::CopyPyListElemsToNumArray(oM, 'd', arM, len3, lenIsSmall);
			if((arM == 0) || (len3 != 3) || lenIsSmall) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of magnetization vector");
		}
		if(odMdr != 0)
		{
			if(!(CPyParse::CopyPyNestedListElemsToNumAr(odMdr, 'd', ar_dMdr, len9))) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of magnetization linear coefficients");
			if((ar_dMdr == 0) || (len9 != 9)) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of magnetization linear coefficients");
		}
		if(oJ != 0)
		{
			CPyParse::CopyPyListElemsToNumArray(oJ, 'd', arJ, len3, lenIsSmall);
			if((arJ == 0) || (len3 != 3) || lenIsSmall) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of magnetization vector");
		}
		if(odJdr != 0)
		{
			if(!(CPyParse::CopyPyNestedListElemsToNumAr(odJdr, 'd', ar_dJdr, len9))) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of current density linear coefficients");
			if((ar_dJdr == 0) || (len9 != 9)) throw CombErStr(strEr_BadFuncArg, ": ObjPolyhdr, incorrect definition of current density linear coefficients");
		}

		char sRelAbs[1024]; 
		strcpy(sRelAbs, "rel");
		if(oRelAbs != 0)
		{
			CPyParse::CopyPyStringToC(oRelAbs, sRelAbs, 1024);
		}

		int ind = 0;
		g_pyParse.ProcRes(RadObjPolyhdr(&ind, arCrd, nVerts, arFaceInds, arFaceLens, nFaces, arM, ar_dMdr, arJ, ar_dJdr));
		//to add sRelAbs to the end of RadObjPolyhdr!

 		oResInd = Py_BuildValue("i", ind);
 		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arCrd != 0) delete[] arCrd;
	if(arFaceInds != 0) delete[] arFaceInds;
	if(arFaceLens != 0) delete[] arFaceLens;
	if(arM != 0) delete[] arM;
	if(ar_dMdr != 0) delete[] ar_dMdr;
	if(arJ != 0) delete[] arJ;
	if(ar_dJdr != 0) delete[] ar_dJdr;
	return oResInd;
}

/************************************************************************//**
* Creates a uniformly magnetized finite-length arc of polygonal cross-section.
***************************************************************************/
static PyObject* radia_ObjArcPgnMag(PyObject *self, PyObject *args)
{
	PyObject *oP=0, *oOrnt=0, *oPgn=0, *oPhi=0, *oSymNo=0, *oM=0, *oResInd=0;
	double *arCrd=0;
	try
	{	
		int nseg = 0;
		if(!PyArg_ParseTuple(args, "OOOOi|OO:ObjArcPgnMag", &oP, &oOrnt, &oPgn, &oPhi, &nseg, &oSymNo, &oM)) throw CombErStr(strEr_BadFuncArg, ": ObjArcPgnMag");
		if((oP == 0) || (oPgn == 0) || (oPhi == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjArcPgnMag");
		if(nseg <= 0) throw CombErStr(strEr_BadFuncArg, ": ObjArcPgnMag, number of segments should be positive");

		char a = ParseOrnt(oOrnt, 'x', "ObjArcPgnMag");

		double arP[2], arPhi[2];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 2, CombErStr(strEr_BadFuncArg, ": ObjArcPgnMag, incorrect definition of initial center point"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oPhi, 'd', arPhi, 2, CombErStr(strEr_BadFuncArg, ": ObjArcPgnMag, incorrect definition of initial and final angles"));

		char sSymNo[1024]; *sSymNo = '\0';
		if(oSymNo != 0) CPyParse::CopyPyStringToC(oSymNo, sSymNo, 1024);

		int nCrd=0;
		if(!(CPyParse::CopyPyNestedListElemsToNumAr(oPgn, 'd', arCrd, nCrd))) throw CombErStr(strEr_BadFuncArg, ": ObjArcPgnMag, incorrect definition of cross-section polygon");
		int nv = nCrd >> 1;
		if(nCrd != (nv << 1)) throw CombErStr(strEr_BadFuncArg, ": ObjArcPgnMag, incorrect definition of cross-section polygon (each vertex point should be defined by 2 cartesian coordinates)");

		double arM[] = {0.,0.,0.};
		ParseM(arM, oM, "ObjArcPgnMag");

		int ind = 0;
		g_pyParse.ProcRes(RadObjArcPgnMag(&ind, arP, a, arCrd, nv, arPhi, nseg, sSymNo[0], arM));

		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arCrd != 0) delete[] arCrd;
	return oResInd;
}

/************************************************************************//**
* Attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on slices.
***************************************************************************/
static PyObject* radia_ObjMltExtPgn(PyObject *self, PyObject *args)
{
	PyObject *oSlices=0, *oM=0, *oResInd=0;
	int *arSliceLens=0;
	double *arAlt=0, *arCrd=0;

	try
	{
		if(!PyArg_ParseTuple(args, "O|O:ObjMltExtPgn", &oSlices, &oM)) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtPgn");
		if(oSlices == 0) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtPgn");

		//Parsing slices structure
		if(!PyList_Check(oSlices)) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtPgn");

		int ns = (int)PyList_Size(oSlices);
		if(ns < 2) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtPgn, number of slices / layers should be >= 2");

		arSliceLens = new int[ns];
		arAlt = new double[ns];

		vector<double> vCrd;
		int prevSize_vCrd=0, size_vCrd;
		const char strErCom[] = ": ObjMltExtPgn, incorrect slices / layers definition";
		for(int i=0; i<ns; i++)
		{
			PyObject *o = PyList_GetItem(oSlices, (Py_ssize_t)i);
			if(o == 0) throw CombErStr(strEr_BadFuncArg, strErCom);
			if(!PyList_Check(o)) throw CombErStr(strEr_BadFuncArg, strErCom);

			PyObject *oPgn = PyList_GetItem(o, 0);
			if(oPgn == 0) throw CombErStr(strEr_BadFuncArg, strErCom);
		
			if(!(CPyParse::CopyPyNestedListElemsToNumVect(oPgn, 'd', &vCrd))) throw CombErStr(strEr_BadFuncArg, strErCom);
			size_vCrd =  (int)vCrd.size();
			arSliceLens[i] = (size_vCrd - prevSize_vCrd) >> 1;
			prevSize_vCrd = size_vCrd;
		
			PyObject *oAlt = PyList_GetItem(o, 1);
			if(oAlt == 0) throw CombErStr(strEr_BadFuncArg, strErCom);

			if(!PyNumber_Check(oAlt)) throw CombErStr(strEr_BadFuncArg, strErCom);
			arAlt[i] = PyFloat_AsDouble(oAlt);
		}

		CAuxParse::DoubleVect2Arr(vCrd, arCrd);

		double arM[] = {0.,0.,0.};
		ParseM(arM, oM, "ObjMltExtPgn");

		int ind = 0;
		g_pyParse.ProcRes(RadObjMltExtPgn(&ind, arCrd, arSliceLens, arAlt, ns, arM));

		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}

	if(arSliceLens != 0) delete[] arSliceLens;
	if(arAlt != 0) delete[] arAlt;
	if(arCrd != 0) delete[] arCrd;
	return oResInd;
}

/************************************************************************//**
* Attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on rectangular slices
***************************************************************************/
static PyObject* radia_ObjMltExtRtg(PyObject *self, PyObject *args)
{
	PyObject *oSlices=0, *oM=0, *oResInd=0;
	double *arDims=0, *arCrd=0;

	try
	{
		if(!PyArg_ParseTuple(args, "O|O:ObjMltExtRtg", &oSlices, &oM)) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtRtg");
		if(oSlices == 0) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtRtg");

		//Parsing slices structure
		if(!PyList_Check(oSlices)) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtPgn");

		int ns = (int)PyList_Size(oSlices);
		if(ns < 2) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtRtg, number of slices / layers should be >= 2");

		const char strErCom[] = ": ObjMltExtRtg, incorrect slices / layers definition";
		const char strErComCenPt[] = ": ObjMltExtRtg, incorrect slices / layers definition (rentangle center point)";
		const char strErComDims[] = ": ObjMltExtRtg, incorrect slices / layers definition (rentangle dimensons)";
		vector<double> vCrd, vDims;

		double arAuxP[3], arAuxD[2];
		for(int i=0; i<ns; i++)
		{
			PyObject *o = PyList_GetItem(oSlices, (Py_ssize_t)i);
			if(o == 0) throw CombErStr(strEr_BadFuncArg, strErCom);
			if(!PyList_Check(o)) throw CombErStr(strEr_BadFuncArg, strErCom);

			PyObject *oP = PyList_GetItem(o, 0);
			if(oP == 0) throw CombErStr(strEr_BadFuncArg, strErComCenPt);
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arAuxP, 3, CombErStr(strEr_BadFuncArg, strErComCenPt));
			vCrd.push_back(arAuxP[0]); vCrd.push_back(arAuxP[1]); vCrd.push_back(arAuxP[2]);

			PyObject *oD = PyList_GetItem(o, 1);
			if(oD == 0) throw CombErStr(strEr_BadFuncArg, strErComDims);
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oD, 'd', arAuxD, 2, CombErStr(strEr_BadFuncArg, strErComDims));
			vDims.push_back(arAuxD[0]); vDims.push_back(arAuxD[1]);
		}

		int nCrd = (int)vCrd.size();
		if(nCrd != ns*3) throw CombErStr(strEr_BadFuncArg, strErComCenPt);
		int nDims = (int)vDims.size();
		if(nDims != (ns << 1)) throw CombErStr(strEr_BadFuncArg, strErComDims);
		CAuxParse::DoubleVect2Arr(vCrd, arCrd);
		CAuxParse::DoubleVect2Arr(vDims, arDims);

		double arM[] = {0.,0.,0.};
		ParseM(arM, oM, "ObjMltExtRtg");

		int ind = 0;
		g_pyParse.ProcRes(RadObjMltExtRtg(&ind, arCrd, arDims, ns, arM));
		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arCrd != 0) delete[] arCrd;
	if(arDims != 0) delete[] arDims;
	return oResInd;
}

/************************************************************************//**
* Creates triangulated extruded polygon block, i.e. an extruded polygon with its bases subdivided by triangulation.
***************************************************************************/
static PyObject* radia_ObjMltExtTri(PyObject *self, PyObject *args)
{
	PyObject *oPgn=0, *oSbd=0, *oOrnt=0, *oM=0, *oOpt=0, *oResInd=0;
	double *arCrd=0, *arSbd=0;
	try
	{	
		double xc = 0, lx = 0;
		if(!PyArg_ParseTuple(args, "ddOO|OOO:ObjMltExtTri", &xc, &lx, &oPgn, &oSbd, &oOrnt, &oM, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtTri");
		if((lx <= 0) || (oPgn == 0) || (oSbd == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtTri");

		int nCrd=0;
		if(!(CPyParse::CopyPyNestedListElemsToNumAr(oPgn, 'd', arCrd, nCrd))) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtTri, incorrect polygon definition");
		int nv = nCrd >> 1;
		if(nCrd != (nv << 1)) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtTri, incorrect polygon definition (total number of coordinates should be even)");

		int nSbd=0;
		if(!(CPyParse::CopyPyNestedListElemsToNumAr(oSbd, 'd', arSbd, nSbd))) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtTri, incorrect definition of subdivision parameters");
		int nSegm = nSbd >> 1;
		if(nSbd != (nSegm << 1)) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtTri, incorrect definition of subdivision parameters");
		if(nv != nSegm) throw CombErStr(strEr_BadFuncArg, ": ObjMltExtTri, inconsistent definition of polygon and subdivision parameters for its segments");

		char a = ParseOrnt(oOrnt, 'x', "ObjMltExtTri");

		double arM[] = {0.,0.,0.};
		ParseM(arM, oM, "ObjMltExtTri");

		char sOpt[1024]; *sOpt = '\0';
		CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		int ind = 0;
		g_pyParse.ProcRes(RadObjMltExtTri(&ind, xc, lx, arCrd, arSbd, nv, a, arM, sOpt));
		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arCrd != 0) delete[] arCrd;
	if(arSbd != 0) delete[] arSbd;
	return oResInd;
}

/************************************************************************//**
* Creates triangulated extruded polygon block, i.e. an extruded polygon with its bases subdivided by triangulation.
***************************************************************************/
static PyObject* radia_ObjCylMag(PyObject *self, PyObject *args)
{
	PyObject *oP=0, *oOrnt=0, *oM=0, *oResInd=0;

	try
	{	
		double r = 0, h = 0;
		int nseg = 0;

		if(!PyArg_ParseTuple(args, "Oddi|OO:ObjCylMag", &oP, &r, &h, &nseg, &oOrnt, &oM)) throw CombErStr(strEr_BadFuncArg, ": ObjCylMag");
		if((oP == 0) || (r <= 0) || (h <= 0) || (nseg <= 0)) throw CombErStr(strEr_BadFuncArg, ": ObjCylMag");

		double arP[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjCylMag, incorrect definition of center point"));

		char a = ParseOrnt(oOrnt, 'z', "ObjCylMag");

		double arM[] = {0.,0.,0.};
		ParseM(arM, oM, "ObjCylMag");

		int ind = 0;
		g_pyParse.ProcRes(RadObjCylMag(&ind, arP, r, h, nseg, a, arM));

		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Field Source: Rectangular Parallelepiped with color {RGB[0],RGB[1],RGB[2]}.
 * The block is magnetized according to M[3] then subdivided according to K[3] and added into the container grp (grp should be defined in advance by calling RadObjCnt()).
 ***************************************************************************/
static PyObject* radia_ObjFullMag(PyObject *self, PyObject *args)
{
	PyObject *oP=0, *oL=0, *oM=0, *oK=0, *oRGB=0, *oResInd=0;

	try
	{
		int indGrp=0, indMat=0;
		if(!PyArg_ParseTuple(args, "OOOOiiO:ObjFullMag", &oP, &oL, &oM, &oK, &indGrp, &indMat, &oRGB)) throw CombErStr(strEr_BadFuncArg, ": ObjFullMag");
		if((oP == 0) || (oL == 0) || (oM == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjFullMag");

		double arP[3], arL[3], arM[3], arK[9], arRGB[3];
		bool lenIsSmall = false;
		int lenK = 9;
		double *p = arK;
		CPyParse::CopyPyListElemsToNumArray(oK, 'd', p, lenK, lenIsSmall);
		if(lenIsSmall) throw CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of subdivision parameters");
	
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of center point"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oL, 'd', arL, 3, CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of dimensions"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oM, 'd', arM, 3, CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of magnetization vector"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oRGB, 'd', arRGB, 3, CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of RGB color"));

		int ind = 0;
		g_pyParse.ProcRes(RadObjFullMag(&ind, arP, arL, arM, arK, lenK, indGrp, indMat, arRGB));

		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
* Creates a current carrying rectangular parallelepiped block.
***************************************************************************/
static PyObject* radia_ObjRecCur(PyObject *self, PyObject *args)
{
	PyObject *oP=0, *oL=0, *oJ=0, *oResInd=0;

	try
	{	
		if(!PyArg_ParseTuple(args, "OOO:ObjRecCur", &oP, &oL, &oJ)) throw CombErStr(strEr_BadFuncArg, ": ObjRecCur");
		if((oP == 0) || (oL == 0) || (oJ == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjRecCur");

		double arP[3], arL[3], arJ[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjRecCur, incorrect definition of center point"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oL, 'd', arL, 3, CombErStr(strEr_BadFuncArg, ": ObjRecCur, incorrect definition of dimensions"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oJ, 'd', arJ, 3, CombErStr(strEr_BadFuncArg, ": ObjRecCur, incorrect definition of current density"));

		int ind = 0;
		g_pyParse.ProcRes(RadObjRecCur(&ind, arP, arL, arJ));
		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
* Creates a current carrying finite-length arc of rectangular cross-section.
***************************************************************************/
static PyObject* radia_ObjArcCur(PyObject *self, PyObject *args)
{
	PyObject *oP=0, *oR=0, *oPhi=0, *oManAuto=0, *oOrnt=0, *oResInd=0;

	try
	{
		double h = 0, j = 0;
		int nseg = 0;
		if(!PyArg_ParseTuple(args, "OOOdid|OO:ObjArcCur", &oP, &oR, &oPhi, &h, &nseg, &j, &oManAuto, &oOrnt)) throw CombErStr(strEr_BadFuncArg, ": ObjArcCur");
		if((oP == 0) || (oR == 0) || (oPhi == 0) || (h <= 0) || (nseg < 1)) throw CombErStr(strEr_BadFuncArg, ": ObjArcCur");

		double arP[3], arR[2], arPhi[2];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjArcCur, incorrect definition of center point"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oR, 'd', arR, 2, CombErStr(strEr_BadFuncArg, ": ObjArcCur, incorrect definition of radii"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oPhi, 'd', arPhi, 2, CombErStr(strEr_BadFuncArg, ": ObjArcCur, incorrect definition of angles"));

		char sManAuto[256];
		sManAuto[0] = 'm';
		if(oManAuto != 0) 
		{
			CPyParse::CopyPyStringToC(oManAuto, sManAuto, 256);
			if((sManAuto[0] != 'm') && (sManAuto[0] != 'a')) throw CombErStr(strEr_BadFuncArg, ": ObjArcCur, incorrect definition of \'man\' or \'auto\' parameter");
		}

		char a = ParseOrnt(oOrnt, 'z', "ObjArcCur");

		int ind = 0;
		g_pyParse.ProcRes(RadObjArcCur(&ind, arP, arR, arPhi, h, nseg, sManAuto[0], a, j));
		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Field Source: Recetrack conductor with rectangular cross-secton and constant current density over the volume
 ***************************************************************************/
static PyObject* radia_ObjRaceTrk(PyObject *self, PyObject *args)
{//The coil consists of four 90-degree bents connected by four straight parts parallel to the XY plane.

	PyObject *oP=0, *oR=0, *oL=0, *oManAuto=0, *oOrnt=0, *oResInd=0;

	try
	{
		double h=0, curDens=0;
		int nseg=0;
		if(!PyArg_ParseTuple(args, "OOOdid|OO:ObjRaceTrk", &oP, &oR, &oL, &h, &nseg, &curDens, &oManAuto, &oOrnt)) throw CombErStr(strEr_BadFuncArg, ": ObjRaceTrk");
		if((oP == 0) || (oR == 0) || (oL == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjRaceTrk");

		double arP[3], arR[2], arL[2];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjRaceTrk, incorrect definition of center point"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oR, 'd', arR, 2, CombErStr(strEr_BadFuncArg, ": ObjRaceTrk, incorrect definition of radii"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oL, 'd', arL, 2, CombErStr(strEr_BadFuncArg, ": ObjRaceTrk, incorrect definition of straight interval lengths"));

		char sManAuto[256];
		sManAuto[0] = 'm';
		if(oManAuto != 0) 
		{
			CPyParse::CopyPyStringToC(oManAuto, sManAuto, 256);
			if((sManAuto[0] != 'm') && (sManAuto[0] != 'a')) throw CombErStr(strEr_BadFuncArg, ": ObjRaceTrk, incorrect definition of \'man\' or \'auto\' parameter");
		}

		char a = ParseOrnt(oOrnt, 'z', "ObjRaceTrk");

		int ind = 0;
		g_pyParse.ProcRes(RadObjRaceTrk(&ind, arP, arR, arL, h, nseg, *sManAuto, a, curDens));

		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
* Creates a filament polygonal line conductor with current.
***************************************************************************/
static PyObject* radia_ObjFlmCur(PyObject *self, PyObject *args)
{
	PyObject *oPts=0, *oResInd=0;
	double *arCrd=0;

	try
	{
		double I = 0;
		if(!PyArg_ParseTuple(args, "Od:ObjFlmCur", &oPts, &I)) throw CombErStr(strEr_BadFuncArg, ": ObjFlmCur");
		if(oPts == 0) throw CombErStr(strEr_BadFuncArg, ": ObjFlmCur, incorrect definition of points");

		int nCrd=0;
		if(!(CPyParse::CopyPyNestedListElemsToNumAr(oPts, 'd', arCrd, nCrd))) throw CombErStr(strEr_BadFuncArg, ": ObjFlmCur, incorrect definition of points");
		int np = (int)round(nCrd/3.);
		if(nCrd != np*3) throw CombErStr(strEr_BadFuncArg, ": ObjFlmCur, incorrect definition of points (each point should be defined by 3 cartesian coordinates)");

		int ind = 0;
		g_pyParse.ProcRes(RadObjFlmCur(&ind, arCrd, np, I));
		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arCrd != 0) delete[] arCrd;
	return oResInd;
}

/************************************************************************//**
* Scales current (density) in a 3D object by multiplying it by a constant.
***************************************************************************/
static PyObject* radia_ObjScaleCur(PyObject *self, PyObject *args)
{
	PyObject *oResInd=0;

	try
	{
		int ind = 0;
		double scaleCoef = 0;
		if(!PyArg_ParseTuple(args, "id:ObjScaleCur", &ind, &scaleCoef)) throw CombErStr(strEr_BadFuncArg, ": ObjScaleCur");
		if(ind <= 0) throw CombErStr(strEr_BadFuncArg, ": ObjScaleCur, incorrect object index");

		g_pyParse.ProcRes(RadObjScaleCur(ind, scaleCoef));
		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
* Creates a source of uniform background magnetic field.
***************************************************************************/
static PyObject* radia_ObjBckg(PyObject *self, PyObject *args)
{
	PyObject *oB=0, *oResInd=0;

	try
	{
		if(!PyArg_ParseTuple(args, "O:ObjBckg", &oB)) throw CombErStr(strEr_BadFuncArg, ": ObjBckg");
		if(oB == 0) throw CombErStr(strEr_BadFuncArg, ": ObjBckg");

		double arB[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oB, 'd', arB, 3, CombErStr(strEr_BadFuncArg, ": ObjBckg, incorrect definition of magnetic field strength"));

		int ind = 0;
		g_pyParse.ProcRes(RadObjBckg(&ind, arB));
		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Field Source: Container of Magnetic Field Sources
 ***************************************************************************/
static PyObject* radia_ObjCnt(PyObject *self, PyObject *args)
{
	PyObject *oInds=0, *oResInd=0;
	int *arInds=0;

	try
	{
		if(!PyArg_ParseTuple(args, "O:ObjCnt", &oInds)) throw CombErStr(strEr_BadFuncArg, ": ObjCnt");
		if(oInds == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCnt");

		int nInds = 0;
		bool lenIsSmall = false;
		CPyParse::CopyPyListElemsToNumArray(oInds, 'i', arInds, nInds, lenIsSmall);

		int ind = 0;
		g_pyParse.ProcRes(RadObjCnt(&ind, arInds, nInds));
		oResInd = Py_BuildValue("i", ind);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arInds != 0) delete[] arInds;
	return oResInd;
}

/************************************************************************//**
 * Adds objects to the container object
 ***************************************************************************/
static PyObject* radia_ObjAddToCnt(PyObject *self, PyObject *args)
{
	PyObject *oInds=0, *oResInd=0;
	int *arInds=0;

	try
	{
		int indCnt = 0;
		if(!PyArg_ParseTuple(args, "iO:ObjAddToCnt", &indCnt, &oInds)) throw CombErStr(strEr_BadFuncArg, ": ObjAddToCnt");
		if((indCnt <= 0) || (oInds == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjAddToCnt");

		int nInds = 0;
		bool lenIsSmall = false;
		CPyParse::CopyPyListElemsToNumArray(oInds, 'i', arInds, nInds, lenIsSmall);

		g_pyParse.ProcRes(RadObjAddToCnt(indCnt, arInds, nInds));
		oResInd = Py_BuildValue("i", indCnt);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arInds != 0) delete[] arInds;
	return oResInd;
}

/************************************************************************//**
 * Calculates the number of objects in the container cnt.
 ***************************************************************************/
//static PyObject* radia_ObjCntSize(PyObject *self, PyObject *args)
//{
//	PyObject *oResInd=0;
//
//	try
//	{
//		int cnt = 0;
//
//		if(!PyArg_ParseTuple(args, "i:ObjCntSize" , &cnt)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjCntSize");
//		if(cnt == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCntSize");
//
//		int ind = 0;
//		g_pyParse.ProcRes(RadObjCntSize(&ind, cnt));
//
//		oResInd = Py_BuildValue("i", ind);
//		Py_XINCREF(oResInd); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return oResInd;
//}

/************************************************************************//**
 * Fills-in an array with the reference numbers of objects present in the container cnt
 ***************************************************************************/
//static PyObject* radia_ObjCntStuf(PyObject *self, PyObject *args)
//{
//	PyObject *oResInd=0;
//
//	try
//	{
//		int cnt = 0;
//
//		if(!PyArg_ParseTuple(args, "i:ObjCntStuf" ,&cnt)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjCntStuf");
//		if(cnt == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCntStuf");
//
//		int ind = 0;
//		g_pyParse.ProcRes(RadObjCntStuf(&ind, cnt));
//
//		oResInd = Py_BuildValue("i", ind);
//		Py_XINCREF(oResInd); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return oResInd;
//}

/************************************************************************//**
 * Duplicates the object obj
 ***************************************************************************/
//static PyObject* radia_ObjDpl(PyObject *self, PyObject *args)
//{
//	PyObject *object_opt=0, *oResInd=0;
//
//	try
//	{
//		int obj = 0;
//
//		if(!PyArg_ParseTuple(args, "iO:ObjDpl" ,&obj ,&object_opt)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjDpl");
//
//		char option[256];
//		CPyParse::CopyPyStringToC(object_opt, option, 256);
//
//		int ind = 0;
//		g_pyParse.ProcRes(RadObjDpl(&ind, obj, option));
//
//		oResInd = Py_BuildValue("i", ind);
//		Py_XINCREF(oResInd); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return oResInd;
//}

/************************************************************************//**
 * Provides coordinates of geometrical center point and magnetization of the object obj
 ***************************************************************************/
//static PyObject* radia_ObjM(PyObject *self, PyObject *args)
//{
//	PyObject *oResInd=0;
//
//	try
//	{
//		int obj = 0;
//
//		if(!PyArg_ParseTuple(args, "i:ObjM", &obj)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjM");
//		if(obj == 0) throw CombErStr(strEr_BadFuncArg, ": ObjM");
//
//		double M = 0;
//		//double *M_ptr = 0;
//		g_pyParse.ProcRes(RadObjM(&M, obj));
//
//		oResInd = Py_BuildValue("d", M);
//		Py_XINCREF(oResInd); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return oResInd;
//}

/************************************************************************//**
 * Provides coordinates of geometrical center point and magnetization of the object obj.
 ***************************************************************************/
//static PyObject* radia_ObjCenFld(PyObject *self, PyObject *args)
//{
//	PyObject *oResInd=0;
//
//	try
//	{
//		int obj = 0;
//		char type;
//
//		if(!PyArg_ParseTuple(args, "ic:ObjCenFld" ,&obj ,&type)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjCenFld");
//		if(obj == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCenFld");
//
//		double B = 0;
//		//double *B_ptr = 0;
//
//		g_pyParse.ProcRes(RadObjCenFld(&B, obj, type));
//
//		oResInd = Py_BuildValue("d", B);
//		Py_XINCREF(oResInd); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return oResInd;
//}

/************************************************************************//**
 * Provides coordinates of geometrical center point and magnetization of the object obj.
 ***************************************************************************/
//static PyObject* radia_ObjSetM(PyObject *self, PyObject *args)
//{
//	PyObject *oResInd=0;
//
//	try
//	{
//		int obj = 0;
//		if(!PyArg_ParseTuple(args, "i:ObjSetM" ,&obj )) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjSetM");
//
//		double M = 0;
//		g_pyParse.ProcRes(RadObjSetM(obj, &M));
//
//		oResInd = Py_BuildValue("d", M);
//		Py_XINCREF(oResInd); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return oResInd;
//}

/************************************************************************//**
 * Cuts the object obj by a plane passing through a given point normally to a given vector.
 ***************************************************************************/
//static PyObject* radia_ObjCutMag(PyObject *self, PyObject *args)
//{
//	PyObject *objectResult_nobj=0, *objectResult_Objs=0;
//	PyObject *object_P=0, *object_N=0, *object_opt=0;
//
//	try
//	{
//		int obj = 0;
//
//		if(!PyArg_ParseTuple(args, "iOOO:ObjCutMag" ,&obj ,&object_P
//			,&object_N ,&object_opt)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjCutMag");
//		if(obj == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCutMag");
//
//		int n_pts = 3;
//		double array_P[n_pts];
//		double array_N[n_pts];
//
//		double *array_P_ptr = array_P;
//		double *array_N_ptr = array_N;
//
//		CPyParse::CopyPyListElemsToNumArray(object_P, 'd', array_P_ptr, n_pts);
//		CPyParse::CopyPyListElemsToNumArray(object_N, 'd', array_N_ptr, n_pts);
//
//		char option[256];
//		CPyParse::CopyPyStringToC(object_opt, option, 256);
//
//		int n_obj = 0;
//		int Objs = 0;
//
//		g_pyParse.ProcRes(RadObjCutMag(&Objs, &n_obj, obj, array_P, array_N, option));
//
//		objectResult_nobj = Py_BuildValue("i", n_obj);
//		objectResult_Objs = Py_BuildValue("i", Objs);
//
//		Py_XINCREF(objectResult_nobj);
//		Py_XINCREF(objectResult_Objs); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return objectResult_Objs; //, objectResult_nobj;
//}

/************************************************************************//**
 * Subdivides (segments) the object obj by 3 sets of parallel planes.
 ***************************************************************************/
//static PyObject* radia_ObjDivMagPln(PyObject *self, PyObject *args)
//{
//	PyObject *object_SbdPar=0, *object_FlatNorm=0;
//	PyObject *object_option=0, *objectResult_n=0;
//
//	try
//	{
//		int obj = 0;
//		int nSbdPar = 0;
//
//		if(!PyArg_ParseTuple(args, "iOiOO:ObjDivMagPln", &obj, &object_SbdPar,
//			&nSbdPar, &object_FlatNorm, &object_option)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjDivMagPln");
//
//		int n_flatnorm = 9;
//		double array_SbdPar[nSbdPar];
//		double array_FlatNorm[n_flatnorm];
//
//		double *array_SbdPar_ptr = 0;
//		double *array_FlatNorm_ptr = 0;
//
//		CPyParse::CopyPyListElemsToNumArray(object_SbdPar, 'd', array_SbdPar_ptr, nSbdPar);
//		CPyParse::CopyPyListElemsToNumArray(object_FlatNorm, 'd', array_FlatNorm_ptr, n_flatnorm);
//		
//		char option[256];
//		CPyParse::CopyPyStringToC(object_option, option, 256);
//
//		int n = 0;
//		g_pyParse.ProcRes(RadObjDivMagPln(&n, obj, array_SbdPar,
//			nSbdPar, array_FlatNorm, option));
//
//		objectResult_n = Py_BuildValue("i", n);
//		Py_XINCREF(objectResult_n); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return objectResult_n;
//}

/************************************************************************//**
 * Subdivides (segments) the object obj by a set of coaxial elliptic cylinders. 
 ***************************************************************************/
//static PyObject* radia_ObjDivMagCyl(PyObject *self, PyObject *args)
//{
//	PyObject *objectResult_index=0;
//	PyObject *object_SbdPar=0, *object_FlatCylPar=0, *object_option=0;
//
//	try
//	{
//		int obj = 0;
//		int nSbdPar = 0;
//		int rat = 0;
//
//		if(!PyArg_ParseTuple(args, "iOiOdO:ObjDivMagCyl", &obj, &object_SbdPar,
//			&nSbdPar, &object_FlatCylPar, &rat, &object_option)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjDivMagCyl");
//
//		int n_flatcyl = 9;
//		double array_SbdPar[nSbdPar];
//		double array_FlatCylPar[n_flatcyl];
//
//		double *ptr_SbdPar = 0;
//		double *ptr_FlatCylPar = 0;
//
//		CPyParse::CopyPyListElemsToNumArray(object_SbdPar, 'd', ptr_SbdPar, nSbdPar);
//		CPyParse::CopyPyListElemsToNumArray(object_FlatCylPar, 'd', ptr_FlatCylPar, n_flatcyl);
//
//		char option[256];
//		CPyParse::CopyPyStringToC(object_option, option, 256);
//
//		int n = 0;
//		g_pyParse.ProcRes(RadObjDivMagCyl(&n, obj, array_SbdPar, nSbdPar,
//			array_FlatCylPar, rat, option));
//
//		objectResult_index = Py_BuildValue("i", n);
//		Py_XINCREF(objectResult_index); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return objectResult_index;
//}

/************************************************************************//**
 * Computes geometrical volume of a 3D object.
 ***************************************************************************/
//static PyObject* radia_ObjGeoVol(PyObject *self, PyObject *args)
//{
//	PyObject *objectResult=0;
//
//	try
//	{
//		int obj = 0;
//
//		if(!PyArg_ParseTuple(args, "i:ObjGeoVol", &obj)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjGeoVol");
//
//
//		double v = 0;
//		g_pyParse.ProcRes(RadObjGeoVol(&v, obj));
//
//		objectResult = Py_BuildValue("d", v);
//		Py_XINCREF(objectResult); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return objectResult;
//}

/************************************************************************//**
 *  Computes coordinates of object extrimities in laboratory frame.
 ***************************************************************************/
//static PyObject* radia_ObjGeoLim(PyObject *self, PyObject *args)
//{
//	PyObject *objectResult_L=0;
//
//	try
//	{
//		int obj = 0;
//
//		if(!PyArg_ParseTuple(args, "i:ObjGeoLim", &obj)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjGeoLim");
//
//		double array_L[6];
//
//		g_pyParse.ProcRes(RadObjGeoLim(array_L, obj));
//
//		objectResult_L = Py_BuildValue("d", array_L);
//		Py_XINCREF(objectResult_L); //?
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		//PyErr_PrintEx(1);
//	}
//
//	return objectResult_L;
//}

/************************************************************************//**
 * Starts an application for viewing of 3D geometry of the object obj
 ***************************************************************************/
//static PyObject* radia_ObjDrwQD3D(PyObject *self, PyObject *args)
//{
//	PyObject *obj_option=0;
//	PyObject *objectError=0;
//	try
//	{
//		int obj = 0;
//
//		if(!PyArg_ParseTuple(args, "iO:ObjDrwQD3D", &obj, &obj_option)) 
//			throw CombErStr(strEr_BadFuncArg, ": ObjDrwQD3D");
//
//		char option[256];
//		CPyParse::CopyPyStringToC(obj_option, option, 256);
//
//		g_pyParse.ProcRes(RadObjDrwQD3D(obj, option));
//		objectError = Py_BuildValue("i", 0);
//
//	}
//	catch(const char* erText) 
//	{
//		PyErr_SetString(PyExc_RuntimeError, erText);
//		objectError = Py_BuildValue("i", 1);
//		//PyErr_PrintEx(1);
//	}
//
//	return objectError;
//}

/************************************************************************//**
 * Magnetic Field Sources: Assigns drawing attributes - RGB color (r,g,b) and line thickness thcn - to a Magnetic Field Source object
 ***************************************************************************/
static PyObject* radia_ObjDrwAtr(PyObject *self, PyObject *args)
{
	PyObject *oInd=0, *oRGB=0;

	try
	{
		double thcn=0.001;
		if(!PyArg_ParseTuple(args, "OO|d:ObjDrwAtr", &oInd, &oRGB, &thcn)) throw CombErStr(strEr_BadFuncArg, ": ObjDrwAtr");
		if((oInd == 0) || (oRGB == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjDrwAtr");

		if(!PyNumber_Check(oInd)) throw CombErStr(strEr_BadFuncArg, ": ObjDrwAtr");
		int ind = (int)PyLong_AsLong(oInd);

		double arRGB[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oRGB, 'd', arRGB, 3, CombErStr(strEr_BadFuncArg, ": ObjDrwAtr, incorrect definition of RGB color"));

		g_pyParse.ProcRes(RadObjDrwAtr(ind, arRGB, thcn));
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oInd;
}

/************************************************************************//**
 * Magnetic Field Sources: Starts an application for viewing of 3D geometry of the object obj (the viewer is based on the GLUT / OpenGL graphics library)
 ***************************************************************************/
static PyObject* radia_ObjDrwOpenGL(PyObject *self, PyObject *args)
{
	PyObject *oInd=0, *oOpt=0;

	try
	{
		if(!PyArg_ParseTuple(args, "O|O:ObjDrwOpenGL", &oInd, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjDrwOpenGL");
		if(oInd == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDrwOpenGL");

		if(!PyNumber_Check(oInd)) throw CombErStr(strEr_BadFuncArg, ": ObjDrwOpenGL");
		int ind = (int)PyLong_AsLong(oInd);

		char sOpt[1024]; *sOpt = '\0';
		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		g_pyParse.ProcRes(RadObjDrwOpenGL(ind, sOpt));
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oInd;
}

/************************************************************************//**
 * Space Transformations: Creates an object mirror with respect to a plane.
 * The object mirror possesses the same geometry as obj, but its magnetization and/or current densities 
 * are modified in such a way that the magnetic field produced by the obj and its mirror in the plane 
 * of mirroring is PERPENDICULAR to this plane.
 ***************************************************************************/
static PyObject* radia_TrfZerPara(PyObject *self, PyObject *args)
{
	PyObject *oP=0, *oN=0, *oResInd=0;

	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "iOO:TrfZerPara", &ind, &oP, &oN)) throw CombErStr(strEr_BadFuncArg, ": TrfZerPara");
		if((oP == 0) || (oN == 0)) throw CombErStr(strEr_BadFuncArg, ": TrfZerPara");

		double arP[3], arN[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjRaceTrk, incorrect definition of point in symmery plane"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oN, 'd', arN, 3, CombErStr(strEr_BadFuncArg, ": ObjRaceTrk, incorrect definition of normal vector to symmery plane"));

		int indRes = 0;
		g_pyParse.ProcRes(RadTrfZerPara(&indRes, ind, arP, arN));
		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Space Transformations: Creates an object mirror with respect to a plane.
 * The object mirror possesses the same geometry as obj, but its magnetization and/or current densities 
 * are modified in such a way that the magnetic field produced by the obj and its mirror in the plane 
 * of mirroring is PARALLEL to this plane.
 ***************************************************************************/
static PyObject* radia_TrfZerPerp(PyObject *self, PyObject *args)
{
	PyObject *oP=0, *oN=0, *oResInd=0;

	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "iOO:TrfZerPerp", &ind, &oP, &oN)) throw CombErStr(strEr_BadFuncArg, ": TrfZerPerp");
		if((oP == 0) || (oN == 0)) throw CombErStr(strEr_BadFuncArg, ": TrfZerPerp");

		double arP[3], arN[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjRaceTrk, incorrect definition of point in symmery plane"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oN, 'd', arN, 3, CombErStr(strEr_BadFuncArg, ": ObjRaceTrk, incorrect definition of normal vector to symmery plane"));

		int indRes = 0;
		g_pyParse.ProcRes(RadTrfZerPerp(&indRes, ind, arP, arN));

		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Materials: a pre-defined magnetic material (the material is identified by its name/formula, e.g. \"NdFeB\").
 ***************************************************************************/
static PyObject* radia_MatStd(PyObject *self, PyObject *args)
{
	PyObject *oMatId=0, *oResInd=0;

	try
	{
		double absM = 0;
		if(!PyArg_ParseTuple(args, "Od:MatStd", &oMatId, &absM)) throw CombErStr(strEr_BadFuncArg, ": MatStd");
		if(oMatId == 0) throw CombErStr(strEr_BadFuncArg, ": MatStd");

		char sMatId[256];
		CPyParse::CopyPyStringToC(oMatId, sMatId, 256);

		int indRes = 0;
		g_pyParse.ProcRes(RadMatStd(&indRes, sMatId, absM));

		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Materials: a nonlinear isotropic magnetic material with the M versus H curve defined by the list of pairs corresponding values of H and M (H1,M1,H2,M2,...)
 ***************************************************************************/
static PyObject* radia_MatSatIsoTab(PyObject *self, PyObject *args)
{
	PyObject *oMatData=0, *oResInd=0;
	double *arMatData=0;

	try
	{
		if(!PyArg_ParseTuple(args, "O:MatSatIsoTab", &oMatData)) throw CombErStr(strEr_BadFuncArg, ": MatSatIsoTab");
		if(oMatData == 0) throw CombErStr(strEr_BadFuncArg, ": MatSatIsoTab");

		int nDataTot = 0;
		char resP = CPyParse::CopyPyNestedListElemsToNumAr(oMatData, 'd', arMatData, nDataTot);
		if(resP == 0) throw CombErStr(strEr_BadFuncArg, ": MatSatIsoTab");
		int nDataP = (int)round(0.5*nDataTot);

		int indRes = 0;
		g_pyParse.ProcRes(RadMatSatIsoTab(&indRes, arMatData, nDataP));

		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Materials: computes magnetization from magnetic field strength vector for a given material
 ***************************************************************************/
static PyObject* radia_MatMvsH(PyObject *self, PyObject *args)
{
	PyObject *oCmpnId=0, *oH=0, *oResM=0;

	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "iOO:MatMvsH", &ind, &oCmpnId, &oH)) throw CombErStr(strEr_BadFuncArg, ": MatMvsH");
		if((oCmpnId == 0) || (oH == 0)) throw CombErStr(strEr_BadFuncArg, ": MatMvsH");

		char sCmpnId[256];
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 256);

		double arH[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oH, 'd', arH, 3, CombErStr(strEr_BadFuncArg, ": MatMvsH, incorrect definition of magnetic field strength vector"));

		double arM[6];
		int lenM = 3;
		g_pyParse.ProcRes(RadMatMvsH(arM, &lenM, ind, sCmpnId, arH));

		if(lenM == 1) oResM = Py_BuildValue("d", *arM);
		else if(lenM > 1) oResM = CPyParse::SetDataListOfLists(arM, lenM, 1);
		if(oResM) Py_XINCREF(oResM);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResM;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Builds an interaction matrix and performs a relaxation procedure
 ***************************************************************************/
static PyObject* radia_Solve(PyObject *self, PyObject *args)
{
	PyObject *oRes=0;

	try
	{
		int ind=0, numIt=1000, meth=0;
		double prec=0.0001;

		if(!PyArg_ParseTuple(args, "idi|i:Solve", &ind, &prec, &numIt, &meth)) throw CombErStr(strEr_BadFuncArg, ": Solve");

		double arResSolve[12];
		int lenResSolve = 4;
		g_pyParse.ProcRes(RadSolve(arResSolve, &lenResSolve, ind, prec, numIt, meth));

		if(lenResSolve == 1) oRes = Py_BuildValue("d", *arResSolve);
		else if(lenResSolve > 1) oRes = CPyParse::SetDataListOfLists(arResSolve, lenResSolve, 1);
		if(oRes) Py_XINCREF(oRes);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	
	return oRes;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes field created by the object obj at one or many points
 ***************************************************************************/
static PyObject* radia_Fld(PyObject *self, PyObject *args)
{
	PyObject *oCmpnId=0, *oP=0, *oResB=0;
	double *arCoord=0, *arB=0;

	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "iOO:Fld", &ind, &oCmpnId, &oP)) throw CombErStr(strEr_BadFuncArg, ": Fld");

		char sCmpnId[256];
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 256);

		int nCoord = 0;
		char resP = CPyParse::CopyPyNestedListElemsToNumAr(oP, 'd', arCoord, nCoord);
		if(resP == 0) throw CombErStr(strEr_BadFuncArg, ": Fld: array / list of points");

		int nP = (int)round(nCoord/3.);
		if(nP*3 != nCoord) throw CombErStr(strEr_BadFuncArg, ": Fld: array / list of points");

		const int maxNumFldCmpn = 14;
		arB = new double[maxNumFldCmpn*nP];
		int nB = 0;
		g_pyParse.ProcRes(RadFld(arB, &nB, ind, sCmpnId, arCoord, nP));

		if(nB == 1) oResB = Py_BuildValue("d", *arB);
		else if(nB > 1) oResB = CPyParse::SetDataListOfLists(arB, nB, nP);
		if(oResB) Py_XINCREF(oResB);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arCoord != 0) delete[] arCoord;
	if(arB != 0) delete[] arB;

	return oResB;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes magnetic field integral produced by the object obj along a straight line specified by points P1 and P2.
 * Depending on the InfOrFin variable value, the integral is infinite ("inf") or finite ("fin"), from P1 to P2; the field integral component is specified by the id input variable. The unit is T*mm.
 ***************************************************************************/
static PyObject* radia_FldInt(PyObject *self, PyObject *args)
{
	PyObject *oInfOrFin=0, *oCmpnId=0, *oP1=0, *oP2=0, *oResIB=0;

	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "iOOOO:FldInt", &ind, &oInfOrFin, &oCmpnId, &oP1, &oP2)) throw CombErStr(strEr_BadFuncArg, ": FldInt");
		if((oInfOrFin == 0) || (oCmpnId == 0) || (oP1 == 0) || (oP2 == 0)) throw CombErStr(strEr_BadFuncArg, ": FldInt");

		char sInfOrFin[256], sCmpnId[256];
		CPyParse::CopyPyStringToC(oInfOrFin, sInfOrFin, 256);
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 256);

		double arP1[3], arP2[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP1, 'd', arP1, 3, CombErStr(strEr_BadFuncArg, ": FldInt, incorrect definition of first point defining integration line"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP2, 'd', arP2, 3, CombErStr(strEr_BadFuncArg, ": FldInt, incorrect definition of second point defining integration line"));

		double arIB[9];
		int nIB = 0;
		g_pyParse.ProcRes(RadFldInt(arIB, &nIB, ind, sInfOrFin, sCmpnId, arP1, arP2));

		if(nIB == 1) oResIB = Py_BuildValue("d", *arIB);
		else if(nIB > 1) oResIB = CPyParse::SetDataListOfLists(arIB, nIB, 1);
		if(oResIB) Py_XINCREF(oResIB);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResIB;
}

/************************************************************************//**
 * Python C API stuff: module & method definition2, etc.
 ***************************************************************************/

static PyMethodDef radia_methods[] = {
	{"ObjRecMag", radia_ObjRecMag, METH_VARARGS, "ObjRecMag() instantiates Rectangular Parallelepiped with Constant Magnetizatiom over volume"},
	{"ObjThckPgn", radia_ObjThckPgn, METH_VARARGS, "ObjThckPgn() instantiates a uniformly magnetized extruded polygon"},
	{"ObjPolyhdr", radia_ObjPolyhdr, METH_VARARGS, "ObjPolyhdr() creates a uniformly magnetized polyhedron (closed volume limited by planes)"},
	{"ObjArcPgnMag", radia_ObjArcPgnMag, METH_VARARGS, "ObjArcPgnMag() creates a uniformly magnetized finite-length arc of polygonal cross-section"},
	{"ObjMltExtPgn", radia_ObjMltExtPgn, METH_VARARGS, "ObjMltExtPgn() attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on slices."},
	{"ObjMltExtRtg", radia_ObjMltExtRtg, METH_VARARGS, "ObjMltExtRtg() attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on rectangular slices"},
	{"ObjMltExtTri", radia_ObjMltExtTri, METH_VARARGS, "ObjMltExtTri() creates triangulated extruded polygon block, i.e. an extruded polygon with its bases subdivided by triangulation"},
	{"ObjCylMag", radia_ObjCylMag, METH_VARARGS, "ObjCylMag() creates a cylindrical magnet"},
	{"ObjFullMag", radia_ObjFullMag, METH_VARARGS, "ObjFullMag() instantiates Rectangular Parallelepiped with Constant Magnetizatiom over volume, subdivided, and added to a group, with some color assigned"},
	
	{"ObjRecCur", radia_ObjRecCur, METH_VARARGS, "ObjRecCur() creates a current carrying rectangular parallelepiped block"},
	{"ObjArcCur", radia_ObjArcCur, METH_VARARGS, "ObjArcCur() creates a current carrying finite-length arc of rectangular cross-section"},
	{"ObjRaceTrk", radia_ObjRaceTrk, METH_VARARGS, "ObjRaceTrk() instantiates Recetrack conductor with rectangular cross-secton and constant current density over the volume"},
	{"ObjFlmCur", radia_ObjFlmCur, METH_VARARGS, "ObjFlmCur() creates a filament polygonal line conductor with current"},
	
	{"ObjScaleCur", radia_ObjScaleCur, METH_VARARGS, "ObjScaleCur() scales current (density) in a 3D object by multiplying it by a constant"},
	{"ObjBckg", radia_ObjBckg, METH_VARARGS, "ObjBckg() creates a source of uniform background magnetic field"},
	{"ObjCnt", radia_ObjCnt, METH_VARARGS, "ObjCnt() instantiates Container of Magnetic Field Sources"},
	{"ObjAddToCnt", radia_ObjAddToCnt, METH_VARARGS, "ObjAddToCnt() adds objects to the container object cnt"},
	//{"ObjCntSize", radia_ObjCntSize, METH_VARARGS, "ObjCntSize() Calculates the number of objects in the container cnt"},
	//{"ObjCntStuf", radia_ObjCntStuf, METH_VARARGS, "ObjCntStuf() Fills-in an array with the reference numbers of objects present in the container cnt"}, 
	//{"ObjDpl", radia_ObjDpl, METH_VARARGS, "ObjDpl() Duplicates the object obj"},
	//{"ObjM", radia_ObjM, METH_VARARGS, "ObjM() Provides coordinates of geometrical center point and magnetization of the object obj"},
	//{"ObjCenFld", radia_ObjCenFld, METH_VARARGS, "ObjCenFld() Provides coordinates of geometrical center point and magnetization of the object obj"},
	//{"ObjSetM", radia_ObjSetM, METH_VARARGS, "ObjSetM() Provides coordinates of geometrical center point and magnetization of the object obj."},

	//{"ObjCutMag", radia_ObjCutMag, METH_VARARGS, "ObjCutMag() Cuts the object obj by a plane passing through a given point normally to a given vector."},
	//{"ObjDivMagPln", radia_ObjDivMagPln, METH_VARARGS, "ObjDivMagPln() Subdivides (segments) the object obj by 3 sets of parallel planes. "},
	//{"ObjDivMagCyl", radia_ObjDivMagCyl, METH_VARARGS, "ObjDivMagCyl() Subdivides (segments) the object obj by a set of coaxial elliptic cylinders. "},
	//{"ObjGeoVol", radia_ObjGeoVol, METH_VARARGS, "ObjGeoVol() Computes geometrical volume of a 3D object."},
	//{"ObjGeoLim", radia_ObjGeoLim, METH_VARARGS, "ObjGeoLim() Computes coordinates of object extrimities in laboratory frame"},
	
	{"ObjDrwAtr", radia_ObjDrwAtr, METH_VARARGS, "ObjDrwAtr() assigns drawing attributes - RGB color (r,g,b) and line thickness thcn - to a Magnetic Field Source object"},
	{"ObjDrwOpenGL", radia_ObjDrwOpenGL, METH_VARARGS, "ObjDrwOpenGL() assigns drawing attributes - RGB color (r,g,b) and line thickness thcn - to a Magnetic Field Source object"},
	//{"ObjDrwQD3D", radia_ObjDrwQD3D, METH_VARARGS, "ObjDrwQD3D() Starts an application for viewing of 3D geometry of the object obj"},

	{"TrfZerPara", radia_TrfZerPara, METH_VARARGS, "TrfZerPara() creates an object mirror with respect to a plane. The object mirror possesses the same geometry as obj, but its magnetization and/or current densities are modified in such a way that the magnetic field produced by the obj and its mirror in the plane of mirroring is perpendicular to this plane"},
	{"TrfZerPerp", radia_TrfZerPerp, METH_VARARGS, "TrfZerPerp() creates an object mirror with respect to a plane. The object mirror possesses the same geometry as obj, but its magnetization and/or current densities are modified in such a way that the magnetic field produced by the obj and its mirror in the plane of mirroring is parallel to this plane"},

	{"MatStd", radia_MatStd, METH_VARARGS, "MatStd() creates a pre-defined magnetic material (the material is identified by its name/formula, e.g. \"NdFeB\")"},
	{"MatSatIsoTab", radia_MatSatIsoTab, METH_VARARGS, "MatSatIsoTab() creates a nonlinear isotropic magnetic material with the M versus H curve defined by the list of pairs corresponding values of H and M [[H1,M1],[H2,M2],...]"},
	{"MatMvsH", radia_MatMvsH, METH_VARARGS, "MatMvsH() computes magnetization from magnetic field strength vector for a given material"},

	{"Solve", radia_Solve, METH_VARARGS, "Solve() solves a magnetostatic problem, i.e. builds an interaction matrix and performs a relaxation procedure"},

	{"Fld", radia_Fld, METH_VARARGS,  "Fld() computes field created by the object obj at one or many points"},
	{"FldInt", radia_FldInt, METH_VARARGS, "FldInt() computes magnetic field integral produced by magnetic field source object along a straight line"},

	{NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef radiamodule = {
    PyModuleDef_HEAD_INIT,
    "radia",
    "radia module is Python binding of Radia 3D Magnetistatics code / library",
    -1,
    radia_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_radia(void)
{ 
	//setting pointer to function to be eventually called from SRWLIB
	//srwlUtiSetWfrModifFunc(&ModifySRWLWfr);

    return PyModule_Create(&radiamodule);
}

#else

PyMODINIT_FUNC initradia(void)
{
	//setting pointer to function to be eventually called from SRWLIB
	//srwlUtiSetWfrModifFunc(&ModifySRWLWfr);

	Py_InitModule("radia", radia_methods);
}

#endif
