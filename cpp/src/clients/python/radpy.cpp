/************************************************************************//**
 * File: radpy.cpp
 * Description: Python binding
 * Project: Radia
 * First release: June 2018
 *
 * @authors O. Chubar (BNL), J. Edelen (RadiaSoft), G. Le Bec (ESRF)
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
 * Auxiliary function for parsing Subdivision Params, including eventual gradients
 ***************************************************************************/
static void ParseSubdPar(double arSbdPar[6], PyObject* oSbdPar, const char* sFuncName=0)
//static void ParseSubdPar(double arSbdPar[6], PyObject* oSbdPar, char* sFuncName=0)
{//OC29022020
	char sErrMes[2000];
	strcpy(sErrMes, ": \0");
	if(sFuncName == 0) strcat(sErrMes, "ObjDivMag\0");
	else strcat(sErrMes, sFuncName);

	strcat(sErrMes, ", incorrect subdivision parameters");

	if(oSbdPar == 0) throw CombErStr(strEr_BadFuncArg, sErrMes);

	double arSbdParLoc[6];
	double *pSbdParLoc = arSbdParLoc;
	int nSbdParLoc = 6;
	char resP = CPyParse::CopyPyNestedListElemsToNumAr(oSbdPar, 'd', pSbdParLoc, nSbdParLoc);
	if(resP == 0) throw CombErStr(strEr_BadFuncArg, sErrMes);
	
	if(nSbdParLoc == 6)
	{
		for(int i=0; i<nSbdParLoc; i++) arSbdPar[i] = arSbdParLoc[i];
		return;
	}

	int arLenSbdPar[6];
	int *pLenSbdPar = arLenSbdPar;
	int nLenSbdPar = 6;
	CPyParse::FindLengthsOfElemListsOrArrays(oSbdPar, pLenSbdPar, nLenSbdPar, true);

	vector<double> vSbdPar;
	for(int i=0; i<nSbdParLoc; i++) vSbdPar.push_back(arSbdParLoc[i]);

	for(int iLen=0; iLen<nLenSbdPar; iLen++)
	{
		if(arLenSbdPar[iLen] == 0) vSbdPar.insert(vSbdPar.begin() + (iLen*2 + 1), 1.);
	}

	for(int j=0; j<6; j++) arSbdPar[j] = vSbdPar[j];
}

/************************************************************************//**
 * Magnetic Field Source: Rectangular Parallelepiped with Constant Magnetizatiom over volume
 ***************************************************************************/
static PyObject* radia_ObjRecMag(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjThckPgn(PyObject* self, PyObject* args)
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

		char a = 'x';
		double arM[] = {0.,0.,0.};
		if((oOrnt != 0) && (oM == 0)) //OC29022020
		{
			if(PyList_Check(oOrnt))
			{
				oM = oOrnt; oOrnt = 0;
			}
		}

		a = ParseOrnt(oOrnt, 'x', "ObjThckPgn"); //OC29022020
		//char a = ParseOrnt(oOrnt, 'x', "ObjThckPgn");

		//double arM[] = {0.,0.,0.};
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
static PyObject* radia_ObjPolyhdr(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjArcPgnMag(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjMltExtPgn(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjMltExtRtg(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjMltExtTri(PyObject* self, PyObject* args)
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

		char a = 'x';
		double arM[] = {0.,0.,0.};
		char sOrnt[1024]; *sOrnt = '\0';
		char sOpt[1024]; *sOpt = '\0';
		bool MagnIsDef = false;
		bool OptIsDef = false;
		//oOrnt can be Orientation or Magnetization or Option
		if(oOrnt != 0)
		{
			bool isListOrAr = PyList_Check(oOrnt);
			if(!isListOrAr) isListOrAr = PyObject_CheckBuffer(oOrnt);
			if(isListOrAr)
			{
				ParseM(arM, oOrnt, "ObjMltExtTri");
				MagnIsDef = true;
			}
			else
			{
				CPyParse::CopyPyStringToC(oOrnt, sOrnt, 1024);
				if(strlen(sOrnt) <= 1) a = ParseOrnt(oOrnt, 'x', "ObjMltExtTri");
				else
				{
					strcpy(sOpt, sOrnt);
					OptIsDef = true;
				}
			}
		}
		//oM can be Magnetization or Option
		if((oM != 0) && ((!MagnIsDef) || (!OptIsDef)))
		{
			bool isListOrAr = PyList_Check(oM);
			if(!isListOrAr) isListOrAr = PyObject_CheckBuffer(oM);
			if(isListOrAr) ParseM(arM, oM, "ObjMltExtTri");
			else
			{
				CPyParse::CopyPyStringToC(oM, sOpt, 1024);
				OptIsDef = true;
			}
		}

		if((oOpt != 0) && (!OptIsDef)) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

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
static PyObject* radia_ObjCylMag(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjFullMag(PyObject* self, PyObject* args)
{
	PyObject *oP=0, *oL=0, *oM=0, *oK=0, *oRGB=0, *oResInd=0;

	try
	{
		int indGrp=0, indMat=0;
		if(!PyArg_ParseTuple(args, "OOOOiiO:ObjFullMag", &oP, &oL, &oM, &oK, &indGrp, &indMat, &oRGB)) throw CombErStr(strEr_BadFuncArg, ": ObjFullMag");
		if((oP == 0) || (oL == 0) || (oM == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjFullMag");

		double arP[3], arL[3], arM[3], arK[9], arRGB[3];
		//bool lenIsSmall = false;
		//int lenK = 9;
		//double *p = arK;
		//CPyParse::CopyPyListElemsToNumArray(oK, 'd', p, lenK, lenIsSmall);
		//if(lenIsSmall) throw CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of subdivision parameters");
		//CPyParse::CopyPyNestedListElemsToNumAr(oK, 'd', p, lenK); //OC27022020
		//if((lenK != 3) && (lenK != 6)) throw CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of subdivision parameters");
		ParseSubdPar(arK, oK, "ObjFullMag\0"); //OC29022020
	
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of center point"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oL, 'd', arL, 3, CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of dimensions"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oM, 'd', arM, 3, CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of magnetization vector"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oRGB, 'd', arRGB, 3, CombErStr(strEr_BadFuncArg, ": ObjFullMag, incorrect definition of RGB color"));

		int ind = 0;
		g_pyParse.ProcRes(RadObjFullMag(&ind, arP, arL, arM, arK, 6, indGrp, indMat, arRGB)); //OC29022020
		//g_pyParse.ProcRes(RadObjFullMag(&ind, arP, arL, arM, arK, lenK, indGrp, indMat, arRGB));

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
static PyObject* radia_ObjRecCur(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjArcCur(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjRaceTrk(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjFlmCur(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjScaleCur(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjBckg(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjCnt(PyObject* self, PyObject* args)
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
static PyObject* radia_ObjAddToCnt(PyObject* self, PyObject* args)
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
 * Calculates the number of objects in the container.
 ***************************************************************************/
static PyObject* radia_ObjCntSize(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	//PyObject *oOpt=0, *oRes=0;
	try
	{
		int indCnt = 0;
		if(!PyArg_ParseTuple(args, "i|O:ObjCntSize" , &indCnt)) throw CombErStr(strEr_BadFuncArg, ": ObjCntSize");
		if(indCnt == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCntSize");

		//char sOpt[1024]; *sOpt = '\0';
		//if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		int sizeCnt = 0;
		g_pyParse.ProcRes(RadObjCntSize(&sizeCnt, indCnt));
		//g_pyParse.ProcRes(RadObjCntSize(&sizeCnt, indCnt, sOpt));

		oRes = Py_BuildValue("i", sizeCnt);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oRes;
}

/************************************************************************//**
 * Returns list with the reference numbers of objects present in the container.
 ***************************************************************************/
static PyObject* radia_ObjCntStuf(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	int *arInds=0;
	try
	{
		int indCnt = 0;
		if(!PyArg_ParseTuple(args, "i:ObjCntStuf" ,&indCnt)) throw CombErStr(strEr_BadFuncArg, ": ObjCntStuf");
		if(indCnt == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCntStuf");

		int sizeCnt = 0;
		g_pyParse.ProcRes(RadObjCntSize(&sizeCnt, indCnt));

		if(sizeCnt <= 0)
		{
			oRes = PyList_New(0);
		}
		else
		{
			arInds = new int[sizeCnt];
			g_pyParse.ProcRes(RadObjCntStuf(arInds, indCnt));
			oRes = CPyParse::SetDataListOfLists(arInds, sizeCnt, 1, (char*)"i");
		}
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arInds != 0) delete[] arInds;
	return oRes;
}

/************************************************************************//**
 * Duplicates the object obj
 ***************************************************************************/
static PyObject* radia_ObjDpl(PyObject* self, PyObject* args)
{
	PyObject *oOpt=0, *oResInd=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "i|O:ObjDpl", &ind, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjDpl");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDpl");

		char sOpt[1024]; *sOpt = '\0';
		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		int indDpl = 0;
		g_pyParse.ProcRes(RadObjDpl(&indDpl, ind, sOpt));

		oResInd = Py_BuildValue("i", indDpl);
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
 * Provides coordinates of geometrical center point and magnetization of the object obj
 ***************************************************************************/
static PyObject* radia_ObjM(PyObject* self, PyObject* args)
{
	PyObject *oResM=0;
	double *arPM=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "i:ObjM", &ind)) throw CombErStr(strEr_BadFuncArg, ": ObjM");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjM");

		int arMesh[21];
		g_pyParse.ProcRes(RadObjM(0, arMesh, ind)); //OC27092018
		//g_pyParse.ProcRes(RadObjM(arMesh, ind));
		
		long nDim = *arMesh;
		if(nDim > 0)
		{
			long long nTot = 1;
			for(int i=1; i<=nDim; i++) nTot *= arMesh[i];
			arPM = new double[nTot];

			g_pyParse.ProcRes(RadUtiDataGet((char*)arPM, (char*)"mad", ind)); //OC27092018
			//g_pyParse.ProcRes(RadUtiDataGet(arPM, ind));

			double *arPMorig = arPM;
			oResM = CPyParse::SetDataListsNested(arPM, arMesh, (char*)"d");
			arPM = arPMorig;
		}
		if(oResM) Py_XINCREF(oResM);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arPM != 0) delete[] arPM;
	return oResM;
}

/************************************************************************//**
 * Provides coordinates of geometrical center point and magnetic field at that point.
 ***************************************************************************/
static PyObject* radia_ObjCenFld(PyObject* self, PyObject* args)
{
	PyObject *oCmpnId=0, *oRes=0;
	double *arPB=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "iO:ObjCenFld", &ind, &oCmpnId)) throw CombErStr(strEr_BadFuncArg, ": ObjCenFld");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCenFld");

		char sCmpnId[256];
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 256);

		int arMesh[21];
		g_pyParse.ProcRes(RadObjCenFld(0, arMesh, ind, *sCmpnId)); //OC27092018
		//g_pyParse.ProcRes(RadObjCenFld(arMesh, ind, *sCmpnId));

		long nDim = *arMesh;
		if(nDim > 0)
		{
			long long nTot = 1;
			for(int i=1; i<=nDim; i++) nTot *= arMesh[i];
			arPB = new double[nTot];

			g_pyParse.ProcRes(RadUtiDataGet((char*)arPB, (char*)"mad", ind)); //OC27092018
			//g_pyParse.ProcRes(RadUtiDataGet(arPB, ind));
			double *arPBorig = arPB;
			oRes = CPyParse::SetDataListsNested(arPB, arMesh, (char*)"d");
			arPB = arPBorig;
		}
		if(oRes) Py_XINCREF(oRes);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arPB != 0) delete[] arPB;
	return oRes;
}

/************************************************************************//**
 * Sets magnetization of the object obj.
 ***************************************************************************/
static PyObject* radia_ObjSetM(PyObject* self, PyObject* args)
{
	PyObject *oM=0, *oResInd=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "iO:ObjSetM", &ind, &oM)) throw CombErStr(strEr_BadFuncArg, ": ObjSetM");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjSetM");

		double arM[] = {0.,0.,0.};
		ParseM(arM, oM, "ObjSetM");

		g_pyParse.ProcRes(RadObjSetM(ind, arM));

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
 * Cuts the object obj by a plane passing through a given point perpendicularly to a given vector.
 ***************************************************************************/
static PyObject* radia_ObjCutMag(PyObject* self, PyObject* args)
{
	PyObject *oP=0, *oN=0, *oOpt=0, *oRes=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "iOO|O:ObjCutMag", &ind, &oP, &oN, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjCutMag");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjCutMag");

		double arP[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": ObjCutMag, incorrect definition of point in the cutting plane"));
		double arN[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oN, 'd', arN, 3, CombErStr(strEr_BadFuncArg, ": ObjCutMag, incorrect definition of cutting plane normal"));

		char sOpt[1024]; *sOpt = '\0';
		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		int arInds[10];
		int nObj = 0;
		g_pyParse.ProcRes(RadObjCutMag(arInds, &nObj, ind, arP, arN, sOpt));

		oRes = CPyParse::SetDataListOfLists(arInds, nObj, 1, (char*)"i");
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oRes;
}

/************************************************************************//**
 * Subdivides (segments) the object obj by 3 sets of parallel planes.
 ***************************************************************************/
static PyObject* radia_ObjDivMagPln(PyObject* self, PyObject* args)
{
	PyObject *oSbdPar=0, *oN1=0, *oN2=0, *oN3=0, *oOpt=0, *oResInd=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "iO|OOOO:ObjDivMagPln", &ind, &oSbdPar, &oN1, &oN2, &oN3, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjDivMagPln");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDivMagPln");

		if((oN1 != 0) && (oN2 == 0) && (oN3 == 0) && (oOpt == 0))
		{//To suport call like: rad.ObjDivMagPln(mag01, [[2,0.5],[3,0.2],[4,0.1]], 'Frame->Lab')
			oOpt = oN1; oN1 = 0;
		}

		double arSbdPar[6];
		//double *pSbdPar = arSbdPar;
		//int nSbdPar = 6;
		//char resP = CPyParse::CopyPyNestedListElemsToNumAr(oSbdPar, 'd', pSbdPar, nSbdPar);
		//if(resP == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDivMagPln, incorrect definition of cutting plane normal vectors");
		ParseSubdPar(arSbdPar, oSbdPar, "ObjDivMagPln\0"); //OC29022020

		double arN1N2N3[] = {1,0,0, 0,1,0, 0,0,1};
		if(oN1 != 0)
		{
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oN1, 'd', arN1N2N3, 3, CombErStr(strEr_BadFuncArg, ": ObjDivMagPln, incorrect definition of first cutting plane normal vector"));
		}
		if(oN2 != 0)
		{
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oN2, 'd', arN1N2N3 + 3, 3, CombErStr(strEr_BadFuncArg, ": ObjDivMagPln, incorrect definition of second cutting plane normal vector"));
		}
		if(oN3 != 0)
		{
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oN3, 'd', arN1N2N3 + 6, 3, CombErStr(strEr_BadFuncArg, ": ObjDivMagPln, incorrect definition of third cutting plane normal vector"));
		}

		char sOpt[1024]; *sOpt = '\0';
		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		int indNew = 0;
		g_pyParse.ProcRes(RadObjDivMagPln(&indNew, ind, arSbdPar, 6, arN1N2N3, sOpt)); //OC29022020
		//g_pyParse.ProcRes(RadObjDivMagPln(&indNew, ind, arSbdPar, nSbdPar, arN1N2N3, sOpt));

		oResInd = Py_BuildValue("i", indNew);
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
 * Subdivides (segments) the object obj by a set of coaxial elliptic cylinders. 
 ***************************************************************************/
static PyObject* radia_ObjDivMagCyl(PyObject* self, PyObject* args)
{
	PyObject *oSbdPar=0, *oA=0, *oV=0, *oP=0, *oOpt=0, *oResInd=0;
	try
	{
		int ind = 0;
		double rat = 0;
		if(!PyArg_ParseTuple(args, "iOOOOd|O:ObjDivMagCyl", &ind, &oSbdPar, &oA, &oV, &oP, &rat, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjDivMagCyl");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDivMagCyl");

		double arSbdPar[6];
		//double *pSbdPar = arSbdPar;
		//int nSbdPar = 6;
		//char resP = CPyParse::CopyPyNestedListElemsToNumAr(oSbdPar, 'd', pSbdPar, nSbdPar);
		//if(resP == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDivMagCyl");
		ParseSubdPar(arSbdPar, oSbdPar, "ObjDivMagCyl\0"); //OC29022020

		double arAVP[] = {0,0,0, 0,0,0, 0,0,0};
		if(oA != 0)
		{
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oA, 'd', arAVP, 3, CombErStr(strEr_BadFuncArg, ": ObjDivMagCyl, incorrect definition of point on cylinder axis"));
		}
		if(oV != 0)
		{
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oV, 'd', arAVP + 3, 3, CombErStr(strEr_BadFuncArg, ": ObjDivMagCyl, incorrect definition of vector of cylinder axis"));
		}
		if(oP != 0)
		{
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arAVP + 6, 3, CombErStr(strEr_BadFuncArg, ": ObjDivMagCyl, incorrect definition of point in elliptical cylinder base"));
		}

		char sOpt[1024]; *sOpt = '\0';
		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		int indNew = 0;
		g_pyParse.ProcRes(RadObjDivMagCyl(&indNew, ind, arSbdPar, 6, arAVP, rat, sOpt)); //OC29022020
		//g_pyParse.ProcRes(RadObjDivMagCyl(&indNew, ind, arSbdPar, nSbdPar, arAVP, rat, sOpt));

		oResInd = Py_BuildValue("i", indNew);
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
 * Subdivides (segments) the object obj by a set of coaxial elliptic cylinders. 
 ***************************************************************************/
static PyObject* radia_ObjDivMag(PyObject* self, PyObject* args)
{
	PyObject *oSbdPar=0, *oType=0, *oDir=0, *oOpt=0, *oResInd=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "iO|OOO:ObjDivMag", &ind, &oSbdPar, &oType, &oDir, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjDivMag");
		//if(!PyArg_ParseTuple(args, "iOOO|O:ObjDivMag", &ind, &oSbdPar, &oType, &oDir, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjDivMag");
		if((ind == 0) || (oSbdPar == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjDivMag");

		double arSbdPar[6];
		//double *pSbdPar = arSbdPar;
		//int nSbdPar = 6;
		//char resP = CPyParse::CopyPyNestedListElemsToNumAr(oSbdPar, 'd', pSbdPar, nSbdPar);
		//if(resP == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDivMag");	
		ParseSubdPar(arSbdPar, oSbdPar, "ObjDivMag\0"); //OC29022020

		//char sType[1024];
		//strcpy(sType, "pln\0");
		char sOpt[1024]; *sOpt = '\0';
		if((oType != 0) && (oDir == 0) && (oOpt == 0)) 
		{//Treating the case like rad.ObjDivMag(u, ndiv, 'Frame->LabTot')
			if(oType != 0) CPyParse::CopyPyStringToC(oType, sOpt, 1024);
			char *pEndOptName = strrchr(sOpt, '>');
			if(pEndOptName != 0)
			{
				if(*(--pEndOptName) == '-') oType = 0;
			}
		}

		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		//char sType[32]; *sType = '\0';
		char sType[] = "pln";
		if(oType != 0) CPyParse::CopyPyStringToC(oType, sType, 32);

		int indNew = 0;
		if((strcmp(sType, "pln") == 0) || (strcmp(sType, "Pln") == 0) || (strcmp(sType, "PLN") == 0))
		{
			double arN1N2N3[] = {1,0,0, 0,1,0, 0,0,1};
			if(oDir != 0)
			{
				double *pN1N2N3 = arN1N2N3;
				int nCrd = 6;
				char resDir = CPyParse::CopyPyNestedListElemsToNumAr(oDir, 'd', pN1N2N3, nCrd);
				if((resDir == 0) || (nCrd != 6)) throw CombErStr(strEr_BadFuncArg, ": ObjDivMag, incorrect definition of cutting plane normal vectors");
			}
			g_pyParse.ProcRes(RadObjDivMagPln(&indNew, ind, arSbdPar, 6, arN1N2N3, sOpt)); //OC29022020
			//g_pyParse.ProcRes(RadObjDivMagPln(&indNew, ind, arSbdPar, nSbdPar, arN1N2N3, sOpt));
		}
		else if((strcmp(sType, "cyl") == 0) || (strcmp(sType, "Cyl") == 0) || (strcmp(sType, "CYL") == 0))
		{
			if(oDir == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDivMag, incorrect definition of parameters for subdivision by elliptical cylinders");

			double arAVPr[] = {0,0,0, 0,0,0, 0,0,0, 0};
			double *pAVPr = arAVPr;
			int nCrdAVPr = 10;
			char resDir = CPyParse::CopyPyNestedListElemsToNumAr(oDir, 'd', pAVPr, nCrdAVPr);
			if((resDir == 0) || (nCrdAVPr != 10)) throw CombErStr(strEr_BadFuncArg, ": ObjDivMag, incorrect definition of parameters for subdivision by elliptical cylinders");

			g_pyParse.ProcRes(RadObjDivMagCyl(&indNew, ind, arSbdPar, 6, arAVPr, arAVPr[9], sOpt)); //OC29022020
			//g_pyParse.ProcRes(RadObjDivMagCyl(&indNew, ind, arSbdPar, nSbdPar, arAVPr, arAVPr[9], sOpt));
		}
		else throw CombErStr(strEr_BadFuncArg, ": ObjDivMag, incorrect definition of subdivision (segmentation) type");

		oResInd = Py_BuildValue("i", indNew);
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
 * Computes geometrical volume of a 3D object.
 ***************************************************************************/
static PyObject* radia_ObjGeoVol(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "i:ObjGeoVol", &ind)) throw CombErStr(strEr_BadFuncArg, ": ObjGeoVol");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjGeoVol");

		double v = 0;
		g_pyParse.ProcRes(RadObjGeoVol(&v, ind));

		oRes = Py_BuildValue("d", v);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oRes;
}

/************************************************************************//**
 *  Computes coordinates of object extrimities in laboratory frame.
 ***************************************************************************/
static PyObject* radia_ObjGeoLim(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "i:ObjGeoLim", &ind)) throw CombErStr(strEr_BadFuncArg, ": ObjGeoLim");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjGeoVol");

		double arLim[6];
		g_pyParse.ProcRes(RadObjGeoLim(arLim, ind));

		oRes = CPyParse::SetDataListOfLists(arLim, 6, 1, (char*)"d");
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oRes;
}

/************************************************************************//**
 *  Computes coordinates of object extrimities in laboratory frame.
 ***************************************************************************/
static PyObject* radia_ObjDegFre(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	try
	{
		int ind = 0;
		if(!PyArg_ParseTuple(args, "i:ObjDegFre", &ind)) throw CombErStr(strEr_BadFuncArg, ": ObjDegFre");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDegFre");

		int numDegFre = 0;
		g_pyParse.ProcRes(RadObjDegFre(&numDegFre, ind));

		oRes = Py_BuildValue("i", numDegFre);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oRes;
}

/************************************************************************//**
 * Magnetic Field Sources: Assigns drawing attributes - RGB color (r,g,b) and line thickness thcn - to a Magnetic Field Source object
 ***************************************************************************/
static PyObject* radia_ObjDrwAtr(PyObject* self, PyObject* args)
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
		Py_XINCREF(oInd); //?
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
static PyObject* radia_ObjDrwOpenGL(PyObject* self, PyObject* args)
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
		Py_XINCREF(oInd); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		PyErr_PrintEx(1); //OC04102018 Need to clear the error in this case, to let Py script runing further if ObjDrwOpenGL is not implemented
	}
	return oInd;
}

/************************************************************************//**
 * Magnetic Field Sources: outputs data for viewing of 3D geometry of an object using the VTK software system
 ***************************************************************************/
static PyObject* ParseGeomDataDrwVTK(double* arCrdV, int nV, int* arLenP, float* arColP, int nP) //OC05112019 (from by R. Nagler's static PyObject* radia_ObjGeometry_helper)
{
	PyObject *lengths = PyList_New(nP);
	PyObject *colors = PyList_New(nP*3);
	PyObject *vertices = PyList_New(nV*3);
	if(!vertices || !lengths || !colors) throw strEr_MAF;

	for(Py_ssize_t i = PyList_Size(lengths); --i >= 0;) 
	{
		PyObject *o = PyLong_FromLong(arLenP[i]);
		if(o == NULL || PyList_SetItem(lengths, i, o) < 0) throw strEr_MAF;
	}

	for(Py_ssize_t i = PyList_Size(colors); --i >= 0;) 
	{
		PyObject *o = PyFloat_FromDouble(arColP[i]); //OC04062020 (from Dan Abell) //why colors are float at input if they are casted to long here?
		//PyObject *o = PyLong_FromDouble(arColP[i]); //why colors are float at input if they are casted to long here?
		if(o == NULL || PyList_SetItem(colors, i, o) < 0) throw strEr_MAF;
	}

	for(Py_ssize_t i = PyList_Size(vertices); --i >= 0;) 
	{
		PyObject *o = PyFloat_FromDouble(arCrdV[i]);
		if(o == NULL || PyList_SetItem(vertices, i, o) < 0) throw strEr_MAF;
	}

	PyObject *oRes = Py_BuildValue(
		"{s:N,s:N,s:N}",
		"colors",
		colors,
		"lengths",
		lengths,
		"vertices",
		vertices
	);
	if(oRes == NULL) throw strEr_MAF;
	return oRes;
}

/************************************************************************//**
 * Magnetic Field Sources: outputs data for viewing of 3D geometry of an object using the VTK software system
 ***************************************************************************/
static PyObject* radia_ObjDrwVTK(PyObject* self, PyObject* args) //OC03112019 (requested by R. Nagler)
//static PyObject* radia_ObjGeometry(PyObject *self, PyObject *args) 
{
	PyObject *oInd=0, *oOpt=0, *oRes=0;
	double *arCrdVP=0, *arCrdVL=0;
	float *arColP=0, *arColL=0;
	int *arLenP=0, *arLenL=0;

	try
	{
		if(!PyArg_ParseTuple(args, "O|O:ObjDrwVTK", &oInd, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjDrwVTK");
		//if(!PyArg_ParseTuple(args, "O|O:ObjGeometry", &oInd, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": ObjGeometry");
		if(oInd == 0) throw CombErStr(strEr_BadFuncArg, ": ObjDrwVTK");

		if(!PyNumber_Check(oInd)) throw CombErStr(strEr_BadFuncArg, ": ObjDrwVTK");
		int ind = (int)PyLong_AsLong(oInd);

		char sOpt[1024]; *sOpt = '\0';
		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		int nVertPgns=0, nPgns=0, nVertLines=0, nLines=0, key=-1;
		g_pyParse.ProcRes(RadObjDrwVTK(&nVertPgns, &nPgns, &nVertLines, &nLines, &key, ind, sOpt));

		if(nVertPgns > 0) arCrdVP = new double[nVertPgns*3];
		if(nPgns > 0)
		{
			arLenP = new int[nPgns];
			arColP = new float[nPgns*3];
		}
		if(nVertLines > 0) arCrdVL = new double[nVertLines*3];
		if(nLines > 0)
		{
			arLenL = new int[nLines];
			arColL = new float[nLines*3];
		}

		g_pyParse.ProcRes(RadObjDrwDataGetVTK(arCrdVP, arLenP, arColP, arCrdVL, arLenL, arColL, key));

		oRes = Py_BuildValue(
			"{s:N,s:N}",
			"polygons",
			ParseGeomDataDrwVTK(arCrdVP, nVertPgns, arLenP, arColP, nPgns),
			"lines",
			ParseGeomDataDrwVTK(arCrdVL, nVertLines, arLenL, arColL, nLines)
		);
		if(oRes == NULL) throw strEr_MAF;
		//TODO(robnagler) ref counts are invalid at this point,
		//  but then this is a malloc error...

		Py_XINCREF(oRes); //? NOTE: seem to have experienced crashes without this
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}

	if(arCrdVP != 0) delete[] arCrdVP;
	if(arCrdVL != 0) delete[] arCrdVL;
	if(arColP != 0) delete[] arColP;
	if(arColL != 0) delete[] arColL;
	if(arLenP != 0) delete[] arLenP;
	if(arLenL != 0) delete[] arLenL;

	return oRes;
}

/************************************************************************//**
 * Space Transformations: Creates a symmetry with respect to plane defined by a point and a normal vector.
 ***************************************************************************/
static PyObject* radia_TrfPlSym(PyObject* self, PyObject* args)
{
	PyObject *oP=0, *oN=0, *oResInd=0;
	try
	{
		if(!PyArg_ParseTuple(args, "OO:TrfPlSym", &oP, &oN)) throw CombErStr(strEr_BadFuncArg, ": TrfPlSym");
		if((oP == 0) || (oN == 0)) throw CombErStr(strEr_BadFuncArg, ": TrfPlSym");

		double arP[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": TrfPlSym, incorrect definition of point in the symmetry plane"));

		double arN[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oN, 'd', arN, 3, CombErStr(strEr_BadFuncArg, ": TrfPlSym, incorrect definition of vector normal to the symmetry plane"));

		int ind = 0;
		g_pyParse.ProcRes(RadTrfPlSym(&ind, arP, arN));

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
 * Space Transformations: Creates a rotation about an axis.
 ***************************************************************************/
static PyObject* radia_TrfRot(PyObject* self, PyObject* args)
{
	PyObject *oP=0, *oV=0, *oResInd=0;
	try
	{
		double phi = 0;
		if(!PyArg_ParseTuple(args, "OOd:TrfRot", &oP, &oV, &phi)) throw CombErStr(strEr_BadFuncArg, ": TrfRot");
		if((oP == 0) || (oV == 0)) throw CombErStr(strEr_BadFuncArg, ": TrfRot");

		double arP[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": TrfRot, incorrect definition of point on axis of rotation"));

		double arV[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oV, 'd', arV, 3, CombErStr(strEr_BadFuncArg, ": TrfRot, incorrect definition of rotation axis vector"));

		int ind = 0;
		g_pyParse.ProcRes(RadTrfRot(&ind, arP, arV, phi));

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
 * Space Transformations: Creates a translation.
 ***************************************************************************/
static PyObject* radia_TrfTrsl(PyObject* self, PyObject* args)
{
	PyObject *oV=0, *oResInd=0;
	try
	{
		if(!PyArg_ParseTuple(args, "O:TrfTrsl", &oV)) throw CombErStr(strEr_BadFuncArg, ": TrfTrsl");
		if(oV == 0) throw CombErStr(strEr_BadFuncArg, ": TrfTrsl");

		double arV[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oV, 'd', arV, 3, CombErStr(strEr_BadFuncArg, ": TrfTrsl, incorrect definition of translation vector"));

		int ind = 0;
		g_pyParse.ProcRes(RadTrfTrsl(&ind, arV));

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
 * Space Transformations: Creates a field inversion.
 ***************************************************************************/
static PyObject* radia_TrfInv(PyObject* self, PyObject* args)
{
	PyObject *oResInd=0;
	try
	{
		int ind = 0;
		g_pyParse.ProcRes(RadTrfInv(&ind));

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
 * Space Transformations: Multiplies original space transformation by another transformation from left.
 ***************************************************************************/
static PyObject* radia_TrfCmbL(PyObject* self, PyObject* args)
{
	PyObject *oResInd=0;
	try
	{
		int indTrfOrig = 0, indTrf = 0;
		if(!PyArg_ParseTuple(args, "ii:TrfCmbL", &indTrfOrig, &indTrf)) throw CombErStr(strEr_BadFuncArg, ": TrfCmbL");
		if((indTrfOrig == 0) || (indTrf == 0)) throw CombErStr(strEr_BadFuncArg, ": TrfCmbL");

		int ind = 0;
		g_pyParse.ProcRes(RadTrfCmbL(&ind, indTrfOrig, indTrf));

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
 * Space Transformations: Multiplies original space transformation by another transformation from right.
 ***************************************************************************/
static PyObject* radia_TrfCmbR(PyObject* self, PyObject* args)
{
	PyObject *oResInd=0;
	try
	{
		int indTrfOrig = 0, indTrf = 0;
		if(!PyArg_ParseTuple(args, "ii:TrfCmbR", &indTrfOrig, &indTrf)) throw CombErStr(strEr_BadFuncArg, ": TrfCmbR");
		if((indTrfOrig == 0) || (indTrf == 0)) throw CombErStr(strEr_BadFuncArg, ": TrfCmbR");

		int ind = 0;
		g_pyParse.ProcRes(RadTrfCmbR(&ind, indTrfOrig, indTrf));

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
 * Space Transformations: Creates mlt-1 symmetry objects of the object obj. 
 * Each symmetry object is derived from the previous one by applying the transformation trf to it. 
 * Following this, the object obj becomes equivalent to mlt different objects.
 ***************************************************************************/
static PyObject* radia_TrfMlt(PyObject* self, PyObject* args)
{
	PyObject *oResInd=0;
	try
	{
		int indObj = 0, indTrf = 0, mlt = 0;
		if(!PyArg_ParseTuple(args, "iii:TrfMlt", &indObj, &indTrf, &mlt)) throw CombErStr(strEr_BadFuncArg, ": TrfMlt");
		if((indObj == 0) || (indTrf == 0) || (mlt == 0)) throw CombErStr(strEr_BadFuncArg, ": TrfMlt");

		int ind = 0;
		g_pyParse.ProcRes(RadTrfMlt(&ind, indObj, indTrf, mlt));

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
 * Space Transformations: Orients object obj by applying a space transformation to it once. 
 ***************************************************************************/
static PyObject* radia_TrfOrnt(PyObject* self, PyObject* args)
{
	PyObject *oResInd=0;
	try
	{
		int indObj = 0, indTrf = 0;
		if(!PyArg_ParseTuple(args, "ii:TrfOrnt", &indObj, &indTrf)) throw CombErStr(strEr_BadFuncArg, ": TrfOrnt");
		if((indObj == 0) || (indTrf == 0)) throw CombErStr(strEr_BadFuncArg, ": TrfOrnt");

		int ind = 0;
		g_pyParse.ProcRes(RadTrfOrnt(&ind, indObj, indTrf));

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
 * Space Transformations: Creates an object mirror with respect to a plane.
 * The object mirror possesses the same geometry as obj, but its magnetization and/or current densities 
 * are modified in such a way that the magnetic field produced by the obj and its mirror in the plane 
 * of mirroring is PERPENDICULAR to this plane.
 ***************************************************************************/
static PyObject* radia_TrfZerPara(PyObject* self, PyObject* args)
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
static PyObject* radia_TrfZerPerp(PyObject* self, PyObject* args)
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
 * Applies magnetic material to a 3D object.
 ***************************************************************************/
static PyObject* radia_MatApl(PyObject* self, PyObject* args)
{
	PyObject *oResInd=0;
	try
	{
		int indObj = 0, indMat = 0;
		if(!PyArg_ParseTuple(args, "ii:MatApl", &indObj, &indMat)) throw CombErStr(strEr_BadFuncArg, ": MatApl");
		if((indObj == 0) || (indMat == 0)) throw CombErStr(strEr_BadFuncArg, ": MatApl");

		int indRes = 0;
		g_pyParse.ProcRes(RadMatApl(&indRes, indObj, indMat));

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
static PyObject* radia_MatStd(PyObject* self, PyObject* args)
{
	PyObject *oMatId=0, *oResInd=0;
	try
	{
		double absM = -1; //OC02112019
		//double absM = 0;
		if(!PyArg_ParseTuple(args, "O|d:MatStd", &oMatId, &absM)) throw CombErStr(strEr_BadFuncArg, ": MatStd"); //OC02112019
		//if(!PyArg_ParseTuple(args, "Od:MatStd", &oMatId, &absM)) throw CombErStr(strEr_BadFuncArg, ": MatStd");
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
 * Magnetic Materials: a linear anisotropic magnetic material.
 ***************************************************************************/
static PyObject* radia_MatLin(PyObject* self, PyObject* args)
{
	PyObject *oKsi=0, *oMr=0, *oResInd=0;
	try
	{
		if(!PyArg_ParseTuple(args, "OO:MatLin", &oKsi, &oMr)) throw CombErStr(strEr_BadFuncArg, ": MatLin");
		if((oKsi == 0) || (oMr == 0)) throw CombErStr(strEr_BadFuncArg, ": MatLin");

		double arKsi[2];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oKsi, 'd', arKsi, 2, CombErStr(strEr_BadFuncArg, ": MatLin, incorrect definition of susceptibility parameters"));

		int nMr=0;
		double arMr[3];
		if(PyNumber_Check(oMr))
		{
			*arMr = PyFloat_AsDouble(oMr);
			nMr = 1;
		}
		else
		{
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oMr, 'd', arMr, 3, CombErStr(strEr_BadFuncArg, ": MatLin, incorrect definition of remanent magnetization"));
			nMr = 3;
		}

		int indRes=0;
		g_pyParse.ProcRes(RadMatLin(&indRes, arKsi, arMr, nMr));

		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch (const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Materials: a nonlinear isotropic magnetic material with the magnetization magnitude equal M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3), where H is the magnitude of the magnetic field strength vector (in Tesla).
 ***************************************************************************/
static PyObject* radia_MatSatIsoFrm(PyObject* self, PyObject* args)
{
	PyObject *oPair1=0, *oPair2=0, *oPair3=0, *oResInd=0;
	double *arKsiMs1=0, *arKsiMs2=0, *arKsiMs3=0;
	try
	{
		if(!PyArg_ParseTuple(args, "O|OO:MatSatIsoFrm", &oPair1, &oPair2, &oPair3)) throw CombErStr(strEr_BadFuncArg, ": MatSatIsoFrm");
		if(oPair1 == 0) throw CombErStr(strEr_BadFuncArg, ": MatSatIsoFrm");
		
		int nElem = 0;
		if((!CPyParse::CopyPyNestedListElemsToNumAr(oPair1, 'd', arKsiMs1, nElem)) || (nElem != 2)) throw CombErStr(strEr_BadFuncArg, ": MatSatIsoFrm, pair of Ksi and Ms is expected");

		if(oPair2 != 0)
		{
			if((!CPyParse::CopyPyNestedListElemsToNumAr(oPair2, 'd', arKsiMs2, nElem)) || (nElem != 2)) throw CombErStr(strEr_BadFuncArg, ": MatSatIsoFrm, pair of Ksi and Ms is expected");
		}
		if(oPair3 != 0)
		{
			if((!CPyParse::CopyPyNestedListElemsToNumAr(oPair3, 'd', arKsiMs3, nElem)) || (nElem != 2)) throw CombErStr(strEr_BadFuncArg, ": MatSatIsoFrm, pair of Ksi and Ms is expected");
		}

		int indRes = 0;
		g_pyParse.ProcRes(RadMatSatIsoFrm(&indRes, arKsiMs1, arKsiMs2, arKsiMs3));

		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arKsiMs1 != 0) delete[] arKsiMs1;
	if(arKsiMs2 != 0) delete[] arKsiMs2;
	if(arKsiMs3 != 0) delete[] arKsiMs3;
	return oResInd;
}

/************************************************************************//**
 * Magnetic Materials: a nonlinear isotropic magnetic material with the M versus H curve defined by the list of pairs corresponding values of H and M (H1,M1,H2,M2,...)
 ***************************************************************************/
static PyObject* radia_MatSatIsoTab(PyObject* self, PyObject* args)
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
		int nDataP = (int)round(0.5 * nDataTot);

		int indRes = 0;
		g_pyParse.ProcRes(RadMatSatIsoTab(&indRes, arMatData, nDataP));

		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch (const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arMatData != 0) delete[] arMatData;
	return oResInd;
}

/************************************************************************//**
 * Magnetic Materials: laminated nonlinear anisotropic magnetic material with packing factor p and the lamination planes perpendicular to the vector N. The magnetization magnitude vs magnetic field strength for the corresponding isotropic material is defined by the formula M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3), where H is the magnitude of the magnetic field strength vector (in Tesla); ksi1, ms1, ksi2, ms2, ksi3, ms3 constants are given by parameters KsiMs1, KsiMs2, KsiMs3.
 ***************************************************************************/
static PyObject* radia_MatSatLamFrm(PyObject* self, PyObject* args)
{
	PyObject *oPair1=0, *o2=0, *o3=0, *o4=0, *o5=0, *oResInd=0;
	double *arKsiMs1=0, *arKsiMs2=0, *arKsiMs3=0, *arN=0;
	try
	{
		if(!PyArg_ParseTuple(args, "OO|OOO:MatSatLamFrm", &oPair1, &o2, &o3, &o4, &o5)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm");
		if(oPair1 == 0) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm");

		int nElem=0;
		if((!CPyParse::CopyPyNestedListElemsToNumAr(oPair1, 'd', arKsiMs1, nElem)) || (nElem != 2)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm, pair of Ksi and Ms is expected");

		double p=0;
		if(o2 != 0)
		{//Can be another pair or staking factor 
			if(PyNumber_Check(o2)) p = PyFloat_AsDouble(o2);
			else
			{
				if((!CPyParse::CopyPyNestedListElemsToNumAr(o2, 'd', arKsiMs2, nElem)) || (nElem != 2)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm, pair of Ksi and Ms is expected");
			}
		}
		if(o3 != 0)
		{
			if(PyNumber_Check(o2)) p = PyFloat_AsDouble(o3);
			else
			{
				if(p == 0)
				{
					if((!CPyParse::CopyPyNestedListElemsToNumAr(o3, 'd', arKsiMs2, nElem)) || (nElem != 2)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm, pair of Ksi and Ms is expected");
				}
				else
				{
					if((!CPyParse::CopyPyNestedListElemsToNumAr(o3, 'd', arN, nElem)) || (nElem != 3)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm, 3D vector coordinate list / array is expected");
				}
			}
		}
		if(o4 != 0)
		{
			if(p == 0)
			{
				if(PyNumber_Check(o2)) p = PyFloat_AsDouble(o4);
			}
			else
			{
				if((!CPyParse::CopyPyNestedListElemsToNumAr(o4, 'd', arN, nElem)) || (nElem != 3)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm, 3D vector coordinate list / array is expected");
			}
		}
		if(o5 != 0)
		{
			if((!CPyParse::CopyPyNestedListElemsToNumAr(o5, 'd', arN, nElem)) || (nElem != 3)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm, 3D vector coordinate list / array is expected");
		}

		if(p == 0) throw CombErStr(strEr_BadFuncArg, ": MatSatLamFrm, positive lamination stacking factor is not defined");

		int indRes=0;
		g_pyParse.ProcRes(RadMatSatLamFrm(&indRes, arKsiMs1, arKsiMs2, arKsiMs3, p, arN));

		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch (const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arKsiMs1 != 0) delete[] arKsiMs1;
	if(arKsiMs2 != 0) delete[] arKsiMs2;
	if(arKsiMs3 != 0) delete[] arKsiMs3;
	if(arN != 0) delete[] arN;
	return oResInd;
}

/************************************************************************//**
 * Magnetic Materials: laminated nonlinear anisotropic magnetic material with packing factor p and the lamination planes perpendicular to the vector N. The magnetization magnitude vs magnetic field strength for the corresponding isotropic material is defined by pairs of values H, M in Tesla.
 ***************************************************************************/
static PyObject* radia_MatSatLamTab(PyObject* self, PyObject* args)
{
	PyObject *oMatData=0, *oN=0, *oResInd=0;
	double *arMatData=0, *arN=0;
	try
	{
		double p=0;
		if(!PyArg_ParseTuple(args, "Od|O:MatSatLamTab", &oMatData, &p, oN)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamTab");
		if((oMatData == 0) || (p == 0)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamTab");

		int nDataTot=0;
		char resP = CPyParse::CopyPyNestedListElemsToNumAr(oMatData, 'd', arMatData, nDataTot);
		if(resP == 0) throw CombErStr(strEr_BadFuncArg, ": MatSatLamTab");
		int nDataP = (int)round(0.5*nDataTot);

		if(oN != 0)
		{
			int nElem=0;
			if((!CPyParse::CopyPyNestedListElemsToNumAr(oN, 'd', arN, nElem)) || (nElem != 3)) throw CombErStr(strEr_BadFuncArg, ": MatSatLamTab, 3D vector coordinate list / array is expected");
		}

		int indRes=0;
		g_pyParse.ProcRes(RadMatSatLamTab(&indRes, arMatData, nDataP, p, arN));

		oResInd = Py_BuildValue("i", indRes);
		Py_XINCREF(oResInd); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arMatData != 0) delete[] arMatData;
	if(arN != 0) delete[] arN;
	return oResInd;
}

/************************************************************************//**
 * Magnetic Materials: a nonlinear anisotropic magnetic material. The magnetization vector component parallel to the easy axis is computed either by the formula: ms1*tanh(ksi1*(hpa-hc1)/ms1) + ms2*tanh(ksi2*(hpa-hc2)/ms2) + ms3*tanh(ksi3*(hpa-hc3)/ms3) + ksi0*(hpa-hc0), where hpa is the field strength vector component parallel to the easy axis, or by ksi0*hpa. The magnetization vector component perpendicular to the easy axis is computed either by the formula: ms1*tanh(ksi1*hpe/ms1) + ms2*tanh(ksi2*hpe/ms2) + ms3*tanh(ksi3*hpe/ms3) + ksi0*hpe, where hpe is the field strength vector component perpendicular to the easy axis, or by ksi0*hpe. At least one of the magnetization components should non-linearly depend on the field strength. The direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is later applied.
 ***************************************************************************/
static PyObject* radia_MatSatAniso(PyObject* self, PyObject* args)
{
	PyObject *oPar=0, *oPer=0, *oResInd=0;
	try
	{
		if(!PyArg_ParseTuple(args, "OO:MatSatAniso", &oPar, &oPer)) throw CombErStr(strEr_BadFuncArg, ": MatSatAniso");
		if((oPar == 0) || (oPer == 0)) throw CombErStr(strEr_BadFuncArg, ": MatSatAniso");

		double arPar[11], arPer[11];
		int nPar=0, nPer=0;
		if(PyNumber_Check(oPar))
		{
			*arPar = PyFloat_AsDouble(oPar);
			nPar = 1;
		}
		else
		{
			nPar = 11;
			double *p = arPar;
			if(!CPyParse::CopyPyNestedListElemsToNumAr(oPar, 'd', p, nPar)) throw CombErStr(strEr_BadFuncArg, ": MatSatAniso, incorrect definition of magnetic material constants in the direction parallel to the easy axis");
		}
		if(PyNumber_Check(oPer))
		{
			*arPer = PyFloat_AsDouble(oPer);
			nPer = 1;
		}
		else
		{
			nPer = 11;
			double *p = arPer;
			if(!CPyParse::CopyPyNestedListElemsToNumAr(oPer, 'd', p, nPer)) throw CombErStr(strEr_BadFuncArg, ": MatSatAniso, incorrect definition of magnetic material constants in the direction perpendicular to the easy axis");
		}

		int indRes=0;
		g_pyParse.ProcRes(RadMatSatAniso(&indRes, arPar, nPar, arPer, nPer));

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
static PyObject* radia_MatMvsH(PyObject* self, PyObject* args)
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
 * Magnetic Field Calculation Methods: Builds interaction matrix for an object.
 ***************************************************************************/
static PyObject* radia_RlxPre(PyObject* self, PyObject* args)
{
	PyObject *oResInd=0;
	try
	{
		int indObj = 0, indSrc = 0;
		//if(!PyArg_ParseTuple(args, "ii:RlxPre", &indObj, &indSrc)) throw CombErStr(strEr_BadFuncArg, ": RlxPre");
		if(!PyArg_ParseTuple(args, "i|i:RlxPre", &indObj, &indSrc)) throw CombErStr(strEr_BadFuncArg, ": RlxPre"); //AB22112019
		if(indObj == 0) throw CombErStr(strEr_BadFuncArg, ": RlxPre");

		int ind = 0;
		g_pyParse.ProcRes(RadRlxPre(&ind, indObj, indSrc));

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
 * Magnetic Field Calculation Methods: Executes manual relaxation procedure for the interaction matrix intrc.
 ***************************************************************************/
static PyObject* radia_RlxMan(PyObject* self, PyObject* args)
{
	PyObject *oResInd=0;
	try
	{
		int ind = 0, meth = 0, numIt = 0;
		double rlxPar = 0;
		if(!PyArg_ParseTuple(args, "iiid:RlxMan", &ind, &meth, &numIt, &rlxPar)) throw CombErStr(strEr_BadFuncArg, ": RlxMan");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": RlxMan");

		double arRes[12];
		int lenRes = 0;
		g_pyParse.ProcRes(RadRlxMan(arRes, &lenRes, ind, meth, numIt, rlxPar));

		if(lenRes == 1) oResInd = Py_BuildValue("d", *arRes);
		else if(lenRes > 1) oResInd = CPyParse::SetDataListOfLists(arRes, lenRes, 1);
		if(oResInd) Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Executes automatic relaxation procedure for the interaction matrix intrc.
 ***************************************************************************/
static PyObject* radia_RlxAuto(PyObject* self, PyObject* args)
{
	PyObject *oOpt=0, *oResInd=0;
	try
	{
		int ind=0, numIt=0, meth=4; //OC30122019
		//int ind=0, numIt=0, meth=0;
		double prec = 0;
		if(!PyArg_ParseTuple(args, "idi|iO:RlxAuto", &ind, &prec, &numIt, &meth, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": RlxAuto");
		//if(!PyArg_ParseTuple(args, "idii|O:RlxAuto", &ind, &prec, &numIt, &meth, &oOpt)) throw CombErStr(strEr_BadFuncArg, ": RlxAuto");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": RlxAuto");

		char sOpt[1024]; *sOpt = '\0';
		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		double arRes[12];
		int lenRes = 0;
		g_pyParse.ProcRes(RadRlxAuto(arRes, &lenRes, ind, prec, numIt, meth, sOpt));

		if(lenRes == 1) oResInd = Py_BuildValue("d", *arRes);
		else if(lenRes > 1) oResInd = CPyParse::SetDataListOfLists(arRes, lenRes, 1);
		if(oResInd) Py_XINCREF(oResInd);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oResInd;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Updates external field data for the relaxation (to take into account e.g. modification of currents in coils, if any) without rebuilding the interaction matrix.
 ***************************************************************************/
static PyObject* radia_RlxUpdSrc(PyObject* self, PyObject* args)
{
	PyObject *oResInd = 0;

	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "i:RlxUpdSrc", &ind)) throw CombErStr(strEr_BadFuncArg, ": RlxUpdSrc");
		if(ind == 0) throw CombErStr(strEr_BadFuncArg, ": RlxUpdSrc");

		g_pyParse.ProcRes(RadRlxUpdSrc(ind));

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
 * Magnetic Field Calculation Methods: Builds an interaction matrix and performs a relaxation procedure.
 ***************************************************************************/
static PyObject* radia_Solve(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	try
	{
		int ind=0, numIt=1000, meth=4; //OC02112019
		//int ind=0, numIt=1000, meth=0;
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
static PyObject* radia_Fld(PyObject* self, PyObject* args)
{
	PyObject *oCmpnId=0, *oP=0, *oResB=0;
	double *arCoord=0, *arB=0;
	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "iOO:Fld", &ind, &oCmpnId, &oP)) throw CombErStr(strEr_BadFuncArg, ": Fld");
		if((ind == 0) || (oCmpnId == 0) || (oP == 0)) throw CombErStr(strEr_BadFuncArg, ": Fld");

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

		if(nB <= 0) { nB = 1; *arB = 0;} //OC10012020 (to ensure some dummy return, e.g. in the case of MPI, rank > 0)

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
 * Magnetic Field Calculation Methods: Computes magnetic field created by object obj in np equidistant points along a line segment from P1 to P2; the field component is specified by the id input variable.
 ***************************************************************************/
static PyObject* radia_FldLst(PyObject* self, PyObject* args)
{
	PyObject *oCmpnId=0, *oP1=0, *oP2=0, *oOpt=0, *oResB=0;
	//PyObject *oCmpnId=0, *oP1=0, *oP2=0, *oOpt=0, *oResB=0;
	double *arB=0;
	try
	{
		int ind=0, nP=0;
		double start=0.;
		if(!PyArg_ParseTuple(args, "iOOOi|Od:FldLst", &ind, &oCmpnId, &oP1, &oP2, &nP, &oOpt, &start)) throw CombErStr(strEr_BadFuncArg, ": FldLst");
		//if(!PyArg_ParseTuple(args, "iOOOO|OO:FldLst", &ind, &oCmpnId, &oP1, &oP2, nP, &oOpt, &start)) throw CombErStr(strEr_BadFuncArg, ": FldLst");
		if((ind <= 0) || (oCmpnId == 0) || (oP1 == 0) || (oP2 == 0) || (nP <= 0)) throw CombErStr(strEr_BadFuncArg, ": FldLst");

		char sCmpnId[256];
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 256);

		double arP1[3], arP2[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP1, 'd', arP1, 3, CombErStr(strEr_BadFuncArg, ": FldLst, incorrect definition of first point defining line segment along which field should be calculated"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP2, 'd', arP2, 3, CombErStr(strEr_BadFuncArg, ": FldLst, incorrect definition of second point defining line segment along which field should be calculated"));

		char sOpt[1024]; *sOpt = '\0';
		if(oOpt != 0) CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		const int maxNumFldCmpn = 14;
		arB = new double[maxNumFldCmpn*nP];
		int nB = 0;
		g_pyParse.ProcRes(RadFldLst(arB, &nB, ind, sCmpnId, arP1, arP2, nP, sOpt, start));

		if(nB == 1) oResB = Py_BuildValue("d", *arB);
		else if(nB > 1) oResB = CPyParse::SetDataListOfLists(arB, nB, nP);
		if(oResB) Py_XINCREF(oResB);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arB != 0) delete[] arB;
	return oResB;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes magnetic field integral produced by the object obj along a straight line specified by points P1 and P2.
 * Depending on the InfOrFin variable value, the integral is infinite ("inf") or finite ("fin"), from P1 to P2; the field integral component is specified by the id input variable. The unit is T*mm.
 ***************************************************************************/
static PyObject* radia_FldInt(PyObject* self, PyObject* args)
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

/***************************************************************************
 * Magnetic Field Calculation Methods: Computes transverse coordinates and their derivatives (angles) of a relativistic charged trajectory in the magnetic field produced by object obj, 
 * using a Runge-Kutta integration. The charge of the particle is that of electron. All positions are in millimeters and angles in radians.
 ***************************************************************************/
static PyObject* radia_FldPtcTrj(PyObject* self, PyObject* args) 
{
	PyObject *oInitCond=0, *oLongLim=0, *oResTrj=0;
	double *arTrj=0; //OC

	try
	{
		int ind=0, Np;
		double E;

		if(!PyArg_ParseTuple(args, "idOOi:FldPtcTrj", &ind, &E, &oInitCond, &oLongLim, &Np)) throw CombErStr(strEr_BadFuncArg, ": FldPtcTrj");
		if((oInitCond == 0) || (oLongLim == 0)) throw CombErStr(strEr_BadFuncArg, ": FldPtcTrj");

		double arInitCond[4], arLongLim[2];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oInitCond, 'd', arInitCond, 4, CombErStr(strEr_BadFuncArg, ": FldPtcTrj, incorrect definition of initial conditions"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oLongLim, 'd', arLongLim, 2, CombErStr(strEr_BadFuncArg, ": FldPtcTrj, incorrect definition of Y range"));

		//double *arTrj = new double[Np*5];
		arTrj = new double[Np*5]; //OC
		int nTrjPts = 0;
		g_pyParse.ProcRes(RadFldPtcTrj(arTrj, &nTrjPts, ind, E, arInitCond, arLongLim, Np));

		if(nTrjPts == 1) oResTrj = Py_BuildValue("d", *arTrj);
		else if (nTrjPts > 1) oResTrj = CPyParse::SetDataListOfLists(arTrj, nTrjPts, Np);
		if(oResTrj) Py_XINCREF(oResTrj);

		//if (arTrj != 0) delete[] arTrj; //OC (commented-out)
	}
	catch (const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	if(arTrj != 0) delete[] arTrj; //OC
	return oResTrj;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes potential energy (in Joule) of the object objdst in the field created by the object objsrc. If SbdPar = 0, the function performes the computation based on absolute accuracy value for the energy (by default 10 Joule; can be modified by the function FldCmpPrc). Otherwise, the computation is performed based on the destination object subdivision numbers (kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]).
 ***************************************************************************/
static PyObject* radia_FldEnr(PyObject* self, PyObject* args)
{
	PyObject *oSbdPar=0, *oResE=0;
	try
	{
		int indDst=0, indSrc=0;
		if(!PyArg_ParseTuple(args, "ii|O:FldEnr", &indDst, &indSrc, &oSbdPar)) throw CombErStr(strEr_BadFuncArg, ": FldEnr");
		if((indDst <= 0) || (indSrc <= 0)) throw CombErStr(strEr_BadFuncArg, ": FldEnr");

		int arSbdPar[] = {1,1,1};
		if(oSbdPar != 0) CPyParse::CopyPyListElemsToNumArrayKnownLen(oSbdPar, 'i', arSbdPar, 3, CombErStr(strEr_BadFuncArg, ": FldEnr, incorrect definition of subdivision numbers"));

		double resE=0;
		g_pyParse.ProcRes(RadFldEnr(&resE, indDst, indSrc, arSbdPar));

		oResE = Py_BuildValue("d", resE);
		Py_XINCREF(oResE); //?
	}
	catch (const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oResE;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes force (in Newton) acting on the object objdst in the field produced by the object objsrc. If SbdPar = 0, the function performes the computation based on absolute accuracy value for the force (by default 10 Newton; can be modified by the function FldCmpPrc). Otherwise, the computation is performed based on the destination object subdivision numbers {kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]}.
 ***************************************************************************/
static PyObject* radia_FldEnrFrc(PyObject* self, PyObject* args)
{
	PyObject *oCmpnId=0, *oSbdPar=0, *oResF=0;
	try
	{
		int indDst=0, indSrc=0;
		if(!PyArg_ParseTuple(args, "iiO|O:FldEnrFrc", &indDst, &indSrc, &oCmpnId, &oSbdPar)) throw CombErStr(strEr_BadFuncArg, ": FldEnrFrc");
		if((indDst <= 0) || (indSrc <= 0) || (oCmpnId == 0)) throw CombErStr(strEr_BadFuncArg, ": FldEnrFrc");

		char sCmpnId[256];
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 256);

		int arSbdPar[] = {1,1,1};
		if(oSbdPar != 0) CPyParse::CopyPyListElemsToNumArrayKnownLen(oSbdPar, 'i', arSbdPar, 3, CombErStr(strEr_BadFuncArg, ": FldEnrFrc, incorrect definition of subdivision numbers"));

		double arF[9];
		int nF=0;
		g_pyParse.ProcRes(RadFldEnrFrc(arF, &nF, indDst, indSrc, sCmpnId, arSbdPar));

		if(nF == 1) oResF = Py_BuildValue("d", *arF);
		else if(nF > 1) oResF = CPyParse::SetDataListOfLists(arF, nF, 1);
		if(oResF) Py_XINCREF(oResF);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oResF;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes torque (in Newton*mm) with respect to point P, acting on the object objdst in the field produced by the object objsrc. If SbdPar = 0, the function performes the computation based on absolute accuracy value for the torque (by default 10 Newton*mm; can be modified by the function FldCmpPrc). Otherwise, the computation is performed based on the destination object subdivision numbers {kx=SbdPar[0],ky=SbdPar[1],kz=SbdPar[2]}.
 ***************************************************************************/
static PyObject* radia_FldEnrTrq(PyObject* self, PyObject* args)
{
	PyObject *oCmpnId=0, *oP=0, *oSbdPar=0, *oResT=0;
	try
	{
		int indDst=0, indSrc=0;
		if(!PyArg_ParseTuple(args, "iiOO|O:FldEnrTrq", &indDst, &indSrc, &oCmpnId, &oP, &oSbdPar)) throw CombErStr(strEr_BadFuncArg, ": FldEnrTrq");
		if((indDst <= 0) || (indSrc <= 0) || (oCmpnId == 0) || (oP == 0)) throw CombErStr(strEr_BadFuncArg, ": FldEnrTrq"); //OC14042020
		//if((indDst <= 0) || (indSrc <= 0) || (oCmpnId == 0) || (oP)) throw CombErStr(strEr_BadFuncArg, ": FldEnrTrq");

		char sCmpnId[256];
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 255);

		int arSbdPar[]={1,1,1};
		if(oSbdPar != 0) CPyParse::CopyPyListElemsToNumArrayKnownLen(oSbdPar, 'i', arSbdPar, 3, CombErStr(strEr_BadFuncArg, ": FldEnrTrq, incorrect definition of subdivision numbers"));

		double arP[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": FldEnrTrq, incorrect definition of point the torque should be calculated with respect to"));

		double arT[9];
		int nT=0;
		g_pyParse.ProcRes(RadFldEnrTrq(arT, &nT, indDst, indSrc, sCmpnId, arP, arSbdPar));

		if(nT == 1) oResT = Py_BuildValue("d", *arT);
		else if(nT > 1) oResT = CPyParse::SetDataListOfLists(arT, nT, 1);
		if(oResT) Py_XINCREF(oResT);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oResT;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes force of the field produced by the object obj onto a shape defined by shape. shape can be the result of RadObjRecMag (parallelepiped) or RadFldFrcShpRtg (rectangular surface). This function uses the algorithm based on Maxwell tensor, which may not always provide high efficiency. We suggest to use the function RadFldEnrFrc instead of this function.
 ***************************************************************************/
static PyObject* radia_FldFrc(PyObject* self, PyObject* args)
{
	PyObject *oResF=0;
	try
	{
		int indObj=0, indShp=0;
		if(!PyArg_ParseTuple(args, "ii:FldFrc", &indObj, &indShp)) throw CombErStr(strEr_BadFuncArg, ": FldFrc");
		if((indObj <= 0) || (indShp <= 0)) throw CombErStr(strEr_BadFuncArg, ": FldFrc");

		double arF[9];
		int nF=0;
		g_pyParse.ProcRes(RadFldFrc(arF, &nF, indObj, indShp));

		if(nF == 1) oResF = Py_BuildValue("d", *arF);
		else if(nF > 1) oResF = CPyParse::SetDataListOfLists(arF, nF, 1);
		if(oResF) Py_XINCREF(oResF);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oResF;
}

/************************************************************************//**
 * Auxiliary Shape Object: Creates a rectangle with central point P and dimensions W (to be used for force computation via Maxwell tensor)
 ***************************************************************************/
static PyObject* radia_FldFrcShpRtg(PyObject* self, PyObject* args)
{
	PyObject *oP=0, *oW=0, *oResInd=0;
	try
	{
		if(!PyArg_ParseTuple(args, "OO:FldFrcShpRtg", &oP, &oW)) throw CombErStr(strEr_BadFuncArg, ": FldFrcShpRtg");
		if((oP == 0) || (oW == 0)) throw CombErStr(strEr_BadFuncArg, ": FldFrcShpRtg");

		double arP[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP, 'd', arP, 3, CombErStr(strEr_BadFuncArg, ": FldFrcShpRtg, incorrect definition of center point"));

		double arW[2];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oW, 'd', arW, 2, CombErStr(strEr_BadFuncArg, ": FldFrcShpRtg, incorrect definition of dimensions"));

		int ind=0;
		g_pyParse.ProcRes(RadFldFrcShpRtg(&ind, arP, arW));

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
 * Magnetic Field Calculation Methods: Computes "focusing potential" for trajectory of relativistic charged particle in magnetic field produced by the object obj. The integration is made from P1 to P2 with np equidistant points.
 ***************************************************************************/
static PyObject* radia_FldFocPot(PyObject* self, PyObject* args)
{
	PyObject *oP1=0, *oP2=0, *oResFP=0;
	try
	{
		int indObj=0, np=0;
		if(!PyArg_ParseTuple(args, "iOOi:FldFocPot", &indObj, &oP1, &oP2)) throw CombErStr(strEr_BadFuncArg, ": FldFocPot");
		if((indObj <= 0) || (np <= 0)) throw CombErStr(strEr_BadFuncArg, ": FldFocPot");

		double arP1[3], arP2[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP1, 'd', arP1, 3, CombErStr(strEr_BadFuncArg, ": FldFocPot, incorrect definition of first end point"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP2, 'd', arP2, 3, CombErStr(strEr_BadFuncArg, ": FldFocPot, incorrect definition of second end point"));

		double resFP=0;
		g_pyParse.ProcRes(RadFldFocPot(&resFP, indObj, arP1, arP2, np));

		oResFP = Py_BuildValue("d", resFP);
		Py_XINCREF(oResFP); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oResFP;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes matrices of 2nd order kicks for trajectory of relativistic charged particle in periodic magnetic field produced by the object obj. The computed kick matrices can be used in charged particle tracking codes.  The longitudinal integration along one period starts at point P1 and is done along direction pointed by vector Ns; one direction of the transverse grid is pointed by vector Ntr, the other transverse direction is given by vector product of Ntr and Ns.
 ***************************************************************************/
static PyObject* radia_FldFocKickPer(PyObject* self, PyObject* args)
{
	PyObject *oP1=0, *oNs=0, *oN1=0, *oCom=0, *oPrecPar=0, *oUnits=0, *oFrm=0, *oRes=0;
	double *arM1=0, *arM2=0, *arIntBtrE2=0, *arArg1=0, *arArg2=0;
	char *sOut=0;
	try
	{
		int ind=0, nPer=0, np1=0, np2=0;
		double per=0., r1=0., r2=0., E=1.;
		if(!PyArg_ParseTuple(args, "iOOdiOdidi|OOOdO:FldFocKickPer", &ind, &oP1, &oNs, &per, &nPer, &oN1, &r1, &np1, &r2, &np2, &oCom, &oPrecPar, &oUnits, &E, &oFrm)) throw CombErStr(strEr_BadFuncArg, ": FldFocKickPer");
		if((oP1 == 0) || (oNs == 0) || (oN1 == 0) || (np1 <= 0) || (np2 <= 0)) throw CombErStr(strEr_BadFuncArg, ": FldFocKickPer");

		double arP1[3], arNs[3], arN1[3];
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP1, 'd', arP1, 3, CombErStr(strEr_BadFuncArg, ": FldFocKickPer, incorrect definition of the integration start point"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oNs, 'd', arNs, 3, CombErStr(strEr_BadFuncArg, ": FldFocKickPer, incorrect definition of a vector defining the longitudinal direction for the periodic magnetic field"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oN1, 'd', arN1, 3, CombErStr(strEr_BadFuncArg, ": FldFocKickPer, incorrect definition of a vector defining first transverse direction"));

		char sCom[1025];
		*sCom = '\0';
		if(oCom != 0) CPyParse::CopyPyStringToC(oCom, sCom, 1024);

		int nh=1, nps=8;
		double d1=0., d2=0.;
		if(oPrecPar != 0)
		{
			double arPrecPar[4];
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oPrecPar, 'd', arPrecPar, 4, CombErStr(strEr_BadFuncArg, ": FldFocKickPer, incorrect definition of the integration precision parameters"));
			nh = (int)arPrecPar[0]; 
			nps = (int)arPrecPar[1];
			d1 = arPrecPar[2];
			d2 = arPrecPar[3];
		}

		char sUnits[256];
		strcpy(sUnits, "T2m2\0");
		if(oUnits != 0) CPyParse::CopyPyStringToC(oUnits, sUnits, 255);

		char sFrm[256];
		strcpy(sFrm, "fix\0");
		if(oFrm != 0) CPyParse::CopyPyStringToC(oFrm, sFrm, 255);

		long np1_np2 = np1*np2;
		arM1 = new double[np1_np2];
		arM2 = new double[np1_np2];
		arIntBtrE2 = new double[np1_np2];
		arArg1 = new double[np1];
		arArg2 = new double[np2];
		int sLen=0;
		g_pyParse.ProcRes(RadFldFocKickPer(arM1, arM2, arIntBtrE2, arArg1, arArg2, &sLen, ind, arP1, arNs, per, nPer, nps, arN1, r1, np1, d1, r2, np2, d2, nh, sCom, sUnits, E, sFrm)); //OC03112019
		//g_pyParse.ProcRes(RadFldFocKickPer(arM1, arM2, arIntBtrE2, arArg1, arArg2, &sLen, ind, arP1, arNs, per, nPer, nps, arN1, r1, np1, d1, r2, np2, d2, nh, sCom));

		sOut = new char[sLen + 1]; //OC19052020
		//sOut = new char[sLen];
		g_pyParse.ProcRes(RadFldFocKickPerFormStr(sOut, arM1, arM2, arIntBtrE2, arArg1, arArg2, np1, np2, per, nPer, sCom));

		oRes = PyTuple_New(6);
		PyObject *oM1 = CPyParse::SetDataListOfLists(arM1, np1_np2, (long)np2); //OC19052020 //check dims; make it a flat array?
		//PyObject *oM1 = CPyParse::SetDataListOfLists(arM1, np1, np2); //check dims; make it a flat array?
		PyTuple_SET_ITEM(oRes, 0, oM1);
		PyObject *oM2 = CPyParse::SetDataListOfLists(arM2, np1_np2, (long)np2); //OC19052020 //check dims; make it a flat array?
		//PyObject *oM2 = CPyParse::SetDataListOfLists(arM2, np1, np2); //check dims; make it a flat array?
		PyTuple_SET_ITEM(oRes, 1, oM2);
		PyObject *oIntBtrE2 = CPyParse::SetDataListOfLists(arIntBtrE2, np1_np2, (long)np2); //OC19052020 //check dims; make it a flat array?
		//PyObject *oIntBtrE2 = CPyParse::SetDataListOfLists(arIntBtrE2, np1, np2); //check dims; make it a flat array?
		PyTuple_SET_ITEM(oRes, 2, oIntBtrE2);
		PyObject *oArg1 = CPyParse::SetDataListOfLists(arArg1, np1, 1);
		PyTuple_SET_ITEM(oRes, 3, oArg1);
		PyObject *oArg2 = CPyParse::SetDataListOfLists(arArg2, np2, 1);
		PyTuple_SET_ITEM(oRes, 4, oArg2);
		PyObject *oOutStr = Py_BuildValue("s", sOut); //to check
		PyTuple_SET_ITEM(oRes, 5, oOutStr);

		Py_XINCREF(oRes); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}

	if(arM1 != 0) delete[] arM1;
	if(arM2 != 0) delete[] arM2;
	if(arIntBtrE2 != 0) delete[] arIntBtrE2;
	if(arArg1 != 0) delete[] arArg1;
	if(arArg2 != 0) delete[] arArg2;
	if(sOut != 0) delete[] sOut;
	return oRes;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Sets general absolute accuracy levels for computation of magnetic field induction (prcB), vector potential (prcA), induction integrals along straight line (prcBint), field force (prcFrc), relativistic particle trajectory coordinates (prcTrjCrd) and angles (prcTrjAng)
 ***************************************************************************/
static PyObject* radia_FldCmpCrt(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	try
	{
		double prcB=0, prcA=0, prcBInt=0, prcFrc=0, prcTrjCrd=0, prcTrjAng=0;
		if(!PyArg_ParseTuple(args, "dddddd:FldCmpCrt", &prcB, &prcA, &prcBInt, &prcFrc, &prcTrjCrd, &prcTrjAng)) throw CombErStr(strEr_BadFuncArg, ": FldCmpCrt");
		if((prcB <= 0) || (prcA <= 0) || (prcBInt <= 0) || (prcFrc <= 0) || (prcTrjCrd <= 0) || (prcTrjAng <= 0)) throw CombErStr(strEr_BadFuncArg, ": FldCmpCrt");

		int dummy;
		g_pyParse.ProcRes(RadFldCmpCrt(&dummy, prcB, prcA, prcBInt, prcFrc, prcTrjCrd, prcTrjAng));

		oRes = Py_BuildValue("i", dummy);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oRes;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Sets general absolute accuracy levels for computation of magnetic field induction (PrcB), vector potential (PrcA), induction integral along straight line (PrcBInt), field force (PrcForce), torque (PrcTorque), energy (PrcEnergy); relativistic charged particle trajectory coordinates (PrcCoord) and angles (PrcAngle). The function works according to the mechanism of string options. The name(s) of the option(s) should be: PrcB, PrcA, PrcBInt, PrcForce, PrcTorque, PrcEnergy, PrcCoord, PrcAngle.
 ***************************************************************************/
static PyObject* radia_FldCmpPrc(PyObject* self, PyObject* args)
{
	PyObject *oOpt=0, *oRes=0;
	try
	{
		if(!PyArg_ParseTuple(args, "O:FldCmpPrc", &oOpt)) throw CombErStr(strEr_BadFuncArg, ": FldCmpPrc");
		if(oOpt == 0) throw CombErStr(strEr_BadFuncArg, ": FldCmpPrc");

		char sOpt[1024]; *sOpt = '\0';
		CPyParse::CopyPyStringToC(oOpt, sOpt, 1024);

		int dummy;
		g_pyParse.ProcRes(RadFldCmpPrc(&dummy, sOpt));

		oRes = Py_BuildValue("i", dummy);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oRes;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Shows the physical units currently in use.
 ***************************************************************************/
static PyObject* radia_FldUnits(PyObject* self, PyObject* args)
{
	PyObject* oResUnits = 0;
	try
	{
		char sUnits[2048];
		g_pyParse.ProcRes(RadFldUnits(sUnits));

		oResUnits = Py_BuildValue("s", sUnits); //to check
		Py_XINCREF(oResUnits); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oResUnits;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Switches on or off the randomization of all the length values. The randomization magnitude can be set by the function radFldLenTol.
 ***************************************************************************/
static PyObject* radia_FldLenRndSw(PyObject* self, PyObject* args)
{
	PyObject *oOnOff=0, *oRes=0;
	try
	{
		if(!PyArg_ParseTuple(args, "O:FldLenRndSw", &oOnOff)) throw CombErStr(strEr_BadFuncArg, ": FldLenRndSw");
		if(oOnOff == 0) throw CombErStr(strEr_BadFuncArg, ": FldLenRndSw");

		char sOnOff[256]; *sOnOff = '\0';
		CPyParse::CopyPyStringToC(oOnOff, sOnOff, 255);

		int res=0;
		g_pyParse.ProcRes(RadFldLenRndSw(&res, sOnOff));

		oRes = Py_BuildValue("i", res);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oRes;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Sets absolute and relative randomization magnitudes for all the length values, including coordinates and dimensions of the objects producing magnetic field, and coordinates of points where the field is computed. Optimal values of the variables can be: RelVal=10^(-11), AbsVal=L*RelVal, ZeroVal=AbsVal, where L is the distance scale value (in mm) for the problem to be solved. Too small randomization magnitudes can result in run-time code errors.
 ***************************************************************************/
static PyObject* radia_FldLenTol(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	try
	{
		double AbsVal=0, RelVal=0, ZeroVal=0;
		if(!PyArg_ParseTuple(args, "dd|d:FldLenTol", &AbsVal, &RelVal, &ZeroVal)) throw CombErStr(strEr_BadFuncArg, ": FldLenTol");
		//if((AbsVal <= 0) || (RelVal <= 0) || (ZeroVal <= 0)) throw CombErStr(strEr_BadFuncArg, ": FldLenTol");
		if((AbsVal <= 0) || (RelVal <= 0) || (ZeroVal < 0)) throw CombErStr(strEr_BadFuncArg, ": FldLenTol"); //AB21112019: ZeroVal=0 is allowed

		int dummy;
		g_pyParse.ProcRes(RadFldLenTol(&dummy, AbsVal, RelVal, ZeroVal));

		oRes = Py_BuildValue("i", dummy);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oRes;
}

/************************************************************************//**
 * Magnetic Field Calculation Methods: Computes a virtual "shim signature", i.e. variation of a given magnetic field component introduced by given displacement of magnetic field source object.
 ***************************************************************************/
static PyObject* radia_FldShimSig(PyObject* self, PyObject* args)
{
	PyObject *oCmpnId=0, *oD=0, *oP1=0, *oP2=0, *oV=0, *oResShimSig=0;
	double *arB=0;
	try
	{
		int ind=0, np=0;
		if(!PyArg_ParseTuple(args, "iOOOOi|O:FldShimSig", &ind, &oCmpnId, &oD, &oP1, &oP2, &np, &oV)) throw CombErStr(strEr_BadFuncArg, ": FldShimSig");
		if((oCmpnId == 0) || (oD == 0) || (oP1 == 0) || (oP2 == 0)) throw CombErStr(strEr_BadFuncArg, ": FldShimSig");

		char sCmpnId[256];
		CPyParse::CopyPyStringToC(oCmpnId, sCmpnId, 255);

		double arD[3], arP1[3], arP2[3], arV[]={0,0,0};
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oD, 'd', arD, 3, CombErStr(strEr_BadFuncArg, ": FldShimSig, incorrect definition of displacement vector"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP1, 'd', arP1, 3, CombErStr(strEr_BadFuncArg, ": FldShimSig, incorrect definition of the first point of the segment along which the shim signature should be calculated"));
		CPyParse::CopyPyListElemsToNumArrayKnownLen(oP2, 'd', arP2, 3, CombErStr(strEr_BadFuncArg, ": FldShimSig, incorrect definition of the last point of the segment along which the shim signature should be calculated"));
		if(oV != 0)
		{
			CPyParse::CopyPyListElemsToNumArrayKnownLen(oV, 'd', arV, 3, CombErStr(strEr_BadFuncArg, ": FldShimSig, incorrect definition of the vector specifying integration direction (for field integral calculation)"));
		}

		const int maxNumFldCmpn = 9;
		arB = new double[maxNumFldCmpn*np];
		int nTot = 0; //?
		g_pyParse.ProcRes(RadFldShimSig(arB, &nTot, ind, sCmpnId, arD, arP1, arP2, np, arV));

		if(nTot == 1) oResShimSig = Py_BuildValue("d", *arB);
		else if(nTot > 1) oResShimSig = CPyParse::SetDataListOfLists(arB, nTot, np);
		if(oResShimSig) Py_XINCREF(oResShimSig);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}
	if(arB != 0) delete[] arB;
	return oResShimSig;
}

/************************************************************************//**
 * Utilities: Outputs information about an object or list of objects.
 ***************************************************************************/
static PyObject* radia_UtiDmp(PyObject* self, PyObject* args)
{
	PyObject *oInds=0, *oType=0, *oRes=0;
	int *arInds=0;
	char *sDmp=0;
	try
	{
		if(!PyArg_ParseTuple(args, "O|O:UtiDmp", &oInds, &oType)) throw CombErStr(strEr_BadFuncArg, ": UtiDmp"); //OC20122019
		//if(!PyArg_ParseTuple(args, "OO:UtiDmp", &oInds, &oType)) throw CombErStr(strEr_BadFuncArg, ": UtiDmp");
		if(oInds == 0) throw CombErStr(strEr_BadFuncArg, ": UtiDmp");
		//if((oInds == 0) || (oType == 0)) throw CombErStr(strEr_BadFuncArg, ": UtiDmp");

		char sType[] = "asc\0\0"; //OC20122019
		//char sType[256];
		if(oType != 0) CPyParse::CopyPyStringToC(oType, sType, 4);

		int nInds = 0;
		if(PyList_Check(oInds))
		{
			bool lenIsSmall = false;
			CPyParse::CopyPyListElemsToNumArray(oInds, 'i', arInds, nInds, lenIsSmall);
			if((arInds == 0) || (nInds <= 0)) throw CombErStr(strEr_BadFuncArg, ": UtiDmp");
		}
		else if(PyNumber_Check(oInds))
		{
			int ind = PyLong_AsLong(oInds);
			if(ind <= 0) throw CombErStr(strEr_BadFuncArg, ": UtiDmp");
		
			arInds = new int[1]; //OC25022020
			//arInds = new int;
			*arInds = ind;
			nInds = 1;
		}
		else throw CombErStr(strEr_BadFuncArg, ": UtiDmp");

		int nBytes = 0;
		g_pyParse.ProcRes(RadUtiDmp(0, &nBytes, arInds, nInds, sType));

		if(nBytes > 0) 
		{
			sDmp = new char[nBytes];
			g_pyParse.ProcRes(RadUtiDataGet(sDmp, sType));

			if((strcmp(sType, "asc") == 0) || (strcmp(sType, "ASC") == 0))
			{
				oRes = Py_BuildValue("s", sDmp);
			}
			else if((strcmp(sType, "bin") == 0) || (strcmp(sType, "BIN") == 0))
			{
				oRes = CPyParse::Py_BuildValueByteStr(sDmp, nBytes);
			}

			if(oRes) Py_XINCREF(oRes);
		}
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arInds != 0) delete[] arInds;
	if(sDmp != 0) delete[] sDmp;
	return oRes;
}

/************************************************************************//**
 * Utilities: Parses byte-string produced by UtiDmp(elem,"bin") and attempts to instantiate 
 * elem objects(s); returns either index of one instantiated object 
 * (if elem was an index of one object) or a list of indexes of instantiated 
 * objects (if elem was a list of objects).
 ***************************************************************************/
static PyObject* radia_UtiDmpPrs(PyObject* self, PyObject* args)
{
	PyObject *oBytes=0, *oRes=0;
	int *arInds=0;
	try
	{
		if(!PyArg_ParseTuple(args, "O:UtiDmpPrs", &oBytes)) throw CombErStr(strEr_BadFuncArg, ": UtiDmpPrs");
		if(oBytes == 0) throw CombErStr(strEr_BadFuncArg, ": UtiDmpPrs");

		char *sBytes=0;
		Py_ssize_t nBytes=0;
		//int nBytes=0;
		//CPyParse::CopyPyByteArrToC(oBytes, sBytes, nBytes);
		if(PyBytes_AsStringAndSize(oBytes, &sBytes, &nBytes) == -1) throw strEr_BadStr; //OC04102018
		//No deallocation of sBytes is required after this!
		if((sBytes == 0) || (nBytes <= 0)) throw CombErStr(strEr_BadFuncArg, ": UtiDmpPrs, object(s) can not be instantiated from this string/byte array");

		int nElem = 0;
		g_pyParse.ProcRes(RadUtiDmpPrs(0, &nElem, (unsigned char*)sBytes, (int)nBytes));
		if(nElem <= 0) throw CombErStr(strEr_BadFuncArg, ": UtiDmpPrs, object(s) can not be instantiated from this string/byte array");

		bool resIsList = (bool)sBytes[0];
		if(resIsList)
		{
			arInds = new int[nElem];
			g_pyParse.ProcRes(RadUtiDataGet((char*)arInds, "mai", 0));
		
			if(nElem == 1) oRes = Py_BuildValue("i", *arInds); //if list has only one element, returning this element (not list)
			else oRes = CPyParse::SetDataListOfLists(arInds, nElem, 1, "i");
			delete[] arInds;
		}
		else
		{
			int auxInd = 0; 
			g_pyParse.ProcRes(RadUtiDataGet((char*)(&auxInd), "i", 0));

			oRes = Py_BuildValue("i", auxInd);
		}

		if(oRes!= 0) Py_XINCREF(oRes); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oRes;
}

/************************************************************************//**
 * Utilities: Deletes an object.
 ***************************************************************************/
static PyObject* radia_UtiDel(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	try
	{
		int ind=0;
		if(!PyArg_ParseTuple(args, "i:UtiDel", &ind)) throw CombErStr(strEr_BadFuncArg, ": UtiDel");

		int nDummy = 0;
		g_pyParse.ProcRes(RadUtiDel(&nDummy, ind));

		oRes = Py_BuildValue("i", nDummy);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oRes;
}

/************************************************************************//**
 * Deletes all previously created objects.
 ***************************************************************************/
static PyObject* radia_UtiDelAll(PyObject* self, PyObject* args)
{
	PyObject *oRes=0;
	try
	{
		int nDummy = 0;
		g_pyParse.ProcRes(RadUtiDelAll(&nDummy));

		oRes = Py_BuildValue("i", nDummy);
		Py_XINCREF(oRes); //?
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oRes;
}

/************************************************************************//**
 * Utilities: Returns Radia library version number.
 ***************************************************************************/
static PyObject* radia_UtiVer(PyObject* self, PyObject* args)
{
	PyObject *oVerNum = 0;
	try
	{
		double verNum = 0.;
		g_pyParse.ProcRes(RadUtiVer(&verNum));

		oVerNum = Py_BuildValue("d", verNum);
		if(oVerNum) Py_XINCREF(oVerNum);
	}
	catch (const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	return oVerNum;
}

/************************************************************************//**
 * Utilities: Initializes or finishes MPI (to support parallel computation)
 ***************************************************************************/
static PyObject* radia_UtiMPI(PyObject* self, PyObject* args)
{
	PyObject *oOnOff=0, *oParMPI=0, *oData=0;
	double *arData=0, *arDataAux=0;
	try
	{
		long rankFrom=0, rankTo=-1;
		if(!PyArg_ParseTuple(args, "O|Oii:UtiMPI", &oOnOff, &oData, &rankFrom, &rankTo)) throw CombErStr(strEr_BadFuncArg, ": UtiMPI"); //OC19032020
		//if(!PyArg_ParseTuple(args, "O:UtiMPI", &oOnOff)) throw CombErStr(strEr_BadFuncArg, ": UtiMPI");
		if(oOnOff == 0) throw CombErStr(strEr_BadFuncArg, ": UtiMPI");

		char sOnOff[256]; *sOnOff = '\0';
		long nData = 0; //OC19032020
		CPyParse::CopyPyStringToC(oOnOff, sOnOff, 255);
		if(oData != 0)
		{
			if((strcmp(sOnOff, "share") == 0) || (strcmp(sOnOff, "Share") == 0) || (strcmp(sOnOff, "SHARE") == 0))
			{
				if(PyList_Check(oData) || PyObject_CheckBuffer(oData))
				{
					//bool lenIsSmall = false;
					//CPyParse::CopyPyListElemsToNumArray(oData, 'd', arData, nData, lenIsSmall);
					CPyParse::CopyPyNestedListElemsToNumAr(oData, 'd', arData, nData);
				}
			}
		}

		int arParMPI[] = {-1,0}; //to obtain world rank and size
		g_pyParse.ProcRes(RadUtiMPI(arParMPI, sOnOff, arData, &nData, &rankFrom, &rankTo));
		//g_pyParse.ProcRes(RadUtiMPI(arParMPI, sOnOff));

		if((oData != 0) && (nData > 0)) //OC19032020
		{//Extract data that was shared from buffer
			arDataAux = new double[nData];
			g_pyParse.ProcRes(RadUtiDataGet((char*)arDataAux, "mad")); //OC22032020
			//g_pyParse.ProcRes(RadUtiDataGet((char*)arDataAux, "d", rankTo)); //rankTo is used as key for Data here
			oParMPI = CPyParse::SetDataListOfLists(arDataAux, nData, (long)1, "d");
		}
		else
		{
			oParMPI = Py_BuildValue("i", arParMPI[0]); //returning just rank to simplify its use in scripts
		}

		//oParMPI = Py_BuildValue("i", arParMPI[0]); //returning just rank to simplify its use in scripts
		//oParMPI = CPyParse::SetDataTuple(arParMPI, 2, (char*)"i");
		//oParMPI = CPyParse::SetDataListOfLists(arParMPI, 2, 1, (char*)"i");
		if(oParMPI) Py_XINCREF(oParMPI);
	}
	catch(const char* erText)
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arData != 0) delete[] arData;
	if(arDataAux != 0) delete[] arDataAux;
	return oParMPI;
}

/************************************************************************//**
 * Python C API stuff: module & method definition2, etc.
 ***************************************************************************/

static PyMethodDef radia_methods[] = {
	{"ObjRecMag", radia_ObjRecMag, METH_VARARGS, "ObjRecMag([x,y,z],[wx,wy,wz],[mx,my,mz]:[0,0,0]) creates a rectangular parallelepiped block with center point [x,y,z], dimensions [wx,wy,wz], and magnetization [mx,my,mz]."},
	{"ObjThckPgn", radia_ObjThckPgn, METH_VARARGS, "ObjThckPgn(x,lx,[[y1,z1],[y2,z2],...],a:'x',[mx,my,mz]:[0,0,0]) creates an extruded polygon block; x is the position of the block's center of gravity in the extrusion direction, lx is the thickness, [[y1,z1],[y2,z2],...] is a list of points describing the polygon in 2D; the extrusion direction is defined by the character a (which can be 'x', 'y' or 'z'), [mx,my,mz] is the block magnetization."},
	{"ObjPolyhdr", radia_ObjPolyhdr, METH_VARARGS, "ObjPolyhdr([[x1,y1,z1],[x2,y2,z2],...],[[f1i1,f1i2,...],[f2i1,f2i2,...],...],[mx,my,mz]:[0,0,0],J:[jx,jy,jz]|[[jx,jy,jz],[[djxdy,djxdy,djxdz],[djydy,djydy,djydz],[djzdy,djzdy,djzdz]]],Lin:'Rel') creates a uniformly magnetized polyhedron (closed volume limited by planes). [[x1,y1,z1],[x2,y2,z2],...] is a list of the polyhedron vertex points, [[f1n1,f1n2,...],[f2n1,f2n2,...],...] is a list of lists of indexes of vertex points defining the polyhedron faces, [mx,my,mz] is magnetization inside the polyhedron. The optional parameter J can be used to define constant [jx,jy,jz] or linearly-varying current density vector inside the polyhedron; the linear dependence can be defined through 3x3 matrix of coefficients [[djxdy,djxdy,djxdz],[djydy,djydy,djydz],[djzdy,djzdy,djzdz]]; depending on the value of the optional parameter Lin this linear dependence is treated with respect to the polyhedron center (Lin='Rel', default) or with respect to the origin of the Cartesian frame (Lin='Abs')."},
	{"ObjMltExtPgn", radia_ObjMltExtPgn, METH_VARARGS, "ObjMltExtPgn([[[[x11,y11],[x12,y12],...],z1],[[[x21,y21],[x22,y22],...],z2],...],[mx,my,mz]:[0,0,0]) attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on slices. The slice polygons are defined by the nested list [[[[x11,y11],[x12,y12],...],z1],[[[x21,y21],[x22,y22],...],z2],...], with [[x11,y11],[x12,y12],...],... describing the polygons in 2D, and z1, z2,... giving their attitudes (vertical coordinates). [mx,my,mz] is the magnetization inside the polyhedron(s) created."},
	//{"ObjMltExtPgnCur", radia_ObjMltExtPgnCur, METH_VARARGS, "ObjMltExtPgnCur(z:0,a:\"z\",{{{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...}},I,Frame->Loc|Lab) attempts to create a set of current-carrying convex polyhedron objects by applying a generalized extrusion to the initial planar convex polygon. The initial polygon is defined for the \"attitude\" z (0 by default) by the list of 2D points {{x1,y1},{x2,y2},...}, with the  a  character specifying orientation of this polygon normal in 3D space: if a = \"z\" (default orientation), the polygon is assumed to be parallel to XY plane of the laboratory frame (\"y\" for ZX plane, \"x\" for YZ plane). The extrusion can consist of a number of \"steps\", with each step creating one convex polyhedron defined optionally by one (or combination of) rotation(s), and/or translation(s), and/or one homothety: {Rk,Tk,Hk}, k = 1,2,..., applied to the base polygon (i.e. either the initial base polygon, or the polygon obtained by previous extrusion step). In case if k-th extrusion step contains one rotation Rk about an axis, it is defined as {{xRk,yRk,zRk},{vxRk,vyRk,vzRk},phRk}}, where {xRk,yRk,zRk} and {vxRk,vyRk,vzRk} are respectively 3D coordinates of a point and a vector difining the rotation axis, and phRk the rotation angle in radians; in case if Rk is a combination of \"atomic\" rotations about different axes, it should be defined as list: {Rk1,Rk2,...}. If k-th extrusion step includes translation Tk, it must be defined by vector {vxTk,vyTk,vzTk}; optional homothety with respect to the base polygon center of gravity should be defined either by two different coefficients {pxHk,pyHk} with respect to two orthogonal axes of the base polygon local frame, or by nested list {{pxHk,pyHk},phHk}, where phHk is rotation angle of the two homothety axes in radians. A real number I defines average current in Amperes along the extrusion path. The Frame->Loc or Frame->Lab option specifies whether the transformations at each step of the extrusion path are defined in the frame of the previous base polygon (Frame->Loc, default), or all the transformations are defined in the laboratory frame (Frame->Lab)."},
	//{"ObjMltExtPgnMag", radia_ObjMltExtPgnMag, METH_VARARGS, "ObjMltExtPgnMag(z:0,a:\"z\",{{{x1,y1},{x2,y2},...},{{k1,q1},{k2,q2},...}:{{1,1},{1,1},...},{{R1,T1,H1},{R2,T2,H2},...}},{{mx1,my1,mz1},{mx2,my2,mz2},...}:{{0,0,0},{0,0,0},...},Frame->Loc|Lab,ki->Numb|Size,TriAngMin->...,TriAreaMax->...,TriExtOpt->\"...\") attempts to create a set of uniformly magnetized polyhedron objects by applying a generalized extrusion to the initial planar convex polygon. The initial polygon is defined for the \"altitude\" z (0 by default) by the list of 2D points {{x1,y1},{x2,y2},...}, with the  a  character specifying orientation of this polygon normal in 3D space: if a = \"z\" (default orientation), the polygon is assumed to be parallel to XY plane of the laboratory frame (\"y\" for ZX plane, \"x\" for YZ plane). The extrusion can consist of a number of \"steps\", with each step creating one convex polyhedron defined optionally by one (or combination of) rotation(s), and/or translation(s), and/or one homothety: {Rk,Tk,Hk}, k = 1,2,..., applied to the base polygon (i.e. either the initial base polygon, or the polygon obtained by previous extrusion step). In case if k-th extrusion step contains one rotation Rk about an axis, it is defined as {{xRk,yRk,zRk},{vxRk,vyRk,vzRk},phRk}}, where {xRk,yRk,zRk} and {vxRk,vyRk,vzRk} are respectively 3D coordinates of a point and a vector difining the rotation axis, and phRk the rotation angle in radians; in case if Rk is a combination of \"atomic\" rotations about different axes, it should be defined as a list: {Rk1,Rk2,...}. If k-th extrusion step includes translation Tk, it must be defined by vector {vxTk,vyTk,vzTk}; optional homothety with respect to the base polygon center of gravity should be defined either by two different coefficients {pxHk,pyHk} with respect to two orthogonal axes of the base polygon local frame, or by nested list {{pxHk,pyHk},phHk}, where phHk is rotation angle of the two homothety axes in radians. Optional list {{mx1,my1,mz1},{mx2,my2,mz2},...} defines magnetization vectors in each of polyhedrons to be created. The Frame->Loc or Frame->Lab option specifies whether the transformations at each step of the extrusion path are defined in the frame of the previous base polygon (Frame->Loc, default), or all the transformations are defined in the laboratory frame (Frame->Lab). Optionally, the object can be subdivided by (extruded) triangulation at its creation; this occurs if {{k1,q1},{k2,q2},...} subdivision (triangulation) parameters for each segment of the base polygon border are defined; the meaning of k1, k2,... depends on value of the option ki: if ki->Numb (default), then k1, k2,... are subdivision numbers; if ki->Size, they are average sizes of sub-segments to be produced; q1, q2,... are ratios of the last-to-first sub-segment lengths; the TriAngMin option defines minimal angle of triangles to be produced (in degrees, default is 20); the TriAreaMax option defines maximal area of traingles to be produced (in mm^2, not defined by default); the ExtOpt option allows to specify additional parameters for triangulation function in a string." },
	{"ObjMltExtRtg", radia_ObjMltExtRtg, METH_VARARGS, "ObjMltExtRtg([[[x1,y1,z1],[wx1,wy1]],[[x2,y2,z2],[wx2,wy2]],...],[mx,my,mz]:[0,0,0]) attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on rectangular slices. The slice rectangles are defined by the nested list [[[x1,y1,z1],[wx1,wy1]],[[x2,y2,z2],[wx2,wy2]],...], with [x1,y1,z1], [x2,y2,z2],... being center points of the rectangles, and [wx1,wy1], [wx2,wy2],... their dimensions. [mx,my,mz] is the magnetization inside the polyhedron(s) created."},
	{"ObjMltExtTri", radia_ObjMltExtTri, METH_VARARGS, "ObjMltExtTri(x,lx,[[y1,z1],[y2,z2],...],[[k1,q1],[k2,q2],...],a:'x',[mx,my,mz]:[0,0,0],Opt:'ki->Numb|Size,TriAngMin->...,TriAreaMax->...') creates triangulated extruded polygon block, i.e. a straight prism with its bases being (possibly non-convex) polygons subdivided by triangulation. x is the position of the block's center of gravity in the extrusion direction, lx is the thickness, [[y1,z1],[y2,z2],...] is a list of points describing the polygon in 2D; [[k1,q1],[k2,q2],...] are subdivision (triangulation) parameters for each segment of the base polygon border; the meaning of k1, k2,... depends on value of the option ki (to be defined in the string Opt): if ki->Numb (default), then k1, k2,... are subdivision numbers; if ki->Size, they are average sizes of sub-segments to be produced; q1, q2,... are ratios of the last-to-first sub-segment lengths; the extrusion direction is defined by the character a (which can be 'x', 'y' or 'z'); [mx,my,mz] is magnetization inside the block. The TriAngMin option defines minimal angle of triangles to be produced (in degrees, default is 20); the TriAreaMax option defines maximal area of traingles to be produced (in mm^2, not defined by default)."},
	{"ObjArcPgnMag", radia_ObjArcPgnMag, METH_VARARGS, "ObjArcPgnMag([x,y],a,[[r1,z1],[r2,z2],...],[phimin,phimax],nseg,'sym|nosym':'nosym',[mx,my,mz]:[0,0,0]) creates a uniformly magnetized finite-length arc of polygonal cross-section with the position and orientation of the rotation axis defined by pair of coordinates {x,y} and character a (which can be 'x', 'y' or 'z'), the cross-section 2D polygon [[r1,z1],[r2,z2],...], initial and final rotation angles [phimin,phimax], number of sectors (segments) vs azimuth nseg, and magnetization vector [mx,my,mz]. Depending on the value of the 'sym|nosym' switch, the magnetization vectors in nseg sector polyhedrons are either assumed to satisfy rotational symmetry conditions ('sym'), or are assumed independent (i.e. will be allowed to vary independently at further relaxation)."},
	{"ObjCylMag", radia_ObjCylMag, METH_VARARGS, "ObjCylMag([x,y,z],r,h,nseg,a:'z',[mx,my,mz]:[0,0,0]) creates a cylindrical magnet approximated by a straight prism with a right polygon in base with center point [x,y,z], base radius r, height h, number of segments nseg, orientation of the rotation axis defined by character a (which can be 'x', 'y' or 'z'), and magnetization vector [mx,my,mz]."},
	{"ObjFullMag", radia_ObjFullMag, METH_VARARGS, "ObjFullMag([x,y,z],[wx,wy,wz],[mx,my,mz],[kx,ky,kz],cnt,mat,[r,g,b]) creates rectangular parallelepiped block with constant magnetizatiom over volume, center point [x,y,z], dimensions [wx,wy,wz] and color [r,g,b]. The block is magnetized according to [mx,my,mz], subdivided according to [kx,ky,kz] and added into the container cnt, that should be defined in advance by calling ObjCnt()."},
	{"ObjRecCur", radia_ObjRecCur, METH_VARARGS, "ObjRecCur([x,y,z],[wx,wy,wz],[jx,jy,jz]) creates a current carrying rectangular parallelepiped block with center point [x,y,z], dimensions [wx,wy,wz], and current density [jx,jy,jz]."},
	{"ObjArcCur", radia_ObjArcCur, METH_VARARGS, "ObjArcCur([x,y,z],[rmin,rmax],[phimin,phimax],h,nseg,j,'man|auto':'man',a:'z') creates a current carrying finite-length arc of rectangular cross-section, with center point [x,y,z], inner and outer radii [rmin,rmax], initial and final angles [phimin,phimax], height h, number of segments nseg, and azimuthal current density j. According to the value of the 'man|auto' switch, the field from the arc is computed based on the number of segments nseg ('man'), or on the general absolute precision level specified by the function FldCmpCrt ('auto'). The orientation of the rotation axis is defined by the character a (which can be either 'x', 'y' or 'z')."},
	{"ObjRaceTrk", radia_ObjRaceTrk, METH_VARARGS, "ObjRaceTrk([x,y,z],[rmin,rmax],[lx,ly],h,nseg,j,'man|auto':'man',a:'z') creates a current carrying racetrack coil consisting of four 90-degree bents connected by four straight parts of rectangular straight section, center point [x,y,z], inner and outer bent radii [rmin,rmax], straight section lengths [lx,ly], height h, number of segments in bents nseg, and azimuthal current density j. According to the value of the 'man|auto' switch, the field from the bents is computed based on the number of segments nseg ('man'), or on the general absolute precision level specified by the function FldCmpCrt ('auto'). The orientation of the racetrack axis is defined by the character a (which can be either 'x', 'y' or 'z')."},
	{"ObjFlmCur", radia_ObjFlmCur, METH_VARARGS, "ObjFlmCur([[x1,y1,z1],[x2,y2,z2],...],i) creates a filament polygonal line conductor defined by the sequence of points [[x1,y1,z1],[x2,y2,z2],...] with current i."},
	{"ObjBckg", radia_ObjBckg, METH_VARARGS, "ObjBckg([bx,by,bz]) creates a source of uniform background magnetic field [bx,by,bz]."},
	{"ObjCnt", radia_ObjCnt, METH_VARARGS, "ObjCnt([obj1,obj2,...]) creates a container object for magnetic field source objects [obj1,obj2,...]."},
	{"ObjAddToCnt", radia_ObjAddToCnt, METH_VARARGS, "ObjAddToCnt(cnt,[obj1,obj2,...]) adds objects [obj1,obj2,...] to the container object cnt."},
	{"ObjCntStuf", radia_ObjCntStuf, METH_VARARGS, "ObjCntStuf(obj) returns list of general indexes of the objects present in container if obj is a container; or returns [obj] if obj is not a container."}, 
	{"ObjCntSize", radia_ObjCntSize, METH_VARARGS, "ObjCntSize(cnt) calculates the number of objects in the container cnt."},
	{"ObjCutMag", radia_ObjCutMag, METH_VARARGS, "ObjCutMag(obj,[x,y,z],[nx,ny,nz],'Frame->Loc|Lab|LabTot') cuts 3D object by a plane normal to the vector [nx,ny,nz] and passing through the point [x,y,z]. The 'Frame->Loc', 'Frame->Lab' or 'Frame->LabTot' option specifies whether the cuting plane is defined in the local frame of the object obj or in the laboratory frame (default). The actions of 'Frame->Lab' and 'Frame->LabTot' differ for containers only: 'Frame->Lab' means that each of the objects in the container is cut separately; 'Frame->LabTot' means that the objects in the container are cut as one object, by the same plane. The function returns a list of indexes of the objects produced by the cutting."},
	{"ObjDivMag", radia_ObjDivMag, METH_VARARGS, "ObjDivMag(obj,[[k1,q1],[k2,q2],[k3,q3]],'pln|cyl',[[n1x,n1y,n1z],[n2x,n2y,n2z],[n3x,n3y,n3z]]|[[ax,ay,az],[vx,vy,vz],[px,py,pz],rat],'kxkykz->Numb|Size,Frame->Loc|Lab|LabTot') subdivides (segments) a 3D object obj. The main subdivision parameters are defined by the list [[k1,q1],[k2,q2],[k3,q3]] or [k1,k2,k3]. The meaning of k1, k2 and k3 depends on the value of the option kxkykz: if kxkykz->Numb (default), then k1, k2 and k3 are subdivision numbers; if kxkykz->Size, these are average sizes of the sub-objects to be produced. q1, q2 and q3 in any case are ratios of the last-to-first sub-object sizes (if these parameters are omitted, they are assumed to be equal to 1). The third string variable defines type of the subdivision. If it is not used, the subdivision is performed in directions X, Y and Z (by a set of parallel planes). If it is equal to 'pln', the function performs subdivision by three sets of parallel planes normal to the vectors [n1x,n1y,n1z], [n2x,n2y,n2z] and [n3x,n3y,n3z], that should be defined by the next list variable [[n1x,n1y,n1z],[n2x,n2y,n2z],[n3x,n3y,n3z]]. The distances between the parallel planes are defined by the parameters [k1,q1],[k2,q2] and [k3,q3]. If the third variable is equal to 'cyl', the function performs subdivision by a system of coaxial elliptic cylinders, with their parameters defined by the next list variable, [[ax,ay,az],[vx,vy,vz],[px,py,pz],rat]. The cylinder axis is defined by the point [ax,ay,az] and vector [vx,vy,vz]. One of two axes of the cylinder base ellipses is exactly the perpendicular from the point [px,py,pz] to the cylinder axis; rat is the ratio of the ellipse axes lengths. In the case of the subdivision by elliptic cylinders, the parameters [k1,q1],[k2,q2] and [k3,q3] correspond to radial, azimuthal, and axial directions respectively. If the option Frame is set to 'Frame->Loc' (default), the subdivision is performed in local frame(s) of the object(s); if it is set to 'Frame->Lab' or 'Frame->LabTot', the subdivision is performed in the laboratory frame. The actions of 'Frame->Lab' and 'Frame->LabTot' differ for containers only: 'Frame->Lab' means that each of the objects in the container is subdivided separately; 'Frame->LabTot' means that the objects in the container are subdivided as one object, by the same planes."},
	{"ObjDivMagPln", radia_ObjDivMagPln, METH_VARARGS, "ObjDivMagPln(obj,[[k1,q1],[k2,q2],[k3,q3]],[n1x,n1y,n1z],[n2x,n2y,n2z],[n3x,n3y,n3z],'kxkykz->Numb|Size,Frame->Loc|Lab|LabTot') subdivides (segments) a 3D object by 3 sets of parallel planes. The main subdivision parameters are defined by the list [[k1,q1],[k2,q2],[k3,q3]] or [k1,k2,k3]. The meaning of k1, k2 and k3 depends on the value of the option kxkykz: if kxkykz->Numb (default), then k1, k2 and k3 are subdivision numbers; if kxkykz->Size, these are average sizes of the sub-objects to be produced. q1, q2 and q3 in any case are ratios of the last-to-first sub-object sizes (if these parameters are omitted, they are assumed to be equal to 1). The orientation of the three sets of parallel planes normal to the vectors [n1x,n1y,n1z], [n2x,n2y,n2z] and [n3x,n3y,n3z]. If these variables are not submitted, the subdivision is performed in directions X, Y and Z. The distances between the parallel planes are defined by the parameters [k1,q1],[k2,q2] and [k3,q3]. If the option Frame is set to 'Frame->Loc' (default), the subdivision is performed in local frame(s) of the object(s); if it is set to 'Frame->Lab' or 'Frame->LabTot', the subdivision is performed in the laboratory frame. The actions of 'Frame->Lab' and 'Frame->LabTot' differ for containers only: 'Frame->Lab' means that each of the objects in the container is subdivided separately; 'Frame->LabTot' means that the objects in the container are subdivided as one object, by the same planes."},
	{"ObjDivMagCyl", radia_ObjDivMagCyl, METH_VARARGS, "ObjDivMagCyl(obj,[[k1,q1],[k2,q2],[k3,q3]],[ax,ay,az],[vx,vy,vz],[px,py,pz],rat,'kxkykz->Numb|Size,Frame->Loc|Lab|LabTot') subdivides (segments) a 3D object obj by a set of coaxial elliptical cylinders. The main subdivision parameters are defined by the list [[k1,q1],[k2,q2],[k3,q3]] or [k1,k2,k3]. The meaning of k1, k2 and k3 depends on the value of the option kxkykz: if kxkykz->Numb (default), then k1, k2 and k3 are subdivision numbers; if kxkykz->Size, these are average sizes of the sub-objects to be produced. q1, q2 and q3 in any case are ratios of the last-to-first sub-object sizes (if these parameters are omitted, they are assumed to be equal to 1). The cylinder axis is defined by the point [ax,ay,az] and vector [vx,vy,vz]. One of two axes of the cylinder base ellipses is exactly the perpendicular from the point [px,py,pz] to the cylinder axis; rat is the ratio of the ellipse axes lengths. The parameters [k1,q1],[k2,q2] and [k3,q3] correspond to radial, azimuthal, and axial directions respectively. If the option Frame is set to 'Frame->Loc' (default), the subdivision is performed in local frame(s) of the object(s); if it is set to 'Frame->Lab' or 'Frame->LabTot', the subdivision is performed in the laboratory frame. The actions of 'Frame->Lab' and 'Frame->LabTot' differ for containers only: 'Frame->Lab' means that each of the objects in the container is subdivided separately; 'Frame->LabTot' means that the objects in the container are subdivided as one object, by the same planes."},
	{"ObjDpl", radia_ObjDpl, METH_VARARGS, "ObjDpl(obj,'FreeSym->False|True') duplicates 3D object obj. The option 'FreeSym->False|True' specifies whether the symmetries (transformations with multiplicity more than one) previously applied to the object obj should be simply copied at the duplication ('FreeSym->False', default), or a container of new independent objects should be created in place of any symmetry previously applied to the object obj. In both cases the final object created by the duplication has exactly the same geometry as the initial object obj."},
	{"ObjGeoVol", radia_ObjGeoVol, METH_VARARGS, "ObjGeoVol(obj) computes geometrical volume of 3D object obj."},
	{"ObjGeoLim", radia_ObjGeoLim, METH_VARARGS, "ObjGeoLim(obj) computes geometrical limits of 3D object obj in the laboratory frame. Returns [xmin, xmax, ymin, ymax, zmin, zmax]."},
	{"ObjDegFre", radia_ObjDegFre, METH_VARARGS, "ObjDegFre(obj) gives number of degrees of freedom for the relaxation of 3D object obj."},
	{"ObjM", radia_ObjM, METH_VARARGS, "ObjM(obj) returns coordinates of geometrical center point and magnetization vector of 3D object obj at that point. If obj is a container, a list of the container members' center points and their magnetizations is returned."},
	{"ObjCenFld", radia_ObjCenFld, METH_VARARGS, "ObjCenFld(obj,'A|B|H|J|M') provides coordinates of geometrical center point of the object obj and a field characteristic vector at that point. The type of field characteristic is defined by the second parameter (character); it can be one of the following: 'A' for vector potential, 'B' for magnetic field induction, 'H' for magnetic field strength, 'J' for current density, 'M' for magnetization. If obj is a container, a list of the container members' center points and their field characteristics is returned."},
	{"ObjSetM", radia_ObjSetM, METH_VARARGS, "ObjSetM(obj,[mx,my,mz]) sets magnetization [mx,my,mz] in 3D object obj."},
	{"ObjScaleCur", radia_ObjScaleCur, METH_VARARGS, "ObjScaleCur(obj,k) scales current (density) in 3D object obj by multiplying it by constant k (if obj is a current-carying object). If obj is a container, the current (density) scaling applies to all its members."},
	{"ObjDrwAtr", radia_ObjDrwAtr, METH_VARARGS, "ObjDrwAtr(obj,[r,g,b],thcn) assigns drawing attributes - RGB color [r,g,b] and line thickness thcn - to 3D object obj."},
	{"ObjDrwOpenGL", radia_ObjDrwOpenGL, METH_VARARGS, "ObjDrwOpenGL(obj,'EdgeLines->True|False,Faces->True|False,Axes->True|False') starts an application for viewing 3D geometry of the object obj. The viewer is based on the GLUT / OpenGL graphics library. The option 'EdgeLines->True|False' (default 'EdgeLines->True') highlights the edge lines of objects; the option 'Faces->True|False' (default 'Faces->True') shows faces of the objects; the option 'Axes->True|False' (default 'Axes->True') shows the Cartesian frame axes."},
	{"ObjDrwVTK", radia_ObjDrwVTK, METH_VARARGS, "ObjDrwVTK(obj,'EdgeLines->True|False,Faces->True|False,Axes->True|False') exports data for viewing 3D geometry of the object obj. The data is in the format compatible with VTK graphics library. The option 'EdgeLines->True|False' (default 'EdgeLines->True') highlights the edge lines of objects; the option 'Faces->True|False' (default 'Faces->True') shows faces of the objects; the option 'Axes->True|False' (default 'Axes->True') shows the Cartesian frame axes."},

	{"TrfTrsl", radia_TrfTrsl, METH_VARARGS, "TrfTrsl([vx,vy,vz]) creates a translation by vector [vx,vy,vz]."},
	{"TrfRot", radia_TrfRot, METH_VARARGS, "TrfRot([x,y,z],[vx,vy,vz],phi) creates a rotation of angle phi around the axis defined by the point [x,y,z] and the vector [vx,vy,vz]."},
	{"TrfPlSym", radia_TrfPlSym, METH_VARARGS, "TrfPlSym([x,y,z],[nx,ny,nz]) creates a symmetry with respect to plane defined by point [x,y,z] and normal vector [nx,ny,nz]."},
	{"TrfInv", radia_TrfInv, METH_VARARGS, "TrfInv() creates a field inversion."},
	{"TrfCmbL", radia_TrfCmbL, METH_VARARGS, "TrfCmbL(trfOrig,trf) multiplies original transformation trfOrig by another transformation trf from left."},
	{"TrfCmbR", radia_TrfCmbR, METH_VARARGS, "TrfCmbR(trfOrig,trf) multiplies original transformation trfOrig by another transformation trf from right."},
	{"TrfMlt", radia_TrfMlt, METH_VARARGS, "TrfMlt(obj,trf,mlt) creates mlt-1 objects. Each object is derived from the previous one by applying the transformation trf. Following this, the object obj becomes equivalent to mlt different objects."},
	{"TrfOrnt", radia_TrfOrnt, METH_VARARGS, "TrfOrnt(obj,trf) orients object obj by applying transformation trf to it once."},
	{"TrfZerPara", radia_TrfZerPara, METH_VARARGS, "TrfZerPara(obj,[x,y,z],[nx,ny,nz]) creates an object mirror of obj with respect to the plane with normal [nx,ny,nz] and passing by the point [x,y,z]. The object mirror presents the same geometry as obj, but its magnetization and/or current densities are modified in such a way that the magnetic field produced by the obj and its mirror in the plane of mirroring is PERPENDICULAR to this plane."},
	{"TrfZerPerp", radia_TrfZerPerp, METH_VARARGS, "TrfZerPerp(obj,[x,y,z],[nx,ny,nz]) creates an object mirror of obj with respect to the plane with normal [nx,ny,nz] and passing by the point [x,y,z]. The object mirror presents the same geometry as obj, but its magnetization and/or current densities are modified in such a way that the magnetic field produced by the obj and its mirror in the plane of mirroring is PARALLEL to this plane."},

	{"MatLin", radia_MatLin, METH_VARARGS, "MatLin([ksipar,ksiper],mr) or MatLin([ksipar,ksiper],[mrx,mry,mrz]) creates a linear anisotropic magnetic material with susceptibilities parallel (perpendicular) to the easy magnetization axis given by ksipar (ksiper). In the first form of the function, mr is the magnitude of the remanent magnetization vector; the direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is applied (the magnetization vector is specified at the object creation). In the second form, [mrx,mry,mrz] is the remanent magnetization vector explicitly defining the direction of the easy magnetization axis in any object to which the material is later applied."},
	{"MatSatIsoFrm", radia_MatSatIsoFrm, METH_VARARGS, "MatSatIsoFrm([ksi1,ms1],[ksi2,ms2],[ksi3,ms3]) creates a nonlinear isotropic magnetic material with the M versus H curve defined by the formula M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3), where H is the magnitude of the magnetic field strength vector (in Tesla). The parameters [ksi3,ms3] and [ksi2,ms2] may be omitted; in such a case the corresponding terms in the formula will be omitted too."},
	{"MatSatIsoTab", radia_MatSatIsoTab, METH_VARARGS, "MatSatIsoTab([[H1,M1],[H2,M2],...]) creates a nonlinear isotropic magnetic material with the M versus H curve defined by the list of pairs [[H1,M1],[H2,M2],...] in Tesla."},
	{"MatSatLamFrm", radia_MatSatLamFrm, METH_VARARGS, "MatSatLamFrm([ksi1,ms1],[ksi2,ms2],[ksi3,ms3],p,[nx,ny,nz]) creates laminated nonlinear anisotropic magnetic material with packing factor p and the lamination planes perpendicular to the vector [nx,ny,nz]. The magnetization magnitude vs magnetic field strength for the corresponding isotropic material is defined by the formula M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3), where H is the magnitude of the magnetic field strength vector (in Tesla). The parameters [ksi3,ms3] and [ksi2,ms2] may be omitted; in such a case the corresponding terms in the formula will be omitted too."},
	{"MatSatLamTab", radia_MatSatLamTab, METH_VARARGS, "MatSatLamTab([[H1,M1],[H2,M2],...],p,[nx,ny,nz]) creates laminated nonlinear anisotropic magnetic material with packing factor p and the lamination planes perpendicular to the vector [nx,ny,nz]. The magnetization magnitude vs magnetic field strength for the corresponding isotropic material is defined by pairs [[H1,M1],[H2,M2],...] in Tesla."},
	{"MatSatAniso", radia_MatSatAniso, METH_VARARGS, "MatSatAniso(datapar,dataper) where datapar can be [[ksi1,ms1,hc1],[ksi2,ms2,hc2],[ksi3,ms3,hc3],[ksi0,hc0]] or ksi0, and dataper can be [[ksi1,ms1],[ksi2,ms2],[ksi3,ms3],ksi0] or ksi0 - creates a nonlinear anisotropic magnetic material. If the first argument is set to [[ksi1,ms1,hc1],[ksi2,ms2,hc2],[ksi3,ms3,hc3],[ksi0,hc0]], the magnetization vector component parallel to the easy axis is computed as ms1*tanh(ksi1*(hpa-hc1)/ms1) + ms2*tanh(ksi2*(hpa-hc2)/ms2) + ms3*tanh(ksi3*(hpa-hc3)/ms3) + ksi0*(hpa-hc0), where hpa is the field strength vector component parallel to the easy axis. If the second argument is set to [[ksi1,ms1],[ksi2,ms2],[ksi3,ms3],ksi0], the magnetization vector component perpendicular to the easy axis is computed as ms1*tanh(ksi1*hpe/ms1) + ms2*tanh(ksi2*hpe/ms2) + ms3*tanh(ksi3*hpe/ms3) + ksi0*hpe, where hpe is the field strength vector component perpendicular to the easy axis. If the first or second argument is set to ksi0, the magnetization component parallel (perpendicular) to the easy axis is computed by ksi0*hp, where hp is the corresponding component of field strength vector. At least one of the magnetization vector components should non-linearly depend on the field strength. The direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is later applied."},
	{"MatStd", radia_MatStd, METH_VARARGS, "MatStd(name,mr) creates a pre-defined magnetic material with an optional absolute remanent magnetization mr. The material is identified by its name/formula, which can be one of the following:\n'NdFeB' for NdFeB permanent magnet material (default mr = 1.2 T);\n'SmCo5' for SmCo5 permanent magnet material (default mr = 0.85 T);\n'Sm2Co17' for Sm2Co17 permanent magnet material (default mr = 1.05 T);\n'Ferrite' for Ferrite permanent magnet material (default mr = 0.35 T);\n'Xc06' for an inexpensive Low Carbon Steel with C<0.06% (AFNOR);\n'Steel37' for an inexpensive Steel with C<0.13%;\n'Steel42' for an inexpensive Steel with C<0.19%;\n'AFK502' for a Vanadium Permendur type material from MetalImphy (Fe:49%, Co:49%, V:2%) similar to Vacoflux50 from VacuumSchmelze;\n'AFK1' for an inexpensive FeCo alloy from MetalImphy (Fe:74.2%, Co:25%, Cr:0.3%, Mn:0.5%)."},
	{"MatApl", radia_MatApl, METH_VARARGS, "MatApl(obj,mat) applies material mat to object obj."},
	{"MatMvsH", radia_MatMvsH, METH_VARARGS, "MatMvsH(obj,'mx|my|mz'|'',[hx,hy,hz]) computes magnetization from magnetic field strength vector [hx,hy,hz] for the material of the object obj; the magnetization components are specified by the second argument."},

	{"RlxPre", radia_RlxPre, METH_VARARGS, "RlxPre(obj,srcobj:0) builds an interaction matrix for the object obj, treating the object srcobj as additional external field source."},
	{"RlxMan", radia_RlxMan, METH_VARARGS, "RlxMan(intrc,meth,iternum,rlxpar) executes manual relaxation procedure for interaction matrix intrc using method number meth, by making iternum iterations with relaxation parameter value rlxpar."},
	{"RlxAuto", radia_RlxAuto, METH_VARARGS, "RlxAuto(intrc,prec,maxiter,meth:4,'ZeroM->True|False') executes automatic relaxation procedure with the interaction matrix intrc using the method number meth. Relaxation stops whenever the change in magnetization (averaged over all sub-elements) between two successive iterations is smaller than prec or the number of iterations is larger than maxiter. The option value 'ZeroM->True' (default) starts the relaxation by setting the magnetization values in all paricipating objects to zero; 'ZeroM->False' starts the relaxation with the existing magnetization values in the sub-volumes."},
	{"RlxUpdSrc", radia_RlxUpdSrc, METH_VARARGS, "RlxUpdSrc(intrc) updates external field data for the relaxation (to take into account e.g. modification of currents in coils, if any) without rebuilding the interaction matrix."},
	{"Solve", radia_Solve, METH_VARARGS, "Solve(obj,prec,maxiter,meth:4) solves a magnetostatic problem, i.e. builds an interaction matrix for the object obj and performs a relaxation procedure using the method number meth (default is 4). The relaxation stops whenever the change in magnetization (averaged over all sub-elements) between two successive iterations is smaller than prec or the number of iterations is larger than maxiter."},

	{"Fld", radia_Fld, METH_VARARGS,  "Fld(obj,'bx|by|bz|hx|hy|hz|ax|ay|az|mx|my|mz'|'',[x,y,z]|[[x1,y1,z1],[x2,y2,z2],...]) computes magnetic field created by the object obj in point(s) {x,y,z} ({x1,y1,z1},{x2,y2,z2},...). The field component is specified by the second input variable. The function accepts a list of 3D points of arbitrary nestness: in this case it returns the corresponding list of magnetic field values."},
	{"FldLst", radia_FldLst, METH_VARARGS,  "FldLst(obj,'bx|by|bz|hx|hy|hz|ax|ay|az|mx|my|mz'|'',[x1,y1,z1],[x2,y2,z2],np,'arg|noarg':'noarg',strt:0.) computes magnetic field created by object obj in np equidistant points along a line segment from [x1,y1,z1] to [x2,y2,z2]; the field component is specified by the second input variable; the 'arg|noarg' string variable specifies whether to output a longitudinal position for each point where the field is computed, and strt gives the start-value for the longitudinal position."},
	{"FldInt", radia_FldInt, METH_VARARGS, "FldInt(obj,'inf|fin','ibx|iby|ibz'|'',[x1,y1,z1],[x2,y2,z2]) computes magnetic field induction integral produced by the object obj along a straight line specified by points [x1,y1,z1] and [x2,y2,z2]; depending on the second variable value, the integral is infinite ('inf') or finite from [x1,y1,z1] to [x2,y2,z2] ('fin'); the field integral component is specified by the third input variable. The units are Tesla x millimeters."},
	{"FldPtcTrj", radia_FldPtcTrj, METH_VARARGS, "FldPtcTrj(obj,E,[x0,dxdy0,z0,dzdy0],[y0,y1],np) computes transverse coordinates and its derivatives (angles) of a relativistic charged particle trajectory in 3D magnetic field produced by the object obj, using the 4th order Runge-Kutta integration. The particle energy is E [GeV], initial transverse coordinates and derivatives are {x0,dxdy0,z0,dzdy0}; the longitudinal coordinate y is varied from y0 to y1 in np steps. All positions are in millimeters and angles in radians."},
	{"FldEnr", radia_FldEnr, METH_VARARGS, "FldEnr(objdst,objsrc) or radFldEnr(objdst,objsrc,[kx,ky,kz]) computes potential energy (in Joule) of the object objdst in the field created by the object objsrc. The first form of the function performes the computation based on absolute accuracy value for the energy (by default 10 Joule; can be modified by the function radFldCmpPrc). The second form performs the computation based on the destination object subdivision numbers [kx,ky,kz]."},
	{"FldEnrFrc", radia_FldEnrFrc, METH_VARARGS, "FldEnrFrc(objdst,objsrc,'fx|fy|fz'|'') or radFldEnrFrc(objdst,objsrc,'fx|fy|fz'|'',[kx,ky,kz]) computes force (in Newton) acting on the object objdst in the field produced by the object objsrc. The first form of the function performes the computation based on absolute accuracy value for the force (by default 10 Newton; can be modified by the function FldCmpPrc). The second form performs the computation based on the destination object subdivision numbers [kx,ky,kz]."},
	{"FldEnrTrq", radia_FldEnrTrq, METH_VARARGS, "radFldEnrTrq(objdst,objsrc,'tx|ty|tz'|'',[x,y,z]) or radFldEnrTrq[objdst,objsrc,'tx|ty|tz'|'',[x,y,z],[kx,ky,kz]] computes torque (in Newton*mm) with respect to point [x,y,z], acting on the object objdst in the field produced by the object objsrc. The first form of the function performes the computation based on absolute accuracy value for the torque (by default 10 Newton*mm; can be modified by the function FldCmpPrc). The second form performs the computation based on the destination object subdivision numbers [kx,ky,kz]."},
	{"FldFrc", radia_FldFrc, METH_VARARGS, "FldFrc(obj,shape) computes force of the field produced by the object obj acting onto objects within / after a shape defined by shape. shape can be the result of ObjRecMag() (parallelepiped) or FldFrcShpRtg() (rectangular surface). This function uses the algorithm based on Maxwell tensor, which may not always provide high efficiency. We suggest to use the function FldEnrFrc instead of this function."},
	{"FldFrcShpRtg", radia_FldFrcShpRtg, METH_VARARGS, "FldFrcShpRtg([x,y,z],[wx,wy]) creates a rectangle with central point [x,y,z] and dimensions [wx,wy], to be used for force computation via Maxwell tensor."},
	{"FldFocPot", radia_FldFocPot, METH_VARARGS, "FldFocPot(obj,[x1,y1,z1],[x2,y2,z2],np) computes \"focusing potential\" for trajectory of relativistic charged particle in magnetic field produced by the object obj. The integration is made from [x1,y1,z1] to [x2,y2,z2] with np equidistant points."},
	{"FldFocKickPer", radia_FldFocKickPer, METH_VARARGS, "FldFocKickPer(obj,[x1,y1,z1],[nsx,nsy,nsz],per,nper,[n1x,n1y,n1z],r1,np1,r2,np2,com:'',[nh:1,nps:8,d1:0,d2:0],'T2m2|rad|microrad':'T2m2',en:1,'fix|tab':'fix') computes matrices of 2nd order kicks of trajectory of relativistic charged particle in periodic magnetic field produced by the object obj. The longitudinal integration along one period starts at point [x1,y1,z1] and is done along direction pointed by vector [nsx,nsy,nsz]; per is period length, nper is number of full periods; one direction of the transverse grid is pointed by vector [n1x,n1y,n1z], the other transverse direction is given by vector product of [n1x,n1y,n1z] and [nsx,nsy,nsz]; r1 and r2 are ranges of the transverse grid, np1 and np2 are corresponding numbers of points; com is arbitrary string comment; nh is maximum number of magnetic field harmonics to treat (default 1), nps is number of longitudinal points (default 8), d1 and d2 are steps of transverse differentiation (by default equal to the steps of the transverse grid); the 'T2m2|rad|microrad' string variable specifies the units for the resulting 2nd order kick values (default 'T2m2'); en is electron elergy in GeV (optional, required only if units are 'rad' or 'microrad'); the 'fix|tab' string variable specifies the format of the output data string (i.e. element [5] of the output list), 'fix' for fixed-width (default), 'tab' for tab-delimited. Returns list containing: [0]- matrix of kick values in the first transverse direction, [1]- matrix of kick values in the second transverse direction, [2]- matrix of longitudinally-integrated squared transverse magnetic field calculated on same transverse mesh as kicks, [3],[4]- lists of positions defining the transverse grid, [5]- formatted string containing the computed results (for saving into a text file)."},
	{"FldCmpCrt", radia_FldCmpCrt, METH_VARARGS, "FldCmpCrt(prcB,prcA,prcBint,prcFrc,prcTrjCrd,prcTrjAng) sets general absolute accuracy levels for computation of field induction (prcB), vector potential (prcA), induction integrals along straight line (prcBint), field force (prcFrc), relativistic particle trajectory coordinates (prcTrjCrd) and angles (prcTrjAng)."},
	{"FldCmpPrc", radia_FldCmpPrc, METH_VARARGS, "FldCmpPrc('PrcB->prb,PrcA->pra,PrcBInt->prbint,PrcForce->prfrc,PrcTorque->prtrq,PrcEnergy->pre,PrcCoord->prcrd,PrcAngle->prang') sets general absolute accuracy levels for computation of magnetic field induction, vector potential, induction integral along straight line, field force, torque, energy; relativistic charged particle trajectory coordinates and angles. The function works in line with the Mathematica mechanism of Options. PrcB, PrcA, PrcBInt, PrcForce, PrcTorque, PrcEnergy, PrcCoord, PrcAngle are names of the options; prb, pra, prbint, prfrc, prtrq, pre, prcrd, prang are the corresponding values (real numbers specifying the accuracy levels)."},
	{"FldUnits", radia_FldUnits, METH_VARARGS, "FldUnits() shows the physical units currently in use."},
	{"FldLenRndSw", radia_FldLenRndSw, METH_VARARGS, "FldLenRndSw('on|off') switches on or off the randomization of all the length values. The randomization magnitude can be set by the function FldLenTol."},
	{"FldLenTol", radia_FldLenTol, METH_VARARGS, "FldLenTol(abs,rel,zero:0) sets absolute and relative randomization magnitudes for all the length values, including coordinates and dimensions of the objects producing magnetic field, and coordinates of points where the field is computed. Optimal values of the variables can be: rel=10^(-11), abs=L*rel, zero=abs, where L is the distance scale value (in mm) for the problem to be solved. Too small randomization magnitudes can result in run-time code errors."},
	{"FldShimSig", radia_FldShimSig, METH_VARARGS, "FldShimSig(obj,'bx|by|bz|hx|hy|hz|ibx|iby|ibz'|'',[dx,dy,dz],[x1,y1,z1],[x2,y2,z2],np,[vix,viy,viz]:[0,0,0]) computes virtual 'shim signature', i.e. variation of magnetic field component defined by the second variable, introduced by displacement [dx,dy,dz] of magnetic field source object obj. The field variation is computed at np equidistant points along a line segment from [x1,y1,z1] to [x2,y2,z2]; the vector [vix,viy,viz] is taken into account if a field integral variation should be computed: in this case, it defines orientation of the integration line."},

	{"UtiDmp", radia_UtiDmp, METH_VARARGS, "UtiDmp(elem,'asc|bin':'asc') outputs information about elem, which can be either one element or list of elements; second argument specifies whether the output should be in ASCII ('asc', default) or in Binary ('bin') format."},
	{"UtiDmpPrs", radia_UtiDmpPrs, METH_VARARGS, "UtiDmpPrs(bstr) parses byte-string bstr produced by UtiDmp(elem,'bin') and attempts to instantiate elem objects(s); returns either index of one instantiated object (if elem was an index of one object) or a list of indexes of instantiated objects (if elem was a list of objects)."},
	{"UtiDel", radia_UtiDel, METH_VARARGS, "UtiDel(elem) deletes element elem."},
	{"UtiDelAll", radia_UtiDelAll, METH_VARARGS, "UtiDelAll() deletes all previously created elements."},
	{"UtiVer", radia_UtiVer, METH_VARARGS, "UtiVer() returns version number of the Radia library."},
	{"UtiMPI", radia_UtiMPI, METH_VARARGS, "UtiMPI('on|off|share',data,rankFrom,rankTo) initializes (if argument is 'on') or finalizes (in argument is 'off') the Message Passing Inteface (MPI) for parallel calculations and returns list of basic MPI process parameters (in the case of initialization): rank of a process and total number of processes. In the case of first argument is 'share', the function will send data (list or array) from rankFrom (by default 0) to all processes (by default) of to rankTo."},

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
