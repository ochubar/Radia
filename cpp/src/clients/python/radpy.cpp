/************************************************************************//**
 * File: radpy.cpp
 * Description: Python binding
 * Project: Radia
 * First release: June 2018
 *
 * @author O.Chubar
 * @version 0.01
 ***************************************************************************/

#include "radentry.h"
#include "pyparse.h"

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
 * Magnetic Field Source: Rectangular Parallelepiped with Constant Magnetizatiom over volume
 ***************************************************************************/
static PyObject* radia_ObjRecMag(PyObject *self, PyObject *args)
{//The parallelepiped block is defined through its center point P[3], dimensions L[3], and magnetization M[3].

	PyObject *oP=0, *oL=0, *oM=0, *oResInd=0;

	try
	{
		if(!PyArg_ParseTuple(args, "OOO:ObjRecMag", &oP, &oL, &oM)) throw CombErStr(strEr_BadFuncArg, ": ObjRecMag");
		if((oP == 0) || (oL == 0) || (oM == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjRecMag");

		double arP[3], arL[3], arM[3];
		int lenP = 3;
		double *p = arP;
		CPyParse::CopyPyListElemsToNumArray(oP, 'd', p, lenP);
		p = arL;
		CPyParse::CopyPyListElemsToNumArray(oL, 'd', p, lenP);
		p = arM;
		CPyParse::CopyPyListElemsToNumArray(oM, 'd', p, lenP);
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

		double arP[3], arL[3], arM[3], arK[10], arRGB[3];
		int lenP = 3, lenK = 10;
		double *p = arP;
		CPyParse::CopyPyListElemsToNumArray(oP, 'd', p, lenP);
		p = arL;
		CPyParse::CopyPyListElemsToNumArray(oL, 'd', p, lenP);
		p = arM;
		CPyParse::CopyPyListElemsToNumArray(oM, 'd', p, lenP);
		p = arK;
		CPyParse::CopyPyListElemsToNumArray(oK, 'd', p, lenK);
		p = arRGB;
		CPyParse::CopyPyListElemsToNumArray(oRGB, 'd', p, lenP);

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
 * Magnetic Field Source: Recetrack conductor with rectangular cross-secton and constant current density over the volume
 ***************************************************************************/
static PyObject* radia_ObjRaceTrk(PyObject *self, PyObject *args)
{//The coil consists of four 90-degree bents connected by four straight parts parallel to the XY plane.

	PyObject *oP=0, *oR=0, *oL=0, *o_man_auto=0, *o_a=0, *oResInd=0;

	try
	{
		double h=0, curDens=0;
		int nseg=0;
		if(!PyArg_ParseTuple(args, "OOOdid|OO:ObjRaceTrk", &oP, &oR, &oL, &h, &nseg, &curDens, &o_man_auto, &o_a)) throw CombErStr(strEr_BadFuncArg, ": ObjRaceTrk");
		if((oP == 0) || (oR == 0) || (oL == 0)) throw CombErStr(strEr_BadFuncArg, ": ObjRaceTrk");

		double arP[3], arR[2], arL[2];
		int lenP = 3, lenR = 2, lenL = 2;
		double *p = arP;
		CPyParse::CopyPyListElemsToNumArray(oP, 'd', p, lenP);
		p = arR;
		CPyParse::CopyPyListElemsToNumArray(oR, 'd', p, lenR);
		p = arL;
		CPyParse::CopyPyListElemsToNumArray(oL, 'd', p, lenL);

		char sManAuto[256], sOrnt[256];
		sManAuto[0] = 'm'; sOrnt[0] = 'z';
		if(o_man_auto != 0) CPyParse::CopyPyStringToC(o_man_auto, sManAuto, 256);
		if(o_a != 0) CPyParse::CopyPyStringToC(o_a, sOrnt, 256);

		int ind = 0;
		g_pyParse.ProcRes(RadObjRaceTrk(&ind, arP, arR, arL, h, nseg, *sManAuto, *sOrnt, curDens));

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
		CPyParse::CopyPyListElemsToNumArray(oInds, 'i', arInds, nInds);

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
		int lenRGB = 3;
		double *p = arRGB;
		CPyParse::CopyPyListElemsToNumArray(oRGB, 'd', p, lenRGB);

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
 * Magnetic Field Sources:  Starts an application for viewing of 3D geometry of the object obj (the viewer is based on the GLUT / OpenGL graphics library)
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
		int lenP = 3;
		double *p = arP;
		CPyParse::CopyPyListElemsToNumArray(oP, 'd', p, lenP);
		p = arN;
		CPyParse::CopyPyListElemsToNumArray(oN, 'd', p, lenP);

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
		int lenP = 3;
		double *p = arP;
		CPyParse::CopyPyListElemsToNumArray(oP, 'd', p, lenP);
		p = arN;
		CPyParse::CopyPyListElemsToNumArray(oN, 'd', p, lenP);

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
		int lenH = 3;
		double *p = arH;
		CPyParse::CopyPyListElemsToNumArray(oH, 'd', p, lenH);

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
		int ind=0, numIt=0, meth=0;
		double prec=0;
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
		int lenP = 3;
		double *p = arP1;
		CPyParse::CopyPyListElemsToNumArray(oP1, 'd', p, lenP);
		p = arP2;
		CPyParse::CopyPyListElemsToNumArray(oP2, 'd', p, lenP);

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
	{"ObjFullMag", radia_ObjFullMag, METH_VARARGS, "ObjFullMag() instantiates Rectangular Parallelepiped with Constant Magnetizatiom over volume, subdivided, and added to a group, with some color assigned"},
	{"ObjRaceTrk", radia_ObjRaceTrk, METH_VARARGS, "ObjRaceTrk() instantiates Recetrack conductor with rectangular cross-secton and constant current density over the volume"},
	{"ObjCnt", radia_ObjCnt, METH_VARARGS, "ObjCnt() instantiates Container of Magnetic Field Sources"},

	{"ObjDrwAtr", radia_ObjDrwAtr, METH_VARARGS, "ObjDrwAtr() assigns drawing attributes - RGB color (r,g,b) and line thickness thcn - to a Magnetic Field Source object"},
	{"ObjDrwOpenGL", radia_ObjDrwOpenGL, METH_VARARGS, "ObjDrwOpenGL() assigns drawing attributes - RGB color (r,g,b) and line thickness thcn - to a Magnetic Field Source object"},

	{"TrfZerPara", radia_TrfZerPara, METH_VARARGS, "TrfZerPara() creates an object mirror with respect to a plane. The object mirror possesses the same geometry as obj, but its magnetization and/or current densities are modified in such a way that the magnetic field produced by the obj and its mirror in the plane of mirroring is perpendicular to this plane"},
	{"TrfZerPerp", radia_TrfZerPerp, METH_VARARGS, "TrfZerPerp() creates an object mirror with respect to a plane. The object mirror possesses the same geometry as obj, but its magnetization and/or current densities are modified in such a way that the magnetic field produced by the obj and its mirror in the plane of mirroring is parallel to this plane"},

	{"MatStd", radia_MatStd, METH_VARARGS, "MatStd() creates a pre-defined magnetic material (the material is identified by its name/formula, e.g. \"NdFeB\")"},
	{"MatSatIsoTab", radia_MatSatIsoTab, METH_VARARGS, "MatSatIsoTab() creates a nonlinear isotropic magnetic material with the M versus H curve defined by the list of pairs corresponding values of H and M [[H1,M1],[H2,M2],...]"},
	{"MatMvsH", radia_MatMvsH, METH_VARARGS, "MatMvsH() computes magnetization from magnetic field strength vector for a given material"},

	{"Solve", radia_Solve, METH_VARARGS, "Solve() solves a magnetostatic problem, i.e. builds an interaction matrix and performs a relaxation procedure"},

	{"Fld", radia_Fld, METH_VARARGS, "Fld() computes field created by the object obj at one or many points"},
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
