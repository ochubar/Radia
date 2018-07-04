/*-------------------------------------------------------------------------
*
* File name:      radsbdvp.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 subdivided polyhedron with constant magnetization
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radsbdvp.h"
#include "radg3dgr.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTg3dGraphPresent* radTSubdividedPolyhedron::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTSubdivPolyhedronGraphPresent((radTGroup*)this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

void radTSubdividedPolyhedron::Dump(std::ostream& o, int ShortSign)
{
	((radTg3d*)((radTGroup*)this))->radTg3d::Dump(o, ShortSign);

	o << "Subdivided Polyhedron";

	if(ShortSign==1) return;

	o << endl;
	radTPolyhedron::DumpPureObjInfo(o, ShortSign);
	DumpMaterApplied(o);

	o << endl;
	radTGroup::DumpPureObjInfo(o, ShortSign);

	o << endl;
	((radTg3d*)((radTGroup*)this))->radTg3d::DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied (incl. the content): " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTSubdividedPolyhedron::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//Dumping objects that may be used by this object
	vector<pair<int, int> > vTrfKeys;
	//DumpBin_g3d_TreatTrfs(oStr, mEl, vTrfKeys);
	radTGroup::DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys); //all G3D members should be accessed via radTGroup
	
	vector<int> vGroupMemKeys;
	DumpBin_Group_TreatMembers(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vGroupMemKeys);

	int matKey=0;
	DumpBin_g3dRelax_TreatMat(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, matKey);

	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)radTGroup::Type_g();
	oStr << (char)radTGroup::Type_g3d();
	oStr << (char)Type_Group();
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTg3d
	radTGroup::DumpBin_g3d(oStr, vTrfKeys);

	//Members of radTGroup
	DumpBin_Group_OutMemKeys(oStr, vGroupMemKeys);

	//Members of radTg3dRelax
	DumpBin_g3dRelax(oStr, matKey);

	//Members of radTPolyhedron
	DumpBin_Polyhedron(oStr);

	//Members of radTSubdividedPolyhedron
	//int AmOfSubElem;
	oStr << AmOfSubElem;
}

//-------------------------------------------------------------------------
