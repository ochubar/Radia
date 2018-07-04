/*-------------------------------------------------------------------------
*
* File name:      radsbdep.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 subdivided extruded polygon (prism)
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radsbdep.h"
#include "radg3dgr.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTg3dGraphPresent* radTSubdividedExtrPolygon::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = 0;
	if(!AlgsBasedOnKsQsMayNotWork) g3dGraphPresentPtr = new radTSubdivExtrPolygGraphPresent((radTGroup*)this);
	else g3dGraphPresentPtr = new radTGroupGraphPresent((radTGroup*)this);
	g3dGraphPresentPtr->ShowInternalFacesAfterCut = false;
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

void radTSubdividedExtrPolygon::Dump(std::ostream& o, int ShortSign)
{
	((radTg3d*)((radTGroup*)this))->radTg3d::Dump(o, ShortSign);

	o << "Subdivided ThckPgn";

	if(ShortSign==1) return;

	o << endl;
	radTExtrPolygon::DumpPureObjInfo(o, ShortSign);
	DumpMaterApplied(o);

	o << endl;
	radTGroup::DumpPureObjInfo(o, ShortSign);

	o << endl;
	((radTg3d*)((radTGroup*)this))->radTg3d::DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied (incl. the content): " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTSubdividedExtrPolygon::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
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

	//Members of radTExtrPolygon
	DumpBin_ExtrPolygon(oStr);

	//Members of radTSubdividedExtrPolygon
	//int kx, ky, kz;
	oStr << kx << ky << kz;

	//double qx, qy, qz;
	oStr << qx << qy << qz;

	//int AmOfSubElem;
	oStr << AmOfSubElem;

	//bool AlgsBasedOnKsQsMayNotWork;
	oStr << AlgsBasedOnKsQsMayNotWork;
}

//-------------------------------------------------------------------------

radTSubdividedExtrPolygon::radTSubdividedExtrPolygon(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{
	//Members of radTg3d 
	//((radTRecMag*)this)->DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
	radTGroup::DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers); //all G3D members should be accessed via radTGroup

	//Members of radTGroup 
	DumpBinParse_Group(inStr, mKeysOldNew, gMapOfHandlers);

	//Members of radTg3dRelax 
	DumpBinParse_g3dRelax(inStr, mKeysOldNew, gMapOfHandlers);

	//Members of radTExtrPolygon 
	DumpBinParse_ExtrPolygon(inStr);

	//Members of radTSubdividedExtrPolygon
	//int kx, ky, kz;
	inStr >> kx;
	inStr >> ky;
	inStr >> kz;

	//double qx, qy, qz;
	inStr >> qx;
	inStr >> qy;
	inStr >> qz;

	//int AmOfSubElem;
	inStr >> AmOfSubElem;

	//bool AlgsBasedOnKsQsMayNotWork;
	inStr >> AlgsBasedOnKsQsMayNotWork;
}

//-------------------------------------------------------------------------
