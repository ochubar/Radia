/*-------------------------------------------------------------------------
*
* File name:      radcast.h
*
* Project:        RADIA
*
* Description:    "Dynamic cast" for RADIA classes
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADCAST_H
#define __RADCAST_H

#include "radsend.h"
#include "radrec.h"
#include "radexpgn.h"
#include "radvlpgn.h"
#include "radgroup.h"
#include "radplnr.h"
#include "radflm.h"
#include "radtrans.h"
#include "radmater.h"
#include "radsbdrc.h"
#include "radsbdep.h"
#include "radsbdvp.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTInteraction;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTCast {
public:
	inline static radTg3d* g3dCast(radTg* gPtr);
	inline static radTGroup* GroupCast(radTg3d* g3dPtr);
	inline static radTg3dRelax* g3dRelaxCast(radTg3d* g3dPtr);
	inline static radTRecMag* RecMagCast(radTg3dRelax* g3dRelaxPtr);
	inline static radTSubdividedRecMag* SubdividedRecMagCast(radTGroup* GroupPtr);
	inline static radTSubdividedRecMag* SubdividedRecMagCastFromRelax(radTg3dRelax* g3dRelaxPtr);
	inline static radTExtrPolygon* ExtrPolygonCast(radTg3dRelax* g3dRelaxPtr);
	inline static radTSubdividedExtrPolygon* SubdExtrPolygonCastFromGroup(radTGroup* GroupPtr);
	inline static radTSubdividedExtrPolygon* SubdExtrPolygonCastFromRelax(radTg3dRelax* g3dRelaxPtr);
	inline static radTPolyhedron* PolyhedronCast(radTg3dRelax* g3dRelaxPtr);
	inline static radTSubdividedPolyhedron* SubdPolyhedronCastFromGroup(radTGroup* GroupPtr);

	inline static radTRectangle* RectangleCast(radTg3d* g3dPtr);
	inline static radTFlmLinCur* FlmLinCurCast(radTg3d* g3dPtr);
	inline static radTrans* TransCast(radTg* gPtr);
	inline static radTrans* IdentTransCast(radTrans* TransPtr);
	inline static radTMaterial* MaterCast(radTg* gPtr);
	inline static radTLinearAnisotropMaterial* LinAnisoMaterCast(radTMaterial* MaterPtr);
	inline static radTLinearIsotropMaterial* LinIsoMaterCast(radTMaterial* MaterPtr);
	inline static radTNonlinearIsotropMaterial* NonlinIsoMaterCast(radTMaterial* MaterPtr);
	static radTInteraction* InteractCast(radTg* gPtr);
};

//-------------------------------------------------------------------------

inline radTg3d* radTCast::g3dCast(radTg* gPtr)
{
	radTg3d g3d;
	if(gPtr->Type_g()==g3d.Type_g()) return (radTg3d*)gPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTGroup* radTCast::GroupCast(radTg3d* g3dPtr)
{
	radTGroup Group;
	if(g3dPtr->Type_g3d()==Group.Type_g3d()) return (radTGroup*)g3dPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTg3dRelax* radTCast::g3dRelaxCast(radTg3d* g3dPtr)
{
	radTg3dRelax g3dRelax;
	if(g3dPtr->Type_g3d()==g3dRelax.Type_g3d()) return (radTg3dRelax*)g3dPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTRecMag* radTCast::RecMagCast(radTg3dRelax* g3dRelaxPtr)
{
	radTRecMag RecMag;
	if(g3dRelaxPtr->Type_g3dRelax()==RecMag.Type_g3dRelax()) return (radTRecMag*)g3dRelaxPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTSubdividedRecMag* radTCast::SubdividedRecMagCast(radTGroup* GroupPtr)
{
	radTSubdividedRecMag SubdividedRecMag;
	if(GroupPtr->Type_Group()==SubdividedRecMag.Type_Group()) return (radTSubdividedRecMag*)GroupPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTSubdividedRecMag* radTCast::SubdividedRecMagCastFromRelax(radTg3dRelax* g3dRelaxPtr)
{
	radTSubdividedRecMag SubdividedRecMag;
	if(g3dRelaxPtr->Type_g3dRelax()==SubdividedRecMag.Type_g3dRelax()) return (radTSubdividedRecMag*)g3dRelaxPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTExtrPolygon* radTCast::ExtrPolygonCast(radTg3dRelax* g3dRelaxPtr)
{
	radTExtrPolygon ExtrPolygon;
	if(g3dRelaxPtr->Type_g3dRelax()==ExtrPolygon.Type_g3dRelax()) return (radTExtrPolygon*)g3dRelaxPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTSubdividedExtrPolygon* radTCast::SubdExtrPolygonCastFromGroup(radTGroup* GroupPtr)
{
	radTSubdividedExtrPolygon SubdividedExtrPolygon;
	if(GroupPtr->Type_Group()==SubdividedExtrPolygon.Type_Group()) return (radTSubdividedExtrPolygon*)GroupPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTSubdividedExtrPolygon* radTCast::SubdExtrPolygonCastFromRelax(radTg3dRelax* g3dRelaxPtr)
{
	radTSubdividedExtrPolygon SubdividedExtrPolygon;
	if(g3dRelaxPtr->Type_g3dRelax()==SubdividedExtrPolygon.Type_g3dRelax()) return (radTSubdividedExtrPolygon*)g3dRelaxPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTPolyhedron* radTCast::PolyhedronCast(radTg3dRelax* g3dRelaxPtr)
{
	radTPolyhedron Polyhedron;
	if(g3dRelaxPtr->Type_g3dRelax()==Polyhedron.Type_g3dRelax()) return (radTPolyhedron*)g3dRelaxPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTSubdividedPolyhedron* radTCast::SubdPolyhedronCastFromGroup(radTGroup* GroupPtr)
{
	radTSubdividedPolyhedron SubdividedPolyhedron;
	if(GroupPtr->Type_Group()==SubdividedPolyhedron.Type_Group()) return (radTSubdividedPolyhedron*)GroupPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTRectangle* radTCast::RectangleCast(radTg3d* g3dPtr)
{
	radTRectangle Rectangle;
	if(g3dPtr->Type_g3d()==Rectangle.Type_g3d()) return (radTRectangle*)g3dPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTFlmLinCur* radTCast::FlmLinCurCast(radTg3d* g3dPtr)
{
	radTFlmLinCur FlmLinCur;
	if(g3dPtr->Type_g3d()==FlmLinCur.Type_g3d()) return (radTFlmLinCur*)g3dPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTrans* radTCast::TransCast(radTg* gPtr)
{
	radTrans Trans;
	if(gPtr->Type_g()==Trans.Type_g()) return (radTrans*)gPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTrans* radTCast::IdentTransCast(radTrans* TransPtr)
{
	radIdentTrans IdentTrans;
	if(TransPtr->Type_Trans()==IdentTrans.Type_Trans()) return (radIdentTrans*)TransPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTMaterial* radTCast::MaterCast(radTg* gPtr)
{
	radTMaterial Mater;
	if(gPtr->Type_g()==Mater.Type_g()) return (radTMaterial*)gPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTLinearAnisotropMaterial* radTCast::LinAnisoMaterCast(radTMaterial* MaterPtr)
{
	radTLinearAnisotropMaterial LinAnisoMater;
	if(MaterPtr->Type_Material()==LinAnisoMater.Type_Material()) return (radTLinearAnisotropMaterial*)MaterPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTLinearIsotropMaterial* radTCast::LinIsoMaterCast(radTMaterial* MaterPtr)
{
	radTLinearIsotropMaterial LinIsoMater;
	if(MaterPtr->Type_Material()==LinIsoMater.Type_Material()) return (radTLinearIsotropMaterial*)MaterPtr;
	else return 0;
}

//-------------------------------------------------------------------------

inline radTNonlinearIsotropMaterial* radTCast::NonlinIsoMaterCast(radTMaterial* MaterPtr)
{
	radTNonlinearIsotropMaterial NonlinIsoMater;
	if(MaterPtr->Type_Material()==NonlinIsoMater.Type_Material()) return (radTNonlinearIsotropMaterial*)MaterPtr;
	else return 0;
}

//-------------------------------------------------------------------------

#endif
