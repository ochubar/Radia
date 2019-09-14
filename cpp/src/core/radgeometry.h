/*-------------------------------------------------------------------------
*
* File name:      radgeometry.h
*
* Project:        RADIA
*
* Description:    Geometry Results
*
* Author(s):      Robert Nagler
*
* First release:  2019
* 
* Copyright (C):  2019 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADGEOMETRY_H
#define __RADGEOMETRY_H

struct radGeometryFlattened {
    int VerticesCount = 0;
    double *Vertices = 0;
    int Count = 0;
    int *Lengths = 0;
    float *Colors = 0;
    ~radGeometryFlattened() {
        // allocation is in PrepareGeomPolygDataForViewing
        delete[] Vertices;
        delete[] Lengths;
        delete[] Colors;
    }
};

struct radGeometry {
	struct radGeometryFlattened Polygons;
	struct radGeometryFlattened Lines;
};

//-------------------------------------------------------------------------

#endif // __RADGEOMETRY_H
