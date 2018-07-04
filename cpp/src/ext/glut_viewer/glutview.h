
#ifndef __GLUTVIEW_H
#define __GLUTVIEW_H

#include "AlphaDll.h"

//-------------------------------------------------------------------------

#ifdef __cplusplus  
extern "C" {
#endif

/** Shows a 3D viewer of the scenes composed of polygons and lines.
@param VertexCoord [in] 2-element array of 3D coordinates of vertex points of polygons [0] and lines [1]
@param VertexNumbers [in] 2-element array of numbers of vertex points in polygons [0] and lines [1]
@param VertexIndexes [in] 2-element array of arrays of indexes of vertex points composing polygons [0] and lines [1]
@param PgnAndLineLenghs [in] 2-element array of arrays of numbers of vertex points in polygons [0] and lines [1]
@param PgnAndLineColors [in] 2-element array of arrays of RGB colors of polygons [0] and lines [1]
@param PgnAndLineNumbers [in] 2-element array of numbers of polygons [0] and lines [1]
@param WinTitle [in] string containing window title
@param StartMode [in] wiewer start mode: 0- single-thread mode (function returns when viewer window is closed), 1- multi-thread mode (function starts viewer window and returnes, the window may stay open)
@param ErrWarnText [out] error or warning text (buffer should be allocated by calling application)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL ViewScene3D(double** VertexCoord, int* VertexNumbers, int** VertexIndexes, int** PgnAndLineLenghs, float** PgnAndLineColors, int* PgnAndLineNumbers, char* WinTitle, char StartMode, char* ErrWarnText);

/** Shows a 3D viewer of the scenes composed of polygons and lines.
@param VertexCoord [in] array of 3D coordinates of vertex points of polygons
@param VertexNumbers [in] number of vertex points in polygons
@param VertexIndexes [in] arrays of indexes of vertex points composing the polygons
@param PgnLenghs [in] array of numbers of vertex points in polygons
@param PgnColors [in] array of RGB colors of polygons
@param PgnNumber [in] number of polygons
@param WinTitle [in] string containing window title
@param StartMode [in] wiewer start mode: 0- single-thread mode (function returns when viewer window is closed), 1- multi-thread mode (function starts viewer window and returnes, the window may stay open)
@param ErrWarnText [out] error or warning text (buffer should be allocated by calling application)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL ViewPolygons3D(double* VertexCoord, int VertexNumber, int* VertexIndexes, int* PgnLenghs, float* PgnColors, int PgnNumber, char* WinTitle, char StartMode, char* ErrWarnText);

/** Shows a 2D graph with one or several curves.
@param FuncValues [in] array of arrays of function values
@param ArgValues [in] array of arrays of function argument values (ignored if 0)
@param ArgStart [in] array of function argument start values (ignored if 0)
@param ArgStep [in] array of function argument step values (ignored if 0)
@param Size [in] array of lengths of the function (and argument) arrays
@param CurveOptions [in] array of arrays of numbers describing curve options:
		[0]- red color intensity (bw 0 and 1)
		[1]- green color intensity (bw 0 and 1)
		[2]- blue color intensity (bw 0 and 1)
		[3]- line thickness
		[4]- line style
@param CurveNumber [in] number of curves to plot
@param Units [in] 4-element array of c-strings describing data units
@param Labels [in] 2-element array of c-strings describing data labels
@param GraphOptions [in] array of numbers describing the graph options:
		[0]- font size of abscissa numbers
		[1]- font size of ordinate numbers
		[2]- font size of the bottom label
		[3]- font size of the left label
		[4]- font size of the top label
		[5]- font size of the right label
		[6]- show or not ascissa grid lines (1 or 0)
		[7]- show or not ordinate grid lines (1 or 0)
		[8]- relative tick length (bw 0 and 1)
@param WinTitle [in] string containing window title
@param StartMode [in] wiewer start mode: 0- single-thread mode (function returns when viewer window is closed), 1- multi-thread mode (function starts viewer window and returnes, the window may stay open)
@param ErrWarnText [out] error or warning text (buffer should be allocated by calling application)
@author O.C.
*/
EXP int CALL ViewPlot2D(double** FuncValues, double** ArgValues, double* ArgStart, double* ArgStep, long* Size, double** CurveOptions, int CurveNumber, char** Units, char** Labels, double* GraphOptions, char* WinTitle, char StartMode, char* ErrWarnText);

#ifdef __cplusplus  
}
#endif

//-------------------------------------------------------------------------

#endif
