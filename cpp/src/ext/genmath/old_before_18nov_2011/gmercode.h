
#ifndef __GMERCODE_H
#define __GMERCODE_H

//-------------------------------------------------------------------------
// Error codes that are eventually thown by GenMath routines
//-------------------------------------------------------------------------

#define GM_FIRST_ERR_CODE -5000000
#define GM_RK_MAX_NUM_STEPS_REACHED GM_FIRST_ERR_CODE + 1 // "Maximum number of step subdivisions reached at automatic Runge-Kutta integration.\0",
#define GM_RK_STEP_SIZE_TOO_SMALL GM_FIRST_ERR_CODE + 2 // "Step size is too small in automatic Runge-Kutta integration routine.\0",

#endif


