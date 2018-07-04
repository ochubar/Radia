/*-------------------------------------------------------------------------
*
* File name:      radopnam.h
*
* Project:        RADIA
*
* Description:    Option names for Radia / Mathematica functions
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADOPNAM_H
#define __RADOPNAM_H

#include <string.h>
#include <vector>
#include <map>

//-------------------------------------------------------------------------

struct radTOptionNames {
	char B[25], A[25], BInt[25], Force[25], Torque[25], Energy[25], Coord[25], Angle[25]; // Precisions

	char Frame[25], FrameValues[3][25];
	char SubdParamCode[25], SubdParamBorderCode[25], SubdParamCodeValues[2][25];
	char SubdCoils[25], SubdCoilsValues[4][25];
	char TriAngMin[25], TriAreaMax[25], TriExtOpt[25]; //ki->Numb|Size,TriAngMin->...,TriAreaMax->...,TriExtOpt->\"...\"

	char ShowLines[25], ShowLinesValues[4][25];
	char ShowFaces[25], ShowFacesValues[4][25];
	char ShowFrameAxes[25], ShowFrameAxesValues[4][25];
	char ShowSymChilds[25], ShowSymChildsValues[4][25];

	char FreeSym[25], FreeSymValues[4][25];
	char ZeroM[25], ZeroM_Values[4][25];

	char LinTreat[25]; //, LinCoefValues[4][25];
	char Debug[25];

	//map<string, vector<pair<string, int> > > mOptData;
	map<string, map<string, int> > mOptData;

	radTOptionNames()
	{
		//partially obsolete definitions:
		strcpy(B, "PrcB");
		strcpy(A, "PrcA");
		strcpy(BInt, "PrcBInt");
		strcpy(Force, "PrcForce");
		strcpy(Torque, "PrcTorque");
		strcpy(Energy, "PrcEnergy");
		strcpy(Coord, "PrcCoord");
		strcpy(Angle, "PrcAngle");

		strcpy(Frame, "Frame");
		strcpy(FrameValues[0], "Loc");
		strcpy(FrameValues[1], "LabTot");
		strcpy(FrameValues[2], "Lab");

		strcpy(SubdParamCode, "kxkykz");
		strcpy(SubdParamBorderCode, "ki");
		strcpy(SubdParamCodeValues[0], "Numb");
		strcpy(SubdParamCodeValues[1], "Size");

		strcpy(SubdCoils, "DivCoils");
		strcpy(SubdCoilsValues[0], "No");
		strcpy(SubdCoilsValues[1], "Yes");
		strcpy(SubdCoilsValues[2], "False");
		strcpy(SubdCoilsValues[3], "True");

		strcpy(ShowLines, "EdgeLines");
		strcpy(ShowLinesValues[0], "No");
		strcpy(ShowLinesValues[1], "Yes");
		strcpy(ShowLinesValues[2], "False");
		strcpy(ShowLinesValues[3], "True");

		strcpy(ShowFaces, "Faces");
		strcpy(ShowFacesValues[0], "No");
		strcpy(ShowFacesValues[1], "Yes");
		strcpy(ShowFacesValues[2], "False");
		strcpy(ShowFacesValues[3], "True");

		strcpy(ShowFrameAxes, "Axes");
		strcpy(ShowFrameAxesValues[0], "No");
		strcpy(ShowFrameAxesValues[1], "Yes");
		strcpy(ShowFrameAxesValues[2], "False");
		strcpy(ShowFrameAxesValues[3], "True");

		strcpy(ShowSymChilds, "ShowSym");
		strcpy(ShowSymChildsValues[0], "No");
		strcpy(ShowSymChildsValues[1], "Yes");
		strcpy(ShowSymChildsValues[2], "False");
		strcpy(ShowSymChildsValues[3], "True");

		strcpy(FreeSym, "FreeSym");
		strcpy(FreeSymValues[0], "No");
		strcpy(FreeSymValues[1], "Yes");
		strcpy(FreeSymValues[2], "False");
		strcpy(FreeSymValues[3], "True");

		strcpy(ZeroM, "ZeroM");
		strcpy(ZeroM_Values[0], "No");
		strcpy(ZeroM_Values[1], "Yes");
		strcpy(ZeroM_Values[2], "False");
		strcpy(ZeroM_Values[3], "True");

		strcpy(TriAngMin, "TriAngMin");
		strcpy(TriAreaMax, "TriAreaMax");
		strcpy(TriExtOpt, "TriExtOpt");

		strcpy(LinTreat, "Lin");
		//strcpy(LinCoefValues[0], "Rel");
		//strcpy(LinCoefValues[1], "Loc");
		//strcpy(LinCoefValues[2], "Abs");
		//strcpy(LinCoefValues[3], "Lab");

		strcpy(Debug, "Debug");

		//newer definitions
		map<string, int> vRealVal;
		vRealVal["d"] = 0;
		mOptData[B] = vRealVal;
		mOptData[A] = vRealVal;
		mOptData[BInt] = vRealVal;
		mOptData[Force] = vRealVal;
		mOptData[Torque] = vRealVal;
		mOptData[Energy] = vRealVal;
		mOptData[Coord] = vRealVal;
		mOptData[Angle] = vRealVal;
		mOptData[TriAngMin] = vRealVal;
		mOptData[TriAreaMax] = vRealVal;

		map<string, int> vStringVal;
		vStringVal["s"] = 0;
		mOptData[TriExtOpt] = vStringVal;

		map<string, int> vFrameVal;
		vFrameVal["Loc"] = 0;
		vFrameVal["LabTot"] = 1;
		vFrameVal["Lab"] = 2;
		mOptData[Frame] = vFrameVal;

		map<string, int> vSubdPar;
		vSubdPar["Numb"] = 0;
		vSubdPar["Size"] = 1;
		mOptData[SubdParamCode] = vSubdPar;
		mOptData[SubdParamBorderCode] = vSubdPar;

		map<string, int> vNoYes;
		vNoYes["No"] = 0;
		vNoYes["Yes"] = 1;
		vNoYes["False"] = 0;
		vNoYes["True"] = 1;
		mOptData[SubdCoils] = vNoYes;
		mOptData[ShowLines] = vNoYes;
		mOptData[ShowFaces] = vNoYes;
		mOptData[ShowFrameAxes] = vNoYes;
		mOptData[ShowSymChilds] = vNoYes;
		mOptData[FreeSym] = vNoYes;
		mOptData[ZeroM] = vNoYes;
		mOptData[Debug] = vNoYes;

		map<string, int> vRelAbs;
		vRelAbs["Rel"] = 0;
		vRelAbs["Abs"] = 1;
		vRelAbs["Loc"] = 0;
		vRelAbs["Lab"] = 1;
		mOptData[LinTreat] = vRelAbs;
	}

	bool parseOption(const char* optName, const char* optValueIn, char& cOptValueOut, double& dOptValueOut, char* sOptValueOut, char& outCode)
	{
		map<string, map<string, int> >::const_iterator it = mOptData.find(optName);
		if(it == mOptData.end()) return false;

		const map<string, int> &mapCurValues = it->second;
		int numCurValues = (int)mapCurValues.size();
		if(numCurValues <= 0) return false;

		outCode = 0;

		if(numCurValues == 1)
		{
			const char* c_strKey = mapCurValues.begin()->first.c_str();
			if((*c_strKey == 'd') || (*c_strKey == 'r'))
			{
				if(optValueIn != 0) dOptValueOut = atof(optValueIn);

				outCode = 'd';
				return true;
			}
			else if(*c_strKey == 's')
			{//return without parsing optValueIn
				if(optValueIn != 0)
				{
					if(sOptValueOut != 0) strcpy(sOptValueOut, optValueIn);
				}

				outCode = 's';
				return true;
			}
		}

		if(optValueIn != 0)
		{
			map<string, int>::const_iterator itVal = mapCurValues.find(optValueIn);
			if(itVal == mapCurValues.end()) return false;
			cOptValueOut = (char)(itVal->second);
		}

		outCode = 'c';
		return true;
	}

	bool findParseOptionValues(const char** arAllOptionNamesIn, const char** arAllOptionValuesIn, int numAllOptIn, const char** arOptNamesToFind, int numOptToFind, char* cArOptValsFoundParsed, double* dArOptValsFoundParsed, char** sArOptValsFoundParsed)
	{
		if((numAllOptIn == 0) || (arAllOptionNamesIn == 0) || (arAllOptionValuesIn == 0)) return true;
		if((arOptNamesToFind == 0) || (numOptToFind == 0)) return true;
		if((cArOptValsFoundParsed == 0) && (dArOptValsFoundParsed == 0) && (sArOptValsFoundParsed)) return true;

		char cAux, sAux[1024], cOutCode;
		double dAux;
		sAux[0] = '\0';

		for(int i=0; i<numAllOptIn; i++)
		{
			//if(cArOptValsFoundParsed != 0) p_c = cArOptValsFoundParsed + countOptFound;
			//if(dArOptValsFoundParsed != 0) p_d = dArOptValsFoundParsed + countOptFound;
			//if(sArOptValsFoundParsed != 0) p_s = sArOptValsFoundParsed[countOptFound];

			const char *sCurOptionNameIn = arAllOptionNamesIn[i];
			cOutCode = 0;
			if(!parseOption(sCurOptionNameIn, arAllOptionValuesIn[i], cAux, dAux, sAux, cOutCode)) return false;

			int cOptCount=-1, dOptCount=-1, sOptCount=-1;

			for(int j=0; j<numOptToFind; j++)
			{
				char auxOutCode=0;
				if(!parseOption(arOptNamesToFind[j], 0, cAux, dAux, sAux, auxOutCode)) return false;
				if(auxOutCode == 'c') cOptCount++;
				else if(auxOutCode == 'd') dOptCount++;
				else if(auxOutCode == 's') sOptCount++;

				if(strcmp(sCurOptionNameIn, arOptNamesToFind[j]) == 0)
				{
					if(cArOptValsFoundParsed != 0) 
					{
						//if(cOutCode == 'c') cArOptValsFoundParsed[j] = cAux;
						if(cOutCode == 'c') cArOptValsFoundParsed[cOptCount] = cAux;
					}
					if(dArOptValsFoundParsed != 0) 
					{
						//if(cOutCode == 'd') dArOptValsFoundParsed[j] = dAux;
						if(cOutCode == 'd') dArOptValsFoundParsed[dOptCount] = dAux;
					}
					if(sArOptValsFoundParsed != 0) 
					{
						//if(cOutCode == 's') strcpy(sArOptValsFoundParsed[j], sAux);
						if(cOutCode == 's') strcpy(sArOptValsFoundParsed[sOptCount], sAux);
					}
				}
			}
		}
		return true;
	}
};

//-------------------------------------------------------------------------

#endif
