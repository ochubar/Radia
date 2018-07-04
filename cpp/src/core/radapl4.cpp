/*-------------------------------------------------------------------------
*
* File name:      radapl4.cpp
*
* Project:        RADIA
*
* Description:    Wrapping RADIA application function calls
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radappl.h"
#include "radg3da1.h"
#include "radrec.h"
#include "radplnr.h"
#include "radvlpgn.h"
#include "radcnvrg.h"
#include "radopnam.h"

#include <math.h>

//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------

//int radTApplication::SetUpPolyhedronsFromBaseFacePolygons(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, char frame, radThg& hgOut)
int radTApplication::SetUpPolyhedronsFromBaseFacePolygons(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInExtrSteps, char frame, radThg& hgOut)
{
	if((arPoints2d == 0) || (lenArPoints2d <= 0)) return 0;
	if((arPtrTrfParInExtrSteps == 0) || (arStrTrfOrderInExtrSteps == 0) || (arNumTrfInExtrSteps == 0) || (NumSteps <= 0)) return 0;

	try
	{
		radTGroup* pGroup = new radTGroup();
		if(pGroup == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(pGroup);

		radTPolygon *pBaseFacePgn1 = new radTPolygon(arPoints2d, lenArPoints2d), *pBaseFacePgn2=0;
		//pBaseFacePgn1->CoordZ = zc;
		pBaseFacePgn1->CoordZ = 0;

		TVector3d vCenPointRot1(0,0,0); //, vEZ(0,0,1), vBaseNorm(0,0,0); //OC160209 Test
		//if((*strOrient == 'x') || (*strOrient == 'X')) vBaseNorm.x = 1.;
		//else if((*strOrient == 'y') || (*strOrient == 'Y')) vBaseNorm.y = 1.;
		//else if((*strOrient == 'z') || (*strOrient == 'Z')) vBaseNorm.z = 1.;

		radTrans trfInitBaseOrient, trfInitBaseTransl;
		trfInitBaseOrient.SetupRotationToPermutAxes(vCenPointRot1, 'Z', *strOrient);
		//to find appropriate rotation taking into account direction of polygon normal !!!
		//trfInitBaseOrient.SetupRotation(vCenPointRot1, vEZ, vBaseNorm);

		TVector3d vCenPointTransl(0,0,zc);
		vCenPointTransl = trfInitBaseOrient.TrPoint(vCenPointTransl);
		trfInitBaseTransl.SetupTranslation(vCenPointTransl);

		//radTrans *pBaseFaceTrf1 = new radTrans();
		//TVector3d vCenPointRot1(0,0,zc);
		//pBaseFaceTrf1->SetupRotationToPermutAxes(vCenPointRot1, 'Z', *strOrient);
		//pBaseFaceTrf1->SetupIdent();

		radTrans auxIntermedTrans, auxTotTrans, auxTrans0;
		auxTrans0.SetupIdent();
		if(frame > 0)
		{//Lab frame
			auxTrans0.TrMult(trfInitBaseOrient, 'R');
			auxTrans0.TrMult(trfInitBaseTransl, 'L');
		}

		TVector3d vNormFace1(0,0,1), vCenPointBaseFace1(pBaseFacePgn1->CentrPoint.x, pBaseFacePgn1->CentrPoint.y, zc);
		int NumSteps_mi_1 = NumSteps - 1;

		for(int j=0; j<NumSteps; j++)
		{
			double **arPtrTrfParInCurExtrStep = arPtrTrfParInExtrSteps[j];
			char *arStrTrfOrderInCurExtrStep = arStrTrfOrderInExtrSteps[j];
			int numTrfInCurExtrStep = arNumTrfInExtrSteps[j];
			if((arPtrTrfParInCurExtrStep == 0) || (arStrTrfOrderInCurExtrStep == 0) || (numTrfInCurExtrStep <= 0)) continue;

			pBaseFacePgn2 = new radTPolygon(*pBaseFacePgn1);
			auxTotTrans.SetupIdent();

			for(int k=0; k<numTrfInCurExtrStep; k++)
			{
				double *pCurTrfPar = arPtrTrfParInCurExtrStep[k];
				if(pCurTrfPar == 0) continue;

				char cTypeCurTrf = arStrTrfOrderInCurExtrStep[k];
				if((cTypeCurTrf == 'r') || (cTypeCurTrf == 'R'))
				{
					TVector3d vRotCen(*pCurTrfPar, *(pCurTrfPar+1), *(pCurTrfPar+2));
					TVector3d vRotAxis(*(pCurTrfPar+3), *(pCurTrfPar+4), *(pCurTrfPar+5));
					auxIntermedTrans.SetupRotation(vRotCen, vRotAxis, *(pCurTrfPar+6));
					auxTotTrans.TrMult(auxIntermedTrans, 'L');
				}
				else if((cTypeCurTrf == 't') || (cTypeCurTrf == 'T'))
				{
					TVector3d vTransl(*pCurTrfPar, *(pCurTrfPar+1), *(pCurTrfPar+2));
					auxIntermedTrans.SetupTranslation(vTransl);
					auxTotTrans.TrMult(auxIntermedTrans, 'L');
				}
				else if((cTypeCurTrf == 'h') || (cTypeCurTrf == 'H'))
				{
					//any homothety is applied to a face before space transformations
					//to check if this should be improved
					pBaseFacePgn2->ApplyHomothety(*pCurTrfPar, *(pCurTrfPar+1), *(pCurTrfPar+2));
				}
			}

			radTrans *pBaseFaceTrf1 = new radTrans(auxTrans0);

			if(frame == 0) auxTrans0.TrMult(auxTotTrans, 'R'); //local frame
			else auxTrans0.TrMult(auxTotTrans, 'L'); //laboratory frame

			radTrans *pBaseFaceTrf2 = new radTrans(auxTrans0); //check for memory leak !!!

			//if(frame == 0)
			//{
			//	pBaseFaceTrf1->TrMult(trfInitBaseOrient, 'R');
			//	pBaseFaceTrf1->TrMult(trfInitBaseTransl, 'L');

			//	pBaseFaceTrf2->TrMult(trfInitBaseOrient, 'R');
			//	pBaseFaceTrf2->TrMult(trfInitBaseTransl, 'L');
			//}

			////else
			////{
			////	pBaseFaceTrf1->TrMult(trfInitBaseOrient, 'L');
			////	pBaseFaceTrf1->TrMult(trfInitBaseTransl, 'R');

			////	pBaseFaceTrf2->TrMult(trfInitBaseOrient, 'L');
			////	pBaseFaceTrf2->TrMult(trfInitBaseTransl, 'R');
			////}

			//this doesn't take into account direction of normals and eventual intersection of the faces
			radTHandle<radTPolygon> hPgn1(pBaseFacePgn1), hPgn2(pBaseFacePgn2);
			radTHandle<radTrans> hTrf1(pBaseFaceTrf1), hTrf2(pBaseFaceTrf2);
			radTHandlePgnAndTrans hPgnTrf1(hPgn1, hTrf1), hPgnTrf2(hPgn2, hTrf2);

			//create polyhedron here from pBaseFacePgn1, pBaseFaceTrf1, pBaseFacePgn2, pBaseFaceTrf2
			//radTPolyhedron *pCurStepPlhdr = new radTPolyhedron(hPgnTrf1, hPgnTrf2, avgCur);

			double *pMagComp = 0;
			if(arMagnCompInExtrSteps != 0) pMagComp = arMagnCompInExtrSteps + j*3;
			radTPolyhedron *pCurStepPlhdr = new radTPolyhedron(hPgnTrf1, hPgnTrf2, avgCur, pMagComp);

			if(pCurStepPlhdr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			if(pCurStepPlhdr->SomethingIsWrong) { delete pCurStepPlhdr; return 0;}
			radThg hgCurStepPlhdr(pCurStepPlhdr);
			pGroup->AddElement(AddElementToContainer(hgCurStepPlhdr), hgCurStepPlhdr);

			if(j < NumSteps_mi_1)
			{
				pBaseFacePgn1 = new radTPolygon(*pBaseFacePgn2);

			//	//- treat frame (0- local, otherwise "laboratory");
			//	//if(frame == 0) pBaseFaceTrf1 = new radTrans(*pBaseFaceTrf1); //leave base trf always the same if Frame->Loc ??
			//	//else pBaseFaceTrf1 = new radTrans(*pBaseFaceTrf2); //change base trf if Frame->Lab

			//	if(frame == 0) pBaseFaceTrf1 = new radTrans(*pBaseFaceTrf2); //change base trf if Frame->Loc
			//	else pBaseFaceTrf1 = new radTrans(*pBaseFaceTrf1); //leave base trf always the same if Frame->Lab
			}

/**
			int numRotAtStep = 0;
			if(RotationsAreDefined) 
			{
				numRotAtStep = arNumRotInExtrSteps[j];
				pRotParInExtrSteps = arPtrRotParInExtrSteps[j];
			}
			int j3 = j*3;
			auxTotTrans.SetupIdent();

			//to modify:
			//- update homothety (kx,ky,phi)

			//Calculate transformation for new base face, and new vertex points in case of homothety
			char *pFlatTrfOrderInExtrSteps = arFlatTrfOrderInExtrSteps + j3;
			for(int i=0; i<3; i++) //"3" means trying to find rotations ('r'), translations ('t') and homothety ('h')
			{
				if(*pFlatTrfOrderInExtrSteps == 'r')
				{
					for(int k=0; k<numRotAtStep; k++)
					{
						TVector3d vRotCen(pRotParInExtrSteps[0], pRotParInExtrSteps[1], pRotParInExtrSteps[2]);
						TVector3d vRotAxis(pRotParInExtrSteps[3], pRotParInExtrSteps[4], pRotParInExtrSteps[5]);
						auxIntermedTrans.SetupRotation(vRotCen, vRotAxis, pRotParInExtrSteps[6]);

						auxTotTrans.TrMult(auxIntermedTrans, 'L');
						pRotParInExtrSteps += 7;
					}
				}
				else if(*pFlatTrfOrderInExtrSteps == 't')
				{
					TVector3d vTransl(arTrslParInExtrSteps[j3], arTrslParInExtrSteps[j3 + 1], arTrslParInExtrSteps[j3 + 2]);
					auxIntermedTrans.SetupTranslation(vTransl);
					auxTotTrans.TrMult(auxIntermedTrans, 'L');
				}
				else if(*pFlatTrfOrderInExtrSteps == 'h')
				{
					//pBaseFacePgn2->ApplyHomothety(arHomParInExtrSteps[j]);
					pBaseFacePgn2->ApplyHomothety(arHomParInExtrSteps[j3], arHomParInExtrSteps[j3 + 1], arHomParInExtrSteps[j3 + 2]);
				}
				pFlatTrfOrderInExtrSteps++;
			}
**/
		}

		if(frame == 0)
		{
			auxTrans0.SetupIdent();
			auxTrans0.TrMult(trfInitBaseOrient, 'R');
			auxTrans0.TrMult(trfInitBaseTransl, 'L');
			
			double relTol = 1.E-12; //to ensure correct operation on 32-bit OS
			if(!auxTrans0.IsIdent(relTol))
			{
				radTrans* TransPtr = new radTrans(auxTrans0);
				if(TransPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
				radThg hgTrf(TransPtr);
				//int trfElemKey = AddElementToContainer(hgTrf);
				//if(trfElemKey != 0) hgLoc.rep->AddTransform(1, hgTrf);
				((radTg3d*)(hgLoc.rep))->AddTransform(1, hgTrf);
			}
		}

		hgOut = hgLoc;
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetUpPolyhedronsFromBaseFacePolygonsTri(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, TVector2d* arTriVertPt, int numTriVertPt, int* arTriVertInd, int numTri, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInExtrSteps, char frame, radThg& hgOut)
{
	if((arTriVertPt == 0) || (numTriVertPt <= 0) || (arTriVertInd == 0) || (numTri <= 0)) return 0;
	if((arPtrTrfParInExtrSteps == 0) || (arStrTrfOrderInExtrSteps == 0) || (arNumTrfInExtrSteps == 0) || (NumSteps <= 0)) return 0;

	TVector2d *arTriVertPt1=0, *arTriVertPt2=0;
	int resOK = 0;

	try
	{
		radTGroup* pGroup = new radTGroup();
		if(pGroup == 0) { Send.ErrorMessage("Radia::Error900"); throw 0;}
		radThg hgLoc(pGroup);

		arTriVertPt1 = new TVector2d[numTriVertPt];
		TVector2d *t_arTriVertPt1 = arTriVertPt1, *t_arTriVertPt = arTriVertPt;
		for(int i=0; i<numTriVertPt; i++) *(t_arTriVertPt1++) = *(t_arTriVertPt++);

		arTriVertPt2 = new TVector2d[numTriVertPt];

		//radTPolygon *pBaseFacePgn1 = new radTPolygon(arPoints2d, lenArPoints2d), *pBaseFacePgn2=0;
		//pBaseFacePgn1->CoordZ = 0;

		//radTPolygon genBaseFacePgn1(arPoints2d, lenArPoints2d); //, *pBaseFacePgn2=0;
		//genBaseFacePgn1.CoordZ = 0;
		//radTPolygon genBaseFacePgn2(arPoints2d, lenArPoints2d); //, *pBaseFacePgn2=0;
		//genBaseFacePgn2.CoordZ = 0;

		TVector3d vCenPointRot1(0,0,0);

		radTrans trfInitBaseOrient, trfInitBaseTransl;
		trfInitBaseOrient.SetupRotationToPermutAxes(vCenPointRot1, 'Z', *strOrient);
		//to find appropriate rotation taking into account direction of polygon normal !!!
		//trfInitBaseOrient.SetupRotation(vCenPointRot1, vEZ, vBaseNorm);

		TVector3d vCenPointTransl(0,0,zc);
		vCenPointTransl = trfInitBaseOrient.TrPoint(vCenPointTransl);
		trfInitBaseTransl.SetupTranslation(vCenPointTransl);

		radTrans auxIntermedTrans, auxTotTrans, auxTrans0;
		auxTrans0.SetupIdent();
		if(frame > 0)
		{//Lab frame
			auxTrans0.TrMult(trfInitBaseOrient, 'R');
			auxTrans0.TrMult(trfInitBaseTransl, 'L');
		}

		//TVector3d vNormFace1(0,0,1), vCenPointBaseFace1(pBaseFacePgn1->CentrPoint.x, pBaseFacePgn1->CentrPoint.y, zc);
		TVector3d vNormFace1(0,0,1); //, vCenPointBaseFace1(genBaseFacePgn1.CentrPoint.x, genBaseFacePgn1.CentrPoint.y, zc);
		int NumSteps_mi_1 = NumSteps - 1;

		for(int j=0; j<NumSteps; j++)
		{
			double **arPtrTrfParInCurExtrStep = arPtrTrfParInExtrSteps[j];
			char *arStrTrfOrderInCurExtrStep = arStrTrfOrderInExtrSteps[j];
			int numTrfInCurExtrStep = arNumTrfInExtrSteps[j];
			if((arPtrTrfParInCurExtrStep == 0) || (arStrTrfOrderInCurExtrStep == 0) || (numTrfInCurExtrStep <= 0)) continue;

			//pBaseFacePgn2 = new radTPolygon(*pBaseFacePgn1);
			//radTPolygon genBaseFacePgn2(genBaseFacePgn1);

			TVector2d *t_arTriVertPt1 = arTriVertPt1, *t_arTriVertPt2 = arTriVertPt2;
			for(int i=0; i<numTriVertPt; i++) *(t_arTriVertPt2++) = *(t_arTriVertPt1++);

			auxTotTrans.SetupIdent();

			for(int k=0; k<numTrfInCurExtrStep; k++)
			{
				double *pCurTrfPar = arPtrTrfParInCurExtrStep[k];
				if(pCurTrfPar == 0) continue;

				char cTypeCurTrf = arStrTrfOrderInCurExtrStep[k];
				if((cTypeCurTrf == 'r') || (cTypeCurTrf == 'R'))
				{
					TVector3d vRotCen(*pCurTrfPar, *(pCurTrfPar+1), *(pCurTrfPar+2));
					TVector3d vRotAxis(*(pCurTrfPar+3), *(pCurTrfPar+4), *(pCurTrfPar+5));
					auxIntermedTrans.SetupRotation(vRotCen, vRotAxis, *(pCurTrfPar+6));
					auxTotTrans.TrMult(auxIntermedTrans, 'L');
				}
				else if((cTypeCurTrf == 't') || (cTypeCurTrf == 'T'))
				{
					TVector3d vTransl(*pCurTrfPar, *(pCurTrfPar+1), *(pCurTrfPar+2));
					auxIntermedTrans.SetupTranslation(vTransl);
					auxTotTrans.TrMult(auxIntermedTrans, 'L');
				}
				else if((cTypeCurTrf == 'h') || (cTypeCurTrf == 'H'))
				{
					//any homothety is applied to a face before space transformations
					//to check if this should be improved
					//pBaseFacePgn2->ApplyHomothety(*pCurTrfPar, *(pCurTrfPar+1), *(pCurTrfPar+2)); //ApplyHomothety(double kxH, double kyH, double phi)

					double kxH = *pCurTrfPar, kyH = *(pCurTrfPar+1), phi = *(pCurTrfPar+2);
					
					//genBaseFacePgn2.ApplyHomothety(kxH, kyH, phi); //ApplyHomothety(double kxH, double kyH, double phi)

					//TVector2d &CenPointBase2 = genBaseFacePgn2.CentrPoint;
					TVector2d vSum(0,0);
					for(int jjj=0; jjj<lenArPoints2d; jjj++) vSum += arPoints2d[jjj];
					TVector2d CenPointBase2 = (1./lenArPoints2d)*vSum;

					if((kxH != 0) && (kyH != 0))
					{
						double cosPhi = cos(phi), sinPhi = sin(phi);
						TVector2d vHomLocX(cosPhi, sinPhi), vHomLocY(-sinPhi, cosPhi);
						for(int ip=0; ip<numTriVertPt; ip++)
						{
							TVector2d vRelVertex = arTriVertPt2[ip] - CenPointBase2;
							arTriVertPt2[ip] = CenPointBase2 + ((kxH*(vHomLocX*vRelVertex))*vHomLocX) + ((kyH*(vHomLocY*vRelVertex))*vHomLocY);
						}
					}
				}
			}
			radTrans *pBaseFaceTrf1 = new radTrans(auxTrans0);

			if(frame == 0) auxTrans0.TrMult(auxTotTrans, 'R'); //local frame
			else auxTrans0.TrMult(auxTotTrans, 'L'); //laboratory frame

			radTrans *pBaseFaceTrf2 = new radTrans(auxTrans0); //check for memory leak !!!
			radTHandle<radTrans> hTrf1(pBaseFaceTrf1), hTrf2(pBaseFaceTrf2);

			double *pMagComp = 0;
			if(arMagnCompInExtrSteps != 0) pMagComp = arMagnCompInExtrSteps + j*3;

			TVector2d arTreePt1[3], arTreePt2[3];
			int *t_arTriVertInd = arTriVertInd;
			for(int itr=0; itr<numTri; itr++)
			{
				arTreePt1[0] = arTriVertPt1[*t_arTriVertInd]; arTreePt2[0] = arTriVertPt2[*(t_arTriVertInd++)];
				arTreePt1[1] = arTriVertPt1[*t_arTriVertInd]; arTreePt2[1] = arTriVertPt2[*(t_arTriVertInd++)];
				arTreePt1[2] = arTriVertPt1[*t_arTriVertInd]; arTreePt2[2] = arTriVertPt2[*(t_arTriVertInd++)]; 
				radTPolygon *pBaseFaceTri1 = new radTPolygon(arTreePt1, 3);
				//pBaseFaceTri1->CoordZ = genBaseFacePgn1.CoordZ;
				radTPolygon *pBaseFaceTri2 = new radTPolygon(arTreePt2, 3);
				//pBaseFaceTri2->CoordZ = genBaseFacePgn2.CoordZ;

				//this doesn't take into account direction of normals and eventual intersection of the faces
				radTHandle<radTPolygon> hPgn1(pBaseFaceTri1), hPgn2(pBaseFaceTri2);
				radTHandlePgnAndTrans hPgnTrf1(hPgn1, hTrf1), hPgnTrf2(hPgn2, hTrf2);
				//create polyhedron here from pBaseFacePgn1, pBaseFaceTrf1, pBaseFacePgn2, pBaseFaceTrf2
				radTPolyhedron *pCurStepPlhdr = new radTPolyhedron(hPgnTrf1, hPgnTrf2, avgCur, pMagComp);
				if(pCurStepPlhdr == 0) { Send.ErrorMessage("Radia::Error900"); throw 0;}
				if(pCurStepPlhdr->SomethingIsWrong) { delete pCurStepPlhdr; throw 0;}

				radThg hgCurStepPlhdr(pCurStepPlhdr);
				pGroup->AddElement(AddElementToContainer(hgCurStepPlhdr), hgCurStepPlhdr);
			}

			//if(j < NumSteps_mi_1)
			//{
			//	//pBaseFacePgn1 = new radTPolygon(*pBaseFacePgn2);
			//	genBaseFacePgn1 = genBaseFacePgn2;
			////	//- treat frame (0- local, otherwise "laboratory");
			////	//if(frame == 0) pBaseFaceTrf1 = new radTrans(*pBaseFaceTrf1); //leave base trf always the same if Frame->Loc ??
			////	//else pBaseFaceTrf1 = new radTrans(*pBaseFaceTrf2); //change base trf if Frame->Lab

			////	if(frame == 0) pBaseFaceTrf1 = new radTrans(*pBaseFaceTrf2); //change base trf if Frame->Loc
			////	else pBaseFaceTrf1 = new radTrans(*pBaseFaceTrf1); //leave base trf always the same if Frame->Lab
			//}
		}
		if(frame == 0)
		{
			auxTrans0.SetupIdent();
			auxTrans0.TrMult(trfInitBaseOrient, 'R');
			auxTrans0.TrMult(trfInitBaseTransl, 'L');
			
			double relTol = 1.E-12; //to ensure correct operation on 32-bit OS
			if(!auxTrans0.IsIdent(relTol))
			{
				radTrans* TransPtr = new radTrans(auxTrans0);
				if(TransPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
				radThg hgTrf(TransPtr);
				((radTg3d*)(hgLoc.rep))->AddTransform(1, hgTrf);
			}
		}

		hgOut = hgLoc;
		resOK = 1;
	}
	catch(...)
	{
		Initialize(); //return 0;
		resOK = 0;
	}
	if(arTriVertPt1 != 0) delete[] arTriVertPt1;
	if(arTriVertPt2 != 0) delete[] arTriVertPt2;
	return resOK;
}

//-------------------------------------------------------------------------

int radTApplication::SetUpPolyhedronsFromLayerPolygons(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons, TVector3d& Magn, radThg& hgOut)
{
	radTVectVect3d VertexPointsVect;
	radTVectIntPtrAndInt FacesVect;
	radTPtrsToPgnAndVect2d PtrsToPgnAndVect2d[2];

	double RelZeroToler = 1.E-09;
	RelZeroToler = 10.*((RelZeroToler>radCR.RelRand)? RelZeroToler : radCR.RelRand);
	double TypSize = fabs(CoordsZ[AmOfLayerPolygons - 1] - *CoordsZ);
	double RelAbsTol[] = { RelZeroToler, RelZeroToler*TypSize };

	try
	{
		PtrsToPgnAndVect2d->AmOfPoints = *PtsNumbersInLayerPgns;
		if(*PtsNumbersInLayerPgns < 3)
		{
			PtrsToPgnAndVect2d->pPgn = 0;
			PtrsToPgnAndVect2d->pVect2d = LayerPolygons[0];
		}
		else
		{
			PtrsToPgnAndVect2d->pPgn = new radTPolygon(LayerPolygons[0], PtsNumbersInLayerPgns[0]);
			if(PtrsToPgnAndVect2d->pPgn==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			if(PtrsToPgnAndVect2d->pPgn->SomethingIsWrong) 
			{
				delete PtrsToPgnAndVect2d->pPgn; return 0;
			}
			if(!PtrsToPgnAndVect2d->pPgn->IsConvex) { Send.ErrorMessage("Radia::Error111"); return 0;}
			PtrsToPgnAndVect2d->pVect2d = 0;
		}

		radTGroup* pGroup = new radTGroup();
		if(pGroup == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(pGroup);

		double z1z2[2];
		*z1z2 = *CoordsZ;

		radTPtrsToPgnAndVect2d& PtrsToPgnAndVect2d_1 = PtrsToPgnAndVect2d[1];
		int AmOfLayerPolygons_m_1 = AmOfLayerPolygons - 1;
		for(int i=1; i<AmOfLayerPolygons; i++)
		{
			PtrsToPgnAndVect2d_1.AmOfPoints = PtsNumbersInLayerPgns[i];
			if(PtsNumbersInLayerPgns[i] < 3)
			{
				PtrsToPgnAndVect2d_1.pPgn = 0;
				PtrsToPgnAndVect2d_1.pVect2d = LayerPolygons[i];
			}
			else
			{
				PtrsToPgnAndVect2d_1.pPgn = new radTPolygon(LayerPolygons[i], PtsNumbersInLayerPgns[i]);
				if(PtrsToPgnAndVect2d_1.pPgn==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
				if(PtrsToPgnAndVect2d_1.pPgn->SomethingIsWrong) 
				{
					delete PtrsToPgnAndVect2d_1.pPgn; return 0;
				}
				if(!PtrsToPgnAndVect2d_1.pPgn->IsConvex) { Send.ErrorMessage("Radia::Error111"); return 0;}
				PtrsToPgnAndVect2d_1.pVect2d = 0;
			}

			z1z2[1] = CoordsZ[i];
			char StageChar = (i == 1)? ((i == AmOfLayerPolygons_m_1)? 'T' : 'S') : ((i == AmOfLayerPolygons_m_1)? 'F' : 'O');

			if(!SetUpOnePolyhedronSegment(PtrsToPgnAndVect2d, z1z2, StageChar, RelAbsTol, Magn, &VertexPointsVect, &FacesVect, hgLoc)) return 0;

			if(PtrsToPgnAndVect2d->pPgn != 0) delete PtrsToPgnAndVect2d->pPgn;
			*PtrsToPgnAndVect2d = PtrsToPgnAndVect2d_1;
			*z1z2 = z1z2[1];
		}
		if(PtrsToPgnAndVect2d->pPgn != 0) delete PtrsToPgnAndVect2d->pPgn;

		if(!CheckIfGroupIsNeeded(hgLoc)) return 0;
		hgOut = hgLoc;

		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetUpOnePolyhedronSegment(radTPtrsToPgnAndVect2d* PtrsToPgnAndVect2d, double* z1z2, char StageChar, double* RelAbsTol, TVector3d& Magn, radTVectVect3d* pVertexPointsVect, radTVectIntPtrAndInt* pFacesVect, radThg& hgGroup)
{
	radTVect2dVect LocVect2dVect1, LocVect2dVect2;
	radTVect2dVect *pFirstVect2dVect, *pSecondVect2dVect;
	radTVect2dVect *pBaseVect2dVect, *pNextVect2dVect;
	radTPolygon* pBasePgn;

	radTPtrsToPgnAndVect2d* PtrsToPgnAndVect2d_p_1 = PtrsToPgnAndVect2d + 1;
	if(PtrsToPgnAndVect2d->AmOfPoints < 3)
	{
		TVector2d* tVect2d = PtrsToPgnAndVect2d->pVect2d;
		for(int kk=0; kk<PtrsToPgnAndVect2d->AmOfPoints; kk++) LocVect2dVect1.push_back(*(tVect2d++));
		pFirstVect2dVect = &LocVect2dVect1;
	}
	else pFirstVect2dVect = &(PtrsToPgnAndVect2d->pPgn->EdgePointsVector);
	if(PtrsToPgnAndVect2d_p_1->AmOfPoints < 3)
	{
		TVector2d* tVect2d = PtrsToPgnAndVect2d_p_1->pVect2d;
		for(int kk=0; kk<PtrsToPgnAndVect2d_p_1->AmOfPoints; kk++) LocVect2dVect2.push_back(*(tVect2d++));
		pSecondVect2dVect = &LocVect2dVect2;
	}
	else pSecondVect2dVect = &(PtrsToPgnAndVect2d_p_1->pPgn->EdgePointsVector);

	double aCoordsZ;
	if((StageChar == 'S') || (StageChar == 'T'))
	{
		aCoordsZ = *z1z2;
		char MakeBaseFace = (PtrsToPgnAndVect2d->AmOfPoints > 2);
		int* NewFace = 0;
		if(MakeBaseFace)
		{
			NewFace = new int[PtrsToPgnAndVect2d->AmOfPoints];
			if(NewFace == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

			radTIntPtrAndInt IntPtrAndInt(NewFace, PtrsToPgnAndVect2d->AmOfPoints);
			pFacesVect->push_back(IntPtrAndInt);
		}
		for(int kk=0; kk<PtrsToPgnAndVect2d->AmOfPoints; kk++)
		{
			TVector2d& p2d = (*pFirstVect2dVect)[kk];
			TVector3d aVertexPo(p2d.x, p2d.y, aCoordsZ);
			pVertexPointsVect->push_back(aVertexPo);

			if(MakeBaseFace) NewFace[kk] = kk + 1;
		}
	}

	int PrevVertexPoVectSize = (int)(pVertexPointsVect->size());
	int PrevFacesVectSize = (int)(pFacesVect->size());

	aCoordsZ = z1z2[1];
	for(int kkk=0; kkk<PtrsToPgnAndVect2d_p_1->AmOfPoints; kkk++)
	{
		TVector2d& p2d = (*pSecondVect2dVect)[kkk];
		TVector3d aVertexPo(p2d.x, p2d.y, aCoordsZ);
		pVertexPointsVect->push_back(aVertexPo);
	}

	int BaseVertexPoVectOffset, NextVertexPoVectOffset;
	if(PtrsToPgnAndVect2d->AmOfPoints >= PtrsToPgnAndVect2d_p_1->AmOfPoints)
	{
		pBaseVect2dVect = pFirstVect2dVect; pNextVect2dVect = pSecondVect2dVect;
		pBasePgn = PtrsToPgnAndVect2d->pPgn;
		BaseVertexPoVectOffset = PrevVertexPoVectSize - (int)(pFirstVect2dVect->size()) + 1;
		NextVertexPoVectOffset = PrevVertexPoVectSize + 1;
	}
	else
	{
		pBaseVect2dVect = pSecondVect2dVect; pNextVect2dVect = pFirstVect2dVect;
		pBasePgn = PtrsToPgnAndVect2d_p_1->pPgn;
		BaseVertexPoVectOffset = PrevVertexPoVectSize + 1;
		NextVertexPoVectOffset = PrevVertexPoVectSize - (int)(pFirstVect2dVect->size()) + 1;
	}
	
	char CuttingIsNeeded = 0;
	char CurrentPartCanOnlyBeTetrahedron = 0;
	if((PtrsToPgnAndVect2d->AmOfPoints < 3) && (PtrsToPgnAndVect2d_p_1->AmOfPoints < 3)) 
	{
		if((PtrsToPgnAndVect2d->AmOfPoints == 2) && (PtrsToPgnAndVect2d_p_1->AmOfPoints == 2))
		{
			CurrentPartCanOnlyBeTetrahedron = 1;
			CuttingIsNeeded = 1;
		}
		else { Send.ErrorMessage("Radia::Error063"); return 0;}
	}

	char OnlyOnePointInUpperLayer = (PtrsToPgnAndVect2d_p_1->AmOfPoints == 1);
	char OnlyOnePointInLowerLayer = (PtrsToPgnAndVect2d->AmOfPoints == 1);

	if(!CurrentPartCanOnlyBeTetrahedron)
	{
		radTVertexPointLiaison* NextLayerPointLinks = new radTVertexPointLiaison[pNextVect2dVect->size()];
		if(NextLayerPointLinks == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		int AmOfBasePoints = (int)(pBaseVect2dVect->size());
		int AmOfBasePoints_m_1 = AmOfBasePoints - 1;
		int AmOfNextPoints = (int)(pNextVect2dVect->size());
		int AmOfNextPoints_m_1 = AmOfNextPoints - 1;

		for(int k=0; k<AmOfBasePoints; k++)
		{
			int NextInd = (k != AmOfBasePoints_m_1)? (k+1) : 0;

			TVector2d SegmVect = (*pBaseVect2dVect)[NextInd] - (*pBaseVect2dVect)[k];

			char AmOfLwstPo;
			int IndOfPointFromNextLayer;
			if(!FindLowestPoint(pNextVect2dVect, SegmVect, RelAbsTol, IndOfPointFromNextLayer, AmOfLwstPo)) return 0;

			radTVertexPointLiaison* pNextLayerPointLinks = NextLayerPointLinks + IndOfPointFromNextLayer;
			if(AmOfLwstPo == 2) pNextLayerPointLinks->AdjSegmentUsed = 1;
			pNextLayerPointLinks->FirstIndVect.push_back(k);
			pNextLayerPointLinks->SecondIndVect.push_back(NextInd);

			int AmOfPoInFace = AmOfLwstPo + 2;
			int* NewFace = new int[AmOfPoInFace];
			if(NewFace == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			*NewFace = NextVertexPoVectOffset + IndOfPointFromNextLayer;
			NewFace[1] = BaseVertexPoVectOffset + k;
			NewFace[2] = BaseVertexPoVectOffset + NextInd;
			if(AmOfLwstPo == 2) 
			{
				int NextIndOfPointFromNextLayer = (IndOfPointFromNextLayer != AmOfNextPoints_m_1)? (IndOfPointFromNextLayer+1) : 0;
				NewFace[3] = NextVertexPoVectOffset + NextIndOfPointFromNextLayer;

				radTVertexPointLiaison* pNextNextLayerPointLinks = NextLayerPointLinks + NextIndOfPointFromNextLayer;
				pNextNextLayerPointLinks->SecondIndVect.push_back(NextInd);
			}
			radTIntPtrAndInt IntPtrAndInt(NewFace, AmOfPoInFace);
			pFacesVect->push_back(IntPtrAndInt);
		}

		if(AmOfNextPoints > 1)
		{
			for(int i=0; i<AmOfNextPoints; i++)
			{
				int NextInd = (i != AmOfNextPoints_m_1)? (i+1) : 0;

				radTVertexPointLiaison* pThisPoLinks = NextLayerPointLinks + i;
				radTVertexPointLiaison* pNextPoLinks = NextLayerPointLinks + NextInd;

				if(!pThisPoLinks->AdjSegmentUsed)
				{
					if(!(pThisPoLinks->SecondIndVect.empty()))
					{
						for(int j=0; j<(int)(pThisPoLinks->SecondIndVect.size()); j++)
						{
							int IndOfPoFromBase = (pThisPoLinks->SecondIndVect)[j];
							char LetsMakeFace = 0;
							for(int k=0; k<(int)(pNextPoLinks->FirstIndVect.size()); k++)
							{
								if(IndOfPoFromBase == (pNextPoLinks->FirstIndVect)[k])
								{
									LetsMakeFace = 1; break;
								}
							}
							if(LetsMakeFace)
							{
								int* NewFace = new int[3];
								if(NewFace == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
								*NewFace = NextVertexPoVectOffset + i;
								NewFace[1] = BaseVertexPoVectOffset + IndOfPoFromBase;
								NewFace[2] = NextVertexPoVectOffset + NextInd;
								radTIntPtrAndInt IntPtrAndInt(NewFace, 3);
								pFacesVect->push_back(IntPtrAndInt);
							}
						}
					}
					else if(pThisPoLinks->FirstIndVect.empty())
					{
						int PrevInd = (i != 0)? (i-1) : AmOfNextPoints_m_1;

						TVector2d SegmVect = (*pNextVect2dVect)[NextInd] - (*pNextVect2dVect)[i];

						char AmOfLwstPo;
						int IndOfPointFromBaseLayer;
						if(!FindLowestPoint(pBaseVect2dVect, SegmVect, RelAbsTol, IndOfPointFromBaseLayer, AmOfLwstPo)) return 0;

						int* NewFace = new int[3];
						if(NewFace == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
						*NewFace = NextVertexPoVectOffset + i;
						NewFace[1] = NextVertexPoVectOffset + PrevInd;
						NewFace[2] = BaseVertexPoVectOffset + IndOfPointFromBaseLayer;
						radTIntPtrAndInt IntPtrAndInt1(NewFace, 3);
						pFacesVect->push_back(IntPtrAndInt1);

						NewFace = new int[3];
						if(NewFace == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
						*NewFace = NextVertexPoVectOffset + i;
						NewFace[1] = BaseVertexPoVectOffset + IndOfPointFromBaseLayer;
						NewFace[2] = NextVertexPoVectOffset + NextInd;
						radTIntPtrAndInt IntPtrAndInt2(NewFace, 3);
						pFacesVect->push_back(IntPtrAndInt2);
					}
				}
			}
		}

		if(!(OnlyOnePointInUpperLayer || OnlyOnePointInLowerLayer))
		{
			if((StageChar != 'S') && (StageChar != 'T'))
			{// Check for convexity. If not convex: CuttingIsNeeded = 1;
				if(PtrsToPgnAndVect2d->AmOfPoints < 3) { CuttingIsNeeded = 1;}
				else
				{
					int AmOfPointsOnFirstLayer = (int)(pFirstVect2dVect->size());
					int AmOfPointsOnFirstLayer_m_1 = AmOfPointsOnFirstLayer - 1;
					int FirstLayerVertexPoIndOffset = PrevVertexPoVectSize - (int)(pFirstVect2dVect->size()) + 1;

					for(int k=0; k<AmOfPointsOnFirstLayer; k++)
					{
						int NextInd = (k != AmOfPointsOnFirstLayer_m_1)? (k+1) : 0;
						int ThisVertPoInd = FirstLayerVertexPoIndOffset + k;
						int NextVertPoInd = FirstLayerVertexPoIndOffset + NextInd;

						int OneFaceInd, AnotherFaceInd;
						int aPoOnOneFaceInd, aPoOnAnotherFaceInd;
						if(!FindTwoAdjacentFaces(ThisVertPoInd, NextVertPoInd, pFacesVect, OneFaceInd, aPoOnOneFaceInd, AnotherFaceInd, aPoOnAnotherFaceInd)) { Send.ErrorMessage("Radia::Error110"); return 0;}

						TVector3d& P0 = (*pVertexPointsVect)[ThisVertPoInd - 1];
						TVector3d V = (*pVertexPointsVect)[NextVertPoInd - 1] - P0;
						TVector3d& aP1 = (*pVertexPointsVect)[((*pFacesVect)[OneFaceInd].pInt)[aPoOnOneFaceInd] - 1];
						TVector3d& aP2 = (*pVertexPointsVect)[((*pFacesVect)[AnotherFaceInd].pInt)[aPoOnAnotherFaceInd] - 1];
						TVector3d P1, P2;
						if(aP1.z < aP2.z) { P1 = aP1; P2 = aP2;}
						else { P2 = aP1; P1 = aP2;}
						TVector3d P1mP0 = P1 - P0, P2mP0 = P2 - P0;
						TVector3d VP1 = P1mP0 - (P1mP0*V)*V;
						TVector3d VP2 = P2mP0 - (P2mP0*V)*V;
						TVector3d VP2_VectBy_VP1 = VP2^VP1;

						double ScalProd = VP2_VectBy_VP1*V;
						if(ScalProd > -RelAbsTol[1]) { CuttingIsNeeded = 1; break;}
					}
				}
			}
		}
		else 
		{
			CuttingIsNeeded = 0;
		}

		if(NextLayerPointLinks != 0) delete[] NextLayerPointLinks;
	}

	char FinalPassAfterCut = 0;

	char FinalizeAll = ((StageChar == 'F') || (StageChar == 'T'));
	if(FinalizeAll || CuttingIsNeeded || OnlyOnePointInUpperLayer)
	{
CreatingPolyhedronPiece:

		int AmOfVertexPoints, AmOfPoToDelete, AmOfPoInLastFace, StartPoIndForNewFace;
		if(CuttingIsNeeded)
		{
			AmOfVertexPoints = PrevVertexPoVectSize;
			AmOfPoToDelete = PrevVertexPoVectSize - PtrsToPgnAndVect2d->AmOfPoints;
			AmOfPoInLastFace = PtrsToPgnAndVect2d->AmOfPoints;
			StartPoIndForNewFace = PrevVertexPoVectSize - PtrsToPgnAndVect2d->AmOfPoints + 1;
		}
		else
		{
			AmOfVertexPoints = (int)(pVertexPointsVect->size());
			AmOfPoToDelete = (!FinalPassAfterCut)? PrevVertexPoVectSize : AmOfVertexPoints;
			AmOfPoInLastFace = PtrsToPgnAndVect2d_p_1->AmOfPoints;
			StartPoIndForNewFace = (!FinalPassAfterCut)? (PrevVertexPoVectSize + 1) : ((int)(pVertexPointsVect->size()) - PtrsToPgnAndVect2d_p_1->AmOfPoints + 1);
		}

		TVector3d* ArrayOfVertexPoints = new TVector3d[AmOfVertexPoints];
		if(ArrayOfVertexPoints == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		TVector3d* tArrayOfVertexPoints = ArrayOfVertexPoints;
		for(int p=0; p<AmOfVertexPoints; p++)
		{
			*(tArrayOfVertexPoints++) = (*pVertexPointsVect)[p];
		}
		radTVectVect3d::iterator VertPoVectStartIter = pVertexPointsVect->begin();
		pVertexPointsVect->erase(VertPoVectStartIter, VertPoVectStartIter + AmOfPoToDelete);

		if(AmOfPoInLastFace > 2)
		{
			int* NewFace = new int[AmOfPoInLastFace];
			if(NewFace == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			int* tNewFace = NewFace;
			for(int s=0; s<AmOfPoInLastFace; s++) *(tNewFace++) = StartPoIndForNewFace + s;
			radTIntPtrAndInt IntPtrAndInt(NewFace, AmOfPoInLastFace);

			if(CuttingIsNeeded) 
			{
				radTVectIntPtrAndInt::iterator NewFaceInsertPos = pFacesVect->begin() + PrevFacesVectSize;
				pFacesVect->insert(NewFaceInsertPos, IntPtrAndInt);
			}
			else
			{
				pFacesVect->push_back(IntPtrAndInt);
			}
		}

		int AmOfFaces = CuttingIsNeeded? (PrevFacesVectSize + ((AmOfPoInLastFace > 2)? 1 : 0)) : (int)(pFacesVect->size());
		int* ArrayOfFacesLengths = new int[AmOfFaces];
		if(ArrayOfFacesLengths == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		int** ArrayOfFaces = new int*[AmOfFaces];
		if(ArrayOfFaces == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		int* tArrayOfFacesLengths = ArrayOfFacesLengths;
		int** tArrayOfFaces = ArrayOfFaces;
		for(int m=0; m<AmOfFaces; m++)
		{
			*(tArrayOfFaces++) = (*pFacesVect)[m].pInt;
			*(tArrayOfFacesLengths++) = (*pFacesVect)[m].AnInt;
		}
		radTVectIntPtrAndInt::iterator FacesVectStartIter = pFacesVect->begin();
		int AmOfFacesToDelete = (AmOfPoInLastFace > 2)? (AmOfFaces - 1) : AmOfFaces;
		if(FinalizeAll && (!CuttingIsNeeded)) AmOfFacesToDelete = AmOfFaces;
		pFacesVect->erase(FacesVectStartIter, FacesVectStartIter + AmOfFacesToDelete);

		if(AmOfFaces > 3)
		{
			radTPolyhedron* pNewPlhdr = new radTPolyhedron(ArrayOfVertexPoints, AmOfVertexPoints, ArrayOfFaces, ArrayOfFacesLengths, AmOfFaces, Magn);
			
			if(pNewPlhdr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			if(pNewPlhdr->SomethingIsWrong) { delete pNewPlhdr; return 0;}
			radThg hgNew(pNewPlhdr);
			((radTGroup*)((radTg3d*)(hgGroup.rep)))->AddElement(AddElementToContainer(hgNew), hgNew);
		}

		if(CuttingIsNeeded) 
		{
			if(CurrentPartCanOnlyBeTetrahedron)
			{
				if(!SetUpTetrahedronBasedOnTwoLinSegm(pFirstVect2dVect, z1z2[0], pSecondVect2dVect, z1z2[1], RelAbsTol, Magn, hgGroup)) return 0;

				radTVectVect3d::iterator Iter = pVertexPointsVect->end();
				--Iter; --Iter;
				radTVectVect3d::iterator FinEraseIter = Iter;
				--Iter; --Iter;
				pVertexPointsVect->erase(Iter, FinEraseIter);

				CuttingIsNeeded = 0;
			}
			else
			{
				if(AmOfPoToDelete > 0)
				{
					char DistinguishFirstFace = (AmOfPoInLastFace > 2);
					if(!ShiftVertexPointNumbersInFaces(pFacesVect, AmOfPoToDelete, DistinguishFirstFace)) return 0;
				}
			}
		}

		if(ArrayOfVertexPoints != 0) delete[] ArrayOfVertexPoints;
		if(ArrayOfFacesLengths != 0) delete[] ArrayOfFacesLengths;
		if(ArrayOfFaces != 0) 
		{
			for(int j=0; j<AmOfFacesToDelete; j++) delete[] (ArrayOfFaces[j]);
			delete[] ArrayOfFaces;
		}
		if(FinalizeAll && CuttingIsNeeded)
		{
			CuttingIsNeeded = 0; FinalPassAfterCut = 1;
			goto CreatingPolyhedronPiece;
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::ShiftVertexPointNumbersInFaces(radTVectIntPtrAndInt* pFacesVect, int InAmOfPoToDelete, char DistinguishFirstFace)
{
	char FirstFace = 1;
	if(!DistinguishFirstFace) FirstFace = 0;
	for(radTVectIntPtrAndInt::iterator FaceIt = pFacesVect->begin(); FaceIt != pFacesVect->end(); ++FaceIt)
	{
		int AmOfPoToDelete = FirstFace? InAmOfPoToDelete - 1 : InAmOfPoToDelete;
		int* tPoInd = FaceIt->pInt;
		for(int k=0; k<FaceIt->AnInt; k++)
		{
			*tPoInd -= AmOfPoToDelete;
			if(*tPoInd < 1) { Send.ErrorMessage("Radia::Error110"); return 0;}
			tPoInd++;
		}
		FirstFace = 0;
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::FindTwoAdjacentFaces(int OneVertPoInd, int AnotherVertPoInd, radTVectIntPtrAndInt* pFacesVect, int& OneFaceInd, int& IndOfPoOnOneFace, int& AnotherFaceInd, int& IndOfPoOnAnotherFace)
{
	int LocOneFaceInd, LocAnotherFaceInd;
	int LocIndOfPoOnOneFace, LocIndOfPoOnAnotherFace;
	char OneFaceFound = 0, AnotherFaceFound = 0;
	char SeparatePoOnOneFaceFound = 0, SeparatePoOnAnotherFaceFound = 0;
	int FaceNo = 0;
	for(radTVectIntPtrAndInt::iterator FaceIt = pFacesVect->begin(); FaceIt != pFacesVect->end(); ++FaceIt)
	{
		char OnePoFound = 0, AnotherPoFound = 0;
		int BufPoInd;
		char BufPoIndIsSet = 0;
		int* tPoInd = FaceIt->pInt;
		for(int k=0; k<FaceIt->AnInt; k++)
		{
			if(*tPoInd == OneVertPoInd) OnePoFound = 1;
			else if(*tPoInd == AnotherVertPoInd) AnotherPoFound = 1;
			else { BufPoInd = k; BufPoIndIsSet = 1;}

			if(OnePoFound && AnotherPoFound && BufPoIndIsSet) break;
			tPoInd++;
		}
		if(OnePoFound && AnotherPoFound && BufPoIndIsSet)
		{
			if(!OneFaceFound) 
			{ 
				LocOneFaceInd = FaceNo; OneFaceFound = 1;
				LocIndOfPoOnOneFace = BufPoInd;
			}
			else if(!AnotherFaceFound) 
			{
				LocAnotherFaceInd = FaceNo; AnotherFaceFound = 1;
				LocIndOfPoOnAnotherFace = BufPoInd; break;
			}
		}
		FaceNo++;
	}
	if(!(OneFaceFound && AnotherFaceFound)) return 0;

	OneFaceInd = LocOneFaceInd; AnotherFaceInd = LocAnotherFaceInd;
	IndOfPoOnOneFace = LocIndOfPoOnOneFace; IndOfPoOnAnotherFace = LocIndOfPoOnAnotherFace;
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::FindLowestPoint(radTVect2dVect* pVectP2d, TVector2d& SegmVect, double* RelAbsTol, int& LowestPointInd, char& AmOfPo)
{
	TVector2d V(-SegmVect.y, SegmVect.x);

	double InvLenV = 1./sqrt(V.x*V.x + V.y*V.y);
	TVector2d LocV = InvLenV*V;
	double AbsTol = RelAbsTol[1];

	int pNo = 0;
	int LocLowestPointInd, AnotherLocLowestPointInd = -1, ThirdLocLowestPointInd = -1;
	double MinScalProd = 1.E+23;

	for(radTVect2dVect::iterator It = pVectP2d->begin(); It != pVectP2d->end(); ++It)
	{
		double ScalProd = (*It)*LocV;
		if(fabs(ScalProd - MinScalProd) < AbsTol)
		{
			if(AnotherLocLowestPointInd > -1) ThirdLocLowestPointInd = pNo;
			else AnotherLocLowestPointInd = pNo;
		}
		else if(ScalProd < MinScalProd - AbsTol)
		{
			LocLowestPointInd = pNo;
			MinScalProd = ScalProd;
			if(AnotherLocLowestPointInd > -1) AnotherLocLowestPointInd = -1;
			if(ThirdLocLowestPointInd > -1) ThirdLocLowestPointInd = -1;
		}
		pNo++;
	}

	if(ThirdLocLowestPointInd > -1) { Send.ErrorMessage("Radia::Error112"); return 0;}

	AmOfPo = 1;
	LowestPointInd = LocLowestPointInd;

	if(AnotherLocLowestPointInd > -1)
	{
		AmOfPo = 2;
		
		TVector2d TestV = (*pVectP2d)[AnotherLocLowestPointInd] - (*pVectP2d)[LocLowestPointInd];
		double TestScalProd = TestV*SegmVect;
		if(TestScalProd < 0) LowestPointInd = AnotherLocLowestPointInd;
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::SetUpTetrahedronBasedOnTwoLinSegm(radTVect2dVect* pFirstVect2dVect, double z1, radTVect2dVect* pSecondVect2dVect, double z2, double* RelAbsTol, TVector3d& Magn, radThg& hgGroup)
{
	TVector2d &p11 = (*pFirstVect2dVect)[0], &p12 = (*pFirstVect2dVect)[1];
	TVector2d &p21 = (*pSecondVect2dVect)[0], &p22 = (*pSecondVect2dVect)[1];

	double AbsTol = RelAbsTol[1];
	TVector2d v1 = p12 - p11, v2 = p22 - p21;
	double InvLen_v1 = 1./sqrt(v1.x*v1.x + v1.y*v1.y);
	v1 = InvLen_v1*v1;
	if(fabs(v1.x*v2.y - v2.x*v1.y) < AbsTol) { Send.ErrorMessage("Radia::Error063"); return 0;}

	TVector3d VertexPoints[4];
	TVector3d* tVertexPoints = VertexPoints;
	tVertexPoints->x = p11.x; tVertexPoints->y = p11.y; tVertexPoints->z = z1; tVertexPoints++;
	tVertexPoints->x = p12.x; tVertexPoints->y = p12.y; tVertexPoints->z = z1; tVertexPoints++;
	tVertexPoints->x = p21.x; tVertexPoints->y = p21.y; tVertexPoints->z = z2; tVertexPoints++;
	tVertexPoints->x = p22.x; tVertexPoints->y = p22.y; tVertexPoints->z = z2;
	int* Faces[4];
	int Face0[] = { 1,2,3 }; Faces[0] = Face0;
	int Face1[] = { 1,3,4 }; Faces[1] = Face1;
	int Face2[] = { 3,2,4 }; Faces[2] = Face2;
	int Face3[] = { 1,4,2 }; Faces[3] = Face3;
	int Lengths[] = { 3,3,3,3 };

	radTPolyhedron* pTetr = new radTPolyhedron(VertexPoints, 4, Faces, Lengths, 4, Magn);
	if(pTetr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
	if(pTetr->SomethingIsWrong) { delete pTetr; return 0;}
	radThg hgNew(pTetr);
	((radTGroup*)((radTg3d*)(hgGroup.rep)))->AddElement(AddElementToContainer(hgNew), hgNew);
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::CheckLayerPolygonStructures(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons)
{
	if(AmOfLayerPolygons < 2) { Send.ErrorMessage("Radia::Error063"); return 0;}

	double CurrentZ = *CoordsZ;
	char GoesDown = 0, GoesUp = 0, NoLayersWithPtsNumMoreThan1 = 1;

	for(int i = 0; i < AmOfLayerPolygons; i++)
	{
		if(i > 0)
		{
			if(CoordsZ[i] > CurrentZ)
			{
				if(i == 1) GoesUp = 1;
				else if(GoesDown || (!GoesUp)) { Send.ErrorMessage("Radia::Error063"); return 0;}
			}
			else if(CoordsZ[i] < CurrentZ)
			{
				if(i == 1) GoesDown = 1;
				else if(GoesUp || (!GoesDown)) { Send.ErrorMessage("Radia::Error063"); return 0;}
			}
			else { Send.ErrorMessage("Radia::Error063"); return 0;}
			CurrentZ = CoordsZ[i];
		}
		if(PtsNumbersInLayerPgns[i] > 1) 
		{
			NoLayersWithPtsNumMoreThan1 = 0;
		}
	}
	if(NoLayersWithPtsNumMoreThan1) { Send.ErrorMessage("Radia::Error063"); return 0;}
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::CheckIfGroupIsNeeded(radThg& In_hg)
{
	radThg hgOld = In_hg;
	radTGroup* pGroup = Cast.GroupCast((radTg3d*)(In_hg.rep));
	if(pGroup == 0) return 0;

	int Size = (int)(pGroup->GroupMapOfHandlers.size());
	if(Size == 0) return 0;
	else if(Size == 1)
	{
		In_hg = (*(pGroup->GroupMapOfHandlers.begin())).second;
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::SetUpPolyhedronsFromLayerRectangles(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect, TVector3d& MagnVect, radThg& hg)
{
	TVector2d** LayerPolygons = 0;
	int* PtsNumbersInLayerPgns = 0;
	double* CoordsZ = 0;
	try
	{
		LayerPolygons = new TVector2d*[AmOfLayerRect];
		if(LayerPolygons == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		PtsNumbersInLayerPgns = new int[AmOfLayerRect];
		if(PtsNumbersInLayerPgns == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		CoordsZ = new double[AmOfLayerRect];
		if(CoordsZ == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		TVector3d* tRectCenPoints = RectCenPoints;
		TVector2d* tRectDims = RectDims;
		for(int i=0; i<AmOfLayerRect; i++)
		{
			TVector2d* RectanglePoints = new TVector2d[4];
			if(RectanglePoints == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

			double HalfWx = 0.5*tRectDims->x, HalfWy = 0.5*tRectDims->y;
			double Xmin = tRectCenPoints->x - HalfWx, Xmax = tRectCenPoints->x + HalfWx;
			double Ymin = tRectCenPoints->y - HalfWy, Ymax = tRectCenPoints->y + HalfWy;
			RectanglePoints[0].x = Xmin; RectanglePoints[0].y = Ymin;
			RectanglePoints[1].x = Xmax; RectanglePoints[1].y = Ymin;
			RectanglePoints[2].x = Xmax; RectanglePoints[2].y = Ymax;
			RectanglePoints[3].x = Xmin; RectanglePoints[3].y = Ymax;

			LayerPolygons[i] = RectanglePoints;
			PtsNumbersInLayerPgns[i] = 4;
			CoordsZ[i] = tRectCenPoints->z;

			tRectCenPoints++; tRectDims++;
		}
		char SetUpOK = SetUpPolyhedronsFromLayerPolygons(LayerPolygons, PtsNumbersInLayerPgns, CoordsZ, AmOfLayerRect, MagnVect, hg);

		if(LayerPolygons != 0)
		{
			for(int k=0; k<AmOfLayerRect; k++) 
			{
				delete[] (LayerPolygons[k]);
				LayerPolygons[k] = 0;
			}
			delete[] LayerPolygons;
		}
		if(PtsNumbersInLayerPgns != 0) delete[] PtsNumbersInLayerPgns;
		if(CoordsZ != 0) delete[] CoordsZ;

		if(!SetUpOK) return 0;
		return 1;
	}
	catch(...)
	{
		if(LayerPolygons != 0)
		{
			for(int k=0; k<AmOfLayerRect; k++) 
			{
				if(LayerPolygons[k] != 0) delete[] (LayerPolygons[k]);
			}
			delete[] LayerPolygons;
		}
		if(PtsNumbersInLayerPgns != 0) delete[] PtsNumbersInLayerPgns;
		if(CoordsZ != 0) delete[] CoordsZ;

		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::CheckLayerRectangleStructures(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect)
{
	if(AmOfLayerRect < 2) { Send.ErrorMessage("Radia::Error063"); return 0;}

	TVector3d* tRectCenPoints = RectCenPoints;
	TVector2d* tRectDims = RectDims;

	double CurrentZ = tRectCenPoints->z;
	char GoesDown = 0, GoesUp = 0;

	for(int i = 0; i < AmOfLayerRect; i++)
	{
		if(i > 0)
		{
			if(tRectCenPoints->z > CurrentZ)
			{
				if(i == 1) GoesUp = 1;
				else if(GoesDown || (!GoesUp)) { Send.ErrorMessage("Radia::Error063"); return 0;}
			}
			else if(tRectCenPoints->z < CurrentZ)
			{
				if(i == 1) GoesDown = 1;
				else if(GoesUp || (!GoesDown)) { Send.ErrorMessage("Radia::Error063"); return 0;}
			}
			else { Send.ErrorMessage("Radia::Error063"); return 0;}
			CurrentZ = tRectCenPoints->z;
		}
		if((tRectDims->x <= 0.) && (tRectDims->y <= 0.)) { Send.ErrorMessage("Radia::Error064"); return 0;}
		tRectCenPoints++; tRectDims++;
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::DecodeViewingOptions(const char** OptionNames, const char** OptionValues, int OptionCount, char* OptBits)
{
	radTOptionNames OptNam;
	const char** BufNameString = OptionNames;
	const char** BufValString = OptionValues;

	char& DoShowLines = OptBits[0];
	char& DoShowFaces = OptBits[1];
	char& DoShowFrameAxes = OptBits[2];
	char& DoShowSymChilds = OptBits[3];

	//Defaults
	DoShowLines = 1;
	DoShowFaces = 1;
	DoShowFrameAxes = 1;
	DoShowSymChilds = 1;

	for(int i=0; i<OptionCount; i++)
	{
		if(!strcmp(*BufNameString, OptNam.ShowLines))
		{
			if(!strcmp(*BufValString, (OptNam.ShowLinesValues)[0])) DoShowLines = 0;
			else if(!strcmp(*BufValString, (OptNam.ShowLinesValues)[1])) DoShowLines = 1;
			else if(!strcmp(*BufValString, (OptNam.ShowLinesValues)[2])) DoShowLines = 0;
			else if(!strcmp(*BufValString, (OptNam.ShowLinesValues)[3])) DoShowLines = 1;
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
		}
		else if(!strcmp(*BufNameString, OptNam.ShowFaces))
		{
			if(!strcmp(*BufValString, (OptNam.ShowFacesValues)[0])) DoShowFaces = 0;
			else if(!strcmp(*BufValString, (OptNam.ShowFacesValues)[1])) DoShowFaces = 1;
			else if(!strcmp(*BufValString, (OptNam.ShowFacesValues)[2])) DoShowFaces = 0;
			else if(!strcmp(*BufValString, (OptNam.ShowFacesValues)[3])) DoShowFaces = 1;
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
		}
		else if(!strcmp(*BufNameString, OptNam.ShowFrameAxes))
		{
			if(!strcmp(*BufValString, (OptNam.ShowFrameAxesValues)[0])) DoShowFrameAxes = 0;
			else if(!strcmp(*BufValString, (OptNam.ShowFrameAxesValues)[1])) DoShowFrameAxes = 1;
			else if(!strcmp(*BufValString, (OptNam.ShowFrameAxesValues)[2])) DoShowFrameAxes = 0;
			else if(!strcmp(*BufValString, (OptNam.ShowFrameAxesValues)[3])) DoShowFrameAxes = 1;
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
		}
		else if(!strcmp(*BufNameString, OptNam.ShowSymChilds))
		{
			if(!strcmp(*BufValString, (OptNam.ShowSymChildsValues)[0])) DoShowSymChilds = 0;
			else if(!strcmp(*BufValString, (OptNam.ShowSymChildsValues)[1])) DoShowSymChilds = 1;
			else if(!strcmp(*BufValString, (OptNam.ShowSymChildsValues)[2])) DoShowSymChilds = 0;
			else if(!strcmp(*BufValString, (OptNam.ShowSymChildsValues)[3])) DoShowSymChilds = 1;
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
		}
		else { Send.ErrorMessage("Radia::Error062"); return 0;}
		BufNameString++; BufValString++;
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTApplication::PrepareGeomPolygDataForViewing(radTVectGeomPolygon& GeomPolygons, double*& VertCoord, int& Nv, int*& VertInd, int*& PgLen, float*& PgColors, int& Npg)
{
	Npg = (int)(GeomPolygons.size());
	if(Npg <= 0) return;

	Nv=0;
	for(int i=0; i<Npg; i++) Nv += GeomPolygons[i].Nv;
	if(Nv <= 0) return;

	VertCoord = new double[3*Nv];
	VertInd = new int[Nv];
	PgLen = new int[Npg];
	PgColors = new float[3*Npg];

	double *tVertCoord = VertCoord;
	int *tVertInd = VertInd;
	int *tPgLen = PgLen;
	float *tPgColors = PgColors;

	int AbsPointCount = 0; //1; //OC240804
	for(int j=0; j<Npg; j++)
	{
		radTGeomPolygon& CurPolygon = GeomPolygons[j];
		double *tCurVert = CurPolygon.VertCoords;
		int CurNv = CurPolygon.Nv;

		*(tPgLen++) = CurNv;

		float *tColor = CurPolygon.ColRGB;
		*(tPgColors++) = *(tColor++);
		*(tPgColors++) = *(tColor++);
		*(tPgColors++) = *(tColor++);

		for(int k=0; k<CurNv; k++)
		{
			*(tVertCoord++) = *(tCurVert++);
			*(tVertCoord++) = *(tCurVert++);
			*(tVertCoord++) = *(tCurVert++);

			*(tVertInd++) = (AbsPointCount++); //correct if necessary
		}
	}
}

//-------------------------------------------------------------------------

void radTApplication::DeallocateAuxPgnViewData(double** dArr, int** iArr1, int** iArr2, float** fArr)
{
	for(int i=0; i<2; i++)
	{
		if(dArr != 0)
		{
			if(dArr[i] != 0) { delete[] dArr[i]; dArr[i] = 0;}
		}
		if(iArr1 != 0)
		{
			if(iArr1[i] != 0) { delete[] iArr1[i]; iArr1[i] = 0;}
		}
		if(iArr2 != 0)
		{
			if(iArr2[i] != 0) { delete[] iArr2[i]; iArr2[i] = 0;}
		}
		if(fArr != 0)
		{
			if(fArr[i] != 0) { delete[] fArr[i]; fArr[i] = 0;}
		}
	}
}

//-------------------------------------------------------------------------
