
#ifndef __IGINTERFSTR_H
#define __IGINTERFSTR_H

//*************************************************************************

struct TDataWaveMD {
	char* pData;
	char DataType[2]; // 'f'|'d'|'cf'|'cd'
	long AmOfDims;
	long DimSizes[10];
	double DimStartValues[10];
	double DimSteps[10];
	char DimUnits[10][255];
	char DataUnits[255];
	char DataName[255];

	double* OutDoubleDataPtr(bool& ArrWasAlloc)
	{
		ArrWasAlloc = false;
		if(DataType[0] == 'd') return (double*)pData;
		else if(DataType[0] == 'f')
		{
			long TotAmOfData = AmOfDataValues();
			if(TotAmOfData <= 0) return 0;
			double *pOutData = new double[TotAmOfData];

			double *tOutData = pOutData;
			float *tfData = (float*)pData;
			for(long i=0; i<TotAmOfData; i++) *(tOutData++) = *(tfData++);
			ArrWasAlloc = true;
			return pOutData;
		}
		else return 0;
	}

	long AmOfDataValues()
	{
		if(AmOfDims <= 0) return 0;

		long OutAmOfData = DimSizes[0];
		if(OutAmOfData <= 0) return 0;
		for(int i=1; i<AmOfDims; i++) 
		{
			if(DimSizes[i] <= 0) break;
			OutAmOfData *= DimSizes[i];
		}
        return OutAmOfData;
	}

	void UpdateDataValues(double* pdInData)
	{
		if(pdInData == 0) return;
		long TotAmOfData = AmOfDataValues();
		if(TotAmOfData <= 0) return;

		double *tdInData = pdInData;
		if(DataType[0] == 'd')
		{
            double *tLocData = (double*)pData;
			for(long i=0; i<TotAmOfData; i++) *(tLocData++) = *(tdInData++);
		}
		else if(DataType[0] == 'f')
		{
            float *tLocData = (float*)pData;
			//for(long i=0; i<TotAmOfData; i++) *(tLocData++) = (double)(*(tdInData++));
			for(long i=0; i<TotAmOfData; i++) *(tLocData++) = (float)(*(tdInData++)); //OC080613
		}
	}
};

//*************************************************************************

#endif