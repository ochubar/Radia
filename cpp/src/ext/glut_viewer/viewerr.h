
#ifndef __VIEWERROR_H
#define __VIEWERROR_H

#include <string>
#include <vector>

//#if defined _MSC_VER 
using namespace std;
#pragma warning(disable : 4786) // to suppress annoying warning from STL
//#endif

//-------------------------------------------------------------------------
//Error numbers are supposed to be thrown by application (catched and returned by API functions)
//Warning numbers are supposed to be inserted into internal container m_WarnNosVect via AddWarning function
class CErrWarn {

	static string m_error[];
	static string m_warning[];
	static vector<int> m_WarnNosVect;

	static int m_BaseErrNo, m_BaseWarNo;

public:
	
	static int GetErrorSize(int ErrNo)
	{
		int i = ErrNo - m_BaseErrNo;
		if(i < 0) return (int)(m_error[0].size());
		try { return (int)(m_error[i].size());}
		catch(exception e) { return (int)(m_error[0].size());}
	}
	static const char* GetError(int ErrNo)
	{
		int i = ErrNo - m_BaseErrNo;
		if(i < 0) return m_error[0].c_str();
		try { return m_error[i].c_str();}
		catch(exception e) { return m_error[0].c_str();}
	}
	static int GetWarningSize(int WarnNo)
	{
		int i = WarnNo - m_BaseWarNo;
		if(i < 0) return (int)(m_warning[0].size());
		try { return (int)(m_warning[i].size());}
		catch(exception e) { return (int)(m_warning[0].size());}
	}
	static const char* GetWarning(int WarnNo)
	{
		int i = WarnNo - m_BaseWarNo;
		if(i < 0) return m_warning[0].c_str();
		try { return m_warning[i].c_str();}
		catch(exception e) { return m_warning[0].c_str();}
	}

	static void CopyErrorText(int ErrNo, char* cBuf)
	{
		if(cBuf == 0) return;
		int i = ErrNo - m_BaseErrNo;
        const char* cTmp0 = m_error[0].c_str();
		if(i < 0) { strcpy(cBuf, cTmp0); return;}

		try 
		{ 
			const char* cTmp = m_error[i].c_str();
			strcpy(cBuf, cTmp);
		}
		catch(exception e) { strcpy(cBuf, cTmp0);}
	}
	static void CopyWarningText(int WarnNo, char* cBuf)
	{
		if(cBuf == 0) return;
		int i = WarnNo - m_BaseWarNo;
        const char* cTmp0 = m_warning[0].c_str();
		if(i < 0) { strcpy(cBuf, cTmp0); return;}

		try 
		{
			const char* cTmp = m_warning[i].c_str();
			strcpy(cBuf, cTmp);
		}
		catch(exception e) { strcpy(cBuf, cTmp0);}
	}
	static void CopyAllWarningText(char* cBuf)
	{
		if(cBuf == 0) return;
		int AmOfWarn = (int)(m_WarnNosVect.size());
		if(AmOfWarn <= 0) return;

		for(int i=0; i<AmOfWarn; i++)
		{
			int WarnNo = m_WarnNosVect[i];
			if(WarnNo < 0) continue;
            const char* cTmp0 = m_warning[WarnNo].c_str();
			if(i != 0) strcat(cBuf, "\r\n"); //to check
			strcat(cBuf, cTmp0);
		}
		m_WarnNosVect.clear();
	}

	static void AddWarning(int WarnNo)
	{
		for(vector<int>::iterator iter = m_WarnNosVect.begin(); iter != m_WarnNosVect.end(); ++iter)
		{
			if(*iter == WarnNo) return; //Don't insert more than one same warning
		}
		m_WarnNosVect.push_back(WarnNo);
	}

	static int OutFirstWarningNo()
	{//returns the number of the first warning in the stack vector, but does not erase it from the stack
		if(m_WarnNosVect.empty()) return 0;
		return m_WarnNosVect[0];
	}
};

//-------------------------------------------------------------------------

#endif
