#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


#include "Param.hh"

using namespace std;

namespace Analysis
{
	Param::Param(string name, bool trackvals):
		fName(name),
		fUpperLimSet(false),
		fLowerLimSet(false),
		fUpperLimit(0),
		fLowerLimit(0),
		fHist(NULL),
		fFixed(false),
		fValue(0),
		fLogProb(0),
		fNbins(1000),
		fHistoFlag(false),
		fMCMCstepSizeFlag(false),
		fTrackVals(trackvals)
	{
	
		//Create 2 empty vectors (the likelihood value and the parameters values)
		fValsChain = vector<vector<double> >(2);
	}


	Param::Param(string name, double lowLimit, double upLimit, bool trackvals):
		fName(name),
		fUpperLimSet(true),
		fLowerLimSet(true),
		fUpperLimit(upLimit),
		fLowerLimit(lowLimit),
		fHist(NULL),
		fFixed(false),
		fValue(0),
		fLogProb(0),
		fNbins(1000),
		fHistoFlag(false),
		fMCMCstepSizeFlag(false),
		fTrackVals(trackvals)
	{
	
		//Create 2 empty vectors (the parameter value and the likelihood value)
		fValsChain = vector<vector<double> >(2);
	
	}


	void Param::SetLimits(double lowerlimit, double upperlimit)
	{
		if(lowerlimit>upperlimit){
			cerr << "\nParam::SetLimits(...) --> ERROR: Tryng to set lower limit higher than upper limit for parameter \"" << fName << "\"" << endl << endl;
			return;
		}
	
		fLowerLimit=lowerlimit;
		fLowerLimSet=true;
		fUpperLimit=upperlimit;
		fLowerLimSet=true;
	}

	void Param::SetUpperLimit(double limit)
	{
	
		if(fLowerLimSet){
			if(fLowerLimit>limit){
				cerr << "\nParam::SetUpperLimit(...) --> ERROR: Tryng to set lower limit higher than upper limit for parameter \"" << fName << "\"" << endl << endl;
				return;
			}
		}
	
	
		fUpperLimit=limit;
		fUpperLimSet=true;
	}

	void Param::SetLowerLimit(double limit)
	{
	
		if(fUpperLimSet){
			if(limit>fUpperLimit){
				cerr << "\nParam::SetLowerLimit(...) --> ERROR: Tryng to set lower limit higher than upper limit for parameter \"" << fName << "\"" << endl << endl;
				return;
			}
		}
	
	
		fLowerLimit=limit;
		fLowerLimSet=true;
	}

	bool Param::SetValue(double value, double logprob)
	{
	
		bool flag=true;
	
		if(fLowerLimSet){
			if(value < fLowerLimit){
				fValue = fLowerLimit;
				flag = false;
			}
		}
		if(fUpperLimSet){
			if(value > fUpperLimit){
				fValue = fUpperLimit;
				flag = false;
			}
		}
		if(flag){
			fValue = value;
			fLogProb = logprob;
		}
	
	
		if(fTrackVals){
			fValsChain.at(0).push_back(fValue);
			fValsChain.at(1).push_back(fLogProb);
		}
	
		return flag;
	
	}


	void Param::ForceValue(double value, double logprob)
	{
	
	
		fValue = value;
		fLogProb = logprob;
	
	
		if(fTrackVals){
			fValsChain.at(0).push_back(fValue);
			fValsChain.at(1).push_back(fLogProb);
		}
	
		return;
	
	}


	int Param::IsOutOfRange(double value) const
	{
	
		if(fLowerLimSet){
			if( (fLowerLimit >= value) ) return -1;//I am below the allowed range
		}
	
		if(fUpperLimSet){
			if( (fUpperLimit <= value) ) return 1;//I am above the allowed range
		}
	
		return 0;
	}


	void Param::Fix(double value)
	{
		if(fTrackVals&&(fValsChain.size()>0)){
			//cout << "\nGatorParam::Fix(...) cannot fix the parameter:  " << fName << endl;
			fTrackVals = false; //When I fix the parameter I don't track it any more!
		}
	
		fFixed = true;
		fValue = value;
	
		//cout << '\nFixing "'  << fName << '" parameter to ' << fValue;
		return;
	}

}//End of Analysis namespace

