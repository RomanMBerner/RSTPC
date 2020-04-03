#ifndef __PARAM_HH__
#define __PARAM_HH__

#include <string>
#include <vector>
#include "TH1D.h"

using namespace std;

#if (not defined(__ROOTCLING__)) && (not defined(__CLING__))
namespace Analysis
{
#endif
	
	class Param{

	public:
		Param(string name, bool trackvals=false);
		Param(string name, double lowEdge, double upEdge, bool trackvals=false);
		~Param(){;};
	
		//Getters
		const string & GetName() const { return fName; };
		inline double GetUpperLimit() const { return fUpperLimit; };
		inline double GetLowerLimit() const { return fLowerLimit; };
		inline bool GetUpperLimit(double& val) const { val= fUpperLimit; return IsUpperLimit(); };
		inline bool GetLowerLimit(double& val) const { val= fLowerLimit;  return IsLowerLimit(); };
		//bool FillHistogram(); //This must be implemented later ()
		inline double GetValue() const { return fValue; };
		inline double GetLogProb() const { return fValue; };
		unsigned GetNbins() const { return fNbins; };
		TH1D* GetHisto() const { return fHist;};
		vector<vector<double> > GetTrackingVals() const { return fValsChain; };
		double GetMCMCStepSize(){return fMCMCstepSize;}
	
		//Setters
		void SetName(string name) { fName = name; };
		bool SetValue(double value, double likelihood=0);
		void SetUpperLimit(double limit);
		void SetLowerLimit(double limit);
		void SetLimits(double lowerlimit, double upperlimit);
		void SetHist(bool flag=true) { fHistoFlag=flag; };
		void SetNbins(unsigned nbins) { fNbins = nbins; };
		void SetMCMCStepSize(double step){ fMCMCstepSize=step;  fMCMCstepSizeFlag=true; };
		void ForceValue(double value, double logprob);
	
		//Misc
		inline bool FillHisto() const { return fHistoFlag; };
		inline bool IsFixed() const { return fFixed; };
		inline bool IsUpperLimit() const {return fUpperLimSet;}
		inline bool IsLowerLimit() const {return fLowerLimit;}
		inline bool IsTracking() const { return fTrackVals; };
		void TrackVals(bool trackval=true) { fTrackVals=trackval; };
		void Fix(double value);
		void Unfix(){ fFixed = false; };
		//void grams(bool flag) { fFillHistograms = flag; };
		bool IsAtLimit(double value) const { return ((fLowerLimit == value) || (value == fUpperLimit)); } ;
		int IsOutOfRange(double value) const;
		bool IsMCMCextStepSize(){return fMCMCstepSizeFlag;};
	
	private:
	
		string fName;
	
		double fUpperLimit;
		double fLowerLimit;
		
		bool fUpperLimSet;
		bool fLowerLimSet;
		
		bool fFixed; //If the value is fixed
		double fValue; //The current value of the parameter
		double fLogProb;
		bool fHistoFlag; //If the histogram is filled
		unsigned fNbins; //Number of bins for the histogram
		bool fTrackVals;
		
		double fMCMCstepSize;//Not initialized. It is used only if this parameter is set from outside with the public method.
		bool fMCMCstepSizeFlag;//It is true when the step size is set from the public setter method.
	
		//Here I could dump the profiled likelihood or the marginal posterior pdf for the parameter.
		TH1D* fHist;
		//In the first vector are saved the values of the parameter and in the second the values of the likelihood (or the log-likelihood or the posterior pdf values)
		vector<vector<double> > fValsChain;
	
	};
	
#if (not defined(__ROOTCLING__)) && (not defined(__CLING__))
}//End of Analysis namespace
#endif



#endif