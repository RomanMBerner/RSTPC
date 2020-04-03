#ifndef __LIKELIHOOD_CLASS_HH__
#define __LIKELIHOOD_CLASS_HH__

#include <vector>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMinuit.h"

#include "Param.hh"

#if (not defined(__ROOTCLING__)) && (not defined(__CLING__))
namespace Analysis
{
#endif
	
	//This a container to store the results of the log-likelihood and of the parameter values
	class LikelihoodResults
	{
	public:
		double logprob;
		vector<double> par;
	};
	
	
	
	//This class is made to load from a tree the spectra saved in arrays. Devoloped mostly for analysis of simulated data.
	class LikelihoodClass
	{
	
	public:
	
		LikelihoodClass();
		
		virtual ~LikelihoodClass(){;};
		
		
		enum PropFuctionType{
			kGauss,
			kCauchy
		};
		
		
		//From here I can take all the values for the LogProb given a fixed value for the signal (or snr)
		//virtual vector<double> TreePLs(double signal){return vector<double>(0);};
	
		//Getters
		inline vector<Analysis::Param*>& GetParameters() { return fParams; };
		Analysis::Param* GetParameter(unsigned iPar);
		unsigned GetNParameters() const { return fParsN;};
		virtual double GetMaxLProb() const { return fMaxLogProb; };
		virtual vector<double> GexMaxPars() const {return fMaxPar; };
		TMinuit *GetMinuit(){return fMinuit;}//Obsolete will go away at some point
		
		
		//Setters
		virtual void SetInitPars(vector<double>* pPars=NULL);
		void SetAbsPrec(double prec){ fAbsPrec=prec;};
		void SetMaxSteps(unsigned maxSteps){ fMaxIter=maxSteps; };
		void SetTrackVals(bool trackvals=true){ fTrackVals=trackvals; };
		void SetVerbosity(unsigned verb){ fVerbosity=verb; };
		void SetEngineVerb(unsigned verb){ fEngineVerb=verb; };
		void SetStepProposalFunction(LikelihoodClass::PropFuctionType prop){fPropFncType=prop;};
		
		//Miscellaneous
		void ForceParVal(bool flag=true){ fParForced=flag; };
		void ForceMCMC(bool flag=true){ fMCMCforced=flag; };
		
		
	
	private:
		
		LikelihoodClass::PropFuctionType fPropFncType;
		
		
	protected:
		
		virtual void MetropolisMLE(vector<double> *initpar);
		virtual void MetropolisMLE(){MetropolisMLE(NULL);};
		
		//This should be implemented for a reasonable parameter start when the starting parameters are not given
		virtual vector<double> SelfParInit()=0;
		
		//This method is responsible to prepare the data, create a new tree and make a full analysis calling other methods.
		virtual int DefineParameters()=0;
	
		//This computes the real likelihood must be implemented in a concrete class
		virtual double LogProb(const vector<double>& par)=0;
		
		//This is the method that returns the new proposed parameters
		virtual void MetrProp(const vector<double>& par, vector<double>& par_new);
		
		//This is the Minuit FCN function (Obsolete will go away at some point)
		static void MinuitFCN(int& npar, double* grad, double& fval, double* par, int flag);
		
		unsigned m_nPar;
		TRandom3 *fRndGen;
	
		//Maximization stuff
		TMinuit *fMinuit; //Obsolete will go away at some point
		vector<double> fMaxPar;
		double fMaxLogProb;
	
		double fMinuitArglist[2];
	
		bool fParInit; //Flag to check if the parameters are properly initialized
		
		vector<double> *fInitPars;
		
		vector<Param*> fParams; //Here are stored ALL the parameters
		unsigned fParsN;
		
		vector<double> *fParSteps;
		
		bool fParForced;
		
		bool fMCMCforced;
		
		double fAbsPrec;
		unsigned fMaxIter;
		bool fTrackVals;
		
		unsigned fVerbosity;
		unsigned fEngineVerb;
	};


	//I need this in order to use properly Minuit
	class LikelihoodClassHolder
	{
	private:
		Analysis::LikelihoodClass * global_this;

	public:
	
		LikelihoodClassHolder():global_this(NULL){}//This could also be private
		
		static Analysis::LikelihoodClass* instance(Analysis::LikelihoodClass * obj = NULL){
			static LikelihoodClassHolder inst;
			if (obj) inst.global_this = obj;
			return inst.global_this;
		}
	};
	
#if (not defined(__ROOTCLING__)) && (not defined(__CLING__))
}//End of analysis namespace
#endif



#endif