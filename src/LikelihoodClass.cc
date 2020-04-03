#include <cmath>
#include <sstream>

#include <vector>
#include <string>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMinuit.h"


#include "LikelihoodClass.hh"



using namespace std;

namespace Analysis{
	
	LikelihoodClass::LikelihoodClass(  ):
		fParForced(false),
		fMCMCforced(false),
		fAbsPrec(1e-6),
		fMaxIter(10000),
		fTrackVals(false),
		fParsN(0),
		fParInit(false),
		fInitPars(NULL),
		fParSteps(NULL),
		fMaxLogProb(0.),
		fVerbosity(1),
		fEngineVerb(0),
		fPropFncType(kGauss)
	{
	
		fRndGen = new TRandom3(0);
	
		fMinuit = NULL;
	
		//Initialization of the arglist used for the MIGRAD routine.
		//To have a better understanding make an interactive session in root and use TMinuit::mnhelp("MIGrad")
		//The first argument for MIGRAD is the maximum number of function calls
		fMinuitArglist[0] = 1000; //Optional
		//This is the tollerance
		fMinuitArglist[1] = 0.01; //Optional
	
	
	}


	Param* LikelihoodClass::GetParameter(unsigned iPar)
	{
		if(iPar<2*m_nPar+1){
			return fParams.at(iPar);
		}else{
			cerr << "\nOut of range: parameter index " << iPar << " doesn't exist. Abort program!\n" << endl;
			return NULL;
		}
	};


	void LikelihoodClass::SetInitPars(vector<double>* pPars)
	{
		if(fInitPars){
			delete fInitPars;
			fInitPars=NULL;
		}
		if(!pPars){
			fInitPars=NULL;
			return;
		}
	
		if(pPars->size()==fParams.size()) fInitPars=pPars;
		return;
	}


	void LikelihoodClass::MetropolisMLE(vector<double> *initpar)
	{
	
		cout << "Entering in LikelihoodClass::MetropolisMLE() function" << endl;
	
		//HERE starts the Metropilis search of MLE 
		vector<double> par(fParsN);
		vector<double> par_new(fParsN);
	
		vector<bool> fixedpars(fParsN);
	
		if(initpar){
			SetInitPars(initpar);
		}
	
		//Initialize the parameters (to something reasonable)
		if(!fInitPars){
			if(fEngineVerb>=1){
				cout << "=> Self initialization of the parameters." << endl;
			}
			par =SelfParInit();
		}
	
	
	
		if(fEngineVerb>=1){
			cout << "=> Checking which are the fixed parameters." << endl;
		}
		unsigned nFixed = 0;
		for(unsigned iPar=0; iPar<fParsN; iPar++){
			if(fParams.at(iPar)->IsFixed()){
				fixedpars.at(iPar) = true;
				par.at(iPar) = fParams.at(iPar)->GetValue();
				par_new.at(iPar) = fParams.at(iPar)->GetValue();
				nFixed += 1;
			}else{
				fixedpars.at(iPar) = false;
			}
		}
		if(fEngineVerb>=1){
			cout << "=> " << nFixed << " fixed parameters over " << fParsN << " total parameters." << endl;
		}
	
		if(nFixed>=fParsN){
			fMaxLogProb = LogProb(par);
			fMaxPar = par;
			if(fEngineVerb>=2){
				cout << "==> All Parameter fixed." << endl;
				cout << "==> Max LogProb: " << fMaxLogProb << endl;
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					cout << "==>Parameter <" << fParams.at(iPar)->GetName() << ">: " << par.at(iPar) << endl;
				}
			}
			return;
		}
	
	
		double logprob = LogProb(par);
		double logprob_new;
	
		unsigned iStep = 0;
	
		fMaxLogProb = logprob;
		fMaxPar = par;
	
		fRndGen->SetSeed();
	
		if(fEngineVerb>=1){
			cout << "=> Starting the Metropolis loop for the MLE finding with " << fMaxIter << " steps." << endl;
		}
		do{
			if(/*iStep%1000==0 &&*/ fEngineVerb>=2){
				cout << "==>Step: "<< iStep << endl;
				if(fEngineVerb>=3){
					for(unsigned iPar=0; iPar<fParsN; iPar++){
						cout << "===> Parameter <" << fParams.at(iPar)->GetName() << ">: " << par.at(iPar) << endl;
					}
				}
				cout << "===> Logprob: " << logprob << endl;
			}
		
			bool accepted = false;
			MetrProp(par,par_new);
			logprob_new = LogProb(par_new);
		
			if(/*iStep%1000==0 &&*/ fEngineVerb>=3){
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					cout << "===> Proposed parameter <" << fParams.at(iPar)->GetName() << ">: " << par_new.at(iPar) << endl;
				}
				cout << "===> New LogProb:  " << logprob_new << endl;
			}
		
			if(logprob_new>=logprob){
				accepted = true;
			}else{
				double r = fRndGen->Rndm();
				if( (logprob_new-logprob)>log(r) ){
					accepted = true;
					if(fEngineVerb>=3){
						cout << "===> Accepted with odds = " << r << endl;
					}
				}
			
			}
		
			//Check if the parameter is inside the allowed region if it is not in the range
			if(fEngineVerb>=3){
				cout << "===> Checking if new parameters are in range." << endl;
			}
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				if(fParams.at(iPar)->IsOutOfRange(par_new.at(iPar)) != 0 ) accepted = false;
			}
		
			//Update step
			if(accepted){
				if(/*iStep%1000==0 &&*/ fEngineVerb>=3){
					cout << "===> Accepted (in range)." << endl;
				}
				par = par_new;
				logprob = logprob_new;
			
				if(logprob>fMaxLogProb){
					fMaxLogProb = logprob;
					fMaxPar = par;
				}
			
				if(fTrackVals){
					for(unsigned iPar=0; iPar<fParsN; iPar++){
						if(fParForced){
							fParams.at(iPar)->ForceValue(par.at(iPar),logprob);
						}else{
							fParams.at(iPar)->SetValue(par.at(iPar),logprob);
						}
					}
				}
				iStep++;
			}
		}while(iStep<fMaxIter);
		if(fEngineVerb>=1){
			cout << "=> End of Metropolis loop for the MLE finding." << endl;
		}
	
		if(fEngineVerb>=1){
			cout << "\n=> Max LogProb after MCMC: " << fMaxLogProb << endl;
			//cout << "\n=> Maximizing with Minuit" << endl;
			//cout << "  --------------------------" << endl;
		}
	
		/*
		//As the "MinuitFNC" is a static function and uses variables of this class it would be good to have always a real instance of this class, otherwise I might use variables that are not actually instanced (look also the comments in the "MuinuitFNC").
		XurichAnalysis::LikelihoodClassHolder::instance(this);
	
		//cout << "\nTMinuit pointer: " << fMinuit << endl;
	
		if(fMinuit) delete fMinuit;
	
		fMinuit = new TMinuit(fParsN);
	
		//Setverbose in testing phase
		if(fEngineVerb<=1){
			fMinuit -> SetPrintLevel(-1);
		}
		else if(fEngineVerb==2){
			fMinuit -> SetPrintLevel(0);
		}
		else if(fEngineVerb>2){
			fMinuit -> SetPrintLevel(1);
		}
	
		fMinuit -> SetFCN(&XurichAnalysis::LikelihoodClass::MinuitFCN);
	
		fMinuit -> SetErrorDef(0.5);
	
		int minuitFlag, minuitierr;
		for(unsigned iPar=0; iPar<fParsN; iPar++){
			double width = fParams.at(iPar)->GetUpperLimit()-fParams.at(iPar)->GetLowerLimit();
			fMinuit -> mnparm(iPar, (fParams.at(iPar)->GetName()).c_str(), fMaxPar.at(iPar), width/1000., fParams.at(iPar)->GetLowerLimit(), fParams.at(iPar)->GetUpperLimit(), minuitFlag);
			if(fParams.at(iPar)->IsFixed()){
				//cout << "\nPar " << fParams.at(iPar)->GetName() << " is fixed" << endl << endl;
				fMinuit->FixParameter(iPar);
			}
		}
	
	
		fMinuit -> mnexcm("MIGRAD", fMinuitArglist, 2, minuitFlag);
	
		for(unsigned iPar=0; iPar<fParsN; iPar++){
			double errors;
			fMinuit -> GetParameter(iPar, par_new.at(iPar), errors);
		}
	
	
	
		logprob_new = LogProb(par_new);
		if(fMaxLogProb < logprob_new){
			fMaxLogProb = logprob_new;
			fMaxPar = par_new;
		
		}
		if(fEngineVerb>=2){
			cout << "==> Max LogProb after Minuit: " << fMaxLogProb << endl;
		}
	
	
		if(!fTrackVals){
			for(unsigned iPar=0; iPar<m_nPar; iPar++){
				if(fParForced){
					fParams.at(iPar)->ForceValue(fMaxPar.at(iPar),fMaxLogProb);
				}else{
					fParams.at(iPar)->SetValue(fMaxPar.at(iPar),fMaxLogProb);
				}
				//fParams.at(iPar)->SetValueForce( par.at(iPar),logprob );
			}
		}
	
		if(fVerbosity>=2){
			cout << "==> Max LogProb: " << fMaxLogProb << endl;
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				cout << "==> Parameter <" << fParams.at(iPar) << ">: " << par.at(iPar) << endl;
			}
		}
		*/
	
		cout << "Exiting from LikelihoodClass::MetropolisMLE() function." << endl;
	
		return;
	}


	void LikelihoodClass::MetrProp(const vector<double>& par, vector<double>& par_new)
	{
		LikelihoodClass::PropFuctionType lPropFncType = fPropFncType;
		if( (fPropFncType!=LikelihoodClass::kGauss) && (fPropFncType!=LikelihoodClass::kCauchy) ) lPropFncType = LikelihoodClass::kGauss;
		for(unsigned iPar=0; iPar<fParsN; iPar++){
			if(lPropFncType==LikelihoodClass::kGauss){
				//Using the Gaussian distribution
				if(!fParams.at(iPar)->IsFixed()){
					par_new.at(iPar) = par.at(iPar) + fRndGen->Gaus(0., fParSteps->at(iPar));
				}
			}else if(lPropFncType==LikelihoodClass::kCauchy){
				//Using the Cauchy distribution
				double x = fRndGen->Rndm();
				par_new.at(iPar) = par.at(iPar) + (fParSteps->at(iPar))*(TMath::Tan((x-0.5)*TMath::Pi()));
			}
		}
	
		return;
	
	}


	void LikelihoodClass::MinuitFCN(int& npar, double* grad, double& fval, double* par, int flag)
	{
	
		//REMEMBER THAT THIS IS A STATIC METHOD
		static vector<double>* parameters= NULL;
	
		//The method GetNParameters() is not static, hence if there isn't any instance of the class this method cannot be called
	
		LikelihoodClass *pThisObj = (LikelihoodClass*)LikelihoodClassHolder::instance();
	
		int nparameters = (int)pThisObj->GetNParameters();
	
		if(!parameters){
			parameters = new vector<double>(nparameters);
		}
	
		//parameters.resize(nparameters, 0.0);
	
		//copy the parameters from the input array to the vector
		std::copy(par, par + nparameters, parameters->begin());
	
		//To evaluate change the sign as Minuit searches the minimum
		fval = - pThisObj->LogProb((*parameters));
		//Also here LogProb() is not a static method and uses non static members of the class. All this stuff can be used only from an instanced object.
	
	
		return;
	}
	
}//End of Analysis namespace

