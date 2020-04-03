#ifndef __BCHISTOFITTERFAST_HH__
#define __BCHISTOFITTERFAST_HH__

#include "TF1.h"
#include "TH1D.h"

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCHistogramFitter.h>

#include <vector>


using namespace std;

class BcHistoFitterFast: public BCHistogramFitter
{
private:

	vector<double> fHistVec;
	vector<double> fLowEdges;
	vector<double> fUpEdges;

	TF1 *fFitFunction;
	TH1D* fHistogram;
	int nBins;
	
	bool fFlagIntegration;

public:
	BcHistoFitterFast(TH1D * hist, TF1 * func);
	virtual ~BcHistoFitterFast(){;};

	virtual double LogLikelihood(const std::vector<double> & parameters);
	
	void SetFlagIntegration(bool flag){ fFlagIntegration = flag; };


};



#endif