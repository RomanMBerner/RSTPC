#ifndef HISTO_MANIPULATORS_CC
#define HISTO_MANIPULATORS_CC

#include "HistoManipulators.hh"


vector<double> GaussianFilter(vector<double> histo, int KernelLen, double sigma){
	
	int HistoLenght = histo.size();
	
	//REMEMBER ALL THE QUANTITIES ARE IN BIN UNITS
	vector<double> filteredHisto(HistoLenght);
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	
	TF1 *gausFunc = new TF1("gausFunc","exp(-pow(x/[0],2)/2)",-KernelLen,KernelLen);
	gausFunc->SetParameter(0,sigma);
	
	//Build the mask
	vector<double> gaussMask(KernelLen);
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		gaussMask.at(iSamp) = gausFunc->Eval((double)(iSamp-centralPos));
	}
	
	//Make the convolution
	for(int nSamp=0; nSamp<HistoLenght; nSamp++){
		
		double sum = 0;
		
		for(int kSamp=0; kSamp<KernelLen; kSamp++){
			if( (nSamp - centralPos + kSamp >= 0) && (nSamp - centralPos + kSamp < HistoLenght) ){
				sum += histo.at(nSamp - centralPos + kSamp) * gaussMask.at(kSamp);
			}
		}
		
		filteredHisto.at(nSamp) = sum;
	}
	
	if(gausFunc) delete gausFunc;
	
	return filteredHisto;
}


vector<double> DerivativeFilter(vector<double> histo){
	
	int HistoLenght = histo.size();
	
	//REMEMBER ALL THE QUANTITIES ARE IN BIN UNITS
	vector<double> derivativeHisto(HistoLenght);
	
	for(int iSamp=0; iSamp<HistoLenght; iSamp++){
		
		if((iSamp-1 < 0) || (iSamp+1 >= HistoLenght)){
			derivativeHisto.at(iSamp) = 0;
		}else{
			derivativeHisto.at(iSamp) = (histo.at(iSamp+1)-histo.at(iSamp-1))/2;
		}
		
	}
	
	return derivativeHisto;
}


void ScaleHistoByAbsMax(TH1D* histo){
	
	double max = TMath::Abs( histo->GetBinContent( histo->GetMaximumBin() ) );
	double min = TMath::Abs( histo->GetBinContent( histo->GetMinimumBin() ) );
	
	if(max>=min){
		histo->Scale(1./max);
		cout << "\nHistogram \"" << histo->GetName() << "\" absolute maximum: " << max << endl;
	}else{
		histo->Scale(1./min);
		cout << "\nHistogram \"" << histo->GetName() << "\" absolute maximum: " << min << endl;
	}
	
	return;
}



//////////////////////////////////////////////////
///Find the maximum positive derivative abscissa
EdgeResult FindMaxInRange(TH1D *histo, double xmin, double xmax, double thrs){
	
	EdgeResult res;
	
	int lowBin = histo->FindBin(xmin);
	int upBin = histo->FindBin(xmax);
	
	int maxBin = lowBin;
	double maxval = histo->GetBinContent(maxBin);
	for(int iBin=lowBin; iBin<=upBin; iBin++){
		if(maxval<histo->GetBinContent(iBin)){
			maxval=histo->GetBinContent(iBin);
			maxBin = iBin;
		}
	}
	
	if(maxval<=0){
		cout << endl << "There is no positive derivative maximum in the specified range." << endl;
		return res;
	}
	
	//Get low edge of the maximum
	int leftBin = maxBin;
	while(true){
		leftBin--;
		double val = histo->GetBinContent(leftBin);
		if(val<=thrs*maxval){
			break;
		}
	}
	
	//Get up edge of the maximum
	int rightBin = maxBin;
	while(true){
		rightBin++;
		double val = histo->GetBinContent(rightBin);
		if(val<=thrs*maxval){
			break;
		}
	}
	
	cout << endl << "Maximum at: " << histo->GetBinCenter(maxBin) << endl;
	cout << "Maximum low limit: " << histo->GetBinLowEdge(leftBin) << endl;
	cout << "Maximum up limit: " << histo->GetBinLowEdge(rightBin+1) << endl << endl;
	
	res.lowEdge = histo->GetBinLowEdge(leftBin);
	res.midEdge = histo->GetBinCenter(maxBin);
	res.upEdge = histo->GetBinLowEdge(rightBin+1);
	
	return res;
}


//////////////////////////////////////////////////
///Find the minimum negative derivative abscissa
EdgeResult FindMinInRange(TH1D *histo, double xmin, double xmax, double thrs){
	
	EdgeResult res;
	
	int lowBin = histo->FindBin(xmin);
	int upBin = histo->FindBin(xmax);
	
	
	int minBin = lowBin;
	double minval = histo->GetBinContent(lowBin);
	for(int iBin=lowBin; iBin<=upBin; iBin++){
		//cout << "minBin: " << iBin << "\tminval: " << minval << endl;
		if(minval>histo->GetBinContent(iBin)){
			//cout << " New minval: " << histo->GetBinContent(iBin) << "\tat bin: " << iBin << endl;
			minval=histo->GetBinContent(iBin);
			minBin = iBin;
		}
	}
	
	if(minval>=0){
		cout << endl << "There is no negative derivative minimum in the specified range." << endl;
		return res;
	}
	
	//Get low edge of the minimum
	int leftBin = minBin;
	while(true){
		leftBin--;
		double val = histo->GetBinContent(leftBin);
		if(val>=thrs*minval){
			break;
		}
	}
	
	//Get up edge of the minimum
	int rightBin = minBin;
	while(true){
		rightBin++;
		double val = histo->GetBinContent(rightBin);
		if(val>=thrs*minval){
			break;
		}
	}
	
	cout << endl << "Minimum at: " << histo->GetBinCenter(minBin) << endl;
	cout << "Minimum low limit: " << histo->GetBinLowEdge(leftBin) << endl;
	cout << "Minimum up limit: " << histo->GetBinLowEdge(rightBin+1) << endl << endl;
	
	res.lowEdge = histo->GetBinLowEdge(leftBin);
	res.midEdge = histo->GetBinCenter(minBin);
	res.upEdge = histo->GetBinLowEdge(rightBin+1);
	
	return res;
}



vector<double> GetValuesFromTree(TTree* T1, string variable, TCut cut){
	
	vector<double> returnVec;
	
	if(!T1) return returnVec;
	
	cout << "\nGetValuesFromTree(...):\n" << "Variable: " << variable << endl << "Used cut: " << cut << endl;
	
	T1->SetEstimate(T1->GetEntries());
	
	T1->Draw(variable.c_str(), cut, "goff");
	Long64_t nPnts = T1->GetSelectedRows();
	
	double *p_V1 = T1->GetV1();
	
	returnVec.assign(p_V1, p_V1+nPnts);
	
	cout << "\nGetValuesFromTree(...): exiting." << endl;
	
	return returnVec;
}


void DetermineProbInt(vector<double> &values, double prob, double &low, double &up){
	DetermineProbInt(values, prob, low, up, false);
}


void DetermineProbInt(vector<double> &values, double prob, double &low, double &up, bool sorted){
	
	bool debug = false;
	
	if( !(values.size()>1) ) return;
	
	if(!sorted) std::sort(values.begin(), values.end());
	
	Long64_t nPnts = values.size();
	Long64_t nVals = (Long64_t)(nPnts*prob+0.5);
	Long64_t nSteps = nPnts-nVals+1;
	
	if(debug){
		cout << "\nTot events: " << nPnts << endl;
		cout << "Prob: " << prob << endl;
		cout << "nVals: " << nVals << endl;
		cout << "nSteps: " << nSteps << endl;
	}
	
	low = values.at(0);
	up = values.at(nVals-1);
	
	for(Long64_t iStep=0; iStep<nSteps; iStep++){
		double interval = values.at(iStep+nVals-1)-values.at(iStep);
		if(interval<(up-low)){
			low = values.at(iStep);
			up = values.at(iStep+nVals-1);
		}
	}
	
	if(debug){
		cout << "LowLim: " << low << endl;
		cout << "UpLim: " << up << endl;
	}
	
	
	return;
}


void DetermineProbInt(TTree* T1, string variable, TCut cut, double prob, double &low, double &up){
	
	bool debug = false;
	
	if(!T1) return;
	
	T1->SetEstimate(T1->GetEntries());
	
	T1->Draw(variable.c_str(), cut, "goff");
	Long64_t nPnts = T1->GetSelectedRows();
	
	double *p_V1 = T1->GetV1();
	
	vector<double> *v_points = new vector<double>(p_V1, p_V1+nPnts);
	
	std::sort(v_points->begin(), v_points->end());
	
	Long64_t nVals = (Long64_t)(nPnts*prob+0.5);
	
	Long64_t nSteps = nPnts-nVals+1;
	
	if(debug){
		cout << "\nCut used: " << cut << endl;
		cout << "Tot events: " << v_points->size() << endl;
		cout << "Prob: " << prob << endl;
		cout << "nVals: " << nVals << endl;
		cout << "nSteps: " << nSteps << endl;
	}
	
	
	low = v_points->at(0);
	up = v_points->at(nVals-1);
	for(Long64_t iStep=0; iStep<nSteps; iStep++){
		double interval = v_points->at(iStep+nVals-1)-v_points->at(iStep);
		if(interval<(up-low)){
			low = v_points->at(iStep);
			up = v_points->at(iStep+nVals-1);
		}
	}
	
	if(debug){
		cout << "LowLim: " << low << endl;
		cout << "UpLim: " << up << endl;
	}
	
	
	return;
}


TEventList* IntersectLists(vector<string> ListsName, map<string, TEventList*> *listsMap, string outname){
	
	TEventList *outlist = NULL;
	
	if(!(ListsName.size()>0)) return NULL;
	
	vector<TEventList*> listsvector;
	
	for(unsigned iList=0; iList<ListsName.size(); iList++){
		if(listsMap->find(ListsName.at(iList)) != listsMap->end()){
			TEventList* evList = (*listsMap)[ListsName.at(iList)];
			listsvector.push_back( evList );
		}else{
			cout << "\nWARNING:::: TEventList <" << ListsName.at(iList) << "> not present in the eventlist map. Skipping it." << endl << endl;
		}
	}
	
	for(unsigned iList=0; iList<listsvector.size(); iList++){
		
		if(iList==0){
			outlist = new TEventList();
			(*outlist) = (*listsvector.at(iList));
			outlist->SetName(outname.c_str());
		}
		outlist->Intersect(listsvector.at(iList));
		
	}
	
	cout << "\nCreated List \""<< outlist->GetName() <<"\". gDirectory " << gDirectory->GetName() << endl;
	
	outlist->SetDirectory(gDirectory);
	
	return outlist;
}


TEventList* IntersectLists(TEventList* ListA, TEventList* ListB, string outname){
	
	TEventList *outlist = NULL;
	
	if(!(ListA&&ListB)) return NULL;
	
	outlist = (TEventList*)ListA->Clone(outname.c_str());
	outlist->Intersect(ListB);
	
	//outlist->SetDirectory(gDirectory);
	
	return outlist;
}


TEventList* SubtractLists(TEventList* ListA, TEventList* ListB, string outname){
	TEventList *outlist = NULL;
	
	if(!(ListA&&ListB)) return NULL;
	
	outlist = new TEventList();
	(*outlist) = (*ListA);
	outlist->SetName(outname.c_str());
	outlist->Subtract(ListB);
	
	outlist->SetDirectory(gDirectory);
	
	return outlist;
}


TEventList* MakeList(TTree* T1, vector<TCut>* cutsVec, string listname){
	
	if(!T1) return NULL;
	
	TEventList *outlist = NULL;
	
	string cmd = string(">>") + listname;
	
	if(!(cutsVec->size()>0)){
		
		outlist = new TEventList(listname.c_str());
		
		T1->Draw(cmd.c_str(), "");
		
		return outlist;
	}
	
	
	TEventList *tmpList = new TEventList("tmpList");
	for(int iCut=0; iCut<cutsVec->size(); iCut++){
		T1->Draw(">>tmpList",cutsVec->at(iCut));
		if(iCut==0){
			outlist = (TEventList*)tmpList->Clone(listname.c_str());
		}else{
			outlist->Intersect(tmpList);
		}
		T1->SetEventList(outlist);
	}
	delete tmpList;
	
	return outlist;
	
}


TEventList* MakeListFromFile(TFile* listfile, vector<string> &ListsNames, string outname){
	
	if(!listfile) return NULL;
	if(!listfile->IsOpen()) return NULL;
	if(ListsNames.size()<=0) return NULL;
	
	TEventList *outlist = NULL;
	
	int counter = 0;
	
	for(unsigned iList=0; iList<ListsNames.size(); iList++){
		TEventList *tmpList = (TEventList*)listfile->Get( ListsNames.at(iList).c_str() );
		if(tmpList){
			if(counter==0){
				outlist = (TEventList*)tmpList->Clone( outname.c_str() );
			}else{
				outlist->Intersect(tmpList);
			}
			counter++;
		}
	}
	return outlist;
}


TH1D* MakeAndDraw1DHisto(TTree* T1, string histname, string title, string variable, int nbins, double xmin, double xmax, TCut _cut, string opt, Color_t color){
	
	TH1D *h = new TH1D(histname.c_str(), title.c_str(), nbins, xmin, xmax);
	h->SetLineColor(color);
	
	string cmd = variable + string(">>") + histname;
	
	T1->Draw(cmd.c_str(), _cut, opt.c_str());
	
	return h;
}


TH1D* MakeAndDraw1DHisto(vector<double> vec, string histname, string title, int nbins, double xmin, double xmax, string _opt, Color_t color){
	
	int nEnt = vec.size();
	if(nEnt<=0) return NULL;
	
	TH1D *h = new TH1D(histname.c_str(), title.c_str(), nbins, xmin, xmax);
	h->SetLineColor(color);
	
	for(int iEnt=0; iEnt<nEnt; iEnt++){
		h->Fill(vec.at(iEnt));
	}
	
	h->Draw(_opt.c_str());
	
	return h;
}


TH2D *MakeAndDraw2DHisto(TTree* T1, string histname, string title, string Xvariable, string Yvariable, int nXbins, double xmin, double xmax, int nYbins, double ymin, double ymax, TCut _cut, string opt){
	
	TH2D *h = new TH2D(histname.c_str(), title.c_str(), nXbins, xmin, xmax, nYbins, ymin, ymax);
	
	string cmd = Yvariable + string(":") + Xvariable + string(">>") + histname;
	
	T1->Draw(cmd.c_str(), _cut, opt.c_str());
	
	return h;
	
}


TH2D* MakeHistoFromTree(TTree *tree, string varX, string varY, double centerX, double widthX, double centerY, double widthY, TCut cut, int nBins){
	
	if(!tree) return NULL;
	
	TH2D* h = new TH2D("h", "", nBins, centerX-widthX, centerX+widthX, nBins, centerY-widthY, centerY+widthY);
	
	stringstream ss_tmp; ss_tmp.str(""); ss_tmp << varY << ":" << varX << ">>h";
	
	tree->Draw(ss_tmp.str().c_str(), cut, "goff");
	
	return h;
}


TGraph* MakeGraphFromTree(TTree* tree, string graphname, string VarX, string VarY, TCut cut){
	
	if(!tree) return NULL;
	
	int nEntries = tree->GetEntries();
	
	tree->SetEstimate(nEntries);
	
	string drawCmd = VarY + string(":") + VarX;
	
	tree->Draw(drawCmd.c_str(), cut, "goff");
	
	nEntries = tree->GetSelectedRows();
	
	double *pY = tree->GetV1();
	double *pX = tree->GetV2();
	
	TGraph *gr = new TGraph(nEntries, pX, pY);
	gr->SetName(graphname.c_str());
	
	return gr;
}


TH2D* DrawGraphAsHisto(TGraph *gr, string histname, string histtitle, int nXbins, double xMin, double xMax, int nYbins, double yMin, double yMax, string opt){
	
	if(!gr) return NULL;
	
	TH2D *h = new TH2D(histname.c_str(), histtitle.c_str(), nXbins, xMin, xMax, nYbins, yMin, yMax);
	
	int nEntries = gr->GetN();
	double *pX = gr->GetX();
	double *pY = gr->GetY();
	
	for(int iEnt=0; iEnt<nEntries; iEnt++){
		h->Fill(pX[iEnt],pY[iEnt]);
	}
	
	h->Draw(opt.c_str());
	
	return h;
}


string formatdigits(double var, double err, int dig){
	
//This function is used to format the output of value with its error with a given number of digits
	
	ostringstream s_var, s_err;
	ostringstream s_var_sc, s_err_sc;
	string outstring;
	
	Int_t log_v, log_e;
	
	log_e = (int)(log10(err));
	if(log10(err)<0.) log_e--;//In this way works fine always
	
	log_v = (int)(log10(var));
	if(log10(var)<0.) log_v--;//In this way works fine always
	
	Int_t nse = (dig-1) - (int)(log_e);
	if(nse < 0) nse = 0; //I don't need digits at the right of the point
	
	
	s_err_sc.precision(dig-1);
	s_err_sc << scientific << err;
	
	err = atof(s_err_sc.str().c_str());
	
	s_err.precision(nse);
	s_err  << fixed << err;
	//cout << "s_err = " << s_err.str() << endl;
	
	
	//Int_t nse = dig - (int)(log10(err)+1);
	
	//cout << "nse = " << nse << endl;
	
	Int_t nse_sc = log_v - log_e + 1;
	s_var_sc.precision(nse_sc);
	//s_var << setw(nse);
	s_var_sc << scientific << var;
	
	
	var = atof(s_var_sc.str().c_str());
	//err = atof(s_err.str().c_str());
	
	
	//s_err << setprecision(dig) << err;
	s_var.precision(nse);
	s_var << fixed << var;
	//cout << "s_var = " << s_var.str() << endl;
	
	
	
	outstring = s_var.str();
	outstring += " +- ";
	outstring += s_err.str();
	
	return outstring;
}

string formatdigits(double var, int dig){//By default the digits used for the precision are 2
	
//This function is used to format the output of a value with a given number of digits (for example to be used in the upper limits)
	
	ostringstream s_var, s_var_sc;
	
	string outstring;
	
	Int_t log_v;
	
	log_v = (int)(log10(var));
	if(log10(var)<0.) log_v--;//In this way works fine always
	
	Int_t nse = (dig-1) - (int)(log_v);
	if(nse < 0) nse = 0; //I don't need digits at the right of the point
	
	s_var_sc.precision(dig - 1);
	s_var_sc << scientific << var;
	
	var = atof(s_var_sc.str().c_str());
	
	s_var.precision(nse);
	s_var << fixed << var;
	//cout << "s_var = " << s_var.str() << endl;
	
	
	outstring = s_var.str();
	
	return outstring;
}


TF1* GaussFit1Dhist(TH1D* hist, string funcName, string opt, double xsigmaLow, double xsigmaUp){
	
	if(!hist){
		cerr << "\nError in Fit1Dhist: histogram pointer not valid!" << endl;
		return NULL;
	}
	
	
	int maxbin = hist->GetMaximumBin();
	double maxval = hist->GetBinContent(maxbin);
	double mean = hist->GetBinCenter(maxbin);
	int lowbin = hist->FindFirstBinAbove( maxval*exp(-0.5) );
	int upbin = hist->FindLastBinAbove( maxval*exp(-0.5) );
	double sigma = ( hist->GetBinCenter(upbin)-hist->GetBinCenter(lowbin) )/2;
	
	TF1 *fitFunc = new TF1(funcName.c_str(), "pow([0],2) + [1]*exp( -pow((x-[2])/[3],2)/2 )", mean-xsigmaLow*sigma, mean+xsigmaUp*sigma);
	fitFunc->SetLineWidth(2);
	fitFunc->SetParNames( "Const","Ampl","Mean","Sigma" );
	fitFunc->SetParameters(0, maxval, mean, sigma);
	
	hist->Fit(fitFunc, opt.c_str());
	
	return fitFunc;
}


//New stuff

TH1D* MakeAmplitudeHisto(string name, string title, TH1D* h, Int_t nBins, Double_t win_scaling, Color_t col)
{
	if(!h) return NULL;
	if(name==string("")) return NULL;
	
	Double_t minX = h->GetMinimum();
	Double_t maxX = h->GetMaximum();
	Double_t delta = maxX - minX;
	
	minX = minX - ((win_scaling-1.)/2)*delta;
	maxX = maxX + ((win_scaling-1.)/2)*delta;
	
	TH1D* hout = new TH1D(name.c_str(), title.c_str(), nBins, minX, maxX);
	for(int iBin=1; iBin<=h->GetNbinsX(); iBin++)
	{
		hout->Fill( h->GetBinContent(iBin) );
	}
	
	hout->SetLineColor(col);
	
	return hout;
}

TH1D* MakeAmplitudeHisto(string name, TH1D* h, Int_t nBins, Double_t win_scaling, Color_t col)
{
	return MakeAmplitudeHisto(name, "", h, nBins, win_scaling, col);
};

TH1D* MakeAmplitudeHisto(string name, string title, Int_t nEnt, Double_t *array, Int_t nBins, Double_t win_scaling, Color_t col)
{
	if(!array) return NULL;
	if(name==string("")) return NULL;
	
	Double_t minX = TMath::MinElement(nEnt, array);
	Double_t maxX = TMath::MaxElement(nEnt, array);
	Double_t delta = maxX - minX;
	
	minX = minX - ((win_scaling-1.)/2)*delta;
	maxX = maxX + ((win_scaling-1.)/2)*delta;
	
	TH1D* hout = new TH1D(name.c_str(), title.c_str(), nBins, minX, maxX);
	for(int iEnt=0; iEnt<nEnt; iEnt++)
	{
		hout->Fill( array[iEnt] );
	}
	
	hout->SetLineColor(col);
	
	return hout;
}

TH1D* MakeAmplitudeHisto(string name, Int_t nEnt, Double_t *array, Int_t nBins, Double_t win_scaling, Color_t col)
{
	return MakeAmplitudeHisto(name, nEnt, array, nBins, win_scaling, col);
}


#endif