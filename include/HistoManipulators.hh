#ifndef HISTO_FILTERS_HH
#define HISTO_FILTERS_HH

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TMath.h"
#include "TRandom.h"
#include "TChain.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "THStack.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCut.h"
#include "TPaletteAxis.h"
#include "TList.h"
#include "TColor.h"
#include "TParameter.h"
#include "TEventList.h"
#include "TEllipse.h"
#include "TProof.h" 

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>


class EdgeResult{
public:
	EdgeResult(){
		lowEdge = -1;
		midEdge = -1;
		upEdge = -1;
	};
	
	double lowEdge, midEdge, upEdge;
};


vector<double> GaussianFilter(vector<double> histo, int KernelLen, double sigma);

vector<double> DerivativeFilter(vector<double> histo);

EdgeResult FindMaxInRange(TH1D *histo, double xmin, double xmax, double thrs=0.9);

EdgeResult FindMinInRange(TH1D *histo, double xmin, double xmax, double thrs=0.9);

void ScaleHistoByAbsMax(TH1D* histo);



//This stuff was in "Functions.hh"
vector<double> GetValuesFromTree(TTree* T1, string variable, TCut cut);

void DetermineProbInt(vector<double> &values, double prob, double &low, double &up, bool sorted);

void DetermineProbInt(vector<double> &values, double prob, double &low, double &up);

void DetermineProbInt(TTree* T1, string variable, TCut cut, double prob, double &low, double &up);

TH1D* MakeAndDraw1DHisto(TTree* T1, string histname, string title, string variable, int nbins, double xmin, double xmax, TCut _cut="", string opt="", Color_t color=kBlack);

TH1D* MakeAndDraw1DHisto(vector<double> vec, string histname, string title, int nbins, double xmin, double xmax, string _opt="", Color_t color=kBlack);

TH2D* MakeAndDraw2DHisto(TTree* T1, string histname, string title, string Xvariable, string Yvariable, int nXbins, double xmin, double xmax, int nYbins, double ymin, double ymax, TCut _cut=TCut(""), string opt="colz");

TEventList* MakeListFromFile(TFile* listfile, vector<string> &ListsNames, string outname="");

TEventList* IntersectLists(vector<string> ListsName, map<string, TEventList*> *listsMap, string outname="");

TEventList* IntersectLists(TEventList* ListA, TEventList* ListB, string outname="");

TEventList* SubtractLists(TEventList* ListA, TEventList* ListB, string outname="");

TEventList* MakeList(TTree* T1, vector<TCut>* cutsVec, string listname);

TH2D* MakeHistoFromTree(TTree *tree, string varX, string varY, double centerX, double widthX, double centerY, double widthY, TCut cut="", int nBins=200);

TGraph* MakeGraphFromTree(TTree* tree, string graphname, string VarX, string VarY, TCut cut="");

TH2D* DrawGraphAsHisto(TGraph *gr, string histname, string histtitle, int nXbins, double xMin, double xMax, int nYbins, double yMin, double yMax, string opt="colz");

string formatdigits(double var, double err, int dig=1);

string formatdigits(double var, int dig=2);

TF1* GaussFit1Dhist(TH1D* hist, string funcName, string opt="LL M R N", double xsigmaLow=1.5, double xsigmaUp=1.5);



//Additional stuff
TH1D* MakeAmplitudeHisto(string name, string title, TH1D* h, Int_t nBins=100, Double_t win_scaling=1., Color_t col=kBlack);
TH1D* MakeAmplitudeHisto(string name, TH1D* h, Int_t nBins=100, Double_t win_scaling=1., Color_t col=kBlack);

TH1D* MakeAmplitudeHisto(string name, string title, Int_t nEnt, Double_t *Array, Int_t nBins=100, Double_t win_scaling=1., Color_t col=kBlue);
TH1D* MakeAmplitudeHisto(string name, Int_t nEnt, Double_t *Array, Int_t nBins=100, Double_t win_scaling=1., Color_t col=kBlue);

#endif