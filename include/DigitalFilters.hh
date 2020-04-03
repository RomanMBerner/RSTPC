#ifndef DIGITAL_FILTERS_HH
#define DIGITAL_FILTERS_HH

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


vector<double>* HistoToVector(TH1D* hist);

TH1D* VectorToHisto(vector<double>* vec, double xMin, double xMax);

//Smoothers
vector<double>* CustomFilter(vector<double> *waveform, TF1* func, int KernelLen);

double* CustomFilter(int nSamples, double *waveform, TF1* func, int KernelLen);

void CustomFilter(int nSamples, double *waveform, double *filter, TF1* func, int KernelLen);


vector<double>* GaussianFilter(vector<double> *waveform, int KernelLen, double sigma);

double* GaussianFilter(int nSamples, double *waveform, int KernelLen, double sigma);

void GaussianFilter(int nSamples, double *waveform, double *filter, int KernelLen, double sigma);


//Smoothing derivatives
vector<double>* CustomDerivative(vector<double> *waveform, TF1* func, int KernelLen);

double* CustomDerivative(int nSamples, double *waveform, TF1* func, int KernelLen);

void CustomDerivative(int nSamples, double *waveform, double *filter, TF1* func, int KernelLen);


vector<double>* GaussianDerivative(vector<double> *waveform, int KernelLen, double sigma);

double* GaussianDerivative(int nSamples, double *waveform, int KernelLen, double sigma);

void GaussianDerivative(int nSamples, double *waveform, double *filter, int KernelLen, double sigma);


#endif