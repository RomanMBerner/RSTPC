#include "DigitalFilters.hh"



vector<double>* HistoToVector(TH1D* hist){
	
	if(!hist) return NULL;
	
	int nSamps = hist->GetNbinsX();
	
	if(nSamps<=0) return NULL;
	
	vector<double>* result = new vector<double>(nSamps);
	
	for(int iSamp=1; iSamp<=nSamps; iSamp++){
		result->at(iSamp-1) = hist->GetBinContent(iSamp);
	}
	
	return result;
}

TH1D* VectorToHisto(vector<double>* vec, double xMin, double xMax){
	
	if(!vec) return NULL;
	
	int nSamps = vec->size();
	
	if(nSamps<=0) return NULL;
	
	TH1D *hist = new TH1D("hist","", nSamps, xMin, xMax);
	
	for(int iSamp=0; iSamp<nSamps; iSamp++){
		hist->SetBinContent( iSamp+1, vec->at(iSamp) );
	}
	
	return hist;
}


vector<double>* CustomFilter(vector<double> *waveform, TF1* func, int KernelLen){
	
	if(!waveform) return NULL;
	if(!func) return NULL;
	
	int nSamps = waveform->size();
	
	if(nSamps<=0) return NULL;
	
	
	//REMEMBER ALL THE QUANTITIES ARE IN BIN UNITS
	vector<double> *result = new vector<double>(nSamps);
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	//Build the mask
	vector<double> KernelMask(KernelLen);
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	double sumkern = 0;
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) = func->Eval((double)(iSamp-centralPos));
		sumkern += KernelMask.at(iSamp);
	}
	
	//Normalize the kernel
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) /= sumkern;
	}
	
	
	//Make the convolution
	for(int iSamp=0; iSamp<nSamps; iSamp++){
		
		double sum = 0;
		
		for(int kSamp=0; kSamp<KernelLen; kSamp++){
			if( (iSamp - centralPos + kSamp >= 0) && (iSamp - centralPos + kSamp < nSamps) ){
				sum += waveform->at(iSamp - centralPos + kSamp) * KernelMask.at(kSamp);
			}
		}
		
		result->at(iSamp) = sum;
	}
	
	return result;
}

double* CustomFilter(int nSamples, double *waveform, TF1* func, int KernelLen){
	
	if(nSamples<=0) return NULL;
	if(!waveform) return NULL;
	if(!func) return NULL;
	
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	double *filter = new double[nSamples];
	
	
	//Build the mask
	vector<double> KernelMask(KernelLen);
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	double sumkern = 0;
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) = func->Eval((double)(iSamp-centralPos));
		sumkern += KernelMask.at(iSamp);
	}
	
	//Normalize the kernel
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) /= sumkern;
	}
	
	
	//Make the convolution
	for(int iSamp=0; iSamp<nSamples; iSamp++){
		
		double sum = 0;
		
		for(int kSamp=0; kSamp<KernelLen; kSamp++){
			if( (iSamp - centralPos + kSamp >= 0) && (iSamp - centralPos + kSamp < nSamples) ){
				sum += waveform[iSamp - centralPos + kSamp] * KernelMask.at(kSamp);
			}
		}
		
		filter[iSamp] = sum;
	}
	
	return filter;
	
}

void CustomFilter(int nSamples, double *waveform, double *filter, TF1* func, int KernelLen){
	
	if(nSamples<=0) return;
	if(!waveform) return;
	if(!filter) return;
	if(!func) return;
	
	//REMEMBER ALL THE QUANTITIES ARE IN BIN UNITS
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	//Build the mask
	vector<double> KernelMask(KernelLen);
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	double sumkern = 0;
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) = func->Eval((double)(iSamp-centralPos));
		sumkern += KernelMask.at(iSamp);
	}
	
	//Normalize the kernel
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) /= sumkern;
	}
	
	
	//Make the convolution
	for(int iSamp=0; iSamp<nSamples; iSamp++){
		
		double sum = 0;
		
		for(int kSamp=0; kSamp<KernelLen; kSamp++){
			if( (iSamp - centralPos + kSamp >= 0) && (iSamp - centralPos + kSamp < nSamples) ){
				sum += waveform[iSamp - centralPos + kSamp] * KernelMask.at(kSamp);
			}
		}
		
		filter[iSamp] = sum;
	}
	
	return;
}



vector<double>* GaussianFilter(vector<double> *waveform, int KernelLen, double sigma){
	
	TF1 *gausFunc = new TF1("gausFunc","exp(-pow(x/[0],2)/2)",-KernelLen,KernelLen);
	gausFunc->SetParameter(0,sigma);
	
	vector<double>* filter = CustomFilter(waveform, gausFunc, KernelLen);
	
	if(gausFunc) delete gausFunc;
	
	return filter;
}

double* GaussianFilter(int nSamples, double *waveform, int KernelLen, double sigma){
	
	if(nSamples<=0) return NULL;
	if(!waveform) return NULL;
	
	
	TF1 *gausFunc = new TF1("gausFunc","exp(-pow(x/[0],2)/2)",-KernelLen,KernelLen);
	gausFunc->SetParameter(0,sigma);
	
	double* filter = CustomFilter(nSamples, waveform, gausFunc, KernelLen);
	
	if(gausFunc) delete gausFunc;
	
	return filter;
}

void GaussianFilter(int nSamples, double *waveform, double *filter, int KernelLen, double sigma){
	
	if(nSamples<=0) return;
	if(!waveform) return;
	if(!filter) return;
	
	TF1 *gausFunc = new TF1("gausFunc","exp(-pow(x/[0],2)/2)",-KernelLen,KernelLen);
	gausFunc->SetParameter(0,sigma);
	
	CustomFilter(nSamples, waveform, filter, gausFunc, KernelLen);
	
	if(gausFunc) delete gausFunc;
	
	return;
}




vector<double>* CustomDerivative(vector<double> *waveform, TF1* func, int KernelLen){
	
	//Here the function must already be the derivative!
	
	if(!waveform) return NULL;
	if(!func) return NULL;
	
	int nSamps = waveform->size();
	
	if(nSamps<=0) return NULL;
	
	
	//REMEMBER ALL THE QUANTITIES ARE IN BIN UNITS
	vector<double> *result = new vector<double>(nSamps);
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	//Build the mask
	vector<double> KernelMask(KernelLen);
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	double kerMean = 0;
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) = func->Eval((double)(iSamp-centralPos));
		kerMean += KernelMask.at(iSamp)/KernelLen;
	}
	
	//Shift the kernel to a zero mean kernel (for the derivative)
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) -= kerMean;
	}
	
	
	//Make the convolution
	for(int iSamp=0; iSamp<nSamps; iSamp++){
		
		double sum = 0;
		
		for(int kSamp=0; kSamp<KernelLen; kSamp++){
			if( (iSamp - centralPos + kSamp >= 0) && (iSamp - centralPos + kSamp < nSamps) ){
				sum += waveform->at(iSamp - centralPos + kSamp) * KernelMask.at(kSamp);
			}
		}
		
		result->at(iSamp) = sum;
	}
	
	return result;
}

double* CustomDerivative(int nSamples, double *waveform, TF1* func, int KernelLen){
	
	//Here the function must already be the derivative!
	
	if(nSamples<=0) return NULL;
	if(!waveform) return NULL;
	if(!func) return NULL;
	
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	double *filter = new double[nSamples];
	
	
	//Build the mask
	vector<double> KernelMask(KernelLen);
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	double kerMean = 0;
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) = func->Eval((double)(iSamp-centralPos));
		kerMean += KernelMask.at(iSamp)/KernelLen;
	}
	
	//Normalize the kernel
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) -= kerMean;
	}
	
	
	//Make the convolution
	for(int iSamp=0; iSamp<nSamples; iSamp++){
		
		double sum = 0;
		
		for(int kSamp=0; kSamp<KernelLen; kSamp++){
			if( (iSamp - centralPos + kSamp >= 0) && (iSamp - centralPos + kSamp < nSamples) ){
				sum += waveform[iSamp - centralPos + kSamp] * KernelMask.at(kSamp);
			}
		}
		
		filter[iSamp] = sum;
	}
	
	return filter;
	
}

void CustomDerivative(int nSamples, double *waveform, double *filter, TF1* func, int KernelLen){
	
	//Here the function must already be the derivative!
	
	if(nSamples<=0) return;
	if(!waveform) return;
	if(!filter) return;
	if(!func) return;
	
	//REMEMBER ALL THE QUANTITIES ARE IN BIN UNITS
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	//Build the mask
	vector<double> KernelMask(KernelLen);
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	double kerMean = 0;
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) = func->Eval((double)(iSamp-centralPos));
		kerMean += KernelMask.at(iSamp)/KernelLen;
	}
	
	//Normalize the kernel
	for(int iSamp=0; iSamp<KernelLen; iSamp++){
		KernelMask.at(iSamp) -= kerMean;
	}
	
	
	//Make the convolution
	for(int iSamp=0; iSamp<nSamples; iSamp++){
		
		double sum = 0;
		
		for(int kSamp=0; kSamp<KernelLen; kSamp++){
			if( (iSamp - centralPos + kSamp >= 0) && (iSamp - centralPos + kSamp < nSamples) ){
				sum += waveform[iSamp - centralPos + kSamp] * KernelMask.at(kSamp);
			}
		}
		
		filter[iSamp] = sum;
	}
	
	return;
}



vector<double>* GaussianDerivative(vector<double> *waveform, int KernelLen, double sigma){
	
	if(!waveform) return NULL;
	
	vector<double>* gaussfilt = GaussianFilter(waveform, KernelLen, sigma);
	
	if(!gaussfilt) return NULL;
	
	int nSamp = gaussfilt->size();
	
	vector<double>* filter = new vector<double>(nSamp);
	for(int iSamp=0; iSamp<nSamp; iSamp++){
		
		if(iSamp>0 && iSamp<nSamp-1){
			filter->at(iSamp) = gaussfilt->at(iSamp+1)-gaussfilt->at(iSamp-1);
		}else{
			filter->at(iSamp) = 0;
		}
		
	}
	
	if(gaussfilt) delete gaussfilt;
	
	return filter;
}


double* GaussianDerivative(int nSamples, double *waveform, int KernelLen, double sigma){
	
	if(nSamples<=0) return NULL;
	if(!waveform) return NULL;
	
	double* gaussfilt = GaussianFilter(nSamples, waveform, KernelLen, sigma);
	
	if(!gaussfilt) return NULL;
	
	double* filter = new double[nSamples];
	for(int iSamp=0; iSamp<nSamples; iSamp++){
		
		if(iSamp>0 && iSamp<nSamples-1){
			filter[iSamp] = gaussfilt[iSamp+1] - gaussfilt[iSamp-1];
		}else{
			filter[iSamp] = 0;
		}
	}
	
	
	if(gaussfilt) delete [] gaussfilt;
	
	return filter;
}


void GaussianDerivative(int nSamples, double *waveform, double *filter, int KernelLen, double sigma){
	
	if(nSamples<=0) return;
	if(!waveform) return;
	if(!filter) return;
	
	
	double* gaussfilt = GaussianFilter(nSamples, waveform, KernelLen, sigma);
	
	if(!gaussfilt) return;
	
	for(int iSamp=0; iSamp<nSamples; iSamp++){
		
		if(iSamp>0 && iSamp<nSamples-1){
			filter[iSamp] = gaussfilt[iSamp+1] - gaussfilt[iSamp-1];
		}else{
			filter[iSamp] = 0;
		}
	}
	
	if(gaussfilt) delete [] gaussfilt;
	
	return;
}

/*
vector<double>* GaussianDerivative(vector<double> *waveform, int KernelLen, double sigma){
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	
	//As a first thing i need a normalized gaussian kernel
	double N =0;
	for(int i=0; i<KernelLen; i++){
		N += exp( -pow((i-centralPos)/sigma,2)/2 );
	}
	
	TF1 *gausDerivFunc = new TF1("gausDerivFunc","-[0]*x*exp(-pow(x/[1],2)/2)/pow([1],2)",-KernelLen,KernelLen);
	gausDerivFunc->SetParameters(1./N,sigma);
	
	vector<double>* filter = CustomDerivative(waveform, gausDerivFunc, KernelLen);
	
	if(gausDerivFunc) delete gausDerivFunc;
	
	return filter;
	
}

double* GaussianDerivative(int nSamples, double *waveform, int KernelLen, double sigma){
	
	if(nSamples<=0) return NULL;
	if(!waveform) return NULL;
	
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	
	//As a first thing i need a normalized gaussian kernel
	double N =0;
	for(int i=0; i<KernelLen; i++){
		N += exp( -pow((i-centralPos)/sigma,2)/2 );
	}
	
	
	TF1 *gausDerivFunc = new TF1("gausDerivFunc","-[0]*x*exp(-pow(x/[1],2)/2)/pow([1],2)",-KernelLen,KernelLen);
	gausDerivFunc->SetParameters(1./N,sigma);
	
	double* filter = CustomDerivative(nSamples, waveform, gausDerivFunc, KernelLen);
	
	if(gausDerivFunc) delete gausDerivFunc;
	
	return filter;
	
}

void GaussianDerivative(int nSamples, double *waveform, double *filter, int KernelLen, double sigma){
	
	if(nSamples<=0) return;
	if(!waveform) return;
	if(!filter) return;
	
	
	//Make the kernel always with odd lenght
	if(KernelLen % 2 == 0) KernelLen += 1;
	if(KernelLen<3) KernelLen=3;//Minimum required lenght
	
	int centralPos = (KernelLen-1)/2; //This is always exact since the kernal lenght is always odd
	
	//As a first thing i need a normalized gaussian kernel
	double N =0;
	for(int i=0; i<KernelLen; i++){
		N += exp( -pow((i-centralPos)/sigma,2)/2 );
	}
	
	
	TF1 *gausDerivFunc = new TF1("gausDerivFunc","-[0]*x*exp(-pow(x/[1],2)/2)/pow([1],2)",-KernelLen,KernelLen);
	gausDerivFunc->SetParameters(1./N,sigma);
	
	
	CustomFilter(nSamples, waveform, filter, gausDerivFunc, KernelLen);
	
	if(gausDerivFunc) delete gausDerivFunc;
	
	return;
}
*/