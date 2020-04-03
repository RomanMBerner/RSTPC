#include <iostream>
#include <fstream>
#include <string>
//using namespace std;


// ***************************************************
void plot_MSC_angles() {
// ***************************************************

    gStyle->SetPadRightMargin(0.02);
    //gStyle->SetPadLeftMargin(0.02);
    gStyle->SetPadTopMargin(0.02);
    //gStyle->SetPadBottomMargin(0.02);

	// Batch mode for plots (disable single canvases with canvas->SetBatch(kFALSE);
	//gROOT->SetBatch(kTRUE);
  

    // Get number of lines/hits in the file 'Angles_and_lengths.txt'
    ULong_t lines = 0;
    std::string line;
    ifstream inputfile;
    inputfile.open("Angles_and_lengths.txt");
    if(inputfile.is_open()) {
        while(getline(inputfile,line)) {
          lines++;
        }
        inputfile.close();
    }
    std::cout << " Number of lines: " << lines << std::endl;


    // Define histogram
    TH2F * angle_vs_length = new TH2F("angle_vs_length","angle_vs_length",60,0,12,90,0,180);


    // Fill histograms with hits and distances
    double angle, length;
    inputfile.open("Angles_and_lengths.txt");
    if(inputfile.is_open()) {
        while( std::getline(inputfile,line) ) {
            istringstream strStream(line);
            strStream >> angle >> length;
            angle_vs_length->Fill(sqrt(length),angle);
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'Angles_and_lengths.txt'" << std::endl;


	double x_axis_min = 0.;
	double x_axis_max = 30.;
	double y_axis_min = 0.;
	double y_axis_max = 180.;


    // Plot 1D histogram
	// ====================================

    TCanvas * canv_00 = new TCanvas("canv_00","canv_00");
    gStyle->SetOptStat(0);
    angle_vs_length->SetTitle("");
	angle_vs_length->GetXaxis()->SetTitle("sqrt(length [mm])");
	angle_vs_length->GetYaxis()->SetTitle("#theta [deg]"); // (fMeanTime/20 * drift_vel)
	angle_vs_length->GetXaxis()->SetTitleOffset(1.4);
	angle_vs_length->GetYaxis()->SetTitleOffset(1.4);
	angle_vs_length->SetMarkerStyle(2);
	angle_vs_length->Draw("SCAT");

    gPad->RedrawAxis();
	//canv_00->SaveAs("plots/distances/ColWires_integrated_distance_x.png");
	//canv_00->SaveAs("plots/distances/ColWires_integrated_distance_x.pdf");
    //canv_00->SaveAs("plots/distances/ColWires_integrated_distance_x.root");


    return;
}
