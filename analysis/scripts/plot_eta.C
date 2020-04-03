#include <iostream>
#include <fstream>
#include <string>
//using namespace std;


// ***************************************************
void plot_eta() {
// ***************************************************

    gStyle->SetPadRightMargin(0.17);
    //gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetPadBottomMargin(0.14);

	// Batch mode for plots (disable single canvases with canvas->SetBatch(kFALSE);
	gROOT->SetBatch(kTRUE);
  


    // Get number of lines/hits in the file 'Hits_Residuals_Eta_l_region.txt'
    ULong_t lines = 0;
    std::string line;
    ifstream inputfile;
    inputfile.open("Hits_Residuals_Eta_l_region.txt");
    if(inputfile.is_open()) {
        while(getline(inputfile,line)) {
          lines++;
        }
        inputfile.close();
    }
    std::cout << " Number of lines: " << lines << std::endl;


	// Define binning
	int nBins_eta   =  10; // defined via average_hits_per_good_event
    double eta_min  = 0.-1./(2*nBins_eta);
    double eta_max  = 1.+1./(2*nBins_eta);
    int nBins_resid = 100;
    int resid_min   = -5;
    int resid_max   =  5;
    double resid_min_z =  -2.8;
    double resid_max_z =   2.8;
    int nBins_dist  = 100;
    int dist_min    =   0;
    int dist_max    =  19;


    // Create histograms
	// ======================
    TH2D * hResidual_x_to_PC = new TH2D("hResidual_x_to_PC","hResidual_x_to_PC",nBins_eta,eta_min,eta_max,nBins_resid,resid_min,resid_max);
    TH2D * hResidual_y_to_PC = new TH2D("hResidual_y_to_PC","hResidual_y_to_PC",nBins_eta,eta_min,eta_max,nBins_resid,resid_min,resid_max);
    TH2D * hResidual_z_to_PC = new TH2D("hResidual_z_to_PC","hResidual_z_to_PC",nBins_eta,eta_min,eta_max,nBins_resid,resid_min_z,resid_max_z);
    TH2D * hDist_to_PC       = new TH2D("hDist_to_PC","hDist_to_PC",nBins_eta,eta_min,eta_max,nBins_dist,dist_min,dist_max);


    // Fill histograms with distances and residuals
    // =============================================
    double x, y, z;
    double r_x, r_y, r_z;
    double eta, l;
    int region;
    inputfile.open("Hits_Residuals_Eta_l_region.txt");
    if(inputfile.is_open()) {
        while( std::getline(inputfile,line) ) {
            istringstream strStream(line);
            strStream >> x >> y >> z >> r_x >> r_y >> r_z >> eta >> l >> region;
            hResidual_x_to_PC->Fill(eta,r_x);
            hResidual_y_to_PC->Fill(eta,r_y);
            hResidual_z_to_PC->Fill(eta,r_z);
            hDist_to_PC->Fill(eta,sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)));
            //std::cout << " r_x: " << r_x << " \tr_y: " << r_y << " \tr_z: " << r_z << std::endl;
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'Hits_Residuals_Eta_l_region.txt'" << std::endl;



    // Plot 2D histograms for the distances and residuals
	// ===================================================
	TCanvas * distances_vs_eta = new TCanvas("distances_vs_eta","distances_vs_eta");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    hDist_to_PC->SetTitle("");
    hDist_to_PC->GetZaxis()->SetTitleOffset(1.3);
	hDist_to_PC->GetZaxis()->SetTitle("entries");
	//hDist_to_PC->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    hDist_to_PC->GetXaxis()->SetTitleOffset(1.3);
    hDist_to_PC->GetYaxis()->SetTitleOffset(1.0);
    hDist_to_PC->GetZaxis()->SetTitleOffset(1.25);
    hDist_to_PC->GetXaxis()->SetTitleSize(0.05);
    hDist_to_PC->GetYaxis()->SetTitleSize(0.05);
    hDist_to_PC->GetZaxis()->SetTitleSize(0.05);
    hDist_to_PC->GetXaxis()->SetLabelSize(0.05);
    hDist_to_PC->GetYaxis()->SetLabelSize(0.05);
    hDist_to_PC->GetZaxis()->SetLabelSize(0.05);
	hDist_to_PC->GetXaxis()->SetTitle("#eta [-]");
	hDist_to_PC->GetYaxis()->SetTitle("distance to PC [mm]");
    hDist_to_PC->Draw("colz");
    gPad->RedrawAxis();
	distances_vs_eta->SaveAs("plots/eta/distance_to_PC.png");
	distances_vs_eta->SaveAs("plots/eta/distance_to_PC.pdf");
    distances_vs_eta->SaveAs("plots/eta/distance_to_PC.root");


    TCanvas * residual_x_vs_eta = new TCanvas("residual_x_vs_eta","residual_x_vs_eta");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    hResidual_x_to_PC->SetTitle("");
    hResidual_x_to_PC->GetZaxis()->SetTitleOffset(1.3);
	hResidual_x_to_PC->GetZaxis()->SetTitle("entries");
	//hResidual_x_to_PC->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    hResidual_x_to_PC->GetXaxis()->SetTitleOffset(1.3);
    hResidual_x_to_PC->GetYaxis()->SetTitleOffset(1.0);
    hResidual_x_to_PC->GetZaxis()->SetTitleOffset(1.25);
    hResidual_x_to_PC->GetXaxis()->SetTitleSize(0.05);
    hResidual_x_to_PC->GetYaxis()->SetTitleSize(0.05);
    hResidual_x_to_PC->GetZaxis()->SetTitleSize(0.05);
    hResidual_x_to_PC->GetXaxis()->SetLabelSize(0.05);
    hResidual_x_to_PC->GetYaxis()->SetLabelSize(0.05);
    hResidual_x_to_PC->GetZaxis()->SetLabelSize(0.05);
	hResidual_x_to_PC->GetXaxis()->SetTitle("#eta [-]");
	hResidual_x_to_PC->GetYaxis()->SetTitle("residual x to PC [mm]");
    hResidual_x_to_PC->Draw("colz");
    gPad->RedrawAxis();
	residual_x_vs_eta->SaveAs("plots/eta/residual_x_to_PC.png");
	residual_x_vs_eta->SaveAs("plots/eta/residual_x_to_PC.pdf");
    residual_x_vs_eta->SaveAs("plots/eta/residual_x_to_PC.root");


    TCanvas * residual_y_vs_eta = new TCanvas("residual_y_vs_eta","residual_y_vs_eta");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    hResidual_y_to_PC->SetTitle("");
    hResidual_y_to_PC->GetZaxis()->SetTitleOffset(1.3);
	hResidual_y_to_PC->GetZaxis()->SetTitle("entries");
	//hResidual_y_to_PC->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    hResidual_y_to_PC->GetXaxis()->SetTitleOffset(1.3);
    hResidual_y_to_PC->GetYaxis()->SetTitleOffset(1.0);
    hResidual_y_to_PC->GetZaxis()->SetTitleOffset(1.25);
    hResidual_y_to_PC->GetXaxis()->SetTitleSize(0.05);
    hResidual_y_to_PC->GetYaxis()->SetTitleSize(0.05);
    hResidual_y_to_PC->GetZaxis()->SetTitleSize(0.05);
    hResidual_y_to_PC->GetXaxis()->SetLabelSize(0.05);
    hResidual_y_to_PC->GetYaxis()->SetLabelSize(0.05);
    hResidual_y_to_PC->GetZaxis()->SetLabelSize(0.05);
	hResidual_y_to_PC->GetXaxis()->SetTitle("#eta [-]");
	hResidual_y_to_PC->GetYaxis()->SetTitle("residual y to PC [mm]");
    hResidual_y_to_PC->Draw("colz");
    gPad->RedrawAxis();
	residual_y_vs_eta->SaveAs("plots/eta/residual_y_to_PC.png");
	residual_y_vs_eta->SaveAs("plots/eta/residual_y_to_PC.pdf");
    residual_y_vs_eta->SaveAs("plots/eta/residual_y_to_PC.root");


    TCanvas * residual_z_vs_eta = new TCanvas("residual_z_vs_eta","residual_z_vs_eta");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    hResidual_z_to_PC->SetTitle("");
    hResidual_z_to_PC->GetZaxis()->SetTitleOffset(1.3);
	hResidual_z_to_PC->GetZaxis()->SetTitle("entries");
	//hResidual_z_to_PC->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    hResidual_z_to_PC->GetXaxis()->SetTitleOffset(1.3);
    hResidual_z_to_PC->GetYaxis()->SetTitleOffset(1.0);
    hResidual_z_to_PC->GetZaxis()->SetTitleOffset(1.25);
    hResidual_z_to_PC->GetXaxis()->SetTitleSize(0.05);
    hResidual_z_to_PC->GetYaxis()->SetTitleSize(0.05);
    hResidual_z_to_PC->GetZaxis()->SetTitleSize(0.05);
    hResidual_z_to_PC->GetXaxis()->SetLabelSize(0.05);
    hResidual_z_to_PC->GetYaxis()->SetLabelSize(0.05);
    hResidual_z_to_PC->GetZaxis()->SetLabelSize(0.05);
	hResidual_z_to_PC->GetXaxis()->SetTitle("#eta [-]");
	hResidual_z_to_PC->GetYaxis()->SetTitle("residual z to PC [mm]");
    hResidual_z_to_PC->Draw("colz");
    gPad->RedrawAxis();
	residual_z_vs_eta->SaveAs("plots/eta/residual_z_to_PC.png");
	residual_z_vs_eta->SaveAs("plots/eta/residual_z_to_PC.pdf");
    residual_z_vs_eta->SaveAs("plots/eta/residual_z_to_PC.root");


    return;
}
