#include <iostream>
#include <fstream>
#include <string>
//using namespace std;


// ***************************************************
void plot_dist_and_resid_vs_l() {
// ***************************************************

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
	int nBins_l        =  38;
    double l_min       =   0.;
    double l_max       =  38.;
    int nBins_dist     = 100;
    double dist_min    =   0.;
    double dist_max    =   5.5;
    int nBins_resid    = 100;
    double resid_min   =  -5.5;
    double resid_max   =   5.5;
    double resid_min_z =  -2.8;
    double resid_max_z =   2.8;


    // Create histograms
	// ======================
    TH2D * Dist_vs_l       = new TH2D("Dist_vs_l","Dist_vs_l",nBins_l,l_min,l_max,nBins_dist,dist_min,dist_max);
    TH2D * Residual_x_vs_l = new TH2D("Residual_x_vs_l","Residual_x_vs_l",nBins_l,l_min,l_max,nBins_resid,resid_min,resid_max);
    TH2D * Residual_y_vs_l = new TH2D("Residual_y_vs_l","Residual_y_vs_l",nBins_l,l_min,l_max,nBins_resid,resid_min,resid_max);
    TH2D * Residual_z_vs_l = new TH2D("Residual_z_vs_l","Residual_z_vs_l",nBins_l,l_min,l_max,nBins_resid,resid_min_z,resid_max_z);
    


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
            Dist_vs_l->Fill(l*52.5/31.,sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)));
            Residual_x_vs_l->Fill(l*52.5/31.,r_x);
            Residual_y_vs_l->Fill(l*52.5/31.,r_y);
            Residual_z_vs_l->Fill(l*52.5/31.,r_z);
            //std::cout << " r_x: " << r_x << " \tr_y: " << r_y << " \tr_z: " << r_z << std::endl;
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'Hits_Residuals_Eta_l_region.txt'" << std::endl;



    // Use TProfile to project all distances from one lBin onto the ordinate
    TProfile * profile_Dist_vs_l = Dist_vs_l->ProfileX();
    //profile_Dist_vs_l->Fit("pol1");
    TProfile * profile_residual_x_vs_l = Residual_x_vs_l->ProfileX();
    TProfile * profile_residual_y_vs_l = Residual_y_vs_l->ProfileX();
    TProfile * profile_residual_z_vs_l = Residual_z_vs_l->ProfileX();



    // For plots
    // ==========
	gStyle->SetPadRightMargin(0.17);
    //gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetPadBottomMargin(0.14);



    // Plot 2D histograms for the distances and residuals
	// ===================================================
	TCanvas * distances_vs_l = new TCanvas("distances_vs_l","distances_vs_l");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    Dist_vs_l->SetTitle("");
    Dist_vs_l->GetZaxis()->SetTitleOffset(1.3);
	Dist_vs_l->GetZaxis()->SetTitle("entries");
	//Dist_vs_l->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    Dist_vs_l->GetXaxis()->SetTitleOffset(1.3);
    Dist_vs_l->GetYaxis()->SetTitleOffset(1.0);
    Dist_vs_l->GetZaxis()->SetTitleOffset(1.25);
    Dist_vs_l->GetXaxis()->SetTitleSize(0.05);
    Dist_vs_l->GetYaxis()->SetTitleSize(0.05);
    Dist_vs_l->GetZaxis()->SetTitleSize(0.05);
    Dist_vs_l->GetXaxis()->SetLabelSize(0.05);
    Dist_vs_l->GetYaxis()->SetLabelSize(0.05);
    Dist_vs_l->GetZaxis()->SetLabelSize(0.05);
	Dist_vs_l->GetXaxis()->SetTitle("l [mm]");
	Dist_vs_l->GetYaxis()->SetTitle("distance to PC [mm]");
    Dist_vs_l->Draw("colz");
    gPad->RedrawAxis();
	distances_vs_l->SaveAs("plots/l/distance_vs_l.png");
	distances_vs_l->SaveAs("plots/l/distance_vs_l.pdf");
    distances_vs_l->SaveAs("plots/l/distance_vs_l.root");




    TCanvas * residual_x_vs_l = new TCanvas("residual_x_vs_l","residual_x_vs_l");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    Residual_x_vs_l->SetTitle("");
    Residual_x_vs_l->GetZaxis()->SetTitleOffset(1.3);
	Residual_x_vs_l->GetZaxis()->SetTitle("entries");
	//Residual_x_vs_l->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    Residual_x_vs_l->GetXaxis()->SetTitleOffset(1.3);
    Residual_x_vs_l->GetYaxis()->SetTitleOffset(1.0);
    Residual_x_vs_l->GetZaxis()->SetTitleOffset(1.25);
    Residual_x_vs_l->GetXaxis()->SetTitleSize(0.05);
    Residual_x_vs_l->GetYaxis()->SetTitleSize(0.05);
    Residual_x_vs_l->GetZaxis()->SetTitleSize(0.05);
    Residual_x_vs_l->GetXaxis()->SetLabelSize(0.05);
    Residual_x_vs_l->GetYaxis()->SetLabelSize(0.05);
    Residual_x_vs_l->GetZaxis()->SetLabelSize(0.05);
	Residual_x_vs_l->GetXaxis()->SetTitle("l [mm]");
	Residual_x_vs_l->GetYaxis()->SetTitle("residual x to PC [mm]");
    Residual_x_vs_l->Draw("colz");
    gPad->RedrawAxis();
	residual_x_vs_l->SaveAs("plots/l/residual_x_vs_l.png");
	residual_x_vs_l->SaveAs("plots/l/residual_x_vs_l.pdf");
    residual_x_vs_l->SaveAs("plots/l/residual_x_vs_l.root");




    TCanvas * residual_y_vs_l = new TCanvas("residual_y_vs_l","residual_y_vs_l");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    Residual_y_vs_l->SetTitle("");
    Residual_y_vs_l->GetZaxis()->SetTitleOffset(1.3);
	Residual_y_vs_l->GetZaxis()->SetTitle("entries");
	//Residual_y_vs_l->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    Residual_y_vs_l->GetXaxis()->SetTitleOffset(1.3);
    Residual_y_vs_l->GetYaxis()->SetTitleOffset(1.0);
    Residual_y_vs_l->GetZaxis()->SetTitleOffset(1.25);
    Residual_y_vs_l->GetXaxis()->SetTitleSize(0.05);
    Residual_y_vs_l->GetYaxis()->SetTitleSize(0.05);
    Residual_y_vs_l->GetZaxis()->SetTitleSize(0.05);
    Residual_y_vs_l->GetXaxis()->SetLabelSize(0.05);
    Residual_y_vs_l->GetYaxis()->SetLabelSize(0.05);
    Residual_y_vs_l->GetZaxis()->SetLabelSize(0.05);
	Residual_y_vs_l->GetXaxis()->SetTitle("l [mm]");
	Residual_y_vs_l->GetYaxis()->SetTitle("residual y to PC [mm]");
    Residual_y_vs_l->Draw("colz");
    gPad->RedrawAxis();
	residual_y_vs_l->SaveAs("plots/l/residual_y_vs_l.png");
	residual_y_vs_l->SaveAs("plots/l/residual_y_vs_l.pdf");
    residual_y_vs_l->SaveAs("plots/l/residual_y_vs_l.root");




    TCanvas * residual_z_vs_l = new TCanvas("residual_z_vs_l","residual_z_vs_l");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    Residual_z_vs_l->SetTitle("");
    Residual_z_vs_l->GetZaxis()->SetTitleOffset(1.3);
	Residual_z_vs_l->GetZaxis()->SetTitle("entries");
	//Residual_z_vs_l->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    Residual_z_vs_l->GetXaxis()->SetTitleOffset(1.3);
    Residual_z_vs_l->GetYaxis()->SetTitleOffset(1.0);
    Residual_z_vs_l->GetZaxis()->SetTitleOffset(1.25);
    Residual_z_vs_l->GetXaxis()->SetTitleSize(0.05);
    Residual_z_vs_l->GetYaxis()->SetTitleSize(0.05);
    Residual_z_vs_l->GetZaxis()->SetTitleSize(0.05);
    Residual_z_vs_l->GetXaxis()->SetLabelSize(0.05);
    Residual_z_vs_l->GetYaxis()->SetLabelSize(0.05);
    Residual_z_vs_l->GetZaxis()->SetLabelSize(0.05);
	Residual_z_vs_l->GetXaxis()->SetTitle("l [mm]");
	Residual_z_vs_l->GetYaxis()->SetTitle("residual z to PC [mm]");
    Residual_z_vs_l->Draw("colz");
    gPad->RedrawAxis();
	residual_z_vs_l->SaveAs("plots/l/residual_z_vs_l.png");
	residual_z_vs_l->SaveAs("plots/l/residual_z_vs_l.pdf");
    residual_z_vs_l->SaveAs("plots/l/residual_z_vs_l.root");



    // Plot 2D and 1D histograms for the distances and residuals
	// ==========================================================
    gStyle->SetPadRightMargin(0.17);
    //gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.19);



	TCanvas * distances_vs_l_with_widths = new TCanvas("distances_vs_l_with_widths","distances_vs_l_with_widths");
	distances_vs_l_with_widths->Divide(1,2);
    gStyle->SetOptStat(0);
    distances_vs_l_with_widths->cd(1);
    gPad->SetTicks();
    Dist_vs_l->SetTitle("");
	Dist_vs_l->GetXaxis()->SetTitle("l [mm]");
	Dist_vs_l->GetYaxis()->SetTitle("distance to PC [mm]");
    Dist_vs_l->GetZaxis()->SetTitleOffset(1.3);
	Dist_vs_l->GetZaxis()->SetTitle("entries");
	//Dist_vs_l->GetZaxis()->SetRangeUser(-0.5,19);
    Dist_vs_l->GetXaxis()->SetTitleOffset(1.1);
    Dist_vs_l->GetYaxis()->SetTitleOffset(0.7);
    Dist_vs_l->GetZaxis()->SetTitleOffset(0.5);
    Dist_vs_l->GetXaxis()->SetTitleSize(0.08);
    Dist_vs_l->GetYaxis()->SetTitleSize(0.08);
    Dist_vs_l->GetZaxis()->SetTitleSize(0.08);
    Dist_vs_l->GetXaxis()->SetLabelSize(0.08);
    Dist_vs_l->GetYaxis()->SetLabelSize(0.08);
    Dist_vs_l->GetZaxis()->SetLabelSize(0.08);
    Dist_vs_l->Draw("colz");
    TLine * hline_Dist_vs_l = new TLine(15.5*52.5/31,dist_min,15.5*52.5/31,0.995*dist_max); // (xbegin,ybegin,xend,yend)
    hline_Dist_vs_l->SetLineColor(kRed);
    hline_Dist_vs_l->Draw();
    //gPad->RedrawAxis();
    distances_vs_l_with_widths->cd(2);
    gPad->SetTicks();
    profile_Dist_vs_l->SetTitle("");
    profile_Dist_vs_l->GetXaxis()->SetTitle("l [mm]");
	profile_Dist_vs_l->GetYaxis()->SetTitle("distance to PC [mm]");
	double range_Dist_vs_l_min = 0.;
    double range_Dist_vs_l_max = 1.3;
	profile_Dist_vs_l->GetYaxis()->SetRangeUser(range_Dist_vs_l_min,range_Dist_vs_l_max);
    profile_Dist_vs_l->GetXaxis()->SetTitleOffset(1.1);
    profile_Dist_vs_l->GetYaxis()->SetTitleOffset(0.7);
    profile_Dist_vs_l->GetXaxis()->SetTitleSize(0.08);
    profile_Dist_vs_l->GetYaxis()->SetTitleSize(0.08);
    profile_Dist_vs_l->GetXaxis()->SetLabelSize(0.08);
    profile_Dist_vs_l->GetYaxis()->SetLabelSize(0.08);
    profile_Dist_vs_l->Draw();
    TLine * hline_profile_Dist_vs_l = new TLine(15.5*52.5/31,range_Dist_vs_l_min,15.5*52.5/31,range_Dist_vs_l_max); // (xbegin,ybegin,xend,yend)
    hline_profile_Dist_vs_l->SetLineColor(kRed);
    hline_profile_Dist_vs_l->Draw();
    //gPad->RedrawAxis();
	distances_vs_l_with_widths->SaveAs("plots/l/distance_vs_l_with_widths.png");
	distances_vs_l_with_widths->SaveAs("plots/l/distance_vs_l_with_widths.pdf");
    distances_vs_l_with_widths->SaveAs("plots/l/distance_vs_l_with_widths.root");
    
    
    
    
    TCanvas * residual_x_vs_l_with_widths = new TCanvas("residual_x_vs_l_with_widths","residual_x_vs_l_with_widths");
	residual_x_vs_l_with_widths->Divide(1,2);
    gStyle->SetOptStat(0);
    residual_x_vs_l_with_widths->cd(1);
    gPad->SetTicks();
    Residual_x_vs_l->SetTitle("");
	Residual_x_vs_l->GetXaxis()->SetTitle("l [mm]");
	Residual_x_vs_l->GetYaxis()->SetTitle("residual x to PC [mm]");
    Residual_x_vs_l->GetZaxis()->SetTitleOffset(1.3);
	Residual_x_vs_l->GetZaxis()->SetTitle("entries");
	//Residual_x_vs_l->GetZaxis()->SetRangeUser(-0.5,19);
    Residual_x_vs_l->GetXaxis()->SetTitleOffset(1.1);
    Residual_x_vs_l->GetYaxis()->SetTitleOffset(0.7);
    Residual_x_vs_l->GetZaxis()->SetTitleOffset(0.5);
    Residual_x_vs_l->GetXaxis()->SetTitleSize(0.08);
    Residual_x_vs_l->GetYaxis()->SetTitleSize(0.08);
    Residual_x_vs_l->GetZaxis()->SetTitleSize(0.08);
    Residual_x_vs_l->GetXaxis()->SetLabelSize(0.08);
    Residual_x_vs_l->GetYaxis()->SetLabelSize(0.08);
    Residual_x_vs_l->GetZaxis()->SetLabelSize(0.08);
    Residual_x_vs_l->Draw("colz");
    TLine * hline_Residual_x_vs_l = new TLine(15.5*52.5/31,resid_min,15.5*52.5/31,0.995*resid_max); // (xbegin,ybegin,xend,yend)
    hline_Residual_x_vs_l->SetLineColor(kRed);
    hline_Residual_x_vs_l->Draw();
    //gPad->RedrawAxis();
    residual_x_vs_l_with_widths->cd(2);
    gPad->SetTicks();
    profile_residual_x_vs_l->SetTitle("");
    profile_residual_x_vs_l->GetXaxis()->SetTitle("l [mm]");
	profile_residual_x_vs_l->GetYaxis()->SetTitle("residual x to PC [mm]");
	double range_residual_x_vs_l_min = -0.65;
    double range_residual_x_vs_l_max = 0.65;
	profile_residual_x_vs_l->GetYaxis()->SetRangeUser(range_residual_x_vs_l_min,range_residual_x_vs_l_max);
    profile_residual_x_vs_l->GetXaxis()->SetTitleOffset(1.1);
    profile_residual_x_vs_l->GetYaxis()->SetTitleOffset(0.7);
    profile_residual_x_vs_l->GetXaxis()->SetTitleSize(0.08);
    profile_residual_x_vs_l->GetYaxis()->SetTitleSize(0.08);
    profile_residual_x_vs_l->GetXaxis()->SetLabelSize(0.08);
    profile_residual_x_vs_l->GetYaxis()->SetLabelSize(0.08);
    profile_residual_x_vs_l->Draw();
    TLine * hline_profile_residual_x_vs_l = new TLine(15.5*52.5/31,range_residual_x_vs_l_min,15.5*52.5/31,range_residual_x_vs_l_max); // (xbegin,ybegin,xend,yend)
    hline_profile_residual_x_vs_l->SetLineColor(kRed);
    hline_profile_residual_x_vs_l->Draw();
    //gPad->RedrawAxis();
	residual_x_vs_l_with_widths->SaveAs("plots/l/residual_x_vs_l_with_widths.png");
	residual_x_vs_l_with_widths->SaveAs("plots/l/residual_x_vs_l_with_widths.pdf");
    residual_x_vs_l_with_widths->SaveAs("plots/l/residual_x_vs_l_with_widths.root");
    
    
    
    
    TCanvas * residual_y_vs_l_with_widths = new TCanvas("residual_y_vs_l_with_widths","residual_y_vs_l_with_widths");
	residual_y_vs_l_with_widths->Divide(1,2);
    gStyle->SetOptStat(0);
    residual_y_vs_l_with_widths->cd(1);
    gPad->SetTicks();
    Residual_y_vs_l->SetTitle("");
	Residual_y_vs_l->GetXaxis()->SetTitle("l [mm]");
	Residual_y_vs_l->GetYaxis()->SetTitle("residual y to PC [mm]");
    Residual_y_vs_l->GetZaxis()->SetTitleOffset(1.3);
	Residual_y_vs_l->GetZaxis()->SetTitle("entries");
	//Residual_y_vs_l->GetZaxis()->SetRangeUser(-0.5,19);
    Residual_y_vs_l->GetXaxis()->SetTitleOffset(1.1);
    Residual_y_vs_l->GetYaxis()->SetTitleOffset(0.7);
    Residual_y_vs_l->GetZaxis()->SetTitleOffset(0.5);
    Residual_y_vs_l->GetXaxis()->SetTitleSize(0.08);
    Residual_y_vs_l->GetYaxis()->SetTitleSize(0.08);
    Residual_y_vs_l->GetZaxis()->SetTitleSize(0.08);
    Residual_y_vs_l->GetXaxis()->SetLabelSize(0.08);
    Residual_y_vs_l->GetYaxis()->SetLabelSize(0.08);
    Residual_y_vs_l->GetZaxis()->SetLabelSize(0.08);
    Residual_y_vs_l->Draw("colz");
    TLine * hline_Residual_y_vs_l = new TLine(15.5*52.5/31,resid_min,15.5*52.5/31,0.995*resid_max); // (xbegin,ybegin,xend,yend)
    hline_Residual_y_vs_l->SetLineColor(kRed);
    hline_Residual_y_vs_l->Draw();
    //gPad->RedrawAxis();
    residual_y_vs_l_with_widths->cd(2);
    gPad->SetTicks();
    profile_residual_y_vs_l->SetTitle("");
    profile_residual_y_vs_l->GetXaxis()->SetTitle("l [mm]");
	profile_residual_y_vs_l->GetYaxis()->SetTitle("residual y to PC [mm]");
	double range_residual_y_vs_l_min = -0.65;
    double range_residual_y_vs_l_max = 0.65;
	profile_residual_y_vs_l->GetYaxis()->SetRangeUser(range_residual_y_vs_l_min,range_residual_y_vs_l_max);
    profile_residual_y_vs_l->GetXaxis()->SetTitleOffset(1.1);
    profile_residual_y_vs_l->GetYaxis()->SetTitleOffset(0.7);
    profile_residual_y_vs_l->GetXaxis()->SetTitleSize(0.08);
    profile_residual_y_vs_l->GetYaxis()->SetTitleSize(0.08);
    profile_residual_y_vs_l->GetXaxis()->SetLabelSize(0.08);
    profile_residual_y_vs_l->GetYaxis()->SetLabelSize(0.08);
    profile_residual_y_vs_l->Draw();
    TLine * hline_profile_residual_y_vs_l = new TLine(15.5*52.5/31,range_residual_y_vs_l_min,15.5*52.5/31,range_residual_y_vs_l_max); // (xbegin,ybegin,xend,yend)
    hline_profile_residual_y_vs_l->SetLineColor(kRed);
    hline_profile_residual_y_vs_l->Draw();
    //gPad->RedrawAxis();
	residual_y_vs_l_with_widths->SaveAs("plots/l/residual_y_vs_l_with_widths.png");
	residual_y_vs_l_with_widths->SaveAs("plots/l/residual_y_vs_l_with_widths.pdf");
    residual_y_vs_l_with_widths->SaveAs("plots/l/residual_y_vs_l_with_widths.root");




    TCanvas * residual_z_vs_l_with_widths = new TCanvas("residual_z_vs_l_with_widths","residual_z_vs_l_with_widths");
	residual_z_vs_l_with_widths->Divide(1,2);
    gStyle->SetOptStat(0);
    residual_z_vs_l_with_widths->cd(1);
    gPad->SetTicks();
    Residual_z_vs_l->SetTitle("");
	Residual_z_vs_l->GetXaxis()->SetTitle("l [mm]");
	Residual_z_vs_l->GetYaxis()->SetTitle("residual z to PC [mm]");
    Residual_z_vs_l->GetZaxis()->SetTitleOffset(1.3);
	Residual_z_vs_l->GetZaxis()->SetTitle("entries");
	//Residual_z_vs_l->GetZaxis()->SetRangeUser(-0.5,19);
    Residual_z_vs_l->GetXaxis()->SetTitleOffset(1.1);
    Residual_z_vs_l->GetYaxis()->SetTitleOffset(0.7);
    Residual_z_vs_l->GetZaxis()->SetTitleOffset(0.5);
    Residual_z_vs_l->GetXaxis()->SetTitleSize(0.08);
    Residual_z_vs_l->GetYaxis()->SetTitleSize(0.08);
    Residual_z_vs_l->GetZaxis()->SetTitleSize(0.08);
    Residual_z_vs_l->GetXaxis()->SetLabelSize(0.08);
    Residual_z_vs_l->GetYaxis()->SetLabelSize(0.08);
    Residual_z_vs_l->GetZaxis()->SetLabelSize(0.08);
    Residual_z_vs_l->Draw("colz");
    TLine * hline_Residual_z_vs_l = new TLine(15.5*52.5/31,resid_min_z,15.5*52.5/31,0.995*resid_max_z); // (xbegin,ybegin,xend,yend)
    hline_Residual_z_vs_l->SetLineColor(kRed);
    hline_Residual_z_vs_l->Draw();
    //gPad->RedrawAxis();
    residual_z_vs_l_with_widths->cd(2);
    gPad->SetTicks();
    profile_residual_z_vs_l->SetTitle("");
    profile_residual_z_vs_l->GetXaxis()->SetTitle("l [mm]");
	profile_residual_z_vs_l->GetYaxis()->SetTitle("residual z to PC [mm]");
	double range_residual_z_vs_l_min = -0.28;
    double range_residual_z_vs_l_max = 0.28;
	profile_residual_z_vs_l->GetYaxis()->SetRangeUser(range_residual_z_vs_l_min,range_residual_z_vs_l_max);
    profile_residual_z_vs_l->GetXaxis()->SetTitleOffset(1.1);
    profile_residual_z_vs_l->GetYaxis()->SetTitleOffset(0.7);
    profile_residual_z_vs_l->GetXaxis()->SetTitleSize(0.08);
    profile_residual_z_vs_l->GetYaxis()->SetTitleSize(0.08);
    profile_residual_z_vs_l->GetXaxis()->SetLabelSize(0.08);
    profile_residual_z_vs_l->GetYaxis()->SetLabelSize(0.08);
    profile_residual_z_vs_l->Draw();
    TLine * hline_profile_residual_z_vs_l = new TLine(15.5*52.5/31,range_residual_z_vs_l_min,15.5*52.5/31,range_residual_z_vs_l_max); // (xbegin,ybegin,xend,yend)
    hline_profile_residual_z_vs_l->SetLineColor(kRed);
    hline_profile_residual_z_vs_l->Draw();
    //gPad->RedrawAxis();
	residual_z_vs_l_with_widths->SaveAs("plots/l/residual_z_vs_l_with_widths.png");
	residual_z_vs_l_with_widths->SaveAs("plots/l/residual_z_vs_l_with_widths.pdf");
    residual_z_vs_l_with_widths->SaveAs("plots/l/residual_z_vs_l_with_widths.root");
    
    
    
    return;
}
