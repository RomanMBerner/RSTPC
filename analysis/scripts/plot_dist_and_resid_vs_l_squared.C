#include <iostream>
#include <fstream>
#include <string>
//using namespace std;


// ***************************************************
void plot_dist_and_resid_vs_l_squared() {
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
    double l_max       =  1300; // draw as function of l*l such that the geometrical effect cancels (in [mm]) // 1300 ~ nBins_l*nBins_l
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
    TH2D * Dist_vs_l_squared       = new TH2D("Dist_vs_l_squared","Dist_vs_l_squared",nBins_l,l_min,l_max,nBins_dist,dist_min,dist_max);
    TH2D * Residual_x_vs_l_squared = new TH2D("Residual_x_vs_l_squared","Residual_x_vs_l_squared",nBins_l,l_min,l_max,nBins_resid,resid_min,resid_max);
    TH2D * Residual_y_vs_l_squared = new TH2D("Residual_y_vs_l_squared","Residual_y_vs_l_squared",nBins_l,l_min,l_max,nBins_resid,resid_min,resid_max);
    TH2D * Residual_z_vs_l_squared = new TH2D("Residual_z_vs_l_squared","Residual_z_vs_l_squared",nBins_l,l_min,l_max,nBins_resid,resid_min_z,resid_max_z);
    


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
            Dist_vs_l_squared->Fill(pow(l*52.5/31.,2),sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)));
            Residual_x_vs_l_squared->Fill(pow(l*52.5/31.,2),r_x);
            Residual_y_vs_l_squared->Fill(pow(l*52.5/31.,2),r_y);
            Residual_z_vs_l_squared->Fill(pow(l*52.5/31.,2),r_z);
            //std::cout << " r_x: " << r_x << " \tr_y: " << r_y << " \tr_z: " << r_z << std::endl;
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'Hits_Residuals_Eta_l_region.txt'" << std::endl;



    // Use TProfile to project all distances from one lBin onto the ordinate
    TProfile * profile_Dist_vs_l_squared = Dist_vs_l_squared->ProfileX();
    //profile_Dist_vs_l_squared->Fit("pol1");
    TProfile * profile_residual_x_vs_l_squared = Residual_x_vs_l_squared->ProfileX();
    TProfile * profile_residual_y_vs_l_squared = Residual_y_vs_l_squared->ProfileX();
    TProfile * profile_residual_z_vs_l_squared = Residual_z_vs_l_squared->ProfileX();



    // For plots
    // ==========
	gStyle->SetPadRightMargin(0.17);
    //gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetPadBottomMargin(0.14);



    // Plot 2D histograms for the distances and residuals
	// ===================================================
	TCanvas * distances_vs_l_squared = new TCanvas("distances_vs_l_squared","distances_vs_l_squared");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    Dist_vs_l_squared->SetTitle("");
    Dist_vs_l_squared->GetZaxis()->SetTitleOffset(1.3);
	Dist_vs_l_squared->GetZaxis()->SetTitle("entries");
	//Dist_vs_l_squared->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    Dist_vs_l_squared->GetXaxis()->SetTitleOffset(1.3);
    Dist_vs_l_squared->GetYaxis()->SetTitleOffset(1.0);
    Dist_vs_l_squared->GetZaxis()->SetTitleOffset(1.25);
    Dist_vs_l_squared->GetXaxis()->SetTitleSize(0.05);
    Dist_vs_l_squared->GetYaxis()->SetTitleSize(0.05);
    Dist_vs_l_squared->GetZaxis()->SetTitleSize(0.05);
    Dist_vs_l_squared->GetXaxis()->SetLabelSize(0.05);
    Dist_vs_l_squared->GetYaxis()->SetLabelSize(0.05);
    Dist_vs_l_squared->GetZaxis()->SetLabelSize(0.05);
	Dist_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	Dist_vs_l_squared->GetYaxis()->SetTitle("distance to PC [mm]");
    Dist_vs_l_squared->Draw("colz");
    gPad->RedrawAxis();
	distances_vs_l_squared->SaveAs("plots/l_squared/distance_vs_l_squared.png");
	distances_vs_l_squared->SaveAs("plots/l_squared/distance_vs_l_squared.pdf");
    distances_vs_l_squared->SaveAs("plots/l_squared/distance_vs_l_squared.root");




    TCanvas * residual_x_vs_l_squared = new TCanvas("residual_x_vs_l_squared","residual_x_vs_l_squared");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    Residual_x_vs_l_squared->SetTitle("");
    Residual_x_vs_l_squared->GetZaxis()->SetTitleOffset(1.3);
	Residual_x_vs_l_squared->GetZaxis()->SetTitle("entries");
	//Residual_x_vs_l_squared->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    Residual_x_vs_l_squared->GetXaxis()->SetTitleOffset(1.3);
    Residual_x_vs_l_squared->GetYaxis()->SetTitleOffset(1.0);
    Residual_x_vs_l_squared->GetZaxis()->SetTitleOffset(1.25);
    Residual_x_vs_l_squared->GetXaxis()->SetTitleSize(0.05);
    Residual_x_vs_l_squared->GetYaxis()->SetTitleSize(0.05);
    Residual_x_vs_l_squared->GetZaxis()->SetTitleSize(0.05);
    Residual_x_vs_l_squared->GetXaxis()->SetLabelSize(0.05);
    Residual_x_vs_l_squared->GetYaxis()->SetLabelSize(0.05);
    Residual_x_vs_l_squared->GetZaxis()->SetLabelSize(0.05);
	Residual_x_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	Residual_x_vs_l_squared->GetYaxis()->SetTitle("residual x to PC [mm]");
    Residual_x_vs_l_squared->Draw("colz");
    gPad->RedrawAxis();
	residual_x_vs_l_squared->SaveAs("plots/l_squared/residual_x_vs_l_squared.png");
	residual_x_vs_l_squared->SaveAs("plots/l_squared/residual_x_vs_l_squared.pdf");
    residual_x_vs_l_squared->SaveAs("plots/l_squared/residual_x_vs_l_squared.root");




    TCanvas * residual_y_vs_l_squared = new TCanvas("residual_y_vs_l_squared","residual_y_vs_l_squared");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    Residual_y_vs_l_squared->SetTitle("");
    Residual_y_vs_l_squared->GetZaxis()->SetTitleOffset(1.3);
	Residual_y_vs_l_squared->GetZaxis()->SetTitle("entries");
	//Residual_y_vs_l_squared->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    Residual_y_vs_l_squared->GetXaxis()->SetTitleOffset(1.3);
    Residual_y_vs_l_squared->GetYaxis()->SetTitleOffset(1.0);
    Residual_y_vs_l_squared->GetZaxis()->SetTitleOffset(1.25);
    Residual_y_vs_l_squared->GetXaxis()->SetTitleSize(0.05);
    Residual_y_vs_l_squared->GetYaxis()->SetTitleSize(0.05);
    Residual_y_vs_l_squared->GetZaxis()->SetTitleSize(0.05);
    Residual_y_vs_l_squared->GetXaxis()->SetLabelSize(0.05);
    Residual_y_vs_l_squared->GetYaxis()->SetLabelSize(0.05);
    Residual_y_vs_l_squared->GetZaxis()->SetLabelSize(0.05);
	Residual_y_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	Residual_y_vs_l_squared->GetYaxis()->SetTitle("residual y to PC [mm]");
    Residual_y_vs_l_squared->Draw("colz");
    gPad->RedrawAxis();
	residual_y_vs_l_squared->SaveAs("plots/l_squared/residual_y_vs_l_squared.png");
	residual_y_vs_l_squared->SaveAs("plots/l_squared/residual_y_vs_l_squared.pdf");
    residual_y_vs_l_squared->SaveAs("plots/l_squared/residual_y_vs_l_squared.root");




    TCanvas * residual_z_vs_l_squared = new TCanvas("residual_z_vs_l_squared","residual_z_vs_l_squared");
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    Residual_z_vs_l_squared->SetTitle("");
    Residual_z_vs_l_squared->GetZaxis()->SetTitleOffset(1.3);
	Residual_z_vs_l_squared->GetZaxis()->SetTitle("entries");
	//Residual_z_vs_l_squared->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
    Residual_z_vs_l_squared->GetXaxis()->SetTitleOffset(1.3);
    Residual_z_vs_l_squared->GetYaxis()->SetTitleOffset(1.0);
    Residual_z_vs_l_squared->GetZaxis()->SetTitleOffset(1.25);
    Residual_z_vs_l_squared->GetXaxis()->SetTitleSize(0.05);
    Residual_z_vs_l_squared->GetYaxis()->SetTitleSize(0.05);
    Residual_z_vs_l_squared->GetZaxis()->SetTitleSize(0.05);
    Residual_z_vs_l_squared->GetXaxis()->SetLabelSize(0.05);
    Residual_z_vs_l_squared->GetYaxis()->SetLabelSize(0.05);
    Residual_z_vs_l_squared->GetZaxis()->SetLabelSize(0.05);
	Residual_z_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	Residual_z_vs_l_squared->GetYaxis()->SetTitle("residual z to PC [mm]");
    Residual_z_vs_l_squared->Draw("colz");
    gPad->RedrawAxis();
	residual_z_vs_l_squared->SaveAs("plots/l_squared/residual_z_vs_l_squared.png");
	residual_z_vs_l_squared->SaveAs("plots/l_squared/residual_z_vs_l_squared.pdf");
    residual_z_vs_l_squared->SaveAs("plots/l_squared/residual_z_vs_l_squared.root");



    // Plot 2D and 1D histograms for the distances and residuals
	// ==========================================================
    gStyle->SetPadRightMargin(0.17);
    //gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.19);



	TCanvas * distances_vs_l_squared_with_widths = new TCanvas("distances_vs_l_squared_with_widths","distances_vs_l_squared_with_widths");
	distances_vs_l_squared_with_widths->Divide(1,2);
    gStyle->SetOptStat(0);
    distances_vs_l_squared_with_widths->cd(1);
    gPad->SetTicks();
    Dist_vs_l_squared->SetTitle("");
	Dist_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	Dist_vs_l_squared->GetYaxis()->SetTitle("distance to PC [mm]");
    Dist_vs_l_squared->GetZaxis()->SetTitleOffset(1.3);
	Dist_vs_l_squared->GetZaxis()->SetTitle("entries");
	//Dist_vs_l_squared->GetZaxis()->SetRangeUser(-0.5,19);
    Dist_vs_l_squared->GetXaxis()->SetTitleOffset(1.1);
    Dist_vs_l_squared->GetYaxis()->SetTitleOffset(0.7);
    Dist_vs_l_squared->GetZaxis()->SetTitleOffset(0.5);
    Dist_vs_l_squared->GetXaxis()->SetTitleSize(0.08);
    Dist_vs_l_squared->GetYaxis()->SetTitleSize(0.08);
    Dist_vs_l_squared->GetZaxis()->SetTitleSize(0.08);
    Dist_vs_l_squared->GetXaxis()->SetLabelSize(0.08);
    Dist_vs_l_squared->GetYaxis()->SetLabelSize(0.08);
    Dist_vs_l_squared->GetZaxis()->SetLabelSize(0.08);
    Dist_vs_l_squared->Draw("colz");
    TLine * hline_Dist_vs_l_squared = new TLine(pow(15.5*52.5/31,2),dist_min,pow(15.5*52.5/31,2),0.995*dist_max); // (xbegin,ybegin,xend,yend)
    hline_Dist_vs_l_squared->SetLineColor(kRed);
    hline_Dist_vs_l_squared->Draw();
    //gPad->RedrawAxis();
    distances_vs_l_squared_with_widths->cd(2);
    gPad->SetTicks();
    profile_Dist_vs_l_squared->SetTitle("");
    profile_Dist_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	profile_Dist_vs_l_squared->GetYaxis()->SetTitle("distance to PC [mm]");
	double range_Dist_vs_l_squared_min = 0.;
    double range_Dist_vs_l_squared_max = 1.7;
	profile_Dist_vs_l_squared->GetYaxis()->SetRangeUser(range_Dist_vs_l_squared_min,range_Dist_vs_l_squared_max);
    profile_Dist_vs_l_squared->GetXaxis()->SetTitleOffset(1.1);
    profile_Dist_vs_l_squared->GetYaxis()->SetTitleOffset(0.7);
    profile_Dist_vs_l_squared->GetXaxis()->SetTitleSize(0.08);
    profile_Dist_vs_l_squared->GetYaxis()->SetTitleSize(0.08);
    profile_Dist_vs_l_squared->GetXaxis()->SetLabelSize(0.08);
    profile_Dist_vs_l_squared->GetYaxis()->SetLabelSize(0.08);
    profile_Dist_vs_l_squared->Draw();
    TLine * hline_profile_Dist_vs_l_squared = new TLine(pow(15.5*52.5/31,2),range_Dist_vs_l_squared_min,pow(15.5*52.5/31,2),range_Dist_vs_l_squared_max); // (xbegin,ybegin,xend,yend)
    hline_profile_Dist_vs_l_squared->SetLineColor(kRed);
    hline_profile_Dist_vs_l_squared->Draw();
    //gPad->RedrawAxis();
	distances_vs_l_squared_with_widths->SaveAs("plots/l_squared/distance_vs_l_squared_with_widths.png");
	distances_vs_l_squared_with_widths->SaveAs("plots/l_squared/distance_vs_l_squared_with_widths.pdf");
    distances_vs_l_squared_with_widths->SaveAs("plots/l_squared/distance_vs_l_squared_with_widths.root");
    
    
    
    
    TCanvas * residual_x_vs_l_squared_with_widths = new TCanvas("residual_x_vs_l_squared_with_widths","residual_x_vs_l_squared_with_widths");
	residual_x_vs_l_squared_with_widths->Divide(1,2);
    gStyle->SetOptStat(0);
    residual_x_vs_l_squared_with_widths->cd(1);
    gPad->SetTicks();
    Residual_x_vs_l_squared->SetTitle("");
	Residual_x_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	Residual_x_vs_l_squared->GetYaxis()->SetTitle("residual x to PC [mm]");
    Residual_x_vs_l_squared->GetZaxis()->SetTitleOffset(1.3);
	Residual_x_vs_l_squared->GetZaxis()->SetTitle("entries");
	//Residual_x_vs_l_squared->GetZaxis()->SetRangeUser(-0.5,19);
    Residual_x_vs_l_squared->GetXaxis()->SetTitleOffset(1.1);
    Residual_x_vs_l_squared->GetYaxis()->SetTitleOffset(0.7);
    Residual_x_vs_l_squared->GetZaxis()->SetTitleOffset(0.5);
    Residual_x_vs_l_squared->GetXaxis()->SetTitleSize(0.08);
    Residual_x_vs_l_squared->GetYaxis()->SetTitleSize(0.08);
    Residual_x_vs_l_squared->GetZaxis()->SetTitleSize(0.08);
    Residual_x_vs_l_squared->GetXaxis()->SetLabelSize(0.08);
    Residual_x_vs_l_squared->GetYaxis()->SetLabelSize(0.08);
    Residual_x_vs_l_squared->GetZaxis()->SetLabelSize(0.08);
    Residual_x_vs_l_squared->Draw("colz");
    TLine * hline_Residual_x_vs_l_squared = new TLine(pow(15.5*52.5/31,2),resid_min,pow(15.5*52.5/31,2),0.995*resid_max); // (xbegin,ybegin,xend,yend)
    hline_Residual_x_vs_l_squared->SetLineColor(kRed);
    hline_Residual_x_vs_l_squared->Draw();
    //gPad->RedrawAxis();
    residual_x_vs_l_squared_with_widths->cd(2);
    gPad->SetTicks();
    profile_residual_x_vs_l_squared->SetTitle("");
    profile_residual_x_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	profile_residual_x_vs_l_squared->GetYaxis()->SetTitle("residual x to PC [mm]");
	double range_residual_x_vs_l_squared_min = -0.85;
    double range_residual_x_vs_l_squared_max = 0.85;
	profile_residual_x_vs_l_squared->GetYaxis()->SetRangeUser(range_residual_x_vs_l_squared_min,range_residual_x_vs_l_squared_max);
    profile_residual_x_vs_l_squared->GetXaxis()->SetTitleOffset(1.1);
    profile_residual_x_vs_l_squared->GetYaxis()->SetTitleOffset(0.7);
    profile_residual_x_vs_l_squared->GetXaxis()->SetTitleSize(0.08);
    profile_residual_x_vs_l_squared->GetYaxis()->SetTitleSize(0.08);
    profile_residual_x_vs_l_squared->GetXaxis()->SetLabelSize(0.08);
    profile_residual_x_vs_l_squared->GetYaxis()->SetLabelSize(0.08);
    profile_residual_x_vs_l_squared->Draw();
    TLine * hline_profile_residual_x_vs_l_squared = new TLine(pow(15.5*52.5/31,2),range_residual_x_vs_l_squared_min,pow(15.5*52.5/31,2),range_residual_x_vs_l_squared_max); // (xbegin,ybegin,xend,yend)
    hline_profile_residual_x_vs_l_squared->SetLineColor(kRed);
    hline_profile_residual_x_vs_l_squared->Draw();
    //gPad->RedrawAxis();
	residual_x_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_x_vs_l_squared_with_widths.png");
	residual_x_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_x_vs_l_squared_with_widths.pdf");
    residual_x_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_x_vs_l_squared_with_widths.root");
    
    
    
    
    TCanvas * residual_y_vs_l_squared_with_widths = new TCanvas("residual_y_vs_l_squared_with_widths","residual_y_vs_l_squared_with_widths");
	residual_y_vs_l_squared_with_widths->Divide(1,2);
    gStyle->SetOptStat(0);
    residual_y_vs_l_squared_with_widths->cd(1);
    gPad->SetTicks();
    Residual_y_vs_l_squared->SetTitle("");
	Residual_y_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	Residual_y_vs_l_squared->GetYaxis()->SetTitle("residual y to PC [mm]");
    Residual_y_vs_l_squared->GetZaxis()->SetTitleOffset(1.3);
	Residual_y_vs_l_squared->GetZaxis()->SetTitle("entries");
	//Residual_y_vs_l_squared->GetZaxis()->SetRangeUser(-0.5,19);
    Residual_y_vs_l_squared->GetXaxis()->SetTitleOffset(1.1);
    Residual_y_vs_l_squared->GetYaxis()->SetTitleOffset(0.7);
    Residual_y_vs_l_squared->GetZaxis()->SetTitleOffset(0.5);
    Residual_y_vs_l_squared->GetXaxis()->SetTitleSize(0.08);
    Residual_y_vs_l_squared->GetYaxis()->SetTitleSize(0.08);
    Residual_y_vs_l_squared->GetZaxis()->SetTitleSize(0.08);
    Residual_y_vs_l_squared->GetXaxis()->SetLabelSize(0.08);
    Residual_y_vs_l_squared->GetYaxis()->SetLabelSize(0.08);
    Residual_y_vs_l_squared->GetZaxis()->SetLabelSize(0.08);
    Residual_y_vs_l_squared->Draw("colz");
    TLine * hline_Residual_y_vs_l_squared = new TLine(pow(15.5*52.5/31,2),resid_min,pow(15.5*52.5/31,2),0.995*resid_max); // (xbegin,ybegin,xend,yend)
    hline_Residual_y_vs_l_squared->SetLineColor(kRed);
    hline_Residual_y_vs_l_squared->Draw();
    //gPad->RedrawAxis();
    residual_y_vs_l_squared_with_widths->cd(2);
    gPad->SetTicks();
    profile_residual_y_vs_l_squared->SetTitle("");
    profile_residual_y_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	profile_residual_y_vs_l_squared->GetYaxis()->SetTitle("residual y to PC [mm]");
	double range_residual_y_vs_l_squared_min = -0.85;
    double range_residual_y_vs_l_squared_max = 0.85;
	profile_residual_y_vs_l_squared->GetYaxis()->SetRangeUser(range_residual_y_vs_l_squared_min,range_residual_y_vs_l_squared_max);
    profile_residual_y_vs_l_squared->GetXaxis()->SetTitleOffset(1.1);
    profile_residual_y_vs_l_squared->GetYaxis()->SetTitleOffset(0.7);
    profile_residual_y_vs_l_squared->GetXaxis()->SetTitleSize(0.08);
    profile_residual_y_vs_l_squared->GetYaxis()->SetTitleSize(0.08);
    profile_residual_y_vs_l_squared->GetXaxis()->SetLabelSize(0.08);
    profile_residual_y_vs_l_squared->GetYaxis()->SetLabelSize(0.08);
    profile_residual_y_vs_l_squared->Draw();
    TLine * hline_profile_residual_y_vs_l_squared = new TLine(pow(15.5*52.5/31,2),range_residual_y_vs_l_squared_min,pow(15.5*52.5/31,2),range_residual_y_vs_l_squared_max); // (xbegin,ybegin,xend,yend)
    hline_profile_residual_y_vs_l_squared->SetLineColor(kRed);
    hline_profile_residual_y_vs_l_squared->Draw();
    //gPad->RedrawAxis();
	residual_y_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_y_vs_l_squared_with_widths.png");
	residual_y_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_y_vs_l_squared_with_widths.pdf");
    residual_y_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_y_vs_l_squared_with_widths.root");




    TCanvas * residual_z_vs_l_squared_with_widths = new TCanvas("residual_z_vs_l_squared_with_widths","residual_z_vs_l_squared_with_widths");
	residual_z_vs_l_squared_with_widths->Divide(1,2);
    gStyle->SetOptStat(0);
    residual_z_vs_l_squared_with_widths->cd(1);
    gPad->SetTicks();
    Residual_z_vs_l_squared->SetTitle("");
	Residual_z_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	Residual_z_vs_l_squared->GetYaxis()->SetTitle("residual z to PC [mm]");
    Residual_z_vs_l_squared->GetZaxis()->SetTitleOffset(1.3);
	Residual_z_vs_l_squared->GetZaxis()->SetTitle("entries");
	//Residual_z_vs_l_squared->GetZaxis()->SetRangeUser(-0.5,19);
    Residual_z_vs_l_squared->GetXaxis()->SetTitleOffset(1.1);
    Residual_z_vs_l_squared->GetYaxis()->SetTitleOffset(0.7);
    Residual_z_vs_l_squared->GetZaxis()->SetTitleOffset(0.5);
    Residual_z_vs_l_squared->GetXaxis()->SetTitleSize(0.08);
    Residual_z_vs_l_squared->GetYaxis()->SetTitleSize(0.08);
    Residual_z_vs_l_squared->GetZaxis()->SetTitleSize(0.08);
    Residual_z_vs_l_squared->GetXaxis()->SetLabelSize(0.08);
    Residual_z_vs_l_squared->GetYaxis()->SetLabelSize(0.08);
    Residual_z_vs_l_squared->GetZaxis()->SetLabelSize(0.08);
    Residual_z_vs_l_squared->Draw("colz");
    TLine * hline_Residual_z_vs_l_squared = new TLine(pow(15.5*52.5/31,2),resid_min_z,pow(15.5*52.5/31,2),0.995*resid_max_z); // (xbegin,ybegin,xend,yend)
    hline_Residual_z_vs_l_squared->SetLineColor(kRed);
    hline_Residual_z_vs_l_squared->Draw();
    //gPad->RedrawAxis();
    residual_z_vs_l_squared_with_widths->cd(2);
    gPad->SetTicks();
    profile_residual_z_vs_l_squared->SetTitle("");
    profile_residual_z_vs_l_squared->GetXaxis()->SetTitle("l^{2} [mm^{2}]");
	profile_residual_z_vs_l_squared->GetYaxis()->SetTitle("residual z to PC [mm]");
	double range_residual_z_vs_l_squared_min = -0.28;
    double range_residual_z_vs_l_squared_max = 0.28;
	profile_residual_z_vs_l_squared->GetYaxis()->SetRangeUser(range_residual_z_vs_l_squared_min,range_residual_z_vs_l_squared_max);
    profile_residual_z_vs_l_squared->GetXaxis()->SetTitleOffset(1.1);
    profile_residual_z_vs_l_squared->GetYaxis()->SetTitleOffset(0.7);
    profile_residual_z_vs_l_squared->GetXaxis()->SetTitleSize(0.08);
    profile_residual_z_vs_l_squared->GetYaxis()->SetTitleSize(0.08);
    profile_residual_z_vs_l_squared->GetXaxis()->SetLabelSize(0.08);
    profile_residual_z_vs_l_squared->GetYaxis()->SetLabelSize(0.08);
    profile_residual_z_vs_l_squared->Draw();
    TLine * hline_profile_residual_z_vs_l_squared = new TLine(pow(15.5*52.5/31,2),range_residual_z_vs_l_squared_min,pow(15.5*52.5/31,2),range_residual_z_vs_l_squared_max); // (xbegin,ybegin,xend,yend)
    hline_profile_residual_z_vs_l_squared->SetLineColor(kRed);
    hline_profile_residual_z_vs_l_squared->Draw();
    //gPad->RedrawAxis();
	residual_z_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_z_vs_l_squared_with_widths.png");
	residual_z_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_z_vs_l_squared_with_widths.pdf");
    residual_z_vs_l_squared_with_widths->SaveAs("plots/l_squared/residual_z_vs_l_squared_with_widths.root");
    
    
    
    return;
}
