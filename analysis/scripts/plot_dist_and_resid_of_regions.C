#include <iostream>
#include <fstream>
#include <string>
//using namespace std;


// ***************************************************
void plot_dist_and_resid_of_regions() {
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
    
    
    // Define number of regions (shown in file 'regions.png')
    int n_regions = 25;


	// Define binning
	int nBins_z   =  38;
    double z_min       =   0.;
    double z_max       = 170.;
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
	TH2D * Dist_vs_z[n_regions];
	TH2D * Residual_x_vs_z[n_regions];
	TH2D * Residual_y_vs_z[n_regions];
	TH2D * Residual_z_vs_z[n_regions];
	
	char * name = new char[50];
	for(int reg=0; reg<n_regions; reg++) {
	    if(reg<10) { sprintf(name,"Dist_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(name,"Dist_vs_z_%d",reg); }
	    Dist_vs_z[reg] = new TH2D(name,"",nBins_z,z_min,z_max,nBins_dist,dist_min,dist_max);
	    if(reg<10) { sprintf(name,"Residual_x_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(name,"Residual_x_vs_z_%d",reg); }
	    Residual_x_vs_z[reg] = new TH2D(name,"",nBins_z,z_min,z_max,nBins_resid,resid_min,resid_max);
	    if(reg<10) { sprintf(name,"Residual_y_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(name,"Residual_y_vs_z_%d",reg); }
	    Residual_y_vs_z[reg] = new TH2D(name,"",nBins_z,z_min,z_max,nBins_resid,resid_min,resid_max);
	    if(reg<10) { sprintf(name,"Residual_z_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(name,"Residual_z_vs_z_%d",reg); }
	    Residual_z_vs_z[reg] = new TH2D(name,"",nBins_z,z_min,z_max,nBins_resid,resid_min,resid_max);
	}



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
            Dist_vs_z[region]->Fill(z+r_z,sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)));
            Residual_x_vs_z[region]->Fill(z+r_z,r_x);
            Residual_y_vs_z[region]->Fill(z+r_z,r_y);
            Residual_z_vs_z[region]->Fill(z+r_z,r_z);
            std::cout << " r_x: " << r_x << " \tr_y: " << r_y << " \tr_z: " << r_z << " \tz: " << z << std::endl;
            if( region!=0 && region!=6 && region!=12 && region!=18 && region!=24 ) { std::cout << " ===================== " << region << std::endl; }
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'Hits_Residuals_Eta_l_region.txt'" << std::endl;



    // Use TProfile to project all distances from one zBin onto the ordinate
    TProfile * profile_Dist_vs_z[n_regions];
    TProfile * profile_residual_x_vs_z[n_regions];
    TProfile * profile_residual_y_vs_z[n_regions];
    TProfile * profile_residual_z_vs_z[n_regions];
    
    for(int reg=0; reg<n_regions; reg++) {
        profile_Dist_vs_z[reg] = Dist_vs_z[reg]->ProfileX();
        // profile_Dist_vs_z[reg]->Fit("pol1");
        profile_residual_x_vs_z[reg] = Residual_x_vs_z[reg]->ProfileX();
        profile_residual_y_vs_z[reg] = Residual_y_vs_z[reg]->ProfileX();
        profile_residual_z_vs_z[reg] = Residual_z_vs_z[reg]->ProfileX();
    }



    // For plots
    // ==========
	gStyle->SetPadRightMargin(0.17);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.21);
    gStyle->SetNumberContours(5); // Number of colors in z-axis
    gStyle->SetStripDecimals(kFALSE); // To force the labels to have the same digits in their decimal numbers
    //TGaxis::SetMaxDigits(0);



    // Plot 2D histograms for the distances and residuals
	// ===================================================

	// DISTANCES
	// +++++++++++++++++++
	for(int reg=0; reg<n_regions; reg++) {
	    char * plot_name = new char[200];
	    if(reg<10) { sprintf(plot_name,"distances_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(plot_name,"distances_vs_z_%d",reg); }
	    char * save_name_png = new char[200];
	    if(reg<10) { sprintf(save_name_png,"plots/region/distances/distances_vs_z_0%d.png",reg); }
	    if(reg>=10) { sprintf(save_name_png,"plots/region/distances/distances_vs_z_%d.png",reg); }
	    char * save_name_pdf = new char[200];
	    if(reg<10) { sprintf(save_name_pdf,"plots/region/distances/distances_vs_z_0%d.pdf",reg); }
	    if(reg>=10) { sprintf(save_name_pdf,"plots/region/distances/distances_vs_z_%d.pdf",reg); }
	    
	    TCanvas * canvas = new TCanvas("distances_vs_z",plot_name);
	    canvas->SetGrid(1);
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        Dist_vs_z[reg]->SetTitle("");
        Dist_vs_z[reg]->GetZaxis()->SetTitleOffset(1.5);
	    Dist_vs_z[reg]->GetZaxis()->SetTitle("entries");
	    //Dist_vs_z[reg]->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
        Dist_vs_z[reg]->GetXaxis()->SetTitleOffset(0.001);
        Dist_vs_z[reg]->GetYaxis()->SetTitleOffset(1.0);
        Dist_vs_z[reg]->GetZaxis()->SetTitleOffset(1.25);
        Dist_vs_z[reg]->GetXaxis()->SetTitleSize(0.035);
        Dist_vs_z[reg]->GetYaxis()->SetTitleSize(0.05);
        Dist_vs_z[reg]->GetZaxis()->SetTitleSize(0.05);
        Dist_vs_z[reg]->GetXaxis()->SetLabelSize(0.05);
        Dist_vs_z[reg]->GetYaxis()->SetLabelSize(0.05);
        Dist_vs_z[reg]->GetZaxis()->SetLabelSize(0.05);
	    Dist_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    Dist_vs_z[reg]->GetYaxis()->SetTitle("distance to PC [mm]");
        Dist_vs_z[reg]->Draw("colz");
        
        gPad->RedrawAxis();
        canvas->Update();
	    
	    canvas->SaveAs(save_name_png);
	    canvas->SaveAs(save_name_pdf);
	}
	
	// RESIDUAL_X
	// +++++++++++++++++++
	for(int reg=0; reg<n_regions; reg++) {
	    char * plot_name = new char[200];
	    if(reg<10) { sprintf(plot_name,"residual_x_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(plot_name,"residual_x_vs_z_%d",reg); }
	    char * save_name_png = new char[200];
	    if(reg<10) { sprintf(save_name_png,"plots/region/residual_x/residual_x_vs_z_0%d.png",reg); }
	    if(reg>=10) { sprintf(save_name_png,"plots/region/residual_x/residual_x_vs_z_%d.png",reg); }
	    char * save_name_pdf = new char[200];
	    if(reg<10) { sprintf(save_name_pdf,"plots/region/residual_x/residual_x_vs_z_0%d.pdf",reg); }
	    if(reg>=10) { sprintf(save_name_pdf,"plots/region/residual_x/residual_x_vs_z_%d.pdf",reg); }
	    
	    TCanvas * canvas = new TCanvas("residual_x_vs_z",plot_name);
	    canvas->SetGrid(1);
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        Residual_x_vs_z[reg]->SetTitle("");
        Residual_x_vs_z[reg]->GetZaxis()->SetTitleOffset(1.5);
	    Residual_x_vs_z[reg]->GetZaxis()->SetTitle("entries");
	    //Residual_x_vs_z[reg]->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
        Residual_x_vs_z[reg]->GetXaxis()->SetTitleOffset(0.001);
        Residual_x_vs_z[reg]->GetYaxis()->SetTitleOffset(1.0);
        Residual_x_vs_z[reg]->GetZaxis()->SetTitleOffset(1.25);
        Residual_x_vs_z[reg]->GetXaxis()->SetTitleSize(0.035);
        Residual_x_vs_z[reg]->GetYaxis()->SetTitleSize(0.05);
        Residual_x_vs_z[reg]->GetZaxis()->SetTitleSize(0.05);
        Residual_x_vs_z[reg]->GetXaxis()->SetLabelSize(0.05);
        Residual_x_vs_z[reg]->GetYaxis()->SetLabelSize(0.05);
        Residual_x_vs_z[reg]->GetZaxis()->SetLabelSize(0.05);
	    Residual_x_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    Residual_x_vs_z[reg]->GetYaxis()->SetTitle("residual x [mm]");
        Residual_x_vs_z[reg]->Draw("colz");
        
        gPad->RedrawAxis();
        canvas->Update();
	    
	    canvas->SaveAs(save_name_png);
	    canvas->SaveAs(save_name_pdf);
	}

	// RESIDUAL_Y
	// +++++++++++++++++++
	for(int reg=0; reg<n_regions; reg++) {
	    char * plot_name = new char[200];
	    if(reg<10) { sprintf(plot_name,"residual_y_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(plot_name,"residual_y_vs_z_%d",reg); }
	    char * save_name_png = new char[200];
	    if(reg<10) { sprintf(save_name_png,"plots/region/residual_y/residual_y_vs_z_0%d.png",reg); }
	    if(reg>=10) { sprintf(save_name_png,"plots/region/residual_y/residual_y_vs_z_%d.png",reg); }
	    char * save_name_pdf = new char[200];
	    if(reg<10) { sprintf(save_name_pdf,"plots/region/residual_y/residual_y_vs_z_0%d.pdf",reg); }
	    if(reg>=10) { sprintf(save_name_pdf,"plots/region/residual_y/residual_y_vs_z_%d.pdf",reg); }
	    
	    TCanvas * canvas = new TCanvas("residual_y_vs_z",plot_name);
	    canvas->SetGrid(1);
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        Residual_y_vs_z[reg]->SetTitle("");
        Residual_y_vs_z[reg]->GetZaxis()->SetTitleOffset(1.5);
	    Residual_y_vs_z[reg]->GetZaxis()->SetTitle("entries");
	    //Residual_y_vs_z[reg]->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
        Residual_y_vs_z[reg]->GetXaxis()->SetTitleOffset(0.001);
        Residual_y_vs_z[reg]->GetYaxis()->SetTitleOffset(1.0);
        Residual_y_vs_z[reg]->GetZaxis()->SetTitleOffset(1.25);
        Residual_y_vs_z[reg]->GetXaxis()->SetTitleSize(0.035);
        Residual_y_vs_z[reg]->GetYaxis()->SetTitleSize(0.05);
        Residual_y_vs_z[reg]->GetZaxis()->SetTitleSize(0.05);
        Residual_y_vs_z[reg]->GetXaxis()->SetLabelSize(0.05);
        Residual_y_vs_z[reg]->GetYaxis()->SetLabelSize(0.05);
        Residual_y_vs_z[reg]->GetZaxis()->SetLabelSize(0.05);
	    Residual_y_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    Residual_y_vs_z[reg]->GetYaxis()->SetTitle("residual y [mm]");
        Residual_y_vs_z[reg]->Draw("colz");
        
        gPad->RedrawAxis();
        canvas->Update();
	    
	    canvas->SaveAs(save_name_png);
	    canvas->SaveAs(save_name_pdf);
	}

	// RESIDUAL_Z
	// +++++++++++++++++++
	for(int reg=0; reg<n_regions; reg++) {
	    char * plot_name = new char[200];
	    if(reg<10) { sprintf(plot_name,"residual_z_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(plot_name,"residual_z_vs_z_%d",reg); }
	    char * save_name_png = new char[200];
	    if(reg<10) { sprintf(save_name_png,"plots/region/residual_z/residual_z_vs_z_0%d.png",reg); }
	    if(reg>=10) { sprintf(save_name_png,"plots/region/residual_z/residual_z_vs_z_%d.png",reg); }
	    char * save_name_pdf = new char[200];
	    if(reg<10) { sprintf(save_name_pdf,"plots/region/residual_z/residual_z_vs_z_0%d.pdf",reg); }
	    if(reg>=10) { sprintf(save_name_pdf,"plots/region/residual_z/residual_z_vs_z_%d.pdf",reg); }
	    
	    TCanvas * canvas = new TCanvas("residual_z_vs_z",plot_name);
	    canvas->SetGrid(1);
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        Residual_z_vs_z[reg]->SetTitle("");
        Residual_z_vs_z[reg]->GetZaxis()->SetTitleOffset(1.5);
	    Residual_z_vs_z[reg]->GetZaxis()->SetTitle("entries");
	    //Residual_z_vs_z[reg]->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
        Residual_z_vs_z[reg]->GetXaxis()->SetTitleOffset(0.001);
        Residual_z_vs_z[reg]->GetYaxis()->SetTitleOffset(1.0);
        Residual_z_vs_z[reg]->GetZaxis()->SetTitleOffset(1.25);
        Residual_z_vs_z[reg]->GetXaxis()->SetTitleSize(0.035);
        Residual_z_vs_z[reg]->GetYaxis()->SetTitleSize(0.05);
        Residual_z_vs_z[reg]->GetZaxis()->SetTitleSize(0.05);
        Residual_z_vs_z[reg]->GetXaxis()->SetLabelSize(0.05);
        Residual_z_vs_z[reg]->GetYaxis()->SetLabelSize(0.05);
        Residual_z_vs_z[reg]->GetZaxis()->SetLabelSize(0.05);
	    Residual_z_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    Residual_z_vs_z[reg]->GetYaxis()->SetTitle("residual z [mm]");
        Residual_z_vs_z[reg]->Draw("colz");
        
        gPad->RedrawAxis();
        canvas->Update();
	    
	    canvas->SaveAs(save_name_png);
	    canvas->SaveAs(save_name_pdf);
	}
	
	

    // Plot 2D and 1D histograms for the distances and residuals
	// ==========================================================
    //gStyle->SetPadRightMargin(0.17);
    //gStyle->SetPadLeftMargin(0.12);
    //gStyle->SetPadTopMargin(0.03);
    //gStyle->SetPadBottomMargin(0.21);



	// DISTANCES WITH WIDTHS
	// +++++++++++++++++++++++++
	for(int reg=0; reg<n_regions; reg++) {
	    char * plot_name = new char[200];
	    if(reg<10) { sprintf(plot_name,"distances_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(plot_name,"distances_vs_z_%d",reg); }
	    char * save_name_png = new char[200];
	    if(reg<10) { sprintf(save_name_png,"plots/region/distances_widths/distance_vs_z_0%d.png",reg); }
	    if(reg>=10) { sprintf(save_name_png,"plots/region/distances_widths/distance_vs_z_%d.png",reg); }
	    char * save_name_pdf = new char[200];
	    if(reg<10) { sprintf(save_name_pdf,"plots/region/distances_widths/distance_vs_z_0%d.pdf",reg); }
	    if(reg>=10) { sprintf(save_name_pdf,"plots/region/distances_widths/distance_vs_z_%d.pdf",reg); }
	    

	    TCanvas * distances_vs_z_with_widths = new TCanvas("distances_vs_z_with_widths",plot_name);
	    distances_vs_z_with_widths->Divide(1,2);
        gStyle->SetOptStat(0);
        distances_vs_z_with_widths->cd(1);
        gPad->SetGrid(1);
        gPad->SetTicks();
        Dist_vs_z[reg]->SetTitle("");
        double range_Dist_vs_z_min = 0.;
        double range_Dist_vs_z_max = 4.;
        Dist_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    Dist_vs_z[reg]->GetYaxis()->SetTitle("distance to PC [mm]");
	    Dist_vs_z[reg]->GetZaxis()->RotateTitle(1);
        Dist_vs_z[reg]->GetZaxis()->SetTitleOffset(1.5);
	    Dist_vs_z[reg]->GetZaxis()->SetTitle("entries");
	    Dist_vs_z[reg]->GetZaxis()->SetRangeUser(4,5);
	    Dist_vs_z[reg]->GetZaxis()->SetRangeUser(0.5,5.5);
	    Dist_vs_z[reg]->GetYaxis()->SetRangeUser(range_Dist_vs_z_min,range_Dist_vs_z_max);
        Dist_vs_z[reg]->GetXaxis()->SetTitleOffset(1.00);
        Dist_vs_z[reg]->GetYaxis()->SetTitleOffset(0.6);
        Dist_vs_z[reg]->GetZaxis()->SetTitleOffset(0.6);
        Dist_vs_z[reg]->GetXaxis()->SetTitleSize(0.11);
        Dist_vs_z[reg]->GetYaxis()->SetTitleSize(0.11);
        Dist_vs_z[reg]->GetZaxis()->SetTitleSize(0.11);
        Dist_vs_z[reg]->GetXaxis()->SetLabelSize(0.1);
        Dist_vs_z[reg]->GetYaxis()->SetLabelSize(0.1);
        Dist_vs_z[reg]->GetZaxis()->SetLabelSize(0.1);
        Dist_vs_z[reg]->Draw("colz");
        //gPad->RedrawAxis();
        distances_vs_z_with_widths->cd(2);
        gPad->SetGrid(1);
        gPad->SetTicks();
        gPad->SetLineWidth(5);
        profile_Dist_vs_z[reg]->SetTitle("");
        profile_Dist_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    profile_Dist_vs_z[reg]->GetYaxis()->SetTitle("distance to PC [mm]");
	    profile_Dist_vs_z[reg]->GetYaxis()->SetRangeUser(range_Dist_vs_z_min,range_Dist_vs_z_max);
        profile_Dist_vs_z[reg]->GetXaxis()->SetTitleOffset(1.00);
        profile_Dist_vs_z[reg]->GetYaxis()->SetTitleOffset(0.6);
        profile_Dist_vs_z[reg]->GetXaxis()->SetTitleSize(0.11);
        profile_Dist_vs_z[reg]->GetYaxis()->SetTitleSize(0.11);
        profile_Dist_vs_z[reg]->GetXaxis()->SetLabelSize(0.1);
        profile_Dist_vs_z[reg]->GetYaxis()->SetLabelSize(0.1);
        profile_Dist_vs_z[reg]->Draw();

        gPad->RedrawAxis();
        distances_vs_z_with_widths->Update();

	    distances_vs_z_with_widths->SaveAs(save_name_png);
	    distances_vs_z_with_widths->SaveAs(save_name_pdf);
	}


	// RESIDUAL_X WITH WIDTHS
	// +++++++++++++++++++++++++
	for(int reg=0; reg<n_regions; reg++) {
	    char * plot_name = new char[200];
	    if(reg<10) { sprintf(plot_name,"residuals_x_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(plot_name,"residuals_x_vs_z_%d",reg); }
	    char * save_name_png = new char[200];
	    if(reg<10) { sprintf(save_name_png,"plots/region/residuals_x_widths/residuals_x_vs_z_0%d.png",reg); }
	    if(reg>=10) { sprintf(save_name_png,"plots/region/residuals_x_widths/residuals_x_vs_z_%d.png",reg); }
	    char * save_name_pdf = new char[200];
	    if(reg<10) { sprintf(save_name_pdf,"plots/region/residuals_x_widths/residuals_x_vs_z_0%d.pdf",reg); }
	    if(reg>=10) { sprintf(save_name_pdf,"plots/region/residuals_x_widths/residuals_x_vs_z_%d.pdf",reg); }
	    

	    TCanvas * residuals_x_vs_z_with_widths = new TCanvas("residuals_x_vs_z_with_widths",plot_name);
	    residuals_x_vs_z_with_widths->SetGrid(1);
	    residuals_x_vs_z_with_widths->Divide(1,2);
        gStyle->SetOptStat(0);
        residuals_x_vs_z_with_widths->cd(1);
        gPad->SetGrid(1);
        gPad->SetTicks();
        double range_residual_x_vs_z_min = -2.;
        double range_residual_x_vs_z_max = 2.;
        Residual_x_vs_z[reg]->SetTitle("");
	    Residual_x_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    Residual_x_vs_z[reg]->GetYaxis()->SetTitle("residual x [mm]");
        Residual_x_vs_z[reg]->GetZaxis()->SetTitleOffset(1.5);
	    Residual_x_vs_z[reg]->GetZaxis()->SetTitle("entries");
	    Residual_x_vs_z[reg]->GetZaxis()->RotateTitle(1);
	    Residual_x_vs_z[reg]->GetZaxis()->SetRangeUser(0.5,5.5);
	    Residual_x_vs_z[reg]->GetYaxis()->SetRangeUser(range_residual_x_vs_z_min,range_residual_x_vs_z_max);
        Residual_x_vs_z[reg]->GetXaxis()->SetTitleOffset(1.00);
        Residual_x_vs_z[reg]->GetYaxis()->SetTitleOffset(0.6);
        Residual_x_vs_z[reg]->GetZaxis()->SetTitleOffset(0.6);
        Residual_x_vs_z[reg]->GetXaxis()->SetTitleSize(0.11);
        Residual_x_vs_z[reg]->GetYaxis()->SetTitleSize(0.11);
        Residual_x_vs_z[reg]->GetZaxis()->SetTitleSize(0.11);
        Residual_x_vs_z[reg]->GetXaxis()->SetLabelSize(0.1);
        Residual_x_vs_z[reg]->GetYaxis()->SetLabelSize(0.1);
        Residual_x_vs_z[reg]->GetZaxis()->SetLabelSize(0.1);
        Residual_x_vs_z[reg]->Draw("colz");
        //gPad->RedrawAxis();
        residuals_x_vs_z_with_widths->cd(2);
        gPad->SetGrid(1);
        gPad->SetTicks();
        profile_residual_x_vs_z[reg]->SetLineWidth(5);
        profile_residual_x_vs_z[reg]->SetTitle("");
        profile_residual_x_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    profile_residual_x_vs_z[reg]->GetYaxis()->SetTitle("residual x [mm]");
	    profile_residual_x_vs_z[reg]->GetYaxis()->SetRangeUser(range_residual_x_vs_z_min,range_residual_x_vs_z_max);
        profile_residual_x_vs_z[reg]->GetXaxis()->SetTitleOffset(1.00);
        profile_residual_x_vs_z[reg]->GetYaxis()->SetTitleOffset(0.6);
        profile_residual_x_vs_z[reg]->GetXaxis()->SetTitleSize(0.11);
        profile_residual_x_vs_z[reg]->GetYaxis()->SetTitleSize(0.11);
        profile_residual_x_vs_z[reg]->GetXaxis()->SetLabelSize(0.1);
        profile_residual_x_vs_z[reg]->GetYaxis()->SetLabelSize(0.1);
        profile_residual_x_vs_z[reg]->Draw();

        gPad->RedrawAxis();
        residuals_x_vs_z_with_widths->Update();

	    residuals_x_vs_z_with_widths->SaveAs(save_name_png);
	    residuals_x_vs_z_with_widths->SaveAs(save_name_pdf);
	}


	// RESIDUAL_Y WITH WIDTHS
	// +++++++++++++++++++++++++
	for(int reg=0; reg<n_regions; reg++) {
	    char * plot_name = new char[200];
	    if(reg<10) { sprintf(plot_name,"residuals_y_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(plot_name,"residuals_y_vs_z_%d",reg); }
	    char * save_name_png = new char[200];
	    if(reg<10) { sprintf(save_name_png,"plots/region/residuals_y_widths/residuals_y_vs_z_0%d.png",reg); }
	    if(reg>=10) { sprintf(save_name_png,"plots/region/residuals_y_widths/residuals_y_vs_z_%d.png",reg); }
	    char * save_name_pdf = new char[200];
	    if(reg<10) { sprintf(save_name_pdf,"plots/region/residuals_y_widths/residuals_y_vs_z_0%d.pdf",reg); }
	    if(reg>=10) { sprintf(save_name_pdf,"plots/region/residuals_y_widths/residuals_y_vs_z_%d.pdf",reg); }
	    

	    TCanvas * residuals_y_vs_z_with_widths = new TCanvas("residuals_y_vs_z_with_widths",plot_name);
	    residuals_y_vs_z_with_widths->SetGrid(1);
	    residuals_y_vs_z_with_widths->Divide(1,2);
        gStyle->SetOptStat(0);
        residuals_y_vs_z_with_widths->cd(1);
        gPad->SetGrid(1);
        gPad->SetTicks();
        double range_residual_y_vs_z_min = -2.;
        double range_residual_y_vs_z_max = 2.;
        Residual_y_vs_z[reg]->SetTitle("");
	    Residual_y_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    Residual_y_vs_z[reg]->GetYaxis()->SetTitle("residual y [mm]");
        Residual_y_vs_z[reg]->GetZaxis()->SetTitleOffset(1.5);
	    Residual_y_vs_z[reg]->GetZaxis()->SetTitle("entries");
	    Residual_y_vs_z[reg]->GetZaxis()->RotateTitle(1);
	    Residual_y_vs_z[reg]->GetZaxis()->SetRangeUser(0.5,5.5);
	    Residual_y_vs_z[reg]->GetYaxis()->SetRangeUser(range_residual_y_vs_z_min,range_residual_y_vs_z_max);
        Residual_y_vs_z[reg]->GetXaxis()->SetTitleOffset(1.00);
        Residual_y_vs_z[reg]->GetYaxis()->SetTitleOffset(0.6);
        Residual_y_vs_z[reg]->GetZaxis()->SetTitleOffset(0.6);
        Residual_y_vs_z[reg]->GetXaxis()->SetTitleSize(0.11);
        Residual_y_vs_z[reg]->GetYaxis()->SetTitleSize(0.11);
        Residual_y_vs_z[reg]->GetZaxis()->SetTitleSize(0.11);
        Residual_y_vs_z[reg]->GetXaxis()->SetLabelSize(0.1);
        Residual_y_vs_z[reg]->GetYaxis()->SetLabelSize(0.1);
        Residual_y_vs_z[reg]->GetZaxis()->SetLabelSize(0.1);
        Residual_y_vs_z[reg]->Draw("colz");
        //gPad->RedrawAxis();
        residuals_y_vs_z_with_widths->cd(2);
        gPad->SetGrid(1);
        gPad->SetTicks();
        profile_residual_y_vs_z[reg]->SetLineWidth(5);
        profile_residual_y_vs_z[reg]->SetTitle("");
        profile_residual_y_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    profile_residual_y_vs_z[reg]->GetYaxis()->SetTitle("residual y [mm]");
	    profile_residual_y_vs_z[reg]->GetYaxis()->SetRangeUser(range_residual_y_vs_z_min,range_residual_y_vs_z_max);
        profile_residual_y_vs_z[reg]->GetXaxis()->SetTitleOffset(1.00);
        profile_residual_y_vs_z[reg]->GetYaxis()->SetTitleOffset(0.6);
        profile_residual_y_vs_z[reg]->GetXaxis()->SetTitleSize(0.11);
        profile_residual_y_vs_z[reg]->GetYaxis()->SetTitleSize(0.11);
        profile_residual_y_vs_z[reg]->GetXaxis()->SetLabelSize(0.1);
        profile_residual_y_vs_z[reg]->GetYaxis()->SetLabelSize(0.1);
        profile_residual_y_vs_z[reg]->Draw();

        gPad->RedrawAxis();
        residuals_y_vs_z_with_widths->Update();

	    residuals_y_vs_z_with_widths->SaveAs(save_name_png);
	    residuals_y_vs_z_with_widths->SaveAs(save_name_pdf);
	}


	// RESIDUAL_Z WITH WIDTHS
	// +++++++++++++++++++++++++
	for(int reg=0; reg<n_regions; reg++) {
	    char * plot_name = new char[200];
	    if(reg<10) { sprintf(plot_name,"residuals_z_vs_z_0%d",reg); }
	    if(reg>=10) { sprintf(plot_name,"residuals_z_vs_z_%d",reg); }
	    char * save_name_png = new char[200];
	    if(reg<10) { sprintf(save_name_png,"plots/region/residuals_z_widths/residuals_z_vs_z_0%d.png",reg); }
	    if(reg>=10) { sprintf(save_name_png,"plots/region/residuals_z_widths/residuals_z_vs_z_%d.png",reg); }
	    char * save_name_pdf = new char[200];
	    if(reg<10) { sprintf(save_name_pdf,"plots/region/residuals_z_widths/residuals_z_vs_z_0%d.pdf",reg); }
	    if(reg>=10) { sprintf(save_name_pdf,"plots/region/residuals_z_widths/residuals_z_vs_z_%d.pdf",reg); }
	    

	    TCanvas * residuals_z_vs_z_with_widths = new TCanvas("residuals_z_vs_z_with_widths",plot_name);
	    residuals_z_vs_z_with_widths->SetGrid(1);
	    residuals_z_vs_z_with_widths->Divide(1,2);
        gStyle->SetOptStat(0);
        residuals_z_vs_z_with_widths->cd(1);
        gPad->SetGrid(1);
        gPad->SetTicks();
	    double range_residual_z_vs_z_min = -2.;
        double range_residual_z_vs_z_max = 2.;
        Residual_z_vs_z[reg]->SetTitle("");
	    Residual_z_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    Residual_z_vs_z[reg]->GetYaxis()->SetTitle("residual z [mm]");
        Residual_z_vs_z[reg]->GetZaxis()->SetTitleOffset(1.5);
	    Residual_z_vs_z[reg]->GetZaxis()->SetTitle("entries");
	    Residual_z_vs_z[reg]->GetZaxis()->RotateTitle(1);
	    Residual_z_vs_z[reg]->GetZaxis()->SetRangeUser(0.5,5.5);
	    Residual_z_vs_z[reg]->GetYaxis()->SetRangeUser(range_residual_z_vs_z_min,range_residual_z_vs_z_max);
        Residual_z_vs_z[reg]->GetXaxis()->SetTitleOffset(1.00);
        Residual_z_vs_z[reg]->GetYaxis()->SetTitleOffset(0.6);
        Residual_z_vs_z[reg]->GetZaxis()->SetTitleOffset(0.6);
        Residual_z_vs_z[reg]->GetXaxis()->SetTitleSize(0.11);
        Residual_z_vs_z[reg]->GetYaxis()->SetTitleSize(0.11);
        Residual_z_vs_z[reg]->GetZaxis()->SetTitleSize(0.11);
        Residual_z_vs_z[reg]->GetXaxis()->SetLabelSize(0.1);
        Residual_z_vs_z[reg]->GetYaxis()->SetLabelSize(0.1);
        Residual_z_vs_z[reg]->GetZaxis()->SetLabelSize(0.1);
        Residual_z_vs_z[reg]->Draw("colz");
        //gPad->RedrawAxis();
        residuals_z_vs_z_with_widths->cd(2);
        gPad->SetGrid(1);
        gPad->SetTicks();
        profile_residual_z_vs_z[reg]->SetLineWidth(5);
        profile_residual_z_vs_z[reg]->SetTitle("");
        profile_residual_z_vs_z[reg]->GetXaxis()->SetTitle("z [mm]");
	    profile_residual_z_vs_z[reg]->GetYaxis()->SetTitle("residual z [mm]");
	    profile_residual_z_vs_z[reg]->GetYaxis()->SetRangeUser(range_residual_z_vs_z_min,range_residual_z_vs_z_max);
        profile_residual_z_vs_z[reg]->GetXaxis()->SetTitleOffset(1.00);
        profile_residual_z_vs_z[reg]->GetYaxis()->SetTitleOffset(0.6);
        profile_residual_z_vs_z[reg]->GetXaxis()->SetTitleSize(0.11);
        profile_residual_z_vs_z[reg]->GetYaxis()->SetTitleSize(0.11);
        profile_residual_z_vs_z[reg]->GetXaxis()->SetLabelSize(0.1);
        profile_residual_z_vs_z[reg]->GetYaxis()->SetLabelSize(0.1);
        profile_residual_z_vs_z[reg]->Draw();

        gPad->RedrawAxis();
        residuals_z_vs_z_with_widths->Update();

	    residuals_z_vs_z_with_widths->SaveAs(save_name_png);
	    residuals_z_vs_z_with_widths->SaveAs(save_name_pdf);
	}

    return;
}
