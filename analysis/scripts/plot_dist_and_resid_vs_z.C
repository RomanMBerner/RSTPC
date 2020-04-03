#include <iostream>
#include <fstream>
#include <string>
//using namespace std;


// ***************************************************
void plot_dist_and_resid_vs_z() {
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
	int nBins_z        =   8;
    double z_min       =   0.;
    double z_max       = 170.;
    int nBins_x        = 4;
    double x_min       = -0.5;
    double x_max       = 31.5;
    int nBins_y        = 4;
    double y_min       = -0.5;
    double y_max       = 31.5;
    int nBins_dist     = 100;
    double dist_min    =   0.;
    double dist_max    =  19.;
    int nBins_resid    = 100;
    double resid_min   =  -5.5;
    double resid_max   =   5.5;
    double resid_min_z =  -2.8;
    double resid_max_z =   2.8;


    // Create histograms
	// ======================
    TH2F * Dist_vs_xy_plane[nBins_z];
    TH2F * Residual_x_vs_xy_plane[nBins_z];
    TH2F * Residual_y_vs_xy_plane[nBins_z];
    TH2F * Residual_z_vs_xy_plane[nBins_z];
    char * Dist_name = new char[100];
    char * Res_x_name = new char[100];
    char * Res_y_name = new char[100];
    char * Res_z_name = new char[100];
    for(int z_slice=0; z_slice<nBins_z; z_slice++) {
        sprintf(Dist_name,"Dist_vs_xy_plane_z_slice_%d",z_slice);
        Dist_vs_xy_plane[z_slice] = new TH2F(Dist_name,"",nBins_x,x_min,x_max,nBins_y,y_min,y_max);
        sprintf(Res_x_name,"Residual_x_vs_xy_plane_z_slice_%d",z_slice);
        Residual_x_vs_xy_plane[z_slice] = new TH2F(Res_x_name,"",nBins_x,x_min,x_max,nBins_y,y_min,y_max);
        sprintf(Res_y_name,"Residual_y_vs_xy_plane_z_slice_%d",z_slice);
        Residual_y_vs_xy_plane[z_slice] = new TH2F(Res_y_name,"",nBins_x,x_min,x_max,nBins_y,y_min,y_max);
        sprintf(Res_z_name,"Residual_z_vs_xy_plane_z_slice_%d",z_slice);
        Residual_z_vs_xy_plane[z_slice] = new TH2F(Res_z_name,"",nBins_x,x_min,x_max,nBins_y,y_min,y_max);
    }
    
    // To count the number of entries in each bin
    int entries[nBins_x][nBins_y][nBins_z];
    for(int x=0; x<nBins_x; x++) {
        for(int y=0; y<nBins_y; y++) {
            for(int z=0; z<nBins_z; z++) {
                entries[x][y][z] = 0;
            }
        }
    }

	

    // Fill histograms with distances and residuals
    // =============================================
    double x, y, z;
    double r_x, r_y, r_z;
    double eta, l;
    double sum = 0.;
    int region;
    int x_bin_number = -999;
    int y_bin_number = -999;
    int z_bin_number = -999;
    double x_bin_width = (x_max-x_min)/nBins_x;
    double y_bin_width = (y_max-y_min)/nBins_y;
    double z_bin_width = (z_max-z_min)/nBins_z;
    std::cout << z_min << " \t" << z_max << std::endl;
    inputfile.open("Hits_Residuals_Eta_l_region.txt");
    if(inputfile.is_open()) {
        while( std::getline(inputfile,line) ) {
            istringstream strStream(line);
            strStream >> x >> y >> z >> r_x >> r_y >> r_z >> eta >> l;
            x_bin_number = floor(x/x_bin_width);
            y_bin_number = floor(y/y_bin_width);
            z_bin_number = floor(z/z_bin_width);
            sum += sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2));
            std::cout << x << " \t" << y << " \t" << sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)) << " \t" << sum << std::endl;
            Dist_vs_xy_plane[z_bin_number]->Fill(x,y,sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)));
            Residual_x_vs_xy_plane[z_bin_number]->Fill(x,y,r_x);
            Residual_y_vs_xy_plane[z_bin_number]->Fill(x,y,r_y);
            Residual_z_vs_xy_plane[z_bin_number]->Fill(x,y,r_z);
            entries[x_bin_number][y_bin_number][z_bin_number] += 1;

        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'Hits_Residuals_Eta_l_region.txt'" << std::endl;

    
    // Normalize the bin contents of each voxel to the number of entries in it
    for(int x=0; x<nBins_x; x++) {
        for(int y=0; y<nBins_y; y++) {
            for(int z=0; z<nBins_z; z++) {
                if(entries[x][y][z]!=0) {
                          Dist_vs_xy_plane[z]->SetBinContent(x+1,y+1,      Dist_vs_xy_plane[z]->GetBinContent(x+1,y+1)/entries[x][y][z]);
                    Residual_x_vs_xy_plane[z]->SetBinContent(x+1,y+1,Residual_x_vs_xy_plane[z]->GetBinContent(x+1,y+1)/entries[x][y][z]);
                    Residual_y_vs_xy_plane[z]->SetBinContent(x+1,y+1,Residual_y_vs_xy_plane[z]->GetBinContent(x+1,y+1)/entries[x][y][z]);
                    Residual_z_vs_xy_plane[z]->SetBinContent(x+1,y+1,Residual_z_vs_xy_plane[z]->GetBinContent(x+1,y+1)/entries[x][y][z]);
                }
            }
        }
    }



    // For plots
    // ==========
	gStyle->SetPadRightMargin(0.17);
    //gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetPadBottomMargin(0.14);
    
    double Distances_third_axis_min = -2.5;
    double Distances_third_axis_max =  2.5;
    double Residual_third_axis_min  = -1.5;
    double Residual_third_axis_max  =  1.5;



    // Plot 2D histograms for the distances and residuals
	// ===================================================
    char * Dist_save_name_png = new char[60];
    char * Dist_save_name_pdf = new char[60];
    char * Residual_x_save_name_png = new char[60];
    char * Residual_x_save_name_pdf = new char[60];
    char * Residual_y_save_name_png = new char[60];
    char * Residual_y_save_name_pdf = new char[60];
    char * Residual_z_save_name_png = new char[60];
    char * Residual_z_save_name_pdf = new char[60];
    for(int z_slice=0; z_slice<nBins_z; z_slice++) {
        if(z_slice<10) {
            sprintf(Dist_save_name_png,"plots/z/Dist_vs_xy_plane_z_slice_0%d.png",z_slice);
            sprintf(Dist_save_name_pdf,"plots/z/Dist_vs_xy_plane_z_slice_0%d.pdf",z_slice);
            sprintf(Residual_x_save_name_png,"plots/z/Residual_x_vs_xy_plane_z_slice_0%d.png",z_slice);
            sprintf(Residual_x_save_name_pdf,"plots/z/Residual_x_vs_xy_plane_z_slice_0%d.pdf",z_slice);
            sprintf(Residual_y_save_name_png,"plots/z/Residual_y_vs_xy_plane_z_slice_0%d.png",z_slice);
            sprintf(Residual_y_save_name_pdf,"plots/z/Residual_y_vs_xy_plane_z_slice_0%d.pdf",z_slice);
            sprintf(Residual_z_save_name_png,"plots/z/Residual_z_vs_xy_plane_z_slice_0%d.png",z_slice);
            sprintf(Residual_z_save_name_pdf,"plots/z/Residual_z_vs_xy_plane_z_slice_0%d.pdf",z_slice);
        }
        if(z_slice>9) {
            sprintf(Dist_save_name_png,"plots/z/Dist_vs_xy_plane_z_slice_%d.png",z_slice);
            sprintf(Dist_save_name_pdf,"plots/z/Dist_vs_xy_plane_z_slice_%d.pdf",z_slice);
            sprintf(Residual_x_save_name_png,"plots/z/Residual_x_vs_xy_plane_z_slice_%d.png",z_slice);
            sprintf(Residual_x_save_name_pdf,"plots/z/Residual_x_vs_xy_plane_z_slice_%d.pdf",z_slice);
            sprintf(Residual_y_save_name_png,"plots/z/Residual_y_vs_xy_plane_z_slice_%d.png",z_slice);
            sprintf(Residual_y_save_name_pdf,"plots/z/Residual_y_vs_xy_plane_z_slice_%d.pdf",z_slice);
            sprintf(Residual_z_save_name_png,"plots/z/Residual_z_vs_xy_plane_z_slice_%d.png",z_slice);
            sprintf(Residual_z_save_name_pdf,"plots/z/Residual_z_vs_xy_plane_z_slice_%d.pdf",z_slice);
        }
	
	    TCanvas * canvas_Dist_vs_xy_plane = new TCanvas("canvas_Dist_vs_xy_plane","canvas_Dist_vs_xy_plane");
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        Dist_vs_xy_plane[z_slice]->SetTitle("");
        Dist_vs_xy_plane[z_slice]->GetZaxis()->SetTitleOffset(1.3);
	    Dist_vs_xy_plane[z_slice]->GetZaxis()->SetTitle("mean distance to PC [mm]");
	    Dist_vs_xy_plane[z_slice]->GetZaxis()->SetRangeUser(Distances_third_axis_min,Distances_third_axis_max);
        Dist_vs_xy_plane[z_slice]->GetXaxis()->SetTitleOffset(1.3);
        Dist_vs_xy_plane[z_slice]->GetYaxis()->SetTitleOffset(1.0);
        Dist_vs_xy_plane[z_slice]->GetZaxis()->SetTitleOffset(1.25);
        Dist_vs_xy_plane[z_slice]->GetXaxis()->SetTitleSize(0.05);
        Dist_vs_xy_plane[z_slice]->GetYaxis()->SetTitleSize(0.05);
        Dist_vs_xy_plane[z_slice]->GetZaxis()->SetTitleSize(0.05);
        Dist_vs_xy_plane[z_slice]->GetXaxis()->SetLabelSize(0.05);
        Dist_vs_xy_plane[z_slice]->GetYaxis()->SetLabelSize(0.05);
        Dist_vs_xy_plane[z_slice]->GetZaxis()->SetLabelSize(0.05);
	    Dist_vs_xy_plane[z_slice]->GetXaxis()->SetTitle("IndWireNum [-]");
	    Dist_vs_xy_plane[z_slice]->GetYaxis()->SetTitle("ColWireNum [-]");
        Dist_vs_xy_plane[z_slice]->Draw("colz");
        gPad->RedrawAxis();
	    canvas_Dist_vs_xy_plane->SaveAs(Dist_save_name_png);
	    canvas_Dist_vs_xy_plane->SaveAs(Dist_save_name_pdf);
	    
	    
	    TCanvas * canvas_Residual_x_vs_xy_plane = new TCanvas("canvas_Residual_x_vs_xy_plane","canvas_Residual_x_vs_xy_plane");
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        Residual_x_vs_xy_plane[z_slice]->SetTitle("");
        Residual_x_vs_xy_plane[z_slice]->GetZaxis()->SetTitleOffset(1.3);
	    Residual_x_vs_xy_plane[z_slice]->GetZaxis()->SetTitle("mean residual x to PC [mm]");
	    Residual_x_vs_xy_plane[z_slice]->GetZaxis()->SetRangeUser(Residual_third_axis_min,Residual_third_axis_max);
        Residual_x_vs_xy_plane[z_slice]->GetXaxis()->SetTitleOffset(1.3);
        Residual_x_vs_xy_plane[z_slice]->GetYaxis()->SetTitleOffset(1.0);
        Residual_x_vs_xy_plane[z_slice]->GetZaxis()->SetTitleOffset(1.25);
        Residual_x_vs_xy_plane[z_slice]->GetXaxis()->SetTitleSize(0.05);
        Residual_x_vs_xy_plane[z_slice]->GetYaxis()->SetTitleSize(0.05);
        Residual_x_vs_xy_plane[z_slice]->GetZaxis()->SetTitleSize(0.05);
        Residual_x_vs_xy_plane[z_slice]->GetXaxis()->SetLabelSize(0.05);
        Residual_x_vs_xy_plane[z_slice]->GetYaxis()->SetLabelSize(0.05);
        Residual_x_vs_xy_plane[z_slice]->GetZaxis()->SetLabelSize(0.05);
	    Residual_x_vs_xy_plane[z_slice]->GetXaxis()->SetTitle("IndWireNum [-]");
	    Residual_x_vs_xy_plane[z_slice]->GetYaxis()->SetTitle("ColWireNum [-]");
        Residual_x_vs_xy_plane[z_slice]->Draw("colz");
        gPad->RedrawAxis();
	    canvas_Residual_x_vs_xy_plane->SaveAs(Residual_x_save_name_png);
	    canvas_Residual_x_vs_xy_plane->SaveAs(Residual_x_save_name_pdf);


	    TCanvas * canvas_Residual_y_vs_xy_plane = new TCanvas("canvas_Residual_y_vs_xy_plane","canvas_Residual_y_vs_xy_plane");
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        Residual_y_vs_xy_plane[z_slice]->SetTitle("");
        Residual_y_vs_xy_plane[z_slice]->GetZaxis()->SetTitleOffset(1.3);
	    Residual_y_vs_xy_plane[z_slice]->GetZaxis()->SetTitle("mean residual y to PC [mm]");
	    Residual_y_vs_xy_plane[z_slice]->GetZaxis()->SetRangeUser(Residual_third_axis_min,Residual_third_axis_max);
        Residual_y_vs_xy_plane[z_slice]->GetXaxis()->SetTitleOffset(1.3);
        Residual_y_vs_xy_plane[z_slice]->GetYaxis()->SetTitleOffset(1.0);
        Residual_y_vs_xy_plane[z_slice]->GetZaxis()->SetTitleOffset(1.25);
        Residual_y_vs_xy_plane[z_slice]->GetXaxis()->SetTitleSize(0.05);
        Residual_y_vs_xy_plane[z_slice]->GetYaxis()->SetTitleSize(0.05);
        Residual_y_vs_xy_plane[z_slice]->GetZaxis()->SetTitleSize(0.05);
        Residual_y_vs_xy_plane[z_slice]->GetXaxis()->SetLabelSize(0.05);
        Residual_y_vs_xy_plane[z_slice]->GetYaxis()->SetLabelSize(0.05);
        Residual_y_vs_xy_plane[z_slice]->GetZaxis()->SetLabelSize(0.05);
	    Residual_y_vs_xy_plane[z_slice]->GetXaxis()->SetTitle("IndWireNum [-]");
	    Residual_y_vs_xy_plane[z_slice]->GetYaxis()->SetTitle("ColWireNum [-]");
        Residual_y_vs_xy_plane[z_slice]->Draw("colz");
        gPad->RedrawAxis();
	    canvas_Residual_y_vs_xy_plane->SaveAs(Residual_y_save_name_png);
	    canvas_Residual_y_vs_xy_plane->SaveAs(Residual_y_save_name_pdf);


	    TCanvas * canvas_Residual_z_vs_xy_plane = new TCanvas("canvas_Residual_z_vs_xy_plane","canvas_Residual_z_vs_xy_plane");
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        Residual_z_vs_xy_plane[z_slice]->SetTitle("");
        Residual_z_vs_xy_plane[z_slice]->GetZaxis()->SetTitleOffset(1.3);
	    Residual_z_vs_xy_plane[z_slice]->GetZaxis()->SetTitle("mean residual z to PC [mm]");
	    Residual_z_vs_xy_plane[z_slice]->GetZaxis()->SetRangeUser(Residual_third_axis_min,Residual_third_axis_max);
        Residual_z_vs_xy_plane[z_slice]->GetXaxis()->SetTitleOffset(1.3);
        Residual_z_vs_xy_plane[z_slice]->GetYaxis()->SetTitleOffset(1.0);
        Residual_z_vs_xy_plane[z_slice]->GetZaxis()->SetTitleOffset(1.25);
        Residual_z_vs_xy_plane[z_slice]->GetXaxis()->SetTitleSize(0.05);
        Residual_z_vs_xy_plane[z_slice]->GetYaxis()->SetTitleSize(0.05);
        Residual_z_vs_xy_plane[z_slice]->GetZaxis()->SetTitleSize(0.05);
        Residual_z_vs_xy_plane[z_slice]->GetXaxis()->SetLabelSize(0.05);
        Residual_z_vs_xy_plane[z_slice]->GetYaxis()->SetLabelSize(0.05);
        Residual_z_vs_xy_plane[z_slice]->GetZaxis()->SetLabelSize(0.05);
	    Residual_z_vs_xy_plane[z_slice]->GetXaxis()->SetTitle("IndWireNum [-]");
	    Residual_z_vs_xy_plane[z_slice]->GetYaxis()->SetTitle("ColWireNum [-]");
        Residual_z_vs_xy_plane[z_slice]->Draw("colz");
        gPad->RedrawAxis();
	    canvas_Residual_z_vs_xy_plane->SaveAs(Residual_z_save_name_png);
	    canvas_Residual_z_vs_xy_plane->SaveAs(Residual_z_save_name_pdf);
    }

    return;
}
