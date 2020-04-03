#include <iostream>
#include <fstream>
#include <string>
//using namespace std;


// ***************************************************
void plot_distances() {
// ***************************************************

    gStyle->SetPadRightMargin(0.2);

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


	// Number of wires and time slices
	int n_wires  = 32;
	int n_bins_z = 170; // 150: one for each mm, plus some extra to see if the hit is reconstructed to be outside of the active volume
	int z_min    = 0;
	int z_max    = 170;


    // Create histograms
	// ======================
	// Sliced along the ColWire plane. Have three distance directions (x,y and z)
    TProfile2D * hSliced_along_IndWire[n_wires][3];
    TProfile2D * hSliced_along_ColWire[n_wires][3];

	char * IndWireHistograms_res_x = new char[60];
	char * IndWireHistograms_res_y = new char[60];
	char * IndWireHistograms_res_z = new char[60];
	char * ColWireHistograms_res_x = new char[60];
	char * ColWireHistograms_res_y = new char[60];
	char * ColWireHistograms_res_z = new char[60];

	for(int wire=0; wire<n_wires; wire++) {
		if(wire<10) {
			sprintf(IndWireHistograms_res_x,"Slice_along_IndWire_0%d_distance_x",wire);
			sprintf(IndWireHistograms_res_y,"Slice_along_IndWire_0%d_distance_y",wire);
			sprintf(IndWireHistograms_res_z,"Slice_along_IndWire_0%d_distance_z",wire);
			sprintf(ColWireHistograms_res_x,"Slice_along_ColWire_0%d_distance_x",wire);
			sprintf(ColWireHistograms_res_y,"Slice_along_ColWire_0%d_distance_y",wire);
			sprintf(ColWireHistograms_res_z,"Slice_along_ColWire_0%d_distance_z",wire);
		}
		if(wire>=9) {
			sprintf(IndWireHistograms_res_x,"Slice_along_IndWire_%d_distance_x",wire);
			sprintf(IndWireHistograms_res_y,"Slice_along_IndWire_%d_distance_y",wire);
			sprintf(IndWireHistograms_res_z,"Slice_along_IndWire_%d_distance_z",wire);
			sprintf(ColWireHistograms_res_x,"Slice_along_ColWire_%d_distance_x",wire);
			sprintf(ColWireHistograms_res_y,"Slice_along_ColWire_%d_distance_y",wire);
			sprintf(ColWireHistograms_res_z,"Slice_along_ColWire_%d_distance_z",wire);
		}
		hSliced_along_IndWire[wire][0] = new TProfile2D(IndWireHistograms_res_x,"Slices along IndWires: Distance component x",n_wires+2,-1.5,n_wires+0.5,n_bins_z,z_min,z_max,-20,20);
		hSliced_along_IndWire[wire][1] = new TProfile2D(IndWireHistograms_res_y,"Slices along IndWires: Distance component y",n_wires+2,-1.5,n_wires+0.5,n_bins_z,z_min,z_max,-20,20);
		hSliced_along_IndWire[wire][2] = new TProfile2D(IndWireHistograms_res_z,"Slices along IndWires: Distance component z",n_wires+2,-1.5,n_wires+0.5,n_bins_z,z_min,z_max,-20,20);
		hSliced_along_ColWire[wire][0] = new TProfile2D(ColWireHistograms_res_x,"Slices along ColWires: Distance component x",n_wires+2,-1.5,n_wires+0.5,n_bins_z,z_min,z_max,-20,20);
		hSliced_along_ColWire[wire][1] = new TProfile2D(ColWireHistograms_res_y,"Slices along ColWires: Distance component y",n_wires+2,-1.5,n_wires+0.5,n_bins_z,z_min,z_max,-20,20);
		hSliced_along_ColWire[wire][2] = new TProfile2D(ColWireHistograms_res_z,"Slices along ColWires: Distance component z",n_wires+2,-1.5,n_wires+0.5,n_bins_z,z_min,z_max,-20,20);
    }


    // Fill histograms with hits and distances
    double x, y, z;
    double r_x, r_y, r_z;
    double eta, l;
    int region;
    inputfile.open("Hits_Residuals_Eta_l_region.txt");
    if(inputfile.is_open()) {
        while( std::getline(inputfile,line) ) {
            istringstream strStream(line);
            strStream >> x >> y >> z >> r_x >> r_y >> r_z >> eta >> l >> region;
			for(int wire=0; wire<n_wires; wire++) {
				// Slices along IndWire planes
		        if(x==wire) {
					hSliced_along_IndWire[wire][0]->Fill(x,z,fabs(r_x));
					hSliced_along_IndWire[wire][1]->Fill(x,z,fabs(r_y));
					hSliced_along_IndWire[wire][2]->Fill(x,z,fabs(r_z));
				}
				// Slices along ColWire planes
		        if(y==wire) {
					hSliced_along_ColWire[wire][0]->Fill(y,z,fabs(r_x));
					hSliced_along_ColWire[wire][1]->Fill(y,z,fabs(r_y));
					hSliced_along_ColWire[wire][2]->Fill(y,z,fabs(r_z));
				}
			}
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'Hits_Residuals_Eta_l_region.txt'" << std::endl;


	double x_axis_min = 0.;
	double x_axis_max = 1.;
	double y_axis_min = 0.;
	double y_axis_max = 1.;
	double z_axis_min = 0.;
	double z_axis_max = 0.3;


    // Plot 2D histograms for the distances
	// ====================================

	// Integrated over all ColWire planes
	// -----------------------------------
    TCanvas * ColWires_integrated_00 = new TCanvas("ColWires_integrated_00","ColWires_integrated_00");
    gStyle->SetOptStat(0);
	for(int wire=0; wire<n_wires; wire++) {
		hSliced_along_ColWire[wire][0]->GetZaxis()->SetTitleOffset(1.3);
		hSliced_along_ColWire[wire][0]->GetZaxis()->SetTitle("distance to principal component [wire pitches]");
		hSliced_along_ColWire[wire][0]->GetZaxis()->SetRangeUser(x_axis_min,x_axis_max);
		if(wire==0) {
			hSliced_along_ColWire[wire][0]->GetXaxis()->SetTitleOffset(1.4);
			hSliced_along_ColWire[wire][0]->GetYaxis()->SetTitleOffset(1.4);
			hSliced_along_ColWire[wire][0]->GetXaxis()->SetTitle("ColWireNum [-]");
			hSliced_along_ColWire[wire][0]->GetYaxis()->SetTitle("z [mm]"); // (fMeanTime/20 * drift_vel)
			hSliced_along_ColWire[wire][0]->Draw("COLZ");
		}
		else hSliced_along_ColWire[wire][0]->Draw("COLZ same");
	}
    gPad->RedrawAxis();
	ColWires_integrated_00->SaveAs("plots/distances/ColWires_integrated_distance_x.png");
	ColWires_integrated_00->SaveAs("plots/distances/ColWires_integrated_distance_x.pdf");
    ColWires_integrated_00->SaveAs("plots/distances/ColWires_integrated_distance_x.root");

    TCanvas * ColWires_integrated_01 = new TCanvas("ColWires_integrated_01","ColWires_integrated_01");
    gStyle->SetOptStat(0);
	for(int wire=0; wire<n_wires; wire++) {
		hSliced_along_ColWire[wire][1]->GetZaxis()->SetTitleOffset(1.3);
		hSliced_along_ColWire[wire][1]->GetZaxis()->SetTitle("distance to principal component [wire pitches]");
		hSliced_along_ColWire[wire][1]->GetZaxis()->SetRangeUser(y_axis_min,y_axis_max);
		if(wire==0) {
			hSliced_along_ColWire[wire][1]->GetXaxis()->SetTitleOffset(1.4);
			hSliced_along_ColWire[wire][1]->GetYaxis()->SetTitleOffset(1.4);
			hSliced_along_ColWire[wire][1]->GetXaxis()->SetTitle("ColWireNum [-]");
			hSliced_along_ColWire[wire][1]->GetYaxis()->SetTitle("z [mm]"); // (fMeanTime/20 * drift_vel)
			hSliced_along_ColWire[wire][1]->Draw("COLZ");
		}
		else hSliced_along_ColWire[wire][1]->Draw("COLZ same");
	}
    gPad->RedrawAxis();
	ColWires_integrated_01->SaveAs("plots/distances/ColWires_integrated_distance_y.png");
	ColWires_integrated_01->SaveAs("plots/distances/ColWires_integrated_distance_y.pdf");
    ColWires_integrated_01->SaveAs("plots/distances/ColWires_integrated_distance_y.root");

    TCanvas * ColWires_integrated_02 = new TCanvas("ColWires_integrated_02","ColWires_integrated_02");
    gStyle->SetOptStat(0);
	for(int wire=0; wire<n_wires; wire++) {
		hSliced_along_ColWire[wire][2]->GetZaxis()->SetTitleOffset(1.3);
		hSliced_along_ColWire[wire][2]->GetZaxis()->SetTitle("distance to principal component [wire pitches]");
		hSliced_along_ColWire[wire][2]->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
		if(wire==0) {
			hSliced_along_ColWire[wire][2]->GetXaxis()->SetTitleOffset(1.4);
			hSliced_along_ColWire[wire][2]->GetYaxis()->SetTitleOffset(1.4);
			hSliced_along_ColWire[wire][2]->GetXaxis()->SetTitle("ColWireNum [-]");
			hSliced_along_ColWire[wire][2]->GetYaxis()->SetTitle("z [mm]"); // (fMeanTime/20 * drift_vel)
			hSliced_along_ColWire[wire][2]->Draw("COLZ");
		}
		else hSliced_along_ColWire[wire][2]->Draw("COLZ same");
	}
    gPad->RedrawAxis();
	ColWires_integrated_02->SaveAs("plots/distances/ColWires_integrated_distance_z.png");
	ColWires_integrated_02->SaveAs("plots/distances/ColWires_integrated_distance_z.pdf");
    ColWires_integrated_02->SaveAs("plots/distances/ColWires_integrated_distance_z.root");


	// Integrated over all IndWire planes
	// -----------------------------------
    TCanvas * IndWires_integrated_00 = new TCanvas("IndWires_integrated_00","IndWires_integrated_00");
    gStyle->SetOptStat(0);
	for(int wire=0; wire<n_wires; wire++) {
		hSliced_along_IndWire[wire][0]->GetZaxis()->SetTitleOffset(1.3);
		hSliced_along_IndWire[wire][0]->GetZaxis()->SetTitle("distance to principal component [wire pitches]");
		hSliced_along_IndWire[wire][0]->GetZaxis()->SetRangeUser(x_axis_min,x_axis_max);
		if(wire==0) {
			hSliced_along_IndWire[wire][0]->GetXaxis()->SetTitleOffset(1.4);
			hSliced_along_IndWire[wire][0]->GetYaxis()->SetTitleOffset(1.4);
			hSliced_along_IndWire[wire][0]->GetXaxis()->SetTitle("IndWireNum [-]");
			hSliced_along_IndWire[wire][0]->GetYaxis()->SetTitle("z [mm]"); // (fMeanTime/20 * drift_vel)
			hSliced_along_IndWire[wire][0]->Draw("COLZ");
		}
		else hSliced_along_IndWire[wire][0]->Draw("COLZ same");
	}
    gPad->RedrawAxis();
	IndWires_integrated_00->SaveAs("plots/distances/IndWires_integrated_distance_x.png");
	IndWires_integrated_00->SaveAs("plots/distances/IndWires_integrated_distance_x.pdf");
    IndWires_integrated_00->SaveAs("plots/distances/IndWires_integrated_distance_x.root");

    TCanvas * IndWires_integrated_01 = new TCanvas("IndWires_integrated_01","IndWires_integrated_01");
    gStyle->SetOptStat(0);
	for(int wire=0; wire<n_wires; wire++) {
		hSliced_along_IndWire[wire][1]->GetZaxis()->SetTitleOffset(1.3);
		hSliced_along_IndWire[wire][1]->GetZaxis()->SetTitle("distance to principal component [wire pitches]");
		hSliced_along_IndWire[wire][1]->GetZaxis()->SetRangeUser(y_axis_min,y_axis_max);
		if(wire==0) {
			hSliced_along_IndWire[wire][1]->GetXaxis()->SetTitleOffset(1.4);
			hSliced_along_IndWire[wire][1]->GetYaxis()->SetTitleOffset(1.4);
			hSliced_along_IndWire[wire][1]->GetXaxis()->SetTitle("IndWireNum [-]");
			hSliced_along_IndWire[wire][1]->GetYaxis()->SetTitle("z [mm]"); // (fMeanTime/20 * drift_vel)
			hSliced_along_IndWire[wire][1]->Draw("COLZ");
		}
		else hSliced_along_IndWire[wire][1]->Draw("COLZ same");
	}
    gPad->RedrawAxis();
	IndWires_integrated_01->SaveAs("plots/distances/IndWires_integrated_distance_y.png");
	IndWires_integrated_01->SaveAs("plots/distances/IndWires_integrated_distance_y.pdf");
    IndWires_integrated_01->SaveAs("plots/distances/IndWires_integrated_distance_y.root");

    TCanvas * IndWires_integrated_02 = new TCanvas("IndWires_integrated_02","IndWires_integrated_02");
    gStyle->SetOptStat(0);
	for(int wire=0; wire<n_wires; wire++) {
		hSliced_along_IndWire[wire][2]->GetZaxis()->SetTitleOffset(1.3);
		hSliced_along_IndWire[wire][2]->GetZaxis()->SetTitle("distance to principal component [wire pitches]");
		hSliced_along_IndWire[wire][2]->GetZaxis()->SetRangeUser(z_axis_min,z_axis_max);
		if(wire==0) {
			hSliced_along_IndWire[wire][2]->GetXaxis()->SetTitleOffset(1.4);
			hSliced_along_IndWire[wire][2]->GetYaxis()->SetTitleOffset(1.4);
			hSliced_along_IndWire[wire][2]->GetXaxis()->SetTitle("IndWireNum [-]");
			hSliced_along_IndWire[wire][2]->GetYaxis()->SetTitle("z [mm]"); // (fMeanTime/20 * drift_vel)
			hSliced_along_IndWire[wire][2]->Draw("COLZ");
		}
		else hSliced_along_IndWire[wire][2]->Draw("COLZ same");
	}
    gPad->RedrawAxis();
	IndWires_integrated_02->SaveAs("plots/distances/IndWires_integrated_distance_z.png");
	IndWires_integrated_02->SaveAs("plots/distances/IndWires_integrated_distance_z.pdf");
    IndWires_integrated_02->SaveAs("plots/distances/IndWires_integrated_distance_z.root");


    return;
}
