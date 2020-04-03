#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>

// ROOT includes
#include "TCanvas.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TStyle.h"


// **************************************************************
//std::vector<double> E_field_correction(float hit_x, float hit_y, float hit_z) {
// **************************************************************

//int main() {
void E_field_correction() {

	// Initialize the vector with the E-field corrected hit position
	std::vector<double> hit_position(3,-999.);

    // ==================================================== //
    // ==================================================== //
    // HISTOGRAMS DO NOT ALLOW TO HAVE NEGATIVE ENTRIES     //
    // THIS METHOD THUS CANNOT BE USED
    // ==================================================== //
    // ==================================================== //
/*
    // Should the .root file of the electric field be created?
    bool create_root_file = true;

	// Divide TPC into .. x .. x .. voxels (voxel size: .. x .. x ..)
	// Produce TH3D with voxels which have the electric field information in each voxel
	double x, y, z;
	double E_x, E_y, E_z;
	int n_bins_x =  80;
	double x_min =  -0.005;
	double x_max =   0.075;
	int n_bins_y =  80;
	double y_min =  -0.005;
	double y_max =  0.075;
	int n_bins_z = 150;
	double z_min =  0.0; //-0.01;
	double z_max =  0.15; // 0.16;
	//TH3D * hist_E_x = new TH3D("hist_E_x","hist_E_x",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
	//TH3D * hist_E_y = new TH3D("hist_E_y","hist_E_y",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
	//TH3D * hist_E_z = new TH3D("hist_E_z","hist_E_z",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
	//TProfile2D * projection_xy = new TProfile2D("projection_xy","projection_xy",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,-100,100);

	// Create the .root file with the information about the electric field
	// --------------------------------------------------------------------
    if(create_root_file) {
        gBenchmark->Start("hsimple");

        // Define some histograms to be filled
        TH3D * hist_E_x = new TH3D("hist_E_x","hist_E_x",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
	    TH3D * hist_E_y = new TH3D("hist_E_y","hist_E_y",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
	    TH3D * hist_E_z = new TH3D("hist_E_z","hist_E_z",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
        TH3D * hist_E_x_entries = new TH3D("hist_E_x_entries","hist_E_x_entries",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
	    TH3D * hist_E_y_entries = new TH3D("hist_E_y_entries","hist_E_y_entries",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
	    TH3D * hist_E_z_entries = new TH3D("hist_E_z_entries","hist_E_z_entries",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
        
        // Read the .txt file with the information about the electric field and save data in histograms
	    std::string line;
	    std::ifstream inputfile;
        inputfile.open("E_field_without_header.txt");
        int linenumber = 0;
        if(inputfile.is_open()) {
            while( std::getline(inputfile,line) ) {
                std::istringstream strStream(line);
                strStream >> x >> y >> z >> E_x >> E_y >> E_z;
                //std::cout << " " << x << " " << y << " " << z << " " << E_x << " \t " << E_y << " \t " << E_z << std::endl;
                // Note that the drift volume in COMSOL has been simulated to be in the coordinates x:[-0.035,0.035]; y:[-0.035,0.035]; z:[0.000(cathode),0.150(anode)]
                // -> Transform those coordinates to x,y,z of the measurements:
                //    shift x and y to [0.000,0.07]
                //    For z: invert and add 0.150 such that 0 mm corresponds to the anode and 150 mm to the cathode
                hist_E_x->Fill(x+0.035,y+0.035,0.15-z,E_x); // Note: this only adds the E_x value to the already stored value in the histogram.
                                                            // Will have to normalize (divide by the number of entries) to obtain the correct E_x field value
                hist_E_y->Fill(x+0.035,y+0.035,0.15-z,E_y);
                hist_E_z->Fill(x+0.035,y+0.035,0.15-z,E_z);
                hist_E_x_entries->Fill(x+0.035,y+0.035,0.15-z); // To count the entries in each bin (for the normalization)
                hist_E_y_entries->Fill(x+0.035,y+0.035,0.15-z);
                hist_E_z_entries->Fill(x+0.035,y+0.035,0.15-z);

                std::cout << " linenumber: " << linenumber << std::endl;
                linenumber ++;
                //if(linenumber > 200) break;
            }
            inputfile.close();
        }
        else std::cout << " Could not open file 'E_field_without_header.txt'" << std::endl;

        // Define the .root file
        TFile * root_file = new TFile("E_field.root","recreate");

        // Define the tree in the file
        //auto Tree_Efield = new TTree("Tree_Efield","Tree_Efield");

        // Define the branch
        //Tree_Efield->Branch("hist_E_x","hist_E_x",hist_E_x,32000,0);
        //Tree_Efield->Branch("hist_E_y","hist_E_y",hist_E_y,32000,0);
        //Tree_Efield->Branch("hist_E_z","hist_E_z",hist_E_z,32000,0);
        //auto hpx   = new TH1F("hpx","This is the px distribution",100,-4,4);
        //auto hpxpy = new TH2F("hpxpy","py vs px",40,-4,4,40,-4,4);
        //auto hprof = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);
        //T->Branch("hpx","TH1F",&hpx,32000,0);
        //T->Branch("hpxpy","TH2F",&hpxpy,32000,0);
        //T->Branch("hprof","TProfile",&hprof,32000,0);

        // Define the histogram which should be saved in the .root file
        TH3D * hist_E_x_normalized = new TH3D("hist_E_x_normalized","hist_E_x_normalized",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
        TH3D * hist_E_y_normalized = new TH3D("hist_E_y_normalized","hist_E_y_normalized",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);
        TH3D * hist_E_z_normalized = new TH3D("hist_E_z_normalized","hist_E_z_normalized",n_bins_x,x_min,x_max,n_bins_y,y_min,y_max,n_bins_z,z_min,z_max);

        //projection_xy->Fill(x,y,E_x,1.);
        //Tree_Efield->Fill();

        // Loop over all bins and apply the normalization of the E_x_field
        for(int bin_x=1; bin_x<=n_bins_x; bin_x++) {
            for(int bin_y=1; bin_y<=n_bins_y; bin_y++) {
                for(int bin_z=1; bin_z<=n_bins_z; bin_z++) {
                    Int_t bin_E_x = hist_E_x->GetBin(bin_x,bin_y,bin_z);
                    Int_t bin_E_y = hist_E_y->GetBin(bin_x,bin_y,bin_z);
                    Int_t bin_E_z = hist_E_z->GetBin(bin_x,bin_y,bin_z);
                    if(hist_E_x->GetBinContent(bin_E_x)!=0) {
                        hist_E_x_normalized->Fill(  x_min+(bin_x-1)*(x_max-x_min)/n_bins_x,
                                                    y_min+(bin_y-1)*(y_max-y_min)/n_bins_y,
                                                    z_min+(bin_z-1)*(z_max-z_min)/n_bins_z,
                                                    hist_E_x->GetBinContent(bin_E_x)/hist_E_x_entries->GetBinContent(bin_E_x) );
                    }
                    if(hist_E_y->GetBinContent(bin_E_y)!=0) {
                        hist_E_y_normalized->Fill(  x_min+(bin_x-1)*(x_max-x_min)/n_bins_x,
                                                    y_min+(bin_y-1)*(y_max-y_min)/n_bins_y,
                                                    z_min+(bin_z-1)*(z_max-z_min)/n_bins_z,
                                                    hist_E_y->GetBinContent(bin_E_y)/hist_E_y_entries->GetBinContent(bin_E_y) );
                    }
                    if(hist_E_z->GetBinContent(bin_E_z)!=0) {
                        hist_E_z_normalized->Fill(  x_min+(bin_x-1)*(x_max-x_min)/n_bins_x,
                                                    y_min+(bin_y-1)*(y_max-y_min)/n_bins_y,
                                                    z_min+(bin_z-1)*(z_max-z_min)/n_bins_z,
                                                    hist_E_z->GetBinContent(bin_E_z)/hist_E_z_entries->GetBinContent(bin_E_z) );
                    }
                }
            }
        }

        //Tree_Efield->Print();
        root_file->Write();

        gBenchmark->Show("hsimple");

    } // end if(create_root_file)
*/
    // ==================================================== //
    // ==================================================== //
    // END HISTOGRAM                                        //
    // ==================================================== //
    // ==================================================== //

    // Define array to store the electric field information
    int n_bins_x = 7;
    int n_bins_y = 7;
    int n_bins_z = 15;
    
    double x_min = 0.; //-0.005;
	double x_max = 0.07; // 0.075;
	double y_min = 0.; //-0.005;
	double y_max = 0.07; // 0.075;
	double z_min = 0.; //-0.005;
	double z_max = 0.15; // 0.155;
	
    double wire_pitch = 52.5 / 31; // [mm]
    
    double voxel_size_x = (x_max - x_min) / n_bins_x;
    double voxel_size_y = (y_max - y_min) / n_bins_y;
    double voxel_size_z = (z_max - z_min) / n_bins_z;
    
    std::cout << " Voxel size: \t x: " << voxel_size_x << " \t y: " << voxel_size_y << " \t z: " << voxel_size_z << std::endl;

    // Allocate memory for all voxels:
    // arr_E_x for the field vectors x component and
    // arr_E_x_entries to count the entries for a normalization later
	unsigned int ***arr_E_x = (unsigned int ***)malloc(n_bins_x * sizeof(float **));
	int ***arr_E_y = (int ***)malloc(n_bins_x * sizeof(float **));
	int ***arr_E_z = (int ***)malloc(n_bins_x * sizeof(float **));
    int ***arr_E_x_entries = (int ***)malloc(n_bins_x * sizeof(int **));
    int ***arr_E_y_entries = (int ***)malloc(n_bins_x * sizeof(int **));
    int ***arr_E_z_entries = (int ***)malloc(n_bins_x * sizeof(int **));
	for(int i = 0; i < n_bins_x; i++) {
		arr_E_x[i] = (int **)malloc(n_bins_y * sizeof(float *));
		arr_E_y[i] = (int **)malloc(n_bins_y * sizeof(float *));
		arr_E_z[i] = (int **)malloc(n_bins_y * sizeof(float *));
		arr_E_x_entries[i] = (int **)malloc(n_bins_y * sizeof(int *));
		arr_E_y_entries[i] = (int **)malloc(n_bins_y * sizeof(int *));
		arr_E_z_entries[i] = (int **)malloc(n_bins_y * sizeof(int *));
		for(int j = 0; j < n_bins_y; j++) {
			arr_E_x[i][j] = (int *)malloc(n_bins_z * sizeof(float));
			arr_E_y[i][j] = (int *)malloc(n_bins_z * sizeof(float));
			arr_E_z[i][j] = (int *)malloc(n_bins_z * sizeof(float));
			arr_E_x_entries[i][j] = (int *)malloc(n_bins_z * sizeof(int));
			arr_E_y_entries[i][j] = (int *)malloc(n_bins_z * sizeof(int));
			arr_E_z_entries[i][j] = (int *)malloc(n_bins_z * sizeof(int));
		}
	}

    // Read the .txt file with the information about the electric field and save data in histograms
	std::string line;
	std::ifstream inputfile;
    inputfile.open("E_field_without_header.txt");
    int linenumber = 0;
    double x, y, z;
	double E_x, E_y, E_z;
	int index_x, index_y, index_z;

    if(inputfile.is_open()) {
        while( std::getline(inputfile,line) ) {
            std::istringstream strStream(line);
            strStream >> x >> y >> z >> E_x >> E_y >> E_z;
            //std::cout << " " << x << " " << y << " " << z << " " << E_x << " \t " << E_y << " \t " << E_z << std::endl;
            // Note that the drift volume in COMSOL has been simulated to be in the coordinates x:[-0.035,0.035]; y:[-0.035,0.035]; z:[0.000(cathode),0.150(anode)]
            // -> Transform those coordinates to x,y,z of the measurements:
            //    shift x and y to [0.000,0.07]
            //    For z: invert and add 0.150 such that 0 mm corresponds to the anode and 150 mm to the cathode
            x = x + 0.035;
            y = y + 0.035;
            z = 0.15 - z;

            // Get the index of the array:
            //      All values <x_min (y_min,z_min) go to bin -1 (underflow)
            //      All values >=x_max (y_max,z_max) go to bin n_bins_x,y,z (overflow)
            for(int bin_index=0; bin_index<n_bins_x; bin_index++) {
                if( x<x_min ) index_x = -1;
                if( x>=(x_min+bin_index*(x_max-x_min)/n_bins_x) && x<(x_min+(bin_index+1)*(x_max-x_min)/n_bins_x)) {
                    index_x = bin_index;
                }
                if(x>=x_max) index_x = n_bins_x;
            }
            //std::cout << " x: " << x << " \t index_x: " << index_x << std::endl;

            for(int bin_index=0; bin_index<n_bins_y; bin_index++) {
                if( y<y_min ) index_y = -1;
                if( (y>=(y_min+bin_index*(y_max-y_min)/n_bins_y)) && (y<(y_min+(bin_index+1)*(y_max-y_min)/n_bins_y)) ) {
                    index_y = bin_index;
                }
                if(y>=y_max) index_y = n_bins_y;
            }
            //std::cout << " y: " << y << " \t index_y: " << index_y << std::endl;

            for(int bin_index=0; bin_index<n_bins_z; bin_index++) {
                if( z<z_min ) index_z = -1;
                if( z>=(z_min+bin_index*(z_max-z_min)/n_bins_z) && z<(z_min+(bin_index+1)*(z_max-z_min)/n_bins_z)) {
                    index_z = bin_index;
                }
                if(z>=z_max) index_z = n_bins_z;
            }
            //std::cout << " z: " << z << " \t index_z: " << index_z << std::endl;



            if( index_x<n_bins_x && index_x>=0 && index_y<n_bins_y && index_y>=0 && index_z<n_bins_z && index_z>=0 ) {
                //std::cout << " > 0 0 0 " << std::endl;
                // Fill the array (array = array + value)
                arr_E_x[index_x][index_y][index_z] = arr_E_x[index_x][index_y][index_z] + E_x;
                arr_E_y[index_x][index_y][index_z] = arr_E_x[index_x][index_y][index_z] + E_y;
                arr_E_z[index_x][index_y][index_z] = arr_E_x[index_x][index_y][index_z] + E_z;
                // Fill also an array to count the number of entries (for the normalization)
                arr_E_x_entries[index_x][index_y][index_z] += 1;
                arr_E_y_entries[index_x][index_y][index_z] += 1;
                arr_E_z_entries[index_x][index_y][index_z] += 1;
            }

            if(linenumber%10000==0) std::cout << " linenumber: " << linenumber << std::endl;
            linenumber ++;
            //if(linenumber > 200) break;
        }
        inputfile.close();
    }
    else std::cout << " Could not open file 'E_field_without_header.txt'" << std::endl;


    /*
    // Print some entries for testing purpose
    for(int bin_index_x=0; bin_index_x<n_bins_x; bin_index_x++) {
        for(int bin_index_y=0; bin_index_y<n_bins_y; bin_index_y++) {
            for(int bin_index_z=0; bin_index_z<n_bins_z; bin_index_z++) {
                if(arr_E_x[bin_index_x][bin_index_y][bin_index_z]>0) std::cout << " x: " << bin_index_x << " \t y: " << bin_index_y << " \t z: " << bin_index_z << " \t arr_E_x: " << arr_E_x[bin_index_x][bin_index_y][bin_index_z] << std::endl;
            }
        }
    }
    */


    // Normalize the arr_E_x,y,z entries
    for(int bin_index_x=0; bin_index_x<n_bins_x; bin_index_x++) {
        for(int bin_index_y=0; bin_index_y<n_bins_y; bin_index_y++) {
            for(int bin_index_z=0; bin_index_z<n_bins_z; bin_index_z++) {
                if(arr_E_x_entries[bin_index_x][bin_index_y][bin_index_z]!=0) {
                    arr_E_x[bin_index_x][bin_index_y][bin_index_z] = arr_E_x[bin_index_x][bin_index_y][bin_index_z] / arr_E_x_entries[bin_index_x][bin_index_y][bin_index_z];
                }
                if(arr_E_y_entries[bin_index_x][bin_index_y][bin_index_z]!=0) {
                    arr_E_y[bin_index_x][bin_index_y][bin_index_z] = arr_E_y[bin_index_x][bin_index_y][bin_index_z] / arr_E_y_entries[bin_index_x][bin_index_y][bin_index_z];
                }
                if(arr_E_z_entries[bin_index_x][bin_index_y][bin_index_z]!=0) {
                    arr_E_z[bin_index_x][bin_index_y][bin_index_z] = arr_E_z[bin_index_x][bin_index_y][bin_index_z] / arr_E_z_entries[bin_index_x][bin_index_y][bin_index_z];
                }
                // Set empty bins to the default value for the E_z field:
                // If there is x and y and z information in the pairwise neighboring bins, take the mean value of them
                //if(arr_E_z[bin_index_x][bin_index_y][bin_index_z]==0) {
                    //arr_E_z[bin_index_x][bin_index_y][bin_index_z]=-150000.;
                    //if(bin_index_z==0) { arr_E_z[bin_index_x][bin_index_y][bin_index_z] = arr_E_z[bin_index_x][bin_index_y][bin_index_z+1]; }
                    //if(bin_index_z!=0 && bin_index_z!=n_bins_z) { arr_E_z[bin_index_x][bin_index_y][bin_index_z] = (arr_E_z[bin_index_x][bin_index_y][bin_index_z-1] + arr_E_z[bin_index_x][bin_index_y][bin_index_z+1]) / 2; }
                    //if(bin_index_z==150) { arr_E_z[bin_index_x][bin_index_y][bin_index_z] = arr_E_z[bin_index_x][bin_index_y][bin_index_z-1]; }
                //}
                //if(arr_E_x[bin_index_x][bin_index_y][bin_index_z]==0) {
                    std::cout << " x: " << bin_index_x <<
                                 " \t y: " << bin_index_y <<
                                 " \t z: " << bin_index_z <<
                                 " \t arr_E_x: " << arr_E_x[bin_index_x][bin_index_y][bin_index_z] <<
                                 " \t arr_E_y: " << arr_E_y[bin_index_x][bin_index_y][bin_index_z] <<
                                 " \t arr_E_z: " << arr_E_z[bin_index_x][bin_index_y][bin_index_z] << std::endl;
                //}*/
                /*if( arr_E_x_entries[bin_index_x][bin_index_y][bin_index_z]!=0 ||
                    arr_E_y_entries[bin_index_x][bin_index_y][bin_index_z]!=0 ||
                    arr_E_z_entries[bin_index_x][bin_index_y][bin_index_z]!=0 ) {
                    std::cout << " Zero E-field at: x=" << bin_index_x << " \t y=" << bin_index_y << " \t z=" << bin_index_z << std::endl;
                    */
                //}
            }
        }
    }

    for(int bin_index_x=0; bin_index_x<n_bins_x; bin_index_x++) {
        for(int bin_index_y=0; bin_index_y<n_bins_y; bin_index_y++) {
            for(int bin_index_z=0; bin_index_z<n_bins_z; bin_index_z++) {
                std::cout << " x: " << bin_index_x <<
                " \t y: " << bin_index_y <<
                " \t z: " << bin_index_z <<
                " \t arr_E_x_entries: " << arr_E_x_entries[bin_index_x][bin_index_y][bin_index_z] <<
                " \t arr_E_y_entries: " << arr_E_y_entries[bin_index_x][bin_index_y][bin_index_z] <<
                " \t arr_E_z_entries: " << arr_E_z_entries[bin_index_x][bin_index_y][bin_index_z] << std::endl;
            }
        }
    }



    // +++++++++++++++++++++++++++++++++++++++++++++++++
    // NOW THERE ARE MANY VOXELS WHICH HAVE ZERO FIELD!!
    // +++++++++++++++++++++++++++++++++++++++++++++++++



    // Project 3D histogram to 2D to check plausability
    // ------------------------------------------------
    //gStyle->SetOptStat(0);
    //hist_E_x->Draw();
    //projection_xy->Draw("COLZ");
    //TH1D * projX = projection_xy->ProjectionX();
    //projX->Draw();
    //TH2D * projXZ = hist_E_x->Project3DProfile("zx");
    //projXZ->Draw("COLZ");


    // Propagate the hit along the electric field lines until it reaches it's hit_z position
    // --------------------------------------------------------------------------------------
    // hit_x and hit_y define the start voxel
    float hit_x = 1;        // unit: wire number
    float hit_y = 1;        // unit: wire number
    float hit_z = 149.9;  // unit: mm (0: at cathode; 150: at anode)

    // Convert the hits from wire numbers to absolute (x,y) coordinates, units: m
    hit_x = (hit_x * 52.5/31 + (70-52.5)/2)/1000;
    hit_y = (hit_y * 52.5/31 + (70-52.5)/2)/1000;
    hit_z = hit_z/1000;

    // Information for unphysical events
    if(hit_z<0) {
        std::cout << " hit_z at position z = " << hit_z << " [mm], outside of TPC drift volume. Set hit_z to 0 mm. " << std::endl;
        hit_z = 0.;
    }
    if(hit_z>150) {
        std::cout << " hit_z at position z = " << hit_z << " [mm], outside of TPC drift volume. Set hit_z to 150 mm. " << std::endl;
        hit_z = 150.;
    }

    std::cout << " Input coordinates: \t x: " << hit_x << " \t y: " << hit_y << " \t z: " << hit_z << std::endl;

    // Get the bin number:
    for(int bin_index=0; bin_index<n_bins_x; bin_index++) {
        if( hit_x<x_min ) index_x = -1;
        if( hit_x>=(x_min+bin_index*(x_max-x_min)/n_bins_x) && hit_x<(x_min+(bin_index+1)*(x_max-x_min)/n_bins_x) ) {
            index_x = bin_index;
        }
        if( hit_x>=x_max ) index_x = n_bins_x;
    }
    for(int bin_index=0; bin_index<n_bins_y; bin_index++) {
        if( hit_y<y_min ) index_y = -1;
        if( hit_y>=(y_min+bin_index*(y_max-y_min)/n_bins_y) && hit_y<(y_min+(bin_index+1)*(y_max-y_min)/n_bins_y) ) {
            index_y = bin_index;
        }
        if( hit_y>=y_max ) index_y = n_bins_y;
    }
    for(int bin_index=0; bin_index<n_bins_z; bin_index++) {
        if( hit_z<z_min ) index_z = -1;
        if( hit_z>=(z_min+bin_index*(z_max-z_min)/n_bins_z) && hit_z<(z_min+(bin_index+1)*(z_max-z_min)/n_bins_z) ) {
            index_z = bin_index;
        }
        if( hit_z>=z_max) index_z = n_bins_z;
    }

    std::cout << " Input bin indices: \t x: " << index_x << " \t \t y: " << index_y << " \t \t z: " << index_z << std::endl;

    // Define vector with spatial coordinates of the E-field corrected hit
    std::vector<double> position_corrected(3,-999.9);

    double actual_x_position = hit_x;
    double actual_y_position = hit_y;
    double actual_z_position = 0.5 * voxel_size_z;
    
    std::cout << " Voxel size: \t\t x: " << voxel_size_x << " \t y: " << voxel_size_y << " \t z: " << voxel_size_z << std::endl;
    std::cout << " Starting position: \t x: " << actual_x_position << " \t y: " << actual_y_position << " \t z: " << actual_z_position << std::endl;

    // Propagate the hit along the electric field line to its end position

    //while(actual_z_position<hit_z) {
    for(int z_slice=0; z_slice<n_bins_z; z_slice++) {
        std::cout << " --------------------------------- " << std::endl;
        std::cout << " z slice number: " << z_slice << " of " << n_bins_z-1 << std::endl;
        std::cout << " --------------------------------- " << std::endl;
        
        
    
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // ADD HERE AN ABORTION CRITERION TO ABORT WHEN THE Z_RECO > HIT_Z
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        // Information of the current voxel:
        std::cout << " Position: x: " << actual_x_position << "\t  y: " << actual_y_position << " \t  z: " << actual_z_position << std::endl;
        std::cout << " Bin ind.: x: " << index_x << " \t  y: " << index_y << " \t  z: " << index_z << std::endl;
        std::cout << " E-field: Ex: " << arr_E_x[index_x][index_y][index_z] << " \t Ey: " << arr_E_y[index_x][index_y][index_z] << " \t Ez: " << arr_E_z[index_x][index_y][index_z] << std::endl;

        if(arr_E_z[index_x][index_y][index_z]!=0) hit_x = hit_x + voxel_size_z/(arr_E_z[index_x][index_y][index_z])*arr_E_x[index_x][index_y][index_z];
        //if(arr_E_z[index_x][index_y][index_z]==0) hit_x = hit_x + arr_E_x[index_x][index_y][index_z]; // VORZEICHEN BEI VOXEL_Z/E_Z positiv oder negativ?

        if(arr_E_z[index_x][index_y][index_z]!=0) hit_y = hit_y + voxel_size_z/(arr_E_z[index_x][index_y][index_z])*arr_E_y[index_x][index_y][index_z];
        //if(arr_E_z[index_x][index_y][index_z]==0) hit_x = hit_y + arr_E_y[index_x][index_y][index_z]; // VORZEICHEN BEI VOXEL_Z/E_Z positiv oder negativ?

        actual_z_position = actual_z_position + voxel_size_z;
        std::cout << " Actual z: " << actual_z_position << std::endl;


        // Propagate the old coordinates to the next z slice
        



        // Get the new bin number:
        std::cout << " Old bin indices: x: " << index_x << " \t y: " << index_y << " \t z: " << index_z << std::endl;

        for(int bin_index=0; bin_index<n_bins_x; bin_index++) {
            if( hit_x<x_min ) index_x = -1;
            if( hit_x>=(x_min+bin_index*(x_max-x_min)/n_bins_x) && hit_x<(x_min+(bin_index+1)*(x_max-x_min)/n_bins_x) ) {
                index_x = bin_index;
            }
            if( hit_x>=x_max ) index_x = n_bins_x;
        }
        for(int bin_index=0; bin_index<n_bins_y; bin_index++) {
            if( hit_y<y_min ) index_y = -1;
            if( hit_y>=(y_min+bin_index*(y_max-y_min)/n_bins_y) && hit_y<(y_min+(bin_index+1)*(y_max-y_min)/n_bins_y) ) {
                index_y = bin_index;
            }
            if( hit_y>=y_max ) index_y = n_bins_y;
        }
        for(int bin_index=0; bin_index<n_bins_z; bin_index++) {
            if( actual_z_position<z_min ) index_z = -1;
            if( actual_z_position>=(z_min+bin_index*(z_max-z_min)/n_bins_z) && actual_z_position<(z_min+(bin_index+1)*(z_max-z_min)/n_bins_z) ) {
                //index_z = bin_index;
            }
            if( actual_z_position>=z_max) index_z = n_bins_z;
        }

        std::cout << " New bin indices: x: " << index_x << " \t y: " << index_y << " \t z: " << index_z << std::endl;

    }






    // Read the .root file with the electric field information
/*    auto f = new TFile("E_field.root");
    auto T = (TTree*)f->Get("Tree_Efield");
    //TH1F *hpx = nullptr;
    //TH2F *hpxpy = nullptr;
    //TProfile *hprof = nullptr;
    TH3D * hist_E_x = nullptr;
	TH3D * hist_E_y = nullptr;
	TH3D * hist_E_z = nullptr;
    
    T->SetBranchAddress("hist_E_x",&hist_E_x);
    T->SetBranchAddress("hist_E_y",&hist_E_y);
    T->SetBranchAddress("hist_E_z",&hist_E_z);
    T->GetEntry(10);
    std::cout << " Number of entries in tree: " << T->GetEntries() << std::endl;
    auto c1 = new TCanvas("c1","test",10,10,600,1000);
    c1->Divide(1,3);
    c1->cd(1);
    hist_E_x->Draw();
    c1->cd(2);
    hist_E_y->Draw();
    c1->cd(3);
    hist_E_z->Draw();
    //c1->Print("htr1.png");
*/
    // Actual propagation of hit (at readout plane) along the electric field lines
    
    
    
    
    
    
    
    position_corrected[0] = hit_x * wire_pitch;
    
    // Start at coordinates (hit_x, hit_y, VOXEL_Z_SIZE/2)    


	//return hit_position;
	
	return 0;
}
