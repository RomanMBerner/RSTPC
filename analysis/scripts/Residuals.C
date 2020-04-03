#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include <iostream>
#include <algorithm>
//#include "math.h"



// **************************************************************
std::vector<double> calculate_distance(float hit_x, float hit_y, float hit_z, std::vector<double> principal_component, std::vector<double> barycentre) {
// **************************************************************

	// Initialize the vector (length = number of hits, which is 1) of vectors (with length = 3 for the residuals in each spatial direction)
	std::vector<double> residuals(3,-999.);


	// For each hit, calculate the shortest vector from the hit to the principal component's line
	double lambda = ( ( (hit_x-barycentre[0])*principal_component[0] +
						(hit_y-barycentre[1])*principal_component[1] +
						(hit_z-barycentre[2])*principal_component[2] ) / (sqrt( principal_component[0]*principal_component[0] + principal_component[1]*principal_component[1] + principal_component[2]*principal_component[2] )) );

	residuals[0] = barycentre[0] + lambda * principal_component[0] - hit_x;
	residuals[1] = barycentre[1] + lambda * principal_component[1] - hit_y;
	residuals[2] = barycentre[2] + lambda * principal_component[2] - hit_z;

	//std::cout << " Residuals: " << residuals[0] << " \t" << residuals[1] << " \t" << residuals[2] << std::endl;

	return residuals;
}



// **************************************************************
std::vector<std::vector<double>> calculate_residuals(std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z, std::vector<double> principal_component, std::vector<double> barycentre) {
// **************************************************************

	// Initialize a vector (length = number of hits) of vectors (with length = 3 for the residuals in each spatial direction)
	std::vector<std::vector<double>> residuals(hits_x.size(),std::vector<double>(3,-1.));


    // First, check that all input vectors have the same size
    //double error[hit_x.size()][3];
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() ) {
        std::cout << " ERROR IN FUNCTION 'calculate_residuals': size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return residuals;
    }


	// For each hit, calculate the shortest vector from the hit to the principal component's line
	double lambda;

	for(int hit=0; hit<hits_x.size(); hit++) {
		lambda = ( ( 	(hits_x[hit]-barycentre[0])*principal_component[0] +
						(hits_y[hit]-barycentre[1])*principal_component[1] +
						(hits_z[hit]-barycentre[2])*principal_component[2] ) / (sqrt( principal_component[0]*principal_component[0] + principal_component[1]*principal_component[1] + principal_component[2]*principal_component[2] )) );

		residuals[hit][0] = barycentre[0] + lambda * principal_component[0] - hits_x[hit];
		residuals[hit][1] = barycentre[1] + lambda * principal_component[1] - hits_y[hit];
		residuals[hit][2] = barycentre[2] + lambda * principal_component[2] - hits_z[hit];
	}

	return residuals;
}


// **************************************************************
std::vector<double> calculate_eta(std::vector<float> hits_z, std::vector<std::vector<double>> residuals) {
// **************************************************************

    // eta is a linear parameter which is defined from 0 to 1.
    // Project all hits onto the principal component axis.
    // The hit with the smallest z has eta=0, the hit with the largest z has eta=1.

    std::vector<double> eta(hits_z.size(),-1.);

    // Calculate the z coordinate of each hit projected onto the principal components axis
    // Note: If residual is positive (negative), then the hit is below (above) the principal component.
    for(int hit=0; hit<hits_z.size(); hit++) {
        eta[hit] = hits_z[hit] + residuals[hit][2];
    }


    // Normalize eta such that eta_min=0 and eta_max=1
    double eta_min = *std::min_element(eta.begin(),eta.end());
    double eta_max = *std::max_element(eta.begin(),eta.end());

    for(int hit=0; hit<hits_z.size(); hit++) {
        eta[hit] = (eta[hit]-eta_min) / (eta_max-eta_min);
    }
    
	return eta;
}



// **************************************************************
std::vector<double> calculate_l(std::vector<float> hits_x, std::vector<float> hits_y) {
// **************************************************************

    // l is a linear parameter which is >= 0.
    // Define an axis through the middle of the TPC active volume (in drift direction).
    // l is the smallest distance of a hit to this axis.
    std::vector<double> l(hits_x.size(),-1.);


    // First, check that the input vectors have the same size
    //double error[hit_x.size()][3];
    if( hits_x.size() != hits_y.size() ) {
        std::cout << " ERROR IN FUNCTION 'calculate_l': size of hits_x and hits_y do not match !!" << std::endl;
        return l;
    }


    // Calculate the distance of each hit to l, which has (x,y) coordinates of:
    // hits given in units of [wire pitches] -> calculate l in units of [wire pitches]
    for(int hit=0; hit<hits_x.size(); hit++) {
        l[hit] = sqrt( pow( fabs(31./2.-hits_x[hit]),2 ) + pow( fabs(31./2.-hits_y[hit]),2 ) );
    }

	return l;
}


// **************************************************************
std::vector<int> calculate_region(std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z, std::vector<double> principal_component, std::vector<double> barycentre) {
// **************************************************************

	// Initialize a vector (length = number of hits) of vectors (with length = 3 for the residuals in each spatial direction)
	std::vector<int> region(hits_x.size());


    // First, check that all input vectors have the same size
    //double error[hit_x.size()][3];
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() ) {
        std::cout << " ERROR IN FUNCTION 'calculate_region': size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return region;
    }


	// For each hit, calculate the region it belongs to (the regions are shown in file 'regions.png')
    // The region is chosen as the point on the PC nearest to the hit
	
	double lambda;
	double x,y;
	
	for(int hit=0; hit<hits_x.size(); hit++) {
		lambda = ( ( (hits_x[hit]-barycentre[0])*principal_component[0] +
					 (hits_y[hit]-barycentre[1])*principal_component[1] +
					 (hits_z[hit]-barycentre[2])*principal_component[2] ) / (sqrt( principal_component[0]*principal_component[0] + principal_component[1]*principal_component[1] + principal_component[2]*principal_component[2] )) );
		
		x = barycentre[0] + lambda * principal_component[0];
		y = barycentre[1] + lambda * principal_component[1];
		z = barycentre[2] + lambda * principal_component[2];

        /*
	    if(x>= 0 && x< 6) {
	        if(y>= 0 && y< 6) { region[hit] = 0; }
	        if(y>= 6 && y<12) { region[hit] = 1; }
	        if(y>=12 && y<20) { region[hit] = 2; }
	        if(y>=20 && y<26) { region[hit] = 3; }
	        if(y>=26 && y<32) { region[hit] = 4; }
	    }
	    if(x>= 6 && x<12) {
	        if(y>= 0 && y< 6) { region[hit] = 5; }
	        if(y>= 6 && y<12) { region[hit] = 6; }
	        if(y>=12 && y<20) { region[hit] = 7; }
	        if(y>=20 && y<26) { region[hit] = 8; }
	        if(y>=26 && y<32) { region[hit] = 9; }
	    }
	    if(x>=12 && x<20) {
	        if(y>= 0 && y< 6) { region[hit] = 10; }
	        if(y>= 6 && y<12) { region[hit] = 11; }
	        if(y>=12 && y<20) { region[hit] = 12; }
	        if(y>=20 && y<26) { region[hit] = 13; }
	        if(y>=26 && y<32) { region[hit] = 14; }
	    }
	    if(x>=20 && x<26) {
	        if(y>= 0 && y< 6) { region[hit] = 15; }
	        if(y>= 6 && y<12) { region[hit] = 16; }
	        if(y>=12 && y<20) { region[hit] = 17; }
	        if(y>=20 && y<26) { region[hit] = 18; }
	        if(y>=26 && y<32) { region[hit] = 19; }
	    }
	    if(x>=26 && x<32) {
	        if(y>= 0 && y< 6) { region[hit] = 20; }
	        if(y>= 6 && y<12) { region[hit] = 21; }
	        if(y>=12 && y<20) { region[hit] = 22; }
	        if(y>=20 && y<26) { region[hit] = 23; }
	        if(y>=26 && y<32) { region[hit] = 24; }
	    }
	    */

	    if(x>= 0 && x< 5) {
	        if(y>= 0 && y< 5) { region[hit] = 0; }
	        if(y>= 5 && y<12) { region[hit] = 1; }
	        if(y>=12 && y<20) { region[hit] = 2; }
	        if(y>=20 && y<27) { region[hit] = 3; }
	        if(y>=27 && y<32) { region[hit] = 4; }
	    }
	    if(x>= 5 && x<12) {
	        if(y>= 0 && y< 5) { region[hit] = 5; }
	        if(y>= 5 && y<12) { region[hit] = 6; }
	        if(y>=12 && y<20) { region[hit] = 7; }
	        if(y>=20 && y<27) { region[hit] = 8; }
	        if(y>=27 && y<32) { region[hit] = 9; }
	    }
	    if(x>=12 && x<20) {
	        if(y>= 0 && y< 5) { region[hit] = 10; }
	        if(y>= 5 && y<12) { region[hit] = 11; }
	        if(y>=12 && y<20) { region[hit] = 12; }
	        if(y>=20 && y<27) { region[hit] = 13; }
	        if(y>=27 && y<32) { region[hit] = 14; }
	    }
	    if(x>=20 && x<27) {
	        if(y>= 0 && y< 5) { region[hit] = 15; }
	        if(y>= 5 && y<12) { region[hit] = 16; }
	        if(y>=12 && y<20) { region[hit] = 17; }
	        if(y>=20 && y<27) { region[hit] = 18; }
	        if(y>=27 && y<32) { region[hit] = 19; }
	    }
	    if(x>=27 && x<32) {
	        if(y>= 0 && y< 5) { region[hit] = 20; }
	        if(y>= 5 && y<12) { region[hit] = 21; }
	        if(y>=12 && y<20) { region[hit] = 22; }
	        if(y>=20 && y<27) { region[hit] = 23; }
	        if(y>=27 && y<32) { region[hit] = 24; }
	    }
	}

	return region;
}


// **************************************************************
std::vector<std::vector<double>> calculate_theta_method1(std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z) {
// **************************************************************

	// Initialize a vector (length = number of hits) of vectors (with length = number of hits)
	std::vector<std::vector<double>> theta(hits_x.size(),std::vector<double>(hits_x.size(),-111.));


    // First, check that all input vectors have the same size and at least three hits are available
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() || hits_x.size()<3 ) {
        std::cout << " ERROR IN FUNCTION 'calculate_theta': size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return theta;
    }


	// For each combination of four (three in the fist and last cases) hits, calculate the angle between the direction vectors (those vectors span a 2D plane)
	// ANSATZ: vec{a}*vec{b} = abs{vec{a}}*abs{vec{b}}*cos(theta)
	std::vector<double> a(3,-999.);
	std::vector<double> b(3,-999.);
	//std::cout << " Number of hits: " << hits_x.size() << std::endl;
    //for(int hit=0; hit<hits_x.size(); hit++) {
        //std::cout << " hit " << hit << "  \tx: " << hits_x[hit] << "  \ty: " << hits_y[hit] << "  \tz: " << hits_z[hit] << std::endl;
    //}
	for(int hit1=0; hit1<(hits_x.size()-2); hit1++) {
	    for(int hit2=hit1+1; hit2<hits_x.size()-1; hit2++) {
	        a[0] = hits_x[hit1+1]-hits_x[hit1]; // vectors pointing from hit at t_i to hit with t_(i+1)
	        a[1] = hits_y[hit1+1]-hits_y[hit1];
	        a[2] = hits_z[hit1+1]-hits_z[hit1];
	        b[0] = hits_x[hit2+1]-hits_x[hit2];
	        b[1] = hits_y[hit2+1]-hits_y[hit2];
	        b[2] = hits_z[hit2+1]-hits_z[hit2];
    		theta[hit1][hit2-hit1-1] = acos( (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]) / ( sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2)) * sqrt(pow(b[0],2)+pow(b[1],2)+pow(b[2],2)) ) ) * 360./(2*M_PI); // in degrees
            //std::cout << " " << hit1 << "-->" << hit1+1 << "  \t " << hit2 << "-->" << hit2+1 << "  \tTheta [deg]: " << theta[hit1][hit2-hit1-1] << std::endl;
        }
	}
	
	return theta;
}



// **************************************************************
std::vector<std::vector<double>> calculate_length_method1(std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z) {
// **************************************************************

    // Length of the track: From second point to last point

    // Initialize a vector (length = number of hits) of vectors (with length = number of hits)
    std::vector<std::vector<double>> length(hits_x.size(),std::vector<double>(hits_x.size(),-111.));


    // First, check that all input vectors have the same size and at least three hits are available
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() || hits_x.size()<3  ) {
        std::cout << " ERROR IN FUNCTION 'calculate_length': size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return length;
    }


    // For each combination of four (three in the fist and last cases) hits, calculate the length of the track (distance traversed)
    //std::vector<double> a(3,-999.);
    //std::vector<double> b(3,-999.);
    //std::cout << " Number of hits: " << hits_x.size() << std::endl;
    //for(int hit=0; hit<hits_x.size(); hit++) {
    //    std::cout << " hit " << hit << "  \tx: " << hits_x[hit] << "  \ty: " << hits_y[hit] << "  \tz: " << hits_z[hit] << std::endl;
    //}
    double wire_pitch = 52.5/31; // [mm]
    for(int hit1=0; hit1<(hits_x.size()-2); hit1++) {
        for(int hit2=hit1+1; hit2<hits_x.size()-1; hit2++) {
            length[hit1][hit2-hit1-1] = sqrt( pow((hits_x[hit1+1]-hits_x[hit2+1])*wire_pitch,2) + pow((hits_y[hit1+1]-hits_y[hit2+1])*wire_pitch,2) + pow(hits_z[hit1+1]-hits_z[hit2+1],2) ); // [mm] 
            //std::cout << " " << hit1 << "-->" << hit1+1 << "  \t " << hit2 << "-->" << hit2+1 << "  \tlength [mm]: " << length[hit1][hit2-hit1-1] << std::endl;
        }
    }

    return length;
}


// **************************************************************
std::vector<std::vector<double>> calculate_theta_method2(std::vector<double> principal_component, std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z) {
// **************************************************************

	// Initialize a vector (length = number of hits) of vectors (with length = number of hits)
	std::vector<std::vector<double>> theta(hits_x.size(),std::vector<double>(hits_x.size(),-111.));


    // First, check that all input vectors have the same size and at least three hits are available
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() || hits_x.size()<3 ) {
        std::cout << " ERROR IN FUNCTION 'calculate_theta_method2': size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return theta;
    }


	// For each combination of four (three in the fist and last cases) hits, calculate the angle between the direction vectors (those vectors span a 2D plane)
	// ANSATZ: vec{a}*vec{b} = abs{vec{a}}*abs{vec{b}}*cos(theta)
	std::vector<double> b(3,-999.);
	//std::cout << " Number of hits: " << hits_x.size() << std::endl;
    //for(int hit=0; hit<hits_x.size(); hit++) {
        //std::cout << " hit " << hit << "  \tx: " << hits_x[hit] << "  \ty: " << hits_y[hit] << "  \tz: " << hits_z[hit] << std::endl;
    //}
    int distance = 1; // number of hits to go further (if 1, we're probably limited by the lateral spread of the hit distribution)
	for(int hit1=(4-1); hit1<(hits_x.size()-2); hit1++) {
	    for(int hit2=hit1; hit2<(hits_x.size()-distance); hit2++) {
	        b[0] = -hits_x[hit1+distance]+hits_x[hit1]; // vectors pointing from hit at t_i to hit with t_(i+1)
	        b[1] = -hits_y[hit1+distance]+hits_y[hit1];
	        b[2] = -hits_z[hit1+distance]+hits_z[hit1];
    		theta[hit1-(4-1)][hit2-hit1] = acos( (principal_component[0]*b[0]+principal_component[1]*b[1]+principal_component[2]*b[2]) / ( sqrt(pow(principal_component[0],2)+pow(principal_component[1],2)+pow(principal_component[2],2)) * sqrt(pow(b[0],2)+pow(b[1],2)+pow(b[2],2)) ) ) * 360./(2*M_PI); // in degrees
            //std::cout << " PC of first 4 hits \t " << hit2 << "-->" << hit2+1 << " \tb[0]: " << b[0] << " \tb[1]: " << b[1] << " \tb[2]: " << b[2] << " \tTheta [deg]: " << theta[hit1-(4-1)][hit2-hit1] << std::endl;
        }
	}
	
	return theta;
}



// **************************************************************
std::vector<std::vector<double>> calculate_length_method2(std::vector<double> principal_component, std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z) {
// **************************************************************

    // Length of the track: From second point to last point

    // Initialize a vector (length = number of hits) of vectors (with length = number of hits)
    std::vector<std::vector<double>> length(hits_x.size(),std::vector<double>(hits_x.size(),-111.));


    // First, check that all input vectors have the same size and at least three hits are available
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() || hits_x.size()<3  ) {
        std::cout << " ERROR IN FUNCTION 'calculate_length_method2': size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return length;
    }


    // For each combination of four (three in the fist and last cases) hits, calculate the length of the track (distance traversed)
    //std::vector<double> a(3,-999.);
    //std::vector<double> b(3,-999.);
    //std::cout << " Number of hits: " << hits_x.size() << std::endl;
    //for(int hit=0; hit<hits_x.size(); hit++) {
    //    std::cout << " hit " << hit << "  \tx: " << hits_x[hit] << "  \ty: " << hits_y[hit] << "  \tz: " << hits_z[hit] << std::endl;
    //}
    double wire_pitch = 52.5/31; // [mm]
    int distance = 1; // number of hits to go further (if 1, we're probably limited by the lateral spread of the hit distribution)
    for(int hit1=(4-1); hit1<(hits_x.size()-2); hit1++) {
        for(int hit2=hit1; hit2<(hits_x.size()-1); hit2++) {
            length[hit1-(4-1)][hit2-hit1] = sqrt( pow((hits_x[hit1]-hits_x[hit2+distance])*wire_pitch,2) + pow((hits_y[hit1]-hits_y[hit2+distance])*wire_pitch,2) + pow(hits_z[hit1]-hits_z[hit2+distance],2) ); // [mm] 
            //std::cout << " 4th hit \t " << hit2 << "-->" << hit2+1 << "  \tlength [mm]: " << length[hit1-(4-1)][hit2-hit1] << std::endl;
        }
    }

    return length;
}


// **************************************************************
std::vector<std::vector<double>> calculate_theta_method3(std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z) {
// **************************************************************

	// Initialize a vector (length = number of hits) of vectors (with length = number of hits)
	std::vector<std::vector<double>> theta(hits_x.size(),std::vector<double>(hits_x.size(),-111.));


    // First, check that all input vectors have the same size and at least three hits are available
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() || hits_x.size()<3 ) {
        std::cout << " ERROR IN FUNCTION 'calculate_theta': size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return theta;
    }


	// For each combination of four (three in the fist and last cases) hits, calculate the angle between the direction vectors (those vectors span a 2D plane)
	// ANSATZ: vec{a}*vec{b} = abs{vec{a}}*abs{vec{b}}*cos(theta)
	std::vector<double> a(3,-999.);
	std::vector<double> b(3,-999.);
	//std::cout << " Number of hits: " << hits_x.size() << std::endl;
    //for(int hit=0; hit<hits_x.size(); hit++) {
        //std::cout << " hit " << hit << "  \tx: " << hits_x[hit] << "  \ty: " << hits_y[hit] << "  \tz: " << hits_z[hit] << std::endl;
    //}
	for(int hit1=0; hit1<(hits_x.size()-2); hit1++) {
	    for(int hit2=hit1+1; hit2<hits_x.size()-1; hit2++) {
	        a[0] = hits_x[hit1+1]-hits_x[hit1]; // vectors pointing from hit at t_i to hit with t_(i+1)
	        a[1] = hits_y[hit1+1]-hits_y[hit1];
	        a[2] = hits_z[hit1+1]-hits_z[hit1];
	        b[0] = hits_x[hit2+1]-hits_x[hit2];
	        b[1] = hits_y[hit2+1]-hits_y[hit2];
	        b[2] = hits_z[hit2+1]-hits_z[hit2];
    		theta[hit1][hit2-hit1-1] = acos( (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]) / ( sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2)) * sqrt(pow(b[0],2)+pow(b[1],2)+pow(b[2],2)) ) ) * 360./(2*M_PI); // in degrees
            //std::cout << " " << hit1 << "-->" << hit1+1 << "  \t " << hit2 << "-->" << hit2+1 << "  \tTheta [deg]: " << theta[hit1][hit2-hit1-1] << std::endl;
        }
	}
	
	return theta;
}



// **************************************************************
std::vector<std::vector<double>> calculate_length_method3(std::vector<float> hits_x, std::vector<float> hits_y, std::vector<float> hits_z) {
// **************************************************************

    // Length of the track: From second point to last point

    // Initialize a vector (length = number of hits) of vectors (with length = number of hits)
    std::vector<std::vector<double>> length(hits_x.size(),std::vector<double>(hits_x.size(),-111.));


    // First, check that all input vectors have the same size and at least three hits are available
    if( hits_x.size() != hits_y.size() || hits_x.size() != hits_z.size() || hits_y.size() != hits_z.size() || hits_x.size()<3  ) {
        std::cout << " ERROR IN FUNCTION 'calculate_length': size of hits_x, hits_y, hits_z do not match !!" << std::endl;
        return length;
    }


    // For each combination of four (three in the fist and last cases) hits, calculate the length of the track (distance traversed)
    //std::vector<double> a(3,-999.);
    //std::vector<double> b(3,-999.);
    //std::cout << " Number of hits: " << hits_x.size() << std::endl;
    //for(int hit=0; hit<hits_x.size(); hit++) {
    //    std::cout << " hit " << hit << "  \tx: " << hits_x[hit] << "  \ty: " << hits_y[hit] << "  \tz: " << hits_z[hit] << std::endl;
    //}
    double wire_pitch = 52.5/31; // [mm]
    for(int hit1=0; hit1<(hits_x.size()-2); hit1++) {
        for(int hit2=hit1+1; hit2<hits_x.size()-1; hit2++) {
            length[hit1][hit2-hit1-1] = sqrt( pow((hits_x[hit1+1]-hits_x[hit2+1])*wire_pitch,2) + pow((hits_y[hit1+1]-hits_y[hit2+1])*wire_pitch,2) + pow(hits_z[hit1+1]-hits_z[hit2+1],2) ); // [mm] 
            //std::cout << " " << hit1 << "-->" << hit1+1 << "  \t " << hit2 << "-->" << hit2+1 << "  \tlength [mm]: " << length[hit1][hit2-hit1-1] << std::endl;
        }
    }

    return length;
}
