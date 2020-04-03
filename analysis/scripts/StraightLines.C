#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TObjArray.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <TTree.h>
#include <TFile.h>
#include <vector>

#include "HistoManipulators.hh"
#include "DigitalFilters.hh"
#include "RSTPC_RunProcessor.hh"
#include "RSTPC_Analyser.hh"
#include "MppcTreeWrapper.hh"
#include "RSTPC_T1wrapper.hh"
#include "RSTPC_T2wrapper.hh"
#include "RSTPC_Hits.hh"
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TClassTable.h"

//#include "RSTPC_T1wrapper.hh"
//#include "RSTPC_T2wrapper.hh"
//#include "TestPulses.C"
//#include "../RSTPC_T1wrapper.hh" // for the RSTPC_T1wrapper class; DOESN'T WORK
//#include "../RSTPC_Hits.hh" // for the RSTPC_Hit and RSTPC_Pulse classes

//#include "PCA.C"



// ***************************************************
void StraightLines() {
// ***************************************************

    //gROOT->ProcessLine(".x LoadLibs.C");
    //gROOT->ProcessLine(".L RSTPC_T2wrapper.cc");



    // Define data files
    // ===================================================================
    TFile * tfile_run_2032 = new TFile("/home/francescop/data/ResistiveShell/merged/test/RSTPC_Run000002032_Merged.root");
                                     // /home/rberner   /data/ResistiveShell/merged/test/RSTPC_Run000002032_Merged.root");
                                     
    int run_number = 2032;


	// Produce some plots for initial data inspection
	// (-> which ADC cut values should be chosen, etc.)
	// Functions defined in TestPulses.C
	// ===================================================================
    //PlotColPulses(); // Amplitudes of all ColPulses from the entire run plotted vs. width (= LeftEdge - RightEdge)
    //Plot_ColPulses_Ampl_vs_Width(); // Amplitudes of each single Col Wire plotted vs. width
    //Plot_ColPulses_Ampl_vs_FWHM(); // Amplitudes of each single Col Wire plotted vs. FWHM
    //Plot_ColPulses_Ampl_vs_FWTM(); // Amplitudes of each single Col Wire plotted vs. FWTM
    //Plot_ColPulses_Ampl_vs_sigma(); // Amplitudes of each single Col Wire plotted vs. sigma
    //Plot_IndPulses_Ampl_vs_sigma(); // Amplitudes of each single Ind Wire plotted vs. sigma



    // Access the T1 data
    // ===================================================================
    // Get T1 tree and link variables
    TTree * T1 = (TTree*)tfile_run_2032->Get("T1");

    Int_t       TpcEvent;
    Double_t    TpcTime;
    Double_t    RmsColWires[32];
    Double_t    RmsIndWires[32];
    Int_t       FebEvent;
    Double_t    FebTime;
    UShort_t    FebTopAmp[3];
    Double_t    FebTopTotAmp;
    UShort_t    FebBotAmp[3];
    Double_t    FebBotTotAmp;
    Double_t    FebTotAmp;

    T1->SetBranchAddress("TpcEvent",     &TpcEvent);
    T1->SetBranchAddress("TpcTime",      &TpcTime);
    T1->SetBranchAddress("RmsColWires",  &RmsColWires);
    T1->SetBranchAddress("RmsIndWires",  &RmsIndWires);
    T1->SetBranchAddress("FebEvent",     &FebEvent);
    T1->SetBranchAddress("FebTime",      &FebTime);
    T1->SetBranchAddress("FebTopAmp",    &FebTopAmp);
    T1->SetBranchAddress("FebTopTotAmp", &FebTopTotAmp);
    T1->SetBranchAddress("FebBotAmp",    &FebBotAmp);
    T1->SetBranchAddress("FebBotTotAmp", &FebBotTotAmp);
    T1->SetBranchAddress("FebTotAmp",    &FebTotAmp);

    /*
    std::cout << " ---------------------------------------- "            << std::endl;
    std::cout << " Reading T1 with " << T1->GetEntries() << " entries: " << std::endl;
    std::cout << " ---------------------------------------- "            << std::endl;
    //std::cout << " Tpc \t Tpc \t\t Rms \t Rms \t\t Feb \t Feb \t\t FebT \t FebT \t FebB \t FebB \t Feb " << std::endl;
    //std::cout << " Evnt: \t Time: \t\t ColW: \t IndW: \t\t Evnt: \t Time: \t\t Amp: \t TotAmp: Amp: \t TotAmp: TotAmp: " << std::endl << std::endl;
    for(int entry=2; entry<3; entry++) {
        T1->GetEntry(entry);
      	std::cout << " " << TpcEvent       << " \t " << TpcTime      << " \t "     << RmsColWires[0] << " \t "
                         << RmsIndWires[5] << " \t " << FebEvent     << " \t "     << FebTime        << " \t "
                         << FebTopAmp[2]   << " \t " << FebTopTotAmp << " \t "     << FebBotAmp[0]   << " \t "
                         << FebBotTotAmp   << " \t " << FebTotAmp    << std::endl;
        //T1->Show(entry);
    }
    std::cout << std::endl;
    */



    // Access the T2 data
    // ===================================================================
    // Get T2 tree and link variables
    TTree * T2 = (TTree*)tfile_run_2032->Get("T2");
    TTree * testtree = (TTree*)tfile_run_2032->Get("T2");

    TClonesArray * tclonesarray_T2_ColPulses = new TClonesArray("RSTPC_Pulse");
    TClonesArray * tclonesarray_T2_IndPulses = new TClonesArray("RSTPC_Pulse");
    TClonesArray * tclonesarray_T2_Hits 	 = new TClonesArray("RSTPC_Hit");

    //TClonesArray &ColPulses_array = *tclonesarray_T2_ColPulses;
    //TClonesArray &IndPulses_array = *tclonesarray_T2_IndPulses;

    T2->SetBranchAddress("ColPulses", &tclonesarray_T2_ColPulses);
    T2->SetBranchAddress("IndPulses", &tclonesarray_T2_IndPulses);
    T2->SetBranchAddress("Hits", 	  &tclonesarray_T2_Hits);

    std::cout << " ---------------------------------------- "            << std::endl;
    std::cout << " Reading T2 with " << T2->GetEntries() << " entries: " << std::endl;
    std::cout << " ---------------------------------------- "            << std::endl;



    // Check that T1 as well as T2 have the same number of entries
    if( T1->GetEntries() != T2->GetEntries() ) {
        std::cout << " WARNING: T1 AND T2 HAVE DIFFERENT NUMBER OF ENTRIES !! " << std::endl;
    }

    
    int n_hits = 0;
    int n_good_events = 0;
    double average_hits_per_good_event = 0;


    // Clear file with 3Dhits and residuals & file with angles and lengths
    ofstream clear_file;
    clear_file.open("Hits_Residuals_Eta_l_region.txt");
    if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }
    clear_file.open("Angles_and_lengths.txt");
    if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }


    // Create 1D histogram to show the maximal drift distances (should not be larger than 150 mm, otherwise change e_drift_velocity)
    double e_drift_velocity = 0.212; // cm/us @ 16.4/15 kV/cm for Run2032
    char * hist_name = new char[60];
    sprintf(hist_name,"drift distance first - last hit, at %.4f kV/cm",e_drift_velocity);
    TH1F * max_drift_distances = new TH1F(hist_name,hist_name,50,0,200);

    
    // Create 1D histogram to show the number of hits per event
    int minimum_hits = 6; // will be used as cut
    TH1I * n_hits_per_event = new TH1I("n_hits_per_event","hits per event",25,0,25);
    
    
    // Create 1D histogram to show the mean distance between the hits (distance here is the distance between the hit's projections to the PC)
    TH1F * distance_between_hits = new TH1F("distance_between_hits","distance between hits (on PC)",50,0,100);


    // Create 1D histogram to show the z distribution of the hits
    TH1F * z_distribution_of_hits = new TH1F("z_distribution_of_hits","z_distribution_of_hits",200,-20,180);


        // Loop over all events in T2
        // -------------------------------------------------------------------
        for(int event=0; event<T2->GetEntries(); event++) { //T2->GetEntries(); event++) { // 2: reference track // 4: muon + delta // 7: muon only // 931: long muon track
            T2->GetEntry(event);
            //T2->Show(event);
        //std::cout << "======> EVENT:" << event << std::endl;



        // Get the number of ColPulses, IndPulses and Hits in this event
        UInt_t NColPulses = tclonesarray_T2_ColPulses->GetEntries();
        UInt_t NIndPulses = tclonesarray_T2_IndPulses->GetEntries();
        UInt_t NHits = tclonesarray_T2_Hits->GetEntries();
        //std::cout << " NColPulses: " << NColPulses << " \tNIndPulses: " << NIndPulses << " \tNHits: " << NHits << std::endl;



        // Loop over all ColPulses and put those with a small ADC amplitude in a set 'badColPulseIDs' (containing bad ColPulseIDs)
        Double_t ColPulse_ADC_threshold = 100.;
        std::set<ULong_t> badColPulseIDs;

        for(UInt_t pulse=0; pulse<NColPulses; pulse++) {
            RSTPC_Pulse * t2event_ColPulse = (RSTPC_Pulse *)tclonesarray_T2_ColPulses->At(pulse);
            if( t2event_ColPulse->fMax<ColPulse_ADC_threshold ) {
                badColPulseIDs.insert( t2event_ColPulse->fPulseID );
            }
            /*else {
                std::cout << "  fColPulseID: "      << t2event_ColPulse->fPulseID
                          << "  \tfColCoinNum: "    << t2event_ColPulse->fColCoinNum
                          << "  \tfMax: "           << t2event_ColPulse->fMax
                          << "  \twidth_ColPulse: " << t2event_ColPulse->fRedge - t2event_ColPulse->fLedge;
                          << "  \tfSigma: "         << t2event_ColPulse->fSigma
                          << "  \tfLedge: "         << t2event_ColPulse->fLedge
                          << "  \tfRedge: "         << t2event_ColPulse->fRedge      << std::endl;
            }*/
        }
        //std::cout << std::endl;


        // After the threshold cut, only a few hits should be remaining.
        // Loop over all remaining (those have high ADC value) ColPulses and put those ...
                // ... which are isolated (no neighboring hit widhin 4 wire pitches in each x,y and z direction) in the set badColPulseIDs
                // ... occurring at the same time (within +/- 1 us) in a set 'CoincidentColPulseIDs'
        // First, have to produce vector-map with the ColPulses (the index corresponds to the ColPulseID)
        vector<RSTPC_Pulse*> ColPulsesVec(NColPulses);
        map<UInt_t, RSTPC_Pulse*> ColPulsesMap;
        for(UInt_t pulse=0; pulse<NColPulses; pulse++) {
            RSTPC_Pulse * t2event_ColPulse = (RSTPC_Pulse *)tclonesarray_T2_ColPulses->At(pulse);
            ColPulsesVec.at(pulse) = t2event_ColPulse;
            ColPulsesMap[t2event_ColPulse->fPulseID] = t2event_ColPulse;
        }



        // ============================================================================================================================
        // THIS CUT SHOULD BE AFTER THE E-FIELD CORRECTION
        // ============================================================================================================================
        // Put all isolated hits (those which do not have a neighboring hit within the distance of 4 wire pitches in each x,y and z direction) in the set badColPulseIDs
        // A much nicer algorithm would use the convex envelope which has to be smaller than a certain threshold
        double min_distance_to_other_hits = 999.9;
        for(int hit1=0; hit1<NHits; hit1++) {
            RSTPC_Hit * t2event_Hit1 = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit1);
            RSTPC_Pulse* CollPulse1 = ColPulsesMap[t2event_Hit1->fColPulseID];
            if(CollPulse1->fMax<ColPulse_ADC_threshold) {
                continue;
            }
            for(int hit2=hit1; hit2<NHits; hit2++) {
                if(hit1==hit2) {
                    continue;
                }
                RSTPC_Hit * t2event_Hit2 = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit2);
                RSTPC_Pulse* CollPulse2 = ColPulsesMap[t2event_Hit2->fColPulseID];
                if(CollPulse2->fMax<ColPulse_ADC_threshold) {
                    continue;
                }
                double distance_2D = sqrt( pow( ((double)t2event_Hit1->fColWireNum - (double)t2event_Hit2->fColWireNum),2 ) +
                                           pow( ((double)t2event_Hit1->fIndWireNum - (double)t2event_Hit2->fIndWireNum),2 ) );
                //double distance_3D = sqrt( pow( ((double)t2event_Hit1->fColWireNum - (double)t2event_Hit2->fColWireNum),2 ) +
                //                           pow( ((double)t2event_Hit1->fIndWireNum - (double)t2event_Hit2->fIndWireNum),2 ) +
                //                           pow( (t2event_Hit1->fCentreTime - t2event_Hit2->fCentreTime)/20*e_drift_velocity*10*31/52.5,2) );

                if(distance_2D < min_distance_to_other_hits) { min_distance_to_other_hits = distance_2D; }
            }            
            if(min_distance_to_other_hits>8) {
                badColPulseIDs.insert( CollPulse1->fPulseID ); // see e.g event 177
            }
        }
        // ===========================================================================================================================

        
        // Put the coincident ColPulses in the set 'CoincidentColPulseID'
        std::set<ULong_t> CoincidentColPulseIDs;
        for(int hit1=0; hit1<NHits; hit1++) {
            RSTPC_Hit * t2event_Hit1 = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit1);
            RSTPC_Pulse* CollPulse1 = ColPulsesMap[t2event_Hit1->fColPulseID];
            if( ! (badColPulseIDs.find(t2event_Hit1->fColPulseID) == badColPulseIDs.end() ) ) { // if true: Pulse is already in set 'badColPulseIDs'
                continue;
            }
            for(int hit2=0; hit2<NHits; hit2++ ) {
                if(hit1==hit2) continue;
                RSTPC_Hit * t2event_Hit2 = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit2);
                RSTPC_Pulse* CollPulse2 = ColPulsesMap[t2event_Hit2->fColPulseID];
                if( ! (badColPulseIDs.find(t2event_Hit2->fColPulseID) == badColPulseIDs.end() ) ) { // if true: Pulse is already in set 'badColPulseIDs'
                    continue;
                }
                if( abs(CollPulse1->fMaxPos-CollPulse2->fMaxPos)<50 ) { // 50 time samples = 50*50 ns = 2.5 us
                    CoincidentColPulseIDs.insert( t2event_Hit1->fColPulseID );
                    CoincidentColPulseIDs.insert( t2event_Hit2->fColPulseID );
                }
            }
        }


        // Loop over all IndPulses and put those with a small ADC amplitude in a set 'badIndPulseIDs' (containing bad IndPulseIDs)
        Double_t IndPulse_ADC_threshold = 150.;
        std::set<ULong_t> badIndPulseIDs;
        for(UInt_t pulse=0; pulse<NIndPulses; pulse++) {
            RSTPC_Pulse	* t2event_IndPulse	= (RSTPC_Pulse *)tclonesarray_T2_IndPulses->At(pulse);
            if( t2event_IndPulse->fMax<IndPulse_ADC_threshold ) {
                badIndPulseIDs.insert( t2event_IndPulse->fPulseID );
            }
            /*else {
                std::cout   << "  fIndPulseID: "      << t2event_IndPulse->fPulseID
                            << "  \tfIndCoinNum: "    << t2event_IndPulse->fIndCoinNum
                            << "  \tfMax: "           << t2event_IndPulse->fMax
                            << "  \twidth_IndPulse: " << t2event_IndPulse->fRedge - t2event_IndPulse->fLedge;
                            << "  \tfSigma: "         << t2event_IndPulse->fSigma
                            << "  \tfLedge: "         << t2event_IndPulse->fLedge
                            << "  \tfRedge: "         << t2event_IndPulse->fRedge      << std::endl;
            }*/
        }
        //std::cout << std::endl;


        // Loop over all hits and put the ColPulseIDs of those hits, which correspond to delta electrons (or other, selected by eye) in the set 'badHits'
        std::set<UInt_t> badHits;
        for(int hit=0; hit<NHits; hit++) {
            RSTPC_Hit * t2event_Hit = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit);
            if(event== 17 && t2event_Hit->fColWireNum== 2 && t2event_Hit->fIndWireNum==11) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 20 && t2event_Hit->fColWireNum== 9 && t2event_Hit->fIndWireNum==23) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 35 && t2event_Hit->fColWireNum== 5 && t2event_Hit->fIndWireNum==19) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 44 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==27) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 56 && t2event_Hit->fColWireNum== 8 && t2event_Hit->fIndWireNum== 1) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 75 && t2event_Hit->fColWireNum==10 && t2event_Hit->fIndWireNum== 9) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 75 && t2event_Hit->fColWireNum==11 && t2event_Hit->fIndWireNum== 9) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 75 && t2event_Hit->fColWireNum==12 && t2event_Hit->fIndWireNum== 9) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 75 && t2event_Hit->fColWireNum==15 && t2event_Hit->fIndWireNum== 5) { badHits.insert( t2event_Hit->fHitID ); }
            if(event== 75 && t2event_Hit->fColWireNum==17 && t2event_Hit->fIndWireNum== 4) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==114 && t2event_Hit->fColWireNum==15 && t2event_Hit->fIndWireNum== 2) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==171 && t2event_Hit->fColWireNum==21 && t2event_Hit->fIndWireNum== 6) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==177 && t2event_Hit->fColWireNum==31 && t2event_Hit->fIndWireNum== 9) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==177 && t2event_Hit->fColWireNum==28 && t2event_Hit->fIndWireNum==10) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==177 && t2event_Hit->fColWireNum==30 && t2event_Hit->fIndWireNum==10) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==202 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==13) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==204 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==21) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==206 && t2event_Hit->fColWireNum== 9 && t2event_Hit->fIndWireNum==23) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==206 && t2event_Hit->fColWireNum== 8 && t2event_Hit->fIndWireNum==23) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==262 && t2event_Hit->fColWireNum==13 && t2event_Hit->fIndWireNum==27) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==291 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==14) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==341 && t2event_Hit->fColWireNum==23 && t2event_Hit->fIndWireNum== 6) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==358 && t2event_Hit->fColWireNum==11 && t2event_Hit->fIndWireNum==27) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==366 && t2event_Hit->fColWireNum==25 && t2event_Hit->fIndWireNum==16) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==371 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum== 4) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==402 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==15) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==420 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==18) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==427 && t2event_Hit->fColWireNum== 1 && t2event_Hit->fIndWireNum==17) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==429 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==18) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==477 && t2event_Hit->fColWireNum==23 && t2event_Hit->fIndWireNum==16) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==483 && t2event_Hit->fColWireNum==12 && t2event_Hit->fIndWireNum== 2) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==483 && t2event_Hit->fColWireNum==13 && t2event_Hit->fIndWireNum== 2) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==489 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==23) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==491 && t2event_Hit->fColWireNum==15 && t2event_Hit->fIndWireNum==13) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==521 && t2event_Hit->fColWireNum==28 && t2event_Hit->fIndWireNum==23) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==540 && t2event_Hit->fColWireNum== 2 && t2event_Hit->fIndWireNum==15) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==541 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==15) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==545 && t2event_Hit->fColWireNum==12 && t2event_Hit->fIndWireNum==17) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==545 && t2event_Hit->fColWireNum==13 && t2event_Hit->fIndWireNum==17) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==562 && t2event_Hit->fColWireNum==22 && t2event_Hit->fIndWireNum==26) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==622 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==15) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==634 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==23) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==640 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==27) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==681 && t2event_Hit->fColWireNum== 4 && t2event_Hit->fIndWireNum== 7) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==683 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==27) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==710 && t2event_Hit->fColWireNum==14 && t2event_Hit->fIndWireNum==16) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==710 && t2event_Hit->fColWireNum==13 && t2event_Hit->fIndWireNum==18) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==743 && t2event_Hit->fColWireNum==28 && t2event_Hit->fIndWireNum== 4) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==758 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum== 6) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==765 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum== 7) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==787 && t2event_Hit->fColWireNum==27 && t2event_Hit->fIndWireNum==21) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==815 && t2event_Hit->fColWireNum==20 && t2event_Hit->fIndWireNum==25) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==824 && t2event_Hit->fColWireNum==14 && t2event_Hit->fIndWireNum==15) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==828 && t2event_Hit->fColWireNum==27 && t2event_Hit->fIndWireNum==25) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==829 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==19) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==871 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum== 9) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==883 && t2event_Hit->fColWireNum==22 && t2event_Hit->fIndWireNum== 6) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==900 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==14) { badHits.insert( t2event_Hit->fHitID ); }
            if(event==973 && t2event_Hit->fColWireNum==29 && t2event_Hit->fIndWireNum==21) { badHits.insert( t2event_Hit->fHitID ); }
        }


        // Some events show multiple tracks or very bad tracks (selected by eye). Reject them.
        if(event== 55) continue;
        //if(event== 56) continue;
        if(event== 57) continue;
        if(event== 82) continue;
        if(event==177) continue;
        if(event==215) continue;
        if(event==260) continue;
        if(event==291) continue;
        if(event==311) continue;
        if(event==369) continue;
        if(event==370) continue;
        if(event==399) continue;
        if(event==412) continue;
        //if(event==421) continue;
        if(event==481) continue;
        if(event==451) continue;
        if(event==502) continue;
        //if(event==571) continue; // this one is not very bad...
        if(event==579) continue;
        if(event==595) continue;
        //if(event==683) continue; // this one is not very bad...
        if(event==796) continue;
        if(event==823) continue;
        if(event==885) continue;
        if(event==961) continue;
        if(event==985) continue;
        if(event==988) continue;
        
        
        // Loop over all hits and put the ColPulseIDs of those hits, which occur at the boundary of the drift volume in the set 'boundary_hits'
        // boundaries are defined as: x>65mm (where 5mm ~ 3 wire pitches), x<5mm, y>65mm, y<5mm, z>130mm, z<20mm
        //std::set<UInt_t> boundary_hits;
        //for(int hit=0; hit<NHits; hit++) {
        //    RSTPC_Hit * t2event_Hit = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit);
        //    if( t2event_Hit->fIndWireNum<3 || t2event_Hit->fIndWireNum>28 || t2event_Hit->fColWireNum<3 || t2event_Hit->fColWireNum>28 || (t2event_Hit->fCentreTime/20*e_drift_velocity*10-15)<20 || (t2event_Hit->fCentreTime/20*e_drift_velocity*10-15)>130 ) boundary_hits.insert( t2event_Hit->fHitID ); // -15 in z axis is due to the observation that the first hit occurs ~15 mm after anode and the last hit occurs ~15 mm after cathode
        //}


        // Clear files with 3Dhits, Principal components and barycentre
        clear_file.open("3Dhits.txt");
        if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }
        clear_file.open("PrincipalComponents.txt");
        if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }
        clear_file.open("Barycentre.txt");
        if(clear_file.is_open()) { clear_file << ""; clear_file.close(); }


        // Loop over all hits and look for those with good ColPulseIDs and good IndPulseIDs
        std::vector<float> hits_x, hits_y, hits_z, weights;
        for(int hit=0; hit<NHits; hit++) {
            RSTPC_Hit * t2event_Hit = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit);
            UInt_t ColPulseID = t2event_Hit->fColPulseID;
            UInt_t HitID = t2event_Hit->fHitID;
            if( badColPulseIDs.find(ColPulseID)==badColPulseIDs.end() && badHits.find(HitID)==badHits.end() && CoincidentColPulseIDs.find(ColPulseID)==CoincidentColPulseIDs.end() ) { // if true: good ColPulse is found
                UInt_t IndPulseID = t2event_Hit->fIndPulseID;
                if( badIndPulseIDs.find(IndPulseID) == badIndPulseIDs.end() ) { // if true: good IndPulse is found
                    //          << " \tfCentreTime: "   << t2event_Hit->fCentreTime << std::endl;
                    //t2event_Hit->fX,->fY,->fZ,
                    //           ->fMeanTime,->fCentreTime,
                    //           ->fColWireNum,->fIndWireNum,->fColPulseID,->fIndPulseID,
                    //           ->fMeanHeight,->fLedge,->fRedge

                    // Reject hits at the boundary of the TPC
                    //if ( boundary_hits.find(HitID)==boundary_hits.end() ) {
                    
                    // Only proceed if fX, fY and fMeanTime are > than certain value (or are not -nan)
                    if( t2event_Hit->fX>-999. && t2event_Hit->fY>-999. && t2event_Hit->fCentreTime>-999. ) {
                        // Save the 3D hit coordinates (in units of wire pitches for x and y coordinates)
                        hits_x.push_back(t2event_Hit->fIndWireNum);
                        hits_y.push_back(t2event_Hit->fColWireNum);
                        hits_z.push_back(t2event_Hit->fCentreTime/20 * e_drift_velocity * 10); // [mm]. Factor 1/20 because: 1 sample is 50 ns, so 1/20 us = 50 ns. Factor 10: cm to mm conversion

                        // Give each hit a weight (corresponding to the ADC value, corresponding to ColPulse->fMax)
                        RSTPC_Pulse* CollPulse = ColPulsesMap[t2event_Hit->fColPulseID];
                        weights.push_back(CollPulse->fMax);

                        // Write spatial coordinates to file for offline event display
                        ofstream textfile;
                        textfile.open("3Dhits.txt", ios::out | ios::app);
                        if(textfile.is_open()) {
                           textfile << t2event_Hit->fIndWireNum                            << " "
                                    << t2event_Hit->fColWireNum                            << " "
                                    << t2event_Hit->fCentreTime/20*e_drift_velocity*10 -15 << " " // -15 because the earliest hit is ~15 mm after anode and last hit is ~15 mm after cathode (see histogram z_distribution_of_hits)
                                    << CollPulse->fMax                                     << std::endl;
                           textfile.close();
                        }
                        z_distribution_of_hits->Fill(t2event_Hit->fCentreTime/20*e_drift_velocity*10);
                    }
                    //}
                }
            }
        }

        //std::cout << " Found " << hits_x.size() << " good hits, rejecting those having > 1 ColPulses within 2.5 us (50samples)." << std::endl;


        // Erase the sets
        badColPulseIDs.erase( badColPulseIDs.begin(), badColPulseIDs.end() );
        badIndPulseIDs.erase( badIndPulseIDs.begin(), badIndPulseIDs.end() );


        // Require at least 2 hits to perform the PCA
        n_hits_per_event->Fill(hits_x.size());
        if( hits_x.size()<2 || hits_y.size()<2 || hits_z.size()<2 ) continue;
        

        // Initial Principal Components Analysis (PCA)
        // ===========================================
        //std::cout << " ----------------------------------------------- " << std::endl;
        //std::cout << " Initial Principal Components Analysis: "          << std::endl;
        //std::cout << " ----------------------------------------------- " << std::endl;

        std::vector<double> Lambda_PC_and_barycentre = principal_components(hits_x,hits_y,hits_z,weights);

        std::vector<double> eigenvalues(3,0.);
        std::vector<double> eigenvector1(3,0.);
        std::vector<double> eigenvector2(3,0.);
        std::vector<double> eigenvector3(3,0.);
        std::vector<double> barycentre(3,0.);

        for(int i=0; i<3; i++) {
            eigenvalues[i]  = Lambda_PC_and_barycentre[i];
            eigenvector1[i] = Lambda_PC_and_barycentre[i+3];
            eigenvector2[i] = Lambda_PC_and_barycentre[i+6];
            eigenvector3[i] = Lambda_PC_and_barycentre[i+9];
            barycentre[i]   = Lambda_PC_and_barycentre[i+12];
        }



        // Reject tracks which have a big principal component orthogonal to the main principal component
        // ==============================================================================================
        double eigenvalue_threshold = 4.5;
        if(pow(eigenvalues[0],2)+pow(eigenvalues[1],2)>eigenvalue_threshold) {
            n_hits -= hits_x.size();
            n_good_events -= 1;
            //std::cout << " Removed event nr. " << event << " from the analysis because the track is too broad. " << std::endl;
            continue;
        }

        //std::cout << " Eigenvalues:             (" << eigenvalues[0]  << " , " << eigenvalues[1]  << " , " << eigenvalues[2]  << ")" << std::endl;
        //std::cout << " 1st eigenvector (x,y,z): (" << eigenvector1[0] << " , " << eigenvector1[1] << " , " << eigenvector1[2] << ")" << std::endl;
        //std::cout << " 2nd eigenvector (x,y,z): (" << eigenvector2[0] << " , " << eigenvector2[1] << " , " << eigenvector2[2] << ")" << std::endl;
        //std::cout << " 3rd eigenvector (x,y,z): (" << eigenvector3[0] << " , " << eigenvector3[1] << " , " << eigenvector3[2] << ")" << std::endl;
        //std::cout << " Barycentre (x,y,z):      (" << barycentre[0]   << " , " << barycentre[1]   << " , " << barycentre[2]   << ")" << std::endl;

        /*ofstream textfile2;
        textfile2.open("PrincipalComponents.txt", ios::out | ios::app);
        if(textfile2.is_open()) {
            // Only draw the main principal component (which is the third axis: mu travelling in z direction)
            textfile2 << eigenvector3[0] << " " << eigenvector3[1] << " " << eigenvector3[2] << std::endl;
            textfile2.close();
        }

        textfile2.open("Barycentre.txt", ios::out | ios::app);
        if(textfile2.is_open()) {
            textfile2 << barycentre[0] << " " << barycentre[1] << " " << barycentre[2] << std::endl;
            textfile2.close();
        }*/



        // Inserting back those hits which have ColPulse coincidences (stored in the set 'CoincidentColPulseIDs')
        // ======================================================================================================
        // If hit has a distance of < 1 wire pitch (ADJUST VALUE!), reinsert the hit
        std::vector<double> dist_hit_to_PC(3,-999.);
        for(int hit=0; hit<NHits; hit++) {
            RSTPC_Hit * t2event_Hit = (RSTPC_Hit *)tclonesarray_T2_Hits->At(hit);
            UInt_t ColPulseID = t2event_Hit->fColPulseID;
            if( CoincidentColPulseIDs.find(ColPulseID) != CoincidentColPulseIDs.end()) { // if true: Hit with ColPulse element set 'CoincidentColPulseIDs' is found
                if( t2event_Hit->fX>-999. && t2event_Hit->fY>-999. && t2event_Hit->fCentreTime>-999. ) {
                    dist_hit_to_PC = calculate_distance(t2event_Hit->fIndWireNum,t2event_Hit->fColWireNum,t2event_Hit->fCentreTime/20*e_drift_velocity*10,eigenvector3,barycentre);
                    //std::cout << " Res_x: " << dist_hit_to_PC[0] << " Res_y: " << dist_hit_to_PC[1] << " Res_z: " << dist_hit_to_PC[2] << std::endl;

                    // Reinsert the hits which have residuals smaller than 1 wire pitch (=52.5/31mm) in all three coordinates (x,y,z)
                    if( fabs(dist_hit_to_PC[0])<1 && fabs(dist_hit_to_PC[1])<1 && fabs(dist_hit_to_PC[2])<52.5/31. ) {
                        hits_x.push_back(t2event_Hit->fIndWireNum);
                        hits_y.push_back(t2event_Hit->fColWireNum);
                        hits_z.push_back(t2event_Hit->fCentreTime/20*e_drift_velocity*10);
                        RSTPC_Pulse* CollPulse = ColPulsesMap[t2event_Hit->fColPulseID];
                        weights.push_back(CollPulse->fMax);
                        //std::cout << " Found another good hit: " << hit << std::endl;

                        // Write spatial coordinates to file for offline event display
                        ofstream textfile;
                        textfile.open("3Dhits.txt", ios::out | ios::app);
                        if(textfile.is_open()) {
		                	textfile << t2event_Hit->fIndWireNum                        << " "
                                     << t2event_Hit->fColWireNum                        << " "
                                     << t2event_Hit->fCentreTime/20*e_drift_velocity*10 << " "
                                     << CollPulse->fMax                                 << std::endl;
                            textfile.close();
                        }
                    }
                }
            }
        }


        // Only consider events with at least minimum_hits
        if( hits_x.size()<=minimum_hits || hits_y.size()<=minimum_hits || hits_z.size()<=minimum_hits ) continue;


        //Counting the number of hits to calculate the average number of hits per good event
        n_hits += hits_x.size();
        n_good_events += 1;


        // Second Principal Components Analysis (PCA)
        // ===========================================
        //std::cout << " ----------------------------------------------- " << std::endl;
        //std::cout << " 2nd Principal Components Analysis: "              << std::endl;
        //std::cout << " ----------------------------------------------- " << std::endl;
        Lambda_PC_and_barycentre = principal_components(hits_x,hits_y,hits_z,weights);

        for(int i=0; i<3; i++) {
            eigenvalues[i] 	= Lambda_PC_and_barycentre[i];
            eigenvector1[i] = Lambda_PC_and_barycentre[i+3];
            eigenvector2[i] = Lambda_PC_and_barycentre[i+6];
            eigenvector3[i] = Lambda_PC_and_barycentre[i+9];
            barycentre[i] 	= Lambda_PC_and_barycentre[i+12];
        }

        //std::cout << " Eigenvalues:             (" << eigenvalues[0]  << " , " << eigenvalues[1]  << " , " << eigenvalues[2]  << ")" << std::endl;
        //std::cout << " 1st eigenvector (x,y,z): (" << eigenvector1[0] << " , " << eigenvector1[1] << " , " << eigenvector1[2] << ")" << std::endl;
        //std::cout << " 2nd eigenvector (x,y,z): (" << eigenvector2[0] << " , " << eigenvector2[1] << " , " << eigenvector2[2] << ")" << std::endl;
        //std::cout << " 3rd eigenvector (x,y,z): (" << eigenvector3[0] << " , " << eigenvector3[1] << " , " << eigenvector3[2] << ")" << std::endl;
        //std::cout << " Barycentre (x,y,z):      (" << barycentre[0]   << " , " << barycentre[1]   << " , " << barycentre[2]   << ")" << std::endl;

        ofstream textfile;
        textfile.open("PrincipalComponents.txt", ios::out | ios::app);
        if(textfile.is_open()) {
            // Only draw the main principal component (which is the third axis: mu travelling in z direction)
            textfile << eigenvector3[0] << " "
                     << eigenvector3[1] << " "
                     << eigenvector3[2] << std::endl;
            textfile.close();
        }

        textfile.open("Barycentre.txt", ios::out | ios::app);
        if(textfile.is_open()) {
            textfile << barycentre[0] << " " << barycentre[1] << " " << barycentre[2] << std::endl;
            textfile.close();
        }


        // Sort all hits according to the time (z coordinate), earliest_time = smallest_z first
        // ====================================================================================
        std::vector<float> hits_x_sorted(hits_x.size(),-999.9);
        std::vector<float> hits_y_sorted(hits_y.size(),-999.9);
        std::vector<float> hits_z_sorted(hits_z.size(),-999.9);
        std::vector<float> hits_z_temp(hits_z.size(),-999.9);
        std::vector<float> hits_z_backup(hits_z.size(),-999.9);
        for(int hit=0; hit<hits_z.size(); hit++) {
            hits_z_temp[hit] = hits_z[hit];
            hits_z_backup[hit] = hits_z[hit];
        }
        // sort hits_z_temp, largest first:
        std::sort(hits_z_temp.begin(), hits_z_temp.end());
        std::reverse(hits_z_temp.begin(), hits_z_temp.end());
        // Fill vectors hits_x_sorted, hits_y_sorted and hits_z_sorted which are sorted according to the time (first hit first)
        for(int hit_sorted=0; hit_sorted<hits_z.size(); hit_sorted++) {
            for(int hit_unsorted=0; hit_unsorted<hits_z_temp.size(); hit_unsorted++) {
                if(hits_z_temp[hit_sorted]==hits_z[hit_unsorted]) {
                    hits_x_sorted[hit_sorted] = hits_x[hit_unsorted]; // unit: wire number
                    hits_y_sorted[hit_sorted] = hits_y[hit_unsorted]; // unit: wire number
                    hits_z_sorted[hit_sorted] = hits_z[hit_unsorted]; // unit: mm (0: at anode; 150: at cathode)
                    // Set value of hits_z_temp[hit_sorted] to -888.8 and hits_z[hit_unsorted] to -777.7 such that this hit is not used a second time (if there are >1 hits happening at the same time!)
                    hits_z_temp[hit_sorted] = -888.8;
                    hits_z[hit_unsorted] = -777.7;
                    break;
                }
            }
        }
        // Restore hits_z vector
        for(int hit=0; hit<hits_z.size(); hit++) {
            hits_z[hit] = hits_z_backup[hit];
        }

        // Print hit's coordinates to the screen
        //for(int hit=0; hit<hits_x_sorted.size(); hit++){
        //    std::cout << " hit: " << hit << " \tx: " << hits_x_sorted[hit] << " \ty: " << hits_y_sorted[hit] << " \tz: " << hits_z_sorted[hit] << std::endl;
        //}


        // Computing residual vector for each 3D hit to the straight line (the principal component)
        // ========================================================================================
        //std::cout << " ----------------------------------------------- " << std::endl;
        //std::cout << " Computing the residuals: "                        << std::endl;
        //std::cout << " ----------------------------------------------- " << std::endl;

        // Initialize a vector (length = number of hits) of vectors (with length = 3 for the residuals in each spatial direction)
        std::vector<std::vector<double>> residuals(hits_x_sorted.size(),std::vector<double>(3,0.));

        // Calculate the residuals for all good hits in the event
        // Note that the residual vector points from the hit to the closest point on the principal component's line
        residuals = calculate_residuals(hits_x_sorted,hits_y_sorted,hits_z_sorted,eigenvector3,barycentre); // output: in units of wire pitches (=52.5/31mm)
        
        
        // Computing for each hit the variables eta and l
        // ========================================================================================
        // (which is defined as 0 for the first hits z coord. and 1 for the last hits z coord.)
        //std::cout << " ----------------------------------------------- " << std::endl;
        //std::cout << " Computing eta and l: "                            << std::endl;
        //std::cout << " ----------------------------------------------- " << std::endl;
        std::vector<double> eta(hits_x_sorted.size(),-999.);
      	std::vector<double> l(hits_x_sorted.size(),-999.);

        eta = calculate_eta(hits_z_sorted,residuals);
        l   = calculate_l(hits_x_sorted,hits_y_sorted); // output: in units of wire pitches (=52.5/31mm)
        
        //for(int hit=0; hit<hits_x_soreted.size(); hit++) {
            //std::cout << " x: " << hits_x_sorted[hit] << " \ty: " << hits_y_sorted[hit] << " \tl: " << l[hit] << std::endl;
        //}
        
        
        // Computing for each hit the region (shown in file 'regions.png') the residuals
        // ========================================================================================
        std::vector<int> region(hits_x_sorted.size()); // Give each hit an index between 0 and 24 for each region

        region = calculate_region(hits_x_sorted,hits_y_sorted,hits_z_sorted,eigenvector3,barycentre);
        //for(int hit=0; hit<region.size(); hit++) std::cout << region[hit] << std::endl;


        // Write Hits, residuals and eta in .txt file
        // ========================================================================================
        textfile.open("Hits_Residuals_Eta_l_region.txt", ios::out | ios::app);
        if(textfile.is_open()) {
            for(int hit=0; hit<hits_x_sorted.size(); hit++) {
                //std::cout << " Hit "       << hit
                //          << " \t \tr_x: " << residuals[hit][0]
                //          << " \tr_y: "    << residuals[hit][1]
                //          << " \tr_z: "    << residuals[hit][2] << std::endl;
                textfile << hits_x_sorted[hit]       << " "
                         << hits_y_sorted[hit]       << " "
                         << hits_z_sorted[hit]       << " "
                         << residuals[hit][0] << " "
                         << residuals[hit][1] << " "
                         << residuals[hit][2] << " "
                         << eta[hit]          << " "
                         << l[hit]            << " "
                         << region[hit]       << std::endl;
            }
            //std::cout << " Hits, residuals and eta have been written to file 'Hits_Residuals_Eta_l_region.txt " << std::endl;
            textfile.close();
        }



        // Make 2D histograms to show the hits
        // ========================================================================================
        //plot_good_hits_of_event(hits_x_sorted,hits_y_sorted,hits_z_sorted,weights,run_number,event); // Function defined in TestPulses.C



        // Fill histogram with maximal drift distances
        // ============================================
        max_drift_distances->Fill(hits_z_sorted[0]-hits_z_sorted[hits_z_sorted.size()-1]);
        //std::cout << " Max. hit's z position: " << hits_z_sorted[0] << std::endl;
        if( (hits_z_sorted[0]-hits_z_sorted[hits_z_sorted.size()-1])>140) {
            std::cout << " Event " << event << " has max. hit's z pos. separation of " << hits_z_sorted[0]-hits_z_sorted[hits_z_sorted.size()-1] << " > 140 mm: " << std::endl;
        }


        // Fill the histogram to show the distance between the hits (measured as projections on the PC)
        // ============================================================================================
        double delta_x, delta_y, delta_z; // differences of two neighboring hits on the PC
        for(int hit=0; hit<(hits_x_sorted.size()-1); hit++) {
            delta_x = hits_x_sorted[hit] + residuals[hit][0] - hits_x_sorted[hit+1] - residuals[hit+1][0];
            delta_y = hits_y_sorted[hit] + residuals[hit][1] - hits_y_sorted[hit+1] - residuals[hit+1][1];
            delta_z = hits_z_sorted[hit] + residuals[hit][2] - hits_z_sorted[hit+1] - residuals[hit+1][2];
            distance_between_hits->Fill( sqrt( pow(delta_x,2) + pow(delta_y,2) + pow(delta_z,2) ) );
        }



        // Computing theta, the scattering angle of the Multiple Coulomb Scattering (MSC)
        // ========================================================================================
        //std::cout << " ----------------------------------------------- " << std::endl;
        //std::cout << " Computing theta: "                                << std::endl;
        //std::cout << " ----------------------------------------------- " << std::endl;

        // Initialize vectors (length = num of hits) of vectors (lengths = num of hits for the angles)
        std::vector<std::vector<double>> MSC_angles(hits_x_sorted.size(),std::vector<double>(hits_x_sorted.size()-1,-999.));
        std::vector<std::vector<double>> MSC_lengths(hits_x_sorted.size(),std::vector<double>(hits_x_sorted.size()-1,-999.));


        // Calculate the scattering angles of all combinations of hits (with time increasing!)
        ////// MSC_angles = calculate_theta_method1(hits_x_sorted,hits_y_sorted,hits_z_sorted);
        ////// MSC_lengths = calculate_length_method1(hits_x_sorted,hits_y_sorted,hits_z_sorted);
        
        
        // Calculate the scattering angles of all hits occurring later than the first four hits (initial direction equals the PC of the first four hits)
        // First, calculate the principal component of the first four hits
         std::vector<float> first_hits_x(4,-999.9);
         std::vector<float> first_hits_y(4,-999.9);
         std::vector<float> first_hits_z(4,-999.9);
         std::vector<float> first_hits_weigths(4,1.);
         for(int hit=0; hit<4; hit++) {
             first_hits_x[hit] = hits_x_sorted[hit];
             first_hits_y[hit] = hits_y_sorted[hit];
             first_hits_z[hit] = hits_z_sorted[hit];
         }
         std::vector<double> PC_of_first_four_hits = principal_components(first_hits_x,first_hits_y,first_hits_z,first_hits_weigths); // For the moment: give all hits the same weight
         std::vector<double> first_hits_eigenvector(3,0.);
         for(int i=0; i<3; i++) { first_hits_eigenvector[i] = PC_of_first_four_hits[i+9]; }
         //std::cout << " PC of first " <<  first_hits_x.size() << " hits: " << first_hits_eigenvector[0] << " \t" << first_hits_eigenvector[1] << " \t" << first_hits_eigenvector[2] << std::endl;
         MSC_angles = calculate_theta_method2(first_hits_eigenvector,hits_x_sorted,hits_y_sorted,hits_z_sorted);
         MSC_lengths = calculate_length_method2(first_hits_eigenvector,hits_x_sorted,hits_y_sorted,hits_z_sorted);
        
        
        // Calculate the scattering angles of all hits combinations of even hits as well as all combinations of odd hits
        //MSC_angles = calculate_theta_method3(hits_x_sorted,hits_y_sorted,hits_z_sorted);
        //MSC_lengths = calculate_length_method3(hits_x_sorted,hits_y_sorted,hits_z_sorted);



        // Write angles and lengths in .txt file
        // ========================================================================================
        for(int entry1=0; entry1<(hits_x_sorted.size()); entry1++) {
            //std::cout << " ====== entry1: " << entry1 << " ====== " << std::endl;
            for(int entry2=0; entry2<hits_x_sorted.size(); entry2++) {
                if(MSC_angles[entry1][entry2]!=-111) {
                    //std::cout << " entry2: " << entry2 << "  \tangle: " << MSC_angles[entry1][entry2] << "\t\tlength: " << MSC_lengths[entry1][entry2] << std::endl;
                    textfile.open("Angles_and_lengths.txt", ios::out | ios::app);
                    if(textfile.is_open()) {
                        textfile << MSC_angles[entry1][entry2]  << " " <<
                        MSC_lengths[entry1][entry2] << std::endl;         
                    }
                    textfile.close();
                }
            }
        }
        //std::cout << " Angles and lengths of MSC have been written to file 'Angles_and_lengths.txt " << std::endl;


    } // End loop over all events in T2


    //calculate_theta(vec_x,vec_y,vec_z);


    // Plot drift distances
    average_hits_per_good_event = (double)n_hits / n_good_events;
    std::cout << std::endl;
    std::cout << " ------------------------------------"                         << std::endl;
    std::cout << " Good events: "                 << n_good_events               << std::endl;
    std::cout << " Total number of hits: "        << n_hits                      << std::endl;
    std::cout << " Average hits per good event: " << average_hits_per_good_event << std::endl;
    std::cout << " ------------------------------------"                         << std::endl;

    gROOT->SetBatch(0);
    TCanvas * canvas_drift_distances = new TCanvas("canvas_drift_distances","canvas_drift_distances");
    gStyle->SetOptStat(1);
    max_drift_distances->GetXaxis()->SetTitle("drift distance [mm]");
    max_drift_distances->GetYaxis()->SetTitle("entries [-]");
    max_drift_distances->GetXaxis()->SetTitleOffset(1.4);
    max_drift_distances->GetYaxis()->SetTitleOffset(1.2);
    double y_min = 0.;
    double y_max = 50.;
    double max_drift_length = 150.;
    max_drift_distances->GetYaxis()->SetRangeUser(y_min,y_max);
    max_drift_distances->Draw();
    TLine * hline_max_drift_distances = new TLine(max_drift_length,y_min,max_drift_length,0.995*y_max); // (xbegin,ybegin,xend,yend)
    hline_max_drift_distances->SetLineColor(kRed);
    hline_max_drift_distances->Draw();
    canvas_drift_distances->SaveAs("plots/drift_distances.png");


    // Plot hits per events
    gROOT->SetBatch(1);
    TCanvas * canvas_hits_per_event = new TCanvas("canvas_hits_per_event","canvas_hits_per_event");
    gStyle->SetOptStat(1);
    n_hits_per_event->GetXaxis()->SetTitle("hits per event [-]");
    n_hits_per_event->GetYaxis()->SetTitle("entries [-]");
    n_hits_per_event->GetXaxis()->SetTitleOffset(1.4);
    n_hits_per_event->GetYaxis()->SetTitleOffset(1.2);
    y_min = 0.;
    y_max = 270.;
    n_hits_per_event->GetYaxis()->SetRangeUser(y_min,y_max);
    n_hits_per_event->Draw();
    TLine * hline_n_hits_per_event = new TLine(minimum_hits,y_min,minimum_hits,0.995*y_max); // (xbegin,ybegin,xend,yend)
    hline_n_hits_per_event->SetLineColor(kRed);
    hline_n_hits_per_event->Draw();
    canvas_hits_per_event->SaveAs("plots/hits_per_event.png");
    
    
    // Plot distance between two neighboring hits
    gROOT->SetBatch(0);
    TCanvas * canvas_distance_between_hits = new TCanvas("canvas_distance_between_hits","canvas_distance_between_hits");
    gStyle->SetOptStat(1);
    distance_between_hits->GetXaxis()->SetTitle("distance between hits [mm]");
    distance_between_hits->GetYaxis()->SetTitle("entries [-]");
    distance_between_hits->GetXaxis()->SetTitleOffset(1.4);
    distance_between_hits->GetYaxis()->SetTitleOffset(1.2);
    y_min = 0.;
    y_max = 600.;
    distance_between_hits->GetYaxis()->SetRangeUser(y_min,y_max);
    distance_between_hits->Draw();
    canvas_distance_between_hits->SaveAs("plots/distance_between_hits.png");

    std::cout << " ------------------------------------ " << std::endl;
    std::cout << " Mean distance (measured with hits projected to PC) between hits: " << distance_between_hits->GetMean() << " mm." << std::endl;
    std::cout << " ------------------------------------ " << std::endl;


    // Plot z distribution of the hits
    gROOT->SetBatch(0);
    gStyle->SetOptStat(1);
    TCanvas * canvas_z_distribution_of_hits = new TCanvas("canvas_z_distribution_of_hits","canvas_z_distribution_of_hits");
    gStyle->SetOptStat(0);
    z_distribution_of_hits->GetXaxis()->SetTitle("z [mm]");
    z_distribution_of_hits->GetYaxis()->SetTitle("entries [-]");
    z_distribution_of_hits->GetXaxis()->SetTitleOffset(1.4);
    z_distribution_of_hits->GetYaxis()->SetTitleOffset(1.2);
    y_min = 0.;
    y_max = 270.;
    z_distribution_of_hits->GetYaxis()->SetRangeUser(y_min,y_max);
    z_distribution_of_hits->Draw();
    canvas_z_distribution_of_hits->SaveAs("plots/z_distribution_of_hits.png");



    // Access the T1 and T2 data (via RSTPC_T1wrapper and RSTPC_T2wrapper)
    // ===================================================================
    /*
    RSTPC_T1wrapper * t1w = new RSTPC_T1wrapper("/home/rberner/data/ResistiveShell/merged/RSTPC_Run000002032_Merged.root");
    RSTPC_T2wrapper * t2w = new RSTPC_T2wrapper( t1w );
    //Alternatively: RSTPC_T2wrapper * t2w = new RSTPC_T2wrapper( (TTree*)(t1w->fInfile->Get("T2")) );

    // Check if t1w is initialised:
    if( !t1w->IsInit() )  {
        std::cout << " ERROR: t1w is not initialized: " << std::endl;
        return;
    }
	
    // Read T1 (t1w->fChain is tree)
    std::cout << " Total entries in T1: " << t1w->fChain->GetEntries() << std::endl;
    std::cout << " Tpc \t Tpc \t\t Rms \t Rms \t\t Feb \t Feb \t\t FebT \t FebT \t FebB \t FebB \t Feb " << std::endl;
    std::cout << " Evnt: \t Time: \t\t ColW: \t IndW: \t\t Evnt: \t Time: \t\t Amp: \t TotAmp: Amp: \t TotAmp: TotAmp: " << std::endl << std::endl;
    for(int entry=0; entry<5; entry++) { // }t1w->fChain->GetEntries(); entry++) {
    t1w->GetEntry(entry);
    //t1w->Show(entry);
    std::cout << " " << t1w->TpcEvent       << " \t " << t1w->TpcTime      << " \t " << t1w->RmsColWires[0] << " \t "
                     << t1w->RmsIndWires[5] << " \t " << t1w->FebEvent     << " \t " << t1w->FebTime        << " \t "
                     << t1w->FebTopAmp[2]   << " \t " << t1w->FebTopTotAmp << " \t " << t1w->FebBotAmp[0]   << " \t "
                     << t1w->FebBotTotAmp   << " \t " << t1w->FebTotAmp    << std::endl;
    }

    // Read T2
    std::cout << " Total entries in T2: " << t2w->fChain->GetEntries() << std::endl;
    for(int entry=0; entry<t2w->fChain->GetEntries(); entry++) {
        t2w->GetEntry(entry);
        //t2w->Show(entry);
        std::cout   << " Entry: "             << entry
                    << " \t GoodEvent: "      << t2w->GoodEvent
                    << " \t First ColPulse: " << t2w->ColPulses->First()
                    << " \t n ColPulses: "    << t2w->ColPulses->GetEntries()
                    << std::endl;
    }
    */


}
