/*

	TAU_h115_a15.C

	Michael Manighalam, Brenna Mockler, Daniel Coelho, Elodie Resseguie
	August 7, 2013
	
	Programs that use TAU_h115_a15.C: macroCopy.sh

	This program is a root macro that is programmed to draw and write the following variables to histograms:

	Energy and PT for all the particles
	Phi and Eta of the taus
	DeltaPhi, DeltaEta, and DeltaR between the taus 
	Invariant Mass of some of the tau pairs (12, 14, 32, 34)
	Invariant Mass of all 4 Taus
	The Largest Tau Pt in each event (as well as the 2nd, 3rd, and 4th largest)

	It can also save the histograms as pngs if you uncomment the lines beginning with **SAVEAS** at the bottom. Remember to delete the **SAVEAS**

	This is the template macro, so if you want to change variable plots that are created in singleEventGenerator or MultiEventGenerator, all you have to do is change the code here.


*/

#define TAU_h115_a15_cxx
#include "TAU_h115_a15.h" //this file does not exist, the line is here so that macroCopy.sh can use it to create the macros for the simulations.
#include "deltaPhi.h" //header file for computing delta phi
#include "deltaR.h"   // header file for computing delta r
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <TMath.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <algorithm>

#include <sstream>
#include <vector>
#include <cmath>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TParameter.h>
#include <TH1.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TClonesArray.h>

#define NCANVAS 47 //number of canvases used in drawing the histograms
#define PI 3.14159265358979323846264338327950288419716939937510

//Do not change the following two lines. They are changed in singleEventGenerator and MultiEventGenerator and both of these programs depend on the code remaining the same.
const string mh01("115"); //Set the Higgs Mass
const string ma01("15"); //Set the a01 mass
using namespace reco;


void TAU_h115_a15::Loop()
{
	

	TFile *file = new TFile(("Histograms_h01_" + mh01 + "_a01_" + ma01 + ".root").c_str(), "RECREATE");
	file->cd();

	gStyle->SetOptStat(001111111);
	gStyle->SetStatH(0.05);
	gStyle->SetStatW(0.25);
	  
	stringstream ss;
	TCanvas *can[NCANVAS];
	string canvas_name;

	for(int i = 0; i < NCANVAS; i++)
		{
		ss << i+1;
		canvas_name = "c" + ss.str();
		can[i] = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 600, 450);
		ss.str("");
		}
//the following lines create arrays used in the names of the histograms. kMaxParticle is the number of particles in the simulation process, and it is defined in the header file for each macro (the .h file with the same name as the .C file)
	string part_energy_name[kMaxParticle] = {"Energy_Quark_1","Energy_Quark_2","Energy_W","Energy_h01", "Energy_First_a01", "Energy_Second_a01", "Energy_Positron", "Energy_Neutrino", "Energy_Tau_1", "Energy_Tau_2",  			"Energy_Tau_3", "Energy_Tau_4"};
	string part_pt_name[kMaxParticle] = {"PT_Quark_1","PT_Quark_2","PT_W","PT_h01", "PT_First_a01", "PT_Second_a01", "PT_Positron", "PT_Neutrino", "PT_Tau_1", "PT_Tau_2",  "PT_Tau_3", "PT_Tau_4"};
	string part_theta_name[kMaxParticle] = {"Theta_Quark_1","Theta_Quark_2","Theta_W","Theta_h01", "Theta_First_a01", "Theta_Second_a01", "Theta_Positron", "Theta_Neutrino", "Theta_Tau_1", "Theta_Tau_2", "Theta_Tau_3", 			"Theta_Tau_4"};
	string part_Et_name[kMaxParticle] = {"Transverse_Energy_Quark_1","Transverse_Energy_Quark_2","Transverse_Energy_W","Transverse_Energy_h01", "Transverse_Energy_First_a01", "Transverse_Energy_Second_a01", 			"Transverse_Energy_Positron", "Transverse_Energy_Neutrino", "Transverse_Energy_Tau_1", "Transverse_Energy_Tau_2", "Transverse_Energy_Tau_3", "Transverse_Energy_Tau_4"};
	string part_rapidity_name[kMaxParticle] = {"Rapidity_Quark_1","Rapidity_Quark_2","Rapidity_W","Rapidity_h01", "Rapidity_First_a01", "Rapidity_Second_a01", "Rapidity_Positron", "Rapidity_Neutrino", 			"Rapidity_Tau_1","Rapidity_Tau_2", "Rapidity_Tau_3", "Rapidity_Tau_4"};
	string part_phi_name[kMaxParticle] = {"Phi_Quark_1","Phi_Quark_2","Phi_W","Phi_h01", "Phi_First_a01", "Phi_Second_a01", "Phi_Positron", "Phi_Neutrino", "Phi_Tau_1", "Phi_Tau_2", "Phi_Tau_3", "Phi_Tau_4"};
	string part_eta_name[kMaxParticle] = {"Eta_Quark_1","Eta_Quark_2","Eta_W","Eta_h01", "Eta_First_a01", "Eta_Second_a01", "Eta_Positron", "Eta_Neutrino", "Eta_Tau_1", "Eta_Tau_2", "Eta_Tau_3", "Eta_Tau_4"};

	TH1F *ParticleEnergy[kMaxParticle]; //initializes histogram arrays that will later be filled and drawn
	TH1F *ParticlePT[kMaxParticle];
	TH1F *ParticleTheta[kMaxParticle];
	TH1F *ParticleTransverse_Energy[kMaxParticle];
	TH1F *ParticleRapidity[kMaxParticle];
	TH1F *ParticlePhi[kMaxParticle];
	TH1F *ParticleEta[kMaxParticle];

//creates histograms in the histogram arrays that were initialized above
	for (int i = 0; i < kMaxParticle; i++)
		{	
		ParticleEnergy[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_energy_name[i]).c_str(), ("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_energy_name[i]).c_str(), 200, 0., 500.);
		ParticlePT[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_pt_name[i]).c_str(), ("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_pt_name[i]).c_str(), 200, 0., 200.);
		ParticleTheta[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_theta_name[i]).c_str(), ("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_theta_name[i]).c_str(), 200, -20., 20.);
		ParticleTransverse_Energy[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_Et_name[i]).c_str(), ("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_Et_name[i]).c_str(), 200, 0, 200.);
		ParticleRapidity[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_rapidity_name[i]).c_str(), ("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_rapidity_name[i]).c_str(), 100, -10, 10.);
		ParticlePhi[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_phi_name[i]).c_str(), ("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_phi_name[i]).c_str(), 200, -4, 4.);
		ParticleEta[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_eta_name[i]).c_str(), ("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_eta_name[i]).c_str(), 200, -4, 4.);
		}

//strings used in names of DeltaPhi, DeltaEta, and DeltaR histograms
	string part_DeltaPhi_Tau[6] = {"DeltaPhi_Tau_12","DeltaPhi_Tau_13","DeltaPhi_Tau_14","DeltaPhi_Tau_23","DeltaPhi_Tau_24","DeltaPhi_Tau_34"};
	string part_DeltaEta_Tau[6] = {"DeltaEta_Tau_12","DeltaEta_Tau_13","DeltaEta_Tau_14","DeltaEta_Tau_23","DeltaEta_Tau_24","DeltaEta_Tau_34"};
	string part_DeltaR_Tau[6] = {"DeltaR_Tau_12","DeltaR_Tau_13","DeltaR_Tau_14","DeltaR_Tau_23","DeltaR_Tau_24","DeltaR_Tau_34"};

	TH1F *DeltaPhi_Tau[6]; //initializes DeltaPhi, DeltaEta and DeltaR histogram arrays
	TH1F *DeltaEta_Tau[6];
	TH1F *DeltaR_Tau[6];

//creates histograms in the histogram arrays that were initialized above
	for(int i = 0; i < 6; i++)
		{
		DeltaPhi_Tau[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_DeltaPhi_Tau[i]).c_str(),("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_DeltaPhi_Tau[i]).c_str(),200,-5,5);
		DeltaEta_Tau[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_DeltaEta_Tau[i]).c_str(),("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_DeltaEta_Tau[i]).c_str(),200,-5,5);
		DeltaR_Tau[i] = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_DeltaR_Tau[i]).c_str(),("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + part_DeltaR_Tau[i]).c_str(),200,-10,10);
		}

//initializes and creates histograms for the ordered tau PTs 
     int n_Taus = 4; //define n_Taus to be the total number of Taus in the process generated from Madgraph
     TH1F *Ordered_Tau_PT_hist[n_Taus];
     string Tau_pt_name;
	for(int i = 0; i < n_Taus; i++)
		{
		ss << i+1;
        	Tau_pt_name = "h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + ss.str() + "_Largest_Tau_Pt_in_Each_Event";
        	Ordered_Tau_PT_hist[i] = new TH1F(Tau_pt_name.c_str(), Tau_pt_name.c_str(), 200,0, 200.);
        	ss.str("");
		}


/*	TH1F *DeltaR_a01s  = new TH1F("DeltaR_a01s ","Delta R Between the two a01s",200, -10, 10.);
   	TH1F *Delta_Phi_a01s   = new TH1F("Delta_Phi_a01s","Delta Phi Between the Two a01",200, -20, 20.);
*/
	TH1F *Invariant_Mass_Pair12  = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant_Mass_Pair12").c_str(),("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant Mass of Taus 12").c_str(),200, 0, 20.);
	TH1F *Invariant_Mass_Pair14  = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant_Mass_Pair14").c_str(),("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant Mass of Taus 14").c_str(),200, 0, 200.);
	TH1F *Invariant_Mass_Pair32  = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant_Mass_Pair32").c_str(),("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant Mass of Taus 32").c_str(),200, 0, 200.);
	TH1F *Invariant_Mass_Pair34  = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant_Mass_Pair34").c_str(),("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant Mass of Taus 34").c_str(),200, 0, 20.);
	TH1F *Invariant_Mass_4_Taus  = new TH1F(("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant_Mass_4_Taus").c_str(),("h01 " + mh01 + "GeV, a01 " + ma01 + "GeV: " + "Invariant Mass of the Four Taus").c_str(),200, 100, 	150.);
	  

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {//Bracket that opens the big loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
//	if (Cut(ientry) < 0) continue;


	double PT[kMaxParticle], theta[kMaxParticle], transverse_energy[kMaxParticle], rapidity[kMaxParticle], eta[kMaxParticle], d_phi[kMaxParticle][kMaxParticle], d_eta[kMaxParticle][kMaxParticle], d_r[kMaxParticle][kMaxParticle], Imass[kMaxParticle][kMaxParticle];
		
	for( int i =0; i<kMaxParticle; i++)
		{
		PT[i] = pow (Particle_Px[i]*Particle_Px[i] + Particle_Py[i]*Particle_Py[i], .5);	
		theta[i] = 2*atan(exp(-Particle_Eta[i]));
		transverse_energy[i] = abs(Particle_E[i] * sin(theta[i]));
		rapidity[i] = .5*log((Particle_E[i] + Particle_Pz[i]) / (Particle_E[i] - Particle_Pz[i]));
		for( int j = 0; j < kMaxParticle; j++) 
			{		
			d_eta[i][j] = Particle_Eta[i] - Particle_Eta[j];//delta Eta between the ith and the jth particles
			d_phi[i][j] = deltaPhi(Particle_Phi[i], Particle_Phi[j]);//delta phi between the ith and the jth particles
			d_r[i][j] = sqrt( d_phi[i][j]*d_phi[i][j] + d_eta[i][j]*d_eta[i][j]);//delta R between the ith and the jth particles
			Imass[i][j] = sqrt( pow(Particle_E[i] + Particle_E[j], 2) - pow(Particle_Px[i] + Particle_Px[j],2) - pow(Particle_Py[i] + Particle_Py[j],2) - pow(Particle_Pz[i] + Particle_Pz[j],2));//Invariant Mass between the ith and the jth particles	
			} 
		ParticlePT[i]->Fill(PT[i]);
		ParticleTheta[i]->Fill(theta[i]);
		ParticleTransverse_Energy[i]->Fill(transverse_energy[i]);	
		ParticleRapidity[i]->Fill(rapidity[i]);	
		ParticleEnergy[i]->Fill(Particle_E[i]);
		ParticlePhi[i]->Fill(Particle_Phi[i]);
		ParticleEta[i]->Fill(Particle_Eta[i]);
      		}

	DeltaEta_Tau[0]->Fill(d_eta[8][9]);
	DeltaEta_Tau[1]->Fill(d_eta[8][10]);
	DeltaEta_Tau[2]->Fill(d_eta[8][11]);
	DeltaEta_Tau[3]->Fill(d_eta[9][10]);
	DeltaEta_Tau[4]->Fill(d_eta[9][11]);
	DeltaEta_Tau[5]->Fill(d_eta[10][11]);

	//Delta_Phi_a01s->Fill(d_phi[4][5]);
	DeltaPhi_Tau[0]->Fill(d_phi[8][9]);
	DeltaPhi_Tau[1]->Fill(d_phi[8][10]);
	DeltaPhi_Tau[2]->Fill(d_phi[8][11]);
	DeltaPhi_Tau[3]->Fill(d_phi[9][10]);
	DeltaPhi_Tau[4]->Fill(d_phi[9][11]);
	DeltaPhi_Tau[5]->Fill(d_phi[10][11]);

	//DeltaR_a01s->Fill(d_r[4][5]);
	DeltaR_Tau[0]->Fill(d_r[8][9]);
	DeltaR_Tau[1]->Fill(d_r[8][10]);
	DeltaR_Tau[2]->Fill(d_r[8][11]);
	DeltaR_Tau[3]->Fill(d_r[9][10]);
	DeltaR_Tau[4]->Fill(d_r[9][11]);
	DeltaR_Tau[5]->Fill(d_r[10][11]);

        Invariant_Mass_Pair12->Fill(Imass[8][9]);
        Invariant_Mass_Pair14->Fill(Imass[8][11]);
        Invariant_Mass_Pair32->Fill(Imass[10][9]);
        Invariant_Mass_Pair34->Fill(Imass[10][11]);
	

	double Imass_taus_tot;
	Imass_taus_tot = sqrt( pow(Particle_E[8] + Particle_E[9] + Particle_E[10] + Particle_E[11] , 2) - pow( Particle_Px[8] + Particle_Px[9] + Particle_Px[10] + Particle_Px[11] , 2) - pow( Particle_Py[8] + Particle_Py[9] + Particle_Py[10] + Particle_Py[11],2) - pow(Particle_Pz[8] + Particle_Pz[9] + Particle_Pz[10] + Particle_Pz[11],2)); //sums up the invariant mass of the 4 taus
        Invariant_Mass_4_Taus->Fill(Imass_taus_tot);
	
	double Ordered_Taus_PT[n_Taus];
	int First_Tau=8; //array number of the first Tau (taus are numbered 8-11)	
	for(int i=0; i < n_Taus;i++)
		{
		Ordered_Taus_PT[i] = PT[First_Tau+i];
		}		  
	std::sort(Ordered_Taus_PT,Ordered_Taus_PT+n_Taus); //sorts the tau pts from smallest to largest
	for (int i=0; i< n_Taus; i++)
		{
		Ordered_Tau_PT_hist[i]->Fill(Ordered_Taus_PT[n_Taus-1-i]); //fills the histograms and switches the ordering to go from largest to smallest
		}
	}//bracket that closes the big loop
	
// Now we draw all the histograms onto the canvases and save them as .png files
	
	/* **SAVEAS**
	string directory = "Histograms_h01_" + mh01 + "_a01_" + ma01;

	gSystem->MakeDirectory(directory.c_str());

	directory = directory + "/";*/


	file->Write();

     string Tau_Eta_png;
	string transverse_png;
	string Energy_png;
	string DeltaEta_Tau_png;
	string DeltaPhi_Tau_png;
	string DeltaR_Tau_png;


	int first_Tau = 8;
	int canvas_number = 0; //canvas_number will be the incremented canvas number
	for (int i = 0; i < n_Taus; i++)
		{
		Tau_Eta_png = "h01_" + mh01 + "GeV" + "_a01_" + ma01 + "GeV_" + part_eta_name[first_Tau + i] + ".png";
		can[i]->cd(); //starting at can[0] and going until can[3]
		ParticleEta[first_Tau+i]->Draw(); 
		// **SAVEAS** can[i]->SaveAs((directory + Tau_Eta_png).c_str());//can[0]-[3] saved 
		canvas_number++; //canvas number goes from 0-4
		}

	int ffs_particle = 4; //First Final State Particle (for our process, the first final state particle is the a01)
	int tfs_particles = 8; //Total Final State Particles
	int tmp1 = canvas_number;//so we can have canvas_number be constant throughout each loop iteration
	for (int i= 0; i < tfs_particles; i++)
		{
		transverse_png = "h01_" + mh01 + "_a01_" + ma01 + "_" + part_pt_name[ffs_particle + i] + ".png";
		Energy_png = "h01_" + mh01 + "_a01_" + ma01 + "_" + part_energy_name[ffs_particle + i] + ".png";
		
		can[tmp1 + i]->cd(); //starting at can[4] and going until can[11]
		ParticlePT[ffs_particle + i]->Draw();
		// **SAVEAS** can[tmp1 + i]->SaveAs((directory + transverse_png).c_str());//can[4]-[11] saved
		canvas_number++;

		can[tmp1 + i + tfs_particles]->cd(); //starting at can[12] and going until can[19]
		ParticleEnergy[ffs_particle + i]->Draw();
		// **SAVEAS** can[tmp1 + i + tfs_particles]->SaveAs((directory + Energy_png).c_str());//can[12]-[19] saved
		canvas_number++;

		}

	int Delta_Max = 6; //the number of Delta particle pairs in the process
	int tmp2 = canvas_number;//so we can have canvas_number be constant throughout each loop iteration
	for (int i = 0; i < Delta_Max ; i++)
		{	
		DeltaEta_Tau_png = "h01_" + mh01 + "_a01_" + ma01 + "_" + part_DeltaEta_Tau[i] + ".png";
		DeltaPhi_Tau_png = "h01_" + mh01 + "_a01_" + ma01 + "_" + part_DeltaPhi_Tau[i] + ".png";
		DeltaR_Tau_png = "h01_" + mh01 + "_a01_" + ma01 + "_" + part_DeltaR_Tau[i] + ".png";
		
		can[tmp2 + i]->cd(); //starting at can[20] and going until can[25]
		DeltaEta_Tau[i]->Draw();
		// **SAVEAS** can[tmp2 + i]->SaveAs((directory + DeltaEta_Tau_png).c_str());//can[20]-[25] saved
		canvas_number++;

		can[tmp2 + i + Delta_Max]->cd(); //starting at can[26] and going until can[31]
		DeltaPhi_Tau[i]->Draw();
		// **SAVEAS** can[tmp2 + i + Delta_Max]->SaveAs((directory + DeltaPhi_Tau_png).c_str());//can[26]-[31] saved
		canvas_number++;

		can[tmp2 + i + 2*Delta_Max]->cd(); //starting at can[32] and going until can[37]
		DeltaR_Tau[i]->Draw();
		// **SAVEAS** can[tmp2 + i + 2*Delta_Max]->SaveAs((directory + DeltaR_Tau_png).c_str());//can[32]-[37] saved
		canvas_number++;

		}
	
	int tmp3 =  canvas_number;
	for (int i = 0; i < n_Taus; i++)
		{
		can[tmp3 + i]->cd();
		Ordered_Tau_PT_hist[i]->Draw();
		canvas_number++;
		}

	can[canvas_number]->cd();
	Invariant_Mass_Pair12->Draw(); 
	// **SAVEAS** can[canvas_number]->SaveAs((directory + "h01_" + mh01 + "_a01_" + ma01 + "_" + "Invariant Mass Tau 1 and 2.png").c_str());
	
	can[++canvas_number]->cd();   
	Invariant_Mass_Pair14->Draw();   
	// **SAVEAS** can[canvas_number]->SaveAs((directory + "h01_" + mh01 + "_a01_" + ma01 + "_" + "Invariant Mass Tau 1 and 4.png").c_str());

	can[++canvas_number]->cd();
	Invariant_Mass_Pair32->Draw();  
	// **SAVEAS** can[canvas_number]->SaveAs((directory + "h01_" + mh01 + "_a01_" + ma01 + "_" + "Invariant Mass Tau 2 and 3.png").c_str());	

	can[++canvas_number]->cd();
	Invariant_Mass_Pair34->Draw();  
	// **SAVEAS** can[canvas_number]->SaveAs((directory + "h01_" + mh01 + "_a01_" + ma01 + "_" + "Invariant Mass Tau 3 and 4.png").c_str());

	can[++canvas_number]->cd();
	Invariant_Mass_4_Taus->Draw(); 
	// **SAVEAS** can[canvas_number]->SaveAs((directory + "h01_" + mh01 + "_a01_" + ma01 + "_" + "Invariant Mass of the four Taus.png").c_str());



}

