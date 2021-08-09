
  
#define SMORAN 1


#include <numeric>      // std::iota
#include <vector>
#include <algorithm>
#include "TMath.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCut.h"
#include "TChain.h"
#include "TApplication.h"
#include "math.h"
#include "TLegendEntry.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TLatex.h"

#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TVirtualFitter.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#include <stdio.h>
#include <stdlib.h>
using namespace std;

//------PARAMETERS-----//
//const int N_Nu = 5;  // number of Nu bins
const double E_max = 2.0; 			// Minimum Energy
const double E_min = 0.5;  			// Maximum Energy
const double limit_xf = 0.1;		// xF cut
const int nshift_E = 99;			// Number of shifts in Energy
const double step_E = 2.5/1000.0;	// Size of Shifts in Energy
//const int nbins = 100;				// Number of Energy bins
const double zcut = 0.7;

const Double_t Q2_min = 1.;
const Double_t Q2_max = 4.;
const Int_t N_Q2= 6.;//6

const Double_t Phi_min = -180.;
const Double_t Phi_max = 180.;
const Int_t N_Phi = 12.; //12

const Double_t Pt2_min = 0.;
const Double_t Pt2_max = 1.5;
const Int_t N_Pt2 = 6.;//6

TString Nuclei_Type;
Double_t energy_shift;


int main(int argc, char *argv[]){

	Nuclei_Type = (TString) argv[1];
  	Int_t N_Nu = atoi(argv[2]);
  	Int_t nbins = atoi(argv[3]);
  	Int_t n = atoi(argv[4]);  // ORDER OF THE CHEBYSHEV FUNC

	//double Nu_min = 3.2 + Nu_bin*((4.2-3.2)/N_Nu); 
	//double Nu_max = Nu_min + (4.2-3.2)/N_Nu;
	Double_t delta_Q2 = (Q2_max-Q2_min)/N_Q2;
	Double_t delta_Phi = (Phi_max-Phi_min)/N_Phi;

	cout << "OK" << endl;

	// Parameters arrays for Fits
	Double_t parD[n+1];
	Double_t parS[n+1];

	TFile *fileFit = new TFile(Form("fits/FIT4D_"+Nuclei_Type+"_%dnubins_cheb%d_Ebins%d.root", N_Nu, n, nbins), "RECREATE");

	//--------Starting Loop Over Nu bins---------//
	for(Int_t Nu_bin = 0; Nu_bin < N_Nu; Nu_bin++){

		double Nu_min = 3.2 + Nu_bin*((4.2-3.2)/N_Nu); 
		double Nu_max = Nu_min + (4.2-3.2)/N_Nu;

		cout << "Nu interval studied : " << Nu_min << " - " << Nu_max << endl;

		TFile *acc = new TFile(Form("acc/SM2Dfout_"+Nuclei_Type+"_%dnubin%d_Ebins%d.root", N_Nu, Nu_bin, nbins)); // Acceptance files

		for (int Q2_bin = 0; Q2_bin < N_Q2; Q2_bin++){

			for (int Phi_bin = 0; Phi_bin < N_Phi; Phi_bin++){

      			for (int Pt_bin = 0; Pt_bin < N_Pt2; Pt_bin++){	
				  
				  	// Fitting Functions
					TF1 * funcD = (TF1*) gROOT->GetFunction(Form("chebyshev%d", n));
					funcD->SetRange(E_min,E_max);
					TF1 * funcS = (TF1*) gROOT->GetFunction(Form("chebyshev%d", n));
					funcS->SetRange(E_min,E_max);

					//-------Fitting Deuterium:--------//
					cout << "Fitting Deuterium" << endl;

					TVirtualFitter::SetMaxIterations(100000);

					TH1F *DAcc = (TH1F*)acc->Get(Form("acc_D_nubin%d_q2bin%d_phi%d_ptbin%d", Nu_bin, Q2_bin, Phi_bin, Pt_bin)); // Getting the right acc histogram

					if (Nu_bin==0){
						for (int i = 0; i <=n; ++i) funcD->SetParameter(i,1); 
					}else{
						funcD->SetParameters(parD);	
					}


					DAcc->Fit(funcD, "R");

					funcD->GetParameters(&parD[0]);  // Saving Fit Parameters

					Double_t chi = funcD->GetChisquare()/funcD->GetNDF();

					while (chi>5){
						funcD->SetParameters(parD);
						DAcc->Fit(funcD, "R");
						funcD->GetParameters(&parD[0]);  // Saving Fit Parameters
						chi = funcD->GetChisquare()/funcD->GetNDF();
					}

					gStyle->SetOptFit(1);
					DAcc->Draw();
				
					// Saving fit
					fileFit->cd();
					DAcc->SetName(Form("fitD_"+Nuclei_Type+"_nubin%d_q2bin%d_phibin%d_ptbin%d", Nu_bin, Q2_bin, Phi_bin, Pt_bin));
					DAcc->Write();


					//--------LOOP OVER THE SHIFTS---------//
					cout << "Starting loop over Shifts" << endl;
					for (int shift = 0; shift <= nshift_E; ++shift){ 



						// Starting Loop over Q2 for Solid Target

			   			cout << "Fitting the Solid Target Shift: " << shift << "  Nu bin: " << Nu_bin << endl;

						TH1F *h = (TH1F*)acc->Get(Form("acc_"+Nuclei_Type+"_shift%d_nubin%d_q2bin%d_phibin%d_ptbin%d", shift, Nu_bin, Q2_bin, Phi_bin, Pt_bin)); // Getting the right acc histogram

						if (shift==0){
							for (int i = 0; i <=n; ++i) funcS->SetParameter(i,1); 
						}else{
							funcD->SetParameters(parS);
						}

			   			TVirtualFitter::SetMaxIterations(100000);
			    		h->Fit(funcS, "R");

			    		funcS->GetParameters(&parS[0]);	  //saving parametes  		
					    Double_t sol_chi = funcS->GetChisquare()/funcS->GetNDF();
					    int count=0;
						while (sol_chi>5){
							funcS->SetParameters(parS);
							h->Fit(funcS, "R");
							funcS->GetParameters(&parS[0]);	
							sol_chi = funcS->GetChisquare()/funcS->GetNDF();
						}

						gStyle->SetOptFit(1);
						h->Draw();
					
						// Saving fit
						fileFit->cd();
						h->SetName(Form("fit_"+Nuclei_Type+"_shift%d_nubin%d_q2bin%d_phibin%d_ptbin%d", shift, Nu_bin, Q2_bin, Phi_bin, Pt_bin));
						h->Write();
					} //END OF LOOP OVER SHIFTS
				}//END OF LOOP OVER PT2
			}//END OF LOOP OVER PHIPQ
		}//END OF LOOP OVER Q2
	}//END OF LOOP OVER NU
	cout << "END, BYE BYE" << endl;
	fileFit->Close();
	return 0;
}