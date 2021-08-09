//CODE FOR THE ACCEPTANCE CALCULATION//
#define SMORAN 1

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TCut.h"
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TTree.h"
#include "TNtuple.h"
#include "math.h"
#include "TMath.h"
#include "TLegendEntry.h"
#include "TGraph.h"

using namespace std;

Double_t Nu_min;
Double_t Nu_max;
Int_t N_Nu;

Double_t delta_Q2;
Double_t delta_Nu;
Double_t delta_Phi;
Double_t delta_Pt2;

Double_t  energy_shift;

TString Target;
TString Nuclei_Type;

Int_t Nu_bin;
Int_t Bundle;
Int_t bundle_size;

const Double_t Q2_min = 1.;
const Double_t Q2_max = 4.;
const Int_t N_Q2= 6.;//6

const Double_t Phi_min = -180.;
const Double_t Phi_max = 180.;
const Int_t N_Phi = 12.; //12

const Double_t Pt2_min = 0.;
const Double_t Pt2_max = 1.5;
const Int_t N_Pt2 = 6.;//6
  
const Double_t E_min = 0.;  
const Double_t E_max = 2.5;  

const double limit_xf = 0.1;         // Cut in xF
const int nshift_E = 99;             // total number of shifts in Energy
const double step_E = 1.0/1000.0;    // size of energy shift
//const Int_t nbins = 100;             // number of energy bins

const Int_t nSimuFiles = 4;

#ifdef SMORAN
  const TString f_location = "/user/b/brooksw/bruno/";
#else
  const TString f_location = "/user/h/hamilton/ThesisProj/data/"; // The location of data and simulation files with ntuples inside
#endif

const TString fd_ext = "_data.root";
const TString fs_ext = "_simul.root";


/*
string toString(int i)
{
	stringstream ss;
	ss.str("");
	ss<<i;
	return ss.str();
}
*/


int main(int argc, char **argv){


	//Input variables

	Nuclei_Type = (TString) argv[1];  // C for Carbon, Fe for Iron and Pb for Lead
	N_Nu = atoi(argv[2]);             // Number of Nu bins
	Nu_bin = atoi(argv[3]);           // Nu bin Index
	Int_t nbins = atoi(argv[4]);

	Double_t Max = 4.2;
	Double_t Min = 3.2;
	Double_t inc = (Max-Min)/N_Nu;

	Double_t Nu_min = Min + Nu_bin*inc;
	Double_t Nu_max = Nu_min + inc;

	Bundle = 0;//(Int_t) std::stoi(argv[4]); // in number of bundle
	bundle_size = 100;//(Int_t) std::stoi(argv[5]);
	int shift_min = Bundle*bundle_size;
	int shift_max = Bundle*bundle_size + bundle_size - 1;
	  
	delta_Q2 = (Q2_max-Q2_min)/N_Q2;
	delta_Phi = (Phi_max-Phi_min)/N_Phi;
	delta_Pt2 = (Pt2_max-Pt2_min)/N_Pt2;
	 
	  
	/*---------------------------CUTS-----------------------------------------*/

	TCut Target_cutD = "TargType==1"; //cut for Deuterium
	TCut Target_cutS = "TargType==2"; //cut for Solid Target
	  
	  
	//Simulation Cuts
	TCut Nu_cut_S = Form("Nu>%f && Nu<%f", Nu_min, Nu_max);
	TCut Q2_cut_S = Form("Q2>%f && Q2<%f", Q2_min, Q2_max);                 
	TCut Phi_cut_S = Form("PhiPQ>%f && PhiPQ<%f", Phi_min, Phi_max);
	#ifdef SMORAN
	  TCut Pt2_cut_S = Form("Pt2>%f && Pt2<%f", Pt2_min, Pt2_max);
	#else
	  TCut Pt2_cut_S = Form("Pt*Pt>%f && Pt*Pt<%f", Pt2_min, Pt2_max);
	#endif

	TCut cuts_simulD = Q2_cut_S&&Nu_cut_S&&Phi_cut_S&&Pt2_cut_S;
	TCut cuts_simulS = Q2_cut_S&&Nu_cut_S&&Phi_cut_S&&Pt2_cut_S;
	  

	TCut xf_cut = "Xf>0.1"; //  Typical xF cut
	TCut Q2_cut, Nu_cut, Phi_cut, Pt2_cut, cuts_loop, xf_mod;  //Loops Cuts                           


	//-----------------------OBTAINING THE NTUPLES-------------------------//

	// DATA NTUPLE

	TChain *data = new TChain("ntuple_data");
	data->Add(Form(f_location + Nuclei_Type + fd_ext));
	data->SetBranchStatus("*",0);
	data->SetBranchStatus("Q2",1);
	data->SetBranchStatus("Nu",1);
	data->SetBranchStatus("PhiPQ",1);
	data->SetBranchStatus("Xf",1);
	data->SetBranchStatus("Zh",1);
	data->SetBranchStatus("TargType",1);
	data->SetBranchStatus("P",1);
	#ifdef SMORAN
	  data->SetBranchStatus("Pt2",1);
	#else
	  data->SetBranchStatus("Pt",1);
	  data->SetBranchStatus("W",1);
	#endif

	// RECONSTRUCTED NTUPLE FOR DEUTERIUM
	TChain *reconstructed_D = new TChain("ntuple_accept");
	#ifdef SMORAN
	  reconstructed_D->Add(Form(f_location + "D" + fs_ext));
	#else
	  for(Int_t q = 0; q < nSimuFiles; q++){
	    reconstructed_D->Add(Form(f_location + "D%d"+ fs_ext, q+1));
	  }
	#endif

	reconstructed_D->SetBranchStatus("*",0);
	reconstructed_D->SetBranchStatus("Q2",1);
	reconstructed_D->SetBranchStatus("Nu",1);
	reconstructed_D->SetBranchStatus("PhiPQ",1);
	reconstructed_D->SetBranchStatus("Xf",1);
	reconstructed_D->SetBranchStatus("Zh",1);
	reconstructed_D->SetBranchStatus("TargType",1);
	reconstructed_D->SetBranchStatus("P",1);

	#ifdef SMORAN
	  reconstructed_D->SetBranchStatus("Pt2",1);
	#else
	  reconstructed_D->SetBranchStatus("Pt",1);
	  reconstructed_D->SetBranchStatus("W",1);
	#endif

	reconstructed_D->Draw(">>list_accD",cuts_simulD,"goff");
	reconstructed_D->SetEventList((TEventList*)gDirectory->Get("list_accD"));

	// RECONSTRUCTED NTUPLE FOR SOLID TARGET
	TChain *reconstructed_S = new TChain("ntuple_accept");
	#ifdef SMORAN
	  reconstructed_S->Add(Form(f_location + Nuclei_Type + fs_ext));
	#else
	    for(Int_t q = 0; q < nSimuFiles; q++){
	      reconstructed_S->Add(Form(f_location + Nuclei_Type + "%d"+ fs_ext, q+1));
	    }
	#endif
	    
	reconstructed_S->SetBranchStatus("*",0);
	reconstructed_S->SetBranchStatus("Q2",1);
	reconstructed_S->SetBranchStatus("Nu",1);
	reconstructed_S->SetBranchStatus("PhiPQ",1);
	reconstructed_S->SetBranchStatus("Xf",1);
	reconstructed_S->SetBranchStatus("Zh",1);
	reconstructed_S->SetBranchStatus("TargType",1);
	reconstructed_S->SetBranchStatus("P",1);

	#ifdef SMORAN
	  reconstructed_S->SetBranchStatus("Pt2",1);
	#else
	  reconstructed_S->SetBranchStatus("Pt",1);
	  reconstructed_S->SetBranchStatus("W",1);
	#endif

	reconstructed_S->Draw(">>list_accS",cuts_simulS,"goff");
	reconstructed_S->SetEventList((TEventList*)gDirectory->Get("list_accS"));


	// THROWN FOR DEUTERIUM
	TChain *thrown_D = new TChain("ntuple_thrown");
	#ifdef SMORAN
	  thrown_D->Add(Form(f_location + "D" + fs_ext));
	#else
	  for(Int_t w = 0; w < nSimuFiles; w++){
	    thrown_D->Add(Form(f_location+"D%d" + fs_ext, w+1));
	  }
	#endif

	thrown_D->SetBranchStatus("*",0);
	thrown_D->SetBranchStatus("Q2",1);
	thrown_D->SetBranchStatus("Nu",1);
	thrown_D->SetBranchStatus("PhiPQ",1);
	thrown_D->SetBranchStatus("Xf",1);
	thrown_D->SetBranchStatus("Zh",1);
	thrown_D->SetBranchStatus("P",1);
	#ifdef SMORAN
	  thrown_D->SetBranchStatus("Pt2",1);
	#else
	  thrown_D->SetBranchStatus("W",1);
	  //thrown_D->SetBranchStatus("TargType",1);
	  thrown_D->SetBranchStatus("Pt",1);
	#endif

	thrown_D->Draw(">>list_thrD",cuts_simulD,"goff");
	thrown_D->SetEventList((TEventList*)gDirectory->Get("list_thrD"));
	  
	// THROWN FOR SOLID TARGET
	TChain *thrown_S = new TChain("ntuple_thrown");
	#ifdef SMORAN
	  thrown_S->Add(Form(f_location + Nuclei_Type + fs_ext));
	#else
	  for(Int_t r = 0; r < nSimuFiles; r++){
	    thrown_S->Add(Form(f_location + Nuclei_Type + "%d" + fs_ext, r+1));
	  }
	#endif

	thrown_S->SetBranchStatus("*",0);
	thrown_S->SetBranchStatus("Q2",1);
	thrown_S->SetBranchStatus("Nu",1);
	thrown_S->SetBranchStatus("PhiPQ",1);
	thrown_S->SetBranchStatus("Xf",1);
	thrown_S->SetBranchStatus("Zh",1);
	thrown_S->SetBranchStatus("P",1);
	#ifdef SMORAN
	  thrown_S->SetBranchStatus("Pt2",1);
	#else
	  thrown_D->SetBranchStatus("W",1);
	  //thrown_S->SetBranchStatus("TargType",1);
	  thrown_S->SetBranchStatus("Pt",1);
	#endif

	thrown_S->Draw(">>list_thrS",cuts_simulS,"goff");
	thrown_S->SetEventList((TEventList*)gDirectory->Get("list_thrS"));


	//  CREATING THE OUTPUT FILE
	#ifdef SMORAN
	  TFile *plots = new TFile(Form("output/SM3Dfout_"+Nuclei_Type+"_%dnubin%d_Ebins%d_step%f.root", N_Nu, Nu_bin, nbins, step_E*1000),"RECREATE");
	#else
	  TFile *plots = new TFile(Form("output/HH2Dfout_"+Nuclei_Type+"_%dnubin%d_Ebins%d_step%f.root", N_Nu, Nu_bin, nbins, step_E*1000),"RECREATE");
	#endif

	//--------CREATING HISTOGRAMS--------//

	TH1F *D = new TH1F("D", "D", nbins, E_min, E_max);           // Deuterium DATA
	TH1F *accD = new TH1F("accD", "accD", nbins, E_min, E_max);  // Deuterium Acceptance

	std::map<int,TH1F*> histograms, acceptances;
	for(int i = shift_min; i<=shift_max; i++){
	  histograms[i] = new TH1F(Form("Nuclei_shift%d",i),Form("Nuclei_shift%d",i),nbins,E_min,E_max);
	  histograms[i]->Sumw2();
	  acceptances[i] = new TH1F(Form("accNuclei_shift%d",i), Form("accNuclei_shift%d",i), nbins, E_min, E_max);
	  acceptances[i]->Sumw2();
	}



	//------------ACCEPTANCE CORRECTION------------//

	Nu_cut = Form("Nu>%f && Nu<%f", Nu_min, Nu_max);  //  Cut for Nu bin

  	for (int j = 0; j < N_Q2; j++){

	  	Q2_cut = Form("Q2>%f && Q2<%f", Q2_min + j*delta_Q2 , Q2_min + (j+1)*delta_Q2);  //  Cut for Q2 bin		

	  	for (int k = 0; k < N_Phi; k++){
	
	      	Phi_cut = Form("PhiPQ>%f && PhiPQ<%f", Phi_min + k*delta_Phi , Phi_min + (k+1)*delta_Phi);	

		  	//----ACCEPTANCE FOR DEUTERIUM-----/
		  	cuts_loop=Phi_cut_S&&Pt2_cut_S&&Q2_cut_S&&Nu_cut&&xf_cut;

		  	TH1F *data_histo = new TH1F("data_histo","",nbins,E_min,E_max);
		  	TH1F *thrown_histo = new TH1F("thrown_histo","",nbins,E_min,E_max);
		  	TH1F *reconstructed_histo = new TH1F("reconstructed_histo","",nbins,E_min,E_max);

		  	data->Draw("Zh*Nu>>data_histo",cuts_loop&&Target_cutD,"goff");
		  	data_histo->Sumw2();

		  	thrown_D->Draw("Zh*Nu>>thrown_histo",cuts_loop,"goff");
		  	thrown_histo->Sumw2();
		  	/*TH1F *hDT = (TH1F*) thrown_histo->Clone();
		  	hDT->SetName(Form("thrown_histoD%d%d", Nu_bin, j));
		  	hDT->Write();*/
		   
		  	reconstructed_D->Draw("Zh*Nu>>reconstructed_histo",cuts_loop,"goff");
		  	reconstructed_histo->Sumw2();
		  	/*TH1F *hDR = (TH1F*) reconstructed_histo->Clone();
		  	hDR->SetName(Form("reconstructed_histoD%d%d", Nu_bin, j));
		  	hDR->Write();*/
		      
		  	TH1F *acceptance_histo = new TH1F("acceptance_histo","",nbins,E_min,E_max);
		  	acceptance_histo->Divide(reconstructed_histo,thrown_histo,1,1,"B");
		  	// Saving ACC histograms to file
		  	TH1F *hDA = (TH1F*) acceptance_histo->Clone();
		  	hDA->SetName(Form("acc_D_nubin%d_q2bin%d_phibin%d", Nu_bin, j, k));
		  	hDA->Write();
		  
		  	TH1F *acceptance_correction_histo = new TH1F("acceptance_correction_histo","",nbins,E_min,E_max);
		  	acceptance_correction_histo->Divide(data_histo,acceptance_histo,1,1);
		  	/*TH1F *hDC = (TH1F*) acceptance_correction_histo->Clone();
		  	hDC->SetName(Form("acc_corr_D%d", Nu_bin));
		  	hDC->Write();
		  	*/
		  	// Filling the Corrected Deuterium histogram
		  	D->Add(acceptance_correction_histo, 1);
		  	accD->Add(acceptance_histo, 1);

		  	delete acceptance_correction_histo;
		  	delete acceptance_histo;
		  	delete data_histo;
		  	delete thrown_histo;
		  	delete reconstructed_histo;
		  	//delete hDT;
		  	//delete hDR;
		  	delete hDA;
		 	//delete hDC;

		  	//---ACCEPTANCE FOR SOLID TARGET::::ENERGY SHIFT LOOP----//
		  	for(int shift = shift_min; shift <= shift_max; shift++){

			    energy_shift = step_E*shift;

			    //Double_t W = TMath::Sqrt(0.938272*0.938272 +2*0.938272*Nu -Q2);
			    //Double_t Pt = TMath::Sqrt(Pt2);


			    // XF MODIFIED CUT
			    #ifdef SMORAN
			      TCut xf_mod = Form("((Nu + 0.9385)*(TMath::Sqrt((P+%f)*(P+%f)-(P+%f)*(TMath::Sqrt(Pt2)/P)*(P+%f)*(TMath::Sqrt(Pt2)/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+%f/Nu)*Nu/(Nu+0.9385))/TMath::Sqrt(0.938272*0.938272 +2*0.938272*Nu -Q2))/((TMath::Sqrt(TMath::Power(TMath::Sqrt(0.938272*0.938272 +2*0.938272*Nu -Q2)*TMath::Sqrt(0.938272*0.938272 +2*0.938272*Nu -Q2)-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*TMath::Sqrt(0.938272*0.938272 +2*0.938272*Nu -Q2)*TMath::Sqrt(0.938272*0.938272 +2*0.938272*Nu -Q2))/2./TMath::Sqrt(0.938272*0.938272 +2*0.938272*Nu -Q2)))>0.1", energy_shift, energy_shift, energy_shift, energy_shift, energy_shift);
			    #else
			      TCut xf_mod = Form("((Nu + 0.9385)*(TMath::Sqrt((P+%f)*(P+%f)-(P+%f)*(Pt/P)*(P+%f)*(Pt/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+%f/Nu)*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*W*W)/2./W))>0.1", energy_shift, energy_shift, energy_shift, energy_shift, energy_shift);
			    #endif

			    cuts_loop=Q2_cut&&Nu_cut&&Phi_cut&&Pt2_cut&&xf_mod;   //NO OLVIDARSE EL XF MOD CUT
			    //TCut cut = Q2_cut&&Nu_cut&&Target_cutS;

			    TH1F *data_histo = new TH1F("data_histo","",nbins,E_min,E_max);
			    TH1F *thrown_histo = new TH1F("thrown_histo","",nbins,E_min,E_max);
			    TH1F *reconstructed_histo = new TH1F("reconstructed_histo","",nbins,E_min,E_max);

			    //TH1F *Data_xF = new TH1F("Data_xF", "", nbins, -1.5,1.5);
			    //data->Draw(Form("((Nu + 0.9385)*(TMath::Sqrt((P+%f)*(P+%f)-(P+%f)*(Pt/P)*(P+%f)*(Pt/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+%f/Nu)*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*W*W)/2./W))>>Data_xF", energy_shift, energy_shift, energy_shift, energy_shift, energy_shift), cut, "goff");
			    //Data_xF->Sumw2();
			    //TH1F *Xf = (TH1F*) Data_xF->Clone();
			    //Xf->SetName(Form("Xf_shift%d_%d%d", shift, Nu_bin, j));
			    //Xf->Write();
			     
			    data->Draw(Form("(Nu*Zh)+%f>>data_histo", energy_shift), cuts_loop&&Target_cutS, "goff");
			    data_histo->Sumw2();

			    thrown_S->Draw(Form("(Nu*Zh)+%f>>thrown_histo", energy_shift), cuts_loop, "goff");
			    thrown_histo->Sumw2();
			    /*TH1F *hT = (TH1F*) thrown_histo->Clone();
			    hT->SetName(Form("thrown_histo"+Nuclei_Type+"_shift%d_%d%d", shift, Nu_bin, j));
			    hT->Write();*/
			     
			    reconstructed_S->Draw(Form("(Nu*Zh)+%f>>reconstructed_histo", energy_shift), cuts_loop, "goff");
			    reconstructed_histo->Sumw2();
			    /*TH1F *hR = (TH1F*) reconstructed_histo->Clone();
			    hR->SetName(Form("reconstructed_histo"+Nuclei_Type+"_shift%d_%d%d", shift, Nu_bin, j));
			    hR->Write();*/
			     

			    // HISTOGRAMAS WITH ACCEPTANCE

			    TH1F *acc_histo = new TH1F("acc_histo", "", nbins, E_min, E_max);
			    acc_histo->Divide(reconstructed_histo, thrown_histo, 1, 1, "B");
			    TH1F *hA = (TH1F*) acc_histo->Clone();
			    hA->SetName(Form("acc_"+Nuclei_Type+"_shift%d_nubin%d_q2bin%d_phibin%d", shift, Nu_bin, j, k));
			    hA->Write();
			    
			    TH1F *acc_corr_histo = new TH1F("acc_corr_histo", "", nbins, E_min, E_max);
			    acc_corr_histo->Divide(data_histo, acc_histo, 1, 1);
			    /*TH1F *hC = (TH1F*) acc_corr_histo->Clone();
			    hC->SetName(Form("acc_corrected"+Nuclei_Type+"_shift%d_nubin%d", shift,  Nu_bin));
			    hC->Write();
			    */

			    // Filling corrected histograms
			    histograms[shift]->Add(acc_corr_histo, 1);
			    acceptances[shift]->Add(acc_histo, 1);


			    delete data_histo;
			    delete thrown_histo;
			    delete reconstructed_histo;
			    delete acc_histo;
			    delete acc_corr_histo;
			    //delete hT;
			    //delete hR;
			    delete hA;
			    //delete Xf;
			    //delete Data_xF;
			    //delete hC;
		  	} // Closing Shifts
		}//Closing the PhiPQ
  	}//Closing the Q2 




	//SAVING CORRECTED HISTOS
	D->Write();
	accD->Write();

	for(int i = shift_min; i<= shift_max; i++){
		histograms[i]->Write();
	   	acceptances[i]->Write();
	}


	//-----------------CLOSING THE FILE-----------------------//
	plots->Close();

	return 0;
}
