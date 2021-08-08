
  
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

TString Nuclei_Type;
Double_t energy_shift;

//-------MODIFIED FEYNMAN X--------//
double Calculate_Modified_Xf(Float_t Shift, Float_t Nu, Float_t P, Float_t Pt, Float_t Q2, Float_t W, Float_t Zh)
{
	//cout << "Shift " << Shift << endl;
  	//double xf = ((Nu + 0.9385)*(TMath::Sqrt(P*P-Pt*Pt) - TMath::Sqrt(Q2+Nu*Nu)*Zh*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0  .1395*0.1395*W*W)/2./W));
  	double xf = ((Nu + 0.9385)*(TMath::Sqrt((P+Shift)*(P+Shift)-(P+Shift)*(Pt/P)*(P+Shift)*(Pt/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+Shift/Nu)*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*W*W)/2./W));
  	//cout<< "Xf " << xf << endl;
  	return xf;
}


int main(int argc, char *argv[]){

	Nuclei_Type = (TString) argv[1];
  	Int_t N_Nu = atoi(argv[2]);
  	Int_t nbins = atoi(argv[3]);
  	Int_t n = atoi(argv[4]);  // ORDER OF THE CHEBYSHEV FUNC

	//double Nu_min = 3.2 + Nu_bin*((4.2-3.2)/N_Nu); 
	//double Nu_max = Nu_min + (4.2-3.2)/N_Nu;


	cout << "The Nuclei type studied is " << Nuclei_Type << endl;
	cout<< "The cut on Xf is " << limit_xf << endl;

	//------Opening data files-----//
	//TFile *file = new TFile(Form("/Users/gbibim/Documents/ElossAnalysis/chargedPions/" + Nuclei_Type + "_data.root"));
	TFile *file = new TFile(Form("/user/b/brooksw/bruno/" + Nuclei_Type + "_data.root"));

	//-----Opening TTree----//

	//------ Hayk's Tuples -----//
	TTree* tree = (TTree*)file->Get("ntuple_data");
	//Reading Branches with appropiate variables. 
	Float_t PID;
	Float_t Nphe;
	Float_t TargType;
	Float_t Q2;
	Float_t Nu;
	Float_t Xb;
	//Float_t W;
	//Float_t SectorEl;
	//Float_t ThetaPQ;
	Float_t PhiPQ;
	Float_t Zh;
	Float_t Pt2;  //Hayk's are Pt
	//Float_t W2p;
	Float_t Xf;
	//Float_t T;
	Float_t P;
	Float_t T4;
	Float_t deltaZ;
	//Float_t NmbPion;
	tree->SetBranchAddress("PID",&PID);
	tree->SetBranchAddress("Nphe",&Nphe);
	tree->SetBranchAddress("TargType",&TargType);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("Nu",&Nu);
	tree->SetBranchAddress("Xb",&Xb);
	//tree->SetBranchAddress("W",&W);
	//tree->SetBranchAddress("SectorEl",&SectorEl);
	//tree->SetBranchAddress("ThetaPQ",&ThetaPQ);
	tree->SetBranchAddress("PhiPQ",&PhiPQ);
	tree->SetBranchAddress("Zh",&Zh);
	tree->SetBranchAddress("Pt2",&Pt2);
	//tree->SetBranchAddress("W2p",&W2p);
	tree->SetBranchAddress("Xf",&Xf);
	//tree->SetBranchAddress("T",&T);
	tree->SetBranchAddress("P",&P);
	tree->SetBranchAddress("T4",&T4);
	tree->SetBranchAddress("deltaZ",&deltaZ);
	//tree->SetBranchAddress("NmbPion",&NmbPion);

	Int_t nentries = tree->GetEntries();
	//Int_t nentries = 100000;

	//-----Creating output file-----//	
	TFile *fout = new TFile(Form("output/Z1D_"+Nuclei_Type+"_%dnubins_cheb%d_Ebins%d.root", N_Nu, n, nbins), "RECREATE");

	//-----Creating the Graphs for the Eloss Shift Values------//
	TGraph *gElossKS = new TGraph();  //  Graph for Eloss values for the KS test
	TGraph *gElossWKS = new TGraph(); //  Graph for Eloss values fot the Weighted KS test
	TGraph *gElossKSb = new TGraph(); //  Graph for Eloss values for the Binned KS test
	TGraph *gElossWKSb = new TGraph(); // Graph for Eloss values for the Binned Weighted KS test
	TGraph *gElossAKSb = new TGraph(); // Graph for Eloss values for the Acceptance corrected KS test

	vector<double> dataS;
	vector<double> dataD;
	vector<double> weightD;
	vector<double> weightS;

	// Parameters arrays for Fits
	Double_t parD[n+1];
	Double_t parS[n+1];

	TFile *fileFit = new TFile(Form("output/OLDFIT1D_"+Nuclei_Type+"_%dnubins_cheb%d_Ebins%d.root", N_Nu, n, nbins), "RECREATE");

	//--------Starting Loop Over Nu bins---------//
	for(Int_t Nu_bin = 0; Nu_bin < N_Nu; Nu_bin++){

		double Nu_min = 3.2 + Nu_bin*((4.2-3.2)/N_Nu); 
		double Nu_max = Nu_min + (4.2-3.2)/N_Nu;

		cout << "Nu interval studied : " << Nu_min << " - " << Nu_max << endl;

		//--------Starting Histos and Graphs-------//

		TGraph *gpKS = new TGraph();      	//Graph for the Unbinned KS Test
		TGraph *gpWKS = new TGraph();     	//Graph for the Unbinned Weighted KS Test
		TGraph *gpKSb = new TGraph();     	// Graph for the Binned KS Test
		TGraph *gpWKSb = new TGraph();    	// Graph for the Binned Weighted KS Test
		TGraph *gpAKSb = new TGraph();      // Graph for the Accepted Binned KS Test

/*
		// Parameters arrays for Fits
		Double_t parD[n+1];
		Double_t parS[n+1];
*/
		// Fitting Functions
		TF1 * funcD = (TF1*) gROOT->GetFunction(Form("chebyshev%d", n));
		funcD->SetRange(E_min,E_max);
		TF1 * funcS = (TF1*) gROOT->GetFunction(Form("chebyshev%d", n));
		funcS->SetRange(E_min,E_max);

		TFile *acc = new TFile(Form("output/SM1Dfout_"+Nuclei_Type+"_%dnubin%d_Ebins%d.root", N_Nu, Nu_bin, nbins)); // Acceptance files

		//-----Histograms with Energy distribution-----//
		TH1F *D = new TH1F("D","D",nbins,E_min,E_max);
		D->Sumw2();
		std::map<int,TH1F*> histograms;
		for(int i = 0; i<= nshift_E; i++){
			histograms[i] = new TH1F(Form("Nuclei_shift%d_Nubin%d",i, Nu_bin),Form("Nuclei_bin%d_Nubin%d",i, Nu_bin),nbins,E_min,E_max);
			histograms[i]->Sumw2();
		}

		//-----Histograms for Weighted Energy Distributions-----//
		TH1F *DW = new TH1F("DW","DW",nbins,E_min,E_max);
		DW->Sumw2();
		std::map<int,TH1F*> histogramsW;
		for(int i = 0; i<= nshift_E; i++){
			histogramsW[i] = new TH1F(Form("WNuclei_shift%d_Nubin%d",i, Nu_bin),Form("WNuclei_bin%d_Nubin%d",i, Nu_bin),nbins,E_min,E_max);
			histogramsW[i]->Sumw2();
		}


		//-------Fitting Deuterium:--------//
		cout << "Fitting Deuterium" << endl;

		TVirtualFitter::SetMaxIterations(100000);

		TH1F *DAcc = (TH1F*)acc->Get("accD"); // Getting the right acc histogram

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
		DAcc->SetName(Form("fitD_"+Nuclei_Type+"_nubin%d", Nu_bin));
		DAcc->Write();


		//--------LOOP OVER THE SHIFTS---------//
		cout << "Starting loop over Shifts" << endl;
		for (int shift = 0; shift <= nshift_E; ++shift){ 



			// Starting Loop over Q2 for Solid Target

   			cout << "Fitting the Solid Target Shift: " << shift << "  Nu bin: " << Nu_bin << endl;

			TH1F *h = (TH1F*)acc->Get(Form("accNuclei_shift%d", shift)); // Getting the right acc histogram

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
			h->SetName(Form("fit_"+Nuclei_Type+"_nubin%d_shift%d", Nu_bin, shift));
			h->Write();

		   	cout << "Entering the loop over entries " << nentries << endl;

		   	// Vectors for data and weights

		   	double w1sum=0., w2sum=0.;

     		energy_shift = step_E*shift;      // is doing the right shift

		    //-------Loop over the ENTRIES: Filling the histos and vectors--------//
		    for(Int_t j=1; j<= nentries; j++){
		     	if(j%5000000==0) cout<< "Processing event number " << j/1000000.0 << "M "<< endl;
		      	tree->GetEntry(j);
		      	//Apply Cuts bin in Nu
		      	if(Nu > Nu_max || Nu < Nu_min) continue; 

		      	if(PID!=211) continue; // && P<2.7 && T4<-0.6 && Nphe<15) continue;

		      	//if ( (T4>0.45 && P>3.3) || (T4>0.5 && P<3.3)) continue;

		    	// Deuterium
		      	if(TargType==1 && Xf>limit_xf && Zh*Nu<E_max && Zh*Nu>E_min && Zh<zcut){

		        	funcD->SetParameters(parD);
		        	double w = 1./(funcD->Eval(Zh*Nu, 0, 0));

					dataD.push_back(Zh*Nu);    	// Unbinned KS
		      		weightD.push_back(w);		// Weighted Unbinned
		        	D->Fill(Zh*Nu);				// Binned KS
		      		DW->Fill(Zh*Nu, w);			// Weighted Binned
		      		w1sum = w1sum + w;
		      	}

		      	// Solid Target
		      	if(TargType==2 && Zh*Nu+energy_shift<E_max && Zh*Nu+energy_shift>E_min){ 
		      		// xF modified cut
		      		Double_t W = TMath::Sqrt(0.938272*0.938272 +2*0.938272*Nu -Q2);
		      		Double_t Pt = TMath::Sqrt(Pt2);
		        	double Xf_Nuclei = Calculate_Modified_Xf(energy_shift, Nu, P, Pt,  Q2,  W, Zh);
		        	if(Xf_Nuclei>limit_xf){

		          		funcS->SetParameters(parS);
		        		double w = 1./(funcS->Eval(Zh*Nu+energy_shift, 0, 0));  //weight for E value

		          		dataS.push_back((Zh*Nu)+energy_shift); 	 	 	// Unbinned KS
		      			weightS.push_back(w);						// Weighted Unbinned
		          		histograms[shift]->Fill(Zh*Nu+energy_shift); 		// Binned KS
		      			histogramsW[shift]->Fill(Zh*Nu+energy_shift, w);	// Weighted Binned
		      			w2sum = w2sum + w;
		        	}                                                
		      	}
	    	}

	    	cout<< "END LOOP OVER ENTRIES" << endl;

		    //-----Binned Kolmogorov-Smirnov Test-----//
		    D->Scale(1.0/D->Integral());
		    histograms[shift]->Scale(1.0/histograms[shift]->Integral());
		    //double pCSbinned = D->Chi2Test(histograms[i], "NORM");
		    double pKSbinned = D->KolmogorovTest(histograms[shift], "D");
		    //gpKSb->SetPoint(i, i, -1*TMath::Log10(pKSbinned));
		    gpKSb->SetPoint(shift, energy_shift, pKSbinned);
/*
		    //-----Accepted Corrected Binned Kolmogorov-Smirnov Test-----//
		    TH1F *Acc = (TH1F*)acc->Get("accD"); // Getting the right acc histogram
		    Acc->Scale(1.0/Acc->Integral());
		    TH1F *Target = (TH1F*)acc->Get(Form("accNuclei_shift%d", shift));
		    Target->Scale(1.0/Target->Integral());
		    //double pCSbinned = D->Chi2Test(histograms[i], "NORM");
		    double ApKSbinned = Acc->KolmogorovTest(Target, "D");
		    //gpKSb->SetPoint(i, i, -1*TMath::Log10(pKSbinned));
		    gpAKSb->SetPoint(shift, shift, ApKSbinned);
*/
		    //-----Binned Weighted KS Test-----//
		    DW->Scale(1.0/DW->Integral());
		    histogramsW[shift]->Scale(1.0/histogramsW[shift]->Integral());
		    double pWKSbinned = DW->KolmogorovTest(histogramsW[shift], "D");
		    //gWpKSb->SetPoint(i, i, -1*TMath::Log10(WpKSbinned));
		    gpWKSb->SetPoint(shift, energy_shift, pWKSbinned);	


		    //--------Unbinned Tests--------//

	    	//--------Sorting vectors and filling energy Graphs-------//
	    	int nD = dataD.size();
	    	int nS = dataS.size();

 			vector<int> idxD(nD);
 			vector<int> idxS(nS);
 			int x=0, y=0;
 			std::iota(idxD.begin(), idxD.end(),x++); //Initializing
 			std::iota(idxS.begin(), idxS.end(),y++); //Initializing
 			sort(idxD.begin(), idxD.end(), [&](int i,int j){return dataD[i]<dataD[j];} );
 			sort(idxS.begin(), idxS.end(), [&](int i,int j){return dataS[i]<dataS[j];} );

	    	sort(dataD.begin(), dataD.end());
	    	sort(dataS.begin(), dataS.end());

		    //----------UNBINNED KOLMOGOROV TEST----------//
		    double pKS = TMath::KolmogorovTest(nD, &dataD[0], nS, &dataS[0], "D");
		    //gpKS->SetPoint(i, i, -1*TMath::Log10(pKS));
		    gpKS->SetPoint(shift, energy_shift, pKS);


		    //Double_t w1sum = std::accumulate(weightD.begin(), weightD.end(), 0);
		    //Double_t w2sum = std::accumulate(weightS.begin(), weightS.end(), 0);

		    //cout << w1sum << "   " << nD << endl;

		    //WEIGHTED KOLMOGOROV TEST :: ROOT IMPLEMENTATION
	        Double_t rdiff = 0;
    		Double_t rdmax = 0;
		    Bool_t ok = kFALSE;
    		int j1 = 0;
    		int j2 = 0;
		    for (int l = 0; l < nD+nS; ++l){
		    	int a = idxD[j1];
		    	int b = idxS[j2];
		    	Double_t w1 = weightD[a]/w1sum;
		    	Double_t w2 = weightS[b]/w2sum;

				if (dataD[j1] < dataS[j2]) {
				    rdiff -= w1;
				    j1++;
				    if (j1 >= nD) {ok = kTRUE; break;}
				} else if (dataD[j1] > dataS[j2]) {
				    rdiff += w2;
				    j2++;
				    if (j2 >= nS) {ok = kTRUE; break;}
				} else {
				    // special cases for the ties
				    double x = dataD[j1];
				    while(j1 < nD && dataD[j1] == x) {
				        rdiff -= w1;
				        j1++;
				    }
				    while(j2 < nS && dataS[j2] == x) {
				        rdiff += w2;
				        j2++;
				    }
				    if (j1 >= nD) {ok = kTRUE; break;}
				    if (j2 >= nS) {ok = kTRUE; break;}
				}
				rdmax = TMath::Max(rdmax,TMath::Abs(rdiff));
			}
		    
		    R__ASSERT(ok);
		    if (ok){
		    	Double_t z = rdmax * TMath::Sqrt(nD*nS/(nD+nS));
       			double pWKS = TMath::KolmogorovProb(z);
       			//gpKS->SetPoint(i, i, -1*TMath::Log10(pKS));
		    	gpWKS->SetPoint(shift, energy_shift, pWKS);
		    }
			dataS.clear();
		   	dataD.clear();
		   	weightD.clear();
		   	weightS.clear();
	  	}//-------End of the Loop over shifts------//

	  	std::cout<< "END LOOP OVER SHIFTS" << std::endl;

		//--------Statistical Tests Graphs--------//
		gpKS->SetName(Form("pKS_"+Nuclei_Type+"_nubin%d", Nu_bin));
		gpWKS->SetName(Form("pWKS_"+Nuclei_Type+"_nubin%d", Nu_bin));
		gpKSb->SetName(Form("pKSb_"+Nuclei_Type+"_nubin%d", Nu_bin));
		gpWKSb->SetName(Form("pWKSb_"+Nuclei_Type+"_nubin%d", Nu_bin));

		gpKS->SetLineColor(4);   
		gpWKS->SetLineColor(1);
		gpKSb->SetLineColor(6); 
		gpWKSb->SetLineColor(2);

		gpKS->SetLineWidth(3);
		gpWKS->SetLineWidth(3);
		gpWKS->SetLineStyle(9);
		gpKSb->SetLineWidth(3);
		gpWKSb->SetLineWidth(3);
		gpWKSb->SetLineStyle(9);

		gpKS->SetTitle("Unbinned KS");
		gpWKS->SetTitle("Weighted Unbinned KS");
		gpKSb->SetTitle("Binned KS");
		gpWKSb->SetTitle("Weighted Binned KS");

		fout->cd();
		gpKS->Write();
		gpWKS->Write();
		gpKSb->Write();
		gpWKSb->Write();

		TCanvas *canvas = new TCanvas();
		TMultiGraph *multi = new TMultiGraph();
		multi->Add(gpKS, "L");
		multi->Add(gpWKS,"L");
		multi->Add(gpKSb,"L");
		multi->Add(gpWKSb,"L");

		multi->Draw("AL");
		multi->SetMaximum(1.2);
		multi->SetTitle(Form("Probability curve for %f<Nu<%f, Target: " + Nuclei_Type + "", float(Nu_min), float(Nu_max)));
		multi->GetYaxis()->SetNdivisions(2);
		multi->GetXaxis()->SetTitle("dE [MeV]");
		multi->GetYaxis()->SetTitle("p_{0}"); //"-Log(p_{0})"
		
		canvas->BuildLegend();
		canvas->SaveAs(Form("output/ZProb_"+Nuclei_Type+"_%dnubin%d_%dentries_%dEcut%d_cheb%d_Ebins%d.pdf", N_Nu, Nu_bin, nentries, int(E_min), int(E_max*100), n, nbins));
		canvas->SaveAs(Form("output/ZProb_"+Nuclei_Type+"_%dnubin%d_%dentries_%dEcut%d_cheb%d_Ebins%d.C", N_Nu, Nu_bin, nentries, int(E_min), int(E_max*100), n, nbins));
	
		//-----ELOSS HISTOGRAMS PLOTS-----//
		
		Double_t elossKS=0;
		Double_t elossWKS=0;
		Double_t elossKSb=0;
		Double_t elossWKSb=0;
		Int_t i_KS=0;
		Int_t i_WKS=0;
		Int_t i_KSb=0;
		Int_t i_WKSb=0;

		for (int i = 0; i < nshift_E; ++i){
			Double_t x, y;

			gpKS->GetPoint(i, x, y);
			if(y>elossKS){
				elossKS=y;
				i_KS = i;
			}

			gpWKS->GetPoint(i, x, y);
			if(y>elossWKS){
				elossWKS=y;
				i_WKS = i;
			}

			gpKSb->GetPoint(i, x, y);
			if(y>elossKSb){
				elossKSb=y;
				i_KSb = i;
			}

			gpWKSb->GetPoint(i, x, y); 
			if(y>elossWKSb){
				elossWKSb=y;
				i_WKSb = i;
			}

		}
		cout << "ELOSS VALUE FOR KS: " << i_KS << "   PROB: " << elossKS << endl;
		cout << "ELOSS VALUE FOR WKS: " << i_WKS << "   PROB: " << elossWKS << endl;
		cout << "ELOSS VALUE FOR KSb: " << i_KSb << "   PROB: " << elossKSb << endl;
		cout << "ELOSS VALUE FOR WKSb: " << i_WKSb << "   PROB: " << elossWKSb << endl;


		//-----Energy spectra distributions-----//

		fout->cd();
		for (int d = 0; d < nshift_E; ++d)
		{
			histograms[d]->Write();
			histogramsW[d]->Write();
		}

		TCanvas *Chu = new TCanvas();

		D->SetName(Form("D_%d", Nu_bin));
		D->SetTitle(Form("Energy for Target "+Nuclei_Type+" Nu bin %d", Nu_bin));
		D->SetLineColor(9);
		D->Draw("HIST L P");
		//D->SetMarkerColor(9);
		//D->SetMarkerStyle(33);
		D->Write();

		histograms[i_KS]->SetName(Form("KS_%d", Nu_bin));
		TH1F *KS = (TH1F*) histograms[i_KS]->Clone();

		histogramsW[i_WKS]->SetName(Form("WKS_%d", Nu_bin));
		TH1F *WKS = (TH1F*) histograms[i_WKS]->Clone();

		histograms[i_KSb]->SetName(Form("KSb_%d", Nu_bin));
		TH1F *KSb = (TH1F*) histograms[i_KSb]->Clone();

		histogramsW[i_WKSb]->SetName(Form("WKSb_%d", Nu_bin));
		TH1F *WKSb = (TH1F*) histograms[i_WKSb]->Clone();

		KS->SetLineColor(4);
		KS->SetMarkerColor(4);
		KS->SetMarkerStyle(20);
		KS->Draw("HIST same L P");

		WKS->SetName(Form("WKS_%d", Nu_bin));
		WKS->SetLineColor(1);
		WKS->SetMarkerColor(1);
		WKS->SetMarkerStyle(30);
		WKS->Draw("HIST same L P");

		KSb->SetName(Form("KSb_%d", Nu_bin));
		KSb->SetLineColor(6);
		KSb->SetMarkerColor(6);
		KSb->SetMarkerStyle(22);
		KSb->Draw("HIST same L P");

		WKSb->SetName(Form("WKSb_%d", Nu_bin));
		WKSb->SetLineColor(2);
		WKSb->SetMarkerColor(2);
		WKSb->SetMarkerStyle(23);
		WKSb->Draw("HIST same L P");

		fout->cd();

		DW->SetName(Form("Weighted_D_%d", Nu_bin));
		DW->Write();
		D->Write();
		KS->Write();
		KSb->Write();
		WKS->Write();
		WKSb->Write();

		auto legend = new TLegend(0.3, 0.1, .5, .3, Form(Nuclei_Type+" Target, Nu bin %d", Nu_bin), "brNDC");
	 
	   	legend->SetNColumns(1);
	 
		legend->AddEntry(D, "Deuterium"); 
	   	legend->AddEntry(KS, Form("KS shift %d", i_KS), "l");
	   	legend->AddEntry(WKS, Form("WKS shift %d", i_WKS), "l");
	   	legend->AddEntry(KSb, Form("KSb shift %d", i_KSb), "l");
	   	legend->AddEntry(WKSb, Form("WKSb shift %d", i_WKSb), "l");
	 
	   	legend->Draw();

		//Chu->BuildLegend();
		Chu->SaveAs(Form("output/ZEnergy_"+Nuclei_Type+"_%dnubin%d_%dentries_%dEcut%d_cheb%d_Ebins%d.pdf", N_Nu, Nu_bin, nentries, int(E_min), int(E_max), n, nbins));
		Chu->SaveAs(Form("output/ZEnergy_"+Nuclei_Type+"_%dnubin%d_%dentries_%dEcut%d_cheb%d_Ebins%d.C", N_Nu, Nu_bin, nentries, int(E_min), int(E_max), n, nbins));
		

		//DAcc->SetName(Form("Acc_Corr_D_%d", Nu_bin));
		//DAcc->Write();

	    //------- ELOSS GRAPHS -------//
	    gElossKS->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_KS);
	    gElossWKS->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_WKS);
	    gElossKSb->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_KSb);
	    gElossWKSb->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_WKSb);
	}//  END OF LOOP OVER NU BINS


	fout->cd();
	
	//  ELOSS PLOTS  //
	gElossKS->SetName(Form("ElossKS_"+Nuclei_Type));
	gElossWKS->SetName(Form("ElossWKS_"+Nuclei_Type));
	gElossKSb->SetName(Form("ElossKSb_"+Nuclei_Type));
	gElossWKSb->SetName(Form("ElossWKS_"+Nuclei_Type));

	gElossKS->SetMarkerColor(4);    
	gElossWKS->SetMarkerColor(1);
	gElossKSb->SetMarkerColor(6);   
	gElossWKSb->SetMarkerColor(2);

	gElossKS->SetMarkerStyle(20);
	gElossWKS->SetMarkerStyle(30);
	gElossKSb->SetMarkerStyle(22);
	gElossWKSb->SetMarkerStyle(23);

	gElossKS->SetTitle("Unbinned KS");
	gElossWKS->SetTitle("Weighted Unbinned KS");
	gElossKSb->SetTitle("Binned KS");
	gElossWKSb->SetTitle("Weighted Binned KS");

	gElossKS->Write();
	gElossWKS->Write();
	gElossKSb->Write();
	gElossWKSb->Write();

	TCanvas *canvas = new TCanvas();
	TMultiGraph *multi = new TMultiGraph();
	multi->Add(gElossKS, "");
	multi->Add(gElossWKS,"");
	multi->Add(gElossKSb,"");
	multi->Add(gElossWKSb, "");

	multi->Draw("AP");
	multi->Write();
	multi->SetMaximum(100);
	multi->SetTitle(Form("Energy Loss for %d Nu bins, Target: " + Nuclei_Type + "", N_Nu));
	multi->GetYaxis()->SetNdivisions(2);
	multi->GetXaxis()->SetTitle("Nu [GeV]");
	multi->GetYaxis()->SetTitle("dE [MeV]"); 

	canvas->BuildLegend();
	canvas->SaveAs(Form("output/ZEloss_"+Nuclei_Type+"_%dnubins_%dentries_%dEcut%d_cheb%d_Ebins%d.pdf", N_Nu, nentries, int(E_min), int(E_max), n, nbins));
	canvas->SaveAs(Form("output/ZEloss_"+Nuclei_Type+"_%dnubins_%dentries_%dEcut%d_cheb%d_Ebins%d.C", N_Nu, nentries, int(E_min), int(E_max), n, nbins));
	//canvas->Write();

	std::cout<<" ABOUT TO CLOSE " << std::endl;

	fout->Close();
	file->Close();
	std::cout<< " BYE BYE " << std::endl;
	return 0;
}
