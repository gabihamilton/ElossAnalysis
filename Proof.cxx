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
const double E_max = 100;//2.0;       // Minimum Energy
const double E_min = 0;//0.5;       // Maximum Energy
const double limit_xf = 0.1;    // xF cut
const int nshift_E = 10;      // Number of shifts in Energy
const double step_E = 1.0;//1000.0; // Size of Shifts in Energy
//const int nbins = 100;        // Number of Energy bins

TRandom3* r0 = new TRandom3();
TRandom3* r1 = new TRandom3();

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

  
  Int_t nbins = atoi(argv[1]);
  Int_t N = atoi(argv[2]); // number of events

  //double Nu_min = 3.2 + Nu_bin*((4.2-3.2)/N_Nu); 
  //double Nu_max = Nu_min + (4.2-3.2)/N_Nu;

  Int_t N_EXPS = 5;
  //Int_t N = 10000;
  double x, y, w;

  std::vector<double> dataD, dataS, weightD, weightS; // vectors to store our generated data
  dataD.reserve(N);          // allocate vector beforehand, just for speedup
  dataS.reserve(N);          // allocate vector beforehand, just for speedup
  weightD.reserve(N);
  weightS.reserve(N);


  //-----Creating output file-----//  
  TFile *fout = new TFile("output/KS1D_Proof.root", "RECREATE");

  //-----Histograms with Energy distribution-----//
  TH1F *D = new TH1F("D","D",nbins,E_min,E_max);
  D->Sumw2();
  std::map<int,TH1F*> histograms;
  for(int i = 0; i<= nshift_E; i++){
    histograms[i] = new TH1F(Form("Shift%d",i),Form("Shift%d",i),nbins,E_min,E_max);
    histograms[i]->Sumw2();
  }

  
  //-----Histograms for Weighted Energy Distributions-----//
  TH1F *DW = new TH1F("D","D",nbins,E_min,E_max);
  DW->Sumw2();
  std::map<int,TH1F*> histogramsW;
  for(int i = 0; i<= nshift_E; i++){
    histogramsW[i] = new TH1F(Form("ShiftW%d",i),Form("ShiftW%d",i),nbins,E_min,E_max);
    histogramsW[i]->Sumw2();
  }


  //-----Creating the Graphs for the Eloss Shift Values------//
  TGraph *gElossKS = new TGraph();  //  Graph for Eloss values for the KS test
  TGraph *gElossWKS = new TGraph(); //  Graph for Eloss values fot the Weighted KS test
  TGraph *gElossKSb = new TGraph(); //  Graph for Eloss values for the Binned KS test
  TGraph *gElossWKSb = new TGraph(); // Graph for Eloss values for the Binned Weighted KS test
  TGraph *gElossP = new TGraph();   
  TGraph *gElossA = new TGraph();

  double pKSunbinned, pWKSunbinned, pWKSunbinnedPy;
  for (int i = 0; i < N_EXPS; i++) { // Loop over different pairs of histograms – Should not be the same as N

    // ------- Generating the Distributions
    for (int j = 0; j < N; j++) {  // Loop to fill one pair of histograms – Over N events
      x = r0->Gaus(50, 5);       // Define the first Gaussian distribution with (mean, sigma)
      y = r1->Gaus(50, 5) - 5;     // Define a second Gaussian distribution
      w = 1;
      dataD.push_back(x);           // push the random number into the first dataset
      weightD.push_back(w);
      dataS.push_back(y);           // push the other random number into the second dataset
      weightS.push_back(w);
      D->Fill(x);
      DW->Fill(x,w);

      for (int shift = 0; shift <= nshift_E; ++shift){ 

        energy_shift = step_E*shift;      // is doing the right shift

        histograms[shift]->Fill(y+energy_shift);    // Binned KS
        histogramsW[shift]->Fill(y+energy_shift, w);  // Weighted Binned
      }  
    }
    /*
    TCanvas *c1 = new TCanvas();
    DW->Draw("AL");
    histograms[0]->Draw("same");
    c1->BuildLegend();
    c1->SaveAs("output/D_histo_proof.pdf");
    */

    TGraph *gpKS = new TGraph();        //Graph for the Unbinned KS Test
    TGraph *gpWKS = new TGraph();       //Graph for the Unbinned Weighted KS Test
    TGraph *gpKSb = new TGraph();       // Graph for the Binned KS Test
    TGraph *gpWKSb = new TGraph();      // Graph for the Binned Weighted KS Test
    TGraph *gP_WKS = new TGraph();    // Graph for the Unbinned Weighted KS PYTHON
    TGraph *gpKSAcc = new TGraph();

    // ------ Startng each shift for tests ------ //
    for (int shift = 0; shift < nshift_E; ++shift)
    {
      std::vector<double> shifted_data;
      energy_shift = step_E*shift;      // is doing the right shift

      // ------- Binned KS tests -------- //

      //-----Binned Kolmogorov-Smirnov Test-----//
      D->Scale(1.0/D->Integral());
      histograms[shift]->Scale(1.0/histograms[shift]->Integral());
      //double pCSbinned = D->Chi2Test(histograms[i], "NORM");
      double pKSbinned = D->KolmogorovTest(histograms[shift], "D");
      //gpKSb->SetPoint(i, i, -1*TMath::Log10(pKSbinned));
      gpKSb->SetPoint(shift, shift, pKSbinned);

      //-----Binned Weighted KS Test-----//
      DW->Scale(1.0/DW->Integral());
      histogramsW[shift]->Scale(1.0/histogramsW[shift]->Integral());
      double pWKSbinned = DW->KolmogorovTest(histogramsW[shift], "D");
      //gWpKSb->SetPoint(i, i, -1*TMath::Log10(WpKSbinned));
      gpWKSb->SetPoint(shift, shift, pWKSbinned); 


      // ------- Unbinned KS tests --------- //

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


      for (int s = 0; s < nS; ++s)
      {
        double E = dataS[s] + energy_shift;
        // cout << dataS[s] << "   " << E << endl;
        shifted_data.push_back(E);
      }
      
      //cout << dataS[0] << "     " << shifted_data[0] << endl; 

      //----------UNBINNED KOLMOGOROV TEST----------//
      double pKS = TMath::KolmogorovTest(nD, &dataD[0], nS, &shifted_data[0], "D");
      //gpKS->SetPoint(i, i, -1*TMath::Log10(pKS));
      gpKS->SetPoint(shift, shift, pKS);


      Double_t w1sum = std::accumulate(weightD.begin(), weightD.end(), 0);
      Double_t w2sum = std::accumulate(weightS.begin(), weightS.end(), 0);


      //WEIGHTED KOLMOGOROV TEST :: PYTHON IMPLEMENTATION
      int j1=0, j2=0;
      Double_t d=0;
      Double_t j1w=0., j2w=0., fn1=0., fn2=0.;
      
      while (j1<nD && j2<nS){

        Double_t d1 = dataD[j1];
        Double_t d2 = shifted_data[j2];
        Double_t w1 = weightD[idxD[j1]];
        Double_t w2 = weightS[idxS[j2]];
        
        //cout << nD << "   " << w1sum << endl;

        if (d1<=d2){
          j1+=1;
          j1w+=w1;
          fn1=(j1w)/w1sum;
        }
        if (d2<=d1){
          j2+=1;
          j2w+=w2;
          fn2=(j2w)/w2sum;
        }
        if (TMath::Abs(fn2-fn1)>d){
          d=TMath::Abs(fn2-fn1);
        }
      }

      Double_t Pz = d * TMath::Sqrt(nD*nS/(nD+nS));
      double P_WKS = TMath::KolmogorovProb(Pz);
      //gpKS->SetPoint(i, i, -1*TMath::Log10(pKS));
      gP_WKS->SetPoint(shift, shift, P_WKS);


      //WEIGHTED KOLMOGOROV TEST :: ROOT IMPLEMENTATION
      Double_t rdiff = 0;
      Double_t rdmax = 0;
      Bool_t ok = kFALSE;
      j1 = 0;
      j2 = 0;
      for (int l = 0; l < nD+nS; ++l){
        Double_t w1 = weightD[idxD[j1]]/w1sum;

        Double_t w2 = weightS[idxS[j2]]/w2sum;

        if (dataD[j1] < shifted_data[j2]) {
          rdiff -= w1;
          j1++;
          if (j1 >= nD) {ok = kTRUE; break;}
        } else if (dataD[j1] > shifted_data[j2]) {
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
          while(j2 < nS && shifted_data[j2] == x) {
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
        gpWKS->SetPoint(shift, shift, pWKS);
      }

    }

    std::cout<< "END LOOP OVER SHIFTS" << std::endl;

    //--------Statistical Tests Graphs--------//
    gpKS->SetName("pKS");
    gpWKS->SetName("pWKS");
    gpKSb->SetName("pKSb_");
    gpWKSb->SetName("pWKSb");
    gP_WKS->SetName("P_WKS");
    //gpKSAcc->SetName(Form("pKSAcc_"+Nuclei_Type+"_nubin%d", Nu_bin));


    gpKS->SetLineColor(4);   
    gpWKS->SetLineColor(1);
    gpKSb->SetLineColor(6); 
    gpWKSb->SetLineColor(2);
    gP_WKS->SetLineColor(3); //green
    //gpKSAcc->SetLineColor(28);

    gpKS->SetLineWidth(3);
    gpWKS->SetLineWidth(3);
    gpKSb->SetLineWidth(3);
    gpWKSb->SetLineWidth(3);
    gP_WKS->SetLineWidth(3);
    //gpKSAcc->SetLineWidth(3);

    gpKS->SetLineStyle(1);
    gpWKS->SetLineStyle(2);
    gpKSb->SetLineStyle(3);
    gpWKSb->SetLineStyle(4);
    gP_WKS->SetLineStyle(5);
    //gpKSAcc->SetLineWidth(3);

    gpKS->SetTitle("Unbinned KS");
    gpWKS->SetTitle("Weighted Unbinned KS");
    gpKSb->SetTitle("Binned KS");
    gpWKSb->SetTitle("Weighted Binned KS"); 
    gP_WKS->SetTitle("Python Weighted KS");
    //gpKSAcc->SetTitle("Acc. Corrected KS");

    TCanvas *canvas = new TCanvas();
    TMultiGraph *multi = new TMultiGraph();
    multi->Add(gpKS, "L");
    multi->Add(gpWKS,"L");
    multi->Add(gpKSb,"L");
    multi->Add(gpWKSb,"L");
    multi->Add(gP_WKS,"L");
    //multi->Add(gpKSAcc, "L");

    multi->Draw("AL");
    multi->SetMaximum(1.2);
    multi->SetTitle("Probability curve");
    multi->GetYaxis()->SetNdivisions(2);
    multi->GetXaxis()->SetTitle("dE [MeV]");
    multi->GetYaxis()->SetTitle("p_{0}"); //"-Log(p_{0})"
    
    canvas->BuildLegend();
    canvas->SaveAs(Form("output/Prob_Proof%d.pdf", i));

  
    //-----ELOSS HISTOGRAMS PLOTS-----//
    
    Double_t elossKS=0;
    Double_t elossWKS=0;
    Double_t elossKSb=0;
    Double_t elossWKSb=0;
    Double_t elossP=0;
    Double_t elossA=0;
    Int_t i_KS=0;
    Int_t i_WKS=0;
    Int_t i_KSb=0;
    Int_t i_WKSb=0;
    Int_t i_P=0;
    Int_t i_A=0;

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

      gP_WKS->GetPoint(i, x, y); 
      if(y>elossP){
        elossP=y;
        i_P = i;
      }

      gpKSAcc->GetPoint(i, x, y); 
      if(y>elossA){
        elossA=y;
        i_A = i;
      }
    }
    cout << "ELOSS VALUE FOR KS: " << i_KS << "   PROB: " << elossKS << endl;
    cout << "ELOSS VALUE FOR WKS: " << i_WKS << "   PROB: " << elossWKS << endl;
    cout << "ELOSS VALUE FOR KSb: " << i_KSb << "   PROB: " << elossKSb << endl;
    cout << "ELOSS VALUE FOR WKSb: " << i_WKSb << "   PROB: " << elossWKSb << endl;
    cout << "ELOSS VALUE FOR Python: " << i_P << "   PROB: " << elossP << endl;
    //cout << "ELOSS VALUE FOR Acc: " << i_A << "   PROB: " << elossA << endl;


    //-----Energy spectra distributions-----//
    
    fout->cd();

    histograms[i_KS]->SetName("KS_");
    histograms[i_KS]->Write();

    histogramsW[i_WKS]->SetName("WKS_");
    histogramsW[i_WKS]->Write();

    histograms[i_KSb]->SetName("KSb_");
    histograms[i_KSb]->Write();

    histogramsW[i_WKSb]->SetName("WKSb_");
    histogramsW[i_WKSb]->Write();

    histogramsW[i_P]->SetName("Python_");
    histogramsW[i_P]->Write();

    D->SetName("D_");
    D->Write();

    DW->SetName("Weighted_D_");
    DW->Write();
    
    //------- ELOSS GRAPHS -------//
    gElossKS->SetPoint(i, i, i_KS);
    gElossWKS->SetPoint(i, i, i_WKS);
    gElossKSb->SetPoint(i, i, i_KSb);
    gElossWKSb->SetPoint(i, i, i_WKSb);
    gElossP->SetPoint(i, i, i_P);
    //gElossA->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_A);
  }//  END OF LOOP OVER NU BINS


  fout->cd();
  
  //  ELOSS PLOTS  //
  gElossKS->SetName(Form("ElossKS_"+Nuclei_Type));
  gElossWKS->SetName(Form("ElossWKS_"+Nuclei_Type));
  gElossKSb->SetName(Form("ElossKSb_"+Nuclei_Type));
  gElossWKSb->SetName(Form("ElossWKS_"+Nuclei_Type));
  gElossP->SetName(Form("ElossP_"+Nuclei_Type));
  //gElossA->SetName(Form("ElossP_"+Nuclei_Type));

  gElossKS->Write();
  gElossWKS->Write();
  gElossKSb->Write();
  gElossWKSb->Write();
  gElossP->Write();
  //gElossA->Write();

  gElossKS->SetMarkerColor(4);    
  gElossWKS->SetMarkerColor(1);
  gElossKSb->SetMarkerColor(6);   
  gElossWKSb->SetMarkerColor(2); 
  gElossP->SetMarkerColor(3);
  //gElossA->SetMarkerColor(28); 

  gElossKS->SetMarkerStyle(20);
  gElossWKS->SetMarkerStyle(30);
  gElossKSb->SetMarkerStyle(22);
  gElossWKSb->SetMarkerStyle(23);
  gElossP->SetMarkerStyle(42);
  //gElossA->SetMarkerStyle(21);

  gElossKS->SetTitle("Unbinned KS");
  gElossWKS->SetTitle("Weighted Unbinned KS");
  gElossKSb->SetTitle("Binned KS");
  gElossWKSb->SetTitle("Weighted Binned KS");
  gElossP->SetTitle("Python Weighted");
  //gElossA->SetTitle("Acc Corr KS");

  TCanvas *canvas = new TCanvas();
  TMultiGraph *multi = new TMultiGraph();
  multi->Add(gElossKS, "");
  multi->Add(gElossWKS,"");
  multi->Add(gElossKSb,"");
  multi->Add(gElossWKSb, "");
  multi->Add(gElossP, "");
  //multi->Add(gElossA, "");

  multi->Draw("AP");
  multi->Write();
  multi->SetMaximum(100);
  multi->SetTitle("Energy Loss Proof");
  multi->GetYaxis()->SetNdivisions(2);
  multi->GetXaxis()->SetTitle("Nu [GeV]");
  multi->GetYaxis()->SetTitle("dE [MeV]"); 

  canvas->BuildLegend();
  canvas->SaveAs("output/Eloss_Proof.pdf");

  std::cout<<" ABOUT TO CLOSE " << std::endl;
  fout->Close();
  std::cout<< " BYE BYE " << std::endl;
  return 1;
}
