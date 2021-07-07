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
const double E_max = 2500;//2.0;       // Minimum Energy     MeV
const double E_min = 0;//0.5;       // Maximum Energy
const double limit_xf = 0.1;    // xF cut
const int nshift_E = 100;      // Number of shifts in Energy
const double step_E = 2.5;//1000.0; // Size of Shifts in Energy in MeV
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

double Efficiency(double x, double x0, double a, double b)
{
  return a / (1 + TMath::Exp(-b*(x-x0)));
}

int main(int argc, char *argv[]){

  
  Int_t nbins = atoi(argv[1]);
  Int_t N = atoi(argv[2]); // number of events

  //double Nu_min = 3.2 + Nu_bin*((4.2-3.2)/N_Nu); 
  //double Nu_max = Nu_min + (4.2-3.2)/N_Nu;

  Int_t N_EXPS = 10;
  //Int_t N = 10000;
  double x, y, w, p;

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
  TH1F *DW = new TH1F("DW","DW",nbins,E_min,E_max);
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



  double x0 = 300;
  double a = 1;
  double b = 0.02;


  double pKSunbinned, pWKSunbinned, pWKSunbinnedPy;
  for (int i = 0; i < N_EXPS; i++) { // Loop over different pairs of histograms – Should not be the same as N
    r0->SetSeed(624+i);
    r1->SetSeed(624*2+i);
    dataD.clear();
    dataS.clear();
    weightD.clear();
    weightS.clear();


    // ------- Generating the Distributions
    for (int j = 0; j < N; j++) {  // Loop to fill one pair of histograms – Over N events
      //x = r0->Gaus(750, 250);       // Define the first Gaussian distribution with (mean, sigma)
      //y = r1->Gaus(750, 250) - 20;     // Define a second Gaussian distribution
      x = r0->Landau(500, 200);       // Define the first Beta distribution with (a, b)
      y = r1->Landau(500-20, 200);     // Defie a second Beta distribution
      //w = 1;
      dataD.push_back(x);           // push the random number into the first dataset
      weightD.push_back(1/(a / (1 + TMath::Exp(-b*(x-x0)))));
      //weightD.push_back(1);
      //cout << a << "   " << b << "  " << "   " << x << "   " << a/(1 + TMath::Exp(-b*(x-x0))) << endl;
      dataS.push_back(y);           // push the other random number into the second dataset
      weightS.push_back(1/(a / (1 + TMath::Exp(-b*(y-x0)))));
      //weightS.push_back(1);
      D->Fill(x);
      DW->Fill(x,1/(a/(1 + TMath::Exp(-b*(x-x0)))));
      //DW->Fill(x, 1);

      for (int shift = 0; shift <= nshift_E; ++shift){ 

        energy_shift = step_E*shift;      // is doing the right shift

        histograms[shift]->Fill(y+energy_shift);    // Binned KS
        histogramsW[shift]->Fill(y+energy_shift, 1/(a / (1 + TMath::Exp(-b*(y+energy_shift-x0)))));  // Weighted Binned
        //histogramsW[shift]->Fill(y+energy_shift, 1);
      }  
      energy_shift = 0;
    }
    
    TCanvas *c1 = new TCanvas();
    DW->Draw("AL");
    D->SetLineColor(2);
    D->Draw("same");
    //histograms[0]->Draw("same");
    c1->BuildLegend();
    c1->SaveAs("output/D_histo_proof.pdf");
    

    TGraph *gpKS = new TGraph();        //Graph for the Unbinned KS Test
    TGraph *gpWKS = new TGraph();       //Graph for the Unbinned Weighted KS Test
    TGraph *gpKSb = new TGraph();       // Graph for the Binned KS Test
    TGraph *gpWKSb = new TGraph();      // Graph for the Binned Weighted KS Test
    TGraph *gP_WKS = new TGraph();    // Graph for the Unbinned Weighted KS PYTHON
    TGraph *gpKSAcc = new TGraph();

    TGraph *lpKS = new TGraph();        //Graph for the Unbinned KS Test
    TGraph *lpWKS = new TGraph();       //Graph for the Unbinned Weighted KS Test
    TGraph *lpKSb = new TGraph();       // Graph for the Binned KS Test
    TGraph *lpWKSb = new TGraph();      // Graph for the Binned Weighted KS Test
    TGraph *lP_WKS = new TGraph();    // Graph for the Unbinned Weighted KS PYTHON

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
      lpKSb->SetPoint(shift, energy_shift, -1*TMath::Log10(pKSbinned));
      gpKSb->SetPoint(shift, energy_shift, pKSbinned);

      //-----Binned Weighted KS Test-----//
      DW->Scale(1.0/DW->Integral());
      histogramsW[shift]->Scale(1.0/histogramsW[shift]->Integral());
      double pWKSbinned = DW->KolmogorovTest(histogramsW[shift], "D");
      lpWKSb->SetPoint(shift, energy_shift, -1*TMath::Log10(pWKSbinned));
      gpWKSb->SetPoint(shift, energy_shift, pWKSbinned); 


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
      lpKS->SetPoint(shift, energy_shift, -1*TMath::Log10(pKS));
      gpKS->SetPoint(shift, energy_shift, pKS);


      //Double_t w1sum = std::accumulate(weightD.begin(), weightD.end(), 0);
      //Double_t w2sum = std::accumulate(weightS.begin(), weightS.end(), 0);
      double w1sum=0, w2sum=0;

      for (int m = 0; m < nS; ++m)
      {
        w1sum = w1sum + weightD[m];
        w2sum = w2sum + weightS[m];
      }


      //WEIGHTED KOLMOGOROV TEST :: PYTHON IMPLEMENTATION
      int j1=0, j2=0;
      Double_t d=0;
      Double_t j1w=0., j2w=0., fn1=0., fn2=0.;
/*      
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
      lP_WKS->SetPoint(shift, energy_shift, -1*TMath::Log10(P_WKS));
      gP_WKS->SetPoint(shift, energy_shift, P_WKS);
*/

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
        lpWKS->SetPoint(shift, energy_shift, -1*TMath::Log10(pWKS));
        gpWKS->SetPoint(shift, energy_shift, pWKS);
      }
    }

    std::cout<< "END LOOP OVER SHIFTS" << std::endl;

    //--------Statistical Tests Graph--------//
    gpKS->SetName("pKS");
    gpWKS->SetName("pWKS");
    gpKSb->SetName("pKSb_");
    gpWKSb->SetName("pWKSb");

    gpKS->SetLineColor(4);   
    gpWKS->SetLineColor(1);
    gpKSb->SetLineColor(6); 
    gpWKSb->SetLineColor(2);

    gpKS->SetLineWidth(3);
    gpWKS->SetLineWidth(3);
    gpKSb->SetLineWidth(3);
    gpWKSb->SetLineWidth(3);

    gpKS->SetLineStyle(1);
    gpWKS->SetLineStyle(2);
    gpKSb->SetLineStyle(3);
    gpWKSb->SetLineStyle(4);

    gpKS->SetTitle("Unbinned KS");
    gpWKS->SetTitle("Weighted Unbinned KS");
    gpKSb->SetTitle("Binned KS");
    gpWKSb->SetTitle("Weighted Binned KS");

    TCanvas *canvas = new TCanvas();
    TMultiGraph *multi = new TMultiGraph();
    multi->Add(gpKS, "L");
    multi->Add(gpWKS,"L");
    multi->Add(gpKSb,"L");
    multi->Add(gpWKSb,"L");

    multi->Draw("AL");
    multi->SetMaximum(1.2);
    multi->SetTitle("Probability curve");
    multi->GetYaxis()->SetNdivisions(2);
    multi->GetXaxis()->SetTitle("Shift");
    multi->GetYaxis()->SetTitle("p_{0}"); //"-Log(p_{0})"
    
    canvas->BuildLegend();
    canvas->SaveAs(Form("output/Prob_Proof%d_step%d_%dEbins_%devents.pdf", i, int(step_E*100), nbins, N));


    //--------Statistical Tests Graph LOG--------//
    lpKS->SetName("pKS");
    lpWKS->SetName("pWKS");
    lpKSb->SetName("pKSb_");
    lpWKSb->SetName("pWKSb");

    lpKS->SetLineColor(4);   
    lpWKS->SetLineColor(1);
    lpKSb->SetLineColor(6); 
    lpWKSb->SetLineColor(2);

    lpKS->SetLineWidth(3);
    lpWKS->SetLineWidth(3);
    lpKSb->SetLineWidth(3);
    lpWKSb->SetLineWidth(3);

    lpKS->SetLineStyle(1);
    lpWKS->SetLineStyle(2);
    lpKSb->SetLineStyle(3);
    lpWKSb->SetLineStyle(4);

    lpKS->SetTitle("Unbinned KS");
    lpWKS->SetTitle("Weighted Unbinned KS");
    lpKSb->SetTitle("Binned KS");
    lpWKSb->SetTitle("Weighted Binned KS");

    TCanvas *vas = new TCanvas();
    TMultiGraph *ulti = new TMultiGraph();
    ulti->Add(lpKS, "L");
    ulti->Add(lpWKS,"L");
    ulti->Add(lpKSb,"L");
    ulti->Add(lpWKSb,"L");

    ulti->Draw("AL");
    ulti->SetMaximum(10);
    ulti->SetTitle("Probability curve");
    ulti->GetYaxis()->SetNdivisions(2);
    ulti->GetXaxis()->SetTitle("Shift");
    ulti->GetYaxis()->SetTitle("-Log(p_{0})"); 
    
    vas->BuildLegend();
    vas->SaveAs(Form("output/LOGProb_Proof%d_step%d_%dEbins_%devents.pdf", i, int(step_E*100), nbins, N));



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
    cout << "ELOSS VALUE FOR KS: " << i_KS*step_E << "   PROB: " << elossKS << endl;
    cout << "ELOSS VALUE FOR WKS: " << i_WKS*step_E << "   PROB: " << elossWKS << endl;
    cout << "ELOSS VALUE FOR KSb: " << i_KSb*step_E << "   PROB: " << elossKSb << endl;
    cout << "ELOSS VALUE FOR WKSb: " << i_WKSb*step_E << "   PROB: " << elossWKSb << endl;


    //-----Energy spectra distributions-----//
    
    fout->cd();

    histograms[i_KS]->SetName(Form("KS_"));
    histograms[i_KS]->Write();

    histogramsW[i_WKS]->SetName("WKS_");
    histogramsW[i_WKS]->Write();

    histograms[i_KSb]->SetName("KSb_");
    histograms[i_KSb]->Write();

    histogramsW[i_WKSb]->SetName("WKSb_");
    histogramsW[i_WKSb]->Write();

    D->SetName(Form("D_%d", i));
    D->Write();

    histograms[0]->SetName(Form("S_%d", i));
    histograms[0]->Write();

    DW->SetName(Form("Weighted_D_%d", i));
    DW->Write();
    
    //------- ELOSS GRAPHS -------//
    gElossKS->SetPoint(i, i, i_KS*step_E);
    gElossWKS->SetPoint(i, i, i_WKS*step_E);
    gElossKSb->SetPoint(i, i, i_KSb*step_E);
    gElossWKSb->SetPoint(i, i, i_WKSb*step_E);
  }//  END OF LOOP OVER NU BINS


  fout->cd();
  
  //  ELOSS PLOTS  //
  gElossKS->SetName(Form("ElossKS_"+Nuclei_Type));
  gElossWKS->SetName(Form("ElossWKS_"+Nuclei_Type));
  gElossKSb->SetName(Form("ElossKSb_"+Nuclei_Type));
  gElossWKSb->SetName(Form("ElossWKS_"+Nuclei_Type));

  gElossKS->Write();
  gElossWKS->Write();
  gElossKSb->Write();
  gElossWKSb->Write();

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

  TCanvas *canvas = new TCanvas();
  TMultiGraph *multi = new TMultiGraph();
  multi->Add(gElossKS, "");
  multi->Add(gElossWKS,"");
  multi->Add(gElossKSb,"");
  multi->Add(gElossWKSb, "");

  multi->Draw("AP");
  multi->Write();
  //multi->SetMaximum(100);
  multi->SetTitle("Energy Loss Proof");
  multi->GetYaxis()->SetNdivisions(2);
  multi->GetXaxis()->SetTitle("Experiment");
  multi->GetYaxis()->SetTitle("Selected Shift"); 

  canvas->BuildLegend();
  canvas->SaveAs(Form("output/Eloss_Proof%d_%dEbins_%devents.pdf", int(step_E*100), nbins, N) );

  std::cout<<" ABOUT TO CLOSE " << std::endl;
  fout->Close();
  std::cout<< " BYE BYE " << std::endl;
  return 1;
}
