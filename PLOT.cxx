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


#define N_EXPS 10 // the number of experiments we make to get the KS test distribution \
                    // I tested this, and I think 1000 is an ok number to consider

// Set this as global
TH1D* Hist_ProbKSb = new TH1D("Hist_ProbKSb", "Hist_ProbKSb", 100, 0.0, 1.0);
TH1D* Hist_ProbWKSb = new TH1D("Hist_ProbWKSb", "Hist_ProbWKSb", 100, 0.0, 1.0);
TH1D* Hist_ProbWKSbPy = new TH1D("Hist_ProbWKSbPy", "Hist_ProbWKSbPy", 100, 0.0, 1.0);

const double step_E = 2.5;//1000.0; // Size of Shifts in Energy in MeV

int main(int argc, char *argv[])
{
    std::cout << "PLOTS" << std::endl;

    Int_t nbins = atoi(argv[1]);
    Int_t N = atoi(argv[2]); // number of events

    TFile *file = new TFile(Form("output/KS1D_POC_%dkevents_%dEbins_%fstep.root", N/1000, nbins, step_E)); // Acceptance files


    for (int i = 0; i < N_EXPS; ++i)
    {
        // PROB GRAPHS
        TGraph *gpKS = (TGraph*) file->Get(Form("pKS_exp%d", i));         // Graph for the Unbinned KS Test
        TGraph *gpWKS = (TGraph*) file->Get(Form("pWKS_exp%d", i));       // Graph for the Unbinned Weighted KS Test
        TGraph *gpKSb = (TGraph*) file->Get(Form("pKSb_exp%d", i));       // Graph for the Binned KS Test
        TGraph *gpWKSb = (TGraph*) file->Get(Form("pWKSb_exp%d", i));     // Graph for the Binned Weighted KS Test

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
        canvas->SaveAs(Form("output/Prob_exp%d_step%d_%dEbins_%devents.pdf", i, int(step_E*100), nbins, N));

        // LOG GRAPHS
        TGraph *lpKS = (TGraph*) file->Get(Form("log_pKS_exp%d", i));         // Graph for the Unbinned KS Test
        TGraph *lpWKS = (TGraph*) file->Get(Form("log_pWKS_exp%d", i));       // Graph for the Unbinned Weighted KS Test
        TGraph *lpKSb = (TGraph*) file->Get(Form("log_pKSb_exp%d", i));       // Graph for the Binned KS Test
        TGraph *lpWKSb = (TGraph*) file->Get(Form("log_pWKSb_exp%d", i));     // Graph for the Binned Weighted KS Test

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


        //Distribution Plots

        TGraph *A_D = (TGraph*) file->Get(Form("D_distribution_exp%d", i));
        TGraph *W_D = (TGraph*) file->Get(Form("D_weight_exp%d", i));

        A_D->SetLineWidth(3);
        W_D->SetLineWidth(3);

        TCanvas *c1 = new TCanvas();
        A_D->Draw("AL");
        c1->BuildLegend();
        c1->SaveAs(Form("output/Array_D_exp%d.pdf", i));

        TCanvas *c2 = new TCanvas();
        W_D->Draw("AL");
        c2->BuildLegend();
        c2->SaveAs(Form("output/Weight_D_exp%d.pdf", i));

    }
    cout << "END OF LOOP OVER EXP" << endl;

    //ELoss PLOT

    TGraph *gElossKS = (TGraph*) file->Get("ElossKS_");  //  Graph for Eloss values for the KS test
    TGraph *gElossWKS = (TGraph*) file->Get("ElossWKS_"); //  Graph for Eloss values fot the Weighted KS test
    TGraph *gElossKSb = (TGraph*) file->Get("ElossKSb_"); //  Graph for Eloss values for the Binned KS test
    TGraph *gElossWKSb = (TGraph*) file->Get("ElossWKSb_"); // Graph for Eloss values for the Binned Weighted KS test

    cout << "OK" << endl;

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

    return 1;

}
