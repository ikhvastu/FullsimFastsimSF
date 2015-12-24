#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "TVector.h"
#include "TText.h"
#include "TCOlor.h"

#include "fstream"
#include "sstream"  
//using namespace std;

static double AVAL = 1;

const double pt_min = 10.;
    const double pt_max = 200.;

    const double eta_min = -2.5;
    const double eta_max = 2.5;

    const double HT_min = 0.;
    const double HT_max = 1000.;

    const double Njets_min = 0.;
    const double Njets_max = 9.;

    const double npv_min = 0.;
    const double npv_max = 40.;

    const double Activb_min = 0.002;
    const double Activb_max = 10.;
    const double Activb_bins = 100;

    const double ptratiobins = 50;
    const double ptratioMin = 0.;
    const double ptratioMax = 1.5;

namespace{
  const int bands = 999;
  int rainbow[bands];
  int patriotic[bands];
}

void showHistUP(double, double, TGraphAsymmErrors *[2], string, string, double, double, string, double, double, double, double);
void showHistDOWN(double, double, TGraphAsymmErrors *[2], string, string, double, double, double, double, string, double, double, double, double, double &p1 = AVAL);
TH2D TranslateHisto(const TH2 &);

void Print1D(TGraphAsymmErrors *gh[6][2], const TString &ext, const string &TitleY, const string &lepton, TLegend *);
void Print2D(TH2 const * const h_data_in, TH2 const * const h_mc_in, const TString &ext);
void Print3D(TH3 const * const h_full_in, TH3 const * const h_fast_in, const TString &ext);

void PrintTable(TH2 const * const histo, const TString &ext);
void PrintLine(ofstream &file, TH2 const * const histo, int bin, const TString &label);
void GetPatrioticPalette();
void GetRainbowPalette();
//void FillHisto(TH2D* [4], TH3D*, TH2D*, TH2D*, TH2D*, TH2D*, TH2D*, TH1D *[6], double, double, double, double, double, double);
//void FillHisto(TH2D* [4], TH3D*, TH2D*, TH1D *[6], double, double, double, double, double, double);
void FillHisto(TH2D* [4], TH3D*, TH2D*, TH2D*, TH2D*, TH1D *[6], double, double, double, double, double, double);

double weightCalc(double);

void tree()
{
    GetPatrioticPalette();
    GetRainbowPalette();
    gStyle->SetPaintTextFormat("4.3f");
    const int nLeptonsMax = 6;
    int fontToUse = 42;
    gStyle->SetOptFit(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetNdivisions(505,"XY");
    
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
    
    gStyle->SetLabelFont(fontToUse,"XYZ");
    gStyle->SetLabelSize(0.05,"XYZ");
    gStyle->SetLabelOffset(0.001,"X");
    
    gStyle->SetTitleFont(fontToUse,"");
    gStyle->SetTitleFont(fontToUse,"XYZ");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetTitleSize(0.06,"XY");
    //gStyle->SetTitleXOffset(1.0);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.2);
    
    gStyle->SetErrorX(0.);
    
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadLeftMargin(0.15);
    
    gStyle->SetStatFont(fontToUse);
    gStyle->SetStatColor(10);
    gStyle->SetStatFontSize(0.08);
    gStyle->SetTitleFont(fontToUse);
    gStyle->SetTitleFontSize(0.08);
    
    gStyle->SetMarkerSize(1.);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerColor(gStyle->GetHistLineColor());
    
    gStyle->SetPalette(1,0);
    
    gStyle->SetFuncColor(kRed);
    
    
    int _n_bJets;
    int _n_Jets;
    
    double _jetEta[20];
    double _jetPhi[20];
    double _jetPt[20];
    bool _bTagged[20];
    double _csv[20];
    
    
    //TClonesArray* _leptonP4;
    
    bool _ossf[3];
    double _mll[3];
    
    int _indeces[nLeptonsMax];
    int _flavors[nLeptonsMax];
    double _charges[nLeptonsMax];
    double _isolation[nLeptonsMax];
    double _miniisolation[nLeptonsMax][2];
    double _multiisolation[nLeptonsMax][5];
    double _isolationComponents[nLeptonsMax][4];
    double _isolationMC[nLeptonsMax][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    
    int _origin[nLeptonsMax];
    int _originReduced[nLeptonsMax];
    double _PVchi2;
    double _PVerr[3];
    double _ipPV[nLeptonsMax];
    double _ipPVerr[nLeptonsMax];
    
    double _ipZPV[nLeptonsMax];
    double _ipZPVerr[nLeptonsMax];
    
    double _3dIP[nLeptonsMax];
    double _3dIPerr[nLeptonsMax];
    double _3dIPsig[nLeptonsMax];
    
    double _mt[nLeptonsMax];
    
    double _closeJetPtAll[nLeptonsMax];
    
    double _closeJetAngAll[nLeptonsMax];
    double _ptRelAll[nLeptonsMax];
    double _ptrel[nLeptonsMax];
    double _closeJetPtAllMC[nLeptonsMax];
    double _closeJetPtAllstatus[nLeptonsMax];
    int _partonIdMatched[nLeptonsMax];
    bool _sameParton[nLeptonsMax];
    
    bool _isloose[nLeptonsMax];
    bool _ismedium[nLeptonsMax];
    bool _istight[nLeptonsMax];
    bool _istightID[nLeptonsMax];
    bool _isvetoIDCutBased[nLeptonsMax];
    bool _islooseIDCutBased[nLeptonsMax];
    bool _ismediumIDCutBased[nLeptonsMax];
    bool _istightIDCutBased[nLeptonsMax];
    
    double _mompt[nLeptonsMax];
    double _momphi[nLeptonsMax];
    double _mometa[nLeptonsMax];
    int _mompdg[nLeptonsMax];
    
    int _n_PV;

    ULong64_t _eventNb;
    ULong64_t _runNb;
    ULong64_t _lumiBlock;

    double RhoCN;
    double RhoAll;

    double _mvaValue[nLeptonsMax];
    bool _trigEmulator[nLeptonsMax];
    bool _isotrigEmulator[nLeptonsMax];
    bool _chargeConst[nLeptonsMax];

    int _hitsNumber[nLeptonsMax];

    double _activity[nLeptonsMax];

    bool _vtxFitConversion[nLeptonsMax];

    bool _lgenParticle[nLeptonsMax];

    bool _fromHardProcessFinalState[nLeptonsMax];
    bool _isPromptFinalState[nLeptonsMax];
    
    double _lPt[nLeptonsMax], _lEta[nLeptonsMax], _lPhi[nLeptonsMax], _lE[nLeptonsMax];
    double _lPtmc[nLeptonsMax], _lEtamc[nLeptonsMax], _lPhimc[nLeptonsMax], _lEmc[nLeptonsMax];
    double _lpdgmc[nLeptonsMax];
    bool _flag_lepton[nLeptonsMax];
    double _dR15[nLeptonsMax];
    double _nuPtmc[nLeptonsMax], _nuEtamc[nLeptonsMax], _nuPhimc[nLeptonsMax], _nuEmc[nLeptonsMax];
    
    double _mtmc[nLeptonsMax];

    
    double _met;
    double _met_phi;
    double HT;
    
    int _nLeptons;
    int _nEle;
    int _nMu;
    int _nTau;
    
    TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
   
    TString folderName = "data/";

    const int nSamples = 2;

    TFile *hfile[nSamples];
    TTree* inputTreeFake[nSamples];

    //TString fileList[nSamples] = {"dy50_fastsim.root", "dy50_fullsim_DB.root"};
    //TString fileList[nSamples] = {"ttbar_fastsim.root", "ttbar_fullsim.root"};
    //TString fileList[nSamples] = {"dy50_fastsim.root", "dy50_fullsim.root"};
    //TString fileList[nSamples] = {"ttbar_fullsim_cor.root", "dy50_fullsim_cor.root"};
    TString fileList[nSamples] = {"ttbar_fullsim_cor.root", "test.root"};

    TString sampleName = "test2";
    //TString sampleName = "DY_dR_onZreq";
        
    TString names[nSamples] = {"Fastsim", "Fullsim"};
    
    Color_t colorH[nSamples] = {kBlack, kBlue};
    
    const int numberBinspt2Del = 6;
    const int numberBinspt2Dmu = 6;
    double ptb2Del[numberBinspt2Del] = {10., 20., 30.,40., 50., 200.}; 
    double ptb2Dmu[numberBinspt2Dmu] = {10., 20., 30.,40., 50., 200.}; 

    const int numberBinseta2D = 4;
    const int numberBinseta2Dmu = 5;
    double etab6mu[numberBinseta2Dmu] = {0., 0.9  , 1.2  , 2.1, 2.4};
    double etab6el[numberBinseta2D] = {0., 1.442, 1.556, 2.5};

    const int numberBinsNPV2D = 5;
    double npvb6[numberBinsNPV2D] = {1., 9., 13., 17., 40.};


    TH1D* h_ele_npv_pt[5][nSamples];
    TH1D* h_ele_npv_ptc[5][nSamples];

    TH1D* h_mu_npv_pt[5][nSamples];
    TH1D* h_mu_npv_ptc[5][nSamples];

    /*
    TH1D* h_ptratio_good[nSamples];
    TH1D* h_ptrel_good[nSamples];
    TH1D* h_ptratio_bad[nSamples];
    TH1D* h_ptrel_bad[nSamples];

    TH1D* h_ptratio_good_mu[nSamples];
    TH1D* h_ptrel_good_mu[nSamples];
    TH1D* h_ptratio_bad_mu[nSamples];
    TH1D* h_ptrel_bad_mu[nSamples];

    TH1D* h_RhoNC[nSamples];
    TH1D* h_RhoAll[nSamples];
    */

    TH2D* h_ele_activHT[15][nSamples];
    TH2D* h_mu_activHT[8][nSamples];

    TH2D* h_ele_activnpv[15][nSamples];
    TH2D* h_mu_activnpv[8][nSamples];

    TH2D* h_ele_HTnpv[15][nSamples];
    TH2D* h_mu_HTnpv[8][nSamples];

    const int numberBinsHT2D = 7;
    double HTb2D[numberBinsHT2D] = {0., 50., 100., 150., 200., 400., 1000.};

    const int numberBinsactiv2D = 6;
    double activb2D[numberBinsactiv2D] = {0.002, 0.02, 0.05, 0.15, 1., 10.};

    /*
    const int numberBinsnpv2D = 5;
    double npvb2D[numberBinsnpv2D] = {1., 13., 17., 21., 40.};
    */

    TH2D* h_ele_pteta[15][nSamples];
    TH2D* h_mu_pteta[8][nSamples];

    TH2D* h_ele_pteta_2Dnpv[15][nSamples][4];
    TH2D* h_mu_pteta_2Dnpv[8][nSamples][4];

    TH3D* h_ele_pt_eta_npv[15][nSamples];
    TH3D* h_mu_pt_eta_npv[8][nSamples];
    
    TH2D* h_ele_ptactiv[15][nSamples];
    TH2D* h_mu_ptactiv[8][nSamples];
    
    TH2D* h_ele_etaactiv[15][nSamples];
    TH2D* h_mu_etaactiv[8][nSamples];
    
   
    /*
    
    TH2D* h_ele_ptactivc[nSamples];

    
    TH2D* h_mu_ptactivc[nSamples];
    */
    TH1D* h_test[4][nSamples];

    TH1D* histoele[15][nSamples][6];
    TH1D* histomu[8][nSamples][6];

    TString nameLep[3] = {"electron ", "muon ", "tau "};
    TString xaxisvar[6] = {"p_{T}, [GeV]", "#eta", "H_{T}, [GeV]", "N_{jets}", "# PV", "Activity"};
    const int numberBins[6] = {11, 11, 7, 7, 11, 6};
   
    double ptb[11] = {10., 20., 30.,40.,50., 60., 70., 80., 90., 100., 200.}; 
    double etab[11] = {-2.5, -2., -1.5, -1., -0.5, 0., 0.5, 1.0, 1.5, 2., 2.5};
    double htb[7] = {0., 50., 100., 150., 200., 400., 1000.};
    double njetsb[7] = {0., 1., 2., 3., 4., 6., 9.};
    double npvb[11] = {1., 5., 9., 13., 17., 21., 25., 29., 33., 37., 40.};
    double activb[6] = {0.002, 0.02, 0.05, 0.15, 1., 10.};

    double *Bins[6] = {ptb, etab, htb, njetsb, npvb, activb};

    TGraphAsymmErrors *graph[21][6][nSamples];
       
    TGraphAsymmErrors *gh_ele_npv_ptc[5][nSamples];
    TGraphAsymmErrors *gh_mu_npv_ptc[5][nSamples];

    TGraphAsymmErrors *gh_mu_HT_ptc[5][nSamples];
    TGraphAsymmErrors *gh_mu_activ_ptc[5][nSamples];

    TGraphAsymmErrors *g_test[nSamples];

    TH1D* h_npv_distr[nSamples];
    TH1D* h_activ_distr[3][2][nSamples];

    TH1D* h_mu_HT_pt[5][nSamples];
    TH1D* h_mu_HT_ptc[5][nSamples];

    TH1D* h_mu_activ_pt[5][nSamples];
    TH1D* h_mu_activ_ptc[5][nSamples];

    TH1D* h_activityHT[5][2][nSamples];
    double activ_max = 4.;
    double activ_min = 0.;
    int activbins = 40;

    TH1D* h_deltaR[3][nSamples];

    TH1D* h_deltaPT[4][nSamples];

    TH2D* h_ptpt[nSamples];
    TH1D* h_etaeta[nSamples];
    TH2D* h_activactiv[nSamples];


 for(int sample = 0;sample < nSamples;sample++){
        TString name;
        /*
        name = Form("pt_%d", sample);
        h_ptpt[sample] = new TH2D(name, "ptpt "+names[sample], 19, 10., 200., 19, 10., 200.);
        h_ptpt[sample]->SetLineColor(colorH[sample]);
        h_ptpt[sample]->SetFillColor(colorH[sample]);
        h_ptpt[sample]->SetMarkerColor(colorH[sample]);
        h_ptpt[sample]->SetFillStyle(3000);
        h_ptpt[sample]->Sumw2();

        name = Form("pt_%d", sample);
        h_etaeta[sample] = new TH2D(name, "ptpt "+names[sample], 100, 0., 10., 100, 0., 10.);
        h_etaeta[sample]->SetLineColor(colorH[sample]);
        h_etaeta[sample]->SetFillColor(colorH[sample]);
        h_etaeta[sample]->SetMarkerColor(colorH[sample]);
        h_etaeta[sample]->SetFillStyle(3000);
        h_etaeta[sample]->Sumw2();

        name = Form("activactiv_%d", sample);
        h_activactiv[sample] = new TH2D(name, "activactiv "+names[sample], 100, 0., 1., 100, 0., 1.);
        h_activactiv[sample]->SetLineColor(colorH[sample]);
        h_activactiv[sample]->SetFillColor(colorH[sample]);
        h_activactiv[sample]->SetMarkerColor(colorH[sample]);
        h_activactiv[sample]->SetFillStyle(3000);
        h_activactiv[sample]->Sumw2();
        */
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 3; j++){
                name = Form("activ_distr_%d_%d_%d", sample,i,j);
                h_activ_distr[j][i][sample] = new TH1D(name, "activ_distr "+names[sample]+";activity; a.u.", 50, 0., 0.5);
                h_activ_distr[j][i][sample]->SetLineColor(colorH[sample]);
                h_activ_distr[j][i][sample]->SetFillColor(colorH[sample]);
                h_activ_distr[j][i][sample]->SetMarkerColor(colorH[sample]);
                h_activ_distr[j][i][sample]->SetFillStyle(3000);
                h_activ_distr[j][i][sample]->Sumw2();
            }

            
        }

        for(int i = 0; i < 3; i++){
            name = Form("deltaR_%d_%d", sample,i);
            h_deltaR[i][sample] = new TH1D(name, "deltaR min "+names[sample]+";deltaR; a.u.", 50, 0., 1.);
            h_deltaR[i][sample]->SetLineColor(colorH[sample]);
            h_deltaR[i][sample]->SetFillColor(colorH[sample]);
            h_deltaR[i][sample]->SetMarkerColor(colorH[sample]);
            h_deltaR[i][sample]->SetFillStyle(3000);
            h_deltaR[i][sample]->Sumw2();
        }

        for(int i = 0; i < 4; i++){
            name = Form("deltaPT_%d_%d", sample,i);
            h_deltaPT[i][sample] = new TH1D(name, "deltaPT min "+names[sample]+";deltaPT; a.u.", 100, 0., 10.);
            h_deltaPT[i][sample]->SetLineColor(colorH[sample]);
            h_deltaPT[i][sample]->SetFillColor(colorH[sample]);
            h_deltaPT[i][sample]->SetMarkerColor(colorH[sample]);
            h_deltaPT[i][sample]->SetFillStyle(3000);
            h_deltaPT[i][sample]->Sumw2();
        }

        for(int i = 0; i < 4; i++){
            name = Form("test_%d_%d", sample, i);
            h_test[i][sample] = new TH1D(name, "test "+names[sample]+";#eta; events", numberBins[1]-1, Bins[1]);
            h_test[i][sample]->SetLineColor(colorH[sample]);
            h_test[i][sample]->SetFillColor(colorH[sample]);
            h_test[i][sample]->SetMarkerColor(colorH[sample]);
            h_test[i][sample]->SetFillStyle(3000);
            h_test[i][sample]->Sumw2();
        }
        
        for(int i = 0; i < 5; i++){
            /*
            name = Form("npv_pt_el_%d_%d", sample, i);
            h_ele_npv_pt[i][sample] = new TH1D(name, "npv_pt_el "+names[sample]+";NPV; events", numberBins[4], npvb[0], npvb[numberBins[4]-1]);
            h_ele_npv_pt[i][sample]->SetLineColor(colorH[sample]);
            h_ele_npv_pt[i][sample]->SetFillColor(colorH[sample]);
            h_ele_npv_pt[i][sample]->SetMarkerColor(colorH[sample]);
            h_ele_npv_pt[i][sample]->SetFillStyle(3000);
            h_ele_npv_pt[i][sample]->Sumw2();

            name = Form("npv_pt_elc_%d_%d", sample, i);
            h_ele_npv_ptc[i][sample] = new TH1D(name, "npv_pt_elc "+names[sample]+";NPV; events", numberBins[4], npvb[0], npvb[numberBins[4]-1]);
            h_ele_npv_ptc[i][sample]->SetLineColor(colorH[sample]);
            h_ele_npv_ptc[i][sample]->SetFillColor(colorH[sample]);
            h_ele_npv_ptc[i][sample]->SetMarkerColor(colorH[sample]);
            h_ele_npv_ptc[i][sample]->SetFillStyle(3000);
            h_ele_npv_ptc[i][sample]->Sumw2();
            

            name = Form("npv_pt_mu_%d_%d", sample, i);
            h_mu_npv_pt[i][sample] = new TH1D(name, "npv_pt_mu "+names[sample]+";NPV; events", numberBins[4], npvb[0], npvb[numberBins[4]-1]);
            h_mu_npv_pt[i][sample]->SetLineColor(colorH[sample]);
            h_mu_npv_pt[i][sample]->SetFillColor(colorH[sample]);
            h_mu_npv_pt[i][sample]->SetMarkerColor(colorH[sample]);
            h_mu_npv_pt[i][sample]->SetFillStyle(3000);
            h_mu_npv_pt[i][sample]->Sumw2();

            name = Form("npv_pt_muc_%d_%d", sample, i);
            h_mu_npv_ptc[i][sample] = new TH1D(name, "npv_pt_muc "+names[sample]+";NPV; events", numberBins[4], npvb[0], npvb[numberBins[4]-1]);
            h_mu_npv_ptc[i][sample]->SetLineColor(colorH[sample]);
            h_mu_npv_ptc[i][sample]->SetFillColor(colorH[sample]);
            h_mu_npv_ptc[i][sample]->SetMarkerColor(colorH[sample]);
            h_mu_npv_ptc[i][sample]->SetFillStyle(3000);
            h_mu_npv_ptc[i][sample]->Sumw2();
            */

            name = Form("HT_pt_mu_%d_%d", sample, i);
            h_mu_HT_pt[i][sample] = new TH1D(name, "HT_pt_mu "+names[sample]+";H_{T}; events", numberBins[2], htb[0], htb[numberBins[2]-1]);
            h_mu_HT_pt[i][sample]->SetLineColor(colorH[sample]);
            h_mu_HT_pt[i][sample]->SetFillColor(colorH[sample]);
            h_mu_HT_pt[i][sample]->SetMarkerColor(colorH[sample]);
            h_mu_HT_pt[i][sample]->SetFillStyle(3000);
            h_mu_HT_pt[i][sample]->Sumw2();

            name = Form("HT_pt_muc_%d_%d", sample, i);
            h_mu_HT_ptc[i][sample] = new TH1D(name, "HT_pt_muc "+names[sample]+";H_{T}; events", numberBins[2], htb[0], htb[numberBins[2]-1]);
            h_mu_HT_ptc[i][sample]->SetLineColor(colorH[sample]);
            h_mu_HT_ptc[i][sample]->SetFillColor(colorH[sample]);
            h_mu_HT_ptc[i][sample]->SetMarkerColor(colorH[sample]);
            h_mu_HT_ptc[i][sample]->SetFillStyle(3000);
            h_mu_HT_ptc[i][sample]->Sumw2();

            name = Form("activ_pt_mu_%d_%d", sample, i);
            h_mu_activ_pt[i][sample] = new TH1D(name, "activ_pt_mu "+names[sample]+";activity; events", numberBins[5]-1, Bins[5]);
            h_mu_activ_pt[i][sample]->SetLineColor(colorH[sample]);
            h_mu_activ_pt[i][sample]->SetFillColor(colorH[sample]);
            h_mu_activ_pt[i][sample]->SetMarkerColor(colorH[sample]);
            h_mu_activ_pt[i][sample]->SetFillStyle(3000);
            h_mu_activ_pt[i][sample]->Sumw2();

            name = Form("activ_pt_muc_%d_%d", sample, i);
            h_mu_activ_ptc[i][sample] = new TH1D(name, "activ_pt_muc "+names[sample]+";activity; events", numberBins[5]-1, Bins[5]);
            h_mu_activ_ptc[i][sample]->SetLineColor(colorH[sample]);
            h_mu_activ_ptc[i][sample]->SetFillColor(colorH[sample]);
            h_mu_activ_ptc[i][sample]->SetMarkerColor(colorH[sample]);
            h_mu_activ_ptc[i][sample]->SetFillStyle(3000);
            h_mu_activ_ptc[i][sample]->Sumw2();

            for(int j = 0; j < 2; j++){

                name = Form("h_activityHT_%d_%d_%d", sample, i, j);
                h_activityHT[i][j][sample] = new TH1D(name, "h_activityHT "+names[sample]+";activity; events", activbins, activ_min, activ_max);
                h_activityHT[i][j][sample]->SetLineColor(colorH[sample]);
                h_activityHT[i][j][sample]->SetFillColor(colorH[sample]);
                h_activityHT[i][j][sample]->SetMarkerColor(colorH[sample]);
                h_activityHT[i][j][sample]->SetFillStyle(3000);
                h_activityHT[i][j][sample]->Sumw2();

            }

        }   
        /*
        name = Form("h_RhoNC_%d", sample);
        h_RhoNC[sample] = new TH1D(name, "h_RhoNC "+names[sample]+";h_RhoNC; events", 25, 0., 50.);
        h_RhoNC[sample]->SetLineColor(colorH[sample]);
        h_RhoNC[sample]->SetFillColor(colorH[sample]);
        h_RhoNC[sample]->SetMarkerColor(colorH[sample]);
        h_RhoNC[sample]->SetFillStyle(3000);
        h_RhoNC[sample]->Sumw2();

        name = Form("h_RhoAll_%d", sample);
        h_RhoAll[sample] = new TH1D(name, "h_RhoAll "+names[sample]+";h_RhoAll; events", 25, 0., 50.);
        h_RhoAll[sample]->SetLineColor(colorH[sample]);
        h_RhoAll[sample]->SetFillColor(colorH[sample]);
        h_RhoAll[sample]->SetMarkerColor(colorH[sample]);
        h_RhoAll[sample]->SetFillStyle(3000);
        h_RhoAll[sample]->Sumw2();

        name = Form("ptratio_good_%d", sample);
        h_ptratio_good[sample] = new TH1D(name, "ptratio_good "+names[sample]+";ptratio; events", ptratiobins, ptratioMin, ptratioMax);
        h_ptratio_good[sample]->SetLineColor(colorH[sample]);
        h_ptratio_good[sample]->SetFillColor(colorH[sample]);
        h_ptratio_good[sample]->SetMarkerColor(colorH[sample]);
        h_ptratio_good[sample]->SetFillStyle(3000);
        h_ptratio_good[sample]->Sumw2();
        
        name = Form("ptrel_good_%d", sample);
        h_ptrel_good[sample] = new TH1D(name, "ptrel_good "+names[sample]+";ptrel[GeV]; events", 50, 0., 50.);
        h_ptrel_good[sample]->SetLineColor(colorH[sample]);
        h_ptrel_good[sample]->SetFillColor(colorH[sample]);
        h_ptrel_good[sample]->SetMarkerColor(colorH[sample]);
        h_ptrel_good[sample]->SetFillStyle(3000);
        h_ptrel_good[sample]->Sumw2();
        
        name = Form("ptratio_bad_%d", sample);
        h_ptratio_bad[sample] = new TH1D(name, "ptratio_bad "+names[sample]+";ptratio; events", ptratiobins, ptratioMin, ptratioMax);
        h_ptratio_bad[sample]->SetLineColor(colorH[sample]);
        h_ptratio_bad[sample]->SetFillColor(colorH[sample]);
        h_ptratio_bad[sample]->SetMarkerColor(colorH[sample]);
        h_ptratio_bad[sample]->SetFillStyle(3000);
        h_ptratio_bad[sample]->Sumw2();
        
        name = Form("ptrel_bad_%d", sample);
        h_ptrel_bad[sample] = new TH1D(name, "ptrel_bad "+names[sample]+";ptrel[GeV]; events", 50, 0., 50.);
        h_ptrel_bad[sample]->SetLineColor(colorH[sample]);
        h_ptrel_bad[sample]->SetFillColor(colorH[sample]);
        h_ptrel_bad[sample]->SetMarkerColor(colorH[sample]);
        h_ptrel_bad[sample]->SetFillStyle(3000);
        h_ptrel_bad[sample]->Sumw2();

        name = Form("ptratio_good_mu_%d", sample);
        h_ptratio_good_mu[sample] = new TH1D(name, "ptratio_good_mu "+names[sample]+";ptratio; events", ptratiobins, ptratioMin, ptratioMax);
        h_ptratio_good_mu[sample]->SetLineColor(colorH[sample]);
        h_ptratio_good_mu[sample]->SetFillColor(colorH[sample]);
        h_ptratio_good_mu[sample]->SetMarkerColor(colorH[sample]);
        h_ptratio_good_mu[sample]->SetFillStyle(3000);
        h_ptratio_good_mu[sample]->Sumw2();
        
        name = Form("ptrel_good_mu_%d", sample);
        h_ptrel_good_mu[sample] = new TH1D(name, "ptrel_good_mu "+names[sample]+";ptrel[GeV]; events", 50, 0., 50.);
        h_ptrel_good_mu[sample]->SetLineColor(colorH[sample]);
        h_ptrel_good_mu[sample]->SetFillColor(colorH[sample]);
        h_ptrel_good_mu[sample]->SetMarkerColor(colorH[sample]);
        h_ptrel_good_mu[sample]->SetFillStyle(3000);
        h_ptrel_good_mu[sample]->Sumw2();
        
        name = Form("ptratio_bad_mu_%d", sample);
        h_ptratio_bad_mu[sample] = new TH1D(name, "ptratio_bad_mu "+names[sample]+";ptratio; events", ptratiobins, ptratioMin, ptratioMax);
        h_ptratio_bad_mu[sample]->SetLineColor(colorH[sample]);
        h_ptratio_bad_mu[sample]->SetFillColor(colorH[sample]);
        h_ptratio_bad_mu[sample]->SetMarkerColor(colorH[sample]);
        h_ptratio_bad_mu[sample]->SetFillStyle(3000);
        h_ptratio_bad_mu[sample]->Sumw2();
        
        name = Form("ptrel_bad_mu_%d", sample);
        h_ptrel_bad_mu[sample] = new TH1D(name, "ptrel_bad_mu "+names[sample]+";ptrel[GeV]; events", 50, 0., 50.);
        h_ptrel_bad_mu[sample]->SetLineColor(colorH[sample]);
        h_ptrel_bad_mu[sample]->SetFillColor(colorH[sample]);
        h_ptrel_bad_mu[sample]->SetMarkerColor(colorH[sample]);
        h_ptrel_bad_mu[sample]->SetFillStyle(3000);
        h_ptrel_bad_mu[sample]->Sumw2();
        */


        
        for(int var = 0; var < 6; var++){
            for(int i = 0; i < 15; i++){
                name = Form("histoel_%d_%d_%d", sample, var, i);
                histoele[i][sample][var] = new TH1D(name, ";" + nameLep[0] + xaxisvar[var] + "; events", numberBins[var]-1, Bins[var]);
                histoele[i][sample][var]->SetLineColor(colorH[sample]);
                histoele[i][sample][var]->SetFillColor(colorH[sample]);
                histoele[i][sample][var]->SetMarkerColor(colorH[sample]);
                histoele[i][sample][var]->Sumw2();
            }
        }

        for(int var = 0; var < 6; var++){
            for(int i = 0; i < 8; i++){
                name = Form("histomu_%d_%d_%d", sample, var, i);
                histomu[i][sample][var] = new TH1D(name, ";" + nameLep[1] + xaxisvar[var] + "; events", numberBins[var]-1, Bins[var]);
                histomu[i][sample][var]->SetLineColor(colorH[sample]);
                histomu[i][sample][var]->SetFillColor(colorH[sample]);
                histomu[i][sample][var]->SetMarkerColor(colorH[sample]);
                histomu[i][sample][var]->Sumw2();
            }
        }

        for(int i = 0; i < 15; i++){
                
                name = Form("electron_pteta_%d_%d", sample, i);
                h_ele_pteta[i][sample] = new TH2D(name, "Scale Factors ;Electron p_{T} (GeV); #eta; events", numberBinspt2Del - 1,  ptb2Del, numberBinseta2D - 1, etab6el);
                h_ele_pteta[i][sample]->SetLineColor(colorH[sample]);
                h_ele_pteta[i][sample]->SetFillColor(colorH[sample]);
                h_ele_pteta[i][sample]->SetFillStyle(3000);
                h_ele_pteta[i][sample]->SetMarkerColor(colorH[sample]);
                h_ele_pteta[i][sample]->Sumw2();

                name = Form("electron_ptactiv_%d_%d", sample, i);
                h_ele_ptactiv[i][sample] = new TH2D(name, "Scale Factors ;Electron p_{T} (GeV); Activity; events", numberBinspt2Del - 1,  ptb2Del, numberBins[5]-1, activb);
                h_ele_ptactiv[i][sample]->SetLineColor(colorH[sample]);
                h_ele_ptactiv[i][sample]->SetFillColor(colorH[sample]);
                h_ele_ptactiv[i][sample]->SetFillStyle(3000);
                h_ele_ptactiv[i][sample]->SetMarkerColor(colorH[sample]);
                h_ele_ptactiv[i][sample]->Sumw2();

                name = Form("electron_etaactiv_%d_%d", sample, i);
                h_ele_etaactiv[i][sample] = new TH2D(name, "Scale Factors ;Electron #eta; Activity; events", numberBinseta2D - 1,  etab6el, numberBins[5]-1, activb);
                h_ele_etaactiv[i][sample]->SetLineColor(colorH[sample]);
                h_ele_etaactiv[i][sample]->SetFillColor(colorH[sample]);
                h_ele_etaactiv[i][sample]->SetFillStyle(3000);
                h_ele_etaactiv[i][sample]->SetMarkerColor(colorH[sample]);
                h_ele_etaactiv[i][sample]->Sumw2();

                name = Form("electron_activHT_%d_%d", sample, i);
                h_ele_activHT[i][sample] = new TH2D(name, "Scale Factors ;Activity; H_{T}; events", numberBinsactiv2D - 1,  activb2D, numberBinsHT2D - 1, HTb2D);
                h_ele_activHT[i][sample]->SetLineColor(colorH[sample]);
                h_ele_activHT[i][sample]->SetFillColor(colorH[sample]);
                h_ele_activHT[i][sample]->SetFillStyle(3000);
                h_ele_activHT[i][sample]->SetMarkerColor(colorH[sample]);
                h_ele_activHT[i][sample]->Sumw2();
                
                name = Form("electron_activnpv_%d_%d", sample, i);
                h_ele_activnpv[i][sample] = new TH2D(name, "Scale Factors ;Activity; #PV; events", numberBinsactiv2D - 1,  activb2D, numberBinsNPV2D - 1, npvb6);
                h_ele_activnpv[i][sample]->SetLineColor(colorH[sample]);
                h_ele_activnpv[i][sample]->SetFillColor(colorH[sample]);
                h_ele_activnpv[i][sample]->SetFillStyle(3000);
                h_ele_activnpv[i][sample]->SetMarkerColor(colorH[sample]);
                h_ele_activnpv[i][sample]->Sumw2();

                name = Form("electron_HTnpv_%d_%d", sample, i);
                h_ele_HTnpv[i][sample] = new TH2D(name, "Scale Factors ; H_{T}; #PV; events", numberBinsHT2D - 1,  HTb2D, numberBinsNPV2D - 1, npvb6);
                h_ele_HTnpv[i][sample]->SetLineColor(colorH[sample]);
                h_ele_HTnpv[i][sample]->SetFillColor(colorH[sample]);
                h_ele_HTnpv[i][sample]->SetFillStyle(3000);
                h_ele_HTnpv[i][sample]->SetMarkerColor(colorH[sample]);
                h_ele_HTnpv[i][sample]->Sumw2();


                for(int k = 0; k < 4; k++){
                    name = Form("electron_pteta_2Dnpv%d_%d_%d", sample, i, k);
                    h_ele_pteta_2Dnpv[i][sample][k] = new TH2D(name, "Scale Factors ;Electron p_{T} (GeV); #eta; events", numberBinspt2Del - 1,  ptb2Del, numberBinseta2D - 1, etab6el);
                    h_ele_pteta_2Dnpv[i][sample][k]->SetLineColor(colorH[sample]);
                    h_ele_pteta_2Dnpv[i][sample][k]->SetFillColor(colorH[sample]);
                    h_ele_pteta_2Dnpv[i][sample][k]->SetFillStyle(3000);
                    h_ele_pteta_2Dnpv[i][sample][k]->SetMarkerColor(colorH[sample]);
                    h_ele_pteta_2Dnpv[i][sample][k]->Sumw2();
                }

                name = Form("electron_pt_eta_npv_%d_%d", sample, i);
                h_ele_pt_eta_npv[i][sample] = new TH3D(name, "Scale Factors ;Electron p_{T} (GeV); #eta; # PV", numberBinspt2Del - 1,  ptb2Del, numberBinseta2D - 1, etab6el, numberBinsNPV2D - 1, npvb6);
                h_ele_pt_eta_npv[i][sample]->SetLineColor(colorH[sample]);
                h_ele_pt_eta_npv[i][sample]->SetFillColor(colorH[sample]);
                h_ele_pt_eta_npv[i][sample]->SetFillStyle(3000);
                h_ele_pt_eta_npv[i][sample]->SetMarkerColor(colorH[sample]);
                h_ele_pt_eta_npv[i][sample]->Sumw2();


        }

        for(int i = 0; i < 8; i++){
            
                name = Form("muon_pteta_%d_%d", sample, i);
                h_mu_pteta[i][sample] = new TH2D(name, "Scale Factors ;Muon p_{T} (GeV); #eta; events", numberBinspt2Dmu - 1,  ptb2Dmu,  numberBinseta2Dmu - 1, etab6mu);
                h_mu_pteta[i][sample]->SetLineColor(colorH[sample]);
                h_mu_pteta[i][sample]->SetFillColor(colorH[sample]);
                h_mu_pteta[i][sample]->SetFillStyle(3000);
                h_mu_pteta[i][sample]->SetMarkerColor(colorH[sample]);
                h_mu_pteta[i][sample]->Sumw2();

                name = Form("muon_ptactiv_%d_%d", sample, i);
                h_mu_ptactiv[i][sample] = new TH2D(name, "Scale Factors ;Muon p_{T} (GeV); Activity; events", numberBinspt2Dmu - 1,  ptb2Dmu, numberBins[5]-1, activb);
                h_mu_ptactiv[i][sample]->SetLineColor(colorH[sample]);
                h_mu_ptactiv[i][sample]->SetFillColor(colorH[sample]);
                h_mu_ptactiv[i][sample]->SetFillStyle(3000);
                h_mu_ptactiv[i][sample]->SetMarkerColor(colorH[sample]);
                h_mu_ptactiv[i][sample]->Sumw2();

                name = Form("muon_etaactiv_%d_%d", sample, i);
                h_mu_etaactiv[i][sample] = new TH2D(name, "Scale Factors ;Muon #eta; Activity; events", numberBinseta2Dmu - 1,  etab6mu, numberBins[5]-1, activb);
                h_mu_etaactiv[i][sample]->SetLineColor(colorH[sample]);
                h_mu_etaactiv[i][sample]->SetFillColor(colorH[sample]);
                h_mu_etaactiv[i][sample]->SetFillStyle(3000);
                h_mu_etaactiv[i][sample]->SetMarkerColor(colorH[sample]);
                h_mu_etaactiv[i][sample]->Sumw2();

                name = Form("muon_activHT_%d_%d", sample, i);
                h_mu_activHT[i][sample] = new TH2D(name, "Scale Factors ;Activity; H_{T}; events", numberBinsactiv2D - 1,  activb2D, numberBinsHT2D - 1, HTb2D);
                h_mu_activHT[i][sample]->SetLineColor(colorH[sample]);
                h_mu_activHT[i][sample]->SetFillColor(colorH[sample]);
                h_mu_activHT[i][sample]->SetFillStyle(3000);
                h_mu_activHT[i][sample]->SetMarkerColor(colorH[sample]);
                h_mu_activHT[i][sample]->Sumw2();

                name = Form("muon_activnpv_%d_%d", sample, i);
                h_mu_activnpv[i][sample] = new TH2D(name, "Scale Factors ;Activity; #PV; events", numberBinsactiv2D - 1,  activb2D, numberBinsNPV2D - 1, npvb6);
                h_mu_activnpv[i][sample]->SetLineColor(colorH[sample]);
                h_mu_activnpv[i][sample]->SetFillColor(colorH[sample]);
                h_mu_activnpv[i][sample]->SetFillStyle(3000);
                h_mu_activnpv[i][sample]->SetMarkerColor(colorH[sample]);
                h_mu_activnpv[i][sample]->Sumw2();

                name = Form("muon_HTnpv_%d_%d", sample, i);
                h_mu_HTnpv[i][sample] = new TH2D(name, "Scale Factors ; H_{T}; #PV; events", numberBinsHT2D - 1,  HTb2D, numberBinsNPV2D - 1, npvb6);
                h_mu_HTnpv[i][sample]->SetLineColor(colorH[sample]);
                h_mu_HTnpv[i][sample]->SetFillColor(colorH[sample]);
                h_mu_HTnpv[i][sample]->SetFillStyle(3000);
                h_mu_HTnpv[i][sample]->SetMarkerColor(colorH[sample]);
                h_mu_HTnpv[i][sample]->Sumw2();


                for(int k = 0; k < 4; k++){
                    name = Form("muon_pteta_2Dnpv%d_%d_%d", sample, i, k);
                    h_mu_pteta_2Dnpv[i][sample][k] = new TH2D(name, "Scale Factors ;Muon p_{T} (GeV); #eta; events", numberBinspt2Dmu - 1,  ptb2Dmu, numberBinseta2Dmu - 1, etab6mu);
                    h_mu_pteta_2Dnpv[i][sample][k]->SetLineColor(colorH[sample]);
                    h_mu_pteta_2Dnpv[i][sample][k]->SetFillColor(colorH[sample]);
                    h_mu_pteta_2Dnpv[i][sample][k]->SetFillStyle(3000);
                    h_mu_pteta_2Dnpv[i][sample][k]->SetMarkerColor(colorH[sample]);
                    h_mu_pteta_2Dnpv[i][sample][k]->Sumw2();
                }

                name = Form("muon_pt_eta_npv_%d_%d", sample, i);
                h_mu_pt_eta_npv[i][sample] = new TH3D(name, "Scale Factors ;Muon p_{T} (GeV); #eta; # PV", numberBinspt2Dmu - 1,  ptb2Dmu, numberBinseta2Dmu - 1, etab6mu, numberBinsNPV2D - 1, npvb6);
                h_mu_pt_eta_npv[i][sample]->SetLineColor(colorH[sample]);
                h_mu_pt_eta_npv[i][sample]->SetFillColor(colorH[sample]);
                h_mu_pt_eta_npv[i][sample]->SetFillStyle(3000);
                h_mu_pt_eta_npv[i][sample]->SetMarkerColor(colorH[sample]);
                h_mu_pt_eta_npv[i][sample]->Sumw2();

            
        }

        /*
        name = Form("electron_ptactiv_%d", sample);
        h_ele_ptactiv[sample] = new TH2D(name, "Scale Factors ;Electron p_{T} (GeV); Activity; events", numberBinspt2D - 1,  ptb2D, numberBins[5]-1, activb);
        h_ele_ptactiv[sample]->SetLineColor(colorH[sample]);
        h_ele_ptactiv[sample]->SetFillColor(colorH[sample]);
        h_ele_ptactiv[sample]->SetFillStyle(3000);
        h_ele_ptactiv[sample]->SetMarkerColor(colorH[sample]);
        h_ele_ptactiv[sample]->Sumw2();

        name = Form("electron_ptactiv_cut_%d", sample);
        h_ele_ptactivc[sample] = new TH2D(name, "Scale Factors ;Electron p_{T} (GeV); Activity; events", numberBinspt2D - 1,  ptb2D, numberBins[5]-1, activb);
        h_ele_ptactivc[sample]->SetLineColor(colorH[sample]);
        h_ele_ptactivc[sample]->SetFillColor(colorH[sample]);
        h_ele_ptactivc[sample]->SetFillStyle(3000);
        h_ele_ptactivc[sample]->SetMarkerColor(colorH[sample]);
        h_ele_ptactivc[sample]->Sumw2();

        name = Form("muon_ptactiv_%d", sample);
        h_mu_ptactiv[sample] = new TH2D(name, "Scale Factors ;Muon p_{T} (GeV); Activity; events", numberBinspt2D - 1,  ptb2D, numberBins[5]-1, activb);
        h_mu_ptactiv[sample]->SetLineColor(colorH[sample]);
        h_mu_ptactiv[sample]->SetFillColor(colorH[sample]);
        h_mu_ptactiv[sample]->SetFillStyle(3000);
        h_mu_ptactiv[sample]->SetMarkerColor(colorH[sample]);
        h_mu_ptactiv[sample]->Sumw2();

        name = Form("muon_ptactiv_cut_%d", sample);
        h_mu_ptactivc[sample] = new TH2D(name, "Scale Factors ;Muon p_{T} (GeV); Activity; events", numberBinspt2D - 1,  ptb2D, numberBins[5]-1, activb);
        h_mu_ptactivc[sample]->SetLineColor(colorH[sample]);
        h_mu_ptactivc[sample]->SetFillColor(colorH[sample]);
        h_mu_ptactivc[sample]->SetFillStyle(3000);
        h_mu_ptactivc[sample]->SetMarkerColor(colorH[sample]);
        h_mu_ptactivc[sample]->Sumw2();
        */

        name = Form("h_npv_distr_%d", sample);
        h_npv_distr[sample] = new TH1D(name, "h_npv_distr "+names[sample]+";h_npv_distr; events", 40, 0.5, 40.5);
        h_npv_distr[sample]->SetLineColor(colorH[sample]);
        h_npv_distr[sample]->SetFillColor(colorH[sample]);
        h_npv_distr[sample]->SetMarkerColor(colorH[sample]);
        h_npv_distr[sample]->SetFillStyle(3000);
        h_npv_distr[sample]->Sumw2();

    }


    int total_taus_number;

    int total_muons_number_1;
    int total_muons_number_2;
    int total_electron_number_1;
    int total_electron_number_2;

    int muons_number_HT200_1;
    int muons_number_HT200_2;
    int electron_number_HT200_1;
    int electron_number_HT200_2;

    //0 - vloose, 1 - vloose FOIDEmu, 2 - vlooseFOIDISOEMU, 3 - tight
    double valuesMVA[3][4];
    valuesMVA[0][0] = -0.16;
    valuesMVA[1][0] = -0.65;
    valuesMVA[2][0] = -0.74;

    valuesMVA[0][1] = -0.70;
    valuesMVA[1][1] = -0.83; 
    valuesMVA[2][1] = -0.92; 

    valuesMVA[0][2] = -0.155;
    valuesMVA[1][2] = -0.56;
    valuesMVA[2][2] = -0.76; 

    valuesMVA[0][3] = 0.87;
    valuesMVA[1][3] = 0.60;
    valuesMVA[2][3] = 0.17;

    int activNumber[6][5] = {0};

    TLorentzVector l0p4, l1p4;

    //TFile* file0 = TFile::Open("files/sf_mu_mediumID_multi.root");
    //TFile* file0 = TFile::Open("files/sf_el_tight_IDEmu_ISOEMu_ra5.root");
    //TH3D *hcorr = (TH3D*)file0->Get("histo3D");

    for (int sam = 0; sam < nSamples; ++sam) {
        if(sam == 0) continue;
        
        hfile[sam] = new TFile("/Users/ikhvastu/Desktop/CERN/MCsamples/" + fileList[sam],"read");

        hfile[sam]->cd("FakeElectrons");
        inputTreeFake[sam] = static_cast<TTree*>(hfile[sam]->Get("FakeElectrons/fakeTree"));
        
        inputTreeFake[sam]->SetBranchAddress("_eventNb",   &_eventNb);
        inputTreeFake[sam]->SetBranchAddress("_runNb",     &_runNb);
        inputTreeFake[sam]->SetBranchAddress("_lumiBlock", &_lumiBlock);
        
        inputTreeFake[sam]->SetBranchAddress("_nLeptons", &_nLeptons);
        inputTreeFake[sam]->SetBranchAddress("_nEle", &_nEle);
        inputTreeFake[sam]->SetBranchAddress("_nMu", &_nMu);
        inputTreeFake[sam]->SetBranchAddress("_nTau", &_nTau);

        inputTreeFake[sam]->SetBranchAddress("_lPt", &_lPt);
        inputTreeFake[sam]->SetBranchAddress("_lPhi", &_lPhi);
        inputTreeFake[sam]->SetBranchAddress("_lEta", &_lEta);
        inputTreeFake[sam]->SetBranchAddress("_lE", &_lE);


        inputTreeFake[sam]->SetBranchAddress("_nuPtmc", &_nuPtmc);
        inputTreeFake[sam]->SetBranchAddress("_lPtmc", &_lPtmc);
        inputTreeFake[sam]->SetBranchAddress("_lEtamc", &_lEtamc);
        inputTreeFake[sam]->SetBranchAddress("_lPhimc", &_lPhimc);
        inputTreeFake[sam]->SetBranchAddress("_lEmc", &_lEmc);
        inputTreeFake[sam]->SetBranchAddress("_lpdgmc", &_lpdgmc);

        inputTreeFake[sam]->SetBranchAddress("_mtmc", &_mtmc);

        inputTreeFake[sam]->SetBranchAddress("_ossf", &_ossf);
        inputTreeFake[sam]->SetBranchAddress("_mll", &_mll);
        
        inputTreeFake[sam]->SetBranchAddress("_indeces", &_indeces);
        inputTreeFake[sam]->SetBranchAddress("_flavors", &_flavors); // flavor: 0 - e, 1 - mu, 2 - tau
        inputTreeFake[sam]->SetBranchAddress("_charges", &_charges); // charge (+/-)
        inputTreeFake[sam]->SetBranchAddress("_isolation", &_isolation); // lepton relative isolation
        inputTreeFake[sam]->SetBranchAddress("_miniisolation", &_miniisolation); // lepton miniisolation
        inputTreeFake[sam]->SetBranchAddress("_multiisolation", &_multiisolation); // lepton miniisolation
       
        inputTreeFake[sam]->SetBranchAddress("_origin", &_origin); //lepton detailed origin
        inputTreeFake[sam]->SetBranchAddress("_originReduced", &_originReduced); // lepton origin: 0-prompt, 1-b-quark, 2-c, 3-uds, 4-unmatched
        
        inputTreeFake[sam]->SetBranchAddress("_PVchi2", &_PVchi2); //primary vertex chi2
        inputTreeFake[sam]->SetBranchAddress("_PVerr", &_PVerr); //primary vertex uncertainty
        
        inputTreeFake[sam]->SetBranchAddress("_ipPV", &_ipPV); //lepton impact parameter to collision vertext (primary vertex) in the transverse plane
        inputTreeFake[sam]->SetBranchAddress("_ipPVerr", &_ipPVerr); //error on ipPV
        inputTreeFake[sam]->SetBranchAddress("_ipZPV", &_ipZPV); //lepton impact parameter in z-direction
        inputTreeFake[sam]->SetBranchAddress("_ipZPVerr", &_ipZPVerr);//error on ipZ
        
        inputTreeFake[sam]->SetBranchAddress("_3dIP", &_3dIP); // lepton 3D impact parameter
        inputTreeFake[sam]->SetBranchAddress("_3dIPerr", &_3dIPerr); //error on 3dIP
        inputTreeFake[sam]->SetBranchAddress("_3dIPsig", &_3dIPsig); //significance of 3dIP
        
        
        inputTreeFake[sam]->SetBranchAddress("_mt", &_mt); // transverse mass of lepton and MET (missing transverse energy = neutrino)
        inputTreeFake[sam]->SetBranchAddress("_isloose", &_isloose); //if lepton satisfies loose selection
        inputTreeFake[sam]->SetBranchAddress("_ismedium", &_ismedium); //if lepton satisfies loose selection
        inputTreeFake[sam]->SetBranchAddress("_istight", &_istight); //if lepton satisfies tight selection
        //inputTreeFake[sam]->SetBranchAddress("_istightID", &_istightID); //if lepton satisfies tight selection
        inputTreeFake[sam]->SetBranchAddress("_isvetoIDCutBased", &_isvetoIDCutBased);
        inputTreeFake[sam]->SetBranchAddress("_islooseIDCutBased", &_islooseIDCutBased); //if lepton satisfies tight selection
        inputTreeFake[sam]->SetBranchAddress("_ismediumIDCutBased", &_ismediumIDCutBased); //if lepton satisfies tight selection
        inputTreeFake[sam]->SetBranchAddress("_istightIDCutBased", &_istightIDCutBased); //if lepton satisfies tight selection
        
        inputTreeFake[sam]->SetBranchAddress("_closeJetPtAll", &_closeJetPtAll); //pt of the jet nesting lepton
        inputTreeFake[sam]->SetBranchAddress("_closeJetAngAll", &_closeJetAngAll); //delta R (=sqrt[(delta phi)^2 + (delta eta)^2]) of lepton to jet
        inputTreeFake[sam]->SetBranchAddress("_ptRelAll", &_ptRelAll); //pt of lepton relative to jet = pt of lepton in jet rest frame
        inputTreeFake[sam]->SetBranchAddress("_ptrel", &_ptrel); //pt of lepton relative to jet = pt of lepton in jet rest frame
        
        inputTreeFake[sam]->SetBranchAddress("_closeJetPtAllMC", &_closeJetPtAllMC); //MC-truth pt og the quark matched to lepton
        inputTreeFake[sam]->SetBranchAddress("_closeJetPtAllstatus", &_closeJetPtAllstatus); // MC status of this quark: 3 = hard scattering, 2 = transitional quark
        inputTreeFake[sam]->SetBranchAddress("_partonIdMatched", &_partonIdMatched); //PDG ID of this quark
        
        inputTreeFake[sam]->SetBranchAddress("_n_PV", &_n_PV); // number of primary vertices in the collision
        
        inputTreeFake[sam]->SetBranchAddress("_met", &_met); // missing transverse energy
        inputTreeFake[sam]->SetBranchAddress("_met_phi", &_met_phi); //phi of MET
        inputTreeFake[sam]->SetBranchAddress("HT", &HT); //hadronic activity in the event: sum of jet pt for jets with pt > 40 GeV
        
        inputTreeFake[sam]->SetBranchAddress("_RhoCN", &RhoCN); // missing transverse energy
        inputTreeFake[sam]->SetBranchAddress("_RhoAll", &RhoAll);

        inputTreeFake[sam]->SetBranchAddress("_mompt", &_mompt); //kinematical information of the mother particle of the lepton
        inputTreeFake[sam]->SetBranchAddress("_momphi", &_momphi);
        inputTreeFake[sam]->SetBranchAddress("_mometa", &_mometa);
        inputTreeFake[sam]->SetBranchAddress("_mompdg", &_mompdg);
        
        inputTreeFake[sam]->SetBranchAddress("_n_bJets", &_n_bJets); //number of b-jets with pt > 25 GeV in the event
        inputTreeFake[sam]->SetBranchAddress("_n_Jets", &_n_Jets); //number of jets with pt > 25 GeV in the event
        inputTreeFake[sam]->SetBranchAddress("_bTagged", &_bTagged); //b-tagging info on each jet
        inputTreeFake[sam]->SetBranchAddress("_jetEta", &_jetEta);
        inputTreeFake[sam]->SetBranchAddress("_jetPhi", &_jetPhi);
        inputTreeFake[sam]->SetBranchAddress("_jetPt", &_jetPt);
        inputTreeFake[sam]->SetBranchAddress("_csv", &_csv); //b-tagging discriminator of each jet

        inputTreeFake[sam]->SetBranchAddress("_mvaValue", &_mvaValue);
        inputTreeFake[sam]->SetBranchAddress("_trigEmulator", &_trigEmulator);
        inputTreeFake[sam]->SetBranchAddress("_isotrigEmulator", &_isotrigEmulator);
        inputTreeFake[sam]->SetBranchAddress("_chargeConst", &_chargeConst);

        inputTreeFake[sam]->SetBranchAddress("_hitsNumber", &_hitsNumber);

        inputTreeFake[sam]->SetBranchAddress("_activity", &_activity);

        inputTreeFake[sam]->SetBranchAddress("_vtxFitConversion", &_vtxFitConversion);

        inputTreeFake[sam]->SetBranchAddress("_lgenParticle", &_lgenParticle);

        //if(sam == 1){

            inputTreeFake[sam]->SetBranchAddress("_isPromptFinalState", &_isPromptFinalState);
            inputTreeFake[sam]->SetBranchAddress("_fromHardProcessFinalState", &_fromHardProcessFinalState);

        //}

        int counterEle = 0;

        long nEntries = inputTreeFake[sam]->GetEntries();
        std::cout<<"Sample "<<fileList[sam]<<std::endl;
        std::cout<<"Entries "<<nEntries<<std::endl;
        
        for (long it=0; it!=nEntries; ++it) {
            //if(it > 100000) break;
            inputTreeFake[sam]->GetEntry(it);

            if (it%100000 == 0)
                cout<<'.'<<flush;

            h_npv_distr[sam]->Fill(_n_PV);
            /*
            h_RhoAll[sam]->Fill(RhoAll);
            h_RhoNC[sam]->Fill(RhoCN);
            */
            //if(_eventNb != 884281) continue;
            

            int posCharge = 0;
            int negCharge = 0;
            int nLoc = 0;
            int leptInd[10];

            int index1 = -1;
            int index2 = -1;
            
            
            for(int num = 0; num < _nEle + _nMu; num++){
                //std::cout << "PT: " << _lPt[num] << " " << _flavors[num] << " " << _charges[num] << " " << _activity[num] << std::endl;
                //std::cout << "ORIGIN: " << _origin[num] << " " << _isPromptFinalState[num] << " " << _fromHardProcessFinalState[num] << std::endl;
                if (_lPt[num] < 10.) continue;
                
                //if(!(_lgenParticle[num] && _isPromptFinalState[num])) continue;
                
                if (!(_origin[num] == 0 || _isPromptFinalState[num] || _fromHardProcessFinalState[num])) continue;
                if (_flavors[num] == 0 && fabs(_lpdgmc[num]) != 11) continue;
                if (_flavors[num] == 1 && fabs(_lpdgmc[num]) != 13) continue;
                
                l0p4.SetPtEtaPhiE(_lPt[num],_lEta[num],_lPhi[num],_lE[num]);
                l1p4.SetPtEtaPhiE(_lPtmc[num],_lEtamc[num],_lPhimc[num],_lEmc[num]);
                
                double deltaR_2l = l0p4.DeltaR(l1p4);
                
                //if(deltaR_2l > 0.2) continue;
                
                    
                leptInd[nLoc] = num;
                nLoc++;
            }

            //if(nLoc != 2) continue;

            //if(sam == 1){
            /*
            if(nLoc > 2) 
                std::cout << "more than 2 prompt leptons" << std::endl;                
                */
                
                l0p4.SetPtEtaPhiE(_lPt[leptInd[0]],_lEta[leptInd[0]],_lPhi[leptInd[0]],_lE[leptInd[0]]);
                l1p4.SetPtEtaPhiE(_lPt[leptInd[1]],_lEta[leptInd[1]],_lPhi[leptInd[1]],_lE[leptInd[1]]);

                double deltaR_2l = l0p4.DeltaR(l1p4);
                
                
                //if(deltaR_2l < 0.2) continue;
                

                double deltaMZ = 9999.;
                
                for (int l0 = 0; l0<nLoc; ++l0) {
                    l0p4.SetPtEtaPhiE(_lPt[leptInd[l0]],_lEta[leptInd[l0]],_lPhi[leptInd[l0]],_lE[leptInd[l0]]);
                    for (int l1 = l0+1; l1<nLoc; ++l1) {
                        if (_charges[leptInd[l0]] != _charges[leptInd[l1]]) {
                            l1p4.SetPtEtaPhiE(_lPt[leptInd[l1]],_lEta[leptInd[l1]],_lPhi[leptInd[l1]],_lE[leptInd[l1]]);
                            l1p4+=l0p4;
                            double mdiL = l1p4.M();
                            if (_flavors[leptInd[l0]] == _flavors[leptInd[l1]] ) {
                                //if (mdiL < minll) minll = mdiL;
                                if (fabs(mdiL - 91) < deltaMZ) {
                                    deltaMZ = fabs(mdiL - 91);
                                    index1 = leptInd[l0];
                                    index2 = leptInd[l1];
                                }
                            }
                        }
                    }
                }
                

                //if(deltaMZ == 9999.) continue;
                //if(_flavors[index1] != _flavors[index2]) continue;
                //if(_charges[index1] == _charges[index2]) continue;


                //std::cout << nLoc << std::endl;
                //std::cout << deltaMZ << std::endl;
                //std::cout << index1 << " " << index2 << std::endl;
                /*
                if(deltaR_2l < 0.2)
                    std::cout << deltaR_2l << " " << deltaMZ << std::endl;
                */
                //if(nLoc != 2) continue;
                //if(deltaMZ > 15) continue;
                
                //if(nLoc > 2)
                //    std::cout << "Bulsit" << std::endl;
                
                
            //}
                
            if(_flavors[leptInd[0]] == 0 && _flavors[leptInd[1]] == 0)
                h_deltaR[0][sam]->Fill(deltaR_2l);
            if(_flavors[leptInd[0]] == 1 && _flavors[leptInd[1]] == 1)
                h_deltaR[1][sam]->Fill(deltaR_2l);
            if((_flavors[leptInd[0]] == 0 && _flavors[leptInd[1]] == 1) || (_flavors[leptInd[0]] == 1 && _flavors[leptInd[1]] == 0))
                h_deltaR[2][sam]->Fill(deltaR_2l);

            if(_flavors[leptInd[0]] != _flavors[leptInd[1]]) continue;

            for(int i = 0; i < nLoc; i++){
                if(deltaR_2l < 0.2){
                    if(_flavors[leptInd[i]] == 0)
                        h_deltaPT[0][sam]->Fill(TMath::Abs(_lPt[leptInd[i]] - _lPtmc[leptInd[i]]) / _lPtmc[leptInd[i]]);
                    if(_flavors[leptInd[i]] == 1)
                        h_deltaPT[1][sam]->Fill(TMath::Abs(_lPt[leptInd[i]] - _lPtmc[leptInd[i]]) / _lPtmc[leptInd[i]]);
                }

                if(deltaR_2l > 0.2){
                    if(_flavors[leptInd[i]] == 0)
                        h_deltaPT[2][sam]->Fill(TMath::Abs(_lPt[leptInd[i]] - _lPtmc[leptInd[i]]) / _lPtmc[leptInd[i]]);
                    if(_flavors[leptInd[i]] == 1)
                        h_deltaPT[3][sam]->Fill(TMath::Abs(_lPt[leptInd[i]] - _lPtmc[leptInd[i]]) / _lPtmc[leptInd[i]]);
                }
            }
            

            int goodLep[2] = {index1, index2};


            /*
            int goodLep[10];
            for(int i = 0; i < nLoc; i++)
                goodLep[i] = leptInd[i];
            */
            
            //for(int num = 0; num < _nEle + _nMu; num++){
            /*
            int goodLeptons = 0;
            for(int i = 0; i < nLoc; i++){
                //int num = goodLep[i];
                int num = leptInd[i];
                if (_lPt[num] < 10.) continue;
                
                if (_flavors[num] == 0){
                    //if(fabs(_lpdgmc[num]) != 11) continue;
                    goodLeptons++;
            
                        FillHisto(h_ele_pteta_2Dnpv[0][sam], h_ele_pt_eta_npv[0][sam], h_ele_pteta[0][sam], h_ele_ptactiv[0][sam], histoele[0][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);

                        
                        if(_isvetoIDCutBased[num] && _miniisolation[num][0] < 0.1)
                            FillHisto(h_ele_pteta_2Dnpv[1][sam], h_ele_pt_eta_npv[1][sam], h_ele_pteta[1][sam], h_ele_ptactiv[1][sam], histoele[1][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(_islooseIDCutBased[num] && _miniisolation[num][0] < 0.1)
                            FillHisto(h_ele_pteta_2Dnpv[2][sam], h_ele_pt_eta_npv[2][sam], h_ele_pteta[2][sam], h_ele_ptactiv[2][sam], histoele[2][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(_islooseIDCutBased[num] && _miniisolation[num][0] < 0.4)
                            FillHisto(h_ele_pteta_2Dnpv[3][sam], h_ele_pt_eta_npv[3][sam], h_ele_pteta[3][sam], h_ele_ptactiv[3][sam], histoele[3][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(_ismediumIDCutBased[num] && _miniisolation[num][0] < 0.1)
                            FillHisto(h_ele_pteta_2Dnpv[4][sam], h_ele_pt_eta_npv[4][sam], h_ele_pteta[4][sam], h_ele_ptactiv[4][sam], histoele[4][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(_istightIDCutBased[num] && _miniisolation[num][0] < 0.1){
                            FillHisto(h_ele_pteta_2Dnpv[5][sam], h_ele_pt_eta_npv[5][sam], h_ele_pteta[5][sam], h_ele_ptactiv[5][sam], histoele[5][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            
                           

                        }
                        bool passedMVA = false;
                        if (TMath::Abs(_lEta[num]) < 0.8 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[0][0];
                        } else if (TMath::Abs(_lEta[num]) < 1.479 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[1][0];
                        } else {
                            passedMVA = _mvaValue[num]> valuesMVA[2][0];
                        }

                        if(passedMVA && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _miniisolation[num][0] < 0.1)
                            FillHisto(h_ele_pteta_2Dnpv[6][sam], h_ele_pt_eta_npv[6][sam], h_ele_pteta[6][sam], h_ele_ptactiv[6][sam], histoele[6][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            
                        if (TMath::Abs(_lEta[num]) < 0.8 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[0][3];
                        } else if (TMath::Abs(_lEta[num]) < 1.479 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[1][3];
                        } else {
                            passedMVA = _mvaValue[num]> valuesMVA[2][3];
                        }
                        
                        if(passedMVA && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _3dIPsig[num] < 4 && _miniisolation[num][0] < 0.1 && _hitsNumber[num] == 0 && !_vtxFitConversion[num]){
                            FillHisto(h_ele_pteta_2Dnpv[7][sam], h_ele_pt_eta_npv[7][sam], h_ele_pteta[7][sam], h_ele_ptactiv[7][sam], histoele[7][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            
                            
                        }
                        if(passedMVA && _trigEmulator[num] && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _3dIPsig[num] < 4 && ((_miniisolation[num][0] < 0.12) && ((_lPt[num]/_closeJetPtAll[num] > 0.80) || (_ptrel[num] > 7.2))) && _hitsNumber[num] == 0 && !_vtxFitConversion[num] && _chargeConst[num])
                            FillHisto(h_ele_pteta_2Dnpv[8][sam], h_ele_pt_eta_npv[8][sam], h_ele_pteta[8][sam], h_ele_ptactiv[8][sam], histoele[8][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(passedMVA && _trigEmulator[num] && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _3dIPsig[num] < 4 && ((_miniisolation[num][0] < 0.12) && ((_lPt[num]/_closeJetPtAll[num] > 0.80) || (_ptrel[num] > 7.2))) && _isotrigEmulator[num] && _hitsNumber[num] == 0 && !_vtxFitConversion[num] && _chargeConst[num]){
                            FillHisto(h_ele_pteta_2Dnpv[9][sam], h_ele_pt_eta_npv[9][sam], h_ele_pteta[9][sam], h_ele_ptactiv[9][sam], histoele[9][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);

                        }
                        
    

                        if (TMath::Abs(_lEta[num]) < 0.8 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[0][1];
                        } else if (TMath::Abs(_lEta[num]) < 1.479 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[1][1];
                        } else {
                            passedMVA = _mvaValue[num]> valuesMVA[2][1];
                        }

                        if(passedMVA && _trigEmulator[num] && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _miniisolation[num][0] < 0.4 && _hitsNumber[num] <= 1 && !_vtxFitConversion[num])
                            FillHisto(h_ele_pteta_2Dnpv[10][sam], h_ele_pt_eta_npv[10][sam], h_ele_pteta[10][sam], h_ele_ptactiv[10][sam], histoele[10][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
 
                    //}
                }
                else if(_flavors[num] == 1){
                    //if(fabs(_lpdgmc[num]) != 13) continue;
                    goodLeptons++;
                    //if (_originReduced[num] == 0){
                    //if(_origin[num] == 0){


                        FillHisto(h_mu_pteta_2Dnpv[0][sam], h_mu_pt_eta_npv[0][sam], h_mu_pteta[0][sam], h_mu_ptactiv[0][sam], histomu[0][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        
                        if(_isloose[num] && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _3dIPsig[num] < 4 && _miniisolation[num][0] < 0.2)
                            FillHisto(h_mu_pteta_2Dnpv[1][sam], h_mu_pt_eta_npv[1][sam], h_mu_pteta[1][sam], h_mu_ptactiv[1][sam], histomu[1][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(_ismedium[num] && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _3dIPsig[num] < 4 && _miniisolation[num][0] < 0.2){
                            FillHisto(h_mu_pteta_2Dnpv[2][sam], h_mu_pt_eta_npv[2][sam], h_mu_pteta[2][sam], h_mu_ptactiv[2][sam], histomu[2][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            
                        }
                        if(_ismedium[num] && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _3dIPsig[num] < 4 && (_miniisolation[num][0] < 0.16) && ((_lPt[num]/_closeJetPtAll[num] > 0.76) || (_ptrel[num] > 7.2))){
                            
                            FillHisto(h_mu_pteta_2Dnpv[3][sam], h_mu_pt_eta_npv[3][sam], h_mu_pteta[3][sam], h_mu_ptactiv[3][sam], histomu[3][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        }
                        
                        if(_isloose[num] && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _3dIPsig[num] < 4 && _miniisolation[num][0] < 0.4)
                            FillHisto(h_mu_pteta_2Dnpv[4][sam], h_mu_pt_eta_npv[4][sam], h_mu_pteta[4][sam], h_mu_ptactiv[4][sam], histomu[4][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        
                    //}
                }
            }
            
            */
            
            
            //for(int num = 0; num < _nEle + _nMu; num++){

            //nLoc = 2;
            
            //std::cout << "New event: " << _lumiBlock << " " << _eventNb << std::endl;
            for(int i = 0; i < nLoc; i++){
              
                int num = leptInd[i];


                
                //_activity[num] = _isolation[num] - _miniisolation[num][0];
                
                //std::cout << _flavors[i] << " " << _lPt[i] << " " << _lPtmc[i] << " " << deltaR_2l << " " << _isPromptFinalState[i] << " " << _fromHardProcessFinalState[i] << std::endl;
                
                //if(deltaR_2l > 0.2) continue;
                /*
                if(deltaR_2l > 0.2 && TMath::Abs(_lPt[num] - _lPtmc[num]) / _lPtmc[num] > 0.2 && deltaR_2l_mc > 0.2){
                    std::cout << _flavors[i] << " " << _lPt[i] << " " << _lPtmc[i] << " " << _activity[i] << " " << deltaR_2l_mc << std::endl;
                    return;
                }
                */

                //if(TMath::Abs(_lPt[num] - _lPtmc[num]) / _lPtmc[num] > 0.1) continue;

                

                if (_flavors[num] == 0){
                    
                    if(_lPt[num] < 20){
                        std::cout << std::setw(5) << counterEle + 1 << std::setw(10) << _eventNb << std::setw(15) << _lPt[num] << std::setw(15) << _lEta[num] << std::setw(15) << _lPhi[num] << std::setw(15) << _activity[num] << std::setw(15) << _miniisolation[num][0] << std::setw(15) << _isolation[num] << std::endl;
                        counterEle++;
                    }

                    if(counterEle == 100)
                        return;
                        //double weight = 1.;

                        //if(sam == 1)
                        //    weight = weightCalc(_n_PV);
                        //if(_lPt[num] > 10 && _lPt[num] < 20)
                        if(deltaR_2l < 0.2)
                            h_activ_distr[0][0][sam]->Fill(TMath::Min(_activity[num], 0.499));
                        if(deltaR_2l > 0.2)
                            h_activ_distr[1][0][sam]->Fill(TMath::Min(_activity[num], 0.499));
                        
                        FillHisto(h_ele_pteta_2Dnpv[0][sam], h_ele_pt_eta_npv[0][sam], h_ele_pteta[0][sam], h_ele_ptactiv[0][sam], h_ele_etaactiv[0][sam], histoele[0][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        
                        //if(_lPt[num] > 10 && _lPt[num] < 20 && _activity[num] > 1. && _activity[num] < 10.)
                            
                        
                        if(_istightIDCutBased[num])
                            FillHisto(h_ele_pteta_2Dnpv[1][sam], h_ele_pt_eta_npv[1][sam], h_ele_pteta[1][sam], h_ele_ptactiv[1][sam], h_ele_etaactiv[1][sam], histoele[1][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(_ismediumIDCutBased[num])
                            FillHisto(h_ele_pteta_2Dnpv[2][sam], h_ele_pt_eta_npv[2][sam], h_ele_pteta[2][sam], h_ele_ptactiv[2][sam], h_ele_etaactiv[2][sam], histoele[2][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(_islooseIDCutBased[num])
                            FillHisto(h_ele_pteta_2Dnpv[3][sam], h_ele_pt_eta_npv[3][sam], h_ele_pteta[3][sam], h_ele_ptactiv[3][sam], h_ele_etaactiv[3][sam], histoele[3][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        if(_isvetoIDCutBased[num]){
                          
                            FillHisto(h_ele_pteta_2Dnpv[4][sam], h_ele_pt_eta_npv[4][sam], h_ele_pteta[4][sam], h_ele_ptactiv[4][sam], h_ele_etaactiv[4][sam], histoele[4][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        }
                        
                        bool passedMVA = false;
                        if (TMath::Abs(_lEta[num]) < 0.8 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[0][0];
                        } else if (TMath::Abs(_lEta[num]) < 1.479 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[1][0];
                        } else {
                            passedMVA = _mvaValue[num]> valuesMVA[2][0];
                        }

                        if(passedMVA && _ipPV[num] < 0.05 && _ipZPV[num] < 0.1){
                            FillHisto(h_ele_pteta_2Dnpv[5][sam], h_ele_pt_eta_npv[5][sam], h_ele_pteta[5][sam], h_ele_ptactiv[5][sam], h_ele_etaactiv[5][sam], histoele[5][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            if(_miniisolation[num][0] < 0.1)
                                FillHisto(h_ele_pteta_2Dnpv[6][sam], h_ele_pt_eta_npv[6][sam], h_ele_pteta[6][sam], h_ele_ptactiv[6][sam], h_ele_etaactiv[6][sam], histoele[6][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            if(_miniisolation[num][0] < 0.4)
                                FillHisto(h_ele_pteta_2Dnpv[7][sam], h_ele_pt_eta_npv[7][sam], h_ele_pteta[7][sam], h_ele_ptactiv[7][sam], h_ele_etaactiv[7][sam], histoele[7][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            if(_hitsNumber[num] <= 1 && !_vtxFitConversion[num])
                                FillHisto(h_ele_pteta_2Dnpv[8][sam], h_ele_pt_eta_npv[8][sam], h_ele_pteta[8][sam], h_ele_ptactiv[8][sam], h_ele_etaactiv[8][sam], histoele[8][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        }
                        

                        passedMVA = false;
                        if (TMath::Abs(_lEta[num]) < 0.8 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[0][3];
                        } else if (TMath::Abs(_lEta[num]) < 1.479 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[1][3];
                        } else {
                            passedMVA = _mvaValue[num]> valuesMVA[2][3];
                        }
                        
                        if(passedMVA && TMath::Abs(_ipPV[num]) < 0.05 && TMath::Abs(_ipZPV[num]) < 0.1 && _3dIPsig[num] < 4){
                            FillHisto(h_ele_pteta_2Dnpv[9][sam], h_ele_pt_eta_npv[9][sam], h_ele_pteta[9][sam], h_ele_ptactiv[9][sam], h_ele_etaactiv[9][sam], histoele[9][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            
                            if(_trigEmulator[num]){
                                FillHisto(h_ele_pteta_2Dnpv[10][sam], h_ele_pt_eta_npv[10][sam], h_ele_pteta[10][sam], h_ele_ptactiv[10][sam], h_ele_etaactiv[10][sam], histoele[10][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                                if((_miniisolation[num][0] < 0.12) && ((_lPt[num]/_closeJetPtAll[num] > 0.80) || (_ptrel[num] > 7.2))){
                                    FillHisto(h_ele_pteta_2Dnpv[11][sam], h_ele_pt_eta_npv[11][sam], h_ele_pteta[11][sam], h_ele_ptactiv[11][sam], h_ele_etaactiv[11][sam], histoele[11][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                                    if(_isotrigEmulator[num])
                                        FillHisto(h_ele_pteta_2Dnpv[12][sam], h_ele_pt_eta_npv[12][sam], h_ele_pteta[12][sam], h_ele_ptactiv[12][sam], h_ele_etaactiv[12][sam], histoele[12][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                                }
                                
                            }
                            
                            if(_hitsNumber[num] == 0 && !_vtxFitConversion[num] && _chargeConst[num])
                                FillHisto(h_ele_pteta_2Dnpv[13][sam], h_ele_pt_eta_npv[13][sam], h_ele_pteta[13][sam], h_ele_ptactiv[13][sam], h_ele_etaactiv[13][sam], histoele[13][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            
                        }

                        passedMVA = false;
                        if (TMath::Abs(_lEta[num]) < 0.8 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[0][1];
                        } else if (TMath::Abs(_lEta[num]) < 1.479 ) {
                            passedMVA = _mvaValue[num]> valuesMVA[1][1];
                        } else {
                            passedMVA = _mvaValue[num]> valuesMVA[2][1];
                        }

                        if(passedMVA && _ipPV[num] < 0.05 && _ipZPV[num] < 0.1)
                            FillHisto(h_ele_pteta_2Dnpv[14][sam], h_ele_pt_eta_npv[14][sam], h_ele_pteta[14][sam], h_ele_ptactiv[14][sam], h_ele_etaactiv[14][sam], histoele[14][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        
                    //}
                }
                
                else if(_flavors[num] == 1){
                    /*
                    if(TMath::Abs(_lPt[num] - _lPtmc[num]) / _lPtmc[num] > 0.1)
                        std::cout << "Problem here: " << TMath::Abs(_lPt[num] - _lPtmc[num]) / _lPtmc[num] << std::endl;
                        */
                    //if(TMath::Abs(_lPt[num] - _lPtmc[num]) / _lPtmc[num] > 0.2) continue;
                    //h_test[2][sam]->Fill(_lEta[num]);
                    //if(fabs(_lpdgmc[num]) != 13) continue;
                    //h_test[3][sam]->Fill(_lEta[num]);
                    //goodLeptons++;
                   
                    //if (_originReduced[num] == 0){
                    //if (_origin[num] == 0 || _isPromptFinalState[num] || _fromHardProcessFinalState[num]){
                        if(deltaR_2l < 0.2)
                            h_activ_distr[0][1][sam]->Fill(TMath::Min(_activity[num], 0.499));
                        if(deltaR_2l > 0.2)
                            h_activ_distr[1][1][sam]->Fill(TMath::Min(_activity[num], 0.499));
                        
                        FillHisto(h_mu_pteta_2Dnpv[0][sam], h_mu_pt_eta_npv[0][sam], h_mu_pteta[0][sam], h_mu_ptactiv[0][sam], h_mu_etaactiv[0][sam], histomu[0][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        
                        if(_isloose[num]){
                            FillHisto(h_mu_pteta_2Dnpv[1][sam], h_mu_pt_eta_npv[1][sam], h_mu_pteta[1][sam], h_mu_ptactiv[1][sam], h_mu_etaactiv[1][sam], histomu[1][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);

                            if(_ipPV[num] < 0.05 && _ipZPV[num] < 0.1)
                                FillHisto(h_mu_pteta_2Dnpv[2][sam], h_mu_pt_eta_npv[2][sam], h_mu_pteta[2][sam], h_mu_ptactiv[2][sam], h_mu_etaactiv[2][sam], histomu[2][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            if(_3dIPsig[num] < 4)
                                FillHisto(h_mu_pteta_2Dnpv[3][sam], h_mu_pt_eta_npv[3][sam], h_mu_pteta[3][sam], h_mu_ptactiv[3][sam], h_mu_etaactiv[3][sam], histomu[3][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            if(_miniisolation[num][0] < 0.2)
                                FillHisto(h_mu_pteta_2Dnpv[4][sam], h_mu_pt_eta_npv[4][sam], h_mu_pteta[4][sam], h_mu_ptactiv[4][sam], h_mu_etaactiv[4][sam], histomu[4][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            if(_miniisolation[num][0] < 0.4)
                                FillHisto(h_mu_pteta_2Dnpv[5][sam], h_mu_pt_eta_npv[5][sam], h_mu_pteta[5][sam], h_mu_ptactiv[5][sam], h_mu_etaactiv[5][sam], histomu[5][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        }
                        if(_ismedium[num]){
                            FillHisto(h_mu_pteta_2Dnpv[6][sam], h_mu_pt_eta_npv[6][sam], h_mu_pteta[6][sam], h_mu_ptactiv[6][sam], h_mu_etaactiv[6][sam], histomu[6][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                            if((_miniisolation[num][0] < 0.16) && ((_lPt[num]/_closeJetPtAll[num] > 0.76) || (_ptrel[num] > 7.2)))
                                FillHisto(h_mu_pteta_2Dnpv[7][sam], h_mu_pt_eta_npv[7][sam], h_mu_pteta[7][sam], h_mu_ptactiv[7][sam], h_mu_etaactiv[7][sam], histomu[7][sam], _lPt[num], _lEta[num], HT, _n_Jets, _n_PV, _activity[num]);
                        }
                    //}
                }
                
            }
            
            /*
            if(goodLeptons > 2)
                std::cout << "Very Bad: " << goodLeptons << std::endl;
            */
            
            
            
        }
        std::cout << std::endl;

        //g_test[sam] = new TGraphAsymmErrors(h_test[1][sam], h_test[0][sam]);

        // Get 1D
        
        for(int var = 0; var < 6; var++){
            /*
            //electrons
            for(int i = 0; i < 10; i++)
                graph[i][var][sam] = new TGraphAsymmErrors(histoele[i+1][sam][var], histoele[0][sam][var]);
            
            
            //muons
            for(int i = 0; i < 4; i++)
                graph[i+10][var][sam] = new TGraphAsymmErrors(histomu[i+1][sam][var], histomu[0][sam][var]);
            
            */

            //______________________________________
            
            //electrons
            
            graph[0][var][sam] = new TGraphAsymmErrors(histoele[1][sam][var], histoele[0][sam][var]);
            graph[1][var][sam] = new TGraphAsymmErrors(histoele[2][sam][var], histoele[0][sam][var]);
            graph[2][var][sam] = new TGraphAsymmErrors(histoele[3][sam][var], histoele[0][sam][var]);
            graph[3][var][sam] = new TGraphAsymmErrors(histoele[4][sam][var], histoele[0][sam][var]);
            graph[13][var][sam] = new TGraphAsymmErrors(histoele[14][sam][var], histoele[0][sam][var]);
            
            graph[5][var][sam] = new TGraphAsymmErrors(histoele[6][sam][var], histoele[5][sam][var]);
            graph[6][var][sam] = new TGraphAsymmErrors(histoele[7][sam][var], histoele[5][sam][var]);
            graph[7][var][sam] = new TGraphAsymmErrors(histoele[8][sam][var], histoele[5][sam][var]);
            
            graph[10][var][sam] = new TGraphAsymmErrors(histoele[11][sam][var], histoele[10][sam][var]);
            graph[11][var][sam] = new TGraphAsymmErrors(histoele[12][sam][var], histoele[10][sam][var]);
            graph[12][var][sam] = new TGraphAsymmErrors(histoele[13][sam][var], histoele[9][sam][var]);
            

            graph[8][var][sam] = new TGraphAsymmErrors(histoele[9][sam][var], histoele[0][sam][var]);
            
            graph[4][var][sam] = new TGraphAsymmErrors(histoele[5][sam][var], histoele[0][sam][var]);    
            graph[9][var][sam] = new TGraphAsymmErrors(histoele[10][sam][var], histoele[0][sam][var]);
            
            

            //muons
            
           
            graph[17][var][sam] = new TGraphAsymmErrors(histomu[4][sam][var], histomu[1][sam][var]);
            graph[18][var][sam] = new TGraphAsymmErrors(histomu[5][sam][var], histomu[1][sam][var]);
            graph[20][var][sam] = new TGraphAsymmErrors(histomu[7][sam][var], histomu[6][sam][var]);

            graph[15][var][sam] = new TGraphAsymmErrors(histomu[2][sam][var], histomu[1][sam][var]);
            graph[16][var][sam] = new TGraphAsymmErrors(histomu[3][sam][var], histomu[1][sam][var]);

            graph[14][var][sam] = new TGraphAsymmErrors(histomu[1][sam][var], histomu[0][sam][var]);
            graph[19][var][sam] = new TGraphAsymmErrors(histomu[6][sam][var], histomu[0][sam][var]);
            
        
            

        }

        //____ end code for 1D

        // Get 2D
        /*
        for(int i = 0; i < 10; i++){
            h_ele_pteta[i+1][sam]->Divide(h_ele_pteta[0][sam]);
            h_ele_ptactiv[i+1][sam]->Divide(h_ele_ptactiv[0][sam]);
            //h_ele_activHT[i+1][sam]->Divide(h_ele_activHT[0][sam]);
            //h_ele_activnpv[i+1][sam]->Divide(h_ele_activnpv[0][sam]);
            //h_ele_HTnpv[i+1][sam]->Divide(h_ele_HTnpv[0][sam]);
        }

        for(int i = 0; i < 4; i++){
            h_mu_pteta[i+1][sam]->Divide(h_mu_pteta[0][sam]);
            h_mu_ptactiv[i+1][sam]->Divide(h_mu_ptactiv[0][sam]);
            //h_mu_activHT[i+1][sam]->Divide(h_mu_activHT[0][sam]);
            //h_mu_activnpv[i+1][sam]->Divide(h_mu_activnpv[0][sam]);
            //h_mu_HTnpv[i+1][sam]->Divide(h_mu_HTnpv[0][sam]);
        }
        
        
        for(int j = 0; j < 4 ; j++){
            for(int i = 0; i < 10; i++)
                h_ele_pteta_2Dnpv[i+1][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);

            for(int i =0; i < 4; i++)
                h_mu_pteta_2Dnpv[i+1][sam][j]->Divide(h_mu_pteta_2Dnpv[0][sam][j]);
        }
        */
        
        h_ele_pteta[4][sam]->Divide(h_ele_pteta[0][sam]);
        h_ele_ptactiv[4][sam]->Divide(h_ele_ptactiv[0][sam]);
        h_ele_etaactiv[4][sam]->Divide(h_ele_etaactiv[0][sam]);

        h_ele_pteta[9][sam]->Divide(h_ele_pteta[0][sam]);
        h_ele_ptactiv[9][sam]->Divide(h_ele_ptactiv[0][sam]);
        h_ele_etaactiv[9][sam]->Divide(h_ele_etaactiv[0][sam]);
        

        h_mu_pteta[1][sam]->Divide(h_mu_pteta[0][sam]);
        h_mu_ptactiv[1][sam]->Divide(h_mu_ptactiv[0][sam]);
        h_mu_etaactiv[1][sam]->Divide(h_mu_etaactiv[0][sam]);

        h_mu_pteta[6][sam]->Divide(h_mu_pteta[0][sam]);
        h_mu_ptactiv[6][sam]->Divide(h_mu_ptactiv[0][sam]);
        h_mu_etaactiv[6][sam]->Divide(h_mu_etaactiv[0][sam]);
        
        /*
        h_ele_pteta[1][sam]->Divide(h_ele_pteta[0][sam]);
        h_ele_pteta[2][sam]->Divide(h_ele_pteta[0][sam]);
        h_ele_pteta[3][sam]->Divide(h_ele_pteta[0][sam]);
        h_ele_pteta[4][sam]->Divide(h_ele_pteta[0][sam]);
        h_ele_pteta[14][sam]->Divide(h_ele_pteta[0][sam]);
        
        h_ele_pteta[6][sam]->Divide(h_ele_pteta[5][sam]);
        h_ele_pteta[7][sam]->Divide(h_ele_pteta[5][sam]);
        h_ele_pteta[8][sam]->Divide(h_ele_pteta[5][sam]);

        h_ele_pteta[11][sam]->Divide(h_ele_pteta[10][sam]);
        h_ele_pteta[12][sam]->Divide(h_ele_pteta[10][sam]);
        h_ele_pteta[13][sam]->Divide(h_ele_pteta[9][sam]);
        
        h_ele_pteta[9][sam]->Divide(h_ele_pteta[0][sam]);
        h_ele_pteta[10][sam]->Divide(h_ele_pteta[0][sam]);
        h_ele_pteta[5][sam]->Divide(h_ele_pteta[0][sam]);
        

        
        h_mu_pteta[4][sam]->Divide(h_mu_pteta[1][sam]);
        h_mu_pteta[5][sam]->Divide(h_mu_pteta[1][sam]);
        h_mu_pteta[7][sam]->Divide(h_mu_pteta[6][sam]);

        h_mu_pteta[2][sam]->Divide(h_mu_pteta[1][sam]);
        h_mu_pteta[3][sam]->Divide(h_mu_pteta[1][sam]);

        h_mu_pteta[1][sam]->Divide(h_mu_pteta[0][sam]);
        h_mu_pteta[6][sam]->Divide(h_mu_pteta[0][sam]);
        */
        /*
        for(int j = 0; j < 4 ; j++){
            h_ele_pteta_2Dnpv[1][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);
            h_ele_pteta_2Dnpv[2][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);
            h_ele_pteta_2Dnpv[3][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);
            h_ele_pteta_2Dnpv[4][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);
            h_ele_pteta_2Dnpv[14][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);
        
            h_ele_pteta_2Dnpv[6][sam][j]->Divide(h_ele_pteta_2Dnpv[5][sam][j]);
            h_ele_pteta_2Dnpv[7][sam][j]->Divide(h_ele_pteta_2Dnpv[5][sam][j]);
            h_ele_pteta_2Dnpv[8][sam][j]->Divide(h_ele_pteta_2Dnpv[5][sam][j]);

            h_ele_pteta_2Dnpv[11][sam][j]->Divide(h_ele_pteta_2Dnpv[10][sam][j]);
            h_ele_pteta_2Dnpv[12][sam][j]->Divide(h_ele_pteta_2Dnpv[10][sam][j]);
            h_ele_pteta_2Dnpv[13][sam][j]->Divide(h_ele_pteta_2Dnpv[9][sam][j]);
        
            h_ele_pteta_2Dnpv[9][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);
            h_ele_pteta_2Dnpv[10][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);
            h_ele_pteta_2Dnpv[5][sam][j]->Divide(h_ele_pteta_2Dnpv[0][sam][j]);


            h_mu_pteta_2Dnpv[4][sam][j]->Divide(h_mu_pteta_2Dnpv[1][sam][j]);
            h_mu_pteta_2Dnpv[5][sam][j]->Divide(h_mu_pteta_2Dnpv[1][sam][j]);
            h_mu_pteta_2Dnpv[7][sam][j]->Divide(h_mu_pteta_2Dnpv[6][sam][j]);

            h_mu_pteta_2Dnpv[2][sam][j]->Divide(h_mu_pteta_2Dnpv[1][sam][j]);
            h_mu_pteta_2Dnpv[3][sam][j]->Divide(h_mu_pteta_2Dnpv[1][sam][j]);

            h_mu_pteta_2Dnpv[1][sam][j]->Divide(h_mu_pteta_2Dnpv[0][sam][j]);
            h_mu_pteta_2Dnpv[6][sam][j]->Divide(h_mu_pteta_2Dnpv[0][sam][j]);
        }
        */
        

        // end code for 2D

        // Get 3D _____ start
        /*
        for(int i = 0; i < 10; i++)
            h_ele_pt_eta_npv[i+1][sam]->Divide(h_ele_pt_eta_npv[0][sam]);

        for(int i = 0; i < 4; i++)
            h_mu_pt_eta_npv[i+1][sam]->Divide(h_mu_pt_eta_npv[0][sam]);
        */
        
        /*
        h_ele_pt_eta_npv[1][sam]->Divide(h_ele_pt_eta_npv[0][sam]);
        h_ele_pt_eta_npv[2][sam]->Divide(h_ele_pt_eta_npv[0][sam]);
        h_ele_pt_eta_npv[3][sam]->Divide(h_ele_pt_eta_npv[0][sam]);
        h_ele_pt_eta_npv[4][sam]->Divide(h_ele_pt_eta_npv[0][sam]);
        h_ele_pt_eta_npv[14][sam]->Divide(h_ele_pt_eta_npv[0][sam]);
        
        h_ele_pt_eta_npv[6][sam]->Divide(h_ele_pt_eta_npv[5][sam]);
        h_ele_pt_eta_npv[7][sam]->Divide(h_ele_pt_eta_npv[5][sam]);
        h_ele_pt_eta_npv[8][sam]->Divide(h_ele_pt_eta_npv[5][sam]);

        h_ele_pt_eta_npv[11][sam]->Divide(h_ele_pt_eta_npv[10][sam]);
        h_ele_pt_eta_npv[12][sam]->Divide(h_ele_pt_eta_npv[10][sam]);
        h_ele_pt_eta_npv[13][sam]->Divide(h_ele_pt_eta_npv[9][sam]);
        
        h_ele_pt_eta_npv[9][sam]->Divide(h_ele_pt_eta_npv[0][sam]);
        h_ele_pt_eta_npv[10][sam]->Divide(h_ele_pt_eta_npv[0][sam]);
        h_ele_pt_eta_npv[5][sam]->Divide(h_ele_pt_eta_npv[0][sam]);
        

        
        h_mu_pt_eta_npv[4][sam]->Divide(h_mu_pt_eta_npv[1][sam]);
        h_mu_pt_eta_npv[5][sam]->Divide(h_mu_pt_eta_npv[1][sam]);
        h_mu_pt_eta_npv[7][sam]->Divide(h_mu_pt_eta_npv[6][sam]);

        h_mu_pt_eta_npv[2][sam]->Divide(h_mu_pt_eta_npv[1][sam]);
        h_mu_pt_eta_npv[3][sam]->Divide(h_mu_pt_eta_npv[1][sam]);

        h_mu_pt_eta_npv[1][sam]->Divide(h_mu_pt_eta_npv[0][sam]);
        h_mu_pt_eta_npv[6][sam]->Divide(h_mu_pt_eta_npv[0][sam]);
        */

        //_____ end code for 3D _________________________________


        
        for(int i = 0; i < 5; i++){
            //gh_ele_npv_ptc[i][sam] = new TGraphAsymmErrors(h_ele_npv_ptc[i][sam], h_ele_npv_pt[i][sam]);
            //gh_mu_npv_ptc[i][sam] = new TGraphAsymmErrors(h_mu_npv_ptc[i][sam], h_mu_npv_pt[i][sam]);
            //gh_mu_HT_ptc[i][sam] = new TGraphAsymmErrors(h_mu_HT_ptc[i][sam], h_mu_HT_pt[i][sam]);
            gh_mu_activ_ptc[i][sam] = new TGraphAsymmErrors(h_mu_activ_ptc[i][sam], h_mu_activ_pt[i][sam]);
            
        }
        
    }
    
    std::cout << std::endl;
    /*
    TString title[5] = {"Medium ID + miniiso02", "Medium ID + multiiso", "Electron CB + miniiso01", "Electron MVA Tight + Miniiso01 + convVeto + MissInnerHits = 0", "Electron MVA T TrigIDEmu + MultiIsoTight + ISOEmu + convVeto + missing inner hits = 0 + charge consistency"};

    for(int i = 0; i < 5; i++){
        std::cout << title[i] << std::endl;
        std::cout << "HT > 80 GeV: " << 1.*activNumber[1][i]/activNumber[0][i] << std::endl;
        std::cout << "HT > 300 GeV: " << 1.*activNumber[3][i]/activNumber[2][i] << std::endl;
        std::cout << "HT > 1000 GeV: " << 1.*activNumber[5][i]/activNumber[4][i] << std::endl;
        std::cout << std::endl;
    }

    TString titlePDF[5] = {"MediumIDminiiso02", "MediumIDmultiiso", "ElectronCBminiiso01", "ElectronMVATightMiniiso01", "ElectronMVAT_TrigIDEmuMultiIsoTight"};

    TString pdf = ".pdf";

    TCanvas *c1[5];
    const char * buffer[5] = {"c1","c2","c3","c4","c5"};
    for(int i = 0; i < 5; i++){
        c1[i] = new TCanvas(buffer[i],buffer[i]);
        c1[i]->Divide(2,1);
        c1[i]->cd(1);
        c1[i]->cd(1)->SetLogy();
        h_activityHT[i][0][0]->SetMarkerColor(kBlack);
        h_activityHT[i][0][0]->SetLineColor(kBlack);
        h_activityHT[i][0][0]->DrawNormalized();
        h_activityHT[i][1][0]->SetMarkerColor(kBlue);
        h_activityHT[i][1][0]->SetLineColor(kBlue);
        h_activityHT[i][1][0]->DrawNormalized("same");

        c1[i]->cd(2);
        c1[i]->cd(2)->SetLogy();
        h_activityHT[i][0][1]->SetMarkerColor(kBlack);
        h_activityHT[i][0][1]->SetLineColor(kBlack);
        h_activityHT[i][0][1]->DrawNormalized();
        h_activityHT[i][1][1]->SetMarkerColor(kBlue);
        h_activityHT[i][1][1]->SetLineColor(kBlue);
        h_activityHT[i][1][1]->DrawNormalized("same");
        c1[i]->SaveAs(titlePDF[i] + pdf);
    }
    */

    
    TString idiso[21] = {"el_tightCB", "el_mediumCB", "el_looseCB", "el_vetoCB", 
                            "el_vloose", "el_mini01", "el_mini04", "el_inhit_less2", 
                            "el_tight2d3d","el_tight2d3dIDEmu", "el_multi", "el_multiISO", 
                            "el_inhit_eq0",
                            "el_vlooseFO",
                            "mu_looseID","mu_tightIP2D", "mu_tightIP3D",
                            "mu_mini02", "mu_mini04", 
                            "mu_mediumID", "mu_multi"};
              
                            
     /*
    TString idiso[14] = {"el_vetoCB_mini01", "el_looseCB_mini01", "el_looseCB_mini04", "el_mediumCB_mini01", "el_tightCB_mini01",
                            "el_vloose_mini01", "el_tight_mini01", "el_tight_IDEmu_ra5", "el_tight_IDEmu_ISOEMu_ra5", "el_vlooseFO_ra5", 
                            "mu_looseID_mini02", "mu_mediumID_mini02", "mu_mediumID_multi", "mu_loose_mini04"
                            };
    
       */                     

    TString npv_reg[4] = {"_npv_1_9", "_npv_10_13", "_npv_14_17", "_npv_18_40"};
    
    
    TLegend *legel = new TLegend(0.6, 0.15, 0.9, 0.45, "Electron");
    legel->SetBorderSize(0);
    legel->SetFillColor(kWhite);
    legel->SetFillStyle(0);
    legel->AddEntry(graph[0][0][0], "Fastsim", "l");
    legel->AddEntry(graph[0][0][1], "Fullsim", "l");

    TLegend *legmu = new TLegend(0.6, 0.15, 0.9, 0.45, "Muon");
    legmu->SetBorderSize(0);
    legmu->SetFillColor(kWhite);
    legmu->SetFillStyle(0);
    legmu->AddEntry(graph[0][0][0], "Fastsim", "l");
    legmu->AddEntry(graph[0][0][1], "Fullsim", "l");

    TLegend *legtau = new TLegend(0.2, 0.15, 0.5, 0.45, "Tau");
    legtau->SetBorderSize(0);
    legtau->SetFillColor(kWhite);
    legtau->AddEntry(graph[0][0][0], "Fastsim", "l");
    legtau->AddEntry(graph[0][0][1], "Fullsim", "l");

    TLegend *legeta = new TLegend(0.6, 0.55, 0.75, 0.75, "activity");
    legeta->SetBorderSize(0);
    legeta->SetFillColor(kWhite);
    legeta->AddEntry(graph[0][0][0], "ttbar", "l");
    legeta->AddEntry(graph[0][0][1], "DY", "l");

    //TCanvas *c_test = new TCanvas("c_test", "c_test");
    /*
    g_test[0]->Draw();
    g_test[1]->SetLineColor(kRed);
    g_test[1]->Draw("same");
    */
    //h_test[1][0]->DrawNormalized();
    //h_test[1][1]->SetLineColor(kRed);
    //h_test[1][1]->DrawNormalized("same");

    /*
    TCanvas *c_test1 = new TCanvas("c_test1", "c_test1");
    c_test1->Divide(2,1);
    c_test1->cd(1);
    h_test[0][1]->SetTitle("#eta distribution");
    h_test[0][1]->SetMarkerColor(kBlack);
    h_test[0][1]->Draw();
    h_test[1][1]->Draw("same");

    std::cout << h_test[0][1]->Integral(0,10) << " " << h_test[1][1]->Integral(0,10) << std::endl;

    c_test1->cd(2);
    h_test[2][1]->SetTitle("#eta distribution");
    h_test[2][1]->SetMarkerColor(kBlack);
    h_test[2][1]->Draw();
    h_test[3][1]->Draw("same");

    std::cout << h_test[2][1]->Integral(0,10) << " " << h_test[3][1]->Integral(0,10) << std::endl;
    */
    //legeta->Draw("same");
    
    /*
    for(int i = 1; i < 15; i++)
        Print3D(h_ele_pt_eta_npv[i][1], h_ele_pt_eta_npv[i][0], idiso[i-1]);

    for(int i = 1; i < 8; i++)
        Print3D(h_mu_pt_eta_npv[i][1], h_mu_pt_eta_npv[i][0], idiso[i+14-1]);
    */
    
    //for(int i = 1; i < 15; i++){
    
   //for(int i = 1; i < 11; i++){

        //Print2D(h_ele_ptactiv[i][1], h_ele_ptactiv[i][0], "ptActiv_" + idiso[i-1]);
        //Print2D(h_ele_pteta[i][1], h_ele_pteta[i][0], idiso[i-1]);
        /*
        Print2D(h_ele_activHT[i][1], h_ele_activHT[i][0], "activHT_" + idiso[i-1]);
        Print2D(h_ele_activnpv[i][1], h_ele_activnpv[i][0], "activNPV_" + idiso[i-1]);
        Print2D(h_ele_HTnpv[i][1], h_ele_HTnpv[i][0], "HTnpv_" + idiso[i-1]);
        */
        /*
        for(int j = 0; j < 4; j++)
            Print2D(h_ele_pteta_2Dnpv[i][1][j], h_ele_pteta_2Dnpv[i][0][j], idiso[i-1] + npv_reg[j]);
        
        
        TFile *f = TFile::Open("files/sf_"+idiso[i-1]+".root","recreate");
        h_ele_pt_eta_npv[i][1]->Divide(h_ele_pt_eta_npv[i][0]);
        h_ele_pt_eta_npv[i][1]->SetName("histo3D");
        h_ele_pt_eta_npv[i][1]->Write();
        h_ele_pteta[i][1]->Divide(h_ele_pteta[i][0]);
        h_ele_pteta[i][1]->SetName("histo2D");
        h_ele_pteta[i][1]->Write();
        f->Close();
        */
        
        
        //Print3D(h_ele_pt_eta_npv[i][1], h_ele_pt_eta_npv[i][0], idiso[i-1]);
        
    //}
    
    //for(int i = 1; i < 8; i++){
    //for(int i = 1; i < 5; i++){
    
    int i = 6;
        Print2D(h_mu_ptactiv[i][1], h_mu_ptactiv[i][0], "ptActiv_" + idiso[i+14-1] + "_" + sampleName);
        Print2D(h_mu_pteta[i][1], h_mu_pteta[i][0], "pteta_" + idiso[i+14-1] + "_" +  sampleName);
        Print2D(h_mu_etaactiv[i][1], h_mu_etaactiv[i][0], "etaActiv_" + idiso[i+14-1] + "_" +  sampleName);

    i = 1;

        Print2D(h_mu_ptactiv[i][1], h_mu_ptactiv[i][0], "ptActiv_" + idiso[i+14-1] + "_" + sampleName);
        Print2D(h_mu_pteta[i][1], h_mu_pteta[i][0], "pteta_" + idiso[i+14-1] + "_" +  sampleName);
        Print2D(h_mu_etaactiv[i][1], h_mu_etaactiv[i][0], "etaActiv_" + idiso[i+14-1] + "_" +  sampleName);
    
        
    i = 4;
        Print2D(h_ele_ptactiv[i][1], h_ele_ptactiv[i][0], "ptActiv_" + idiso[i-1] + "_" +  sampleName);
        Print2D(h_ele_pteta[i][1], h_ele_pteta[i][0], "pteta_" + idiso[i-1] + "_" +  sampleName);
        Print2D(h_ele_etaactiv[i][1], h_ele_etaactiv[i][0], "etaActiv_" + idiso[i-1] + "_" +  sampleName);
    
        
    i = 9;
        Print2D(h_ele_ptactiv[i][1], h_ele_ptactiv[i][0], "ptActiv_" + idiso[i-1] + "_" +  sampleName);
        Print2D(h_ele_pteta[i][1], h_ele_pteta[i][0], "pteta_" + idiso[i-1] + "_" +  sampleName);
        Print2D(h_ele_etaactiv[i][1], h_ele_etaactiv[i][0], "etaActiv_" + idiso[i-1] + "_" +  sampleName);
    

    
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->Divide(2,1);
    c1->cd(1);
    h_deltaR[0][1]->Draw();

    c1->cd(2);
    h_deltaR[1][1]->Draw();

    //c1->cd(3);
    //h_deltaR[2][1]->Draw();
    
    
    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->Divide(2,2);
    c2->cd(1);
    h_deltaPT[0][1]->SetTitle("Electron, dR < 0.2");
    h_deltaPT[0][1]->Draw();

    c2->cd(2);
    h_deltaPT[1][1]->SetTitle("Muon, dR < 0.2");
    h_deltaPT[1][1]->Draw();

    c2->cd(3);
    h_deltaPT[2][1]->SetTitle("Electron, dR > 0.2");
    h_deltaPT[2][1]->Draw();

    c2->cd(4);
    h_deltaPT[3][1]->SetTitle("Muon, dR > 0.2");
    h_deltaPT[3][1]->Draw();
    
    

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 400);

    c3->Divide(2,2);
    c3->cd(1);

    h_activ_distr[0][0][1]->SetTitle("Electron, dR < 0.2");
    h_activ_distr[0][0][1]->Draw();

    c3->cd(2);

    h_activ_distr[1][0][1]->SetTitle("Electron, dR > 0.2");
    h_activ_distr[1][0][1]->Draw();

    c3->cd(3);

    h_activ_distr[0][1][1]->SetTitle("Muon, dR < 0.2");
    h_activ_distr[0][1][1]->Draw();

    c3->cd(4);

    h_activ_distr[1][1][1]->SetTitle("Muon, dR > 0.2");
    h_activ_distr[1][1][1]->Draw();
    

    /*
    TCanvas *cactiv_distr[2];// = new TCanvas("cactiv_distr", "cactiv_distr",600,900); 
    
    
    for(int i = 0; i < 2; i++){
        TString name = "cactiv_distr" + std::to_string(i);
        cactiv_distr[i] = new TCanvas(name, name,600,400);
        cactiv_distr[i]->Divide(2,1);
        TLine *line = new TLine(0., 1., 0.5, 1.);
        //electrons
        cactiv_distr[i]->cd(1);
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        
        h_activ_distr[i][0][1]->Scale(1./h_activ_distr[i][0][1]->Integral());
        h_activ_distr[i][0][1]->SetTitle("electrons");
        h_activ_distr[i][0][1]->Draw();

        h_activ_distr[i][0][0]->Scale(1./h_activ_distr[i][0][0]->Integral());
        h_activ_distr[i][0][0]->Draw("same");
        legeta->Draw("same");
        
        cactiv_distr[i]->cd(1);        // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.2);
        //pad2->SetGridx(); // vertical grid
        pad2->Draw();
        pad2->cd();   

        TH1F *h = (TH1F*)h_activ_distr[i][0][1]->Clone("h");
        h->SetLineColor(kBlack);
        //h->SetMinimum(0.8);  // Define Y ..
        //h->SetMaximum(1.35); // .. range
        h->Sumw2();
        h->SetStats(0);      // No statistics on lower plot
        h->Divide(h_activ_distr[i][0][0]);
        
        h->SetTitle("");
        h->SetMaximum(2.);
        h->SetMinimum(0.);
        h->SetMarkerStyle(21);
        h->SetMarkerSize(0.3);
        h->GetXaxis()->SetLabelSize(0.15);
        h->GetXaxis()->SetTitleSize(0.2);
        h->GetXaxis()->SetTitleOffset(0.5);

        h->GetYaxis()->SetLabelSize(0.15);
        h->GetYaxis()->SetTitleSize(0.2);
        h->GetYaxis()->SetTitleOffset(0.2);

        h->Draw("ep"); 
        line->Draw("same"); 

        //muons
        cactiv_distr[i]->cd(2);

        TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
        pad11->SetBottomMargin(0); // Upper and lower plot are joined
        pad11->Draw();             // Draw the upper pad: pad1
        pad11->cd();               // pad1 becomes the current pad

        h_activ_distr[i][1][1]->Scale(1./h_activ_distr[i][1][1]->Integral());
        h_activ_distr[i][1][1]->SetTitle("muons");
        h_activ_distr[i][1][1]->Draw();

        h_activ_distr[i][1][0]->Scale(1./h_activ_distr[i][1][0]->Integral());
        h_activ_distr[i][1][0]->DrawNormalized("same");

        legeta->Draw("same");

        cactiv_distr[i]->cd(2);        // Go back to the main canvas before defining pad2
        TPad *pad22 = new TPad("pad22", "pad22", 0, 0.05, 1, 0.3);
        pad22->SetTopMargin(0);
        pad22->SetBottomMargin(0.2);
        //pad22->SetGridx(); // vertical grid
        pad22->Draw();
        pad22->cd();   

        TH1F *hh = (TH1F*)h_activ_distr[i][1][1]->Clone("h");
        hh->SetLineColor(kBlack);
        //hh->SetMinimum(0.8);  // Define Y ..
        //hh->SetMaximum(1.35); // .. range
        hh->Sumw2();
        hh->SetStats(0);      // No statistics on lower plot
        hh->Divide(h_activ_distr[i][1][0]);
        
        hh->SetTitle("");
        hh->SetMaximum(2.);
        hh->SetMinimum(0.);
        hh->SetMarkerStyle(21);
        hh->SetMarkerSize(0.3);
        hh->GetXaxis()->SetLabelSize(0.15);
        hh->GetXaxis()->SetTitleSize(0.2);
        hh->GetXaxis()->SetTitleOffset(0.5);

        hh->GetYaxis()->SetLabelSize(0.15);
        hh->GetYaxis()->SetTitleSize(0.2);
        hh->GetYaxis()->SetTitleOffset(0.2);
        hh->Draw("ep");  
        line->Draw("same"); 

    }
    */
        /*
        Print2D(h_mu_pteta[i][1], h_mu_pteta[i][0], idiso[i+10-1]);
        for(int j = 0; j < 4; j++)
            Print2D(h_mu_pteta_2Dnpv[i][1][j], h_mu_pteta_2Dnpv[i][0][j], idiso[i+10-1] + npv_reg[j]);
        */
        /*
        Print2D(h_mu_activHT[i][1], h_mu_activHT[i][0], "activHT_" + idiso[i+10-1]);
        Print2D(h_mu_activnpv[i][1], h_mu_activnpv[i][0], "activNPV_" + idiso[i+10-1]);
        Print2D(h_mu_HTnpv[i][1], h_mu_HTnpv[i][0], "HTnpv_" + idiso[i+10-1]);
        */
        //TFile *f = TFile::Open("files/sf_"+idiso[i+10-1]+".root","recreate");
        /*
        TFile *f = TFile::Open("files/sf_"+idiso[i+14-1]+".root","recreate");
        h_mu_pt_eta_npv[i][1]->Divide(h_mu_pt_eta_npv[i][0]);
        h_mu_pt_eta_npv[i][1]->SetName("histo3D");
        h_mu_pt_eta_npv[i][1]->Write();
        h_mu_pteta[i][1]->Divide(h_mu_pteta[i][0]);
        h_mu_pteta[i][1]->SetName("histo2D");
        h_mu_pteta[i][1]->Write();
        f->Close();
        */
        
        
        /*
        Print2D(h_mu_pteta[i][1], h_mu_pteta[i][0], idiso[i+14-1]);
        
        for(int j = 0; j < 4; j++)
            Print2D(h_mu_pteta_2Dnpv[i][1][j], h_mu_pteta_2Dnpv[i][0][j], idiso[i+14-1] + npv_reg[j]);
         */

    //}
    

    
    //for(int i = 1; i < 11; i++){
    //for(int i = 1; i < 15; i++){
        //if(i != 8 && i != 9 && i != 10) continue;
        
        //Print2D(h_ele_pteta[i][1], h_ele_pteta[i][0], idiso[i-1]);
        //Print2D(h_ele_ptactiv[i][1], h_ele_ptactiv[i][0], "ptActiv_" + idiso[i-1]);
        /*
        for(int j = 0; j < 4; j++)
            Print2D(h_ele_pteta_2Dnpv[i][1][j], h_ele_pteta_2Dnpv[i][0][j], idiso[i-1] + npv_reg[j]);
            */
        /*
        TFile *f = TFile::Open("FinalResult/files/sf_"+idiso[i-1]+".root","recreate");
        h_ele_pt_eta_npv[i][1]->Divide(h_ele_pt_eta_npv[i][0]);
        h_ele_pt_eta_npv[i][1]->SetName("histo3D");
        h_ele_pt_eta_npv[i][1]->Write();
        h_ele_pteta[i][1]->Divide(h_ele_pteta[i][0]);
        h_ele_pteta[i][1]->SetName("histo2D");
        h_ele_pteta[i][1]->Write();
        f->Close();
        */
    //}

        
    
    //for(int i = 1; i < 8; i++){
        
        //Print2D(h_mu_pteta[i][1], h_mu_pteta[i][0], idiso[i+14-1]);
        //Print2D(h_mu_ptactiv[i][1], h_mu_ptactiv[i][0], "ptActiv_" + idiso[i+10-1]);
        
        /*
        for(int j = 0; j < 4; j++)
            Print2D(h_mu_pteta_2Dnpv[i][1][j], h_mu_pteta_2Dnpv[i][0][j], idiso[i+10-1] + npv_reg[j]);
         
        TFile *f = TFile::Open("FinalResult/files/sf_"+idiso[i+10-1]+".root","recreate");
        h_mu_pt_eta_npv[i][1]->Divide(h_mu_pt_eta_npv[i][0]);
        h_mu_pt_eta_npv[i][1]->SetName("histo3D");
        h_mu_pt_eta_npv[i][1]->Write();
        h_mu_pteta[i][1]->Divide(h_mu_pteta[i][0]);
        h_mu_pteta[i][1]->SetName("histo2D");
        h_mu_pteta[i][1]->Write();
        f->Close();
        */
    //}
    
    /*
    for(int i = 10; i < 14; i++){
        Print1D(graph[i], idiso[i], "muon ID/iso efficiency", "muon", legmu);
    }

    for(int i = 0; i < 10; i++){
        Print1D(graph[i], idiso[i], "electron ID/iso efficiency", "electron", legel);
    }
    */
    

    //Print1D(graph[8], "test", "electron ID efficiency", "electron", legel);
    

    /*
    for(int i = 14; i < 21; i++){
        if(i == 14 || i == 19)
            Print1D(graph[i], idiso[i], "muon ID efficiency", "muon", legmu);
        else
            Print1D(graph[i], idiso[i], "muon iso efficiency", "muon", legmu);
    }

    for(int i = 0; i < 14; i++){
        if((i > -1 && i < 5) || i == 7 || i == 8 || i == 9 || i == 12 || i == 13)
            Print1D(graph[i], idiso[i], "electron ID efficiency", "electron", legel);
        else
            Print1D(graph[i], idiso[i], "electron iso efficiency", "electron", legel);
    }
    */

//__________________________________________________________________________________________________________________
 
    string TitleYhistEl = "electron ID efficiency";
    string TitleYhistMu = "muon ID efficiency";
    //string TitleYhistTau = "tau ID efficiency";
    //string TitleYhistEl = "electron iso efficiency";
    //string TitleYhistMu = "muon iso efficiency";
    //string TitleYhistTau = "tau iso efficiency";
    string TitleYSF = "Fullsim/Fastsim";

    /*
    TCanvas *c1 = new TCanvas("c1", "electron kin", 0, 0, 1000, 700);
    
    showHistUP(pt_min, pt_max, graph[0][0], "Electron p_{T} (GeV)", TitleYhistEl, 0., 0.5, "", 0.06, 1.2, 0.05, 0.005);
    legel->Draw("same");
    c1->cd();         
    showHistDOWN(pt_min, pt_max,graph[0][0], "Electron p_{T} (GeV)",TitleYSF ,0.,0.5, 0.7, 1.2, "", 26, 1.85, 22., 0.005);

    c1->cd();
    showHistUP(eta_min, eta_max, graph[0][1], "Electron #eta", TitleYhistEl, 0.5, 1.0, "", 0.06, 1.2, 0.05, 0.005);
    legel->Draw("same");
    c1->cd();         
    showHistDOWN(eta_min, eta_max,graph[0][1], "Electron #eta",TitleYSF ,0.5,1.0, 0.7, 1.2, "", 26, 1.85, 22., 0.005);
    
    c1->SaveAs("fastfull/electron.pdf");
    */
//__________________________________________________________________________________________________________________

    /*
    TCanvas *c2 = new TCanvas("c2", "muon kin", 0, 0, 1000, 700);
    showHistUP(pt_min, pt_max, graph[1][0], "Muon p_{T} (GeV)", TitleYhistMu, 0., 0.5, "", 0.06, 1.2, 0.05, 0.005);
    legmu->Draw("same");
    c2->cd();         
    showHistDOWN(pt_min, pt_max,graph[1][0], "Muon p_{T} (GeV)",TitleYSF ,0.,0.5, 0.8, 1.2, "", 26, 1.85, 22., 0.005);

    c2->cd();
    showHistUP(eta_min, eta_max, graph[1][1], "Muon #eta", TitleYhistMu, 0.5, 1.0, "", 0.06, 1.2, 0.05, 0.005);
    legmu->Draw("same");
    c2->cd();          // Go back to the main canvas before defining pad2
    showHistDOWN(eta_min, eta_max,graph[1][1], "Muon #eta",TitleYSF ,0.5,1.0, 0.8, 1.2, "", 26, 1.85, 22., 0.005);

    c2->SaveAs("fastfull/muon.pdf");
    */
    
//__________________________________________________________________________________________________________________

/*
    TCanvas *c22 = new TCanvas("c22", "tau kin", 0, 0, 1000, 700);
    showHistUP(pt_min, pt_max, gh_tau_ptc, "Tau p_{T} (GeV)", TitleYhistTau, 0., 0.5);
    legtau->Draw("same");
    c22->cd();          // Go back to the main canvas before defining pad2
    showHistDOWN(numberBinspt, gh_tau_ptc, pt_min, pt_max,"Tau p_{T} (GeV)",TitleYSF ,0.,0.5, 0, 4.5);

    c22->cd();
    showHistUP(eta_minfull, eta_max, gh_tau_etac, "Tau #eta", TitleYhistTau, 0.5, 1.0);
    legtau->Draw("same");
    c22->cd();          // Go back to the main canvas before defining pad2
    showHistDOWN(numberBinsetafull+1, gh_tau_etac, eta_minfull, eta_max,"Tau #eta",TitleYSF ,0.5,1.0, 0, 4.5);

    //c22->SaveAs("fastfull/tau.pdf");
    */

//__________________________________________________________________________________________________________________

    /*
    TCanvas *c3 = new TCanvas("c3", "ele # PV Activity HT Njets", 0, 0, 1000, 700);
    c3->Divide(2,2);
    c3->cd(1);
    showHistUP(npv_min, npv_max, graph[0][4], "# PV", TitleYhistEl, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    legel->Draw("same");
    c3->cd(1);          
    showHistDOWN(npv_min, npv_max,graph[0][4], "# PV",TitleYSF ,0.,1., 0.95, 1.05, "", 16., 2.4, 16, 0.005);

    c3->cd(2);
    showHistUP(Activb_min, Activb_max, graph[0][5], "Activity", TitleYhistEl, 0, 1., "log", 0.09, 0.7, 0.09, 0.005);
    legel->Draw("same");
    c3->cd(2);          
    showHistDOWN(Activb_min, Activb_max,graph[0][5], "Activity",TitleYSF ,0.,1., 0.9, 1.1, "log", 16., 2.4, 16, 0.005);
    */

//---------
    /*
    c3->cd(3);          
    showHistUP(HT_min, HT_max, graph[0][2], "HT (GeV)", TitleYhistEl, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    legel->Draw("same");
    c3->cd(3);          // Go back to the main canvas before defining pad2
    showHistDOWN(HT_min, HT_max,graph[0][2], "HT (GeV)",TitleYSF ,0.,1., 0.95, 1.05, "", 16., 2.4, 16, 0.005);

    c3->cd(4);          
    showHistUP(Njets_min, Njets_max, graph[0][3], "Njets", TitleYhistEl, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    legel->Draw("same");
    c3->cd(4);         
    showHistDOWN(Njets_min, Njets_max, graph[0][3], "Njets",TitleYSF ,0.,1., 0.95, 1.05, "", 16., 2.4, 16, 0.005);
    
    c3->SaveAs("fastfull/ele_npv_activ_HT_njets.pdf");
    */
    
 //__________________________________________________________________________________________________________________
      
    /*
    TCanvas *c4 = new TCanvas("c4", "mu # PV Activity HT Njets", 0, 0, 1000, 700);
    c4->Divide(2,2);
    c4->cd(1);
    showHistUP(npv_min, npv_max, graph[1][4], "# PV", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    legmu->Draw("same");
    c4->cd(1);          
    showHistDOWN(npv_min, npv_max,graph[1][4], "# PV",TitleYSF ,0.,1., 0.95, 1.05, "", 16., 2.4, 16, 0.005);

    c4->cd(2);
    showHistUP(Activb_min, Activb_max, graph[1][5], "Activity", TitleYhistMu, 0., 1., "log", 0.09, 0.7, 0.09, 0.005);
    legmu->Draw("same");
    c4->cd(2);          
    showHistDOWN(Activb_min, Activb_max,graph[1][5], "Activity",TitleYSF ,0.,1., 0.9, 1.1, "log", 16., 2.4, 16, 0.005);
    

    c4->cd(3);          
    showHistUP(HT_min, HT_max, graph[1][2], "HT (GeV)", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    legmu->Draw("same");
    c4->cd(3);          // Go back to the main canvas before defining pad2
    showHistDOWN(HT_min, HT_max,graph[1][2], "HT (GeV)",TitleYSF ,0.,1., 0.95, 1.05, "", 16., 2.4, 16, 0.005);

    c4->cd(4);          // Go back to the main canvas before defining pad2
    showHistUP(Njets_min, Njets_max, graph[1][3], "Njets", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    legmu->Draw("same");
    c4->cd(4);          // Go back to the main canvas before defining pad2
    showHistDOWN(Njets_min, Njets_max, graph[1][3], "Njets",TitleYSF ,0.,1., 0.95, 1.05, "", 16., 2.4, 16, 0.005);

    c4->SaveAs("fastfull/mu_npv_activ_HT_njets.pdf");
    */
//__________________________________________________________________________________________________________________
    
    /*
    TLine *line1 = new TLine(0.65,0.15,0.65,0.9);
    line1->SetNDC();
    line1->SetLineColor(kRed);
    line1->SetLineWidth(2);
    TLine *line2 = new TLine(0.63,0.15,0.63,0.9);
    line2->SetNDC();
    line2->SetLineColor(kRed);
    line2->SetLineWidth(2);
    TLine *line3 = new TLine(0.27,0.15,0.27,0.9);
    line3->SetNDC();
    line3->SetLineColor(kRed);
    line3->SetLineWidth(2);
    TLine *line4 = new TLine(0.27,0.15,0.27,0.9);
    line4->SetNDC();
    line4->SetLineColor(kRed);
    line4->SetLineWidth(2);
    TText *xlabel1 = new TText(.5, .5, "ptratio < 0.8");
    xlabel1->SetTextFont(43);
    xlabel1->SetTextSize(40);
    xlabel1->SetNDC();
    
    TText *xlabel2 = new TText(.5, .5, "ptratio < 0.76");
    xlabel2->SetTextFont(43);
    xlabel2->SetTextSize(40);
    xlabel2->SetNDC();

    TText *xlabel_fabe = new TText(.2, .2, "Fastsim eff: 0.83");
    xlabel_fabe->SetNDC();
    TText *xlabel_fube = new TText(.2, .3, "Fullsim eff: 0.68");
    xlabel_fube->SetNDC();

    TText *xlabel_fage = new TText(.2, .2, "Fastsim eff: 0.95");
    xlabel_fage->SetNDC();
    TText *xlabel_fuge = new TText(.2, .3, "Fullsim eff: 0.95");
    xlabel_fuge->SetNDC();

    TText *xlabel_fabm = new TText(.2, .2, "Fastsim eff: 0.84");
    xlabel_fabm->SetNDC();
    TText *xlabel_fubm = new TText(.2, .3, "Fullsim eff: 0.77");
    xlabel_fubm->SetNDC();

    TText *xlabel_fagm = new TText(.2, .2, "Fastsim eff: 0.96");
    xlabel_fagm->SetNDC();
    TText *xlabel_fugm = new TText(.2, .3, "Fullsim eff: 0.96");
    xlabel_fugm->SetNDC();
    */

/*
    TCanvas *c5el = new TCanvas("c5el", "el", 0, 0, 1000, 700);
    c5el->Divide(2,2);
    c5el->cd(1);
    h_ptratio_bad[0]->SetTitle("ptratio bad");
    h_ptratio_bad[0]->DrawNormalized();
    h_ptratio_bad[2]->DrawNormalized("same");
    line1->Draw("same");
    legel->Draw("same");
    xlabel_fabe->Draw("same");
    xlabel_fube->Draw("same");
    c5el->cd(2);
    h_ptrel_bad[0]->SetTitle("ptrel bad");
    h_ptrel_bad[0]->DrawNormalized();
    h_ptrel_bad[2]->DrawNormalized("same");
    line3->Draw("same");
    xlabel1->Draw("same");
    c5el->cd(3);
    h_ptratio_good[0]->SetTitle("ptratio good");
    h_ptratio_good[0]->DrawNormalized();
    h_ptratio_good[2]->DrawNormalized("same");
    xlabel_fage->Draw("same");
    xlabel_fuge->Draw("same");
    line1->Draw("same");
    c5el->cd(4);
    h_ptrel_good[0]->SetTitle("ptrel good");
    h_ptrel_good[0]->DrawNormalized();
    h_ptrel_good[2]->DrawNormalized("same");
    line4->Draw("same");
    xlabel1->Draw("same");
    */


     
//__________________________________________________________________________________________________________________
    /*
    TCanvas *c5mu = new TCanvas("c5mu", "mu", 0, 0, 1000, 700);
    c5mu->Divide(2,2);
    c5mu->cd(1);
    h_ptratio_bad_mu[0]->SetTitle("ptratio bad");
    h_ptratio_bad_mu[0]->DrawNormalized();
    h_ptratio_bad_mu[2]->DrawNormalized("same");
    line2->Draw("same");
    legmu->Draw("same");
    xlabel_fabm->Draw("same");
    xlabel_fubm->Draw("same");
    c5mu->cd(2);
    h_ptrel_bad_mu[0]->SetTitle("ptrel bad");
    h_ptrel_bad_mu[0]->DrawNormalized();
    h_ptrel_bad_mu[2]->DrawNormalized("same");
    line3->Draw("same");
    xlabel2->Draw("same");
    c5mu->cd(3);
    h_ptratio_good_mu[0]->SetTitle("ptratio good");
    h_ptratio_good_mu[0]->DrawNormalized();
    h_ptratio_good_mu[2]->DrawNormalized("same");
    xlabel_fagm->Draw("same");
    xlabel_fugm->Draw("same");
    line2->Draw("same");
    c5mu->cd(4);
    h_ptrel_good_mu[0]->SetTitle("ptrel good");
    h_ptrel_good_mu[0]->DrawNormalized();
    h_ptrel_good_mu[2]->DrawNormalized("same");
    line4->Draw("same");
    xlabel2->Draw("same");
    */
    

//__________________________________________________________________________________________________________________

    
    TText *xlabel_el_pt11 = new TText(.5, .2, "[10 20]");
    xlabel_el_pt11->SetTextSize(0.1);
    xlabel_el_pt11->SetNDC();

    TText *xlabel_el_pt22 = new TText(.5, .2, "[20 30]");
    xlabel_el_pt22->SetTextSize(0.1);
    xlabel_el_pt22->SetNDC();

    TText *xlabel_el_pt33 = new TText(.5, .2, "[30 40]");
    xlabel_el_pt33->SetTextSize(0.1);
    xlabel_el_pt33->SetNDC();

    TText *xlabel_el_pt44 = new TText(.5, .2, "[40 50]");
    xlabel_el_pt44->SetTextSize(0.1);
    xlabel_el_pt44->SetNDC();

    TText *xlabel_el_pt55 = new TText(.5, .2, "[50 200]");
    xlabel_el_pt55->SetTextSize(0.1);
    xlabel_el_pt55->SetNDC();


    double graphpt[5] = {15., 25., 35., 45., 125.};
    double graphpter[5] = {5., 5., 5., 5., 75.};
    double graphyer[5] = {0., 0., 0., 0., 0.};
    double elygraph[5] = {0.};

/*
    TCanvas *c3_npv_pt = new TCanvas("c3_npv_pt", "ele # PV # pt", 0, 0, 1000, 700);
    c3_npv_pt->Divide(3,2);
    c3_npv_pt->cd(1); 
    showHistUP(npv_min, npv_max, gh_ele_npv_ptc[0], "# PV", TitleYhistEl, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt11->Draw("same");
    //legel->Draw("same");
    c3_npv_pt->cd(1);          
    showHistDOWN(npv_min, npv_max,gh_ele_npv_ptc[0], "# PV",TitleYSF ,0.,1., 0.7, 1.1, "", 16., 2.4, 16, 0.005, elygraph[0]);

    c3_npv_pt->cd(2); 
    showHistUP(npv_min, npv_max, gh_ele_npv_ptc[1], "# PV", TitleYhistEl, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt22->Draw("same");
    //legel->Draw("same");
    c3_npv_pt->cd(2);          
    showHistDOWN(npv_min, npv_max,gh_ele_npv_ptc[1], "# PV",TitleYSF ,0.,1., 0.8, 1.1, "", 16., 2.4, 16, 0.005, elygraph[1]);

    c3_npv_pt->cd(3); 
    showHistUP(npv_min, npv_max, gh_ele_npv_ptc[2], "# PV", TitleYhistEl, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt33->Draw("same");
    //legel->Draw("same");
    c3_npv_pt->cd(3);          
    showHistDOWN(npv_min, npv_max,gh_ele_npv_ptc[2], "# PV",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, elygraph[2]);


    c3_npv_pt->cd(4);  
    showHistUP(npv_min, npv_max, gh_ele_npv_ptc[3], "# PV", TitleYhistEl, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt44->Draw("same");
    //legel->Draw("same");
    c3_npv_pt->cd(4);          
    showHistDOWN(npv_min, npv_max,gh_ele_npv_ptc[3], "# PV",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, elygraph[3]);

    c3_npv_pt->cd(5); 
    showHistUP(npv_min, npv_max, gh_ele_npv_ptc[4], "# PV", TitleYhistEl, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt55->Draw("same");
    //legel->Draw("same");
    c3_npv_pt->cd(5);          
    showHistDOWN(npv_min, npv_max,gh_ele_npv_ptc[4], "# PV",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, elygraph[4]);

    
    c3_npv_pt->cd(6);
    TGraphErrors *grel = new TGraphErrors(5, graphpt, elygraph, graphpter, graphyer);
    grel->SetTitle("");
    grel->GetXaxis()->SetTitle("p_{T}, [GeV]");
    grel->GetYaxis()->SetTitle("slope");
    grel->Draw("ap");
    legel->Draw("same");
    c3_npv_pt->SaveAs("fastfull/ele_npv_pt_5range.pdf");
    */
/*
    TCanvas *c3_HT_pt = new TCanvas("c3_HT_pt", "ele HT # pt", 0, 0, 1000, 700);
    c3_HT_pt->Divide(3,2);
    c3_HT_pt->cd(1); 
    showHistUP(HT_min, HT_max, gh_mu_HT_ptc[0], "H_{T}", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt11->Draw("same");
    //legel->Draw("same");
    c3_HT_pt->cd(1);          
    showHistDOWN(HT_min, HT_max,gh_mu_HT_ptc[0], "H_{T}",TitleYSF ,0.,1., 0.7, 1.1, "", 16., 2.4, 16, 0.005, elygraph[0]);

    c3_HT_pt->cd(2); 
    showHistUP(HT_min, HT_max, gh_mu_HT_ptc[1], "H_{T}", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt22->Draw("same");
    //legel->Draw("same");
    c3_HT_pt->cd(2);          
    showHistDOWN(HT_min, HT_max,gh_mu_HT_ptc[1], "H_{T}",TitleYSF ,0.,1., 0.8, 1.1, "", 16., 2.4, 16, 0.005, elygraph[1]);

    c3_HT_pt->cd(3); 
    showHistUP(HT_min, HT_max, gh_mu_HT_ptc[2], "H_{T}", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt33->Draw("same");
    //legel->Draw("same");
    c3_HT_pt->cd(3);          
    showHistDOWN(HT_min, HT_max,gh_mu_HT_ptc[2], "H_{T}",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, elygraph[2]);


    c3_HT_pt->cd(4);  
    showHistUP(HT_min, HT_max, gh_mu_HT_ptc[3], "H_{T}", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt44->Draw("same");
    //legel->Draw("same");
    c3_HT_pt->cd(4);          
    showHistDOWN(HT_min, HT_max,gh_mu_HT_ptc[3], "H_{T}",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, elygraph[3]);

    c3_HT_pt->cd(5); 
    
    showHistUP(HT_min, HT_max, gh_mu_HT_ptc[4], "H_{T}", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt55->Draw("same");
    //legel->Draw("same");
    c3_HT_pt->cd(5);          
    showHistDOWN(HT_min, HT_max,gh_mu_HT_ptc[4], "H_{T}",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, elygraph[4]);

    
    c3_HT_pt->cd(6);
    
    TGraphErrors *grel = new TGraphErrors(5, graphpt, elygraph, graphpter, graphyer);
    grel->SetTitle("");
    grel->GetXaxis()->SetTitle("p_{T}, [GeV]");
    grel->GetYaxis()->SetTitle("slope");
    grel->Draw("ap");
    legel->Draw("same");
    c3_HT_pt->SaveAs("fastfull/ele_HT_pt_5range.pdf");
    */

    /*
    TCanvas *c3_activ_pt = new TCanvas("c3_activ_pt", "ele activ # pt", 0, 0, 1000, 700);
    c3_activ_pt->Divide(3,2);
    c3_activ_pt->cd(1); 
    showHistUP(Activb_min, Activb_max, gh_mu_activ_ptc[0], "activity", TitleYhistEl, 0., 1., "log", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt11->Draw("same");
    //legel->Draw("same");
    c3_activ_pt->cd(1);          
    showHistDOWN(Activb_min, Activb_max,gh_mu_activ_ptc[0], "activity",TitleYSF ,0.,1., 0.7, 1.1, "log", 16., 2.4, 16, 0.005, elygraph[0]);

    c3_activ_pt->cd(2); 
    showHistUP(Activb_min, Activb_max, gh_mu_activ_ptc[1], "activity", TitleYhistEl, 0., 1., "log", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt22->Draw("same");
    //legel->Draw("same");
    c3_activ_pt->cd(2);          
    showHistDOWN(Activb_min, Activb_max,gh_mu_activ_ptc[1], "activity",TitleYSF ,0.,1., 0.8, 1.1, "log", 16., 2.4, 16, 0.005, elygraph[1]);

    c3_activ_pt->cd(3); 
    showHistUP(Activb_min, Activb_max, gh_mu_activ_ptc[2], "activity", TitleYhistEl, 0., 1., "log", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt33->Draw("same");
    //legel->Draw("same");
    c3_activ_pt->cd(3);          
    showHistDOWN(Activb_min, Activb_max,gh_mu_activ_ptc[2], "activity",TitleYSF ,0.,1., 0.9, 1.1, "log", 16., 2.4, 16, 0.005, elygraph[2]);


    c3_activ_pt->cd(4);  
    showHistUP(Activb_min, Activb_max, gh_mu_activ_ptc[3], "activity", TitleYhistEl, 0., 1., "log", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt44->Draw("same");
    //legel->Draw("same");
    c3_activ_pt->cd(4);          
    showHistDOWN(Activb_min, Activb_max,gh_mu_activ_ptc[3], "activity",TitleYSF ,0.,1., 0.9, 1.1, "log", 16., 2.4, 16, 0.005, elygraph[3]);

    c3_activ_pt->cd(5); 
    
    showHistUP(Activb_min, Activb_max, gh_mu_activ_ptc[4], "activity", TitleYhistEl, 0., 1., "log", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt55->Draw("same");
    //legel->Draw("same");
    c3_activ_pt->cd(5);          
    showHistDOWN(Activb_min, Activb_max, gh_mu_activ_ptc[4], "activity",TitleYSF ,0.,1., 0.9, 1.1, "log", 16., 2.4, 16, 0.005, elygraph[4]);

    
    c3_activ_pt->cd(6);
    
    TGraphErrors *grel = new TGraphErrors(5, graphpt, elygraph, graphpter, graphyer);
    grel->SetTitle("");
    grel->GetXaxis()->SetTitle("p_{T}, [GeV]");
    grel->GetYaxis()->SetTitle("slope");
    grel->Draw("ap");
    legel->Draw("same");
    c3_activ_pt->SaveAs("fastfull/ele_activ_pt_5range.pdf");
    */
    
    

 //__________________________________________________________________________________________________________________

    /*
    double muygraph[5] = {0.};

    TCanvas *c3_npv_pt_mu = new TCanvas("c3_npv_pt_mu", "mu # PV # pt", 0, 0, 1000, 700);
    c3_npv_pt_mu->Divide(3,2);
    c3_npv_pt_mu->cd(1);   
    showHistUP(npv_min, npv_max, gh_mu_npv_ptc[0], "# PV", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt11->Draw("same");
    //legmu->Draw("same");
    c3_npv_pt_mu->cd(1);          // Go back to the main canvas before defining pad2
    showHistDOWN(npv_min, npv_max,gh_mu_npv_ptc[0], "# PV",TitleYSF ,0.,1., 0.7, 1.1, "", 16., 2.4, 16, 0.005, muygraph[0]);

    c3_npv_pt_mu->cd(2); 
    showHistUP(npv_min, npv_max, gh_mu_npv_ptc[1], "# PV", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt22->Draw("same");
    //legmu->Draw("same");
    c3_npv_pt_mu->cd(2);          // Go back to the main canvas before defining pad2
    showHistDOWN(npv_min, npv_max,gh_mu_npv_ptc[1], "# PV",TitleYSF ,0.,1., 0.8, 1.1, "", 16., 2.4, 16, 0.005, muygraph[1]);

    c3_npv_pt_mu->cd(3); 
    showHistUP(npv_min, npv_max, gh_mu_npv_ptc[2], "# PV", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt33->Draw("same");
    //legmu->Draw("same");
    c3_npv_pt_mu->cd(3);          // Go back to the main canvas before defining pad2
    showHistDOWN(npv_min, npv_max,gh_mu_npv_ptc[2], "# PV",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, muygraph[2]);

    c3_npv_pt_mu->cd(4);
    showHistUP(npv_min, npv_max, gh_mu_npv_ptc[3], "# PV", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt44->Draw("same");
    //legmu->Draw("same");
    c3_npv_pt_mu->cd(4);          // Go back to the main canvas before defining pad2
    showHistDOWN(npv_min, npv_max,gh_mu_npv_ptc[3], "# PV",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, muygraph[3]);

    c3_npv_pt_mu->cd(5); 
    showHistUP(npv_min, npv_max, gh_mu_npv_ptc[4], "# PV", TitleYhistMu, 0., 1., "", 0.09, 0.7, 0.09, 0.005);
    xlabel_el_pt55->Draw("same");
    //legmu->Draw("same");
    c3_npv_pt_mu->cd(5);          // Go back to the main canvas before defining pad2
    showHistDOWN(npv_min, npv_max,gh_mu_npv_ptc[4], "# PV",TitleYSF ,0.,1., 0.9, 1.1, "", 16., 2.4, 16, 0.005, muygraph[4]);

    
    c3_npv_pt_mu->cd(6);
    TGraphErrors *grmu = new TGraphErrors(5, graphpt, muygraph, graphpter, graphyer);
    grmu->SetTitle("");
    grmu->GetXaxis()->SetTitle("p_{T}, [GeV]");
    grmu->GetYaxis()->SetTitle("slope");
    grmu->Draw("ap");
    legmu->Draw("same");
    c3_npv_pt_mu->SaveAs("fastfull/mu_npv_pt_5range.pdf");
    */

//______
    
    /*
    TLegend *leg = new TLegend(0.2, 0.15, 0.5, 0.45, "");
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->AddEntry(h_RhoAll[0], "Fastsim", "l");
    leg->AddEntry(h_RhoAll[1], "Fullsim", "l");


    TCanvas *cRhoAll1 = new TCanvas("cRhoAll1","cRhoAll1");
    h_RhoAll[0]->SetTitle("Rho All");
    h_RhoAll[0]->GetXaxis()->SetTitle("Rho All");
    h_RhoAll[0]->GetYaxis()->SetTitle("a.u.");
    h_RhoAll[0]->DrawNormalized();
    h_RhoAll[1]->SetLineColor(kBlue);
    h_RhoAll[1]->DrawNormalized("same");
    leg->Draw("same");

    TCanvas *cRhoCN1 = new TCanvas("cRhoCN1","cRhoCN1");
    h_RhoNC[0]->SetTitle("Rho Neutral Central");
    h_RhoNC[0]->GetXaxis()->SetTitle("Rho NC");
    h_RhoNC[0]->GetYaxis()->SetTitle("a.u.");
    h_RhoNC[0]->DrawNormalized();
    h_RhoNC[1]->SetLineColor(kBlue);
    h_RhoNC[1]->DrawNormalized("same");
    leg->Draw("same");
    */

    /*
    TCanvas *eleactivity = new TCanvas("eleactivity","eleactivity");
    h_ele_activ[0]->SetTitle("electron activity");
    h_ele_activ[0]->GetXaxis()->SetTitle("ele activity");
    h_ele_activ[0]->GetYaxis()->SetTitle("a.u.");
    h_ele_activ[0]->DrawNormalized("E");
    h_ele_activ[1]->DrawNormalized("Esame");
    */
   
//__________________________________________________________________________________________________________________

    
    /*TCanvas *c55 = new TCanvas("c55", "tau # PV Njets HT", 0, 0, 1000, 700);
    showHistUP(npv_min, npv_max, gh_tau_npvc, "# PV", TitleYhistTau, 0., 0.33);
    legtau->Draw("same");
    c55->cd();          // Go back to the main canvas before defining pad2
    showHistDOWN(numberBinsnpv+1, gh_tau_npvc, npv_min, npv_max,"# PV",TitleYSF ,0.,0.33, 0, 4.5);

    c55->cd();          // Go back to the main canvas before defining pad2
    showHistUP(HT_min, HT_max, gh_tau_HTc, "HT (GeV)", TitleYhistTau, 0.33, 0.66);
    legtau->Draw("same");
    c55->cd();          // Go back to the main canvas before defining pad2
    showHistDOWN(numberBinsHT+1, gh_tau_HTc, HT_min, HT_max,"HT (GeV)",TitleYSF ,0.33,0.66, 0, 4.5);

    c55->cd();          // Go back to the main canvas before defining pad2
    showHistUP(Njets_min, Njets_max, gh_tau_Njetsc, "Njets", TitleYhistTau, 0.66, 1.0);
    legtau->Draw("same");
    c55->cd();          // Go back to the main canvas before defining pad2
    showHistDOWN(numberBinsNjets, gh_tau_Njetsc, Njets_min, Njets_max, "Njets",TitleYSF ,0.66,1.0, 0, 4.5);

    //c55->SaveAs("fastfull/tau_npv_HT_njets.pdf");
    */

//__________________________________________________________________________________________________________________
    
    /*
    h_ele_ptetac[0]->Divide(h_ele_pteta[0]);
    h_ele_ptetac[1]->Divide(h_ele_pteta[1]);
    //h_ele_ptetac[1]->Divide(h_ele_ptetac[0]);
    //TCanvas *c6 = new TCanvas("c6", "2D electron", 0, 0, 1200, 900);
    //h_ele_ptetac[1]->Draw("colz");
    //h_ele_ptetac[1]->Draw("text e same");
    //h_ele_ptetac[1]->SetMarkerColor(kBlack);
    Print2D(h_ele_ptetac[1], h_ele_ptetac[0], "pteta");
   
    //c6->SaveAs("fastfull/2del.pdf");
    
    std::cout << h_mu_ptetac[0]->GetBinContent(5, 3) << " " << h_mu_pteta[0]->GetBinContent(5, 3) << std::endl;
    std::cout << h_mu_ptetac[1]->GetBinContent(5, 3) << " " << h_mu_pteta[1]->GetBinContent(5, 3) << std::endl;
    */
    /*
    h_tau_ptetac[0]->Divide(h_tau_pteta[0]);
    h_tau_ptetac[1]->Divide(h_tau_pteta[1]);
    h_tau_ptetac[1]->Divide(h_tau_ptetac[0]);
    TCanvas *c66 = new TCanvas("c66", "2D tau", 0, 0, 1200, 900);
    h_tau_ptetac[1]->Draw("colz");
    h_tau_ptetac[1]->Draw("text e same");
    h_tau_ptetac[1]->SetMarkerColor(kBlack);
    */
    //c66->SaveAs("fastfull/2dtau.pdf");

    /*
    h_mu_ptetac[0]->Divide(h_mu_pteta[0]);
    h_mu_ptetac[1]->Divide(h_mu_pteta[1]);
    h_mu_ptetac[1]->Divide(h_mu_ptetac[0]);
    TCanvas *c7 = new TCanvas("c7", "2D muon", 0, 0, 1200, 900);
    h_mu_ptetac[1]->Draw("colz");
    h_mu_ptetac[1]->Draw("text e same");
    h_mu_ptetac[1]->SetMarkerColor(kBlack);
    c7->SaveAs("fastfull/2dmu.pdf");
    
    std::cout << h_mu_ptetac[1]->GetBinContent(5, 3) << " " << h_mu_ptetac[1]->GetBinError(5,3)<< std::endl;\

    h_ele_ptactivc[0]->Divide(h_ele_ptactiv[0]);
    h_ele_ptactivc[1]->Divide(h_ele_ptactiv[1]);
    h_ele_ptactivc[1]->Divide(h_ele_ptactivc[0]);
    TCanvas *cptactiv_ele = new TCanvas("cptactiv_ele", "2D electron pt activ", 0, 0, 1200, 900);
    h_ele_ptactivc[1]->Draw("colz");
    h_ele_ptactivc[1]->Draw("text e same");
    h_ele_ptactivc[1]->SetMarkerColor(kBlack);
    cptactiv_ele->SaveAs("fastfull/2delptactiv.pdf");

    h_mu_ptactivc[0]->Divide(h_mu_ptactiv[0]);
    h_mu_ptactivc[1]->Divide(h_mu_ptactiv[1]);
    h_mu_ptactivc[1]->Divide(h_mu_ptactivc[0]);
    TCanvas *cptactiv_mu = new TCanvas("cptactiv_mu", "2D muon pt activ", 0, 0, 1200, 900);
    h_mu_ptactivc[1]->Draw("colz");
    h_mu_ptactivc[1]->Draw("text e same");
    h_mu_ptactivc[1]->SetMarkerColor(kBlack);
    cptactiv_mu->SaveAs("fastfull/2dmuptactiv.pdf");
    */

    
//__________________________________________________________________________________________________________________
   
    /*
    TCanvas *c8 = new TCanvas("c8", "Vertices number", 0, 0, 1000, 700);
    //h_ele_npv[1]->SetLineColor(kRed);
    //h_ele_npv[1]->SetMarkerColor(kRed);
    h_ele_npv[1]->SetLineColor(kBlue);
    h_ele_npv[1]->SetMarkerColor(kBlue);

    h_ele_npv[1]->SetTitle("Primary vertices");
    h_ele_npv[1]->GetXaxis()->SetTitle("# PV");
    h_ele_npv[1]->GetYaxis()->SetTitle("");
    h_ele_npv[1]->SetLineWidth(2);
    h_ele_npv[1]->SetFillStyle(3000);
    h_ele_npv[1]->DrawNormalized("hist"); 
    h_ele_npv[0]->SetLineWidth(2);
    h_ele_npv[0]->SetFillStyle(3000);
    h_ele_npv[0]->DrawNormalized("histsame");
      
    legel->Draw("same");
    */
  
  //__________________________________________________________________________________________________________________
  
    /*
    h_RhoAll[0]->Scale(1);
    h_RhoAll[2]->Scale(1);

    h_RhoNC[0]->Scale(1);
    h_RhoNC[2]->Scale(1);
    */

    /*
    TGraphAsymmErrors *gRhoAll = new TGraphAsymmErrors(h_RhoAll[2], h_RhoAll[0]);
    TGraphAsymmErrors *gRhoNC = new TGraphAsymmErrors(h_RhoNC[2], h_RhoNC[0]);
    */
    /*
    h_RhoAll[2]->Divide(h_RhoAll[0]);
    h_RhoNC[2]->Divide(h_RhoNC[0]);

    int rhoAllbins = h_RhoAll[2]->GetNbinsX();
    int rhoNCbins = h_RhoNC[2]->GetNbinsX();

    for(int i = 0; i < rhoAllbins; i++)
        std::cout << h_RhoAll[2]->GetBinContent(i) << ", " ;
    std::cout << std::endl;

    for(int i = 0; i < rhoNCbins; i++)
        std::cout << h_RhoNC[2]->GetBinContent(i) << ", " ;
    std::cout << std::endl;


    TCanvas *cRhoNC = new TCanvas("cRhoNC","cRhoNC");
    h_RhoNC[2]->Draw();

    TCanvas *cRhoAll = new TCanvas("cRhoAll","cRhoAll");
    h_RhoAll[2]->Draw();
    */
    
}

double weightCalc(double npv){

    return 7.28675e-1 - 6.34788e-3 * npv + 2.39709e-3 * npv * npv - 3.27811e-5 * npv * npv * npv;
}

void readTreeSync()
{
    tree();
}

void showHistUP(double min, double max, TGraphAsymmErrors *gh[2], string TitleX, string TitleY,double xposl,double xposr, string log, double titlesize, double titleoffset, double labelsize, double labeloffset){

    TPad *pad1 = new TPad("pad1", "pad1", xposl, 0.3, xposr, 1.0);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(log == "log")
        pad1->SetLogx();
    pad1->DrawFrame(min, -0.1, max, 1.1);
    TMultiGraph *multi = new TMultiGraph();

    multi->Add(gh[0]);
    multi->Add(gh[1]);
    multi->Draw("P");
    multi->SetTitle("");
    multi->GetXaxis()->SetTitle(TitleX.c_str());
    multi->GetYaxis()->SetTitle(TitleY.c_str());
    multi->GetXaxis()->SetRangeUser(min, max);
    multi->GetYaxis()->SetTitleSize(titlesize);
    multi->GetYaxis()->SetTitleOffset(titleoffset);
    multi->GetYaxis()->SetLabelSize(labelsize);
    multi->GetYaxis()->SetLabelOffset(labeloffset);

}

void showHistDOWN(double min,double max, TGraphAsymmErrors *gh[2],string TitleX, string TitleY, double xposl, double xposr, double miny, double maxy, string log, double titlesize, double titleoffset, double labelsize, double labeloffset, double &p1){

    vector<double> x;
    vector<double> y;
    vector<double> xerh;
    vector<double> xerl;
    vector<double> yerh;
    vector<double> yerl;

    vector<double> xSpr;
    vector<double> ySpr;
    vector<double> xerhSpr;
    vector<double> xerlSpr;
    vector<double> yerhSpr;
    vector<double> yerlSpr;

    int numPoints = gh[0]->GetN();

    for(int i = 0; i < numPoints;i++){

        double xx1, yy1, yyh1, yyl1;
        gh[0]->GetPoint(i, xx1, yy1);
        yyh1 = gh[0]->GetErrorYhigh(i);
        yyl1 = gh[0]->GetErrorYlow(i);
        
        int numPointsSpr = gh[1]->GetN();
        bool flagSpr = false;
        double xx3, yy3, yyh3, yyl3;
        for(int j = 0; j < numPointsSpr;j++){
            gh[1]->GetPoint(j, xx3, yy3);
            yyh3 = gh[1]->GetErrorYhigh(j);
            yyl3 = gh[1]->GetErrorYlow(j);
            if(xx3 == xx1){
                flagSpr = true;
                break;
            }
        }
        
        if(flagSpr){
            ySpr.push_back(yy3 / yy1);
            double a, b; 
            a = yy3/yy1 * TMath::Sqrt(TMath::Power(yyh1/yy1, 2) + TMath::Power(yyh3/yy3, 2));
            b = yy3/yy1 * TMath::Sqrt(TMath::Power(yyl1/yy1, 2) + TMath::Power(yyl3/yy3, 2));
            yerhSpr.push_back(TMath::Max(a,b));
            yerlSpr.push_back(TMath::Max(a,b));
            //std::cout << "xSpr: " << xx3 << " ySpr: " << yy3/yy1 << std::endl;

            xSpr.push_back(xx3);
            xerhSpr.push_back(gh[0]->GetErrorXhigh(i));
            xerlSpr.push_back(gh[0]->GetErrorXlow(i));
        }
    }

    TPad *pad = new TPad("pad", "pad", xposl, 0., xposr, 0.3);
    pad->SetTopMargin(0);
    pad->SetBottomMargin(0.3);
    pad->SetGridx(); // vertical grid
    if(log == "log")
        pad->SetLogx();
    pad->Draw();
    pad->cd();       // pad2 becomes the current pad
    pad->DrawFrame(min, 0.5, max, 1.1);

    TVectorT<double> tx(x.size(), &x[0]);
    TVectorT<double> ty(y.size(), &y[0]);
    TVectorT<double> txerl(xerl.size(), &xerl[0]);
    TVectorT<double> txerh(xerh.size(), &xerh[0]);
    TVectorT<double> tyerl(yerl.size(), &yerl[0]);
    TVectorT<double> tyerh(yerh.size(), &yerh[0]);

    TVectorT<double> txSpr(xSpr.size(), &xSpr[0]);
    TVectorT<double> tySpr(ySpr.size(), &ySpr[0]);
    TVectorT<double> txerlSpr(xerlSpr.size(), &xerlSpr[0]);
    TVectorT<double> txerhSpr(xerhSpr.size(), &xerhSpr[0]);
    TVectorT<double> tyerlSpr(yerlSpr.size(), &yerlSpr[0]);
    TVectorT<double> tyerhSpr(yerhSpr.size(), &yerhSpr[0]);

    TGraphAsymmErrors *gSpr = new TGraphAsymmErrors(txSpr, tySpr, txerlSpr, txerhSpr, tyerlSpr, tyerhSpr);

    TH1F *hgr = pad->DrawFrame(min, miny, max, maxy);

    TLine *line = new TLine(min, 1., max, 1.);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);

    gSpr->SetMarkerColor(kBlue);
    gSpr->SetLineColor(kBlue);
    gSpr->SetMarkerStyle(21);
    line->Draw("");
    //g->Draw("Psame");       // Draw the ratio plot
    gSpr->Draw("Psame");       // Draw the ratio plot
   
    TF1 *myfit = new TF1("myfit","pol1", 0, 30);
    if(p1 != 1.){
        myfit->SetParameter(0, 1.);
        myfit->SetParameter(1, -0.001);
        gSpr->Fit("myfit");
        p1 = myfit->GetParameter(1);
    }


    // Y axis ratio plot settings
    hgr->GetYaxis()->SetTitle(TitleY.c_str());
    hgr->GetXaxis()->SetTitle(TitleX.c_str());
    hgr->GetYaxis()->SetNdivisions(505);
    hgr->GetYaxis()->SetTitleSize(titlesize);
    hgr->GetYaxis()->SetTitleFont(43);
    hgr->GetYaxis()->SetTitleOffset(titleoffset);
    hgr->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hgr->GetYaxis()->SetLabelSize(labelsize);
    hgr->GetYaxis()->SetLabelOffset(labeloffset);
    // X axis ratio plot settings
    hgr->GetXaxis()->SetTitleSize(22);
    hgr->GetXaxis()->SetTitleFont(43);
    hgr->GetXaxis()->SetTitleOffset(4.);
    hgr->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hgr->GetXaxis()->SetLabelSize(22);
  
}


TH2D TranslateHisto(const TH2 &input){
  int nx = input.GetNbinsX();
  int ny = input.GetNbinsY();
  TH2D output(input.GetName(), input.GetTitle(), nx, 0.5, nx+0.5, ny, 0.5, ny+0.5);
  output.Sumw2();
  output.SetStats(false);
  output.SetMarkerSize(2);
  output.SetLabelSize(0.05,"XYZ");
  output.SetTitleSize(0.05,"XYZ");
  
  //std::cout << input.GetName() << " ";

  for(int ix = 0; ix <= nx+1; ++ix){
    //double minActiv = 999.;
    //double maxActiv = -999.;
    for(int iy = 0; iy <= ny+1; ++iy){
      /*
      if(input.GetBinContent(ix, iy) > maxActiv)
        maxActiv = input.GetBinContent(ix, iy);
      if(input.GetBinContent(ix, iy) < minActiv)
        minActiv = input.GetBinContent(ix, iy);
        */
      output.SetBinContent(ix, iy, input.GetBinContent(ix, iy));
      output.SetBinError(ix, iy, input.GetBinError(ix, iy));
    }
    //std::cout << (maxActiv - minActiv) / 2 << " ";
  }

  std::cout << std::endl;
  
  for(int ix = 1; ix <= nx; ++ix){
    const TAxis *iaxis = input.GetXaxis();
    TAxis *oaxis = output.GetXaxis();
    if(iaxis == NULL || oaxis == NULL) continue;
    oaxis->SetTitle(iaxis->GetTitle());
    ostringstream oss;
    oss << iaxis->GetBinLowEdge(ix) << "-" << iaxis->GetBinUpEdge(ix) << flush;
    oaxis->SetBinLabel(ix, oss.str().c_str());
  }

  for(int iy = 1; iy <= ny; ++iy){
    const TAxis *iaxis = input.GetYaxis();
    TAxis *oaxis = output.GetYaxis();
    if(iaxis == NULL || oaxis == NULL) continue;
    oaxis->SetTitle(iaxis->GetTitle());
    ostringstream oss;
    oss << iaxis->GetBinLowEdge(iy) << "-" << iaxis->GetBinUpEdge(iy) << flush;
    oaxis->SetBinLabel(iy, oss.str().c_str());
  }

  return output;
}

void Print2D(TH2 const * const h_data_in, TH2 const * const h_mc_in, const TString &ext){
  if(h_data_in == NULL || h_mc_in == NULL) return;

  TH2D h_data = TranslateHisto(*h_data_in);
  TH2D h_mc = TranslateHisto(*h_mc_in);

  TCanvas canvas;
  gStyle->SetPalette(bands, rainbow);
  h_data.SetMarkerSize(2);
  h_mc.SetMarkerSize(2);

  h_data.SetTitle("Fullsim Eff");
  h_data.GetZaxis()->SetRangeUser(0.,1.);
  h_data.Draw("colz");
  h_data.Draw("textesame");
  canvas.Print("FinalResult/plots/2d_full_"+ext+".pdf");
  //canvas.Print("plots/2d_data_"+ext+".png");
  //PrintTable(&h_data, "data_"+ext);
  h_mc.SetTitle("Fastsim Eff");
  h_mc.GetZaxis()->SetRangeUser(0.,1.);
  h_mc.Draw("colz");
  h_mc.Draw("textesame");
  canvas.Print("FinalResult/plots/2d_fast_"+ext+".pdf");
  //canvas.Print("plots/2d_mc_"+ext+".png");
  //PrintTable(&h_mc, "mc_"+ext);
  
  gStyle->SetPalette(bands, patriotic);
  h_data.Divide(&h_mc);
/*
  int nx = h_data.GetNbinsX();
  int ny = h_data.GetNbinsY();

  std::cout << h_data.GetName() << " ";
  
  for(int ix = 1; ix <= nx; ++ix){
    double minActiv = 999.;
    double maxActiv = -999.;
    for(int iy = 1; iy <= ny; ++iy){
      
      if(h_data.GetBinContent(ix, iy) > maxActiv)
        maxActiv = h_data.GetBinContent(ix, iy);
      if(h_data.GetBinContent(ix, iy) < minActiv)
        minActiv = h_data.GetBinContent(ix, iy);
        

      //std::cout << h_data.GetBinContent(ix, iy) << " "; 
    }
    std::cout << fixed << setprecision(3) << (maxActiv - minActiv) / 2 << " ";
    //std::cout<<std::endl;
  }

  std::cout << std::endl;
  */

  h_data.SetMinimum(0.8);
  h_data.SetMaximum(1.25);
  canvas.SetLogz();
  h_data.GetZaxis()->SetRangeUser(0.5,2);
  h_data.SetTitle("Scale Factors");
  h_data.Draw("colz");
  h_data.Draw("textesame");
  canvas.Print("FinalResult/plots/sf_"+ext+".pdf");
  

  //h_data_in->Divide(h_mc_in);

  
  //TH2D *hnew = (TH2D*)h_data_in->Clone();
  //hnew->SetName("histo");

  //hnew->Divide(h_mc_in);


  //TFile *f = TFile::Open("files/sf_"+ext+".root","recreate");

  //hnew->Write();

  //f->Close();


  //canvas.Print("plots/sf_"+ext+".png");
  //PrintTable(&h_data, "sf_"+ext);
}

void Print3D(TH3 const * const h_full_in, TH3 const * const h_fast_in, const TString &ext){
  if(h_full_in == NULL || h_fast_in == NULL) return;

  TCanvas canvas;

  //h_full_in->Divide(&h_fast_in);
  //h_data.SetMinimum(0.8);
  //h_data.SetMaximum(1.25);
  //canvas.SetLogz();
  //h_full_in->SetTitle("Scale Factors");
  //h_full_in->Draw();
  //canvas.Print("plots/3Dsf_"+ext+".pdf");

  //h_data_in->Divide(h_mc_in);

  TH3D *hnew = (TH3D*)h_full_in->Clone();
  hnew->SetName("histo");

  hnew->Divide(h_fast_in);

  //TFile *f = TFile::Open("files/3Dsf_"+ext+".root","recreate");

  hnew->Write();

  //f->Close();

  //canvas.Print("plots/sf_"+ext+".png");
  //PrintTable(&h_data, "sf_"+ext);
}

void PrintTable(TH2 const * const histo, const TString &ext){
  int eta = ext.Index("abseta");
  int pt = ext.Index("_et_");
  if(eta<0 || pt<0 || pt<eta) return;

  ofstream file("tables/"+ext+".tex");
  file << "\\documentclass{article}\n\n";

  file << "\\begin{document}\n";
  file << "\\begin{table}\n";
  file << "  \\begin{tabular}{r|rrrrr}\n";
  file << "    \\hline\\hline\n";
  file << "    & \\multicolumn{5}{c}{$p_T$ [GeV]}\\\\\n";
  file << "    $|\\eta|$ & 10-20 & 20-30 & 30-40 & 40-50 & 50-200\\\\\n";
  file << "    \\hline\n";
  PrintLine(file, histo, 1, "0-1.442");
  PrintLine(file, histo, 2, "1.442-1.556");
  PrintLine(file, histo, 3, "1.556-2.5");
  file << "    \\hline\\hline\n";
  file << "  \\end{tabular}\n";
  file << "\\end{table}\n";
  file << "\\end{document}\n";

  file.flush();
  file.close();
}

void PrintLine(ofstream &file, TH2 const * const histo, int bin, const TString &label){
  if(!file.is_open()) return;

  if(histo == NULL || histo->GetNbinsX() < bin || histo->GetNbinsY()<5) return;
  file << "    " << label;
  for(int y = 1; y <= 5; ++y){
    file << " & $"
         << fixed << setprecision(3) << histo->GetBinContent(bin, y)
         << "\\pm"
         << fixed << setprecision(3) << histo->GetBinError(bin, y)
         << "$";
  }
  file << "\\\\\n";
}

void GetPatrioticPalette(){
  const unsigned num = 3;   
  double red[num] = {0., 1., 1.};
  double green[num] = {0., 1., 0.};
  double blue[num] = {1., 1., 0.};
  double stops[num] = {0., 0.5, 1.};
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    patriotic[i] = fi + i;
  }
}

void GetRainbowPalette(){
  const unsigned num = 6;
  double red[num] =   {1.,0.,0.,0.,1.,1.};
  double green[num] = {0.,0.,1.,1.,1.,0.};
  double blue[num] =  {1.,1.,1.,0.,0.,0.};
  double stops[num] = {0.,0.2,0.4,0.6,0.8,1.};
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    rainbow[i] = fi+i;
  }
}

//void FillHisto(TH2D* hpteta_2dnpv[4], TH3D* hptetanpv, TH2D* hpteta, TH2D* hptactiv, TH2D* hactivHT, TH2D* hactivnpv, TH2D* hHTnpv, TH1D *hist[6], double pt, double eta, double ht, double njets, double npv, double activity){
//void FillHisto(TH2D* hpteta_2dnpv[4], TH3D* hptetanpv, TH2D* hpteta, TH1D *hist[6], double pt, double eta, double ht, double njets, double npv, double activity){
void FillHisto(TH2D* hpteta_2dnpv[4], TH3D* hptetanpv, TH2D* hpteta, TH2D* hptactiv, TH2D* hetaactiv, TH1D *hist[6], double pt, double eta, double ht, double njets, double npv, double activity){

    if(npv > 0 && npv < 10)
        hpteta_2dnpv[0]->Fill(pt, TMath::Abs(eta));
    if(npv > 9 && npv < 14)
        hpteta_2dnpv[1]->Fill(pt, TMath::Abs(eta));
    if(npv > 13 && npv < 18)
        hpteta_2dnpv[2]->Fill(pt, TMath::Abs(eta));
    if(npv > 17 && npv < 41)
        hpteta_2dnpv[3]->Fill(pt, TMath::Abs(eta)); 
    hptetanpv->Fill(pt, TMath::Abs(eta), npv);
    //if(ht > 200 && ht < 300)
    //if(activity > 0.15 && activity < 1)
    hpteta->Fill(pt, TMath::Abs(eta));
    hptactiv->Fill(pt, TMath::Min(activity, 9.9));
    //if(pt > 30 && pt < 40)
    //    hetaactiv->Fill(TMath::Abs(eta), TMath::Min(activity, 9.9));
    //hactivHT->Fill(TMath::Min(activity, 9.9), ht);
    //hactivnpv->Fill(TMath::Min(activity, 9.9), npv);
    //hHTnpv->Fill(ht, npv);
    hist[0]->Fill(pt);
    hist[1]->Fill(eta);
    hist[4]->Fill(npv);
    hist[2]->Fill(ht);
    hist[3]->Fill(njets);
    hist[5]->Fill(TMath::Min(activity, 9.9));
                           
}

void Print1D(TGraphAsymmErrors *gh[6][2], const TString &ext, const string &TitleY, const string &lepton, TLegend *leg){

    TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 1000, 700);
    
    showHistUP(pt_min, pt_max, gh[0], lepton + " p_{T} (GeV)", TitleY, 0., 1.0, "", 0.06, 1.2, 0.05, 0.005);
    
    TString folder = "FinalResult/plots/";

    leg->Draw("same");
    c1->cd();         
    showHistDOWN(pt_min, pt_max, gh[0], lepton + " p_{T} (GeV)","Fullsim/Fastsim" ,0.,1.0, 0.7, 1.2, "", 26, 1.85, 22., 0.005);
    c1->SaveAs(folder + "pt_" + ext + ".pdf");

    //c1->cd();
    TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 1000, 700);
    
    showHistUP(eta_min, eta_max, gh[1], lepton + " #eta", TitleY, 0., 1.0, "", 0.06, 1.2, 0.05, 0.005);
    leg->Draw("same");
    c2->cd();         
    showHistDOWN(eta_min, eta_max,gh[1], lepton + " #eta", "Fullsim/Fastsim",0.,1.0, 0.7, 1.2, "", 26, 1.85, 22., 0.005);
    
    c2->SaveAs(folder + "eta_" + ext + ".pdf");


    TCanvas *c3 = new TCanvas("c3", "c3", 0, 0, 1000, 700);
    //c2->Divide(2,2);
    //c3->cd(1);
    showHistUP(npv_min, npv_max, gh[4], "# PV", TitleY, 0., 1., "", 0.06, 1.2, 0.05, 0.005);
    leg->Draw("same");
    c3->cd();          
    showHistDOWN(npv_min, npv_max,gh[4], "# PV", "Fullsim/Fastsim",0.,1., 0.8, 1.2, "", 26, 1.85, 22., 0.005);
    c3->SaveAs(folder + "npv_" + ext + ".pdf");

    TCanvas *c4 = new TCanvas("c4", "c4", 0, 0, 1000, 700);
    //c2->cd(2);
    showHistUP(Activb_min, Activb_max, gh[5], "Activity", TitleY, 0, 1., "log", 0.06, 1.2, 0.05, 0.005);
    leg->Draw("same");
    c4->cd();          
    showHistDOWN(Activb_min, Activb_max,gh[5], "Activity","Fullsim/Fastsim" ,0.,1., 0.8, 1.2, "log", 26, 1.85, 22., 0.005);
    c4->SaveAs(folder + "activ_" + ext + ".pdf");
    

//---------
    
    TCanvas *c5 = new TCanvas("c5", "c5", 0, 0, 1000, 700);
    //c2->cd(3);          
    showHistUP(HT_min, HT_max, gh[2], "HT (GeV)", TitleY, 0., 1., "", 0.06, 1.2, 0.05, 0.005);
    leg->Draw("same");
    c5->cd();          // Go back to the main canvas before defining pad2
    showHistDOWN(HT_min, HT_max,gh[2], "HT (GeV)","Fullsim/Fastsim" ,0.,1., 0.8, 1.2, "", 26, 1.85, 22., 0.005);
    c5->SaveAs(folder + "HT_" + ext + ".pdf");

    //c2->cd(4);   
    TCanvas *c6 = new TCanvas("c6", "c6", 0, 0, 1000, 700);       
    showHistUP(Njets_min, Njets_max, gh[3], "Njets", TitleY, 0., 1., "", 0.06, 1.2, 0.05, 0.005);
    leg->Draw("same");
    c6->cd();         
    showHistDOWN(Njets_min, Njets_max, gh[3], "Njets","Fullsim/Fastsim" ,0.,1., 0.8, 1.2, "", 26, 1.85, 22., 0.005);
    
    c6->SaveAs(folder + "njets_" + ext + ".pdf");
    

}


/*
                            if(_lPt[num] > 10 && _lPt[num] < 20)
                                h_mu_npv_pt[0][sam]->Fill(_n_PV); 
                            if(_lPt[num] > 20 && _lPt[num] < 30)
                                h_mu_npv_pt[1][sam]->Fill(_n_PV); 
                            if(_lPt[num] > 30 && _lPt[num] < 40)
                                h_mu_npv_pt[2][sam]->Fill(_n_PV); 
                            if(_lPt[num] > 40 && _lPt[num] < 50)
                                h_mu_npv_pt[3][sam]->Fill(_n_PV); 
                            if(_lPt[num] > 50)
                                h_mu_npv_pt[4][sam]->Fill(_n_PV);
                                */
                            
                            /*
                            if(_lPt[num] > 10 && _lPt[num] < 20 && _lEta[num] > -2.4 && _lEta[num] < -2){
                                total_muons_number_1++;
                                h_ptratio_bad_mu[sam]->Fill(_lPt[num]/_closeJetPtAll[num]);
                                if(_lPt[num]/_closeJetPtAll[num] < 0.76){
                                    h_ptrel_bad_mu[sam]->Fill(_ptrel[num]);
                                    muons_number_HT200_1++;
                                }
                            }
                            if(_lPt[num] > 100 && _lEta[num] > 0. && _lEta[num] < 0.4){
                                total_muons_number_2++;
                                h_ptratio_good_mu[sam]->Fill(_lPt[num]/_closeJetPtAll[num]);
                                if(_lPt[num]/_closeJetPtAll[num] < 0.76){
                                    h_ptrel_good_mu[sam]->Fill(_ptrel[num]);
                                    muons_number_HT200_2++;
                                }
                            }
                            */
                                /*
                                if(_lPt[num] > 10 && _lPt[num] < 20)
                                    h_mu_npv_ptc[0][sam]->Fill(_n_PV); 
                                if(_lPt[num] > 20 && _lPt[num] < 30)
                                    h_mu_npv_ptc[1][sam]->Fill(_n_PV); 
                                if(_lPt[num] > 30 && _lPt[num] < 40)
                                    h_mu_npv_ptc[2][sam]->Fill(_n_PV); 
                                if(_lPt[num] > 40 && _lPt[num] < 50)
                                    h_mu_npv_ptc[3][sam]->Fill(_n_PV); 
                                if(_lPt[num] > 50)
                                    h_mu_npv_ptc[4][sam]->Fill(_n_PV);
                                */