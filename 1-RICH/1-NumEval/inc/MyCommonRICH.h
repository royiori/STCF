#ifndef _MyCommonRICH_h_
#define _MyCommonRICH_h_

#include <map>
#include <string>
#include "TF1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "MyRICHDetector.h"
#include "MyStyle.h"

using namespace std;

const int NHYPO = 4;
const int NMoment = 10;
const int NThetaC = 10;

const Double_t gkMassProton = 0.93827203;      // GeV
const Double_t gkMassNeutron = 0.93956536;     // GeV
const Double_t gkMassElectron = 0.00051099892; // GeV
const Double_t gkMassPion = 0.13957061;        // GeV
const Double_t gkMassKaon = 0.493677;          // GeV
const Double_t gkMassMuon = 0.1134289257;      // GeV
const Double_t pi = TMath::Pi();

class MyDatabaseClass;

class MyCommonRICH
{
public:
    MyCommonRICH();
    ~MyCommonRICH();

    //----------------------------
    /// Mass & Beta function
    int GetHypoID(TString p);
    int GetMomID(double mom);
    int GetThetaID(double the);
    double GetMass(TString p);
    double Beta(Double_t p, Double_t m); //p in GeV, m in GeV
    double BetaByTheta(Double_t theta, Double_t n);
    double Momentum(Double_t beta, Double_t m);
    double Momentum(Double_t theta, Double_t n, Double_t m);

    /// photon function [nm]
    double PhotonEng(Double_t lambda); //nm
    double PhotonWavelength(Double_t energy);
    double CherenkovCosAng(Double_t beta, Double_t n);
    double CherenkovCosAng(Double_t p, Double_t m, Double_t n);
    double CherenkovAng(Double_t beta, Double_t n);
    double CherenkovAng(Double_t p, Double_t m, Double_t n);

    /// Full reflective angle
    double FullReflectiveSinAng(Double_t n1, Double_t n2); //from n1 to n2, (n1>n2)
    double FullReflectiveCosAng(Double_t n1, Double_t n2); //from n1 to n2
    double FullReflectiveAng(Double_t n1, Double_t n2);    //from n1 to n2
    double CherenkovMaxMoment(Double_t n1, Double_t n2, Double_t m);

    double dNdLdLambda(Double_t lambda, Double_t beta, Double_t n);

    //----------------------------
    /// get detector
    MyRICHDetector *GetDetector() { return gDet; }
    MyRICHDetector *GetDetList(int ihypo) { return (0 <= ihypo && ihypo < (int)gDetList.size()) ? gDetList[ihypo] : NULL; }
    MyRICHDetector *GetDetScan(int imom, int itheta0, int ihypo)
    {
        if (0 <= imom && imom < (int)gScanDetList.size())
            if (0 <= itheta0 && itheta0 < (int)gScanDetList[imom].size())
                if (0 <= ihypo && ihypo < (int)gScanDetList[imom][itheta0].size())
                    return gScanDetList[imom][itheta0][ihypo];
        return NULL;
    }

    //----------------------------
    /// get detector hit map histograms
    TH1F *GetDetWLMap() { return gDet->GetWaveLengthHist(); }
    TH1F *GetDetListWLMap() { return gDet->GetWaveLengthHist(); }

    TH2F *GetDetHitMap() { return gDet->GetHitMap(); }
    TH2F *GetDetHitMap(int irad) { return gDet->GetDetHitMap(irad); }
    TH2F *GetDetListHitMap(int ihypo) { return gDetList[ihypo]->GetHitMap(); }
    TH2F *GetDetScanHitMap(int imom, int itheta0, int ihypo) { return (GetDetScan(imom, itheta0, ihypo) == NULL) ? NULL : GetDetScan(imom, itheta0, ihypo)->GetHitMap(); }
    TH2F *GetDetScanNPhMap(TString shypo) { return fNPhMap[GetHypoID(shypo)]; }

    TH1D *GetDetScanNphMapAtMom(TString shypo, double mom) { return MyProjectY(fNPhMap[GetHypoID(shypo)], mom); }
    TH1D *GetDetScanNphMapAtTheta(TString shypo, double the) { return MyProjectX(fNPhMap[GetHypoID(shypo)], the); }
    TH1D *GetDetScanNphMapAtMom(int ihypo, double mom) { return MyProjectY(fNPhMap[ihypo], mom); }
    TH1D *GetDetScanNphMapAtTheta(int ihypo, double the) { return MyProjectX(fNPhMap[ihypo], the); }
    TH1D *GetDetScanNphEachRadMapAtMom(TString shypo, double mom, int irad) { return MyProjectY(fNPhMapEachRad[GetHypoID(shypo)][irad], mom); }
    TH1D *GetDetScanNphEachRadMapAtTheta(TString shypo, double the, int irad) { return MyProjectX(fNPhMapEachRad[GetHypoID(shypo)][irad], the); }

    TH1D *MyProjectX(TH2F *f, double yval)
    {
        int ybin = 0;
        for (int i = 0; i <= f->GetYaxis()->GetNbins(); i++)
            if (f->GetYaxis()->GetBinLowEdge(i) <= yval && yval < f->GetYaxis()->GetBinLowEdge(i + 1))
                ybin = i;

        if (yval == f->GetYaxis()->GetBinLowEdge(f->GetYaxis()->GetNbins() + 1)) //maximum
            ybin = f->GetYaxis()->GetNbins();

        return f->ProjectionX("", ybin, ybin);
    }

    TH1D *MyProjectY(TH2F *f, double xval)
    {
        int xbin = 0;
        for (int i = 0; i <= f->GetXaxis()->GetNbins() + 1; i++)
            if (f->GetXaxis()->GetBinLowEdge(i) <= xval && xval < f->GetXaxis()->GetBinLowEdge(i + 1))
                xbin = i;

        if (xval == f->GetXaxis()->GetBinLowEdge(f->GetXaxis()->GetNbins() + 1))
            xbin = f->GetXaxis()->GetNbins();

        return f->ProjectionY("", xbin, xbin);
    }

    //----------------------------
    /// get detector reconstruction histograms
    TH2F *GetDetRecOffsetMap() { return fOffsetMap; }
    TH1F *GetDetRecOffsetVsNphPlot(int ihypo) { return fOffsetVsNphPlot[ihypo]; }
    TH1F *GetDetRecOffsetVsMomPlot(int ihypo) { return fOffsetVsMomPlot[ihypo]; }
    TH1F *GetDetRecOffsetVsThetaPlot(int ihypo) { return fOffsetVsThetaPlot[ihypo]; }

    TH2F *GetDetRecSigmaMap() { return fSigmaMap; }
    TH1F *GetDetRecSigmaVsNphPlot(int ihypo) { return fSigmaVsNphPlot[ihypo]; }
    TH1F *GetDetRecSigmaVsMomPlot(int ihypo) { return fSigmaVsMomPlot[ihypo]; }
    TH1F *GetDetRecSigmaVsThetaPlot(int ihypo) { return fSigmaVsThetaPlot[ihypo]; }

    TH2F *GetDetRecRing(int imom, int itheta0, int ihypo) { return (GetDetScan(imom, itheta0, ihypo) == NULL) ? NULL : GetDetScan(imom, itheta0, ihypo)->GetRecMap(); }
    TH1D *GetRecMap(TString shypo, int irad, int imom, int itheta0, int iph) { return GenRecMap(GetHypoID(shypo), irad, imom, itheta0, iph); }
    TH1D *GetRecMap(int ihypo, int irad, int imom, int itheta0, int iph) { return GenRecMap(ihypo, irad, imom, itheta0, iph); }

    int GetDetListNumber() { return gDetList.size(); }     //nhypo
    int GetDetScanNumber() { return gScanDetList.size(); } //[动量][角度][粒子种类]
    int GetDetRectNumber() { return fRecMap.size(); }      //[粒子种类][辐射体][动量][角度][光子数]

    //----------------------------
    /// get detector pid map histograms
    TH2F *GetPIDMap(int ihypo) { return fPIDMap[ihypo]; }
    TH1F *GetPIDMapVsMom(int ihypo) { return fPIDVsMomPlot[ihypo];}
    TH1F *GetPIDMapVsTheta(int ihypo) { return fPIDVsThetaPlot[ihypo];}

    //----------------------------
    // draw histograms
    void DrawDetHitMap(TString opt)
    {
        if (gDet->GetHitMap() != NULL)
            gDet->GetHitMap()->Draw(opt);
    }
    void DrawDetHitMap(int irad, TString opt)
    {
        if (gDet->GetDetHitMap(irad) != NULL)
            gDet->GetDetHitMap(irad)->Draw(opt);
    }

    void DrawDetListHitMap(int ihypo, TString opt)
    {
        if (gDetList[ihypo]->GetHitMap() != NULL)
            gDetList[ihypo]->GetHitMap()->Draw(opt);
    }
    void DrawDetListHitMap(int ihypo, int irad, TString opt)
    {
        if (gDetList[ihypo]->GetHitMap(irad) != NULL)
            gDetList[ihypo]->GetHitMap(irad)->Draw(opt);
    }

    void DrawScanDetListHitMap(int imom, int ithe, int ihypo, TString opt)
    {
        MyRICHDetector *det = GetDetScan(imom, ithe, ihypo);
        if (det == NULL)
            return;
        if (det->GetHitMap() == NULL)
            return;
        det->GetHitMap()->Draw(opt);
    }
    void DrawScanDetListHitMap(int imom, int ithe, int ihypo, int irad, TString opt)
    {
        MyRICHDetector *det = GetDetScan(imom, ithe, ihypo);
        if (det == NULL)
            return;
        if (det->fHitMapEachRad[irad] == NULL)
            return;
        det->fHitMapEachRad[irad]->Draw(opt);
    }
    void DrawScanDetListHitMap(TString shypo, int imom, int ithe, int irad, TString opt) { DrawScanDetListHitMap(imom, ithe, GetHypoID(shypo), irad, opt); }

    //----------------------------
    // set functions
    void SetNEvent(int ne) { nEvent = ne; }
    void SetPrecision(double ep) { epsilon = ep; }
    void SetDatabase(MyDatabaseClass *d) { gDb = d; }

    void SaveRings(const char *fname);
    void LoadRings(const char *fname);

    void SaveRecFile(const char *fname);
    void LoadRecFile(const char *fname);

    //----------------------------
    // calculations
    void GenerateDetRing(MyRICHDetector *det = 0);
    void GenerateMultiParticleRICHRings();
    void GenerateTheScanHitMapsForEachDetector();
    void GenerateTheNPhotonMap();

    bool ReconstructForEachDetector();
    void GenerateRecOffsetSigmaMap();
    void GenerateRecHistograms(TString particle, int irad, int imom, int ithe, int iph);
    void GeneratePIDHistograms(TString particle, int imom, int ithe);

    // reconstruction
    void ReconstructRICHBySolver(MyRICHDetector *det = 0);
    double ReconstructRICHBySolver(MyRICHDetector *det, double irad, double Xc, double Yc);
    double ReconstructRICHByBeta(MyRICHDetector *det, double irad, double Xc, double Yc);
    double ReconstructRICHByBetaSolver(MyRICHDetector *det, double irad, double Xc, double Yc);
    double Findz0(MyRICHDetector *det, int irad, double Xc, double Yc);
    double FindPhi(MyRICHDetector *det, double Xc, double Yc, double X0, double Y0);

    // pid efficiency 
    bool CalPIDEfficiency();
    int PIDProb(vector<MyRICHDetector *> detlist, vector<pair<double, double>> hit);
    double CalPIDProb(MyRICHDetector *det, vector<pair<double, double>> hit);


    void DrawFCN(int irad);

private:
    double ProjectToPixel(MyRICHDetector *, double xmin, double xmax, double ymin, double ymax, double lambda);
    double ProjectToPixel(MyRICHDetector *, double xmin, double xmax, double ymin, double ymax);

    double PhotonGenFromRad(MyRICHDetector *det, int irad, int AbsFlag = 1);
    double PhotonGenFromRadWithAbs(MyRICHDetector *det, int irad) { return PhotonGenFromRad(det, irad, 1); }
    double PhotonGenFromRadWithOutAbs(MyRICHDetector *det, int irad) { return PhotonGenFromRad(det, irad, 0); }

    //----------------------------
    // support function
    vector<double> GetThickList(MyRICHDetector *, int);
    vector<double> GetRefIndList(MyRICHDetector *, int, double);
    vector<double> GetAbsLenList(MyRICHDetector *, int, double);

    void UpdateThetaCExp(MyRICHDetector *det, double hypolambda);

    // variables
    int nEvent;
    double epsilon;

    TF1 *fRhoFcn;
    TF1 *fNphFcn;
    MyDatabaseClass *gDb;                                  //
    MyRICHDetector *gDet;                                  //id=0, 用来保存/读取探测器设置
    vector<MyRICHDetector *> gDetList;                     //id=1~NHYPO,[粒子种类]
    vector<vector<vector<MyRICHDetector *>>> gScanDetList; //id=10.....,[动量][角度][粒子种类]

    // histograms to show
    vector<TH2F *> fNPhMap;                //X为动量SCAN范围，Y为入射角度SCAN范围，Z为光子数, [粒子种类]
    vector<vector<TH2F *>> fNPhMapEachRad; //[粒子种类][辐射体]

    vector<vector<vector<vector<vector<TH1D *>>>>> fRecMap;        //为重建的角度填图， 用高斯拟合得offset和sigma, [粒子种类][辐射体][动量][角度][光子数]
    vector<vector<vector<vector<vector<double>>>>> fRecOffList;    //为重建的角度填图， 用高斯拟合得offset和sigma, [粒子种类][辐射体][动量][角度][光子数]
    vector<vector<vector<vector<vector<double>>>>> fRecSigList;    //为重建的角度填图， 用高斯拟合得offset和sigma, [粒子种类][辐射体][动量][角度][光子数]
    vector<vector<vector<vector<vector<double>>>>> fRecOffErrList; //为重建的角度填图， 用高斯拟合得offset和sigma, [粒子种类][辐射体][动量][角度][光子数]
    vector<vector<vector<vector<vector<double>>>>> fRecSigErrList; //为重建的角度填图， 用高斯拟合得offset和sigma, [粒子种类][辐射体][动量][角度][光子数]

    vector<vector<vector<vector<double>>>> fPidEffList; //[imom][ithe][ihypo][pidhypo]

    TH2F *fOffsetMap, *fSigmaMap;
    vector<TH1F *> fOffsetVsNphPlot;   //[hypo]
    vector<TH1F *> fOffsetVsMomPlot;   //[hypo]
    vector<TH1F *> fOffsetVsThetaPlot; //[hypo]
    vector<TH1F *> fSigmaVsNphPlot;
    vector<TH1F *> fSigmaVsMomPlot;
    vector<TH1F *> fSigmaVsThetaPlot;

    vector<TH2F *> fPIDMap;
    vector<TH1F *> fPIDVsMomPlot;
    vector<TH1F *> fPIDVsThetaPlot;

    void ResizeDetList(int n)
    {
        for (int i = 0; i < (int)gDetList.size(); i++)
            if (gDetList[i] != 0)
                delete gDetList[i];

        gDetList.clear();
        gDetList.resize(n);
    }

    void ResizeScanDetList(int i, int j, int k) //[动量][角度][粒子种类]
    {
        for (int i = 0; i < (int)gScanDetList.size(); i++)
            for (int j = 0; j < (int)gScanDetList[i].size(); j++)
                for (int k = 0; k < (int)gScanDetList[i][j].size(); k++)
                {
                    if (gScanDetList[i][j][k] != 0)
                        delete gScanDetList[i][j][k];
                }

        gScanDetList.resize(i);
        for (int imom = 0; imom < i; imom++)
        {
            gScanDetList[imom].resize(j);
            for (int ithe = 0; ithe < j; ithe++)
                gScanDetList[imom][ithe].resize(k);
        }
    }

    void ResizeNPhMap(int n, int r) //[粒子种类][辐射体数]
    {
        for (int i = 0; i < (int)fNPhMap.size(); i++)
            if (fNPhMap[i] != 0)
                delete fNPhMap[i];
        for (int i = 0; i < (int)fNPhMapEachRad.size(); i++)
            for (int j = 0; j < (int)fNPhMapEachRad[i].size(); j++)
                delete fNPhMapEachRad[i][j];

        fNPhMap.resize(n);
        fNPhMapEachRad.resize(n);
        for (int i = 0; i < n; i++)
            fNPhMapEachRad[i].resize(r);
    }

    void ResizeRecMap(int nhypo, int rad, int mom, int the, int nph) //[粒子种类][辐射体][动量][角度]
    {
        fRecOffList.clear();
        fRecSigList.clear();
        fRecOffErrList.clear();
        fRecSigErrList.clear();
        for (int i = 0; i < (int)fRecMap.size(); i++)
            for (int j = 0; j < (int)fRecMap[i].size(); j++)
                for (int k = 0; k < (int)fRecMap[i][j].size(); k++)
                    for (int l = 0; l < (int)fRecMap[i][j][k].size(); l++)
                        for (int m = 0; m < (int)fRecMap[i][j][k][l].size(); m++)
                            if (fRecMap[i][j][k][l][m] != NULL)
                                delete fRecMap[i][j][k][l][m];

        fRecMap.resize(nhypo);
        fRecOffList.resize(nhypo);
        fRecSigList.resize(nhypo);
        fRecOffErrList.resize(nhypo);
        fRecSigErrList.resize(nhypo);
        for (int i = 0; i < nhypo; i++)
        {
            fRecMap[i].resize(rad);
            fRecOffList[i].resize(rad);
            fRecSigList[i].resize(rad);
            fRecOffErrList[i].resize(rad);
            fRecSigErrList[i].resize(rad);
            for (int j = 0; j < rad; j++)
            {
                fRecMap[i][j].resize(mom);
                fRecOffList[i][j].resize(mom);
                fRecSigList[i][j].resize(mom);
                fRecOffErrList[i][j].resize(mom);
                fRecSigErrList[i][j].resize(mom);
                for (int k = 0; k < mom; k++)
                {
                    fRecMap[i][j][k].resize(the);
                    fRecOffList[i][j][k].resize(the);
                    fRecSigList[i][j][k].resize(the);
                    fRecOffErrList[i][j][k].resize(the);
                    fRecSigErrList[i][j][k].resize(the);
                    for (int l = 0; l < the; l++)
                    {
                        fRecMap[i][j][k][l].resize(nph);
                        fRecOffList[i][j][k][l].resize(nph);
                        fRecSigList[i][j][k][l].resize(nph);
                        fRecOffErrList[i][j][k][l].resize(nph);
                        fRecSigErrList[i][j][k][l].resize(nph);
                    }
                }
            }
        }
    }

    TH1D *GenRecMap(int ihypo, int irad, int imom, int ithe, int iph)
    {
        MyRICHDetector *det = GetDetScan(imom, ithe, ihypo);
        if (det == NULL)
            return NULL;
        TH2F *fhitmap = det->fRecMapEachRad[irad];
        TString particle = det->particle;
        TString radiator = det->sRadLayer[irad];
        double momentum = det->momentum;
        double theta0 = det->Theta0;
        TString fname = Form("fRecMap%d_%d_%d_%d_%d", ihypo, irad, imom, ithe, iph);
        if (gDirectory->Get(fname) != NULL)
            delete (TH1D *)gDirectory->Get(fname);
        TH1D *fproj = fhitmap->ProjectionY(fname, iph + 1, iph + 1);
        fproj->SetTitle(Form("Particle: %s [%.1f GeV/c, %.1f#circ]. %d photons from radiator %s", particle.Data(), momentum, theta0, iph, radiator.Data()));
        if (fproj->Integral() > 0)
            fproj->Fit("gaus", "Q", "", fproj->GetMean() - 4 * fproj->GetRMS(), fproj->GetMean() + 4 * fproj->GetRMS());
        return fproj;
    }

    void ResizePIDEffMap(int nmom, int nthe, int hypo)
    {
        fPidEffList.clear();

        fPidEffList.resize(nmom);
        for (int imom = 0; imom < nmom; imom++)
        {
            fPidEffList[imom].resize(nthe);
            for (int ithe = 0; ithe < nthe; ithe++)
            {
                fPidEffList[imom][ithe].resize(hypo);
                for (int ihypo = 0; ihypo < hypo; ihypo++)
                    fPidEffList[imom][ithe][ihypo].resize(hypo);
            }
        }
    }
};

#endif