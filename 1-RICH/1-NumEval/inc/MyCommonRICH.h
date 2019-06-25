#ifndef _MyCommonRICH_h_
#define _MyCommonRICH_h_

#include <map>
#include <string>
#include "TF1.h"
#include "TH2.h"
#include "TRandom.h"
#include "MyRICHDetector.h"

using namespace std;

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
    MyRICHDetector *GetDetList(int ihypo) { return gDetList[ihypo]; }
    MyRICHDetector *GetDetScan(int imom, int itheta0, int ihypo) { return gScanDetList[imom][itheta0][ihypo]; }

    TH1F *GetDetWLMap() { return gDet->GetWaveLengthHist(); }
    TH1F *GetDetListWLMap() { return gDet->GetWaveLengthHist(); }

    TH2F *GetDetHitMap() { return gDet->GetHitMap(); }
    TH2F *GetDetListHitMap(int ihypo) { return gDetList[ihypo]->GetHitMap(); }
    TH2F *GetDetScanHitMap(int imom, int itheta0, int ihypo) { return gScanDetList[imom][itheta0][ihypo]->GetHitMap(); }

    int GetDetListNumber() { return (int)gDetList.size(); }

    //----------------------------
    // set functions
    void SetPrecision(double ep) { epsilon = ep; }
    void SetDatabase(MyDatabaseClass *d) { gDb = d; }

    void SaveRings(const char *fname);
    void LoadRings(const char *fname);

    //----------------------------
    // calculations
    void GenerateDetRing(MyRICHDetector *det = 0);
    void GenerateMultiParticleRICHRings();
    void GenerateTheScanHitMaps();
    void GenerateTheOffsetAndResolutionMaps();
    
    // reconstruction
    vector<double> ReconstructRICHBySolver(MyRICHDetector *det, double xc, double yc);
    void ReconstructCerekovAngleDist(MyRICHDetector *det);

private:
    double ProjectToPixel(MyRICHDetector *, double xmin, double xmax, double ymin, double ymax, double lambda);
    double ProjectToPixel(MyRICHDetector *, double xmin, double xmax, double ymin, double ymax);

    double PhotonGenFromRad(MyRICHDetector *det, int idet, int AbsFlag = 1);
    double PhotonGenFromRadWithAbs(MyRICHDetector *det, int idet) { return PhotonGenFromRad(det, idet, 1); }
    double PhotonGenFromRadWithOutAbs(MyRICHDetector *det, int idet) { return PhotonGenFromRad(det, idet, 0); }

    //----------------------------
    // support function
    vector<double> GetThickList(MyRICHDetector *, int);
    vector<double> GetRefIndList(MyRICHDetector *, int, double);
    vector<double> GetAbsLenList(MyRICHDetector *, int, double);

    // variables
    double epsilon;

    TF1 *fRhoFcn;
    MyDatabaseClass *gDb;                                  //
    MyRICHDetector *gDet;                                  //id=0, 用来保存/读取探测器设置
    vector<MyRICHDetector *> gDetList;                     //id=1~NHYPO,[粒子种类]
    vector<vector<vector<MyRICHDetector *>>> gScanDetList; //id=10.....,[动量][角度][粒子种类]

    void ClearScanList()
    {
        for (int i = 0; i < (int)gScanDetList.size(); i++)
            for (int j = 0; j < (int)gScanDetList[i].size(); j++)
                for (int k = 0; k < (int)gScanDetList[i][j].size(); k++)
                    if (gScanDetList[i][j][k] != 0)
                        delete gScanDetList[i][j][k];
        gScanDetList.clear();
    }
};

#endif