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
    MyRICHDetector *GetDetector(int id = -1) { return (id == -1) ? gDet : gScanDetList[id]; }

    TH2F *Get2DViewer(int id = -1) { return (id == -1) ? gDet->Get2DViewer() : ((0 <= id && id < int(gScanDetList.size())) ? gScanDetList[id]->Get2DViewer() : NULL); }
    TH1D **Get1DViewer(int id = -1) { return (id == -1) ? gDet->Get1DViewer() : ((0 <= id && id < int(gScanDetList.size())) ? gScanDetList[id]->Get1DViewer() : NULL); }

    //----------------------------
    // set functions
    void SetPrecision(double ep) { epsilon = ep; }
    void SetTrackStepSize(double sp) { trkStep = sp; }
    void SetNumberOfPhiStep(int _np) { nphi = _np; }
    void SetDatabase(MyDatabaseClass *d) { gDb = d; }

    void SaveRings(const char *fname);
    void LoadRings(const char *fname);

    //----------------------------
    // calculations
    void GenerateRICHRing(MyRICHDetector *det = 0);
    void GenerateMultiParticleRICHRings(vector<TString> parList);

    void ReconstructRICH(MyRICHDetector *det = 0, TString hypo="pi");
private:
    double ProjectToPixel(MyRICHDetector *, double xmin, double xmax, double ymin, double ymax, double lambda);
    double ProjectToPixel(MyRICHDetector *, double xmin, double xmax, double ymin, double ymax);

    double PhotonGenFromRad(MyRICHDetector *det, int idet, double lambda, TH2F *fRing, double scal = 1., int AbsFlag = 1);
    double PhotonGenFromRadWithAbs(MyRICHDetector *det, int idet, double lambda, TH2F *fRing, double scal) { return PhotonGenFromRad(det, idet, lambda, fRing, scal, 1); }
    double PhotonGenFromRadWithOutAbs(MyRICHDetector *det, int idet, double lambda, TH2F *fRing, double scal) { return PhotonGenFromRad(det, idet, lambda, fRing, scal, 0); }

    //----------------------------
    // support function
    vector<double> GetThickList(MyRICHDetector *, int);
    vector<double> GetRefIndList(MyRICHDetector *, int, double);
    vector<double> GetAbsLenList(MyRICHDetector *, int, double);

    // variables
    double epsilon;
    double trkStep;
    int nphi;

    TF1 *fRhoFcn;

    MyRICHDetector *gDet;
    MyDatabaseClass *gDb;
    vector<MyRICHDetector *> gScanDetList;

    // for calculation
    double fun(double l0);
    void ClearScanList()
    {
        for (int i = 0; i < (int)gScanDetList.size(); i++)
            if (gScanDetList[i] != 0)
                delete gScanDetList[i];
        gScanDetList.clear();
    }
};

#endif