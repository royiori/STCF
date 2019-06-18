#ifndef _MyRICHDetector_h_
#define _MyRICHDetector_h_

#include <map>
#include <string>
#include "TF1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TObject.h"

using namespace std;

class MyRICHDetector : public TObject
{
public:
    MyRICHDetector()
    {
        Init();
    }

    MyRICHDetector(const MyRICHDetector &det)
    {
        ((MyRICHDetector &)det).MyCopy(*this);
    }

    ~MyRICHDetector()
    {
        if (fRing != 0)
            delete fRing;
        for (int i = 0; i < 10; i++)
            if (fR1F[i] != 0)
                delete fR1F[i];
    };

    //定义探测器结构
    int nRadLayer;            //可以有多层辐射体层, 最上层为起始层0
    vector<double> tRadLayer; //辐射体层厚度，thickness
    vector<string> sRadLayer;
    void SetRadiator(int n, vector<double> t, vector<string> s)
    {
        nRadLayer = n;
        tRadLayer = t;
        sRadLayer = s;
    }
    double GetRadTotLeng()
    {
        double len = 0;
        for (int i = 0; i < nRadLayer; i++)
            len += tRadLayer[i];
        return len;
    }
    string GetLayerName(double l)
    {
        for (int i = 0; i < nRadLayer; i++)
        {
            l -= tRadLayer[i];
            if (l <= 0)
                return sRadLayer[i];
        }
        return "";
    }

    int nTransLayer; // 只有一层传输层，但是这层里可以是多个气体的混合
    double tTransLayer;
    vector<double> pTransLayer; //气体混合百分比，percentage
    vector<string> sTransLayer;
    void SetTransLayer(double t, int n, vector<double> p, vector<string> s)
    {
        tTransLayer = t;
        nTransLayer = n;
        pTransLayer = p;
        sTransLayer = s;
    }

    int nImpurities;            // 传输层里的杂质成分
    vector<double> pImpurities; //ppm
    vector<string> sImpurities;
    void SetImpurities(int n, vector<double> p, vector<string> s)
    {
        nImpurities = n;
        pImpurities = p;
        sImpurities = s;
    }

    string PCDetector; // 读出的光电器件类型
    double pixel;
    void SetDetector(string det, double pix)
    {
        PCDetector = det;
        pixel = pix;
        NXBin = (XBinMax - XBinMin) / pixel;
        NYBin = (YBinMax - YBinMin) / pixel;
    }

    //计算分析所需要用的参数
    int nLambda;
    double lambdaMin, lambdaMax, lambdaStep;
    void SetLambdaRange(int n, double min, double max)
    {
        nLambda = n;
        lambdaMin = min;
        lambdaMax = max;
        lambdaStep = (max - min) / n;
    }

    TString particle;
    double mass, momentum, theta0;
    void SetParticleGun(TString p, double m, double mom, double t0)
    {
        particle = p;
        mass = m;
        momentum = mom;
        theta0 = t0;
    }

    int NXBin, NYBin;
    int XBinMin, XBinMax, YBinMin, YBinMax;
    void SetDetectorViewSize(double xmin, double xmax, double ymin, double ymax)
    {
        XBinMin = xmin;
        XBinMax = xmax;
        YBinMin = ymin;
        YBinMax = ymax;
        NXBin = (XBinMax - XBinMin) / pixel;
        NYBin = (YBinMax - YBinMin) / pixel;
    }

    double epsilon;
    double trkStep;
    int nphi;
    void SetPrecision(double ep, double tp, int np)
    {
        epsilon = ep;
        trkStep = tp;
        nphi = np;
    }

    TH2F *Get2DViewer(int flag = 0)
    {
        if (flag == 0 && fRing != 0)
            return fRing;

        if (fRing != 0)
            delete fRing;
        fRing = new TH2F(Form("fRing-%d", id), Form("Cherenkov Ring for %s @ %.1f GeV", particle.Data(), momentum), NXBin, XBinMin, XBinMax, NYBin, YBinMin, YBinMax);
        fRing->GetXaxis()->SetTitle("X[mm]");
        fRing->GetYaxis()->SetTitle("Y[mm]");
        return fRing;
    }

    TH1D **Get1DViewer()
    {
        for (int i = 0; i < 10; i++)
            if (fR1F[i] != 0)
                delete fR1F[i];

        fR1F[0] = new TH1D(Form("fRingX-%d", id), Form("Cherenkov Ring for %s @ %.1f GeV, X - profile", particle.Data(), momentum), NXBin, XBinMin, XBinMax);
        fR1F[1] = new TH1D(Form("fRingY-%d", id), Form("Cherenkov Ring for %s @ %.1f GeV, Y - profile", particle.Data(), momentum), NYBin, YBinMin, YBinMax);
        fR1F[2] = new TH1D(Form("fRingX-all-%d", id), Form("Cherenkov Ring for %s @ %.1f GeV, X - projected", particle.Data(), momentum), NXBin, XBinMin, XBinMax);
        fR1F[3] = new TH1D(Form("fRingY-all-%d", id), Form("Cherenkov Ring for %s @ %.1f GeV, Y - projected", particle.Data(), momentum), NYBin, YBinMin, YBinMax);

        for (int i = 0; i < 4; i++)
            fR1F[i]->GetYaxis()->SetTitle("Entries");
        for (int i = 0; i < 4; i += 2)
            fR1F[i]->GetXaxis()->SetTitle("X[mm]");
        for (int i = 1; i < 4; i += 2)
            fR1F[i]->GetXaxis()->SetTitle("Y[mm]");

        if (fRing == NULL)
            return fR1F;

        fR1F[0] = fRing->ProjectionX(fR1F[0]->GetName());
        fR1F[1] = fRing->ProjectionY(fR1F[1]->GetName());
        fR1F[2] = fRing->ProjectionX(fR1F[2]->GetName(), NXBin/2, NXBin/2);
        fR1F[3] = fRing->ProjectionY(fR1F[3]->GetName(), NYBin/2, NYBin/2);
        return fR1F;
    }

    int id;
    TH1D *fR1F[10]; // xprofile, yprofile, x-all-project, y-all-project
    TH2F *fRing;

    void Init()
    {
        fRing = 0;
        for (int i = 0; i < 10; i++)
            fR1F[i] = 0;
        id = (int)gRandom->Uniform(10000);
    };

    void BuildFrom(MyRICHDetector &det)
    {
        ((MyRICHDetector &)det).MyCopy(*this);
    }

private:
    void MyCopy(MyRICHDetector &det)
    {
        det.SetRadiator(nRadLayer, tRadLayer, sRadLayer);
        det.SetTransLayer(tTransLayer, nTransLayer, pTransLayer, sTransLayer);

        det.SetImpurities(nImpurities, pImpurities, sImpurities);
        det.SetDetector(PCDetector, pixel);
        det.SetParticleGun(particle, mass, momentum, theta0);
        det.SetLambdaRange(nLambda, lambdaMin, lambdaMax);
        det.SetDetectorViewSize(XBinMin, XBinMax, YBinMin, YBinMax);

        det.SetPrecision(epsilon, trkStep, nphi);
        det.Init();
    }

    ClassDef(MyRICHDetector, 1) //Jet class
};

#endif