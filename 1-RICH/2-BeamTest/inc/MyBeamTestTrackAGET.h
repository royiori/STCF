#ifndef MyBeamTestTrackAGET_h
#define MyBeamTestTrackAGET_h

#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "MyBeamTestDetector.h"

using namespace std;

class MyBeamTestTrackAGET : public MyBeamTestDetector
{
public:
    MyBeamTestTrackAGET(int, int, int, double, double);
    ~MyBeamTestTrackAGET();

    //读取DST数据
    void ReadDSTRoot(TString fPath);
    void CloseDSTRoot();

    //读取map文件
    virtual void ReadMapFile();

    //作图
    virtual void DrawHits(TGeoManager *geom, TGeoVolume *stop, TGeoMedium *med, MyBeamTestHitData *hit);
    virtual void DrawDetector(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med);
    virtual pair<double, double> MapPosition(int bd, int chip, int channel);
    virtual void AbsPosition(pair<double, double> epos, double *absPos);
    virtual void GenPedMap() {;}

    //分布
    virtual bool FillDistribution(int i, MyBeamTestHitData *hit);
    TH2F *GetAllHit() { return fAllHit; }
    TH2F *GetAllHit2() { return fAllHit2; }
    TH1F **GetHitXY() { return fHit; }
    TH1F *GetCharge() { return fCharge; }

    //生成DST-root时，将hit合并为cluster
    void AnalysisCluster(MyBeamTestData* fEventList);

private:
    TH2F *fAllHit;
    TH2F *fAllHit2;
    TH1F *fHit[2];
    TH1F *fPed[2];
    TH1F *fCharge;
};

#endif