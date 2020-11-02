#ifndef MyBeamTestRICH_h
#define MyBeamTestRICH_h

#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "MyBeamTestDetector.h"

using namespace std;

class MyBeamTestRICH : public MyBeamTestDetector
{
public:
    MyBeamTestRICH(int);
    virtual ~MyBeamTestRICH();

    //读取 BIN 数据并生成 RAW.root 文件
    bool ReadData2RawRoot(vector<TString> datList, TString fRawName, int force = 1);
    void AnalysisPedestal(TString fRawName, TString fPedName, vector<MyBeamTestRICH *> vRICH);
    void AnalysisPedestal(TString fRawName, TString fPedName, vector<int> boardName, vector<int> chipName, int force = 1);
    bool ReadPedestal(TString fPedName, vector<MyBeamTestRICH *> vRICH);

    //读取 RAW 数据并生成 DST.root 文件
    bool ReadRaw2DstRoot(TString fRawName, TString fPedName, TString fDstPath, vector<MyBeamTestRICH *> vRICH, bool SaveWaveFlag);

    //读取DST数据
    void ReadDSTRoot(TString fPath);
    void CloseDSTRoot();

    //读取map文件
    virtual void ResizeMap();
    virtual void ReadMapFile();

    //作图
    virtual void DrawDetector(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med);
    virtual void DrawHits(TGeoManager *geom, TGeoVolume *stop, TGeoMedium *med, MyBeamTestHitData *hit);
    virtual pair<double, double> MapPosition(int bd, int chip, int channel);
    virtual void AbsPosition(pair<double, double> epos, double *absPos);
    virtual void GenPedMap() {};

    //分布
    virtual bool FillDistribution(int i, MyBeamTestHitData *hit);
    TH2F *GetAllHit() { return fAllHit; }
    TH2F *GetAllHit2() { return fAllHit2; }
    TH1F *GetCharge() { return fCharge; }

    //生成DST-root时，将hit合并为cluster
    void AnalysisCluster(MyBeamTestData *fEventList);

private:
    TH2F *fAllHit;
    TH2F *fAllHit2;
    TH2F *fPed;
    TH1F *fCharge;
};

#endif