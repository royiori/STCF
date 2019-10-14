#ifndef MyBeamTestTrackVMM_h
#define MyBeamTestTrackVMM_h

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

class MyBeamTestTrackVMM : public MyBeamTestDetector
{
public:
    MyBeamTestTrackVMM(int, int, int, double, double);
    ~MyBeamTestTrackVMM(); 

    //读取DST数据
    void ReadDSTRoot(TString fPath);
    void CloseDSTRoot();

    //读取map文件
    virtual void ReadMapFile();

    //作图
    virtual void DrawHits(TGeoManager *geom, TGeoVolume *stop, TGeoMedium *med, MyBeamTestHitData *hit);
    virtual void DrawDetector(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med);
    virtual pair<double, double> MapPosition(int bd, int chip, int channel);

    virtual void AbsPosition(pair<double, double> epos, double *absPos)//epos 是电子学的x/y的channel号定义，也就是MapPosition的输出
    {
        double x0 = (epos.first == -999) ? 0 : epos.first;
        double y0 = (epos.second == -999) ? 0 : epos.second;
        double posX = ((x0 - 1) - Nstrp[0] / 2. + 1 / 2.) * strip[0]; //电子学mapping里的X/Y是这里的Y/X，从1到32
        double posY = ((y0 - 1) - Nstrp[1] / 2. + 1 / 2.) * strip[1]; //中心为0点
        posX = (epos.first == -999) ? 0 : posX;
        posY = (epos.second == -999) ? 0 : posY;


        //cout<<"Strip: ["<<strip[0]<<", "<<strip[1]<<"], mid:["<<posX<<", "<<posY<<"]"<<endl;

        //旋转
        absPos[0] = posX;
        absPos[1] = posY * sin(ROT[0]); //RICH与水平面（-Z方向）的夹角
        absPos[2] = -posY * cos(ROT[0]);

        //平移到RICH的轴心
        absPos[0] = pos[0] + absPos[0];
        absPos[1] = pos[1] + absPos[1] + LtoAnode * cos(ROT[0]);
        absPos[2] = pos[2] + absPos[2] + LtoAnode * sin(ROT[0]);
    }

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
    TH1F *fCharge;
};

#endif