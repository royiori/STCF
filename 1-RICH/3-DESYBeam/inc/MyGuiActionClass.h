#ifndef _MyGuiActionClass_h_
#define _MyGuiActionClass_h_

#include <map>
#include <string>
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBox.h"
#include "TText.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TEnv.h"

#include "MyMainFrameGui.h"
#include "MyBeamTest.h"

#include "MyStyle.h"

using namespace std;

enum ButtonAction
{
    //...
    DrawConfig,
    LoadTextBuf,
    DrawBeamTestConfig, //for beamtest
    LoadBeamTextBuf,

    //material list
    MatList = 10,

    //detector list
    DetList = 50,

    //analysis action
    AnalysisAction = 100,
    DrawSelectedMat,
    DrawSelectedDet,

    LoadDetFile,
    ShowFCN,
    ShowSpecRICH,
    ShowMulParRICH,
    ShowScanRICHList,
    ShowRecRICHList,
    ShowPIDEff,
    SEP1,       // SEP用来显示一个间隔条，在cpp里定义名字时，用SEP开头即可
    SEP2,
    GenSpecRICH,
    GenMulParRICH,
    SEP3,
    SaveDetFile,
    GenScanRICHList,
    GenRecRICHList,
    GenPIDEff,

    //beam-test data analysis
    SetDSPath = 200,
    ReadRawData,
    CheckRawRoot,
    GenDSTRoot,
    CheckDSTRoot,
    CombineDSTRoot,
    LoadDSTRoot,
    //..
    SlidDSTHit,   // Slid & Prev & Next 三个需要一起出现，用来显示一个控制条、 在cpp里定义名字时，用Slider开头即可
    PrevDSTHit,
    NextDSTHit,
    //..
    ZoomToDefault,
    SEP4,
    SEP5,
    AnalysisCMDRoot,
};

const int BATCH = 1;
const int GUI = 2;

//-----------------------------------------------------------------------------------------------
//
class MyGuiActionClass
{
public:
    MyGuiActionClass(const char *path, int flag = GUI); //默认用GUI模式
    ~MyGuiActionClass();

    int GetNTabPage() { return nTabPage; }
    TString GetTabPageName(int i) { return sTabPage[i]; }

    int GetNPageButton(int i) { return nButton[i]; }
    int GetIndexButton(int i) { return iButton[i]; }
    TString GetPageButtonName(int i, int j) { return sButton[i].at(j); }

    TString GenerateSIMSettingsText();
    int ReadSIMSettingsText(const char *fname = NULL);

    void DoReadBatchFile(const char *fname, const char *epoch);
    void ExecButtonClick(Long_t bid, const char *cmd);
    vector<TString> ReadContent(TString LINE);

    //for beamtest
    TString GenerateBEAMSettingsText();
    int ReadBEAMSettingsText();
    void SetRayTracing(bool tr) { gMyBeamTest->SetRayTracing(tr); }

private:
    int BatGuiFlag;
    TString realPath;

    //gui related action
    int nTabPage;
    TString sTabPage[10];

    int nButton[10];
    int iButton[10];
    vector<TString> sButton[10];

    int nSelectedMat, nSelectedDet;
    vector<TString> selectedMat;
    vector<TString> selectedDet;

    //定义探测器结构
    int nRadLayer;            //可以有多层传输层
    vector<double> tRadLayer; //传输层厚度，thickness
    vector<string> sRadLayer;

    int nTransLayer; // 只有一层传输层，但是这层里可以是多个气体的混合
    double tTransLayer;
    vector<double> pTransLayer; //气体混合百分比，percentage
    vector<string> sTransLayer;
    map<string, double> mapTransLayers;

    int nImpurities;            // 传输层里的杂质成分
    vector<double> pImpurities; //ppm
    vector<string> sImpurities;
    map<string, double> mapImpurities;

    string PCDetector; // 读出的光电器件类型
    double pixel;

    //计算分析所需要用的参数
    int nLambda;
    double lambdaMin, lambdaMax;

    TString particle;
    double mass, momentum, Theta0;

    int NXBin, NYBin;
    double XBinMin, XBinMax, YBinMin, YBinMax;

    //扫描的参数范围
    int np, nthe0;
    double pMin, pMax;
    double The0Min, The0Max;

    void SetDetectorParameters();
    void GetDetectorParameters();

    //分析精度要求
    double epsilon;
    double trkStep;
    int nphi;

    //以下参数仅仅用于GUI显示，只保存到env，不会保存到root
    int iph, irad1, irad2; //photon & radiator to show
    int nEvent;
    double pList, thList; //指定显示范围

    //..
    TEnv *env;
    TString filePath;
    MyStyle *gMyStyle;

    //.. button action
    //for beamtest
    void DoLoadBeamTextBuf();
    void DoDrawBeamTestConfig();

    void DoSetDSPath(const char *fname);
    void DoReadRawData(TString cmdStr);
    void DoCheckRawRoot(TString cmdStr);
    void DoGenDSTRoot(TString cmdStr);
    void DoCheckDSTRoot(TString cmdStr);
    void DoCombineDSTRoot(TString cmdStr);
    void DoLoadDSTRoot(TString cmdStr);
    void DoAnalysisCMDRoot();
    void DoSlidDSTHit(TString cmdStr)
    {
        gMyBeamTest->DrawCurtHit(cmdStr.Atoi());
        DoZoomToDefault();
    }
    void DoPrevDSTHit()
    {
        gMyBeamTest->DrawPrevHit();
        DoZoomToDefault();
    }
    void DoNextDSTHit()
    {
        gMyBeamTest->DrawNextHit();
        DoZoomToDefault();
    }

    void DoZoomToDefault();

private:
    bool CheckIdelDatFile(TString cmdStr)
    {
        if (cmdStr == "no")
            return false;

        FileStat_t fs;
        if (gSystem->GetPathInfo(gMyBeamTest->GetDSPath(), fs))
        {
            cout << "#### Error: " << gMyBeamTest->GetDSPath() << " doesn't exist." << endl;
            return false;
        }

        if (!gMyBeamTest->GetDSPath().EndsWith("idle.dat"))
        {
            cout << "#### Error, please select the 'idle.dat' file." << endl;
            return false;
        }
        return true;
    }
};

extern MyGuiActionClass *gMyGuiActionClass;
#endif
