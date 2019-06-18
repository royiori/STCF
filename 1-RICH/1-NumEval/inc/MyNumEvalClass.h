#ifndef _MyNumEvalClass_h_
#define _MyNumEvalClass_h_

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
#include "MyDatabaseClass.h"
#include "MyCommonRICH.h"
#include "MyStyle.h"

using namespace std;

enum ButtonAction
{
    //...
    DrawConfig,
    LoadTextBuf,

    //material list
    MatList = 10,

    //detector list
    DetList = 50,

    //analysis action
    AnalysisAction = 100,
    DrawSelectedMat,
    DrawSelectedDet,
    DrawRICH,
    Draw3Rings,
    Save3Rings,
    Load3Rings,
};

//-----------------------------------------------------------------------------------------------
//
class MyNumEvalClass
{
public:
    MyNumEvalClass();
    ~MyNumEvalClass();

    int GetNTabPage() { return nTabPage; }
    TString GetTabPageName(int i) { return sTabPage[i]; }

    int GetNPageButton(int i) { return nButton[i]; }
    int GetIndexButton(int i) { return iButton[i]; }
    TString GetPageButtonName(int i, int j) { return sButton[i].at(j); }

    TString GenerateSettingsText();
    void ReadSettingsText();

    void ExecButtonClick(Long_t bid, const char *cmd);

private:
    //gui related action
    int nTabPage;
    TString sTabPage[10];

    int nButton[10];
    int iButton[10];
    vector<TString> sButton[10];

    int nSelectedMat, nSelectedDet;
    vector<TString> selectedMat;
    vector<TString> selectedDet;

    bool ValicationCheck();
    void DoDrawConfig(TString str);
    void DoLoadTextBuf();
    void ShowMaterialInfo(TString);
    void ShowDetectorInfo(TString);
    void DoDrawSelectedMat();
    void DoDrawSelectedDet();
    void DoDrawRICH();
    void DoDrawRings();
    void DoSaveRings(const char *fname);
    void DoLoadRings(const char *fname);


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
    double momentum, theta0;
    int NXBin, NYBin;
    double XBinMin, XBinMax, YBinMin, YBinMax;

    void SetDetectorParameters(int id = -1);

    //分析精度要求
    double epsilon;
    double trkStep;
    int nphi;

    //..
    TEnv *env;
    TString filePath;
    MyStyle *gMyStyle;
    MyCommonRICH *gMyCommonRICH;
    MyDatabaseClass *gMyDatabaseClass;
};

extern MyNumEvalClass *gMyNumEvalClass;

#endif
