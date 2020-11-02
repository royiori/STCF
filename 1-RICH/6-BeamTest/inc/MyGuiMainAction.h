#ifndef _MyGuiMainAction_h_
#define _MyGuiMainAction_h_

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

#include "MyGuiMainFrame.h"
#include "MyGuiBeamTest.h"

#include "MyStyle.h"

using namespace std;

enum ButtonAction
{
    //...
    DrawConfig,
    LoadTextBuf,
    DrawBeamTestConfig, //for beamtest
    LoadBeamTextBuf,

    //beam-test data analysis
    SetDSPath = 200,
    ReadRawData,
    CheckRawRoot,
    GenDSTRoot,
    CheckDSTRoot,
    CombineDSTRoot,
    LoadDSTRoot,
    //..
    SlidDSTHit, // Slid & Prev & Next 三个需要一起出现，用来显示一个控制条、 在cpp里定义名字时，用Slider开头即可
    PrevDSTHit,
    NextDSTHit,
    //..
    ZoomToDefault,
    SEP4,
    SEP5,
    AnalysisCMDRoot,
};

//-----------------------------------------------------------------------------------------------
//
class MyGuiMainAction
{
public:
    MyGuiMainAction(const char *path); //默认用GUI模式
    ~MyGuiMainAction();

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

private:
    TString realPath;

    //gui related action
    int nTabPage;
    TString sTabPage[10];

    int nButton[10];
    int iButton[10];
    vector<TString> sButton[10];

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
    void DoSlidDSTHit(TString cmdStr);
    void DoPrevDSTHit();
    void DoNextDSTHit();
    void DoZoomToDefault();

private:
    bool CheckIdelDatFile(TString cmdStr);
};

extern MyGuiMainAction *gMyGuiMainAction;
#endif
