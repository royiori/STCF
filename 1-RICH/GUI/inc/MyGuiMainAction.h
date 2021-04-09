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

#include "MyStyle.h"

using namespace std;

//----------------------------------------------------
// button action structure
class ButtonAct
{
public:
    ButtonAct(TString tit, int ig = -1)
    {
        title = tit;
        iPG = ig;
    }
    ~ButtonAct() {};

    void setID(int id) { uID = id; }
    bool isMe(int id) { return (id == uID); }

    int uID;
    int iPG; //page id
    TString title;
};

//-----------------------------------------------------------------------------------------------
//
class MyGuiMainAction
{
public:
    MyGuiMainAction(const char *path); //默认用GUI模式
    ~MyGuiMainAction();

    //for GUI
    vector<ButtonAct> GetButAct() { return ButAct; }

    // settings-text related actions
    TString GenerateSettingsText();
    int ReadSettingsText();
    vector<TString> ReadContent(TString LINE);

    // button-page related actions 
    void ExecButtonClick(Long_t bid, const char *cmd);

private:
    TString realPath;

    //gui related action
    vector<ButtonAct> ButAct; //button action

    //..
    TEnv *env;
    TString filePath;
    MyStyle *gMyStyle;

private:
    bool CheckIdelDatFile(TString cmdStr);
};

extern MyGuiMainAction *gMyGuiMainAction;
#endif
