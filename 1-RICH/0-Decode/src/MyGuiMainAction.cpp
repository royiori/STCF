//___ 画图 & GUI接口 ___

#include "TPad.h"
#include "TView.h"
#include "TGraph2D.h"
#include "MyGuiMainFrame.h"
#include "MyGuiMainAction.h"
#include "MyDataDecoder.h"

MyGuiMainAction *gMyGuiMainAction = (MyGuiMainAction *)0;

//____构造函数________________________________________________________________
//
MyGuiMainAction::MyGuiMainAction(const char *realpath)
{
    env = new TEnv(gSystem->WorkingDirectory() + TString("/.env"));
    env->SaveLevel(kEnvLocal);

    gMyStyle = new MyStyle();
    gMyStyle->SetDrawOption("c");
    gMyStyle->SetColorPattern(MATHEMATIC_STYLE);

    realPath = TString(realpath);

    // pages
    //#GUI. Add button content
    ButAct.push_back(ButtonAct("Data Decoder", 1)); //the second parm is 1, then add a new page
    ButAct.push_back(ButtonAct("SEP"));
    ButAct.push_back(ButtonAct("Analysis Pedestal Data"));
    ButAct.push_back(ButtonAct("Read Pedestal Analysis Result"));
    ButAct.push_back(ButtonAct("SEP"));
    ButAct.push_back(ButtonAct("Analysis AGET Data"));
    ButAct.push_back(ButtonAct("Read AGET Data"));
    ButAct.push_back(ButtonAct("Slider")); //slider must follower bby LR-L/LR-R
    ButAct.push_back(ButtonAct("LR-L"));
    ButAct.push_back(ButtonAct("LR-R"));
    ButAct.push_back(ButtonAct("Show All Events"));
    ButAct.push_back(ButtonAct("SEP"));
};

MyGuiMainAction::~MyGuiMainAction()
{
    // Destructor.
    // Store env
    env->SaveLevel(kEnvLocal);
}

//____________________________________________________________________
//
TString MyGuiMainAction::GenerateSettingsText()
{
    return "";
}

vector<TString> MyGuiMainAction::ReadContent(TString LINE)
{
    TString line = LINE;
    vector<TString> text;

    line.Remove(TString::kBoth, ' ');
    line.ReplaceAll("  ", " ");
    line.ReplaceAll("  ", " ");
    line.ReplaceAll("  ", " ");

    while (line.Length() > 0 && line.Index(' ') > 0)
    {
        TString head = line;
        line.Remove(0, line.Index(' ') + 1);
        head.Remove(head.Index(' '), head.Length());
        head.Remove(TString::kBoth, ' ');
        if (head.Length() > 0)
            text.push_back(head);
        if (line.Index(' ') == -1)
            text.push_back(line);
    }

    return text;
}

//______________________________________________________________________________
// basic button actions
void MyGuiMainAction::ExecButtonClick(Long_t bid, const char *cmdStr)
{
    //#GUI. Add button actions
    switch (bid)
    {
    case 1: //LoadBeamTextBuf
        break;
    case 2: //DrawBeamTestConfig
        break;
    case 100:
        gMyDataDecoder->AnalysisPedestalData(cmdStr);
        break;
    case 101:
        gMyDataDecoder->ReadPedestalData(cmdStr);
        break;
    case 102:
        gMyDataDecoder->AnalysisAGETData(cmdStr);
        break;
    case 103:
        gMyDataDecoder->ReadAGETData(cmdStr);
        break;
    default:
        break;
    }
}
