#ifndef _MyGuiMainFrame_h_
#define _MyGuiMainFrame_h_

#include "TGFrame.h"
#include "TGFileDialog.h"
#include "TGCanvas.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGMsgBox.h"
#include "TGSlider.h"
#include "TGTab.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGTextEdit.h"
#include "TGComboBox.h"
#include "TG3DLine.h"
#include "TGClient.h"
#include "TGResourcePool.h"
#include "TStyle.h"
#include "TStyleManager.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "Riostream.h"

#include "MyGuiMainAction.h"

using namespace std;

enum ShowButtonAction
{
    InitStatus,
};

class MyGuiMainFrame : public TGMainFrame
{
public:
    MyGuiMainFrame(const TGWindow *p, int w, int h);
    virtual ~MyGuiMainFrame();

    void ClearCanvas(int id) { fCA[id]->Clear(); }
    void ClearAllCanvas()
    {
        for (int i = 1; i < 100; i++)
            if (fCA[i] != NULL)
                fCA[i]->Clear();
    }
    void UpdateCanvas(int id)
    {
        fCA[id]->Modified();
        fCA[id]->Update();
    }
    void SwitchCanvas(int id)
    {
        fCTab->SetTab(id);
        fCA[id]->cd();
    }
    void DivideCanvas(int id, int x, int y) { fCA[id]->Divide(x, y); }
    void SwitchDivCanvas(int id, int cid) { fCA[id]->cd(cid); }
    TCanvas *GetCanvas(int id) { return fCA[id]; }

    void LoadText1(TString txt) { fSettingText1->LoadBuffer(txt); }
    TGText *GetText1() { return fSettingText1->GetText(); }

    void HandleCommand();

    void ToggleButList(bool, ShowButtonAction);
    void InitShowButtonList() { ToggleButList(0, ShowButtonAction::InitStatus); }

    void SetSliderRange(int, ShowButtonAction);

private:
    TGTextEdit *fSettingText1;
    TGTab *fCTab;
    TCanvas *fCA[100];

    map<int, TGTextButton *> butList;
    map<int, TGSlider *> sldList;

    virtual void CloseWindow();
    virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

    TGComboBox *fComboCmd;
    TGTextEntry *fCommand;
    TGTextBuffer *fCommandBuf;

    const int NResPAGE = 11;

    ClassDef(MyGuiMainFrame, 1)
};

extern MyGuiMainFrame *gMyGuiMainFrame;

#endif