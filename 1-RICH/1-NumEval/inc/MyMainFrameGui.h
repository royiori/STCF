#ifndef _MyMainFrameGui_h_
#define _MyMainFrameGui_h_

#include "TGFrame.h"
#include "TGFileDialog.h"
#include "TGCanvas.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGMsgBox.h"
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

#include "MyGuiActionClass.h"

using namespace std;

enum ShowButtonAction
{
    InitStatus,
    EnableShowSpecButton,
    EnableShowMulParButton,
    EnableShowScanRICHButton,
    EnableShowRecRICHButton,
    EnableShowPIDEffButton,
    EnableGenScanRICHButton,
    EnableGenRecRICHButton,
    EnableGenPIDEffButton,
};

class MyMainFrameGui : public TGMainFrame
{
public:
    MyMainFrameGui(const TGWindow *p, int w, int h);
    virtual ~MyMainFrameGui();

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

    void LoadText(TString txt) { fSettingText->LoadBuffer(txt); }
    TGText *GetText() { return fSettingText->GetText(); }

    void HandleCommand();

    void ToggleButList(bool, ShowButtonAction);
    void InitShowButtonList() { ToggleButList(0, ShowButtonAction::InitStatus); }
    void EnableShowSpecButton() { ToggleButList(1, ShowButtonAction::EnableShowSpecButton); }
    void EnableShowMulParButton() { ToggleButList(1, ShowButtonAction::EnableShowMulParButton); }
    void EnableShowScanRICHButton() { ToggleButList(1, ShowButtonAction::EnableShowScanRICHButton); }
    void EnableShowRecRICHButton() { ToggleButList(1, ShowButtonAction::EnableShowRecRICHButton); }
    void EnableShowPIDEffButton() { ToggleButList(1, ShowButtonAction::EnableShowPIDEffButton); }
    void EnableGenScanRICHButton() { ToggleButList(1, ShowButtonAction::EnableGenScanRICHButton); }
    void EnableGenRecRICHButton() { ToggleButList(1, ShowButtonAction::EnableGenRecRICHButton); }
    void EnableGenPIDEffButton() { ToggleButList(1, ShowButtonAction::EnableGenPIDEffButton); }

    int GetNPage() { return NPage; }

private:
    TGTextEdit *fSettingText;
    TGTab *fCTab;
    TCanvas *fCA[100];

    //vector<TGTextButton *> butList;
    map<int, TGTextButton *> butList;

    virtual void CloseWindow();
    virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

    TGComboBox *fComboCmd;
    TGTextEntry *fCommand;
    TGTextBuffer *fCommandBuf;

    int NPage;

    ClassDef(MyMainFrameGui, 1)
};

extern MyMainFrameGui *gMyMainFrameGui;

#endif