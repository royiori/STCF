#include "MyMainFrameGui.h"

MyMainFrameGui *gMyMainFrameGui = (MyMainFrameGui *)0;

using namespace std;

const char *filetypes[] = {"ROOT files", "*.root",
                           "All files", "*",
                           0, 0};

const int NPAGE = 10;
//______________________________________________________________________________
//
MyMainFrameGui::MyMainFrameGui(const TGWindow *p, int w, int h) : TGMainFrame(p, w, h)
{
    TGGC myGC = *gClient->GetResourcePool()->GetFrameGC();
    TGFont *myfont = gClient->GetFont("-adobe-helveticaa-n-r-*-*-16-*-*-*-*-*-iso8859-1");
    if (myfont)
        myGC.SetFont(myfont->GetFontHandle());

    for (int i = 0; i < 100; i++)
        fCA[i] = 0;

    TGLayoutHints *LayoutC = new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2);
    TGLayoutHints *LayoutX = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2);
    TGLayoutHints *LayoutY = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY, 2, 2, 2, 2);
    TGLayoutHints *LayoutXY = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);

    Pixel_t green, yellow;
    fClient->GetColorByName("green", green);
    fClient->GetColorByName("yellow", yellow);

    TGHorizontalFrame *fHFrame0 = new TGHorizontalFrame(this, 10, 10);
    AddFrame(fHFrame0, LayoutXY);

    TGVerticalFrame *fVFrame1 = new TGVerticalFrame(fHFrame0, 10, 10);
    TGVerticalFrame *fVFrame2 = new TGVerticalFrame(fHFrame0, 10, 10);
    TGVerticalFrame *fVFrame3 = new TGVerticalFrame(fHFrame0, 10, 10);
    fHFrame0->AddFrame(fVFrame1, LayoutY);
    fHFrame0->AddFrame(fVFrame2, LayoutY);
    fHFrame0->AddFrame(fVFrame3, LayoutXY);

    //1. Setting Text
    {
        fSettingText = new TGTextEdit(fVFrame1, 400, 580);
        fSettingText->LoadBuffer(gMyGuiActionClass->GenerateSettingsText());
        fSettingText->SetFont(myfont->GetFontStruct());
        fVFrame1->AddFrame(fSettingText, LayoutXY);

        TGHorizontalFrame *fHFrame10 = new TGHorizontalFrame(fVFrame1, 10, 10, kHorizontalFrame);
        fVFrame1->AddFrame(fHFrame10, LayoutX);

        TGTextButton *fButtonTmp1 = new TGTextButton(fHFrame10, "Load Set-text", LoadTextBuf);
        fButtonTmp1->Associate(this);
        fHFrame10->AddFrame(fButtonTmp1, LayoutX);

        TGTextButton *fButtonTmp2 = new TGTextButton(fHFrame10, "Draw Configure", DrawConfig);
        fButtonTmp2->Associate(this);
        fHFrame10->AddFrame(fButtonTmp2, LayoutX);
    }

    //2. Buttons Tab
    {
        TGTab *fTab0 = new TGTab(fVFrame2, 400, 580);
        fTab0->SetTab(0);
        fTab0->Resize(fTab0->GetDefaultSize());
        fVFrame2->AddFrame(fTab0, LayoutXY);

        // 2.1 page0
        for (int i = 0; i < gMyGuiActionClass->GetNTabPage(); i++)
        {
            TGCompositeFrame *fTabPage = fTab0->AddTab(gMyGuiActionClass->GetTabPageName(i));
            fTabPage->SetLayoutManager(new TGVerticalLayout(fTabPage));

            for (int j = 0; j < gMyGuiActionClass->GetNPageButton(i); j++)
            {
                if (gMyGuiActionClass->GetPageButtonName(i, j) == "SEP")
                {
                    TGHorizontal3DLine *fHorizontal3DLine = new TGHorizontal3DLine(fTabPage, 102, 8);
                    fTabPage->AddFrame(fHorizontal3DLine, LayoutX);
                }
                else
                {
                    TGTextButton *fButtonTmp = new TGTextButton(fTabPage, gMyGuiActionClass->GetPageButtonName(i, j), gMyGuiActionClass->GetIndexButton(i) + j);
                    fButtonTmp->Associate(this);
                    fTabPage->AddFrame(fButtonTmp, LayoutX);
                    if (gMyGuiActionClass->GetPageButtonName(i, j).BeginsWith("Load"))
                        fButtonTmp->ChangeBackground(green);
                    if (gMyGuiActionClass->GetPageButtonName(i, j).BeginsWith("Set"))
                        fButtonTmp->ChangeBackground(yellow);
                    if (gMyGuiActionClass->GetPageButtonName(i, j).BeginsWith("Draw"))
                        fButtonTmp->ChangeBackground(yellow);
                    butList.insert(pair<int, TGTextButton *>(gMyGuiActionClass->GetIndexButton(i) + j, fButtonTmp));
                }
            }
        }
    }

    //3. Results tab
    {
        fCTab = new TGTab(fVFrame3, 580, 580);
        NPage = NPAGE;

        for (int i = 0; i < NPage; i++)
        {
            TGCompositeFrame *fTabPage = fCTab->AddTab((i == 0) ? "Config" : Form("Res%d", i));
            TRootEmbeddedCanvas *fEmbeddedCanvas = new TRootEmbeddedCanvas(0, fTabPage, 580, 360);
            fTabPage->AddFrame(fEmbeddedCanvas, LayoutXY);
            fCA[i] = fEmbeddedCanvas->GetCanvas();
            fCA[i]->SetFillColor(0);
            fCA[i]->SetBorderMode(0);
        }

        TGCompositeFrame *fTabPage = fCTab->AddTab("About");
        TRootEmbeddedCanvas *fEmbeddedCanvas = new TRootEmbeddedCanvas(0, fTabPage, 580, 360);
        fTabPage->AddFrame(fEmbeddedCanvas, LayoutXY);
        TPaveText *pt = new TPaveText(0.2, 0.2, 0.8, 0.8);
        pt->AddText("STCF RICH-PID software package");
        pt->AddText("Author: LIUQ");
        pt->AddText("You can follow me on https://github.com/royiori");
        pt->Draw();

        fCTab->SetTab(0);
        fCTab->Resize(fCTab->GetDefaultSize());
        fVFrame3->AddFrame(fCTab, LayoutXY);

        TGHorizontalFrame *fHFrame31 = new TGHorizontalFrame(fVFrame3, 10, 10);
        fComboCmd = new TGComboBox(fHFrame31, "", 1);
        fCommand = fComboCmd->GetTextEntry();
        fCommandBuf = fCommand->GetBuffer();
        fCommand->Connect("ReturnPressed()", "MyMainFrameGui", this, "HandleCommand()");
        fComboCmd->Resize(200, fCommand->GetDefaultHeight());
        fHFrame31->AddFrame(new TGLabel(fHFrame31, "Command:"), LayoutC);
        fHFrame31->AddFrame(fComboCmd, LayoutX);
        fVFrame3->AddFrame(fHFrame31, LayoutX);
    }

    ToggleButList(0, 0);

    SetWindowName("STCF-RICH");
    SetIconName("STCF-RICH");
    MapSubwindows();
    Resize(GetDefaultSize()); // this is used here to init layout algoritme
    MapWindow();
}

MyMainFrameGui::~MyMainFrameGui()
{
    // Destructor.
    CloseWindow();
}

void MyMainFrameGui::CloseWindow()
{
    // Destructor.
    delete gMyGuiActionClass;
    DeleteWindow();
    gApplication->Terminate(0);
}

//______________________________________________________________________________
//

Bool_t MyMainFrameGui::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
    switch (GET_MSG(msg))
    {
    case kC_COMMAND:

        switch (GET_SUBMSG(msg))
        {

        case kCM_BUTTON:
        case kCM_MENU:
            static TString dir(".");
            TGFileInfo fi;
            fi.fFilename = NULL;
            fi.fFileTypes = filetypes;
            fi.fIniDir = StrDup(dir);
            TString cmdStr;

            if (parm1 == SaveDetFile || parm1 == LoadDetFile || parm1 == SaveRecFile || parm1 == LoadRecFile || parm1 == ReadCosmicData || parm1 == LoadCosmicRes || parm1 == SaveCosmicRes)
            {
                new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
                if (fi.fFilename == NULL)
                    return kTRUE;
                cmdStr = fi.fFilename;
                if (!cmdStr.EndsWith("root"))
                    cmdStr += ".root";
            }

            if (parm1 == GenSpecRICH || parm1 == GenMulParRICH || parm1 == GenScanRICHList || parm1 == GenRecRICHList || parm1 == GenPIDEff || parm1 == AnalysisCosmicData)
            {
                int retval;
                new TGMsgBox(gClient->GetRoot(), this, "Message-RICH", "Do you want RE-generate the histograms?",
                             kMBIconExclamation, kMBYes | kMBNo, &retval);
                cmdStr = (retval == 1) ? "yes" : "no";
            }

            if (parm1 == GenSpecRICH || parm1 == SaveDetFile)
                ToggleButList(1, 1);
            if (parm1 == GenMulParRICH || parm1 == SaveDetFile)
                ToggleButList(1, 2);
            if (parm1 == LoadDetFile || parm1 == GenScanRICHList)
                ToggleButList(1, 3);
            if (parm1 == LoadRecFile || parm1 == GenRecRICHList)
                ToggleButList(1, 4);
            if (parm1 == LoadRecFile || parm1 == GenPIDEff)
                ToggleButList(1, 5);
            if (parm1 == SaveDetFile)
                ToggleButList(1, 6);
            if (parm1 == SaveRecFile)
                ToggleButList(1, 7);

            gMyGuiActionClass->ExecButtonClick(parm1, cmdStr.Data());
            break;
        }
        break;
    }
    parm2 = 0.;
    return kTRUE;
}

//______________________________________________________________________________
//

void MyMainFrameGui::HandleCommand()
{
    const char *string = fCommandBuf->GetString();
    if (strlen(string) > 1)
    {
        gROOT->ProcessLine(string);
        fComboCmd->InsertEntry(string, 0, -1);
        fCommand->Clear();
    }
}

void MyMainFrameGui::ToggleButList(bool on, int flag)
{
    EButtonState stat = (on) ? kButtonUp : kButtonDisabled;

    map<int, TGTextButton *>::iterator mIter;

    if (flag == 0) //初始状态
    {
        butList.find(ShowSpecRICH)->second->SetState(stat);
        butList.find(ShowMulParRICH)->second->SetState(stat);
        butList.find(ShowScanRICHList)->second->SetState(stat);
        butList.find(ShowRecRICHList)->second->SetState(stat);
        butList.find(ShowPIDEff)->second->SetState(stat);
        butList.find(GenScanRICHList)->second->SetState(stat);
        butList.find(GenRecRICHList)->second->SetState(stat);
        butList.find(GenPIDEff)->second->SetState(stat);
    }
    if (flag == 1)
        butList.find(ShowSpecRICH)->second->SetState(stat);
    if (flag == 2)
        butList.find(ShowMulParRICH)->second->SetState(stat);
    if (flag == 3)
        butList.find(ShowScanRICHList)->second->SetState(stat);
    if (flag == 4)
        butList.find(ShowRecRICHList)->second->SetState(stat);
    if (flag == 5)
        butList.find(ShowPIDEff)->second->SetState(stat);
    if (flag == 6)
        butList.find(GenScanRICHList)->second->SetState(stat);
    if (flag == 7)
        butList.find(GenRecRICHList)->second->SetState(stat);
    if (flag == 7)
        butList.find(GenPIDEff)->second->SetState(stat);
}
