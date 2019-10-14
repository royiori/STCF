#include "MyMainFrameGui.h"

MyMainFrameGui *gMyMainFrameGui = (MyMainFrameGui *)0;

using namespace std;

const char *filetypes[] = {"ROOT files", "*.root",
                           "Dat files", "*.dat",
                           "All files", "*",
                           0, 0};

const char *filetypes2[] = {"Dat files", "*.dat",
                            "All files", "*",
                            0, 0};
const int NPAGE = 11;
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

    //1. Setting Text Tab
    {
        TGTab *fTab1 = new TGTab(fVFrame1, 300, 580);
        fTab1->SetTab(0);
        fTab1->Resize(fTab1->GetDefaultSize());
        fVFrame1->AddFrame(fTab1, LayoutXY);

        {
            TGCompositeFrame *fTabPage = fTab1->AddTab("Beamtest");
            fSettingText2 = new TGTextEdit(fTabPage, 400, 580);
            fSettingText2->LoadBuffer(gMyGuiActionClass->GenerateBEAMSettingsText());
            fTabPage->AddFrame(fSettingText2, LayoutXY);

            TGHorizontalFrame *fHFrame11 = new TGHorizontalFrame(fTabPage, 10, 10, kHorizontalFrame);
            fTabPage->AddFrame(fHFrame11, LayoutX);

            TGTextButton *fButtonTmp3 = new TGTextButton(fHFrame11, "Load BeamTest Set-text", LoadBeamTextBuf);
            fButtonTmp3->Associate(this);
            fHFrame11->AddFrame(fButtonTmp3, LayoutX);

            TGTextButton *fButtonTmp4 = new TGTextButton(fHFrame11, "Draw BeamTest Configure", DrawBeamTestConfig);
            fButtonTmp4->Associate(this);
            fHFrame11->AddFrame(fButtonTmp4, LayoutX);
        }

        {
            TGCompositeFrame *fTabPage = fTab1->AddTab("RICH-sim");
            fTabPage->SetLayoutManager(new TGVerticalLayout(fTabPage));

            fSettingText1 = new TGTextEdit(fTabPage, 400, 580);
            fSettingText1->LoadBuffer(gMyGuiActionClass->GenerateSIMSettingsText());
            fSettingText1->SetFont(myfont->GetFontStruct());
            fTabPage->AddFrame(fSettingText1, LayoutXY);

            TGHorizontalFrame *fHFrame10 = new TGHorizontalFrame(fTabPage, 10, 10, kHorizontalFrame);
            fTabPage->AddFrame(fHFrame10, LayoutX);

            TGTextButton *fButtonTmp1 = new TGTextButton(fHFrame10, "Load SIM Set-text", LoadTextBuf);
            fButtonTmp1->Associate(this);
            fHFrame10->AddFrame(fButtonTmp1, LayoutX);

            TGTextButton *fButtonTmp2 = new TGTextButton(fHFrame10, "Draw SIM Configure", DrawConfig);
            fButtonTmp2->Associate(this);
            fHFrame10->AddFrame(fButtonTmp2, LayoutX);
        }
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
                else if (gMyGuiActionClass->GetPageButtonName(i, j) == "Slider")
                {
                    TGHSlider *fHSlider = new TGHSlider(fTabPage, 200, kSlider1 | kScaleBoth, gMyGuiActionClass->GetIndexButton(i) + j, kHorizontalFrame);
                    fHSlider->Associate(this);
                    fHSlider->SetRange(0, 1024);
                    fHSlider->SetPosition(512);
                    fTabPage->AddFrame(fHSlider, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                    sldList.insert(pair<int, TGSlider *>(gMyGuiActionClass->GetIndexButton(i) + j, fHSlider));

                    TGHorizontalFrame *fSliderFrame = new TGHorizontalFrame(fTabPage, 220, 100, kHorizontalFrame);
                    fTabPage->AddFrame(fSliderFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                    {

                        TGTextButton *fButtonTmp1 = new TGTextButton(fSliderFrame, "Pre Channel <", gMyGuiActionClass->GetIndexButton(i) + j + 1);
                        fButtonTmp1->Associate(this);
                        fSliderFrame->AddFrame(fButtonTmp1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                        butList.insert(pair<int, TGTextButton *>(gMyGuiActionClass->GetIndexButton(i) + j + 1, fButtonTmp1));

                        TGTextButton *fButtonTmp2 = new TGTextButton(fSliderFrame, "Next Channel >", gMyGuiActionClass->GetIndexButton(i) + j + 2);
                        fButtonTmp2->Associate(this);
                        fSliderFrame->AddFrame(fButtonTmp2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                        butList.insert(pair<int, TGTextButton *>(gMyGuiActionClass->GetIndexButton(i) + j + 2, fButtonTmp2));
                    }
                    j = j + 2;
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
        fCTab = new TGTab(fVFrame3, 900, 580);
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
        TRootEmbeddedCanvas *fEmbeddedCanvas = new TRootEmbeddedCanvas(0, fTabPage, 900, 360);
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

    InitShowButtonList();

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
    case kC_HSLIDER:
        switch (GET_SUBMSG(msg))
        {
        case kSL_RELEASE: //slider
            gMyGuiActionClass->ExecButtonClick(parm1, TString(Form("%d", sldList.find(parm1)->second->GetPosition())));
            break;
        }
        break;
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

            if (parm1 == SaveDetFile || parm1 == LoadDetFile)
            {
                new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
                if (fi.fFilename == NULL)
                    return kTRUE;
                cmdStr = fi.fFilename;
                if (!cmdStr.EndsWith("root"))
                    cmdStr += ".root";
            }
            if (parm1 == SetDSPath)
            {
                fi.fFileTypes = filetypes2;
                new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
                if (fi.fFilename == NULL)
                    return kTRUE;
                cmdStr = fi.fFilename;
            }

            if (parm1 == GenSpecRICH || parm1 == GenMulParRICH || parm1 == GenScanRICHList || parm1 == GenRecRICHList || parm1 == GenPIDEff || parm1 == ReadRawData || parm1 == GenDSTRoot || parm1 == CombineDSTRoot || parm1 == LoadDSTRoot)
            {
                int retval;
                TString msg("Do you want RE-generate the histograms?");
                switch (parm1)
                {
                case ReadRawData:
                    msg = TString("Do you want read the raw data \"*.dat\" from the DS-Path?");
                    break;
                case GenDSTRoot:
                    msg = TString("Do you want read the raw root \"*-raw.root\" and generate the DST-root \"*-dst.root\" from the DS-Path?");
                    break;
                case CombineDSTRoot:
                    msg = TString("Do you want read the dst root data \"*-dst.root\" and combine to \"*-combined-dst.root\" from the DS-Path?");
                    break;
                case LoadDSTRoot:
                    msg = TString("Do you want load the combined root \"*-combine-dst.root\" from the DS-Path?");
                    break;
                }
                new TGMsgBox(gClient->GetRoot(), this, "Message-RICH", msg, kMBIconExclamation, kMBYes | kMBNo, &retval);
                cmdStr = (retval == 1) ? "yes" : "no";
            }

            if (parm1 == CheckRawRoot)
            {
                int retval;
                TString msg("Which one do you want to check? RICH ? TrackAGET ? TrackVMM");
                new TGMsgBox(gClient->GetRoot(), this, "Message-RICH", msg, kMBIconExclamation, kMBYes | kMBNo | kMBNoAll, &retval);
                if (retval == 1)
                    cmdStr = "RICH";
                if (retval == 2)
                    cmdStr = "TrackerAGET";
                if (retval == 1024)
                    cmdStr = "TrackerVMM";
            }

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

void MyMainFrameGui::SetSliderRange(int rmax, ShowButtonAction flag)
{
    if (flag == ShowButtonAction::EnableShowGeoSliderButton)
    {
        sldList.find(SlidDSTHit)->second->SetRange(0, rmax);
        sldList.find(SlidDSTHit)->second->SetPosition(0);
    }
}

void MyMainFrameGui::ToggleButList(bool on, ShowButtonAction flag)
{
    EButtonState stat = (on) ? kButtonUp : kButtonDisabled;

    map<int, TGTextButton *>::iterator mIter;

    if (flag == InitStatus) //初始状态
    {
        butList.find(ShowSpecRICH)->second->SetState(stat);
        butList.find(ShowMulParRICH)->second->SetState(stat);
        butList.find(ShowScanRICHList)->second->SetState(stat);
        butList.find(ShowRecRICHList)->second->SetState(stat);
        butList.find(ShowPIDEff)->second->SetState(stat);
        butList.find(GenScanRICHList)->second->SetState(stat);
        butList.find(GenRecRICHList)->second->SetState(stat);
        butList.find(GenPIDEff)->second->SetState(stat);
        butList.find(PrevDSTHit)->second->SetState(stat);
        butList.find(NextDSTHit)->second->SetState(stat);
        sldList.find(SlidDSTHit)->second->SetEnabled(on);
    }

    if (flag == ShowButtonAction::EnableShowSpecButton)
        butList.find(ShowSpecRICH)->second->SetState(stat);
    if (flag == ShowButtonAction::EnableShowMulParButton)
        butList.find(ShowMulParRICH)->second->SetState(stat);
    if (flag == ShowButtonAction::EnableShowScanRICHButton)
        butList.find(ShowScanRICHList)->second->SetState(stat);
    if (flag == ShowButtonAction::EnableShowRecRICHButton)
        butList.find(ShowRecRICHList)->second->SetState(stat);
    if (flag == ShowButtonAction::EnableShowPIDEffButton)
        butList.find(ShowPIDEff)->second->SetState(stat);

    if (flag == ShowButtonAction::EnableGenScanRICHButton)
        butList.find(GenScanRICHList)->second->SetState(stat);
    if (flag == ShowButtonAction::EnableGenRecRICHButton)
        butList.find(GenRecRICHList)->second->SetState(stat);
    if (flag == ShowButtonAction::EnableGenPIDEffButton)
        butList.find(GenPIDEff)->second->SetState(stat);

    if (flag == ShowButtonAction::EnableShowGeoSliderButton)
    {
        butList.find(PrevDSTHit)->second->SetState(stat);
        butList.find(NextDSTHit)->second->SetState(stat);
        sldList.find(SlidDSTHit)->second->SetEnabled(on);
    }
}
