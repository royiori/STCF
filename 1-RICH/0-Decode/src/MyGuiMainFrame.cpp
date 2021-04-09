#include "MyGuiMainFrame.h"

MyGuiMainFrame *gMyGuiMainFrame = (MyGuiMainFrame *)0;

using namespace std;

const char *filetypes[] = {"ROOT files", "*.root",
                           "Dat files", "*.dat",
                           "All files", "*",
                           0, 0};

const char *filetypes2[] = {"Dat files", "*.dat",
                            "All files", "*",
                            0, 0};
//______________________________________________________________________________
//
MyGuiMainFrame::MyGuiMainFrame(const TGWindow *p, int w, int h) : TGMainFrame(p, w, h)
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
            fSettingText1 = new TGTextEdit(fTabPage, 400, 580);
            fSettingText1->LoadBuffer(gMyGuiMainAction->GenerateSettingsText());
            fTabPage->AddFrame(fSettingText1, LayoutXY);

            TGHorizontalFrame *fHFrame11 = new TGHorizontalFrame(fTabPage, 10, 10, kHorizontalFrame);
            fTabPage->AddFrame(fHFrame11, LayoutX);

            //#GUI. Define setting-text buttons
            TGTextButton *fButtonTmp3 = new TGTextButton(fHFrame11, "Load BeamTest Set-text", 1);
            fButtonTmp3->Associate(this); 
            fHFrame11->AddFrame(fButtonTmp3, LayoutX);

            TGTextButton *fButtonTmp4 = new TGTextButton(fHFrame11, "Draw BeamTest Configure", 2);
            fButtonTmp4->Associate(this);
            fHFrame11->AddFrame(fButtonTmp4, LayoutX);
        }
    }

    //2. Buttons Tab
    {
        TGTab *fTab0 = new TGTab(fVFrame2, 400, 580);
        fTab0->SetTab(0);
        fTab0->Resize(fTab0->GetDefaultSize());
        fVFrame2->AddFrame(fTab0, LayoutXY);

        // 2.1 page0
        vector<ButtonAct> butAct = gMyGuiMainAction->GetButAct();

        int id = 100;
        for (int i = 0; i < (int)butAct.size(); i++)
        {
            if (butAct[i].iPG != -1)
                continue;
            if (butAct[i].title == "SEP")
                continue;
            butAct[i].setID(id);
            cout << "--> Setting " << id++ << " to " << butAct[i].title << endl;
        }

        TGCompositeFrame *fTabPage;
        for (int i = 0; i < (int)butAct.size(); i++)
        {
            if (butAct[i].iPG != -1)
            {
                fTabPage = fTab0->AddTab(butAct[i].title);
                fTabPage->SetLayoutManager(new TGVerticalLayout(fTabPage));
                continue;
            }

            if (butAct[i].title == "SEP")
            {
                TGHorizontal3DLine *fHorizontal3DLine = new TGHorizontal3DLine(fTabPage, 102, 8);
                fTabPage->AddFrame(fHorizontal3DLine, LayoutX);
            }
            else if (butAct[i].title == "Slider")
            {
                TGHSlider *fHSlider = new TGHSlider(fTabPage, 200, kSlider1 | kScaleBoth, butAct[i].uID, kHorizontalFrame);
                fHSlider->Associate(this);
                fHSlider->SetRange(0, 1024);
                fHSlider->SetPosition(512);
                fTabPage->AddFrame(fHSlider, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                sldList.insert(pair<int, TGSlider *>(butAct[i].uID, fHSlider));

                TGHorizontalFrame *fSliderFrame = new TGHorizontalFrame(fTabPage, 220, 100, kHorizontalFrame);
                fTabPage->AddFrame(fSliderFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                {
                    i++;
                    TGTextButton *fButtonTmp1 = new TGTextButton(fSliderFrame, "Pre Channel <", butAct[i].uID);
                    fButtonTmp1->Associate(this);
                    fSliderFrame->AddFrame(fButtonTmp1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                    butList.insert(pair<int, TGTextButton *>(butAct[i].uID, fButtonTmp1));

                    i++;
                    TGTextButton *fButtonTmp2 = new TGTextButton(fSliderFrame, "Next Channel >", butAct[i].uID);
                    fButtonTmp2->Associate(this);
                    fSliderFrame->AddFrame(fButtonTmp2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                    butList.insert(pair<int, TGTextButton *>(butAct[i].uID, fButtonTmp2));
                }
            }
            else
            {
                TGTextButton *fButtonTmp = new TGTextButton(fTabPage, butAct[i].title, butAct[i].uID);
                fButtonTmp->Associate(this);
                fTabPage->AddFrame(fButtonTmp, LayoutX);
                if (butAct[i].title.BeginsWith("Analysis"))
                    fButtonTmp->ChangeBackground(green);
                if (butAct[i].title.BeginsWith("Set"))
                    fButtonTmp->ChangeBackground(yellow);
                if (butAct[i].title.BeginsWith("Draw"))
                    fButtonTmp->ChangeBackground(yellow);
                butList.insert(pair<int, TGTextButton *>(butAct[i].uID, fButtonTmp));
            }
        }
    }

    //3. Results tab
    {
        fCTab = new TGTab(fVFrame3, 900, 580);

        for (int i = 0; i < NResPAGE; i++)
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
        pt->AddText("STCF RICH beamtest software package");
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
        fCommand->Connect("ReturnPressed()", "MyGuiMainFrame", this, "HandleCommand()");
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

MyGuiMainFrame::~MyGuiMainFrame()
{
    // Destructor.
    CloseWindow();
}

void MyGuiMainFrame::CloseWindow()
{
    // Destructor.
    delete gMyGuiMainAction;
    DeleteWindow();
    gApplication->Terminate(0);
}

//______________________________________________________________________________
//

Bool_t MyGuiMainFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
    switch (GET_MSG(msg))
    {
    case kC_HSLIDER:
        switch (GET_SUBMSG(msg))
        {
        case kSL_RELEASE: //slider
            gMyGuiMainAction->ExecButtonClick(parm1, TString(Form("%d", sldList.find(parm1)->second->GetPosition())));
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

            ///#GUI. Open file dialog
            //open file dialog
            if (parm1 == 100 || parm1 == 102) // dat file dialog
            {
                fi.fFileTypes = filetypes2;
                new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
                if (fi.fFilename == NULL)
                    return kTRUE;
                cmdStr = fi.fFilename;
            }

            if (parm1 == 101 || parm1 == 103) //root file dialog
            {
                fi.fFileTypes = filetypes;
                new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
                if (fi.fFilename == NULL)
                    return kTRUE;
                cmdStr = fi.fFilename;
            }

            gMyGuiMainAction->ExecButtonClick(parm1, cmdStr.Data());
            break;
        }
        break;
    }
    parm2 = 0.;
    return kTRUE;
}

//______________________________________________________________________________
//

void MyGuiMainFrame::HandleCommand()
{
    const char *string = fCommandBuf->GetString();
    if (strlen(string) > 1)
    {
        gROOT->ProcessLine(string);
        fComboCmd->InsertEntry(string, 0, -1);
        fCommand->Clear();
    }
}

void MyGuiMainFrame::SetSliderRange(int rmax, ShowButtonAction flag)
{
    /*
    if (flag == ShowButtonAction::EnableShowGeoSliderButton)
    {
        sldList.find(SlidDSTHit)->second->SetRange(0, rmax);
        sldList.find(SlidDSTHit)->second->SetPosition(0);
    }
    */
}

void MyGuiMainFrame::ToggleButList(bool on, ShowButtonAction flag)
{
    EButtonState stat = (on) ? kButtonUp : kButtonDisabled;

    map<int, TGTextButton *>::iterator mIter;

    //#GUI. Set button status
    if (flag == InitStatus) //初始状态
    {
        sldList.find(104)->second->SetEnabled(on);
        butList.find(105)->second->SetState(stat);
        butList.find(106)->second->SetState(stat);
        butList.find(107)->second->SetState(stat);
    }
    else
    {
        sldList.find(104)->second->SetEnabled(on);
        butList.find(105)->second->SetState(stat);
        butList.find(106)->second->SetState(stat);
        butList.find(107)->second->SetState(stat);
    }
}

