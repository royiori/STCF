//___ 画图 & GUI接口 ___

#include "TPad.h"
#include "TView.h"
#include "TGraph2D.h"
#include "MyGuiMainFrame.h"
#include "MyGuiMainAction.h"
#include "MyGuiBeamTest.h"

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

    //for beamtest
    gMyGuiBeamTest = new MyGuiBeamTest(env);

    // pages
    nTabPage = 1;
    sTabPage[0] = TString("Beamtest");

    // page
    int ipage = 0;
    iButton[ipage] = SetDSPath;
    sButton[ipage].push_back("Set Data-Structure Path");
    sButton[ipage].push_back("Read binary data to Raw-root file");
    sButton[ipage].push_back("Check the Raw-Root file");
    sButton[ipage].push_back("Generate DST-root file");
    sButton[ipage].push_back("Check the DST-root file");
    sButton[ipage].push_back("Combine All Detectors' DST-root file");
    sButton[ipage].push_back("Load the Combined DST-root file");
    sButton[ipage].push_back("Slider");
    sButton[ipage].push_back("LR-L");
    sButton[ipage].push_back("LR-R");
    sButton[ipage].push_back("Zoom to default view");
    sButton[ipage].push_back("SEP");
    sButton[ipage].push_back("SEP");
    sButton[ipage].push_back("Analysis the Combined DST-root file");
    nButton[ipage] = (int)sButton[ipage].size();
};

MyGuiMainAction::~MyGuiMainAction()
{
    // Destructor.
    // Store env
    gMyGuiBeamTest->StoreEnv();
    env->SaveLevel(kEnvLocal);
}

//____________________________________________________________________
//
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

bool MyGuiMainAction::CheckIdelDatFile(TString cmdStr)
{
    if (cmdStr == "no")
        return false;

    FileStat_t fs;
    if (gSystem->GetPathInfo(gMyGuiBeamTest->GetDSPath(), fs))
    {
        cout << "#### Error: " << gMyGuiBeamTest->GetDSPath() << " doesn't exist." << endl;
        return false;
    }

    if (!gMyGuiBeamTest->GetDSPath().EndsWith("idle.dat"))
    {
        cout << "#### Error, please select the 'idle.dat' file." << endl;
        return false;
    }
    return true;
}

//______________________________________________________________________________
// basic button actions
void MyGuiMainAction::ExecButtonClick(Long_t bid, const char *cmdStr)
{
    //for beamtest
    if (bid == LoadBeamTextBuf)
        return DoLoadBeamTextBuf();

    if (bid == DrawBeamTestConfig)
        DoDrawBeamTestConfig();

    if (bid == SetDSPath)
        DoSetDSPath(cmdStr);
    if (bid == ReadRawData)
        DoReadRawData(cmdStr);
    if (bid == CheckRawRoot)
        DoCheckRawRoot(cmdStr);

    if (bid == GenDSTRoot)
        DoGenDSTRoot(cmdStr);
    if (bid == CheckDSTRoot)
        DoCheckDSTRoot(cmdStr);

    if (bid == CombineDSTRoot)
        DoCombineDSTRoot(cmdStr);

    if (bid == LoadDSTRoot)
        DoLoadDSTRoot(cmdStr);
    if (bid == SlidDSTHit)
        DoSlidDSTHit(cmdStr);
    if (bid == PrevDSTHit)
        DoPrevDSTHit();
    if (bid == NextDSTHit)
        DoNextDSTHit();
    if (bid == ZoomToDefault)
        DoZoomToDefault();

    if (bid == AnalysisCMDRoot)
        DoAnalysisCMDRoot();
}

//______________________________________________________________________________
// for beamtest data analysis
TString MyGuiMainAction::GenerateBEAMSettingsText()
{
    return gMyGuiBeamTest->GenerateSettingText();
}

void MyGuiMainAction::DoLoadBeamTextBuf()
{
    gMyGuiMainFrame->LoadText1(GenerateBEAMSettingsText());
}
void MyGuiMainAction::DoDrawBeamTestConfig()
{
    gMyGuiMainFrame->SwitchCanvas(0);
    gMyGuiMainFrame->ClearCanvas(0);
    gMyGuiBeamTest->DrawConfig();
    gMyGuiMainFrame->UpdateCanvas(0);
}

void MyGuiMainAction::DoSlidDSTHit(TString cmdStr)
{
    gMyGuiBeamTest->DrawCurtHit(cmdStr.Atoi());
    DoZoomToDefault();
}
void MyGuiMainAction::DoPrevDSTHit()
{
    gMyGuiBeamTest->DrawPrevHit();
    DoZoomToDefault();
}
void MyGuiMainAction::DoNextDSTHit()
{
    gMyGuiBeamTest->DrawNextHit();
    DoZoomToDefault();
}
void MyGuiMainAction::DoZoomToDefault()
{
    TView *view = gPad->GetView();
    if (view == NULL)
        return;
    view->RotateView(90, 0);
    view->ZoomIn();
    view->ZoomIn();
    view->ZoomIn();
    view->ZoomIn();
    view->ZoomIn();
    view->ZoomIn();
    gPad->Modified();
    gPad->Update();
}

//___1. 生成&检查 raw.root______________________
//
void MyGuiMainAction::DoSetDSPath(const char *fname)
{
    gMyGuiBeamTest->SetDSPath(fname);
    gMyGuiBeamTest->ReadSettingsText();
    DoLoadBeamTextBuf();
}

void MyGuiMainAction::DoReadRawData(TString cmdStr)
{
    if (!CheckIdelDatFile(cmdStr))
        return;

    gMyGuiBeamTest->ConvtBinaryToRawRoot();
}

void MyGuiMainAction::DoCheckRawRoot(TString cmdStr)
{
    int ID = RICH;
    if (cmdStr == "TrackerAGET")
        ID = TrackerAGET;

    { // 读入pedestal的分布并展示
        TH1::AddDirectory(kFALSE);
        TH1F *fNoiseMean;
        TH1F *fNoiseRMS;
        TGraph2D *fWaveAll[10];

        TFile fPedFile(gMyGuiBeamTest->GenPath(ID, PED));
        if (!fPedFile.IsOpen())
            return;
        cout << "--> Opening pedestal file: " << gMyGuiBeamTest->GenPath(ID, PED) << endl;

        fNoiseMean = (TH1F *)fPedFile.Get("fNoiseMean");
        fNoiseRMS = (TH1F *)fPedFile.Get("fNoiseRMS");
        fWaveAll[0] = (TGraph2D *)fPedFile.Get("gWaveAll0");
        fWaveAll[1] = (TGraph2D *)fPedFile.Get("gWaveAll1");
        fWaveAll[2] = (TGraph2D *)fPedFile.Get("gWaveAll2");
        fWaveAll[3] = (TGraph2D *)fPedFile.Get("gWaveAll3");

        int cid = 1;
        gMyGuiMainFrame->ClearAllCanvas();
        gMyGuiMainFrame->SwitchCanvas(cid);
        if (fNoiseMean != NULL)
            fNoiseMean->Draw();
        gMyGuiMainFrame->UpdateCanvas(cid++);

        gMyGuiMainFrame->SwitchCanvas(cid);
        if (fNoiseRMS != NULL)
            fNoiseRMS->Draw();
        gMyGuiMainFrame->UpdateCanvas(cid++);

        gMyGuiMainFrame->SwitchCanvas(cid);
        if (fWaveAll[0] != NULL)
            fWaveAll[0]->Draw("pcol");
        gMyGuiMainFrame->UpdateCanvas(cid++);

        gMyGuiMainFrame->SwitchCanvas(cid);
        if (fWaveAll[1] != NULL)
            fWaveAll[1]->Draw("pcol");
        gMyGuiMainFrame->UpdateCanvas(cid++);

        gMyGuiMainFrame->SwitchCanvas(cid);
        if (fWaveAll[2] != NULL)
            fWaveAll[2]->Draw("pcol");
        gMyGuiMainFrame->UpdateCanvas(cid++);

        gMyGuiMainFrame->SwitchCanvas(cid);
        if (fWaveAll[3] != NULL)
            fWaveAll[3]->Draw("pcol");
        gMyGuiMainFrame->UpdateCanvas(cid++);

        fPedFile.Close();
    }

    env->SetValue("checkRawRoot", gMyGuiBeamTest->GenPath(ID, RAW));
    env->SetValue("checkPedRoot", gMyGuiBeamTest->GenPath(ID, PED));
    env->SetValue("checkRawType", ID);
    env->Save();
    gSystem->Exec(Form("(root %s/data/check1RawRoot.C)", realPath.Data()));
}

//____2. 生成&检查 dst.root_____________________
// convert beam-test raw root file to a new data structure
void MyGuiMainAction::DoGenDSTRoot(TString cmdStr)
{
    if (!CheckIdelDatFile(cmdStr))
        return;

    gMyGuiBeamTest->ConvtRawToDstRoot();
}

void MyGuiMainAction::DoCheckDSTRoot(TString cmdStr)
{
    if (!CheckIdelDatFile(cmdStr))
        return
        
    env->SetValue("checkDSTRoot", gMyGuiBeamTest->GenPath(RICH, DST));
    env->Save();
    gSystem->Exec(Form("(root %s/data/check2DSTRoot.C)", realPath.Data()));
}

//____3. 生成&检查 combin.root_____________________
// Combine Results
void MyGuiMainAction::DoCombineDSTRoot(TString cmdStr)
{
    if (!CheckIdelDatFile(cmdStr))
        return;

    gMyGuiBeamTest->CombineDSTRoot(gMyGuiBeamTest->GetDSPath().Data());
}

void MyGuiMainAction::DoLoadDSTRoot(TString cmdStr)
{
    if (cmdStr == "no")
        return;

    DoDrawBeamTestConfig();

    if (!gMyGuiBeamTest->GetDSPath().EndsWith("idle.dat"))
    {
        cout << "#### Error, please select the 'idle.dat' file." << endl;
        return;
    }

    gMyGuiMainFrame->SetWindowName(gMyGuiBeamTest->GetDSPath().ReplaceAll("idle.dat", ""));
    gMyGuiMainFrame->EnableShowGeoSliderButton();

    //读入文件
    int cid = 0;
    gMyGuiMainFrame->SwitchCanvas(cid);
    gMyGuiMainFrame->ClearCanvas(cid);
    gMyGuiBeamTest->LoadDSTRoot(gMyGuiBeamTest->GetDSPath().Data());
    DoZoomToDefault();
    gMyGuiMainFrame->UpdateCanvas(cid++);
    gMyGuiMainFrame->SetGeoSliderRange(gMyGuiBeamTest->GetDSTEntries()); //根据文件内的事例数，设置slider范围

    //画分布图-RICH
    gStyle->SetPalette(57);
    gMyGuiMainFrame->SwitchCanvas(cid);
    gMyGuiMainFrame->ClearCanvas(cid);
    TCanvas *can1 = gMyGuiMainFrame->GetCanvas(cid);
    can1->Divide(3, 2);
    can1->cd(1);
    gMyGuiBeamTest->GetRICHDet()[0]->GetAllHit()->Draw("colz");
    can1->cd(2);
    gMyGuiBeamTest->GetRICHDet()[0]->GetCharge()->Draw();
    can1->cd(4);
    gMyGuiBeamTest->GetRICHDet()[0]->GetAllHit2()->Draw("colz");
    can1->cd(3);
    gMyGuiBeamTest->GetRICHDet()[0]->GetAllHit()->ProjectionX()->Draw();
    can1->cd(6);
    gMyGuiBeamTest->GetRICHDet()[0]->GetAllHit()->ProjectionY()->Draw();
    gMyGuiMainFrame->UpdateCanvas(cid++);

    //Track02
    gMyGuiMainFrame->SwitchCanvas(cid);
    gMyGuiMainFrame->ClearCanvas(cid);
    TCanvas *can2 = gMyGuiMainFrame->GetCanvas(cid);
    can2->Divide(3, 2);
    can2->cd(1);
    gMyGuiBeamTest->GetTrackAGET()[0]->GetAllHit()->Draw("colz");
    can2->cd(4);
    gMyGuiBeamTest->GetTrackAGET()[0]->GetAllHit2()->Draw("colz");
    can2->cd(2);
    gMyGuiBeamTest->GetTrackAGET()[0]->GetCharge()->Draw();
    can2->cd(3);
    gMyGuiBeamTest->GetTrackAGET()[0]->GetHitXY()[0]->Draw();
    can2->cd(6);
    gMyGuiBeamTest->GetTrackAGET()[0]->GetHitXY()[1]->Draw();
    gMyGuiMainFrame->UpdateCanvas(cid++);

    //Track03
    gMyGuiMainFrame->SwitchCanvas(cid);
    gMyGuiMainFrame->ClearCanvas(cid);
    TCanvas *can3 = gMyGuiMainFrame->GetCanvas(cid);
    can3->Divide(3, 2);
    can3->cd(1);
    gMyGuiBeamTest->GetTrackAGET()[1]->GetAllHit()->Draw("colz");
    can3->cd(4);
    gMyGuiBeamTest->GetTrackAGET()[1]->GetAllHit2()->Draw("colz");
    can3->cd(2);
    gMyGuiBeamTest->GetTrackAGET()[1]->GetCharge()->Draw();
    can3->cd(3);
    gMyGuiBeamTest->GetTrackAGET()[1]->GetHitXY()[0]->Draw();
    can3->cd(6);
    gMyGuiBeamTest->GetTrackAGET()[1]->GetHitXY()[1]->Draw();
    gMyGuiMainFrame->UpdateCanvas(cid++);

    //Track06
    gMyGuiMainFrame->SwitchCanvas(cid);
    gMyGuiMainFrame->ClearCanvas(cid);
    TCanvas *can4 = gMyGuiMainFrame->GetCanvas(cid);
    can4->Divide(3, 2);
    can4->cd(1);
    gMyGuiBeamTest->GetTrackAGET()[2]->GetAllHit()->Draw("colz");
    can4->cd(4);
    gMyGuiBeamTest->GetTrackAGET()[2]->GetAllHit2()->Draw("colz");
    can4->cd(2);
    gMyGuiBeamTest->GetTrackAGET()[2]->GetCharge()->Draw();
    can4->cd(3);
    gMyGuiBeamTest->GetTrackAGET()[2]->GetHitXY()[0]->Draw();
    can4->cd(6);
    gMyGuiBeamTest->GetTrackAGET()[2]->GetHitXY()[1]->Draw();
    gMyGuiMainFrame->UpdateCanvas(cid++);

    /*
    TH2F *ftmp = gMyGuiBeamTest->GetTrackAGET()[2]->GetAllHit();
    double meanX = ftmp->GetMean(1);
    double meanY = ftmp->GetMean(2);
    cout<<"Get events from ("<<meanX<<", "<<meanY<<")"<<endl;
    double entry = 0;
    for(int i=0; i<128; i++)
        for(int j=0; j<128; j++)
        {
            if(fabs(ftmp->GetXaxis()->GetBinCenter(i)-meanX)<5/0.4 && 
               fabs(ftmp->GetYaxis()->GetBinCenter(j)-meanY)<5/0.4)
                entry+=ftmp->GetBinContent(i, j);
        }
    cout<<"In area: "<<entry<<" / Total entry: "<<ftmp->GetEntries()<<endl;
    */

    gMyGuiMainFrame->SwitchCanvas(0);
}

void MyGuiMainAction::DoAnalysisCMDRoot()
{
    gMyGuiBeamTest->AnalysisDSTRoot(gMyGuiBeamTest->GetDSPath().Data());
}
