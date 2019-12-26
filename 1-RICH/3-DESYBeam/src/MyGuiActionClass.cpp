#include "TPad.h"
#include "TView.h"
#include "MyMainFrameGui.h"
#include "MyGuiActionClass.h"
#include "MyBeamTest.h"
#include "CombineData.h"

MyGuiActionClass *gMyGuiActionClass = (MyGuiActionClass *)0;

//______________________________________________________________________________
//
MyGuiActionClass::MyGuiActionClass(const char *realpath, int flag)
{
    BatGuiFlag = flag;
    env = new TEnv(gSystem->WorkingDirectory() + TString("/.env"));
    env->SaveLevel(kEnvLocal);

    gMyStyle = new MyStyle();
    gMyStyle->SetDrawOption("c");
    gMyStyle->SetColorPattern(MATHEMATIC_STYLE);

    realPath = TString(realpath);

    //for beamtest
    gMyBeamTest = new MyBeamTest(env);

    // pages
    nTabPage = 4;
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

MyGuiActionClass::~MyGuiActionClass()
{
    // Destructor.
    // Store env
    ReadBEAMSettingsText();
    gMyBeamTest->StoreEnv();
    env->SaveLevel(kEnvLocal);
}

vector<TString> MyGuiActionClass::ReadContent(TString LINE)
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
void MyGuiActionClass::ExecButtonClick(Long_t bid, const char *cmdStr)
{
    //for beamtest
    if (bid == LoadBeamTextBuf)
        return DoLoadBeamTextBuf();

    ReadBEAMSettingsText();
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
void MyGuiActionClass::DoLoadBeamTextBuf()
{
    gMyMainFrameGui->LoadText2(GenerateBEAMSettingsText());
}

TString MyGuiActionClass::GenerateBEAMSettingsText()
{
    return gMyBeamTest->GenerateSettingText();
}

int MyGuiActionClass::ReadBEAMSettingsText()
{
    return gMyBeamTest->ReadSettingsText(gMyMainFrameGui->GetText2());
}

void MyGuiActionClass::DoDrawBeamTestConfig()
{
    gMyMainFrameGui->SwitchCanvas(0);
    gMyMainFrameGui->ClearCanvas(0);
    gMyBeamTest->DrawConfig();
    gMyMainFrameGui->UpdateCanvas(0);
}

void MyGuiActionClass::DoZoomToDefault()
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

//_________________________
// read beam-test dat file and save to -raw.root
void MyGuiActionClass::DoSetDSPath(const char *fname)
{
    gMyBeamTest->SetDSPath(fname);
    DoLoadBeamTextBuf();
}

void MyGuiActionClass::DoReadRawData(TString cmdStr)
{
    if (!CheckIdelDatFile(cmdStr))
        return;

    TString fileDir = gMyBeamTest->GetDSPath();
    fileDir.ReplaceAll("idle.dat", "");

    TString fileDir1(fileDir + "/RICH/");
    TString fileDir2(fileDir + "/TrackAGET/");
    TString fileDir3(fileDir + "/TrackVMM/");
    TString fileDir4(fileDir + "/Combine/");

    vector<TString> datList1;
    vector<TString> datList2;
    vector<TString> datList3;

    GetFileList(fileDir1, ".dat", datList1);
    GetFileList(fileDir2, ".dat", datList2);
    GetFileList(fileDir3, ".bin", datList3);

    if (datList1.size() > 0)
    {
        ReadRICHData2Root(datList1, gMyBeamTest->GenPath(RICH, RAW), 0);
        GenerateRICHPed(gMyBeamTest->GenPath(RICH, RAW), gMyBeamTest->GenPath(RICH, PED), gMyBeamTest->GetRICHDet());
    }

    if (datList2.size() > 0)
    {
        ReadTrackAGTData2Root(datList2, gMyBeamTest->GenPath(TrackerAGET, RAW), 0);
        GenerateAGETPed(gMyBeamTest->GenPath(TrackerAGET, RAW), gMyBeamTest->GenPath(TrackerAGET, PED), gMyBeamTest->GetTrackAGET());
    }

    if (datList3.size() > 0)
        ReadTrackVMMData2Root(datList3, gMyBeamTest->GenPath(TrackerVMM, RAW), 1);
}

void MyGuiActionClass::DoCheckRawRoot(TString cmdStr)
{
    int ID = RICH;
    if (cmdStr == "TrackerAGET")
        ID = TrackerAGET;
    if (cmdStr == "TrackerVMM")
        ID = TrackerVMM;

    if (cmdStr != "TrackerVMM")
    {
        TH1::AddDirectory(kFALSE);
        TH1F *fNoiseMean;
        TH1F *fNoiseRMS;
        TGraph2D *fWaveAll[10];

        TFile fPedFile(gMyBeamTest->GenPath(ID, PED));
        if (!fPedFile.IsOpen())
            return;
        cout << "--> Opening pedestal file: " << gMyBeamTest->GenPath(ID, PED) << endl;

        fNoiseMean = (TH1F *)fPedFile.Get("fNoiseMean");
        fNoiseRMS = (TH1F *)fPedFile.Get("fNoiseRMS");
        fWaveAll[0] = (TGraph2D *)fPedFile.Get("gWaveAll0");
        fWaveAll[1] = (TGraph2D *)fPedFile.Get("gWaveAll1");
        fWaveAll[2] = (TGraph2D *)fPedFile.Get("gWaveAll2");
        fWaveAll[3] = (TGraph2D *)fPedFile.Get("gWaveAll3");

        int cid = 1;
        gMyMainFrameGui->ClearAllCanvas();
        gMyMainFrameGui->SwitchCanvas(cid);
        if (fNoiseMean != NULL)
            fNoiseMean->Draw();
        gMyMainFrameGui->UpdateCanvas(cid++);

        gMyMainFrameGui->SwitchCanvas(cid);
        if (fNoiseRMS != NULL)
            fNoiseRMS->Draw();
        gMyMainFrameGui->UpdateCanvas(cid++);

        gMyMainFrameGui->SwitchCanvas(cid);
        if (fWaveAll[0] != NULL)
            fWaveAll[0]->Draw("pcol");
        gMyMainFrameGui->UpdateCanvas(cid++);

        gMyMainFrameGui->SwitchCanvas(cid);
        if (fWaveAll[1] != NULL)
            fWaveAll[1]->Draw("pcol");
        gMyMainFrameGui->UpdateCanvas(cid++);

        gMyMainFrameGui->SwitchCanvas(cid);
        if (fWaveAll[2] != NULL)
            fWaveAll[2]->Draw("pcol");
        gMyMainFrameGui->UpdateCanvas(cid++);

        gMyMainFrameGui->SwitchCanvas(cid);
        if (fWaveAll[3] != NULL)
            fWaveAll[3]->Draw("pcol");
        gMyMainFrameGui->UpdateCanvas(cid++);

        fPedFile.Close();
    }

    env->SetValue("checkRawRoot", gMyBeamTest->GenPath(ID, RAW));
    env->SetValue("checkPedRoot", gMyBeamTest->GenPath(ID, PED));
    env->SetValue("checkRawType", ID);
    env->Save();
    gSystem->Exec(Form("(root %s/data/checkRawRoot.C)", realPath.Data()));
}

//_________________________
// convert beam-test raw root file to a new data structure
void MyGuiActionClass::DoGenDSTRoot(TString cmdStr)
{
    if (!CheckIdelDatFile(cmdStr))
        return;

    bool SaveWaveFlag = kFALSE;
    gMyBeamTest->ConvtRICHRoot(gMyBeamTest->GetDSPath().Data(), SaveWaveFlag);
    gMyBeamTest->ConvtTrackAGTRoot(gMyBeamTest->GetDSPath().Data(), SaveWaveFlag);
    gMyBeamTest->ConvtTrackVMMRoot(gMyBeamTest->GetDSPath().Data());
}

void MyGuiActionClass::DoCheckDSTRoot(TString cmdStr)
{
    //显示dst数据的噪声在实际位置上的分布
    if (cmdStr != "TrackerVMM")
    {
        gMyBeamTest->ReadRICHPed(gMyBeamTest->GenPath(RICH, PED));
        gMyBeamTest->ReadTrackAGTPed(gMyBeamTest->GenPath(TrackerAGET, PED));
    }

    env->SetValue("checkDSTRoot", gMyBeamTest->GenPath(RICH, DST));
    env->Save();
    gSystem->Exec(Form("(root %s/data/checkDSTRoot.C)", realPath.Data()));
}

//_________________________
// Combine Results
void MyGuiActionClass::DoCombineDSTRoot(TString cmdStr)
{
    if (!CheckIdelDatFile(cmdStr))
        return;

    gMyBeamTest->CombineDSTRoot(gMyBeamTest->GetDSPath().Data());
}

//_________________________
// Load Results
void MyGuiActionClass::DoLoadDSTRoot(TString cmdStr)
{
    if (cmdStr == "no")
        return;

    DoDrawBeamTestConfig();

    if (!gMyBeamTest->GetDSPath().EndsWith("idle.dat"))
    {
        cout << "#### Error, please select the 'idle.dat' file." << endl;
        return;
    }

    gMyMainFrameGui->SetWindowName(gMyBeamTest->GetDSPath().ReplaceAll("idle.dat", ""));
    gMyMainFrameGui->EnableShowGeoSliderButton();

    //读入文件
    int cid = 0;
    gMyMainFrameGui->SwitchCanvas(cid);
    gMyMainFrameGui->ClearCanvas(cid);
    gMyBeamTest->LoadDSTRoot(gMyBeamTest->GetDSPath().Data());
    DoZoomToDefault();
    gMyMainFrameGui->UpdateCanvas(cid++);
    gMyMainFrameGui->SetGeoSliderRange(gMyBeamTest->GetDSTEntries()); //根据文件内的事例数，设置slider范围

    //画分布图-RICH
    gStyle->SetPalette(57);
    gMyMainFrameGui->SwitchCanvas(cid);
    gMyMainFrameGui->ClearCanvas(cid);
    TCanvas *can1 = gMyMainFrameGui->GetCanvas(cid);
    can1->Divide(3, 2);
    can1->cd(1);
    gMyBeamTest->GetRICHDet()[0]->GetAllHit()->Draw("colz");
    can1->cd(2);
    gMyBeamTest->GetRICHDet()[0]->GetCharge()->Draw();
    can1->cd(4);
    gMyBeamTest->GetRICHDet()[0]->GetAllHit2()->Draw("colz");
    can1->cd(3);
    gMyBeamTest->GetRICHDet()[0]->GetAllHit()->ProjectionX()->Draw();
    can1->cd(6);
    gMyBeamTest->GetRICHDet()[0]->GetAllHit()->ProjectionY()->Draw();
    gMyMainFrameGui->UpdateCanvas(cid++);

    //Track02
    gMyMainFrameGui->SwitchCanvas(cid);
    gMyMainFrameGui->ClearCanvas(cid);
    TCanvas *can2 = gMyMainFrameGui->GetCanvas(cid);
    can2->Divide(3, 2);
    can2->cd(1);
    gMyBeamTest->GetTrackAGET()[0]->GetAllHit()->Draw("colz");
    can2->cd(4);
    gMyBeamTest->GetTrackAGET()[0]->GetAllHit2()->Draw("colz");
    can2->cd(2);
    gMyBeamTest->GetTrackAGET()[0]->GetCharge()->Draw();
    can2->cd(3);
    gMyBeamTest->GetTrackAGET()[0]->GetHitXY()[0]->Draw();
    can2->cd(6);
    gMyBeamTest->GetTrackAGET()[0]->GetHitXY()[1]->Draw();
    gMyMainFrameGui->UpdateCanvas(cid++);

    //Track03
    gMyMainFrameGui->SwitchCanvas(cid);
    gMyMainFrameGui->ClearCanvas(cid);
    TCanvas *can3 = gMyMainFrameGui->GetCanvas(cid);
    can3->Divide(3, 2);
    can3->cd(1);
    gMyBeamTest->GetTrackAGET()[1]->GetAllHit()->Draw("colz");
    can3->cd(4);
    gMyBeamTest->GetTrackAGET()[1]->GetAllHit2()->Draw("colz");
    can3->cd(2);
    gMyBeamTest->GetTrackAGET()[1]->GetCharge()->Draw();
    can3->cd(3);
    gMyBeamTest->GetTrackAGET()[1]->GetHitXY()[0]->Draw();
    can3->cd(6);
    gMyBeamTest->GetTrackAGET()[1]->GetHitXY()[1]->Draw();
    gMyMainFrameGui->UpdateCanvas(cid++);

    //Track06
    gMyMainFrameGui->SwitchCanvas(cid);
    gMyMainFrameGui->ClearCanvas(cid);
    TCanvas *can4 = gMyMainFrameGui->GetCanvas(cid);
    can4->Divide(3, 2);
    can4->cd(1);
    gMyBeamTest->GetTrackAGET()[2]->GetAllHit()->Draw("colz");
    can4->cd(4);
    gMyBeamTest->GetTrackAGET()[2]->GetAllHit2()->Draw("colz");
    can4->cd(2);
    gMyBeamTest->GetTrackAGET()[2]->GetCharge()->Draw();
    can4->cd(3);
    gMyBeamTest->GetTrackAGET()[2]->GetHitXY()[0]->Draw();
    can4->cd(6);
    gMyBeamTest->GetTrackAGET()[2]->GetHitXY()[1]->Draw();
    gMyMainFrameGui->UpdateCanvas(cid++);

    //Track04
    gMyMainFrameGui->SwitchCanvas(cid);
    gMyMainFrameGui->ClearCanvas(cid);
    TCanvas *can5 = gMyMainFrameGui->GetCanvas(cid);
    can5->Divide(3, 2);
    can5->cd(1);
    gMyBeamTest->GetTrackVMM()[0]->GetAllHit()->Draw("colz");
    can5->cd(4);
    gMyBeamTest->GetTrackVMM()[0]->GetAllHit2()->Draw("colz");
    can5->cd(2);
    gMyBeamTest->GetTrackVMM()[0]->GetCharge()->Draw();
    can5->cd(3);
    gMyBeamTest->GetTrackVMM()[0]->GetHitXY()[0]->Draw();
    can5->cd(6);
    gMyBeamTest->GetTrackVMM()[0]->GetHitXY()[1]->Draw();
    gMyMainFrameGui->UpdateCanvas(cid++);

    //Track05
    gMyMainFrameGui->SwitchCanvas(cid);
    gMyMainFrameGui->ClearCanvas(cid);
    TCanvas *can6 = gMyMainFrameGui->GetCanvas(cid);
    can6->Divide(3, 2);
    can6->cd(1);
    gMyBeamTest->GetTrackVMM()[1]->GetAllHit()->Draw("colz");
    can6->cd(4);
    gMyBeamTest->GetTrackVMM()[1]->GetAllHit2()->Draw("colz");
    can6->cd(2);
    gMyBeamTest->GetTrackVMM()[1]->GetCharge()->Draw();
    can6->cd(3);
    gMyBeamTest->GetTrackVMM()[1]->GetHitXY()[0]->Draw();
    can6->cd(6);
    gMyBeamTest->GetTrackVMM()[1]->GetHitXY()[1]->Draw();
    gMyMainFrameGui->UpdateCanvas(cid++);

    /*
    TH2F *ftmp = gMyBeamTest->GetTrackAGET()[2]->GetAllHit();
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

    gMyMainFrameGui->SwitchCanvas(0);
}

void MyGuiActionClass::DoAnalysisCMDRoot()
{
    gMyBeamTest->AnalysisDSTRoot(gMyBeamTest->GetDSPath().Data());
}