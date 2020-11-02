#include <iostream>
#include <vector>
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGSlider.h"
#include "TGButton.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGFileDialog.h"
#include "TApplication.h"
#include "TRootEmbeddedCanvas.h"
#include "TSystem.h"
#include "../inc/MyBeamTestDetector.h"

using namespace std;

const char *filetypes[] = {
    "ROOT files", "*-dst.root",
    "All files", "*",
    0, 0};

enum
{
    RICHType,
    T02Type,
    T03Type,
    T06Type,
    T04Type,
    T05Type
};

const char *TYPE[6] = {"RICH", "T2", "T3", "T6", "T4", "T5"};

int Color[7] = {kRed, kBlue, kGreen, kOrange, kMagenta, kYellow, kBlack};

class MyReadRootGUI : public TGMainFrame
{
private:
    int ip = 0;

    TGHSlider *fHSlider1;
    TGTextButton *fPreButton;
    TGTextButton *fNextButton;
    TGTextButton *fShowButton;
    TGTextButton *fCheckButton;
    TGTextButton *fFileButton;

    void DrawEvent(int id);
    void DrawRICHEvent(int id);
    void DrawAGETandVMMEvent(int id);

    TFile *fFile = 0;
    TTree *fTree = 0;
    TTree *fTree2 = 0;
    MyBeamTestData *fEvent = 0;
    long nentries;
    int type;

    TH2F *fHit = 0;
    vector<TGraph *> fCluster;
    TCanvas *c1;
    TH2F *hAllHit = 0;
    TH2F *hQvT = 0;
    TH1F *hAllHitX = 0;
    TH1F *hAllHitY = 0;
    TH1F *hCluster = 0; //cluster size
    TH1F *hBranch = 0;  //branch size = number of cluster
    TH1F *hCharge = 0;

    TString path;

public:
    MyReadRootGUI();

    void CloseWindow();
    void OpenFile();
    void ReadFile(TString fName);
    void PrintData();
    void CheckData();
    void Created();
    void ShowData();

    void DrawSelected()
    {
        ip = fHSlider1->GetPosition();
        DrawEvent(ip);
    }
    void DrawPre()
    {
        ip = (ip - 1 < 0) ? 0 : ip - 1;
        DrawEvent(ip);
    }
    void DrawNext()
    {
        ip = (ip + 1 >= nentries) ? ip : ip + 1;
        DrawEvent(ip);
    }

    ClassDef(MyReadRootGUI, 0)
};

void MyReadRootGUI::CloseWindow()
{
    gApplication->Terminate();
}

void MyReadRootGUI::Created()
{
    TString envPath = gSystem->WorkingDirectory() + TString("/.env");
    TEnv *env = new TEnv(envPath);
    path = env->GetValue("checkDSTRoot", "NOTSET");
    cout << "Read env from " << envPath << " and try to read the DST file from : " << path << endl;
    if (path == "NOTSET")
    {
        path = ".";
        return;
    }
    SetWindowName(path);
}

void MyReadRootGUI::OpenFile()
{
    static TString dir(path);
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    fi.fIniDir = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (fi.fFilename == NULL)
        return;
    ReadFile(fi.fFilename);
}

void MyReadRootGUI::ReadFile(TString fName)
{
    cout << "Opening: " << fName << endl;

    fFileButton->SetEnabled(kFALSE);
    fHSlider1->SetEnabled(kTRUE);
    fPreButton->SetEnabled(kTRUE);
    fNextButton->SetEnabled(kTRUE);
    fShowButton->SetEnabled(kTRUE);
    fCheckButton->SetEnabled(kTRUE);

    fFile = new TFile(fName);
    if (!fFile->IsOpen())
        return;

    fTree = (TTree *)fFile->Get("tree");
    fTree->SetBranchAddress("event", &fEvent);
    nentries = fTree->GetEntriesFast();
    fHSlider1->SetRange(-1, nentries);
    //fTree->Print();

    fTree2 = (TTree *)fFile->Get("tree2");
    fTree2->SetBranchAddress("type", &type);
    fTree2->GetEntry(0);
    cout << "This file belongs to " << TYPE[type] << endl;

    int NXbin = 0, Xbinlw = 0, Xbinup = 0;
    if (type == RICHType)
    {
        NXbin = 32;
        Xbinup = 32;
    }
    else
    {
        NXbin = 129;
        Xbinup = 129;
    }

    fHit = new TH2F("fHit", Form("Hit map for %s", TYPE[type]), NXbin, Xbinlw, Xbinup, NXbin, Xbinlw, Xbinup);
    fHit->SetXTitle("X");
    fHit->SetYTitle("Y");

    hAllHit = new TH2F("fAllHit", Form("Hit map for %s", TYPE[type]), NXbin, Xbinlw, Xbinup, NXbin, Xbinlw, Xbinup);
    hAllHit->SetXTitle("X");
    hAllHit->SetYTitle("Y");

    hAllHitX = new TH1F("fAllHitX", Form("Hit map for %s-X", TYPE[type]), NXbin, Xbinlw, Xbinup);
    hAllHitX->SetXTitle("X");
    hAllHitX->SetYTitle("Entries");

    hAllHitY = new TH1F("fAllHitY", Form("Hit map for %s-Y", TYPE[type]), NXbin, Xbinlw, Xbinup);
    hAllHitY->SetXTitle("Y");
    hAllHitY->SetYTitle("Entries");

    hQvT = new TH2F("hQvT", Form("Charge vs. time for %s", TYPE[type]), 128, 0, 512, 300, 100, 400);
    hQvT->SetXTitle("Charge");
    hQvT->SetYTitle("Time");

    hCharge = new TH1F("hCharge", Form("charge distribution for %s", TYPE[type]), 600, 0, 6000);
    hCharge->SetXTitle("Charge");
    hCharge->SetYTitle("Entries");

    hBranch = new TH1F("hBranch", Form("Number of clusters for %s", TYPE[type]), 20, 0, 20);
    hBranch->SetXTitle("number of Cluster");
    hBranch->SetYTitle("Time");

    hCluster = new TH1F("hCluster", Form("Cluster size for %s", TYPE[type]), 20, 0, 20);
    hCluster->SetXTitle("Cluster size");
    hCluster->SetYTitle("Entries");

    ShowData();
}

void MyReadRootGUI::DrawEvent(int ip)
{
    if (type == RICHType)
        DrawRICHEvent(ip);
    else if (type == T02Type || type == T03Type || type == T06Type)
        DrawAGETandVMMEvent(ip);
    else
        DrawAGETandVMMEvent(ip);
}

void MyReadRootGUI::ShowData()
{
    if (hAllHit->GetEntries() == 0)
    {
        if (type == RICHType)
        {
            cout << "Total entries: " << nentries << endl;

            for (int i = 0; i < nentries; i++)
            {
                Long64_t ii = fTree->LoadTree(i);
                if (ii < 0)
                    continue;
                fTree->GetEntry(ii);

                for (int j = 0; j < (int)fEvent->hit.size(); j++)
                {
                    //if (fEvent->hit[j].second < 25)
                    //    continue;
                    if (fEvent->time[j] < 230 || fEvent->time[j] > 245)
                        continue;
                    hAllHit->Fill(fEvent->hit[j].first, fEvent->hit[j].second); //, fEvent->charge[j]);
                    hAllHitX->Fill(fEvent->hit[j].first);
                    hAllHitY->Fill(fEvent->hit[j].second);
                    hQvT->Fill(fEvent->charge[j], fEvent->time[j]);
                    //hCharge->Fill(fEvent->charge[j]);
                }

                hBranch->Fill((int)fEvent->branch.size());
                for (int j = 0; j < (int)fEvent->branch.size(); j++)
                    hCluster->Fill(fEvent->branch[j].size());

                if (fEvent->branch.size() > 10)
                    continue;
                for (int j = 0; j < (int)fEvent->branch.size(); j++)
                {
                    int sum = 0;
                    int save = 1;
                    //if (fEvent->branch[j].size() != 2) // && fEvent->charge[j]<100)
                    //    continue;
                    if (fEvent->branch[j].size() >4 || fEvent->branch[j].size()<2)
                        continue;
                    for (int k = 0; k < (int)fEvent->branch[j].size(); k++)
                    {
                        if (fEvent->branch[j][k].hit.second < 26 || fEvent->branch[j][k].hit.second > 30)
                            save = 0;
                        if (fEvent->branch[j][k].t < 230 || fEvent->branch[j][k].t > 245)
                            save = 0;
                        sum += fEvent->branch[j][k].q;
                    }
                    if (save)
                        hCharge->Fill(sum);
                }
            }
        }
        else
        {
            for (int i = 0; i < nentries; i++)
            {
                Long64_t ii = fTree->LoadTree(i);
                if (ii < 0)
                    continue;
                fTree->GetEntry(ii);

                for (int j = 0; j < (int)fEvent->hit.size(); j++)
                {
                    double xhit = fEvent->hit[j].first;
                    double yhit = fEvent->hit[j].second;

                    if (xhit != -999)
                        hAllHitX->Fill(xhit);
                    if (yhit != -999)
                        hAllHitY->Fill(yhit);

                    for (int k = 0; k < hAllHit->GetXaxis()->GetNbins(); k++)
                    {
                        if (xhit == -999)
                        {
                            //yhit = (yhit % 2 == 0) ? yhit - 1 : yhit + 1;
                            hAllHit->Fill(hAllHit->GetXaxis()->GetBinCenter(k), yhit, fEvent->charge[j]);
                        }
                        if (yhit == -999)
                        {
                            //xhit = (xhit % 2 == 0) ? xhit - 1 : xhit + 1;
                            hAllHit->Fill(xhit, hAllHit->GetYaxis()->GetBinCenter(k), fEvent->charge[j]);
                        }
                    }
                    hQvT->Fill(fEvent->charge[j], fEvent->time[j]);
                    hCharge->Fill(fEvent->charge[j]);

                    hBranch->Fill(fEvent->XYbranch[0].size());
                    hBranch->Fill(fEvent->XYbranch[1].size());
                }
            }
        }
    }
    c1->Clear();
    c1->Divide(3, 2);
    c1->cd(1);
    hAllHit->Draw("colz");
    c1->cd(2);
    hCharge->Draw();
    c1->cd(3);
    hBranch->Draw();
    c1->cd(4);
    hAllHitX->Draw("colz");
    c1->cd(5);
    hAllHitY->Draw("colz");
    c1->cd(6);
    hCluster->Draw();
    c1->Modified();
    c1->Update();
}

void MyReadRootGUI::DrawRICHEvent(int ip)
{
    cout << "--> Draw Event: " << ip << endl;
    Long64_t ii = fTree->LoadTree(ip);
    if (ii < 0)
        return;
    fTree->GetEntry(ii);

    for (int i = 0; i < (int)fCluster.size(); i++)
        delete fCluster[i];
    fCluster.clear();

    TGraph *ff = new TGraph();
    ff->SetTitle(Form("Cluster Map for ID: %d", fEvent->event));
    ff->GetXaxis()->SetTitle("X");
    ff->GetYaxis()->SetTitle("Y");
    ff->SetPoint(0, fHit->GetXaxis()->GetXmin(), fHit->GetYaxis()->GetXmin());
    ff->SetPoint(1, fHit->GetXaxis()->GetXmin(), fHit->GetYaxis()->GetXmax());
    ff->SetPoint(2, fHit->GetXaxis()->GetXmax(), fHit->GetYaxis()->GetXmin());
    ff->SetPoint(3, fHit->GetXaxis()->GetXmax(), fHit->GetYaxis()->GetXmax());
    fCluster.push_back(ff);

    fHit->Reset();
    for (int j = 0; j < (int)fEvent->hit.size(); j++)
        fHit->Fill(fEvent->hit[j].first, fEvent->hit[j].second, fEvent->time[j]); // fEvent->charge[j]);

    for (int j = 0; j < (int)fEvent->branch.size(); j++)
    {
        TGraph *ftmp = new TGraph();
        fCluster.push_back(ftmp);
        for (int k = 0; k < (int)fEvent->branch[j].size(); k++)
            ftmp->SetPoint(k, fEvent->branch[j][k].hit.first, fEvent->branch[j][k].hit.second);
        ftmp->SetMarkerColor(Color[j % 7]);
        ftmp->SetMarkerSize(1.4);
        ftmp->SetMarkerStyle(21);
    }

    c1->Clear();
    c1->Divide(2, 1);
    c1->cd(1);
    fHit->Draw("lego2");
    c1->cd(2);
    fCluster[0]->Draw("ap");
    for (int i = 1; i < (int)fCluster.size(); i++)
        fCluster[i]->Draw("psame");
    c1->Modified();
    c1->Update();

    PrintData();
}

void MyReadRootGUI::DrawAGETandVMMEvent(int ip)
{
    cout << "--> Draw Event: " << ip << endl;
    Long64_t ii = fTree->LoadTree(ip);
    if (ii < 0)
        return;
    fTree->GetEntry(ii);

    for (int i = 0; i < (int)fCluster.size(); i++)
        delete fCluster[i];
    fCluster.clear();

    TGraph *ff = new TGraph();
    ff->SetTitle(Form("Cluster Map for ID: %d", fEvent->event));
    ff->GetXaxis()->SetTitle("X");
    ff->GetYaxis()->SetTitle("Y");
    ff->SetPoint(0, fHit->GetXaxis()->GetXmin(), fHit->GetYaxis()->GetXmin());
    ff->SetPoint(1, fHit->GetXaxis()->GetXmin(), fHit->GetYaxis()->GetXmax());
    ff->SetPoint(2, fHit->GetXaxis()->GetXmax(), fHit->GetYaxis()->GetXmin());
    ff->SetPoint(3, fHit->GetXaxis()->GetXmax(), fHit->GetYaxis()->GetXmax());
    fCluster.push_back(ff);

    fHit->Reset();
    for (int j = 0; j < (int)fEvent->hit.size(); j++)
    {
        double xhit = fEvent->hit[j].first;
        double yhit = fEvent->hit[j].second;

        for (int k = 0; k < fHit->GetXaxis()->GetNbins(); k++)
        {
            if (xhit == -999)
                fHit->Fill(fHit->GetXaxis()->GetBinCenter(k), yhit, fEvent->charge[j]);
            if (yhit == -999)
                fHit->Fill(xhit, fHit->GetYaxis()->GetBinCenter(k), fEvent->charge[j]);
        }
    }

    //XYbranch[0]
    for (int j = 0; j < (int)fEvent->XYbranch[0].size(); j++)
    {
        TGraph *ftmp = new TGraph();
        fCluster.push_back(ftmp);

        int ip = 0;
        for (int k = 0; k < (int)fEvent->XYbranch[0][j].size(); k++)
        {
            double xhit = fEvent->XYbranch[0][j][k].hit.first;
            double yhit = fEvent->XYbranch[0][j][k].hit.second;

            for (int m = 0; m < fHit->GetXaxis()->GetNbins(); m++)
            {
                if (xhit == -999)
                    ftmp->SetPoint(ip++, fHit->GetXaxis()->GetBinCenter(m), yhit);
                if (yhit == -999)
                    ftmp->SetPoint(ip++, xhit, fHit->GetYaxis()->GetBinCenter(m));
            }
        }
        ftmp->SetMarkerColor(Color[j % 7]);
        ftmp->SetMarkerSize(0.3);
        ftmp->SetMarkerStyle(21);
    }

    //ybranch
    for (int j = 0; j < (int)fEvent->XYbranch[1].size(); j++)
    {
        TGraph *ftmp = new TGraph();
        fCluster.push_back(ftmp);

        int ip = 0;
        for (int k = 0; k < (int)fEvent->XYbranch[1][j].size(); k++)
        {
            double xhit = fEvent->XYbranch[1][j][k].hit.first;
            double yhit = fEvent->XYbranch[1][j][k].hit.second;

            for (int m = 0; m < fHit->GetXaxis()->GetNbins(); m++)
            {
                if (xhit == -999)
                    ftmp->SetPoint(ip++, fHit->GetXaxis()->GetBinCenter(m), yhit);
                if (yhit == -999)
                    ftmp->SetPoint(ip++, xhit, fHit->GetYaxis()->GetBinCenter(m));
            }
        }
        ftmp->SetMarkerColor(Color[j % 7]);
        ftmp->SetMarkerSize(0.3);
        ftmp->SetMarkerStyle(21);
    }

    c1->Clear();
    c1->Divide(2, 1);
    c1->cd(1);
    fHit->Draw("colz");
    c1->cd(2);
    fCluster[0]->Draw("ap");
    for (int i = 1; i < (int)fCluster.size(); i++)
        fCluster[i]->Draw("psame");
    c1->Modified();
    c1->Update();

    PrintData();
}

void MyReadRootGUI::PrintData()
{
    cout << "Event: " << (int)fEvent->event << " has entries = " << fEvent->board.size() << " " << fEvent->chip.size() << " " << fEvent->channel.size() << " " << endl;
    for (int i = 0; i < fEvent->board.size(); i++)
    {
        cout << "->Entry " << i << ", board: " << fEvent->board[i] << " chip: " << fEvent->chip[i] << " channel: " << fEvent->channel[i] << " charge: " << fEvent->charge[i] << " @ (" << fEvent->hit[i].first << " " << fEvent->hit[i].second << ")" << endl;
    }
    cout << endl;
}

void MyReadRootGUI::CheckData()
{
    //...
    cout << "checking data..." << endl;
    int tmp = -1;
    fTree->GetEntry(0);
    int trgstart = fEvent->event;
    fTree->GetEntry(nentries - 1);
    int trgstop = fEvent->event;

    int total = 0;
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = fTree->LoadTree(i);
        if (ii < 0)
            return;
        fTree->GetEntry(ii);
        //cout << Event << endl;
        tmp = (tmp == -1) ? fEvent->event : tmp;
        if (fabs(tmp - fEvent->event) > 1)
            cout << "Warning: Event id jumps from " << tmp << " to " << fEvent->event << endl;
        if (fabs(tmp - fEvent->event) != 0)
            total++;
        tmp = (tmp == fEvent->event) ? tmp : fEvent->event;
    }
    cout << "--> check is done" << endl;
    cout << "--> Trig start from " << trgstart << " to " << trgstop << ", total recorded trig number: " << total << endl;
}

MyReadRootGUI::MyReadRootGUI() : TGMainFrame(gClient->GetRoot(), 10, 10, kHorizontalFrame)
{
    //c1 = new TCanvas();
    TGLayoutHints *LayoutC = new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2);
    TGLayoutHints *LayoutX = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2);
    TGLayoutHints *LayoutY = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY, 2, 2, 2, 2);
    TGLayoutHints *LayoutXY = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);

    Connect("CloseWindow()", "MyReadRootGUI", this, "CloseWindow()");

    TGVerticalFrame *fVFrame1 = new TGVerticalFrame(this, 220, 500, kVerticalFrame);
    AddFrame(fVFrame1, LayoutY);
    {
        TGGroupFrame *fGroupFrame1 = new TGGroupFrame(this, "Display single event:");
        fVFrame1->AddFrame(fGroupFrame1, LayoutC);
        fGroupFrame1->SetLayoutManager(new TGVerticalLayout(fGroupFrame1));
        {
            fFileButton = new TGTextButton(fGroupFrame1, "Open file", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
            fFileButton->Connect("Clicked()", "MyReadRootGUI", this, "OpenFile()");
            fGroupFrame1->AddFrame(fFileButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));

            fHSlider1 = new TGHSlider(fGroupFrame1, 200, kSlider1 | kScaleBoth, -1, kHorizontalFrame);
            fHSlider1->Connect("Released()", "MyReadRootGUI", this, "DrawSelected()");
            fHSlider1->SetRange(0, 40);
            fHSlider1->SetPosition(20);
            fGroupFrame1->AddFrame(fHSlider1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            fHSlider1->SetEnabled(kFALSE);

            TGHorizontalFrame *fSliderFrame = new TGHorizontalFrame(fGroupFrame1, 220, 100, kHorizontalFrame);
            fGroupFrame1->AddFrame(fSliderFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            {
                fPreButton = new TGTextButton(fSliderFrame, "<", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fPreButton->Connect("Clicked()", "MyReadRootGUI", this, "DrawPre()");
                fSliderFrame->AddFrame(fPreButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fPreButton->SetEnabled(kFALSE);

                fNextButton = new TGTextButton(fSliderFrame, ">", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fNextButton->Connect("Clicked()", "MyReadRootGUI", this, "DrawNext()");
                fSliderFrame->AddFrame(fNextButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fNextButton->SetEnabled(kFALSE);
            }

            fShowButton = new TGTextButton(fGroupFrame1, "Show Data", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
            fShowButton->Connect("Clicked()", "MyReadRootGUI", this, "ShowData()");
            fGroupFrame1->AddFrame(fShowButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
            fShowButton->SetEnabled(kFALSE);

            fCheckButton = new TGTextButton(fGroupFrame1, "Check Data", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
            fCheckButton->Connect("Clicked()", "MyReadRootGUI", this, "CheckData()");
            fGroupFrame1->AddFrame(fCheckButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
            fCheckButton->SetEnabled(kFALSE);
        }
    }

    TRootEmbeddedCanvas *fEmbeddedCanvas = new TRootEmbeddedCanvas(0, this, 1200, 580);
    c1 = fEmbeddedCanvas->GetCanvas();
    AddFrame(fEmbeddedCanvas, LayoutXY);
    SetWindowName("Check DST-ROOT");

    MapSubwindows();
    Resize();
    MapWindow();

    Connect("Created()", "MyReadRootGUI", this, "Created()");
    Created();
}

void check2DSTRoot()
{
    new MyReadRootGUI();
}
