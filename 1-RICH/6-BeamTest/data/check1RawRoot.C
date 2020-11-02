#include <iostream>
#include <vector>
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGSlider.h"
#include "TGButton.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGFileDialog.h"
#include "TApplication.h"
#include "TRootEmbeddedCanvas.h"
#include "TSystem.h"
#include "../inc/MyBeamTestDetector.h"
using namespace std;

const char *filetypes[] = {
    "ROOT files", "*-raw.root",
    "All files", "*",
    0, 0};

enum
{
    RICHType,
    AGETType,
    VMMType
};
const char *TYPE[3] = {"RICH", "Trk-AGET", "Trk_VMM"};

class MyReadRootGUI : public TGMainFrame
{
private:
    int ip = 0; //ped指针
    int ic = 0; //channel指针
    int ie = 0; //event指针
    int iboard = 0, ichip = 0, ichannel = 0;
    TH1F *fped, *fwav;

    TGHSlider *fHSlider1;
    TGTextButton *fPreButton1;
    TGTextButton *fNextButton1;
    TGHSlider *fHSlider2;
    TGTextButton *fPreButton2;
    TGTextButton *fNextButton2;
    TGTextButton *fPreButton3;
    TGTextButton *fNextButton3;
    TGTextButton *fPrintButton;
    TGTextButton *fCheckButton;
    TGTextButton *fFileButton;

    void DrawPed(int id);
    void DrawChannel(int id1, int id2);
    void DrawEvent(int id);
    void ReadFile(TString);

    TFile *fFile = 0;
    TTree *fTree = 0;
    TTree *fTree2 = 0;
    int type;

    TFile *fFile2 = 0;

    //AGET
    int Event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    UShort_t bad;
    long nentries;
    TBranch *b_Event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_data;
    TBranch *b_bad;
    vector<double> *data = 0;

    //VMM
    UShort_t PDO;
    UShort_t BCID;
    UShort_t TDO;
    TBranch *b_PDO;
    TBranch *b_BCID;
    TBranch *b_TDO;

    TGraph *gWave = 0;
    TCanvas *c1;

    vector<int> boardname;
    vector<int> chipname;

public:
    MyReadRootGUI();

    void CloseWindow();
    void OpenFile();
    void PrintData();
    void CheckData();
    void Created();

    //Draw pedstal
    void DrawSelected()
    {
        ip = fHSlider1->GetPosition();
        DrawPed(ip);
    }
    void DrawPre()
    {
        ip = (ip - 1 < 0) ? 0 : ip - 1;
        DrawPed(ip);
    }
    void DrawNext()
    {
        ip = (ip + 1 >= nentries) ? ip : ip + 1;
        DrawPed(ip);
    }

    //Draw channel
    void DrawSelectedChannel()
    {
        ic = fHSlider2->GetPosition();
        DrawChannel(ic, 0);
    }
    void DrawPreChannel()
    {
        ic = (ic - 1 < 0) ? 0 : ic - 1;
        DrawChannel(ic, 0);
    }
    void DrawNextChannel()
    {
        ic = (ic + 1 >= 1024) ? ic : ic + 1;
        DrawChannel(ic, 0);
    }

    //Draw Event
    void DrawPreEvent() { DrawEvent(0); }
    void DrawNextEvent() { DrawEvent(++ie); }

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
    int type = env->GetValue("checkRawType", -1);
    cout << "Read env from " << envPath << " and the check raw type: " << type << endl;
    if (type == -1)
        return;

    TString path = env->GetValue("checkRawRoot", "NOTSET");
    SetWindowName(path);
    ReadFile(path);

    /*
    if (type == RICH)
    {
        for (int i = 0; i < env->GetValue("NRICH", 1); i++)
        {
            for (int j = 0; j < env->GetValue(Form("RICH%d-Nboard", i), 4); j++)
                boardname.push_back(env->GetValue(Form("RICH%d-board%d", i, j), 1));

            for (int j = 0; j < env->GetValue(Form("RICH%d-NChip", i), 4); j++)
                chipname.push_back(env->GetValue(Form("RICH%d-chip%d", i, j), 1));
        }
    }

    if (type == TrackerAGET)
    {
        for (int i = 0; i < env->GetValue("NTrkAGT", 3); i++)
        {
            for (int j = 0; j < env->GetValue(Form("TrkAGET%d-Nboard", i), 4); j++)
                boardname.push_back(env->GetValue(Form("TrkAGET%d-board%d", i, j), 1));

            for (int j = 0; j < env->GetValue(Form("TrkAGET%d-NChip", i), 4); j++)
                chipname.push_back(env->GetValue(Form("TrkAGET%d-chip%d", i, j), 1));
        }
    }

    if (type == TrackerVMM)
    {
        for (int i = 0; i < env->GetValue("NTrkVMM", 3); i++)
        {
            for (int j = 0; j < env->GetValue(Form("TrkVMM%d-Nboard", i), 4); j++)
                boardname.push_back(env->GetValue(Form("TrkVMM%d-board%d", i, j), 1));

            for (int j = 0; j < env->GetValue(Form("TrkVMM%d-NChip", i), 4); j++)
                chipname.push_back(env->GetValue(Form("TrkVMM%d-chip%d", i, j), 1));
        }
    }
    */

    if (type != TrackerVMM)
    {
        path = env->GetValue("checkPedRoot", "NOTSET");
        fFile2 = new TFile(path);
    }
}

void MyReadRootGUI::OpenFile()
{
    static TString dir(".");
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
    fFileButton->SetEnabled(kFALSE);
    fHSlider1->SetEnabled(kTRUE);
    fPreButton1->SetEnabled(kTRUE);
    fNextButton1->SetEnabled(kTRUE);
    fHSlider2->SetEnabled(kTRUE);
    fPreButton2->SetEnabled(kTRUE);
    fNextButton2->SetEnabled(kTRUE);
    fPrintButton->SetEnabled(kTRUE);
    fCheckButton->SetEnabled(kTRUE);
    fPreButton3->SetEnabled(kFALSE);
    fNextButton3->SetEnabled(kFALSE);

    fFile = new TFile(fName);
    if (!fFile->IsOpen())
        return;

    boardname.clear();
    chipname.clear();
    vector<int> *bdlist = 0;
    vector<int> *chlist = 0;
    fTree = (TTree *)fFile->Get("tree");
    fTree2 = (TTree *)fFile->Get("tree2");
    fTree2->SetBranchAddress("type", &type);
    fTree2->SetBranchAddress("bdlist", &bdlist);
    fTree2->SetBranchAddress("chlist", &chlist);
    fTree2->GetEntry(0);

    for(int i=0; i<bdlist->size(); i++)
        boardname.push_back(bdlist->at(i));
    for(int i=0; i<chlist->size(); i++)
        chipname.push_back(chlist->at(i));

    cout << "This file belongs to " << TYPE[type] << "." << boardname.size() << " " << chipname.size() << endl;;
    cout << "The board: " << boardname[0] << " " << boardname[1] << "...";
    cout << "The chip: " << chipname[0] << " " << chipname[1] << "..." << endl;

    if (type == RICHType || type == AGETType)
    {
        fTree->SetBranchAddress("event", &Event, &b_Event);
        fTree->SetBranchAddress("board", &board, &b_board);
        fTree->SetBranchAddress("chip", &chip, &b_chip);
        fTree->SetBranchAddress("channel", &channel, &b_channel);
        fTree->SetBranchAddress("wave", &data, &b_data);
    }

    if (type == VMMType)
    {
        fTree->SetBranchAddress("event", &Event, &b_Event);
        fTree->SetBranchAddress("board", &board, &b_board);
        fTree->SetBranchAddress("chip", &chip, &b_chip);
        fTree->SetBranchAddress("channel", &channel, &b_channel);
        fTree->SetBranchAddress("PDO", &PDO, &b_PDO);
        fTree->SetBranchAddress("BCID", &BCID, &b_BCID);
        fTree->SetBranchAddress("TDO", &TDO, &b_TDO);
    }

    nentries = fTree->GetEntriesFast();
    fHSlider1->SetRange(0, nentries);
    //fTree->Print();

    if (type != VMMType)
    {
        //读ped-root
        fName.ReplaceAll("raw.root", "ped.root");
        fFile2 = new TFile(fName);
        if (!fFile2->IsOpen())
        {
            fPrintButton->SetEnabled(kFALSE);
            fCheckButton->SetEnabled(kFALSE);
            fHSlider2->SetEnabled(kFALSE);
            fPreButton2->SetEnabled(kFALSE);
            fNextButton2->SetEnabled(kFALSE);
            fPreButton3->SetEnabled(kFALSE);
            fNextButton3->SetEnabled(kFALSE);
        }
        cout << "Pedestal file opened: " << fName << endl;
    }
}

void MyReadRootGUI::DrawPed(int ip)
{
    Long64_t ii = fTree->LoadTree(ip);
    if (ii < 0)
        return;
    fTree->GetEntry(ii);

    if (type == VMMType)
    {
        PrintData();
        return;
    }

    if (gWave != 0)
        delete gWave;
    gWave = new TGraph(data->size());
    for (int i = 0; i < data->size(); i++)
        gWave->SetPoint(i + 1, i + 1, data->at(i));
    gWave->SetMarkerStyle(8);
    gWave->SetMarkerSize(1);
    gWave->SetMarkerColor(kRed);
    gWave->SetTitle(Form("Event:%d, Board:%d, Chip:%d, Channel:%d", (int)Event, board, chip, channel));
    c1->Clear();
    gWave->Draw("apl");
    c1->Modified();
    c1->Update();
}

//id1是通道号的编号，id2是event的编号
void MyReadRootGUI::DrawChannel(int id1, int id2)
{
    if (fFile2 == 0)
        return;

    iboard = id1 / (chipname.size() * 64);
    ichip = (id1 - iboard * chipname.size() * 64) / 64;
    ichannel = id1 - iboard * chipname.size() * 64 - ichip * 64;
    fped = (TH1F *)fFile2->Get(Form("Ped_%d_%d_%d", boardname[iboard], chipname[ichip], ichannel));
    fwav = (TH1F *)fFile2->Get(Form("Wave_%d_%d_%d", boardname[iboard], chipname[ichip], ichannel));
    //fwav->Scale(1. / (fwav->GetEntries() / 512.));
    cout << "--> Drawing " << id1 << ", " << ichip << ", " << ichannel << ": " << boardname[iboard] << " " << chipname[ichip] << " " << ichannel << endl;

    c1->Clear();
    c1->Divide(3, 2);
    if (fped != NULL || fwav != NULL)
    {
        c1->cd(1);
        if (fped != NULL)
            fped->Draw();
        c1->cd(4);
        if (fwav != NULL)
            fwav->Draw();
        c1->Modified();
        c1->Update();
    }

    if (type == VMMType)
        return;

    fPreButton3->SetEnabled(kTRUE);
    fNextButton3->SetEnabled(kTRUE);
    DrawEvent(id2);
}

void MyReadRootGUI::DrawEvent(int id)
{
    ie = id;

    int nw = 0;
    for (int i = 0; i < fTree->GetEntries(); i++)
    {
        fTree->GetEntry(i);

        if (Event < ie)
            continue;

        if (board == boardname[iboard] && chip == chipname[ichip] && ichannel == channel)
        {

            ie = Event;
            TGraph *gg = new TGraph(data->size());
            for (int ii = 0; ii < data->size(); ii++)
                gg->SetPoint(ii + 1, ii + 1, data->at(ii));
            gg->SetMarkerStyle(8);
            gg->SetMarkerSize(1);
            gg->SetMarkerColor(kRed);
            gg->SetTitle(Form("Event:%d, Board:%d, Chip:%d, Channel:%d", (int)Event, board, chip, channel));
            nw++;
            if (nw == 1)
                c1->cd(2);
            if (nw == 2)
                c1->cd(3);
            if (nw == 3)
                c1->cd(5);
            if (nw == 4)
                c1->cd(6);
            gg->Draw("apl");
            cout << "--> find " << nw << " event = " << Event << endl;
            if (nw == 4)
                break;
        }
    }

    c1->Modified();
    c1->Update();
    cout << "-->done.\n\n";
}

void MyReadRootGUI::PrintData()
{
    if (type == AGETType)
    {
        cout << "Event: " << (int)Event << " Board:" << board << " Chip:" << chip << " Channel:" << channel << endl
             << "--\t";
        /*for(int i=0; i<data->size(); i++)
    {
        cout<<data->at(i)<<" ";
        if((i+1)%8==0) cout<<" -- ";
        if((i+1)%32==0) cout<<endl<<"--\t";    
    }
    */
        cout << endl;
    }
    if (type == VMMType)
    {
        cout << "Event: " << (int)Event << " Board:" << board << " Chip:" << chip << " Channel:" << channel;
        cout << " --  Q: " << PDO << " T_corse:" << BCID << " T_fine:" << TDO << endl;
        ;
    }
}

void MyReadRootGUI::CheckData()
{
    //...
    cout << "checking data..." << endl;
    int tmp = -1;
    fTree->GetEntry(0);
    int trgstart = Event;
    fTree->GetEntry(nentries - 1);
    int trgstop = Event;

    int total = 0;
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = fTree->LoadTree(i);
        if (ii < 0)
            return;
        fTree->GetEntry(ii);
        //cout << Event << endl;
        tmp = (tmp == -1) ? Event : tmp;
        if (fabs(tmp - Event) > 1)
            cout << "Warning: Event id jumps from " << tmp << " to " << Event << endl;
        if (fabs(tmp - Event) != 0)
            total++;
        tmp = (tmp == Event) ? tmp : Event;
    }
    cout << "--> check is done" << endl;
    cout << "--> Trig start from " << trgstart << " to " << trgstop << ", total recorded trig number: " << total << endl;
}

MyReadRootGUI::MyReadRootGUI() : TGMainFrame(gClient->GetRoot(), 10, 10, kHorizontalFrame)
{
    TGLayoutHints *LayoutC = new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2);
    TGLayoutHints *LayoutX = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2);
    TGLayoutHints *LayoutY = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY, 2, 2, 2, 2);
    TGLayoutHints *LayoutXY = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);

    Connect("CloseWindow()", "MyReadRootGUI", this, "CloseWindow()");

    TGVerticalFrame *fVFrame1 = new TGVerticalFrame(this, 220, 500, kVerticalFrame);
    AddFrame(fVFrame1, LayoutY);
    {
        TGGroupFrame *fGroupFrame1 = new TGGroupFrame(fVFrame1, "Display single event:");
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
                fPreButton1 = new TGTextButton(fSliderFrame, "<", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fPreButton1->Connect("Clicked()", "MyReadRootGUI", this, "DrawPre()");
                fSliderFrame->AddFrame(fPreButton1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fPreButton1->SetEnabled(kFALSE);

                fNextButton1 = new TGTextButton(fSliderFrame, ">", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fNextButton1->Connect("Clicked()", "MyReadRootGUI", this, "DrawNext()");
                fSliderFrame->AddFrame(fNextButton1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fNextButton1->SetEnabled(kFALSE);
            }

            fPrintButton = new TGTextButton(fGroupFrame1, "Print Data", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
            fPrintButton->Connect("Clicked()", "MyReadRootGUI", this, "PrintData()");
            fGroupFrame1->AddFrame(fPrintButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
            fPrintButton->SetEnabled(kFALSE);

            fCheckButton = new TGTextButton(fGroupFrame1, "Check Data", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
            fCheckButton->Connect("Clicked()", "MyReadRootGUI", this, "CheckData()");
            fGroupFrame1->AddFrame(fCheckButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
            fCheckButton->SetEnabled(kFALSE);
        }

        TGGroupFrame *fGroupFrame2 = new TGGroupFrame(fVFrame1, "Display each channel:");
        fVFrame1->AddFrame(fGroupFrame2, LayoutC);
        fGroupFrame2->SetLayoutManager(new TGVerticalLayout(fGroupFrame2));
        {
            fHSlider2 = new TGHSlider(fGroupFrame2, 200, kSlider1 | kScaleBoth, -1, kHorizontalFrame);
            fHSlider2->Connect("Released()", "MyReadRootGUI", this, "DrawSelectedChannel()");
            fHSlider2->SetRange(0, 1024);
            fHSlider2->SetPosition(0);
            fGroupFrame2->AddFrame(fHSlider2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            fHSlider2->SetEnabled(kFALSE);

            TGHorizontalFrame *fSliderFrame2 = new TGHorizontalFrame(fGroupFrame2, 220, 100, kHorizontalFrame);
            fGroupFrame2->AddFrame(fSliderFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            {

                fPreButton2 = new TGTextButton(fSliderFrame2, "Pre Channel <", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fPreButton2->Connect("Clicked()", "MyReadRootGUI", this, "DrawPreChannel()");
                fSliderFrame2->AddFrame(fPreButton2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fPreButton2->SetEnabled(kFALSE);

                fNextButton2 = new TGTextButton(fSliderFrame2, "Next Channel >", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fNextButton2->Connect("Clicked()", "MyReadRootGUI", this, "DrawNextChannel()");
                fSliderFrame2->AddFrame(fNextButton2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fNextButton2->SetEnabled(kFALSE);
            }

            TGHorizontalFrame *fSliderFrame3 = new TGHorizontalFrame(fGroupFrame2, 220, 100, kHorizontalFrame);
            fGroupFrame2->AddFrame(fSliderFrame3, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            {

                fPreButton3 = new TGTextButton(fSliderFrame3, "Pre Event <", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fPreButton3->Connect("Clicked()", "MyReadRootGUI", this, "DrawPreEvent()");
                fSliderFrame3->AddFrame(fPreButton3, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fPreButton3->SetEnabled(kFALSE);

                fNextButton3 = new TGTextButton(fSliderFrame3, "Next Event >", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fNextButton3->Connect("Clicked()", "MyReadRootGUI", this, "DrawNextEvent()");
                fSliderFrame3->AddFrame(fNextButton3, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fNextButton3->SetEnabled(kFALSE);
            }
        }
    }

    TRootEmbeddedCanvas *fEmbeddedCanvas = new TRootEmbeddedCanvas(0, this, 1200, 580);
    c1 = fEmbeddedCanvas->GetCanvas();
    AddFrame(fEmbeddedCanvas, LayoutXY);
    SetWindowName("Check Raw-ROOT");

    MapSubwindows();
    Resize();
    MapWindow();

    Connect("Created()", "MyReadRootGUI", this, "Created()");
    Created();
}

//TBrowser b;

void check1RawRoot()
{
    new MyReadRootGUI();
}
