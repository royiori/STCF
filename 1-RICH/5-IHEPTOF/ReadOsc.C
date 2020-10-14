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

#include "OscIOFuncs.h"
#include "OscData.h"
using namespace std;

const char *filetypes[] = {
    "bin files", "*.bin",
    "root files", "*.root",
    "All files", "*",
    0, 0};

enum
{
    READFILE
};

class MyReadRootGUI : public TGMainFrame
{
private:
    int ie = 0;       //event指针
    TGraph *gWave[4]; //最多4个channel的wave
    TCanvas *c1;

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

    void DrawEvent(int inc);

    TFile *fRtFile = 0;
    TTree *fRtTree = 0;
    long nEntry = 0; //nEntries;

    double t0[4], ts[4];
    TBranch *b_data1;
    TBranch *b_data2;
    TBranch *b_data3;
    TBranch *b_data4;
    vector<float> *data1 = 0;
    vector<float> *data2 = 0;
    vector<float> *data3 = 0;
    vector<float> *data4 = 0;

public:
    MyReadRootGUI();

    void CloseWindow();
    void OpenFile();

    //Draw Event
    void DrawPreEvent() { DrawEvent(-1); }
    void DrawNextEvent() { DrawEvent(1); }

private:
    void GUIStatus(int flag);
    void GetFileList(TString filePath, TString filePattern, vector<TString> &fList);

    void ReadBinFiles(TString fDir);
    void ReadBinToRoot(vector<TString> binLists, TString fName);
    void ReadRootFile(TString fName);

    void FillGraph(TGraph *g1, double x0, double xs, vector<float> *data)
    {
        if (data->size() == 0)
            return;

        for (int i = 0; i < data->size(); i++)
            g1->SetPoint(i, x0 + i * xs, data->at(i));
    }

    ClassDef(MyReadRootGUI, 0)
};

void MyReadRootGUI::CloseWindow()
{
    gApplication->Terminate();
}

//-----------------------------------------------------
// 改变控件状态
void MyReadRootGUI::GUIStatus(int flag)
{
    if (flag == READFILE)
    {
        fFileButton->SetEnabled(kFALSE);
        fHSlider1->SetEnabled(kFALSE);
        fPreButton1->SetEnabled(kFALSE);
        fNextButton1->SetEnabled(kFALSE);
        fHSlider2->SetEnabled(kFALSE);
        fPreButton2->SetEnabled(kFALSE);
        fNextButton2->SetEnabled(kFALSE);
        fPrintButton->SetEnabled(kFALSE);
        fCheckButton->SetEnabled(kFALSE);
        fPreButton3->SetEnabled(kTRUE);
        fNextButton3->SetEnabled(kTRUE);
    }
}

//-----------------------------------------------------
// 读取文件列表
void MyReadRootGUI::GetFileList(TString filePath, TString filePattern, vector<TString> &fList)
{
    char line[1000];
    fList.clear();

    FILE *fp = gSystem->OpenPipe("ls -l " + filePath + "/*" + filePattern, "r");
    if (!fp)
    {
        cout << "----> NO data files exists in " << filePath << "!" << endl;
        return;
    }

    while (fgets(line, sizeof(line), fp))
    {
        TString s(line);
        if (s.Index(filePattern) == -1)
            continue;

        int ind = -1;
        while (1)
        {
            if (s.Index(" ", ind) == -1)
                break;
            else
                ind = s.Index(" ", ind) + 1;
        }
        if (ind < 0)
            continue;

        s.Remove(0, ind);
        fList.push_back(s.ReplaceAll("\n", ""));
    }
}

//-----------------------------------------------------
// 选择文件所在目录
void MyReadRootGUI::OpenFile()
{
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    fi.fIniDir = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (fi.fFilename == NULL)
        return;

    if (TString(fi.fFilename).EndsWith("root"))
        ReadRootFile(fi.fFilename);
    else
        ReadBinFiles(fi.fIniDir);
}

//-----------------------------------------------------
// 打开文件夹
void MyReadRootGUI::ReadBinFiles(TString fDir)
{
    vector<TString> binLists;
    GetFileList(fDir, ".bin", binLists);

    if (binLists.size() > 0)
    {
        ReadBinToRoot(binLists, fDir + "/combine.root");
    }
    else
    {
        cout << "#### Warning: The folder \"" << fDir << "\" doesn't contain any bin files." << endl;
        return;
    }
}

//-----------------------------------------------------
// 读取一系列bin文件并保存到root
void MyReadRootGUI::ReadBinToRoot(vector<TString> binLists, TString fName)
{
    if (binLists.size() == 0)
        return;

    //打开root文件
    double tOrg[4]; //t0
    double tInc[4]; //tstep
    vector<float> wave1;
    vector<float> wave2;
    vector<float> wave3;
    vector<float> wave4;

    TFile *fFile = new TFile(fName, "RECREATE");
    if (!fFile->IsOpen())
        return;

    TTree *fTree = new TTree("waveform", "WaveForms");
    {
        fTree->Branch("wave1", &wave1);
        fTree->Branch("wave2", &wave2);
        fTree->Branch("wave3", &wave3);
        fTree->Branch("wave4", &wave4);

        fTree->Branch("tOrg1", &tOrg[0], "tOrg1/D");
        fTree->Branch("tOrg2", &tOrg[1], "tOrg2/D");
        fTree->Branch("tOrg3", &tOrg[2], "tOrg3/D");
        fTree->Branch("tOrg4", &tOrg[3], "tOrg4/D");

        fTree->Branch("tInc1", &tInc[0], "tInc1/D");
        fTree->Branch("tInc2", &tInc[1], "tInc2/D");
        fTree->Branch("tInc3", &tInc[2], "tInc3/D");
        fTree->Branch("tInc4", &tInc[3], "tInc4/D");
    }

    // 读取bin文件
    for (int i = 0; i < (int)binLists.size(); i++)
    {
        InputFile = fopen(binLists[i], "rb");
        cout << "--> Reading a new bin file-" << i << " : " << binLists[i] << endl;

        for (int j = 0; j < 4; j++)
        {
            tOrg[j] = 0;
            tInc[j] = 0;
        }

        size_t BytesRead = 0L;
        BytesRead = fread(&fileHeader, 1, 12, InputFile); //sizeof(FileHeader), InputFile);
        if ((fileHeader.Cookie[0] != COOKIE[0]) || (fileHeader.Cookie[1] != COOKIE[1]) || BytesRead != 12)
        {
            cout << "#### Error: The cookie does not fit DSO9254A waveform data, please check your input file! " << BytesRead << endl;
            return;
        }

        while (!feof(InputFile))
        {
            BytesRead = fread(&waveformHeader, 1, 140, InputFile); //sizeof(WaveformHeader), InputFile);
            if (BytesRead != 140)
                continue;

            if (waveformHeader.WaveformType != PB_NORMAL && waveformHeader.WaveformType != PB_AVERAGE)
            {
                cout << "#### Error: This is not the right waveform data type, please check your input file!" << BytesRead << endl;
                return;
            }

            double XIncrement = waveformHeader.XIncrement;
            double XOrigin = waveformHeader.XOrigin;
            TString fileLabel(waveformHeader.WaveformLabel);
            fileLabel.ReplaceAll("Channel ", "");
            int id = fileLabel.Atoi() - 1;
            if (id < 0 || id > 3)
            {
                cout << "#### Error: Channel id is wrong! Read in from " << waveformHeader.WaveformLabel << " to " << id << endl;
                return;
            }

            tOrg[id] = XOrigin;
            tInc[id] = XIncrement;

            cout << "---> * id: " << id << endl;
            cout << "     * T0: " << tOrg[id] << ", Ts: " << tInc[id] << endl;
            //cout << "     * NWaveformBuffers:" << waveformHeader.NWaveformBuffers << endl;
            //cout << "     * Points: " << waveformHeader.Points << " Average_Count: " << waveformHeader.Count << endl;
            //cout << "     * XUnits:" << waveformHeader.XUnits << " YUnits: " << waveformHeader.YUnits << endl;
            //cout << "     * Date: " << waveformHeader.Date << " " << waveformHeader.Time << ";  Frame: " << waveformHeader.Frame << ";  WaveformLabel: " << waveformHeader.WaveformLabel << endl;
            //cout << "     * Timetag:" << waveformHeader.TimeTag << "  SegmentIndex: " << waveformHeader.SegmentIndex << endl;

            BytesRead = fread(&waveformDataHeader, 1, 12, InputFile); //sizeof(WaveformDataHeader), InputFile);
            //cout << "---> * HeaderSize: " << waveformDataHeader.HeaderSize << endl;
            cout << "     * BufferType: " << waveformDataHeader.BufferType << " BytesPerPoint:" << waveformDataHeader.BytesPerPoint << " BufferSize: " << waveformDataHeader.BufferSize << endl;

            BytesRead = fread((char *)Volts, 1, waveformDataHeader.BufferSize, InputFile);
            //cout << "     * Byte to be read: " << BytesRead << endl;

            switch (id)
            {
            case 0:
                wave1.clear();
                break;
            case 1:
                wave2.clear();
                break;
            case 2:
                wave3.clear();
                break;
            default:
                wave4.clear();
                break;
            }

            for (int j = 0; j < (BytesRead) / waveformDataHeader.BytesPerPoint; j++)
            {
                switch (id)
                {
                case 0:
                    wave1.push_back(Volts[j]);
                    break;
                case 1:
                    wave2.push_back(Volts[j]);
                    break;
                case 2:
                    wave3.push_back(Volts[j]);
                    break;
                default:
                    wave4.push_back(Volts[j]);
                    break;
                }
            }
        }

        fTree->Fill();
    }

    cout << "--> Total waveforms = " << fTree->GetEntries() << endl;
    cout << "--> Binary files have been convert to root-file: " << fName << ".\n\n";

    fTree->Write();
    fFile->Flush();
    fFile->Close();

    ReadRootFile(fName);
}

//-----------------------------------------------------
// 读取一系列bin文件并保存到root
void MyReadRootGUI::ReadRootFile(TString fName)
{
    fRtFile = new TFile(fName);
    if (!fRtFile->IsOpen())
        return;

    GUIStatus(READFILE);

    fRtTree = (TTree *)fRtFile->Get("waveform");
    fRtTree->SetBranchAddress("wave1", &data1, &b_data1);
    fRtTree->SetBranchAddress("wave2", &data2, &b_data2);
    fRtTree->SetBranchAddress("wave3", &data3, &b_data3);
    fRtTree->SetBranchAddress("wave4", &data4, &b_data4);

    fRtTree->SetBranchAddress("tOrg1", &t0[0]);
    fRtTree->SetBranchAddress("tOrg2", &t0[1]);
    fRtTree->SetBranchAddress("tOrg3", &t0[2]);
    fRtTree->SetBranchAddress("tOrg4", &t0[3]);

    fRtTree->SetBranchAddress("tInc1", &ts[0]);
    fRtTree->SetBranchAddress("tInc2", &ts[1]);
    fRtTree->SetBranchAddress("tInc3", &ts[2]);
    fRtTree->SetBranchAddress("tInc4", &ts[3]);

    nEntry = fRtTree->GetEntries();

    ie = 0;
    DrawEvent(0);
}

//-----------------------------------------------------
void MyReadRootGUI::DrawEvent(int inc)
{
    if (ie + inc < 0)
        ie = 0;
    else if (ie + inc >= nEntry)
        ie = nEntry - 1;
    else
        ie = ie + inc;

    fRtTree->GetEntry(ie);
    cout << "--> Drawing events " << ie << endl;

    FillGraph(gWave[0], t0[0], ts[0], data1);
    FillGraph(gWave[1], t0[1], ts[1], data2);
    FillGraph(gWave[2], t0[2], ts[2], data3);
    FillGraph(gWave[3], t0[3], ts[3], data4);

    c1->Clear();

    int chstart = -1;
    for (int i = 0; i < 4; i++)
    {
        if (gWave[i]->GetN() == 0)
            continue;

        if (chstart == -1)
        {
            gWave[i]->Draw("apl");
            chstart = i;
        }
        else
        {
            gWave[i]->Draw("aplsame");
        }
    }

    c1->Modified();
    c1->Update();
}

MyReadRootGUI::MyReadRootGUI() : TGMainFrame(gClient->GetRoot(), 10, 10, kHorizontalFrame)
{
    int fColor[4] = {kRed, kOrange - 1, kBlue, kGreen + 2};

    for (int i = 0; i < 4; i++)
    {
        gWave[i] = new TGraph();
        gWave[i]->SetName(Form("gWave%d", i));
        gWave[i]->SetTitle(Form("Waveform for channel %d", i + 1));
        gWave[i]->GetXaxis()->SetTitle("T");
        gWave[i]->GetYaxis()->SetTitle("V");
        gWave[i]->SetLineColor(fColor[i]);
        gWave[i]->SetMarkerColor(fColor[i]);
    }

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
            //fHSlider1->Connect("Released()", "MyReadRootGUI", this, "DrawSelected()");
            fHSlider1->SetRange(0, 40);
            fHSlider1->SetPosition(20);
            fGroupFrame1->AddFrame(fHSlider1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            fHSlider1->SetEnabled(kFALSE);

            TGHorizontalFrame *fSliderFrame = new TGHorizontalFrame(fGroupFrame1, 220, 100, kHorizontalFrame);
            fGroupFrame1->AddFrame(fSliderFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            {
                fPreButton1 = new TGTextButton(fSliderFrame, "<", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                //fPreButton1->Connect("Clicked()", "MyReadRootGUI", this, "DrawPre()");
                fSliderFrame->AddFrame(fPreButton1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fPreButton1->SetEnabled(kFALSE);

                fNextButton1 = new TGTextButton(fSliderFrame, ">", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                //fNextButton1->Connect("Clicked()", "MyReadRootGUI", this, "DrawNext()");
                fSliderFrame->AddFrame(fNextButton1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fNextButton1->SetEnabled(kFALSE);
            }

            fPrintButton = new TGTextButton(fGroupFrame1, "Print Data", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
            //fPrintButton->Connect("Clicked()", "MyReadRootGUI", this, "PrintData()");
            fGroupFrame1->AddFrame(fPrintButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
            fPrintButton->SetEnabled(kFALSE);

            fCheckButton = new TGTextButton(fGroupFrame1, "Check Data", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
            //fCheckButton->Connect("Clicked()", "MyReadRootGUI", this, "CheckData()");
            fGroupFrame1->AddFrame(fCheckButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
            fCheckButton->SetEnabled(kFALSE);
        }

        TGGroupFrame *fGroupFrame2 = new TGGroupFrame(fVFrame1, "Display each channel:");
        fVFrame1->AddFrame(fGroupFrame2, LayoutC);
        fGroupFrame2->SetLayoutManager(new TGVerticalLayout(fGroupFrame2));
        {
            fHSlider2 = new TGHSlider(fGroupFrame2, 200, kSlider1 | kScaleBoth, -1, kHorizontalFrame);
            //fHSlider2->Connect("Released()", "MyReadRootGUI", this, "DrawSelectedChannel()");
            fHSlider2->SetRange(0, 1024);
            fHSlider2->SetPosition(0);
            fGroupFrame2->AddFrame(fHSlider2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            fHSlider2->SetEnabled(kFALSE);

            TGHorizontalFrame *fSliderFrame2 = new TGHorizontalFrame(fGroupFrame2, 220, 100, kHorizontalFrame);
            fGroupFrame2->AddFrame(fSliderFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            {

                fPreButton2 = new TGTextButton(fSliderFrame2, "Pre Channel <", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                //fPreButton2->Connect("Clicked()", "MyReadRootGUI", this, "DrawPreChannel()");
                fSliderFrame2->AddFrame(fPreButton2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                fPreButton2->SetEnabled(kFALSE);

                fNextButton2 = new TGTextButton(fSliderFrame2, "Next Channel >", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                //fNextButton2->Connect("Clicked()", "MyReadRootGUI", this, "DrawNextChannel()");
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
}

void ReadOsc()
{
    new MyReadRootGUI();
}
