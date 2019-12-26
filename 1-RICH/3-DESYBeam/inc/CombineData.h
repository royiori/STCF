#ifndef CombineData_h
#define CombineData_h
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGraph2D.h"
using namespace std;

enum
{
    RICHType,
    AGETType,
    VMMType
};

//---------------------
// support functions
unsigned short AGETExchange(unsigned short memblock)
{
    unsigned short aa = memblock;
    unsigned short low = memblock & 0xff;
    unsigned short high = (memblock >> 8) & 0xFF;
    aa = low;
    aa <<= 8;
    aa |= high;
    return aa;
}

unsigned short VMMExchange(unsigned short memblock)
{
    unsigned short aa = memblock;
    unsigned short t = 0;
    unsigned short rt = 0;
    int NO_OF_BITS = 12; // number of bits to be conversed

    for (int i = 0; i < NO_OF_BITS; i++)
    {
        t = aa & (1 << i);
        if (t)
            rt |= (1 << (NO_OF_BITS - 1 - i));
    }

    string gray = bitset<12>(rt).to_string(); // convert Gray code to string
    string binary;
    binary = gray[0];

    for (int i = 0; i < NO_OF_BITS - 1; i++) // convert string number to string of binary number
    {
        if (binary[i] == gray[i + 1])
            binary = binary + "0";
        else
            binary = binary + "1";
    }

    //cout << "gray = " << gray << endl;
    //cout << "binary = " << binary << endl;
    rt = (unsigned short)stoul(binary, nullptr, 2); // convert string of binary number to decimal number
    //cout << "rt = " << rt << endl;
    return rt;
}

unsigned short Reverse(unsigned short number, int kbits)
{
    unsigned short aa = number;
    unsigned short t = 0;
    unsigned short rt = 0;
    int NO_OF_BITS = kbits; // number of bits to be conversed

    for (int i = 0; i < NO_OF_BITS; i++)
    {
        t = aa & (1 << i);
        if (t)
            rt |= (1 << (NO_OF_BITS - 1 - i));
    }

    return rt;
}

unsigned short BitAttraction(unsigned short memblock, int kbits, int startbit)
{
    return (((1 << kbits) - 1) & (memblock >> (startbit - 1)));
}

unsigned short BitAttraction2(unsigned short memblock, int kbits, int startbit)
{
    unsigned short a = (1 << kbits) - 1;
    unsigned short b = memblock >> (startbit - 1);
    unsigned short c = a & b;
    unsigned short d = c << 2;
    unsigned short e = Reverse(d, 8);
    //cout<<memblock<<"+"<<kbits<<"+"<<startbit<<" -> "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<endl;
    return e;
}

//-------------------------
// 读取文件列表
void GetFileList(TString filePath, TString filePattern, vector<TString> &fList)
{
    char line[1000];
    fList.clear();

    FILE *fp = gSystem->OpenPipe("ls " + filePath + "*" + filePattern, "r");
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
        fList.push_back(s.ReplaceAll("\n", ""));
    }
}

//-------------------------------
// RICH Data to Raw-ROOT
void ReadRICHData2Root(vector<TString> datList, TString fRawName, int force = 1);
void ReadRICHData2Root(vector<TString> datList, TString fRawName, int force)
{
    if (datList.size() == 0)
        return;

    int length = sizeof(unsigned short);
    unsigned short memblock;
    int boardstart = 0;
    int chipstart = 0;
    int type;
    //bool headerFlag;
    //bool trailerFlag;
    double counting = 0;
    //unsigned short triggerType;
    int badNumber = 0;
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> wave;

    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fRawName, fStat);
        if (fStat.fSize != 0)
            return;
    }

    TFile *fFile = new TFile(fRawName, "RECREATE");
    if (!fFile->IsOpen())
        return;
    TTree *fTree = new TTree("tree", "tree");
    fTree->Branch("event", &event);
    fTree->Branch("board", &board);
    fTree->Branch("chip", &chip);
    fTree->Branch("channel", &channel);
    fTree->Branch("wave", &wave);

    TTree *fTree2 = new TTree("tree2", "tree2");
    fTree2->Branch("type", &type);
    type = RICHType;
    fTree2->Fill();
    fTree2->Write();

    vector<int> bdList;
    for (int i = 0; i < (int)datList.size(); i++)
    {
        fstream InDat;
        InDat.open(datList[i], ios::in | ios::binary);
        cout << "--> Reading a new dat file-" << i << " : " << datList[i] << endl;

        while (!InDat.eof())
        {
            InDat.read((char *)(&memblock), length);
            memblock = AGETExchange(memblock);
            if (memblock == 0xeeee)
            {
                //headerFlag = 1;
                counting = 0;
                event = 0; // reset
                           //cout << "header" << endl;
            }              // data packet header
            //else if (memblock&0xe == 0xe && headerFlag == 1 && counting == 0)	// board number, chip number, trigger type
            else if (counting == 1) // board number, chip number, trigger type
            {
                //board = (((memblock&0x0f00)>>8 - boardstart)); //cout << "board number = " << board << endl;
                board = (((memblock & 0x0f00) >> 8) - boardstart); //cout << "board number = " << board << endl;

                if (bdList.size() < 4)
                {
                    int flag = -1;
                    for (size_t ii = 0; ii < bdList.size(); ii++)
                        if (board == bdList[ii])
                            flag = 1;
                    if (flag == -1)
                        bdList.push_back(board);
                }
                chip = (((memblock & 0x00f0) >> 4) - chipstart); //cout << "chip number = " << chip << endl;
                //triggerType = (memblock & 0x000f);               //cout << "trigger type = " << triggerType << endl;
                //                                                 //cout << "fixed number" << endl;// counting = 1;
            }
            else if (counting > 1 && counting < 6)
            {
            }
            else if (counting == 6)
            {
                channel = memblock; /*cout << "channel number = " << channel << endl;*/
            }
            else if (counting > 6 && counting < 10)
            {
                if (counting == 7)
                    event += memblock * (2 ^ 32);
                else if (counting == 8)
                    event += memblock * (2 ^ 16);
                else if (counting == 9)
                {
                    event += memblock; /*cout << "trigger Number = " << event << endl;*/
                }
            }
            //else if (counting > 9 && counting < 521)
            else if (memblock == 0xffff)
            {
                //trailerFlag = 1;
                if (wave.size() == 512)
                {
                    fTree->Fill();
                }
                else if (wave.size() != 0)
                { /*cout << "wave.size() = " << wave.size() << endl;*/
                    badNumber++;
                }
                wave.clear();
                //cout << "clear ..." << endl;
            }
            else
            {
                double buffer = memblock & 0x0fff;
                wave.push_back(buffer);
                //cout << counting << "\t" << buffer << "\t" << wave.size() << endl;
            }

            counting++;
        }
    }
    cout << "--> Total entries = " << fTree->GetEntries() << endl;
    cout << "--> Bad Events = " << badNumber << endl;
    cout << "--> Board list: ";
    for (size_t i = 0; i < bdList.size(); i++)
        cout << bdList[i] << " ";
    cout << endl;
    fTree->Write();
    //fTree->Print();
    fFile->Flush();
    fFile->Close();

    cout << "--> RICH Raw data file has been convert to root-file: " << fRawName << ".\n"
         << endl;
}

void GenerateRICHPed(TString fRawName, TString fPedName, vector<int> boardName, vector<int> chipName, int force = 1)
{
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fPedName, fStat);
        if (fStat.fSize != 0)
            return;
    }

    //init
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> *wave = 0;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;

    int boardNo = 0;
    int chipNo = 0;
    boardNo = boardName.size();
    chipNo = chipName.size();
    int dim = boardNo * chipNo * 64;

    //定义分布图
    vector<TH1D *> hNoise;
    vector<TH1D *> hWave;
    vector<TGraph2D *> gWaveAll;
    hNoise.resize(dim);
    hWave.resize(dim);
    gWaveAll.resize(boardNo);
    for (int i = 0; i < boardNo; i++)
    {
        gWaveAll[i] = new TGraph2D();
        gWaveAll[i]->SetName(Form("gWaveAll%d", i));
        gWaveAll[i]->SetTitle(Form("Wave Sum for Board %d", boardName[i]));
        gWaveAll[i]->GetXaxis()->SetTitle("chip+channel");
        gWaveAll[i]->GetYaxis()->SetTitle("time");
        gWaveAll[i]->GetZaxis()->SetTitle("amplitude");

        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {
                int id = i * (64 * 4) + j * 64 + k;
                TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", boardName[i], chipName[j], k));
                if (tmp != NULL)
                    delete tmp;
                hNoise[id] = new TH1D(Form("Ped_%d_%d_%d", boardName[i], chipName[j], k), Form("Pedestal for Board=%d, Chip=%d, Channel=%d", boardName[i], chipName[j], k), 1000, 0, 1000);
                hNoise[id]->GetXaxis()->SetTitle("amplitude");
                hNoise[id]->GetYaxis()->SetTitle("Entries");

                tmp = (TH1D *)gDirectory->Get(Form("Wave_%d_%d_%d", boardName[i], chipName[j], k));
                if (tmp != NULL)
                    delete tmp;
                hWave[id] = new TH1D(Form("Wave_%d_%d_%d", boardName[i], chipName[j], k), Form("Waveform for Board=%d, Chip=%d, Channel=%d", boardName[i], chipName[j], k), 512, 0, 512);
                hWave[id]->GetXaxis()->SetTitle("time");
                hWave[id]->GetYaxis()->SetTitle("amplitude");
            }
    }

    TH1F *tmp = (TH1F *)gDirectory->Get("fNoiseMean");
    if (tmp != NULL)
        delete tmp;
    tmp = (TH1F *)gDirectory->Get("fNoiseRMS");
    if (tmp != NULL)
        delete tmp;
    TH1F *fNoiseMean = new TH1F("fNoiseMean", "noise mean distribution", dim + 1, 0, dim + 1);
    TH1F *fNoiseRMS = new TH1F("fNoiseRMS", "noise rms distribution", dim + 1, 0, dim + 1);
    fNoiseMean->GetXaxis()->SetTitle("board+chip+channel");
    fNoiseRMS->GetXaxis()->SetTitle("board+chip+channel");
    fNoiseMean->GetYaxis()->SetTitle("Pedestal mean");
    fNoiseRMS->GetYaxis()->SetTitle("Pedestal sigma");

    //读取raw-root文件
    TFile *rawFile = new TFile(fRawName);
    if (!rawFile->IsOpen())
    {
        cout << "#### Error: Can't open " << fRawName << ". Please check the filename." << endl;
        return;
    }
    TFile *pedFile = new TFile(fPedName, "recreate");
    cout << "--> Analyzing pedestal files from " << fRawName << endl;
    TTree *tree = (TTree *)rawFile->Get("tree");
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);

    //填图
    int nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;

        tree->GetEntry(ii);

        //找出board对应的编号
        int iBoardIndex = -1;
        for (int j = 0; j < (int)boardName.size(); j++)
            if (boardName[j] == board)
                iBoardIndex = j;

        if (iBoardIndex == -1)
        {
            cout << "The board " << board << " does not exist in the current RICH list, check your settings." << endl;
            break;
        }

        int index = iBoardIndex * (64 * 4) + (chip - 10) * 64 + channel;

        //填图
        for (int j = 0; j < (int)wave->size(); j++)
        {
            hWave[index]->SetBinContent(j, hWave[index]->GetBinContent(j) + wave->at(j));
            if (j < 150 || j > 320)
                hNoise[index]->Fill(wave->at(j));
        }
    }

    for (int i = 0; i < boardNo; i++)
    {
        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {
                double mean = hNoise[i * (64 * 4) + j * 64 + k]->GetMean();
                double rms = hNoise[i * (64 * 4) + j * 64 + k]->GetRMS();
                int id = i * (64 * 4) + j * 64 + k;
                if (hNoise[id]->GetEntries() > 0)
                    hNoise[id]->Fit("gaus", "RQ", "", mean - 3 * rms, mean + 3 * rms);
                hNoise[id]->Write();
                hWave[id]->Write();
                fNoiseMean->SetBinContent(id, mean);
                fNoiseMean->SetBinError(id, rms);
                fNoiseRMS->SetBinContent(id, rms);

                for (int ii = 0; ii < 512; ii++)
                    gWaveAll[i]->SetPoint((j * 64 + k) * 521 + ii, j * 64 + k, ii, hWave[id]->GetBinContent(ii));
            }
        gWaveAll[i]->Write();
    }
    fNoiseMean->Write();
    fNoiseRMS->Write();

    pedFile->Write();
    pedFile->Close();
    rawFile->Close();
    cout << "--> Pedestal has been stored to " << fPedName << "\n\n";
}

void GenerateRICHPed(TString fRawName, TString fPedName, vector<MyBeamTestRICH *> vRICH)
{
    if (vRICH.size() == 0)
        return;

    vector<int> boardName;
    vector<int> chipName;

    boardName = vRICH[0]->GetBoardList();
    chipName = vRICH[0]->GetChipList();

    GenerateRICHPed(fRawName, fPedName, boardName, chipName, 1);
}

//-------------------------------
// Track-AGET Data to Raw-ROOT
void ReadTrackAGTData2Root(vector<TString> datList, TString fRawName, int force = 1);
void ReadTrackAGTData2Root(vector<TString> datList, TString fRawName, int force)
{
    if (datList.size() == 0)
        return;

    int length = sizeof(unsigned short);
    unsigned short memblock;
    int boardstart = 0;
    int chipstart = 0;
    int type;
    //bool headerFlag;
    //bool trailerFlag;
    double counting = 0;
    //unsigned short triggerType;
    int badNumber = 0;
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> wave;

    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fRawName, fStat);
        if (fStat.fSize != 0)
            return;
    }

    TFile *fFile = new TFile(fRawName, "RECREATE");
    if (!fFile->IsOpen())
        return;
    TTree *fTree = new TTree("tree", "tree");
    fTree->Branch("event", &event);
    fTree->Branch("board", &board);
    fTree->Branch("chip", &chip);
    fTree->Branch("channel", &channel);
    fTree->Branch("wave", &wave);

    TTree *fTree2 = new TTree("tree2", "tree2");
    fTree2->Branch("type", &type);
    type = AGETType;
    fTree2->Fill();
    fTree2->Write();

    vector<int> bdList;
    for (int i = 0; i < (int)datList.size(); i++)
    {
        fstream InDat;
        InDat.open(datList[i], ios::in | ios::binary);
        cout << "--> Reading a new dat file-" << i << " : " << datList[i] << endl;

        while (!InDat.eof())
        {
            InDat.read((char *)(&memblock), length);
            memblock = AGETExchange(memblock);
            if (memblock == 0xeeee)
            {
                //headerFlag = 1;
                counting = 0;
                event = 0; // reset
                           //cout << "header" << endl;
            }              // data packet header
            //else if (memblock&0xe == 0xe && headerFlag == 1 && counting == 0)	// board number, chip number, trigger type
            else if (counting == 1) // board number, chip number, trigger type
            {
                //board = (((memblock&0x0f00)>>8 - boardstart)); //cout << "board number = " << board << endl;
                board = (((memblock & 0x0f00) >> 8) - boardstart); //cout << "board number = " << board << endl;
                //if (board == 7) board = 5;
                //if (board == 9) board = 1;
                if (bdList.size() < 4)
                {
                    int flag = -1;
                    for (size_t ii = 0; ii < bdList.size(); ii++)
                        if (board == bdList[ii])
                            flag = 1;
                    if (flag == -1)
                        bdList.push_back(board);
                }

                chip = (((memblock & 0x00f0) >> 4) - chipstart); //cout << "chip number = " << chip << endl;
                //triggerType = (memblock & 0x000f);               //cout << "trigger type = " << triggerType << endl;
                //                                                 //cout << "fixed number" << endl;// counting = 1;
            }
            else if (counting > 1 && counting < 6)
            {
            }
            else if (counting == 6)
            {
                channel = memblock; /*cout << "channel number = " << channel << endl;*/
            }
            else if (counting > 6 && counting < 10)
            {
                if (counting == 7)
                    event += memblock * (2 ^ 32);
                else if (counting == 8)
                    event += memblock * (2 ^ 16);
                else if (counting == 9)
                {
                    event += memblock; /*cout << "trigger Number = " << event << endl;*/
                }
            }
            //else if (counting > 9 && counting < 521)
            else if (memblock == 0xffff)
            {
                //trailerFlag = 1;
                if (wave.size() == 512)
                {
                    fTree->Fill();
                }
                else if (wave.size() != 0)
                { /*cout << "wave.size() = " << wave.size() << endl;*/
                    badNumber++;
                }
                wave.clear();
                //cout << "clear ..." << endl;
            }
            else
            {
                double buffer = memblock & 0x0fff;
                wave.push_back(buffer);
                //cout << counting << "\t" << buffer << "\t" << wave.size() << endl;
            }

            counting++;
        }
    }
    cout << "--> Total entries = " << fTree->GetEntries() << endl;
    cout << "--> Bad Events = " << badNumber << endl;
    cout << "--> Board list: ";
    for (size_t i = 0; i < bdList.size(); i++)
        cout << bdList[i] << " ";
    cout << endl;
    fTree->Write();
    //fTree->Print();
    fFile->Flush();
    fFile->Close();

    cout << "--> Track-AGET Raw data file has been convert to root-file: " << fRawName << ".\n"
         << endl;
}

void GenerateAGETPed(TString fRawName, TString fPedName, vector<int> boardName, vector<int> chipName, int force = 1)
{
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fPedName, fStat);
        if (fStat.fSize != 0)
            return;
    }

    //init
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> *wave = 0;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;

    int boardNo = 0;
    int chipNo = 0;

    boardNo = boardName.size();
    chipNo = chipName.size();
    int dim = boardNo * chipNo * 64;

    //定义分布图
    vector<TH1D *> hNoise;
    vector<TH1D *> hWave;
    vector<TGraph2D *> gWaveAll;
    hNoise.resize(dim);
    hWave.resize(dim);
    gWaveAll.resize(boardNo);
    for (int i = 0; i < boardNo; i++)
    {
        gWaveAll[i] = new TGraph2D();
        gWaveAll[i]->SetName(Form("gWaveAll%d", i));
        gWaveAll[i]->SetTitle(Form("Wave Sum for Board %d", boardName[i]));
        gWaveAll[i]->GetXaxis()->SetTitle("chip+channel");
        gWaveAll[i]->GetYaxis()->SetTitle("time");
        gWaveAll[i]->GetZaxis()->SetTitle("amplitude");

        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {
                int id = i * (64 * 4) + j * 64 + k;
                TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", boardName[i], chipName[j], k));
                if (tmp != NULL)
                    delete tmp;
                hNoise[id] = new TH1D(Form("Ped_%d_%d_%d", boardName[i], chipName[j], k), Form("Pedestal for Board=%d, Chip=%d, Channel=%d", boardName[i], chipName[j], k), 1000, 0, 1000);
                hNoise[id]->GetXaxis()->SetTitle("amplitude");
                hNoise[id]->GetYaxis()->SetTitle("Entries");

                tmp = (TH1D *)gDirectory->Get(Form("Wave_%d_%d_%d", boardName[i], chipName[j], k));
                if (tmp != NULL)
                    delete tmp;
                hWave[id] = new TH1D(Form("Wave_%d_%d_%d", boardName[i], chipName[j], k), Form("Waveform for Board=%d, Chip=%d, Channel=%d", boardName[i], chipName[j], k), 512, 0, 512);
                hWave[id]->GetXaxis()->SetTitle("time");
                hWave[id]->GetYaxis()->SetTitle("amplitude");
            }
    }

    TH1F *tmp = (TH1F *)gDirectory->Get("fNoiseMean");
    if (tmp != NULL)
        delete tmp;
    tmp = (TH1F *)gDirectory->Get("fNoiseRMS");
    if (tmp != NULL)
        delete tmp;
    TH1F *fNoiseMean = new TH1F("fNoiseMean", "noise mean distribution", dim + 1, 0, dim + 1);
    TH1F *fNoiseRMS = new TH1F("fNoiseRMS", "noise rms distribution", dim + 1, 0, dim + 1);
    fNoiseMean->GetXaxis()->SetTitle("board+chip+channel");
    fNoiseRMS->GetXaxis()->SetTitle("board+chip+channel");
    fNoiseMean->GetYaxis()->SetTitle("Pedestal mean");
    fNoiseRMS->GetYaxis()->SetTitle("Pedestal sigma");

    hNoise.resize(dim);
    for (int i = 0; i < boardNo; i++)
        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {
                TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", boardName[i], chipName[j], k));
                if (tmp != NULL)
                    delete tmp;
                hNoise[i * (64 * 4) + j * 64 + k] = new TH1D(Form("Ped_%d_%d_%d", boardName[i], chipName[j], k), Form("Pedestal for Board=%d, Chip=%d, Channel=%d", boardName[i], chipName[j], k), 1000, 0, 1000);
            }

    //read raw file
    TFile *rawFile = new TFile(fRawName);
    if (!rawFile->IsOpen())
    {
        cout << "#### Error: Can't open " << fRawName << ". Please check the filename." << endl;
        return;
    }
    TFile *pedFile = new TFile(fPedName, "recreate");
    cout << "--> Analyzing pedestal files from " << fRawName << endl;
    TTree *tree = (TTree *)rawFile->Get("tree");
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);

    int nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;

        tree->GetEntry(ii);

        //找出board对应的编号
        int iBoardIndex = -1;
        for (int j = 0; j < (int)boardName.size(); j++)
            if (boardName[j] == board)
                iBoardIndex = j;

        if (iBoardIndex == -1)
        {
            cout << "The board " << board << " does not exist in the current Track-AGET list, check your settings." << endl;
            break;
        }

        int index = iBoardIndex * (64 * 4) + (chip - 10) * 64 + channel;

        for (int j = 0; j < (int)wave->size(); j++)
        {
            hWave[index]->SetBinContent(j, hWave[index]->GetBinContent(j) + wave->at(j));
            if (j < 150 || j > 320)
                hNoise[index]->Fill(wave->at(j));
        }
    }

    for (int i = 0; i < boardNo; i++)
    {
        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {
                double mean = hNoise[i * (64 * 4) + j * 64 + k]->GetMean();
                double rms = hNoise[i * (64 * 4) + j * 64 + k]->GetRMS();
                int id = i * (64 * 4) + j * 64 + k;
                if (hNoise[id]->GetEntries() > 0)
                    hNoise[id]->Fit("gaus", "RQ", "", mean - 3 * rms, mean + 3 * rms);
                hNoise[id]->Write();
                hWave[id]->Write();
                fNoiseMean->SetBinContent(id, mean);
                fNoiseMean->SetBinError(id, rms);
                fNoiseRMS->SetBinContent(id, rms);

                for (int ii = 0; ii < 512; ii++)
                    gWaveAll[i]->SetPoint((j * 64 + k) * 521 + ii, j * 64 + k, ii, hWave[id]->GetBinContent(ii));
            }
        gWaveAll[i]->Write();
    }
    fNoiseMean->Write();
    fNoiseRMS->Write();

    pedFile->Write();
    pedFile->Close();
    rawFile->Close();
    cout << "--> Pedestal has been stored to " << fPedName << endl
         << endl;
}

void GenerateAGETPed(TString fRawName, TString fPedName, vector<MyBeamTestTrackAGET *> vTrkAGET)
{
    if (vTrkAGET.size() == 0)
        return;

    vector<int> boardName;
    vector<int> chipName; //假设所有的board的chiplist都相同

    chipName = vTrkAGET[0]->GetChipList();

    for (int i = 0; i < (int)vTrkAGET.size(); i++)
        for (int j = 0; j < (int)vTrkAGET[i]->GetBoardList().size(); j++)
            boardName.push_back(vTrkAGET[i]->GetBoardList()[j]);

    GenerateAGETPed(fRawName, fPedName, boardName, chipName, 1);
}

//-------------------------------
// Track-VMM Data to Raw-ROOT
class vmmTDCData
{
public:
    bool StoreData(vector<unsigned short> dat)
    {
        if (dat.size() != 16 || dat[0] != 0xAB || (dat[1] & 0xC0) != 0xC0)
            return false;

        trigID = ((dat[1] & 0x0F) << 12) | (dat[2] << 4) | ((dat[3] & 0xF0) >> 4);
        trigBCID = (dat[3] & 0x0F << 8) | dat[4];
        tdcValue = 1;
        status = 1;
        return true;
    }
    void Print()
    {
        if (!status)
            return;
        cout << "--> This is TDC data: triger ID = " << trigID << ", BCID = " << trigBCID << endl;
    }
    int status = 0;
    int trigID = 0;
    int trigBCID = 0;
    int tdcValue = 0;
};

class vmmTLUData
{
public:
    bool StoreData(vector<unsigned short> dat)
    {
        if (dat.size() != 16 || dat[0] != 0xFA || dat[1] != 0xFD || dat[14] != 0xAF || dat[15] != 0xDF)
            return false;

        trigID = (dat[2] << 8 | dat[3]);
        status = 1;
        return true;
    }
    void Print()
    {
        if (!status)
            return;
        cout << "--> This is TLU data: triger ID = " << trigID << endl;
    }
    int status = 0;
    int trigID = 0;
};

class vmmDetData
{
public:
    bool StoreData(vector<unsigned short> dat)
    {
        if (dat.size() != 16 || dat[0] != 0xEF || dat[1] != 0x0D || dat[15] != 0xFE)
            return false;

        board = dat[2] - 48;
        BCID = VMMExchange((dat[3] & 0x0F) << 8 | dat[4]);
        trigID = dat[5] << 8 | dat[6];
        chip = BitAttraction(dat[7], 3, 1);
        channel = BitAttraction2(dat[8], 6, 3);
        PD0 = Reverse((BitAttraction(dat[8], 2, 1) << 8) | dat[9], 10);
        TD0 = Reverse(dat[10], 8);
        status = 1;
        return true;
    }
    void Print()
    {
        if (!status)
            return;
        cout << "--> This is VMM data: triger ID = " << trigID << ", board=" << board << ", chip=" << chip << ", channel=" << channel << ", Q=" << PD0 << ", TD0=" << TD0 << ", BCID=" << BCID << endl;
    }
    int board = 0;
    int chip = 0;
    int channel = 0;
    int trigID = 0;
    int BCID = 0;
    int PD0 = 0;
    int TD0 = 0;
    int status = 0;
};

void ReadTrackVMMData2Root(vector<TString> datList, TString fRawName, int force = 1);
void ReadTrackVMMData2Root(vector<TString> datList, TString fRawName, int force)
{
    if ((int)datList.size() == 0)
        return;

    int type;
    int event;
    UShort_t board = 0;
    UShort_t chip = 0;
    UShort_t channel = 0;
    UShort_t PDO = 0;
    UShort_t BCID = 0;
    UShort_t TDO = 0;

    int length = sizeof(unsigned char);
    unsigned short memblock = 0;
    vector<unsigned short> memlist;

    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fRawName, fStat);
        if (fStat.fSize != 0)
            return;
    }

    TFile *fFile = new TFile(fRawName, "RECREATE");
    TTree *fTree = new TTree("tree", "tree");
    fTree->Branch("event", &event);
    fTree->Branch("board", &board);
    fTree->Branch("chip", &chip);
    fTree->Branch("channel", &channel);
    fTree->Branch("PDO", &PDO);
    fTree->Branch("BCID", &BCID);
    fTree->Branch("TDO", &TDO);

    TTree *fTree2 = new TTree("tree2", "tree2");
    fTree2->Branch("type", &type);
    type = VMMType;
    fTree2->Fill();
    fTree2->Write();

    TString fHead(datList[0]);
    fHead.ReplaceAll("1.bin", "");

    //int bufferFlag = 0;
    vector<vmmTDCData *> tdcList;
    vector<vmmTLUData *> trgList;
    vector<vmmDetData *> detList;

    for (int i = 0; i < (int)datList.size(); i++)
    {
        TString fName = fHead + Form("/%d.bin", i + 1);
        fstream InDat;
        InDat.open(fName, ios::in | ios::binary);
        cout << "--> Reading a new dat file-" << i << " : " << fName << endl;

        //int ntot = 0;
        while (!InDat.eof())
        {
            InDat.read((char *)(&memblock), length);

            //ntot++;
            //cout << dec << ntot << " reading: " << hex << memblock << dec << endl;

            if (memlist.size() < 16)
            {
                //cout << "    -> push" << endl;
                memlist.push_back(memblock);
                if (memlist.size() < 16)
                    continue;
            }
            else
            {
                //cout << "     -> erase" << endl;
                memlist.erase(memlist.begin());
                memlist.push_back(memblock);
            }

            vmmTDCData *tdcData = new vmmTDCData();
            vmmTLUData *trgData = new vmmTLUData();
            vmmDetData *detData = new vmmDetData();

            int stat = 0;
            stat += tdcData->StoreData(memlist);
            stat += trgData->StoreData(memlist);
            stat += detData->StoreData(memlist);

            if (stat != 1)
            {
                //cout << "      ";
                //for (int i = 0; i < memlist.size(); i++)
                //    cout << hex << memlist[i] << " ";
                //cout << "     not a good one: " << stat << endl;
                continue;
            }
            memlist.clear();

            if (tdcData->status)
            {
                tdcList.push_back(tdcData);
                //tdcData->Print();
            }

            if (trgData->status)
            {
                trgList.push_back(trgData);
                //trgData->Print();
            }

            if (detData->status)
            {
                //store data
                detList.push_back(detData);
                //detData->Print();
                event = detData->trigID;
                board = detData->board;
                chip = detData->chip;
                PDO = detData->PD0;
                BCID = detData->BCID;
                TDO = detData->TD0;
                fTree->Fill();
            }

            /*
            if (memblock == 0xEF)
            {
                counting = -2;
                headerFlag = 0;
            }
            else if (counting == -1)
            {
                if (memblock == 0x0D)
                    headerFlag = 1;
            }
            //else if (memblock >= 0x30 && memblock <= 0x37)
            else if (headerFlag == 1 && counting == 0)
            {
                board = memblock - 48;
            }
            else if (counting >= 1 && counting <= 2) // coarse time
            {
                if (counting == 1)
                    BCID = (memblock & 0xF);
                if (counting == 2)
                    BCID += memblock;
                BCID = VMMExchange(BCID);
            }
            else if (counting >= 3 && counting <= 4) // Event ID
            {
                if (counting == 3)
                {
                    event = memblock << 8;
                }
                if (counting == 4)
                {
                    event = event + memblock;
                }
            }
            else if (counting == 5) // VMM ID
            {
                //chip = Reverse(memblock, 8);
                chip = BitAttraction(memblock, 3, 1);
            }
            else if (counting == 6) // hit channel number
            {
                channel = BitAttraction2(memblock, 6, 3);
                PDO = BitAttraction(memblock, 2, 1);
                //cout << "2bit_PDO = " << PDO << "\t";
            }
            else if (counting == 7) // actually, the ADC
            {
                PDO = (PDO << 8) + memblock;
                //cout << "\t 10_bit_PDO = " << bitset<10>(PDO) << "\t ";
                PDO = Reverse(PDO, 10);
                //cout << "\t PDO = " << PDO << endl;
            }
            else if (counting == 8) // fine time -> TDC
            {
                TDO = Reverse(memblock, 8);
            }
            else if (counting >= 9 && counting <= 12) //
            {
                if (counting == 9)
                    buffer = memblock << 24;
                if (counting == 10)
                    buffer += memblock << 16;
                if (counting == 11)
                    buffer += memblock << 8;
                if (counting == 12)
                {
                    buffer += memblock;
                    //if (buffer == 0x80000000)
                    bufferFlag = 1;
                    //else
                    //    bufferFlag = 0;
                }
            }
            else if (counting == 13)
            {
                if (memblock == 0xFE)
                    trailerFlag = 1;
                else
                    trailerFlag = 0;
                //if (headerFlag == 1 && bufferFlag == 1 && trailerFlag == 1)
                //cout << "counting = " << counting << "\t buffer = " << hex << buffer << "\t memblock = "<< memblock << endl;
            }

            if (headerFlag == 1 && bufferFlag == 1 && trailerFlag == 1)
            {
                fTree->Fill();
                headerFlag = 0, bufferFlag = 0, trailerFlag = 0;
            }

            counting++;
            */
        }
    }
    cout << "--> Total entries = " << fTree->GetEntries() << endl;
    fTree->Write();
    //fTree->Print();
    fFile->Flush();
    fFile->Close();

    cout << "--> Track-VMM Raw data file has been convert to root-file: " << fRawName << ".\n"
         << endl;
}

#endif
