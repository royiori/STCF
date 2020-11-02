#include "TSystem.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "TGeoVolume.h"
#include "MyBeamTestTrackAGET.h"

//______________________________________________________________________________
MyBeamTestTrackAGET::MyBeamTestTrackAGET(int _id) : MyBeamTestDetector(_id)
{
    fAllHit = 0;
    fHit[0] = 0;
    fHit[1] = 0;
    fPed[0] = 0;
    fPed[1] = 0;
    fCharge = 0;
}

MyBeamTestTrackAGET::~MyBeamTestTrackAGET()
{
}

//___ 0. 读取map文件 ___
//
void MyBeamTestTrackAGET::ResizeMap()
{
    mapX.resize(boardList.size());
    for (int i = 0; i < (int)boardList.size(); i++)
    {
        mapX[i].resize(chipList.size());
        for (int j = 0; j < (int)chipList.size(); j++)
            mapX[i][j].resize(64);
    }

    mapY.resize(boardList.size());
    for (int i = 0; i < (int)boardList.size(); i++)
    {
        mapY[i].resize(chipList.size());
        for (int j = 0; j < (int)chipList.size(); j++)
            mapY[i][j].resize(64);
    }
}

void MyBeamTestTrackAGET::ReadMapFile()
{
    fstream mapFp;

    TString mapfile = mapping;
    mapFp.open(mapfile.Data(), ios::in);
    if (!mapFp.is_open())
        cout << "Track-AGET Detector Map Open Failed: " << mapfile << "." << endl;

    string a;
    int board, chip, channel, connector, pos;
    getline(mapFp, a);
    ResizeMap();
    while (!mapFp.eof())
    {
        mapFp >> board >> chip >> channel >> connector >> pos;
        if (IndexChip(chip) == -1 || channel < 0 || channel >= 256)
            continue;

        if (IndexBD(board) != -1)
        {
            mapX[IndexBD(board)][IndexChip(chip)][channel - IndexChip(chip) * 64] = (pos - 6);
            //cout<<board<<"("<<IndexBD(board)<<") "<<chip<<"("<<IndexChip(chip)<<") "<<channel<<"("<<channel-IndexChip(chip)*64<<") = "<<pos<<"("<<pos-6<<")"<<endl;
        }
    }

    mapFp.close();
}

//___ 1. 读取bin数据并转换为raw.root文件 ___
//
void MyBeamTestTrackAGET::ReadData2RawRoot(vector<TString> datList, TString fRawName, vector<MyBeamTestTrackAGET *> vTrkAGET, int force)
{
    if (datList.size() == 0)
        return;

    int length = sizeof(unsigned short);
    unsigned short memblock;
    int boardstart = 0;
    int chipstart = 0;
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

    int type;
    vector<int> bdlist;
    vector<int> chlist;
    TTree *fTree2 = new TTree("tree2", "tree2");
    fTree2->Branch("type", &type);
    fTree2->Branch("bdlist", &bdlist);
    fTree2->Branch("chlist", &chlist);
    type = TrackerAGET;
    bdlist.clear();
    chlist.clear();
    for (int i = 0; i < (int)vTrkAGET.size(); i++)
    {
        for (int j = 0; j < (int)vTrkAGET[i]->boardList.size(); j++)
        {
            bdlist.push_back(vTrkAGET[i]->boardList.at(j));
        }
    }
    for (int j = 0; j < (int)vTrkAGET[0]->chipList.size(); j++)
    {
        chlist.push_back(vTrkAGET[0]->chipList.at(j));
    } //#### Tracker的chip都是一样的
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

            //exchange the low 8 & high 8
            unsigned short aa = memblock;
            unsigned short low = memblock & 0xff;
            unsigned short high = (memblock >> 8) & 0xFF;
            aa = low;
            aa <<= 8;
            aa |= high;
            memblock = aa;

            //
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

void MyBeamTestTrackAGET::AnalysisPedestal(TString fRawName, TString fPedName, vector<int> boardName, vector<int> chipName, int force)
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
    cout << "--> Pedestal has been stored to " << fPedName << "\n\n";
}

void MyBeamTestTrackAGET::AnalysisPedestal(TString fRawName, TString fPedName, vector<MyBeamTestTrackAGET *> vTrkAGET)
{
    if (vTrkAGET.size() == 0)
        return;

    vector<int> boardName;
    vector<int> chipName; //假设所有的board的chiplist都相同

    chipName = vTrkAGET[0]->GetChipList();

    for (int i = 0; i < (int)vTrkAGET.size(); i++)
        for (int j = 0; j < (int)vTrkAGET[i]->GetBoardList().size(); j++)
            boardName.push_back(vTrkAGET[i]->GetBoardList()[j]);

    AnalysisPedestal(fRawName, fPedName, boardName, chipName, 1);
}

bool MyBeamTestTrackAGET::ReadPedestal(TString fPedName, vector<MyBeamTestTrackAGET *> vTrkAGT)
{
    TFile *fPedFile = new TFile(fPedName);
    if (!fPedFile->IsOpen())
    {
        cout << "#### Can't open pedestal " << fPedName << " to read. Please check the path." << endl;
        return false;
    }
    cout << "--> Now reading the Track-AGET pedestal file: " << fPedName << endl;

    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        vTrkAGT[i]->SetNPedMean();
        vTrkAGT[i]->SetNPedRMS();

        for (int ii = 0; ii < vTrkAGT[i]->GetNBoard(); ii++)
            for (int jj = 0; jj < vTrkAGT[i]->GetNChip(); jj++)
                for (int kk = 0; kk < 64; kk++)
                {
                    TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", vTrkAGT[i]->GetBoardList()[ii], vTrkAGT[i]->GetChipList()[jj], kk));
                    if (tmp == 0)
                    {
                        cout << "Check the TrackAGET pedestal file: Ped_" << vTrkAGT[i]->GetBoardList()[ii] << "_" << vTrkAGT[i]->GetChipList()[jj] << "_" << kk << " exists or not" << endl;
                        continue;
                    }
                    double mean = 0;
                    double rms = 0;
                    if (tmp->GetEntries() > 0)
                    {
                        mean = tmp->GetFunction("gaus")->GetParameter(1);
                        rms = tmp->GetFunction("gaus")->GetParameter(2);
                    }
                    vTrkAGT[i]->SetPedMean(ii, jj, kk, mean);
                    vTrkAGT[i]->SetPedRMS(ii, jj, kk, rms);
                }
        vTrkAGT[i]->GenPedMap();
    }
    cout << "--> TrkAGET Pedestal file has been read.\n\n";
    return true;
}

//___ 2. 读取raw数据并转换为dst.root文件 ___
//
bool MyBeamTestTrackAGET::ReadRaw2DstRoot(TString fRawName, TString fPedName, TString fDstPath, vector<MyBeamTestTrackAGET *> vTrkAGT, bool SaveWaveFlag)
{
    // 确认Raw-ROOT路径及文件存在
    TFile *fFile = new TFile(fRawName);
    if (!fFile->IsOpen())
        return false;

    // 确认ped-root路径及文件存在
    if (!ReadPedestal(fPedName, vTrkAGT))
        return false;

    cout << "--> Now opening Track-AGET raw root file: " << fRawName << endl;

    // 生成DST的文件路径, 三个track会生成三个-dst.root文件
    vector<TFile *> fFileList;
    vector<TTree *> fTreeList;
    vector<MyBeamTestData *> fEventList;

    fFileList.resize(vTrkAGT.size());
    fTreeList.resize(vTrkAGT.size());
    fEventList.resize(vTrkAGT.size());

    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        TString fDstName = fDstPath;
        fDstName += Form("%s-dst.root", vTrkAGT[i]->GetName().Data());
        fFileList[i] = new TFile(fDstName, "recreate");
        fTreeList[i] = new TTree("tree", "beam test data structure");
        fEventList[i] = new MyBeamTestData;
        fTreeList[i]->Branch("event", "MyBeamTestData", &fEventList[i], 8000, 2);
    }

    // 读入Raw-ROOT文件
    TTree *tree = (TTree *)fFile->Get("tree");
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;
    vector<double> *wave = 0;

    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);

    int id = -1;
    long nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;
        tree->GetEntry(ii);

        //根据board判断这个数据属于哪个探测器，放入id中。
        id = -1;
        for (int j = 0; j < (int)vTrkAGT.size(); j++)
            for (int k = 0; k < vTrkAGT[j]->GetNBoard(); k++)
                if (board == vTrkAGT[j]->GetBoardList()[k])
                    id = j;

        if (id == -1)
        {
            cout << "#### Warning: event id " << event << "'s board: " << board << " doesn't exist in the Track-AGET board list." << endl;
            continue;
        }

        //判断这个事例所属于的探测器
        fEventList[id]->event = (fEventList[id]->event == -1) ? event : fEventList[id]->event;
        if ((fEventList[id]->event) - event > 0)
            cout << "#### Warning: " << vTrkAGT[id]->GetName() << "'s event id changed backwards from " << fEventList[id]->event << " to " << event << "\n\n";

        //发现新的事例，则填tree
        if (fabs(fEventList[id]->event - event) != 0)
        {
            //cout<<"Store a new event: "<<fEventList[id]->event<<endl;
            fFileList[id]->cd();
            fEventList[id]->Analysis(0, 150, 200, 300);   //计算Q/T
            vTrkAGT[id]->AnalysisCluster(fEventList[id]); //分析cluster

            if (!SaveWaveFlag) //如果不保存波形数据，则清空波形数组
                for (int ii = 0; ii < (int)fEventList[id]->board.size(); ii++)
                    fEventList[id]->wave[ii].clear();
            fTreeList[id]->Fill();
            fEventList[id]->Init(event, vTrkAGT[id]->id);
        }

        //数据保存进vector
        fEventList[id]->board.push_back(board);
        fEventList[id]->chip.push_back(chip);
        fEventList[id]->channel.push_back(channel);
        fEventList[id]->pedeMean.push_back(vTrkAGT[id]->GetPedMean(board, chip, channel));
        fEventList[id]->pedeRms.push_back(vTrkAGT[id]->GetPedRMS(board, chip, channel));
        int nwave = fEventList[id]->wave.size();
        fEventList[id]->wave.resize(nwave + 1);
        for (int j = 0; j < (int)wave->size(); j++)
            fEventList[id]->wave[nwave].push_back(wave->at(j));
    }

    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        fTreeList[i]->Fill();
        fFileList[i]->cd();
        fTreeList[i]->Write();

        int type;
        TTree *fTree2 = new TTree("tree2", "tree2");
        fTree2->Branch("type", &type);
        type = vTrkAGT[i]->id;
        fTree2->Fill();
        fTree2->Write();

        fFileList[i]->Flush();
        fFileList[i]->Close();
    }
    fFile->Close();

    cout << "--> Track-AGET Data file converted to a new structure." << endl
         << endl;
    return true;
}

void MyBeamTestTrackAGET::AnalysisCluster(MyBeamTestData *fEvent)
{
    //转换map
    for (int i = 0; i < fEvent->GetHitSize(); i++)
    {
        pair<double, double> hpos = MapPosition(fEvent->board[i], fEvent->chip[i], fEvent->channel[i]);
        fEvent->hit.push_back(hpos);
    }

    //vector<vector<BeamHit>> branch; // [cluster][hit] 按击中hit的channel号分簇
    //先排序
    vector<BeamHit> xBHit;
    vector<BeamHit> yBHit;

    for (int i = 0; i < 128; i++)
    {
        for (int j = 0; j < fEvent->GetHitSize(); j++)
        {
            BeamHit BHit;
            BHit.id = id;
            BHit.hit = fEvent->hit[j];
            BHit.q = fEvent->charge[j];
            BHit.t = fEvent->time[j];

            if (BHit.hit.first == i) //Y读出条
                xBHit.push_back(BHit);

            if (BHit.hit.second == i) //X读出条
                yBHit.push_back(BHit);
        }
    }

    /*
    cout << "-->hitmap: " << endl;
    for (int i = 0; i < fEvent->GetHitSize(); i++)
        cout << "(" << fEvent->hit[i].first << ", " << fEvent->hit[i].second << ")" << endl;

    cout << "Sorted-X = " << xBHit.size() << endl;
    for (int ii = 0; ii < (int)xBHit.size(); ii++)
        cout << "(" << xBHit[ii].hit.first << ", " << xBHit[ii].hit.second << ")" << endl;

    cout << "Sorted-Y = " << yBHit.size() << endl;
    for (int ii = 0; ii < (int)yBHit.size(); ii++)
        cout << "(" << yBHit[ii].hit.first << ", " << yBHit[ii].hit.second << ")" << endl;
    */

    //找X方向的cluster，形状为X长Y短的条，X为-999，给出Y坐标
    for (int i = 0; i < (int)xBHit.size(); i++)
    {
        int checked = false; //是否保存进当前的branch
        int nbranch = fEvent->XYbranch[0].size();
        for (int j = 0; j < nbranch; j++)
        {
            int slen = fEvent->XYbranch[0][j].size();
            for (int jj = 0; jj < slen; jj++)
                if (fabs(xBHit[i].hit.first - fEvent->XYbranch[0][j][jj].hit.first) <= 1)
                {
                    fEvent->XYbranch[0][j].push_back(xBHit[i]);
                    checked = true;
                    break;
                }
            if (checked)
                break;
        }

        if (!checked)
        {
            fEvent->XYbranch[0].resize(fEvent->XYbranch[0].size() + 1);
            fEvent->XYbranch[0][fEvent->XYbranch[0].size() - 1].push_back(xBHit[i]);
        }
    }

    //找Y方向的cluster，形状为X长Y短的条，X为-999，给出Y坐标
    for (int i = 0; i < (int)yBHit.size(); i++)
    {
        int checked = false; //是否保存进当前的branch
        int nbranch = fEvent->XYbranch[1].size();
        for (int j = 0; j < nbranch; j++)
        {
            int slen = fEvent->XYbranch[1][j].size();
            for (int jj = 0; jj < slen; jj++)
                if (fabs(yBHit[i].hit.second - fEvent->XYbranch[1][j][jj].hit.second) <= 1)
                {
                    fEvent->XYbranch[1][j].push_back(yBHit[i]);
                    checked = true;
                    break;
                }
            if (checked)
                break;
        }

        if (!checked)
        {
            fEvent->XYbranch[1].resize(fEvent->XYbranch[1].size() + 1);
            fEvent->XYbranch[1][fEvent->XYbranch[1].size() - 1].push_back(yBHit[i]);
        }
    }

    /*
    cout << "\nafter mearge, XYbranch[0] = " << fEvent->XYbranch[0].size() << endl;
    for (int ii = 0; ii < (int)fEvent->XYbranch[0].size(); ii++)
    {
        cout << ii << ":" << endl;
        for (int jj = 0; jj < (int)fEvent->XYbranch[0][ii].size(); jj++)
        {
            cout << "(" << fEvent->XYbranch[0][ii][jj].hit.first << ", " << fEvent->XYbranch[0][ii][jj].hit.second << ") ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "\nafter mearge, Ybranch = " << fEvent->XYbranch[1].size() << endl;
    for (int ii = 0; ii < (int)fEvent->XYbranch[1].size(); ii++)
    {
        cout << ii << ":" << endl;
        for (int jj = 0; jj < (int)fEvent->XYbranch[1][ii].size(); jj++)
        {
            cout << "(" << fEvent->XYbranch[1][ii][jj].hit.first << ", " << fEvent->XYbranch[1][ii][jj].hit.second << ") ";
        }
        cout << endl;
    }
    cout << endl;
    */

    /*
    TH1F *ftmp = new TH1F("ftmp", "ftmp", 128, 0, 128);
    RealBeamHit hittmp;
    hittmp.id = id;
    for (int j = 0; j < (int)fEvent->XYbranch[0].size(); j++)
    {
        ftmp->Reset();
        double qsum = 0;
        for (int k = 0; k < (int)fEvent->XYbranch[0][j].size(); k++)
        {
            qsum += fEvent->XYbranch[0][j][k].q;
            ftmp->Fill(fEvent->XYbranch[0][j][k].hit.first, fEvent->XYbranch[0][j][k].q);
        }
        hittmp.q = qsum;
        hittmp.nhit = (int)fEvent->XYbranch[0][j].size();
        hittmp.hit[0] = ftmp->GetMean();
        hittmp.hit[1] = -999;
        hittmp.hit[2] = -999;
        hittmp.hiterr[0] = ftmp->GetRMS();
        hittmp.hiterr[1] = -999;
        hittmp.hiterr[2] = -999;
        fEvent->Xcluster.push_back(hittmp);
    }

    for (int j = 0; j < (int)fEvent->Ybranch.size(); j++)
    {
        ftmp->Reset();
        double qsum = 0;
        for (int k = 0; k < (int)fEvent->Ybranch[j].size(); k++)
        {
            qsum += fEvent->Ybranch[j][k].q;
            ftmp->Fill(fEvent->Ybranch[j][k].hit.second, fEvent->Ybranch[j][k].q);
        }
        hittmp.q = qsum;
        hittmp.nhit = (int)fEvent->Ybranch[j].size();
        hittmp.hit[0] = -999;
        hittmp.hit[1] = ftmp->GetMean();
        hittmp.hit[2] = -999;
        hittmp.hiterr[0] = -999;
        hittmp.hiterr[1] = ftmp->GetRMS();
        hittmp.hiterr[2] = -999;
        fEvent->Ycluster.push_back(hittmp);
    }
    delete ftmp;
    */
}

//___ 3. 读取dst数据并合并为combine.root文件 ___
//
void MyBeamTestTrackAGET::ReadDSTRoot(TString fPath)
{
    cout << "\n\n--> Reading Track-AGET DST-root file:" << endl;

    TString fName = fPath;
    fName.ReplaceAll("TrackAGET-raw", name + "-dst");

    fFile = new TFile(fName);
    if (!fFile->IsOpen())
        return;
    cout << "--> Opening " << fName << endl;

    fTree = (TTree *)fFile->Get("tree");
    fTree->SetBranchAddress("event", &fEvent);
    nentries = fTree->GetEntries();
}

void MyBeamTestTrackAGET::CloseDSTRoot()
{
    if (fFile != NULL)
        fFile->Close();
    fTree = 0;
    fEvent = 0;
}

pair<double, double> MyBeamTestTrackAGET::MapPosition(int board, int chip, int channel)
{
    if (IndexBD(board) == -1 || IndexChip(chip) == -1 || channel >= 64 || channel < 0)
    {
        cout << "#### TrackAGET Warning: Hit Board=" << board << " Chip=" << chip << " channel=" << channel << " is not correct." << endl;
        return make_pair(0, 0);
    }

    // 这里的坐标是全局坐标的X和Y。
    // Chip<2是全局坐标形状为Y的条（Y长X短）的条，给出的是X坐标，长方向的长度设为-999, position从1开始
    if (IndexChip(chip) < 2)
        return make_pair(mapX[IndexBD(board)][IndexChip(chip)][channel], -999);

    // Chip>=2是全局坐标形状为X的条（X长Y短），给出的是Y坐标，position从1开始
    return make_pair(-999, mapX[IndexBD(board)][IndexChip(chip)][channel]);
}

void MyBeamTestTrackAGET::AbsPosition(pair<double, double> epos, double *absPos) //epos 是电子学的x/y的channel号定义，也就是MapPosition的输出
{
    double x0 = (epos.first == -999) ? 0 : epos.first;
    double y0 = (epos.second == -999) ? 0 : epos.second;
    double posX = ((x0 - 1) - Nstrp[0] / 2. + 1 / 2.) * strip[0]; //电子学mapping里的X/Y是这里的Y/X，从1到32
    double posY = ((y0 - 1) - Nstrp[1] / 2. + 1 / 2.) * strip[1]; //中心为0点
    posX = (epos.first == -999) ? 0 : posX;
    posY = (epos.second == -999) ? 0 : posY;
    //cout<<"Strip: ["<<strip[0]<<", "<<strip[1]<<"], mid:["<<posX<<", "<<posY<<"]"<<endl;

    //旋转
    absPos[0] = posX;
    absPos[1] = posY * sin(ROT[0]); //RICH与水平面（-Z方向）的夹角
    absPos[2] = -posY * cos(ROT[0]);
    //cout<<"absPos1: "<<absPos[0]<<" "<<absPos[1]<<" "<<absPos[2]<<endl;

    //平移到探测器的轴心
    absPos[0] = pos[0] + absPos[0];
    absPos[1] = pos[1] + absPos[1] + LtoAnode * cos(ROT[0]);
    absPos[2] = pos[2] + absPos[2] + LtoAnode * sin(ROT[0]);
    //cout<<"absPos2: "<<absPos[0]<<" "<<absPos[1]<<" "<<absPos[2]<<endl;
}
//______________________________________________________________________________
// 画三维击中
void MyBeamTestTrackAGET::DrawHits(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med, MyBeamTestHitData *eventList)
{
    //Xpad是平行X方向的读出条, X方向为长方向(对应Geom是Z为长方向)
    //Ypad是平行Y方向的读出条，Y方向为长方形(对应Geom是Y为长方向)
    TGeoVolume *Xpad = geom->MakeBox("pad", med, 1.0, strip[0] * 0.8 / 2, Nstrp[1] * strip[1] * 0.8 / 2);
    TGeoVolume *Ypad = geom->MakeBox("pad", med, 1.0, Nstrp[0] * strip[0] * 0.8 / 2, strip[1] * 0.8 / 2);
    Xpad->SetLineColorAlpha(kBlue, 50);
    Ypad->SetLineColorAlpha(kRed, 50);

    int xcopy = 0;
    int ycopy = 0;
    for (int i = 0; i < (int)eventList->detector.size(); i++)
    {
        if (eventList->detector[i]->id != id)
            continue;

        double *geoPos, absPos[3];
        for (int j = 0; j < (int)eventList->detector[i]->hit.size(); j++)
        {
            AbsPosition(eventList->detector[i]->hit[j], absPos);
            geoPos = GeoPosition(absPos);
            TGeoTranslation *gtr = new TGeoTranslation(geoPos[0], geoPos[1], geoPos[2]);

            if (eventList->detector[i]->hit[j].first == -999) // X-pad, X为长方向，长度设为-999
                top->AddNode(Xpad, xcopy++, gtr);
            else
                top->AddNode(Ypad, ycopy++, gtr);
            //cout << name << "-" << i << ": GEOM:[" << geoPos[0] << "," << geoPos[1] << "," << geoPos[2] << "] for ABS:[" << absPos[0] << "," << absPos[1] << "," << absPos[2] << "]"<<endl;
        }
    }
}

void MyBeamTestTrackAGET::DrawDetector(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med)
{
    TGeoVolume *box = geom->MakeBox(name, med, siz[2] / 2., siz[1] / 2., siz[0] / 2.);
    TGeoTranslation *tr = new TGeoTranslation(pos[2], -pos[1], -pos[0]);
    TGeoRotation *rott = new TGeoRotation(name, rot[0] - 90, rot[1], 0.);
    TGeoCombiTrans *comb = new TGeoCombiTrans(*tr, *rott);
    box->SetLineColorAlpha(kBlack, 50);
    top->AddNode(box, 1, comb);
}

//______________________________________________________________________________
// 生成分布图
bool MyBeamTestTrackAGET::FillDistribution(int i, MyBeamTestHitData *eventList)
{
    bool hitted = false;

    if (i == 0)
    {
        fAllHit = new TH2F(Form("AHitTrk%d", id), Form("All hits for %s", name.Data()), 128, 0, 128, 128, 0, 128);
        fAllHit2 = new TH2F(Form("AHitTrk_2%d", id), Form("All hits for %s", name.Data()), 150, -30 + pos[0], 30 + pos[0], 150, -30 + pos[1], 30 + pos[1]);
        fHit[0] = new TH1F(Form("HitXTrk%d", id), "X Strip hits for single event", 128, 0, 128);
        fHit[1] = new TH1F(Form("HitYTrk%d", id), "Y Strip hits for single event", 128, 0, 128);
        fPed[0] = new TH1F(Form("PedXTrk%d", id), "X Strip pedestal", 128, 0, 128);
        fPed[1] = new TH1F(Form("PedYTrk%d", id), "Y Strip pedestal", 128, 0, 128);
        fCharge = new TH1F(Form("ChargTrk%d", id), Form("Charge distribution of %s", name.Data()), 200, 0, 10000);

        fAllHit->SetXTitle("channel");
        fAllHit->SetYTitle("channel");
        fAllHit2->SetXTitle("mm");
        fAllHit2->SetYTitle("mm");
        fCharge->SetXTitle("q");
        fCharge->SetYTitle("Entries");
        fHit[0]->SetXTitle("channel");
        fHit[1]->SetXTitle("channel");
        fHit[0]->SetYTitle("Entries");
        fHit[1]->SetYTitle("Entries");
    }

    vector<pair<double, double>> x;
    vector<pair<double, double>> y;

    bool hittedX = false;
    bool hittedY = false;
    double sum = 0;
    for (int i = 0; i < (int)eventList->detector.size(); i++)
    {
        if (eventList->detector[i]->id != id)
            continue;

        for (int j = 0; j < (int)eventList->detector[i]->hit.size(); j++)
        {
            if (eventList->detector[i]->hit[j].first == -999) // X-pad, X为长方向，长度设为-999
            {
                y.push_back(eventList->detector[i]->hit[j]);
                fHit[1]->Fill(eventList->detector[i]->hit[j].second);
                hittedY = true;
            }
            else // Y-pad, Y为长方向，长度设为-999, 提供X坐标
            {
                x.push_back(eventList->detector[i]->hit[j]);
                fHit[0]->Fill(eventList->detector[i]->hit[j].first);
                hittedX = true;
            }
            sum += eventList->detector[i]->charge[i];
        }
        fCharge->Fill(sum);

        if (hittedX && hittedY)
            hitted = true;

        double absPos[3];
        for (int i = 0; i < (int)x.size(); i++)
            for (int j = 0; j < (int)y.size(); j++)
            {
                fAllHit->Fill(x[i].first, y[j].second);

                AbsPosition(x[i], absPos);
                double _x = absPos[0];
                AbsPosition(y[j], absPos);
                double _y = absPos[1];

                //cout << name << " " << _x << " " << _y << endl;
                fAllHit2->Fill(_x, _y); //gRandom->Gaus(0, 3), gRandom->Gaus(0, 3));
            }
    }
    return hitted;
}
