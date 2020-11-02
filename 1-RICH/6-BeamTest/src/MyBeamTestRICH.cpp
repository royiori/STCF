#include "TSystem.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TGeoVolume.h"
#include "MyBeamTestRICH.h"

//______________________________________________________________________________
//
MyBeamTestRICH::MyBeamTestRICH(int _id) : MyBeamTestDetector(_id)
{
    fPed = 0;
    fAllHit = 0;
    fCharge = 0;
}

MyBeamTestRICH::~MyBeamTestRICH()
{
}

//___ 0. 读取map文件 ___
//
void MyBeamTestRICH::ResizeMap()
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

void MyBeamTestRICH::ReadMapFile()
{
    fstream mapFp;
    TString mapfile = gSystem->WorkingDirectory() + TString("/") + mapping;
    mapFp.open(mapfile.Data(), ios::in);
    if (!mapFp.is_open())
        cout << "RICH Detector Map Open Failed: " << mapfile << "." << endl;

    string a;
    int board, chip, channel, posx, posy;
    getline(mapFp, a);
    ResizeMap();
    while (!mapFp.eof())
    {
        mapFp >> board >> chip >> channel >> posy >> posx;
        if (IndexChip(chip) == -1 || channel < 0 || channel >= 256)
            continue;

        if (IndexBD(board) != -1)
        {
            mapX[IndexBD(board)][IndexChip(chip)][channel] = 32 - posx;
            mapY[IndexBD(board)][IndexChip(chip)][channel] = posy;
            //cout<<board<<"("<<IndexBD(board)<<") "<<chip<<"("<<IndexChip(chip)<<") "<<channel<<" = ["<<posx<<", "<<posy<<"]"<<endl;
        }
    }

    mapFp.close();
}

//___ 1. 读取bin数据并转换为raw.root文件 ___
//
bool MyBeamTestRICH::ReadData2RawRoot(vector<TString> datList, TString fRawName, vector<MyBeamTestRICH *> vRICH, int force)
{
    if (datList.size() == 0)
        return false;

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
            return false;
    }

    TFile *fFile = new TFile(fRawName, "RECREATE");
    if (!fFile->IsOpen())
        return false;
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
    type = RICH;
    bdlist.clear();
    chlist.clear();
    for(int i=0; i<(int)vRICH.size(); i++)
    {
        for(int j=0; j<(int)vRICH[i]->boardList.size(); j++)
        {
            bdlist.push_back(vRICH[i]->boardList.at(j));
        }

        for(int j=0; j<(int)vRICH[i]->chipList.size(); j++)
        {
            chlist.push_back(vRICH[i]->chipList.at(j));
        }
    }
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
            }              // data packet header
            //else if (memblock&0xe == 0xe && headerFlag == 1 && counting == 0)	// board number, chip number, trigger type
            else if (counting == 1) // board number, chip number, trigger type
            {
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
            }
            else if (counting > 1 && counting < 6)
            {
            }
            else if (counting == 6)
            {
                channel = memblock; //cout << "channel number = " << channel << endl;
            }
            else if (counting > 6 && counting < 10)
            {
                if (counting == 7)
                    event += memblock * (2 ^ 32);
                else if (counting == 8)
                    event += memblock * (2 ^ 16);
                else if (counting == 9)
                {
                    event += memblock; //cout << "trigger Number = " << event << endl;
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
                { //cout << "wave.size() = " << wave.size() << endl;
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
    return true;
}

void MyBeamTestRICH::AnalysisPedestal(TString fRawName, TString fPedName, vector<int> boardName, vector<int> chipName, int force)
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

void MyBeamTestRICH::AnalysisPedestal(TString fRawName, TString fPedName, vector<MyBeamTestRICH *> vRICH)
{
    if (vRICH.size() == 0)
        return;

    vector<int> boardName;
    vector<int> chipName;

    boardName = vRICH[0]->GetBoardList();
    chipName = vRICH[0]->GetChipList();

    AnalysisPedestal(fRawName, fPedName, boardName, chipName, 1);
}

bool MyBeamTestRICH::ReadPedestal(TString fPedName, vector<MyBeamTestRICH *> vRICH)
{
    TFile *fPedFile = new TFile(fPedName);
    if (!fPedFile->IsOpen())
    {
        cout << "#### Can't open pedestal " << fPedName << " to read. Please check the path." << endl;
        return false;
    }
    cout << "--> Now reading the RICH pedestal file: " << fPedName << endl;

    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        vRICH[i]->SetNPedMean();
        vRICH[i]->SetNPedRMS();

        for (int ii = 0; ii < vRICH[i]->GetNBoard(); ii++)
            for (int jj = 0; jj < vRICH[i]->GetNChip(); jj++)
                for (int kk = 0; kk < 64; kk++)
                {
                    TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", vRICH[i]->GetBoardList()[ii], vRICH[i]->GetChipList()[jj], kk));
                    if (tmp == 0)
                    {
                        cout << "Can't not find the RICH pedestal file for: Ped_" << vRICH[i]->GetBoardList()[ii] << "_" << vRICH[i]->GetChipList()[jj] << "_" << kk << ", please check!" << endl;
                        continue;
                    }

                    double mean = 0;
                    double rms = 0;
                    if (tmp->GetEntries() > 0)
                    {
                        mean = tmp->GetFunction("gaus")->GetParameter(1);
                        rms = tmp->GetFunction("gaus")->GetParameter(2);
                    }
                    vRICH[i]->SetPedMean(ii, jj, kk, mean);
                    vRICH[i]->SetPedRMS(ii, jj, kk, rms);
                }
        vRICH[i]->GenPedMap();
    }
    cout << "--> RICH Pedestal file has been read.\n\n";
    return true;
}

//___ 2. 读取raw数据并转换为dst.root文件 ___
//
bool MyBeamTestRICH::ReadRaw2DstRoot(TString fRawName, TString fPedName, TString fDstPath, vector<MyBeamTestRICH *> vRICH, bool SaveWaveFlag)
{
    // 确认Raw-ROOT路径及文件存在
    TFile *fFile = new TFile(fRawName);
    if (!fFile->IsOpen())
        return false;

    // 确认ped-root路径及文件存在
    if (!ReadPedestal(fPedName, vRICH))
        return false;

    cout << "--> Now opening RICH raw root file: " << fRawName << endl;

    // 生成DST的文件路径，定义DST文件的ROOT格式
    vector<TFile *> fFileList;
    vector<TTree *> fTreeList;
    vector<MyBeamTestData *> fEventList;

    fFileList.resize(vRICH.size());
    fTreeList.resize(vRICH.size());
    fEventList.resize(vRICH.size());

    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        TString fDstName = fDstPath;
        fDstName += Form("%s-dst.root", vRICH[i]->GetName().Data());
        fFileList[i] = new TFile(fDstName, "recreate");
        fTreeList[i] = new TTree("tree", "beam test data structure");
        fEventList[i] = new MyBeamTestData;
        fTreeList[i]->Branch("event", "MyBeamTestData", &fEventList[i], 8000, 2);
    }

    // 读入Raw-ROOT文件，定义读取的格式
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

    // 读入Raw-ROOT文件
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
        for (int j = 0; j < (int)vRICH.size(); j++)
            for (int k = 0; k < vRICH[j]->GetNBoard(); k++)
                if (board == vRICH[j]->GetBoardList()[k])
                    id = j;

        if (id == -1)
        {
            cout << "#### Warning: event id " << event << "'s board: " << board << " doesn't exist in the RICH board list." << endl;
            continue;
        }

        //发现新的事例，则填tree
        if (fEventList[id]->event == -1 || fabs(fEventList[id]->event - event) != 0)
        {
            if (fEventList[id]->event >= 0 && fEventList[id]->event - event > 0)
                cout << "#### Warning: " << vRICH[id]->GetName() << "'s event id changed backwards from " << fEventList[id]->event << " to " << event << "\n\n";
            fEventList[id]->event = (fEventList[id]->event == -1) ? event : fEventList[id]->event;

            fFileList[id]->cd();
            fEventList[id]->Analysis(0, 150, 200, 300); //计算Q/T, 参数为ped和charge的计算所用的范围
            vRICH[id]->AnalysisCluster(fEventList[id]); //分析cluster

            if (!SaveWaveFlag) //如果不保存波形数据，则清空波形数组
                for (int ii = 0; ii < (int)fEventList[id]->board.size(); ii++)
                    fEventList[id]->wave[ii].clear();

            fTreeList[id]->Fill();
            fEventList[id]->Init(event, vRICH[id]->id);
        }

        //数据保存进vector
        fEventList[id]->board.push_back(board);
        fEventList[id]->chip.push_back(chip);
        fEventList[id]->channel.push_back(channel);
        fEventList[id]->pedeMean.push_back(vRICH[id]->GetPedMean(board, chip, channel));
        fEventList[id]->pedeRms.push_back(vRICH[id]->GetPedRMS(board, chip, channel));
        int nwave = fEventList[id]->wave.size();
        fEventList[id]->wave.resize(nwave + 1);
        for (int j = 0; j < (int)wave->size(); j++)
            fEventList[id]->wave[nwave].push_back(wave->at(j));
    }

    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        fFileList[i]->cd();
        fTreeList[i]->Write();

        int type;
        TTree *fTree2 = new TTree("tree2", "tree2");
        fTree2->Branch("type", &type);
        type = vRICH[i]->id;
        fTree2->Fill();
        fTree2->Write();

        fFileList[i]->Flush();
        fFileList[i]->Close();
    }
    fFile->Close();

    cout << "--> RICH Data file converted to a new structure." << endl
         << endl;
    return true;
}

void MyBeamTestRICH::AnalysisCluster(MyBeamTestData *fEvent)
{
    //转换map
    for (int i = 0; i < fEvent->GetHitSize(); i++)
    {
        pair<double, double> hpos = MapPosition(fEvent->board[i], fEvent->chip[i], fEvent->channel[i]);
        fEvent->hit.push_back(hpos);
    }

    //vector<vector<BeamHit>> branch; // [cluster][hit] 按击中hit的channel号分簇
    //找cluster-seed
    for (int i = 0; i < fEvent->GetHitSize(); i++)
    {
        double xhit = fEvent->hit[i].first;
        double yhit = fEvent->hit[i].second;

        BeamHit BHit;
        BHit.id = id;
        BHit.hit = fEvent->hit[i];
        BHit.q = fEvent->charge[i];
        BHit.t = fEvent->time[i];

        int checked = false; //是否保存进当前的branch
        int nbranch = fEvent->branch.size();
        for (int ii = 0; ii < nbranch; ii++)
        {
            int slen = fEvent->branch[ii].size();
            for (int jj = 0; jj < slen; jj++)
                if (fabs(xhit - fEvent->branch[ii][jj].hit.first) <= 1 && fabs(yhit - fEvent->branch[ii][jj].hit.second) <= 1)
                {
                    fEvent->branch[ii].push_back(BHit);
                    checked = true;
                    break;
                }
            if (checked)
                break;
        }

        if (!checked) //与所有branch都不匹配，则存一个新的branch
        {
            fEvent->branch.resize(fEvent->branch.size() + 1);
            fEvent->branch[fEvent->branch.size() - 1].push_back(BHit);
        }
    }

    /*
    cout << "-->hitmap: " << endl;
    for (int i = 0; i < fEvent->GetHitSize(); i++)
        cout << "(" << fEvent->hit[i].first << ", " << fEvent->hit[i].second << ")" << endl;

    cout << "Branch = " << fEvent->branch.size() << endl;
    for (int ii = 0; ii < (int)fEvent->branch.size(); ii++)
    {
        cout << ii << ":" << endl;
        for (int jj = 0; jj < (int)fEvent->branch[ii].size(); jj++)
        {
            cout << "(" << fEvent->branch[ii][jj].hit.first << ", " << fEvent->branch[ii][jj].hit.second << ") ";
        }
        cout << endl;
    }
    */

    //合并cluster
    if (fEvent->GetClusterSize() <= 1) //cluster数太少，不需要合并
        return;

    int nold = 0;
    int nnew = 1;
    while (nold != nnew)
    {
        int nbranch = fEvent->GetClusterSize();
        nold = nbranch;
        for (int i = 0; i < nbranch - 1; i++)
            for (int j = i + 1; j < nbranch; j++)
            {
                int mergeFlag = 0;
                int nbranch1 = fEvent->branch[i].size();
                int nbranch2 = fEvent->branch[j].size();
                for (int ii = 0; ii < nbranch1; ii++)
                    for (int jj = 0; jj < nbranch2; jj++)
                        if (fabs(fEvent->branch[i][ii].hit.first - fEvent->branch[j][jj].hit.first) <= 1 &&
                            fabs(fEvent->branch[i][ii].hit.second - fEvent->branch[j][jj].hit.second) <= 1)
                        {
                            mergeFlag = 1;
                            break;
                        }

                if (mergeFlag) //cluster_i和_j需要合并
                {
                    //cout << "cluster " << i << " " << j << " needs to be merged." << endl;
                    int nbranch1 = fEvent->branch[i].size();
                    int nbranch2 = fEvent->branch[j].size();
                    for (int jj = 0; jj < nbranch2; jj++)
                    {
                        int duplicateFlag = 0;
                        for (int ii = 0; ii < nbranch1; ii++)
                            if (fEvent->branch[i][ii].hit.first == fEvent->branch[j][jj].hit.first &&
                                fEvent->branch[i][ii].hit.second == fEvent->branch[j][jj].hit.second)
                            {
                                duplicateFlag = 1;
                                break;
                            }

                        if (duplicateFlag == 0)
                            fEvent->branch[i].push_back(fEvent->branch[j][jj]);
                    }

                    fEvent->branch[j].resize(0);

                    //cout << "\nremove the empty items." << endl;
                    int nbranch3 = (int)fEvent->branch.size() - 1;
                    for (int ii = nbranch3; ii >= 0; ii--)
                        if (fEvent->branch[ii].size() == 0)
                            fEvent->branch.erase(fEvent->branch.begin() + ii);
                }
            }
        nnew = fEvent->GetClusterSize();
    }

    /*
    cout << "\nafter mearge, branch = " << fEvent->branch.size() << endl;
    for (int ii = 0; ii < (int)fEvent->branch.size(); ii++)
    {
        cout << ii << ":" << endl;
        for (int jj = 0; jj < (int)fEvent->branch[ii].size(); jj++)
        {
            cout << "(" << fEvent->branch[ii][jj].hit.first << ", " << fEvent->branch[ii][jj].hit.second << ") ";
        }
        cout << endl;
    }
    cout << endl;
    */

    //每个cluster用重心法给出一个击中位置和q-sum的值
    /*
    double absPos[3];
    for (int ii = 0; ii < (int)branch.size(); ii++)
    {
        RealBeamHit bhit;
        bhit.id = id;
        bhit.hit[0] = 0;
        bhit.hit[1] = 0;
        bhit.hit[2] = 0;
        bhit.q = 0;

        for (int jj = 0; jj < (int)branch[ii].size(); jj++)
        {
            AbsPosition(branch[ii][jj].hit, absPos);
            bhit.hit[0] += absPos[0] * branch[ii][jj].q;
            bhit.hit[1] += absPos[1] * branch[ii][jj].q;
            bhit.hit[2] += absPos[2] * branch[ii][jj].q;
            bhit.q += branch[ii][jj].q;
        }

        bhit.hit[0] /= bhit.q;
        bhit.hit[1] /= bhit.q;
        bhit.hit[2] /= bhit.q;
        cluster.push_back(bhit);
    }
    */
}

//___ 3. 读取dst数据并合并为combine.root文件 ___
//
void MyBeamTestRICH::ReadDSTRoot(TString fPath)
{
    cout << "\n\n--> Reading RICH DST-root file:" << endl;

    TString fName = fPath;
    fName.ReplaceAll("RICH-raw", name + "-dst");

    fFile = new TFile(fName);
    if (!fFile->IsOpen())
        return;
    cout << "--> Opening " << fName << endl;

    fTree = (TTree *)fFile->Get("tree");
    fTree->SetBranchAddress("event", &fEvent);
    nentries = fTree->GetEntries();
}

void MyBeamTestRICH::CloseDSTRoot()
{
    if (fFile != NULL)
        fFile->Close();
    fTree = 0;
    fEvent = 0;
}

pair<double, double> MyBeamTestRICH::MapPosition(int board, int chip, int channel)
{
    if (IndexBD(board) == -1 || IndexChip(chip) == -1)
    {
        cout << "#### Warning: RICH Hit Board=" << board << " Chip=" << chip << " channel=" << channel << " is not correct." << endl;
        return make_pair(0, 0);
    }

    return make_pair(mapX[IndexBD(board)][IndexChip(chip)][channel], mapY[IndexBD(board)][IndexChip(chip)][channel]);
}

void MyBeamTestRICH::AbsPosition(pair<double, double> epos, double *absPos) //epos 是电子学的x/y的channel号定义，也就是MapPosition的输出
{
    double x0 = epos.first;
    double y0 = epos.second;
    double posX = -1 * ((x0 - 1) - Nstrp[0] / 2. + 1 / 2.) * strip[0]; //电子学mapping里的X/Y是这里的Y/X，从1到32
    double posY = ((y0 - 1) - Nstrp[1] / 2. + 1 / 2.) * strip[1];      //中心为0点

    //cout<<"Strip: ["<<strip[0]<<", "<<strip[1]<<"], mid:["<<posX<<", "<<posY<<"]"<<endl;

    //旋转
    absPos[0] = posX;
    absPos[1] = posY * sin(ROT[0]); //RICH与水平面（-Z方向）的夹角
    absPos[2] = -posY * cos(ROT[0]);

    //平移到RICH的轴心
    absPos[0] = pos[0] + absPos[0];
    absPos[1] = pos[1] + absPos[1] + LtoAnode * cos(ROT[0]);
    absPos[2] = pos[2] + absPos[2] + LtoAnode * sin(ROT[0]);
}
//______________________________________________________________________________
// 画三维击中
void MyBeamTestRICH::DrawHits(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med, MyBeamTestHitData *eventList)
{
    //TGeoVolume *anode = geom->MakeBox("anode", med, 1.0, siz[1], siz[2]);
    //anode->SetVisibility(kFALSE);
    TGeoVolume *pad = geom->MakeBox("pad", med, 1.0, strip[0] * 0.8 / 2, strip[0] * 0.8 / 2);
    pad->SetLineColorAlpha(kRed, 50);

    int ncopy = 1;
    for (int i = 0; i < (int)eventList->detector.size(); i++)
    {
        if (eventList->detector[i]->id != id)
            continue;

        double *geoPos, absPos[3];

        for (int j = 0; j < (int)eventList->detector[i]->hit.size(); j++)
        {
            AbsPosition(eventList->detector[i]->hit[j], absPos);
            geoPos = GeoPosition(absPos);
            //cout << name << "-" << i << ": GEOM:[" << geoPos[0] << "," << geoPos[1] << "," << geoPos[2] << "] for ABS:[" << absPos[0] << "," << absPos[1] << "," << absPos[2] << "]"<<endl;
            TGeoTranslation *gtr = new TGeoTranslation(geoPos[0], geoPos[1], geoPos[2]);
            top->AddNode(pad, ncopy++, gtr);
        }
    }

    /*
    TGeoTranslation *gtr = new TGeoTranslation(pos[0], pos[1], pos[2]);
    TGeoRotation *grot = new TGeoRotation("rot", rot[0], rot[1], 0.);
    TGeoCombiTrans *comb = new TGeoCombiTrans(*gtr, *grot);
    top->AddNode(anode, 1, 0); //comb);
    */
}

void MyBeamTestRICH::DrawDetector(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med)
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
bool MyBeamTestRICH::FillDistribution(int i, MyBeamTestHitData *eventList)
{
    int hitted = false;
    if (i == 0)
    {
        fAllHit = new TH2F(Form("AHitRH%d", id), Form("All hits for %s", name.Data()), 32, 0, 32, 32, 0, 32);
        fAllHit2 = new TH2F(Form("AHitRH%d_2", id), Form("All hits for %s", name.Data()), 36, -90 + pos[0], 90 + pos[0], 36, -90 + pos[1], 90 + pos[1]);
        fPed = new TH2F(Form("fPedRH%d_2", id), Form("Pedestal for %s", name.Data()), 32, 0, 32, 32, 0, 32);
        fCharge = new TH1F(Form("ChargRH%d", id), Form("Charge distribution of %s", name.Data()), 200, 0, 10000);
        fAllHit->SetXTitle("channel");
        fAllHit->SetYTitle("channel");
        fAllHit2->SetXTitle("mm");
        fAllHit2->SetYTitle("mm");
        fCharge->SetXTitle("q");
        fCharge->SetYTitle("Entries");
    }

    double absPos[3];
    double sum = 0;
    for (int i = 0; i < (int)eventList->detector.size(); i++)
    {
        if (eventList->detector[i]->id != id)
            continue;

        for (int j = 0; j < (int)eventList->detector[i]->hit.size(); j++)
        {
            AbsPosition(eventList->detector[i]->hit[j], absPos);
            fAllHit->Fill(eventList->detector[i]->hit[j].first, eventList->detector[i]->hit[j].second);
            fAllHit2->Fill(absPos[0], absPos[1]);
            if (absPos[1] > 0)
                sum += eventList->detector[i]->charge[j];
            hitted = true;
        }
        fCharge->Fill(sum);
    }
    return hitted;
}
