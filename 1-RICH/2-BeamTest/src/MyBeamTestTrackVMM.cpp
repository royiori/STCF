#include "TGeoVolume.h"
#include "MyBeamTestTrackVMM.h"

//______________________________________________________________________________
MyBeamTestTrackVMM::MyBeamTestTrackVMM(int _id, int nx, int ny, double wx, double wy) : MyBeamTestDetector(_id, nx, ny, wx, wy)
{
}

MyBeamTestTrackVMM::~MyBeamTestTrackVMM()
{
}

//______________________________________________________________________________
// 读取map文件
void MyBeamTestTrackVMM::ReadMapFile()
{
    fstream mapFp;
    mapFp.open("./data/Map/Tracker_VMM_Map.txt", ios::in);
    if (!mapFp.is_open())
        cout << "Track-VMM Detector Map Open Failed! " << endl;

    string a;
    int board, chip, channel, pos;
    getline(mapFp, a);
    ResizeMap();
    while (!mapFp.eof())
    {
        mapFp >> board >> chip >> channel >> pos;
        if (IndexChip(chip) == -1 || channel < 0 || channel >= 64)
            continue;

        if (IndexBD(board) != -1)
        {
            pos -= 6;
            if (pos < 64)
                pos += 64;
            else
                pos -= 64;
            mapX[IndexBD(board)][IndexChip(chip)][channel] = (pos - 6);
            //cout << name << " " << board << "(" << IndexBD(board) << ") " << chip << "(" << IndexChip(chip) << ") " << channel << "(" << channel << ") = " << pos << "(" << pos - 6 << ")" << endl;
        }
    }

    mapFp.close();
}

//______________________________________________________________________________
// Combine数据
void MyBeamTestTrackVMM::ReadDSTRoot(TString fPath)
{
    cout << "\n\n--> Reading Track-VMM DST-root file:" << endl;

    TString fName = fPath;
    fName.ReplaceAll("TrackVMM-raw", name + "-dst");

    fFile = new TFile(fName);
    if (!fFile->IsOpen())
        return;
    cout << "--> Opening " << fName << endl;

    fTree = (TTree *)fFile->Get("tree");
    fTree->SetBranchAddress("event", &fEvent);
    nentries = fTree->GetEntries();
}

void MyBeamTestTrackVMM::CloseDSTRoot()
{
    if (fFile != NULL)
        fFile->Close();
    fTree = 0;
    fEvent = 0;
}

pair<double, double> MyBeamTestTrackVMM::MapPosition(int board, int chip, int channel)
{
    if (IndexBD(board) == -1 || IndexChip(chip) == -1 || channel >= 64 || channel < 0)
    {
        cout << "#### TrackVMM Warning: Hit Board=" << board << " Chip=" << chip << " channel=" << channel << " is not correct." << endl;
        return make_pair(0, 0);
    }

    //cout<<IndexBD(board)<<" "<<IndexChip(chip)<<" "<<channel<<" "<<mapX.size()<<" "<<mapX[0].size()<<" "<<mapX[0][0].size()<<endl;
    // 这里的坐标是全局坐标的X和Y。
    // board=0是全局坐标形状为Y的条（Y长X短）的条，给出的是X坐标，长方向的长度设为-999, position从1开始
    if (IndexBD(board) == 0)
        return make_pair(mapX[IndexBD(board)][IndexChip(chip)][channel], -999);

    // board=1是全局坐标形状为X的条（X长Y短），给出的是Y坐标，position从1开始
    return make_pair(-999, mapX[IndexBD(board)][IndexChip(chip)][channel]);
}

//______________________________________________________________________________
// 画三维击中
void MyBeamTestTrackVMM::DrawHits(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med, MyBeamTestHitData *eventList)
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

void MyBeamTestTrackVMM::DrawDetector(TGeoManager *geom, TGeoVolume *top, TGeoMedium *med)
{
    TGeoVolume *box = geom->MakeBox(name, med, siz[2] / 2., siz[1] / 2., siz[0] / 2.);
    TGeoTranslation *tr = new TGeoTranslation(pos[2], -pos[1], -pos[0]);
    TGeoRotation *rott = new TGeoRotation(name, rot[0] - 90, rot[1], 0.);
    TGeoCombiTrans *comb = new TGeoCombiTrans(*tr, *rott);
    box->SetLineColorAlpha(kAzure - 9, 50);
    top->AddNode(box, 1, comb);
}
//______________________________________________________________________________
// 生成分布图
bool MyBeamTestTrackVMM::FillDistribution(int i, MyBeamTestHitData *eventList)
{
    bool hitted = false;

    if (i == 0)
    {
        fAllHit = new TH2F(Form("AHit%d", id), Form("All hits for %s", name.Data()), 128, 0, 128, 128, 0, 128);
        fAllHit2 = new TH2F(Form("AHit_2%d", id), Form("All hits for %s", name.Data()), 150, -30 + pos[0], 30 + pos[0], 150, -30 + pos[1], 30 + pos[1]);
        fHit[0] = new TH1F(Form("HitX%d", id), "X Strip hits for single event", 128, 0, 128);
        fHit[1] = new TH1F(Form("HitY%d", id), "Y Strip hits for single event", 128, 0, 128);
        fCharge = new TH1F(Form("Charg%d", id), Form("Charge distribution of %s", name.Data()), 200, 0, 10000);

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

//______________________________________________________________________________
//生成DST-root时，将hit合并为cluster
void MyBeamTestTrackVMM::AnalysisCluster(MyBeamTestData *fEvent)
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

    cout << "\nafter mearge, Ybranch = " << fEvent->Ybranch.size() << endl;
    for (int ii = 0; ii < (int)fEvent->Ybranch.size(); ii++)
    {
        cout << ii << ":" << endl;
        for (int jj = 0; jj < (int)fEvent->Ybranch[ii].size(); jj++)
        {
            cout << "(" << fEvent->Ybranch[ii][jj].hit.first << ", " << fEvent->Ybranch[ii][jj].hit.second << ") ";
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
