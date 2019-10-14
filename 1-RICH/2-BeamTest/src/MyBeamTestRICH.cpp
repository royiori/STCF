#include "TH2F.h"
#include "MyBeamTestRICH.h"

//______________________________________________________________________________
//
MyBeamTestRICH::MyBeamTestRICH(int _id, int nx, int ny, double wx, double wy, double la) : MyBeamTestDetector(_id, nx, ny, wx, wy, la)
{
    fAllHit = 0;
    fCharge = 0;
}

MyBeamTestRICH::~MyBeamTestRICH()
{
}

//______________________________________________________________________________
// 读取map文件
void MyBeamTestRICH::ReadMapFile()
{
    fstream mapFp;
    mapFp.open("./data/Map/RICH_Map.txt", ios::in);
    if (!mapFp.is_open())
        cout << "RICH Detector Map Open Failed! " << endl;

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

//______________________________________________________________________________
// Combine数据
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
        fAllHit = new TH2F(Form("AHit%d", id), Form("All hits for %s", name.Data()), 32, 0, 32, 32, 0, 32);
        fAllHit2 = new TH2F(Form("AHit%d_2", id), Form("All hits for %s", name.Data()), 36, -90 + pos[0], 90 + pos[0], 36, -90 + pos[1], 90 + pos[1]);
        fCharge = new TH1F(Form("Charg%d", id), Form("Charge distribution of %s", name.Data()), 200, 0, 10000);
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

//______________________________________________________________________________
//生成DST-root时，将hit合并为cluster
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
