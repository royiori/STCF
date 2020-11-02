#include "TEnv.h"
#include "TView.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TPolyLine3D.h"
#include "TGLViewer.h"
#include "TGraphErrors.h"

#include "MyGuiBeamTest.h"
#include "MyGuiMainAction.h"
#include "MyBeamTestDetector.h"

MyGuiBeamTest *gMyGuiBeamTest = (MyGuiBeamTest *)0;

//______________________________________________________________________________
// 构造和析构函数
MyGuiBeamTest::MyGuiBeamTest(TEnv *ev)
{
    env = ev;
    geom = NULL;

    fDSPath = env->GetValue("DSPath", "NOTSet");
    fSTPath = env->GetValue("STPath", "NOTSet");

    NRICH = 0;
    NTrkAGT = 0;
    ReadSettingsText();
}

MyGuiBeamTest::~MyGuiBeamTest()
{
    // Destructor.
    // Store env
    StoreEnv();
}

void MyGuiBeamTest::StoreEnv()
{
    env->SetValue("DSPath", fDSPath);
    env->SetValue("STPath", fSTPath);

    env->SaveLevel(kEnvLocal);
}

//______________________________________________________________________________
// 
void MyGuiBeamTest::DrawConfig()
{
    if (geom != NULL)
        delete geom;

    //------------------------------
    // geom顶层结构为top
    geom = new TGeoManager("beamgeo", "Beamtest geometry");
    TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
    media = new TGeoMedium("Vacuum", 1, mat);
    top = geom->MakeBox("Top", media, 1000., 500., 500.);
    geom->SetTopVolume(top);
    geom->SetTopVisible(false);
    geom->SetVerboseLevel(0);

    //------------------------------
    // top -> det -> sub-detectors
    TGeoVolume *det = geom->MakeBox("Detector", media, 1000., 500., 500.);
    det->SetVisibility(kFALSE);
    top->AddNode(det, 1, 0);

    for (int i = 0; i < NRICH; i++)
        vRICH[i]->DrawDetector(geom, det, media);
    for (int i = 0; i < NTrkAGT; i++)
        vTrkAGT[i]->DrawDetector(geom, det, media);

    TGeoVolume *ground = geom->MakeBox("Ground", media, 1000., 1., 1000);
    det->AddNode(ground, 1, new TGeoTranslation(500, 500, 0));
    ground->SetFillColor(kGray);

    //make beam
    {
        TGeoVolume *arrow = geom->MakeBox("arrow", media, 1000., 500., 500.);
        arrow->SetVisibility(kFALSE);

        //Z axis - 实验室坐标系
        TGeoVolume *xbar = geom->MakeBox("xbar", media, 550., 1., 1.);
        xbar->SetLineColor(kRed);
        arrow->AddNode(xbar, 1, new TGeoTranslation(450., 0., 0.));
        TGeoVolume *xarr = geom->MakeCone("xarr", media, 20, 0., 0.1, 0., 10.);
        xarr->SetLineColor(kRed);
        arrow->AddNode(xarr, 1, new TGeoCombiTrans(1000., 0., 0., new TGeoRotation("arrowRotX", 90, -90, 0)));

        //Y axis - 实验室坐标系
        TGeoVolume *ybar = geom->MakeBox("ybar", media, 1., 100., 1.);
        ybar->SetLineColor(kGreen);
        arrow->AddNode(ybar, 1, new TGeoTranslation(-100., -90., 0.));
        TGeoVolume *yarr = geom->MakeCone("yarr", media, 20, 0., 0.1, 0., 10.);
        yarr->SetLineColor(kGreen);
        arrow->AddNode(yarr, 1, new TGeoCombiTrans(-100., -200., 0., new TGeoRotation("arrowRotY", 0, -90, -90)));

        //X axis - 实验室坐标系
        TGeoVolume *zbar = geom->MakeBox("zbar", media, 1., 1., 100.);
        zbar->SetLineColor(kBlue);
        arrow->AddNode(zbar, 1, new TGeoTranslation(-100., 0., -90.));
        TGeoVolume *zarr = geom->MakeCone("zarr", media, 20, 0., 0.1, 0., 10.);
        zarr->SetLineColor(kBlue);
        arrow->AddNode(zarr, 1, new TGeoCombiTrans(-100., 0., -240., new TGeoRotation("arrowRotZ", 0, 0, 0)));

        det->AddNode(arrow, 1, new TGeoTranslation(0., 0., 0.));
    }

    //------------------------------
    // top -> hits -> sub-det's hits
    geohits = geom->MakeBox("geohits", media, 1000., 500., 500.);
    geohits->SetVisibility(kFALSE);
    top->AddNode(geohits, 1, 0);

    top->Draw();

    ptext = new TPaveText(.1, .8, .9, .97);
    ptext->AddText("Beam-test event display");
    ptext->Draw();
}

//___1. 生成&检查 raw.root______________________

//___1.1 生成和读取configure文本
//
TString MyGuiBeamTest::GenerateSettingText()
{
    TString settings("");

    settings += "\n#---------------------------";
    settings += "\n#    Beam-test settings";
    settings += "\n#---------------------------";
    settings += "\n\n# root file names for beam-test: \n";
    settings += Form("DSPath: %s\n", fDSPath.Data());
    settings += Form("STPath: %s\n", fSTPath.Data());

    
    settings += "\n\n#---------------------------";
    settings += "\n# RICH parameters: (in mm & degrees) \n";
    settings += Form("NRICH: %d\n", NRICH);
    for (int i = 0; i < (int)vRICH.size(); i++)
        settings += vRICH[i]->GenerateIntro();

    settings += "\n\n#---------------------------";
    settings += "\n# Tracker with AGET parameters: (in mm & degrees) \n";
    settings += Form("NTrkAGT: %d\n", NTrkAGT);
    for (int i = 0; i < (int)vTrkAGT.size(); i++)
        settings += vTrkAGT[i]->GenerateIntro();

    return settings;
}

void MyGuiBeamTest::ReadSettingsText()
{
    ifstream fin(fSTPath.Data());
    if (!fin)
    {
        cout << "####> Error opening " << fSTPath << " to read!" << endl;
        return;
    }

    string s;
    while (getline(fin, s))
    {
        TString line(s);

        //cout << "Read from file: " << s << endl;
        if (line.BeginsWith("#") || line.Length() < 5)
            continue;

        vector<TString> cont = gMyGuiMainAction->ReadContent(line);
        cont.erase(cont.begin());

        if (line.BeginsWith("NRICH"))
        {
            int ndet = cont[0].Atoi();
            if (ndet > NRICH)
                for (int i = NRICH; i < ndet; i++)
                    vRICH.push_back(new MyBeamTestRICH(i));
            NRICH = ndet;
        }

        if (line.BeginsWith("RICH"))
        {
            for (int i = 0; i < NRICH; i++)
            {
                if(line.BeginsWith(Form("RICH%d", i)))
                    vRICH[i]->SetParameters(line, cont);
            }
        }

        if (line.BeginsWith("NTrkAGET"))
        {
            int ndet = cont[0].Atoi();
            if (ndet > NTrkAGT)
                for (int i = NTrkAGT; i < ndet; i++)
                    vTrkAGT.push_back(new MyBeamTestTrackAGET(i));
            NTrkAGT = ndet;
        }

        if (line.BeginsWith("TrkAGET"))
        {
            for (int i = 0; i < NTrkAGT; i++)
            {
                if(line.BeginsWith(Form("TrkAGET%d", i)))
                    vTrkAGT[i]->SetParameters(line, cont);
            }
        }
    }

    //reset id:
    for (int i = 0; i < NRICH; i++)
        vRICH[i]->SetDetID(i);
    for (int i = NRICH; i < NRICH+NTrkAGT; i++)
        vTrkAGT[i]->SetDetID(i);

    //read mapping:
    for (int i = 0; i < NRICH; i++)
        vRICH[i]->ReadMapFile();
    for (int i = 0; i < NTrkAGT; i++)
        vTrkAGT[i]->ReadMapFile();
}

//___1.2 设置路径 & 查找配置文件
// 
void MyGuiBeamTest::SetDSPath(const char *fileName)
{
    fDSPath = TString(fileName);
    if (!fDSPath.EndsWith("idle.dat"))
    {
        cout << "##### you can simply choose 'idle.dat', all files will be generated automatically." << endl;
        fDSPath = TString("NOSET");
    }

    // 遍历此目录及向上三级文件夹下的所有文件，查找setting.txt文件
    // 请在setting.txt里保存此data文件的设置
    TString savdir = gSystem->WorkingDirectory();

    int upfold = 0;
    int found = 0;
    while (upfold < 3)
    {
        void *dirp = gSystem->OpenDirectory(".");
        const char *name;
        while ((name = gSystem->GetDirEntry(dirp)) != 0)
        {
            if (TString(name) == "settings.txt")
            {
                found = 1;
                fSTPath = gSystem->WorkingDirectory() + TString("/") + name;
                break;
            }
        }

        if (found)
            break;
        if (!gSystem->ChangeDirectory("../"))
            break;
        upfold++;
    }

    gSystem->ChangeDirectory(savdir.Data());

    if (!found)
    {
        cout << " #### could not find settings.txt, which contains the electronics & mapping information and is needed." << endl;
        fSTPath = TString("NOSET");
    }
}

TString MyGuiBeamTest::GenPath(int type1, int type2, const char *fileName)
{
    TString fName = (fileName == NULL) ? fDSPath : TString(fileName);

    if (!fName.EndsWith("idle.dat"))
    {
        cout << "##### you can simply choose 'idle.dat', all files will be generated automatically." << endl;
        return fName;
    }

    if (type1 == RICH && type2 == BIN)
        fName.ReplaceAll("idle.dat", "/RICH/");
    if (type1 == RICH && type2 == RAW)
        fName.ReplaceAll("idle.dat", "/Combine/RICH-raw.root");
    if (type1 == RICH && type2 == PED)
        fName.ReplaceAll("idle.dat", "/Combine/RICH-ped.root");

    if (type1 == TrackerAGET && type2 == BIN)
        fName.ReplaceAll("idle.dat", "/TrackAGET/");
    if (type1 == TrackerAGET && type2 == RAW)
        fName.ReplaceAll("idle.dat", "/Combine/TrackAGET-raw.root");
    if (type1 == TrackerAGET && type2 == PED)
        fName.ReplaceAll("idle.dat", "/Combine/TrackAGET-ped.root");

    if (type2 == DST)
        fName.ReplaceAll("idle.dat", "/Combine/");

    if (type1 == ALL && type2 == ALL)
        fName.ReplaceAll("idle.dat", "/Combine/Combined-dst.root");

    return fName;
}

void MyGuiBeamTest::GetFileList(TString filePath, TString filePattern, vector<TString> &fList)
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

//___1.3 读取二进制并转换为raw.root & ped.root 文件
// 
void MyGuiBeamTest::ConvtBinaryToRawRoot()
{
    TString fileDir = fDSPath;
    fileDir.ReplaceAll("idle.dat", "");

    TString fileDir1 = GenPath(RICH, BIN);
    TString fileDir2 = GenPath(TrackerAGET, BIN);

    vector<TString> datList1;
    vector<TString> datList2;

    GetFileList(fileDir1, ".dat", datList1);
    GetFileList(fileDir2, ".dat", datList2);

    if (datList1.size() > 0)
    {
        MyBeamTestRICH *_rich = new MyBeamTestRICH(-1);
        _rich->ReadData2RawRoot(datList1, GenPath(RICH, RAW), GetRICHDet(), 1);
        _rich->AnalysisPedestal(GenPath(RICH, RAW), GenPath(RICH, PED), GetRICHDet());
        delete _rich;
    }

    if (datList2.size() > 0)
    {
        MyBeamTestTrackAGET *_trk1 = new MyBeamTestTrackAGET(-1);
        _trk1->ReadData2RawRoot(datList2, GenPath(TrackerAGET, RAW), GetTrackAGET(), 1);
        _trk1->AnalysisPedestal(GenPath(TrackerAGET, RAW), GenPath(TrackerAGET, PED), GetTrackAGET());
        delete _trk1;
    }

    cout<<"-->> All binary data has been converted to raw.root. <<--"<<endl;
}

//___2. 生成&检查 dst.root______________________
// 转换RICH的raw-root 为 dst-root文件，同时分析波形得到Q/T
void MyGuiBeamTest::ConvtRawToDstRoot()
{
    MyBeamTestRICH *_rich = new MyBeamTestRICH(-1);
    _rich->ReadRaw2DstRoot(GenPath(RICH, RAW), GenPath(RICH, PED), GenPath(RICH, DST), GetRICHDet(), SaveWaveFlag);

    MyBeamTestTrackAGET *_trk1 = new MyBeamTestTrackAGET(-1);
    _trk1->ReadRaw2DstRoot(GenPath(TrackerAGET, RAW), GenPath(TrackerAGET, PED), GenPath(TrackerAGET, DST), GetTrackAGET(), SaveWaveFlag);

    cout<<"-->> All raw.root data has been converted to dst.root. <<--"<<endl;
}

//___3. 生成&检查 combine.root______________________
//
void MyGuiBeamTest::CombineDSTRoot(const char *fileName)
{
    //init
    for (int i = 0; i < (int)vRICH.size(); i++)
        vRICH[i]->ReadDSTRoot(GenPath(RICH, RAW, fileName));
    for (int i = 0; i < (int)vTrkAGT.size(); i++)
        vTrkAGT[i]->ReadDSTRoot(GenPath(TrackerAGET, RAW, fileName));

    //
    //设定扫描trig的区间
    int nMax = 0;
    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        int nEntries = vRICH[i]->GetEntries();
        int tmin = vRICH[i]->GetFirstTrig();
        int tmax = vRICH[i]->GetLastTrig();
        if (nEntries == 0)
            continue;
        nMax = (nMax < nEntries) ? nEntries : nMax;
        cout << "--> " << vRICH[i]->GetName() << " Events: " << nEntries << ", tigger range: " << tmin << " " << tmax << endl;
    }

    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        int nEntries = vTrkAGT[i]->GetEntries();
        int tmin = vTrkAGT[i]->GetFirstTrig();
        int tmax = vTrkAGT[i]->GetLastTrig();
        if (nEntries == 0)
            continue;
        nMax = (nMax < nEntries) ? nEntries : nMax;
        cout << "--> " << vTrkAGT[i]->GetName() << " Events: " << nEntries << ", tigger range: " << tmin << " " << tmax << endl;
    }

    cout << "--> Searching trigger from " << 0 << " to " << nMax << endl;

    //
    //定义全部组合后的数据结构
    TString fName = GenPath(ALL, ALL, fileName);
    TFile *fFile = new TFile(fName, "recreate");
    TTree *fTree = new TTree("tree", "beam test data structure");
    MyBeamTestHitData *fEventList = new MyBeamTestHitData;
    fTree->Branch("event", "MyBeamTestHitData", &fEventList, 8000, 2);

    vector<pair<double, double>> hit;
    vector<double> q;
    vector<double> t;

    //
    //扫描trig来合并事例
    for (int trig = 0; trig < nMax; trig++)
    {
        if (trig % 100 == 0)
            cout << "--> Searching trigger = " << trig << endl;

        fEventList->Init(trig);
        for (int i = 0; i < (int)vRICH.size(); i++)
            vRICH[i]->SearchTrigID(trig, fEventList);

        for (int i = 0; i < (int)vTrkAGT.size(); i++)
            vTrkAGT[i]->SearchTrigID(trig, fEventList);

        fTree->Fill();
    }
    cout << "--> " << fTree->GetEntries() << " events are stored." << endl;
    fTree->Write();
    fFile->Close();

    for (int i = 0; i < (int)vRICH.size(); i++)
        vRICH[i]->CloseDSTRoot();
    for (int i = 0; i < (int)vTrkAGT.size(); i++)
        vTrkAGT[i]->CloseDSTRoot();
}

void MyGuiBeamTest::LoadDSTRoot(const char *fileName)
{
    if (fDSTFile != 0)
    {
        fDSTFile->Close();
        fDSTFile = 0;
        fDSTTree = 0;
        fDSTEvent = 0;
    }

    //读入数据
    TString fName = GenPath(ALL, ALL, fileName);
    fDSTFile = new TFile(fName);
    if (!fDSTFile->IsOpen())
    {
        cout << "#### Error: Can't open " << fName << " to read." << endl;
        return;
    }
    cout << "--> Reading " << fName << "." << endl;

    fDSTTree = (TTree *)fDSTFile->Get("tree");
    fDSTTree->SetBranchAddress("event", &fDSTEvent);
    DSTNEntries = fDSTTree->GetEntries();
    cout << "--> Total hit entries = " << DSTNEntries << endl;

    //生成分布数据，计算效率，拟合并求残差
    int NTOT = NRICH + NTrkAGT;
    vector<int> hitList;  //此探测器是否击中的list
    vector<int> NHitList; //每个探测器有多少击中的list
    hitList.resize(NTOT);
    NHitList.resize(NTOT);
    int NHits = 0;

    for (int ip = 0; ip < DSTNEntries; ip++)
    {
        fDSTTree->GetEntry(ip);

        //生成分布
        for (int i = 0; i < NRICH; i++)
            hitList[i] = vRICH[i]->FillDistribution(ip, fDSTEvent);
        for (int i = 0; i < NTrkAGT; i++)
            hitList[i + NRICH] = vTrkAGT[i]->FillDistribution(ip, fDSTEvent);

        //计算效率
        for (int i = 0; i < NTOT; i++)
        {
            int temp = 1;
            for (int ii = 0; ii < NTOT; ii++)
                temp *= ((ii == i) ? 1 : hitList[ii]);
            if (temp == 1)
                NHitList[i]++;
        }

        int temp = 1;
        for (int i = 0; i < (int)hitList.size(); i++)
            temp *= hitList[i];
        if (temp == 1)
            NHits++;

        //拟合径迹
        /*vector<vector<vector<double>>> hitpos; //Ndetector, Ntrack, Ndim
        hitpos.resize(NTOT);
        for(int i=0; i<NTOT; i++)
            hitpos.resize(4); 
        for (int i = 0; i < NRICH; i++)
            vRICH[i]->FillHitVector(hitpos[i], fDSTEvent);
        for (int i = 0; i < NTrkAGT; i++)
            vTrkAGT[i]->FillHitVector(hitpos[i+NRICH], fDSTEvent);
        */
    }

    cout << "Total Events=" << DSTNEntries << ", Effect Events=" << NHits << ", Ratio=" << NHits * 1.0 / DSTNEntries << endl;
    for (int i = 0; i < NRICH; i++)
        cout << ": " << vRICH[i]->GetName() << " Efficiency = " << ((NHitList[i] > 0) ? NHits * 1.0 / NHitList[i] : 0) << endl;
    for (int i = 0; i < NTrkAGT; i++)
        cout << ": " << vTrkAGT[i]->GetName() << " Efficiency = " << ((NHitList[i + NRICH] > 0) ? NHits * 1.0 / NHitList[i + NRICH] : 0) << endl;

    DrawCurtHit(0);
}

//___4. 事例显示 ___________________________________
//
void MyGuiBeamTest::DrawDSTHit(int entry)
{
    if (entry < 0 || entry > DSTNEntries - 1)
        return;

    ptext->Clear();
    ptext->AddText(Form("Event: %d", entry));

    //remove existing hits
    for (int i = 0; i < geohits->GetNdaughters(); i++)
    {
        geohits->RemoveNode(geohits->GetNode(0));
    }

    fDSTTree->GetEntry(entry);
    for (int i = 0; i < NRICH; i++)
        vRICH[i]->DrawHits(geom, geohits, media, fDSTEvent);
    for (int i = 0; i < NTrkAGT; i++)
        vTrkAGT[i]->DrawHits(geom, geohits, media, fDSTEvent);

    top->Draw();
    ptext->Draw("same");
}

void MyGuiBeamTest::DrawPrevHit()
{
    iEntry = (iEntry - 1 < 0) ? 0 : iEntry - 1;
    DrawDSTHit(iEntry);
}
void MyGuiBeamTest::DrawNextHit()
{
    iEntry = (iEntry + 1 >= DSTNEntries) ? iEntry : iEntry + 1;
    DrawDSTHit(iEntry);
}
void MyGuiBeamTest::DrawCurtHit(int entry)
{
    if (entry >= 0 && entry < DSTNEntries)
    {
        iEntry = entry;
        DrawDSTHit(iEntry);
    }
}

// for testing function only
void MyGuiBeamTest::AnalysisDSTRoot(const char *fileName)
{
    if (fDSTFile != 0)
    {
        fDSTFile->Close();
        fDSTFile = 0;
        fDSTTree = 0;
        fDSTEvent = 0;
    }

    //读入数据
    TString fName = GenPath(ALL, ALL, fileName);
    fDSTFile = new TFile(fName);
    if (!fDSTFile->IsOpen())
    {
        cout << "#### Error: Can't open " << fName << " to read." << endl;
        return;
    }
    cout << "--> Reading " << fName << "." << endl;

    fDSTTree = (TTree *)fDSTFile->Get("tree");
    fDSTTree->SetBranchAddress("event", &fDSTEvent);
    DSTNEntries = fDSTTree->GetEntries();
    cout << "--> Total hit entries = " << DSTNEntries << endl;

    for (int ip = 0; ip < DSTNEntries; ip++)
    {
        if (ip % 1000 == 0 || ip < 10)
            cout << "--> event: " << ip << endl;

        fDSTTree->GetEntry(ip);
    }

    return;
}