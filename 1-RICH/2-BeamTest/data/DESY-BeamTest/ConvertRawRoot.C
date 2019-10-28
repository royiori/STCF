//--------------------------------
//
// 批量将二进制数据转换为raw-root格式
// raw-root格式只包含二进制解码的数据，不做任何处理
//
//---------------------------------
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
#include "../../inc/MyBeamTestDetector.h"
#include "../../inc/MyBeamTestRICH.h"
#include "../../inc/MyBeamTestTrackAGET.h"
#include "../../inc/MyBeamTestTrackVMM.h"
#include "../../inc/CombineData.h"
using namespace std;

TString GenPath(int type1, int type2, TString fFolder)
{
    TString fName = fFolder;

    if (type1 == RICH && type2 == RAW)
        fName += "/Combine/RICH-raw.root";
    if (type1 == RICH && type2 == PED)
        fName += "/Combine/RICH-ped.root";

    if (type1 == TrackerAGET && type2 == RAW)
        fName += "/Combine/TrackAGET-raw.root";
    if (type1 == TrackerAGET && type2 == PED)
        fName += "/Combine/TrackAGET-ped.root";

    if (type1 == TrackerVMM && type2 == RAW)
        fName += "/Combine/TrackVMM-raw.root";
    if (type1 == TrackerVMM && type2 == PED)
        fName += "/Combine/TrackVMM-ped.root";

    if (type1 == ALL && type2 == ALL)
        fName += "/Combine/Combined-dst.root";

    return fName;
}

void ConvertRawRoot()
{
    //RICH
    vector<int> boardList1 = {2, 10, 5, 7};
    vector<int> chipList1 = {10, 11, 12, 13};
    //Track-AGET
    vector<int> boardList2 = {6, 9, 1};
    vector<int> chipList2 = {10, 11, 12, 13};

    for (int i = 86; i < 87; i++)
    {
        FileStat_t fStat;
        TString fFolder(Form("./RUN%02d", i));
        gSystem->GetPathInfo(fFolder, fStat);
        if (fStat.fSize == 0)
            continue;

        cout << "\n==========================================\n";
        cout << "-->Opening " << fFolder << endl;

        TString fileDir1(fFolder + "/RICH/");
        TString fileDir2(fFolder + "/TrackAGET/");
        TString fileDir3(fFolder + "/TrackVMM/");

        vector<TString> datList1;
        vector<TString> datList2;
        vector<TString> datList3;

        GetFileList(fileDir1, ".dat", datList1);
        GetFileList(fileDir2, ".dat", datList2);
        GetFileList(fileDir3, ".bin", datList3);

        /*
        //快速检查事例数
        //...
        cout << "checking data..." << endl;
        TFile *fFile = new TFile(GenPath(RICH, RAW, fFolder));
        if (!fFile->IsOpen())
            continue;
        TTree *fTree = (TTree *)fFile->Get("tree");
        int Event;
        TBranch *b_Event;
        fTree->SetBranchAddress("event", &Event, &b_Event);
        int nentries = fTree->GetEntriesFast();

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
            if (fabs(tmp - Event) != 0)
                total++;
            tmp = (tmp == Event) ? tmp : Event;
        }
        cout << "--> Trig start from " << trgstart << " to " << trgstop << " , total recorded trig number: " << total << endl;

        continue;
        */

        if (datList1.size() > 0)
        {
            //ReadRICHData2Root(datList1, GenPath(RICH, RAW, fFolder), 0);
            //GenerateRICHPed(GenPath(RICH, RAW, fFolder), GenPath(RICH, PED, fFolder), boardList1, chipList1, 0);
        }

        if (datList2.size() > 0)
        {
            //ReadTrackAGTData2Root(datList2, GenPath(TrackerAGET, RAW, fFolder), 0);
            //GenerateAGETPed(GenPath(TrackerAGET, RAW, fFolder), GenPath(TrackerAGET, PED, fFolder), boardList2, chipList2, 0);
        }

        if (datList3.size() > 0)
            ReadTrackVMMData2Root(datList3, GenPath(TrackerVMM, RAW, fFolder), 1);
    }
}
