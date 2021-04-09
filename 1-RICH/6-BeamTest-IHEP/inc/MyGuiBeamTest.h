#ifndef MyGuiBeamTest_h
#define MyGuiBeamTest_h

#include "TGText.h"
#include "TPaveText.h"
#include "TGeoManager.h"

#include "MyBeamTestDatStruct.h"
#include "MyBeamTestRICH.h"
#include "MyBeamTestTrackAGET.h"

class MyGuiBeamTest
{
public:
    MyGuiBeamTest(TEnv *e);
    ~MyGuiBeamTest();

    void StoreEnv();
    void DrawConfig();

    void SetDSPath(const char *fname);
    TString GetDSPath() { return fDSPath; }

    void ReadSettingsText();
    TString GenerateSettingText();

    TString GenPath(int type1, int type2, const char *fileName = NULL);
    void GetFileList(TString filePath, TString filePattern, vector<TString> &fList);

    //文件类型转换
    void ConvtBinaryToRawRoot();
    void ConvtRawToDstRoot();
    void CombineDSTRoot(const char *fileName);
    void LoadDSTRoot(const char *fileName);
    void AnalysisDSTRoot(const char *fileName);
    int GetDSTEntries() { return DSTNEntries; }

    void DrawDSTHit(int entry);
    void DrawCurtHit(int entry);
    void DrawPrevHit();
    void DrawNextHit();

    vector<MyBeamTestRICH *> GetRICHDet() { return vRICH; }
    vector<MyBeamTestTrackAGET *> GetTrackAGET() { return vTrkAGT; }

private:
    TEnv *env;
    TGeoManager *geom;
    TGeoVolume *top;
    TGeoVolume *geohits;
    TGeoMedium *media;

    TPaveText *ptext;

    TFile *fDSTFile = 0;
    TTree *fDSTTree = 0;
    MyBeamTestHitData *fDSTEvent = 0;
    int DSTNEntries = -1;
    int iEntry = 0;

    TString fDSPath; //data path
    TString fSTPath; //setting.txt path

    bool SaveWaveFlag = kFALSE; //dst.root里是否保存wave波形，默认不保存，减少文件大小

    int NRICH, NTrkAGT;
    vector<MyBeamTestRICH *> vRICH;
    vector<MyBeamTestTrackAGET *> vTrkAGT;

    void UpdateDetID()
    {
        NRICH = vRICH.size();
        NTrkAGT = vTrkAGT.size();

        for (int i = 0; i < NRICH; i++)
            vRICH[i]->SetDetID(i);
        for (int i = 0; i < NTrkAGT; i++)
            vTrkAGT[i]->SetDetID(i + NRICH);
    }
};

extern MyGuiBeamTest *gMyGuiBeamTest;

#endif