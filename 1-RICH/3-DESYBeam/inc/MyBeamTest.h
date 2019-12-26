#ifndef MyBeamTest_h
#define MyBeamTest_h

#include "TGText.h"
#include "TGeoManager.h"
#include "MyBeamTestRICH.h"
#include "MyBeamTestTrackAGET.h"
#include "MyBeamTestTrackVMM.h"

class MyBeamTest
{
public:
    MyBeamTest(TEnv *e);
    ~MyBeamTest();

    void StoreEnv();
    void DrawConfig();

    void SetDSPath(const char *fname);
    TString GetDSPath() { return fDSPath; }

    int ReadSettingsText(TGText *text);
    TString GenerateSettingText();
    TString GenPath(int type1, int type2, const char *fileName = NULL);

    //bool ConvtTrackVMMRoot(const char *fileName);
    bool ConvtRICHRoot(const char *fileName, int SaveWaveFlag);
    bool ConvtTrackAGTRoot(const char *fileName, int SaveWaveFlag);
    bool ConvtTrackVMMRoot(const char *fileName);
    bool ReadRICHPed(TString);
    bool ReadTrackAGTPed(TString);

    void CombineDSTRoot(const char *fileName);
    void LoadDSTRoot(const char *fileName);
    void AnalysisDSTRoot(const char *fileName);
    void SetRayTracing(bool tr) { tracing = tr; }
    int GetDSTEntries() { return DSTNEntries; }

    void DrawDSTHit(int entry);
    void DrawCurtHit(int entry);
    void DrawPrevHit();
    void DrawNextHit();

    vector<MyBeamTestRICH *> GetRICHDet() { return vRICH; }
    vector<MyBeamTestTrackAGET *> GetTrackAGET() { return vTrkAGT; }
    vector<MyBeamTestTrackVMM *> GetTrackVMM() { return vTrkVMM; }

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

    TString fDSPath; //

    bool tracing;

    int NRICH, NTrkAGT, NTrkVMM;
    vector<MyBeamTestRICH *> vRICH;
    vector<MyBeamTestTrackAGET *> vTrkAGT;
    vector<MyBeamTestTrackVMM *> vTrkVMM;

    void UpdateDetID()
    {
        NRICH = vRICH.size();
        NTrkAGT = vTrkAGT.size();
        NTrkVMM = vTrkVMM.size();

        for (int i = 0; i < NRICH; i++)
            vRICH[i]->SetDetID(i);
        for (int i = 0; i < NTrkAGT; i++)
            vTrkAGT[i]->SetDetID(i + NRICH);
        for (int i = 0; i < NTrkVMM; i++)
            vTrkVMM[i]->SetDetID(i + NRICH + NTrkAGT);
    }
};

extern MyBeamTest *gMyBeamTest;

#endif