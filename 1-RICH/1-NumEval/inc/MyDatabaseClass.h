#ifndef _MyDatabaseClass_h_
#define _MyDatabaseClass_h_

#include <map>
#include <vector>
#include "TSystem.h"
#include "TString.h"
#include "TGraph.h"
#include "MyReadData.h"

using namespace std;

//-----------------------------------------------------------------------------------------------
//
class MyDatabaseClass
{
public:
    MyDatabaseClass();
    ~MyDatabaseClass();

    int GetNMaterial() { return nMaterial; }
    vector<TString> GetMaterialList() { return MaterialList; }
    TString GetMaterialName(int id) { return MaterialList[id]; }

    int GetNDetector() { return nDetector; }
    vector<TString> GetDetectorList() { return DetectorList; }
    TString GetDetectorName(int id) { return DetectorList[id]; }

    double GetMatAbsValue(TString matStr, double lambda, double ppm);
    double GetMatRefValue(TString matStr, double lambda);
    double GetDetQEValue(TString detStr, double lambda);

    TGraph *GetMatAbsGraph(TString matStr, int nL, double lmin, double lmax, double ppm);
    TGraph *GetMatTrsGraph(TString matStr, int nL, double lmin, double lmax, double ppm);
    TGraph *GetMatRefGraph(TString matStr, int nL, double lmin, double lmax);
    TGraph *GetDetQEGraph(TString detStr, int nL, double lmin, double lmax);

private:
    int nMaterial;
    vector<TString> MaterialList;

    int nDetector;
    vector<TString> DetectorList;

    vector<MyReadData *> MatAbsData;
    vector<MyReadData *> MatRefData;
    vector<MyReadData *> DetQEData;

    vector<TGraph *> MatAbsGraph;
    vector<TGraph *> MatTrsGraph;
    vector<TGraph *> MatRefGraph;
    vector<TGraph *> DetQEGraph;

    int GetID(TString tstr, vector<TString> tlist)
    {
        int id = -1;
        for (int i = 0; i < (int)tlist.size(); i++)
            if (tlist.at(i) == tstr)
                id = i;
        return id;
    }

};

#endif
