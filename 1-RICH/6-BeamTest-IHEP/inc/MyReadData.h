//-----------------------------------------------------
// to read data file
//
// usage as example:
//
/*
    MyReadData* myReadData = new MyReadData();
    myReadData->ReadData(fileName);

    // print basic info.
    myReadData->PrintInfo();
    myReadData->GetStatus(); // return status flag
    
    // get raw data
    vector<vector<double>> vdata;
    vdata = myReadData->GetData();

    // set seperator
    myReadData->SetSeperator('|'); // default is ','

    // sort data
    myReadData->SetIndex(2);  // use colnum 2 as index， 0 as default
    myReadData->SortData();   // sort data using colnum 2

    // get data
    myReadData->GetData(0, 2);   // return [0][2], 0 as default
    myReadData->GetValue(1.3);   // return y @ index=1.3 by extrapolate

    // draw plot
    TGraph *ftmp1;
    TH1F *ftmp2;
    ftmp1 = myReadData->GetGraph(0, 1);  //plot colnum 0 as x, colnum 1 as y
    ftmp2 = myReadData->GetHistogram(0); //plot colnum 0 as histogram
    
    // save to root
    myReadData->SetBranchName(0, "x");
    myReadData->SetBranchName(1, "y");
    myReadData->SetBranchName(2, "val");
    myReadData->SaveToRootFile(rootfileName);
*/
// data file format:
/*
#this is data format to show how it works
#start with LABEL will be read as label
#LABEL x   y   val
0.1, 1.1, 3.2
0.2, 2.2, 4.2
...
*/
//
// return status meaning
//   -1  -- something went wrong
//    0  -- init
//    1  -- read data correctly
//
//                            Q.LIU
//-----------------------------------------------------

#ifndef MYREADDATA_H
#define MYREADDATA_H

#include "TGraph.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "Riostream.h"

using namespace std;

class MyReadData
{
public:
    MyReadData(const char *fn, const char *nm)
    {
        init();
        name = nm;
        ReadData(fn);
        //if (ncol > 0 && nrec > 0)
        //    PrintInfo();
    };
    MyReadData() { init(); };
    ~MyReadData(){};

    // 1. read data.
    void ReadData(TString fn);
    void ReadData(const char *fn) { ReadData(TString(fn)); }
    void SetSeperator(const char c) { sep = c; }

    // 2. get info.
    void PrintInfo(int printlevel = 1);
    int GetStatus() { return status; }
    TString GetName() { return name; }
    TString GetFileName() { return fileName; }

    // 3. sort data.
    void SortData();
    void SetIndex(int _idx = 0);

    // 4. get data.
    double GetValue(double xval, int yrec=1);
    double GetData(int xrec, int yrec);

    // 5. draw plots.
    TGraph *GetGraph(int x, int y);
    TH1 *GetHistogram(int x);

    // 6. save root.
    void SetBranchName(int id, const char *bname) { SetBranchName(id, TString(bname)); }
    void SetBranchName(int id, TString bname);
    void SaveToRootFile(const char *rname) { SaveToRootFile(TString(rname)); }
    void SaveToRootFile(TString rname);

    // 7. get data
    vector<vector<double>> GetData() { return data; }
    void SetVerbose(int verb) { verbose = verb; }

private:
    void init()
    {
        idx = 0;
        ncol = 0;
        nrec = 0;
        status = 0;
        verbose = 0;
        sep = ',';
        g = NULL;
        h = NULL;
    }
    void saveRecord(TString line);
    void saveLabels(TString line);
    void trimLine(TString &line);

    vector<vector<double>> data;
    vector<TString> label;
    TH1 *h;
    TGraph *g;
    TString name;
    TString fileName;

    int idx;
    int status;
    int verbose;

    int ncol; //column
    int nrec;
    char sep;
};

#endif