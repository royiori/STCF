#ifndef MyDataDecoder_h
#define MyDataDecoder_h

#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"

using namespace std;

class MyDataDecoder
{
public:
    MyDataDecoder() {;}
    ~MyDataDecoder() {;}

    //读取 dat 数据并生成 RAW.root 文件
    void AnalysisPedestalData(TString fRawName);
    void ReadPedestalData(TString fRawName);
   
    void AnalysisAGETData(TString fRawName);
    void ReadAGETData(TString fRawName);

private:
    
};

extern MyDataDecoder *gMyDataDecoder;
#endif