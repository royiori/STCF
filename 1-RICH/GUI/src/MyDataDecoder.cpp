#include "TSystem.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "TGeoVolume.h"
#include "MyDataDecoder.h"

MyDataDecoder *gMyDataDecoder = (MyDataDecoder *)0;

void MyDataDecoder::AnalysisPedestalData(TString fRawName)
{
    cout<<"AnalysisPedestalData: "<<fRawName<<endl;
}

void MyDataDecoder::ReadPedestalData(TString fRawName)
{
    cout<<"ReadPedestalData: "<<fRawName<<endl;
}

void MyDataDecoder::AnalysisAGETData(TString fRawName)
{
    cout<<"AnalysisAGETData: "<<fRawName<<endl;
}

void MyDataDecoder::ReadAGETData(TString fRawName)
{
    cout<<"ReadAGETData: "<<fRawName<<endl;
}