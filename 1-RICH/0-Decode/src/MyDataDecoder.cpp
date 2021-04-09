#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TGeoVolume.h"
#include "MyDataDecoder.h"

#include "MyDatabaseClass.h"

MyDataDecoder *gMyDataDecoder = (MyDataDecoder *)0;

void MyDataDecoder::AnalysisPedestalData(TString fRawName)
{
    cout << "AnalysisPedestalData: " << fRawName << endl;
}

void MyDataDecoder::ReadPedestalData(TString fRawName)
{
    cout << "ReadPedestalData: " << fRawName << endl;
}

void MyDataDecoder::AnalysisAGETData(TString fRawName)
{
    fstream InDat;
    InDat.open(fRawName, ios::in | ios::binary);

    cout << "--> AnalysisData : " << fRawName << endl;
    int length = sizeof(unsigned short);
    unsigned short memblock;

    TString fRootName = fRawName;
    fRootName.ReplaceAll(".dat", ".root");

    TFile *fFile = new TFile(fRootName, "RECREATE");
    MyDatabaseClass *fDat = new MyDatabaseClass();
    TTree *fTree = new TTree("tree", "tree");
    fTree->Branch("event", "MyDatabaseClass", &fDat, 64000, 0);

    fDat->evtID = 1;
    fTree->Fill();
    
    fDat->evtID = 2;
    fTree->Fill();

    fDat->evtID = 3;
    fTree->Fill();
    /*
    while (!InDat.eof())
    {
        //1. read head
        InDat.read((char *)(&memblock), length);

        if (memblock != 0xAC0F)
        {
            cout << "Error!!!:....." << endl;
            return;
        }

        //2. read protocol header
        InDat.read((char *)(&memblock), length);

        int soe = (memblock & 0x4000) >> 14;
        int eoe = (memblock & 0x2000) >> 13;
        int packsize = (memblock & 0x1fff);

        if (soe == 1 && eoe == 0)
        {
            //start of event
            InDat.read((char *)(&memblock), length);

            // time stamp
            int pfx_start = (memblock & 0xFF00) >> 8;
            int sid = (memblock & 0x1F);
            int xy = (sid == 16) ? 0 : 1;

            InDat.read((char *)(&memblock), length);
            long int lsb = memblock;
            InDat.read((char *)(&memblock), length);
            long int msb = memblock;
            InDat.read((char *)(&memblock), length);
            long int hsb = memblock;

            fDat->timeStamp[xy] = hsb << 32 | msb << 16 | lsb;
            // event cnt
            //...
        }
        else if(soe == 0 && eoe == 0)
        {
            // package data

        }
        else if(soe == 0 && eoe == 1)
        {
            // eoe
            //...


            fDat->eoeCnt ++;
            if(fDat->eoeCnt == 2)
            {
                fTree->Fill();
                fDat->Clear();
            }
        }
        else 
        {
            cout << "ERROR"<<endl;
            return;
        }
        //...........

    
    }
    */


    fFile->Flush();
    fFile->Close();

    cout << "--> Raw data file has been convert to root-file: " << fRootName << ".\n"
         << endl;
}

void MyDataDecoder::ReadAGETData(TString fRawName)
{
    cout << "ReadAGETData: " << fRawName << endl;
    //fTree->SetBranchAddress("wave", &gwave);
}