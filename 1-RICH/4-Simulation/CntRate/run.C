//------- const values ------------
const int NType = 8;
int NEvent = 200;
int EvtID[NType + 1] = {0, 7716, 15675, 16275, 18928, 21639, 25049, 25246, 25332}; //
double freqlist[NType] = {79800, 3.9e+07, 5.59e+08, 3.5e+08, 2e+07, 2e+07, 1.03e+09, 1.07e+08};

void readFreq();
void genBatch();
void anaRoot();
//
//--
void run()
{
    //readFreq();
    //genBatch();
    anaRoot();
}

//
//-- read frequency from backgrounddata.root
void readFreq()
{
    double freq;
    TFile froot("./backgrounddata.root");
    TTree *branch = (TTree *)froot.Get("datatree");
    branch->SetBranchAddress("BackgroundFrequency", &freq);

    for (int i = 0; i < 8; i++)
    {
        branch->GetEntry(EvtID[i]);
        cout << i << " " << setprecision(5) << freq << endl;
        freqlist[i] = freq;
    }

    cout << "double freqlist[8] = {";
    for (int i = 0; i < 7; i++)
        cout << freqlist[i] << ", ";
    cout << freqlist[7] << "};\n";


}
//
//--
void genBatch()
{
    gSystem->Exec("mkdir ./batch");

    ofstream runfile;
    runfile.open("./batch.sh");

    for (int i = 0; i < NType; i++)
    {
        gSystem->Exec(Form("mkdir ./batch/type%d", i));
        int entStart = EvtID[i];
        int entStop = EvtID[i + 1];
        int TotEvent = entStop - entStart;

        for (int j = 0; j < TotEvent / NEvent + 1; j++)
        {
            int nevent = ((TotEvent > NEvent * (j + 1)) ? NEvent : (TotEvent - NEvent * j));
            if (nevent == 0)
                continue;

            ofstream outfile;
            TString filename1(Form("./batch/type%d/bg%d.txt", i, j));
            TString filename2(Form("./batch/type%d/result-%d.root", i, j));
            outfile.open(filename1);
            outfile << "#\n";
            outfile << "/control/verbose 0\n";
            outfile << "/run/verbose 0\n";
            outfile << "/process/verbose 0\n";
            outfile << "/tracking/verbose 0\n";
            outfile << "#\n";
            outfile << "/globalField/setValue 0 0 1 tesla\n";
            outfile << "#\n";
            outfile << "/MyRun/SetG4BasedFileName " << filename2 << "\n";
            outfile << "/MyGun/SetGunBGFile backgrounddata.root\n";
            outfile << "/MyGun/SetEntry0 " << entStart + j * NEvent << "\n";
            outfile << "#\n";
            outfile << "/run/initialize\n";
            outfile << "#\n";
            outfile << "/run/beamOn " << nevent << "\n";
            outfile.close();
            //gSystem->Exec(Form("./gdmlDet %s", filename.Data()));
            runfile << "./gdmlDet " << filename1 << endl;
        }
    }
    runfile.close();
    gSystem->Exec("chmod 777 ./batch.sh");
}

TChain charged("Charged");
TChain optph("OptPh");
void anaRoot()
{
    double Ncnt = 0;
    for (int i = 0; i < NType; i++)
    {
        charged.Reset();
        optph.Reset();

        int entStart = EvtID[i];
        int entStop = EvtID[i + 1];
        int TotEvent = entStop - entStart;

        for (int j = 0; j < TotEvent / NEvent + 1; j++)
        {
            TString filename2(Form("./batch/type%d/result-%d.root", i, j));
            charged.Add(filename2);
            optph.Add(filename2);
        }

        double ncharge = charged.GetEntries();
        double nopt = optph.GetEntries();
        double ncharge2 = (ncharge) / 10000 * freqlist[i];
        double nopt2 = (nopt) / 10000 * freqlist[i];
        double ncnt = (ncharge + nopt) / 10000 * freqlist[i];
        Ncnt += ncnt;

        cout << "----------------------------------------------------" << endl;
        cout << i << " Ncharged: " << ncharge << "; Nopt: " << nopt << endl;
        cout << "  Tot: " << ncnt << endl;
    }
    cout << "----------------------------------------------------" << endl;
    cout << "Sum: " << Ncnt << endl;

    charged.Reset();
    optph.Reset();
    for (int i = 0; i < 1; i++)//NType; i++)
    {
        int entStart = EvtID[i];
        int entStop = EvtID[i + 1];
        int TotEvent = entStop - entStart;

        for (int j = 0; j < TotEvent / NEvent + 1; j++)
        {
            TString filename2(Form("./batch/type%d/result-%d.root", i, j));
            charged.Add(filename2);
            optph.Add(filename2);
        }
    }
}