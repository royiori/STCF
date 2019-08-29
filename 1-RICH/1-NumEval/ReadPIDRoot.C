#include "TString.h"
#include "TFile.h"
#include "TTree.h"

//设置文件名
TString HeadName("./sim/STCF-5mm-batch");

//读pid文件
void ReadPIDRoot()
{
    const int NHYPO = 5;
    const char *SHYPO[NHYPO] = {"p", "k", "#pi", "#mu", "e"};

    //1. 读一个pid文件, 取出关键量
    TFile f(HeadName + TString("_pid_0.root"));
    if (!f.IsOpen())
        return;

    TTree *T1 = (TTree *)f.Get("Configure");
    if (T1 == NULL)
        return;

    int nhyp, nrad, nmom, nthe, nph;
    double pMin, pMax, The0Min, The0Max;

    T1->SetBranchAddress("nhyp", &nhyp);
    T1->SetBranchAddress("nrad", &nrad);
    T1->SetBranchAddress("nmom", &nmom);
    T1->SetBranchAddress("nthe", &nthe);
    T1->SetBranchAddress("nph", &nph);
    T1->SetBranchAddress("pMin", &pMin);
    T1->SetBranchAddress("pMax", &pMax);
    T1->SetBranchAddress("The0Min", &The0Min);
    T1->SetBranchAddress("The0Max", &The0Max);
    T1->GetEntry(0);

    cout << "---->Reading: nHypo=" << nhyp << ", nRadiator=" << nrad << ", nMomentum=" << nmom << ", nTheta0=" << nthe << ", nPhoton=" << nph << endl;

    double pstep = (nmom <= 1) ? 0 : (pMax - pMin) / (nmom - 1);
    double thstep = (nthe <= 1) ? 0 : (The0Max - The0Min) / (nthe - 1);

    f.Close();

    //------------------------------------
    //2. 读取PID的效率和误判率

    vector<vector<vector<vector<double>>>> fPidEffList; //[imom][ithe][ihypo][pidhypo]

    fPidEffList.clear();
    fPidEffList.resize(nmom);
    for (int imom = 0; imom < nmom; imom++)
    {
        fPidEffList[imom].resize(nthe);
        for (int ithe = 0; ithe < nthe; ithe++)
        {
            fPidEffList[imom][ithe].resize(nhyp);
            for (int ihypo = 0; ihypo < nhyp; ihypo++)
                fPidEffList[imom][ithe][ihypo].resize(nhyp);
        }
    }

    for (int imom = 0; imom < nmom; imom++)
    {
        TString pidFile = HeadName;
        pidFile += TString(Form("_pid_%d.root", imom));

        TFile f(pidFile);
        if (!f.IsOpen())
        {
            cout << "##### Can't open " << pidFile << "." << endl;
            return;
        }

        TTree *T2 = (TTree *)f.Get("PidTree");
        if (T2 == NULL)
            return;

        double pideff;
        T2->SetBranchAddress("pideff", &pideff);

        int ientry = 0;
        for (int ithe = 0; ithe < nthe; ithe++)
            for (int ihypo = 0; ihypo < nhyp; ihypo++)
                for (int jhypo = 0; jhypo < nhyp; jhypo++)
                {
                    T2->GetEntry(ientry++);
                    fPidEffList[imom][ithe][ihypo][jhypo] = pideff;

                    //根据需要拿对应的中心值和误差的结果：
                    double mom = pMin + imom * pstep;
                    double the = The0Min + ithe * thstep;
                    cout << "Momentum = " << mom << " GeV/c, theta = " << the << " degree, "
                         << "particle = " << SHYPO[ihypo] << " identified as " << SHYPO[jhypo] << " prob = " << pideff << endl;
                }
    }

    cout << "----> Total entries loaded. " << endl;
}
