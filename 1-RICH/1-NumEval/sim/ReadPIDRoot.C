#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

//设置文件名
TString HeadName("./STCF-5mm-batch");

//读pid文件
void ReadPIDRoot()
{
    const int NHYPO = 5;
    const char *SHYPO[NHYPO] = {"p", "k", "#pi", "#mu", "e"};

    //1. 读一个pid文件, 取出关键量
    TFile f(HeadName + TString("_pid.root"));
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

    TTree *T2 = (TTree *)f.Get("PidTree");
    if (T2 == NULL)
        return;

    double pideff;
    T2->SetBranchAddress("pideff", &pideff);

    //画图
    TFile I(HeadName + TString("_out.root"), "RECREATE");
    if (!I.IsOpen())
        return;

    int dMom = 81;
    TH2D * hMom[dMom][NHYPO];
    for (int ihypo = 0; ihypo < nhyp; ihypo++) 
    {
      for(int jhypo = 0; jhypo < nhyp; jhypo++)
      {
        const string name = Form("PID efficiency for %s as %s", SHYPO[ihypo], SHYPO[jhypo]);
        hMom[ihypo][jhypo] = new TH2D(Form("hmom_%d_%d", ihypo, jhypo), name.c_str(), 81, 0, 4, 71, 0, 70 );
      }
        
    }

    int ientry = 0;
    for (int imom = 0; imom < nmom; imom++)
        for (int ithe = 0; ithe < nthe; ithe++)
            //for (int ithe = 0; ithe < nhyp; ithe++)
            for (int ihypo = 0; ihypo < nhyp; ihypo++)
                for (int jhypo = 0; jhypo < nhyp; jhypo++)
                {
                    T2->GetEntry(ientry++);
                    fPidEffList[imom][ithe][ihypo][jhypo] = pideff;

                    //根据需要拿对应的中心值和误差的结果：
                    double mom = pMin + imom * pstep;
                    double the = The0Min + ithe * thstep;
                    //cout << "Momentum = " << mom << " GeV/c, theta = " << the << " degree, "
                         //<< "particle = " << SHYPO[ihypo] << " identified as " << SHYPO[jhypo] << " prob = " << pideff << endl;
                   
                    hMom[ihypo][jhypo] -> Fill(mom,the,pideff);
                }

    cout << "----> Total entries loaded. " << endl;
   
    for (int ihypo = 0; ihypo < nhyp; ihypo++) 
    {
      for(int jhypo = 0; jhypo < nhyp; jhypo++)
      {
        hMom[ihypo][jhypo] -> GetXaxis() -> SetTitle("/#it{GeV}");
        hMom[ihypo][jhypo] -> GetYaxis() -> SetTitle("#it{#theta}");
        hMom[ihypo][jhypo]-> Write();
      }
        
    }
    // 出图检查

    I.Close();
    f.Close();
}
