void ReadPIDRoot()
{
    //可以从hit.root里读扫描的动量/角度/光子数范围min和max，简单起见也可以直接这里输入：
    double pmin = 0;
    double pmax = 4.;
    double themin = 0.;
    double themax = 30.;

    //读pid文件
    TString fname("./simdata/STCF-5mm-pid.root");
    TFile f(fname);
    if (!f.IsOpen())
        return;

    TTree *T1 = (TTree *)f.Get("TTree1");
    TTree *T2 = (TTree *)f.Get("TTree2");
    TTree *T3 = (TTree *)f.Get("TTree3");

    //1. T1是分的bin数
    if (T1 == NULL)
        return;

    int nhyp, nrad, nmom, nthe, nph;
    T1->SetBranchAddress("nhyp", &nhyp);
    T1->SetBranchAddress("nrad", &nrad);
    T1->SetBranchAddress("nmom", &nmom);
    T1->SetBranchAddress("nthe", &nthe);
    T1->SetBranchAddress("nph", &nph);
    T1->GetEntry(0);

    cout << "---->Reading: nHypo=" << nhyp << ", nRadiator=" << nrad << ", nMomentum=" << nmom << ", nTheta0=" << nthe << ", nPhoton=" << nph << endl;
    cout << "---->Reading: nEntries=" << T2->GetEntries() << endl;

    double pstep = (pmax - pmin) / (nmom - 1);
    double thstep = (themax - themin) / (nthe - 1);

    //2. T2是重建的中心值和分辨率及对应的拟合误差
    if (T2 == NULL)
        return;

    double mean, merr;
    double sigm, serr;
    T2->SetBranchAddress("mean", &mean);
    T2->SetBranchAddress("merr", &merr);
    T2->SetBranchAddress("sigm", &sigm);
    T2->SetBranchAddress("serr", &serr);

    int ientry = 0;
    for (int ihypo = 0; ihypo < nhyp; ihypo++)
        for (int irad = 0; irad < nrad; irad++)
            for (int imom = 0; imom < nmom; imom++)
                for (int ithe = 0; ithe < nthe; ithe++)
                    for (int iph = 0; iph < nph; iph++)
                    {
                        T2->GetEntry(ientry++);

                        //根据需要拿对应的中心值和误差的结果：
                        double mom = pmin + imom * pstep;
                        double the = themin + ithe * thstep;

                        //只输出一部分的结果
                        if (iph != 10)
                            continue;
                        if (imom < 10 || imom > 20)
                            continue;
                        if (ithe != 0)
                            continue;
                        if (irad != 0)
                            continue;

                        cout << "Momentum = " << mom << " GeV/c, theta = " << the << " degree, "
                             << "particle = " << ihypo << ", radiator = " << irad << ", photon = " << iph << ", "
                             << "offset = " << mean << "+-" << merr << ", sigma=" << sigm << "+-" << serr << endl;
                    }

    if (T3 == NULL)
        return;

    double pideff;
    T3->SetBranchAddress("pideff", &pideff);

    ientry = 0;
    for (int imom = 0; imom < nmom; imom++)
        for (int ithe = 0; ithe < nthe; ithe++)
            for (int ihypo = 0; ihypo < nhyp; ihypo++)
                for (int jhypo = 0; jhypo < nhyp; jhypo++)
                {
                    T3->GetEntry(ientry++);

                    //根据需要拿对应的中心值和误差的结果：
                    double mom = pmin + imom * pstep;
                    double the = themin + ithe * thstep;

                    cout << "Momentum = " << mom << " GeV/c, theta = " << the << " degree, "
                         << "particle = " << ihypo << " identified as " << jhypo << " prob = " << pideff << endl;
                }

    cout << "----> Total entries loaded. " << endl;
}
