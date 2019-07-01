#include "MyCommonRICH.h"
#include "MyGuiActionClass.h"

vector<double> reflist; //reflective index of each radiator list
vector<double> thklist; //thickness of each radiator list
vector<double> lenlist; //light path length in each radiator list
vector<double> abslist; //absorption length in each radiator list

const char *SHYPO[NHYPO] = {"p", "k", "#pi", "#mu"};

//______________________________________________________________________________
// 1. 全局坐标系定义： 最上层辐射体的上表面为Z = 0 处，往上为正，往下为负，即粒子沿-Z方向运动.
//                  入射带电粒子进入Z = 0处的坐标为原点 X = 0, Y = 0, 且粒子始终在X = 0的平面.
//                  粒子入射与Z轴正向夹角为theta0， 取值范围为[-pi/2, pi/2]， 逆时针为正，与+Y夹角为pi/2
//
// 2. 粒子在当前辐射体（厚度t0）内的发光点的坐标为 z0(循环里的变量), 则在全局坐标系里有：Z0 = -1*t(上面所有层辐射体厚度) - z0
//                  此发光点在X/Y平面的坐标为： X0 = 0, y0 = Z0 * tan(theta0)
//
// 3. findRho是在极坐标系下求解，此时极坐标系下的中心点为发光点位置，因此在计算时用的有效发光长度 l0 = t0 - z0
//                  极坐标与+X夹角为phi，求解的结果转换到全局坐标系有：
//                      Xr = rho * cos(phi); Yr = rho * sin(phi) + Y0
//______________________________________________________________________________

//______________________________________________________________________________
// math func
double findA(double theta0, double thetac, double phi)
{
    double cc = cos(thetac);
    double s0 = sin(theta0);
    double c0 = cos(theta0);
    double t0 = tan(theta0);
    double sp = sin(phi);
    double c2c = cos(2 * thetac);
    double c20 = cos(2 * theta0);
    double c2d = cos(2 * phi);

    if (c0 == 0)
        return 0;

    double a = 1;
    double b = 1 + c20 - 2 * c2c - 2 * c2d * s0 * s0;
    double c = 2 * cc * cc / c0 / c0 - 2 * sp * sp * t0 * t0;

    if (b < 0 || c <= 0)
        return 0;

    a *= cc / c0 / c0 * sqrt(b) + 2 * sp * t0;
    a /= c;
    return a;
}

double findRho(Double_t *x, Double_t *par)
{
    double phi = x[0];

    double tt = thklist[0]; //t0 : 当前辐射体厚度
    double z0 = par[0];     // z0 : 发光点到当前辐射体上表面的距离
    double l0 = tt - z0;
    double theta0 = par[1];
    double thetac = par[2];
    //double verbose = par[3];

    //cout<<"tt="<<tt<<" z0="<<z0<<" l0="<<l0<<endl;
    //cout<<"theta0="<<theta0<<" thetaC="<<thetac<<endl;

    double a = findA(theta0, thetac, phi);
    if (a == 0)
        return 0;

    double nSinTheta = reflist[0] * (a / sqrt(a * a + 1));
    a *= l0;
    lenlist.clear();
    lenlist.push_back(sqrt(a * a + l0 * l0));
    //cout<<"a0="<<a<<" n0="<<reflist[0]<<" n*sin="<<nSinTheta<<" "<<a<<endl;

    for (int i = 1; i < (int)thklist.size(); i++)
    {
        double lt = thklist[i];
        double nt = reflist[i];
        if (nSinTheta / nt > 1)
            return 0; //full reflected
        double thetat = asin(nSinTheta / nt);
        double at = lt * tan(thetat);
        a += at;
        lenlist.push_back(sqrt(at * at + lt * lt));
        //cout<<"a"<<i<<"="<<a<<" n"<<i<<"="<<nt<<" "<<lt<<" "<<nSinTheta<<" "<<nt<<endl;
    }

    return a;
}

//______________________________________________________________________________
// 重建发光点，在知道切伦科夫角+波长时可以严格求解出发光点
// 多层都时候没考虑！！！！y0计算都不对！！！
double findZ0Solver(TF1 *fcn1, double x1, double y1, double z) //计算半径和距离的差, z是到顶层的距离
{
    double the0 = fcn1->GetParameter(1);
    double x0 = 0;
    double y0 = (thklist[0] - z) * tan(the0); //注意z的定义
    double rad = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
    double phi = atan((y1 - y0) / (x1 - x0));
    fcn1->SetParameter(0, z);
    double rho = fcn1->Eval(phi);
    return rho - rad;
}

double findZ0(TF1 *fcn1, double x1, double y1, double zmin, double zmax, double err) //根据阳极上X/Y，在Zmin～Zmax区间求解发光点Z
{
    if (findZ0Solver(fcn1, x1, y1, zmin) * findZ0Solver(fcn1, x1, y1, zmax) > 0)
        return -1;

    int WDog = 0;
    double zmid;

    while (WDog < 1000)
    {
        zmid = (zmin + zmax) / 2.;

        double val1 = findZ0Solver(fcn1, x1, y1, zmin);
        double val2 = findZ0Solver(fcn1, x1, y1, zmid);

        if (fabs(val2) < err)
            return zmid;
        if (val1 * val2 > 0)
            zmin = zmid;
        else
            zmax = zmid;

        WDog++;
    }
    return -1;
}

//______________________________________________________________________________
// 重建theta_c，不知道发光点和波长，假设发光点在辐射体的中心处，且波长为185nm来求解
double findThetaC(TF1 *fcn1, double x1, double y1, double x0, double y0, double thetaExp, double err) //根据阳极上X/Y，求解theta_c
{
    int WDog = 0;
    double thetamid;
    double thetaStep = 0.001;
    double thetamin = thetaStep;
    double thetamax = 1.0;

    double rho;
    double rad = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
    double phi = atan((y1 - y0) / (x1 - x0));

    while (WDog < 1000)
    {
        thetamid = (thetamin + thetamax) / 2.;

        fcn1->SetParameter(2, thetamin);
        rho = fcn1->Eval(phi);
        double val1 = (rho == 0) ? 0 : rho - rad;

        fcn1->SetParameter(2, thetamid);
        rho = fcn1->Eval(phi);
        double val2 = (rho == 0) ? 0 : rho - rad;

        if (val1 == 0) //无解
        {
            thetamin += thetaStep;
            if (thetamin >= thetamax)
                break;
            continue;
        }
        if (val2 == 0) //无解
        {
            thetamax -= thetaStep;
            if (thetamax <= thetamin)
                break;
            continue;
        }

        if (fabs(val2) < err)
            return thetamin;

        if (val1 * val2 > 0)
            thetamin = thetamid;
        else
            thetamax = thetamid;

        WDog++;
    }

    return -999;
}
//______________________________________________________________________________
//
//   MyCommonRICH
//
//______________________________________________________________________________
//
MyCommonRICH::MyCommonRICH()
{
    fRhoFcn = new TF1("fRhoFcn", findRho, 0, 2 * TMath::Pi(), 4);
    fRhoFcn->SetParameters(0, 0, 0, 0);
    gDet = new MyRICHDetector(0);
    gDet->SetNHypothesis(NHYPO);
    epsilon = 1e-5;
    nEvent = 1e4;
}

MyCommonRICH::~MyCommonRICH()
{
    // Destructor.
}

//______________________________________________________________________________
/// Mass & Beta function
int MyCommonRICH::GetHypoID(TString p)
{
    p.ToLower();
    if (p == "mu" || p == "muon" || p == "#mu")
        return 3;
    if (p == "pi" || p == "pion" || p == "#pi")
        return 2;
    if (p == "k" || p == "kaon")
        return 1;
    if (p == "p" || p == "proton")
        return 0;
    return 3;
}

double MyCommonRICH::GetMass(TString p)
{
    p.ToLower();
    if (p == "mu" || p == "muon" || p == "#mu")
        return gkMassMuon;
    if (p == "pi" || p == "pion" || p == "#pi")
        return gkMassPion;
    if (p == "k" || p == "kaon")
        return gkMassKaon;
    if (p == "p" || p == "proton")
        return gkMassProton;
    return gkMassPion;
}

double MyCommonRICH::Beta(double p, double m) //p in GeV, m in GeV
{
    if (p == 0 && m == 0)
        return 0;

    return sqrt(p * p / (p * p + m * m));
}

double MyCommonRICH::BetaByTheta(double theta, double n)
{
    if (n == 0)
        return 0;
    if (cos(theta) == 0)
        return 0;
    return 1 / n / cos(theta);
}

double MyCommonRICH::Momentum(double beta, double m)
{
    if (fabs(beta) == 1.)
        return 0.;
    double gamma = 1 / sqrt(1 - beta * beta);
    return m * beta * gamma;
}

double MyCommonRICH::Momentum(double theta, double n, double m)
{
    return Momentum(BetaByTheta(theta, n), m);
}

//______________________________________________________________________________
/// photon function [nm]
double MyCommonRICH::PhotonEng(double lambda) //nm
{
    if (lambda == 0)
        return 0;

    return 1241.53 / lambda; //eV
}

double MyCommonRICH::PhotonWavelength(double energy) //eV
{
    if (energy == 0)
        return 0;

    return 1241.53 / energy; //nm
}

//______________________________________________________________________________
/// Cherenkov angle
double MyCommonRICH::CherenkovCosAng(double beta, double n)
{
    if (beta == 0 || n == 0)
        return 1;

    double cang = 1 / (beta * n);

    if (cang > 1)
        return 1;
    else
        return cang;
}

double MyCommonRICH::CherenkovCosAng(double p, double m, double n) //MeV
{
    return CherenkovCosAng(Beta(p, m), n);
}

double MyCommonRICH::CherenkovAng(double beta, double n)
{
    return acos(CherenkovCosAng(beta, n));
}

double MyCommonRICH::CherenkovAng(double p, double m, double n)
{
    return acos(CherenkovCosAng(p, m, n));
}

//______________________________________________________________________________
/// Full reflective angle
double MyCommonRICH::FullReflectiveSinAng(double n1, double n2) //from n1 to n2, (n1>n2)
{
    if (n1 < n2 || n2 == 0. || n1 < 0. || n2 < 0.)
        return 1.;
    return n2 / n1;
}

double MyCommonRICH::FullReflectiveCosAng(double n1, double n2) //from n1 to n2
{
    return sqrt(1. - FullReflectiveSinAng(n1, n2) * FullReflectiveSinAng(n1, n2));
}

double MyCommonRICH::FullReflectiveAng(double n1, double n2) //from n1 to n2
{
    return acos(FullReflectiveCosAng(n1, n2));
}

double MyCommonRICH::CherenkovMaxMoment(double n1, double n2, double m) //the maximum momentum particle(with mass as m in [GeV]) allowed to have cherenkov coming from n1 to n2
{
    double fang = FullReflectiveCosAng(n1, n2);
    if (fang == 0. || n1 <= 0.)
        return -1.;

    double beta = 1. / n1 / fang;
    if (beta > 1.)
        return -1;

    double gamma = 1. / sqrt(1 - beta * beta);
    return m * beta * gamma;
}

//______________________________________________________________________________
/// Cherenkov photon numbers per length [/cm] per Lambda [/nm] (default is /m/m)
double MyCommonRICH::dNdLdLambda(double lambda, double beta, double n) // [nm, GeV, nan]
{
    if (lambda == 0 || beta == 0 || n == 0)
        return 0;

    if (1 / (beta * n) > 1)
        return 0; // no Cherenkov light generated.

    return 2. * pi * 1.E6 * (1 / 137.) / (lambda * lambda) * (1 - 1 / beta / beta / n / n);
}

//______________________________________________________________________________
/// support functions
vector<double> MyCommonRICH::GetThickList(MyRICHDetector *det, int id)
{
    vector<double> Tlist;
    for (int i = id; i < det->nRadLayer; i++)
        Tlist.push_back(det->tRadLayer[i]);
    Tlist.push_back(det->tTransLayer);
    return Tlist;
}

vector<double> MyCommonRICH::GetRefIndList(MyRICHDetector *det, int id, double lambda)
{
    vector<double> Nlist;
    for (int i = id; i < det->nRadLayer; i++)
        Nlist.push_back(gDb->GetMatRefValue(det->sRadLayer[i], lambda));
    Nlist.push_back(gDb->GetMatRefValue(det->sTransLayer[0], lambda));
    return Nlist;
}

vector<double> MyCommonRICH::GetAbsLenList(MyRICHDetector *det, int id, double lambda) //[mm]
{
    vector<double> Alist;
    for (int i = id; i < det->nRadLayer; i++)
        Alist.push_back(10 * gDb->GetMatAbsValue(det->sRadLayer[i], lambda, 0));

    for (int i = 0; i < det->nImpurities; i++)
        Alist.push_back(10 * gDb->GetMatAbsValue(det->sImpurities[i], lambda, det->pImpurities[i]));

    return Alist;
}

void MyCommonRICH::UpdateThetaCExp(MyRICHDetector *det, double lambda)
{
    vector<double> rlist = GetRefIndList(det, 0, lambda);

    det->fThetaCReal.clear();
    for (int i = 0; i < det->nRadLayer; i++)
        det->fThetaCReal.push_back(CherenkovAng(Beta(det->momentum, det->mass), rlist[i]));
}

void MyCommonRICH::SaveRings(const char *fname)
{
    TFile f(fname, "recreate");
    if (!f.IsOpen())
        return;

    // save detector infor.
    TTree *T = new TTree("TTree", "NumClass");
    MyRICHDetector *det = new MyRICHDetector();
    T->Branch("det", "MyRICHDetector", &det, 20000, 1);

    //1. gDet
    det = gDet; //->BuildFrom(*gDet, -1, 1);
    T->Fill();

    //2. gDetList
    for (int i = 0; i < (int)gDetList.size(); i++)
    {
        det = gDetList[i]; //->BuildFrom(*gDetList[i]);
        T->Fill();
    }

    //3. gScanDetList
    for (int imom = 0; imom < (int)gScanDetList.size(); imom++)
        for (int ithe = 0; ithe < (int)gScanDetList[imom].size(); ithe++)
            for (int ihypo = 0; ihypo < (int)gScanDetList[imom][ithe].size(); ihypo++)
            {
                det = gScanDetList[imom][ithe][ihypo]; //->BuildFrom(*gScanDetList[imom][ithe][ihypo]);
                T->Fill();
            }

    T->Write();
    cout << "----> Current results are stored to " << fname << "." << endl;
}

void MyCommonRICH::LoadRings(const char *fname)
{
    TH1::AddDirectory(kFALSE);
    TFile f(fname);
    if (!f.IsOpen())
        return;

    MyRICHDetector *det = new MyRICHDetector(0);
    TTree *T = (TTree *)f.Get("TTree");
    T->SetBranchAddress("det", &det);
    cout << "Total entries: " << T->GetEntries() << endl;

    //1. gDet
    int ientry = 0;
    T->GetEntry(ientry++);
    gDirectory->cd();
    gDet->BuildFrom(*det, -1, 1);

    //2. gDetList
    ResizeDetList(gDet->nhypo);
    for (int i = 0; i < gDet->nhypo; i++)
    {
        T->GetEntry(ientry++);
        gDetList[i] = new MyRICHDetector(*det, -1, 1);
    }

    //3. gScanDetList
    ResizeScanDetList(gDet->np, gDet->nthe0, gDet->nhypo);

    for (int imom = 0; imom < gDet->np; imom++)
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
            for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            {
                T->GetEntry(ientry++);
                gScanDetList[imom][ithe][ihypo] = new MyRICHDetector(*det, -1, 1);
            }

    f.Close();
}

/*
double MyCommonRICH::ProjectToPixel(MyRICHDetector *det, double xmin, double xmax, double ymin, double ymax, double lambda)
{
    double the0 = det->theta0;
    double beta = Beta(det->momentum, det->mass);

    int nx = 10, ny = 50;
    double xc = (xmin + xmax) / 2;
    double yc = (ymin + ymax) / 2;
    double xs = (xmax - xmin) / nx;
    double ys = (ymax - ymin) / ny;

    double w = 0;

    // loop along the track in each radiator layer
    for (int i = 0; i < 1; i++) //det->nRadLayer; i++)
    {
        thklist.clear();
        reflist.clear();
        thklist = GetThickList(det, i);
        reflist = GetRefIndList(det, i, lambda);

        double n = reflist[0];
        double thec = CherenkovAng(beta, n);
        double density = dNdLdLambda(lambda, beta, n);
        fRhoFcn->SetParameters(0, the0, thec);

        if (findZ0Solver(fRhoFcn, xc, yc, 0) * findZ0Solver(fRhoFcn, xc, yc, thklist[0]) > 0)
            continue;

        vector<double> l0list;
        for (int ii = 0; ii < nx; ii++)
        {
            for (int jj = 0; jj < ny; jj++)
            {
                double xcc = xmin + ii * xs;
                double ycc = ymin + jj * ys;
                double lc0 = findZ0(fRhoFcn, xcc, ycc, 0, thklist[0], epsilon);
                double phi = atan((ycc - lc0 * tan(the0)) / xcc);
                double xcd = (xcc > 0) ? xcc + xs * cos(phi) : xcc - xs * cos(phi);
                double ycd = (ycc > lc0 * tan(the0)) ? ycc + ys * sin(phi) : ycc - ys * sin(phi);
                double lc1 = findZ0(fRhoFcn, xcd, ycd, 0, thklist[0], epsilon);

                if (lc0 < 0 || lc1 < 0)
                    continue;
                double dphi = sqrt(xs * xs + ys * ys) / sqrt(xcc * xcc + pow(ycc - lc0 * tan(the0), 2));
                double dL = fabs(lc0 - lc1);

                // absorption length

                // total
                w += dL * dphi / 2 / TMath::Pi(); // * density;
                //cout << density << " " << dL << " " << dphi << endl;
            }
        }
    }
    return w;
}

double MyCommonRICH::ProjectToPixel(MyRICHDetector *det, double xmin, double xmax, double ymin, double ymax)
{
    double w1 = ProjectToPixel(det, xmin, xmax, ymin, ymax, det->lambdaMin);
    double w2 = ProjectToPixel(det, xmin, xmax, ymin, ymax, det->lambdaMin);

    if (w1 + w2 == 0)
        return 0;

    // loop in lambda range
    double w = 0;
    for (int i = 0; i < 1; i++) //det->nLambda; i++)
    {
        double lambda = det->lambdaMin + (i + 0.5) * det->lambdaStep;
        double val = ProjectToPixel(det, xmin, xmax, ymin, ymax, lambda);
        w += val; // * det->lambdaStep;
    }
    return w;
}
*/

//______________________________________________________________________________
/// 计算从某个辐射体发出来的光子数
double MyCommonRICH::PhotonGenFromRad(MyRICHDetector *det, int idet, int absFlag)
{
    if (det == 0)
        return 0;

    //1.初始化
    double beta = Beta(det->momentum, det->mass);
    double tbase = det->CalZ0(idet);
    double theta0 = det->theta0;
    thklist = GetThickList(det, idet);

    double detPh = 0;

    //2. 在波长范围循环
    for (int ilamb = 0; ilamb < det->nLambda; ilamb++)
    {
        double lambda = det->lambdaMin + (ilamb + 0.5) * det->lambdaStep;
        double QE = gDb->GetDetQEValue(det->PCDetector, lambda);
        if (QE == 0)
            continue;

        //2.1 根据波长初始化变量
        reflist = GetRefIndList(det, idet, lambda);
        abslist = GetAbsLenList(det, idet, lambda);
        double n0 = reflist[0];
        double t0 = thklist[0];

        double thetac = CherenkovAng(beta, n0);
        double density = dNdLdLambda(lambda, beta, n0);
        fRhoFcn->SetParameters(0, theta0, thetac);

        //2.2 在辐射体内循环发光点位置
        for (int ilen = 0; ilen < (int)t0 / det->trkStep; ilen++)
        {
            double z0 = (ilen + 0.5) * det->trkStep;
            fRhoFcn->SetParameter(0, z0);

            double Z0 = -1 * tbase - z0; //全局坐标系
            double X0 = 0;
            double Y0 = Z0 * tan(theta0);

            //3. 在出射角范围内循环
            for (int iphi = 0; iphi < det->nphi; iphi++)
            {
                // 3. findRho是在极坐标系下求解，此时极坐标系下的中心点为发光点位置，因此在计算时用的有效发光长度 l0 = t0 - z0
                //                  极坐标与+X夹角为phi，求解的结果转换到全局坐标系有：
                //                      Xr = rho * cos(phi); Yr = rho * sin(phi) + Y0
                double phi = iphi * det->phiStep;
                double rho = fRhoFcn->Eval(phi);
                if (rho == 0 || rho != rho)
                    continue;

                // 全局坐标系下的圆心坐标
                double Xr = rho * cos(phi) + X0;
                double Yr = rho * sin(phi) + Y0;

                //  产生的光子数
                double weight = density * det->trkStep * det->phiStep * det->lambdaStep * QE / 2 / TMath::Pi();

                // 考虑辐射体的吸收
                double coeff = 1.;
                for (int ii = 0; ii < (int)lenlist.size(); ii++)
                {
                    if (absFlag == 0)
                        continue;

                    double len = lenlist[ii];
                    double abs = abslist[ii];
                    if (ii == (int)lenlist.size() - 1)
                        for (int jj = ii; jj < (int)abslist.size(); jj++)
                            coeff *= exp(-1 * len / abslist[jj]);
                    else
                        coeff *= exp(-1 * len / abs);
                }
                weight *= coeff;
                detPh += weight;

                det->GetHitMap()->Fill(Xr, Yr, weight);
                det->GetDetHitMap(idet)->Fill(Xr, Yr, weight);
                det->GetWaveLengthHist()->Fill(lambda, weight);
            }
        }
    }
    return detPh;
}

/// 计算所有辐射体发射出来的光子数分布图
void MyCommonRICH::GenerateDetRing(MyRICHDetector *det)
{
    if (det == 0)
        det = gDet;

    //清除已经生成的图像
    det->GetHitMap(1);
    det->GetWaveLengthHist(1);
    det->Gen2DRingListForEachRad();

    double totPh = 0;
    for (int idet = 0; idet < det->nRadLayer; idet++)
    {
        double detPh = PhotonGenFromRadWithAbs(det, idet);
        cout << "-- Detector " << idet << " " << det->sRadLayer[idet] << " generates " << detPh << " ph." << endl;

        totPh += detPh;
    }
    cout << "Total photon from " << det->particle << ", p=" << det->momentum << "GeV/c, theta0=" << det->Theta0 << " : nPh=" << totPh << endl;
}

/// 生成mu/pi/k/p四种粒子的光子数分布图
void MyCommonRICH::GenerateMultiParticleRICHRings()
{
    ResizeDetList(gDet->nhypo);
    cout << "\n---------------------------";
    for (int i = 0; i < gDet->nhypo; i++)
    {
        cout << "\n-- generating: " << endl;
        gDetList[i] = new MyRICHDetector(*gDet, i + 1);
        gDetList[i]->SetParticleGun(SHYPO[i], GetMass(SHYPO[i]));
        GenerateDetRing(gDetList[i]);
    }
}

/// 根据momentum/theta范围生成四种粒子的光子数分布图
void MyCommonRICH::GenerateTheScanHitMapsForEachDetector()
{
    ResizeScanDetList(gDet->np, gDet->nthe0, gDet->nhypo);

    cout << "\n---------------------------";
    int id = 10;
    for (int imom = 0; imom < gDet->np; imom++)
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
            for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            {
                double mom = gDet->pMin + imom * gDet->pStep;
                double Theta0 = gDet->The0Min + ithe * gDet->The0Step;
                cout << "\n-- generating: [" << imom << " " << ithe << " " << ihypo << "]" << endl;
                gScanDetList[imom][ithe][ihypo] = new MyRICHDetector(*gDet, id++);
                gScanDetList[imom][ithe][ihypo]->SetParticleGun(SHYPO[ihypo], GetMass(SHYPO[ihypo]));
                gScanDetList[imom][ithe][ihypo]->SetParticleGun(mom, Theta0);
                GenerateDetRing(gScanDetList[imom][ithe][ihypo]);
            }
}

void MyCommonRICH::GenerateTheNPhotonMap()
{
    ResizeNPhMap(gDet->nhypo, gDet->nRadLayer);

    for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
    {
        fNPhMap[ihypo] = new TH2F(Form("fNPhMap%d", ihypo), Form("Number of photon distribution for %s", SHYPO[ihypo]), gDet->np, gDet->pMin, gDet->pMax, gDet->nthe0, gDet->The0Min, gDet->The0Max);
        fNPhMap[ihypo]->GetXaxis()->SetTitle("momentum[GeV/c]");
        fNPhMap[ihypo]->GetYaxis()->SetTitle("#theta[degree]");

        for (int irad = 0; irad < gDet->nRadLayer; irad++)
        {
            fNPhMapEachRad[ihypo][irad] = new TH2F(Form("fNPhMap%d_%d", ihypo, irad), Form("Number of photon distribution for %s from %s", SHYPO[ihypo], gDet->sRadLayer[irad].c_str()), gDet->np, gDet->pMin, gDet->pMax, gDet->nthe0, gDet->The0Min, gDet->The0Max);
            fNPhMapEachRad[ihypo][irad]->GetXaxis()->SetTitle("momentum[GeV/c]");
            fNPhMapEachRad[ihypo][irad]->GetYaxis()->SetTitle("#theta[degree]");
        }

        for (int imom = 0; imom < gDet->np; imom++)
            for (int ithe = 0; ithe < gDet->nthe0; ithe++)
            {
                fNPhMap[ihypo]->SetBinContent(imom + 1, ithe + 1, gScanDetList[imom][ithe][ihypo]->GetHitMap()->Integral());
                for (int irad = 0; irad < gDet->nRadLayer; irad++)
                {
                    double nph = 0;
                    if (gScanDetList[imom][ithe][ihypo]->GetDetHitMap(irad) != NULL)
                        nph = gScanDetList[imom][ithe][ihypo]->GetDetHitMap(irad)->Integral();
                    fNPhMapEachRad[ihypo][irad]->SetBinContent(imom + 1, ithe + 1, nph);
                }
            }
    }
}

//______________________________________________________________________________
/// reconstruct functions

//-----------ALICE的beta重建方法-------------------------------------------
//  要调用这个函数，需要将下面三个vector都赋值
//    double hypoLambda = 185;
//    UpdateThetaCExp(det, hypoLambda);
//    thklist = GetThickList(gDet, irad);
//    reflist = GetRefIndList(gDet, irad, hypoLambda); //折射率选用假设的波长处的折射率
//    abslist = GetAbsLenList(gDet, irad, hypoLambda);
//    double recThe = ReconstructRICHByBeta(gDet, irad, Xr, Yr);
//
double MyCommonRICH::ReconstructRICHByBeta(MyRICHDetector *det, double irad, double Xc, double Yc)
{
    if (det == 0)
        return -999;

    //1. 初始化
    double theta0 = det->theta0;
    double thetac = det->fThetaCReal[irad];
    double beta = Beta(det->momentum, det->mass);
    double z0 = thklist[0] / 2;
    double Xep = thklist[0] - z0;
    double Tg = det->tTransLayer;

    //2. 若theta0不为0，则将击中点进行变换到theta0=0下进行求解
    double R = sqrt(Xc * Xc + Yc * Yc);

    if (theta0 != 0)
    {
        double tbase = det->CalZ0(irad);
        double Z0 = -1 * tbase - z0;
        double X0 = 0;
        double Y0 = Z0 * tan(theta0);
        double phi = atan((Yc - Y0) / (Xc - X0));
        fRhoFcn->SetParameters(z0, theta0, thetac);
        double a1 = fRhoFcn->Eval(phi);
        fRhoFcn->SetParameters(z0, 0, thetac);
        double a2 = fRhoFcn->Eval(phi);
        R = R * a2 / a1;
        //cout<<"In rec: phi="<<phi<<", theta0="<<theta0<<" thetac="<<thetac<<endl;
        //cout<<"   cal: a2="<<a2<<"; a1="<<a1<<endl;
        //cout<<"   Yc="<<Yc<<" Xc="<<Xc<<" Y0="<<Y0<<" X0="<<X0<<endl;
    }

    //3. 求R0
    double Rrad = Xep * tan(thetac);
    double Rqz = 0;
    double nSinTheta = reflist[0] * sin(thetac);

    //！！注意！！全反射的情况就是简单将Rqz设为0，并没有返回-999
    //！！注意！！这样直接导致R0变大，从而导致重建的thetac变小，极限为acos(1/beta)。
    for (int i = 1; i < (int)thklist.size() - 1; i++)
    {
        if (nSinTheta / reflist[i] > 1)
            continue;
        Rqz += thklist[i] * tan(asin(nSinTheta / reflist[i]));
    }

    double R0 = R - Rrad - Rqz;

    //4. 求重建的rthec
    double rthec;
    rthec = acos(1 / sqrt(1 + beta * beta / (1 + Tg * Tg / R0 / R0)));

    return rthec;
}

//-----------利用求解投影a来重建的算法-------------------------------------------
//  对每一个[X, Y]的击中进行重建，假设：发光点在辐射体中心处 & 波长185nm
//  返回为假设从该辐射体中心处发185nm的光对应的切伦科夫辐射角。无解则为-999
//  要调用这个函数，需要将下面三个vector都赋值
//    double hypoLambda = 185;
//    thklist = GetThickList(gDet, irad);
//    reflist = GetRefIndList(gDet, irad, hypoLambda); //折射率选用假设的波长处的折射率
//    abslist = GetAbsLenList(gDet, irad, hypoLambda);
//    double recThe = ReconstructRICHBySolver(gDet, irad, Xr, Yr, thetac);
//
double MyCommonRICH::ReconstructRICHBySolver(MyRICHDetector *det, double irad, double Xc, double Yc)
{
    if (det == 0)
        return -999;

    //1. 假设发光点在中心处，计算phi
    double L = thklist[0];
    double tbase = det->CalZ0(irad);
    double thetaExp = det->fThetaCReal[irad];

    double z0, z0Prim;
    double X0, Y0, Z0;
    double phi1, phi2;

    z0 = L / 2.;
    Z0 = -1 * tbase - z0;
    X0 = 0;
    Y0 = Z0 * tan(det->theta0);
    phi1 = atan((Yc - Y0) / (Xc - X0));
    phi2 = 0;

    while (fabs(phi1 - phi2) > 0.001)
    {
        //2. 根据phi计算出射光的夹角Theta
        double a = findA(det->theta0, thetaExp, phi1);
        double Theta = atan(a / z0);

        //3. 计算发光点的期望值
        double absl = abslist[0] * cos(Theta);
        z0Prim = L - (absl - absl * exp(L / absl) + L) / (1 - exp(L / absl));

        //4. 更新期望的发光点位置Z0，这里phi也发生了变化
        Z0 = -1 * tbase - z0Prim; //update  Z0
        Y0 = Z0 * tan(det->theta0);

        phi2 = atan((Yc - Y0) / (Xc - X0));
    }
    cout << "z0 is now: " << z0Prim << endl;

    fRhoFcn->SetParameters(z0, det->theta0);

    return findThetaC(fRhoFcn, Xc, Yc, X0, Y0, thetaExp, epsilon);
}

// 根据辐射体的光击中分布，求解中心值偏差及展宽的分布图
void MyCommonRICH::ReconstructRICHBySolver(MyRICHDetector *det)
{
    if (det == 0)
        det = gDet;

    if (det->fHitMap == 0)
        return;

    //1. 初始化
    double hypoLambda = 185;
    det->GetRecMap(1);
    det->Gen2DRecRingListForEachRad();
    UpdateThetaCExp(det, hypoLambda);

    //2. 在辐射体产生的光子分布里循环
    double xc, yc;
    for (int irad = 0; irad < det->nRadLayer; irad++)
    {
        thklist.clear();
        reflist.clear();
        thklist = GetThickList(det, irad);
        reflist = GetRefIndList(det, irad, hypoLambda); //折射率选用假设的波长处的折射率
        abslist = GetAbsLenList(det, irad, hypoLambda); //吸收系数用假设的波长处的折射率
        double muPhoton = det->fHitMapEachRad[irad]->Integral();
        if (muPhoton == 0)
            continue;
        //cout << "-->" << det->particle << " " << det->momentum << " " << det->theta0 << " rad:" << irad << " ph=" << muPhoton << " realTh=" << det->fThetaCReal[irad] << endl;

        //2.1 光子期望值循环
        for (int iph = 1; iph < 50; iph++)
        {
            //2.2 进行NEvent个事例数
            for (int i = 0; i < nEvent; i++)
            {
                //2.3 单个Event里会有iph个光子
                double avgth = 0;
                for (int j = 0; j < iph; j++)
                {
                    det->fHitMapEachRad[irad]->GetRandom2(xc, yc);
                    xc = int(xc / det->pixel) * det->pixel + det->pixel / 2;
                    yc = int(yc / det->pixel) * det->pixel + det->pixel / 2;
                    double theta = ReconstructRICHByBeta(det, irad, xc, yc);
                    //double theta = ReconstructRICHBySolver(det, irad, xc, yc);
                    avgth += theta;
                }
                avgth /= iph;

                det->fRecMap->Fill(iph, avgth);
                det->fRecMapEachRad[irad]->Fill(iph, avgth);
                //if(i%50==0) cout << "   evt: " << i << " ph = " << iph << " theta = " << avgth << endl;
            }
        }
    }
}

/// 生成momentum/theta范围的重建中心值偏差及展宽的分布图, 必须先完成GenerateTheScanHitMapsForEachDetector
bool MyCommonRICH::ReconstructForEachDetector()
{
    if (GetDetScanNumber() == 0)
    {
        cout << "--> Must generate the hitmap first before performing the reconstruction." << endl;
        return false;
    }

    for (int imom = 0; imom < (int)gDet->np; imom++)
        for (int ithe = 0; ithe < (int)gDet->nthe0; ithe++)
            for (int ihypo = 0; ihypo < (int)gDet->nhypo; ihypo++)
                ReconstructRICHBySolver(gScanDetList[imom][ithe][ihypo]);

    return true;
}

void MyCommonRICH::GenerateRecOffsetSigmaMap()
{
    //
    ResizeRecMap(gDet->nhypo, gDet->nRadLayer, gDet->np, gDet->nthe0, gDet->NPhoton);

    // 1. 拟合重建结果
    for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
        for (int irad = 0; irad < gDet->nRadLayer; irad++)
            for (int imom = 0; imom < gDet->np; imom++)
                for (int ithe = 0; ithe < gDet->nthe0; ithe++)
                {
                    TH2F *fhitmap = gScanDetList[imom][ithe][ihypo]->fRecMapEachRad[irad];
                    TString particle = gScanDetList[imom][ithe][ihypo]->particle;
                    TString radiator = gScanDetList[imom][ithe][ihypo]->sRadLayer[irad];
                    double momentum = gScanDetList[imom][ithe][ihypo]->momentum;
                    double theta0 = gScanDetList[imom][ithe][ihypo]->Theta0;
                    for (int iph = 1; iph < gDet->NPhoton; iph++)
                    {
                        TH1D *fproj = fhitmap->ProjectionY(Form("fRecMap%d_%d_%d_%d_%d", ihypo, irad, imom, ithe, iph), iph + 1, iph + 1);
                        fproj->SetTitle(Form("Particle: %s [%.1f GeV/c, %.1f#circ]. %d photons from radiator %s", particle.Data(), momentum, theta0, iph, radiator.Data()));
                        if (fproj->Integral() > 0)
                            fproj->Fit("gaus", "Q");
                        double mean = (fproj->Integral() == 0) ? 0 : fproj->GetFunction("gaus")->GetParameter(1);
                        double sigm = (fproj->Integral() == 0) ? 0 : fproj->GetFunction("gaus")->GetParameter(2);
                        double merr = (fproj->Integral() == 0) ? 0 : fproj->GetFunction("gaus")->GetParError(1);
                        double serr = (fproj->Integral() == 0) ? 0 : fproj->GetFunction("gaus")->GetParError(2);
                        fRecOffList[ihypo][irad][imom][ithe][iph] = mean;
                        fRecSigList[ihypo][irad][imom][ithe][iph] = sigm;
                        fRecOffErrList[ihypo][irad][imom][ithe][iph] = merr;
                        fRecSigErrList[ihypo][irad][imom][ithe][iph] = serr;
                        fRecMap[ihypo][irad][imom][ithe][iph] = fproj;
                    }
                }
}

void MyCommonRICH::GenerateRecHistograms(TString particle, int irad, int imom, int ithe, int iph)
{
    //1. clear
    fOffsetMap = (TH2F *)gDirectory->Get("_offset_");
    if (fOffsetMap != 0)
        delete fOffsetMap;

    fSigmaMap = (TH2F *)gDirectory->Get("_sigma_");
    if (fSigmaMap != 0)
        delete fSigmaMap;

    for (int i = 0; i < (int)fOffsetVsNphPlot.size(); i++)
        delete fOffsetVsNphPlot[i];
    for (int i = 0; i < (int)fOffsetVsMomPlot.size(); i++)
        delete fOffsetVsMomPlot[i];
    for (int i = 0; i < (int)fOffsetVsThetaPlot.size(); i++)
        delete fOffsetVsThetaPlot[i];

    for (int i = 0; i < (int)fSigmaVsNphPlot.size(); i++)
        delete fSigmaVsNphPlot[i];
    for (int i = 0; i < (int)fSigmaVsMomPlot.size(); i++)
        delete fSigmaVsMomPlot[i];
    for (int i = 0; i < (int)fSigmaVsThetaPlot.size(); i++)
        delete fSigmaVsThetaPlot[i];

    fOffsetVsNphPlot.resize(gDet->nhypo);
    fOffsetVsMomPlot.resize(gDet->nhypo);
    fOffsetVsThetaPlot.resize(gDet->nhypo);
    fSigmaVsNphPlot.resize(gDet->nhypo);
    fSigmaVsMomPlot.resize(gDet->nhypo);
    fSigmaVsThetaPlot.resize(gDet->nhypo);

    //2. define histograms
    int ihypo = GetHypoID(particle);
    double mom = gScanDetList[imom][ithe][0]->momentum;
    double the = gScanDetList[imom][ithe][0]->Theta0;
    TString rad = gDet->sRadLayer[irad];

    fOffsetMap = new TH2F("_offset_", Form("offset map for %s from radiator %s with %d photons", particle.Data(), rad.Data(), iph), gDet->np, gDet->pMin, gDet->pMax, gDet->nthe0, gDet->The0Min, gDet->The0Max);
    fSigmaMap = new TH2F("_sigma_", Form("sigma map for %s from radiator %s with %d photons", particle.Data(), rad.Data(), iph), gDet->np, gDet->pMin, gDet->pMax, gDet->nthe0, gDet->The0Min, gDet->The0Max);

    for (int i = 0; i < gDet->nhypo; i++)
    {
        fOffsetVsNphPlot[i] = new TH1F(Form("_offvnph_%d", i), Form("Offset vs. Nphoton for %s[%.1fGeV/c, %.1f#circ] from radiator %s", particle.Data(), mom, the, rad.Data()), gDet->NPhoton, 0, gDet->NPhoton);
        fSigmaVsNphPlot[i] = new TH1F(Form("_sigvnph_%d", i), Form("Sigma vs. Nphoton for %s[%.1fGeV/c, %.1f#circ] from radiator %s", particle.Data(), mom, the, rad.Data()), gDet->NPhoton, 0, gDet->NPhoton);

        fOffsetVsMomPlot[i] = new TH1F(Form("_offvmom_%d", i), Form("Offset vs. Momentum for %s[%.1f#circ] from radiator %s using %d photons", particle.Data(), the, rad.Data(), iph), gDet->np, gDet->pMin, gDet->pMax);
        fSigmaVsMomPlot[i] = new TH1F(Form("_sigvmom_%d", i), Form("Sigma vs. Momentum for %s[%.1f#circ] from radiator %s using %d photons", particle.Data(), the, rad.Data(), iph), gDet->np, gDet->pMin, gDet->pMax);

        fOffsetVsThetaPlot[i] = new TH1F(Form("_offvthe_%d", i), Form("Offset vs. Momentum for %s[%.1fGeV/c] from radiator %s using %d photons", particle.Data(), mom, rad.Data(), iph), gDet->nthe0, gDet->The0Min, gDet->The0Max);
        fSigmaVsThetaPlot[i] = new TH1F(Form("_sigvthe_%d", i), Form("Sigma vs. Momentum for %s[%.1fGeV/c] from radiator %s using %d photons", particle.Data(), mom, rad.Data(), iph), gDet->nthe0, gDet->The0Min, gDet->The0Max);
    }

    //3. fill histograms
    for (int imom = 0; imom < gDet->np; imom++)
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
        {
            fOffsetMap->SetBinContent(imom + 1, ithe + 1, fRecOffList[ihypo][irad][imom][ithe][iph]);
            fOffsetMap->SetBinError(imom + 1, ithe + 1, fRecOffErrList[ihypo][irad][imom][ithe][iph]);
            fSigmaMap->SetBinContent(imom + 1, ithe + 1, fRecSigList[ihypo][irad][imom][ithe][iph]);
            fSigmaMap->SetBinError(imom + 1, ithe + 1, fRecSigErrList[ihypo][irad][imom][ithe][iph]);
        }

    for (int iihypo = 0; iihypo < gDet->nhypo; iihypo++)
    {
        for (int iiph = 0; iiph < gDet->NPhoton; iiph++)
        {
            fOffsetVsNphPlot[iihypo]->SetBinContent(iiph+1, fRecOffList[iihypo][irad][imom][ithe][iiph]);
            fOffsetVsNphPlot[iihypo]->SetBinError(iiph+1, fRecOffErrList[iihypo][irad][imom][ithe][iiph]);
            fSigmaVsNphPlot[iihypo]->SetBinContent(iiph+1, fRecSigList[iihypo][irad][imom][ithe][iiph]);
            fSigmaVsNphPlot[iihypo]->SetBinError(iiph+1, fRecSigErrList[iihypo][irad][imom][ithe][iiph]);
        }
        for (int iimom = 0; iimom < gDet->np; iimom++)
        {
            fOffsetVsMomPlot[iihypo]->SetBinContent(iimom+1, fRecOffList[iihypo][irad][iimom][ithe][iph]);
            fOffsetVsMomPlot[iihypo]->SetBinError(iimom+1, fRecOffErrList[iihypo][irad][iimom][ithe][iph]);
            fSigmaVsMomPlot[iihypo]->SetBinContent(iimom+1, fRecSigList[iihypo][irad][iimom][ithe][iph]);
            fSigmaVsMomPlot[iihypo]->SetBinError(iimom+1, fRecSigErrList[iihypo][irad][iimom][ithe][iph]);
        }
        for (int iithe = 0; iithe < gDet->nthe0; iithe++)
        {
            fOffsetVsThetaPlot[iihypo]->SetBinContent(iithe+1, fRecOffList[iihypo][irad][imom][iithe][iph]);
            fOffsetVsThetaPlot[iihypo]->SetBinError(iithe+1, fRecOffErrList[iihypo][irad][imom][iithe][iph]);
            fSigmaVsThetaPlot[iihypo]->SetBinContent(iithe+1, fRecSigList[iihypo][irad][imom][iithe][iph]);
            fSigmaVsThetaPlot[iihypo]->SetBinError(iithe+1, fRecSigErrList[iihypo][irad][imom][iithe][iph]);
        }
    }
}

//______________________________________________________________________________
//
void MyCommonRICH::DrawFCN(int irad)
{
    // 1. 出射一个真实光子
    double lambda = 165;
    double z0 = 0.1 * gDet->tRadLayer[irad];
    double phi = 10. * TMath::Pi() / 180.;

    double beta = Beta(gDet->momentum, gDet->mass);
    double tbase = gDet->CalZ0(irad);
    double theta0 = gDet->theta0;

    thklist = GetThickList(gDet, irad);
    reflist = GetRefIndList(gDet, irad, lambda);
    abslist = GetAbsLenList(gDet, irad, lambda);
    double n0 = reflist[0];
    double thetac = CherenkovAng(beta, n0);
    fRhoFcn->SetParameters(z0, theta0, thetac);

    double Z0 = -1 * tbase - z0; //全局坐标系
    double X0 = 0;
    double Y0 = Z0 * tan(theta0);

    double rho = fRhoFcn->Eval(phi);

    // 全局坐标系下的圆心坐标
    double Xr = rho * cos(phi) + X0;
    double Yr = rho * sin(phi) + Y0;

    cout << "光子出射的参数：z0=" << z0 << " phi=" << phi << " theta0=" << theta0 << " thetac=" << thetac << endl;
    cout << "             (X0,Y0,Z0)=(" << X0 << ", " << Y0 << "," << Z0 << ")" << endl;
    cout << "             (Xr,Yr   )=(" << Xr << ", " << Yr << " ), rho=" << rho << endl;

    // 2. 重建一个真实光子
    double hypoLambda = 165;

    thklist = GetThickList(gDet, irad);
    reflist = GetRefIndList(gDet, irad, hypoLambda); //折射率选用假设的波长处的折射率
    abslist = GetAbsLenList(gDet, irad, hypoLambda);
    UpdateThetaCExp(gDet, hypoLambda);

    TGraph *g = new TGraph();
    for (int i = 0; i < 90; i++)
    {
        double thetac = TMath::Pi() / 2 / 90 * i;
        fRhoFcn->SetParameter(2, thetac);
        g->SetPoint(i, thetac, fRhoFcn->Eval(phi));
    }
    g->Draw("apl");

    double recThe1 = ReconstructRICHByBeta(gDet, irad, Xr, Yr);
    double recThe2 = 0; //ReconstructRICHBySolver(gDet, irad, Xr, Yr);
    cout << "光子重建的参数方法-Beta：   thetac1=" << recThe1 << endl;
    cout << "光子重建的参数方法-Solver： thetac2=" << recThe2 << " " << epsilon << endl;
    /* 
    thklist = GetThickList(gDet, irad);
    reflist = GetRefIndList(gDet, irad, hypoLambda); //折射率选用假设的波长处的折射率
    fRhoFcn->SetParameters(thklist[0] / 2., gDet->theta0);
    TH2F *f = new TH2F("ftmp", "fmp", 90, 0, TMath::Pi() / 2, 360, -TMath::Pi(), TMath::Pi());
    for (int i = 0; i < 90; i++)
    {
        double thetac = f->GetXaxis()->GetBinCenter(i);
        fRhoFcn->SetParameter(2, thetac);
        for (int j = 0; j < 360; j++)
        {
            double phi = f->GetYaxis()->GetBinCenter(j);
            f->SetBinContent(i, j, fRhoFcn->Eval(phi));
        }
    }
    f->GetXaxis()->SetTitle("thetac");
    f->GetYaxis()->SetTitle("phi");
    f->Draw("cont");
    */

    /* 
    TH2F *f = new TH2F("ftmp", "fmp", 90, 0, TMath::Pi() / 2, 100, 0, thklist[0]);
    for (int i = 0; i < 90; i++)
    {
        double thetac = f->GetXaxis()->GetBinCenter(i);
        fRhoFcn->SetParameter(2, thetac);
        for (int j = 0; j < 100; j++)
        {
            double z0 = f->GetYaxis()->GetBinCenter(j);
            fRhoFcn->SetParameter(0, z0);
            f->SetBinContent(i, j, fRhoFcn->Eval(0));
        }
    }
    f->GetXaxis()->SetTitle("thetac");
    f->GetYaxis()->SetTitle("L");
    f->Draw("cont");
    */

    /* 
    double thetac = 0.53;
    fRhoFcn->SetParameter(2, thetac);

    TH2F *f = new TH2F("ftmp", "fmp", 360, -TMath::Pi(), TMath::Pi(), 100, 0, thklist[0]);
    for (int i = 0; i < 360; i++)
    {
        double phi = f->GetYaxis()->GetBinCenter(i);
        for (int j = 0; j < 100; j++)
        {
            double z0 = f->GetYaxis()->GetBinCenter(j);
            fRhoFcn->SetParameter(0, z0);
            f->SetBinContent(i, j, fRhoFcn->Eval(phi));
        }
    }
    f->GetXaxis()->SetTitle("phi");
    f->GetYaxis()->SetTitle("L");
    f->Draw("cont");
    */
}