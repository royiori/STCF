#include "MyCommonRICH.h"
#include "MyGuiActionClass.h"

vector<double> reflist; //reflective index of each radiator list
vector<double> thklist; //thickness of each radiator list
vector<double> lenlist; //light path length in each radiator list
vector<double> abslist; //absorption length in each radiator list

const int NHYPO = 4;
const char *SHYPO[NHYPO] = {"mu", "pi", "k", "p"};

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
double findRho(Double_t *x, Double_t *par)
{
    double phi = x[0];

    double tt = thklist[0]; //t0 : 当前辐射体厚度
    double z0 = par[0];     // z0 : 发光点到当前辐射体上表面的距离
    double l0 = tt - z0;
    double theta0 = par[1];
    double thetac = par[2];
    //double verbose = par[3];

    double cc = cos(thetac);
    double s0 = sin(theta0);
    double c0 = cos(theta0);
    double t0 = tan(theta0);
    double sd = sin(phi);
    double c2c = cos(2 * thetac);
    double c20 = cos(2 * theta0);
    double c2d = cos(2 * phi);

    if (c0 == 0)
        return 0;

    double a = 1;
    a *= cc * pow(1 / c0, 2) * sqrt(1 + c20 - 2 * c2c - 2 * c2d * s0 * s0) + 2 * sd * t0;
    a /= (2 * cc * cc / c0 / c0 - 2 * sd * sd * t0 * t0);

    double nSinTheta = reflist[0] * (a / sqrt(a * a + 1));
    a *= l0;
    lenlist.clear();
    lenlist.push_back(sqrt(a * a + l0 * l0));

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
double findThetaC(TF1 *fcn1, double x1, double y1, double x0, double y0, double err) //根据阳极上X/Y，求解theta_c
{
    cout << "---------------" << endl;
    /* 
    double diff = -1;
    double thetamin = -1, thetamax = -1;
    double rad = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
    double phi = atan((y1 - y0) / (x1 - x0));

    for (int i = 1; i <= 18; i++)
    {
        double thec = i * TMath::Pi() / 36.;
        fcn1->SetParameter(2, thec);
        double rho = fcn1->Eval(phi);
        if (rho == 0)
            continue;

        if (i == -1)
        {
            diff = rho - rad;
            continue;
        }

        if(diff * (rho-rad) < 0) 
        {
            thetamin = thec - TMath::Pi()/36.;
            thetamax = thec;
            break;
        }

        diff = rho - rad;
    }
    */

    //if (findThetaCSolver(fcn1, x1, y1, z1, thetamin) * findThetaCSolver(fcn1, x1, y1, z1, thetamax) > 0)
    //    return -999;

    int WDog = 0;
    double thetamid;
    double thetaStep = TMath::Pi() / 360.;
    double thetamin = thetaStep;
    double thetamax = TMath::Pi() / 2.;

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

        if (val1 == 0)
        {
            thetamin += thetaStep;
            if (thetamin >= thetamax)
                break;
            continue;
        }
        if (val2 == 0)
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
}

MyCommonRICH::~MyCommonRICH()
{
    // Destructor.
}

//______________________________________________________________________________
/// Mass & Beta function
double MyCommonRICH::GetMass(TString p)
{
    p.ToLower();
    if (p == "mu" || p == "muon")
        return gkMassMuon;
    if (p == "pi" || p == "pion")
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

void MyCommonRICH::SaveRings(const char *fname)
{

    TFile f(fname, "recreate");
    if (!f.IsOpen())
        return;

    // save detector infor.
    TTree *T = new TTree("TTree", "NumClass");
    MyRICHDetector *det = new MyRICHDetector();
    T->Branch("det", "MyRICHDetector", &det, 8000, 2);

    //1. gDet
    det->BuildFrom(*gDet);
    T->Fill();

    //2. gDetList
    for (int i = 0; i < (int)gDetList.size(); i++)
    {
        det->BuildFrom(*gDetList[i]);
        T->Fill();
    }

    //3. gScanDetList
    for (int imom = 0; imom < (int)gScanDetList.size(); imom++)
        for (int ithe = 0; ithe < (int)gScanDetList[imom].size(); ithe++)
            for (int ihypo = 0; ihypo < (int)gScanDetList[imom][ithe].size(); ihypo++)
            {
                det->BuildFrom(*gScanDetList[imom][ithe][ihypo]);
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

    MyRICHDetector *det = 0;
    TTree *T = (TTree *)f.Get("TTree");
    T->SetBranchAddress("det", &det);
    cout << "Total entries: " << T->GetEntries() << endl;

    //1. gDet
    int ientry = 0;
    T->GetEntry(ientry++);
    gDirectory->cd();
    gDet->BuildFrom(*det);

    //2. gDetList
    for (int i = 0; i < (int)gDetList.size(); i++)
        delete gDetList[i];

    for (int i = 0; i < gDet->nhypo; i++)
    {
        T->GetEntry(ientry++);
        gDetList.push_back(new MyRICHDetector(*det));
    }

    //3. gScanDetList
    for (int imom = 0; imom < (int)gScanDetList.size(); imom++)
        for (int ithe = 0; ithe < (int)gScanDetList[imom].size(); ithe++)
            for (int ihypo = 0; ihypo < (int)gScanDetList[imom][ithe].size(); ihypo++)
                delete gScanDetList[imom][ithe][ihypo];

    gScanDetList.resize(gDet->np);
    for (int imom = 0; imom < gDet->np; imom++)
    {
        gScanDetList[imom].resize(gDet->nthe0);
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
        {
            gScanDetList[imom][ithe].resize(gDet->nhypo);
            for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            {
                T->GetEntry(ientry++);
                gScanDetList[imom][ithe][ihypo] = new MyRICHDetector(*det);
            }
        }
    }
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

/// 生成某个探测器能接收到的光子数分布图
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
    for (int i = 0; i < (int)gDetList.size(); i++)
        if (gDetList[i] != 0)
            delete gDetList[i];

    gDetList.resize(gDet->nhypo);
    for (int i = 0; i < gDet->nhypo; i++)
    {
        gDetList[i] = new MyRICHDetector(*gDet, i + 1);
        gDetList[i]->SetParticleGun(SHYPO[i], GetMass(SHYPO[i]));
        GenerateDetRing(gDetList[i]);
    }
}

/// 根据momentum/theta范围生成四种粒子的光子数分布图
void MyCommonRICH::GenerateTheScanHitMaps()
{
    for (int imom = 0; imom < (int)gScanDetList.size(); imom++)
        for (int ithe = 0; ithe < (int)gScanDetList[imom].size(); ithe++)
            for (int ihypo = 0; ihypo < (int)gScanDetList[imom][ithe].size(); ihypo++)
                delete gScanDetList[imom][ithe][ihypo];

    int id = 10;
    gScanDetList.resize(gDet->np);
    for (int imom = 0; imom < gDet->np; imom++)
    {
        gScanDetList[imom].resize(gDet->nthe0);
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
        {
            gScanDetList[imom][ithe].resize(gDet->nhypo);
            for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            {
                double mom = gDet->pMin + imom * gDet->pStep;
                double Theta0 = gDet->The0Min + ithe * gDet->The0Step;
                cout << "-- generating: " << endl;
                gScanDetList[imom][ithe][ihypo] = new MyRICHDetector(*gDet, id++);
                gScanDetList[imom][ithe][ihypo]->SetParticleGun(SHYPO[ihypo], GetMass(SHYPO[ihypo]));
                gScanDetList[imom][ithe][ihypo]->SetParticleGun(mom, Theta0);
                GenerateDetRing(gScanDetList[imom][ithe][ihypo]);
            }
        }
    }
}

/// 生成momentum/theta范围的重建中心值偏差及展宽的分布图, 必须先完成GenerateTheScanHitMaps
void MyCommonRICH::GenerateTheOffsetAndResolutionMaps()
{
    if (gScanDetList.size() == 0)
    {
        cout<<"--> Must generate the hitmap first before performing the reconstruction."<<endl;
        return;
    }

}

//______________________________________________________________________________
/// reconstruct functions
//     对每一个[X, Y]的击中进行重建，假设：发光点在辐射体中心处 & 波长185nm
//     返回为rthec，里面按辐射体编号排序，保存数据为假设从该辐射体中心处发185nm的光对应的切伦科夫辐射角。无解则为-999
/* 
vector<double> MyCommonRICH::ReconstructRICHBySolver(MyRICHDetector *det, double Xc, double Yc)
{
    vector<double> rthec;
    rthec.resize(det->nRadLayer);

    double the0 = det->theta0;
    double hypoLambda = 185;

    // loop in each radiator layer
    for (int irad = 0; irad < det->nRadLayer; irad++)
    {
        thklist.clear();
        reflist.clear();
        thklist = GetThickList(det, irad);
        reflist = GetRefIndList(det, irad, hypoLambda); //折射率选用假设的波长处的折射率

        double tbase = det->CalZ0(irad);
        double n = reflist[0];
        double z0 = thklist[0] / 2.; //z为辐射体中心处，辐射体坐标系
        double Z0 = -1 * tbase - z0; //全局坐标系
        double X0 = 0;
        double Y0 = Z0 * tan(det->theta0);

        fRhoFcn->SetParameters(z0, the0);

        rthec[irad] = findThetaC(fRhoFcn, Xc, Yc, X0, Y0, epsilon);
        cout << irad << ": real=" << CherenkovAng(Beta(det->momentum, det->mhypo[0]), n) << " rec=" << rthec[irad]<< endl;
    }
    return rthec;
}

// 生成重建theta_c的分布
void MyCommonRICH::ReconstructCerekovAngleDist(MyRICHDetector *det)
{
    GenerateMultiParticleRICHRings(det);

    double x, y;
    vector<double> rthec;
    TH2F **fRing = det->Get2DRingList();
    TH2F **fRecR = det->Get2DRecList();
    for (int ihypo = 0; ihypo < 1; ihypo++) //det->nhypo; ihypo++)
    {
        //double muPhoton = fRing[ihypo]->Integral();
        for (int i = 0; i < 1; i++) //Nevents
        {
            //if (i % 100 == 0)
            for (int nph = 2; nph < 3; nph++) //40; nph++)
            {
                double ravg = 0;
                for (int iph = 1; iph < nph; iph++)
                {
                    fRing[ihypo]->GetRandom2(x, y);
                    rthec.clear();
                    rthec = ReconstructRICHBySolver(det, x, y);
                    ravg += rthec[0]; //need to update
                }
                ravg /= nph;
                fRecR[ihypo]->Fill(nph, ravg);
            }
        }
    }
    cout << "---> Reconstruction is done." << endl;
}
*/

/* 
void MyCommonRICH::ReconstructRICHBySolver(MyRICHDetector *det)
{
    TH2F *fHitPos = det->Get2DViewer();
    TH2F *fHitRec = det->Get2DRecRing();
    double muPhoton = fHitPos->Integral();

    int Nevent = 1;
    double xc, yc;
    vector<vector<double>> dthec; //thetaC_rec - thetaC_hypo [ 第几个光子 ][ 这个光子可能是哪个辐射体的 ]

    fHitRec->Clear();

    for (int i = 0; i < Nevent; i++)
    {
        int nPhoton = gRandom->Poisson(muPhoton);

        dthec.clear();
        dthec.resize(nPhoton);
        for (int j = 0; j < nPhoton; j++)
        {
            fHitPos->GetRandom2(xc, yc);
            ReconstructRICHBySolver(det, hypo, xc, yc, dthec[j]);
        }

        // 对每个光子的重建结果，可能存在这个情况：
        // 1.--  多个辐射体都重建出一个角度。此时选择距离0最近的那个辐射体的重建结果
        // 2.--  将所有的辐射体结果都拿出来，按光子数进行平均，然后填图
        double rad = 0;
        for (int j = 0; j < nPhoton; j++)
        {
            double radmin = TMath::Pi();
            for (int k = 0; k < det->nRadLayer; k++)
                radmin = (fabs(dthec[j][k]) < fabs(radmin)) ? dthec[j][k] : radmin;
            if(radmin == TMath::Pi()) 
            {
                cout<<"This point is not used."<<endl;
                nPhoton--
            }
            rad += radmin;
        }

        fHitRec->Fill();
    }
}

void MyCommonRICH::ReconstructRICHBySolver(MyRICHDetector *det)
{
    if (det->Get2DViewer() == NULL)
        GenerateRICHRing(det);

    TH2F **frec = det->Get2DRecRing();

    for (int i = 0; i < 4; i++) //4 hypothesis
    {
        ReconstructRICHBySolver(det, det->hypo[i]);
    }
}
*/