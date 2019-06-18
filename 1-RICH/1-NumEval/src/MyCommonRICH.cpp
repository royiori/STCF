#include "MyCommonRICH.h"
#include "MyNumEvalClass.h"

vector<double> reflist; //reflective index of each radiator list
vector<double> thklist; //thickness of each radiator list
vector<double> lenlist; //light path length in each radiator list
vector<double> abslist; //absorption length in each radiator list

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
//
MyCommonRICH::MyCommonRICH()
{
    fRhoFcn = new TF1("fRhoFcn", findRho, 0, 2 * TMath::Pi(), 3);
    gDet = new MyRICHDetector();

    nphi = 360;
    trkStep = 0.01;
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

    gDet->Get2DViewer()->Write();
    return;

    // save detector infor.
    TTree *T = new TTree("TTree", "NumClass");
    MyRICHDetector *det = new MyRICHDetector();
    T->Branch("det", "MyRICHDetector", &det, 8000, 2);
    det->BuildFrom(*gDet);
    cout<<det->nphi<<endl;
    T->Fill();
    for (int i = 0; i < (int)gScanDetList.size(); i++)
    {
        det->BuildFrom(*gScanDetList[i]);
        T->Fill();
    }
    T->Write();
    cout << "----> Current results are stored to " << fname << "." << endl;
}

void MyCommonRICH::LoadRings(const char *fname)
{
    TFile f(fname);
    if (!f.IsOpen())
        return;

    MyRICHDetector *det = 0;
    TTree *T = (TTree *)f.Get("TTree");
    T->SetBranchAddress("det", &det);
    Int_t nentries = (Int_t)T->GetEntries();
    cout<<"----> Total "<<nentries<<" entries stored in root file."<<endl;

    T->GetEntry(0);
    gDet->BuildFrom(*det);
    ClearScanList();
    for(int i=1; i<nentries; i++)
    {
        T->GetEntry(i);
        gScanDetList.push_back(new MyRICHDetector(*det));
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

double MyCommonRICH::PhotonGenFromRad(MyRICHDetector *det, int idet, double lambda, TH2F *fRing, double scal, int absFlag)
{
    thklist = GetThickList(det, idet);
    reflist = GetRefIndList(det, idet, lambda);
    abslist = GetAbsLenList(det, idet, lambda);

    double phiStep = 2 * TMath::Pi() / nphi;
    double n0 = reflist[0];
    double t0 = thklist[0];
    double beta = Beta(det->momentum, det->mass);
    double theta0 = det->theta0;
    double thetac = CherenkovAng(beta, n0);
    double density = dNdLdLambda(lambda, beta, n0);
    fRhoFcn->SetParameters(0, theta0, thetac);

    double tbase = 0;
    for (int i = 0; i < idet; i++)
        tbase += det->tRadLayer[i];

    // 2. 粒子在当前辐射体（厚度t0）内的发光点的坐标为 z0(循环里的变量), 则在全局坐标系里有：Z0 = -1*t(上面所有层辐射体厚度) - z0
    //                  此发光点在X/Y平面的坐标为： X0 = 0, y0 = Z0 * tan(theta0)
    double nph = 0;
    for (int i = 0; i < (int)t0 / trkStep; i++)
    {
        double z0 = (i + 0.5) * trkStep;
        fRhoFcn->SetParameter(0, z0);

        double Z0 = -1 * tbase - z0; //全局坐标系
        double X0 = 0;
        double Y0 = Z0 * tan(det->theta0);

        for (int j = 0; j < nphi; j++)
        {
            // 3. findRho是在极坐标系下求解，此时极坐标系下的中心点为发光点位置，因此在计算时用的有效发光长度 l0 = t0 - z0
            //                  极坐标与+X夹角为phi，求解的结果转换到全局坐标系有：
            //                      Xr = rho * cos(phi); Yr = rho * sin(phi) + Y0
            double phi = j * phiStep;
            double rho = fRhoFcn->Eval(phi);
            if (rho == 0 || rho != rho)
                continue;

            double Xr = rho * cos(phi) + X0;
            double Yr = rho * sin(phi) + Y0;

            double weight = density * trkStep * phiStep / 2 / TMath::Pi();

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
            weight *= coeff * scal;
            nph += weight;

            fRing->Fill(Xr, Yr, weight);
        }
    }
    return nph;
}

void MyCommonRICH::GenerateRICHRing(MyRICHDetector *det)
{
    if (det == 0)
        det = gDet;

    cout << "----> Generating..." << endl;

    TH2F *fRing = det->Get2DViewer(1);

    double totPh = 0;
    for (int idet = 0; idet < det->nRadLayer; idet++)
    {
        double detPh = 0;

        for (int i = 0; i < det->nLambda; i++)
        {
            double lamStep = det->lambdaStep;
            double lambda = det->lambdaMin + (i + 0.5) * lamStep;
            double QE = gDb->GetDetQEValue(det->PCDetector, lambda);
            if (QE == 0)
                continue;

            double photon = PhotonGenFromRadWithAbs(det, idet, lambda, fRing, lamStep * QE);
            detPh += photon;
        }
        totPh += detPh;
        cout << "Detector " << idet << " " << det->sRadLayer[idet] << " generates " << detPh << " ph." << endl;
    }
    cout << "Total photon: " << totPh << endl;
}

void MyCommonRICH::GenerateMultiParticleRICHRings(vector<TString> parList)
{
    ClearScanList();
    for (int i = 0; i < (int)parList.size(); i++)
    {
        gScanDetList.push_back(new MyRICHDetector(*gDet));
        gScanDetList[i]->particle = parList[i];
        gScanDetList[i]->mass = GetMass(parList[i]);
        GenerateRICHRing(gScanDetList[i]);
    }
}

//______________________________________________________________________________
/// reconstruct functions

