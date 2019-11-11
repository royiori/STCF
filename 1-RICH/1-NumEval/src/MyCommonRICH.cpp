#include "MyCommonRICH.h"
#include "MyGuiActionClass.h"

//vector<double> reflist; //reflective index of each radiator list
//vector<double> thklist; //thickness of each radiator list
//vector<double> lenlist; //light path length in each radiator list
//vector<double> abslist; //absorption length in each radiator list

const char *SHYPO[NHYPO] = {"p", "k", "#pi", "#mu", "e"};

MyCommonRICH *gMyCommonRICH = (MyCommonRICH *)0;

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
// 调用TF1求解rho
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

    if (b < 0 || c == 0)
        return 0;

    a *= fabs(cc) / c0 / c0 * sqrt(b) + 2 * sp * t0;
    a /= c;
    return a;
}

double findRho(Double_t *x, Double_t *par)
{
    double phi = x[0];
    double z0 = par[0]; // z0 : 发光点到当前辐射体上表面的距离
    double theta0 = par[1];
    double thetac = par[2];

    vector<double> thklist;
    vector<double> reflist;
    vector<double> lenlist;
    int n = par[3];
    for (int i = 0; i < n; i++)
    {
        thklist.push_back(par[4 + 2 * i]);
        reflist.push_back(par[5 + 2 * i]);
    }

    double tt = thklist[0]; //t0 : 当前辐射体厚度
    double l0 = tt - z0;

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

    for (int i = 0; i < (int)lenlist.size(); i++)
        par[4 + 2 * n + i] = lenlist[i];

    return a;
}

void SetFcnParameters(TF1 *fcn, Double_t z0, Double_t theta0, Double_t thetac, vector<Double_t> thklist, vector<Double_t> reflist)
{
    fcn->SetParameter(0, z0);
    fcn->SetParameter(1, theta0);
    fcn->SetParameter(2, thetac);
    fcn->SetParameter(3, thklist.size());
    for (int i = 0; i < (int)thklist.size(); i++)
    {
        fcn->SetParameter(4 + i * 2, thklist[i]);
        fcn->SetParameter(5 + i * 2, reflist[i]);
    }
}

void GetFcnParameters(TF1 *fcn, vector<Double_t> &lenlist)
{
    int n = fcn->GetParameter(3);
    lenlist.clear();
    for (int i = 0; i < n; i++)
        lenlist.push_back(fcn->GetParameter(4 + 2 * n + i));
}

//______________________________________________________________________________
// 直接输入参数求解rho
double CalculateRho(Double_t phi, Double_t z0, Double_t theta0, Double_t thetac, vector<double> thklist, vector<double> reflist, vector<double> &lenlist)
{
    double tt = thklist[0]; //t0 : 当前辐射体厚度
    double l0 = tt - z0;
    double a = findA(theta0, thetac, phi);
    if (a == 0)
        return 0;

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
// Beta方法：重建theta_c，假设发光点在辐射体的中心处，且波长为185nm来求解
double BetaThetaC(double theta0, double theta1, double phi, double err) //根据theta1，求thetaC
{
    int WDog = 0;
    double thetaCmid;
    double thetaCStep = 0.001;
    double thetaCmin = thetaCStep;
    double thetaCmax = TMath::Pi() / 2.;

    while (WDog < 2000)
    {
        thetaCmid = (thetaCmin + thetaCmax) / 2.;

        double a1 = findA(theta0, thetaCmin, phi);
        double val1 = (a1 == 0) ? 0 : atan(a1) - theta1;

        double a2 = findA(theta0, thetaCmid, phi);
        double val2 = (a2 == 0) ? 0 : atan(a2) - theta1;

        if (val1 == 0) //无解
        {
            thetaCmin += thetaCStep;
            if (thetaCmin >= thetaCmax)
                break;
            continue;
        }
        if (val2 == 0) //无解
        {
            thetaCmax -= thetaCStep;
            if (thetaCmax <= thetaCmin)
                break;
            continue;
        }

        if (fabs(val2) < err)
            return thetaCmid;

        if (val1 * val2 > 0)
            thetaCmin = thetaCmid;
        else
            thetaCmax = thetaCmid;

        WDog++;
    }

    return -999;
}

//______________________________________________________________________________
// Solver方法：重建theta_c，假设发光点在辐射体的中心处，且波长为185nm来求解
/* 
double SolverThetaC(TF1 *fcn1, double x1, double y1, double x0, double y0, double thetaExp, double err) //根据阳极上X/Y，求解theta_c
{
    int WDog = 0;
    double thetamid;
    double thetaStep = 0.001;
    double thetamin = thetaStep;
    double thetamax = 1.2 * thetaExp; //1.0;

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
*/

//______________________________________________________________________________
//
//   MyCommonRICH
//
//______________________________________________________________________________
//
MyCommonRICH::MyCommonRICH()
{
    fNphFcn = new TF1("fNphFcn", "sqrt([0]/x+[1])", 1, 49);
    gDet = new MyRICHDetector(0);
    gDet->SetNHypothesis(NHYPO);
    epsilon = 1e-5;
    nEvent = 1e4;
    NThread = 4;
    Epoch = -1; //-1指扫描全范围
    fileName = TString("./tmp.root");
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
    if (p == "e")
        return 4;
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

int MyCommonRICH::GetMomID(double momentum)
{
    return (gDet->pStep == 0) ? 0 : int(round((momentum - gDet->pMin) / gDet->pStep));
}

int MyCommonRICH::GetThetaID(double Theta)
{
    return (gDet->The0Step == 0) ? 0 : int(round((Theta - gDet->The0Min) / gDet->The0Step));
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
    if (p == "e")
        return gkMassElectron;
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

//______________________________________________________________________________
/// 计算从某个辐射体发出来的光子数
double MyCommonRICH::PhotonGenFromRad(MyRICHDetector *det, int irad, int absFlag)
{
    if (det == 0)
        return 0;

    //1.初始化
    TF1 *fRhoFcn = new TF1("fRhoFcn", findRho, 0, 2 * TMath::Pi(), 100);

    vector<double> thklist;
    vector<double> reflist;
    vector<double> abslist;
    vector<double> lenlist;

    double beta = Beta(det->momentum, det->mass);
    double tbase = det->CalZ0(irad);
    double theta0 = det->theta0;
    thklist = GetThickList(det, irad);

    double detPh = 0;

    //2. 在波长范围循环
    for (int ilamb = 0; ilamb < det->nLambda; ilamb++)
    {
        double lambda = det->lambdaMin + (ilamb + 0.5) * det->lambdaStep;
        double QE = gDb->GetDetQEValue(det->PCDetector, lambda);
        if (QE == 0)
            continue;

        //2.1 根据波长初始化变量
        reflist = GetRefIndList(det, irad, lambda);
        abslist = GetAbsLenList(det, irad, lambda);
        double n0 = reflist[0];
        double t0 = thklist[0];

        double thetac = CherenkovAng(beta, n0);
        double density = dNdLdLambda(lambda, beta, n0);
        SetFcnParameters(fRhoFcn, 0, theta0, thetac, thklist, reflist);

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
                GetFcnParameters(fRhoFcn, lenlist);

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
                det->GetDetHitMap(irad)->Fill(Xr, Yr, weight);
                det->GetWaveLengthHist()->Fill(lambda, weight);
            }
        }
    }

    delete fRhoFcn;
    return detPh;
}

//---- 被GuiAction调用的函数
/// 计算所有辐射体发射出来的光子数分布图
void MyCommonRICH::GenerateDetRing(MyRICHDetector *det, int verbose)
{
    if (det == 0)
        det = gDet;

    //清除已经生成的图像
    det->GetHitMap(1);
    det->GetWaveLengthHist(1);
    det->Gen2DRingListForEachRad();

    double totPh = 0;
    for (int irad = 0; irad < det->nRadLayer; irad++)
    {
        double detPh = PhotonGenFromRadWithAbs(det, irad);
        if (verbose)
            cout << "-- Detector " << irad << " " << det->sRadLayer[irad] << " generates " << detPh << " ph." << endl;

        totPh += detPh;
    }

    if (verbose)
        cout << "Total photon from " << det->particle << ", p=" << det->momentum << "GeV/c, theta0=" << det->Theta0 << " : nPh=" << totPh << endl;
}

//---- 被GuiAction调用的函数
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

//---- 被GuiAction调用的函数 -- hitmap --
/// 根据momentum/theta范围生成四种粒子的光子数hitmap分布图
void GetMomentumScanRange(long ip, int &ibegin, int &iend)
{
    MyRICHDetector *gDet = gMyCommonRICH->GetDetector();

    int NTHREAD = gMyCommonRICH->GetNThread();
    int epoch = gMyCommonRICH->GetEpoch();
    int np = gDet->np;

    ibegin = -1;
    iend = -1;

    //0. ip的范围为【0，NTHREAD-1】
    //1. EPoch!=-1: 扫描Epoch指定的动量段
    if (epoch != -1)
    {
        if (epoch > np || epoch + 1 > np)
            return;
        ibegin = epoch;
        iend = epoch + 1;
        return;
    }

    //2. EPoch==-1: NTHREAD > np，如NTHREAD为4， np为2， 则只分2个线程。此时ip取值>=2的时候不扫描
    //              NTHREAD <= np, 分段扫描
    int nthread = (NTHREAD < np) ? NTHREAD : np;
    ibegin = ip * np / nthread;
    iend = (ip + 1) * np / nthread;
    if (ibegin > np || iend > np)
    {
        ibegin = -1;
        iend = -1;
    }
}

void *GenerateTheScanHitMapsHandler(void *ptr)
{
    long ip = (long)ptr;
    int ibegin, iend;

    MyRICHDetector *gDet = gMyCommonRICH->GetDetector();
    vector<vector<vector<MyRICHDetector *>>> gScanDetList = gMyCommonRICH->GetScanDetector();

    GetMomentumScanRange(ip, ibegin, iend);

    if (ibegin == -1 || iend == -1)
        return 0;

    for (int imom = ibegin; imom < iend; imom++)
    {
        // 1. 循环生成hitmap
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
            for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            {
                gMyCommonRICH->GenerateDetRing(imom, ithe, ihypo);
            }

        // 2. 保存root文件
        gMyCommonRICH->SaveHitFile(imom);
    }

    return 0;
}

void MyCommonRICH::GenerateDetRing(int imom, int ithe, int ihypo)
{
    double mom = gDet->pMin + imom * gDet->pStep;
    double Theta0 = gDet->The0Min + ithe * gDet->The0Step;
    cout << "-- generating : [" << imom << " " << ithe << " " << ihypo << "]: "
         << "mom=" << mom << "GeV/c, theta0=" << Theta0 << ", pid=" << SHYPO[ihypo] << endl;
    gScanDetList[imom][ithe][ihypo]->SetParticleGun(SHYPO[ihypo], gMyCommonRICH->GetMass(SHYPO[ihypo]));
    gScanDetList[imom][ithe][ihypo]->SetParticleGun(mom, Theta0);
    GenerateDetRing(gScanDetList[imom][ithe][ihypo], 0);
}

void MyCommonRICH::GenerateTheScanHitMapsForEachDetector()
{
    cout << "\n---------------------------\n";
    ResizeScanDetList(gDet->np, gDet->nthe0, gDet->nhypo);
    int nthread = (NThread < gDet->np) ? NThread : gDet->np;

    TThread *thread[nthread];

    for (int i = 0; i < nthread; i++)
        thread[i] = new TThread(Form("thit%d", i), GenerateTheScanHitMapsHandler, (void *)(size_t)i);

    for (int i = 0; i < nthread; i++)
        thread[i]->Run();

    for (int i = 0; i < nthread; i++)
        thread[i]->Join();

    cout << "----> Hit-map generated." << endl;
}

//---- 被GuiAction调用的函数 -- hitmap --
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
                if (gScanDetList[imom][ithe][ihypo]->GetHitMap() != NULL)
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
// 重建函数
double MyCommonRICH::Findz0(MyRICHDetector *det, int irad, double Xr, double Yr, vector<double> thklist, vector<double> abslist)
{
    double L = thklist[0];
    double tbase = det->CalZ0(irad);
    double theta0 = det->theta0;
    double thetaExp = det->fThetaCReal[irad];

    double z0, z0Prim;
    double X0, Y0, Z0;
    double phi1, phi2;

    z0 = L / 2.;
    z0Prim = z0;
    Z0 = -1 * tbase - z0;
    X0 = 0;
    Y0 = Z0 * tan(theta0);
    phi1 = FindPhi(det, Xr, Yr, X0, Y0);
    phi2 = 0;

    int wdog = 0;
    while (fabs(phi1 - phi2) > 0.001 && wdog < 100)
    {
        wdog++;
        //2. 根据phi计算出射光的夹角Theta
        double a = findA(theta0, thetaExp, phi1);
        double Theta = atan(a / z0);

        //3. 计算发光点的期望值
        double absl = abslist[0] * cos(Theta);
        z0Prim = L - (absl - absl * exp(L / absl) + L) / (1 - exp(L / absl));

        //4. 更新期望的发光点位置Z0，这里phi也发生了变化
        Z0 = -1 * tbase - z0Prim; //update  Z0
        Y0 = Z0 * tan(det->theta0);

        phi2 = FindPhi(det, Xr, Yr, X0, Y0);
    }
    if (wdog >= 100)
    {
        z0Prim = 0.52 * z0;
        //cout<<"###Warning: Can't find z0 for ("<<Xr<<", "<<Yr<<") hit, using 0.5*z instead."<<endl;
    }
    //cout << "z0 is now: " << z0Prim << endl;
    return z0Prim;
}

//求解phi
double MyCommonRICH::FindPhi(MyRICHDetector *det, double Xr, double Yr, double X0, double Y0)
{
    //double R = sqrt((Xr - X0) * (Xr - X0) + (Yr - Y0) * (Yr - Y0));

    double phi = atan((Yr - Y0) / (Xr - X0));

    //if (R > det->tTransLayer)
    {
        phi = (phi < 0) ? phi + 2 * TMath::Pi() : phi;
        if (Yr - Y0 < 0 && Xr - X0 < 0)
            phi += TMath::Pi();
        if (Yr - Y0 > 0 && Xr - X0 < 0)
            phi -= TMath::Pi();
    }
    //else //对于大于30degree的才可能出现
    //     // phi有可能出现差一个TMath::Pi()的情况
    //     // 而且对于入射角40度开始，辐射体上界面和下界面的出光点位置 与 光子击中位置重叠，极有可能误判象限
    //     // 最简单的方式就是两个都求解，然后找与R差别小的作为phi角
    //{
    //    double phi1 = phi;
    //    double phi2 = phi1 + TMath::Pi();
    //    double a11 = fRhoFcn->Eval(phi1);
    //    double a12 = fRhoFcn->Eval(phi2);
    //    phi = (fabs(a11 - R) < fabs(a12 - R)) ? phi1 : phi2;
    //}

    return phi;
}

//-----------ALICE的beta重建方法-------------------------------------------
//  要调用这个函数，需要将下面三个vector都赋值
//    double hypoLambda = 185;
//    UpdateThetaCExp(det, hypoLambda);
//    thklist = GetThickList(gDet, irad);
//    reflist = GetRefIndList(gDet, irad, hypoLambda); //折射率选用假设的波长处的折射率
//    abslist = GetAbsLenList(gDet, irad, hypoLambda); //折射率选用假设的波长处的折射率
//    double recThe = ReconstructRICHByBeta(gDet, irad, Xr, Yr, thklist, reflist, abslist);
//
double MyCommonRICH::ReconstructRICHByBeta(MyRICHDetector *det, double irad, double Xr, double Yr, vector<double> thklist, vector<double> reflist, vector<double> abslist)
{
    if (det == 0)
        return -999;

    //1. 初始化
    double theta0 = det->theta0;
    double thetac = det->fThetaCReal[irad];
    double beta = Beta(det->momentum, det->mass);
    double Tg = det->tTransLayer;
    double z0 = Findz0(det, irad, Xr, Yr, thklist, abslist);
    double Xep = thklist[0] - z0;

    //2. 求解出光角度phi
    double tbase = det->CalZ0(irad);
    double Z0 = -1 * tbase - z0;
    double X0 = 0;
    double Y0 = Z0 * tan(theta0);
    double phi = FindPhi(gDet, Xr, Yr, X0, Y0);
    double R = sqrt((Xr - X0) * (Xr - X0) + (Yr - Y0) * (Yr - Y0));

    //fRhoFcn->SetParameters(z0, theta0, thetac);
    //fRhoFcn->SetParameters(z0, 0, thetac);
    //double a2 = fRhoFcn->Eval(phi);
    //cout << "重建: z0=" << z0 << " phi=" << phi << ", theta0=" << theta0 << " thetac=" << thetac << endl;
    //cout << "     (X0,Y0,Z0)=(" << X0 << ", " << Y0 << "," << Z0 << ")" << endl;
    //cout << "     (Xr,Yr   )=(" << Xr << ", " << Yr << " )" << endl;

    double a1 = findA(theta0, thetac, phi);
    double theta1 = atan(a1);

    //3. 求R0
    double Rrad = Xep * tan(theta1);
    double Rqz = 0;
    double nSinTheta = reflist[0] * sin(theta1);

    for (int i = 1; i < (int)thklist.size() - 1; i++)
    {
        if (nSinTheta / reflist[i] > 1)
            return -999; //全反射

        Rqz += thklist[i] * tan(asin(nSinTheta / reflist[i]));
    }

    double R0 = R - Rrad - Rqz;

    /*
    //4. 求重建的rec-theta1
    double rthec1;
    rthec1 = acos(1 / sqrt(1 + beta * beta / (1 + Tg * Tg / R0 / R0)));

    //5. theta1->thetac
    double rthec;
    rthec = acos(cos(theta0) * cos(rthec1) + sin(theta0) * sin(rthec1) * sin(phi));

    return rthec;
    //return BetaThetaC(theta0, rthec1, phi, epsilon);
    */

    //4. 求重建的n
    double theta2 = atan(R0 / Tg);
    double nquartz2 = sin(theta2) / sin(theta1);
    double rthec = acos(1 / nquartz2 / beta);
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
/* 
double MyCommonRICH::ReconstructRICHBySolver(MyRICHDetector *det, double irad, double Xc, double Yc)
{
    if (det == 0)
        return -999;

    vector<double>thklist;
    vector<double>reflist;
    vector<double>abslist;

           thklist.clear();
        reflist.clear();
        abslist.clear();
        thklist = GetThickList(det, irad);
        reflist = GetRefIndList(det, irad, hypoLambda); //折射率选用假设的波长处的折射率
        abslist = GetAbsLenList(det, irad, hypoLambda); //吸收系数用假设的波长处的折射率


    //1. 假设发光点在中心处，计算phi
    double L = thklist[0];
    double tbase = det->CalZ0(irad);
    double thetaExp = det->fThetaCReal[irad];

    double z0, z0Prim = 0.52 * L;
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
        double Theta = atan(a / 1.);

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

    return SolverThetaC(fRhoFcn, Xc, Yc, X0, Y0, thetaExp, epsilon);
}
*/

// 根据辐射体的光击中分布，求解中心值偏差及展宽的分布图
void MyCommonRICH::ReconstructRICHDetector(MyRICHDetector *det)
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

    vector<double> thklist;
    vector<double> reflist;
    vector<double> abslist;

    //2. 在辐射体产生的光子分布里循环
    double Xr, Yr;
    for (int irad = 0; irad < det->nRadLayer; irad++)
    {
        double muPhoton = det->fHitMapEachRad[irad]->Integral();
        if (muPhoton == 0)
            continue;

        cout << "-->" << det->particle << " " << det->momentum << " " << det->theta0 << " rad:" << irad << " ph=" << muPhoton << " realTh=" << det->fThetaCReal[irad] << endl;

        thklist.clear();
        reflist.clear();
        thklist = GetThickList(det, irad);
        reflist = GetRefIndList(det, irad, hypoLambda);
        abslist = GetAbsLenList(det, irad, hypoLambda);

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
                    det->fHitMapEachRad[irad]->GetRandom2(Xr, Yr);
                    Xr = int(Xr / det->pixel) * det->pixel + det->pixel / 2;
                    Yr = int(Yr / det->pixel) * det->pixel + det->pixel / 2;
                    double theta = ReconstructRICHByBeta(det, irad, Xr, Yr, thklist, reflist, abslist);
                    avgth += theta;
                }
                avgth /= iph;

                det->fRecMap->Fill(iph, avgth);
                det->fRecMapEachRad[irad]->Fill(iph, avgth);
            }
        }
    }
}

//---- 被GuiAction调用的函数 -- reconstruct --
/// 生成momentum/theta范围的重建中心值偏差及展宽的分布图, 必须先完成GenerateTheScanHitMapsForEachDetector
void *ReconstructionHandle(void *ptr)
{
    MyRICHDetector *gDet = gMyCommonRICH->GetDetector();

    long ip = (long)ptr;
    int ibegin, iend;
    GetMomentumScanRange(ip, ibegin, iend);

    if (ibegin == -1 || iend == -1)
        return 0;

    for (int imom = ibegin; imom < iend; imom++)
        for (int ithe = 0; ithe < (int)gDet->nthe0; ithe++)
            for (int ihypo = 0; ihypo < (int)gDet->nhypo; ihypo++)
                gMyCommonRICH->ReconstructRICHDetector(imom, ithe, ihypo);

    return 0;
}

bool MyCommonRICH::ReconstructForEachDetector()
{
    if (GetDetScanNumber() == 0)
    {
        cout << "--> Must generate the hitmap first before performing the reconstruction." << endl;
        return false;
    }

    int nthread = (NThread < gDet->np) ? NThread : gDet->np;

    TThread *thread[nthread];
    cout << "----> Applying " << nthread << " thread to perform reconstruction." << endl;
    for (int i = 0; i < nthread; i++)
        thread[i] = new TThread(Form("trec%d", i), ReconstructionHandle, (void *)(size_t)i);

    for (int i = 0; i < nthread; i++)
        thread[i]->Run();

    for (int i = 0; i < nthread; i++)
        thread[i]->Join();

    cout << "----> Recontruction finished." << endl;
    return true;
}

//---- 被GuiAction调用的函数, 生成相应的histogram
void MyCommonRICH::GenerateRecOffsetSigmaMap()
{
    //
    cout << "-->Generating the offset & sigma map." << endl;
    ResizeRecMap(gDet->nhypo, gDet->nRadLayer, gDet->np, gDet->nthe0, gDet->NPhoton);

    // 1. 拟合重建结果
    int p0 = (Epoch == -1) ? 0 : Epoch;
    int p1 = (Epoch == -1) ? gDet->np : Epoch + 1;
    if (p0 > gDet->np || p1 > gDet->np)
        return;

    for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
    {
        for (int irad = 0; irad < gDet->nRadLayer; irad++)
        {
            for (int imom = p0; imom < p1; imom++)
                for (int ithe = 0; ithe < gDet->nthe0; ithe++)
                    for (int iph = 1; iph < gDet->NPhoton; iph++)
                    {
                        TH1D *fproj = GenRecMap(ihypo, irad, imom, ithe, iph);
                        double mean = (fproj->Integral() == 0) ? 0 : fproj->GetFunction("gaus")->GetParameter(1);
                        double sigm = (fproj->Integral() == 0) ? 0 : fproj->GetFunction("gaus")->GetParameter(2);
                        double merr = (fproj->Integral() == 0) ? 0 : fproj->GetFunction("gaus")->GetParError(1);
                        double serr = (fproj->Integral() == 0) ? 0 : fproj->GetFunction("gaus")->GetParError(2);
                        if (fproj->Integral() != 0 && fproj->GetFunction("gaus")->GetNDF() == 0)
                        {
                            mean = (fproj->Integral() == 0) ? 0 : fproj->GetMean();
                            sigm = (fproj->Integral() == 0) ? 0 : fproj->GetRMS();
                            merr = (fproj->Integral() == 0) ? 0 : fproj->GetMeanError();
                            serr = (fproj->Integral() == 0) ? 0 : fproj->GetRMSError();
                        }
                        fRecOffList[ihypo][irad][imom][ithe][iph] = mean;
                        fRecSigList[ihypo][irad][imom][ithe][iph] = sigm;
                        fRecOffErrList[ihypo][irad][imom][ithe][iph] = merr;
                        fRecSigErrList[ihypo][irad][imom][ithe][iph] = serr;
                        //fRecMap[ihypo][irad][imom][ithe][iph] = fproj;
                    }
        }
    }

    for (int imom = p0; imom < p1; imom++)
    {
        gMyCommonRICH->SaveHitFile(imom);
        SaveRecFile(imom);
    }
    cout << "-->the offset & sigma maps are generated." << endl;
}

//---- 被GuiAction调用的函数, 生成相应的histogram
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
    fOffsetMap->SetXTitle("Momentum(GeV/c)");
    fOffsetMap->SetYTitle("#theta_0(degree)");

    fSigmaMap = new TH2F("_sigma_", Form("sigma map for %s from radiator %s with %d photons", particle.Data(), rad.Data(), iph), gDet->np, gDet->pMin, gDet->pMax, gDet->nthe0, gDet->The0Min, gDet->The0Max);
    fSigmaMap->SetXTitle("Momentum(GeV/c)");
    fSigmaMap->SetYTitle("#theta_0(degree)");
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
    for (int jmom = 0; jmom < gDet->np; jmom++)
        for (int jthe = 0; jthe < gDet->nthe0; jthe++)
        {
            fOffsetMap->SetBinContent(jmom + 1, jthe + 1, fRecOffList[ihypo][irad][jmom][jthe][iph]);
            fOffsetMap->SetBinError(jmom + 1, jthe + 1, fRecOffErrList[ihypo][irad][jmom][jthe][iph]);
            fSigmaMap->SetBinContent(jmom + 1, jthe + 1, fRecSigList[ihypo][irad][jmom][jthe][iph]);
            fSigmaMap->SetBinError(jmom + 1, jthe + 1, fRecSigErrList[ihypo][irad][jmom][jthe][iph]);
        }

    for (int iihypo = 0; iihypo < gDet->nhypo; iihypo++)
    {
        for (int iiph = 0; iiph < gDet->NPhoton; iiph++)
        {
            fOffsetVsNphPlot[iihypo]->SetBinContent(iiph + 1, fRecOffList[iihypo][irad][imom][ithe][iiph]);
            fOffsetVsNphPlot[iihypo]->SetBinError(iiph + 1, fRecOffErrList[iihypo][irad][imom][ithe][iiph]);
            fSigmaVsNphPlot[iihypo]->SetBinContent(iiph + 1, fRecSigList[iihypo][irad][imom][ithe][iiph]);
            fSigmaVsNphPlot[iihypo]->SetBinError(iiph + 1, fRecSigErrList[iihypo][irad][imom][ithe][iiph]);
        }
        fNphFcn->SetParLimits(0, 0, 2 * fRecSigList[iihypo][irad][imom][ithe][1]);
        fNphFcn->SetParLimits(1, 0, 2 * fRecSigList[iihypo][irad][imom][ithe][gDet->NPhoton - 2]);
        fSigmaVsNphPlot[iihypo]->Fit("fNphFcn");
        for (int iimom = 0; iimom < gDet->np; iimom++)
        {
            fOffsetVsMomPlot[iihypo]->SetBinContent(iimom + 1, fRecOffList[iihypo][irad][iimom][ithe][iph]);
            fOffsetVsMomPlot[iihypo]->SetBinError(iimom + 1, fRecOffErrList[iihypo][irad][iimom][ithe][iph]);
            fSigmaVsMomPlot[iihypo]->SetBinContent(iimom + 1, fRecSigList[iihypo][irad][iimom][ithe][iph]);
            fSigmaVsMomPlot[iihypo]->SetBinError(iimom + 1, fRecSigErrList[iihypo][irad][iimom][ithe][iph]);
        }
        for (int iithe = 0; iithe < gDet->nthe0; iithe++)
        {
            fOffsetVsThetaPlot[iihypo]->SetBinContent(iithe + 1, fRecOffList[iihypo][irad][imom][iithe][iph]);
            fOffsetVsThetaPlot[iihypo]->SetBinError(iithe + 1, fRecOffErrList[iihypo][irad][imom][iithe][iph]);
            fSigmaVsThetaPlot[iihypo]->SetBinContent(iithe + 1, fRecSigList[iihypo][irad][imom][iithe][iph]);
            fSigmaVsThetaPlot[iihypo]->SetBinError(iithe + 1, fRecSigErrList[iihypo][irad][imom][iithe][iph]);
        }
    }
}

//______________________________________________________________________________
// 计算在当前det的hypothesis下的prob大小
double MyCommonRICH::CalPIDProb(MyRICHDetector *det, vector<pair<double, double>> hit, double &chi2, double &ndf)
{
    int ihypo = GetHypoID(det->particle);
    int imom = GetMomID(det->momentum);
    int ithe = GetThetaID(det->Theta0);

    double hypoLambda = 185;
    UpdateThetaCExp(det, hypoLambda);

    vector<vector<double>> recList; //[nrad][第n个光子重建的切伦科夫角]
    recList.resize(det->nRadLayer);

    for (int i = 0; i < (int)hit.size(); i++)
    {
        int fromRad = -1;    //光从哪个辐射体发出的
        double recThe = 999; //重建的切伦科夫角
        double recDis = 999; //重建的切伦科夫角和真值的差

        for (int irad = 0; irad < det->nRadLayer; irad++)
        {
            vector<double> thklist = GetThickList(det, irad);
            vector<double> reflist = GetRefIndList(det, irad, hypoLambda); //折射率选用假设的波长处的折射率
            vector<double> abslist = GetAbsLenList(det, irad, hypoLambda);

            double rectmp = ReconstructRICHByBeta(det, irad, hit[i].first, hit[i].second, thklist, reflist, abslist);
            //double rectmp = ReconstructRICHBySolver(det, irad, hit[i].first, hit[i].second);
            if (fabs(rectmp - det->fThetaCReal[irad]) < recDis)
            {
                fromRad = irad;
                recThe = rectmp;
                recDis = fabs(rectmp - det->fThetaCReal[irad]);
            }
        }

        if (fromRad == -1)
            continue;
        //if (fRecSigList[ihypo][fromRad][imom][ithe][1] == 0) continue;
        //recList[fromRad].push_back(pow(recThe - fRecOffList[ihypo][fromRad][imom][ithe][1], 2) / pow(fRecSigList[ihypo][fromRad][imom][ithe][1], 2));
        recList[fromRad].push_back(recThe);
    }

    ndf = 0;
    chi2 = 0;
    for (int irad = 0; irad < det->nRadLayer; irad++)
    {
        double avg = 0;
        int nph = recList[irad].size();
        if (nph <= 0 || fRecSigList[ihypo][irad][imom][ithe][nph] == 0)
            continue;
        ndf++;

        for (int i = 0; i < nph; i++)
            avg += recList[irad][i];
        avg /= nph;

        chi2 += pow((avg - fRecOffList[ihypo][irad][imom][ithe][nph]) / fRecSigList[ihypo][irad][imom][ithe][nph], 2);
    }

    return TMath::Prob(chi2, ndf);
}

//根据击中位置坐标判断粒子种类
int MyCommonRICH::PIDProb(vector<MyRICHDetector *> detlist, vector<pair<double, double>> hit, vector<double> chilist, vector<double> ndflist)
{
    int iprob = -1;
    double mprob = -1;
    double chi2, ndf;


    //计算那种hypothesis具有最大的prob
    for (int ihypo = 0; ihypo < (int)detlist.size(); ihypo++)
    {
        double prob = CalPIDProb(detlist[ihypo], hit, chi2, ndf);
        if (prob > mprob)
        {
            mprob = prob;
            iprob = ihypo;
        }
        chilist[ihypo] = chi2;
        ndflist[ihypo] = ndf;
    }
    return iprob;
}


void MyCommonRICH::CalculatePIDForDetector(MyRICHDetector *det)
{
    if (det == 0)
        det = gDet;

    int ihypo = GetHypoID(det->particle);
    int imom = GetMomID(det->momentum);
    int ithe = GetThetaID(det->Theta0);

    if (det->fHitMap == 0)
        return;

    TH2F *fhm = det->GetHitMap();
    double muPhoton = fhm->Integral();
    if (muPhoton < 1)
    {
        for (int jhypo = 0; jhypo < gDet->nhypo; jhypo++)
            gMyCommonRICH->SetPIDEff(imom, ithe, ihypo, jhypo, 0);
        return;
    }

    vector<double> pideff;
    vector<double> chilist;
    vector<double> ndflist;
    pideff.resize(gDet->nhypo);
    chilist.resize(gDet->nhypo);
    ndflist.resize(gDet->nhypo);
    for (int jhypo = 0; jhypo < gDet->nhypo; jhypo++)
        pideff[jhypo] = 0;

    vector<pair<double, double>> hit;
    int nEvent = gMyCommonRICH->GetNevent();

    for (int i = 0; i < nEvent; i++)
    {
        double nPh = gRandom->Poisson(muPhoton);

        hit.clear();
        for (int iph = 0; iph < nPh; iph++)
        {
            double Xr, Yr;
            fhm->GetRandom2(Xr, Yr);
            Xr = int(Xr / det->pixel) * det->pixel + det->pixel / 2;
            Yr = int(Yr / det->pixel) * det->pixel + det->pixel / 2;
            hit.push_back(make_pair(Xr, Yr));
        }

        int pid = gMyCommonRICH->PIDProb(gScanDetList[imom][ithe], hit, chilist, ndflist);
        gMyCommonRICH->FillChiSquare(imom, ithe, ihypo, chilist, ndflist);
        pideff[pid]++;
    }

    for (int jhypo = 0; jhypo < gDet->nhypo; jhypo++)
        gMyCommonRICH->SetPIDEff(imom, ithe, ihypo, jhypo, pideff[jhypo] / nEvent);

    return;
}

//计算PID效率和误判率
void *CalPIDEfficiencyHandle(void *ptr)
{
    MyRICHDetector *gDet = gMyCommonRICH->GetDetector();

    long ip = (long)ptr;
    int ibegin, iend;
    GetMomentumScanRange(ip, ibegin, iend);

    if (ibegin == -1 || iend == -1)
        return 0;

    for (int imom = ibegin; imom < iend; imom++)
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
            for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
                gMyCommonRICH->CalculatePIDForDetector(imom, ithe, ihypo);

    for (int imom = ibegin; imom < iend; imom++)
        gMyCommonRICH->SavePidFile(imom);

    return 0;
}

bool MyCommonRICH::CalPIDEfficiency()
{
    //if (fRecSigList.size() == 0)
    //    return false;

    ResizePIDEffMap(gDet->np, gDet->nthe0, gDet->nhypo);

    int nthread = (NThread < gDet->np) ? NThread : gDet->np;

    TThread *thread[nthread];
    cout << "----> Applying " << nthread << " thread to perform PID calculation." << endl;
    for (int i = 0; i < nthread; i++)
        thread[i] = new TThread(Form("tpid%d", i), CalPIDEfficiencyHandle, (void *)(size_t)i);

    for (int i = 0; i < nthread; i++)
        thread[i]->Run();

    for (int i = 0; i < nthread; i++)
        thread[i]->Join();

    cout << "----> PID efficiency calculation finished." << endl;
    return true;
}

void MyCommonRICH::GeneratePIDHistograms(TString particle, int imom, int ithe)
{
    //1. clear
    for (int i = 0; i < (int)fPIDMap.size(); i++)
        delete fPIDMap[i];
    for (int i = 0; i < (int)fPIDVsMomPlot.size(); i++)
        delete fPIDVsMomPlot[i];
    for (int i = 0; i < (int)fPIDVsThetaPlot.size(); i++)
        delete fPIDVsThetaPlot[i];

    fPIDMap.resize(gDet->nhypo);
    fPIDVsMomPlot.resize(gDet->nhypo);
    fPIDVsThetaPlot.resize(gDet->nhypo);

    //2. define histograms
    int ihypo = GetHypoID(particle);
    double mom = gScanDetList[imom][ithe][0]->momentum;
    double the = gScanDetList[imom][ithe][0]->Theta0;

    for (int i = 0; i < gDet->nhypo; i++)
    {
        fPIDMap[i] = new TH2F(Form("_pid%d_", i), Form("PID efficiency map for %s identified as %s", particle.Data(), SHYPO[i]), gDet->np, gDet->pMin, gDet->pMax, gDet->nthe0, gDet->The0Min, gDet->The0Max);
        fPIDMap[i]->SetXTitle("Momentum(GeV/c)");
        fPIDMap[i]->SetYTitle("#theta_0(degree)");

        fPIDVsMomPlot[i] = new TH1F(Form("_pidvmom_%d", i), Form("PID efficiency vs. Momentum for %s[%.1f#circ] identified as %s", particle.Data(), the, SHYPO[i]), gDet->np, gDet->pMin, gDet->pMax);
        fPIDVsMomPlot[i]->SetXTitle("Momentum(GeV/c)");
        fPIDVsMomPlot[i]->SetYTitle("PID Efficiency(%)");

        fPIDVsThetaPlot[i] = new TH1F(Form("_pidvthe_%d", i), Form("PID efficiency vs. Momentum for %s[%.1fGeV/c] identified as %s", particle.Data(), mom, SHYPO[i]), gDet->nthe0, gDet->The0Min, gDet->The0Max);
        fPIDVsThetaPlot[i]->SetXTitle("#theta_0(degree)");
        fPIDVsThetaPlot[i]->SetYTitle("PID Efficiency(%)");
    }

    //3. fill histograms
    for (int jmom = 0; jmom < gDet->np; jmom++)
        for (int jthe = 0; jthe < gDet->nthe0; jthe++)
            for (int jhypo = 0; jhypo < gDet->nhypo; jhypo++)
                fPIDMap[jhypo]->SetBinContent(jmom + 1, jthe + 1, fPidEffList[jmom][jthe][ihypo][jhypo]);

    for (int jhypo = 0; jhypo < gDet->nhypo; jhypo++)
    {

        for (int iimom = 0; iimom < gDet->np; iimom++)
        {
            fPIDVsMomPlot[jhypo]->SetBinContent(iimom + 1, fPidEffList[iimom][ithe][ihypo][jhypo]);
            fPIDVsMomPlot[jhypo]->SetBinError(iimom + 1, 1. / sqrt(nEvent));
        }
        for (int iithe = 0; iithe < gDet->nthe0; iithe++)
        {
            fPIDVsThetaPlot[jhypo]->SetBinContent(iithe + 1, fPidEffList[imom][iithe][ihypo][jhypo]);
            fPIDVsThetaPlot[jhypo]->SetBinError(iithe + 1, 1. / sqrt(nEvent));
        }
    }
}

//______________________________________________________________________________
/// FILE IO functions
int MyCommonRICH::LoadRings(const char *fname)
{
    // return: -1 文件不存在
    //          1 gdet文件存在
    //          2 hitmap文件存在
    //          3 recmap文件存在
    //          4 pidmap文件存在

    fileName = TString(fname);
    if (fileName.Index(".root") == -1)
        return -1;

    headName = fileName;
    headName.Remove(headName.Index(".root"), 5);

    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    if (LoadDetFile() < 0)
        return -1;

    if (LoadHitFile() < 0)
        return 1;

    if (LoadRecFile() < 0)
        return 2;

    if (LoadPidFile() < 0)
        return 3;

    cout << "----> All files are loaded. " << endl;
    return 4;
}

void MyCommonRICH::SaveRings(const char *fname)
{
    //hitmap文件规则：
    //随着分析数据范围的变多以及分析时间过长，将文件划分为多个子文件将有助于增加鲁棒性
    //hitmap文件作为分析的第一步，按以下规则来保存：
    //1. 在用户点击SaveRings，仅仅保存用户设置的探测器参数，即gDet，同时保存这个文件名
    //2. 在用户点击Scan Hitmap后，将进行多重循环来生成hitmap
    //   此时以momentum为单位，在文件名后面添加"_mom_%d.root"来保存每一轮的hitmap
    //3. 在用户点击Scan Recmap后，会全部分析完然后保存成"_rec.root"的文件
    //4. 在用户点击Scan PIDmap后，会全部分析完然后保存成"_pid.root"的文件

    fileName = TString(fname);
    headName = fileName;
    headName.Remove(headName.Index(".root"), 5);

    if (Epoch == -1 || Epoch == 0)
        SaveDetFile();
}

int MyCommonRICH::LoadDetFile()
{
    //-----------------------------------
    // 1. 读取基本gDet的root文件
    //
    TFile f(fileName);
    if (!f.IsOpen())
        return -1;

    MyRICHDetector *det = new MyRICHDetector(0);
    TTree *T = (TTree *)f.Get("TTree");
    T->SetBranchAddress("det", &det);
    cout << "----> Total entries: " << T->GetEntries() << endl;

    //1.1 gDet
    int ientry = 0;
    T->GetEntry(ientry++);
    gDirectory->cd();
    gDet->BuildFrom(*det, -1, 1);
    gDet->SetDirectory();

    //1.2 gDetList
    ResizeDetList(gDet->nhypo);
    for (int i = 0; i < gDet->nhypo; i++)
    {
        T->GetEntry(ientry++);
        gDetList[i] = new MyRICHDetector(*det, -1, 1);
        gDetList[i]->SetDirectory();
    }
    f.Close();

    return 0;
}

void MyCommonRICH::SaveDetFile()
{
    //1. fileName保存用户设置的det信息
    TFile f(fileName, "recreate");
    if (!f.IsOpen())
        return;

    // save detector infor.
    TTree *T = new TTree("TTree", "NumClass");
    MyRICHDetector *det = new MyRICHDetector();
    T->Branch("det", "MyRICHDetector", &det, 20000, 1);

    //1. 保存gDet
    det = gDet; //->BuildFrom(*gDet, -1, 1);
    T->Fill();

    //2. gDetList
    for (int i = 0; i < (int)gDetList.size(); i++)
    {
        det = gDetList[i]; //->BuildFrom(*gDetList[i]);
        T->Fill();
    }

    T->Write();
    f.Close();
}

int MyCommonRICH::LoadHitFile()
{
    //-----------------------------------
    // 2. 读取gScanDetList，也就是hitmap文件
    //
    ResizeScanDetList(gDet->np, gDet->nthe0, gDet->nhypo, 0);

    int ibegin, iend;
    GetMomentumScanRange(0, ibegin, iend);

    for (int imom = ibegin; imom < iend; imom++)
    {
        TString hitFile = headName;
        hitFile += TString(Form("_imom_%d.root", imom));
        TFile f(hitFile);

        if (!f.IsOpen())
        {
            cout << "##### Fatal error: Can't open " << hitFile << " to read!" << endl;
            return -1;
        }

        MyRICHDetector *detector = new MyRICHDetector(0);
        TTree *T = (TTree *)f.Get("TTree");
        T->SetBranchAddress("det", &detector);
        int ientry = 0;

        cout << "----> Loading imom=" << imom << " from " << hitFile << endl;
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
            for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            {
                T->GetEntry(ientry++);
                gScanDetList[imom][ithe][ihypo] = new MyRICHDetector(*detector, -1, 1);
                gScanDetList[imom][ithe][ihypo]->SetDirectory();
            }

        f.Close();
    }

    return 0;
}

void MyCommonRICH::SaveHitFile(int imom)
{
    TString hitFile = headName;
    hitFile += TString(Form("_imom_%d.root", imom));
    TFile f(hitFile, "recreate");

    // save detector infor.
    TTree *T = new TTree("TTree", "NumClass");
    MyRICHDetector *det = new MyRICHDetector();
    T->Branch("det", "MyRICHDetector", &det, 20000, 1);

    for (int ithe = 0; ithe < gDet->nthe0; ithe++)
        for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
        {
            det = gScanDetList[imom][ithe][ihypo];
            T->Fill();
        }

    T->Write();
    f.Close();
}

int MyCommonRICH::LoadRecFile()
{
    //------------------
    // 3. 读取recmap文件
    //
    ResizeRecMap(gDet->nhypo, gDet->nRadLayer, gDet->np, gDet->nthe0, gDet->NPhoton);

    int ibegin, iend;
    GetMomentumScanRange(0, ibegin, iend);

    for (int imom = ibegin; imom < iend; imom++)
    {
        TString recFile = headName;
        recFile += TString(Form("_rec_%d.root", imom));

        TFile f(recFile);
        if (!f.IsOpen())
            return -1;

        TTree *T1 = (TTree *)f.Get("Configure");
        TTree *T2 = (TTree *)f.Get("RecTree");

        if (T1 == NULL)
            return -1;

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

        cout << "----> Loading REC file: nHypo=" << nhyp << ", nRadiator=" << nrad << ", nMomentum=" << nmom << ", nTheta0=" << nthe << ", nPhoton=" << nph << endl;
        if ((nhyp != gDet->nhypo) || (nrad != gDet->nRadLayer) || (nmom != gDet->np) ||
            (nthe != gDet->nthe0) || (nph != gDet->NPhoton))
        {
            cout << "##### This map doesn't match with the hitmap, please check this data file is the right one or not." << endl;
            return -1;
        }

        if (pMin != gDet->pMin || pMax != gDet->pMax || The0Min != gDet->The0Min || The0Max != gDet->The0Max)
        {
            cout << "##### This map doesn't math with the global detector settings(pMin/pMax/the0Min/the0Max), please check this data file is the right one or not." << endl;
            cout << "##### Read values are momentum: (" << pMin << ", " << pMax << "), theta: (" << The0Min << ", " << The0Max << ")" << endl;
            cout << "##### gDet values are momentum: (" << gDet->pMin << ", " << gDet->pMax << "), theta: (" << gDet->The0Min << ", " << gDet->The0Max << ")" << endl;
            return -1;
        }

        if (T2 == NULL)
            return -1;

        double mean, merr;
        double sigm, serr;
        T2->SetBranchAddress("mean", &mean);
        T2->SetBranchAddress("merr", &merr);
        T2->SetBranchAddress("sigm", &sigm);
        T2->SetBranchAddress("serr", &serr);

        int ientry = 0;
        for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            for (int irad = 0; irad < gDet->nRadLayer; irad++)
                for (int ithe = 0; ithe < gDet->nthe0; ithe++)
                    for (int iph = 0; iph < gDet->NPhoton; iph++)
                    {
                        T2->GetEntry(ientry++);
                        fRecOffList[ihypo][irad][imom][ithe][iph] = mean;
                        fRecOffErrList[ihypo][irad][imom][ithe][iph] = merr;
                        fRecSigList[ihypo][irad][imom][ithe][iph] = sigm;
                        fRecSigErrList[ihypo][irad][imom][ithe][iph] = serr;
                    }
        f.Close();
    }

    return 0;
}

void MyCommonRICH::SaveRecFile(int imom)
{
    TString recFile = headName;
    recFile += TString(Form("_rec_%d.root", imom));
    TFile f(recFile, "recreate");
    if (!f.IsOpen())
        return;

    TTree *T1 = new TTree("Configure", "Detector Info");
    TTree *T2 = new TTree("RecTree", "Mean/Sigma Map");

    //Save Tree1
    int nhyp, nrad, nmom, nthe, nph;
    double pMin, pMax, The0Min, The0Max;
    T1->Branch("nhyp", &nhyp, "nhyp/I");
    T1->Branch("nrad", &nrad, "nrad/I");
    T1->Branch("nmom", &nmom, "nmom/I");
    T1->Branch("nthe", &nthe, "nthe/I");
    T1->Branch("nph", &nph, "nph/I");
    T1->Branch("pMin", &pMin, "pMin/D");
    T1->Branch("pMax", &pMax, "pMax/D");
    T1->Branch("The0Min", &The0Min, "The0Min/D");
    T1->Branch("The0Max", &The0Max, "The0Max/D");

    nhyp = gDet->nhypo;
    nrad = gDet->nRadLayer;
    nmom = gDet->np;
    nthe = gDet->nthe0;
    nph = gDet->NPhoton;
    pMin = gDet->pMin;
    pMax = gDet->pMax;
    The0Min = gDet->The0Min;
    The0Max = gDet->The0Max;

    T1->Fill();
    T1->Write();

    //Save Tree2
    double mean, merr;
    double sigm, serr;
    T2->Branch("mean", &mean, "mean/D");
    T2->Branch("merr", &merr, "merr/D");
    T2->Branch("sigm", &sigm, "sigm/D");
    T2->Branch("serr", &serr, "serr/D");

    // save detector infor.
    int ientry = 0;
    for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
        for (int irad = 0; irad < gDet->nRadLayer; irad++)
            for (int ithe = 0; ithe < gDet->nthe0; ithe++)
                for (int iph = 0; iph < gDet->NPhoton; iph++)
                {
                    mean = fRecOffList[ihypo][irad][imom][ithe][iph];
                    merr = fRecOffErrList[ihypo][irad][imom][ithe][iph];
                    sigm = fRecSigList[ihypo][irad][imom][ithe][iph];
                    serr = fRecSigErrList[ihypo][irad][imom][ithe][iph];
                    T2->Fill();
                    ientry++;
                }

    T2->Write();
    f.Close();
}

int MyCommonRICH::LoadPidFile()
{
    //------------------
    //4. 读取pidmap文件
    ResizePIDEffMap(gDet->np, gDet->nthe0, gDet->nhypo);

    for (int imom = 0; imom < gDet->np; imom++)
    {
        TString pidFile = headName;
        pidFile += TString(Form("_pid_%d.root", imom));

        TFile f(pidFile);
        if (!f.IsOpen())
            return -1;

        TTree *T1 = (TTree *)f.Get("Configure");
        TTree *T2 = (TTree *)f.Get("PidTree");

        if (T1 == NULL)
            return -1;

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

        cout << "---->Reading PID root: nHypo=" << nhyp << ", nRadiator=" << nrad << ", nMomentum=" << nmom << ", nTheta0=" << nthe << ", nPhoton=" << nph << endl;
        if ((nhyp != gDet->nhypo) || (nrad != gDet->nRadLayer) || (nmom != gDet->np) ||
            (nthe != gDet->nthe0) || (nph != gDet->NPhoton))
        {
            cout << "##### This map doesn't match with the hitmap, please check this data file is the right one or not." << endl;
            return -1;
        }

        if (pMin != gDet->pMin || pMax != gDet->pMax || The0Min != gDet->The0Min || The0Max != gDet->The0Max)
        {
            cout << "##### This map doesn't math with the global detector settings(pMin/pMax/the0Min/the0Max), please check this data file is the right one or not." << endl;
            cout << "##### Read values are momentum: (" << pMin << ", " << pMax << "), theta: (" << The0Min << ", " << The0Max << ")" << endl;
            cout << "##### gDet values are momentum: (" << gDet->pMin << ", " << gDet->pMax << "), theta: (" << gDet->The0Min << ", " << gDet->The0Max << ")" << endl;
            return -1;
        }

        if (T2 == NULL)
            return -1;

        double pideff;
        T2->SetBranchAddress("pideff", &pideff);

        int ientry = 0;
        for (int ithe = 0; ithe < gDet->nthe0; ithe++)
            for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
                for (int jhypo = 0; jhypo < gDet->nhypo; jhypo++)
                {
                    T2->GetEntry(ientry++);
                    fPidEffList[imom][ithe][ihypo][jhypo] = pideff;
                }
        f.Close();
    }
    return 0;
}

void MyCommonRICH::SavePidFile(int imom)
{
    TString pidFile = headName;
    pidFile += TString(Form("_pid_%d.root", imom));
    TFile f(pidFile, "recreate");
    if (!f.IsOpen())
        return;
    if (fPidEffList.size() == 0)
        return;

    TTree *T1 = new TTree("Configure", "Detector Info");
    TTree *T2 = new TTree("PidTree", "PID Efficiency");

    //Save Tree1
    int nhyp, nrad, nmom, nthe, nph;
    double pMin, pMax, The0Min, The0Max;
    T1->Branch("nhyp", &nhyp, "nhyp/I");
    T1->Branch("nrad", &nrad, "nrad/I");
    T1->Branch("nmom", &nmom, "nmom/I");
    T1->Branch("nthe", &nthe, "nthe/I");
    T1->Branch("nph", &nph, "nph/I");
    T1->Branch("pMin", &pMin, "pMin/D");
    T1->Branch("pMax", &pMax, "pMax/D");
    T1->Branch("The0Min", &The0Min, "The0Min/D");
    T1->Branch("The0Max", &The0Max, "The0Max/D");

    nhyp = gDet->nhypo;
    nrad = gDet->nRadLayer;
    nmom = gDet->np;
    nthe = gDet->nthe0;
    nph = gDet->NPhoton;
    pMin = gDet->pMin;
    pMax = gDet->pMax;
    The0Min = gDet->The0Min;
    The0Max = gDet->The0Max;

    T1->Fill();
    T1->Write();

    //Save Tree2
    double pideff;
    T2->Branch("pideff", &pideff, "pideff/D");
    for (int ithe = 0; ithe < gDet->nthe0; ithe++)
        for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            for (int jhypo = 0; jhypo < gDet->nhypo; jhypo++)
            {
                pideff = fPidEffList[imom][ithe][ihypo][jhypo];
                T2->Fill();
            }
    T2->Write();

    f.cd();
    for (int ithe = 0; ithe < gDet->nthe0; ithe++)
        for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            for (int jhypo = 0; jhypo < gDet->nhypo; jhypo++)
            {
                fPidChiHist[imom][ithe][ihypo][jhypo][0]->Write();
                fPidChiHist[imom][ithe][ihypo][jhypo][1]->Write();
            }

    f.Close();
}

//______________________________________________________________________________
//
void MyCommonRICH::DrawFCN(int irad)
{
    /* 
    double beta = Beta(gDet->momentum, gDet->mass);
    double tbase = gDet->CalZ0(irad);
    double theta0 = gDet->theta0;

    // 1. 用beta方法重建光子
    int ip = 0;
    TGraph2D *g2 = new TGraph2D();
    for (int iphi = 0; iphi < 360; iphi++)
        for (int ilambda = 160; ilambda < 220; ilambda++)
        {
            double phi = iphi * TMath::Pi() / 180;
            double lambda = ilambda;

            // 生成击中
            thklist = GetThickList(gDet, irad);
            reflist = GetRefIndList(gDet, irad, lambda);
            abslist = GetAbsLenList(gDet, irad, lambda);
            double z0 = 0.5 * gDet->tRadLayer[irad];
            double n0 = reflist[0];
            double thetac = CherenkovAng(beta, n0);
            fRhoFcn->SetParameters(z0, theta0, thetac);

            double Z0 = -1 * tbase - z0; //全局坐标系
            double X0 = 0;
            double Y0 = Z0 * tan(theta0);

            double rho = fRhoFcn->Eval(phi);
            if (rho == 0 || rho != rho)
                continue;

            double Xr = rho * cos(phi) + X0;
            double Yr = rho * sin(phi) + Y0;

            // 重建
            double hypoLambda = 185;

            thklist = GetThickList(gDet, irad);
            reflist = GetRefIndList(gDet, irad, hypoLambda); //折射率选用假设的波长处的折射率
            abslist = GetAbsLenList(gDet, irad, hypoLambda);
            UpdateThetaCExp(gDet, hypoLambda);

            //double phi2 = FindPhi(gDet, Xr, Yr, X0, Y0);
            double recThe1 = ReconstructRICHByBeta(gDet, irad, Xr, Yr);
            g2->SetPoint(ip++, iphi, lambda, recThe1 - thetac);
            //if (fabs(recThe1 - thetac) > 1)
            {
                //cout << phi << " " << phi2 << endl;
                //cout << "光子出射的参数：z0=" << z0 << " phi=" << phi << " theta0=" << theta0 << " thetac=" << thetac << endl;
                //cout << "             (X0,Y0,Z0)=(" << X0 << ", " << Y0 << "," << Z0 << ")" << endl;
                //cout << "             (Xr,Yr   )=(" << Xr << ", " << Yr << " ), rho=" << rho << endl;
                //cout << "光子重建的参数方法-Beta： thetac2=" << recThe1 << endl;
            }
        }
    g2->Draw("pcol");
    */

    // 2. 用solver方法来重建, 可以看到对于相同的a，其实会有两个解
    /* 
    double beta = Beta(gDet->momentum, gDet->mass);
    double theta0 = gDet->theta0;

    double phi = 30 * TMath::Pi() / 180;
    double lambda = 185;

    // 生成击中
    thklist = GetThickList(gDet, irad);
    reflist = GetRefIndList(gDet, irad, lambda);
    abslist = GetAbsLenList(gDet, irad, lambda);
    double z0 = 0.5 * gDet->tRadLayer[irad];
    double n0 = reflist[0];
    double thetac = CherenkovAng(beta, n0);
    fRhoFcn->SetParameters(z0, theta0, thetac);

    TGraph *g = new TGraph();
    for (int i = 0; i < 90; i++)
    {
        double thetac = TMath::Pi() / 2 / 90 * i;
        fRhoFcn->SetParameter(2, thetac);
        g->SetPoint(i, thetac, fRhoFcn->Eval(phi));
    }
    g->Draw("apl");
    g->GetXaxis()->SetTitle("#theta");
    g->GetYaxis()->SetTitle("S[mm]");
    //double recThe2 = ReconstructRICHBySolver(gDet, irad, Xr, Yr);
    //cout << "光子重建的参数方法-Solver： thetac2=" << recThe2 << " " << epsilon << endl;
    */

    // 3. 用beta+solver方法来重建

    double theta0 = gDet->theta0;
    double beta = Beta(gDet->momentum, gDet->mass);
    double lambda = 185;
    vector<double> thklist = GetThickList(gDet, irad);
    vector<double> reflist = GetRefIndList(gDet, irad, lambda);
    vector<double> abslist = GetAbsLenList(gDet, irad, lambda);
    double z0 = 0.5 * gDet->tRadLayer[irad];
    double n0 = reflist[0];
    double thetacReal = CherenkovAng(beta, n0);
    double phi = 0;
    double a = findA(theta0, thetacReal, phi);
    double theta1 = atan(a);

    cout << "input a=" << a << ", theta1 = " << theta1 << endl;

    TGraph *g = new TGraph();
    for (int i = 0; i < 90; i++)
    {
        double thetac = TMath::Pi() / 2 / 90 * i;
        g->SetPoint(i, thetac, findA(theta0, thetac, phi));
    }
    g->Draw("apl");
    g->GetXaxis()->SetTitle("#theta");
    g->GetYaxis()->SetTitle("S[mm]");

    double calthe = z0 * cos(theta0) + (a * z0) * sin(phi) * sin(theta0);
    calthe /= sqrt(a * a * z0 * z0 + z0 * z0);
    calthe = acos(calthe);

    cout << "output thetac = " << BetaThetaC(theta0, theta1, phi, epsilon) << endl;
    cout << "calculate = " << calthe << endl
         << endl;
}