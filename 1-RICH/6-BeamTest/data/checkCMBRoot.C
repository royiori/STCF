#include <iostream>
#include <vector>
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGSlider.h"
#include "TGButton.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGFileDialog.h"
#include "TApplication.h"
#include "TRootEmbeddedCanvas.h"
#include "TSystem.h"
#include "../inc/MyBeamTestDatStruct.h"
#include "../inc/MyBeamTestDetector.h"
#include "../inc/MyBeamTestRICH.h"
#include "../inc/MyBeamTestTrackAGET.h"
#include "checkCMBRoot.h"

//#include "/usr/local/include/eigen3/Eigen/Eigen"
//#include "/Users/chad/Work/src/source/GeneralBrokenLines/cpp/include/GblTrajectory.h"
//using namespace gbl;
//using namespace Eigen;

using namespace std;
TH2D *fhmap = NULL;

//----------
//常数定义
enum
{
    X = 0,
    Y = 1,
    Z = 2,
    ANGLE = 3,
} XYENUM;
char xyName[3][2] = {"X", "Y", "Z"};

enum
{
    RICHid = 0,
    T02id = 1,
    T03id = 2,
    T06id = 3,
    T04id = 4,
    T05id = 5,
} DetectorID;
char detName[NDET][10] = {"RICH", "T02", "T03", "T06", "T04", "T05"};

//-------------
// histograms
TH1F *ftmp1;
TH2F *ftmp2;
TH1F *fNCluster[NDET][2]; // number of cluster
TH1F *fCluster[NDET][2];  //cluster size
TH2F *fNCluster2D[NDET];
TH1F *fCharge[NDET];
TH1F *fTotCharge[NDET];
TH1F *fPeakingTime[NDET];
TH2F *fQT[NDET];
TH2F *fhitmap[NDET][2];
TH2F *fhitmapReal[NDET][2];
TGraph *fFootmap;
TH2F *fHitmap;
TH1F *fDeviation[NDET][2];
TH1F *fRec;                //单个cluster的重建结果
TH1F *fRec2;               // 一个事例里的cluster重建结果
TH1F *fQphoton, *fTphoton; //重建有效后，确认是光子的信号，然后再填Q/T的分布
TH1F *fTotQphoton;
TH2F *fQTphoton;
TH1F *fClusterPhoton;
TH1F *fNClusterPhoton;
TH1F *fNRec;
TGraph2D *gRec;

//RUN63:
//double Xoff[NDET] = {0, -0.06816, 0.06108, -0.5638, 0, 0};
//double Yoff[NDET] = {0, -0.15200, 0.13520, -1.2460, 0, 0};

//分析用
class MyBeamTestData2
{
public:
    int id = -1;                    //探测器的id号
    int event = -1;                 //触发事例号
    vector<vector<BeamHit>> branch; //RICH击中信息
    RealBeamHit XYcluster[2];       //条读出，分簇后按照重心得到的三维真实坐标的击中信息
};
vector<MyBeamTestData2 *> detector;

class My3DLine
{
public:
    double A, B, C;

public:
    void SetSurface(double theta, pair<double, double> center, double dist)
    {
        A = tan(theta);
        B = -1;
        C = center.second - A * center.first - dist / cos(theta);
    }

    void SetLine(TF1 *func1, TF1 *func2)
    {
        A = func1->GetParameter(1) * func2->GetParameter(0);
        B = func1->GetParameter(0) * func2->GetParameter(1);
        C = func1->GetParameter(0) * func2->GetParameter(0);
    }

    void SetLine(TF1 *func1)
    {
        A = func1->GetParameter(1);
        B = -1;
        C = func1->GetParameter(0);
    }

    double CrossLineZ(My3DLine *l) { return (C - l->C) / (l->A - A); } //两条直线点交点
    double CrossLineY(My3DLine *l) { return (-l->A * C + A * l->C) / (l->A * B - A * l->B); }
    double PedalZ(pair<double, double> point) { return 1 / (A * A + B * B) * (B * B * point.first - A * B * point.second - A * C); } //点到直线点垂足
    double PedalY(pair<double, double> point) { return 1 / (A * A + B * B) * (-A * B * point.first + A * A * point.second - B * C); }
    void Print() { cout << A << ", " << B << ", " << C; }
};

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

    a *= cc / c0 / c0 * sqrt(b) + 2 * sp * t0;
    a /= c;
    return a;
}

double findA2(double theta0, double thetac, double phi)
{
    double cosThetac = cos(thetac);
    double sinTheta0 = sin(theta0);
    double cosTheta0 = cos(theta0);
    double tanTheta0 = tan(theta0);
    double sinPhi = sin(phi);
    double cosPhi = cos(phi);
    double cos2Thetac = cos(2 * thetac);
    double cos2Theta0 = cos(2 * theta0);
    double cos2Phi = cos(2 * phi);

    if (cosTheta0 == 0)
        return 0;

    double secTheta0 = 1 / cosTheta0;

    double a = 1;
    double b = pow(cosThetac, 2) * pow(secTheta0, 4) * (-cos2Thetac + cos2Theta0 * pow(cosPhi, 2) + pow(sinPhi, 2));
    double c = 4 * sinPhi * tanTheta0;
    double d = 4 * pow(cosThetac, 2) * pow(secTheta0, 2) - 4 * pow(sinPhi, 2) * pow(tanTheta0, 2);

    if (b < 0 || d == 0)
        return 0;

    a = (2 * sqrt(2) * sqrt(b) + c) / d;
    return a;
}

double findA3(double theta0, double thetac, double phi)
{
    double cosThetac = cos(thetac);
    double sinTheta0 = sin(theta0);
    double cosTheta0 = cos(theta0);
    double tanTheta0 = tan(theta0);
    double sinPhi = sin(phi);
    double cosPhi = cos(phi);
    double cos2Thetac = cos(2 * thetac);
    double cos2Theta0 = cos(2 * theta0);
    double cos2Phi = cos(2 * phi);

    if (cosTheta0 == 0)
        return 0;

    double secTheta0 = 1 / cosTheta0;

    double a = 1;
    double b = pow(cosThetac, 2) * pow(secTheta0, 4) * (-cos2Thetac + cos2Theta0 * pow(cosPhi, 2) + pow(sinPhi, 2));
    double c = 4 * sinPhi * tanTheta0;
    double d = 4 * pow(cosThetac, 2) * pow(secTheta0, 2) - 4 * pow(sinPhi, 2) * pow(tanTheta0, 2);

    if (b < 0 || d == 0)
        return 0;

    a = (-2 * sqrt(2) * sqrt(b) + c) / d;
    return a;
}

double ReconstructRICHByBeta(double Xr, double Yr)
//在入射点为原点的坐标系下的坐标， 粒子沿着+y和-z方向前进。
//入射角theta为[0,pi/2]范围, Z出发到束流负方向，逆时针为正
//phi为[-pi,pi]范围，从X轴出发到切伦科夫光的投影线的夹角，逆时针为正
{
    //if (Xr < 0 || Xr > 20 || Yr < 60 || Yr > 70)
    //    return 0;

    //Xr=65;
    //Yr=5;
    //1. 初始化
    double nquartz = 1.585;
    double theta0 = angle / 180 * TMath::Pi(); //入射粒子和法线的夹角
    double beta = 1;
    double thetac = acos(1 / nquartz);
    double Tg = 93;       //gas的厚度
    double z0 = 5;        //预期出光点位置到辐射体上表面的平均位置～辐射体厚度/2
    double Xep = 10 - z0; //出光点位置到辐射体下表面的距离

    //2. 求解出光角度phi
    double tbase = 0;
    double Z0 = -1 * tbase - z0; //发光点的坐标
    double X0 = 0;
    double Y0 = fabs(Z0) * tan(theta0);
    double phi = atan((Yr - Y0) / (Xr - X0));
    double R = sqrt((Xr - X0) * (Xr - X0) + (Yr - Y0) * (Yr - Y0));

    //phi取值区间移到[0,2pi]
    if (Yr > 0 && Xr > 0)
        phi = phi;
    if (Yr > 0 && Xr < 0)
        phi = phi + TMath::Pi();
    if (Yr < 0 && Xr < 0)
        phi = phi + TMath::Pi();
    if (Yr < 0 && Xr > 0)
        phi = phi + 2 * TMath::Pi();

    double a[2];
    a[0] = findA2(theta0, thetac, phi); //两个解
    a[1] = findA3(theta0, thetac, phi);

    double rthec[2];

    //cout << "\n(*\n(" << Yr << "-" << Y0 << ") / (" << Xr << ", " << X0 << ")" << endl;
    //cout << "phi = " << phi << "; thetac=" << thetac << "; theta0=" << theta0 << "; Xep=" << Xep << ";" << endl;

    for (int i = 0; i < 2; i++)
    {
        double theta1 = atan(a[i]);

        //3. 求R0
        double Rrad = Xep * tan(theta1);
        double Rqz = 0;
        double nSinTheta = nquartz * sin(theta1);

        //只有quartz一个辐射体，所以不需要计算
        //for (int i = 1; i < (int)thklist.size() - 1; i++)
        //{
        //    if (nSinTheta / reflist[i] > 1)
        //        continue;
        //    Rqz += thklist[i] * tan(asin(nSinTheta / reflist[i]));
        //}

        double R0 = R - Rrad - Rqz;

        //4. 求重建的rec-theta1
        double theta2 = atan(R0 / Tg);
        double nquartz2 = sin(theta2) / sin(theta1);
        rthec[i] = acos(1 / nquartz2 / beta);

        //cout << "a=" << a[i] * Xep << " R0=" << R0 << " Theta2=" << theta2 << " " << sin(theta2) << "/" << sin(theta1) << endl;
        //cout << "n=" << nquartz2 << " rec=" << rthec[i] << endl;
    }

    if (verbose)
    {
        cout << "\n(*\n(" << Yr << "-" << Y0 << ") / (" << Xr << ", " << X0 << ")" << endl;
        cout << "phi = " << phi << "; thetac=" << thetac << "; theta0=" << theta0 << "; Xep=" << Xep << ";" << endl;
        cout << "a1=" << a[0] * Xep << " a2=" << a[1] * Xep << endl;
        cout << "R01=" << R - a[0] * Xep << " R02=" << R - a[1] * Xep << endl;
        //", theta1=" << theta1 << endl;
        //cout << "R0=" << R0 << ", R=" << R << ", Rrad=" << Rrad << endl;
        //cout << "nquartz2 = " << nquartz2 << ", sinTheta2=" << R0 / sqrt(R0 * R0 + Tg * Tg) << ", sinTheta1 = " << sin(theta1) << endl;
        cout << "rthec=" << rthec[0] << ", " << rthec[1] << endl;
        cout << "*)\n\n";
    }

    if (isnan(rthec[0]))
        return rthec[1];
    if (isnan(rthec[1]))
        return rthec[0];
    if (fabs(rthec[0] - thetac) < fabs(rthec[1] - thetac))
        return rthec[0];
    return rthec[1];
}

//---------------------------------------------------------------------------
// 分析Tracker数据
// 计算Tracker在束流实验系下的hit坐标
double CalTrackPosition(int detID, int flag, double channel = 0)
{
    double position = 0;
    if (detID == RICHid)
        return position;

    if (detID == T02id || detID == T03id || detID == T06id || detID==T04id || detID==T05id)
    {
        if (flag == X)
            position = Xoff[detID] + (128 / 2. - channel + 1 / 2.) * 0.4;
        if (flag == Y)
            position = Yoff[detID] + (channel - 128 / 2. + 1 / 2.) * 0.4;
        if (flag == Z)
            position = Zpos[detID];
    }
    return position;
}

// 判断该cluster是否有效，若有效则给出此cluster的mean和rms, control用来控制不同探测器的不同条件
bool IsTrackClusterEffective(int id, vector<BeamHit> cluster, double &mean, double &rms, double &qsum, int control = 0)
{
    mean = 0;
    rms = 0;
    qsum = 0;

    //判断是否是放电信号或噪声信号
    if (id == T02id || id == T03id || id == T06id)
    {
        for (int i = 0; i < (int)cluster.size(); i++)
            if (cluster[i].t < 230 || cluster[i].t > 260 || cluster[i].q > 3600 || cluster[i].q < 10)
                return false;
    }
    else
    {
        for (int i = 0; i < (int)cluster.size(); i++)
            if (cluster[i].q > 3600 || cluster[i].q < 10)
                return false;
    }

    ftmp1->Reset();
    for (int i = 0; i < (int)cluster.size(); i++)
    {
        qsum += cluster[i].q;
        if (cluster[i].hit.first == -999)
            ftmp1->Fill(cluster[i].hit.second, cluster[i].q);
        if (cluster[i].hit.second == -999)
            ftmp1->Fill(cluster[i].hit.first, cluster[i].q);
    }
    mean = ftmp1->GetMean();
    rms = ftmp1->GetRMS();

    return true;
}

void FillTrackHistogram(vector<MyBeamTestData *> detlist, int detID)
{
    for (int i = 0; i < (int)detlist.size(); i++)
    {
        if (detlist[i]->id != detID)
            continue;

        MyBeamTestData *detector = detlist[i];

        // cluster 分布
        for (int xy = 0; xy < 2; xy++)
        {
            fNCluster[detID][xy]->Fill(detector->XYbranch[xy].size());
            for (int j = 0; j < (int)detector->XYbranch[xy].size(); j++)
                fCluster[detID][xy]->Fill(detector->XYbranch[xy][j].size());
        }
        fNCluster2D[detID]->Fill(detector->XYbranch[X].size(), detector->XYbranch[Y].size());

        //电荷Q/达峰时间T分布
        double qsum = 0;
        for (int j = 0; j < (int)detector->charge.size(); j++)
        {
            qsum += detector->charge[j];
            fCharge[detID]->Fill(detector->charge[j]);
        }
        fTotCharge[detID]->Fill(qsum);

        for (int j = 0; j < (int)detector->time.size(); j++)
        {
            fPeakingTime[detID]->Fill(detector->time[j]);
            fQT[detID]->Fill(detector->charge[j], detector->time[j]);
        }

        //fhitmap分布
        double xmean, xrms, xq;
        double ymean, yrms, yq;
        for (int ii = 0; ii < (int)detector->XYbranch[X].size(); ii++)
            for (int jj = 0; jj < (int)detector->XYbranch[Y].size(); jj++)
            {
                for (int iii = 0; iii < (int)detector->XYbranch[X][ii].size(); iii++)
                    for (int jjj = 0; jjj < (int)detector->XYbranch[Y][jj].size(); jjj++)
                        fhitmap[detID][0]->Fill(detector->XYbranch[X][ii][iii].hit.first, detector->XYbranch[Y][jj][jjj].hit.second);

                if (!IsTrackClusterEffective(detID, detector->XYbranch[X][ii], xmean, xrms, xq))
                    continue;
                if (!IsTrackClusterEffective(detID, detector->XYbranch[Y][jj], ymean, yrms, yq))
                    continue;

                for (int iii = 0; iii < (int)detector->XYbranch[X][ii].size(); iii++)
                    for (int jjj = 0; jjj < (int)detector->XYbranch[Y][jj].size(); jjj++)
                    {
                        fhitmap[detID][1]->Fill(detector->XYbranch[X][ii][iii].hit.first, detector->XYbranch[Y][jj][jjj].hit.second);

                        //在束流实验室系下的坐标在这里写出来
                        double x = CalTrackPosition(detID, X, detector->XYbranch[X][ii][iii].hit.first);
                        double y = CalTrackPosition(detID, Y, detector->XYbranch[Y][jj][jjj].hit.second);
                        fhitmapReal[detID][0]->Fill(x, y);
                    }
                double x = CalTrackPosition(detID, X, xmean);
                double y = CalTrackPosition(detID, Y, ymean);
                fhitmapReal[detID][1]->Fill(x, y);
            }
    }
}

void AnalysisTrack(vector<int> uselist, int target, vector<MyBeamTestData *> detlist, int xyflag)
{
    TGraphErrors *g1 = new TGraphErrors(); //方法1: 只要是好的cluster，就使用这个信息
    TGraphErrors *g2 = new TGraphErrors(); //方法2: 只用q最大的那个好的cluster

    int ipoint1 = 0;
    int ipoint2 = 0;
    double mean, rms, qsum;

    //填图
    for (int i = 0; i < (int)uselist.size(); i++)
    {
        //作为use的track必须要有两个以上的击中
        int id = -1;
        int detID = uselist[i];
        for (int j = 0; j < detlist.size(); j++)
            if (detlist[j]->id == detID)
                id = j;
        if (id == -1)
            continue;

        MyBeamTestData *detector = detlist[id];
        vector<vector<BeamHit>> branch = detector->XYbranch[xyflag];

        //作为use的track的number of cluster不能太多
        if (branch.size() >= 2)
            continue;

        //填图以拟合
        double qmax = -1;
        double mmax = -1;
        for (int j = 0; j < (int)branch.size(); j++)
        {
            if (!IsTrackClusterEffective(id, branch[j], mean, rms, qsum))
                continue;

            double xy = CalTrackPosition(detID, xyflag, mean);
            double z = CalTrackPosition(detID, Z);
            g1->SetPoint(ipoint1++, z, xy);

            mmax = (qmax < qsum) ? mean : mmax;
            qmax = (qmax < qsum) ? qsum : qmax;
        }

        if (mmax != -1)
        {
            double xy = CalTrackPosition(detID, xyflag, mmax);
            double z = CalTrackPosition(detID, Z);
            g2->SetPoint(ipoint2++, z, xy);
        }
    }

    //拟合
    if (g1->GetN() < 2)
        return;
    if (g2->GetN() < 2)
        return;
    g1->Fit("pol1", "q");
    g2->Fit("pol1", "q");

    //target上击中的期望值
    double exp1 = g1->GetFunction("pol1")->Eval(CalTrackPosition(target, Z));
    double exp2 = g1->GetFunction("pol1")->Eval(CalTrackPosition(target, Z));

    //找出target在detlist里的编号
    int id = -1;
    for (int j = 0; j < detlist.size(); j++)
        if (detlist[j]->id == target)
            id = j;

    //如果target没有击中，则填默认值
    if (id == -1)
    {
        fDeviation[target][xyflag]->Fill(-25);
        return;
    }

    MyBeamTestData *detector = detlist[id];
    vector<vector<BeamHit>> branch = detector->XYbranch[xyflag];

    //循环target的所有击中，求残差
    for (int j = 0; j < (int)branch.size(); j++)
    {
        if (!IsTrackClusterEffective(id, branch[j], mean, rms, qsum))
            continue;

        double xy = CalTrackPosition(target, xyflag, mean);
        fDeviation[target][xyflag]->Fill(xy - exp2);
    }

    delete g1;
    delete g2;
}

//---------------------------------------------------------------------------
// 分析RICH数据
// 计算RICH在束流实验系下的hit坐标
void CalRICHPosition(int detID, double xch, double ych, double &x, double &y, double &z)
{
    x = 0;
    y = 0;
    z = 0;

    if (detID != RICHid)
        return;

    double posx = (-1 * (xch - 1) + 32 / 2. + 1 / 2.) * 5.0; //在探测器坐标系下的位置
    double posy = ((ych - 1) - 32 / 2. + 1 / 2.) * 5.0;

    double rot = (180 - (90 - angle)) / 180. * TMath::Pi();
    double LtoAnode = 50.;

    x = Xoff[detID] + posx;
    y = Yoff[detID] - LtoAnode * cos(rot) + posy * sin(rot);
    z = Zpos[detID] + LtoAnode * sin(rot) + posy * cos(rot);
}

//RICH-pad读出
bool IsRICHClusterEffective(vector<BeamHit> cluster, double &xmean, double &xrms, double &ymean, double &yrms, double qsum, int control = 0)
{
    xmean = 0;
    ymean = 0;
    xrms = 0;
    yrms = 0;
    qsum = 0;

    //判断是否是放电信号或噪声信号
    for (int i = 0; i < (int)cluster.size(); i++)
        if ((processControl == 3 || processControl == 5 || processControl == 6 || processControl == 7) && (cluster[i].t < 230 || cluster[i].t > 260 || cluster[i].q > 3600 || cluster[i].q < 10))
            return false;

    ftmp2->Reset();
    for (int i = 0; i < (int)cluster.size(); i++)
    {
        qsum += cluster[i].q;
        if ((processControl == 3 || processControl == 5) && (fabs(cluster[i].hit.first - 15) > 3 || fabs(cluster[i].hit.second - 3) > 4))
            return false;
        if ((processControl == 7) && (fabs(cluster[i].hit.first - 15) < 5 && fabs(cluster[i].hit.second - 3) < 5))
            return false;

        ftmp2->Fill(cluster[i].hit.first, cluster[i].hit.second, cluster[i].q);
    }

    xmean = ftmp2->GetMean(1);
    ymean = ftmp2->GetMean(2);
    xrms = ftmp2->GetRMS(1);
    yrms = ftmp2->GetRMS(2);

    return true;
}

void FillRICHHistogram(vector<MyBeamTestData *> detlist, int detID)
{
    for (int i = 0; i < (int)detlist.size(); i++)
    {
        if (detlist[i]->id != detID)
            continue;

        MyBeamTestData *detector = detlist[i];

        // cluster 分布
        fNCluster[detID][0]->Fill(detector->branch.size());
        for (int j = 0; j < (int)detector->branch.size(); j++)
            fCluster[detID][0]->Fill(detector->branch[j].size());

        //电荷Q/达峰时间T分布
        double qsum = 0;
        for (int j = 0; j < (int)detector->charge.size(); j++)
        {
            qsum += detector->charge[j];
            //fCharge[detID]->Fill(detector->charge[j]);
        }
        fTotCharge[detID]->Fill(qsum);

        for (int j = 0; j < (int)detector->time.size(); j++)
        {
            //fPeakingTime[detID]->Fill(detector->time[j]);
            //fQT[detID]->Fill(detector->charge[j], detector->time[j]);
        }

        //fhitmap分布
        vector<vector<BeamHit>> branch = detector->branch;

        double xmean, xrms;
        double ymean, yrms;
        double x, y, z;

        for (int ii = 0; ii < (int)branch.size(); ii++)
        {
            for (int iii = 0; iii < (int)branch[ii].size(); iii++)
                fhitmap[detID][0]->Fill(branch[ii][iii].hit.first, branch[ii][iii].hit.second);

            if (!IsRICHClusterEffective(branch[ii], xmean, xrms, ymean, yrms, qsum))
                continue;

            for (int iii = 0; iii < (int)branch[ii].size(); iii++)
            {
                fCharge[detID]->Fill(branch[ii][iii].q);
                fPeakingTime[detID]->Fill(branch[ii][iii].t);
                fQT[detID]->Fill(branch[ii][iii].q, branch[ii][iii].t);

                fhitmap[detID][1]->Fill(branch[ii][iii].hit.first, branch[ii][iii].hit.second);

                //在束流实验室系下的坐标在这里写出来
                CalRICHPosition(detID, branch[ii][iii].hit.first, branch[ii][iii].hit.second, x, y, z);
                fhitmapReal[detID][0]->Fill(x, y);
            }
            CalRICHPosition(detID, xmean, ymean, x, y, z);
            fhitmapReal[detID][1]->Fill(x, y);
        }
    }
}

void AnalysisRICH(vector<int> uselist, int target, vector<MyBeamTestData *> detlist)
{
    TGraphErrors *g[2];

    double mean, rms, qsum;
    for (int ixy = 0; ixy < 2; ixy++)
    {
        g[ixy] = new TGraphErrors();
        int ipoint = 0;
        for (int i = 0; i < (int)uselist.size(); i++)
        {
            //作为use的track必须要有击中
            int id = -1;
            int detID = uselist[i];
            for (int j = 0; j < detlist.size(); j++)
                if (detlist[j]->id == detID)
                    id = j;
            if (id == -1)
                return;

            MyBeamTestData *detector = detlist[id];
            vector<vector<BeamHit>> branch = detector->XYbranch[ixy];

            //作为use的track的number of cluster不能太多
            if (branch.size() >= 2)
                return;

            //填图以拟合
            double qmax = -1;
            double mmax = -1;
            double rmax = -1;
            for (int j = 0; j < (int)branch.size(); j++)
            {
                if (!IsTrackClusterEffective(id, branch[j], mean, rms, qsum))
                    continue;

                mmax = (qmax < qsum) ? mean : mmax;
                rmax = (qmax < qsum) ? rms : rmax;
                qmax = (qmax < qsum) ? qsum : qmax;
            }

            if (mmax != -1)
            {
                double xy = CalTrackPosition(detID, ixy, mmax);
                double z = CalTrackPosition(detID, Z);
                double zerr = 0.5; //z的测量误差
                g[ixy]->SetPoint(ipoint, z, xy);
                g[ixy]->SetPointError(ipoint, zerr, rmax * 0.4);
                g[ixy]->SetPoint(ipoint++, z, xy);
            }
        }

        if (g[ixy]->GetN() < 2)
            return;
        g[ixy]->Fit("pol1", "q");
    }

    //找RICH的击中
    //在mathematic里画图
    double LtoAnode = 50.;
    double LtoQuartzTop = -53;
    double LtoQuartzBot = -43;
    double rot = (180 - (90 - angle)) / 180 * TMath::Pi();

    //需要求解入射点的坐标，然后将击中信息转换为对应的坐标系下。
    //坐标系为入射点为o点，入射粒子沿着-y和-z的方向入射
    //所以首先需要粒子和辐射体上表面的交点坐标
    My3DLine *eBeamLine = new My3DLine();
    My3DLine *quartzTop = new My3DLine();
    My3DLine *quartzBot = new My3DLine();
    My3DLine *anodePlat = new My3DLine();
    eBeamLine->SetLine(g[1]->GetFunction("pol1"));
    quartzTop->SetSurface(rot, make_pair(Zpos[target], Yoff[target]), LtoQuartzTop);
    quartzBot->SetSurface(rot, make_pair(Zpos[target], Yoff[target]), LtoQuartzBot);
    anodePlat->SetSurface(rot, make_pair(Zpos[target], Yoff[target]), LtoAnode);

    double origin[3]; //原点为辐射体上表面和束流的交点
    origin[2] = quartzTop->CrossLineZ(eBeamLine);
    origin[1] = quartzTop->CrossLineY(eBeamLine);
    origin[0] = g[0]->GetFunction("pol1")->Eval(origin[2]);

    double foot[3]; // 原点到阳极面上的垂足
    foot[2] = anodePlat->PedalZ(make_pair(origin[2], origin[1]));
    foot[1] = anodePlat->PedalY(make_pair(origin[2], origin[1]));
    foot[0] = g[0]->GetFunction("pol1")->Eval(foot[2]);

    double cross[3]; // 电离信号是阳极板和束流的交点
    cross[2] = anodePlat->CrossLineZ(eBeamLine);
    cross[1] = anodePlat->CrossLineY(eBeamLine);
    cross[0] = g[0]->GetFunction("pol1")->Eval(cross[2]);

    if (verbose)
    {
        cout << "--------------------------------" << endl;
        cout << "Show[ContourPlot[{" << endl;
        cout << "y == MyLine[";
        eBeamLine->Print();
        cout << ", z], (*beamline*)" << endl;

        cout << "y == MyLine[";
        quartzTop->Print();
        cout << ", z], (*quartzTop*)" << endl;

        cout << "y == MyLine[";
        quartzBot->Print();
        cout << ", z], (*quartzBot*)" << endl;

        cout << "y == MyLine[";
        anodePlat->Print();
        cout << ", z] (*anode*)" << endl;
        cout << "}, {z, 0, 500}, {y, -250, 250}, ContourStyle -> {{Red, Dashed}, {Blue, Dotted}, {Blue, Dotted}, {Orange, Dotted}, ImageSize -> 800}]," << endl;

        for (int i = 0; i < g[1]->GetN(); i++)
        {
            double z, y;
            g[1]->GetPoint(i, z, y);
            cout << "Graphics[{{PointSize[Large], Blue, Point[{" << z << ", " << y << "}]}}], (*beamline*)" << endl;
        }

        cout << "Graphics[{{PointSize[Large], Red, Point[{" << Zpos[target] << ", " << Yoff[target] << "}]}, {Dotted, Thin, Circle[{" << Zpos[target] << ", " << Yoff[target] << "}, 45]}}], (*center*)" << endl;
        cout << "Graphics[{{PointSize[Large], Blue, Point[{" << origin[2] << ", " << origin[1] << "}]}}], (*origin*)" << endl;
        cout << "Graphics[{{PointSize[Large], Orange, Point[{" << cross[2] << ", " << cross[1] << "}]}}], (*cross*)" << endl;
        cout << "Graphics[{{PointSize[Large], Pink, Point[{" << foot[2] << ", " << foot[1] << "}]}}], (*foot*)" << endl;
    }

    //-----------
    // 填残差图
    //找出target在detlist里的编号
    int id = -1;
    for (int j = 0; j < detlist.size(); j++)
        if (detlist[j]->id == target)
            id = j;

    //如果target没有击中，则填默认值
    if (id == -1)
    {
        fDeviation[target][0]->Fill(-25);
        fDeviation[target][1]->Fill(-25);
        return;
    }

    MyBeamTestData *detector = detlist[id];
    vector<vector<BeamHit>> branch = detector->branch;

    //循环target的所有击中，求残差
    double xmean, xrms;
    double ymean, yrms;
    double x, y, z;
    double recMean = 0, nrec = 0;
    for (int i = 0; i < (int)branch.size(); i++)
    {
        if (!IsRICHClusterEffective(branch[i], xmean, xrms, ymean, yrms, qsum))
            continue;

        CalRICHPosition(target, xmean, ymean, x, y, z);
        if (verbose)
            cout << "Graphics[{{PointSize[Large], Green, Point[{" << z << ", " << y << "}]}}], (*pad*)" << endl;

        //填图残差
        if (processControl == 5)
        {
            double expRICH[2];
            expRICH[0] = g[0]->GetFunction("pol1")->Eval(z);
            expRICH[1] = g[1]->GetFunction("pol1")->Eval(z);

            fDeviation[target][0]->Fill(expRICH[0] - x);
            fDeviation[target][1]->Fill(expRICH[1] - y);
        }

        //角度重建
        if (processControl == 7)
        {
            //给出在探测器坐标系下的击中位置的坐标
            //探测器坐标系的零点为辐射体上表面与束流的交点》在阳极板上的投影
            double Xr = -1 * (x - foot[0]);
            double Yr = sqrt(pow(y - foot[1], 2) + pow(z - foot[2], 2)) * (foot[1] - y) / fabs(y - foot[1]);
            //cout<<(y - foot[1]) / fabs(y - foot[1])<<endl;
            if (Yr > 40)
                continue;
            double rec = ReconstructRICHByBeta(Xr, Yr);
            fFootmap->SetPoint(fFootmap->GetN() + 1, foot[0], foot[1]);
            fHitmap->Fill(Xr, Yr);
            fRec->Fill(rec);
            gRec->SetPoint(gRec->GetN() + 1, Xr, Yr, rec);
            recMean += rec;
            nrec++;
            if (verbose)
                cout << "(*-->Real=(" << Xr << ", " << Yr << ") from (" << x << ", " << y << ", " << z << ") and foot=(" << foot[0] << ", " << foot[1] << ", " << foot[2] << "), rec=" << rec << "*)" << endl;
        }
    }

    if (nrec > 0)
    {
        double rec = recMean / nrec;
        fRec2->Fill(rec);
        if (fabs(rec - 0.9) < 3 * 0.03) //认为是一个有效的光子击中
        {
            int nclu = 0;
            fNRec->Fill(nrec);
            for (int i = 0; i < (int)branch.size(); i++) //光子的Q和T
            {
                if (!IsRICHClusterEffective(branch[i], xmean, xrms, ymean, yrms, qsum))
                    continue;

                CalRICHPosition(target, xmean, ymean, x, y, z);
                double Xr = -1 * (x - foot[0]);
                double Yr = sqrt(pow(y - foot[1], 2) + pow(z - foot[2], 2)) * (foot[1] - y) / fabs(y - foot[1]);
                if (Yr > -40 || y < 40)
                    continue;
                double rec = ReconstructRICHByBeta(Xr, Yr);
                if (fabs(rec - 0.9) > 3 * 0.03)
                    continue;

                nclu++;
                double sum = 0;
                fClusterPhoton->Fill(branch[i].size());
                for (int j = 0; j < (int)branch[i].size(); j++)
                {
                    fQphoton->Fill(branch[i][j].q);
                    fTphoton->Fill(branch[i][j].t);
                    fQTphoton->Fill(branch[i][j].q, branch[i][j].t);
                    sum += branch[i][j].q;
                }
                fTotQphoton->Fill(sum);
            }
            fNClusterPhoton->Fill(branch.size()); //fNClusterPhoton 应该和 fNRec 分布差不多
        }
    }

    if (verbose)
        cout << "Graphics[{{Blue, Dotted, Line[{{" << origin[2] << ", " << origin[1] << "}, {" << foot[2] << ", " << foot[1] << "}}]}}]]" << endl;
}

//---------------------------------------------------------------------------
// 用GBL重建径迹
/*
inline Eigen::Matrix<double, 5, 5> Jac55new(double ds)
{
    Eigen::Matrix<double, 5, 5> jac = Eigen::Matrix<double, 5, 5>::Identity();
    // jac.UnitMatrix();
    jac(3, 1) = ds; // x = x0 + xp * ds
    jac(4, 2) = ds; // y = y0 + yp * ds
    return jac;
}

void GBLAnalysis(vector<int> uselist, vector<MyBeamTestData *> detlist)
{
    TGraphErrors *g[2];

    double mean, rms, qsum;
    for (int ixy = 0; ixy < 2; ixy++)
    {
        g[ixy] = new TGraphErrors();
        int ipoint = 0;
        for (int i = 0; i < (int)uselist.size(); i++)
        {
            //作为use的track必须要有击中
            int id = -1;
            int detID = uselist[i];
            for (int j = 0; j < detlist.size(); j++)
                if (detlist[j]->id == detID)
                    id = j;
            if (id == -1)
                continue;

            MyBeamTestData *detector = detlist[id];
            vector<vector<BeamHit>> branch = detector->XYbranch[ixy];

            //作为use的track的number of cluster不能太多
            if (branch.size() >= 2)
                continue;

            //填图以拟合
            double qmax = -1;
            double mmax = -1;
            double rmax = -1;
            for (int j = 0; j < (int)branch.size(); j++)
            {
                if (!IsTrackClusterEffective(id, branch[j], mean, rms, qsum))
                    continue;

                mmax = (qmax < qsum) ? mean : mmax;
                rmax = (qmax < qsum) ? rms : rmax;
                qmax = (qmax < qsum) ? qsum : qmax;
            }

            if (mmax != -1)
            {
                double xy = CalTrackPosition(detID, ixy, mmax);
                double z = CalTrackPosition(detID, Z);
                double zerr = 0.5; //z的测量误差
                g[ixy]->SetPoint(ipoint, z, xy);
                g[ixy]->SetPointError(ipoint, zerr, rmax * 0.4);
                g[ixy]->SetPoint(ipoint++, z, xy);
            }
        }

        if (g[ixy]->GetN() < 2)
            return;
        g[ixy]->Fit("pol1", "q");
    }
}
*/


//---------------------------------------------------------------------------
// 分析数据的控制流程
void AnalysisEvent(MyBeamTestHitData *fDSTEvent)
{
    if (processControl == 1 || processControl == 999)
    {
        FillTrackHistogram(fDSTEvent->detector, T02id);
        FillTrackHistogram(fDSTEvent->detector, T03id);
        FillTrackHistogram(fDSTEvent->detector, T06id);
        FillTrackHistogram(fDSTEvent->detector, T04id);
        FillTrackHistogram(fDSTEvent->detector, T05id);
    }

    if (processControl == 2 || processControl == 999)
    {
        AnalysisTrack(vector<int>{T02id, T03id}, T06id, fDSTEvent->detector, X);
        AnalysisTrack(vector<int>{T02id, T03id}, T06id, fDSTEvent->detector, Y);
        AnalysisTrack(vector<int>{T02id, T06id}, T03id, fDSTEvent->detector, X);
        AnalysisTrack(vector<int>{T02id, T06id}, T03id, fDSTEvent->detector, Y);
        AnalysisTrack(vector<int>{T03id, T06id}, T02id, fDSTEvent->detector, X);
        AnalysisTrack(vector<int>{T03id, T06id}, T02id, fDSTEvent->detector, Y);
        //AnalysisTrack(vector<int>{T02id, T03id, T05id, T06id}, T04id, fDSTEvent->detector, X);
        //AnalysisTrack(vector<int>{T02id, T03id, T05id, T06id}, T04id, fDSTEvent->detector, Y);
        AnalysisTrack(vector<int>{T02id, T03id, T06id}, T05id, fDSTEvent->detector, X);
        AnalysisTrack(vector<int>{T02id, T03id, T06id}, T05id, fDSTEvent->detector, Y);
    }

    if (processControl == 3 || processControl == 4 || processControl == 6 || processControl == 999)
        FillRICHHistogram(fDSTEvent->detector, RICHid);

    if (processControl == 5 || processControl == 7 || processControl == 999)
        AnalysisRICH(vector<int>{T02id, T03id, T06id}, RICHid, fDSTEvent->detector);
    
//    if (processControl == 8) 
//        GBLAnalysis(vector<int>{T02id, T03id, RICHid, T05id, T06id}, fDSTEvent->detector);
}

//---------------
// draw画图
void DrawTracker(TString cname, int id)
{
    TCanvas *cc = new TCanvas(cname);
    cc->SetTitle(fName + " " + detName[id]);
    cc->Clear();
    cc->Divide(4, 4);
    cc->cd(1);
    fhitmap[id][0]->Draw("colz");
    cc->cd(2);
    fNCluster[id][0]->Draw("");
    cc->cd(3);
    fCluster[id][0]->Draw("");
    cc->cd(4);
    fNCluster2D[id]->Draw("boxtext");
    cc->cd(5);
    fhitmap[id][1]->Draw("colz");
    cc->cd(6);
    fNCluster[id][1]->Draw("");
    cc->cd(7);
    fCluster[id][1]->Draw("");
    cc->cd(8);
    fTotCharge[id]->Draw("");
    cc->cd(9);
    fhitmapReal[id][0]->Draw("colz");
    cc->cd(10);
    fPeakingTime[id]->Draw("");
    cc->cd(11);
    fCharge[id]->Draw("");
    cc->cd(12);
    fQT[id]->Draw("colz");
    cc->cd(13);
    fhitmapReal[id][1]->Draw("colz");
    cc->cd(14);
    fhitmapReal[id][1]->ProjectionX()->Draw();
    cc->cd(15);
    fhitmapReal[id][1]->ProjectionY()->Draw();
}

void DrawRICH(TString cname, int id)
{
    TCanvas *cc = new TCanvas(cname);
    cc->SetTitle(fName + " " + detName[id]);
    cc->Clear();
    cc->Divide(4, 3);
    cc->cd(1);
    cc->GetPad(1)->SetLogz();
    fhitmap[id][0]->Draw("colz");
    cc->cd(2);
    fNCluster[id][0]->Draw("");
    cc->cd(3);
    fCluster[id][0]->Draw("");
    cc->cd(4);
    fTotCharge[id]->Draw("");

    cc->cd(5);
    cc->GetPad(5)->SetLogz();
    fhitmap[id][1]->Draw("colz");
    cc->cd(6);
    fPeakingTime[id]->Draw("");
    cc->cd(7);
    fCharge[id]->Draw("");
    cc->cd(8);
    fQT[id]->Draw("colz");

    cc->cd(9);
    cc->GetPad(9)->SetLogz();
    fhitmapReal[id][0]->Draw("colz");
    cc->cd(10);
    fhitmapReal[id][0]->ProjectionX()->Draw();
    cc->cd(11);
    fhitmapReal[id][0]->ProjectionY()->Draw();

    cc->cd(12);
    cc->GetPad(12)->SetLogz();
    fhitmapReal[id][1]->Draw("colz");
}

//---------------
// 输出信息
void PrintOffset(int xyz)
{
    if (xyz == ANGLE)
    {
        cout << "{" << fName.ReplaceAll("RUN", "") << ", " << angle << "}, // angle" << endl;
        return;
    }

    cout.precision(4);
    cout << "{" << fName.ReplaceAll("RUN", "") << ", ";
    for (int i = 0; i < NDET; i++) //FG125111
    {
        if (xyz == X)
        {
            if (i == RICHid)
            {
                if (processControl == 5 && fDeviation[i][0]->GetFunction("gaus") != 0)
                    cout << fDeviation[i][0]->GetFunction("gaus")->GetParameter(1) + Xoff[i];
                else
                    cout << -1 * fhitmapReal[i][0]->GetMean(1) + Xoff[i];
            }
            else
            {
                if (i != T02id && i != T03id && fDeviation[i][0]->GetFunction("gaus") != 0) //固定前两个tracker去标定后面几个tracker的位置
                    cout << -1 * fDeviation[i][0]->GetFunction("gaus")->GetParameter(1) + Xoff[i];
                else
                    cout << -1 * fhitmapReal[i][0]->GetMean(1) + Xoff[i];
            }
        }

        if (xyz == Y)
        {
            if (i == RICHid)
            {
                if (processControl == 5 && fDeviation[i][0]->GetFunction("gaus") != 0)
                    cout << fDeviation[i][1]->GetFunction("gaus")->GetParameter(1) + Yoff[i];
                else
                    cout << -1 * fhitmapReal[i][0]->GetMean(2) + Yoff[i];
            }
            else
            {
                if (i != T02id && i != T03id && fDeviation[i][1]->GetFunction("gaus") != 0)
                    cout << -1 * fDeviation[i][1]->GetFunction("gaus")->GetParameter(1) + Yoff[i];
                else
                    cout << -1 * fhitmapReal[i][0]->GetMean(2) + Yoff[i];
            }
        }
        if (xyz == Z)
            cout << Zpos[i];
        if (i != NDET - 1)
            cout << ", ";
    }
    cout << "}, // " << xyName[xyz] << endl;
    cout.precision(6);
}

void PrintInfo()
{
    cout << "\n\n----------------Offset for " << fName << "-------------------\n";
    PrintOffset(X);
    PrintOffset(Y);
    PrintOffset(Z);
    PrintOffset(ANGLE);
}

//---------------
// 读取offset
void ReadOffset()
{
    for (int i = 0; i < (int)offset.size(); i++)
    {
        if (fName == TString(Form("RUN%02d", (int)offset[i][0])))
        {
            cout << "==> Reading offset from : " << offset[i][0] << ". " << endl;
            for (int j = 0; j < NDET; j++)
            {
                Xoff[j] = offset[i + 0][j + 1];
                Yoff[j] = offset[i + 1][j + 1];
                Zpos[j] = offset[i + 2][j + 1];
            }
            angle = offset[i + 3][1];
            return;
        }
    }

    cout << "Warning: can't find the offset for " << fName << ". Please check." << endl;
}

void DrawDeviation(TString cname, vector<int> uselist)
{
    TCanvas *cc = new TCanvas(cname);
    cc->SetTitle("Deviation distribution for " + fName);
    cc->Clear();
    cc->Divide(uselist.size(), 2);
    for (int i = 0; i < (int)uselist.size(); i++)
    {
        cc->cd(1 + i);
        fDeviation[uselist[i]][0]->Fit("gaus", "q");
        cc->cd(1 + i + uselist.size());
        fDeviation[uselist[i]][1]->Fit("gaus", "q");

        cout << detName[uselist[i]]
             << " efficiency_x = " << fDeviation[uselist[i]][0]->Integral(1, 100) / fDeviation[uselist[i]][0]->GetEntries()
             << " efficiency_y = " << fDeviation[uselist[i]][1]->Integral(1, 100) / fDeviation[uselist[i]][1]->GetEntries() << endl;
    }
}

void DrawReconstruction(TString cname, int target)
{
    TCanvas *cc = new TCanvas(cname);
    cc->SetTitle("Deviation distribution for " + fName);
    cc->Clear();
    cc->Divide(5, 2);
    cc->cd(1);
    fRec->SetXTitle("Reconstructed Cherenkov angle(rad)");
    fRec->SetYTitle("Entries");
    fRec->Draw();
    cc->cd(6);
    fRec2->SetXTitle("Reconstructed Cherenkov angle(rad)");
    fRec2->SetYTitle("Entries");
    fRec2->Fit("gaus");
    cc->cd(2);
    //gRec->SetTitle("Reconstruction 2D map");
    //gRec->Draw("pcol");
    cc->cd(7);
    fHitmap->Draw("colz");
    cc->cd(3);
    fFootmap->SetTitle("Foot position");
    fFootmap->GetXaxis()->SetTitle("X(mm)");
    fFootmap->GetYaxis()->SetTitle("Y(mm)");
    fFootmap->Draw("ap");
    cc->cd(9);
    fNRec->Draw();
    cc->cd(4);
    fTotQphoton->Draw();
    cc->cd(8);
    fQTphoton->Draw("colz");
    cc->cd(5);
    fClusterPhoton->Draw();
    cc->cd(10);
    fNClusterPhoton->Draw();
}

//---------------
// 主函数
void checkCMBRoot()
{
    //TFile froot("./beamtest_hist.root");
    //fhmap = (TH2D *)froot.Get("fHitMap0");

    //--------------------------
    // 读取combine-dst.root
    TString fileName = "./BES/" + fName + "/Combine/Combined-dst.root";

    TFile *fDSTFile = new TFile(fileName);
    if (!fDSTFile->IsOpen())
        return;
    cout << "--> Reading " << fileName << "." << endl;

    TTree *fDSTTree = (TTree *)fDSTFile->Get("tree");
    MyBeamTestHitData *fDSTEvent = 0;
    fDSTTree->SetBranchAddress("event", &fDSTEvent);
    DSTEntries = (DSTEntries == 10) ? 10 : fDSTTree->GetEntries();
    cout << "==> Total hit entries = " << DSTEntries << endl;

    ReadOffset();

    //--------------------------
    // 定义histogram & graph
    int nbin;

    ftmp1 = new TH1F("ftmp1", "tmp", 128, 0, 128);
    ftmp2 = new TH2F("ftmp2", "tmp", 32, 0, 32, 32, 0, 32);
    //光子重建相关
    fRec = new TH1F("fRec", "reconstructed cherenkov angle by single cluster", 100, 0.5, 1.5);           //3.1415926/2);
    fRec2 = new TH1F("fRec2", "reconstructed mean cherenkov angle by multiple clusters", 100, 0.5, 1.5); //3.1415926/2);
    fNRec = new TH1F("fNRec", "number reconstructed cherenkov photons", 10, 0, 10);                      //3.1415926/2);
    fQphoton = new TH1F("fQphoton", "charge for photon signal", 256, 0, 4096 * 5);
    fTphoton = new TH1F("fTphoton", "time for photon signal", 100, 200, 300);
    fQTphoton = new TH2F("fQTphoton", "Q vs T for photons", 256, 0, 4096, 100, 200, 300);
    fTotQphoton = new TH1F("fTotQphoton", "total photon charges", 256, 0, 4096);
    fNClusterPhoton = new TH1F("fNClusterPhoton", "number of cluster for photons", 10, 0, 10);
    fClusterPhoton = new TH1F("fClusterPhoton", "cluster size for photons", 10, 0, 10);

    gRec = new TGraph2D();
    gRec->GetXaxis()->SetTitle("X");
    gRec->GetYaxis()->SetTitle("Y");
    fFootmap = new TGraph(); //new TH2F("fFootmap", "foot map for RICH",  40, -40 / 2 * 5, 40 / 2 * 5, 40, -40 / 2 * 5, 40 / 2 * 5);
    fHitmap = new TH2F("fHitmap", "Real Position hit map", 40, -40 / 2 * 5, 40 / 2 * 5, 40, -40 / 2 * 5, 40 / 2 * 5);
    fHitmap->GetXaxis()->SetTitle("X");
    fHitmap->GetYaxis()->SetTitle("Y");

    //事例填图
    for (int i = 0; i < NDET; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            nbin = (i == RICHid) ? 32 : 128;
            fhitmap[i][j] = new TH2F(Form("fhitmap%d_%d", i, j), Form("hit map for %s", detName[i]), nbin, 0, nbin, nbin, 0, nbin);
            fhitmap[i][j]->GetXaxis()->SetTitle("X");
            fhitmap[i][j]->GetYaxis()->SetTitle("Y");

            nbin = (i == T06id) ? 25 : 5;
            nbin = (i == RICHid) ? 25 : nbin;
            fDeviation[i][j] = new TH1F(Form("fdev%d_%d", i, j), Form("Deviation of Track %s for %s", detName[i], xyName[j]), 100, -1 * nbin, nbin);
            fDeviation[i][j]->SetXTitle(Form("%s(mm)", xyName[j]));
            fDeviation[i][j]->SetYTitle("Entries");

            if (i == RICHid)
            {
                fNCluster[i][j] = new TH1F(Form("fnclu%d_%d", i, j), Form("Number of cluster of %s", detName[i]), 10, 0, 10);
                fCluster[i][j] = new TH1F(Form("fclu%d_%d", i, j), Form("Cluster size of %s", detName[i]), 10, 0, 10);
                fhitmapReal[i][j] = new TH2F(Form("fhitmapReal%d_%d", i, j), Form("Real Position hit map for %s", detName[i]), 40, -40 / 2 * 5, 40 / 2 * 5, 40, -40 / 2 * 5, 40 / 2 * 5);
            }
            else
            {
                fNCluster[i][j] = new TH1F(Form("fnclu%d_%d", i, j), Form("Number of cluster of %s for %s", detName[i], xyName[j]), 10, 0, 10);
                fCluster[i][j] = new TH1F(Form("fclu%d_%d", i, j), Form("Cluster size of %s for %s", detName[i], xyName[j]), 10, 0, 10);
                fhitmapReal[i][j] = new TH2F(Form("fhitmapReal%d_%d", i, j), Form("Real Position hit map for %s", detName[i]), 160, -160 / 2 * 0.4, 160 / 2 * 0.4, 160, -160 / 2 * 0.4, 160 / 2 * 0.4);
            }
        }
        //nbin = (i == 0) ? 32 : 128;
        fNCluster2D[i] = new TH2F(Form("fnclu2%d", i), Form("Number of cluster of %s", detName[i]), 10, 0, 10, 10, 0, 10);

        fPeakingTime[i] = new TH1F(Form("fT%d", i), Form("Peaking time of %s", detName[i]), 100, 0, 3000);
        fPeakingTime[i]->GetXaxis()->SetTitle("T");
        fPeakingTime[i]->GetYaxis()->SetTitle("Entries");

        fCharge[i] = new TH1F(Form("fQ%d", i), Form("Charge of %s", detName[i]), 256, 0, 4096);
        fCharge[i]->GetXaxis()->SetTitle("Q");
        fCharge[i]->GetYaxis()->SetTitle("Entries");

        fTotCharge[i] = new TH1F(Form("ftotQ%d", i), Form("Total charge of %s", detName[i]), 256, 0, 4096 * 5);
        fTotCharge[i]->GetXaxis()->SetTitle("Q");
        fTotCharge[i]->GetYaxis()->SetTitle("Entries");

        fQT[i] = new TH2F(Form("fQT%d", i), Form("Q vs. T of %s", detName[i]), 256, 0, 4096, 100, 200, 300);
        fQT[i]->GetXaxis()->SetTitle("Q");
        fQT[i]->GetYaxis()->SetTitle("T");
    }

    //--------------------------
    // 分析数据
    //DSTEntries = (DSTEntries > 20000) ? 20000 : DSTEntries;
    for (int ip = 0; ip < DSTEntries; ip++)
    {
        if (ip % 10000 == 0 || ip < 10)
            cout << "--> event: " << ip << endl;

        fDSTTree->GetEntry(ip);
        AnalysisEvent(fDSTEvent);
    }

    //--------------------------
    // 画图
    gStyle->SetOptFit(1);

    if (processControl == 1 || processControl == 999)
    {
        DrawTracker("c1", T02id);
        DrawTracker("c2", T03id);
        DrawTracker("c3", T06id);
        DrawTracker("c4", T04id);
        DrawTracker("c5", T05id);
    }

    if (processControl == 2 || processControl == 999)
    {
        DrawDeviation("c1", vector<int>{T02id, T03id, T05id, T06id});
    }
    
    if (processControl == 3 || processControl == 4 || processControl == 6 || processControl == 999)
        DrawRICH("c5", RICHid);

    if (processControl == 5 || processControl == 999)
        DrawDeviation("c6", vector<int>{RICHid});

    if (processControl == 7 || processControl == 999)
        DrawReconstruction("c7", RICHid);

    PrintInfo();
}
