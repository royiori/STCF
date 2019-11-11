#include "MyStyle.h"
#include "TRandom.h"

using namespace std;

//-------------------------------
// Draw 1D histogram
//
void MyStyle::DrawHistograms(TH1 *f1, double min, double max, TH1 *f2, TH1 *f3, TH1 *f4, TH1 *f5, TH1 *f6, TH1 *f7, TH1 *f8, TH1 *f9, TH1 *f10)
{
    gPresetH1.clear();
    SetHistogram(f1);
    SetHistogram(f2);
    SetHistogram(f3);
    SetHistogram(f4);
    SetHistogram(f5);
    SetHistogram(f6);
    SetHistogram(f7);
    SetHistogram(f8);
    SetHistogram(f9);
    SetHistogram(f10);
    Draw1DHistVector(gPresetH1, min, max);
}

void MyStyle::Draw1DHistVector(vector<TH1 *> hvec, double Min, double Max)
{
    if (hvec[0] == NULL)
        return;

    SetHist1Theme(hvec);

    gH1 = (TH1F *)hvec[0]->Clone(Form("__tmpH1__%d", (int)gRandom->Uniform() * 1000));
    gH1List.push_back(gH1);
    gH1->SetXTitle(Label[0]);
    gH1->SetYTitle(Label[1]);
    gH1->SetTitle(Title);

    double min = gH1->GetMinimum();
    double max = gH1->GetMaximum();

    for (int i = 0; i < (int)hvec.size(); i++)
    {
        min = (min > hvec[i]->GetMinimum()) ? hvec[i]->GetMinimum() : min;
        max = (max < hvec[i]->GetMaximum()) ? hvec[i]->GetMaximum() : max;
    }
    min = (Min != -1) ? Min : min * 0.95;
    max = (Max != -1) ? Max : max * 1.05;

    gH1->GetYaxis()->SetRangeUser(min, max);
    gH1->Draw(Option);
    for (int i = 1; i < (int)hvec.size(); i++)
        MyDraw(hvec[i], Option + TString(";same"));

    if (hvec[1] == NULL)
        return;
    TLegend *leg = new TLegend(0.8, 0.7, 1.0, 0.8);
    for (int i = 0; i < (int)hvec.size(); i++)
        MyAddLegend(leg, hvec[i], Legends[i], Option + TString("l"));
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw("same");
}

//-------------------------------
// Draw 2D histogram
//
void MyStyle::DrawHistograms(TH2 *f1, TH2 *f2, TH2 *f3, TH2 *f4, TH2 *f5, TH2 *f6, TH2 *f7, TH2 *f8, TH2 *f9, TH2 *f10)
{
    gPresetH2.clear();
    SetHistogram(f1);
    SetHistogram(f2);
    SetHistogram(f3);
    SetHistogram(f4);
    SetHistogram(f5);
    SetHistogram(f6);
    SetHistogram(f7);
    SetHistogram(f8);
    SetHistogram(f9);
    SetHistogram(f10);
    Draw2DHistVector(gPresetH2);
}

void MyStyle::Draw2DHistVector(vector<TH2 *> hvec)
{
    if (hvec[0] == NULL)
        return;

    SetHist2Theme(hvec);

    gH2 = (TH2F *)hvec[0]->Clone(Form("__tmpH2__%d", (int)gRandom->Uniform() * 1000));
    gH2List.push_back(gH2);

    gH2->SetXTitle(Label[0]);
    gH2->SetYTitle(Label[1]);
    gH2->SetTitle(Title);

    //// Y axis
    //gH2->GetYaxis()->SetAxisColor(2000);
    gH2->GetYaxis()->SetTickLength(0.01);

    //gH2->GetYaxis()->SetLabelColor(1001);
    //gH2->GetYaxis()->SetLabelFont(55);
    gH2->GetYaxis()->SetLabelSize(0.03);
    //gH2->GetYaxis()->SetLabelOffset(0.005);

    gH2->GetYaxis()->SetTitleFont(62);
    gH2->GetYaxis()->SetTitleSize(0.03);
    //gH2->GetYaxis()->SetTitleColor(1001);
    gH2->GetYaxis()->SetTitleOffset(1.5);

    //// X axis
    //gH2->GetXaxis()->SetAxisColor(2000);
    gH2->GetXaxis()->SetTickLength(0.01);

    //gH2->GetXaxis()->SetLabelColor(1001);
    //gH2->GetXaxis()->SetLabelFont(55);
    gH2->GetXaxis()->SetLabelSize(0.03);
    //gH2->GetXaxis()->SetLabelOffset(0.025);

    gH2->GetXaxis()->SetTitleFont(62);
    gH2->GetXaxis()->SetTitleSize(0.03);
    //gH2->GetXaxis()->SetTitleColor(1001);
    gH2->GetXaxis()->SetTitleOffset(1.20);
    gH2->Draw(Option);

    for (int i = 1; i < (int)hvec.size(); i++)
        MyDraw(hvec[i], Option + TString("same"));

    if ((int)hvec.size() == 1)
        return;
    TLegend *leg = new TLegend(0.8, 0.7, 1.0, 0.8);
    for (int i = 0; i < (int)hvec.size(); i++)
        MyAddLegend(leg, hvec[i], Legends[i], Option + TString("p"));
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw("same");
}

//-------------------------------
// Draw Graphs
//
void MyStyle::DrawGraphs(TGraph *g1, TGraph *g2, TGraph *g3, TGraph *g4, TGraph *g5, TGraph *g6, TGraph *g7, TGraph *g8, TGraph *g9, TGraph *g10)
{
    gPresetGph.clear();
    SetGraph(g1);
    SetGraph(g2);
    SetGraph(g3);
    SetGraph(g4);
    SetGraph(g5);
    SetGraph(g6);
    SetGraph(g7);
    SetGraph(g8);
    SetGraph(g9);
    SetGraph(g10);
    DrawGraphVector(gPresetGph);
}

void MyStyle::DrawGraphVector(vector<TGraph *> gvec)
{
    if (gvec[0] == NULL)
        return;

    SetGraphTheme(gvec);

    gG = (TGraph*)gvec[0]->Clone();
    gGList.push_back(gG);

    if (CP == MATHEMATIC_STYLE)
    {
        //// Y axis
        //gG->GetYaxis()->SetAxisColor(2000);
        //gG->GetYaxis()->SetTickLength(0.01);

        //gG->GetYaxis()->SetLabelFont(55);
        //gG->GetYaxis()->SetLabelSize(0.03);
        //gG->GetYaxis()->SetLabelColor(1001);
        //gG->GetYaxis()->SetLabelOffset(0.005);

        //gG->GetYaxis()->SetTitleFont(62);
        //gG->GetYaxis()->SetTitleSize(0.03);
        //gG->GetYaxis()->SetTitleColor(1001);
        //gG->GetYaxis()->SetTitleOffset(1.45);

        //// X axis
        //gG->GetXaxis()->SetAxisColor(2000);
        //gG->GetXaxis()->SetTickLength(0.01);

        //gG->GetXaxis()->SetLabelColor(1001);
        //gG->GetXaxis()->SetLabelFont(55);
        //gG->GetXaxis()->SetLabelSize(0.03);
        //gG->GetXaxis()->SetLabelOffset(0.025);

        //gG->GetXaxis()->SetTitleFont(62);
        //gG->GetXaxis()->SetTitleSize(0.03);
        //gG->GetXaxis()->SetTitleColor(1001);
        //gG->GetXaxis()->SetTitleOffset(1.20);
    }

    double min = gG->GetYaxis()->GetXmin();
    double max = gG->GetYaxis()->GetXmax();
    for (int i = 0; i < (int)gvec.size(); i++)
        min = (min > gvec[i]->GetYaxis()->GetXmin()) ? gvec[i]->GetYaxis()->GetXmin() : min;
    for (int i = 0; i < (int)gvec.size(); i++)
        max = (max < gvec[i]->GetYaxis()->GetXmax()) ? gvec[i]->GetYaxis()->GetXmax() : max;

    gG->SetTitle(Title);
    gG->GetXaxis()->SetTitle(Label[0]);
    gG->GetYaxis()->SetTitle(Label[1]);
    gG->GetYaxis()->SetRangeUser(min, max);
    gG->Draw(Option + "a");

    for (int i = 1; i < (int)gvec.size(); i++)
        MyDraw(gvec[i], Option + TString("same"));

    if (gvec[1] == NULL)
        return;
    TLegend *leg = new TLegend(0.8, 0.7, 1.0, 0.8);
    for (int i = 0; i < (int)gvec.size(); i++) 
        MyAddLegend(leg, gvec[i], Legends[i], Option + TString("l"));
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw("same");
}

//-------------------------------
// Plot theme
//
void MyStyle::UpdateTheme()
{
    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(0);

    if (cID <= 2002)
    {
        new TColor(cID++, 0x8f / 255., 0xb0 / 255., 0x31 / 255.);
        new TColor(cID++, 0xe2 / 255., 0x9e / 255., 0x29 / 255.);
        new TColor(cID++, 0x5e / 255., 0x81 / 255., 0xb5 / 255.);
        new TColor(cID++, 0xeb / 255., 0x63 / 255., 0x35 / 255.);
        new TColor(cID++, 0x84 / 255., 0x7a / 255., 0xaf / 255.);
        new TColor(cID++, 0xb9 / 255., 0x71 / 255., 0x31 / 255.);
        new TColor(cID++, 0x6b / 255., 0x9d / 255., 0xc4 / 255.);
        new TColor(cID++, 0xf6 / 255., 0xc1 / 255., 0x42 / 255.);
        new TColor(cID++, 0x9c / 255., 0x65 / 255., 0x99 / 255.);
        new TColor(cID++, 0x92 / 255., 0x94 / 255., 0x2e / 255.);
    }

    double lcolor1[__NHIST__] = {2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011};
    double lstyle1[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double lwidth1[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double mstyle1[__NHIST__] = {2, 3, 5, 8, 29, 13, 14, 15, 16, 18};
    double msize1[__NHIST__] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
    double fstyle1[__NHIST__] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    double lcolor2[__NHIST__] = {kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack};
    double lstyle2[__NHIST__] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double lwidth2[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double mstyle2[__NHIST__] = {2, 3, 5, 7, 29, 13, 14, 15, 16, 18};
    double msize2[__NHIST__] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double fstyle2[__NHIST__] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; //{3004,3005,3006,3007,3010,3013,3001,3002,3003,3021};

    double lcolor3[__NHIST__] = {kRed, kBlue, kOrange, kViolet, kTeal, kAzure, kMagenta, kYellow, kPink, kSpring}; //kCyan, kTeal,
    double lstyle3[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double lwidth3[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double mstyle3[__NHIST__] = {2, 3, 5, 7, 29, 13, 14, 15, 16, 18};
    double msize3[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double fstyle3[__NHIST__] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    switch (CP)
    {
    case 1: //mathematica 的 缺省设置
        StoreTheme(lcolor1, lstyle1, lwidth1, lcolor1, mstyle1, msize1, lcolor1, fstyle1);
        break;

    case 10: //黑白设置
        StoreTheme(lcolor2, lstyle2, lwidth2, lcolor2, mstyle2, msize2, lcolor2, fstyle2);
        break;

    default: //缺省设置
        StoreTheme(lcolor3, lstyle3, lwidth3, lcolor3, mstyle3, msize3, lcolor3, fstyle3);
        break;
    }
}

void MyStyle::StoreTheme(double *lcolor, double *lstyle, double *lwidth,
                         double *mcolor, double *mstyle, double *msize,
                         double *fcolor, double *fstyle)
{
    for (int i = 0; i < __NHIST__; i++)
    {
        LineColor[i] = lcolor[i];
        LineStyle[i] = lstyle[i];
        LineWidth[i] = lwidth[i];

        MarkerColor[i] = mcolor[i];
        MarkerStyle[i] = mstyle[i];
        MarkerSize[i] = msize[i];

        FillColor[i] = fcolor[i];
        FillStyle[i] = fstyle[i];
    }
}

void MyStyle::SetTheme(TH1 *hist, int id)
{
    if (hist == NULL)
        return;

    hist->SetLineColor(LineColor[id]);
    hist->SetLineStyle(LineStyle[id]);
    hist->SetLineWidth(LineWidth[id]);

    hist->SetMarkerColor(MarkerColor[id]);
    hist->SetMarkerStyle(MarkerStyle[id]);
    hist->SetMarkerSize(MarkerSize[id]);

    hist->SetFillColor(FillColor[id]);
    ///hist->SetFillStyle(FillStyle[id]);
}

void MyStyle::SetTheme(TH2 *hist, int id)
{
    if (hist == NULL)
        return;
    hist->SetLineColor(LineColor[id]);
    hist->SetLineStyle(LineStyle[id]);
    hist->SetLineWidth(LineWidth[id]);

    hist->SetMarkerColor(MarkerColor[id]);
    hist->SetMarkerStyle(MarkerStyle[id]);
    hist->SetMarkerSize(MarkerSize[id]);

    hist->SetFillColor(FillColor[id]);
    hist->SetFillStyle(FillStyle[id]);
}

void MyStyle::SetTheme(TGraph *graph, int id)
{
    if (graph == NULL)
        return;

    graph->SetLineColor(LineColor[id]);
    graph->SetLineStyle(LineStyle[id]);
    graph->SetLineWidth(LineWidth[id]);

    graph->SetMarkerColor(MarkerColor[id]);
    graph->SetMarkerStyle(MarkerStyle[id]);
    graph->SetMarkerSize(MarkerSize[id]);

    graph->SetFillColor(FillColor[id]);
    graph->SetFillStyle(FillStyle[id]);
}
