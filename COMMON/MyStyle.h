//-----------------------------------------------------
// to draw multiple histograms
//
// usage as example:
//
/*
    MyStyle* myStyle = new MyStyle();
    myStyle->SetColorPattern(MATHEMATIC_STYLE);
    myStyle->SetDrawOption("c");
    myStyle->SetLegends(0, "Legend0");
    myStyle->SetLegends(1, "Legend1");
    myStyle->SetXLabel("XLabel");
    myStyle->SetYLabel("YLabel");
    myStyle->Draw2Graph(g1, g2);
*/
//
//  
//                            Q.LIU
//-----------------------------------------------------

#ifndef MYSTYLE_H
#define MYSTYLE_H

#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TColor.h>
#include <TLegend.h>
#include <TStyle.h>

#define ROOT_STYLE 0
#define MATHEMATIC_STYLE 1

const int __NHIST__ = 10;
using namespace std;

class MyStyle
{
public:

   MyStyle() { CP=0; cID = 2000; UpdateTheme(); gH1=NULL; gH2=NULL; gG=NULL;};
   ~MyStyle() {};

   //缺省的颜色设置
   void SetColorPattern( int _CP = 0 ) {CP=_CP;}

   void SetTitle(const char *title) { Title    = TString(title); }
   void SetXLabel(const char *xlab) { Label[0] = TString(xlab);  }
   void SetYLabel(const char *ylab) { Label[1] = TString(ylab);  }
   void SetLegends(int id, const char *leg) { if(id<0 || id>=__NHIST__) return; Legends[id] = TString(leg); }
   void SetDrawOption(const char *opt) { Option = TString(opt); }

   void Draw1Histograms(TH1 *f1, double min=0, double max=0) {DrawHistograms(min, max, f1, NULL, NULL, NULL, NULL, NULL, NULL);}
   void Draw2Histograms(TH1 *f1, TH1 *f2, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, NULL, NULL, NULL, NULL, NULL);}
   void Draw3Histograms(TH1 *f1, TH1 *f2, TH1 *f3, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, NULL, NULL, NULL, NULL);}
   void Draw4Histograms(TH1 *f1, TH1 *f2, TH1 *f3, TH1 *f4, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, f4, NULL, NULL, NULL);}
   void Draw5Histograms(TH1 *f1, TH1 *f2, TH1 *f3, TH1 *f4, TH1 *f5, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, f4, f5, NULL, NULL);}
   void Draw6Histograms(TH1 *f1, TH1 *f2, TH1 *f3, TH1 *f4, TH1 *f5, TH1 *f6, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, f4, f5, f6, NULL);}
   void Draw7Histograms(TH1 *f1, TH1 *f2, TH1 *f3, TH1 *f4, TH1 *f5, TH1 *f6, TH1 *f7, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, f4, f5, f6, f7);}

   void Draw1Histograms(TH2 *f1, double min=0, double max=0) {DrawHistograms(min, max, f1, NULL, NULL, NULL, NULL, NULL, NULL);}
   void Draw2Histograms(TH2 *f1, TH2 *f2, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, NULL, NULL, NULL, NULL, NULL);}
   void Draw3Histograms(TH2 *f1, TH2 *f2, TH2 *f3, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, NULL, NULL, NULL, NULL);}
   void Draw4Histograms(TH2 *f1, TH2 *f2, TH2 *f3, TH2 *f4, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, f4, NULL, NULL, NULL);}
   void Draw5Histograms(TH2 *f1, TH2 *f2, TH2 *f3, TH2 *f4, TH2 *f5, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, f4, f5, NULL, NULL);}
   void Draw6Histograms(TH2 *f1, TH2 *f2, TH2 *f3, TH2 *f4, TH2 *f5, TH2 *f6, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, f4, f5, f6, NULL);}
   void Draw7Histograms(TH2 *f1, TH2 *f2, TH2 *f3, TH2 *f4, TH2 *f5, TH2 *f6, TH2 *f7, double min=0, double max=0) {DrawHistograms(min, max, f1, f2, f3, f4, f5, f6, f7);}

   void Draw1Graph(TGraph *g1) {DrawGraphs(g1, NULL, NULL, NULL, NULL, NULL, NULL);}
   void Draw2Graph(TGraph *g1, TGraph *g2) {DrawGraphs(g1, g2, NULL, NULL, NULL, NULL, NULL);}
   void Draw3Graph(TGraph *g1, TGraph *g2, TGraph*g3) {DrawGraphs(g1, g2, g3, NULL, NULL, NULL, NULL);}
   void Draw4Graph(TGraph *g1, TGraph *g2, TGraph*g3, TGraph*g4) {DrawGraphs(g1, g2, g3, g4, NULL, NULL, NULL);}
   void Draw5Graph(TGraph *g1, TGraph *g2, TGraph*g3, TGraph*g4, TGraph*g5) {DrawGraphs(g1, g2, g3, g4, g5, NULL, NULL);}
   void Draw6Graph(TGraph *g1, TGraph *g2, TGraph*g3, TGraph*g4, TGraph*g5, TGraph*g6) {DrawGraphs(g1, g2, g3, g4, g5, g6, NULL);}
   void Draw7Graph(TGraph *g1, TGraph *g2, TGraph*g3, TGraph*g4, TGraph*g5, TGraph*g6, TGraph*g7) {DrawGraphs(g1, g2, g3, g4, g5, g6, g7);}
   void DrawGraphs(vector<TGraph*>);

   void SetHistTheme(TH1 *f1, TH1 *f2=NULL, TH1 *f3=NULL, TH1 *f4=NULL, TH1 *f5=NULL, TH1 *f6=NULL, TH1 *f7=NULL);
   void SetHistTheme(TH2 *f1, TH2 *f2=NULL, TH2 *f3=NULL, TH2 *f4=NULL, TH2 *f5=NULL, TH2 *f6=NULL, TH2 *f7=NULL);
   void SetGraphTheme(TGraph *g1, TGraph *g2=NULL, TGraph*g3=NULL, TGraph*g4=NULL, TGraph*g5=NULL, TGraph*g6=NULL, TGraph*g7=NULL);

private:

   void DrawHistograms(double min, double max, TH1 *f1, TH1 *f2, TH1 *f3, TH1 *f4, TH1 *f5, TH1 *f6, TH1 *f7);
   void DrawHistograms(double min, double max, TH2 *f1, TH2 *f2, TH2 *f3, TH2 *f4, TH2 *f5, TH2 *f6, TH2 *f7);

   void DrawGraphs(TGraph *g1, TGraph *g2, TGraph*g3, TGraph*g4, TGraph*g5, TGraph*g6, TGraph*g7);

   void SetTheme(TH1 *hist, int id);
   void SetTheme(TH2 *hist, int id);
   void SetTheme(TGraph *graph, int id);
   void UpdateTheme();
   void StoreTheme(double *lcolor, double *lstyle, double *lwidth, 
                   double *mcolor, double *mstyle, double *msize,  
                   double *fcolor, double *fstyle);

   int CP;
   int cID;  //new defined color ID

   double LineColor[__NHIST__];
   double LineStyle[__NHIST__];
   double LineWidth[__NHIST__];

   double MarkerColor[__NHIST__];
   double MarkerStyle[__NHIST__];
   double MarkerSize [__NHIST__];

   double FillColor[__NHIST__];
   double FillStyle[__NHIST__];

   TH1F *gH1;
   TH2F *gH2;
   TGraph *gG;

   TString Legends[__NHIST__];
   TString Title;
   TString Label[2];
   TString Option;
};

#endif
