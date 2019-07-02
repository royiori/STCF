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

   // for 1D:
   myStyle->Draw2DHistograms(f1, 0, 100, f2);
   // or
   myStyle->ClearPreset1DHist();
   myStyle->SetHistogram(f1);
   myStyle->SetHistogram(f2);
   myStyle->Draw1DHistograms(0, 100);
   // or
   vector<TH1 *>vec;
   vec.push_back(f1);
   vec.push_back(f2);
   myStyle->Draw1DHistVector(vec, 0, 100);

   // for 2D:
   myStyle->Draw2DHistograms(h1, h2);
   //...
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
#include <iostream>

#define ROOT_STYLE 0
#define MATHEMATIC_STYLE 1

const int __NHIST__ = 10;
using namespace std;

class MyStyle
{
public:

   MyStyle() { CP=0; cID = 2002; UpdateTheme(); gH1=NULL; gH2=NULL; gG=NULL; };
   ~MyStyle() {};

   //缺省的颜色设置
   void SetColorPattern( int _CP = 0 ) {CP=_CP;}

   void SetTitle(const char *title) { Title    = TString(title); }
   void SetXLabel(const char *xlab) { Label[0] = TString(xlab);  }
   void SetYLabel(const char *ylab) { Label[1] = TString(ylab);  }
   void SetLegends(int id, const char *leg) { if(id < 0 || id >= __NHIST__) return; Legends[id] = TString(leg); }
   void SetDrawOption(const char *opt) { Option = TString(opt); }

   // 1D histograms
   void ClearPreset1DHist() { gPresetH1.clear(); }
   void SetHistogram(TH1 *f) { if(f!=NULL) gPresetH1.push_back(f);}
   void Draw1DHistograms(double min=-1, double max=-1) { Draw1DHistVector(gPresetH1, min, max); }
   void Draw1DHistVector(vector<TH1 *>hvec, double min=-1, double max=-1);
   void DrawHistograms(TH1 *f1, double min = -1, double max = -1, TH1 *f2=NULL, TH1 *f3=NULL, TH1 *f4=NULL, TH1 *f5=NULL, TH1 *f6=NULL, TH1 *f7=NULL, TH1 *f8=NULL, TH1 *f9=NULL, TH1 *f10=NULL);

   // 2D histograms
   void DrawHistograms(TH2 *f1, TH2 *f2=NULL, TH2 *f3=NULL, TH2 *f4=NULL, TH2 *f5=NULL, TH2 *f6=NULL, TH2 *f7=NULL, TH2 *f8=NULL, TH2 *f9=NULL, TH2 *f10=NULL);
   void SetHistogram(TH2 *f) { if(f!=NULL) gPresetH2.push_back(f);}
   void ClearPreset2DHist() { gPresetH2.clear(); }
   void Draw2DHistograms() { Draw2DHistVector(gPresetH2); }
   void Draw2DHistVector(vector<TH2 *>hvec);

   // Graphs
   void DrawGraphs(TGraph *g1, TGraph *g2=NULL, TGraph *g3=NULL, TGraph *g4=NULL, TGraph *g5=NULL, TGraph *g6=NULL, TGraph *g7=NULL, TGraph *g8=NULL, TGraph *g9=NULL, TGraph *g10=NULL);
   void SetGraph(TGraph *g) { if(g!=NULL) gPresetGph.push_back(g); }
   void ClearPresetGraphs() { gPresetGph.clear(); }
   void DrawPresetGraph();

   void SetGraphTheme(TGraph *g1, TGraph *g2=NULL, TGraph*g3=NULL, TGraph*g4=NULL, TGraph*g5=NULL, TGraph*g6=NULL, TGraph*g7=NULL);

private:




   void MyDraw(TH1 *f, TString opt) { if(f!=NULL) f->Draw(opt); }
   void MyDraw(TH2 *f, TString opt) { if(f!=NULL) f->Draw(opt); }
   void MyDraw(TGraph *g, TString opt) { if(g!=NULL) g->Draw(opt); }

   void MyAddLegend(TLegend* leg, TH1* f, TString text, TString opt) { if(f!=NULL) leg->AddEntry(f, text, opt); }
   void MyAddLegend(TLegend* leg, TH2* f, TString text, TString opt) { if(f!=NULL) leg->AddEntry(f, text, opt); }
   void MyAddLegend(TLegend* leg, TGraph* g, TString text, TString opt) { if(g!=NULL) leg->AddEntry(g, text, opt); }

   void SetHist1Theme() { UpdateTheme(); for (int i = 0; i < (int)gPresetH1.size(); i++) SetTheme(gPresetH1[i], i); }
   void SetHist2Theme() { UpdateTheme(); for (int i = 0; i < (int)gPresetH2.size(); i++) SetTheme(gPresetH2[i], i); }
   void SetGraphTheme() { UpdateTheme(); for (int i = 0; i < (int)gPresetGph.size(); i++) SetTheme(gPresetGph[i], i); }
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
   double MarkerSize[__NHIST__];

   double FillColor[__NHIST__];
   double FillStyle[__NHIST__];

   TH1F *gH1;
   TH2F *gH2;
   TGraph *gG;
   
   vector<TH1 *>gPresetH1;
   vector<TH2 *>gPresetH2;
   vector<TGraph *>gPresetGph;

   vector<TH1 *>gH1List;
   vector<TH2 *>gH2List;
   vector<TGraph *>gGList;

   TString Option;
   TString Title;
   TString Label[2];
   TString Legends[__NHIST__];
};

#endif
