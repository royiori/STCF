#include "MyStyle.h"

using namespace std;

void MyStyle::SetHistTheme(TH1 *f1, TH1 *f2, TH1 *f3, TH1 *f4, TH1 *f5, TH1 *f6, TH1 *f7)
{
    UpdateTheme();

    SetTheme(f1, 0);
    if(f2!=NULL) SetTheme(f2, 1);
    if(f3!=NULL) SetTheme(f3, 2);
    if(f4!=NULL) SetTheme(f4, 3);
    if(f5!=NULL) SetTheme(f5, 4);
    if(f6!=NULL) SetTheme(f6, 5);
    if(f7!=NULL) SetTheme(f7, 6);
}

void MyStyle::SetHistTheme(TH2 *f1, TH2 *f2, TH2 *f3, TH2 *f4, TH2 *f5, TH2 *f6, TH2 *f7)
{
    UpdateTheme();

    SetTheme(f1, 0);
    if(f2!=NULL) SetTheme(f2, 1);
    if(f3!=NULL) SetTheme(f3, 2);
    if(f4!=NULL) SetTheme(f4, 3);
    if(f5!=NULL) SetTheme(f5, 4);
    if(f6!=NULL) SetTheme(f6, 5);
    if(f7!=NULL) SetTheme(f7, 6);
}


void MyStyle::DrawHistograms(double min, double max, TH1 *f1, TH1 *f2, TH1 *f3, TH1 *f4, TH1 *f5, TH1 *f6, TH1 *f7)
{
    if(f1==NULL) return;

    SetHistTheme(f1, f2, f3, f4, f5, f6, f7);

    if(gH1!=NULL) delete gH1;
    gH1 = new TH1F(*(TH1F*)f1);

    gH1->SetXTitle(Label[0]);
    gH1->SetYTitle(Label[1]);
    gH1->SetTitle(Title);

    if(min!=max && min<max) gH1->GetYaxis()->SetRangeUser(min, max);

    gH1->Draw(Option);
    if(f2!=NULL) f2->Draw(Option+TString("same"));
    if(f3!=NULL) f3->Draw(Option+TString("same"));
    if(f4!=NULL) f4->Draw(Option+TString("same"));
    if(f5!=NULL) f5->Draw(Option+TString("same"));
    if(f6!=NULL) f6->Draw(Option+TString("same"));
    if(f7!=NULL) f7->Draw(Option+TString("same"));
    if(gH1!=NULL) gH1->Draw(Option+TString("same"));

    if(f2==NULL) return;
    TLegend *leg = new TLegend(0.8, 0.7, 1.0, 0.8);
    leg->AddEntry(gH1, Legends[0], Option+TString("l"));
    if(f2!=NULL) leg->AddEntry(f2, Legends[1], Option+TString("l"));
    if(f3!=NULL) leg->AddEntry(f3, Legends[2], Option+TString("l"));
    if(f4!=NULL) leg->AddEntry(f4, Legends[3], Option+TString("l"));
    if(f5!=NULL) leg->AddEntry(f5, Legends[4], Option+TString("l"));
    if(f6!=NULL) leg->AddEntry(f6, Legends[5], Option+TString("l"));
    if(f7!=NULL) leg->AddEntry(f7, Legends[6], Option+TString("l"));
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw("same");
}

void MyStyle::DrawHistograms(double min, double max, TH2 *f1, TH2 *f2, TH2 *f3, TH2 *f4, TH2 *f5, TH2 *f6, TH2 *f7)
{
    if(f1==NULL) return;

    SetHistTheme(f1, f2, f3, f4, f5, f6, f7);

    if(gH2!=NULL) delete gH2;
    gH2 = new TH2F(*(TH2F*)f1);

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

    if(min!=max && min<max) {
        gH2->GetXaxis()->SetRangeUser(min, max);
        gH2->GetYaxis()->SetRangeUser(min, max);
    }

    gH2->Draw(Option);
    if(f2!=NULL) f2->Draw(Option+TString("same"));
    if(f3!=NULL) f3->Draw(Option+TString("same"));
    if(f4!=NULL) f4->Draw(Option+TString("same"));
    if(f5!=NULL) f5->Draw(Option+TString("same"));
    if(f6!=NULL) f6->Draw(Option+TString("same"));
    if(f7!=NULL) f7->Draw(Option+TString("same"));
    if(gH2!=NULL) gH2->Draw(Option+TString("same"));

    if(f2==NULL) return;
    TLegend *leg = new TLegend(0.8, 0.7, 1.0, 0.8);
    leg->AddEntry(gH2, Legends[0], Option+TString("p"));
    if(f2!=NULL) leg->AddEntry(f2, Legends[1], Option+TString("p"));
    if(f3!=NULL) leg->AddEntry(f3, Legends[2], Option+TString("p"));
    if(f4!=NULL) leg->AddEntry(f4, Legends[3], Option+TString("p"));
    if(f5!=NULL) leg->AddEntry(f5, Legends[4], Option+TString("p"));
    if(f6!=NULL) leg->AddEntry(f6, Legends[5], Option+TString("p"));
    if(f7!=NULL) leg->AddEntry(f7, Legends[6], Option+TString("p"));
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw("same");
}

void MyStyle::DrawGraphs(TGraph *g1, TGraph *g2, TGraph*g3, TGraph*g4, TGraph*g5, TGraph*g6, TGraph*g7)
{
    if(g1==NULL) return;

    UpdateTheme();

    if(gG!=NULL) delete gG;
    gG = new TGraph(*g1);

    if(CP==MATHEMATIC_STYLE) {
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

    int nG = 1;

    SetTheme(gG, 0);
    if(g2!=NULL) { nG++; SetTheme(g2, 1); } 
    if(g3!=NULL) { nG++; SetTheme(g3, 2); }
    if(g4!=NULL) { nG++; SetTheme(g4, 3); }
    if(g5!=NULL) { nG++; SetTheme(g5, 4); }
    if(g6!=NULL) { nG++; SetTheme(g6, 5); }
    if(g7!=NULL) { nG++; SetTheme(g7, 6); }

    double ymin = gG->GetYaxis()->GetXmin();
    double ymax = gG->GetYaxis()->GetXmax();
    if(g2!=NULL) { ymin = (g2->GetYaxis()->GetXmin()>ymin) ? ymin : g2->GetYaxis()->GetXmin();  ymax = (g2->GetYaxis()->GetXmax()<ymax) ? ymax : g2->GetYaxis()->GetXmax(); }
    if(g3!=NULL) { ymin = (g3->GetYaxis()->GetXmin()>ymin) ? ymin : g3->GetYaxis()->GetXmin();  ymax = (g3->GetYaxis()->GetXmax()<ymax) ? ymax : g3->GetYaxis()->GetXmax(); }
    if(g4!=NULL) { ymin = (g4->GetYaxis()->GetXmin()>ymin) ? ymin : g4->GetYaxis()->GetXmin();  ymax = (g4->GetYaxis()->GetXmax()<ymax) ? ymax : g4->GetYaxis()->GetXmax(); }
    if(g5!=NULL) { ymin = (g5->GetYaxis()->GetXmin()>ymin) ? ymin : g5->GetYaxis()->GetXmin();  ymax = (g5->GetYaxis()->GetXmax()<ymax) ? ymax : g5->GetYaxis()->GetXmax(); }
    if(g6!=NULL) { ymin = (g6->GetYaxis()->GetXmin()>ymin) ? ymin : g6->GetYaxis()->GetXmin();  ymax = (g6->GetYaxis()->GetXmax()<ymax) ? ymax : g6->GetYaxis()->GetXmax(); }
    if(g7!=NULL) { ymin = (g7->GetYaxis()->GetXmin()>ymin) ? ymin : g7->GetYaxis()->GetXmin();  ymax = (g7->GetYaxis()->GetXmax()<ymax) ? ymax : g7->GetYaxis()->GetXmax(); }

    gG->SetTitle(Title);
    gG->GetXaxis()->SetTitle(Label[0]);
    gG->GetYaxis()->SetTitle(Label[1]);
    gG->GetYaxis()->SetRangeUser(ymin, ymax);

    gG->Draw(Option+"a");
    if(g2!=NULL) g2->Draw(Option+TString("same"));
    if(g3!=NULL) g3->Draw(Option+TString("same"));
    if(g4!=NULL) g4->Draw(Option+TString("same"));
    if(g5!=NULL) g5->Draw(Option+TString("same"));
    if(g6!=NULL) g6->Draw(Option+TString("same"));
    if(g7!=NULL) g7->Draw(Option+TString("same"));

    TLegend *leg = new TLegend(0.8, 0.6 , 1.0, 0.6 +nG*0.05);
    leg->SetFillColor(0);
    leg->AddEntry(gG, Legends[0], Option+TString("l"));
    if(g2!=NULL) leg->AddEntry(g2, Legends[1], Option+TString("l"));
    if(g3!=NULL) leg->AddEntry(g3, Legends[2], Option+TString("l"));
    if(g4!=NULL) leg->AddEntry(g4, Legends[3], Option+TString("l"));
    if(g5!=NULL) leg->AddEntry(g5, Legends[4], Option+TString("l"));
    if(g6!=NULL) leg->AddEntry(g6, Legends[5], Option+TString("l"));
    if(g7!=NULL) leg->AddEntry(g7, Legends[6], Option+TString("l"));

    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw("same");
}

void MyStyle::DrawGraphs(vector<TGraph *> gvec)
{
    int ng = gvec.size();
    if(ng==0) return;

    UpdateTheme();

    double ymin = gvec[0]->GetYaxis()->GetXmin();
    double ymax = gvec[0]->GetYaxis()->GetXmax();

    for(int i=0; i<(int)ng; i++)
    {
        SetTheme(gvec[i], i);
        ymin = (gvec[i]->GetYaxis()->GetXmin()>ymin) ? ymin : gvec[i]->GetYaxis()->GetXmin(); 
        ymax = (gvec[i]->GetYaxis()->GetXmax()<ymax) ? ymax : gvec[i]->GetYaxis()->GetXmax();
    }

    gvec[0]->GetYaxis()->SetRangeUser(ymin, ymax);
    for(int i=0; i<(int)ng; i++)
    {
        if(i==0) 
            gvec[i]->Draw(Option+"a");
        else 
            gvec[i]->Draw(Option+TString("same"));
    }

    TLegend *leg = new TLegend(0.8, 0.6 , 1.0, 0.6 + ng*0.05);
    leg->SetFillColor(0);
    for(int i=0; i<(int)ng; i++)
        leg->AddEntry(gvec[i], Legends[i], Option+TString("l"));

    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw("same");
}



void MyStyle::UpdateTheme()
{
    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(0);

    if(cID<=2000)
    {
        new TColor(cID++, 0x99/255., 0x99/255., 0x99/255.);
        new TColor(cID++, 0x6e/255., 0x6e/255., 0x6e/255.);
        new TColor(cID++, 0x8f/255., 0xb0/255., 0x31/255.);
        new TColor(cID++, 0xe2/255., 0x9e/255., 0x29/255.);
        new TColor(cID++, 0x5e/255., 0x81/255., 0xb5/255.);
        new TColor(cID++, 0xeb/255., 0x63/255., 0x35/255.);
        new TColor(cID++, 0x8b/255., 0x7c/255., 0x76/255.);
    }


    double lcolor1[__NHIST__] = {2002, 2003, 2004, 2005, 2006, 2000, 2000, 2000, 2000, 2000};
    double lstyle1[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double lwidth1[__NHIST__] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    double mstyle1[__NHIST__] = {2, 3, 5, 8,29,13,14,15,16,18};
    double msize1 [__NHIST__] = {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4};
    double fstyle1[__NHIST__] = {1,2,3,4,5,6,7,8,9,10};

    double lcolor2[__NHIST__] = {kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack};
    double lstyle2[__NHIST__] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double lwidth2[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double mstyle2[__NHIST__] = {2, 3, 5, 7,29,13,14,15,16,18};
    double msize2 [__NHIST__] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double fstyle2[__NHIST__] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//{3004,3005,3006,3007,3010,3013,3001,3002,3003,3021};

    double lcolor3[__NHIST__] = {kRed, kBlue, kOrange, kViolet, kTeal, kAzure, kMagenta, kYellow, kPink, kSpring};  //kCyan, kTeal,
    double lstyle3[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double lwidth3[__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double mstyle3[__NHIST__] = {2, 3, 5, 7,29,13,14,15,16,18};
    double msize3 [__NHIST__] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double fstyle3[__NHIST__] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
 
    switch(CP)
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
    for(int i=0; i<__NHIST__; i++)
    {
        LineColor[i] = lcolor[i];
        LineStyle[i] = lstyle[i];
        LineWidth[i] = lwidth[i];

        MarkerColor[i] = mcolor[i];
        MarkerStyle[i] = mstyle[i];
        MarkerSize [i] = msize [i];

        FillColor[i] = fcolor[i];
        FillStyle[i] = fstyle[i];
    }
}


void MyStyle::SetTheme(TH1 *hist, int id)
{
   hist->SetLineColor(LineColor[id]);
   hist->SetLineStyle(LineStyle[id]);
   hist->SetLineWidth(LineWidth[id]);

   hist->SetMarkerColor(MarkerColor[id]);
   hist->SetMarkerStyle(MarkerStyle[id]);
   hist->SetMarkerSize (MarkerSize [id]);

   hist->SetFillColor(FillColor[id]);
   hist->SetFillStyle(FillStyle[id]);
}

void MyStyle::SetTheme(TH2 *hist, int id)
{
   hist->SetLineColor(LineColor[id]);
   hist->SetLineStyle(LineStyle[id]);
   hist->SetLineWidth(LineWidth[id]);

   hist->SetMarkerColor(MarkerColor[id]);
   hist->SetMarkerStyle(MarkerStyle[id]);
   hist->SetMarkerSize (MarkerSize [id]);

   hist->SetFillColor(FillColor[id]);
   hist->SetFillStyle(FillStyle[id]);
}

void MyStyle::SetTheme(TGraph *graph, int id)
{
   graph->SetLineColor(LineColor[id]);
   graph->SetLineStyle(LineStyle[id]);
   graph->SetLineWidth(LineWidth[id]);

   graph->SetMarkerColor(MarkerColor[id]);
   graph->SetMarkerStyle(MarkerStyle[id]);
   graph->SetMarkerSize (MarkerSize [id]);

   graph->SetFillColor(FillColor[id]);
   graph->SetFillStyle(FillStyle[id]);
}
