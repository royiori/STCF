/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
///  'SPECIAL PDFS' RooFit tutorial macro #706
///
///  Histogram based p.d.f.s and functions
///
/// \macro_image
/// \macro_output
/// \macro_code
/// \author 07/2008 - Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TRandom.h"
#include "TH2D.h"
#include "RooPlot.h"
using namespace RooFit;

TH2D *MakeRICHRing(double x0, double y0, double R, TString hname);
RooDataSet *makeFakeDataXY(double x0, double y0, double R, int nevent);

void RICHFit()
{
   // C r e a t e   PDF
   // ---------------------------------------------------
   TH2D *fRing = MakeRICHRing(0, 0, 140, "hh");
   TCanvas *cc = new TCanvas();
   fRing->Draw("colz");
   return;

   fRing->SetLineColor(kBlue);

   RooRealVar x("x", "x", -10, 10);
   RooRealVar y("y", "y", -10, 10);
   RooRealVar X0("X0", "X0", 5.0,  4, 6);
   RooRealVar Y0("Y0", "Y0", 5.0,  4, 6);
   RooFormulaVar x1("x1", "x - X0", RooArgSet(x, X0));
   RooFormulaVar y1("y1", "y - Y0", RooArgSet(y, Y0));

   RooDataHist hist1("hist1", "hist1", RooArgSet(x, y), Import(*fRing));
   RooHistPdf histpdf1("histpdf1", "histpdf1", RooArgSet(x1, y1), RooArgSet(x, y), hist1, 2);

   // background
   RooUniform bkpdf("bkpdf", "bkpdf", RooArgSet(x1, y1));
   RooRealVar nsig("nsig", "number of signal events in signalRange", 10, 9, 111);
   RooRealVar nbkg("nbkg", "number of background events", 0.1, 0, 0.2);
   RooExtendPdf esig("esig", "extended signal p.d.f", histpdf1, nsig);
   RooExtendPdf ebkg("ebkg", "extended background p.d.f", bkpdf, nbkg);
   RooAddPdf model("model", "model", RooArgList(ebkg, esig));

   // C r e a t e   Data
   // ---------------------------------------------------
   // Obtain fake external experimental dataset with values for x and y
   RooDataSet *data = makeFakeDataXY(5, 5, 4, 10);

   RooFitResult* res = model.fitTo(*data, Extended(kTRUE), Save());
   
   // Plot unbinned data and histogram pdf overlaid
   // ---------------------------------------------------
   RooPlot *frame1 = x.frame(Title("2D histogram pdf"));
   data->plotOn(frame1);
   model.plotOn(frame1);//, VisualizeError(*res));

   TH1* hh = model.createHistogram("hh",x,Binning(50),YVar(y,Binning(50)));
   hh->SetLineColor(kRed);

   TCanvas *c = new TCanvas();
   //c->Divide(1, 2);
   //c->cd(1);
   //frame1->GetYaxis()->SetTitleOffset(1.4);
   //frame1->Draw();
   //cout<<frame1->chiSquare()<<endl;
    //hfit->Draw("colz");

   //c->cd(2);
   hh->Draw("colz");
   
   cout<<res->edm()<<" "<<res->minNll()<<endl;
}

TH2D *MakeRICHRing(double x0, double y0, double R, TString hname)
{
   TH2D *hh = new TH2D(hname, "hh", 80, -200, 200, 80, -200, 200);
   for (int i = 0; i < 10000000; i++)
   {
      double phi = gRandom->Uniform(0, 2*TMath::Pi());
      double x = R*sin(phi) + gRandom->Gaus(0, 5);
      double y = R*cos(phi) + gRandom->Gaus(0, 5);
      hh->Fill(x - x0, y - y0);
   }

   //for (int i = 0; i < 100000; i++)
   //{
   //   double x = gRandom->Uniform(-200, 200);
   //   double y = gRandom->Uniform(-200, 200);
   //   hh->Fill(x - x0, y - y0);
   //}

   return hh;
}

RooDataSet *makeFakeDataXY(double x0, double y0, double R, int nevent)
{
   RooRealVar x("x", "x", -200, 200);
   RooRealVar y("y", "y", -200, 200);
   RooArgSet coord(x, y);

   RooDataSet *d = new RooDataSet("d", "d", RooArgSet(x, y));

   TH2D *f = MakeRICHRing(x0, y0, R, "htemp");
   for (int i = 0; i < nevent; i++)
   {
      Double_t tmpy;
      Double_t tmpx;
      f->GetRandom2(tmpx, tmpy);
      x = tmpx;
      y = tmpy;
      d->add(coord);
   }
   return d;
}
