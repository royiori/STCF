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

TH1D *MakeHistogram(double x0, double sig, TString hname);
RooDataSet *makeFakeData(double x0, double sigma, int nevent);
RooDataSet *shiftData(RooDataSet *data);

void test()
{
   RooRealVar x("x", "x", -50, 50);

/*
   // C r e a t e   Data
   // ---------------------------------------------------
   RooDataSet *data = makeFakeData(0, 3, 1000);

   RooDataSet *data2 = shiftData(data);

   // ---------------------------------------------------
   RooPlot *frame1 = x.frame(Title("2D histogram pdf"));
   data->plotOn(frame1, MarkerColor(kRed));
   data2->plotOn(frame1, MarkerColor(kBlue));
   TCanvas *c = new TCanvas();
   frame1->Draw();
*/
   TH1D *fRing = MakeHistogram(0, 10, "hh");

   RooRealVar X0("X0", "X0", 1.0,  0, 5);
   RooFormulaVar x1("x1", "x - X0", RooArgSet(x, X0));
   RooDataHist hist1("hist1", "hist1", x, Import(*fRing));
   RooHistPdf histpdf1("histpdf1", "histpdf1", x1, x, hist1, 2);

   RooDataSet *data = makeFakeData(2, 10, 10000);
   RooFitResult* res = histpdf1.fitTo(*data, Save());

   RooPlot *frame1 = x.frame(Title("2D histogram pdf"));
   data->plotOn(frame1);
   histpdf1.plotOn(frame1);
   TCanvas *c = new TCanvas();
   frame1->Draw();
}

TH1D *MakeHistogram(double x0, double sig, TString hname)
{
   TH1D *hh = new TH1D(hname, "hh", 80, -50, 50);
   for (int i = 0; i < 10000000; i++)
   {
      double x = gRandom->Gaus(x0, sig);
      hh->Fill(x);
   }

   return hh;
}


RooDataSet *makeFakeData(double x0, double sigma, int nevent)
{
   RooRealVar x("x", "x", -50, 50);

   RooDataSet *d = new RooDataSet("d", "d", x);

   for (int i = 0; i < nevent; i++)
   {
      x = gRandom->Gaus(x0, sigma);
      d->add(x);
   }
   return d;
}

RooDataSet *shiftData(RooDataSet *data)
{
   RooRealVar x("x", "x", -1000, 1000);  
    
   RooDataSet *d = new RooDataSet("d", "d", x);

   //TTree *tree = data->GetClonedTree();
   //double x0;
   //tree->SetBranchAddress("x", &x0);

   //for(long i=0; i<tree->GetEntries(); i++)
   //{
   //    //tree->GetEntry(i);
   //    x = x0 * 100 + 20;
   //    d->add(x);
   //}
   int entries = data->numEntries();
    for (int i=0; i<entries; i++) {
        double x0 = data->get(i)[0].getRealValue("x");
        cout<<x0<<endl;
        x = 100 * x0 + 20;
        d->add(x);
    }

   return d;
}
