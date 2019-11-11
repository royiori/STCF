TCanvas c1;
TFile f("res.root");

#include "MyStyle.h"
#include "MyStyle.cpp"

MyStyle *myStyle = 0;

void Draw(int mom, int the)
{
if(myStyle==0)
  myStyle = new MyStyle();
myStyle->SetColorPattern(ROOT_STYLE);

c1.Clear();
vector<TH1 *>vec;
for(int i=0; i<5; i++)
{
  TH1F *ftmp1 = (TH1F*)gDirectory->Get(Form("Rec%d_%d_%d_0", i, mom, the));
  TH1F *ftmp2 = (TH1F*)gDirectory->Get(Form("Rec%d_%d_%d_1", i, mom, the));
  ftmp1->Add(ftmp2);
  vec.push_back(ftmp1);
}
myStyle->SetDrawOption("bar2");
myStyle->Draw1DHistVector(vec);
c1.Modified();
c1.Update();
}
