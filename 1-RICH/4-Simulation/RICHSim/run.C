//------- const values ------------
const int Npar = 3;
char particle[Npar][10] = {"pi+", "kaon+", "proton"};
char partFileName[Npar][10] = {"pip_", "kp_", "p_"};

int NEvent = 1000;

double pmin = 2.5;
double pmax = 2.5;
double pstp = 0.1; // GeV

double cos_min = 0.5;
double cos_max = 1.0;
double cos_stp = 0.1; //cos

double theta_min = 48.0;
double theta_max = 48.0;
double theta_stp = 1.; //theta

void genBatch();
void anaRoot();
//
//--
void run()
{
   genBatch();
   //anaRoot();
}

//
//--
void genBatch()
{
   gSystem->Exec("mkdir ./batch");

   ofstream setfile;
   setfile.open("./batch/setting.txt");
   setfile << NEvent << "\n";
   setfile << pmin << " " << pmax << " " << pstp << "\n";
   setfile << theta_min << " " << theta_max << " " << theta_stp << "\n";
   setfile.close();

   ofstream runfile;
   runfile.open("./batch.sh");

   for (int i = 0; i < (int)((pmax - pmin) / pstp + 1); i++)
   {
      //for (int j = 0; j < (int)((cos_max - cos_min) / cos_stp + 1); j++)
      for (int j = 0; j < (int)((theta_max - theta_min) / theta_stp + 1); j++)
      {
         for (int k = 0; k < Npar; k++)
         {
            double p = pmin + pstp * i;
            //double cos = cos_min + cos_stp * j;
            double theta0 = theta_min + theta_stp * j;
            double theta1 = (theta0) / 180. * TMath::Pi();

            ofstream outfile;
            //TString filename1(Form("./batch/%sp%.2fcos%.2f.txt", partFileName[k], p, cos));
            //TString filename2(Form("./batch/%sp%.2fcos%.2f.root", partFileName[k], p, cos));
            TString filename1(Form("./batch/%sp%.2ftheta%.2f.txt", partFileName[k], p, theta0));
            TString filename2(Form("./batch/%sp%.2ftheta%.2f.root", partFileName[k], p, theta0));
            outfile.open(filename1);
            outfile << "#\n";
            outfile << "/control/verbose 0\n";
            outfile << "/run/verbose 0\n";
            outfile << "/process/verbose 0\n";
            outfile << "/tracking/verbose 0\n";
            outfile << "#\n";
            outfile << "/myrun/filename " << filename2 << "\n";
            outfile << "/globalField/setValue 0 0 1 tesla\n";
            outfile << "#\n";
            outfile << "/run/initialize\n";
            outfile << "#\n";
            outfile << "/gun/particle " << particle[k] << "\n";
            //outfile << "/gun/direction " << cos << " 0 " << sqrt(1 - cos * cos) << "\n";
            outfile << "/gun/direction " << cos(theta1) << " 0 " << sin(theta1) << "\n";
            outfile << "/gun/energy " << p << " GeV\n";
            outfile << "/run/beamOn " << NEvent << "\n";
            outfile.close();
            //gSystem->Exec(Form("./gdmlDet %s", filename.Data()));
            runfile << "./gdmlDet " << filename1 << endl;
         }
      }
   }
   runfile.close();
   gSystem->Exec("chmod 777 ./batch.sh");
}

//
//--
double getMeanFromTree(TTree *tree, TString lname)
{
   TH1F *htemp = (TH1F *)gDirectory->Get("htemp");
   if(htemp!=0)
      delete htemp;
   tree->Draw(lname+">>htemp");
   htemp = (TH1F *)gDirectory->Get("htemp");
   return htemp->GetMean();
}

double Reconstruct(double x0, double y0, double z0, double px, double py, double pz, double x, double y, double z)
{
   double 
}

void anaRoot()
{
   gStyle->SetOptFit(1);
   gStyle->SetDrawOption("ap");
   TF1 *func = new TF1("func", "[0]*(x-[1])*(x-[1])+[2]", -200, 200);

   TFile *rootFile = new TFile("result.root", "recreate");

   for (int i = 0; i < (int)((pmax - pmin) / pstp + 1); i++)
   {
      //for (int j = 0; j < (int)((cos_max - cos_min) / cos_stp + 1); j++)
      for (int j = 0; j < (int)((theta_max - theta_min) / theta_stp + 1); j++)
      {
         for (int k = 0; k < Npar; k++)
         {
            //-- open file
            double p = pmin + pstp * i;
            double theta0 = theta_min + theta_stp * j;
            double theta1 = (theta0) / 180. * TMath::Pi();

            TString filename2(Form("./batch/%sp%.2ftheta%.2f.root", partFileName[k], p, theta0));
            TFile rootfp(filename2);
            if (rootfp.IsOpen() != true)
               continue;
            cout << "opening: " << filename2 << endl;

            //-- define constant
            double x0, y0, z0;
            double px, py, pz;            
            double X, Y, Z;

            x0 = getMeanFromTree((TTree *)rootfp.Get("Init"), "X");
            y0 = getMeanFromTree((TTree *)rootfp.Get("Init"), "Y");
            z0 = getMeanFromTree((TTree *)rootfp.Get("Init"), "Z");
            px = getMeanFromTree((TTree *)rootfp.Get("Init"), "PX");
            py = getMeanFromTree((TTree *)rootfp.Get("Init"), "PY");
            pz = getMeanFromTree((TTree *)rootfp.Get("Init"), "PZ");
            cout<<"---->"<<x0<<" "<<y0<<" "<<z0<<" "<<px<<" "<<py<<" "<<pz<<endl;

            TTree *tree = (TTree *)rootfp.Get("SD");
            tree->SetBranchAddress("X", &X);
            tree->SetBranchAddress("Y", &Y);
            tree->SetBranchAddress("Z", &Z);

            //--- define histogram
            TH2F *fHitmap = new TH2F(Form("h2_%sp%.2f_theta%.2f", partFileName[k], p, theta0),
                                     Form("yzHist %s p=%.2f theta=%.2f", partFileName[k], p, theta0),
                                     80, -200, 200, 120, 800, 1400);

            Long64_t nentries = tree->GetEntries();
            for (Long64_t i = 0; i < nentries; i++)
            {
               tree->GetEntry(i);
               fHitmap->Fill(Y, Z);

               //double rec0 = Reconstruct(Y, Z, 0);
               //double rec0 = Reconstruct(Y, Z, 1);
               //double rec0 = Reconstruct(Y, Z, 2);
            }

            fHitmap->Draw("colz");

            rootFile->cd();
            fHitmap->Write();
         }
      }
   }
}