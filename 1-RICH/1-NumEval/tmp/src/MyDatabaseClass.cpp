#include "math.h"
#include "MyDatabaseClass.h"

//______________________________________________________________________________
//
MyDatabaseClass::MyDatabaseClass()
{
    MaterialList.push_back("C6F14");
    MaterialList.push_back("Quartz");
    MaterialList.push_back("Quartz311");
    MaterialList.push_back("Quartz3001");
    MaterialList.push_back("Quartz3301");
    MaterialList.push_back("MgF2");
    MaterialList.push_back("CaF2");
    MaterialList.push_back("Water");
    MaterialList.push_back("Aerogel105");
    MaterialList.push_back("Aerogel106");
    MaterialList.push_back("Aerogel107");
    MaterialList.push_back("Aerogel108");
    MaterialList.push_back("Aerogel109");
    MaterialList.push_back("Aerogel110");
    MaterialList.push_back("Ar");
    MaterialList.push_back("N2");
    MaterialList.push_back("CH4");
    MaterialList.push_back("H2O");
    MaterialList.push_back("O2");
    nMaterial = (int)MaterialList.size();
    MatAbsGraph.resize(nMaterial);
    MatTrsGraph.resize(nMaterial);
    MatRefGraph.resize(nMaterial);

    DetectorList.push_back("CsI");
    DetectorList.push_back("MPPC1025");
    DetectorList.push_back("MPPC1050");
    DetectorList.push_back("MPPC2050");
    DetectorList.push_back("APD");
    nDetector = (int)DetectorList.size();
    DetQEGraph.resize(nDetector);

    // read data
    TString datpath = gSystem->WorkingDirectory() + TString("/database/");
    for (int i = 0; i < nMaterial; i++)
        MatAbsData.push_back(new MyReadData(datpath + MaterialList[i] + "_abs.txt", MaterialList[i] + "_abs"));
    for (int i = 0; i < nMaterial; i++)
        MatRefData.push_back(new MyReadData(datpath + MaterialList[i] + "_ref.txt", MaterialList[i] + "_ref"));
    for (int i = 0; i < nDetector; i++)
        DetQEData.push_back(new MyReadData(datpath + DetectorList[i] + "_qe.txt", DetectorList[i] + "_qe"));
}

//--------------------------------
// get absorption length for each material [in cm]
//
double MyDatabaseClass::GetMatAbsValue(TString matStr, double lambda, double ppm)
{
    int id = GetID(matStr, MaterialList);
    if (id == -1)
        return 0.;
    if (ppm == 0)
        return MatAbsData.at(id)->GetValue(lambda);

    // ppm!=0 so it's cross section, use ppm to calcuate the absorption length
    double A = 1;
    A = (matStr == "H2O") ? 18 : A;
    A = (matStr == "O2") ? 32 : A;
    double Xsec = MatAbsData.at(id)->GetValue(lambda);
    double xsec = Xsec * 1E6 * 1E-28; //m^2
    double density = ppm / A * 6.02214129 * 1E23;
    double absleng = 100. / 1 / xsec / density; // /sqrt(2)  //[cm]
    return (xsec == 0 || density == 0) ? 10000 : absleng;
}

TGraph *MyDatabaseClass::GetMatAbsGraph(TString matStr, int nL, double lmin, double lmax, double ppm)
{
    int id = GetID(matStr, MaterialList);
    if (id == -1 || nL == 0 || lmax <= lmin)
        return 0;
    if (MatAbsGraph[id] != NULL)
        delete MatAbsGraph[id];
    MatAbsGraph[id] = new TGraph();

    double lstep = (lmax - lmin) / nL;
    for (int i = 0; i < nL; i++)
    {
        double lambda = lmin + (i + 0.5) * lstep;
        MatAbsGraph[id]->SetPoint(i, lambda, GetMatAbsValue(matStr, lambda, ppm));
    }
    MatAbsGraph[id]->SetTitle("Absorption length of " + matStr);
    MatAbsGraph[id]->GetXaxis()->SetTitle("lambda");
    MatAbsGraph[id]->GetYaxis()->SetTitle("Absorption length[cm]");
    return MatAbsGraph[id];
}


TGraph *MyDatabaseClass::GetMatTrsGraph(TString matStr, int nL, double lmin, double lmax, double ppm)
{
    int id = GetID(matStr, MaterialList);
    if (id == -1 || nL == 0 || lmax <= lmin)
        return 0;
    if (MatTrsGraph[id] != NULL)
        delete MatTrsGraph[id];
    MatTrsGraph[id] = new TGraph();

    double lstep = (lmax - lmin) / nL;
    for (int i = 0; i < nL; i++)
    {
        double lambda = lmin + (i + 0.5) * lstep;
        double abslen = GetMatAbsValue(matStr, lambda, ppm);
        MatTrsGraph[id]->SetPoint(i, lambda, exp(-1./abslen)); //[cm]
    }
    MatTrsGraph[id]->SetTitle("Transmission for 10mm " + matStr);
    MatTrsGraph[id]->GetXaxis()->SetTitle("lambda");
    MatTrsGraph[id]->GetYaxis()->SetTitle("Transmission[%]");
    return MatTrsGraph[id];
}
//--------------------------------
// get refractive index for each material [in cm]
//
double MyDatabaseClass::GetMatRefValue(TString matStr, double lambda)
{
    int id  = GetID(matStr, MaterialList);
    return (id == -1) ? 1 : MatRefData.at(id)->GetValue(lambda);
}

TGraph *MyDatabaseClass::GetMatRefGraph(TString matStr, int nL, double lmin, double lmax)
{
    int id = GetID(matStr, MaterialList);
    if (id == -1 || nL == 0 || lmax <= lmin)
        return 0;
    if (MatRefGraph[id] != NULL)
        delete MatRefGraph[id];
    MatRefGraph[id] = new TGraph();

    double lstep = (lmax - lmin) / nL;
    for (int i = 0; i < nL; i++)
    {
        double lambda = lmin + (i + 0.5) * lstep;
        MatRefGraph[id]->SetPoint(i, lambda, GetMatRefValue(matStr, lambda));
    }
    MatRefGraph[id]->SetTitle("Reflective index of " + matStr);
    MatRefGraph[id]->GetXaxis()->SetTitle("lambda");
    MatRefGraph[id]->GetYaxis()->SetTitle("Reflective index");
    return MatRefGraph[id];
}

//--------------------------------
// get Q.E. for each detector
//
double MyDatabaseClass::GetDetQEValue(TString detStr, double lambda)
{
    int id = GetID(detStr, DetectorList);
    return (id == -1) ? 0 : DetQEData.at(id)->GetValue(lambda) / 100.;
}

TGraph *MyDatabaseClass::GetDetQEGraph(TString detStr, int nL, double lmin, double lmax)
{
    int id = GetID(detStr, DetectorList);
    if (id == -1 || nL == 0 || lmax <= lmin)
        return 0;
    if (DetQEGraph[id] != NULL)
        delete DetQEGraph[id];
    DetQEGraph[id] = new TGraph();

    double lstep = (lmax - lmin) / nL;
    for (int i = 0; i < nL; i++)
    {
        double lambda = lmin + (i + 0.5) * lstep;
        DetQEGraph[id]->SetPoint(i, lambda, GetDetQEValue(detStr, lambda));
    }
    DetQEGraph[id]->SetTitle("Quantum efficiency of " + detStr);
    DetQEGraph[id]->GetXaxis()->SetTitle("lambda");
    DetQEGraph[id]->GetYaxis()->SetTitle("Q.E.[%]");
    return DetQEGraph[id];
}