#include "MyMainFrameGui.h"
#include "MyNumEvalClass.h"

MyNumEvalClass *gMyNumEvalClass = (MyNumEvalClass *)0;

//______________________________________________________________________________
//
MyNumEvalClass::MyNumEvalClass()
{
    env = new TEnv(gSystem->WorkingDirectory() + TString("/.env"));
    env->SaveLevel(kEnvLocal);

    gMyStyle = new MyStyle();
    gMyStyle->SetDrawOption("c");
    gMyStyle->SetColorPattern(MATHEMATIC_STYLE);

    gMyCommonRICH = new MyCommonRICH();

    gMyDatabaseClass = new MyDatabaseClass();

    // pages
    nTabPage = 4;
    sTabPage[0] = TString("Act");
    sTabPage[1] = TString("Mat");
    sTabPage[2] = TString("Det");
    sTabPage[3] = TString("Mis");

    // page 0
    nButton[0] = 4;
    iButton[0] = DrawRICH;
    sButton[0].push_back("Draw selected C-Ring");
    sButton[0].push_back("Draw mu/pi/K/p C-Rings");
    sButton[0].push_back("Save mu/pi/K/p C-Rings");
    sButton[0].push_back("Load mu/pi/K/p C-Rings");

    // page1
    nButton[1] = gMyDatabaseClass->GetNMaterial();
    iButton[1] = MatList;
    sButton[1] = gMyDatabaseClass->GetMaterialList();

    // page 2
    nButton[2] = gMyDatabaseClass->GetNDetector();
    iButton[2] = DetList;
    sButton[2] = gMyDatabaseClass->GetDetectorList();

    // page 3
    nButton[3] = 2;
    iButton[3] = DrawSelectedMat;
    sButton[3].push_back("Draw Sel-Mat");
    sButton[3].push_back("Draw Sel-Det");

    // init parameters
    nRadLayer = env->GetValue("nRadLayer", 2);
    for (int i = 0; i < nRadLayer; i++)
    {
        sRadLayer.push_back(env->GetValue(Form("sRad%d", i), (i == 0) ? "C6F14" : "Quartz"));
        tRadLayer.push_back(env->GetValue(Form("tRad%d", i), (i == 0) ? 10 : 4));
    }

    tTransLayer = env->GetValue("tTransLayer", 100);
    nTransLayer = env->GetValue("nTransLayer", 2);
    for (int i = 0; i < nTransLayer; i++)
    {
        sTransLayer.push_back(env->GetValue(Form("sTrans%d", i), (i == 0) ? "Ar" : "CH4"));
        pTransLayer.push_back(env->GetValue(Form("tTrans%d", i), 50));
    }

    nImpurities = env->GetValue("nImpurities", 2);
    for (int i = 0; i < nImpurities; i++)
    {
        sImpurities.push_back(env->GetValue(Form("sImp%d", i), (i == 0) ? "H2O" : "O2"));
        pImpurities.push_back(env->GetValue(Form("tImp%d", i), 1));
    }

    PCDetector = env->GetValue("PCDetector", "CsI");
    pixel = env->GetValue("pixel", 5);

    nLambda = env->GetValue("nLambda", 60);
    lambdaMin = env->GetValue("lambdaMin", 160);
    lambdaMax = env->GetValue("lambdaMax", 220);

    particle = env->GetValue("particle", "pion");
    momentum = env->GetValue("momentum", 2.);
    theta0 = env->GetValue("theta0", 0.);

    XBinMin = env->GetValue("XBinMin", -200);
    XBinMax = env->GetValue("XBinMax", 200);
    YBinMin = env->GetValue("YBinMin", -200);
    YBinMax = env->GetValue("YBinMax", 200);
    NXBin = XBinMax - XBinMin;
    NYBin = YBinMax - YBinMin;

    epsilon = env->GetValue("epsilon", 1e-5);
    trkStep = env->GetValue("trkStep", 1e-1);
    nphi = env->GetValue("nphi", 360);

    nSelectedMat = env->GetValue("nSelectedMat", 2);
    for (int i = 0; i < nSelectedMat; i++)
        selectedMat.push_back(env->GetValue(Form("selectedMat%d", i), (i == 0) ? "Quartz" : "MgF2"));

    nSelectedDet = env->GetValue("nSelectedDet", 2);
    for (int i = 0; i < nSelectedDet; i++)
        selectedDet.push_back(env->GetValue(Form("selectedDet%d", i), (i == 0) ? "CsI" : "APD"));
};

MyNumEvalClass::~MyNumEvalClass()
{
    // Destructor.
}

//______________________________________________________________________________
// read / generate settings
TString MyNumEvalClass::GenerateSettingsText()
{
    TString settings("");

    //settings += "\n# The pre-defined material list\nMaterial_List: ";
    //for (int i = 0; i < nButton[0]; i++)
    //    settings += sButton[0].at(i) + " ";

    //settings += "\n\n# The pre-defined detector list\nDetector_List: ";
    //for (int i = 0; i < nButton[1]; i++)
    //    settings += sButton[1].at(i) + " ";
    settings += "\n#---------------------------";
    settings += "\n#    Detector settings";
    settings += "\n#---------------------------";
    settings += "\n\n# The number of radiator layers, thickness in [mm]\n";
    for (int i = 0; i < nRadLayer; i++)
        settings += Form("Radiator %s %.1f\n", sRadLayer[i].c_str(), tRadLayer[i]);

    settings += "\n# The transition layer, content is percentile [%]\n";
    for (int i = 0; i < nTransLayer; i++)
        settings += Form("TransLayer %s %.1f\n", sTransLayer[i].c_str(), pTransLayer[i]);
    settings += "\n# The thickness of the transition layer in [mm]\n";
    settings += Form("TransThick %.1f\n", tTransLayer);

    settings += "\n# The Impurities in the transition layer, content is [ppm]\n";
    for (int i = 0; i < nImpurities; i++)
        settings += Form("Impurities %s %.1f\n", sImpurities[i].c_str(), pImpurities[i]);

    settings += "\n# The photocathode detector, pixel size in [mm]\n";
    settings += Form("PCDetector %s %.1f\n", PCDetector.c_str(), pixel);

    settings += "\n# The n_pixel range for view, [xmin, xmax, ymin, ymax]\n";
    settings += Form("PixelRange %.1f %.1f %.1f %.1f\n", XBinMin, XBinMax, YBinMin, YBinMax);

    settings += "\n# The lambda range for analysis\n";
    settings += Form("Lambda %d %.1f %.1f\n", nLambda, lambdaMin, lambdaMax);

    settings += "\n#---------------------------";
    settings += "\n#    Particle settings";
    settings += "\n#---------------------------";

    settings += "\n\n# Particle (mass[GeV], momentum[GeV], angle[default:rad, or:deg])\n";
    settings += Form("Particle: %s %.1f %.1f\n", particle.Data(), momentum, theta0);

    settings += "\n#---------------------------";
    settings += "\n#    Precision settings";
    settings += "\n#---------------------------";

    settings += "\n\n# Step size for charge particle track[mm]\n";
    settings += Form("trkStep: %f\n", trkStep);

    settings += "\n# Number of Step for phi [0~2pi] * 360\n";
    settings += Form("nPhi: %.1f\n", nphi / 360.);
 
    settings += "\n# Solver precision\n";
    settings += Form("Precision: %f\n", epsilon);

    settings += "\n#---------------------------";
    settings += "\n#    Others";
    settings += "\n#---------------------------";

    settings += "\n\n# Draw the selected material\nSelectedMat: ";
    for (int i = 0; i < nSelectedMat; i++)
        settings += selectedMat.at(i) + " ";

    settings += "\n\n# Draw the selected detector\nSelectedDet: ";
    for (int i = 0; i < nSelectedDet; i++)
        settings += selectedDet.at(i) + " ";

    settings += "\n#---------------------------";
    return settings;
}

vector<TString> ReadContent(TString LINE)
{
    TString line = LINE;
    vector<TString> text;

    line.Remove(TString::kBoth, ' ');
    line.ReplaceAll("  ", " ");
    line.ReplaceAll("  ", " ");
    line.ReplaceAll("  ", " ");

    while (line.Length() > 0 && line.Index(' ') > 0)
    {
        TString head = line;
        line.Remove(0, line.Index(' ') + 1);
        head.Remove(head.Index(' '), head.Length());
        head.Remove(TString::kBoth, ' ');
        if (head.Length() > 0)
            text.push_back(head);
        if (line.Index(' ') == -1)
            text.push_back(line);
    }

    return text;
}

void MyNumEvalClass::ReadSettingsText()
{
    int nline = 0;
    TGLongPosition pos(0, nline);

    tRadLayer.clear();
    sRadLayer.clear();
    pTransLayer.clear();
    sTransLayer.clear();
    mapTransLayers.clear();
    pImpurities.clear();
    sImpurities.clear();
    mapImpurities.clear();

    TGText *text = gMyMainFrameGui->GetText();
    while (text->GetLineLength(nline) != -1)
    {
        if (text->GetChar(pos) != -1)
        {
            TString line(text->GetLine(pos, 100));
            pos.fY = (nline++);

            if (line.BeginsWith("#") || line.Length() < 5)
                continue;

            vector<TString> cont = ReadContent(line);
            cont.erase(cont.begin());

            if (line.BeginsWith("Radiator"))
                sRadLayer.push_back(cont[0].Data());
            if (line.BeginsWith("Radiator"))
                tRadLayer.push_back(cont[1].Atof());
            if (line.BeginsWith("TransLayer"))
                sTransLayer.push_back(cont[0].Data());
            if (line.BeginsWith("TransLayer"))
                pTransLayer.push_back(cont[1].Atof());
            if (line.BeginsWith("TransThick"))
                tTransLayer = TString(cont[0].Data()).Atof();
            if (line.BeginsWith("Impurities"))
                sImpurities.push_back(cont[0].Data());
            if (line.BeginsWith("Impurities"))
                pImpurities.push_back(cont[1].Atof());
            if (line.BeginsWith("PCDetector"))
                PCDetector = cont[0].Data();
            if (line.BeginsWith("PCDetector"))
                pixel = cont[1].Atof();
            if (line.BeginsWith("PixelRange"))
                XBinMin = cont[0].Atof();
            if (line.BeginsWith("PixelRange"))
                XBinMax = cont[1].Atof();
            if (line.BeginsWith("PixelRange"))
                YBinMin = cont[2].Atof();
            if (line.BeginsWith("PixelRange"))
                YBinMax = cont[3].Atof();
            if (line.BeginsWith("Lambda"))
                nLambda = cont[0].Atoi();
            if (line.BeginsWith("Lambda"))
                lambdaMin = cont[1].Atof();
            if (line.BeginsWith("Lambda"))
                lambdaMax = cont[2].Atof();
            if (line.BeginsWith("trkStep"))
                trkStep = cont[0].Atof();
            if (line.BeginsWith("nPhi"))
                nphi = int(cont[0].Atof() * 360);
            if (line.BeginsWith("Precision"))
                epsilon = cont[0].Atof();

            if (line.BeginsWith("SelectedMat"))
                selectedMat = cont;
            if (line.BeginsWith("SelectedDet"))
                selectedDet = cont;

            if (line.BeginsWith("Particle"))
                particle = cont[0];
            if (line.BeginsWith("Particle"))
                momentum = cont[1].Atof();
            if (line.BeginsWith("Particle"))
            {
                if (cont[2].Index("deg") > 0)
                    theta0 = cont[2].Atof() / 180 * TMath::Pi();
                else
                    theta0 = cont[2].Atof(); //default is rad
            }

            continue;
        }
        pos.fY = (nline++);
    }

    nRadLayer = (int)tRadLayer.size();
    nTransLayer = (int)pTransLayer.size();
    nImpurities = (int)pImpurities.size();
    nSelectedMat = (int)selectedMat.size();
    nSelectedDet = (int)selectedDet.size();
    NXBin = (XBinMax - XBinMin) / pixel;
    NYBin = (YBinMax - YBinMin) / pixel;
    for (int i = 0; i < nTransLayer; i++)
        mapTransLayers[sTransLayer[i]] = pTransLayer[i];
    for (int i = 0; i < nImpurities; i++)
        mapImpurities[sImpurities[i]] = pImpurities[i];
}

//______________________________________________________________________________
// basic button actions
void MyNumEvalClass::ExecButtonClick(Long_t bid, const char *cmdStr)
{
    ReadSettingsText();
        
    if (bid == DrawConfig)
        DoDrawConfig("Show the configurations");
    if (bid == LoadTextBuf)
        DoLoadTextBuf();
    if (bid == DrawSelectedMat)
        DoDrawSelectedMat();
    if (bid == DrawSelectedDet)
        DoDrawSelectedDet();
    if (bid == DrawRICH)
        DoDrawRICH();
    if (bid == Draw3Rings)
        DoDrawRings();
    if (bid == Save3Rings)
        DoSaveRings(cmdStr);
    if (bid == Load3Rings)
        DoLoadRings(cmdStr);
    if (MatList <= bid && bid < DetList)
        ShowMaterialInfo(gMyDatabaseClass->GetMaterialName(bid - MatList));
    if (DetList <= bid && bid < AnalysisAction)
        ShowDetectorInfo(gMyDatabaseClass->GetDetectorName(bid - DetList));
}

void MyNumEvalClass::DoDrawConfig(TString scap)
{
    double x0 = 0.02, x1 = 0.98;
    double y0 = 0.10, y1 = 0.15, y2 = 0.50, y3 = 0.98;

    gMyMainFrameGui->ClearCanvas(0);
    gMyMainFrameGui->SwitchCanvas(0);

    //detector
    TBox *bdet = new TBox(x0, y0, x1, y1);
    bdet->SetFillColor(kOrange - 2);
    bdet->Draw();
    TText *tdet = new TText(0.5, (y0 + y1) / 2., Form("%s = %.1f[mm]", PCDetector.c_str(), pixel));
    tdet->SetTextAlign(22);
    tdet->SetTextSize(.04);
    tdet->SetTextFont(42);
    tdet->Draw();

    //trans
    TBox *btrans = new TBox(x0, y1, x1, y2);
    btrans->SetFillColor(kCyan - 10);
    btrans->Draw();

    TString ttext1(""), ttext2("");
    for (int i = 0; i < nTransLayer; i++)
    {
        ttext1 += sTransLayer[i].c_str();
        ttext2 += Form("%.1f", pTransLayer[i]);
        if (i != nTransLayer - 1)
            ttext1 += ":";
        if (i != nTransLayer - 1)
            ttext2 += ":";
    }
    TString ttext3("");
    for (int i = 0; i < nImpurities; i++)
        ttext3 += Form("%s = %.1f[ppm]  ", sImpurities[i].c_str(), pImpurities[i]);
    TText *ttrans = new TText(0.5, (y1 + y2) / 2., ttext1 + " = " + ttext2);
    ttrans->SetTextAlign(22);
    ttrans->SetTextSize(.04);
    ttrans->SetTextFont(42);
    ttrans->Draw();

    TPaveText *pt = new TPaveText(.2, .25, 0.8, .3);
    pt->AddText(ttext3);
    pt->AddText(Form("Thickness = %.1f[mm]", tTransLayer));
    pt->SetFillColor(-1);
    pt->SetBorderSize(0);
    pt->SetTextFont(42);
    pt->Draw();

    //radiators
    int color[2] = {kAzure + 1, kAzure - 9};
    for (int i = 0; i < nRadLayer; i++)
    {
        double _y0 = y2 + (nRadLayer - i) * (y3 - y2) / nRadLayer;
        double _y1 = y2 + (nRadLayer - i - 1) * (y3 - y2) / nRadLayer;
        TBox *brad = new TBox(x0, _y0, x1, _y1);
        brad->SetFillColor(color[i % 2]);
        brad->Draw();
        TText *trad = new TText(0.5, (_y0 + _y1) / 2., Form("%s = %.1f[mm]", sRadLayer[i].c_str(), tRadLayer[i]));
        trad->SetTextAlign(22);
        trad->SetTextSize(.04);
        trad->SetTextFont(42);
        trad->Draw();
    }

    // current running parameters
    TText *tcap = new TText(0.5, 0.05, scap);
    tcap->SetTextFont(132);
    tcap->SetTextAlign(22);
    tcap->SetTextSize(.03);
    tcap->Draw();

    gMyMainFrameGui->UpdateCanvas(0);
}

void MyNumEvalClass::DoLoadTextBuf()
{
    gMyMainFrameGui->LoadText(GenerateSettingsText());
}

//______________________________________________________________________________
// analysis button actions
void MyNumEvalClass::DoDrawRICH()
{
    SetDetectorParameters();
    gMyCommonRICH->GenerateRICHRing();
    TH2F *f2 = gMyCommonRICH->Get2DViewer();
    TH1D **fl = gMyCommonRICH->Get1DViewer();

    gMyMainFrameGui->ClearAllCanvas();
    gMyMainFrameGui->SwitchCanvas(2);
    fl[0]->Draw("pl");
    gMyMainFrameGui->SwitchCanvas(3);
    fl[1]->Draw("pl");
    gMyMainFrameGui->SwitchCanvas(4);
    fl[2]->Draw("pl");
    gMyMainFrameGui->SwitchCanvas(5);
    fl[3]->Draw("pl");
    gMyMainFrameGui->SwitchCanvas(1);
    f2->Draw("colz");
    TPaveText *pt = new TPaveText(XBinMax * 5 / 8, YBinMax * 5 / 8, XBinMax * 9 / 10, YBinMax * 9 / 10);
    pt->AddText(Form("Nph=%.1f", f2->Integral()));
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->Draw();
    gMyMainFrameGui->UpdateCanvas(5);
    gMyMainFrameGui->UpdateCanvas(4);
    gMyMainFrameGui->UpdateCanvas(3);
    gMyMainFrameGui->UpdateCanvas(2);
    gMyMainFrameGui->UpdateCanvas(1);
}

void MyNumEvalClass::DoDrawRings()
{
    SetDetectorParameters();
    vector<TString> parlist;
    parlist.push_back("muon");
    parlist.push_back("pion");
    parlist.push_back("kaon");
    parlist.push_back("proton");

    gMyCommonRICH->GenerateMultiParticleRICHRings(parlist);
    gMyMainFrameGui->ClearAllCanvas();
    gMyMainFrameGui->SwitchCanvas(2);
    gMyCommonRICH->Get2DViewer(0)->Draw("colz");
    gMyMainFrameGui->SwitchCanvas(3);
    gMyCommonRICH->Get2DViewer(1)->Draw("colz");
    gMyMainFrameGui->SwitchCanvas(4);
    gMyCommonRICH->Get2DViewer(2)->Draw("colz");
    gMyMainFrameGui->SwitchCanvas(5);
    gMyCommonRICH->Get2DViewer(3)->Draw("colz");

    gMyMainFrameGui->SwitchCanvas(1);
    gMyStyle->SetTitle("Cherenkov Rings for mu/pi/K/p");
    gMyStyle->SetXLabel("X[mm]");
    gMyStyle->SetYLabel("Y[mm]");
    for (int i = 0; i < (int)parlist.size(); i++)
        gMyStyle->SetLegends(i, Form("%s Nph=%.1f", parlist[i].Data(), gMyCommonRICH->Get2DViewer(i)->Integral()));
    gMyStyle->SetDrawOption("");
    gMyStyle->Draw4Histograms(gMyCommonRICH->Get2DViewer(0), gMyCommonRICH->Get2DViewer(1), gMyCommonRICH->Get2DViewer(2), gMyCommonRICH->Get2DViewer(3));
    gMyMainFrameGui->UpdateCanvas(5);
    gMyMainFrameGui->UpdateCanvas(4);
    gMyMainFrameGui->UpdateCanvas(3);
    gMyMainFrameGui->UpdateCanvas(2);
    gMyMainFrameGui->UpdateCanvas(1);
}

void MyNumEvalClass::DoSaveRings(const char *fname)
{
    if (fname == NULL)
        return;
    gMyCommonRICH->SaveRings(fname);
}

void MyNumEvalClass::DoLoadRings(const char *fname)
{
    if (fname == NULL)
        return;
    gMyCommonRICH->LoadRings(fname);
    MyRICHDetector *gDet = gMyCommonRICH->GetDetector();
    nRadLayer = gDet->nRadLayer;
    tRadLayer = gDet->tRadLayer;
    sRadLayer = gDet->sRadLayer;

    tTransLayer = gDet->tTransLayer;
    nTransLayer = gDet->nTransLayer;
    pTransLayer = gDet->pTransLayer;
    sTransLayer = gDet->sTransLayer;

    nImpurities = gDet->nImpurities;
    pImpurities = gDet->pImpurities;
    sImpurities = gDet->sImpurities;

    PCDetector = gDet->PCDetector;
    pixel = gDet->pixel;

    particle = gDet->particle;
    momentum = gDet->momentum;
    theta0 = gDet->theta0;

    nLambda = gDet->nLambda;
    lambdaMin = gDet->lambdaMin;
    lambdaMax = gDet->lambdaMax;

    XBinMin = gDet->XBinMin;
    XBinMax = gDet->XBinMax;
    YBinMin = gDet->YBinMin;
    YBinMax = gDet->YBinMax;

    epsilon = gDet->epsilon;
    trkStep = gDet->trkStep;
    nphi = gDet->nphi;

    DoLoadTextBuf();
    DoDrawConfig("");
}

//______________________________________________________________________________
// show infor button actions
void MyNumEvalClass::ShowMaterialInfo(TString matName)
{
    if (nLambda == 0 || lambdaMin >= lambdaMax)
        return;

    double ppm = mapImpurities[matName.Data()];

    gMyMainFrameGui->ClearCanvas(1);
    TCanvas *c1 = gMyMainFrameGui->GetCanvas(1);
    c1->Divide(2, 2);
    c1->cd(1);
    gMyDatabaseClass->GetMatAbsGraph(matName, nLambda, lambdaMin, lambdaMax, ppm);

    gMyStyle->Draw1Graph(gMyDatabaseClass->GetMatAbsGraph(matName, nLambda, lambdaMin, lambdaMax, ppm));
    c1->cd(2);
    gMyStyle->Draw1Graph(gMyDatabaseClass->GetMatTrsGraph(matName, nLambda, lambdaMin, lambdaMax, ppm));
    c1->cd(3);
    gMyStyle->Draw1Graph(gMyDatabaseClass->GetMatRefGraph(matName, nLambda, lambdaMin, lambdaMax));
    gMyMainFrameGui->UpdateCanvas(1);
}

void MyNumEvalClass::ShowDetectorInfo(TString detName)
{
    gMyMainFrameGui->ClearCanvas(1);
    gMyStyle->Draw1Graph(gMyDatabaseClass->GetDetQEGraph(detName, nLambda, lambdaMin, lambdaMax));
    gMyMainFrameGui->UpdateCanvas(1);
}

void MyNumEvalClass::DoDrawSelectedMat()
{
    vector<TGraph *> tmp1;
    vector<TGraph *> tmp2;

    for (int i = 0; i < nSelectedMat; i++)
    {
        TString matName = selectedMat[i];
        double ppm = mapImpurities[matName.Data()];
        TGraph *g1 = gMyDatabaseClass->GetMatTrsGraph(matName, nLambda, lambdaMin, lambdaMax, ppm);
        TGraph *g2 = gMyDatabaseClass->GetMatRefGraph(matName, nLambda, lambdaMin, lambdaMax);
        if (g1 != NULL)
            tmp1.push_back(g1);
        if (g2 != NULL)
            tmp2.push_back(g2);
        gMyStyle->SetLegends(i, matName);
    }

    gMyMainFrameGui->ClearCanvas(1);
    gMyMainFrameGui->SwitchCanvas(1);
    gMyStyle->DrawGraphs(tmp1);
    gMyMainFrameGui->UpdateCanvas(1);

    gMyMainFrameGui->ClearCanvas(2);
    gMyMainFrameGui->SwitchCanvas(2);
    gMyStyle->DrawGraphs(tmp2);
    gMyMainFrameGui->UpdateCanvas(2);
}

void MyNumEvalClass::DoDrawSelectedDet()
{
    vector<TGraph *> tmp1;

    for (int i = 0; i < nSelectedDet; i++)
    {
        TString detName = selectedDet[i];
        TGraph *g1 = gMyDatabaseClass->GetDetQEGraph(detName, nLambda, lambdaMin, lambdaMax);
        if (g1 != NULL)
            tmp1.push_back(g1);
        gMyStyle->SetLegends(i, detName);
    }

    gMyMainFrameGui->ClearCanvas(1);
    gMyMainFrameGui->SwitchCanvas(1);
    gMyStyle->DrawGraphs(tmp1);
    gMyMainFrameGui->UpdateCanvas(1);
}

//______________________________________________________________________________
// private func
void MyNumEvalClass::SetDetectorParameters(int id)
{
    gMyCommonRICH->SetPrecision(epsilon);
    gMyCommonRICH->SetTrackStepSize(trkStep);
    gMyCommonRICH->SetNumberOfPhiStep(nphi);
    gMyCommonRICH->SetDatabase(gMyDatabaseClass);

    MyRICHDetector *gDet = gMyCommonRICH->GetDetector(id);
    gDet->SetRadiator(nRadLayer, tRadLayer, sRadLayer);
    gDet->SetTransLayer(tTransLayer, nTransLayer, pTransLayer, sTransLayer);
    gDet->SetImpurities(nImpurities, pImpurities, sImpurities);
    gDet->SetDetector(PCDetector, pixel);

    gDet->SetParticleGun(particle, gMyCommonRICH->GetMass(particle), momentum, theta0);
    gDet->SetLambdaRange(nLambda, lambdaMin, lambdaMax);
    gDet->SetDetectorViewSize(XBinMin, XBinMax, YBinMin, YBinMax);

    gDet->SetPrecision(epsilon, trkStep, nphi);
}