#include "MyMainFrameGui.h"
#include "MyGuiActionClass.h"

MyGuiActionClass *gMyGuiActionClass = (MyGuiActionClass *)0;

//______________________________________________________________________________
//
MyGuiActionClass::MyGuiActionClass()
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
    nButton[0] = 6;
    iButton[0] = LoadDetFile;
    sButton[0].push_back("Load");
    sButton[0].push_back("Save");
    sButton[0].push_back("Show the spec RICH");
    sButton[0].push_back("Show multi-Particles RICH");
    sButton[0].push_back("Generate the XY-hitmap for SCAN");
    sButton[0].push_back("Scan the Rec-RICH");

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
        pTransLayer.push_back(env->GetValue(Form("pTrans%d", i), 50));
    }

    nImpurities = env->GetValue("nImpurities", 2);
    for (int i = 0; i < nImpurities; i++)
    {
        sImpurities.push_back(env->GetValue(Form("sImp%d", i), (i == 0) ? "H2O" : "O2"));
        pImpurities.push_back(env->GetValue(Form("pImp%d", i), 1));
    }

    PCDetector = env->GetValue("PCDetector", "CsI");
    pixel = env->GetValue("pixel", 5);

    nLambda = env->GetValue("nLambda", 60);
    lambdaMin = env->GetValue("lambdaMin", 160);
    lambdaMax = env->GetValue("lambdaMax", 220);

    particle = env->GetValue("particle", "pi");
    momentum = env->GetValue("momentum", 2.);
    Theta0 = env->GetValue("Theta0", 0.);
    mass = gMyCommonRICH->GetMass(particle);

    XBinMin = env->GetValue("XBinMin", -200);
    XBinMax = env->GetValue("XBinMax", 200);
    YBinMin = env->GetValue("YBinMin", -200);
    YBinMax = env->GetValue("YBinMax", 200);
    NXBin = XBinMax - XBinMin;
    NYBin = YBinMax - YBinMin;

    epsilon = env->GetValue("epsilon", 1e-5);
    trkStep = env->GetValue("trkStep", 1e-1);
    nphi = env->GetValue("nphi", 360);

    np = env->GetValue("np", 10);
    pMin = env->GetValue("pMin", 0.1);
    pMax = env->GetValue("pMax", 4.0);

    nthe0 = env->GetValue("nthe0", 10);
    The0Min = env->GetValue("The0Min", 0);
    The0Max = env->GetValue("The0Max", 90);

    for (int i = 0; i < 5; i++)
        pList.push_back(env->GetValue(Form("pList%d", i), i * 0.4 + 0.2));
    for (int i = 0; i < 5; i++)
        thList.push_back(env->GetValue(Form("thList%d", i), i * 10));

    nSelectedMat = env->GetValue("nSelectedMat", 2);
    for (int i = 0; i < nSelectedMat; i++)
        selectedMat.push_back(env->GetValue(Form("selectedMat%d", i), (i == 0) ? "Quartz" : "MgF2"));

    nSelectedDet = env->GetValue("nSelectedDet", 2);
    for (int i = 0; i < nSelectedDet; i++)
        selectedDet.push_back(env->GetValue(Form("selectedDet%d", i), (i == 0) ? "CsI" : "APD"));
};

MyGuiActionClass::~MyGuiActionClass()
{
    // Destructor.
    // Store env
    
    env->SetValue("nRadLayer", nRadLayer);
    for (int i = 0; i < nRadLayer; i++)
    {
        env->SetValue(Form("sRad%d", i), sRadLayer[i].c_str());
        env->SetValue(Form("tRad%d", i), tRadLayer[i]);
    }


    env->SetValue("tTransLayer", tTransLayer);
    env->SetValue("nTransLayer", nTransLayer);
    for (int i = 0; i < nTransLayer; i++)
    {
        env->SetValue(Form("sTrans%d", i), sTransLayer[i].c_str());
        env->SetValue(Form("pTrans%d", i), pTransLayer[i]);
    }

    env->SetValue("nImpurities", nImpurities);
    for (int i = 0; i < nImpurities; i++)
    {
        env->SetValue(Form("sImp%d", i), sImpurities[i].c_str());
        env->SetValue(Form("pImp%d", i), pImpurities[i]);
    }

    env->SetValue("PCDetector", PCDetector.c_str());
    env->SetValue("pixel", pixel);

    env->SetValue("nLambda", nLambda);
    env->SetValue("lambdaMin", lambdaMin);
    env->SetValue("lambdaMax", lambdaMax);

    env->SetValue("particle", particle);
    env->SetValue("momentum", momentum);
    env->SetValue("Theta0", Theta0);

    env->SetValue("XBinMin", XBinMin);
    env->SetValue("XBinMax", XBinMax);
    env->SetValue("YBinMin", YBinMin);
    env->SetValue("YBinMax", YBinMax);

    env->SetValue("epsilon", epsilon);
    env->SetValue("trkStep", trkStep);
    env->SetValue("nphi", nphi);

    env->SetValue("np", np);
    env->SetValue("pMin", pMin);
    env->SetValue("pMax", pMax);

    env->SetValue("nthe0", nthe0);
    env->SetValue("The0Min", The0Min);
    env->SetValue("The0Max", The0Max);

    for (int i = 0; i < 5; i++)
        env->SetValue(Form("pList%d", i), pList[i]);
    for (int i = 0; i < 5; i++)
        env->SetValue(Form("thList%d", i), thList[i]);

    env->SetValue("nSelectedMat", nSelectedMat);
    for (int i = 0; i < nSelectedMat; i++)
        env->SetValue(Form("selectedMat%d", i), selectedMat[i]);

    env->SetValue("nSelectedDet", nSelectedDet);
    for (int i = 0; i < nSelectedDet; i++)
        env->SetValue(Form("selectedDet%d", i), selectedDet[i]);
    env->SaveLevel(kEnvLocal);
}

//______________________________________________________________________________
// read / generate settings
TString MyGuiActionClass::GenerateSettingsText()
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

    settings += "\n\n# Particle to show (momentum[GeV], angle[degree])\n";
    settings += Form("Particle: %s %.1f %.1f\n", particle.Data(), momentum, Theta0);

    settings += "\n\n# Particle to scan (N, Min, Max) for momentum[GeV/c])\n";
    settings += Form("Momentum: %d %.1f %.1f\n", np, pMin, pMax);

    settings += "\n\n# Particle to scan (N, Min, Max) for theta_0[degree])\n";
    settings += Form("Theta0: %d %.1f %.1f\n", nthe0, The0Min, The0Max);

    settings += "\n\n# Particle to show for momentum[GeV/c] and theta_0[degree])\n";
    settings += Form("ShowMom: %.1f %.1f %.1f %.1f %.1f\n", pList[0], pList[1], pList[2], pList[3], pList[4]);
    settings += Form("ShowThe: %.1f %.1f %.1f %.1f %.1f\n", thList[0], thList[1], thList[2], thList[3], thList[4]);

    settings += "\n#---------------------------";
    settings += "\n#    Precision settings";
    settings += "\n#---------------------------";

    settings += "\n\n# Step size for charge particle track[mm]\n";
    settings += Form("trkStep: %f\n", trkStep);

    settings += "\n# Number of Step for phi [0~2pi] * 360\n";
    settings += Form("nPhi: %d\n", int(nphi / 360.));

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

void MyGuiActionClass::ReadSettingsText()
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

            PCDetector = (line.BeginsWith("PCDetector")) ? cont[0].Data() : PCDetector;
            pixel = (line.BeginsWith("PCDetector")) ? cont[1].Atof() : pixel;

            XBinMin = (line.BeginsWith("PixelRange")) ? cont[0].Atof() : XBinMin;
            XBinMax = (line.BeginsWith("PixelRange")) ? cont[1].Atof() : XBinMax;
            YBinMin = (line.BeginsWith("PixelRange")) ? cont[2].Atof() : YBinMin;
            YBinMax = (line.BeginsWith("PixelRange")) ? cont[3].Atof() : YBinMax;

            nLambda = (line.BeginsWith("Lambda")) ? cont[0].Atoi() : nLambda;
            lambdaMin = (line.BeginsWith("Lambda")) ? cont[1].Atof() : lambdaMin;
            lambdaMax = (line.BeginsWith("Lambda")) ? cont[2].Atof() : lambdaMax;

            trkStep = (line.BeginsWith("trkStep")) ? cont[0].Atof() : trkStep;
            nphi = (line.BeginsWith("nPhi")) ? int(cont[0].Atof() * 360) : nphi;
            epsilon = (line.BeginsWith("Precision")) ? cont[0].Atof() : epsilon;

            selectedMat = (line.BeginsWith("SelectedMat")) ? cont : selectedMat;
            selectedDet = (line.BeginsWith("SelectedDet")) ? cont : selectedDet;

            particle = (line.BeginsWith("Particle")) ? cont[0] : particle;
            momentum = (line.BeginsWith("Particle")) ? cont[1].Atof() : momentum;
            Theta0 = (line.BeginsWith("Particle")) ? cont[2].Atof() : Theta0;

            np = (line.BeginsWith("Momentum")) ? cont[0].Atoi() : np;
            pMin = (line.BeginsWith("Momentum")) ? cont[1].Atof() : pMin;
            pMax = (line.BeginsWith("Momentum")) ? cont[2].Atof() : pMax;

            nthe0 = (line.BeginsWith("Theta0")) ? cont[0].Atoi() : nthe0;
            The0Min = (line.BeginsWith("Theta0")) ? cont[1].Atof() : The0Min;
            The0Max = (line.BeginsWith("Theta0")) ? cont[2].Atof() : The0Max;

            pList[0] = (line.BeginsWith("ShowMom")) ? cont[0].Atof() : pList[0];
            pList[1] = (line.BeginsWith("ShowMom")) ? cont[1].Atof() : pList[1];
            pList[2] = (line.BeginsWith("ShowMom")) ? cont[2].Atof() : pList[2];
            pList[3] = (line.BeginsWith("ShowMom")) ? cont[3].Atof() : pList[3];
            pList[4] = (line.BeginsWith("ShowMom")) ? cont[4].Atof() : pList[4];
            thList[0] = (line.BeginsWith("ShowThe")) ? cont[0].Atof() : thList[0];
            thList[1] = (line.BeginsWith("ShowThe")) ? cont[1].Atof() : thList[1];
            thList[2] = (line.BeginsWith("ShowThe")) ? cont[2].Atof() : thList[2];
            thList[3] = (line.BeginsWith("ShowThe")) ? cont[3].Atof() : thList[3];
            thList[4] = (line.BeginsWith("ShowThe")) ? cont[4].Atof() : thList[4];

            continue;
        }
        pos.fY = (nline++);
    }

    mass = gMyCommonRICH->GetMass(particle);
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
void MyGuiActionClass::ExecButtonClick(Long_t bid, const char *cmdStr)
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
    if (bid == LoadDetFile)
        DoLoadDetFile(cmdStr);
    if (bid == SaveDetFile)
        DoSaveDetFile(cmdStr);
    if (bid == ShowSpecRICH)
        DoShowSpecRICH(cmdStr);
    if (bid == ShowMulParRICH)
        DoShowMulParRICH(cmdStr);
    if (bid == GenRICHList)
        DoGenHitMaps(cmdStr);
    if (bid == RecRICHList)
        DoRecRings(cmdStr);
    if (MatList <= bid && bid < DetList)
        ShowMaterialInfo(gMyDatabaseClass->GetMaterialName(bid - MatList));
    if (DetList <= bid && bid < AnalysisAction)
        ShowDetectorInfo(gMyDatabaseClass->GetDetectorName(bid - DetList));
}

void MyGuiActionClass::DoDrawConfig(TString scap)
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

void MyGuiActionClass::DoLoadTextBuf()
{
    gMyMainFrameGui->LoadText(GenerateSettingsText());
}

//______________________________________________________________________________
// analysis button actions
void MyGuiActionClass::DoShowSpecRICH(TString cmdStr)
{
    SetDetectorParameters();
    if (cmdStr == "yes" || gMyCommonRICH->GetDetHitMap()->GetEntries() == 0)
        gMyCommonRICH->GenerateDetRing();

    gMyMainFrameGui->SwitchCanvas(2);
    gMyCommonRICH->GetDetWLMap()->Draw();
    gMyMainFrameGui->SwitchCanvas(1);
    gMyCommonRICH->GetDetHitMap()->Draw("colz");
    gMyMainFrameGui->UpdateCanvas(2);
    gMyMainFrameGui->UpdateCanvas(1);
}

void MyGuiActionClass::DoShowMulParRICH(TString cmdStr)
{
    SetDetectorParameters();
    if (cmdStr == "yes" || gMyCommonRICH->GetDetListNumber()==0)
        gMyCommonRICH->GenerateMultiParticleRICHRings();

    gMyMainFrameGui->ClearAllCanvas();
    gMyMainFrameGui->SwitchCanvas(2);
    gMyCommonRICH->GetDetListHitMap(0)->Draw("colz");
    gMyMainFrameGui->SwitchCanvas(3);
    gMyCommonRICH->GetDetListHitMap(1)->Draw("colz");
    gMyMainFrameGui->SwitchCanvas(4);
    gMyCommonRICH->GetDetListHitMap(2)->Draw("colz");
    gMyMainFrameGui->SwitchCanvas(5);
    gMyCommonRICH->GetDetListHitMap(3)->Draw("colz");

    gMyMainFrameGui->SwitchCanvas(1);
    gMyStyle->SetTitle("Cherenkov Rings for mu/pi/K/p");
    gMyStyle->SetXLabel("X[mm]");
    gMyStyle->SetYLabel("Y[mm]");
    for (int ihypo = 0; ihypo < gMyCommonRICH->GetDetector()->nhypo; ihypo++)
        gMyStyle->SetLegends(ihypo, Form("%s Nph=%.1f", gMyCommonRICH->GetDetList(ihypo)->particle.Data(), gMyCommonRICH->GetDetListHitMap(ihypo)->Integral()));
    gMyStyle->SetDrawOption("");
    gMyStyle->Draw4Histograms(gMyCommonRICH->GetDetListHitMap(0), gMyCommonRICH->GetDetListHitMap(1), gMyCommonRICH->GetDetListHitMap(2), gMyCommonRICH->GetDetListHitMap(3));
    gMyMainFrameGui->UpdateCanvas(5);
    gMyMainFrameGui->UpdateCanvas(4);
    gMyMainFrameGui->UpdateCanvas(3);
    gMyMainFrameGui->UpdateCanvas(2);
    gMyMainFrameGui->UpdateCanvas(1);
}

void MyGuiActionClass::DoGenHitMaps(TString cmdStr)
{
    SetDetectorParameters();
    if (cmdStr == "yes")
        gMyCommonRICH->GenerateTheScanHitMaps();

    for (int ican = 0; ican < 5; ican++)
    {
        MyRICHDetector *gDet = gMyCommonRICH->GetDetector();
        int imom = (pList[ican] - gDet->pMin) / gDet->pStep;
        int ithe = (thList[ican] - gDet->The0Min) / gDet->The0Step;

        if(imom >= gDet->np) continue;
        if(ithe >= gDet->nthe0) continue;

        gMyStyle->SetTitle(Form("Cherenkov Rings for [%.1fGeV/c, %.1f] mu/pi/K/p", gMyCommonRICH->GetDetScan(imom, ithe, 0)->momentum, gMyCommonRICH->GetDetScan(imom, ithe, 0)->Theta0));
        gMyStyle->SetXLabel("X[mm]");
        gMyStyle->SetYLabel("Y[mm]");
        for (int ihypo = 0; ihypo < gDet->nhypo; ihypo++)
            gMyStyle->SetLegends(ihypo, Form("%s Nph=%.1f", gMyCommonRICH->GetDetScan(imom, ithe, ihypo)->particle.Data(), gMyCommonRICH->GetDetScanHitMap(imom, ithe, ihypo)->Integral()));
        gMyStyle->SetDrawOption("");
        gMyMainFrameGui->SwitchCanvas(ican + 1);
        gMyStyle->Draw4Histograms(gMyCommonRICH->GetDetScanHitMap(imom, ithe, 0), gMyCommonRICH->GetDetScanHitMap(imom, ithe, 1), gMyCommonRICH->GetDetScanHitMap(imom, ithe, 2), gMyCommonRICH->GetDetScanHitMap(imom, ithe, 3));
        gMyMainFrameGui->UpdateCanvas(ican + 1);
    }
}

void MyGuiActionClass::DoRecRings(TString cmdStr)
{
    SetDetectorParameters();
    if (cmdStr == "yes")
        gMyCommonRICH->GenerateTheOffsetAndResolutionMaps();
    
    /* 
    SetDetectorParameters();

    MyRICHDetector *gDet = gMyCommonRICH->GetDetector();

    gMyCommonRICH->ReconstructCerekovAngleDist(gDet);
    gMyMainFrameGui->ClearAllCanvas();
    gMyMainFrameGui->SwitchCanvas(2);
    gMyCommonRICH->Get2DRecRing(0)->Draw("colz");
    return;

    gMyMainFrameGui->SwitchCanvas(3);
    gMyCommonRICH->Get2DRing(1)->Draw("colz");
    gMyMainFrameGui->SwitchCanvas(4);
    gMyCommonRICH->Get2DRing(2)->Draw("colz");
    gMyMainFrameGui->SwitchCanvas(5);
    gMyCommonRICH->Get2DRing(3)->Draw("colz");

    gMyMainFrameGui->SwitchCanvas(1);
    gMyStyle->SetTitle("Cherenkov Rings for mu/pi/K/p");
    gMyStyle->SetXLabel("X[mm]");
    gMyStyle->SetYLabel("Y[mm]");
    for (int i = 0; i < gDet->nhypo; i++)
        gMyStyle->SetLegends(i, Form("%s Nph=%.1f", gDet->hypo[i], gMyCommonRICH->Get2DRing(i)->Integral()));
    gMyStyle->SetDrawOption("");
    gMyStyle->Draw4Histograms(gMyCommonRICH->Get2DRing(0), gMyCommonRICH->Get2DRing(1), gMyCommonRICH->Get2DRing(2), gMyCommonRICH->Get2DRing(3));
    gMyMainFrameGui->UpdateCanvas(5);
    gMyMainFrameGui->UpdateCanvas(4);
    gMyMainFrameGui->UpdateCanvas(3);
    gMyMainFrameGui->UpdateCanvas(2);
    gMyMainFrameGui->UpdateCanvas(1);
    */
}

void MyGuiActionClass::DoSaveDetFile(const char *fname)
{
    if (fname == NULL)
        return;
    gMyCommonRICH->SaveRings(fname);
}

void MyGuiActionClass::DoLoadDetFile(const char *fname)
{
    if (fname == NULL)
        return;
    gMyCommonRICH->LoadRings(fname);
    GetDetectorParameters();

    DoLoadTextBuf();
    DoDrawConfig("Show the configurations");
}



//______________________________________________________________________________
// show infor button actions
void MyGuiActionClass::ShowMaterialInfo(TString matName)
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

void MyGuiActionClass::ShowDetectorInfo(TString detName)
{
    gMyMainFrameGui->ClearCanvas(1);
    gMyStyle->Draw1Graph(gMyDatabaseClass->GetDetQEGraph(detName, nLambda, lambdaMin, lambdaMax));
    gMyMainFrameGui->UpdateCanvas(1);
}

void MyGuiActionClass::DoDrawSelectedMat()
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

void MyGuiActionClass::DoDrawSelectedDet()
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
void MyGuiActionClass::SetDetectorParameters()
{
    gMyCommonRICH->SetPrecision(epsilon);
    gMyCommonRICH->SetDatabase(gMyDatabaseClass);

    MyRICHDetector *gDet = gMyCommonRICH->GetDetector();

    gDet->SetRadiator(nRadLayer, tRadLayer, sRadLayer);
    gDet->SetTransLayer(tTransLayer, nTransLayer, pTransLayer, sTransLayer);
    gDet->SetImpurities(nImpurities, pImpurities, sImpurities);
    gDet->SetDetector(PCDetector, pixel);

    gDet->SetParticleGun(particle, mass);
    gDet->SetParticleGun(momentum, Theta0);
    gDet->SetLambdaRange(nLambda, lambdaMin, lambdaMax);
    gDet->SetDetectorViewSize(XBinMin, XBinMax, YBinMin, YBinMax);

    gDet->SetMomentumScanRange(np, pMin, pMax);
    gDet->SetTheta0ScanRange(nthe0, The0Min, The0Max);
    gDet->SetPrecision(epsilon, trkStep, nphi);
}

void MyGuiActionClass::GetDetectorParameters()
{
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
    Theta0 = gDet->Theta0;
    mass = gDet->mass;

    nLambda = gDet->nLambda;
    lambdaMin = gDet->lambdaMin;
    lambdaMax = gDet->lambdaMax;

    XBinMin = gDet->XBinMin;
    XBinMax = gDet->XBinMax;
    YBinMin = gDet->YBinMin;
    YBinMax = gDet->YBinMax;

    np = gDet->np;
    pMin = gDet->pMin;
    pMax = gDet->pMax;
    nthe0 = gDet->nthe0;
    The0Min = gDet->The0Min;
    The0Max = gDet->The0Max;

    epsilon = gDet->epsilon;
    trkStep = gDet->trkStep;
    nphi = gDet->nphi;
}