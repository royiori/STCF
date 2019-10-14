#include "TEnv.h"
#include "TView.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TPolyLine3D.h"
#include "TGLViewer.h"
#include "TGraphErrors.h"

#include "MyGuiActionClass.h"
#include "MyBeamTest.h"
#include "MyBeamTestDetector.h"

MyBeamTest *gMyBeamTest = (MyBeamTest *)0;

//______________________________________________________________________________
// 构造和析构函数
MyBeamTest::MyBeamTest(TEnv *ev)
{
    env = ev;
    geom = NULL;
    tracing = false;

    fDSPath = env->GetValue("DSPath", "NOTSet");

    NRICH = env->GetValue("NRICH", 1);
    NTrkAGT = env->GetValue("NTrkAGT", 3);
    NTrkVMM = env->GetValue("NTrkVMM", 2);

    for (int i = 0; i < NRICH; i++)
        vRICH.push_back(new MyBeamTestRICH(i, 32, 32, 5, 5, 30));
    for (int i = 0; i < NTrkAGT; i++)
        vTrkAGT.push_back(new MyBeamTestTrackAGET(i + NRICH, 128, 128, 0.4, 0.4));
    for (int i = 0; i < NTrkVMM; i++)
        vTrkVMM.push_back(new MyBeamTestTrackVMM(i + NRICH + NTrkAGT, 128, 128, 0.4, 0.4));

    for (int i = 0; i < NRICH; i++)
    {
        vRICH[i]->SetType(RICH);
        vRICH[i]->SetName(env->GetValue(Form("RICH%d-Name", i), "RICH"));
        vRICH[i]->SetFullSize(env->GetValue(Form("RICH%d-FullSize-x", i), 10), env->GetValue(Form("RICH%d-FullSize-y", i), 10), env->GetValue(Form("RICH%d-FullSize-z", i), 10));
        vRICH[i]->SetPosition(env->GetValue(Form("RICH%d-Position-x", i), 10), env->GetValue(Form("RICH%d-Position-y", i), 10), env->GetValue(Form("RICH%d-Position-z", i), 10));
        vRICH[i]->SetRotation(env->GetValue(Form("RICH%d-Rotation-x", i), 0), env->GetValue(Form("RICH%d-Rotation-y", i), 0));
        vRICH[i]->SetName(env->GetValue(Form("RICHName%d", i), Form("RICH%d", i)));

        //board: 2,10,5,7
        vRICH[i]->SetNBoard(env->GetValue(Form("RICH%d-Nboard", i), 4));
        vRICH[i]->SetNChip(env->GetValue(Form("RICH%d-NChip", i), 4));
        for (int j = 0; j < vRICH[i]->GetNBoard(); j++)
            vRICH[i]->SetBoard(j, env->GetValue(Form("RICH%d-board%d", i, j), 1));
        for (int j = 0; j < vRICH[i]->GetNChip(); j++)
            vRICH[i]->SetChip(j, env->GetValue(Form("RICH%d-chip%d", i, j), 1));

        vRICH[i]->ReadMapFile();
    }

    for (int i = 0; i < NTrkAGT; i++)
    {
        vTrkAGT[i]->SetType(TrackerAGET);
        vTrkAGT[i]->SetName(env->GetValue(Form("TrkAGET%d-Name", i), "Tk"));
        vTrkAGT[i]->SetFullSize(env->GetValue(Form("TrkAGET%d-FullSize-x", i), 10), env->GetValue(Form("TrkAGET%d-FullSize-y", i), 10), env->GetValue(Form("TrkAGET%d-FullSize-z", i), 10));
        vTrkAGT[i]->SetPosition(env->GetValue(Form("TrkAGET%d-Position-x", i), 10), env->GetValue(Form("TrkAGET%d-Position-y", i), 10), env->GetValue(Form("TrkAGET%d-Position-z", i), 10));
        vTrkAGT[i]->SetRotation(env->GetValue(Form("TrkAGET%d-Rotation-x", i), 0), env->GetValue(Form("TrkAGET%d-Rotation-y", i), 0));
        vTrkAGT[i]->SetName(env->GetValue(Form("TrkAGETName%d", i), Form("T%d", i)));

        //1(T06) 9(T03) 6(T02) <---(beam)
        vTrkAGT[i]->SetNBoard(env->GetValue(Form("TrkAGET%d-Nboard", i), 4));
        vTrkAGT[i]->SetNChip(env->GetValue(Form("TrkAGET%d-NChip", i), 4));
        for (int j = 0; j < vTrkAGT[i]->GetNBoard(); j++)
            vTrkAGT[i]->SetBoard(j, env->GetValue(Form("TrkAGET%d-board%d", i, j), 1));
        for (int j = 0; j < vTrkAGT[i]->GetNChip(); j++)
            vTrkAGT[i]->SetChip(j, env->GetValue(Form("TrkAGET%d-chip%d", i, j), 1));

        vTrkAGT[i]->ReadMapFile();
    }

    for (int i = 0; i < NTrkVMM; i++)
    {
        vTrkVMM[i]->SetType(TrackerVMM);
        vTrkVMM[i]->SetName(env->GetValue(Form("TrkVMM%d-Name", i), "Tk"));
        vTrkVMM[i]->SetFullSize(env->GetValue(Form("TrkVMM%d-FullSize-x", i), 10), env->GetValue(Form("TrkVMM%d-FullSize-y", i), 10), env->GetValue(Form("TrkVMM%d-FullSize-z", i), 10));
        vTrkVMM[i]->SetPosition(env->GetValue(Form("TrkVMM%d-Position-x", i), 10), env->GetValue(Form("TrkVMM%d-Position-y", i), 10), env->GetValue(Form("TrkVMM%d-Position-z", i), 10));
        vTrkVMM[i]->SetRotation(env->GetValue(Form("TrkVMM%d-Rotation-x", i), 0), env->GetValue(Form("TrkVMM%d-Rotation-y", i), 0));
        vTrkVMM[i]->SetName(env->GetValue(Form("TrkVMMName%d", i), Form("T%d", i)));

        vTrkVMM[i]->SetNBoard(env->GetValue(Form("TrkVMM%d-Nboard", i), 4));
        vTrkVMM[i]->SetNChip(env->GetValue(Form("TrkVMM%d-NChip", i), 4));
        for (int j = 0; j < vTrkVMM[i]->GetNBoard(); j++)
            vTrkVMM[i]->SetBoard(j, env->GetValue(Form("TrkVMM%d-board%d", i, j), 1));
        for (int j = 0; j < vTrkVMM[i]->GetNChip(); j++)
            vTrkVMM[i]->SetChip(j, env->GetValue(Form("TrkVMM%d-chip%d", i, j), 1));

        vTrkVMM[i]->ReadMapFile();
    }
}

MyBeamTest::~MyBeamTest()
{
    // Destructor.
    // Store env
    StoreEnv();
}

void MyBeamTest::StoreEnv()
{
    NRICH = (int)vRICH.size();
    NTrkAGT = (int)vTrkAGT.size();
    NTrkVMM = (int)vTrkVMM.size();

    env->SetValue("DSPath", fDSPath);

    env->SetValue("NRICH", NRICH);
    env->SetValue("NTrkAGT", NTrkAGT);
    env->SetValue("NTrkVMM", NTrkVMM);

    for (int i = 0; i < NRICH; i++)
    {
        env->SetValue(Form("RICH%d-Name", i), vRICH[i]->GetName());
        env->SetValue(Form("RICH%d-FullSize-x", i), vRICH[i]->GetFullSize()[0]);
        env->SetValue(Form("RICH%d-FullSize-y", i), vRICH[i]->GetFullSize()[1]);
        env->SetValue(Form("RICH%d-FullSize-z", i), vRICH[i]->GetFullSize()[2]);
        env->SetValue(Form("RICH%d-Position-x", i), vRICH[i]->GetPosition()[0]);
        env->SetValue(Form("RICH%d-Position-y", i), vRICH[i]->GetPosition()[1]);
        env->SetValue(Form("RICH%d-Position-z", i), vRICH[i]->GetPosition()[2]);
        env->SetValue(Form("RICH%d-Rotation-x", i), vRICH[i]->GetRotation()[0]);
        env->SetValue(Form("RICH%d-Rotation-y", i), vRICH[i]->GetRotation()[1]);

        env->SetValue(Form("RICHName%d", i), vRICH[i]->GetName());
        env->SetValue(Form("RICH%d-Nboard", i), vRICH[i]->GetNBoard());
        env->SetValue(Form("RICH%d-NChip", i), vRICH[i]->GetNChip());
        for (int j = 0; j < vRICH[i]->GetNBoard(); j++)
            env->SetValue(Form("RICH%d-board%d", i, j), vRICH[i]->GetBoardList()[j]);
        for (int j = 0; j < vRICH[i]->GetNChip(); j++)
            env->SetValue(Form("RICH%d-chip%d", i, j), vRICH[i]->GetChipList()[j]);
    }

    for (int i = 0; i < NTrkAGT; i++)
    {
        env->SetValue(Form("TrkAGET%d-Name", i), vTrkAGT[i]->GetName());
        env->SetValue(Form("TrkAGET%d-FullSize-x", i), vTrkAGT[i]->GetFullSize()[0]);
        env->SetValue(Form("TrkAGET%d-FullSize-y", i), vTrkAGT[i]->GetFullSize()[1]);
        env->SetValue(Form("TrkAGET%d-FullSize-z", i), vTrkAGT[i]->GetFullSize()[2]);
        env->SetValue(Form("TrkAGET%d-Position-x", i), vTrkAGT[i]->GetPosition()[0]);
        env->SetValue(Form("TrkAGET%d-Position-y", i), vTrkAGT[i]->GetPosition()[1]);
        env->SetValue(Form("TrkAGET%d-Position-z", i), vTrkAGT[i]->GetPosition()[2]);
        env->SetValue(Form("TrkAGET%d-Rotation-x", i), vTrkAGT[i]->GetRotation()[0]);
        env->SetValue(Form("TrkAGET%d-Rotation-y", i), vTrkAGT[i]->GetRotation()[1]);

        env->SetValue(Form("TrkAGETName%d", i), vTrkAGT[i]->GetName());
        env->SetValue(Form("TrkAGET%d-Nboard", i), vTrkAGT[i]->GetNBoard());
        env->SetValue(Form("TrkAGET%d-NChip", i), vTrkAGT[i]->GetNChip());
        for (int j = 0; j < vTrkAGT[i]->GetNBoard(); j++)
            env->SetValue(Form("TrkAGET%d-board%d", i, j), vTrkAGT[i]->GetBoardList()[j]);
        for (int j = 0; j < vTrkAGT[i]->GetNChip(); j++)
            env->SetValue(Form("TrkAGET%d-chip%d", i, j), vTrkAGT[i]->GetChipList()[j]);
    }

    for (int i = 0; i < NTrkVMM; i++)
    {
        env->SetValue(Form("TrkVMM%d-Name", i), vTrkVMM[i]->GetName());
        env->SetValue(Form("TrkVMM%d-FullSize-x", i), vTrkVMM[i]->GetFullSize()[0]);
        env->SetValue(Form("TrkVMM%d-FullSize-y", i), vTrkVMM[i]->GetFullSize()[1]);
        env->SetValue(Form("TrkVMM%d-FullSize-z", i), vTrkVMM[i]->GetFullSize()[2]);
        env->SetValue(Form("TrkVMM%d-Position-x", i), vTrkVMM[i]->GetPosition()[0]);
        env->SetValue(Form("TrkVMM%d-Position-y", i), vTrkVMM[i]->GetPosition()[1]);
        env->SetValue(Form("TrkVMM%d-Position-z", i), vTrkVMM[i]->GetPosition()[2]);
        env->SetValue(Form("TrkVMM%d-Rotation-x", i), vTrkVMM[i]->GetRotation()[0]);
        env->SetValue(Form("TrkVMM%d-Rotation-y", i), vTrkVMM[i]->GetRotation()[1]);

        env->SetValue(Form("TrkVMMName%d", i), vTrkVMM[i]->GetName());
        env->SetValue(Form("TrkVMM%d-Nboard", i), vTrkVMM[i]->GetNBoard());
        env->SetValue(Form("TrkVMM%d-NChip", i), vTrkVMM[i]->GetNChip());
        for (int j = 0; j < vTrkVMM[i]->GetNBoard(); j++)
            env->SetValue(Form("TrkVMM%d-board%d", i, j), vTrkVMM[i]->GetBoardList()[j]);
        for (int j = 0; j < vTrkVMM[i]->GetNChip(); j++)
            env->SetValue(Form("TrkVMM%d-chip%d", i, j), vTrkVMM[i]->GetChipList()[j]);
    }
    env->SaveLevel(kEnvLocal);
}

//______________________________________________________________________________
// 生成和读取configure文本
TString MyBeamTest::GenerateSettingText()
{
    TString settings("");

    settings += "\n#---------------------------";
    settings += "\n#    Beam-test settings";
    settings += "\n#---------------------------";
    settings += "\n\n# root file names for beam-test: \n";
    settings += Form("DSPath: %s\n", fDSPath.Data());

    settings += "\n\n#---------------------------";
    settings += "\n# RICH parameters: (in mm & degrees) \n";
    settings += Form("NRICH: %d\n", NRICH);
    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        settings += Form("RICH%d-Name: %s\n", i, vRICH[i]->GetName().Data());
        settings += Form("RICH%d-FullSize: %.2f %.2f %.2f\n", i, vRICH[i]->GetFullSize()[0], vRICH[i]->GetFullSize()[1], vRICH[i]->GetFullSize()[2]);
        settings += Form("RICH%d-Position: %.2f %.2f %.2f\n", i, vRICH[i]->GetPosition()[0], vRICH[i]->GetPosition()[1], vRICH[i]->GetPosition()[2]);
        settings += Form("RICH%d-Rotation: %.2f %.2f\n", i, vRICH[i]->GetRotation()[0], vRICH[i]->GetRotation()[1]);
        settings += Form("RICH%d-NBoard: %d\n", i, vRICH[i]->GetNBoard());
        settings += Form("RICH%d-Boards: ", i);
        for (int j = 0; j < vRICH[i]->GetNBoard(); j++)
            settings += Form("%d ", vRICH[i]->GetBoardList()[j]);
        settings += Form("\nRICH%d-NChips: %d\n", i, vRICH[i]->GetNChip());
        settings += Form("RICH%d-Chips: ", i);
        for (int j = 0; j < vRICH[i]->GetNChip(); j++)
            settings += Form("%d ", vRICH[i]->GetChipList()[j]);
        settings += "\n\n";
    }

    settings += "\n\n#---------------------------";
    settings += "\n# Tracker with AGET parameters: (in mm & degrees) \n";
    settings += Form("NTrkAGT: %d\n", NTrkAGT);
    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        settings += Form("TrkAGET%d-Name: %s\n", i, vTrkAGT[i]->GetName().Data());
        settings += Form("TrkAGET%d-FullSize: %.2f %.2f %.2f\n", i, vTrkAGT[i]->GetFullSize()[0], vTrkAGT[i]->GetFullSize()[1], vTrkAGT[i]->GetFullSize()[2]);
        settings += Form("TrkAGET%d-Position: %.2f %.2f %.2f\n", i, vTrkAGT[i]->GetPosition()[0], vTrkAGT[i]->GetPosition()[1], vTrkAGT[i]->GetPosition()[2]);
        settings += Form("TrkAGET%d-Rotation: %.2f %.2f\n", i, vTrkAGT[i]->GetRotation()[0], vTrkAGT[i]->GetRotation()[1]);
        settings += Form("TrkAGET%d-NBoard: %d\n", i, vTrkAGT[i]->GetNBoard());
        settings += Form("TrkAGET%d-Boards: ", i);
        for (int j = 0; j < vTrkAGT[i]->GetNBoard(); j++)
            settings += Form("%d ", vTrkAGT[i]->GetBoardList()[j]);
        settings += Form("\nTrkAGET%d-NChips: %d\n", i, vTrkAGT[i]->GetNChip());
        settings += Form("TrkAGET%d-Chips: ", i);
        for (int j = 0; j < vTrkAGT[i]->GetNChip(); j++)
            settings += Form("%d ", vTrkAGT[i]->GetChipList()[j]);
        settings += "\n\n";
    }

    settings += "\n\n#---------------------------";
    settings += "\n# Tracker with VMM parameters: (in mm & degrees) \n";
    settings += Form("NTrkVMM: %d\n", NTrkVMM);
    for (int i = 0; i < (int)vTrkVMM.size(); i++)
    {
        settings += Form("TrkVMM%d-Name: %s\n", i, vTrkVMM[i]->GetName().Data());
        settings += Form("TrkVMM%d-FullSize: %.2f %.2f %.2f\n", i, vTrkVMM[i]->GetFullSize()[0], vTrkVMM[i]->GetFullSize()[1], vTrkVMM[i]->GetFullSize()[2]);
        settings += Form("TrkVMM%d-Position: %.2f %.2f %.2f\n", i, vTrkVMM[i]->GetPosition()[0], vTrkVMM[i]->GetPosition()[1], vTrkVMM[i]->GetPosition()[2]);
        settings += Form("TrkVMM%d-Rotation: %.2f %.2f\n", i, vTrkVMM[i]->GetRotation()[0], vTrkVMM[i]->GetRotation()[1]);
        settings += Form("TrkVMM%d-NBoard: %d\n", i, vTrkVMM[i]->GetNBoard());
        settings += Form("TrkVMM%d-Boards: ", i);
        for (int j = 0; j < vTrkVMM[i]->GetNBoard(); j++)
            settings += Form("%d ", vTrkVMM[i]->GetBoardList()[j]);
        settings += Form("\nTrkVMM%d-NChips: %d\n", i, vTrkVMM[i]->GetNChip());
        settings += Form("TrkVMM%d-Chips: ", i);
        for (int j = 0; j < vTrkVMM[i]->GetNChip(); j++)
            settings += Form("%d ", vTrkVMM[i]->GetChipList()[j]);
        settings += "\n\n";
    }

    return settings;
}

int MyBeamTest::ReadSettingsText(TGText *text)
{
    int nline = 0;
    TGLongPosition pos(0, nline);

    while (text->GetLineLength(nline) != -1)
    {
        if (text->GetChar(pos) != -1)
        {
            TString line(text->GetLine(pos, 100));
            pos.fY = (nline++);

            if (line.BeginsWith("#") || line.Length() < 5)
                continue;

            vector<TString> cont = gMyGuiActionClass->ReadContent(line);
            cont.erase(cont.begin());

            if (line.BeginsWith("NRICH"))
            {
                int ndet = cont[0].Atoi();
                if (ndet > NRICH)
                    for (int i = NRICH; i < ndet; i++)
                        vRICH.push_back(new MyBeamTestRICH(i, 32, 32, 5, 5, 30));
            }

            if (line.BeginsWith("RICH"))
            {
                for (int i = 0; i < NRICH; i++)
                {
                    if (line.BeginsWith(Form("RICH%d-Name", i)))
                        vRICH[i]->SetName(cont[0]);
                    if (line.BeginsWith(Form("RICH%d-FullSize", i)))
                        vRICH[i]->SetFullSize(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
                    if (line.BeginsWith(Form("RICH%d-Position", i)))
                        vRICH[i]->SetPosition(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
                    if (line.BeginsWith(Form("RICH%d-Rotation", i)))
                        vRICH[i]->SetRotation(cont[0].Atof(), cont[1].Atof());
                    if (line.BeginsWith(Form("RICH%d-NBoard", i)))
                        vRICH[i]->SetNBoard(cont[0].Atoi());
                    if (line.BeginsWith(Form("RICH%d-NChip", i)))
                        vRICH[i]->SetNChip(cont[0].Atoi());
                    if (line.BeginsWith(Form("RICH%d-Boards", i)))
                        for (int j = 0; j < vRICH[i]->GetNBoard(); j++)
                            vRICH[i]->SetBoard(j, cont[j].Atoi());
                    if (line.BeginsWith(Form("RICH%d-Chips", i)))
                        for (int j = 0; j < vRICH[i]->GetNChip(); j++)
                            vRICH[i]->SetChip(j, cont[j].Atoi());
                }
            }

            if (line.BeginsWith("NTrkAGET"))
            {
                int ndet = cont[0].Atoi();
                if (ndet > NTrkAGT)
                    for (int i = NTrkAGT; i < ndet; i++)
                        vTrkAGT.push_back(new MyBeamTestTrackAGET(i, 128, 128, 0.4, 0.4));
            }

            if (line.BeginsWith("TrkAGET"))
            {
                for (int i = 0; i < NTrkAGT; i++)
                {
                    if (line.BeginsWith(Form("TrkAGET%d-Name", i)))
                        vTrkAGT[i]->SetName(cont[0]);
                    if (line.BeginsWith(Form("TrkAGET%d-FullSize", i)))
                        vTrkAGT[i]->SetFullSize(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
                    if (line.BeginsWith(Form("TrkAGET%d-Position", i)))
                        vTrkAGT[i]->SetPosition(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
                    if (line.BeginsWith(Form("TrkAGET%d-Rotation", i)))
                        vTrkAGT[i]->SetRotation(cont[0].Atof(), cont[1].Atof());
                    if (line.BeginsWith(Form("TrkAGET%d-NBoard", i)))
                        vTrkAGT[i]->SetNBoard(cont[0].Atoi());
                    if (line.BeginsWith(Form("TrkAGET%d-NChip", i)))
                        vTrkAGT[i]->SetNChip(cont[0].Atoi());
                    if (line.BeginsWith(Form("TrkAGET%d-Boards", i)))
                        for (int j = 0; j < vTrkAGT[i]->GetNBoard(); j++)
                            vTrkAGT[i]->SetBoard(j, cont[j].Atoi());
                    if (line.BeginsWith(Form("TrkAGET%d-Chips", i)))
                        for (int j = 0; j < vTrkAGT[i]->GetNChip(); j++)
                            vTrkAGT[i]->SetChip(j, cont[j].Atoi());
                }
            }

            if (line.BeginsWith("NTrkVMM"))
            {
                int ndet = cont[0].Atoi();
                if (ndet > NTrkVMM)
                    for (int i = NTrkVMM; i < ndet; i++)
                        vTrkVMM.push_back(new MyBeamTestTrackVMM(i, 128, 128, 0.4, 0.4));
            }

            if (line.BeginsWith("TrkVMM"))
            {
                for (int i = 0; i < NTrkVMM; i++)
                {
                    if (line.BeginsWith(Form("TrkVMM%d-Name", i)))
                        vTrkVMM[i]->SetName(cont[0]);
                    if (line.BeginsWith(Form("TrkVMM%d-FullSize", i)))
                        vTrkVMM[i]->SetFullSize(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
                    if (line.BeginsWith(Form("TrkVMM%d-Position", i)))
                        vTrkVMM[i]->SetPosition(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
                    if (line.BeginsWith(Form("TrkVMM%d-Rotation", i)))
                        vTrkVMM[i]->SetRotation(cont[0].Atof(), cont[1].Atof());
                    if (line.BeginsWith(Form("TrkVMM%d-NBoard", i)))
                        vTrkVMM[i]->SetNBoard(cont[0].Atoi());
                    if (line.BeginsWith(Form("TrkVMM%d-NChip", i)))
                        vTrkVMM[i]->SetNChip(cont[0].Atoi());
                    if (line.BeginsWith(Form("TrkVMM%d-Boards", i)))
                        for (int j = 0; j < vTrkVMM[i]->GetNBoard(); j++)
                            vTrkVMM[i]->SetBoard(j, cont[j].Atoi());
                    if (line.BeginsWith(Form("TrkVMM%d-Chips", i)))
                        for (int j = 0; j < vTrkVMM[i]->GetNChip(); j++)
                            vTrkVMM[i]->SetChip(j, cont[j].Atoi());
                }
            }

            continue;
        }
        pos.fY = (nline++);
    }

    UpdateDetID();
    return 1;
}

//______________________________________________________________________________
// 画出结构示意图
void MyBeamTest::DrawConfig()
{
    if (geom != NULL)
        delete geom;

    //------------------------------
    // geom顶层结构为top
    geom = new TGeoManager("beamgeo", "Beamtest geometry");
    TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
    media = new TGeoMedium("Vacuum", 1, mat);
    top = geom->MakeBox("Top", media, 1000., 500., 500.);
    geom->SetTopVolume(top);
    geom->SetTopVisible(false);
    geom->SetVerboseLevel(0);

    //------------------------------
    // top -> det -> sub-detectors
    TGeoVolume *det = geom->MakeBox("Detector", media, 1000., 500., 500.);
    det->SetVisibility(kFALSE);
    top->AddNode(det, 1, 0);

    for (int i = 0; i < NRICH; i++)
        vRICH[i]->DrawDetector(geom, det, media);
    for (int i = 0; i < NTrkAGT; i++)
        vTrkAGT[i]->DrawDetector(geom, det, media);
    for (int i = 0; i < NTrkVMM; i++)
        vTrkVMM[i]->DrawDetector(geom, det, media);

    TGeoVolume *ground = geom->MakeBox("Ground", media, 1000., 1., 1000);
    det->AddNode(ground, 1, new TGeoTranslation(500, 500, 0));
    ground->SetFillColor(kGray);

    //make beam
    {
        TGeoVolume *arrow = geom->MakeBox("arrow", media, 1000., 500., 500.);
        arrow->SetVisibility(kFALSE);

        //Z axis - 实验室坐标系
        TGeoVolume *xbar = geom->MakeBox("xbar", media, 550., 1., 1.);
        xbar->SetLineColor(kRed);
        arrow->AddNode(xbar, 1, new TGeoTranslation(450., 0., 0.));
        TGeoVolume *xarr = geom->MakeCone("xarr", media, 20, 0., 0.1, 0., 10.);
        xarr->SetLineColor(kRed);
        arrow->AddNode(xarr, 1, new TGeoCombiTrans(1000., 0., 0., new TGeoRotation("arrowRotX", 90, -90, 0)));

        //Y axis - 实验室坐标系
        TGeoVolume *ybar = geom->MakeBox("ybar", media, 1., 100., 1.);
        ybar->SetLineColor(kGreen);
        arrow->AddNode(ybar, 1, new TGeoTranslation(-100., -90., 0.));
        TGeoVolume *yarr = geom->MakeCone("yarr", media, 20, 0., 0.1, 0., 10.);
        yarr->SetLineColor(kGreen);
        arrow->AddNode(yarr, 1, new TGeoCombiTrans(-100., -200., 0., new TGeoRotation("arrowRotY", 0, -90, -90)));

        //X axis - 实验室坐标系
        TGeoVolume *zbar = geom->MakeBox("zbar", media, 1., 1., 100.);
        zbar->SetLineColor(kBlue);
        arrow->AddNode(zbar, 1, new TGeoTranslation(-100., 0., -90.));
        TGeoVolume *zarr = geom->MakeCone("zarr", media, 20, 0., 0.1, 0., 10.);
        zarr->SetLineColor(kBlue);
        arrow->AddNode(zarr, 1, new TGeoCombiTrans(-100., 0., -240., new TGeoRotation("arrowRotZ", 0, 0, 0)));

        det->AddNode(arrow, 1, new TGeoTranslation(0., 0., 0.));
    }

    //------------------------------
    // top -> hits -> sub-det's hits
    geohits = geom->MakeBox("geohits", media, 1000., 500., 500.);
    geohits->SetVisibility(kFALSE);
    top->AddNode(geohits, 1, 0);

    top->Draw();

    ptext = new TPaveText(.1, .8, .9, .97);
    ptext->AddText("Beam-test event display");
    ptext->Draw();
}

void MyBeamTest::SetDSPath(const char *fileName)
{
    fDSPath = TString(fileName);
    if (!fDSPath.EndsWith("idle.dat"))
    {
        cout << "##### you can simply choose 'idle.dat', all files will be generated automatically." << endl;
        fDSPath = TString("NOSET");
    }
}

TString MyBeamTest::GenPath(int type1, int type2, const char *fileName)
{
    TString fName = (fileName == NULL) ? fDSPath : TString(fileName);

    if (!fName.EndsWith("idle.dat"))
    {
        cout << "##### you can simply choose 'idle.dat', all files will be generated automatically." << endl;
        return fName;
    }

    if (type1 == RICH && type2 == RAW)
        fName.ReplaceAll("idle.dat", "/Combine/RICH-raw.root");
    if (type1 == RICH && type2 == PED)
        fName.ReplaceAll("idle.dat", "/Combine/RICH-ped.root");

    if (type1 == TrackerAGET && type2 == RAW)
        fName.ReplaceAll("idle.dat", "/Combine/TrackAGET-raw.root");
    if (type1 == TrackerAGET && type2 == PED)
        fName.ReplaceAll("idle.dat", "/Combine/TrackAGET-ped.root");

    if (type1 == TrackerVMM && type2 == RAW)
        fName.ReplaceAll("idle.dat", "/Combine/TrackVMM-raw.root");
    if (type1 == TrackerVMM && type2 == PED)
        fName.ReplaceAll("idle.dat", "/Combine/TrackVMM-ped.root");

    if (type2 == DST)
        fName.ReplaceAll("idle.dat", "/Combine/");

    if (type1 == ALL && type2 == ALL)
        fName.ReplaceAll("idle.dat", "/Combine/Combined-dst.root");

    return fName;
}

//______________________________________________________________________________
// 转换RICH的raw-root 为 dst-root文件，同时分析波形得到Q/T
bool MyBeamTest::ReadRICHPed(TString fPedName)
{
    TFile *fPedFile = new TFile(fPedName);
    if (!fPedFile->IsOpen())
    {
        cout << "#### Can't open pedestal " << fPedName << " to read. Please check the path." << endl;
        return false;
    }
    cout << "--> Now reading the RICH pedestal file: " << fPedName << endl;

    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        vRICH[i]->SetNPedMean();
        vRICH[i]->SetNPedRMS();

        for (int ii = 0; ii < vRICH[i]->GetNBoard(); ii++)
            for (int jj = 0; jj < vRICH[i]->GetNChip(); jj++)
                for (int kk = 0; kk < 64; kk++)
                {
                    TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", vRICH[i]->GetBoardList()[ii], vRICH[i]->GetChipList()[jj], kk));
                    if (tmp == 0)
                    {
                        cout << "Check the RICH pedestal file: Ped_" << vRICH[i]->GetBoardList()[ii] << "_" << vRICH[i]->GetChipList()[jj] << "_" << kk << " exists or not" << endl;
                        continue;
                    }

                    double mean = 0;
                    double rms = 0;
                    if (tmp->GetEntries() > 0)
                    {
                        mean = tmp->GetFunction("gaus")->GetParameter(1);
                        rms = tmp->GetFunction("gaus")->GetParameter(2);
                    }
                    vRICH[i]->SetPedMean(ii, jj, kk, mean);
                    vRICH[i]->SetPedRMS(ii, jj, kk, rms);
                }
    }
    cout << "--> Pedestal file has been read.\n"
         << endl;
    return true;
}

bool MyBeamTest::ConvtRICHRoot(const char *fileName, int SaveWaveFlag)
{
    // 确认Raw-ROOT路径及文件存在
    TString fName = GenPath(RICH, RAW, fileName);
    TFile *fFile = new TFile(fName);
    if (!fFile->IsOpen())
        return false;

    // 确认ped-root路径及文件存在
    TString fPedName = GenPath(RICH, PED, fileName);
    if (!ReadRICHPed(fPedName))
        return false;

    cout << "--> Now opening RICH raw root file: " << fName << endl;

    // 生成DST的文件路径
    vector<TFile *> fFileList;
    vector<TTree *> fTreeList;
    vector<MyBeamTestData *> fEventList;

    fFileList.resize(vRICH.size());
    fTreeList.resize(vRICH.size());
    fEventList.resize(vRICH.size());

    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        TString fName2 = fName;
        fName2.ReplaceAll("RICH-raw.root", Form("%s-dst.root", vRICH[i]->GetName().Data()));
        fFileList[i] = new TFile(fName2, "recreate");
        fTreeList[i] = new TTree("tree", "beam test data structure");
        fEventList[i] = new MyBeamTestData;
        fTreeList[i]->Branch("event", "MyBeamTestData", &fEventList[i], 8000, 2);
    }

    // 读入Raw-ROOT文件
    TTree *tree = (TTree *)fFile->Get("tree");
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;
    vector<double> *wave = 0;

    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);

    int id = -1;
    long nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;
        tree->GetEntry(ii);

        //根据board判断这个数据属于哪个探测器，放入id中。
        id = -1;
        for (int j = 0; j < (int)vRICH.size(); j++)
            for (int k = 0; k < vRICH[j]->GetNBoard(); k++)
                if (board == vRICH[j]->GetBoardList()[k])
                    id = j;

        if (id == -1)
        {
            cout << "#### Warning: event id " << event << "'s board: " << board << " doesn't exist in the RICH board list." << endl;
            continue;
        }

        //判断这个事例所属于的探测器
        fEventList[id]->event = (fEventList[id]->event == -1) ? event : fEventList[id]->event;
        if ((fEventList[id]->event) - event > 0)
            cout << "#### Warning: " << vRICH[id]->GetName() << "'s event id changed backwards from " << fEventList[id]->event << " to " << event << "\n\n";

        //发现新的事例，则填tree
        if (fabs(fEventList[id]->event - event) != 0)
        {
            //cout<<"Store a new event: "<<fEventList[id]->event<<endl;
            fFileList[id]->cd();
            fEventList[id]->Analysis(0, 150, 200, 300); //计算Q/T, 参数为ped和charge的计算所用的范围
            vRICH[id]->AnalysisCluster(fEventList[id]); //分析cluster

            if (!SaveWaveFlag) //如果不保存波形数据，则清空波形数组
                for (int ii = 0; ii < (int)fEventList[id]->board.size(); ii++)
                    fEventList[id]->wave[ii].clear();
            fTreeList[id]->Fill();
            fEventList[id]->Init(event, vRICH[id]->id);
        }

        //数据保存进vector
        fEventList[id]->board.push_back(board);
        fEventList[id]->chip.push_back(chip);
        fEventList[id]->channel.push_back(channel);
        fEventList[id]->pedeMean.push_back(vRICH[id]->GetPedMean(board, chip, channel));
        fEventList[id]->pedeRms.push_back(vRICH[id]->GetPedRMS(board, chip, channel));
        int nwave = fEventList[id]->wave.size();
        fEventList[id]->wave.resize(nwave + 1);
        for (int j = 0; j < (int)wave->size(); j++)
            fEventList[id]->wave[nwave].push_back(wave->at(j));
    }

    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        fFileList[i]->cd();
        fTreeList[i]->Write();

        int type;
        TTree *fTree2 = new TTree("tree2", "tree2");
        fTree2->Branch("type", &type);
        type = vRICH[i]->id;
        fTree2->Fill();
        fTree2->Write();

        fFileList[i]->Flush();
        fFileList[i]->Close();
    }
    fFile->Close();

    cout << "--> RICH Data file converted to a new structure." << endl
         << endl;
    return true;
}

//______________________________________________________________________________
// 转换TrackAGET的raw-root 为 dst-root文件，同时分析波形得到Q/T
bool MyBeamTest::ReadTrackAGTPed(TString fPedName)
{
    TFile *fPedFile = new TFile(fPedName);
    if (!fPedFile->IsOpen())
    {
        cout << "#### Can't open pedestal " << fPedName << " to read. Please check the path." << endl;
        return false;
    }
    cout << "--> Now reading the Track-AGET pedestal file: " << fPedName << endl;

    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        vTrkAGT[i]->SetNPedMean();
        vTrkAGT[i]->SetNPedRMS();

        for (int ii = 0; ii < vTrkAGT[i]->GetNBoard(); ii++)
            for (int jj = 0; jj < vTrkAGT[i]->GetNChip(); jj++)
                for (int kk = 0; kk < 64; kk++)
                {
                    TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", vTrkAGT[i]->GetBoardList()[ii], vTrkAGT[i]->GetChipList()[jj], kk));
                    if (tmp == 0)
                    {
                        cout << "Check the TrackAGET pedestal file: Ped_" << vRICH[i]->GetBoardList()[ii] << "_" << vRICH[i]->GetChipList()[jj] << "_" << kk << " exists or not" << endl;
                        continue;
                    }
                    double mean = 0;
                    double rms = 0;
                    if (tmp->GetEntries() > 0)
                    {
                        mean = tmp->GetFunction("gaus")->GetParameter(1);
                        rms = tmp->GetFunction("gaus")->GetParameter(2);
                    }
                    vTrkAGT[i]->SetPedMean(ii, jj, kk, mean);
                    vTrkAGT[i]->SetPedRMS(ii, jj, kk, rms);
                }
    }
    cout << "--> Pedestal file has been read.\n"
         << endl;
    return true;
}

bool MyBeamTest::ConvtTrackAGTRoot(const char *fileName, int SaveWaveFlag)
{
    // 确认Raw-ROOT路径及文件存在
    TString fName = GenPath(TrackerAGET, RAW, fileName);
    TFile *fFile = new TFile(fName);
    if (!fFile->IsOpen())
        return false;

    // 确认ped-root路径及文件存在
    TString fPedName = GenPath(TrackerAGET, PED, fileName);
    if (!ReadTrackAGTPed(fPedName))
        return false;

    cout << "--> Now opening Track-AGET raw root file: " << fName << endl;

    // 生成DST的文件路径, 三个track会生成三个-dst.root文件
    vector<TFile *> fFileList;
    vector<TTree *> fTreeList;
    vector<MyBeamTestData *> fEventList;

    fFileList.resize(vTrkAGT.size());
    fTreeList.resize(vTrkAGT.size());
    fEventList.resize(vTrkAGT.size());

    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        TString fName2 = fName;
        fName2.ReplaceAll("TrackAGET-raw.root", Form("%s-dst.root", vTrkAGT[i]->GetName().Data()));
        fFileList[i] = new TFile(fName2, "recreate");
        fTreeList[i] = new TTree("tree", "beam test data structure");
        fEventList[i] = new MyBeamTestData;
        fTreeList[i]->Branch("event", "MyBeamTestData", &fEventList[i], 8000, 2);
    }

    // 读入Raw-ROOT文件
    TTree *tree = (TTree *)fFile->Get("tree");
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;
    vector<double> *wave = 0;

    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);

    int id = -1;
    long nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;
        tree->GetEntry(ii);

        //根据board判断这个数据属于哪个探测器，放入id中。
        id = -1;
        for (int j = 0; j < (int)vTrkAGT.size(); j++)
            for (int k = 0; k < vTrkAGT[j]->GetNBoard(); k++)
                if (board == vTrkAGT[j]->GetBoardList()[k])
                    id = j;

        if (id == -1)
        {
            cout << "#### Warning: event id " << event << "'s board: " << board << " doesn't exist in the Track-AGET board list." << endl;
            continue;
        }

        //判断这个事例所属于的探测器
        fEventList[id]->event = (fEventList[id]->event == -1) ? event : fEventList[id]->event;
        if ((fEventList[id]->event) - event > 0)
            cout << "#### Warning: " << vTrkAGT[id]->GetName() << "'s event id changed backwards from " << fEventList[id]->event << " to " << event << "\n\n";

        //发现新的事例，则填tree
        if (fabs(fEventList[id]->event - event) != 0)
        {
            //cout<<"Store a new event: "<<fEventList[id]->event<<endl;
            fFileList[id]->cd();
            fEventList[id]->Analysis(0, 150, 200, 300);   //计算Q/T
            vTrkAGT[id]->AnalysisCluster(fEventList[id]); //分析cluster

            if (!SaveWaveFlag) //如果不保存波形数据，则清空波形数组
                for (int ii = 0; ii < (int)fEventList[id]->board.size(); ii++)
                    fEventList[id]->wave[ii].clear();
            fTreeList[id]->Fill();
            fEventList[id]->Init(event, vTrkAGT[id]->id);
        }

        //数据保存进vector
        fEventList[id]->board.push_back(board);
        fEventList[id]->chip.push_back(chip);
        fEventList[id]->channel.push_back(channel);
        fEventList[id]->pedeMean.push_back(vTrkAGT[id]->GetPedMean(board, chip, channel));
        fEventList[id]->pedeRms.push_back(vTrkAGT[id]->GetPedRMS(board, chip, channel));
        int nwave = fEventList[id]->wave.size();
        fEventList[id]->wave.resize(nwave + 1);
        for (int j = 0; j < (int)wave->size(); j++)
            fEventList[id]->wave[nwave].push_back(wave->at(j));
    }

    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        fTreeList[i]->Fill();
        fFileList[i]->cd();
        fTreeList[i]->Write();

        int type;
        TTree *fTree2 = new TTree("tree2", "tree2");
        fTree2->Branch("type", &type);
        type = vTrkAGT[i]->id;
        fTree2->Fill();
        fTree2->Write();

        fFileList[i]->Flush();
        fFileList[i]->Close();
    }
    fFile->Close();

    cout << "--> Track-AGET Data file converted to a new structure." << endl
         << endl;
    return true;
}

//______________________________________________________________________________
// 转换TrackVMM的raw-root 为 dst-root文件
bool MyBeamTest::ConvtTrackVMMRoot(const char *fileName)
{
    // 确认Raw-ROOT路径及文件存在
    TString fName = GenPath(TrackerVMM, RAW, fileName);
    TFile *fFile = new TFile(fName);
    if (!fFile->IsOpen())
        return false;

    cout << "--> Now opening Track-VMM raw root file: " << fName << endl;

    // 生成DST的文件路径, 二个track会生成二个-dst.root文件
    vector<TFile *> fFileList;
    vector<TTree *> fTreeList;
    vector<MyBeamTestData *> fEventList;

    fFileList.resize(vTrkVMM.size());
    fTreeList.resize(vTrkVMM.size());
    fEventList.resize(vTrkVMM.size());

    for (int i = 0; i < (int)vTrkVMM.size(); i++)
    {
        TString fName2 = fName;
        fName2.ReplaceAll("TrackVMM-raw.root", Form("%s-dst.root", vTrkVMM[i]->GetName().Data()));
        fFileList[i] = new TFile(fName2, "recreate");
        fTreeList[i] = new TTree("tree", "beam test data structure");
        fEventList[i] = new MyBeamTestData;
        fTreeList[i]->Branch("event", "MyBeamTestData", &fEventList[i], 8000, 2);
    }

    // 读入Raw-ROOT文件
    TTree *tree = (TTree *)fFile->Get("tree");
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    UShort_t PDO;
    UShort_t BCID;
    UShort_t TDO;
    TBranch *b_Event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_PDO;
    TBranch *b_BCID;
    TBranch *b_TDO;

    tree->SetBranchAddress("event", &event, &b_Event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("PDO", &PDO, &b_PDO);
    tree->SetBranchAddress("BCID", &BCID, &b_BCID);
    tree->SetBranchAddress("TDO", &TDO, &b_TDO);

    int id = -1;
    long nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;
        tree->GetEntry(ii);

        //根据board判断这个数据属于哪个探测器，放入id中。
        id = -1;
        for (int j = 0; j < (int)vTrkVMM.size(); j++)
            for (int k = 0; k < vTrkVMM[j]->GetNBoard(); k++)
                if (board == vTrkVMM[j]->GetBoardList()[k])
                    id = j;

        if (id == -1)
        {
            cout << "#### Warning: event id " << event << "'s board: " << board << " doesn't exist in the Track-VMM board list." << endl;
            continue;
        }

        //判断这个事例所属于的探测器
        fEventList[id]->event = (fEventList[id]->event == -1) ? event : fEventList[id]->event;
        if ((fEventList[id]->event) - event > 0)
            cout << "#### Warning: " << vTrkVMM[id]->GetName() << "'s event id changed backwards from " << fEventList[id]->event << " to " << event << "\n\n";

        //发现新的事例，则填tree
        if (fabs(fEventList[id]->event - event) != 0)
        {
            fFileList[id]->cd();
            vTrkVMM[id]->AnalysisCluster(fEventList[id]); //分析cluster
            fTreeList[id]->Fill();
            fEventList[id]->Init(event, vTrkVMM[id]->id);
        }

        //数据保存进vector
        fEventList[id]->board.push_back(board);
        fEventList[id]->chip.push_back(chip);
        fEventList[id]->channel.push_back(channel);
        fEventList[id]->charge.push_back(PDO);
        fEventList[id]->time.push_back(300); //(TDO | BCID);
    }

    for (int i = 0; i < (int)vTrkVMM.size(); i++)
    {
        fTreeList[i]->Fill();
        fFileList[i]->cd();
        fTreeList[i]->Write();

        int type;
        TTree *fTree2 = new TTree("tree2", "tree2");
        fTree2->Branch("type", &type);
        type = vTrkVMM[i]->id;
        fTree2->Fill();
        fTree2->Write();

        fFileList[i]->Flush();
        fFileList[i]->Close();
    }
    fFile->Close();

    cout << "--> Track-VMM Data file converted to a new structure." << endl
         << endl;
    return true;
}

//______________________________________________________________________________
// 三个系统联合分析
void MyBeamTest::CombineDSTRoot(const char *fileName)
{
    //init
    for (int i = 0; i < (int)vRICH.size(); i++)
        vRICH[i]->ReadDSTRoot(GenPath(RICH, RAW, fileName));
    for (int i = 0; i < (int)vTrkAGT.size(); i++)
        vTrkAGT[i]->ReadDSTRoot(GenPath(TrackerAGET, RAW, fileName));
    for (int i = 0; i < (int)vTrkVMM.size(); i++)
        vTrkVMM[i]->ReadDSTRoot(GenPath(TrackerVMM, RAW, fileName));

    //
    //设定扫描trig的区间
    int nMax = 0;
    for (int i = 0; i < (int)vRICH.size(); i++)
    {
        int nEntries = vRICH[i]->GetEntries();
        int tmin = vRICH[i]->GetFirstTrig();
        int tmax = vRICH[i]->GetLastTrig();
        if (nEntries == 0)
            continue;
        nMax = (nMax < nEntries) ? nEntries : nMax;
        cout << "--> " << vRICH[i]->GetName() << " Events: " << nEntries << ", tigger range: " << tmin << " " << tmax << endl;
    }

    for (int i = 0; i < (int)vTrkAGT.size(); i++)
    {
        int nEntries = vTrkAGT[i]->GetEntries();
        int tmin = vTrkAGT[i]->GetFirstTrig();
        int tmax = vTrkAGT[i]->GetLastTrig();
        if (nEntries == 0)
            continue;
        nMax = (nMax < nEntries) ? nEntries : nMax;
        cout << "--> " << vTrkAGT[i]->GetName() << " Events: " << nEntries << ", tigger range: " << tmin << " " << tmax << endl;
    }

    for (int i = 0; i < (int)vTrkVMM.size(); i++)
    {
        int nEntries = vTrkVMM[i]->GetEntries();
        int tmin = vTrkVMM[i]->GetFirstTrig();
        int tmax = vTrkVMM[i]->GetLastTrig();
        if (nEntries == 0)
            continue;
        nMax = (nMax < nEntries) ? nEntries : nMax;
        cout << "--> " << vTrkVMM[i]->GetName() << " Events: " << nEntries << ", tigger range: " << tmin << " " << tmax << endl;
    }
    cout << "--> Searching trigger from " << 0 << " to " << nMax << endl;

    //
    //定义全部组合后的数据结构
    TString fName = GenPath(ALL, ALL, fileName);
    TFile *fFile = new TFile(fName, "recreate");
    TTree *fTree = new TTree("tree", "beam test data structure");
    MyBeamTestHitData *fEventList = new MyBeamTestHitData;
    fTree->Branch("event", "MyBeamTestHitData", &fEventList, 8000, 2);

    vector<pair<double, double>> hit;
    vector<double> q;
    vector<double> t;

    //
    //扫描trig来合并事例
    for (int trig = 0; trig < nMax; trig++)
    {
        if (trig % 100 == 0)
            cout << "--> Searching trigger = " << trig << endl;

        fEventList->Init(trig);
        for (int i = 0; i < (int)vRICH.size(); i++)
            vRICH[i]->SearchTrigID(trig, fEventList);

        for (int i = 0; i < (int)vTrkAGT.size(); i++)
            vTrkAGT[i]->SearchTrigID(trig, fEventList);

        for (int i = 0; i < (int)vTrkVMM.size(); i++)
            vTrkVMM[i]->SearchTrigID(trig, fEventList);

        fTree->Fill();
    }
    cout << "--> " << fTree->GetEntries() << " events are stored." << endl;
    fTree->Write();
    fFile->Close();

    for (int i = 0; i < (int)vRICH.size(); i++)
        vRICH[i]->CloseDSTRoot();
    for (int i = 0; i < (int)vTrkAGT.size(); i++)
        vTrkAGT[i]->CloseDSTRoot();
    for (int i = 0; i < (int)vTrkVMM.size(); i++)
        vTrkAGT[i]->CloseDSTRoot();
}

//______________________________________________________________________________
// 生成击中图
void MyBeamTest::LoadDSTRoot(const char *fileName)
{
    if (fDSTFile != 0)
    {
        fDSTFile->Close();
        fDSTFile = 0;
        fDSTTree = 0;
        fDSTEvent = 0;
    }

    //读入数据
    TString fName = GenPath(ALL, ALL, fileName);
    fDSTFile = new TFile(fName);
    if (!fDSTFile->IsOpen())
    {
        cout << "#### Error: Can't open " << fName << " to read." << endl;
        return;
    }
    cout << "--> Reading " << fName << "." << endl;

    fDSTTree = (TTree *)fDSTFile->Get("tree");
    fDSTTree->SetBranchAddress("event", &fDSTEvent);
    DSTNEntries = fDSTTree->GetEntries();
    cout << "--> Total hit entries = " << DSTNEntries << endl;

    //生成分布数据，计算效率，拟合并求残差
    int NTOT = NRICH + NTrkAGT + NTrkVMM;
    vector<int> hitList;  //此探测器是否击中的list
    vector<int> NHitList; //每个探测器有多少击中的list
    hitList.resize(NTOT);
    NHitList.resize(NTOT);
    int NHits = 0;

    for (int ip = 0; ip < DSTNEntries; ip++)
    {
        fDSTTree->GetEntry(ip);

        //生成分布
        for (int i = 0; i < NRICH; i++)
            hitList[i] = vRICH[i]->FillDistribution(ip, fDSTEvent);
        for (int i = 0; i < NTrkAGT; i++)
            hitList[i + NRICH] = vTrkAGT[i]->FillDistribution(ip, fDSTEvent);
        for (int i = 0; i < NTrkVMM; i++)
            hitList[i + NRICH + NTrkAGT] = vTrkVMM[i]->FillDistribution(ip, fDSTEvent);

        //计算效率
        for (int i = 0; i < NTOT; i++)
        {
            bool temp = 1;
            for (int ii = 0; ii < NTOT; ii++)
                temp *= ((ii == i) ? 1 : hitList[ii]);
            if (temp == 1)
                NHitList[i]++;
        }

        bool temp = 1;
        for (int i = 0; i < (int)hitList.size(); i++)
            temp *= hitList[i];
        if (temp == 1)
            NHits++;

        //拟合径迹
        /*vector<vector<vector<double>>> hitpos; //Ndetector, Ntrack, Ndim
        hitpos.resize(NTOT);
        for(int i=0; i<NTOT; i++)
            hitpos.resize(4); 
        for (int i = 0; i < NRICH; i++)
            vRICH[i]->FillHitVector(hitpos[i], fDSTEvent);
        for (int i = 0; i < NTrkAGT; i++)
            vTrkAGT[i]->FillHitVector(hitpos[i+NRICH], fDSTEvent);
        for (int i = 0; i < NTrkVMM; i++)
            vTrkVMM[i]->FillHitVector(hitops[i+NRICH+NTrkAGT], fDSTEvent);
        */
    }

    cout << "Total Events=" << DSTNEntries << ", Effect Events=" << NHits << ", Ratio=" << NHits * 1.0 / DSTNEntries << endl;
    for (int i = 0; i < NRICH; i++)
        cout << ": " << vRICH[i]->GetName() << " Efficiency = " << ((NHitList[i] > 0) ? NHits * 1.0 / NHitList[i] : 0) << endl;
    for (int i = 0; i < NTrkAGT; i++)
        cout << ": " << vTrkAGT[i]->GetName() << " Efficiency = " << ((NHitList[i + NRICH] > 0) ? NHits * 1.0 / NHitList[i + NRICH] : 0) << endl;
    for (int i = 0; i < NTrkVMM; i++)
        cout << ": " << vTrkVMM[i]->GetName() << " Efficiency = " << ((NHitList[i + NRICH + NTrkAGT] > 0) ? NHits * 1.0 / NHitList[i + NRICH + NTrkAGT] : 0) << endl;

    DrawCurtHit(0);
}

void MyBeamTest::DrawDSTHit(int entry)
{
    if (entry < 0 || entry > DSTNEntries - 1)
        return;

    ptext->Clear();
    ptext->AddText(Form("Event: %d", entry));

    //remove existing hits
    for (int i = 0; i < geohits->GetNdaughters(); i++)
    {
        geohits->RemoveNode(geohits->GetNode(0));
    }

    fDSTTree->GetEntry(entry);
    for (int i = 0; i < NRICH; i++)
        vRICH[i]->DrawHits(geom, geohits, media, fDSTEvent);
    for (int i = 0; i < NTrkAGT; i++)
        vTrkAGT[i]->DrawHits(geom, geohits, media, fDSTEvent);
    for (int i = 0; i < NTrkVMM; i++)
        vTrkVMM[i]->DrawHits(geom, geohits, media, fDSTEvent);

    top->Draw();
    ptext->Draw("same");
}

void MyBeamTest::DrawPrevHit()
{
    iEntry = (iEntry - 1 < 0) ? 0 : iEntry - 1;
    DrawDSTHit(iEntry);
}
void MyBeamTest::DrawNextHit()
{
    iEntry = (iEntry + 1 >= DSTNEntries) ? iEntry : iEntry + 1;
    DrawDSTHit(iEntry);
}

void MyBeamTest::DrawCurtHit(int entry)
{
    if (entry >= 0 && entry < DSTNEntries)
    {
        iEntry = entry;
        DrawDSTHit(iEntry);
    }
}

// for testing function only
void MyBeamTest::AnalysisDSTRoot(const char *fileName)
{
    if (fDSTFile != 0)
    {
        fDSTFile->Close();
        fDSTFile = 0;
        fDSTTree = 0;
        fDSTEvent = 0;
    }

    //读入数据
    TString fName = GenPath(ALL, ALL, fileName);
    fDSTFile = new TFile(fName);
    if (!fDSTFile->IsOpen())
    {
        cout << "#### Error: Can't open " << fName << " to read." << endl;
        return;
    }
    cout << "--> Reading " << fName << "." << endl;

    fDSTTree = (TTree *)fDSTFile->Get("tree");
    fDSTTree->SetBranchAddress("event", &fDSTEvent);
    DSTNEntries = fDSTTree->GetEntries();
    cout << "--> Total hit entries = " << DSTNEntries << endl;

    for (int ip = 0; ip < DSTNEntries; ip++)
    {
        if (ip % 1000 == 0 || ip < 10)
            cout << "--> event: " << ip << endl;

        fDSTTree->GetEntry(ip);
    }

    return;
}