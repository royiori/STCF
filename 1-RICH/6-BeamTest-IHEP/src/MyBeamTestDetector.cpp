#include "MyBeamTestDetector.h"

//______________________________________________________________________________
//
MyBeamTestDetector::MyBeamTestDetector(int _id)
{
    id = _id;
}

MyBeamTestDetector::~MyBeamTestDetector()
{
    if (fFile == 0)
        return;
    fFile->Close();
}

//_____读取&生成“设置文本”_______________________________________________________________________
//
TString MyBeamTestDetector::GenerateIntro()
{
    TString settings("");

    settings += Form("%s-Name: %s\n", surname.Data(), name.Data());
    settings += Form("%s-FullSize: %.2f %.2f %.2f\n", surname.Data(), siz[0], siz[1], siz[2]);
    settings += Form("%s-Position: %.2f %.2f %.2f\n", surname.Data(), pos[0], pos[1], pos[2]);
    settings += Form("%s-Rotation: %.2f %.2f\n", surname.Data(), rot[0], rot[1]);
    settings += Form("%s-NX-NY: %d %d\n", surname.Data(), Nstrp[0], Nstrp[1]);
    settings += Form("%s-WX-WY: %.2f %.2f %.2f\n", surname.Data(), strip[0], strip[1], LtoAnode);
    settings += Form("%s-Mappingg: %s\n", surname.Data(), mapping.Data());
    settings += Form("%s-NBoard: %d\n", surname.Data(), GetNBoard());
    settings += Form("%s-Boards: ", surname.Data());
    for (int j = 0; j < GetNBoard(); j++)
        settings += Form("%d ", boardList[j]);
    settings += Form("\n%s-NChips: %d\n", surname.Data(), GetNChip());
    settings += Form("%s-Chips: ", surname.Data());
    for (int j = 0; j < GetNChip(); j++)
        settings += Form("%d ", chipList[j]);
    settings += "\n\n";

    return settings;
}

void MyBeamTestDetector::SetParameters(TString line, vector<TString> cont)
{
    if (line.Index(Form("%d-Name", id)) > 0)
    {
        SetSurName(line.Remove(line.Index("-Name")));
        SetName(cont[0]);
    }
    if (line.Index(Form("%d-NX-NY", id)) > 0)
        SetNXNY(cont[0].Atoi(), cont[1].Atoi());
    if (line.Index(Form("%d-WX-WY", id)) > 0)
        SetWXWY(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
    if (line.Index(Form("%d-Mapping", id)) > 0)
        SetMapping(cont[0]);
    if (line.Index(Form("%d-FullSize", id)) > 0)
        SetFullSize(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
    if (line.Index(Form("%d-Position", id)) > 0)
        SetPosition(cont[0].Atof(), cont[1].Atof(), cont[2].Atof());
    if (line.Index(Form("%d-Rotation", id)) > 0)
        SetRotation(cont[0].Atof(), cont[1].Atof());
    if (line.Index(Form("%d-NBoard", id)) > 0)
        SetNBoard(cont[0].Atoi());
    if (line.Index(Form("%d-NChip", id)) > 0)
        SetNChip(cont[0].Atoi());
    if (line.Index(Form("%d-Boards", id)) > 0)
        for (int j = 0; j < GetNBoard(); j++)
            SetBoard(j, cont[j].Atoi());
    if (line.Index(Form("%d-Chips", id)) > 0)
        for (int j = 0; j < GetNChip(); j++)
            SetChip(j, cont[j].Atoi());
}

//_____电子学&ID转换_________________________________________________________________________
//
int MyBeamTestDetector::CalID(int ibd, int chip, int chan)
{
    return ibd * chipList.size() * 64 + chip * 64 + chan;
}

int MyBeamTestDetector::IndexBD(int board)
{
    int index = -1;
    for (int i = 0; i < (int)boardList.size(); i++)
        if (board == boardList[i])
            index = i;
    return index;
}

int MyBeamTestDetector::IndexChip(int chip)
{
    int index = -1;
    for (int i = 0; i < (int)chipList.size(); i++)
        if (chip == chipList[i])
            index = i;
    return index;
}

//_____查找探测器的dst.root数据，根据trigger查找event_________________________________________________________________________
//
int MyBeamTestDetector::GetEntries()
{
    return (fTree == 0) ? 0 : fTree->GetEntries();
}

int MyBeamTestDetector::GetFirstTrig()
{
    return GetTrigID(0);
};

int MyBeamTestDetector::GetLastTrig()
{
    return GetTrigID(nentries - 1);
};

int MyBeamTestDetector::GetTrigID(int ientry)
{
    if (fTree == 0 || ientry < 0 || ientry > nentries - 1)
        return -1;

    fTree->GetEntry(ientry);
    return fEvent->event;
}

bool MyBeamTestDetector::SearchTrigID(int trigid, MyBeamTestHitData *fEventList)
{
    if (fFile == NULL || fTree == NULL || nentries == 0)
        return false;

    int nrange = 100;
    int beg = (iEntry - nrange < 0) ? 0 : iEntry - nrange;
    int end = (iEntry + nrange > nentries) ? nentries : iEntry + nrange;
    for (int i = beg; i < end; i++)
    {
        fTree->GetEntry(i);
        if (fEvent->event == trigid)
        {
            iEntry = i;
            //PushBack(fEventList->hit, fEvent);
            fEventList->detector.push_back(fEvent);
            return true;
        }
    }
    return false;
}
