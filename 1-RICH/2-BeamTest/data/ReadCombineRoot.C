//-----------------------------
// hit击中信息合并处理后的 combine-DST-root 的数据结构
class MyBeamTestHitData
{
public:
    int event = -1; //触发事例号
    vector<int> detector;
    vector<double> charge;
    vector<pair<double, double>> hit;

    void Init(int ev)
    {
        event = ev;
        detector.clear();
        charge.clear();
        hit.clear();
    }
};

const char *DetName[4] = {"RICH", "T02", "T03", "T06"};

void ReadCombineRoot()
{
    TFile *fDSTFile = new TFile("./Data09131/Combine/Combined-dst.root");
    TTree *fDSTTree = (TTree *)fDSTFile->Get("tree");
    MyBeamTestHitData *fDSTEvent = 0;
    fDSTTree->SetBranchAddress("event", &fDSTEvent);
    long DSTNEntries = fDSTTree->GetEntries();

    for (int ip = 0; ip < 100; ip++)//DSTNEntries; ip++)
    {
        fDSTTree->GetEntry(ip);
        cout<<"Entry "<<ip<<": Event="<<fDSTEvent->event<<" NHits="<<fDSTEvent->detector.size()<<endl;
        for(int i=0; i<fDSTEvent->detector.size(); i++)
            cout<<"  Hit "<<i<<": @Detector "<<DetName[fDSTEvent->detector[i]]<<" Q="<<fDSTEvent->charge[i]<<" Pos=("<<fDSTEvent->hit[i].first<<", "<<fDSTEvent->hit[i].second<<")"<<endl;
    	cout<<endl;
        //getchar();
    }
}
