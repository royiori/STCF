#ifndef MyBeamTestDetector_h
#define MyBeamTestDetector_h

using namespace std;

#include "TGeoManager.h"

enum
{
    RICH,
    TrackerAGET,
    TrackerVMM,
    ALL
};

enum
{
    RAW,
    DST,
    PED
};

class MyBeamTestRICH;
class MyBeamTestTrackAGET;
class MyBeamTestTrackVMM;

//-----------------------------
// hit击中信息合并处理后的击中数据结构
struct BeamHit //hit的击中channel坐标及q值
{
    int id;
    double q;
    double t;
    pair<double, double> hit;
};

struct RealBeamHit //cluster的中心值及q值
{
    int id;
    int nhit; 
    double q;
    double t;
    double hit[3];
    double hiterr[3];
};

//-----------------------------
// DST-root 的数据结构，包含的是每个hit及cluster的channel坐标，而不是绝对坐标
class MyBeamTestData
{
public:
    int id = -1; //探测器的id号
    int event = -1; //触发事例号
    vector<UShort_t> board;
    vector<UShort_t> chip;
    vector<UShort_t> channel;
    vector<vector<double>> wave; //一次event有多个hit，这里保存每个hit的q/t/ped/wave
    vector<double> charge;       //用Analysis计算这个Q/T
    vector<double> time;
    vector<double> pedeMean;
    vector<double> pedeRms;

    vector<pair<double, double>> hit; //这是mapping后的探测器阳极板对应的channel号（不是FEE电子学的channel号）
    
    //将hit进行分簇branch
    vector<vector<BeamHit>> branch;  //Pad读出：将一次event的hit分簇，branch.size()就是每个的cluster_size
    vector<vector<BeamHit>> XYbranch[2]; //条读出：X/Y的hit分簇，branch.size()就是每个的cluster_size, Xbranch给出X坐标，所以表示Y的second为-999

    //将branch按重心得到击中信息
    vector<RealBeamHit> cluster; //Pad读出，分簇后按照重心得到的三维真实坐标的击中信息，暂时没用
    vector<RealBeamHit> XYcluster[2]; //条读出，分簇后按照重心得到的三维真实坐标的击中信息

    int GetHitSize() { return board.size(); }
    int GetClusterSize() { return branch.size(); }

    void Init(int _ev, int _id)
    {
        id = _id;
        event = _ev;
        board.clear();
        chip.clear();
        channel.clear();
        wave.clear();
        charge.clear();
        time.clear();
        pedeMean.clear();
        pedeRms.clear();
        hit.clear();
        branch.clear();
        XYbranch[0].clear();
        XYbranch[1].clear();
        cluster.clear();
        XYcluster[0].clear();
        XYcluster[1].clear();
    }

    void Analysis(int pedRangeMin, int pedRangeMax, int sigRangeMin, int sigRangeMax) //分析每个hit的Q/T
    {
        charge.clear();
        charge.resize(board.size());
        time.clear();
        time.resize(board.size());

        for (int i = 0; i < (int)board.size(); i++)
        {
            time[i] = 0;

            //扣除pedstal
            //double pedstal = pedeMean[i];
            double pedstal = 0;
            for (int j = pedRangeMin; j < pedRangeMax; j++) //用0～150的时间
                pedstal += wave[i][j];
            pedstal /= (pedRangeMax - pedRangeMin);

            charge[i] = wave[i][0] - pedstal;
            for (int j = sigRangeMin; j < sigRangeMax; j++)
                if (charge[i] < wave[i][j] - pedstal)
                {
                    charge[i] = wave[i][j] - pedstal;
                    time[i] = j;
                }
        }
    }
};

//-----------------------------
// hit击中信息合并处理后的 combine-DST-root 的数据结构
class MyBeamTestHitData
{
public:
    int event = -1; //触发事例号
    vector<RealBeamHit> hit;
    vector<MyBeamTestData *> detector;

    void Init(int ev)
    {
        event = ev;
        hit.clear();
        detector.clear();
    }
};

//-----------------------------
// 探测器的基类
class MyBeamTestDetector
{
public:
    MyBeamTestDetector(int _id, int nsX, int nsY, double wsX, double wsY, double la = 0)
    {
        id = _id;
        Nstrp[0] = nsX;
        Nstrp[1] = nsY;
        strip[0] = wsX;
        strip[1] = wsY;
        LtoAnode = la;
    }
    ~MyBeamTestDetector()
    {
        if (fFile == 0)
            return;
        fFile->Close();
    }

    //基本信息：id，类型，坐标，旋转，尺寸
    void SetType(int tmp) { type = tmp; }
    void SetDetID(int idd) { id = idd; }
    void SetName(TString nm) { name = nm; }
    void SetFullSize(double a, double b, double c)
    {
        siz[0] = a;
        siz[1] = b;
        siz[2] = c;
    }
    void SetPosition(double a, double b, double c)
    {
        pos[0] = a;
        pos[1] = b;
        pos[2] = c;
    }
    void SetRotation(double a, double b)
    {
        rot[0] = a;
        rot[1] = b;
        ROT[0] = a / 180. * 3.1415926;
        ROT[1] = b / 180. * 3.11415926;
    }

    int GetType() { return type; }
    int GetDetID() { return id; }
    TString GetName() { return name; }
    double *GetFullSize() { return siz; }
    double *GetPosition() { return pos; }
    double *GetRotation() { return rot; }
    double *GetROTATION() { return ROT; }

    //
    // 板号和芯片号设置
    vector<int> chipList;
    vector<int> boardList;
    int GetNChip() { return (int)chipList.size(); }
    int GetNBoard() { return (int)boardList.size(); }
    void SetNChip(int nc) { chipList.resize(nc); }
    void SetNBoard(int nb) { boardList.resize(nb); }
    vector<int> GetChipList() { return chipList; }
    vector<int> GetBoardList() { return boardList; }
    void SetChip(int ic, int nc) { chipList[ic] = nc; }
    void SetBoard(int ib, int nb) { boardList[ib] = nb; }

    //
    //电子学mapping
    vector<vector<vector<double>>> mapX;
    vector<vector<vector<double>>> mapY;

    virtual void ReadMapFile() = 0;
    void ResizeMap()
    {
        mapX.resize(boardList.size());
        for (int i = 0; i < (int)boardList.size(); i++)
        {
            mapX[i].resize(chipList.size());
            for (int j = 0; j < (int)chipList.size(); j++)
                mapX[i][j].resize(64);
        }

        mapY.resize(boardList.size());
        for (int i = 0; i < (int)boardList.size(); i++)
        {
            mapY[i].resize(chipList.size());
            for (int j = 0; j < (int)chipList.size(); j++)
                mapY[i][j].resize(64);
        }
    }

    //
    // pedestal的读取和设置
    vector<double> pedMean;
    vector<double> pedRMS;
    void SetNPedMean() { pedMean.resize(boardList.size() * chipList.size() * 64); }
    void SetNPedRMS() { pedRMS.resize(boardList.size() * chipList.size() * 64); }
    void SetPedMean(int ibd, int ichip, int ichan, double nm) { pedMean[CalID(ibd, ichip, ichan)] = nm; }
    void SetPedRMS(int ibd, int ichip, int ichan, double nr) { pedRMS[CalID(ibd, ichip, ichan)] = nr; }
    double GetPedMean(int board, int chip, int chan) { return pedMean[CalID(IndexBD(board), IndexChip(chip), chan)]; }
    double GetPedRMS(int board, int chip, int chan) { return pedRMS[CalID(IndexBD(board), IndexChip(chip), chan)]; }

    //
    // 纯虚函数：作图
    virtual void DrawHits(TGeoManager *geom, TGeoVolume *stop, TGeoMedium *med, MyBeamTestHitData *hit) = 0;
    virtual void DrawDetector(TGeoManager *geom, TGeoVolume *stop, TGeoMedium *med) = 0;
    virtual pair<double, double> MapPosition(int bd, int chip, int channel) = 0;
    //AbsPosition计算实验室系下的绝对坐标。0点为束流点，束流方向为+Z，垂直向上方向为+Y，左手系可得+X
    //                           ^ +Y
    //                           |   /| +X
    //                           |  /
    //                           | /
    //   +Z <--------------------|/
    virtual void AbsPosition(pair<double, double> epos, double *absPos) = 0;

    //GeoPosition计算在GeoMetry下的坐标。0点为束流点，束流方向为+X，垂直向下方向为+Z，右手系可得+Z
    //   +X <--------------------|
    //                         / |
    //                        /  |
    //                    +Z|/    |
    //                          \ / +y
    double *GeoPosition(double *absPos) //将实验室系坐标转换为TGeometry的坐标
    {
        static double geoPos[3];
        geoPos[0] = absPos[2];
        geoPos[1] = -absPos[1];
        geoPos[2] = -absPos[0];
        return geoPos;
    }

    virtual bool FillDistribution(int i, MyBeamTestHitData *hit) = 0;
    //virtual void FillHitVector(vector<vector<double>> hitpos, MyBeamTestHitData *hitDat) = 0;

public:
    int CalID(int ibd, int chip, int chan) { return ibd * chipList.size() * 64 + chip * 64 + chan; }
    int IndexBD(int board)
    {
        int index = -1;
        for (int i = 0; i < (int)boardList.size(); i++)
            if (board == boardList[i])
                index = i;
        return index;
    }
    int IndexChip(int chip)
    {
        int index = -1;
        for (int i = 0; i < (int)chipList.size(); i++)
            if (chip == chipList[i])
                index = i;
        return index;
    }

public:
    int id; //此探测器的id在整个系统里的编号，一些图可能会根据这个id号来进行计算编号。
    int type;
    TString name;
    double siz[3];
    double pos[3];
    double rot[2];
    double ROT[2];
    int Nstrp[2];        //实验室系下X/Y方向的strip条数
    double strip[2];     //实验室系下X/Y方向的strip条宽，[0]为平行X轴的条的宽度，[1]为平行与Y轴的条的宽度
    double LtoAnode = 0; //阳极板与轴心的距离

    //DST数据文件指针
    TFile *fFile = 0;
    TTree *fTree = 0;
    MyBeamTestData *fEvent = 0;
    Long64_t nentries = 0;

    //数据相关
    int iEntry = 0; //current_id for TTree's entry
    int GetEntries() { return (fTree == 0) ? 0 : fTree->GetEntries(); }
    int GetFirstTrig() { return GetTrigID(0); };
    int GetLastTrig() { return GetTrigID(nentries - 1); };
    int GetTrigID(int ientry)
    {
        if (fTree == 0 || ientry < 0 || ientry > nentries - 1)
            return -1;

        fTree->GetEntry(ientry);
        return fEvent->event;
    }

    bool SearchTrigID(int trigid, MyBeamTestHitData *fEventList)
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
};

#endif