#ifndef MyBeamTestDetector_h
#define MyBeamTestDetector_h

#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "MyBeamTestDatStruct.h"

class MyBeamTestRICH;
class MyBeamTestTrackAGET;
class MyBeamTestTrackVMM;

//-----------------------------
// 探测器的基类
class MyBeamTestDetector
{
public:
    MyBeamTestDetector(int);
    ~MyBeamTestDetector();

public:
    //
    //___ 设置基本信息：id，类型，坐标，旋转，尺寸 ___
    TString GenerateIntro();
    void SetParameters(TString, vector<TString>);
    void SetDetID(int idd) { id = idd; }
    void SetNXNY(int nx, int ny) { Nstrp[0] = nx; Nstrp[1] = ny; }
    void SetWXWY(double wx, double wy, double ltoa) { strip[0] = wx; strip[1] = wy; LtoAnode = ltoa; }
    void SetName(TString nm) { name = nm; }
    void SetSurName(TString nm) { surname = nm; }
    void SetMapping(TString mp) { mapping = mp; }
    void SetFullSize(double a, double b, double c) { siz[0] = a; siz[1] = b; siz[2] = c; }
    void SetPosition(double a, double b, double c) { pos[0] = a; pos[1] = b; pos[2] = c; }
    void SetRotation(double a, double b) { rot[0] = a; rot[1] = b; ROT[0] = a / 180. * 3.1415926; ROT[1] = b / 180. * 3.11415926; }

    int GetDetID() { return id; }
    TString GetName() { return name; }
    double *GetFullSize() { return siz; }
    double *GetPosition() { return pos; }
    double *GetRotation() { return rot; }
    double *GetROTATION() { return ROT; }

    //
    //___ 板号和芯片号设置 ___
    int GetNChip() { return (int)chipList.size(); }
    int GetNBoard() { return (int)boardList.size(); }
    void SetNChip(int nc) { chipList.resize(nc); }
    void SetNBoard(int nb) { boardList.resize(nb); }

    vector<int> GetChipList() { return chipList; }
    vector<int> GetBoardList() { return boardList; }
    void SetChip(int ic, int nc) { chipList[ic] = nc; }
    void SetBoard(int ib, int nb) { boardList[ib] = nb; }

    //
    //___ 电子学mapping ___
    virtual void ResizeMap() = 0;
    virtual void ReadMapFile() = 0;

    int CalID(int ibd, int chip, int chan);
    int IndexBD(int board);
    int IndexChip(int chip);

    //
    //___ pedestal的读取和设置 ___
    void SetNPedMean() { pedMean.resize(boardList.size() * chipList.size() * 64); }
    void SetNPedRMS() { pedRMS.resize(boardList.size() * chipList.size() * 64); }
    void SetPedMean(int ibd, int ichip, int ichan, double nm) { pedMean[CalID(ibd, ichip, ichan)] = nm; }
    void SetPedRMS(int ibd, int ichip, int ichan, double nr) { pedRMS[CalID(ibd, ichip, ichan)] = nr; }
    double GetPedMean(int board, int chip, int chan) { return pedMean[CalID(IndexBD(board), IndexChip(chip), chan)]; }
    double GetPedRMS(int board, int chip, int chan) { return pedRMS[CalID(IndexBD(board), IndexChip(chip), chan)]; }
    virtual void GenPedMap() = 0;
    
    //
    //___ 纯虚函数：作图 ___
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

    //
    //___ DST数据操作___
    int GetEntries();
    int GetFirstTrig();
    int GetLastTrig();
    int GetTrigID(int ientry);
    bool SearchTrigID(int trigid, MyBeamTestHitData *fEventList);

public:
    int id = -1; //此探测器的id在整个系统里的编号，一些图可能会根据这个id号来进行计算编号。
    TString name;
    TString surname;     //在log里的名字 
    TString mapping;     //mapping的文件路径
    double siz[3] = {0, 0, 0};
    double pos[3] = {0, 0, 0};
    double rot[2] = {0, 0};       //in degree
    double ROT[2] = {0, 0};       //in rad
    int Nstrp[2]  = {1, 1};       //实验室系下X/Y方向的strip条数或pad数
    double strip[2] = {1, 1};     //实验室系下X/Y方向的strip条宽，[0]为平行X轴的条的宽度，[1]为平行与Y轴的条的宽度
    double LtoAnode = 0;          //阳极板与探测器轴心的距离

    vector<int> chipList;         //电子学芯片
    vector<int> boardList;        //电子学板号
    vector<vector<vector<double>>> mapX;    //X mapping
    vector<vector<vector<double>>> mapY;    //Y mapping

    vector<double> pedMean;
    vector<double> pedRMS;

    TGeoMatrix* geoMatrix;

    //DST数据文件指针
    TFile *fFile = 0;
    TTree *fTree = 0;
    MyBeamTestData *fEvent = 0;
    Long64_t nentries = 0;
    int iEntry = 0; //current_id for TTree's entry

};

#endif