
#ifndef MyBeamTestDatStructure_h
#define MyBeamTestDatStructure_h

using namespace std;

enum
{
    BIN,
    RAW,
    DST,
    PED
};

enum
{
    RICH,
    TrackerAGET,
    TrackerVMM,
    ALL
};

//-----------------------------
// DST-root里hit击中信息合并处理后的击中数据结构
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
// dst.root 的数据结构，包含的是每个hit及cluster的channel坐标，而不是绝对坐标
class MyBeamTestData
{
public:
    int id = -1;    //探测器的id号
    int event = -1; //触发事例号

    //1. 以下为raw.root里的原始数据，vector长度为一个event里的有效信号数
    vector<UShort_t> board;
    vector<UShort_t> chip;
    vector<UShort_t> channel;
    vector<vector<double>> wave; //wave可以通过gMyGuiBeamTest里的SaveWaveFlag来决定是否保存到dst.root里
    //2. 以下为用Analysis计算这个charge/time, vector长度为一个event里的有效信号数
    vector<double> charge; //charge=qmax
    vector<double> time;   //time=charge@qmax
    vector<double> pedeMean;
    vector<double> pedeRms;
    vector<pair<double, double>> hit; //这是mapping后的探测器阳极板对应的channel号（不是FEE电子学的channel号）

    //将hit进行分簇branch
    vector<vector<BeamHit>> branch;      //Pad读出：将一次event的hit分簇，branch.size()就是每个的cluster_size
    vector<vector<BeamHit>> XYbranch[2]; //条读出：X/Y的hit分簇，branch.size()就是每个的cluster_size, Xbranch给出X坐标，所以表示Y的second为-999

    //将branch按重心得到击中信息
    vector<RealBeamHit> cluster;      //Pad读出，分簇后按照重心得到的三维真实坐标的击中信息，暂时没用
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
            double max_charge = charge[i] + pedstal;
            //charge[i] = wave[i][0] - pedstal;
            for (int j = sigRangeMin; j < sigRangeMax; j++)
            {
                if ((wave[i][j] - pedstal) > 0.15 * (max_charge - pedstal))
                {
                    time[i] = j;
                    //printf ("Pedstal = %5.4f, max_charge = %5.4f, j = %d \n", pedstal, max_charge, j);
                    break;
                }
            }
        }
    }
};

//-----------------------------
// hit击中信息合并处理后的 combine.root 的数据结构
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

#endif