//----------
//分析数据控制参数
TString fName("RUN86"); //./DESY-BeamTest/RUN64/Combine/Combined-dst.root
int processControl = 7; //1: 用hitmapReal的直方图去做Offset，运行几遍直到输出的offset不怎么变化，同时检查二维直方图，确保束流信号在0，0位置，！！！需要手动更新offset数据！！！
                        //2: 用T02/T03去看其他tracker，然后更新offset。运行几遍直到输出的offset不怎么编号，同时检查残差分布图， ！！！需要手动更新offset数据！！！
                        //3: 看RICH的hitmap及其他相关的分布, 这里会cut选择束流击中位置，从而来做RICH的offset修正，因此需要调整一下在IsRICHClusterEffective里的参数，并运行几遍 ！！！需要手动更新offset数据！！！
                        //4: 修正完后再看RICH的分布，此时不做任何cut，用来检查原始数据过offset修正后是否正常
                        //5: RICH的offset的残差分布，和3一样，也会cut选择束流击中位置，多运行几遍 ！！！需要手动更新offset数据！！！
                        //6: 修正完后再看RICH的分布，此时对cluster进行cut
                        //7: 重建光子

//这里的offset是探测器坐标系的零点，在束流实验系下的坐标位置。
//束流实验系是束流为0点，方向为z，向上为y，向东为x
//实验记录上的角度是和探测器表面的法线的夹角！！！
vector<vector<double>> offset = {
    {63, -16.92, -5.609, -4.664, 3.735, 0, 0}, // X
    {63, 24.6, -9.753, -9.218, -6.235, 0, 0},  // Y
    {63, 375, 0, 97, 875, 0, 0},               // Z
    {63, 45},                                  // angle

    {65, 0, 6.796, 7.72, 17.33, 0, 0}, // X
    {65, 0, 7.08, 8.197, 14.76, 0, 0}, // Y
    {65, 375, 0, 97, 875, 0, 0},       // Z
    {65, 45},                          // angle

    {86, -14.21, 0.203, 0.86, -2.287, 0, 0},  // X
    {86, 31.6, 0.195, -0.2778, -5.591, 0, 0}, // Y
    {86, 375, 0, 97, 875, 0, 0},              // Z
    {86, 30},                                 // angle

    {999, 0, 0, 0, 0, 0, 0},      //X
    {999, 0, 0, 0, 0, 0, 0},      //Y
    {999, 375, 0, 97, 875, 0, 0}, //Z
    {999, 45}                     //angle
};

const int NDET = 6;
double Xoff[NDET] = {0, 0, 0, 0, 0, 0};
double Yoff[NDET] = {0, 0, 0, 0, 0, 0};
double Zpos[NDET] = {375, 0, 97, 875, 0, 0};
double angle = 45;

long DSTEntries = 1000;                  //设置为小于100的数，就跑DSTEntries个事例，设置为大与100的数，就跑数据里的所有数据。
int verbose = (DSTEntries < 100) ? 1 : 0; //输出信息后可以在mathematicas里检查
