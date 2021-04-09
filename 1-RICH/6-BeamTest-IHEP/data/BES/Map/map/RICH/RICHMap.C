#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

void RICHMap()
{
    fstream mapFp1; // 探测器map
    fstream mapFp2; // FEE的map
    fstream mapFp3; // 探测器和FEE的链接
    string a;

    mapFp1.open("./map-Detector", ios::in);
    mapFp2.open("./map-FEE", ios::in);
    mapFp3.open("./map-Connector", ios::in);

    //读入connector的map
    int board, hiro1, hiro2;
    vector<int> boardlist;
    vector<int> hiro1list;
    vector<int> hiro2list;
    getline(mapFp3, a);
    while (!mapFp3.eof())
    {
        mapFp3 >> board >> hiro1 >> hiro2;
        boardlist.push_back(board);
        hiro1list.push_back(hiro1);
        hiro2list.push_back(hiro2);
    }

    //读入探测器map
    int hiro, pin, posx, posy;
    vector<int> hirolist;
    vector<int> pinlist;
    vector<int> xlist;
    vector<int> ylist;
    getline(mapFp1, a);
    while (!mapFp1.eof())
    {
        mapFp1 >> hiro >> pin >> posx >> posy;
        hirolist.push_back(hiro);
        pinlist.push_back(pin);
        xlist.push_back(posx);
        ylist.push_back(posy);
    }

    //读入FEE的map
    int chip, channel, pin2;
    vector<int> chiplist;
    vector<int> channellist;
    vector<int> pin2list;
    getline(mapFp2, a);
    while (!mapFp2.eof())
    {
        mapFp2 >> chip >> channel >> pin2;
        chiplist.push_back(chip);
        channellist.push_back(channel);
        pin2list.push_back(pin2);
    }

    for (int i = 0; i < (int)xlist.size(); i++)
    {
        //根据hiro的接头编号找board
        int Board = -1;
        int chipFlag = -1;
        for (int j = 0; j < (int)boardlist.size(); j++)
        {
            if (hirolist[i] == hiro1list[j])
            {
                chipFlag = 1;
                Board = boardlist[j];
            }
            if (hirolist[i] == hiro2list[j])
            {
                chipFlag = 2;
                Board = boardlist[j];
            }
        }

        //根据hiro pin的编号找chip 和 channel
        int Chip = -1;
        int Channel = -1;
        for (int j = 0; j < (int)pin2list.size(); j++)
        {
            if (pinlist[i] == pin2list[j])
            {
                Chip = chiplist[j];
                Channel = channellist[j];
                Chip = (chipFlag % 2 == 1) ? Chip + 10 : Chip + 12;
                break;
            }
        }

        if (Board == -1 || Chip == -1 || Channel == -1)
        {
            cout << "failed to find the map for HIRO=" << hirolist[i] << " PIN=" << pinlist[i] << endl;
            return;
        }

        cout << Board << "\t" << Chip << "\t" << Channel << "\t" << xlist[i] << "\t" << ylist[i] << endl;
    }
}
