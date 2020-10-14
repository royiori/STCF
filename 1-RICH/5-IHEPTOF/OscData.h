#include <iostream>
#include <vector>
#define CHMAX 4

class MyBeamTestDetector
{
public:
    MyBeamTestDetector() { ; }

    ~MyBeamTestDetector() { ; }


    void PushData(int nch, vector<double> data)
    {
        if (nch < 0 || nch >= CHMAX)
            return;
        fData[nch] = data;
    }

    void PushTime(vector<double> time) { fTime = time; }

public:
    int NChannel;
    int NPoints;

    vector<double> fTime;
    vector<double> fData[CHMAX];
};