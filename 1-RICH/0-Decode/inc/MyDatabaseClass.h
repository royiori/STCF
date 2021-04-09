#ifndef _MyDatabaseClass_h_
#define _MyDatabaseClass_h_

#include <map>
#include <vector>

using namespace std;

//-----------------------------------------------------------------------------------------------
//
class MyDatabaseClass : public TObject
{
public:
    MyDatabaseClass();
    virtual ~MyDatabaseClass();

    void reset()
    {
        evtID = 0;
        timeStamp[0] = 0;
        timeStamp[1] = 0;
        for(size_t i=0; i<wave.size(); i++)
            wave[i].clear();
            
        t0.clear();
        //...
        eoeCnt = 0;
    }

public:
    int evtID = 0;
    long int timeStamp[2] = {0, 0};

    vector<vector<int>> wave;
    vector<int> t0;
    int t0step = 0;

    vector<int> xyFlag;
    vector<int> xyID;

    int eoeCnt = 0;

    ClassDef(MyDatabaseClass,1)
};

#endif
