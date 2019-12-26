#ifndef MyDetGeoManager_h
#define MyDetGeoManager_h


// C++
#include <array>
#include <map>
#include <memory>
#include <string>

// Eigen
#include <Eigen/Core>

// ROOT
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

static const double PI = TMath::Pi();
static const double DEG = 180. / PI;
static const double RADIAN = PI / 180.;

class MyDetGeoManager
{
public:
    bool isGeoInitialized = false;

    /* 探测器数目 */
    static unsigned ndet;
    /* 探测器种类 */
    std::map<int, TString> detTypeMap;

    std::vector<std::unique_ptr<EUTelLayer>> _telescopeLayers;
    /* 每个平面的法线方向 */
    std::map<int, Eigen::Vector3d> planeNormMap;
    /* 每个平面的X轴 */
    std::map<int, Eigen::Vector3d> planeXmap;
    /* 每个平面的Y轴 */
    std::map<int, Eigen::Vector3d> planeYMap;
    /* 每个平面的辐射长度 */
    std::map<int, double> planeRadMap;

    MyDetGeoManager();
}


#endif