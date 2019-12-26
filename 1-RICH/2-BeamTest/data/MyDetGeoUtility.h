// Modified from MyDetescope - EUTelutility.h
//
#ifndef MyDetGeoUtility_h
#define MyDetGeoUtility_h

// ROOT
#include "TVectorD.h"

// system includes <>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace Utility
{

//FIXME: check doubling with existing "AlignmentMode" enum!
enum class alignMode
{
    XYShiftsRotZ,
    XYShifts,
    XYShiftsAllRot,
    XYZShiftsRotZ,
    XYZShiftsRotXYZ
};

const long double PI = 3.141592653589793238L;

TMatrixD setPrecision(TMatrixD mat, double mod);

class Hit;

/**
     * Vector of Hits in plane
     */
typedef std::vector<Hit> HitsVec;
typedef std::vector<Hit *> HitsPVec;

/**
     * Represents geometric properties of hits
     */
class Hit
{
public:
    Hit() : _PlaneID(-1), _x(0.0), _y(0.0), _z(0.0) {}

    Hit(int id, double x, double y, double z)
        : _PlaneID(id), _x(x), _y(y), _z(z) {}

    bool operator<(const Hit &b) const { return (_z < b.Z()); }

    inline void SetZ(double z) { this->_z = z; }

    inline double Z() const { return _z; }

    inline void SetY(double y) { this->_y = y; }

    inline double Y() const { return _y; }

    inline void SetX(double x) { this->_x = x; }

    inline double X() const { return _x; }

    void SetPlaneID(int _PlaneID) { this->_PlaneID = _PlaneID; }

    int GetPlaneID() const { return _PlaneID; }

protected:
    int _PlaneID;

    double _x;
    double _y;
    double _z;
};

/**
     * Specify a Rectangular in a sensor
     */
class SensorRectangular
{
protected:
    int sensor; // SensorID
    int A;      // lowest pixel in X direction
    int B;      // lowest pixel in Y direction
    int C;      // highest pixel in X direction
    int D;      // highest pixel in Y direction
public:
    SensorRectangular(int s, int a, int b, int c, int d)
        : sensor(s), A(a), B(b), C(c), D(d){};

    SensorRectangular() : sensor(0), A(0), B(0), C(0), D(0){};

    int getSensor() const { return sensor; }

    // look if x and y are inside the foreseen rectangular

    bool isInside(int x, int y) const
    {
        return (x >= A && x <= C && y >= B && y <= D);
    }
};

class RectangularArray
{
protected:
    std::map<int, SensorRectangular> _rect;

public:
    void addRectangular(SensorRectangular &s) { _rect[s.getSensor()] = s; }

    bool isInside(int s, int x, int y)
    {
        std::map<int, SensorRectangular>::iterator it = _rect.find(s);
        if (it == _rect.end())
        {   // not in the map means no limit on this sensor
            // -> always true
            return true;
        }
        SensorRectangular cSensor = _rect[s];
        return cSensor.isInside(x, y);
    }
};

/*!
     * Fills indices of not excluded planes
     */

void FillNotExcludedPlanesIndices(std::vector<int> &,
                                  const std::vector<unsigned int> &,
                                  unsigned int = 0);

/*!
     * Checks if a hit belongs to a hot-pixel-map
     */

//! Called for first event per run
/*! Reads hotpixel information from hotPixelCollection into hotPixelMap
     * to be used in the sensor exclusion area logic
     */

int cantorEncode(int X, int Y);
std::map<int, std::vector<int>>
readNoisyPixelList(LCEvent *event,
                   std::string const &noisyPixelCollectionName);

Eigen::Matrix3d rotationMatrixFromAngles(long double alpha,
                                         long double beta,
                                         long double gamma);
Eigen::Vector3d getRotationAnglesFromMatrix(Eigen::Matrix3d rotMat);

std::unique_ptr<EUTelClusterDataInterfacerBase>
getClusterData(IMPL::TrackerDataImpl *const data, SparsePixelType type);
std::unique_ptr<EUTelClusterDataInterfacerBase>
getClusterData(IMPL::TrackerDataImpl *const data, int type);

std::unique_ptr<EUTelTrackerDataInterfacer>
getSparseData(IMPL::TrackerDataImpl *const data, SparsePixelType type);
std::unique_ptr<EUTelTrackerDataInterfacer>
getSparseData(IMPL::TrackerDataImpl *const data, int type);

std::map<std::string, bool>
FillHotPixelMap(EVENT::LCEvent *event,
                const std::string &hotPixelCollectionName);

bool HitContainsHotPixels(const IMPL::TrackerHitImpl *hit,
                          const std::map<std::string, bool> &hotPixelMap);

std::unique_ptr<EUTelVirtualCluster>
GetClusterFromHit(const IMPL::TrackerHitImpl *);

int getSensorIDfromHit(EVENT::TrackerHit *hit);

/** Highland's formula for multiple scattering */
double getThetaRMSHighland(double, double);

/** Calculate median */
double getMedian(std::vector<double> &);

/** Calculate track's 2D curvature */
double getCurvature(double, double, double);

/** Solve quadratic equation a*x^2 + b*x + c = 0 */
std::vector<double> solveQuadratic(double, double, double);

/** getClusterSize method from TrackerHit object:: assumes known cluster
     * types */
void getClusterSize(const IMPL::TrackerHitImpl *hit, int &sizeX,
                    int &sizeY);

/** */
void copyLCCollectionHitVec(LCCollectionVec *, LCCollectionVec *);

/** */
void copyLCCollectionTrackVec(LCCollectionVec *, LCCollectionVec *);

} // namespace Utility

#endif