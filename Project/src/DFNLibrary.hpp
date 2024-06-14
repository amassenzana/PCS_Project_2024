#pragma once

#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <variant>


namespace DFNLibrary{
struct Point {
    double x, y;
};
struct Frattura{
    int id;
    int numVert;
    std::vector<Eigen::Vector3d> vertices;
    Eigen::Vector3d center;
    double radius;
    Eigen::Vector3d planeC;
    double planeD;
    void computeCenter();
    void computeRadius();
    void computePlane();
    std::vector<Eigen::Vector3d> IntersectionPoints;
    std::vector<Eigen::Vector2d> GoodValsDouble;
    std::vector<std::vector<Eigen::Vector3d>> GoodValsVec; // "Tensore"
    // Eigen::Vector3d PlaneDir;
};

struct Traccia{
    Eigen::Vector3d origin, end;
    Frattura* generator[2];
    bool passante[2];
};

struct DFN{
    int numberFractures;
    std::vector<Frattura> fractures;
    std::vector<Traccia> tracce;
    void computeDFN();
};

bool importDFN(std::string path, DFN& dfn);
bool checkIntersection(Frattura& f1, Frattura& f2, std::vector<Eigen::Vector3d>& v);
int sign(double d);
Eigen::Vector2d ProjXY( Eigen::Vector3d tmpvec_xy);
Eigen::Vector2d ProjXZ( Eigen::Vector3d tmpvec_xz);
Eigen::Vector2d ProjYZ( Eigen::Vector3d tmpvec_yz);
std::pair<Eigen::Vector2d, bool> ProjIntersection(Eigen::Vector2d a1, Eigen::Vector2d a2, Eigen::Vector2d b1, Eigen::Vector2d b2);
Eigen::Vector3d IntersectVerif(Eigen::Vector3d p1, Eigen::Vector3d d1, Eigen::Vector3d p2, Eigen::Vector3d d2);
bool areParallel(DFNLibrary::Point line1_p1, DFNLibrary::Point line1_p2, DFNLibrary::Point line2_p1, DFNLibrary::Point line2_p2);
Eigen::Vector3d IntersecPointCalc(int i, int j);
} // namespace DFNLibrary


