#pragma once

#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <utility>
#include <list>


namespace DFNLibrary{


struct Traccia{
    int id;
    Eigen::Vector3d origin, end;
    int idF1, idF2;
    double length();
};

struct Frattura{
    int id;
    int numVert;
    std::vector<Eigen::Vector3d> vertices;
    Eigen::Vector3d center;
    double radius;
    Eigen::Vector3d planeC;
    double planeD;
    std::vector<std::pair<Traccia*, bool>> tracce;

    void computeCenter();
    void computeRadius();
    void computePlane();
};





struct DFN{
    int numberFractures;
    std::vector<Frattura> fractures;
    std::list<Traccia> tracce;
    void computeDFN();
    void plotFracture();
    void elaborateV(Frattura& f1, Frattura& f2, std::vector<Eigen::Vector3d>& v);
    void output();
};

bool importDFN(std::string path, DFN& dfn);
bool checkIntersection(Frattura& f1, Frattura& f2, std::vector<Eigen::Vector3d>& v);
int sign(double d);
bool pointSort(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);
bool compareTrace(std::pair<Traccia*, bool>& T1, std::pair<Traccia*, bool>& T2);
} // namespace DFNLibrary


