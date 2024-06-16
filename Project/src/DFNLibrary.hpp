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
    std::vector<std::pair<Traccia*, bool>> tracce;

    Eigen::Vector3d computeCenter();
    double computeRadius();
    double computePlane(Eigen::Vector3d& planeC);
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
double angle(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);
double angleReference (const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& center);
} // namespace DFNLibrary


