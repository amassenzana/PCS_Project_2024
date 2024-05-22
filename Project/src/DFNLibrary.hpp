#pragma once

#include "Eigen/Eigen"
#include <iostream>
#include <vector>


namespace DFNLibrary{

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
bool checkIntersection(Frattura& f1, Frattura& f2);
int sign(double d);
} // namespace DFNLibrary


