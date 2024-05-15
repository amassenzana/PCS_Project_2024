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
    void computeCenter();
    void computeRadius();
};

struct Traccia{
    Eigen::Vector3d origin, end;
    Frattura* generator[2];
};

struct DFN{
    int numberFractures;
    std::vector<Frattura> fractures;
    std::vector<Traccia> tracce;
};

bool importDFN(std::string path, DFN& dfn);


} // namespace DFNLibrary


