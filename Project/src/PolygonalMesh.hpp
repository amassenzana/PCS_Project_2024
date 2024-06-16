#ifndef __POLYGONAL__MESH_H
#define __POLYGONAL__MESH_H

#include <map>
#include <vector>
#include <list>
#include "Eigen/Eigen"

namespace PolygonalLibrary{

struct PolygonalMesh{
    unsigned int NumberCell0D = 0;
    std::vector<unsigned int> Cell0DId = {};
    std::vector<Eigen::Vector2d> Cell0DCoordinates = {};
    std::map<unsigned int, std::list<unsigned int>> Cell0DMarkers = {};

    unsigned int NumberCell1D = 0;
    std::vector<unsigned int> Cell1DId = {};
    std::vector<Eigen::Vector2i> Cell1DVertices = {};
    std::map<unsigned int, std::list<unsigned int>> Cell1DMarkers = {};

    unsigned int NumberCell2D = 0;
    std::vector<unsigned int> Cell2DId = {};
    std::vector<std::vector<unsigned int>> Cell2DVertices = {};
    std::vector<std::vector<unsigned int>> Cell2DEdges = {};
};



} // namespace PolygonalMesh











#endif
