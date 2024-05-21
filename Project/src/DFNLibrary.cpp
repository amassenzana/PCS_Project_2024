#include "DFNLibrary.hpp"
#include <fstream>
#include <list>
#define tol 1e-10

using namespace std;

bool DFNLibrary::importDFN(std::string path, DFN &dfn){
    path = "DFN/" + path;
    ifstream file(path);

    if(file.fail()){
        cerr << "Errore nell'apertura del file" << endl;
        return false;
    }

    string line;
    getline(file, line);  // BRUCIO # Number of Fractures
    getline(file, line);
    dfn.numberFractures = stoi(line);
    dfn.fractures.reserve(dfn.numberFractures);

    for(size_t k = 0; k < dfn.numberFractures; k++){
        Frattura f;
        getline(file, line); // BRUCIO # FractureId; NumVertices
        getline(file, line, ';');
        f.id = stoi(line);
        getline(file, line);
        f.numVert = stoi(line);
        f.vertices.resize(f.numVert);
        getline(file, line); // BRUCIO # Vertices

        vector<string> lines(3);

        for(size_t i = 0; i < 3; i++){
            getline(file, lines[i]);
            for(size_t j = 0; j < f.numVert-1; j++){
                int pos = lines[i].find(';',0);
                f.vertices[j](i) = stod(lines[i].substr(0,pos));
                lines[i] = lines[i].substr(pos+1, string::npos);
            }
            f.vertices[f.numVert-1](i) = stod(lines[i]);
        }
        f.computeCenter();
        f.computeRadius();
        f.computePlane();
        dfn.fractures.push_back(f);
    }


    return true;
}

void DFNLibrary::Frattura::computeCenter(){


    Eigen::Vector3d c = {0, 0, 0};
    double sw = 0;

    for(int i = 1; i < numVert+1; i++){
        double w =  (vertices[i%numVert]-vertices[(i+1)%numVert]).norm() +
                    (vertices[i%numVert]-vertices[(i-1)%numVert]).norm();
        c += vertices[i%numVert]*w;
        sw += w;
    }

    center[0] = c[0]/sw;
    center[1] = c[1]/sw;
    center[2] = c[2]/sw;
    // cout << center << endl;
}

void DFNLibrary::Frattura::computeRadius(){
    double max = 0, t = 0;
    for(int i = 0; i < numVert; i++){
        t = (center-vertices[i]).norm();
        if(t>max)
            max = t;
    }

    radius = max;
    // cout << radius << endl;
}







void DFNLibrary::Frattura::computePlane(){

    Eigen::Vector3d p1 = vertices[0];
    Eigen::Vector3d p2 = vertices[1];
    Eigen::Vector3d p3 = center;

    Eigen::Vector3d v1 = p2-p1;
    Eigen::Vector3d v2 = p3-p1;
    Eigen::Vector3d n = v1.cross(v2);

    Eigen::Vector3d plane{{n.x(), n.y(), n.z()}};
    planeC = plane;
    planeD = -n.dot(p1);
}

void DFNLibrary::DFN::computeDFN(){
    for(int i = 0; i < numberFractures; i++){
        for(int j = 0; j < numberFractures; j++){
            if(i!=j && checkIntersection(fractures[i], fractures[j]))
                cout << "Intersezione tra f[" << i << "] e f[" << j << "]" << endl;
        }
    }


}

bool DFNLibrary::checkIntersection(Frattura &f1, Frattura &f2){
    double distance = (f1.center-f2.center).norm();
    if(distance > f1.radius+f2.radius)
        return false;


    int s = sign(f1.planeC.dot(f2.vertices[0]) + f1.planeD);
    bool flag = false;

    for(int i = 1; i < f2.numVert; i++){
        if(!(sign(f1.planeC.dot(f2.vertices[i]) + f1.planeD) == s))
            flag = true;
    }
    if(!flag) return false;

    s = sign(f2.planeC.dot(f1.vertices[0]) + f2.planeD);
    flag = false;
    for(int i = 1; i < f1.numVert; i++){
        double k = f2.planeC.dot(f1.vertices[i]) + f2.planeD;
        if(! (k == s))
            flag = true;
    }
    if(!flag) return false;



    return true;
}

int DFNLibrary::sign(double d){
    if (d>tol) return 1;
    if (d<-tol) return -1;

    cout << "Careful" << endl;

    return 0;
}

