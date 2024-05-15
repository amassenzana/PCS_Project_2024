#include "DFNLibrary.hpp"
#include <fstream>
#include <list>

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
    cout << center << endl;
}

void DFNLibrary::Frattura::computeRadius(){
    double max = 0, t = 0;
    for(int i = 0; i < numVert; i++){
        t = (center-vertices[i]).norm();
        if(t>max)
            max = t;
    }

    radius = max;
    cout << radius << endl;
}
