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
        f.id = stoi(line);      // Salvo FractureId dal file in f.id
        getline(file, line);
        f.numVert = stoi(line);
        f.vertices.resize(f.numVert);   // Salvo NumVertices dal file in f.numVert
        getline(file, line); // BRUCIO # Vertices

        vector<string> lines(3);

        for(size_t i = 0; i < 3; i++){
            getline(file, lines[i]);
            istringstream converter(lines[i]);          // int pos = lines[i].find(';',0);
            for(size_t j = 0; j < f.numVert; j++){
                char tmp;                               // f.vertices[j](i) = stod(lines[i].substr(0,pos));
                converter >> f.vertices[j][i] >> tmp;   // lines[i] = lines[i].substr(pos+1, string::npos);
            }
            // f.vertices[f.numVert-1](i) = stod(lines[i]);
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
        for(int j = i+1; j < numberFractures; j++){
            vector<Eigen::Vector3d> v;
            if(i!=j && checkIntersection(fractures[i], fractures[j], v)){
                cout << "Intersezione tra f[" << i << "] e f[" << j << "]:" << endl;

            }
        }
    }


}

bool DFNLibrary::checkIntersection(Frattura &f1, Frattura &f2, vector<Eigen::Vector3d> &v){

    // TEST 1: La distanza tra f1 e f2 è maggiore della somma dei raggi
    //  delle ipersfere che le contengono
    double distance = (f1.center-f2.center).norm();
    if(distance > f1.radius+f2.radius)
        return false;


    // TEST 2: Uno dei poligoni sta interamente a destra (o a sinistra)
    //          del piano contente l'altro
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
        double k = sign(f2.planeC.dot(f1.vertices[i]) + f2.planeD);
        if(! (k == s))
            flag = true;
    }
    if(!flag) return false;


    Eigen::Matrix2d MAT;
    Eigen::Vector2d TN;

    MAT << f1.planeC[0], f1.planeC[1],
         f2.planeC[0], f2.planeC[1];
    TN << -f1.planeD, -f2.planeD;

    Eigen::Vector2d x = MAT.colPivHouseholderQr().solve(TN);
    Eigen::Vector3d p1{{x(0),x(1),0}};
    // Eigen::Vector3d P{{0.017584598939879986, 0, 0}};
    Eigen::Vector3d v1 = f1.planeC.cross(f2.planeC);


    // vector<Eigen::Vector3d> v;
    v.reserve(4);

    for(int i = 0; i < f1.numVert; i++){
        Eigen::Vector3d B = f1.vertices[i], p2 = f1.vertices[(i+1)%f1.numVert];
        Eigen::Vector3d v2 = B-p2;

        Eigen::Vector3d w0 = p1-p2;
        double a = v1.dot(v1);
        double b = v1.dot(v2);
        double c = v2.dot(v2);
        double d = v1.dot(w0);
        double e = v2.dot(w0);

        double denom = a*c - b*b;
        if (denom == 0) {
            // std::cout << "Le rette sono parallele o coincidenti.\n";
            continue;
        }
        double sc = (b*e - c*d) / denom;
        double tc = (a*e - b*d) / denom;

        Eigen::Vector3d intersection = p2 + tc*v2;
        if(tc > -tol && tc < 1+tol)
            v.push_back(intersection);
    }


    for(int i = 0; i < f2.numVert; i++){
        Eigen::Vector3d B = f2.vertices[i], p2 = f2.vertices[(i+1)%f2.numVert];
        Eigen::Vector3d v2 = B-p2;

        Eigen::Vector3d w0 = p1-p2;
        double a = v1.dot(v1);
        double b = v1.dot(v2);
        double c = v2.dot(v2);
        double d = v1.dot(w0);
        double e = v2.dot(w0);

        double denom = a*c - b*b;
        if (denom == 0) {
            // std::cout << "Le rette sono parallele o coincidenti.\n";
            continue;
        }
        double sc = (b*e - c*d) / denom;
        double tc = (a*e - b*d) / denom;

        Eigen::Vector3d intersection = p2 + tc*v2;
        if(tc > -tol && tc < 1+tol)
            v.push_back(intersection);

    }

    if(v.size()>0)
        return true;




    return false;
}

int DFNLibrary::sign(double d){
    if (d>tol) return 1;
    if (d<-tol) return -1;

    cout << "Careful" << endl;

    return 0;
}




