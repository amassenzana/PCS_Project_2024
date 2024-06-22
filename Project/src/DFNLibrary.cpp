#include "DFNLibrary.hpp"
#include <fstream>
#include <iomanip>
#include <list>
#include <algorithm>
#define tol 1e-6

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

    for(int k = 0; k < dfn.numberFractures; k++){
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
            for(int j = 0; j < f.numVert; j++){
                char tmp;                               // f.vertices[j](i) = stod(lines[i].substr(0,pos));
                converter >> f.vertices[j][i] >> tmp;   // lines[i] = lines[i].substr(pos+1, string::npos);
            }
            // f.vertices[f.numVert-1](i) = stod(lines[i]);
        }

        dfn.fractures.push_back(f);
    }


    return true;
}

Eigen::Vector3d DFNLibrary::Frattura::computeCenter(){


    Eigen::Vector3d c = {0, 0, 0};
    double sw = 0;

    for(int i = 1; i < numVert+1; i++){
        double w =  (vertices[i%numVert]-vertices[(i+1)%numVert]).norm() +
                    (vertices[i%numVert]-vertices[(i-1)%numVert]).norm();
        c += vertices[i%numVert]*w;
        sw += w;
    }

    c[0] = c[0]/sw;
    c[1] = c[1]/sw;
    c[2] = c[2]/sw;

    return c;
}

double DFNLibrary::Frattura::computeRadius(){
    double max = 0, t = 0;
    for(int i = 0; i < numVert; i++){
        t = (computeCenter()-vertices[i]).norm();
        if(t>max)
            max = t;
    }

    return max;
}


double DFNLibrary::Frattura::computePlane(Eigen::Vector3d& planeC){

    Eigen::Vector3d p1 = vertices[0];
    Eigen::Vector3d p2 = vertices[1];
    Eigen::Vector3d p3 = computeCenter();

    Eigen::Vector3d v1 = p2-p1;
    Eigen::Vector3d v2 = p3-p1;
    Eigen::Vector3d n = v1.cross(v2);

    planeC = {n.x(), n.y(), n.z()};
    return -n.dot(p1);
}

void DFNLibrary::DFN::computeDFN(){
    for(int i = 0; i < numberFractures; i++){
        for(int j = i+1; j < numberFractures; j++){
            vector<Eigen::Vector3d> v;
            if(i!=j && checkIntersection(fractures[i], fractures[j], v)){
                // cout << "Intersezione tra f[" << i << "] e f[" << j << "]:" << endl;
                elaborateV(fractures[i], fractures[j], v);
            }
        }
    }


}

void DFNLibrary::DFN::plotFracture(){
    ofstream file("MLPlot.txt");
    char xyz[3] = {'x', 'y', 'z'};
    char color[] = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};
    int k = 0;
    file << "figure\nhold on\n";
    for(const DFNLibrary::Frattura& f : fractures){
        for(int i = 0; i < 3; i++){
            file << xyz[i] << " = [";
            for(int j = 0; j < f.numVert; j++){
                file << f.vertices[j][i] << "; ";
            }
            file << "]';" << endl;
        }
        file << "fill3(x,y,z,'"<< color[(k++)%8] << "');\n\n" << endl;
    }

    k = 0;
    for(const DFNLibrary::Traccia& T : tracce){
        file << "plot3(";
        for(int i = 0; i < 3; i++){
            file << T.origin[i] << ", ";
        }
        file << "'o', 'LineWidth',3);";
        file << "\nplot3(";
        for(int i = 0; i < 3; i++){
            file << T.end[i] << ", ";
        }
        file << "'o', 'LineWidth',3);" << endl;
        if((T.origin-T.end).norm() < tol)
            k++;
        // cout << (T.origin-T.end).norm() << "  " << endl;
    }
    // cout << "Da escludere: " << k << endl;
}

void DFNLibrary::DFN::elaborateV(Frattura &f1, Frattura &f2, std::vector<Eigen::Vector3d> &v){
    int coppie = 0;
    for(int i = 0; i < 2; i++){
        for(int j = 2; j < 4; j++){
            if((v[i] - v[j]).norm() < tol)
                coppie++;
        }
    }

    if(coppie==2){
        // cout << "Traccia passante per entrambe" << endl;
        Traccia T{{int(tracce.size())}, {v[0]}, {v[1]}, f1.id, f2.id};
        for(size_t i = 0; i < v.size(); i++){
            if((T.origin - v[i]).norm() > tol){
                T.end = v[i];
                break;
            }
        }
        tracce.push_back(T);
        fractures[f1.id].tracce.push_back(make_pair(&(tracce.back()), true));
        fractures[f2.id].tracce.push_back(make_pair(&(tracce.back()), true));
    }

    else if(coppie==0){
        Traccia T{{int(tracce.size())}, {}, {}, {f1.id}, {f2.id}};
        Eigen::Vector3d temp, dir{v[1]-v[0]};
        bool passante[2] = {false, false};
        int k;
        vector<Eigen::Vector3d> backup = v;
        // bool swapped = false;
        for(k = 0; k < 3; k++){
            if(fabs(dir[k])>tol)
                break;
        }
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++)
                if( v[i][k] > v[i+1][k] ){
                    temp = v[i];
                    v[i] = v[i+1];
                    v[i+1] = temp;
                }

        T.origin = v[1];
        T.end = v[2];

        if( (T.origin == backup[0] || T.origin == backup[1]) &&
            (T.end == backup[0] || T.end == backup[1]) ){
            passante[0] = true;
            // cout << "Traccia passante per f["<< f1.id <<"] e non passante per f["<< f2.id << "]"<< endl;
        }

        else if( (T.origin == backup[2] || T.origin == backup[3]) &&
                 (T.end == backup[2] || T.end == backup[3]) ){
            passante[1] = true;
            // cout << "Traccia non passante per f["<< f1.id <<"] e passante per f["<< f2.id << "]"<< endl;
        }
        // else
            // cout << "Traccia non passante per entrambe" << endl;

        tracce.push_back(T);
        fractures[f1.id].tracce.push_back(make_pair(&(tracce.back()), passante[0]));
        fractures[f2.id].tracce.push_back(make_pair(&(tracce.back()), passante[1]));
    }
    else if(coppie==1){
        Traccia T{{int(tracce.size())}, {},{},{f1.id}, {f2.id}};
        bool passante[2] = {false, false};
        if((v[0]-v[1]).norm() < (v[2]-v[3]).norm()){
            passante[0] = true;
            T.origin = v[0];
            T.end = v[1];
            // cout << "Traccia passante per f["<< f1.id <<"] e non passante per f["<< f2.id << "]"<< endl;
        }
        else {
            passante[1] = true;
            T.origin = v[2];
            T.end = v[3];
            // cout << "Traccia non passante per f["<< f1.id <<"] e passante per f["<< f2.id << "]"<< endl;
        }
        tracce.push_back(T);
        fractures[f1.id].tracce.push_back(make_pair(&(tracce.back()), passante[0]));
        fractures[f2.id].tracce.push_back(make_pair(&(tracce.back()), passante[1]));
    }

}

void DFNLibrary::DFN::output(){
    ofstream file("Output.txt");
    file << "# Number of Traces\n" << tracce.size() << endl;

    file << scientific << setprecision(16);
    file << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    for(const Traccia& T : tracce){
        file << T.id << "; " << T.idF1 << "; " << T.idF2 << "; "
             << T.origin[0] << "; " << T.origin[1] << "; " << T.origin[2] << "; "
             << T.end[0]    << "; " << T.end[1]    << "; " << T.end[2]   << endl;
    }

    file.close();

    // ORDINAMENTO TRACCE

    file.open("OutputFractures.txt");
    file << scientific << setprecision(16);
    for(Frattura& f : fractures){
        file << "# FractureId; NumTraces" << endl;
        file << f.id << "; " << f.tracce.size() << endl;
        file << "# TraceId; Tips; Length" << endl;


        sort( (f.tracce).begin(), (f.tracce).end(),  compareTrace );


        for(const auto& coppia : f.tracce){
            file << coppia.first->id << "; " << coppia.second
                 << "; " << coppia.first->length() << endl;
        }
    }

}

bool DFNLibrary::checkIntersection(Frattura &f1, Frattura &f2, vector<Eigen::Vector3d>& v){

    // TEST 1: La distanza tra f1 e f2 è maggiore della somma dei raggi
    //  delle ipersfere che le contengono
    double distance = (f1.computeCenter()-f2.computeCenter()).norm();
    if(distance > f1.computeRadius()+f2.computeRadius())
        return false;


    // TEST 2: Uno dei poligoni sta interamente a destra (o a sinistra)
    //          del piano contente l'altro
    Eigen::Vector3d f1pc, f2pc;
    double f1pd = f1.computePlane(f1pc), f2pd = f2.computePlane(f2pc);

    int s = sign(f1pc.dot(f2.vertices[0]) + f1pd);
    bool flag = false;

    for(int i = 1; i < f2.numVert; i++){
        if(!(sign(f1pc.dot(f2.vertices[i]) + f1pd) == s))
            flag = true;
    }
    if(!flag) return false;

    s = sign(f2pc.dot(f1.vertices[0]) + f2pd);
    flag = false;
    for(int i = 1; i < f1.numVert; i++){
        double k = sign(f2pc.dot(f1.vertices[i]) + f2pd);
        if(! (k == s))
            flag = true;
    }
    if(!flag) return false;


    // TEST 3: Calcolo intersezione effettiva


    Eigen::Vector3d v1 = f1pc.cross(f2pc);

    Eigen::MatrixXd MAT = Eigen::MatrixXd::Identity(2,3);
    MAT.row(0) = f1pc;
    MAT.row(1) = f2pc;

    Eigen::VectorXd TN{{-f1pd}, {-f2pd}};
    Eigen::VectorXd p1(3,1);
    p1 = MAT.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(TN);



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
        // double sc = (b*e - c*d) / denom;
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
        // double sc = (b*e - c*d) / denom;
        double tc = (a*e - b*d) / denom;

        Eigen::Vector3d intersection = p2 + tc*v2;
        if(tc > -tol && tc < 1+tol)
            v.push_back(intersection);

    }


    if(v.size()==4){
        vector<Eigen::Vector3d> backup = v;
        int k;
        Eigen::Vector3d temp;
        for(k = 0; k < 3; k++){
            if(fabs(v1[k])>tol)
                break;
        }
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++)
                if( backup[i][k] > backup[i+1][k] ){
                    temp = backup[i];
                    backup[i] = backup[i+1];
                    backup[i+1] = temp;
                }

        if( v[0][k] > v[1][k]  ){
            temp = v[0];
            v[0] = v[1];
            v[1] = temp;
        }

        if( v[2][k] > v[3][k]  ){
            temp = v[2];
            v[2] = v[3];
            v[3] = temp;
        }

        if( (backup[0]==v[0] && backup[1]==v[1]) ||
            (backup[0]==v[2] && backup[1]==v[3]) ) return false;

        else if((backup[0]==v[1] && backup[1]==v[0]) ||
                (backup[0]==v[3] && backup[1]==v[2]) ) return false;

        return true;
    }





    return false;
}

int DFNLibrary::sign(double d){
    if (d>tol) return 1;
    if (d<-tol) return -1;

    // cout << "Careful - Book case" << endl;

    return 0;
}

bool DFNLibrary::pointSort(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2){
    return p1.norm() < p2.norm();
}

double DFNLibrary::Traccia::length(){
    return (origin-end).norm();
}

bool DFNLibrary::compareTrace(std::pair<Traccia*, bool>& T1, std::pair<Traccia*, bool>& T2){
    if(T1.second != T2.second){
        if(T1.second) return true;
        else return false;
    }

    return (T1.first->length() > T2.first->length());
}


double DFNLibrary::angle(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2){
    double d = v1.dot(v2);
    double lv1 = v1.norm();
    double lv2 = v2.norm();
    return std::acos(d / (lv1 * lv2));
}

double DFNLibrary::angleReference(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& center){
    Eigen::Vector3d v1c = center-v1;
    Eigen::Vector3d v2c = center-v2;
    return angle(v1c, v2c);
}
