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
            if(checkIntersection(fractures[i], fractures[j], v)){
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
//    p1 = MAT.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(TN);
    p1 = MAT.colPivHouseholderQr().solve(TN);



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

bool DFNLibrary::lies(Eigen::Vector3d &A, Eigen::Vector3d &B, Eigen::Vector3d &P){
    Eigen::Vector3d AB = B-A;
    Eigen::Vector3d AP = P-A;

    Eigen::Vector3d cross = AP.cross(AB);

    if(fabs(cross.norm()) > tol){
        return false;
    }

    double dot1 = AB.dot(AP);
    double dot2 = AB.dot(AB);

    if(dot1 < 0 || dot1>dot2){
        return false;
    }

    return true;
}

void DFNLibrary::cutDFN(DFN &dfn, std::vector<PolygonalLibrary::PolygonalMesh> &globalMesh){
    for(Frattura& F : dfn.fractures){
        PolygonalLibrary::PolygonalMesh mesh;
        // INIZIALIZZAZIONI

        mesh.Cell0DId.reserve(F.numVert);
        mesh.Cell0DCoordinates.reserve(F.numVert);
        mesh.NumberCell0D = F.numVert;

        mesh.Cell1DId.reserve(F.numVert);
        mesh.Cell1DVertices.reserve(F.numVert);
        mesh.NumberCell1D = F.numVert;

        mesh.Cell2DEdges.resize( F.tracce.size() +1 );
        mesh.Cell2DVertices.resize( F.tracce.size() +1 );
        mesh.Cell2DId.push_back(mesh.NumberCell2D++);


        // int k = 0;

        for(int i = 0; i < F.numVert; i++){
            mesh.Cell0DId.push_back(i);
            mesh.Cell0DCoordinates.push_back(F.vertices[i]);

            mesh.Cell1DId.push_back(i);
            mesh.Cell1DVertices.push_back({i, (i+1)%F.numVert});

            mesh.Cell2DVertices[0].push_back(i);
            mesh.Cell2DEdges[0].push_back(i);
        }

        list<pair<Eigen::Vector3d, Eigen::Vector3d>>  listaTagli;
        for( pair<Traccia*, bool>& T : F.tracce ){
            listaTagli.push_back(make_pair( (T.first)->origin , (T.first)->end ));
        }

        // TAGLIO RICORSIVO
        cutPolygon(listaTagli, mesh, 0);

        globalMesh.push_back(mesh);
    }
}

void DFNLibrary::cutPolygon(std::list<std::pair<Eigen::Vector3d, Eigen::Vector3d> > &listaTagli, PolygonalLibrary::PolygonalMesh &mesh, unsigned int idP){

    if(listaTagli.empty()) return;

    pair<Eigen::Vector3d, Eigen::Vector3d> T = listaTagli.front();
    listaTagli.pop_front();

    unsigned int k = 0, i = 0;
    // bool cutIt[2] = {true,true};
    Eigen::Vector3d inters[2] = {};
    unsigned int lati[2], inizioFinePD[2] = {0,0}, ifPD[2] = {0,0};
    for( const unsigned int& lato : mesh.Cell2DEdges[idP] ){
        Eigen::Vector2i punti = mesh.Cell1DVertices[lato];
        Eigen::Vector3d A = mesh.Cell0DCoordinates[punti[0]];
        Eigen::Vector3d B = mesh.Cell0DCoordinates[punti[1]];
        if( prolungate(A,B,T.first,T.second,inters[k]) ){
            lati[k] = lato;
            if( (A-inters[k]).norm() < tol){
                inters[k] = A;
                // cutIt[k] = false;
            }
            else if( (B-inters[k]).norm() < tol){
                inters[k] = A;
                // cutIt[k] = false;
            }
            k++;
            if(k==2) break;
        }
        i++;
    }

    // Conosco i punti della prolungazione e i lati coinvolti; Procedo a tagliare
    vector<unsigned int> PDEdges, PSEdges;
    vector<unsigned int> PDVert, PSVert;
    mesh.Cell0DCoordinates.push_back(inters[0]);
    mesh.Cell0DId.push_back(mesh.NumberCell0D++);
    mesh.Cell0DCoordinates.push_back(inters[1]);
    mesh.Cell0DId.push_back(mesh.NumberCell0D++);

    mesh.Cell1DVertices.push_back(Eigen::Vector2i{mesh.NumberCell0D -2, mesh.NumberCell0D -1});
    mesh.Cell1DId.push_back(mesh.NumberCell1D++);

    mesh.Cell1DVertices.push_back(Eigen::Vector2i{mesh.Cell1DVertices[lati[0]][0],mesh.NumberCell0D -2});
    mesh.Cell1DVertices.push_back(Eigen::Vector2i{mesh.Cell1DVertices[lati[0]][1],mesh.NumberCell0D -2});
    mesh.Cell1DVertices.push_back(Eigen::Vector2i{mesh.Cell1DVertices[lati[1]][0],mesh.NumberCell0D -1});
    mesh.Cell1DVertices.push_back(Eigen::Vector2i{mesh.Cell1DVertices[lati[1]][1],mesh.NumberCell0D -1});

    mesh.Cell1DId.push_back(mesh.NumberCell1D++);
    mesh.Cell1DId.push_back(mesh.NumberCell1D++);
    mesh.Cell1DId.push_back(mesh.NumberCell1D++);
    mesh.Cell1DId.push_back(mesh.NumberCell1D++);

    for( i = 0; i < mesh.Cell2DVertices[idP].size(); i++ ){
        if( mesh.Cell2DEdges[idP][i] == lati[0]  ){

            PSEdges.push_back(mesh.NumberCell1D-4);
            PSEdges.push_back(mesh.NumberCell1D-5);

            inizioFinePD[0] = 1;

        } else if( mesh.Cell2DEdges[idP][i] == lati[1] ){

            PDEdges.push_back(mesh.NumberCell1D - 2);
            PDEdges.push_back(mesh.NumberCell1D - 5);
            PDEdges.push_back(mesh.NumberCell1D - 3);

            PSEdges.push_back(mesh.NumberCell1D - 1);

            inizioFinePD[1] = 1;

        } else {
            if( inizioFinePD[0] == 0 ){
                PSEdges.push_back(mesh.Cell2DEdges[idP][i]);
            }
            else if(inizioFinePD[0] == 1 && inizioFinePD[1] == 0){
                PDEdges.push_back(mesh.Cell2DEdges[idP][i]);
            }
            else{
                PSEdges.push_back(mesh.Cell2DEdges[idP][i]);
            }
        }

        Eigen::Vector3d A = mesh.Cell0DCoordinates[mesh.Cell2DVertices[idP][i]];
        Eigen::Vector3d B = mesh.Cell0DCoordinates[mesh.Cell2DVertices[idP][(i+1)%mesh.Cell2DVertices[idP].size()]];

        if( lies(A, B, inters[0]) ){
            PSVert.push_back(mesh.Cell2DVertices[idP][i]);
            PDVert.push_back(mesh.NumberCell0D-2);
            PSVert.push_back(mesh.NumberCell0D-2);


            if(ifPD[0] == 1) ifPD[1] = 1;
            ifPD[0] = 1;
        }
        else if( lies(A,B, inters[1]) ){
            PDVert.push_back(mesh.Cell2DVertices[idP][i]);
            PDVert.push_back(mesh.NumberCell0D-1);
            PSVert.push_back(mesh.NumberCell0D-1);

            if(ifPD[0] == 1) ifPD[1] = 1;
            ifPD[0] = 1;
        }
        else if( ifPD[0] == 0 ){
            PSVert.push_back(mesh.Cell2DVertices[idP][i]);
        }
        else if( ifPD[0] == 1 && ifPD[1] == 0 ){
            PDVert.push_back(mesh.Cell2DVertices[idP][i]);
        }
        else {
            PSVert.push_back(mesh.Cell2DVertices[idP][i]);
        }


    }

    // CONTROLLARE TUTTI GLI ALTRI POLIGONI E RIMANEGGIARE PER MANTENERE COERENZA:
    // I LATI TAGLIATI, NEI POLIGONI VICINI VANNO RIAGGIORNATI PER AVERE DUE LATI PARALLELI
    for( i = 0; i < mesh.NumberCell2D; i++ ){
        if( mesh.Cell2DId[i] == idP ) continue;

        bool yes = false;
        unsigned int temp, values;
        for( k = 0; k < mesh.Cell2DEdges[i].size(); k++ ){
            if(yes){
                mesh.Cell2DEdges[i][k] = temp;
                temp = mesh.Cell2DEdges[i][k+1];
            }
            else if( mesh.Cell2DEdges[i][k] == lati[0] ){
                mesh.Cell2DEdges[i].push_back(0);
                yes = true;
                mesh.Cell2DEdges[i][k] = mesh.NumberCell1D - 3;
                temp = mesh.Cell2DEdges[i][k+1];
                mesh.Cell2DEdges[i][k+1] = mesh.NumberCell1D - 4;

                values = 0;

                k++;
            }
            else if( mesh.Cell2DEdges[i][k] == lati[1] ){
                mesh.Cell2DEdges[i].push_back(0);
                yes = true;
                mesh.Cell2DEdges[i][k] = mesh.NumberCell1D - 1;
                temp = mesh.Cell2DEdges[i][k+1];
                mesh.Cell2DEdges[i][k+1] = mesh.NumberCell1D - 2;

                values = 1;

                k++;
            }
        }
        if(yes){
            // mesh.Cell2DEdges[i].push_back(temp);
            unsigned int v1 = mesh.Cell1DVertices[lati[values]][0];
            unsigned int v2 = mesh.Cell1DVertices[lati[values]][1];

            for(auto it = mesh.Cell2DVertices[i].begin(); it != mesh.Cell2DVertices[i].end(); it++){
                auto it2 = it+1;
                if(it2 == mesh.Cell2DVertices[i].end()) it2 = mesh.Cell2DVertices[i].begin();
                if(*it == v1 && *it2 == v2){
                    mesh.Cell2DVertices[i].insert(++it, mesh.NumberCell0D-1);
                    break;
                } else if(*it == v2 && *it2 == v1){
                    mesh.Cell2DVertices[i].insert(++it, mesh.NumberCell0D-2);
                    break;
                }
            }

        }

    }

    // Ridistribuire le tracce su PDX e PSX
    list<pair<Eigen::Vector3d, Eigen::Vector3d>> listaTagliDx, listaTagliSx;
    for( unsigned int j = 0; j < listaTagli.size(); j++ ){
        pair<Eigen::Vector3d, Eigen::Vector3d> taglio = listaTagli.front();
        listaTagli.pop_front();

        if( staDentro(taglio, mesh, PDVert) ){
            listaTagliDx.push_back(taglio);
        } else {
            listaTagliSx.push_back(taglio);
        }
    }

    mesh.Cell2DEdges[idP] = PDEdges;
    mesh.Cell2DVertices[idP] = PDVert;
    mesh.Cell2DEdges[mesh.NumberCell2D] = PSEdges;
    mesh.Cell2DVertices[mesh.NumberCell2D] = PSVert;
    mesh.Cell2DId.push_back(mesh.NumberCell2D++);

    cutPolygon(listaTagliDx, mesh, idP);
    cutPolygon(listaTagliSx, mesh, mesh.NumberCell2D-1);
}

bool DFNLibrary::prolungate(Eigen::Vector3d &A, Eigen::Vector3d &B, Eigen::Vector3d &C, Eigen::Vector3d &D, Eigen::Vector3d& inters){
    if( lies(A,B,C) ){
        inters = C;
        return true;
    }
    if( lies(A,B,D) ){
        inters = D;
        return true;
    }

    Eigen::Vector3d AB = B-A;
    Eigen::Vector3d CD = D-C;
    Eigen::Vector3d AC = C-A;

    Eigen::Matrix2d matrix;
    matrix << AB.dot(AB), -CD.dot(AB),
              AB.dot(CD), -CD.dot(CD);

    if (std::abs(matrix.determinant()) < 1e-10) {
        // Parallele o coincidenti
        return false;
    }

    Eigen::Vector2d rhs;
    rhs << AB.dot(AC),
        CD.dot(AC);

    Eigen::Vector2d ts = matrix.inverse() * rhs;
    double t = ts(0);
    // double s = ts(1);

    if(t < 0-tol || t > 1+tol) {
        // L'intersezione esterna ad AB
        return false;
    }

    inters = A + t * AB;

    return true;
}

bool DFNLibrary::staDentro(std::pair<Eigen::Vector3d, Eigen::Vector3d> &taglio, PolygonalLibrary::PolygonalMesh &mesh, std::vector<unsigned int> &P){

    Eigen::Vector3d normal = (mesh.Cell0DCoordinates[P[1]] - mesh.Cell0DCoordinates[P[0]]).cross(mesh.Cell0DCoordinates[P[2]] - mesh.Cell0DCoordinates[P[0]]).normalized();

    for( unsigned int i = 0; i < P.size(); i++ ){
        const Eigen::Vector3d& A = mesh.Cell0DCoordinates[P[i]];
        const Eigen::Vector3d& B = mesh.Cell0DCoordinates[P[(i + 1) % P.size()]];
        Eigen::Vector3d edgeNormal = (B - A).cross(normal).normalized();

        double d = edgeNormal.dot(taglio.first - A);
        if(d < 0) return false;
    }

    return true;
}
