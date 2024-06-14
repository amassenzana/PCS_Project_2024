#include "DFNLibrary.hpp"
#include <fstream>
#include <list>
#include <iomanip>
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
                if ( fractures[i].IntersectionPoints.size() != 0)
                {
                    for (int k = 0; k < fractures[i].IntersectionPoints.size(); k++)
                    {
                        cout << setprecision(16) << scientific << fractures[i].IntersectionPoints[k] << endl;
                        cout << "" << endl;
                    }
                }
                else
                {
                    cout << "non va bene" << endl;
                }

            }
        }
    }


}

// 1)   Funzioni di proiezione
Eigen::Vector2d DFNLibrary::ProjXY(Eigen::Vector3d tmpvec_xy)
{
    Eigen::Vector2d projvec_xy (tmpvec_xy[0], tmpvec_xy[1]);
    return projvec_xy;
}

Eigen::Vector2d DFNLibrary::ProjXZ(Eigen::Vector3d tmpvec_xz)
{
    Eigen::Vector2d projvec_xz (tmpvec_xz[0], tmpvec_xz[2]);
    return projvec_xz;
}

Eigen::Vector2d DFNLibrary::ProjYZ(Eigen::Vector3d tmpvec_yz)
{
    Eigen::Vector2d projvec_yz (tmpvec_yz[1], tmpvec_yz[2]);
    return projvec_yz;
}


// 2)   Ricerca del punto di intersezione nelle proiezioni

std::pair<Eigen::Vector2d, bool> DFNLibrary::ProjIntersection(Eigen::Vector2d a1, Eigen::Vector2d a2, Eigen::Vector2d b1, Eigen::Vector2d b2)
{
    using Line2 = Eigen::Hyperplane<double,2>;
    using Vec2 = Eigen::Vector2d;
    bool tmpcheck = false;

    Vec2 a(a1[0], a1[1]);
    Vec2 b(a2[0], a2[1]);

    Vec2 c(b1[0], b1[1]);
    Vec2 d(b2[0], b2[1]);

    Line2 ab = Line2::Through(a,b);
    Line2 cd = Line2::Through(c,d);

    Vec2 qp = ab.intersection(cd);

    // Escludo le linee coincidenti
    if (!ab.isApprox(cd))
    {
        // Escludo le linee parallele
        if (qp[0] != std::numeric_limits<double>::infinity() || qp[1] != std::numeric_limits<double>::infinity())
        {
            tmpcheck = true;
        }
    }

    return std::make_pair(qp,tmpcheck);
};

// qp è il punto di partenza della retta d'intersezione tra i due piani generati dalle fratture in uso.

// 3)   Funzione che calcola l'intersezione tra due linee in uno spazio 3D.
// Adattare l'algoritmo per ogni proiezione possibile dai 3 piani.
// qp è il punto ottenuto dalla proiezione, dqp e il vettore direzione della proiezione del punto in una linea infinita, mentre p1 e p2 sono due punti della retta con la
//quale si effettua il confronto.
Eigen::Vector3d DFNLibrary::IntersectVerif(Eigen::Vector3d qp_l, Eigen::Vector3d dqp, Eigen::Vector3d p1, Eigen::Vector3d p2) {
    Eigen::Vector3d d1 = p2-p1;
    Eigen::Vector3d n = dqp.cross(d1);
    Eigen::Vector3d sr = qp_l - p1;

    // Check if the lines are parallel
    if (n.isZero(1e-6)) {
        std::cout << "The lines are parallel, returning NaN vector." << std::endl;
        return Eigen::Vector3d(std::nan(""), std::nan(""), std::nan(""));
    }

    double t = (sr.cross(d1)).dot(n) / n.dot(n);

    // Check if t is infinite
    if (std::isinf(t)) {
        std::cout << "The lines are infinitely long, returning NaN vector." << std::endl;
        return Eigen::Vector3d(std::nan(""), std::nan(""), std::nan(""));
    }

    return qp_l + dqp * t;
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

    // INTERSEZIONE
    Eigen::Vector3d tmpvec1_1 , tmpvec1_2;
    Eigen::Vector3d tmpvec2_1 , tmpvec2_2;
    Eigen::Vector2d prj1_1_xy , prj1_1_xz, prj1_1_yz, prj1_2_xy , prj1_2_xz, prj1_2_yz;
    Eigen::Vector2d prj2_1_xy , prj2_1_xz, prj2_1_yz, prj2_2_xy , prj2_2_xz, prj2_2_yz;
    int tmpcount = 0;
    Eigen::Vector3d PlaneDir = f1.planeC;
    for (int i = 0; i <= f1.numVert; i++)
    {
        Eigen::Vector3d tmpvec1_1 = f1.vertices[i%f1.numVert];
        Eigen::Vector3d tmpvec1_2 = f1.vertices[(i+1)%f1.numVert];
        Eigen::Vector2d prj1_1_xy = ProjXY(tmpvec1_1);
        Eigen::Vector2d prj1_1_xz = ProjXZ(tmpvec1_1);
        Eigen::Vector2d prj1_1_yz = ProjYZ(tmpvec1_1);
        Eigen::Vector2d prj1_2_xy = ProjXY(tmpvec1_2);
        Eigen::Vector2d prj1_2_xz = ProjXZ(tmpvec1_2);
        Eigen::Vector2d prj1_2_yz = ProjYZ(tmpvec1_2);
        for (int j = 0; j < f2.numVert; j++)
        {
            tmpcount = 0;
            tmpvec2_1 = f2.vertices[j%f2.numVert];
            tmpvec2_2 = f2.vertices[(j+1)%f2.numVert];
            prj2_1_xy = ProjXY(tmpvec2_1);
            prj2_1_xz = ProjXZ(tmpvec2_1);
            prj2_1_yz = ProjYZ(tmpvec2_1);
            prj2_2_xy = ProjXY(tmpvec2_2);
            prj2_2_xz = ProjXZ(tmpvec2_2);
            prj2_2_yz = ProjYZ(tmpvec2_2);

            auto result = ProjIntersection(prj1_1_xy, prj1_2_xy, prj2_1_xy, prj2_2_xy);
            bool tmpflag = result.second;
            if (tmpflag == true)
            {
                Eigen::Vector2d qp = result.first;
                Eigen::Vector3d qp_xy (qp[0], qp[1], 0);
                // Eigen::Vector3d dqp (0, 0, 1); -> FALSO
                f1.IntersectionPoints.push_back(IntersectVerif(qp_xy, PlaneDir, tmpvec2_1, tmpvec2_2));
                tmpcount++;
                f1.GoodValsDouble.push_back(Eigen::Vector2d (i,j));
                std::vector<Eigen::Vector3d> tmpvec = {qp_xy, PlaneDir, tmpvec2_1, tmpvec2_2};
                f1.GoodValsVec.push_back(tmpvec);
            }
            else
            {
                auto result = ProjIntersection(prj1_1_xz, prj1_2_xz, prj2_1_xz, prj2_2_xz);
                bool tmpflag = result.second;
                if (tmpflag == true)
                {
                    Eigen::Vector2d qp = result.first;
                    Eigen::Vector3d qp_xz (qp[0], 0, qp[1]);
                    // Eigen::Vector3d dqp (0, 1, 0); -> FALSO
                    f1.IntersectionPoints.push_back(IntersectVerif(qp_xz, PlaneDir, tmpvec2_1, tmpvec2_2));
                    tmpcount++;
                    f1.GoodValsDouble.push_back(Eigen::Vector2d (i,j));
                    std::vector<Eigen::Vector3d> tmpvec = {qp_xz, PlaneDir, tmpvec2_1, tmpvec2_2};
                    f1.GoodValsVec.push_back(tmpvec);
                }
                else
                {
                    auto result = ProjIntersection(prj1_1_yz, prj1_2_yz, prj2_1_yz, prj2_2_yz);
                    bool tmpflag = result.second;
                    if (tmpflag == true)
                    {
                        Eigen::Vector2d qp = result.first;
                        Eigen::Vector3d qp_yz (0, qp[0], qp[1]);
                        // Eigen::Vector3d dqp (1, 0, 0); -> FALSO
                        f1.IntersectionPoints.push_back(IntersectVerif(qp_yz, PlaneDir, tmpvec2_1, tmpvec2_2));
                        tmpcount++;
                        f1.GoodValsDouble.push_back(Eigen::Vector2d (i,j));
                        std::vector<Eigen::Vector3d> tmpvec = {qp_yz, PlaneDir, tmpvec2_1, tmpvec2_2};
                        f1.GoodValsVec.push_back(tmpvec);
                    }
                }
            }
        }
    }

    // for (int m = 0; m < f1.GoodValsDouble.size(); m++)
    // {
    //     for (int n = 0; n < f2.GoodValsDouble.size(); n++)
    //     {
    //         std::vector j1 = {f1.GoodValsDouble[m][0], f1.GoodValsDouble[m][0]+1};
    //         std::vector j2 = {f1.GoodValsDouble[m][1], f1.GoodValsDouble[m][0]+1};
    //         tmpvec2_1 = f2.vertices[int (j1[0])%f2.numVert];
    //         tmpvec2_2 = f2.vertices[int (j1[1])%f2.numVert];
    //         Eigen::Vector3d tmpvecfinal = IntersectVerif(f1.GoodValsVec[m][0], f1.GoodValsVec[m][1], f1.GoodValsVec[m][2], f1.GoodValsVec[m][3]);
    //         // Eigen::Vector3d tmpvecfinal = IntersectVerif(f1.GoodValsVec[m][0], f1.GoodValsVec[m][1], tmpvec2_1, tmpvec2_2);
    //         f1.IntersectionPoints.push_back(tmpvecfinal);
    //     }
    // }

    if (tmpcount > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
    // tmpcount = 0;
}

// Eigen::Vector3d IntersecPointCalc


int DFNLibrary::sign(double d){
    if (d>tol) return 1;
    if (d<-tol) return -1;

    cout << "Careful" << endl;

    return 0;
}




