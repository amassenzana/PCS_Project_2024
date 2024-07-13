#pragma once

#include <gtest/gtest.h>
#include "Eigen/Eigen"
#include "DFNLibrary.hpp"
#include <iostream>

#define testtol 1e-6

using namespace std;
using namespace Eigen;


namespace DFNLibrary{

TEST(DFN_UTILITIES, TestComputeCenter){
    Frattura f;
    f.numVert = 4;
    f.vertices = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}};
    Vector3d center = f.computeCenter(), C{0,0,0};
    EXPECT_LT(fabs(center[0]-C[0]), testtol);
    EXPECT_LT(fabs(center[1]-C[1]), testtol);
    EXPECT_LT(fabs(center[2]-C[2]), testtol);


    Frattura f2{{1},{4},{{1,1,0}, {1,-1,0}, {-1,1,0}, {-1,-1,0}},{}};
    EXPECT_LT(fabs(C[0]-f2.computeCenter()[0]), testtol);
    EXPECT_LT(fabs(C[1]-f2.computeCenter()[1]), testtol);
    EXPECT_LT(fabs(C[2]-f2.computeCenter()[2]), testtol);
}

TEST(DFN_UTILITIES, TestComputeRadius){
    Frattura f;
    f.numVert = 4;
    f.vertices = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}};
    double r = f.computeRadius();
    EXPECT_LT(fabs(r - 1), testtol); // r-1 = 0 --> r-1 < tol LT = LESS THAN

    Frattura f2;
    f2.numVert = 4;
    f2.vertices = {{1,1,0}, {1,-1,0}, {-1,1,0}, {-1,-1,0}};
    EXPECT_LT(fabs(f2.computeRadius() - sqrt(2)), testtol);
}

TEST(DFN_UTILITIES, TestComputePlane){
    Vector3d planeCoeff{0,0,0};
    Frattura f{{1},{4},{{1,1,0}, {1,-1,0}, {-1,1,0}, {-1,-1,0}},{}}; // Mi aspetto vettore normale = (0,0,1)
    double planeConstant = f.computePlane(planeCoeff);
    EXPECT_LT(fabs(planeConstant), testtol);
    EXPECT_LT(fabs(planeCoeff[0]), testtol);
    EXPECT_LT(fabs(planeCoeff[1]), testtol);
    EXPECT_NE(fabs(planeCoeff[2]), 0);          // NOT EQUALITY
}

TEST(DFN_UTILITIES, TestTraceLength){
    Traccia T;
    T.id = 0; T.origin = {0,0,0}; T.end={1,0,0};
    EXPECT_LT(fabs(T.length() - 1), testtol);

    Traccia T2;
    T2.origin = {-1, -1, -1}; T2.end = {1, 1, 1};
    EXPECT_LT(fabs(T2.length() - 2*sqrt(3)), testtol);
}

TEST(DFN_UTILITIES, TestSign){
    double d = 3.16;
    EXPECT_EQ(sign(d), 1);
    EXPECT_EQ(sign(-d), -1);
    EXPECT_EQ(sign(0.00000001), 0);
}

TEST(DFN_UTILITIES, TestCompareTrace){
    Traccia T1{{0}, {0,0,0}, {2,2,2}, {}, {}};
    Traccia T2{{1}, {0,0,0}, {1,1,1}, {}, {}};
    pair<Traccia*, bool> PT1{&T1, false};   // T1 Come NON passante
    pair<Traccia*, bool> PT2{&T2, true};    // T2 passante

    EXPECT_EQ(compareTrace(PT1,PT2), 0);

    PT1 = {&T1, false};
    PT2 = {&T2, false};
    EXPECT_EQ(compareTrace(PT1,PT2), 1);

    PT1 = {&T1, true};
    PT2 = {&T2, true};
    EXPECT_EQ(compareTrace(PT1,PT2), 1);
}

// TEST(DFN_UTILITIES, TestAngle){
//     Vector3d v1 = {1,0,0},  v2 = {0,1,0};
//     EXPECT_LT(fabs( angle(v1,v2) - M_PI/2 ), testtol);

//     Vector3d C = {0.5, 0.5, 0};
//     EXPECT_LT(fabs( angleReference(v1,v2,C) - M_PI ), testtol);
// }

// TEST(DFN_UTILITIES, TestPointSort){
//     Vector3d v1 = {1,0,0}, v2 = {1.2,0,0};
//     EXPECT_EQ(pointSort(v1,v2), 1);
// }

TEST(DFN_UTILITIES, TestLiesOnSegment){
    Eigen::Vector3d P={0,0,0};
    Eigen::Vector3d A={1,1,1};
    Eigen::Vector3d B={3,3,3};

    EXPECT_EQ( lies(A, B, P), 0 );

    A = {-1,-1,-1};
    EXPECT_EQ( lies(B, A, P), 1 );

    P = {-1.1,-1.1,-1.1};
    EXPECT_EQ( lies(B, A, P), 0 );
}


TEST(DFN_TESTS, TestImport){
    DFN dfn3;
    EXPECT_EQ(importDFN("FR3_data.txt", dfn3), 1);
    EXPECT_EQ(dfn3.numberFractures, 3);

    DFN dfn10;
    EXPECT_EQ(importDFN("FR10_data.txt", dfn10), 1);
    EXPECT_EQ(dfn10.numberFractures, 10);

    DFN dfn50;
    EXPECT_EQ(importDFN("FR50_data.txt", dfn50), 1);
    EXPECT_EQ(dfn50.numberFractures, 50);

    DFN dfn82;
    EXPECT_EQ(importDFN("FR82_data.txt", dfn82), 1);
    EXPECT_EQ(dfn82.numberFractures, 82);

    DFN dfn200;
    EXPECT_EQ(importDFN("FR200_data.txt", dfn200), 1);
    EXPECT_EQ(dfn200.numberFractures, 200);

    DFN dfn362;
    EXPECT_EQ(importDFN("FR362_data.txt", dfn362), 1);
    EXPECT_EQ(dfn362.numberFractures, 362);
}

TEST(DFN_TESTS, TestLiesOn){
    DFN dfn3;
    EXPECT_EQ(importDFN("FR3_data.txt", dfn3), 1);
    dfn3.computeDFN();
    dfn3.output();


    Traccia&  T1 = dfn3.tracce.front();
    Frattura& F1 = dfn3.fractures[T1.idF1];
    Frattura& F2 = dfn3.fractures[T1.idF2];

    for(int i = 0; i < F1.numVert; i++){
        if( lies( F1.vertices[i], F1.vertices[(i+1)%F1.numVert], T1.origin ) ){
            EXPECT_EQ(i,0);
        }

        if( lies( F1.vertices[i], F1.vertices[(i+1)%F1.numVert], T1.end ) ){
            EXPECT_EQ(i,2);
        }
    }

    for(int i = 0; i < F2.numVert; i++){
        if( lies( F2.vertices[i], F2.vertices[(i+1)%F2.numVert], T1.origin ) ) {
            EXPECT_EQ(i,0);
        }

        if( lies( F2.vertices[i], F2.vertices[(i+1)%F2.numVert], T1.end ) ) {EXPECT_EQ(i,2);}
    }

    Traccia&  T2 = dfn3.tracce.back();
    F1 = dfn3.fractures[T2.idF1];
    F2 = dfn3.fractures[T2.idF2];

    for(int i = 0; i < F1.numVert; i++){
        if( lies( F1.vertices[i], F1.vertices[(i+1)%F1.numVert], T2.origin ) ) {EXPECT_EQ(i,3);}
        if( lies( F1.vertices[i], F1.vertices[(i+1)%F1.numVert], T2.end ) ) {EXPECT_EQ(i,3);}
    }

    for(int i = 0; i < F2.numVert; i++){
        if( lies( F2.vertices[i], F2.vertices[(i+1)%F2.numVert], T2.origin ) ) {EXPECT_EQ(i,1);}
        if( lies( F2.vertices[i], F2.vertices[(i+1)%F2.numVert], T2.end ) ) {EXPECT_EQ(i,1);}
    }
}

TEST(DFN_TESTS, TestDFN3){
    DFN dfn3;
    EXPECT_EQ(importDFN("FR3_data.txt", dfn3), 1);
    dfn3.computeDFN();
    dfn3.output();

    EXPECT_EQ(dfn3.tracce.size(), 2);

    Vector3d expectOrigin = {0.8, 0, 0}, expectEnd = {0.8, 1, 0};
    Vector3d x = dfn3.tracce.front().origin - expectOrigin;
    EXPECT_LT(fabs(x[0]), testtol);
    EXPECT_LT(fabs(x[1]), testtol);
    EXPECT_LT(fabs(x[2]), testtol);

    x = dfn3.tracce.front().end - expectEnd;
    EXPECT_LT(fabs(x[0]), testtol);
    EXPECT_LT(fabs(x[1]), testtol);
    EXPECT_LT(fabs(x[2]), testtol);

    expectOrigin = {0, 0.5, 0};
    expectEnd = {0.3161837, 0.5, 0};
    x = dfn3.tracce.back().origin - expectOrigin;
    EXPECT_LT(fabs(x[0]), testtol);
    EXPECT_LT(fabs(x[1]), testtol);
    EXPECT_LT(fabs(x[2]), testtol);

    x = dfn3.tracce.back().end - expectEnd;
    EXPECT_LT(fabs(x[0]), testtol);
    EXPECT_LT(fabs(x[1]), testtol);
    EXPECT_LT(fabs(x[2]), testtol);
}


TEST(DFN_TESTS, TestDFN10){
    DFN dfn;
    EXPECT_EQ(importDFN("FR10_data.txt", dfn), 1);
    dfn.computeDFN();

    EXPECT_EQ(dfn.tracce.size(), 25);
}

TEST(DFN_TESTS, TestDFN82){
    DFN dfn;
    EXPECT_EQ(importDFN("FR82_data.txt", dfn), 1);
    dfn.computeDFN();

    EXPECT_EQ(dfn.tracce.size(), 1);
}
/*
TEST(DFN_TESTS, TestDFN362){
    DFN dfn;
    EXPECT_EQ(importDFN("FR362_data.txt", dfn), 1);
    dfn.computeDFN();

    EXPECT_EQ(dfn.tracce.size(), 1);
}
*/


} //namespace DFNLibrary

