//
//  PolyDecomp.hpp
//  CMPS 299
//
//  Created by Ali Jbaily on 3/8/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#ifndef PolyDecomp_hpp
#define PolyDecomp_hpp

#include <iostream>
#include <set>
#include <map>
#include "convex_hull.hpp"
#include <math.h>
#include <algorithm>
#include <vector>
#include <string.h>
#include <climits>
#include "mpi.h"
#include "ArrayX.cpp"

typedef vector<Point> Polygon;

typedef struct{
    int k, i, u, x, y;
} buffer_t;

typedef struct {
    unsigned int u;
    set<pair<int, int>> S;
} IntAndSet;

typedef struct {
    unsigned int u;
    ArrayX<std::pair<int, int>> a;
} IntAndX;

typedef struct {
    map<Point, IntAndSet> *A;
    unsigned int u;
} Alg2Type;

typedef struct {
    IntAndX **A;
    unsigned int u;
} Alg2TypeX;

struct edge {
    int n;
    Point e;
};

int gcd (int, int);
bool PolyDecomp();
bool PolyDecompArray();
bool MPIPolyDecomp();
bool PolyDecompDFS(int, int, Point*);
bool MultiThreadedPolyDecompDFS(int threads);
Alg2TypeX MPIPolyDecompNum();
Alg2TypeX MPIPolyDecompNumX();
bool belongToIP(Point*);
void printer();
void init(Polygon);
void clean();
void printPoint(Point *);

#endif /* PolyDecomp_hpp */
