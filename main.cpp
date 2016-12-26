//
//  main.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 2/13/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "PolyDecomp.hpp"
#include <fstream>


typedef pair<unsigned int, IntAndSet*> summands;

using namespace std;

vector<Point> randomGenerator();
vector<Point> circularGenerator(double, int);

int main(int argc, const char * argv[]) {
    MPI_Init(NULL, NULL);
    
    Polygon ConvexHull;
    
    fstream file;
    file.open(argv[1]);
    if (!file.is_open())
        return 9;
    
    coord_t x, y;
    while (file >> x >> y)
        ConvexHull.push_back({x, y}); // read input
    
    file.close();
    
    int m = (int)ConvexHull.size() - 1, maxN = 0, minN = INT32_MAX;
    
    double start, end;
    start = MPI_Wtime();
    
    init(ConvexHull);
    
    end = MPI_Wtime();
    
    cout << "Preprocessing time: " << end - start << endl;
    
    bool b = MPIPolyDecomp();

    MPI_Barrier(MPI_COMM_WORLD);
    
//    start = MPI_Wtime();
    
//    Alg2TypeX ret = MPIPolyDecompNumX();
    
//    end = MPI_Wtime();
    
//    cout << "MPIPolyDecompNum time: " << end - start << " seconds, " << ret.u << " summands" << endl;

//	if(b) MPI_Abort(MPI_COMM_WORLD, 0);


	MPI_Finalize();
    
    return 0;
}
vector<Point> circularGenerator(double c, int r) {
    vector<Point> v;
    coord_t x, y;
    for (int i = 0; i < c; i++) {
        x = ceil(r + r*cos(2*i*M_PI/c)), y = ceil(r + r*sin(2*i*M_PI/c));
        v.push_back({x,y});
    }
    return v;
}

