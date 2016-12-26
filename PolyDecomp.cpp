//
//  PolyDecomp.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 3/8/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "PolyDecomp.hpp"
#include <stdlib.h>


#define CHUNKSIZE 1000

Point** IP, *v0;
vector<Point> latticePoints;
Polygon ConvexHull;
int **flags;

int ymin = INT32_MAX, ymax = 0, xmin = INT32_MAX, xmax = 0, m, world_rank, world_size;
edge* edge_sequence;
bool B;
IntAndX** newA;

void printPoint(Point *p){
    cout << p->x << ", " << p->y << endl;
}

void info(){
    int minN = INT32_MAX, maxN = 0;
    double avgN = 0, count, ratio;
    for (int i = 0; i < m; i++) {
        int n = edge_sequence[i].n;
        if (n > maxN) maxN = n;
        else if (abs(n) < minN) minN = abs(n);
        avgN += n;
    }
    
    for (coord_t y = 0; y <= ymax; y++)
        for (coord_t x = 0; x <= xmax; x++) {
            Point p = {x,y};
            if (belongToIP(&p))
                count++;
        }
    
    ratio = count/((xmax-xmin+1)*(ymax-ymin+1));
    avgN /= m;
    cout << "m = " << m << ", max n = " << maxN << ", min n = " << minN << ", average n = " << avgN << endl;
    cout << "xmin = " << xmin << ", xmax = " << xmax << ", ymin = " << ymin << ", ymax = " << ymax << endl;
    cout << "dx = " << xmax - xmin << ", dy = " << ymax - ymin << ", of size = " << (ymax-ymin+1)*(xmax-xmin+1) << endl;
    cout << "IP size = " << count << ", ratio = " << ratio << endl;
}

void allocateMemory(int* size, Point* array) {
    *size <<= 1;
    Point *oldptr = array;
    array = new Point[*size];
    memcpy(array, oldptr, sizeof(Point)*(*size >> 1));
    delete [] oldptr;
    //    cout << "new allocation" << endl;
}

bool PolyDecomp() {
    
    Point p;
    set<Point> oldSet, newSet;
    
    for (int k = world_rank+1; k <= edge_sequence[0].n; k+=world_size) {
        p = *v0 + (edge_sequence[0].e * k);
        if (belongToIP(&p)) // 2.1
            newSet.insert(p);
    }
    oldSet = newSet;
    
    for(int i = 1; i < m - 1; i++) {
        for (int k = 1; k <= edge_sequence[i].n; k++) {
            if (world_rank == 0) {
                p = *v0 + (edge_sequence[i].e * k);
                if (belongToIP(&p)) // 2.1
                    newSet.insert(p);
            }
            for (auto u : oldSet) { // 2.2
                p = u + edge_sequence[i].e * k;
                if (belongToIP(&p))
                    newSet.insert(p);
            }
        }
        oldSet = newSet;
    }
    
    for (auto u : oldSet) // 3
        for (int k = 0; k < edge_sequence[m-1].n; k++) {
            p = u + edge_sequence[m-1].e * k;
            if (belongToIP(&p))
                newSet.insert(p);
        }
    
    //     return v0 ∈ Am ? 4
    bool result = newSet.find(*v0) != newSet.end();
    cout << "set size " << newSet.size() << endl;
    return result;
}

bool MPIPolyDecomp() {
    double start, end;
    start = MPI_Wtime();
    bool result = PolyDecomp();
    end = MPI_Wtime();
//    int start, end, length = edge_sequence[0].n, num = length % world_size, range = length/world_size;
//    
//    if (world_rank < num) {
//        start = world_rank * (range + 1) + 1;
//        end = start + range + 1;
//    }
//    else {
//        start = num * (range + 1) + (world_rank - num) * range + 1;
//        end = start + range;
//    }
//    cout << "start: " << start << ", end: " << end << ", n0 = " << edge_sequence[0].n << ", num: " << num << ", range: " << range << ", rank: " << world_rank << endl;

    cout << "rank " << world_rank << " took " << end - start << " seconds with result " << boolalpha << result << endl;
    return false;
}

bool between(Point *p, Point *p1, Point *p2) {
    if(p1->y < p->y && p->y < p2->y) // between
        return true;
    else if((p1->x <= p->x && p1->y == p->y) || (p->x <= p2->x && p->y == p2->y)) // extremities
        return true;
    else if(p1->x <= p->x && p->x <= p2->x && p1->y == p2->y) // on the same line
        return true;
    else
        return false;
}

int getSector(Point *p){
    for (int r = 0; r < world_size; r++) {
        int start = r*((xmax-xmin)/world_size) + xmin, end;
        r == world_size - 1 ? end = xmax : end = ((r+1)*((xmax-xmin)/world_size)) - 1 + xmin;
        if (p->x >= start && p->x <= end)
            return r;
    }
    cout << "-1\n";
    return -1;
}

void getLatticePoints() {
    int start = world_rank*((xmax-xmin)/world_size), end;
    world_rank == world_size - 1 ? end = xmax-xmin: end = ((world_rank+1)*((xmax-xmin)/world_size)) - 1;
    if (end < start) end = start;
    for (coord_t y = ymin; y <= ymax; y++)
        for (coord_t x = xmin + start; x <= xmin + end; x++) {
            Point v = {x, y};
            if (belongToIP(&v))
                latticePoints.push_back(v);
        }
}

Alg2TypeX PolyDecompNumX() {
    
    double startTime, endTime, PreTime = 0, ComputeTime = 0, SendTime = 0, RecvTime = 0;
    startTime = MPI_Wtime();
    
    int start = world_rank*((xmax-xmin)/world_size) + xmin, end;
    world_rank == world_size - 1 ? end = xmax : end = ((world_rank+1)*((xmax-xmin)/world_size)) - 1 + xmin;
    if (end < start) end = start;
    
    newA = new IntAndX*[ymax+1];
    IntAndX** oldA = new IntAndX*[ymax+1];
    
    for (int y = 0; y <= ymax; y++) {
        newA[y] = new IntAndX[end-start+1];
        oldA[y] = new IntAndX[end-start+1];
    }
    
    int rank = getSector(v0);
    
    if(world_rank == rank) {
        oldA[(int)v0->y][((int)v0->x)-start].u = 1;
        newA[(int)v0->y][((int)v0->x)-start].u = 1;
    }
    
    vector<int> *buffer = new vector<int>[world_size];
    
    endTime = MPI_Wtime();
    PreTime = endTime - startTime;
    
    int counter = 0;
    
    for (int i = 0; i < m; i++) {
        startTime = MPI_Wtime();
        for (int y = ymin; y <= ymax; y++)
            for (int x = start; x <= end; x++) {
                Point v = {x, y};
                if (oldA[y][x-start].u > 0) {
                    for (int k = 1; k <= edge_sequence[i].n; k++) {
                        Point vp = v + edge_sequence[i].e*k;
                        if (belongToIP(&vp)) {
                            counter++;
                            int xp = vp.x, yp = vp.y, sector = getSector(&vp);
                            if (sector == world_rank) {
                                newA[yp][xp-start].u += oldA[y][x-start].u;
                                newA[yp][xp-start].a.add({k, i});
                            }
                            else {
                                int elems[5] = {k,i,oldA[y][x-start] .u, xp,yp};
                                buffer[sector].insert(buffer[sector].end(), elems, elems + 5);
                            }
                        }
                    }
                }
            }
        endTime = MPI_Wtime();
        ComputeTime += (endTime - startTime);
        startTime = MPI_Wtime();
        /*****************      Send        *****************/
        for (int j = 0; j < world_size; j++) {
            if(j != world_rank) {
                int count = buffer[j].size();
                MPI_Send(&count, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                
                for (int chunk = 0; chunk < count/CHUNKSIZE; chunk++)
                    MPI_Send(&buffer[j][chunk*CHUNKSIZE], CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD);
                
                if(count%CHUNKSIZE > 0)
                    MPI_Send(&buffer[j][(count/CHUNKSIZE)*CHUNKSIZE], count%CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD);
            }
        }
        /****************************************************/
        endTime = MPI_Wtime();
        SendTime += (endTime - startTime);
        
        startTime = MPI_Wtime();
        /*****************       Receive        *****************/
        for (int j = 0; j < world_size; j++) {
            if(j != world_rank) {
                
                int count = -1, *array;
                MPI_Recv(&count, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                array = new int[count];
                
                for (int chunk = 0; chunk < count/CHUNKSIZE; chunk++)
                    MPI_Recv(&array[chunk*CHUNKSIZE], CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                if(count%CHUNKSIZE > 0)
                    MPI_Recv(&array[(count/CHUNKSIZE)*CHUNKSIZE], count%CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                /******** After Receiving ********/
                
                for (int k = 0; k < count; k+=5) {
                    int x = array[k+3], y = array[k+4];
                    newA[y][x-start].u += array[k+2];
                    newA[y][x-start].a.add({array[k], array[k+1]});
                }
                delete array;
            }
        }
        /*******************************************************/
        endTime = MPI_Wtime();
        RecvTime += (endTime - startTime);
        
        startTime = MPI_Wtime();
        for (int y = ymin; y <= ymax; y++)
            for (int x = 0; x <= end - start; x++)
                oldA[y][x].u = newA[y][x].u;
        
        for (int sec = 0; sec < world_size; sec++)
            if(sec != world_rank)
                buffer[sec].erase(buffer[sec].begin(), buffer[sec].end());
        endTime = MPI_Wtime();
        ComputeTime += (endTime - startTime);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    for (int i = 0; i <= ymax; i++)
        delete [] oldA[i];
    delete [] oldA;

    cout << "Pre time: " << PreTime << ", Compute time: " << ComputeTime << ", Send time: " << SendTime << ", Receive time: " << RecvTime << ", count: " << counter << endl;
    
    return {newA, newA[(int)v0->y][(int)v0->x - start].u};
}

Alg2TypeX MPIPolyDecompNumX() {
    double start, end;
    start = MPI_Wtime();
    Alg2TypeX ret = PolyDecompNumX();
    end = MPI_Wtime();
    
    
    int rank = getSector(v0);
    
    if(world_rank == rank)
        cout << "p" << world_rank << " finished in " << end - start << " with " << ret.u << " summands" << endl;
    else
        cout << "p" << world_rank << " finished in " << end - start << endl;
    return ret;
}

Point** CalcIP() {
    int convexSize = (int)ConvexHull.size();
    ymax = ymin = ConvexHull[0].y;
    for (int i = 1; i < convexSize; i++) {
        if (ConvexHull[i].y > ymax)
            ymax = ConvexHull[i].y;
        if (ConvexHull[i].y < ymin)
            ymin = ConvexHull[i].y;
        if (ConvexHull[i].x > xmax)
            xmax = ConvexHull[i].x;
        if (ConvexHull[i].x < xmin)
            xmin = ConvexHull[i].x;
    }
    
    IP = new Point*[ymax + 1];
    for (int i = 0; i <= ymax; i++)
        IP[i] = new Point[2];
    
    double slope[convexSize];
    for (int i = 0; i < convexSize - 1; i++)
        slope[i]=(ConvexHull[i].y-ConvexHull[i+1].y)/(ConvexHull[i].x - ConvexHull[i+1].x);
    
    for (intmax_t y = ymin; y <= ymax; y++) {
        int count = 0;
        double x[2] = {-1, -1};
        for (int p = 0; p < convexSize && count < 2; p++) {
            const int y1 = ConvexHull[p].y;
            const int y2 = ConvexHull[p+1].y;
            if (y >= min(y1, y2) && y <= max(y1, y2)) {
                if (abs(slope[p]) < 1e-2) {
                    x[0] = ConvexHull[p].x;
                    x[1] = ConvexHull[(p+1)%convexSize].x;
                    count = 2;
                }
                else {
                    double value = ((y-ConvexHull[p].y)/slope[p])+ConvexHull[p].x;
                    if (x[0] != value)
                        x[count++] = value;
                }
            }
        }
        const int x1 = (const int)x[0], x2 = (const int)x[1];
        IP[y][0].x = ceil(min(x1, x2));
        IP[y][1].x = floor(max(x1, x2));
        IP[y][0].y = IP[y][1].y = y;
    }
    
    if (IP[ymin][0].x == -1)
        IP[ymin][0].x = IP[ymin][1].x;
    if (IP[ymax][0].x == -1)
        IP[ymax][0].x = IP[ymax][1].x;
    
    return IP;
}

bool belongToIP(Point* p) {
    return p->y <= ymax && p->y >= ymin && p->x >= IP[(int)p->y][0].x && p->x <= IP[(int)p->y][1].x;
}

int gcd (int a, int b) {
    int c;
    while (a != 0) {
        c = a; a = b%a;  b = c;
    }
    return b;
}

void printer() {
    for (int i = 0; i < ymax; i++) {
        cout << "(" << IP[i][0].x << ", " << IP[i][0].y << ")" << "\t";
        cout << "(" << IP[i][1].x << ", " << IP[i][1].y << ")" << endl;
    }
}

void init(Polygon p) {
    ConvexHull = p;
    m = (int)ConvexHull.size() - 1;
    edge_sequence = new edge[m];
    
    for (int i = 0; i < m; i++) {
        coord_t dx = ConvexHull[i+1].x - ConvexHull[i].x, dy = ConvexHull[i+1].y - ConvexHull[i].y;
        int n = gcd(dx, dy);
        n > 0 ? edge_sequence[i] = {n, dx/n, dy/n} : edge_sequence[i] = {-n, dx/-n, dy/-n};
    }
    
    v0 = &ConvexHull[0];
    CalcIP();
    
    
    MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &world_size);
    
}

void clean() {
    delete [] edge_sequence;
    for (int i = 0; i <= ymax; i++)
        delete [] IP[i];
    delete [] IP;
    delete newA;
}

//Alg2Type PolyDecompNum() {
//    getLatticePoints();
//
//    int start = world_rank*(latticePoints.size()/world_size), end, counter;
//    world_rank == world_size - 1 ? end = latticePoints.size() - 1: end = ((world_rank+1)*(latticePoints.size()/world_size)) - 1;
//    if (end < start) end = start;
//    vector<Point>::const_iterator first = latticePoints.begin() + start;
//    vector<Point>::const_iterator last = latticePoints.begin() + (end+1);
//    vector<Point> newVec(first, last);
//
//    map<Point, IntAndSet> *oldA = new map<Point, IntAndSet>();
//    newA = new map<Point, IntAndSet>();
//
//    int rank = getSector(v0);
//
//    if(world_rank == rank)
//        (*oldA)[*v0].u = 1;
//
//    (*newA)[*v0].u = 1;
//
//    vector<int> *buffer = new vector<int>[world_size];
//
//    for (int i = 0; i < m; i++) {
//        for (auto v : newVec) {
//            if ((*oldA)[v].u > 0) {
//                for (int k = 1; k <= edge_sequence[i].n; k++) {
//                    Point vp = v + edge_sequence[i].e*k;
//                    if (belongToIP(&vp)) {
//                        counter++;
//                        int sector = getSector(&vp);
//                        if (sector == world_rank) {
//                            (*newA)[vp].u += (*oldA)[v].u;
//                            (*newA)[vp].S.insert({k,i});
//                        }
//                        else {
//                            buffer[sector].push_back(k);
//                            buffer[sector].push_back(i);
//                            buffer[sector].push_back((*oldA)[v].u);
//                            buffer[sector].push_back(vp.x);
//                            buffer[sector].push_back(vp.y);
//                        }
//                    }
//                }
//            }
//        }
//        /*****************      Send        *****************/
//        for (int j = 0; j < world_size; j++) {
//            if(j != world_rank) {
//                int count = buffer[j].size();
//                MPI_Send(&count, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
//
//                //                cout << "sent count = " << count << endl;
//
//                for (int chunk = 0; chunk < count/CHUNKSIZE; chunk++)
//                    MPI_Send(&buffer[j][chunk*CHUNKSIZE], CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD);
//
//                if(count%CHUNKSIZE > 0)
//                    MPI_Send(&buffer[j][(count/CHUNKSIZE)*CHUNKSIZE], count%CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD);
//            }
//        }
//        /************************************************/
////        MPI_Barrier(MPI_COMM_WORLD);
//
//        /*****************       Receive        *****************/
//        for (int j = 0; j < world_size; j++) {
//            if(j != world_rank) {
//
//                int count = -1, *array;
//                MPI_Recv(&count, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                //                cout << "receive at " << world_rank << " from " << j << " = " << count << endl;
//                array = new int[count];
//
//                for (int chunk = 0; chunk < count/CHUNKSIZE; chunk++)
//                    MPI_Recv(&array[chunk*CHUNKSIZE], CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//                if(count%CHUNKSIZE > 0)
//                    MPI_Recv(&array[(count/CHUNKSIZE)*CHUNKSIZE], count%CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//                /******** After Receiving ********/
//
//                for (int k = 0; k < count; k+=5) {
//                    Point vx = {(coord_t)array[k+3], (coord_t)array[k+4]};
//                    (*newA)[vx].u += array[k+2];
//                    (*newA)[vx].S.insert({array[k], array[k+1]});
//                    //                    cout << array[k] << ", " << array[k + 1] << ", " << array[k + 2] << ", " << array[k + 3] << ", " << array[k + 4] << ", " << endl;
//                }
//                delete array;
//            }
//        }
////        MPI_Barrier(MPI_COMM_WORLD); // think about it
//        /************************************************/
//        *oldA = *newA;
//        for (int sec = 0; sec < world_size; sec++)
//            if(sec != world_rank)
//                buffer[sec].erase(buffer[sec].begin(), buffer[sec].end());
//    }
//
//    delete oldA;
//
//    cout << "rank " << world_rank << ": " << counter << endl;
//
//    return {newA, (*newA)[*v0].u};
//}