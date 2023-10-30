#pragma once

#include <iostream>
#include <getopt.h>
#include <vector>
#include <string>
#include <bitset>
#include <fstream>
#include <sstream>
#include <random>
#include <time.h>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <iterator>
#include <memory>
#include <tuple>
#include <unistd.h>

using namespace std;

typedef tuple<int, int*, int> grid;

typedef unordered_map<int, vector<grid>> sig_grid_Map;

typedef struct Config
{
    int Pc = 500;
    int generation = 150;
    int max_descendant = 500;
    int ymin = 5;
    int rank = 1;
    bool to_optimize = true;
    string infile;
    string outfile;
} Config;

typedef struct bounding_box
{
    int left;
    int right;
    int bottom;
    int top;

    bounding_box(int left, int right, int bottom, int top): left(left),
                                                            right(right),
                                                            bottom(bottom),
                                                            top(top) {};
} bounding_box;


typedef struct point2D 
{   
    int x; 
    int y;

    point2D(int x, int y) : x(x), y(y) {};
} point2D;


typedef struct point3D 
{   
    int x; 
    int y;
    int z;

    point3D(int x, int y, int z) : x(x), y(y), z(z) {};
} point3D;


typedef struct visit_table
{
    vector<char> t;
    int nrow, ncol;

    visit_table() = default;
    visit_table(int nrow, int ncol): nrow(nrow), ncol(ncol) { t.resize(nrow * ncol * 2, 0); }
    
    void clear() { fill(t.begin(), t.end(), 'x'); }
    void addRows(int n) { t.resize(t.size() + n*ncol*2); nrow+=n; }
    void deleteRows(int n) { t.resize(t.size() - n*ncol*2); nrow-=n; }
    void addCols(int n) { t.resize(t.size() + n*nrow*2); ncol+=n; }
    char& operator()(int x, int y, int z) { return t[y * ncol * 2 + x * 2 + z]; }
} Visit_Table;


#include "parsing.hpp"
#include "io.hpp"