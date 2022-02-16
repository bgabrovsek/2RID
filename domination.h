//
// Created by Bostjan on 05/04/2020.
//
// GENERAL DOMINATION FUNCTIONS
//

#ifndef RAINBOW_DOMINATION_H
#define RAINBOW_DOMINATION_H

#include <iostream>
#include <functional>
#include "defs.h"
#include "graph.h"
#include "matrix.h"

using namespace std;

// if graph (A,B,C) has coloring ca+cb+cc, where |ca| = a and |cb| = b, return the color of v
#define TRICOLOR(ca,cb,cc,a,b,v) (( (v) < (a) ) ? (ca[v]) : ( ( (v) < ((a) + (b)) ) ? ( cb[(v)-(a)]) : ( cc[(v)-(a)-(b)]) ))

// common degree functions
#define CUBIC [](float x){return 3;}
#define QUADRATIC [](float x){return 4;}
#define DVA [](float x){return 2;}

// function templates
typedef const function< u (Graph *, const vector<u> &, u, fuu) > &func_weight; // coloring specific weight function
typedef const function<vector<vector<u>> (Graph *, const vector<u> &, const vector<u> &, u, fuu) > &func_partial_coloring;
typedef const function<vector<vector<u>> (Graph *, const vector<vector<u>> &, const vector<vector<u>> &, u, fuu) > &func_filter_ex;
typedef const function<bool(Graph *, const vu &, u, fuu) > & func_good_coloring;

//typedef const function< vector<vector<u>> (Graph *g, const vector<u> &, const vector<vector<u>> &, u, fu) > &func_partial_coloring_bunch;
// general use functions
u nbh_color(Graph *g, const vector<u> &coloring, int vertex, u64 filter = INF64);
u nbh_tricolor(Graph *g, const vector<u> &col_in, const vector<u> &col_mid, const vector<u> &col_out, int vertex, u64 filter = INF64);



matrix * Coloring_Matrix_Serial(Graph *g, const vector<vector<u>> &color_in, const vector<vector<u>> &color_out, u r, fuu f_deg, func_partial_coloring f_coloring,
Graph *g_ex, const vector<vector<u>> &color_ex, func_filter_ex f_filter);
matrix * Coloring_Matrix_Parallel(Graph *g, const vector<vector<u>> &colorings, u r, fuu f_deg, func_partial_coloring f_coloring, int N_THREAD,
        Graph *g_ex, const vector<vector<u>> &color_ex, func_filter_ex f_filter);

// matix powers
void compute_matrix_powers(vvu &coloring, matrix &start_matrix, matrix &m, fsu f, fuu e, int max_power=64);
void compute_matrix_powers(vvu &coloring, matrix &m, fsu f, fuu e, int max_power = 64);
// extended colorings
vector<vector<u>> Middle_Colorings(Graph *g, const vvu &ex_in_cols, const vvu &ex_out_cols, u r, fuu f_deg, func_partial_coloring f, vvvu & ret_in, vvvu & ret_out);


matrix * MatrixEx(Graph *g, Graph * g_ex,
const vvu & colorings,
const vvvu & ex_in, const vvvu & ex_out,
        u r, fuu f_deg,  func_good_coloring good);

#endif //RAINBOW_DOMINATION_H