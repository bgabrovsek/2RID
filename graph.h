//
// Created by Bostjan on 19/03/2020.
//

#ifndef RAINBOW_GRAPH_H
#define RAINBOW_GRAPH_H

#include <iostream>
#include <set>
#include <vector>
#include "defs.h"


using namespace std;

#define EDGE_BITS 8u
#define EDGE_FIRST 0xFFu

#define EDGE(a, b) ( ((a)<(b)) ? (((a)<<EDGE_BITS)+(b)) : (((b)<<EDGE_BITS)+(a)) ) // edge {a,b}
#define FIRST(e) ((e) >> EDGE_BITS) // returns a from {a,b | a < b}
#define SECOND(e) ((e) & EDGE_FIRST) // returns b from {a,b | a < b}
#define HAS(e, v) ((FIRST(e) == v) or (SECOND(e) == v)) // does edge e have vertex v?

class Graph {
public:
    u n;
    set<u> edges;

    Graph(u n = 0);
    ~Graph() {cout << "DEL GRAPH" << endl;}

    Graph &operator += (u edge);

    set<u> neighbourhood(u vertex, u64 filter = INF64);
    u degree(u vertex, u64 filter = INF64);
    vector<u> degrees();

    set<u> cut_vertices(u64 filter_A, u64 filter_B);

};


std::ostream& operator<<(std::ostream& os, const set<u>& edges);
std::ostream& operator<<(std::ostream& os, const Graph& g);

// SPECIAL GRAPHS

Graph * CycleCycle(u a, u b, u start_index = 0);
Graph * CyclePath(u a, u b, u start_index = 0);
Graph * PathPath(u a, u b, u start_index = 0);

Graph * PetersenMonograph(u k, u n, u start_index = 0);

#endif //RAINBOW_GRAPH_H