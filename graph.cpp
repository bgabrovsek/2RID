//
// Created by Bostjan on 19/03/2020.
//

#include "graph.h"

/**
 * Graph constructor, creates an edgeless graph
 * @param n number of vertices
 */
Graph::Graph(u n) {
    this->n = n;
}

// Add edge to graph
Graph &Graph::operator += (u edge) {
//cout << "add" << edge << endl;
    this->edges.insert(edge);
    return *this;
}

/**
 * Returns a list of neighbourhood vertices, if filter is given, returns vertices {i} that have i'th bit 1 in filter
 * @param vertex
 * @return
 */
set<u> Graph::neighbourhood(u vertex, u64 filter) {
    set<u> nbh;
    bool no_filter = (filter == INF64);
    for(const u &e : edges) {
        if (( FIRST(e) == vertex) and ((no_filter) or (BIT(filter,SECOND(e))))) nbh.insert(SECOND(e));
        if ((SECOND(e) == vertex) and ((no_filter) or (BIT(filter, FIRST(e))))) nbh.insert(FIRST(e));
    }
    return nbh;
}


/**
 * returns degree of vertex, but counts only vertices in filter (v bit 1), if given
 */
u Graph::degree(u vertex, u64 filter) {
    u d = 0;
    bool no_filter = (filter == INF64);
    for(const u &e : edges)
        if (((FIRST(e) == vertex) and ((no_filter) or (BIT(filter, SECOND(e)))))
            or ((SECOND(e) == vertex) and ((no_filter) or (BIT(filter, FIRST(e)))))) d++;
    return d;
}


/**
 * @return list of vertex degrees
 */
vector<u> Graph::degrees() {
    vector<u> d(this->n, 0);
    for(const int &e : edges) {
        d[FIRST(e)]++;
        d[SECOND(e)]++;
    }
    return d;
}

/**
 * For a partition (or subgraphs) {A,B}, returns list of vertices that is connected to the other partition
 * @param filter_A
 * @param filter_B
 * @return
 */
set<u> Graph::cut_vertices(u64 filter_A, u64 filter_B) {
    set<u> v;
    for(const u &e : edges) {
        if ((BIT(filter_A, FIRST(e)) and BIT(filter_B, SECOND(e))) or
            (BIT(filter_B, FIRST(e)) and BIT(filter_A, SECOND(e)))) {
            v.insert(FIRST(e));
            v.insert(SECOND(e));
        }
    }
    return v;
}



// Prints a set of edges
std::ostream& operator<<(std::ostream& os, const set<u>& edges)
{
    for(const u &e : edges)
        os << "(" << FIRST(e) << "," << SECOND(e) << ") ";
    return os;
}

// Prints a graph
std::ostream& operator<<(std::ostream& os, const Graph& g)
{
    os << "Graph " << g.edges << "["<< g.n << "]" <<endl;
    return os;
}

// Returns cylce C_a × C_b
Graph * CycleCycle(u a, u b, u start_index) {
    auto *G = new Graph(a*b);
    for (u bi = 0; bi < b ; bi++) {
        for (u ai = 0; ai < a; ai++) {
            *G += EDGE(a * bi + ai, a * bi + (ai + 1) % a);
            *G += EDGE(a * bi + ai, a * ((bi + 1) % b) + ai);
        }
    }
    return G;
}

// Returns C_a × P_b
Graph * CyclePath(u a, u b, u start_index) {
    auto *G = new Graph(a*b);
    for (u bi = 0; bi < b ; bi++) {
        for (u ai = 0; ai < a; ai++) {
            *G += EDGE(a * bi + ai, a * bi + (ai + 1) % a);
            if (bi+1 < b)
                *G += EDGE(a * bi + ai, a * ((bi + 1) % b) + ai);
        }
    }
    return G;
}

// Returns P_a × P_b
Graph * PathPath(u a, u b, u start_index) {
    auto *G = new Graph(a*b);
    for (u bi = 0; bi < b ; bi++) {
        for (u ai = 0; ai < a; ai++) {
            if (ai+1 < a)
                *G += EDGE(a * bi + ai, a * bi + (ai + 1) % a);
            if (bi+1 < b)
                *G += EDGE(a * bi + ai, a * ((bi + 1) % b) + ai);
        }
    }
    return G;
}

/**
 * Returns the Petersen monograph
 * @param k Petersen graph P(n,k) parameter
 * @param n size of monograph (monograph will have 2n vertices)
 * @param start_index vertices will not start from 0, but start_index
 */
Graph * PetersenMonograph(u k, u n, u start_index) {
    auto *G = new Graph(2*n);
    for (u i=0; i < n; i++) {
        *G += EDGE(2*i,2*i+1);
        if (i+1 < n)    *G += EDGE(2*i+1, 2*i+3);
        if (i + k < n)  *G += EDGE(2*i,2*i+2*k);
    }
    return G;
}
