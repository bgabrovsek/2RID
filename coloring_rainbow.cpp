//
// Created by Bostjan on 19/03/2020.
//

#include <iostream>
#include "coloring_rainbow.h"
#include "defs.h"
#include "domination.h"
#include <math.h>
using namespace std;

/**
 * @return conjectured values for irrdf(P(n,k))
 */
u I2RD_conjecture(u n, u k, u r) {
    if ((k == 2) and (r==2))
        return vector<u>({0,3,2,2,2,2,3,2,2,1})[n%10] + (u)ceil(.8*n);
    if ((k == 3) and (r==2))
        return vector<u>({0,2,0,2,0,1,2,1,1,2,1,2,1,1,0,1})[n%16] + (u)ceil(.875*n);
    if ((k == 0) and (r==3))
        return vector<u>({0,4,3,3,4,3,3,2,3,2})[n%10] +  (u)ceil(1.2*n);
    return 0;
}


u IPC(u n, u r) {
    if (n%2 == 0) return n;
    return n+1;
}


/**
 * is the PIRD fucntion extendable (left and right)
 * @param g extended graph
 * @param colors_ex extension colors
 * @param color PIRD fucntion
 * @param r number of colors
 * @param f_deg returns defree of vertex
 * @return true/false

bool PIRD_extendableQ(Graph *g, const vector<vector<u>> & colors_ex, const vector<u> & color, u r, fuu f_deg) {

    u64 filter = ONES(g->n);
    for (const auto &c_ex_in : colors_ex)
        for (const auto &c_ex_out : colors_ex) {
            bool good = true;
            vector<u> c = concatenate({c_ex_in, color, c_ex_out});
            // check independence
            for (u v = 0; v < g->n; v++) {
                if ((c[v]) and nbh_color(g, c, v, filter)) good = false;
                if ((c[v] == 0) and (g->degree(v, filter) == f_deg(v)) and (nbh_color(g,c,v,filter) < ONES(r))) good = false;
            }
            if (good) return true;
        }
    return false; // no good colorings found
}


// not used
vector<vector<u>> PIRD_filter_extendable(Graph *g, const vector<vector<u>> & colors_ex, const vector<vector<u>> & colors, u r, fuu f_deg) {
    vector<vector<u>> colorings;
    for (const auto & c : colors)
        if (PIRD_extendableQ(g, colors_ex, c, r, f_deg)) colorings.push_back(c);
    return colorings;
}
 */

/**
 * Determines if c is a good Independent Rainbow domination coloring for vertex v
 * @param g full graph (IN, MID, OUT)
 * @param color_in colors for IN monograph
 * @param coloring (partial) colors for MID monograph
 * @param color_out colors for OUT monograph
 * @param vertex new vertex
 * @param color color of new vertex
 * @param filter 64-bit array, which vertices to include (some may not yet have been colored)
 * @return true if coloring is OK, false if not
 */
// TODO: do we really need this function seperately?
bool good_PIRD_color_for_vertex(Graph *g,
                                const vector<u> &color_in, const vector<u> &coloring, const vector<u> &color_out,
                                u n_io, u n_mid, u vertex, u color, u64 filter, u r, fuu f_deg, bool check_nn) {

    // is the new vertex color good to add (i.e. defines a PIRDF)?
    // check independence


    if ((color) & (nbh_tricolor(g, color_in, coloring, color_out, vertex, filter))) { return false;} // NEODVISNOST

    // check if all colors are in neighbourhood
    if ((not color) and (g->degree(vertex, filter) == f_deg(vertex)) and
        (nbh_tricolor(g, color_in, coloring, color_out, vertex, filter) < ONES(r)))
        return false;

    if (not check_nn) return true;
    // check if neighbours are fine
    for (const auto &v : g->neighbourhood(vertex, filter)) { // loop through adjacent vertices
        if ((TRICOLOR(color_in, coloring, color_out, n_io, n_mid, v) == 0) and // vertex is empty
            (g->degree(v, filter) + (check_nn ? 1 : 0) == f_deg(v)) and //
            ((nbh_tricolor(g, color_in, coloring, color_out, v, filter) | color) < ONES(r))) {
            return false;
        }
    }
    return true;
}




/**
 * Computes all partial r-rainbow domination colorings of a graph (IN,MID,OUT), where IN and OUT are given
 * @param g the graph
 * @param color_in coloring of IN part of the graph
 * @param color_out coloring of OUT part of the graph
 * @param r number of colors
 * @return vector of colorings of MID part (? check ?)
 */
// TODO: should be parallelized for large graphs?
vector<vector<u>> Partial_Independent_Rainbow_Domination_Colorings(Graph* g, const vector<u> &color_in, const vector<u> &color_out, u r, fuu f_deg) {
    vector<vector<u>> colorings(1, vector<u>()); // vector of all possible colorings, init with trivial coloring of zero vertices

    u n_io = color_in.size(); // size of IO matrices
    u n_mid = g->n - (n_io<<1u); // size of middle matrix

    u64 filter_in = ONES(n_io);
    u64 filter_out = (ONES(n_io) <<(n_io+n_mid));

    // check if color_in and color_our are compatible (if there are edges between them)
    set<u> cut = g->cut_vertices(filter_in, filter_out);
    for (const u& v : cut) {
        if (not good_PIRD_color_for_vertex(g, color_in, vector<u>(), color_out,
                                           n_io, n_mid, v, (v < n_io) ? color_in[v] : color_out[v-n_io-n_mid],
                                           CLEAR_BIT(filter_in | filter_out,v), r, f_deg, false)) {
            return vector<vector<u>>(); // return empty vector
        }
    }

    u64 filter = filter_in | filter_out; // add all new vrtices in mid graph
    for (u vertex=n_io; vertex < n_io+n_mid; vertex++) { // loop through vertices in middle graph
        vector<vector<u>>new_colorings; // colorings of current (next) depth
        for (auto & coloring : colorings) { // loop through existing colorings of vertices 0,...,n-2.
            //for (u col = 0; col < (1u<<r); col++) {// loop through the powerset (if there are no colors in the neighbourhood)
            for (u col = 0; col <= r; col++) {// loop through the powerset (if there are no colors in the neighbourhood)
                u c = col ? (1u << (col-1)) : 0; // only one color
                if (good_PIRD_color_for_vertex(g, color_in, coloring, color_out,n_io, n_mid, vertex, c, filter, r, f_deg, true)) {
                    vector<u> new_coloring(coloring); // copy old coloring in new one
                    new_coloring.push_back(c); // add current color
                    new_colorings.push_back(new_coloring); // add graph coloring to all colorings
                }
            } // powerset
        } // colorings
        colorings = new_colorings; // replace old n-1 colorings with the new n colorings
        filter = SET_BIT(filter, vertex); // include vertex in new filter
    } // vertices
    return colorings;
}

/**
 * Check if coloring is good
 * @param g graph
 * @param c coloring
 * @param r numbder of colors
 * @param f_deg degree of vertices
 * @return true/false
 */

bool GoodPIRDQ(Graph *g, const vu & c, u r, fuu f_deg) {

    u64 filter = ONES(g->n);
    if (c.size() != g->n) cout << "Coloring lengths mismatch." << endl;

    for (u v = 0; v < g->n; v++) {
        if ((c[v]) &  nbh_color(g, c, v, filter)) return false; // NEODVISNOST
        if ((c[v] == 0) and (g->degree(v, filter) == f_deg(v)) and (nbh_color(g,c,v,filter) < ONES(r))) return false;
    }
    return true;

}

/*
 *
 * P(2,2) -> tr(A*M^0) = 65535 (4)
P(4,2) -> tr(A*M^1) = 6 (6)
P(6,2) -> tr(A*M^2) = 6 (8)
P(8,2) -> tr(A*M^3) = 8 (9)
P(10,2) -> tr(A*M^4) = 8 (8)
P(12,2) -> tr(A*M^5) = 12 (12)
P(14,2) -> tr(A*M^6) = 14 (14)
P(16,2) -> tr(A*M^7) = 14 (16)
P(18,2) -> tr(A*M^8) = 16 (17)
P(20,2) -> tr(A*M^9) = 16 (16)
P(22,2) -> tr(A*M^10) = 20 (20)
P(24,2) -> tr(A*M^11) = 22 (22)
P(26,2) -> tr(A*M^12) = 22 (24) = P(16,2) -> tr(A*M^7) + [8]
It holds: tr(A*M^(n+5)) = tr(A*M^n) + [8] for n >= 7.
P(3,2) -> tr(A*M^0) = 4 (5)
P(5,2) -> tr(A*M^1) = 6 (6)
P(7,2) -> tr(A*M^2) = 7 (8)
P(9,2) -> tr(A*M^3) = 8 (9)
P(11,2) -> tr(A*M^4) = 12 (12)
P(13,2) -> tr(A*M^5) = 13 (13)
P(15,2) -> tr(A*M^6) = 14 (14)
P(17,2) -> tr(A*M^7) = 15 (16)
P(19,2) -> tr(A*M^8) = 16 (17)
P(21,2) -> tr(A*M^9) = 20 (20)
P(23,2) -> tr(A*M^10) = 21 (21)
P(25,2) -> tr(A*M^11) = 22 (22) = P(15,2) -> tr(A*M^6) + [8]
 */