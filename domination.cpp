//
// Created by Bostjan on 05/04/2020.
//

// OK
#include "domination.h"
#include <future>
#include <algorithm>
//#include <thread>         // std::thread
#include <deque>
#include "defs.h"


/**
 * get colors of neighbours
 * @param g graph
 * @param coloring coloring vector
 * @param vertex
 * @param only_smaller_vertices consider only vertices smaller than vertex
 * @return binary colors
 */
u nbh_color(Graph *g, const vector<u> &coloring, int vertex, u64 filter) {
    u color = 0;
    for (const u &v : g->neighbourhood(vertex, filter))
        color |= coloring[v];
    return color;
}


/**
 * get colors of neighbours, where graph has colors (IN, MID +[0], OUT).
 * @param g the graph (IN, MID +[0], OUT).
 * @param col_in coloring for IN part
 * @param col_mid coloring for MID part (can be incomplete)
 * @param col_out coloring for OUT part
 * @param vertex which vertex we are considering
 * @param filter bool array of which vertices are colored
 * @return
 */
//
u nbh_tricolor(Graph *g, const vector<u> &col_in, const vector<u> &col_mid, const vector<u> &col_out, int vertex, u64 filter) {
    u color = 0;
    for (const u &v : g->neighbourhood(vertex, filter))
        color |= TRICOLOR(col_in, col_mid, col_out, col_in.size(), g->n - ((u)col_in.size() << 1), v);
    return color;
}


/**
 * Weight w/o checking if coloring is good
 * @param coloring coloring function
 * @return weight of coloring
 */
// TODO: do we need to pass a custom weight function for general dominations?
u weight(const vector<u> &coloring) {
    u w = 0;
    for (auto& c : coloring) w += POPCOUNT(c);
    return w;
}

/**
 * minimal weight w/o checking if coloring is goof
 * @param colorings vector of colorigs
 * @return minimal weight
 */
u min_weight(const vector<vector<u>> &colorings) {
    u min_w = INF;
    for (auto &c : colorings) {
        u w = weight(c);
        if (w < min_w) min_w = w;
    }
    return min_w;
}

/**
 * returns a good subset of colorings
 * @param g graph
 * @param colorings vector of colorings
 * @param r number of colors
 * @param f_deg degree function
 * @param f_w weight function
 * @return
 */
vector<vector<u>> filter_colorings(Graph *g, const vector<vector<u>> &colorings, u r, fuu f_deg, func_weight f_w) {
    vector<vector<u>> ird;
    for (const auto &coloring : colorings)
        if (f_w(g, coloring, r, f_deg) != INF)
            ird.push_back(coloring);
    return  ird;
}


/**
 * COmputes colorings of middle (missing graphs for mod > 0)
 * @param g graph
 * @param ex_in_cols coloring extension left
 * @param ex_out_cols coloring extension right
 * @param r number of colors
 * @param f_deg degrees of vertices
 * @param f coloring function
 * @param ret_in
 * @param ret_out
 * @return
 */
vector<vector<u>> Middle_Colorings(Graph *g, const vvu &ex_in_cols, const vvu &ex_out_cols, u r, fuu f_deg, func_partial_coloring f, vvvu & ret_in, vvvu & ret_out) {

    vvu colorings;
    // get all colorings for IN graph
    int ind = 0;
    for (const auto & ex_in : ex_in_cols) {
        for (const auto &ex_out : ex_out_cols) {
            //      cout << "color: " << ex_in  << " + " << ex_out << endl;
            auto cs = f(g, ex_in, ex_out, r, f_deg); // partial colorings
            if (cs.size() == 0) {
                //     cout << " not good." << endl;
                continue;
            }


            for (const auto & c : cs) {
                // find if new coloring already in vector
                auto it = std::find(colorings.begin(), colorings.end(), c);
                if (it == colorings.end()) {
                    // new coloring

                    colorings.push_back(c);
                    ret_in.push_back({ex_in});
                    ret_out.push_back({ex_out});
                } else {
                    // coloring already found
                    int c_index = distance(colorings.begin(), it);
                    if (not InQ(ret_in[c_index], ex_in)) ret_in[c_index].push_back(ex_in);
                    if (not InQ(ret_out[c_index], ex_out)) ret_out[c_index].push_back(ex_out);
                }
            }
        }
    }
    return colorings;
}


///// MATRICES /////

/**
 * Returns the coloring matrix of a graph (non parallel version)
 * @param g garph of the form (IN, MID, OUT)
 * @param colorings vector of colorings for graphs IN and OUT
 * @param r number of colors (or parameter)
 * @param f_deg function u->u that returns degree of a vertex
 * @param f_coloring coloring function for MID graph
 * @return returns the matrix as a vector of size colorings[0].size() × colorings[0].size()
 */
matrix * Coloring_Matrix_Serial(Graph *g, const vector<vector<u>> &color_in, const vector<vector<u>> &color_out, u r, fuu f_deg, func_partial_coloring f_coloring,
                                Graph *g_ex, const vector<vector<u>> &color_ex, func_filter_ex f_filter) {

    auto mat = new matrix(); // the matrix
    u w_out, w_mid; // weights
    int ci = 0;

    for (const auto &col_in : color_in) { // loop through rows (IN colors)
        vector<u> row;
        row.reserve(color_in.size());
        for (const auto &col_out : color_out) { // loop through columns (OUT colors)
            w_out = weight(col_out);
            auto middle_colorings = f_coloring(g, col_in, col_out, r, f_deg);

            // filter middle colorings so they are extendable
            if (middle_colorings.size() == 0)
                w_mid = INF;
            else {
                // has middle colorings
                auto ex_middle_colorings = trim(f_filter(g_ex, color_ex, concatenate(col_in, middle_colorings, col_out), r, f_deg),col_in.size(),col_out.size());
                w_mid = min_weight(ex_middle_colorings);
            }
            row.push_back((w_mid == INF) ? INF :  (w_out + w_mid)); //(w_out + w_in + w_mid)); // exclude w_in
        }
        mat->push_row(row);

        ci++;
    }
    return mat;
}


/**
 * Returns the coloring matrix of a graph (parallel w/ multi-threading)
 * @param g garph of the form (IN, MID, OUT)
 * @param colorings vector of colorings for graphs IN and OUT
 * @param r number of colors (or parameter)
 * @param f_deg function u->u that returns degree of a vertex
 * @param f_coloring coloring function for MID graph
 * @param N_THREAD number of threads to use
 * @return returns the matrix as a vector of size colorings[0].size() × colorings[0].size()

matrix * Coloring_Matrix_Parallel(Graph *g, const vector<vector<u>> &colorings, u r, fuu f_deg, func_partial_coloring f_coloring, int N_THREAD,
                                  Graph *g_ex, const vector<vector<u>> &color_ex, func_filter_ex f_filter) {

    auto mat = new matrix(); // the matrix

    u w_in, w_out, w_mid; // weights

    vector<vector<vector<u>>> col_bunch = split_to_bunch(colorings, N_THREAD);
    vector<future<matrix *>> futures(N_THREAD);


    // parallel compute weights
    for (int Q=0; Q<N_THREAD; Q++)
        //futures[Q] = async(launch::async, Bunch_Minimal_Coloring, g, col_in, col_bunch[Q], r, f_deg, f_coloring);
        futures[Q] = async(launch::async, Coloring_Matrix_Serial, g, col_bunch[Q], colorings, r, f_deg, f_coloring, g_ex, color_ex, f_filter);


    // get results from threads

    for (int Q=0; Q < N_THREAD; Q++)
        mat->append(futures[Q].get());
    //for_each(result.begin(), result.end(), [=](u& elt) { if (elt != INF) elt+= w_in;}); // add w_in to each entry of result, use [&] to capture by reference
    //mat.insert(mat.end(), result.begin(), result.end()); // insert results into matrix

    return mat;

}

 */
matrix * MatrixEx(Graph *g, Graph * g_ex,
                  const vvu & colorings,
                  const vvvu & ex_in, const vvvu & ex_out,
                  u r, fuu f_deg, func_good_coloring good) {


    auto mat = new matrix(); // the matrix
    u w;
    int ind_a = 0;

    u m = g_ex->n - 5*g->n/2;

    cout << "Computing matrix (mod " << m/2<< ")" << endl;
    for (const auto &a : colorings) { // loop through rows (IN colors)
        vu row;
        int ind_b = 0;
        for (const auto & b : colorings) {

            int reason = 0;

            auto part1 = vu(a.begin() + a.size()/2+m,a.end());
            auto part2 = vu(b.begin(), b.begin()+ b.size()/2-m);
            if (part1 == part2) {
                // matrices compatible
                bool extendable = false;
                //check if extendable
                for (const auto &aa : ex_in[ind_a])
                    for (const auto &bb : ex_out[ind_b]) {
                        auto penta_coloring = concatenate({aa,a,vu(b.begin()+ b.size()/2-m, b.end()),bb});
                        extendable |= good(g_ex, penta_coloring, r, f_deg);
                    }

                if (extendable) {
                    w = weight( vu(b.begin(), b.begin() + b.size()/2)  ) + weight( vu(a.begin()+a.size()/2, a.begin()+a.size()/2+m));
                }
                else {
                    reason = 1;
                    w = INF;
                }
            }
            else {w = INF; reason = 2;}

            cout << "row in " << a << " col out " << b << " val " << w << " (reason " << reason << ")" << endl;
            row.push_back(w);
            ind_b++;
        }

        mat->push_row(row);
        ind_a++;
    }
    return mat;
}






/**
 * Computes matrix powers start_matrix(×)m^n for n = 0,...,max_power and prints out minimal diagonal element
 * @param start_matrix
 * @param m
 * @param max_power
 */
void compute_matrix_powers(vvu &coloring, matrix &start_matrix, matrix &m, fsu f, fuu e, int max_power) {

    matrix mt = m.transpose(); // transpose the matrix
    matrix result = start_matrix;
    deque<matrix> history;
    for (int i=0; i < max_power; i++) {
        if (i) result = result * mt;

        u tr_value = result.trace();
        u tr_index = result.trace(true);

        cout << f(i)  <<" -> tr(A*M^" << i << ") = " << tr_value << " (" <<e(i)<< ")";



        // BACKTRACK COLORING
        u ind_row = tr_index;
        u ind_col = tr_index;
        u val = tr_value;
        int sz = coloring[tr_index].size() /2;
        vu COL;


        for (int hi=0; hi < history.size(); hi++) {
            // we want to analyse row/col index
            auto mat = history[hi];
            vu row = mat.get_row(ind_row);
            vu col = mt.get_row(ind_col);

            for (int a = 0; a < row.size(); a++) {
                //cout << row[a] + col[a] << " ";
                if (row[a] + col[a] == val) {
                    COL.insert(COL.begin(), coloring[a].begin(), coloring[a].begin()+sz);
                    //cout << "  " << coloring[a] << endl;
                    ind_col = a;
                    val = row[a];
                    break;
                }
                //cout << endl;
            }
        }
        COL.insert(COL.begin(), coloring[tr_index].begin(), coloring[tr_index].begin()+sz);

        // print in nice form
        cout << " coloring: ";
        for (int ind = 0; ind < COL.size(); ind++) {
            cout << COL[ind];
            if (ind % sz == sz-1) cout << " ";
        }
        cout << endl;
        //cout << " coloring: " << COL << endl;
        //cout << "COL " << COL << endl;




        for (int hi=0; hi < history.size(); hi++) {
            s diff = result - history[hi];
            if (diff < SINF) {
                cout << " = " << f(i-hi-1) << " -> tr(A*M^" << (i-hi-1) << ") + [" << diff << "]" << endl;
                cout << "It holds: tr(A*M^(n+" << (hi+1) << ")) = tr(A*M^n) + [" << diff << "] for n >= "<< (i-hi-1) <<"."<<endl;
                return;
            }
        }
        //cout << endl;
        history.push_front(result);
    }
    cout << "No repetition found." << endl;
}





/**
 * Computes matrix powers m^n for n = 0,...,max_power and prints out minimal diagonal element without star matrix
 * @param start_matrix
 * @param m
 * @param max_power
 */
void compute_matrix_powers(vvu &coloring, matrix &m, fsu f, fuu e, int max_power) {

    matrix mt = m.transpose(); // transpose the matrix

    matrix result = *identity(m.size());
    deque<matrix> history;
    for (int i=0; i < max_power; i++) {
        if (i) result = result * mt;

        u tr_value = result.trace();
        u tr_index = result.trace(true);

        cout << f(i)  <<" = \\tr(M^{" << i << "}) = " << tr_value << "$\\\\";



        // BACKTRACK COLORING
        u ind_row = tr_index;
        u ind_col = tr_index;
        u val = tr_value;
        int sz = coloring[tr_index].size() /2;
        vu best_coloring;


        for (int hi=0; hi < history.size(); hi++) {
            // we want to analyse row/col index
            auto mat = history[hi];
            vu row = mat.get_row(ind_row);
            vu col = mt.get_row(ind_col);

            for (int a = 0; a < row.size(); a++) {
                //cout << row[a] + col[a] << " ";
                if (row[a] + col[a] == val) {
                    best_coloring.insert(best_coloring.begin(), coloring[a].begin(), coloring[a].begin()+sz);
                    //cout << "  " << coloring[a] << endl;
                    ind_col = a;
                    val = row[a];
                    break;
                }
                //cout << endl;
            }
        }

        //best_coloring.insert(best_coloring.begin(), coloring[tr_index].begin(), coloring[tr_index].begin()+sz);

        // print in nice form
        cout << endl<<"Minimal coloring: ";
        for (int ind = 0; ind < best_coloring.size(); ind++) {
            cout << best_coloring[ind];
            if (ind % sz == sz-1) cout << " ";
        }
        cout << endl << endl;





        for (int hi=0; hi < history.size(); hi++) {
            s diff = result - history[hi];
            if (diff < SINF) {
                cout << " = " << f(i-hi-1) << " -> tr(A*M^" << (i-hi-1) << ") + [" << diff << "]" << endl;
                cout << "It holds: $\\tr(M^{n+" << (hi+1) << "}) = \\tr(M^n) + [" << diff << "] \\text{ for } n >= "<< (i-hi-1) <<"$."<<endl;
                return;
            }
        }
        //cout << endl;
        history.push_front(result);
    }
    cout << "No repetition found." << endl;
}