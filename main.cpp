#include <iostream>
#include <string>
#include <vector>
//#include <thread>
#include "defs.h"
#include "graph.h"
#include "coloring_rainbow.h"
#include "matrix.h"
#include "domination.h"

using namespace std;

// Check if compiled using 32-bit or 64-bit compiler
// Check windows
#if _WIN32 || _WIN64
#if _WIN64
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif
// Check GCC
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif

#define PETER 0
#define CnxPm 1
#define CnxCm 2

void main_rainbow(u K, u r, u select) {

    Graph *monograph, *duograph, *quatrograph, *pentagraph, *start_pentagraph[K];
    function<u(u)> vert_deg_f;

    vvvu ex_colorings_in, ex_colorings_out; // store coloring extensions (left and right)
    vvu ex_colorings, colorings; // list of coloring of the  monograph

    matrix * start_mat[K]; // starting matrice for PETERSEN, not needed for others

    // GENERATE MONOGRAPHS

    switch(select) {
        case PETER:
            cout << "----" << endl << "Petersen graph (M = "<< K << ")" << endl;
            monograph = PetersenMonograph(K,K); // one-side monograph
            duograph = PetersenMonograph(K,2*K); // double monograph
            quatrograph = PetersenMonograph(K, 4*K); // four monographs
            pentagraph = PetersenMonograph(K, 5*K); // five monographs
            for (int m = 0; m < K; m++)  // generate starting graphs
                start_pentagraph[m] = PetersenMonograph(K, 5*K+m);
            break;
            vert_deg_f = CUBIC; // vertex degree function

        case CnxPm:
            cout << "----" << endl << "Cn×Pm (M = "<< K << ")" << endl;
            monograph = PathPath(K, 1); // one-side monograph
            duograph = PathPath(K, 2); // double monograph
            quatrograph = PathPath(K, 4); // four monographs
            pentagraph = PathPath(K, 5); // five monographs
            vert_deg_f = [K](int a) { if ((a % K == 0) or (a % K == K-1)) return 3; else return 4; };
            break;

        case CnxCm:
            cout << "----" << endl << "Cn×Cm (M = "<< K << ")" << endl;
            monograph = CyclePath(K, 1); // one-side monograph
            duograph = CyclePath(K, 2); // double monograph
            quatrograph = CyclePath(K, 4); // four monographs
            pentagraph = CyclePath(K, 5); // five monographs
            vert_deg_f = QUADRATIC;
    }

    cout << "Monograph:" << *duograph << endl;
    //cout << "Monograph degrees: ";
    //for (int i = 0; i < N; i++) cout << vert_deg_f(i) << " ";
    //cout << endl;

    // compute colorings of (small) monograph
    // compute colorings of graphs (left/right) extensions
    ex_colorings = Partial_Independent_Rainbow_Domination_Colorings(monograph, vector<u>(), vector<u>(), r, vert_deg_f);

    // compute colorings of middle graph
    colorings = Middle_Colorings(quatrograph, ex_colorings, ex_colorings, r, vert_deg_f,
                                 Partial_Independent_Rainbow_Domination_Colorings, ex_colorings_in, ex_colorings_out);


    cout << "Number of colorings: " << colorings.size() << endl;

    // compute matrices
    matrix * mat =  MatrixEx(duograph, pentagraph,colorings, ex_colorings_in, ex_colorings_out, r, vert_deg_f, GoodPIRDQ);

    cout << "Monomatrix: " << *mat << endl;

    // compute start matrices

    switch (select) {
        case PETER:
            start_mat[0] = mat;
            for (int m = 1; m < K; m++)
                start_mat[m] = MatrixEx(duograph, start_pentagraph[m], colorings, ex_colorings_in, ex_colorings_out, r,
                                        CUBIC, GoodPIRDQ);
            break;

        case CnxPm:
        case CnxCm:
            for (int m = 0; m < K; m++)
                start_mat[0] = mat; // starting matrices all equal
    }

    // compute and analyse the matrices

    switch (select) {

        case PETER:
            for (u m = 0; m < K; m++) {
                u n0 = monograph->n/2 + m, nm = monograph->n/2;
                compute_matrix_powers(colorings, *start_mat[m], *mat,
                                      [n0, nm, m, K](int i) { return "P("+to_string(n0 + i*nm)+","+to_string(K)+")";},
                                      [n0, nm, m, K, r](int i) { return I2RD_conjecture(n0 + i*nm, K, r);}
                );
            }
            break;

        case CnxPm:
        case CnxCm:

            int N_START = 1;
            int N = K;

            for (u m = 0; m < N_START; m++) {
                u n0 = monograph->n/N + m,nm = monograph->n/N;
                n0 = 0;
                compute_matrix_powers(colorings, *mat,
                                      [n0, nm, m, N](int i) { return "$i_{ri2}(C_{"+to_string(n0 + i*nm)+"} \\times C_{"+to_string(N)+"})";},
                                      [n0, nm, m, N, r](int i) { return IPC(n0 + i*nm, r);}
                );
            }

    }










}



//int N_THREADS;

/**
 * Computes Independent r-rainbow domination numbers for Petersen graphs P(n,k)
 * @param k parameter of P(n,k)
 * @param r number of colors
 */
void main_rainbow_petersen(u K, u r) {

    // GENERATE MONOGRAPHS
    Graph *monograph = PetersenMonograph(K, K); // one-side monograph
    Graph *duograph = PetersenMonograph(K, 2*K); // double monograph
    Graph *quatrograph = PetersenMonograph(K, 4*K); // four monographs
    Graph *pentagraph = PetersenMonograph(K, 5*K); // five monographs
    Graph *start_pentagraph[K]; // starting monographs (mod)    break;
    // print the monograph
    cout << "Monograph:" << *duograph << endl;

    // generate starting graphs
    for (int m = 0; m < K; m++)
        start_pentagraph[m] = PetersenMonograph(K, 5*K+m);

    vvvu ex_colorings_in, ex_colorings_out; // store coloring extensions (left and right)
    vvu colorings; // list of coloring of the  monograph

    // compute colorings of (small) monograph
    vvu ex_colorings = Partial_Independent_Rainbow_Domination_Colorings(monograph, vector<u>(), vector<u>(), r, CUBIC);

    // compute colorings of middle graph
    colorings = Middle_Colorings(quatrograph, ex_colorings, ex_colorings, r, CUBIC, Partial_Independent_Rainbow_Domination_Colorings, ex_colorings_in, ex_colorings_out);
    cout << "Number of colorings: " << colorings.size() << endl;

    // compute matrices
    matrix * mat =  MatrixEx(duograph, pentagraph,colorings, ex_colorings_in, ex_colorings_out, r, CUBIC, GoodPIRDQ);
    cout << "Monomatrix: " << *mat << endl;

    // compoute starting matrices
    matrix * start_mat[K];

    start_mat[0] = mat;
    for (int m=1; m < K; m++) {
        start_mat[m] = MatrixEx(duograph, start_pentagraph[m], colorings, ex_colorings_in, ex_colorings_out, r, CUBIC, GoodPIRDQ);
        cout << "Matrix mod " << m << ": " << *start_mat[m] << endl;
    }

    // compute matrix powers
    for (u m = 0; m < K; m++) {
        u n0 = monograph->n/2 + m, nm = monograph->n/2;
        compute_matrix_powers(colorings, *start_mat[m], *mat,
                              [n0, nm, m, K](int i) { return "P("+to_string(n0 + i*nm)+","+to_string(K)+")";},
                              [n0, nm, m, K, r](int i) { return I2RD_conjecture(n0 + i*nm, K, r);}
        );
    }







}





/**
 * Computes Independent r-rainbow domination numbers for Petersen graphs P(n,k)
 * @param k parameter of P(n,k)
 * @param r number of colors
 */
void main_rainbow_path_cycle(u N, u r) {

    int N_START = 1;

/*    Graph *monograph = PathPath(N, 1); // one-side monograph
    Graph *duograph = PathPath(N, 2); // double monograph
    Graph *quatrograph = PathPath(N, 4); // four monographs
    Graph *pentagraph = PathPath(N, 5); // five monographs
*/

    Graph *monograph = CyclePath(N, 1); // one-side monograph
    Graph *duograph = CyclePath(N, 2); // double monograph
    Graph *quatrograph = CyclePath(N, 4); // four monographs
    Graph *pentagraph = CyclePath(N, 5); // five monographs


    cout << "Monograph:" << *monograph << endl;
    cout << "Duograph:" << *duograph << endl;
    cout << "Quatrograph:" << *quatrograph << endl;

    vvvu ex_colorings_in, ex_colorings_out; // store coloring extensions
    vvu colorings; // coloring of monograph


    //auto vert_deg_f = [N](int a) { return vert_degree_path(a, N);};
    auto vert_deg_f = [N](int a) { return 4;};
    //auto vdeg = [N](int a) { return 2;};
    //auto vdeg = [N](int a) { return 2;};

    cout << "Monograph degrees: ";
    for (int i = 0; i < N; i++) cout << vert_deg_f(i) << " ";
    cout << endl;

    // compute colorings of graphs (left/right) extensions
    vvu ex_colorings = Partial_Independent_Rainbow_Domination_Colorings(monograph, vector<u>(), vector<u>(), r, vert_deg_f);
    colorings = Middle_Colorings(quatrograph, ex_colorings, ex_colorings, r, vert_deg_f,
            Partial_Independent_Rainbow_Domination_Colorings, ex_colorings_in, ex_colorings_out);

    cout << "Number of colorings: " << colorings.size() << endl;

    // compute matrices
    matrix * mat =  MatrixEx(duograph, pentagraph,colorings, ex_colorings_in, ex_colorings_out, r, vert_deg_f, GoodPIRDQ);
    matrix * start_mat[N_START];
    start_mat[0] = mat;

    // compute matrix powers
    for (u m = 0; m < N_START; m++) {
        u n0 = monograph->n/N + m,nm = monograph->n/N;
        n0 = 0;
        compute_matrix_powers(colorings, *mat,
                              [n0, nm, m, N](int i) { return "$i_{ri2}(C_{"+to_string(n0 + i*nm)+"} \\times C_{"+to_string(N)+"})";},
                              [n0, nm, m, N, r](int i) { return IPC(n0 + i*nm, r);}
        );
    }
}



int main() {

#ifdef ENVIRONMENT64
    cout << "Using 64-bit architecture.\n";
#else
    cout << "Using 32-bit architecture.\n";
#endif

/*    N_THREADS = std::thread::hardware_concurrency()/2;

    cout << "Number of supported threads: " << std::thread::hardware_concurrency() << endl;
    cout << "Number of used threads:      " << N_THREADS << endl;
*/


    timer_start();
    cout << "Hello.";
    main_rainbow_path_cycle(3,2);

    // PETERSEN
/*    for (u M = 1; M <= 4; M++)
        main_rainbow(M, 2, PETER);
*/
    // CnxPm
//    for (u M = 1; M <= 4; M++)
  //      main_rainbow(M, 2, CnxPm);
/*
    // CnxPm
    for (u M = 1; M <= 4; M++)
        main_rainbow(M, 2, CnxCm);

*/
    cout << "Time: " <<timer_value() << endl;
    return 0;
}