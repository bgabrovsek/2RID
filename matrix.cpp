//
// Created by Bostjan on 23/03/2020.
//

#include "matrix.h"
#include <string>
#include <iostream>
#include <algorithm>
#include "defs.h"

using namespace std;

#define PATH_MULTIPLY(a, b) ( (((a) == INF) or ((b) == INF)) ? (INF) : ( (a) + (b) ) )
#define PERCENT(a, n) (0.1*uround(1000.0*(a) / (n)))


/**
 * @param n matrix dimension
 * @return identity matrix
 */
matrix * identity(u n) {
    matrix * m = new matrix();
    for (int i=0;i<n;i++) {
        vu row = vu(n,INF);
        row[i] = 0;
        m->push_row(row);
    }
    return m;
}


/**
 * pushes a row vector at the end of the matrix
 * @param v row vector
 */
void matrix::push_row(vector<u> v) {

    if (this->n == 0) { this->n = v.size(); this->data.reserve(this->n);} // first row added
    if (this->n != v.size()) {exit(1);}
    auto itr = find(this->data.begin(), this->data.end(), v);
    if (itr == this->data.end()) {
        // new row that does not exist in the matrix data
        this->dict.push_back(this->data.size()); // store the index of the new row
        this->data.push_back(v); // store the new row in data
    } else {
        // new row already exists in the matrix data
        this->dict.push_back(distance(this->data.begin(), itr)); // index of the row
    }
}

/**
 * append rows of m to the matrix
 * @param m
 */
void matrix::append(matrix *m) {

    //   for (auto row : *m)
    //       this->push_row(*row);

    for (int i = 0; i < m->size(); i++)
        this->push_row(m->data[m->dict[i]]);

}

int matrix::size() {
    return this->dict.size();
}

int matrix::data_size() {
    return this->data.size();
}

/**
 * @return number of non-infinity values
 */
int matrix::counta() {
    int result = 0;
    for (u i=0; i < this->size(); i++)
        for (const auto &elt : this->data[dict[i]])
            if (elt != INF) result++;
    return result;
}

/**
 * @return minimal non-infinity value
 */
int matrix::min() {
    int result = INF;
    for (u i=0; i < this->size(); i++)
        for (const auto &elt : this->data[dict[i]])
            if (elt < result) result = elt;
    return result;
}

/**
 * @return maximal non-infinity value
 */
int matrix::max() {
    int result = 0;
    for (u i=0; i < this->size(); i++)
        for (const auto &elt : this->data[dict[i]])
            if ((elt != INF)  and (elt > result)) result = elt;
    return result;
}


// MATRIX/MATH OPERATORS

/**
 * @return tranposed matrix
 */
// TODO: make faster w/ C++11
matrix matrix::transpose() const{
    matrix m;
    for (int i=0; i<this->n; i++) {
        vector<u> row;
        for (int j=0; j<this->n; j++)
            row.push_back(this->data[this->dict[j]][i]);
        m.push_row(row);
    }
    return m;
}


// OPERATORS


/**
 * Product of two matrices, note that the second matrix must be transposed
 * @param tm tranposed matrix m
 * @return product this (×) m (m is transposed tm), where (×) is path-algebra multiplication
 */
matrix matrix::operator * (const matrix & mt) {
    vector<u> * row, col;
    // matrix tm = m.transpose(); // TODO: does this copy the tranposed matrix?
    //  cout << "Multiplying" << this->dict.size()<< "and" << tm.dict.size() <<endl;
    matrix result;
    for (auto  row_ind : this->dict) { // loop over first matrix rows
        //  cout << row_ind << endl;
        vector<u> new_row;
        new_row.reserve(this->n);
        for (auto col_ind : mt.dict) { // loop over second matrix columns
            new_row.push_back( dot(this->data[row_ind], mt.data[col_ind]));
        }

        result.push_row(new_row);

    }
    // cout << "result: " << result << endl;
    return result;
}




/**
 *
 * @param tm tranposed matrix m
 * @return product this (×) m (m is transposed tm), where (×) is path-algebra multiplication
 */
s matrix::operator - (const matrix & m) {
    s diff = SINF;
    if (this->dict == m.dict) {

        for (auto row_index : this->dict) {

            for (int i = 0; i < this->n; i++) {
                u val_a  = this->data[row_index][i];
                u val_b = m.data[row_index][i];
                if ((val_a == INF) xor (val_b == INF)) return SINF;
                if ((val_a != INF) and (val_b != INF)) {

                    if (diff == SINF)
                        diff = (s)val_a - (s)val_b ;
                    else
                    if (diff != (s)val_a - (s)val_b) return SINF;
                }
            }
        }
    }
    return diff;
}

/**
 * @return min-plus trace of matrix
 */
u matrix::trace(bool return_index) {
    u min_elt = INF;
    u min_index = 0;
    for (int i=0; i<this->n; i++)
        if (this->data[this->dict[i]][i] < min_elt) {
            min_elt = this->data[this->dict[i]][i];
            min_index = i;
        }

    if (return_index)
        return min_index;
    else
        return min_elt;// + min_index*100;
}

// LOAD/SAVE


vector<uc> matrix::to_bytes() {
    vector<uc> v;
    v.reserve(this->n*this->n);
    for (int i = 0; i < this->n; i++) {

        for (const auto & val : this->data[this->dict[i]]) {

            if ((val != INF) and (val >= 0x7F)) exit(1);
            v.push_back( (val == INF) ? (UCINF) : ((uc)val));
        }
    }
    return v;
}


void matrix::from_bytes(const vector<uc> & v) {
    int n = usqrt(v.size()); // matrix size
    //TODO: empty data?
    for (int row_index=0; row_index < n; row_index++) {
        vector<u> row;
        row.reserve(v.size());
        for (int i  = 0; i < n; i++) {
            row.push_back((v[row_index * n + i] == UCINF) ? INF : (u)v[row_index * n + i]);
        }
        this->push_row(row);
    }
}

void matrix::save(const string &filename) {
    cout << "Saving " << filename <<endl;
    auto vec = this->to_bytes();
    ofstream f(filename, ios::out | ios::binary);
    f.write(reinterpret_cast<const char*>(&vec[0]), vec.size() * sizeof(uc));
    f.close();
}


/**
 * loads binary matrix
 * @param filename
 * @return true if loaded, false otherwise
 */
bool matrix::load(const string &filename) {
    return false; // leave line uncommented to compute each time
    if (not file_existsQ(filename)) return false;
    cout << "Loading " << filename <<endl;
    ifstream f(filename, ios::in | ios::binary);
    const size_t count = filesize(filename) / sizeof(uc);
    vector<uc> vec(count);
    f.read(reinterpret_cast<char*>(&vec[0]), count*sizeof(uc));
    //  cout << "read." << endl;
    f.close();
    this->from_bytes(vec);
    return true;

}

/**
 * path-algebra scalar product ab
 * @param a first vector
 * @param b second vector
 * @return scalar
 */
u dot(const vector<u> &a, const vector<u> &b) {
    u value, min_value = INF;
    auto ita = a.begin();
    auto itb = b.begin();
    auto ena = a.end();
    while (ita != ena) {
//        cout << ita << endl;
        value = PATH_MULTIPLY(*ita, *itb);
        //value = 0;
        if (value < min_value) min_value = value;
        ++ita;
        ++itb;
    }
    return min_value;
}


ostream& operator<<(ostream& os, matrix& m) {
    os << "Matrix " << m.size() << "x" << m.size()
       << " Values: " << m.min() << "-" << m.max() << " "
       << PERCENT(m.counta(), m.size()*m.size()) << "% full "
       << PERCENT(m.data_size(),m.size()) << "% compression"
       << " (trace " << m.trace() << ").";
    //return s.str();

    return os;
}

ostream& operator<<(ostream& os, matrix * m) {
    for (auto it : *m) {
        os << it << endl;
    }
    return os;
}