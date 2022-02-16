//
// Created by Bostjan on 23/03/2020.
//

#ifndef RAINBOW_MATRIX_H
#define RAINBOW_MATRIX_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "defs.h"
using namespace std;

#include "defs.h"


/**
 * class for a square (?) matrix, optimized in such a way that a lot of rows repeat
 */
class matrix {
    vector<vector<u>> data; // rows
    vector<u> dict; // convert to data row, dictionary if row_dict[7] = 2, then matrix[7] = data[2];
    u n = 0;
public:
    //explicit matrix(u n);
    void push_row(vector<u>);
    void append(matrix *m);
    //void analyse();
    // iterator
    //int begin() { return 0; }
    //int end() { return this->n; }

    // size
    int size();
    int data_size();
    int counta();
    int min();
    int max();

    // converting
    vector<uc> to_bytes();
    void from_bytes(const vector<uc> & v);

    // matrix operators
    matrix transpose() const;
    u trace(bool return_index = false);

    vector<u> get_row(u index) {return this->data[this->dict[index]];}

    // operators
    //vector<u> operator [] (int row_index); // &!!
    matrix operator * (const matrix & m);
    s operator - (const matrix &m);
    //int operator - (const matrix & m);
    // save/load
    void save(const string & filename);
    bool load(const string & filename);

    // iterator, copies the rows, do not use for fast computation

    class the_iterator {
    private:
        u row_index_;
        matrix *parent_;
    public:
        explicit the_iterator(u row_index, matrix *parent) : row_index_(row_index), parent_(parent) {}
        vector<u>  operator*() const { return parent_->data[parent_->dict[row_index_]]; } // returns pointer to vector
        bool operator!=(const the_iterator& other) const { return row_index_ != other.row_index_; }
        void operator++() { ++row_index_; }
    };

    the_iterator begin() {return the_iterator(0, this); }
    the_iterator end() {return the_iterator(this->n, this); }

};

matrix * identity(u n);

u dot(const vector<u> &a, const vector<u> &b);

ostream& operator<<(ostream& os, matrix& m);
ostream& operator<<(ostream& os, matrix * m);

#endif //RAINBOW_MATRIX_H