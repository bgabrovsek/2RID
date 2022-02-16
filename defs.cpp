//
// Created by Bostjan on 21/03/2020.
//

#include "defs.h"
#include <algorithm>

#include <sstream>

#if _WIN32 || _WIN64
#define stat64 _stat64
#endif


high_resolution_clock::time_point value_start_timer[4] = {high_resolution_clock::now(), high_resolution_clock::now(), high_resolution_clock::now(), high_resolution_clock::now()};

//u INF = (std::numeric_limits<u>::max()); // infinity is biggest integer




// maximal element of list
u MAX(const vector<u> &v) {
    return *max_element(v.begin(), v.end());
}

// concatenate two lists
vector<u> concatenate(const vector<vector<u>> & V) {
    vector<u> vect;
    for (const auto & v : V)
        vect.insert(vect.end(), v.begin(), v.end());
    return vect;
}

// concatenate lists from B with A and C: A+B[0]+C, A+B[1]+C,...
vector<vector<u>> concatenate(const vector<u> & A, const vector<vector<u>> & B, const vector<u> & C) {
    vector<vector<u>> vects;
    for (const auto & b : B) {
        vector<u> v;
        v.insert(v.end(), A.begin(), A.end());
        v.insert(v.end(), b.begin(), b.end());
        v.insert(v.end(), C.begin(), C.end());
        vects.push_back(v);
    }
    return vects;
}

// returns left/right trimmed sublist
vector<u> trim(const vector<u> & V, u a, u b) {
    return vector<u>(V.begin()+a, V.end()-b);
}

// trim multiple lists
vector<vector<u>> trim(const vector<vector<u>> & V,u a,  u b) {
    vector<vector<u>> T;
    for (const auto & v : V)
        T.push_back( trim(v,a,b));
    return T;
}

// searches if vector x is inside vector v
bool InQ(const vector<vector<u>> & v, const vector<u> & x) {
    return (find(v.begin(), v.end(), x) != v.end());
}

// return number of bits in bin, works for 16 bit unsigned
u POPCOUNT(u x) {
    x = x - ((x >> 1u) & 0x5555u);
    x = (x & 0x3333u) + ((x >> 2u) & 0x3333u);
    return (((x + (x >> 4u)) & 0x0F0Fu) * 0x0101u) >> 8u;
}

ostream& operator<<(ostream& os, const vector<u>& v) {
    os << "[";
    for (const u &i : v) os << i << " ";
    if (v.size()) os <<"\b]";
    return os;
}


void timer_start(int i) {
    value_start_timer[i] = high_resolution_clock::now();
}

double timer_value(int i) {
    auto value_end_timer = high_resolution_clock::now();
    auto difference = std::chrono::duration_cast<std::chrono::milliseconds>(value_end_timer - value_start_timer[i]);
    return difference.count();
}

// unsigned long to binary string
string BINARY(u64 b) {
    stringstream s;
    bool one = false;
    for (u64 i=63;i<=63;i--) {
        if ((one) or (BIT(b, i))) {
            s << BIT(b, i);
            one = true;
        }
    }
    s << "b";
    return s.str();
}

// STRING
string filename_data(string s0, u i0, string s1, u i1, string s2) {
    stringstream s;
    s << "data/" << s0 << "-" << i0 << "-" << s1 << "-" << i1 << "-" << s2 << ".dat";
    return s.str();
}

vector<string> filename_datas(string s0, u i0, string s1, u i1, string s2, u range2) {
    vector<string> v;
    for (u m=0; m<range2; m++) {
        stringstream s;
        s << "data/" << s0 << "-" << i0 << "-" << s1 << "-" << i1 << "-" << s2 << "-" << m << ".dat";
        v.push_back(s.str());
    }
    return v;
}

// return bool if the file exists in system
bool file_existsQ (const std::string& filename) { ifstream f(filename); return bool(f);}
//TODO: change to stat64 for large files

long filesize(std::string filename) {
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

// unsigned square root
u64 usqrt(u64 x) { return (u64)(lround(sqrt(x)));}

// unsgined round
int uround(double x) { return (u)(lround(x)); }

// THREADING
/**
 * split vector into n partitions/bunches of approx. the same size
 * @tparam T any type
 * @param v the input vector
 * @param N number of bunches
 * @return a vector of vectors
 */
//template<class T>
vector<vector<vector<u>>> split_to_bunch(const vector<vector<u>> &v, int N) {
    int BUNCH_SIZE = v.size() / N;
    int MOD = v.size() % N;

    vector<vector<vector<u>>> bunch; //  a batch is a subset of the colorings, which is a vector of colorings (the coloring is also a vector)
    for (size_t i = 0; i < v.size(); i += BUNCH_SIZE + (MOD-- > 0 ? 1 : 0)) {
        bunch.emplace_back(v.begin() + i, v.begin() + i + BUNCH_SIZE + (MOD > 0 ? 1 : 0));
    }
    return bunch;
}

