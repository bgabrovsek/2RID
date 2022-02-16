//
// Created by Bostjan on 21/03/2020.
//

#ifndef RAINBOW_DEFS_H
#define RAINBOW_DEFS_H

#include <iostream>
#include <functional>
#include <vector>
#include <string>
//#include <sstring>
#include <math.h>
#include <chrono>
#include <sys/stat.h>


//#include <sys/stat.h>
//#include <unistd.h>
//#include <string>
#include <fstream>

using namespace std;
using namespace std::chrono;

typedef unsigned short int u; // unsigned 16-bit, weights up to 65536
typedef signed short int s; // signed 16 bit
typedef unsigned char uc;
typedef unsigned long long int u64; // 64-bit

//vector shortcuts
typedef vector<u> vu;
typedef vector<vector<u>> vvu;
typedef vector<vector<vector<u>>> vvvu;

// functions
typedef const function<u(u)> &fuu;
//typedef const function<u(u,u)> &fuuu;
typedef const function<string(u)> &fsu;

// bit operations
//#define BIT(n) (0x1 << (n))
#define BIT(d, n) (((d) >> (n)) & 1u) // get n-th bit, BIT(0x01001,3) = 1
#define ONES(n) (((1u) << (n)) - 1u) // first n-bits to 1, BITS(3) = 0111b
#define SET_BIT(d, n) ((d) | (1u << (n))) // SET n+th bit of d to 1
#define CLEAR_BIT(d, n)  ((d) & ~(1u << (n))) // SET n+th bit of d to 0
#define BIT_RANGE(a, b) (((1u << ((b)-(a))) - 1u) << (a)) // set bits form a to b-1 to 1, BIT_RANGE(2,4) = 00001100b


#define INF 0xFFFFu
#define INF64 0xFFFFFFFFFFFFFFFFu
#define UCINF 0xFFu // unsigned byte
#define SINF 0x7FFF // signed infinity
// LISTS

u MAX(const vector<u> &v);
vector<u> concatenate(const vector<vector<u>> & V);
vector<vector<u>> concatenate(const vector<u> & A, const vector<vector<u>> & B, const vector<u> & C);
vector<u> trim(const vector<u> & V, u a, u b);
vector<vector<u>> trim(const vector<vector<u>> & V, u a, u b);
bool InQ(const vector<vector<u>> & v, const vector<u> & x);

// BINARY

u POPCOUNT(u bin); // returns number of bits
string BINARY(u64 b);

// FILE/STRING

string filename_data(string s0, u i0, string s1, u i1, string s2);
vector<string> filename_datas(string s0, u i0, string s1, u i1, string s2, u range2);
bool file_existsQ (const std::string& filename);

long filesize(std::string filename);
//int fstat64(int fildes, struct stat64 *buf);

// MATH

u64 usqrt(u64 x);
int uround(double x);

// TIMER
void timer_start(int i = 0);
double timer_value(int i = 0);

// PRINTING
ostream& operator<<(ostream& os, const vector<u>& v);

// MULTITHREADING

vector<vector<vector<u>>> split_to_bunch(const vector<vector<u>> &v, int N);

#endif //RAINBOW_DEFS_H