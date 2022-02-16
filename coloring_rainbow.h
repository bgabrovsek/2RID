//
// Created by Bostjan on 19/03/2020.
//

#ifndef RAINBOW_COLORING_RAINBOW_H
#define RAINBOW_COLORING_RAINBOW_H

#include "graph.h"
#include "defs.h"
#include <functional>
#include <iostream>
#include <forward_list>
#include <array>
#include <vector>

u I2RD_conjecture(u n, u k, u r);
u IPC(u n, u r);

//vector<vector<u>> Monograph_Partial_Independent_Rainbow_Domination_Colorings(Graph* g, u r, fu f_deg);
vector<vector<u>> Partial_Independent_Rainbow_Domination_Colorings(Graph* g, const vector<u> &color_in, const vector<u> &color_out, u r, fuu f_deg);

bool PIRD_extendableQ(Graph *g, const vector<vector<u>> & colors_ex, const vector<u> & color, u r, fuu f_deg);
vector<vector<u>> PIRD_filter_extendable(Graph *g, const vector<vector<u>> & colors_ex, const vector<vector<u>> & colors, u r, fuu f_deg);

bool GoodPIRDQ(Graph *g, const vu & c, u r, fuu f_deg);

#endif //RAINBOW_COLORING_RAINBOW_H