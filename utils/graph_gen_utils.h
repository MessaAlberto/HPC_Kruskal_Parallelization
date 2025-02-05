#ifndef GRAPH_GEN_UTILS_H
#define GRAPH_GEN_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void open_file(FILE** fp, const char* filename, const char* mode);
void get_params(int argc, char* argv[], int* V, int* max_wgt);
long long random_num_gen(long long min, long long max);
void shuffle(int* nodes, int v);
double calculate_density(int V, int E);

#endif