#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>

typedef struct Edge {
  int u;
  int v;
  int weight;
} Edge;

Edge* read_graph(int* V, int* E, char* in_file);
int compare_edges(const void* a, const void* b);
void create_output_filename(char* input_filename, char* prefix, char** output_filename);

#endif