#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

// Subset structure
typedef struct Subset {
  int parent;
  int rank;
} Subset;

// Find the root of the set
int find(Subset* sets, int i);

// Union two sets
void union_sets(Subset* sets, int root_u, int root_v);

#endif