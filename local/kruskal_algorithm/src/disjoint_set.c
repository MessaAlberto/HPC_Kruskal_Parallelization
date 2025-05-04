#include "disjoint_set.h"

// Find the root of the set
int find(Subset* sets, int i) {
  if (sets[i].parent != i) {
    sets[i].parent = find(sets, sets[i].parent);
  }
  return sets[i].parent;
}

// Union two sets
void union_sets(Subset* sets, int root_u, int root_v) {
  if (sets[root_u].rank < sets[root_v].rank) {
    sets[root_u].parent = root_v;
  } else if (sets[root_u].rank > sets[root_v].rank) {
    sets[root_v].parent = root_u;
  } else {
    sets[root_v].parent = root_u;
    sets[root_u].rank++;
  }
}