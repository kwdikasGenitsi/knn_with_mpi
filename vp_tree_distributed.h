#pragma once
#include "dataset.h"
#include "vp_tree_local.h"

typedef struct
{
  Dataset dataset;
  Array vps;
  Array medians;
  VPTree local_tree;
} VPTreeDistributed;

VPTreeDistributed vp_tree_dist_from_dataset (Dataset dataset);
void vp_tree_dist_free (VPTreeDistributed tree);

void vp_tree_dist_write_vp (VPTreeDistributed tree, number_t *vp, size_t index);
void vp_tree_dist_read_vp (VPTreeDistributed tree, number_t *vp, size_t index);

/**
 * Performs a search in a distributed VP tree.
 * @param tree The tree to search.
 * @param target The target point (must also have an ID at the start).
 * @param k The number of nearest neighbors to find.
 */
Dataset vp_tree_dist_find_knn (VPTreeDistributed tree, number_t *target,
                               size_t k);

void test_vp_tree_distributed ();
