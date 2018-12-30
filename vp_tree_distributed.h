#pragma once
#include "dataset.h"

typedef struct
{
  Dataset dataset;
  Array vps;
  Array medians;
} VPTreeDistributed;

VPTreeDistributed vp_tree_dist_from_dataset (Dataset dataset);
void vp_tree_dist_free (VPTreeDistributed tree);

void vp_tree_dist_write_vp (VPTreeDistributed tree, number_t *vp, size_t index);
void vp_tree_dist_read_vp (VPTreeDistributed tree, number_t *vp, size_t index);

void test_vp_tree_distributed ();
