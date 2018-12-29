#pragma once

#include "array.h"
#include "dataset.h"
#include <assert.h>

typedef struct
{
  Dataset dataset;
  Array vp_heap;
  Array median_heap;
} VPTree;

VPTree vp_tree_from_dataset (Dataset dataset);
void vp_tree_free (VPTree vp_tree);

void vp_tree_read_vp (VPTree vp_tree, number_t *vp, size_t heap_index);
void vp_tree_write_vp (VPTree vp_tree, number_t *vp, size_t heap_index);

/**
 * Generates an entire VP tree, recursively.
 * The heaps must have size at least dataset.size - 1.
 * @param dataset The dataset to split.
 * @param vp_heap A binary heap with the vantage points.
 * @param median_heap A binary heap with the median distances.
 * @param heap_root The index of within the heap where the next
 *                  vantage point and median should be written.
 *                  Usually, this should be 0.
 * @warning This function is NOT thread-safe.
 */
void vp_tree_local_split (ArraySlice dataset, Array *vp_heap,
                          Array *median_heap, size_t heap_root);

void test_vp_tree_local ();