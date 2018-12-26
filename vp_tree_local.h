#pragma once

#include "array.h"
#include <assert.h>

/**
 * The points are stored in an one-dimensional array, tightly packed
 * with the following format:
 * [id, feature_1, feature_2, feature_3, ...]
 * Therefore, to access the second feature of point 3, the index is
 * (number of features + 1) * 3 + (1 + 2).
 * This system makes copying around data points easier since only
 * continuous chunks of memory need to be moved.
 */

/** Indexing helper functions. */
#define point_offset(feature_count, index) ((index) * (feature_count + 1))
#define feature_offset(feature_count, index, feature)                          \
  ((index) * (feature_count + 1) + (feature) + 1)
#define point_count(feature_count, dataset_size)                               \
  ((dataset_size) / (feature_count + 1))

typedef struct
{
  Array data;
  size_t size;
  size_t feature_count;
} Dataset;

typedef struct
{
  Dataset dataset;
  Array vp_heap;
  Array median_heap;
} VPTree;

VPTree vp_tree_from_dataset (Dataset dataset);
void vp_tree_free (VPTree vp_tree);

number_t dataset_distance (size_t index1, size_t index2);

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