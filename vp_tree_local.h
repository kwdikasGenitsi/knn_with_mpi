#pragma once

#include "array.h"
#include "dataset.h"
#include <assert.h>

/**
 * @brief Structure representing a VP tre.
 *
 * The vantage points and the medians are stored in two binary heaps. For
 * vp_heap, keep in mind that each point needs (feature count + 1) slots. You
 * should therefore use vp_tree_read_vp() and vp_tree_write_vp() in conjunction
 * with heap_parent_of(), heap_left_child_of() and heap_right_child_of().
 */
typedef struct
{
  Dataset dataset;
  Array vp_heap;
  Array median_heap;
} VPTree;

/**
 * Creates a VP tree from the given dataset. The items in the dataset are
 * rearranged.
 * @param dataset The dataset.
 * @return A new VP tree based on the dataset.
 */
VPTree vp_tree_from_dataset (Dataset dataset);

/**
 * Frees the heaps of a VP tree.
 * This does not deallocate the memory used by the dataset.
 */
void vp_tree_free (VPTree vp_tree);

/**
 * Reads a vantage point from the heap.
 * @param vp_tree The tree to read from.
 * @param vp The buffer to write the result into.
 * @param heap_index The index (in points) within the heap.
 */
void vp_tree_read_vp (VPTree vp_tree, number_t *vp, size_t heap_index);

/**
 * @sa vp_tree_read_vp
 */
void vp_tree_write_vp (VPTree vp_tree, number_t *vp, size_t heap_index);

/**
 * @brief Unit test.
 */
void test_vp_tree_local ();

/**
 * Verifies the correctness of a local VP tree.
 */
void vp_tree_local_verify (VPTree tree, size_t offset, size_t size,
                           size_t heap_root);

/**
 * Performs a search in a local VP tree.
 * @param tree The tree to search.
 * @param target The target point (must also have an ID at the start).
 * @param k The number of nearest neighbors to find.
 */
Dataset vp_tree_find_knn (VPTree tree, number_t *target, size_t k);

void queue_insert (Dataset knn, number_t *point, number_t *target);

void test_vp_queue ();