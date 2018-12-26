#include "vp_tree_local.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void
vp_tree_free (VPTree vp_tree)
{
  array_free (vp_tree.vp_heap);
  array_free (vp_tree.median_heap);
}

VPTree
vp_tree_from_dataset (Dataset dataset)
{
  VPTree tree;
  /* @todo */
  return tree;
}

number_t
dataset_distance (size_t index1, size_t index2)
{
}

#if 0
/**
 * Global variable holding the vantage point to be used by the
 * qsort comparison callback, compare_vp_distances().
 */
number_t g_vp;

int
compare_vp_distances (const void *a, const void *b)
{
  number_t fa = *(const number_t *) a;
  number_t fb = *(const number_t *) b;
  number_t da = fabs (fa - g_vp);
  number_t db = fabs (fb - g_vp);
  return (da > db) - (da < db);
}

void
vp_tree_local_split (ArraySlice dataset, Array *vp_heap, Array *median_heap,
                     size_t heap_root)
{
  if (dataset.size <= 1)
    return;

  assert (dataset.size % 2 == 0);

  assert (heap_root < vp_heap->size);
  assert (heap_root < median_heap->size);

  /* Pick a vantage point. */
  number_t vp = dataset.data[rand () % dataset.size];
  vp_heap->data[heap_root] = vp;

  /* Sort by distance from the vp. */
  g_vp = vp;
  qsort (dataset.data, dataset.size, sizeof (number_t), compare_vp_distances);

  /* Find the median. */
  size_t middle = dataset.size / 2 - 1;
  number_t median = 0.5f
                    * (fabs (dataset.data[middle] - vp)
                       + fabs (dataset.data[middle + 1] - vp));
  median_heap->data[heap_root] = median;

  ArraySlice left = {dataset.data, dataset.size / 2};
  ArraySlice right = {dataset.data + dataset.size / 2, dataset.size / 2};
  assert (left.size + right.size == dataset.size);

  vp_tree_local_split (left, vp_heap, median_heap,
                       heap_left_child_of (heap_root));
  vp_tree_local_split (right, vp_heap, median_heap,
                       heap_right_child_of (heap_root));
}

void
verify_tree (ArraySlice dataset, Array *vp_heap, Array *median_heap,
             size_t heap_root)
{
  if (dataset.size <= 1)
    return;

  assert (heap_root < vp_heap->size);
  assert (heap_root < median_heap->size);

  /* Get the vantage point. */
  number_t vp = vp_heap->data[heap_root];
  /* Get the median. */
  number_t median = median_heap->data[heap_root];

  size_t middle = dataset.size / 2;
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t distance = fabs (dataset.data[i] - vp);
      if (i < middle)
        {
          assert (distance <= median);
        }
      else
        {
          assert (distance >= median);
        }
    }

  /* Verify the children. */
  ArraySlice left = {dataset.data, dataset.size / 2};
  ArraySlice right = {dataset.data + dataset.size / 2, dataset.size / 2};
  assert (left.size + right.size == dataset.size);

  verify_tree (left, vp_heap, median_heap, heap_left_child_of (heap_root));
  verify_tree (right, vp_heap, median_heap, heap_right_child_of (heap_root));
}

#endif

/**
 * Verifies the correctness of a VP tree.
 * @param tree The tree to verify.
 * @param offset The offset from the start of the dataset to start at.
 * @param size The size of the dataset (as in the number of points).
 * @param heap_root The location of the vp and the median for this particular
 *                  level within the heaps.
 */
static void
verify_tree (VPTree tree, size_t offset, size_t size, size_t heap_root)
{
  if (size <= 1)
    return;

  size_t fc = tree.dataset.feature_count;

  /* Find the VP and the median. */
  number_t vp[tree.dataset.feature_count + 1];
  memcpy (vp, tree.vp_heap.data + point_offset (fc, heap_root),
          4 * sizeof (number_t));

  number_t median = tree.median_heap.data[heap_root];

  /* For each point, confirm the vp tree properties. */
  size_t middle = size / 2;

  for (size_t i = 0; i < size; i++)
    {
    }

  /* Verify the two children. */
}

void
test_vp_tree_local ()
{
  Array data = array_new (256 * 3);
  array_fill_random (data);

  Dataset dataset;
  dataset.data = data;
  dataset.feature_count = 2;
  dataset.size = 256;

  VPTree tree = vp_tree_from_dataset (dataset);
  verify_tree (tree, 0, dataset.size, 0);

  vp_tree_free (tree);
  array_free (data);
}