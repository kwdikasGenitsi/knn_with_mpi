#include "vp_tree_local.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

void
test_vp_tree_local ()
{
  Array *dataset = array_new (256);
  array_fill_random (dataset);

  /*
   * Our tree will have log2(x) levels and thus
   * have 2^log2(x) - 1 = x - 1 elements, where
   * x the size of the dataset.
   */
  size_t heap_size = dataset->size - 1;
  Array *vp_heap = array_new (heap_size);
  Array *median_heap = array_new (heap_size);

  vp_tree_local_split (array_get_slice (dataset, 0, dataset->size), vp_heap,
                       median_heap, 0);
  verify_tree (array_get_slice (dataset, 0, dataset->size), vp_heap,
               median_heap, 0);

  array_free (dataset);
  array_free (vp_heap);
  array_free (median_heap);
}