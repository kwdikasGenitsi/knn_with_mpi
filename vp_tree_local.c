#include "vp_tree_local.h"
#include "sort_r.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
  number_t *vp;
  size_t feature_count;
} VPComparisonContext;

static int
compare_distances_from_vp (const void *a, const void *b, void *arg)
{
  VPComparisonContext *ctx = (VPComparisonContext *) arg;
  number_t *point_a = (number_t *) a;
  number_t *point_b = (number_t *) b;
  number_t distance_a = point_distance (ctx->vp, point_a, ctx->feature_count);
  number_t distance_b = point_distance (ctx->vp, point_b, ctx->feature_count);
  if (distance_a == distance_b)
    return 0;
  return (distance_a > distance_b) ? 1 : -1;
}

void
vp_tree_free (VPTree vp_tree)
{
  array_free (vp_tree.vp_heap);
  array_free (vp_tree.median_heap);
}

#ifdef PEDANTIC_CHECKS
static void
check_sorted (VPComparisonContext ctx, Dataset dataset, size_t offset,
              size_t size)
{
  assert (ctx.feature_count == dataset.feature_count);
  number_t prev_distance = 0.0f;
  for (size_t i = 0; i < size; i++)
    {
      number_t distance
        = point_distance (ctx.vp, dataset_point (dataset, offset + i),
                          ctx.feature_count);
      assert (prev_distance <= distance);
      prev_distance = distance;
    }
}
#endif

static void
vp_tree_split (VPTree tree, size_t offset, size_t size, size_t heap_root)
{
  if (size <= 1)
    return;

  assert (size % 2 == 0);

  assert (heap_root < tree.vp_heap.size);
  assert (heap_root < tree.median_heap.size);

  /* Pick a vantage point. */
  size_t vp_index = rand () % size;
  number_t vp[tree.dataset.feature_count + 1];
  dataset_read (tree.dataset, vp, vp_index, 1);
  vp_tree_write_vp (tree, vp, heap_root);

  /* Sort by distance from the vp. */
  VPComparisonContext ctx;
  ctx.feature_count = tree.dataset.feature_count;
  ctx.vp = vp;
  sort_r (dataset_point (tree.dataset, offset), size, 3 * sizeof (number_t),
          compare_distances_from_vp, &ctx);

#ifdef PEDANTIC_CHECKS
  check_sorted (ctx, tree.dataset, offset, size);
#endif

  /* Find the median. */
  size_t middle = size / 2 - 1;
  number_t *point1 = dataset_point (tree.dataset, offset + middle);
  number_t *point2 = dataset_point (tree.dataset, offset + middle + 1);
  number_t median
    = 0.5f
      * (point_distance (point1, vp, tree.dataset.feature_count)
         + point_distance (point2, vp, tree.dataset.feature_count));
  tree.median_heap.data[heap_root] = median;

  vp_tree_split (tree, offset, size / 2, heap_left_child_of (heap_root));
  vp_tree_split (tree, offset + size / 2, size / 2,
                 heap_right_child_of (heap_root));
}

VPTree
vp_tree_from_dataset (Dataset dataset)
{
  VPTree tree;
  size_t heap_size = dataset.size - 1;
  tree.dataset = dataset;
  tree.median_heap = array_new (heap_size);
  tree.vp_heap = array_new (heap_size * (dataset.feature_count + 1));

  vp_tree_split (tree, 0, dataset.size, 0);

  return tree;
}

void
vp_tree_read_vp (VPTree vp_tree, number_t *vp, size_t heap_index)
{
  size_t chunk_size = (vp_tree.dataset.feature_count + 1);
  size_t offset = heap_index * (vp_tree.dataset.feature_count + 1);
  memcpy (vp, vp_tree.vp_heap.data + offset, chunk_size * sizeof (number_t));
}

void
vp_tree_write_vp (VPTree vp_tree, number_t *vp, size_t heap_index)
{
  size_t chunk_size = (vp_tree.dataset.feature_count + 1);
  size_t offset = heap_index * (vp_tree.dataset.feature_count + 1);
  memcpy (vp_tree.vp_heap.data + offset, vp, chunk_size * sizeof (number_t));
}

/**
 * Verifies the correctness of a VP tree.
 * @param tree The tree to verify.
 * @param offset The offset from the start of the dataset to start at.
 * @param size The size of the dataset (as in the number of points).
 * @param heap_root The location of the vp and the median for this particular
 *                  level within the heaps.
 */
void
vp_tree_local_verify (VPTree tree, size_t offset, size_t size, size_t heap_root)
{
  if (size <= 1)
    return;

  assert (heap_root < tree.vp_heap.size);
  assert (heap_root < tree.median_heap.size);

  /* Get the vantage point. */
  number_t vp[tree.dataset.feature_count + 1];
  vp_tree_read_vp (tree, vp, heap_root);

  /* Get the median. */
  number_t median = tree.median_heap.data[heap_root];

  size_t middle = size / 2;
  for (size_t i = 0; i < size; i++)
    {
      number_t *point = dataset_point (tree.dataset, i + offset);
      number_t distance
        = point_distance (vp, point, tree.dataset.feature_count);
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
  vp_tree_local_verify (tree, offset, size / 2, heap_left_child_of (heap_root));
  vp_tree_local_verify (tree, offset + size / 2, size / 2,
                        heap_right_child_of (heap_root));
}

void
test_vp_tree_local ()
{
  Dataset dataset = dataset_new (2, 16);
  dataset_fill_random (dataset);

  VPTree tree = vp_tree_from_dataset (dataset);

  vp_tree_local_verify (tree, 0, dataset.size, 0);

  size_t k = 4;
  Dataset knn = dataset_new (dataset.feature_count, k);
  for (size_t i = 0; i < knn.size; i++)
    {
      number_t *point = dataset_point (knn, i);
      point[0] = __FLT_MAX__;
      point[1] = 1.0f;
      point[2] = 666.f;
    }

  number_t target[3] = {0.0f, 0.0f, 0.0f};

  print_dataset (&dataset);
  vp_tree_search (tree, knn, target, 0, dataset.size, 0);
  print_dataset (&knn);
  dataset_free (knn);

  dataset_free (dataset);
  vp_tree_free (tree);
}

void write_distances_from_vp (Array dest, Dataset dataset, number_t *vp);

static void
queue_insert (Dataset knn, number_t *point, number_t *target)
{
  Array distances = array_new (knn.size);
  write_distances_from_vp (distances, knn, target);

  number_t new_distance = point_distance (point, target, knn.feature_count);

  size_t insert_location = 0;

  while ((distances.data[insert_location] < new_distance)
         && (insert_location < distances.size))
    insert_location++;

  if (insert_location == distances.size)
    return;

  size_t stride = knn.feature_count + 1;

  assert (distances.size - 1 >= insert_location);
  for (size_t i = distances.size - 1; i > insert_location; i--)
    {
      number_t *current = dataset_point (knn, i - 1);
      number_t *next = dataset_point (knn, i);
      memcpy (next, current, stride * sizeof (number_t));
    }

  number_t *slot = dataset_point (knn, insert_location);
  memcpy (slot, point, stride * sizeof (number_t));
}

void
test_vp_queue ()
{
  Dataset dataset = dataset_new (1, 10);
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      point[0] = __FLT_MAX__;
    }

  number_t origin[2];
  origin[0] = 0.0f;
  origin[1] = 0.0f;
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t point[2];
      point[0] = 0.0f;
      point[1] = (number_t) ((50 - i) % 10);
      queue_insert (dataset, point, origin);
    }

  for (size_t i = 1; i < dataset.size; i++)
    {
      number_t *prev = dataset_point (dataset, i - 1);
      number_t *current = dataset_point (dataset, i);
      assert (prev[1] <= current[1]);
    }
}

void
vp_tree_search (VPTree tree, Dataset knn, number_t *target, size_t offset,
                size_t size, size_t heap_root)
{
  /* Get the vantage point. */
  number_t vp[tree.dataset.feature_count + 1];
  vp_tree_read_vp (tree, vp, heap_root);

  /* Get the median. */
  number_t median = tree.median_heap.data[heap_root];

  number_t distance = point_distance (vp, target, tree.dataset.feature_count);

  assert (size > 1);

  if (size == 2)
    {
      number_t *first = dataset_point (tree.dataset, offset);
      number_t *second = dataset_point (tree.dataset, offset + 1);
      queue_insert (knn, first, target);
      queue_insert (knn, second, target);
      return;
    }

  if (distance < median)
    {
      vp_tree_search (tree, knn, target, offset, size / 2,
                      heap_left_child_of (heap_root));
      number_t furthest_distance
        = point_distance (target, dataset_point (knn, knn.size - 1),
                          knn.feature_count);
      if (furthest_distance >= (median - distance))
        {
          /* Search the outer child. */
          vp_tree_search (tree, knn, target, offset + size / 2, size / 2,
                          heap_right_child_of (heap_root));
        }
    }
  else
    {
      vp_tree_search (tree, knn, target, offset + size / 2, size / 2,
                      heap_right_child_of (heap_root));
      number_t furthest_distance
        = point_distance (target, dataset_point (knn, knn.size - 1),
                          knn.feature_count);
      if (furthest_distance >= (distance - median))
        {
          /* Search the inner child. */
          vp_tree_search (tree, knn, target, offset, size / 2,
                          heap_left_child_of (heap_root));
        }
    }
}
