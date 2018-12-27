#pragma once

#include "array.h"
#include <assert.h>

typedef struct
{
  Array data;
  size_t feature_count;
} Point;

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

Point point_new (size_t feature_count);
Dataset dataset_new (size_t feature_count, size_t point_count);
void dataset_fill_random (Dataset dataset);

void point_free (Point point);
void dataset_free (Dataset dataset);

Point get_point_from_dataset (Dataset *dataset, size_t index);
void enter_point_to_dataset (Dataset *dataset, Point new_point, size_t index);

void print_dataset (Dataset *dataset);
void print_point (Point point);