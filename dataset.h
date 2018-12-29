#pragma once

#include "array.h"
#include <assert.h>

/** Incredibly inefficient data structure. */
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
  size_t size; /**< The number of points in the dataset. */
  size_t feature_count;
} Dataset;

Dataset dataset_new (size_t feature_count, size_t point_count);
void dataset_fill_random (Dataset dataset);
void dataset_free (Dataset dataset);

/**
 * Writes a point to the dataset, possibly overwriting an existing point.
 * @param dataset The dataset to write into.
 * @param data The index and the coordinates of the point.
 * @param index The offset (in points) from the beginning of the dataset.
 * @param count The number of points to write.
 */
void dataset_write (Dataset dataset, number_t *data, size_t index,
                    size_t count);

/**
 * Reads a point from the dataset, possibly overwriting the buffer given.
 * @param dataset The dataset to read from.
 * @param data The buffer to write into.
 * @param index The offset (in points) from the beginning of the dataset.
 * @param count The number of points to read.
 */
void dataset_read (Dataset dataset, number_t *data, size_t index, size_t count);

void dump_point (number_t *p, size_t feature_count);

number_t point_distance (number_t *p1, number_t *p2, size_t feature_count);

number_t *dataset_point (Dataset dataset, size_t index);

/**
 * @note Leaving these here for compatibility purposes, but they should not
 * be used since they are extremely low-performance given the 4-byte dynamic
 * allocations.
 */
Point point_new (size_t feature_count);
void point_free (Point point);
Point get_point_from_dataset (Dataset *dataset, size_t index);
void enter_point_to_dataset (Dataset *dataset, Point new_point, size_t index);

void print_dataset (Dataset *dataset);
void print_point (Point point);

void test_dataset ();