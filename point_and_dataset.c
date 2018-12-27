#include "point_and_dataset.h"
#include <stdio.h>
#include <stdlib.h>

Dataset
dataset_new (size_t feature_count, size_t point_count)
{
  size_t size = (feature_count + 1) * point_count;
  Dataset dataset;
  dataset.size = size;
  dataset.feature_count = feature_count;
  dataset.data = array_new (size);
  return dataset;
}

Point
point_new (size_t feature_count)
{
  Point point;
  point.feature_count = feature_count;
  point.data = array_new (feature_count + 1);
  return point;
}

void
point_free (Point point)
{
  array_free (point.data);
}

void
dataset_free (Dataset dataset)
{
  array_free (dataset.data);
}

void
dataset_fill_random (Dataset dataset)
{
  array_fill_random (dataset.data);
  for (int i = 0; i < point_count (dataset.feature_count, dataset.size); i++)
    {
      dataset.data.data[point_offset (dataset.feature_count, i)] = i;
    }
  return dataset;
}

Point
get_point_from_dataset (Dataset *dataset, size_t index)
{
  Point point;
  point.feature_count = dataset->feature_count;
  point.data = array_new (point.feature_count + 1);
  point.data.data[0]
    = dataset->data.data[point_offset (dataset->feature_count,
                                       index)]; // get point's index
  for (size_t i = 0; i < dataset->feature_count; i++)
    {
      point.data.data[i + 1]
        = dataset->data.data[feature_offset (dataset->feature_count, index, i)];
    }
  return point;
}

void
enter_point_to_dataset (Dataset *dataset, Point point, size_t index)
{
  assert (dataset->feature_count == point.feature_count);
  dataset->data.data[point_offset (dataset->feature_count, index)]
    = point.data.data[0];
  for (size_t i = 0; i < dataset->feature_count; i++)
    {
      dataset->data.data[feature_offset (dataset->feature_count, index, i)]
        = point.data.data[i + 1];
    }
}

void
print_dataset (Dataset *dataset)
{
  for (size_t i = 0; i < point_count (dataset->feature_count, dataset->size);
       i++)
    {
      printf ("Point with index %d: (",
              (int)
                dataset->data.data[point_offset (dataset->feature_count, i)]);
      for (size_t j = 0; j < dataset->feature_count; j++)
        {
          printf (
            " %f ",
            dataset->data.data[feature_offset (dataset->feature_count, i, j)]);
        }
      printf (")\n");
    }
}

void
print_point (Point point)
{
  printf ("Point with index %d: (", (int) point.data.data[0]);
  for (int i = 1; i <= point.feature_count; i++)
    {
      printf (" %f ", point.data.data[i]);
    }
  printf (")\n");
}