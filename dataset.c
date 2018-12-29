#include "dataset.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
dataset_write (Dataset dataset, number_t *data, size_t index, size_t count)
{
  size_t chunk_size = (dataset.feature_count + 1) * count;
  size_t offset = index * (dataset.feature_count + 1);
  memcpy (dataset.data.data + offset, data, chunk_size * sizeof (number_t));
}

void
dataset_read (Dataset dataset, number_t *data, size_t index, size_t count)
{
  size_t chunk_size = (dataset.feature_count + 1) * count;
  size_t offset = index * (dataset.feature_count + 1);
  memcpy (data, dataset.data.data + offset, chunk_size * sizeof (number_t));
}

number_t
point_distance (number_t *p1, number_t *p2, size_t feature_count)
{
  /* Increase the two pointers by 1 to skip the index. */
  p1++;
  p2++;
  number_t square_sum = 0.0f;
  for (size_t i = 0; i < feature_count; i++)
    {
      number_t distance = p1[i] - p2[i];
      square_sum += distance * distance;
    }
  return sqrtf (square_sum); /**< Replace with sqrt() if number_t is double. */
}

void
dataset_fill_random (Dataset dataset)
{
  array_fill_random (dataset.data);
  for (int i = 0; i < point_count (dataset.feature_count, dataset.size); i++)
    {
      dataset.data.data[point_offset (dataset.feature_count, i)] = i;
    }
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

void
test_dataset ()
{
  Dataset d = dataset_new (2, 10);
  for (size_t i = 0; i < 10; i++)
    {
      for (size_t j = 0; j < 3; j++)
        {
          d.data.data[i * 3 + j] = j;
        }
    }

  number_t buffer[] = {3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  dataset_write (d, buffer, 3, 2);

  for (size_t i = 9; i < 15; i++)
    {
      assert (d.data.data[i] == buffer[i - 9]);
    }

  dataset_read (d, buffer, 1, 2);

  for (size_t i = 0; i < 6; i++)
    {
      assert (buffer[i] == (number_t) (i % 3));
    }

  number_t point1[3];
  number_t point2[3];
  dataset_read (d, point1, 0, 1);
  dataset_read (d, point2, 3, 1);

  number_t distance = point_distance (point1, point2, d.feature_count);
  assert (fabs (distance - sqrtf (18)) <= __FLT_EPSILON__);

  dataset_free (d);
}
