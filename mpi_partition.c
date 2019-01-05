#include "mpi_partition.h"
#include "stack.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

void
mpi_partition_by_value (Dataset dataset, Array values, number_t threshold,
                        MPI_Comm comm)
{
  int rank, comm_size;
  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &comm_size);
  assert (comm_size > 1);

  Array master_values = array_new (values.size * comm_size);
  Dataset master_dataset
    = dataset_new (dataset.feature_count, dataset.size * comm_size);

  MPI_Gather (values.data, values.size, MPI_NUMBER_T, master_values.data,
              values.size, MPI_NUMBER_T, 0, comm);
  MPI_Gather (dataset.data.data, dataset.data.size, MPI_NUMBER_T,
              master_dataset.data.data, dataset.data.size, MPI_NUMBER_T, 0,
              comm);

  if (rank == 0)
    {
      /* Standard Hoare's partition. */
      size_t low = 0;
      size_t high = master_dataset.size - 1;
      while (true)
        {
          while (master_values.data[low] < threshold)
            low++;
          while (master_values.data[high] > threshold)
            high--;

          if (low >= high)
            break;

          /* Swap the values .*/
          number_t temp = master_values.data[high];
          master_values.data[high] = master_values.data[low];
          master_values.data[low] = temp;

          /* Swap the points. */
          number_t *a = dataset_point (master_dataset, low);
          number_t *b = dataset_point (master_dataset, high);
          for (size_t i = 0; i <= dataset.feature_count; i++)
            {
              number_t temp = b[i];
              b[i] = a[i];
              a[i] = temp;
            }
        }
    }

  MPI_Barrier (comm);

  MPI_Scatter (master_dataset.data.data, dataset.data.size, MPI_NUMBER_T,
               dataset.data.data, dataset.data.size, MPI_NUMBER_T, 0, comm);
  MPI_Scatter (master_values.data, values.size, MPI_NUMBER_T, values.data,
               values.size, MPI_NUMBER_T, 0, comm);

  array_free (master_values);
  dataset_free (master_dataset);
}

float get_median_distance (Array distances, MPI_Comm comm, int group_rank,
                           int group_size);

void write_distances_from_vp (Array dest, Dataset dataset, number_t *vp);

void validation (float median, int partLength, int size, float *numberPart,
                 int world_rank, MPI_Comm comm);

void
test_mpi_partition_by_value ()
{
  int rank, comm_size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &comm_size);

  Dataset dataset = dataset_new (1, 64);

  Array values = array_new (dataset.size);

  array_fill_random (values);

  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      point[0] = values.data[i];
      point[1] = values.data[i];
    }

  Array values_to_shuffle = array_copy (values);
  number_t median
    = get_median_distance (values_to_shuffle, MPI_COMM_WORLD, rank, comm_size);
  array_free (values_to_shuffle);

  MPI_Barrier (MPI_COMM_WORLD);
  mpi_partition_by_value (dataset, values, median, MPI_COMM_WORLD);

  MPI_Barrier (MPI_COMM_WORLD);

  bool is_left = rank < (comm_size / 2) ? true : false;

  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      if (is_left)
        {
          assert (point[1] <= median);
        }
      else
        {
          assert (point[1] >= median);
        }
    }

  dataset_free (dataset);
}
