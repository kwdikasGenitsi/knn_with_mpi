#include "median.h"
#include "mpi_partition.h"
#include "stack.h"
#include "vp_tree_local.h"
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define UNUSED(expr)                                                           \
  do                                                                           \
    {                                                                          \
      (void) (expr);                                                           \
    }                                                                          \
  while (0)

int world_size;
int world_rank;

void
test_log2i ()
{
  int n = 1;
  for (int i = 0; i < 15; i++)
    {
      assert (log2i (n) == i);
      n *= 2;
    }
}

float masterPart (int world_size, int world_rank, int size, int partLength,
                  float *numberPart, MPI_Comm comm);

void slavePart (int world_rank, int partLength, float *numberPart, int size,
                MPI_Comm comm);

/**
 * Returns the number of processes in every working group, at step l.
 * For instance at the 2nd stage (l = 1), 4 processes would yield a
 * group size of 2.
 */
int
group_size (int l)
{
  return world_size / (1 << l);
}

/**
 * Should be the same for all processes within the same working group.
 * To be used with functions such as MPI_Comm_split().
 */
int
group_number (int l)
{
  return world_rank / group_size (l);
}

float
get_median_distance (Array distances, MPI_Comm comm, int group_rank,
                     int group_size)
{
  float median_distance;
  if (group_rank == 0)
    {
      median_distance
        = masterPart (group_size, group_rank, distances.size * group_size,
                      distances.size, distances.data, comm);
    }
  else // if i am a cheap slave
    {
      slavePart (group_rank, distances.size, distances.data,
                 distances.size * group_size, comm);
    }
  MPI_Bcast (&median_distance, 1, MPI_FLOAT, 0, comm);
  return median_distance;
}

int
master_rank (int l)
{
  return world_rank - (world_rank % group_size (l));
}

void
test_partitioning ()
{
  world_size = 16;
  world_rank = 12;
  assert (master_rank (0) == 0);
  assert (master_rank (1) == 8);
  assert (master_rank (2) == 12);
}

MPI_Comm
communicator_for_level (int l)
{
  MPI_Comm c;
  MPI_Comm_split (MPI_COMM_WORLD, group_number (l), 0, &c);
  return c;
}

void
write_distances_from_vp (Array dest, Dataset dataset, number_t *vp)
{
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t distance = point_distance (dataset_point (dataset, i), vp,
                                          dataset.feature_count);
      dest.data[i] = distance;
    }
}

int
less_than_median (Array *distances, number_t median_distance)
{
  int count_points
    = 0; // counts how many points are closer to vantage point than median.
  for (size_t i = 0; i < distances->size; i++)
    {
      if (distances->data[i] <= median_distance)
        {
          count_points++;
        }
    }
  return count_points;
}

#if 0
void
split_parallel (Array *dataset, Stack *vp_stack, Stack *median_stack, int l)
{
  MPI_Comm comm = communicator_for_level (l);
  int rank;
  int comm_size;
  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &comm_size);

  int chunk_size = (dataset->size) / (1 << l);

  /* Pick a vantage point and announce it. */
  int vp_index = rand () % chunk_size;
  number_t vp_value = dataset->data[vp_index];
  MPI_Bcast (&vp_value, 1, MPI_NUMBER_T, 0, comm);

  /* Find the distances from the VP. */
  Array *distances = array_distances_from_vp (dataset, vp_value);

  /* Find the median. */
  float median = 0.0f;
  if (rank == 0)
    {
      /** @todo Use masterPart() to find the median. */
    }
  else
    {
      /** @todo Invoke slavePart() here. */
    }

  /* Split the free from the commies. */
  free (distances);

  /* Push the vantage point and the median onto the stack. */
  stack_push (vp_stack, vp_value);
  stack_push (median_stack, median);

  MPI_Comm_free (&comm);
}
#endif

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);

  /* If the tests are invoked, have the root run them. */
  if (argc == 2 && (strcmp (argv[1], "test") == 0))
    {
      if (world_rank == 0)
        {
          if (getenv ("PAUSE"))
            {
              float f = 0.0f;
              while (f != 3.0f)
                scanf (" %f\n", &f);
            }
          printf ("Running the unit tests...\n");
          test_partitioning ();
          test_log2i ();
          test_stack ();
          test_dataset ();
          test_vp_tree_local ();
        }
      MPI_Barrier (MPI_COMM_WORLD);
      test_mpi_partition_by_value ();

      /* For some reason this gets corrupted. */
      MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
      if (world_rank == 0)
        printf ("All tests sucessful.\n");

      MPI_Finalize ();
      return EXIT_SUCCESS;
    }

  /*
  Array dataset = array_new (5);
  array_fill_random (dataset);
  int tree_depth = log2i (dataset.size * world_size);
  Stack *vp_stack = stack_new (tree_depth);
  Stack *median_stack = stack_new (tree_depth);

  for (int l = 0; group_size (l) > 1; l++)
    {
      split_parallel (dataset, vp_stack, median_stack, l);
    }

  array_free (dataset);
  */
  MPI_Finalize ();

  return EXIT_SUCCESS;
}
