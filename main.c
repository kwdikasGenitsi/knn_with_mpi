#include "median.h"
#include "mpi_partition.h"
#include "stack.h"
#include "vp_tree_distributed.h"
#include "vp_tree_local.h"
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

static long
delta_in_millis (struct timeval *start, struct timeval *end)
{
  long start_ms = start->tv_sec * 1000 + start->tv_usec / 1000;
  long end_ms = end->tv_sec * 1000 + end->tv_usec / 1000;
  return end_ms - start_ms;
}

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
          test_log2i ();
          test_stack ();
          test_dataset ();
          test_vp_tree_local ();
          test_vp_queue ();
        }
      MPI_Barrier (MPI_COMM_WORLD);
      test_mpi_partition_by_value ();
      MPI_Barrier (MPI_COMM_WORLD);
      test_vp_tree_distributed ();

      /* For some reason this gets corrupted. */
      MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
      if (world_rank == 0)
        printf ("All tests sucessful.\n");

      MPI_Finalize ();
      return EXIT_SUCCESS;
    }

  /* Parse the arguments. */
  size_t n = (size_t) (1 << atoi (argv[1]));
  size_t k = (size_t) (1 << atoi (argv[2]));

  size_t n_per_process = n / world_size;

  Dataset data = dataset_new (2, n_per_process);
  dataset_fill_random (data);

  struct timeval start, end;
  MPI_Barrier (MPI_COMM_WORLD);
  gettimeofday (&start, NULL);
  VPTreeDistributed tree = vp_tree_dist_from_dataset (data);
  MPI_Barrier (MPI_COMM_WORLD);
  gettimeofday (&end, NULL);

  long construction_time_ms = delta_in_millis (&start, &end);
  if (world_rank == 0)
    printf ("Construction time: %ld ms\n", construction_time_ms);

  vp_tree_dist_free (tree);

  MPI_Barrier (MPI_COMM_WORLD);
  gettimeofday (&start, NULL);
  for (int i = 0; i < world_size; i++)
    {
      for (size_t j = 0; j < data.size; j++)
        {
          number_t *target = dataset_point (data, j);
          Dataset nearest = vp_tree_dist_find_knn (tree, target, k);
          dataset_free (nearest);
        }
    }
  MPI_Barrier (MPI_COMM_WORLD);
  gettimeofday (&end, NULL);

  long all_knn_ms = delta_in_millis (&start, &end);
  if (world_rank == 0)
    printf ("All KNN time: %ld s\n", all_knn_ms);

  dataset_free (data);
  MPI_Finalize ();

  return EXIT_SUCCESS;
}
