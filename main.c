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
      // test_vp_tree_distributed ();

      /* For some reason this gets corrupted. */
      MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
      if (world_rank == 0)
        printf ("All tests sucessful.\n");

      MPI_Finalize ();
      return EXIT_SUCCESS;
    }

  test_vp_tree_distributed ();
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
