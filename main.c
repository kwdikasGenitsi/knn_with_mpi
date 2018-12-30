#include "median.h"
#include "stack.h"
#include "vp_master_buffer.h"
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

size_t g_feature_count;
number_t *dataset;

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

/**
 * Partitions the given dataset across a communicator, so that entries for which
 * values[i] < threshold are in the left processes and vice versa.
 */
void
partition_by_value (Dataset dataset, Array values, number_t threshold,
                    MPI_Comm comm)
{
  int rank, comm_size;
  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &comm_size);
  assert (comm_size > 1);

  size_t stride
    = dataset.feature_count + 1; /* Useful for navigating the stacks. */

  bool is_left = rank < (comm_size / 2) ? true : false;

  /* Create a buffer and fill it with the values that need to be sent. */
  Stack *stack = stack_new (dataset.data.size);

  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t value = values.data[i];
      if ((is_left && (value > threshold)) || (!is_left && (value < threshold)))
        {
          number_t *point = dataset_point (dataset, i);
          stack_push_buffer (stack, point, stride);
          /* Mark the empty slots by setting their first coordinate to
           * the threshold. */
          point[1] = threshold;
        }
    }

  size_t elements_to_swap = stack_get_size (stack);
  /* Stores the number of items each process wants to swap.
   * Only valid for the master.
   */
  size_t elements[comm_size];
  MPI_Gather (&elements_to_swap, 1, MPI_UINT64_T, elements, 1, MPI_UINT64_T, 0,
              comm);

  Array receive_buffer = array_new (elements_to_swap);
  if (rank == 0)
    {
      Stack *buffer_stack = stack_new (dataset.data.size * 2);
      int middle_rank = comm_size / 2;
      int sending_rank = middle_rank;

      for (int i = 0; i < middle_rank; i++)
        {
          /* How many elements does this process need? */
          size_t requested_elements = elements[i];

          while (stack_get_size (buffer_stack) < requested_elements)
            {
              if (sending_rank >= comm_size)
                fprintf (stderr, "Imbalanced request going left.\n");

              size_t provided_elements = elements[sending_rank];
              Array buffer = array_new (provided_elements);

              MPI_Recv (buffer.data, provided_elements, MPI_NUMBER_T,
                        sending_rank, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

              stack_push_array (buffer_stack, buffer);
              sending_rank++;

              array_free (buffer);
            }

          /* Send the elements to the process. */
          Array buffer = array_new (requested_elements);

          stack_pop_array (buffer_stack, buffer, 0, requested_elements);
          if (i == 0)
            memcpy (receive_buffer.data, buffer.data,
                    sizeof (number_t) * requested_elements);
          else
            MPI_Send (buffer.data, requested_elements, MPI_NUMBER_T, i, 0,
                      comm);

          array_free (buffer);
        }

      /* The stack should be empty by this point. */
      assert (stack_get_size (buffer_stack) == 0);

      sending_rank = 0;
      for (int i = middle_rank; i < comm_size; i++)
        {
          /* How many elements does this process need? */
          size_t requested_elements = elements[i];

          while (stack_get_size (buffer_stack) < requested_elements)
            {
              if (sending_rank >= middle_rank)
                fprintf (stderr, "Imbalanced request going right.\n");

              size_t provided_elements = elements[sending_rank];
              Array buffer = array_new (provided_elements);

              if (sending_rank == 0)
                {
                  memcpy (buffer.data, stack->stack,
                          provided_elements * sizeof (number_t));
                  stack->top = 0; /* The data has been sent, empty the stack. */
                }
              else
                MPI_Recv (buffer.data, provided_elements, MPI_NUMBER_T,
                          sending_rank, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

              stack_push_array (buffer_stack, buffer);
              sending_rank++;

              array_free (buffer);
            }

          /* Send the elements to the process. */
          Array buffer = array_new (requested_elements);

          stack_pop_array (buffer_stack, buffer, 0, requested_elements);
          MPI_Send (buffer.data, requested_elements, MPI_NUMBER_T, i, 0, comm);

          array_free (buffer);
        }

      /* The stack should once again be empty. */
      assert (stack_get_size (buffer_stack) == 0);

      stack_free (buffer_stack);
    }
  else
    {
      MPI_Status status;

      /* Do a send/receive so that the master can coordinate the rest.
       */
      MPI_Sendrecv (stack->stack, elements_to_swap, MPI_NUMBER_T, 0, 0,
                    receive_buffer.data, elements_to_swap, MPI_NUMBER_T, 0,
                    MPI_ANY_TAG, comm, &status);

      int actual_elements_received = -1;
      MPI_Get_count (&status, MPI_NUMBER_T, &actual_elements_received);
      assert ((int) elements_to_swap == actual_elements_received);

      /* Empty the stack, since the data has been sent. */
      stack->top = 0;
    }
  stack_push_array (stack, receive_buffer);
  array_free (receive_buffer);

  /* All MPI communication should be finished by this point. */
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      if (point[1] == threshold)
        {
          stack_pop_buffer (stack, point, stride);
        }
    }
  assert (stack_get_size (stack) == 0);
  stack_free (stack);
}

static void
test_partition_by_value ()
{
  Dataset dataset = dataset_new (1, 1024);
  int rank, comm_size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &comm_size);

  Array values = array_new (dataset.size);

  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      point[0] = (number_t) (i + rank * dataset.size); /* Random ID. */
      /* Should produce an equal number of 0s and 1s. */
      point[1] = (number_t) (i % 2);
      values.data[i] = point[1];
    }

  partition_by_value (dataset, values, 0.5f, MPI_COMM_WORLD);

  bool is_left = rank < (comm_size / 2) ? true : false;

  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      if (is_left)
        {
          assert (point[1] <= 0.5f);
        }
      else
        {
          assert (point[1] >= 0.5f);
        }
    }

  dataset_free (dataset);
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
          test_partitioning ();
          test_log2i ();
          test_stack ();
          test_dataset ();
          test_vp_tree_local ();
          test_partition_by_value ();
          printf ("All tests sucessful.\n");
        }
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
