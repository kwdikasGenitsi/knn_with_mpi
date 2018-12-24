#include "median.h"
#include "vp_master_buffer.h"
#include "vp_stack.h"
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int world_size;
int world_rank;

/**
 * To know the depth of our tree, we need to calculate log2(x) accurately.
 * We can't use floating point numbers since we risk getting an off-by-one error
 * due to floating point accuracy.
 * We therefore use this GCC bultin, which returns the number of leading 0s in
 * an integer. On x86 this is done with the BSR instruction.
 */
static inline int
log2i (int x)
{
  assert (x > 0);
  return sizeof (int) * 8 - __builtin_clz (x) - 1;
}

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

Array *
array_new_random (int size)
{
  Array *array = (Array *) malloc (sizeof (*array));
  array->size = size;
  array->data = (number_t *) malloc (sizeof (number_t) * size);
  int cal = 5;
  srand ((cal + 1) * time (NULL));

  for (int i = 0; i < size; i++)
    array->data[i]
      = world_rank * size + i; // (number_t)(rand() - rand()) * 0.05f;
  return array;
}

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

Array *
array_distances_from_vp (Array *dataset, number_t vantage_point)
{
  Array *distances = (Array *) malloc (sizeof (Array));
  distances->size = dataset->size;
  distances->data = (number_t *) malloc (sizeof (number_t) * distances->size);
  for (int i = 0; i < dataset->size; i++)
    {
      distances->data[i] = fabs (
        vantage_point - dataset->data[i]); // we calculate the distance of the
                                           // ith point from the vantage point
    }
  return distances;
}

void
test_array_distances_from_vp ()
{
  Array *dataset = array_new_random (5);
  number_t vp = 4.56f;
  Array *distances = array_distances_from_vp (dataset, vp);
  for (int i = 0; i < dataset->size; i++)
    {
      printf (
        "PROCESS %d: vantage_point = %f, point[%d] = %f and distance = %f\n",
        world_rank, vp, i, dataset->data[i], distances->data[i]);
    }
  array_free (dataset);
  /** @todo Do some actual testing with assert() */
}

int
less_than_median (Array *distances, number_t vantage_point, number_t median)
{
  int count_points
    = 0; // counts how many points are closer to vantage point than median.
  for (int i = 0; i < distances->size; i++)
    {
      if (distances->data[i] <= median)
        {
          count_points++;
        }
    }
  return count_points;
}

Array *
points_for_transfer (Array *dataset, Array *distances, number_t vantage_point,
                     number_t median, const int is_process_left,
                     int less_than_median)
{
  Array *transfer_buffer = (Array *) malloc (sizeof (*transfer_buffer));
  if (is_process_left)
    {
      transfer_buffer->data = (number_t *) malloc (
        sizeof (number_t *) * (dataset->size - less_than_median));
      transfer_buffer->size = dataset->size - less_than_median;
    }
  else
    {
      transfer_buffer->data
        = (number_t *) malloc (sizeof (number_t *) * less_than_median);
      transfer_buffer->size = less_than_median;
    }

  int k = 0;
  for (int i = 0; i < dataset->size; i++)
    {
      if (is_process_left)
        {
          if (distances->data[i] > median)
            {
              transfer_buffer->data[k] = dataset->data[i];
              k++;
            }
        }
      else
        {
          if (distances->data[i] <= median)
            {
              transfer_buffer->data[k] = dataset->data[i];
              k++;
            }
        }
    }
  return transfer_buffer;
}

void
test_master_buffer_functionality ()
{
  Array *input_array = array_new (3);
  MasterBuffer *master_buffer = master_buffer_new (9);
  int k = 0, l = 0;
  for (int i = 0; i < 15; i++)
    {
      input_array->data[k] = i;
      k++;
      if (k == 3)
        {
          l++;
          k = 0;
          master_buffer_fill (master_buffer, input_array);
        }
      if (l == 2)
        {
          l = 0;
          array_free (input_array);
          input_array = master_buffer_throw (master_buffer, 3);
        }
    }
  int test_array[9] = {0, 1, 2, 6, 7, 8, 12, 13, 14};
  for (int i = 0; i < 9; i++)
    {
      assert (test_array[i] == (int) master_buffer->buffer->data[i]);
    }
  array_free (input_array);
  master_buffer_free (master_buffer);
}

void
test_points_for_transfer (Array *dataset, Array *distances, number_t median)
{
  int vantage_point = 4.04f;
  Array *transfer_buffer
    = points_for_transfer (dataset, distances, vantage_point, median, 0,
                           less_than_median (distances, vantage_point, median));
  if (world_rank == 0)
    {
      printf ("transfer_buffer size is %d\n", transfer_buffer->size);
      printf ("median distance is %f\n", median);
      for (int i = 0; i < distances->data[i]; i++)
        {
          printf ("distances[%d] = %f\n", i, distances->data[i]);
        }
      for (int i = 0; i < transfer_buffer->size; i++)
        {
          printf ("transfer_buffer[%d] = %f\n", i, transfer_buffer->data[i]);
        }
    }
}

void
transfer_points_subteam (Array *dataset, Array *distances,
                         number_t vantage_point, number_t median,
                         MPI_Comm subteam_comm, MPI_Comm subteam_rank,
                         MPI_Comm subteam_size)
{
  // -- ------split it on another function later-------------------------------
  int less_than_median_int
    = less_than_median (distances, vantage_point, median);
  int *less_than_median_array = NULL;
  MasterBuffer *master_buffer = NULL;
  if (subteam_rank == 0)
    {
      less_than_median_array = (int *) malloc (sizeof (int) * subteam_size);
      master_buffer = master_buffer_new (2 * subteam_size);
    }
  MPI_Gather (&less_than_median_int, 1, MPI_INT, less_than_median_array, 1,
              MPI_INT, 0, subteam_comm);
  // -------------------------------------------------------------------------
  int is_process_left = 1 - ((subteam_size / 2) + subteam_rank) / subteam_size;
  Array *points_to_send
    = points_for_transfer (dataset, distances, vantage_point, median,
                           is_process_left, less_than_median_int);
  Array *recieved_points = array_new (points_to_send->size);
  if (subteam_rank == 0)
    {
      int recieving_process = 0, sending_process = subteam_size / 2;
      while (recieving_process < subteam_size / 2
             || sending_process < subteam_size)
        {
          request_master_recieve (&sending_process, master_buffer,
                                  less_than_median_array, subteam_comm,
                                  subteam_size);
        }
    }
  else if (subteam_rank >= subteam_size / 2)
    {
      MPI_Send (points_to_send->data, points_to_send->size, MPI_FLOAT, 0, 0,
                subteam_comm);
    }
}

void
request_master_recieve (int *sending_process, MasterBuffer *master_buffer,
                        int *less_than_median_array, MPI_Comm subteam_comm,
                        MPI_Comm subteam_size)
{
  // if (sending_process >= subteam_size
  //  || less_than_median_array[sending_process]) // !! carefull with sending
  // size!
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
split_parallel (Array *dataset, VpStack *vp_stack, int l)
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
  vp_stack_push (vp_stack,
                 (VpNode){.vantage_point = vp_value, .median = median});

  MPI_Comm_free (&comm);
}

/**
 * Generates one level of the VP tree serially.
 * @param dataset The global dataset.
 * @param offset The start of the chunk to be processed.
 * @param l The depth of the node (determines the chunk size).
 */
void
split_serial (Array *dataset, int offset, int l)
{
  /** @todo */
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
          test_partitioning ();
          test_log2i ();
          test_vp_stack ();
          test_master_buffer_functionality ();
          printf ("All tests sucessful.\n");
        }
      return EXIT_SUCCESS;
    }

  Array *dataset = array_new_random (5);
  VpStack *vp_stack = vp_stack_new (log2i (dataset->size * world_size));

  for (int l = 0; group_size (l) > 1; l++)
    {
      split_parallel (dataset, vp_stack, l);
    }

  array_free (dataset);
  MPI_Finalize ();
  return EXIT_SUCCESS;
}
