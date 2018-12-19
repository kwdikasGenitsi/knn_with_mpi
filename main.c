#include "median.h"
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef float number_t;

typedef struct
{
  number_t *data;
  int size;
} Array;

int world_size;
int world_rank;

typedef struct
{
  number_t vantage_point;
  number_t median;
} VPNode;

VPNode *vp_stack;
int vp_stack_length;
int vp_stack_top;

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

void
array_free (Array *array)
{
  free (array->data);
  free (array);
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
split_parallel (int l)
{
  MPI_Comm comm = communicator_for_level (l);

  /* Pick a vantage point and announce it. */
  /* Find the median. */
  /* Split the free from the commies. */
  /* Push the vantage point and the median onto the stack. */

  MPI_Comm_free (&comm);
}

int
main (int argc, char **argv)
{
  test_partitioning ();
  test_log2i ();

  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);

  Array *dataset = array_new_random (5);
  int dataset_size = dataset->size * world_size;

  vp_stack_length = log2i (dataset_size);
  vp_stack_top = 0;
  vp_stack = malloc (sizeof (*vp_stack));

  number_t vantage_point = 5;

  Array *distances = array_distances_from_vp (dataset, vantage_point);

  for (int l = 0; group_size (l) > 1; l++)
    {
      split_parallel (l);
    }

  free (vp_stack);
  array_free (dataset);
  MPI_Finalize ();
  return 0;
}
