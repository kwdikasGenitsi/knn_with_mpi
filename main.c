#include "median.h"
#include "stack.h"
#include "vp_master_buffer.h"
#include "vp_tree_local.h"
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

number_t
eucledian_distance (Point point1, Point point2)
{
  number_t distance, squares_sum = 0;
  assert (point1.feature_count == point2.feature_count);
  for (size_t i = 0; i < point1.feature_count; i++)
    {
      squares_sum += (point1.data.data[i + 1] - point2.data.data[i + 1])
                     * (point1.data.data[i + 1] - point2.data.data[i + 1]);
    }
  // distance = sqrt (squares_sum);
  distance = squares_sum;
  return distance;
}

void
test_point_dataset ()
{
  Point points[10], points_from_dataset[10];
  Dataset dataset = dataset_new (2, 10);
  for (int i = 0; i < 10; i++)
    {
      points[i] = point_new (2);
      point_fill_random (points[i]);
      enter_point_to_dataset (&dataset, points[i], i);
      points_from_dataset[i] = get_point_from_dataset (&dataset, i);
      assert (eucledian_distance (points_from_dataset[i], points[i]) == 0);
      point_free (points[i]);
      point_free (points_from_dataset[i]);
    }
  dataset_free (dataset);
}

Array
array_distances_from_vp (Dataset *dataset, Point vantage_point)
{
  size_t number_of_points = dataset->size;
  Array distances = array_new (number_of_points);
  for (int i = 0; i < number_of_points; i++)
    {
      distances.data[i]
        = eucledian_distance (get_point_from_dataset (dataset, i),
                              vantage_point);
    }
  return distances;
}

int
less_than_median (Array *distances, number_t median_distance)
{
  int count_points
    = 0; // counts how many points are closer to vantage point than median.
  for (int i = 0; i < distances->size; i++)
    {
      if (distances->data[i] <= median_distance)
        {
          count_points++;
        }
    }
  return count_points;
}

Array
data_to_send (Dataset *dataset, Array *distances, number_t median_distance,
              const int is_process_left, size_t less_than_median)
{
  Dataset dataset_to_send;
  size_t number_of_points = dataset->size;
  assert (number_of_points == distances->size);
  if (is_process_left)
    {
      dataset_to_send = dataset_new (dataset->feature_count,
                                     distances->size - less_than_median);
    }
  else
    {
      dataset_to_send = dataset_new (dataset->feature_count, less_than_median);
    }

  int k = 0;
  for (size_t point_index = 0; point_index < number_of_points; point_index++)
    {
      if (is_process_left)
        {
          if (distances->data[point_index] > median_distance)
            {
              Point point_to_enter
                = get_point_from_dataset (dataset, point_index);
              enter_point_to_dataset (&dataset_to_send, point_to_enter, k);
              point_free (point_to_enter);
              k++;
            }
        }
      else
        {
          if (distances->data[point_index] <= median_distance)
            {
              Point point_to_enter
                = get_point_from_dataset (dataset, point_index);
              enter_point_to_dataset (&dataset_to_send, point_to_enter, k);
              point_free (point_to_enter);
              k++;
            }
        }
    }
  return dataset_to_send.data;
}

void
test_exchanged_data_validity (Point vantage_point, Array array,
                              number_t median_distance,
                              size_t should_be_less_than_median,
                              size_t feature_count)
{
  Dataset dataset = array_to_dataset (array, feature_count);
  Array distances = array_distances_from_vp (&dataset, vantage_point);
  if (should_be_less_than_median)
    {
      // All distances from vp should be less than median.
      assert (dataset.size == less_than_median (&distances, median_distance));
    }
  else
    {
      // All distances from vp should be more than median.
      assert (less_than_median (&distances, median_distance) == 0);
    }
}
void
request_master_recieve (int *sending_process, Stack *master_buffer,
                        int *less_than_median_array, size_t feature_count,
                        size_t dataset_points_count, int is_sender_left,
                        int group_comm, int group_size)
{
  if (is_sender_left)
    {
      if (*sending_process >= group_size / 2)
        {
          return;
        }
      size_t sending_points_count
        = dataset_points_count - less_than_median_array[*sending_process];
      if (sending_points_count
          > point_count (feature_count, stack_get_free_space (master_buffer)))
        {
          return;
        }
      else
        {
          Array new_points
            = array_new (dataset_size (feature_count, sending_points_count));
          MPI_Recv (new_points.data, new_points.size, MPI_FLOAT,
                    *sending_process, 0, group_comm, MPI_STATUS_IGNORE);
          stack_push_array (master_buffer, &new_points);
          array_free (new_points);
          (*sending_process)++;
        }
    }
  else
    {
      if (*sending_process >= group_size)
        {
          return;
        }
      size_t sending_points_count = less_than_median_array[*sending_process];
      if (sending_points_count
          > point_count (feature_count, stack_get_free_space (master_buffer)))
        {
          return;
        }
      else
        {
          Array new_points
            = array_new (dataset_size (feature_count, sending_points_count));
          printf ("before master recv! with sending process %d\n",
                  *sending_process);
          MPI_Recv (new_points.data, new_points.size, MPI_FLOAT,
                    *sending_process, 0, group_comm, MPI_STATUS_IGNORE);
          printf ("after master recv!\n");
          stack_push_array (master_buffer, &new_points);
          array_free (new_points);
          (*sending_process)++;
        }
    }
}

void
request_master_send (int *recieving_process, Stack *master_buffer,
                     int *less_than_median_array, size_t feature_count,
                     size_t dataset_points_count, int is_reciever_left,
                     int group_comm, int group_size)
{
  if (is_reciever_left)
    {
      if (*recieving_process >= group_size / 2)
        {
          return;
        }
      size_t recieving_points_count
        = dataset_points_count - less_than_median_array[*recieving_process];
      if (recieving_points_count
          > point_count (feature_count, stack_get_size (master_buffer)))
        {
          return;
        }
      else
        {
          Array old_points
            = array_new (dataset_size (feature_count, recieving_points_count));
          stack_pop_array (master_buffer, &old_points, 0, old_points.size);
          printf ("before master Send with recieving_process: %d\n",
                  *recieving_process);
          MPI_Send (old_points.data, old_points.size, MPI_FLOAT,
                    *recieving_process, 0, group_comm);
          array_free (old_points);
          (*recieving_process)++;
        }
    }
  else
    {
      if (*recieving_process >= group_size)
        {
          return;
        }
      size_t recieving_points_count
        = less_than_median_array[*recieving_process];
      if (recieving_points_count
          > point_count (feature_count, stack_get_size (master_buffer)))
        {
          return;
        }
      else
        {
          Array old_points
            = array_new (dataset_size (feature_count, recieving_points_count));
          stack_pop_array (master_buffer, &old_points, 0, old_points.size);
          MPI_Send (old_points.data, old_points.size, MPI_FLOAT,
                    *recieving_process, 0, group_comm);
          array_free (old_points);
          (*recieving_process)++;
        }
    }
}

void
merge_recieved_points (Dataset *dataset, Array *distances,
                       Dataset *recieved_points, number_t median_distance,
                       size_t is_process_left)
{

  size_t number_of_points = dataset->size;
  assert (number_of_points == distances->size);
  int k = 0;
  for (size_t point_index = 0; point_index < number_of_points; point_index++)
    {
      if (distances->data[point_index] == median_distance)
        {
          printf ("a distance is equal to median!\n");
          return;
        }

      if (is_process_left)
        {

          if (distances->data[point_index] == median_distance)
            {
              printf ("its EQUAL to median from process: %d\n", world_rank);
            }
          if (distances->data[point_index] > median_distance)
            {
              Point point_to_merge
                = get_point_from_dataset (recieved_points, k++);
              enter_point_to_dataset (dataset, point_to_merge, point_index);
              point_free (point_to_merge);
            }
        }
      else
        {
          if (distances->data[point_index] < median_distance)
            {
              Point point_to_merge
                = get_point_from_dataset (recieved_points, k++);
              enter_point_to_dataset (dataset, point_to_merge, point_index);
              point_free (point_to_merge);
            }
        }
    }
}

void
exchange_group_data (Dataset *dataset, Array *distances,
                     number_t median_distance, Point vantage_point,
                     MPI_Comm group_comm, int group_rank, int group_size)
{
  // -- ------split it on another function later-------------------------------
  int less_than_median_int = less_than_median (distances, median_distance);
  int *less_than_median_array = NULL;
  Stack *master_buffer = NULL;
  if (group_rank == 0)
    {
      less_than_median_array = (int *) malloc (sizeof (int) * group_size);
      master_buffer = stack_new (2 * dataset->data.size);
    }
  MPI_Gather (&less_than_median_int, 1, MPI_INT, less_than_median_array, 1,
              MPI_INT, 0, group_comm);
  // -------------------------------------------------------------------------
  int is_process_left = 1 - ((group_size / 2) + group_rank) / group_size;
  int should_be_less_than_median;
  if (group_rank < group_size / 2)
    {
      is_process_left = 1;
      should_be_less_than_median = 0;
    }
  else
    {
      is_process_left = 0;
      should_be_less_than_median = 1;
    }
  Array points_to_send = data_to_send (dataset, distances, median_distance,
                                       is_process_left, less_than_median_int);
  test_exchanged_data_validity (vantage_point, points_to_send, median_distance,
                                should_be_less_than_median,
                                dataset->feature_count);
  Array recieved_points = array_new (points_to_send.size);
  if (group_rank == 0)
    {
      int recieving_process = 1, sending_process = group_size / 2;
      while (recieving_process < group_size / 2 || sending_process < group_size)
        {
          request_master_recieve (&sending_process, master_buffer,
                                  less_than_median_array,
                                  dataset->feature_count, dataset->size, 0,
                                  group_comm, group_size);
          // printf ("master recieved! sending_process: %d\n",
          // (sending_process - 1));
          request_master_send (&recieving_process, master_buffer,
                               less_than_median_array, dataset->feature_count,
                               dataset->size, 1, group_comm, group_size);
          // printf ("master send! recieving_process: %d\n",
          //   (recieving_process - 1));
        }
    }
  else if (group_rank >= group_size / 2)
    {
      // test_exchanged_data_validity (vantage_point, points_to_send,
      //                            median_distance, 1, dataset->feature_count);
      MPI_Send (points_to_send.data, points_to_send.size, MPI_FLOAT, 0, 0,
                group_comm);
    }
  else if (group_rank < group_size / 2)
    {
      MPI_Recv (recieved_points.data, recieved_points.size, MPI_FLOAT, 0, 0,
                group_comm, MPI_STATUS_IGNORE);
    }
  else
    {
      printf ("its imposible to come here!\n");
    }
  MPI_Barrier (group_comm);
  Dataset recieved_points_dataset
    = array_to_dataset (recieved_points, dataset->feature_count);
  if (group_rank < group_size / 2)
    {
      merge_recieved_points (dataset, distances, &recieved_points_dataset,
                             median_distance, 1);
    }
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

Point
get_vantage_point (size_t feature_count, MPI_Comm comm, int group_rank)
{
  Point vantage_point = point_new (feature_count);
  if (group_rank == 0)
    {
      point_fill_random (vantage_point);
    }
  MPI_Bcast (vantage_point.data.data, vantage_point.data.size, MPI_FLOAT, 0,
             comm);
  return vantage_point;
}
void
test_exchange_group_data (MPI_Comm comm, int group_rank, int group_size)
{
  Dataset dataset = dataset_new (1, 5);
  dataset_fill_random (dataset);
  Point vantage_point = get_vantage_point (1, comm, group_rank);
  if (group_rank == 0)
    {
      printf ("vantage_point = \n");
      print_point (vantage_point);
    }
  // printf ("old dataset rank %d\n", group_rank);
  // print_dataset (&dataset);
  Array distances = array_distances_from_vp (&dataset, vantage_point);

  number_t median_distance
    = get_median_distance (distances, comm, group_rank, group_size);
  // printf ("hi im before exchange!\n");
  MPI_Barrier (comm);

  exchange_group_data (&dataset, &distances, median_distance, vantage_point,
                       comm, group_rank, group_size);
  Array new_distances = array_distances_from_vp (&dataset, vantage_point);
  // printf ("new dataset rank %d\n", group_rank);
  // print_dataset (&dataset);
  if (group_rank < group_size / 2)
    {
      for (size_t i = 0; i < new_distances.size; i++)
        {
          // printf ("im in!\n");
          if (new_distances.data[i] > median_distance)
            {
              printf ("ERROR!! LARGER THAN MEDIAN!! with rank %d\n",
                      group_rank);
              return;
            }
        }
    }
  printf ("from process: %d, test_exchange_group_data: PASSED!\n", group_rank);
  return;
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
