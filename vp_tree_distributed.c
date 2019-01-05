#include "vp_tree_distributed.h"
#include "mpi_partition.h"
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
group_size (int world_size, int l)
{
  return world_size / (1 << l);
}

/**
 * Should be the same for all processes within the same working group.
 * To be used with functions such as MPI_Comm_split().
 */
int
group_number (int world_rank, int world_size, int l)
{
  return world_rank / group_size (world_size, l);
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
  else
    {
      slavePart (group_rank, distances.size, distances.data,
                 distances.size * group_size, comm);
    }
  MPI_Bcast (&median_distance, 1, MPI_FLOAT, 0, comm);
  return median_distance;
}

MPI_Comm
communicator_for_level (int l)
{
  MPI_Comm c;
  int world_rank;
  int world_size;
  MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);
  MPI_Comm_split (MPI_COMM_WORLD, group_number (world_rank, world_size, l), 0,
                  &c);
  return c;
}

void
write_distances_from_vp (Array dest, Dataset dataset, number_t *vp)
{
  assert (dest.size >= dataset.size);
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t distance = point_distance (dataset_point (dataset, i), vp,
                                          dataset.feature_count);
      dest.data[i] = distance;
    }
}

VPTreeDistributed
vp_tree_dist_from_dataset (Dataset dataset)
{
  int world_size, world_rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);
  VPTreeDistributed tree;

  size_t tree_depth = (size_t) log2i (world_size);

  tree.vps = array_new (tree_depth * (dataset.feature_count + 1));
  tree.medians = array_new (tree_depth);
  tree.dataset = dataset;

  Array distances = array_new (dataset.size);
  for (size_t i = 0; i < tree_depth; i++)
    {
      /* Create a communicator. */
      MPI_Comm comm = communicator_for_level (i);

      int comm_size, comm_rank;
      MPI_Comm_size (comm, &comm_size);
      MPI_Comm_rank (comm, &comm_rank);

      /* Pick a vantage point. */
      number_t vp[dataset.feature_count + 1];
      if (comm_rank == 0)
        {
          size_t vp_index = 0; // (size_t) (rand () % (dataset.size);
          dataset_read (dataset, vp, vp_index, 1);
        }
      MPI_Bcast (vp, dataset.feature_count + 1, MPI_NUMBER_T, 0, comm);

      /* Find the distances within the local dataset. */
      write_distances_from_vp (distances, dataset, vp);

      /* Find the median distance. */
      Array values_to_shuffle = array_copy (distances);
      number_t median
        = get_median_distance (values_to_shuffle, comm, comm_rank, comm_size);
      array_free (values_to_shuffle);

      /* Partition the data around the median distance. */
      mpi_partition_by_value (dataset, distances, median, comm);

      /* Write the vp into the array. */
      vp_tree_dist_write_vp (tree, vp, i);

      /* Write the median into the array. */
      tree.medians.data[i] = median;
    }
  array_free (distances);

  tree.local_tree = vp_tree_from_dataset (tree.dataset);

  return tree;
}

void
vp_tree_dist_free (VPTreeDistributed tree)
{
  array_free (tree.vps);
  array_free (tree.medians);
  vp_tree_free (tree.local_tree);
}

void
vp_tree_dist_write_vp (VPTreeDistributed tree, number_t *vp, size_t index)
{
  memcpy (&tree.vps.data[index * (tree.dataset.feature_count + 1)], vp,
          sizeof (number_t) * (tree.dataset.feature_count + 1));
}

void
vp_tree_dist_read_vp (VPTreeDistributed tree, number_t *vp, size_t index)
{
  memcpy (vp, &tree.vps.data[index * (tree.dataset.feature_count + 1)],
          sizeof (number_t) * (tree.dataset.feature_count + 1));
}

Dataset
vp_tree_dist_find_knn (VPTreeDistributed tree, number_t *target, size_t k)
{
  int world_size, world_rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);
  MPI_Bcast (target, tree.dataset.feature_count + 1, MPI_NUMBER_T, 0,
             MPI_COMM_WORLD);
  MPI_Bcast (&k, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  /* Find the local nearest neighbors. */
  Dataset knn = vp_tree_find_knn (tree.local_tree, target, k);

  size_t tree_depth = tree.medians.size;

  for (size_t i = 0; i < tree_depth; i++)
    {
      /* Create a communicator. */
      MPI_Comm comm = communicator_for_level (tree_depth - i - 1);

      int comm_size, comm_rank;
      MPI_Comm_size (comm, &comm_size);
      MPI_Comm_rank (comm, &comm_rank);
      int middle_rank = comm_size / 2;

      /* If we are the root, receive the right leaf and merge the results. */
      if (comm_rank == 0)
        {
          Dataset incoming = dataset_new (knn.feature_count, knn.size);
          MPI_Recv (incoming.data.data, k * (knn.feature_count + 1),
                    MPI_NUMBER_T, middle_rank, 0, comm, MPI_STATUS_IGNORE);
          for (size_t j = 0; j < incoming.size; j++)
            {
              number_t *point = dataset_point (incoming, j);
              queue_insert (knn, point, target);
            }
          dataset_free (incoming);
        }
      else if (comm_rank == comm_size / 2)
        {
          /* Send our so-far findings to the master. */
          MPI_Send (knn.data.data, k * (knn.feature_count + 1), MPI_NUMBER_T, 0,
                    0, comm);
        }
      else
        continue;
      MPI_Comm_free (&comm);
    }

  /* Only the world root will return the correct result. */
  return knn;
}

void
verify_vp_tree_dist (VPTreeDistributed tree)
{
  size_t tree_depth = tree.medians.size;
  Array distances = array_new (tree.dataset.size);
  for (size_t i = 0; i < tree_depth; i++)
    {
      number_t vp[tree.dataset.feature_count + 1];
      vp_tree_dist_read_vp (tree, vp, i);
      number_t median = tree.medians.data[i];

      MPI_Comm comm = communicator_for_level (i);

      int comm_size, comm_rank;
      MPI_Comm_size (comm, &comm_size);
      MPI_Comm_rank (comm, &comm_rank);

      bool is_left = comm_rank < (comm_size / 2);

      write_distances_from_vp (distances, tree.dataset, vp);

      for (size_t i = 0; i < tree.dataset.size; i++)
        {
          if (is_left)
            assert (distances.data[i] <= median);
          else
            assert (distances.data[i] >= median);
        }

      MPI_Comm_free (&comm);
    }
  array_free (distances);
  vp_tree_local_verify (tree.local_tree, 0, tree.dataset.size, 0);
}

void
test_vp_tree_distributed ()
{
  Dataset dataset = dataset_new (2, 1024);
  dataset_fill_random (dataset);

  VPTreeDistributed tree = vp_tree_dist_from_dataset (dataset);
  verify_vp_tree_dist (tree);

  number_t target[3] = {0.0f, 0.0f, 0.0f};
  Dataset nearest8 = vp_tree_dist_find_knn (tree, target, 8);

  number_t *last_neighbor = dataset_point (nearest8, nearest8.size - 1);
  MPI_Bcast (last_neighbor, dataset.feature_count + 1, MPI_NUMBER_T, 0,
             MPI_COMM_WORLD);
  number_t furthest_distance
    = point_distance (last_neighbor, target, dataset.feature_count);

  /* Count the number of elements further away from the target than the last
   * neighbor, and make sure it equals dataset size - k.
   */
  int local_count = 0;
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      number_t distance = point_distance (point, target, dataset.feature_count);
      if (distance > furthest_distance)
        local_count++;
    }

  int total_count = 0;
  MPI_Reduce (&local_count, &total_count, 1, MPI_NUMBER_T, MPI_SUM, 0,
              MPI_COMM_WORLD);

  int rank, world_size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);

  if (rank == 0)
    {
      assert ((size_t) total_count
              == dataset.size * world_size - nearest8.size);
    }

  dataset_free (nearest8);
  dataset_free (dataset);
}