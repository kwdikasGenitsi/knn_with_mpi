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
  for (size_t i = 0; i < 1; i++)
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

      if (comm_rank == 0)
        {
          printf ("Before partitioning: ");
          array_dump (distances);
        }

      /* Find the median distance. */
      number_t median
        = get_median_distance (distances, comm, comm_rank, comm_size)
          + 0.0f; // 0.000005f;

      /* Partition the data around the median distance. */
      mpi_partition_by_value (dataset, distances, median, comm);

      if (comm_rank == 0)
        {
          printf ("Median: %f - ", median);
          write_distances_from_vp (distances, dataset, vp);
          array_dump (distances);
        }

      /* Write the vp into the array. */
      vp_tree_dist_write_vp (tree, vp, i);

      /* Write the median into the array. */
      tree.medians.data[i] = median;
    }
  array_free (distances);

  return tree;
}

void
vp_tree_dist_free (VPTreeDistributed tree)
{
  array_free (tree.vps);
  array_free (tree.medians);
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
}

void
test_vp_tree_distributed ()
{
  Dataset dataset = dataset_new (2, 4);
  dataset_fill_random (dataset);

  VPTreeDistributed tree = vp_tree_dist_from_dataset (dataset);
  verify_vp_tree_dist (tree);
  UNUSED (tree);

  dataset_free (dataset);
}